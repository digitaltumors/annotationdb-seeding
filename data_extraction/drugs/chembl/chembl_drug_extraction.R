#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
    library(httr)
    library(future.apply)
    library(AnnotationGx)
})

args <- commandArgs(trailingOnly = TRUE)
prop_out <- "../pubchem/output_data/union/mar-27-2026/union_out.csv"
out_dir <- "../pubchem/output_data/union/mar-27-2026"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Logging --------------------
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logf <- function(fmt, ...) cat(sprintf("[%s] %s\n", timestamp(), sprintf(fmt, ...)))
log_step <- function(label) logf("== %s ==", label)

# -------------------- Config --------------------
PAR <- 0
if (is.na(PAR) || PAR < 0) PAR <- 0
CHUNK <- 200
LIMIT <- 5000
ua <- httr::user_agent("pubchem-chembl-script/1.0")

# Flatten list-columns returned by getChemblMechanism (variant_sequence, mechanism_refs)
flatten_list_cols <- function(dt) {
    list_cols <- c("variant_sequence", "mechanism_refs")
    for (col in intersect(list_cols, names(dt))) {
        set(dt, j = col, value = vapply(dt[[col]], function(x) {
            if (is.null(x) || (is.logical(x) && all(is.na(x)))) {
                return(NA_character_)
            }
            tryCatch(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"),
                error = function(e) NA_character_
            )
        }, character(1)))
    }
    dt
}


if (PAR > 0) {
    logf("Enabling parallel fetching with %d workers (using mclapply for real-time logs)", PAR)
} else {
    logf("Running sequentially (PAR=0)")
}

# -------------------- ID collection --------------------
extract_chembl_ids_from_csv <- function(path) {
    if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
        return(character(0))
    }
    ids <- character(0)
    dt <- tryCatch(fread(path, fill = TRUE, quote = "\"'", showProgress = FALSE), error = function(e) NULL)
    if (!is.null(dt)) {
        vals <- unlist(lapply(dt, function(col) as.character(col)), use.names = FALSE)
        vals <- vals[!is.na(vals) & nzchar(vals)]
        ids <- unique(unlist(regmatches(vals, gregexpr("CHEMBL\\d+", vals, perl = TRUE))))
    }
    if (length(ids)) {
        return(ids)
    }
    txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
    if (!length(txt)) {
        return(character(0))
    }
    unique(unlist(regmatches(txt, gregexpr("CHEMBL\\d+", txt, perl = TRUE))))
}

log_step("Collecting ChEMBL IDs")
chembl_ids <- extract_chembl_ids_from_csv(prop_out)

completed_path <- file.path(out_dir, "chembl_completed_ids.csv")
if (file.exists(completed_path)) {
    done_ids <- as.character(fread(completed_path, header = FALSE)$V1)
    before <- length(chembl_ids)
    chembl_ids <- setdiff(chembl_ids, done_ids)
    logf("Found %d already completed IDs. Reduced from %d total to %d remaining to query.", length(done_ids), before, length(chembl_ids))
}

if (!length(chembl_ids)) {
    logf("No valid ChEMBL IDs found in %s; nothing to do.", prop_out)
    quit(save = "no", status = 0)
}
logf(
    "Collected %d IDs (e.g., %s%s)",
    length(chembl_ids),
    paste(utils::head(chembl_ids, 3), collapse = ", "),
    if (length(chembl_ids) > 3) ", ..." else ""
)

# -------------------- HTTP helpers --------------------
get_json <- function(url, tries = 5, timeout_sec = 60) {
    for (k in seq_len(tries)) {
        resp <- tryCatch(httr::GET(url, httr::timeout(timeout_sec), ua), error = function(e) NULL)
        if (!is.null(resp)) {
            sc <- httr::status_code(resp)
            if (sc == 200) {
                txt <- httr::content(resp, as = "text", encoding = "UTF-8")
                return(tryCatch(jsonlite::fromJSON(txt, flatten = TRUE), error = function(e) NULL))
            } else {
                logf("HTTP %s for %s", sc, url)
            }
        }
        Sys.sleep(10)
    }
    NULL
}

chunk_vec <- function(x, n) split(x, ceiling(seq_along(x) / n))

flatten_lists <- function(df) {
    if (!is.data.frame(df) || !NROW(df)) {
        return(df)
    }
    for (nm in names(df)) {
        if (is.list(df[[nm]])) {
            df[[nm]] <- vapply(
                df[[nm]],
                function(z) {
                    if (is.null(z)) {
                        return(NA_character_)
                    }
                    if (is.atomic(z) && length(z) == 1 && is.null(dim(z))) {
                        return(as.character(z))
                    }
                    jsonlite::toJSON(z, auto_unbox = TRUE, null = "null")
                },
                character(1),
                USE.NAMES = FALSE
            )
        }
    }
    df
}

# Helper: remove empty rows with no meaningful mechanism data
clean_mechanism_rows <- function(df) {
    if (!nrow(df)) {
        return(df)
    }

    # Ensure important columns exist and normalize to trimmed character
    ensure_col <- function(nm) {
        if (!(nm %in% names(df))) df[[nm]] <<- NA_character_
        v <- as.character(df[[nm]])
        v[is.na(v)] <- ""
        trimws(v)
    }

    moa <- ensure_col("mechanism_of_action")
    act <- ensure_col("action_type")
    tgt <- ensure_col("target_chembl_id")
    tnm <- ensure_col("target_name")

    # Keep row if at least one mechanism-descriptive field is non-empty
    keep <- nzchar(moa) | nzchar(act) | nzchar(tgt) | nzchar(tnm)

    # Also drop rows where ALL columns except molecule_chembl_id are empty/NA
    # (extra guard for strange payloads)
    non_mol_cols <- setdiff(names(df), "molecule_chembl_id")
    # Build a per-row "has any value" check across descriptive columns we normalized
    keep2 <- keep

    # Apply filters
    out <- df[keep2, , drop = FALSE]

    # If column types got promoted, return as data.table consistently
    data.table::as.data.table(out)
}


fetch_page <- function(endpoint, param_name, id_vec, limit = LIMIT, offset = 0) {
    base <- "https://www.ebi.ac.uk/chembl/api/data"
    q <- paste(id_vec, collapse = ",")
    url <- sprintf(
        "%s/%s?%s__in=%s&format=json&limit=%d&offset=%d",
        base, endpoint, param_name, utils::URLencode(q, reserved = TRUE),
        limit, offset
    )
    get_json(url)
}

fetch_chunk <- function(endpoint, param_name, id_vec, limit = LIMIT, chunk_idx = NA, n_chunks = NA) {
    key <- switch(endpoint,
        "mechanism"       = "mechanisms",
        "drug"            = "drugs",
        "molecule"        = "molecules",
        "molecule_form"   = "molecule_forms",
        stop(sprintf("Unknown endpoint: %s", endpoint))
    )

    rows <- list()
    offset <- 0L
    page <- 1L
    repeat {
        x <- fetch_page(endpoint, param_name, id_vec, limit = limit, offset = offset)
        if (is.null(x)) {
            logf(
                "[%s] chunk %s/%s page %s: NULL response, stopping",
                endpoint,
                ifelse(is.na(chunk_idx), "?", chunk_idx),
                ifelse(is.na(n_chunks), "?", n_chunks),
                page
            )
            break
        }
        df <- x[[key]]
        total <- suppressWarnings(as.integer(x$page_meta$total_count))
        got <- if (!is.null(df)) nrow(df) else 0

        logf(
            "[%s] chunk %s/%s page %s: got %d rows (offset=%d%s)",
            endpoint,
            ifelse(is.na(chunk_idx), "?", chunk_idx),
            ifelse(is.na(n_chunks), "?", n_chunks),
            page, got, offset,
            if (!is.na(total)) sprintf(", total=%d", total) else ""
        )

        if (is.null(df) || !NCOL(df) || !NROW(df)) break

        if (!("molecule_chembl_id" %in% names(df))) df$molecule_chembl_id <- NA_character_
        df <- flatten_lists(df)

        # Clean empty mechanisms immediately
        if (endpoint == "mechanism") {
            before <- nrow(df)
            df <- clean_mechanism_rows(df)
            after <- nrow(df)
            if (before > after) logf("[%s] filtered out %d empty rows", endpoint, before - after)
        }

        rows[[length(rows) + 1]] <- as.data.table(df)

        if (!is.na(total) && (offset + got) >= total) break
        if (got < limit) break
        offset <- offset + limit
        page <- page + 1L
        Sys.sleep(1.0)
    }

    if (!length(rows)) {
        return(data.table())
    }
    rbindlist(rows, fill = TRUE, use.names = TRUE)
}

fetch_union_chunks <- function(endpoint, ids, fields = c("molecule_chembl_id"),
                               chunk_size = CHUNK, limit = LIMIT, parallel_workers = PAR) {
    chunks <- chunk_vec(ids, chunk_size)
    n_chunks <- length(chunks)
    logf(
        "[%s] %d IDs → %d chunks (size %d), fields: %s",
        endpoint, length(ids), n_chunks, chunk_size, paste(fields, collapse = ", ")
    )

    fetch_field_one <- function(field) {
        logf("[%s] starting field '%s'", endpoint, field)
        if (parallel_workers > 0) {
            parallel::mclapply(seq_along(chunks), function(i) {
                ch <- chunks[[i]]
                fetch_chunk(endpoint, field, ch, limit = limit, chunk_idx = i, n_chunks = n_chunks)
            }, mc.cores = parallel_workers)
        } else {
            lapply(seq_along(chunks), function(i) {
                ch <- chunks[[i]]
                res <- fetch_chunk(endpoint, field, ch, limit = limit, chunk_idx = i, n_chunks = n_chunks)
                Sys.sleep(1.5)
                res
            })
        }
    }

    all_parts <- list()
    for (f in fields) {
        parts <- fetch_field_one(f)
        all_parts <- c(all_parts, parts)
        Sys.sleep(0.1)
    }
    if (!length(all_parts)) {
        return(data.table())
    }
    out <- rbindlist(all_parts, fill = TRUE, use.names = TRUE)
    if (endpoint == "mechanism") out <- clean_mechanism_rows(out)
    logf("[%s] combined rows: %d", endpoint, nrow(out))
    out
}

# -------------------- Outer Loop for Incremental Writes --------------------
OUTER_BATCH <- 1000 # Number of IDs to process before saving to disk
chunks_outer <- split(chembl_ids, ceiling(seq_along(chembl_ids) / OUTER_BATCH))
n_outer <- length(chunks_outer)

write_or_rm <- function(dt, path) {
    if (nrow(dt)) {
        is_append <- file.exists(path)
        fwrite(dt, path, append = is_append)
        logf("%s %d rows → %s", ifelse(is_append, "Appended", "Wrote"), nrow(dt), basename(path))
    }
}

log_step(sprintf("Processing %d total IDs in %d outer batches of %d", length(chembl_ids), n_outer, OUTER_BATCH))

for (b_idx in seq_along(chunks_outer)) {
    batch_ids <- chunks_outer[[b_idx]]
    log_step(sprintf("Outer Batch %d/%d (IDs: %d)", b_idx, n_outer, length(batch_ids)))

    # 1. Fetch molecules
    mol_df <- fetch_union_chunks("molecule", batch_ids, fields = "molecule_chembl_id", chunk_size = CHUNK, limit = LIMIT, parallel_workers = PAR)
    valid_molecule_ids <- character(0)
    if (nrow(mol_df) && "molecule_chembl_id" %in% names(mol_df)) {
        v <- unique(as.character(na.omit(mol_df$molecule_chembl_id)))
        valid_molecule_ids <- v[nzchar(v)]
    }

    if (!length(valid_molecule_ids)) {
        logf("No valid molecules in batch %d. Marking IDs complete.", b_idx)
        fwrite(data.table(batch_ids), completed_path, append = file.exists(completed_path), col.names = FALSE)
        next
    }

    # 2. Fetch molecule_form
    form_df <- fetch_union_chunks("molecule_form", valid_molecule_ids, fields = "parent_chembl_id", chunk_size = CHUNK, limit = LIMIT, parallel_workers = PAR)

    parent_child_map <- data.table()
    if (nrow(form_df)) {
        setDT(form_df)
        parent_child_map <- unique(
            form_df[
                !is.na(parent_chembl_id) & !is.na(molecule_chembl_id),
                .(
                    parent_id = as.character(parent_chembl_id),
                    child_id = as.character(molecule_chembl_id)
                )
            ]
        )
    }

    parents_of_inputs <- unique(parent_child_map[parent_id %chin% valid_molecule_ids, parent_id])
    children_of_inputs <- unique(parent_child_map[parent_id %chin% valid_molecule_ids, child_id])
    parents_from_child <- unique(parent_child_map[child_id %chin% valid_molecule_ids, parent_id])

    all_ids <- unique(c(valid_molecule_ids, parents_of_inputs, children_of_inputs, parents_from_child))

    # 3. Fetch mechanism via AnnotationGx::getChemblMechanism()
    #    Chunk manually so we can log progress and stay polite to the API.
    mech_chunks <- chunk_vec(all_ids, CHUNK)
    mech_parts <- lapply(seq_along(mech_chunks), function(ci) {
        logf("  [mechanism] chunk %d/%d (%d IDs)", ci, length(mech_chunks), length(mech_chunks[[ci]]))
        res <- tryCatch(
            getChemblMechanism(mech_chunks[[ci]]),
            error = function(e) {
                logf("  ERROR: %s", conditionMessage(e))
                NULL
            }
        )
        if (!is.null(res) && nrow(res) > 0) {
            setDT(res)
            res <- flatten_list_cols(res)
            if (endpoint_ok <- "molecule_chembl_id" %in% names(res)) res else NULL
        } else {
            NULL
        }
    })
    mech_parts <- Filter(Negate(is.null), mech_parts)
    mech_direct <- if (length(mech_parts)) rbindlist(mech_parts, fill = TRUE, use.names = TRUE) else data.table()

    remap <- data.table()
    if (nrow(parent_child_map)) {
        remap <- rbindlist(list(
            parent_child_map[, .(from_id = parent_id, to_id = child_id)],
            parent_child_map[, .(from_id = child_id, to_id = parent_id)]
        ), use.names = TRUE, fill = TRUE)
        remap <- unique(remap)
    }

    mech_remapped <- data.table()
    if (nrow(mech_direct) && nrow(remap)) {
        setDT(mech_direct)
        setDT(remap)
        mech_tmp <- merge(
            remap, mech_direct,
            by.x = "from_id", by.y = "molecule_chembl_id",
            all.x = TRUE, allow.cartesian = TRUE
        )
        if (nrow(mech_tmp)) {
            mech_tmp[, molecule_chembl_id := to_id]
            mech_tmp[, c("from_id", "to_id") := NULL]
            mech_remapped <- clean_mechanism_rows(mech_tmp)
        }
    }

    chembl_mechanism <- clean_mechanism_rows(rbindlist(list(mech_direct, mech_remapped), use.names = TRUE, fill = TRUE))
    setDT(chembl_mechanism)
    key_cols <- intersect(c("molecule_chembl_id", "record_id", "mec_id"), names(chembl_mechanism))
    if (length(key_cols) >= 2) {
        chembl_mechanism <- unique(chembl_mechanism, by = key_cols)
    } else {
        alt_cols <- intersect(c("molecule_chembl_id", "target_chembl_id", "action_type"), names(chembl_mechanism))
        if (length(alt_cols) >= 2) chembl_mechanism <- unique(chembl_mechanism, by = alt_cols)
    }

    # 4. Write data incrementally
    write_or_rm(chembl_mechanism, file.path(out_dir, "chembl_mechanism.csv"))

    # 5. Mark IDs completely tracked
    fwrite(data.table(batch_ids), completed_path, append = file.exists(completed_path), col.names = FALSE)
    logf("Batch %d fully parsed and saved to disk.", b_idx)
}

log_step("Done!")
