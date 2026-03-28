#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(httr)
    library(jsonlite)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
printf <- function(...) {
    cat(sprintf(...), "\n")
    flush.console()
}

# ---------------------------
# Args
# ---------------------------
batch_ids_path <- NULL
out_dir_override <- NULL
args <- commandArgs(trailingOnly = TRUE)
i <- 1
while (i <= length(args)) {
    if (args[i] == "--batch_input_ids" && i < length(args)) {
        batch_ids_path <- args[i + 1]
        i <- i + 2
    } else if (args[i] == "--out_dir" && i < length(args)) {
        out_dir_override <- args[i + 1]
        i <- i + 2
    } else {
        i <- i + 1
    }
}
if (is.null(batch_ids_path) || !nzchar(batch_ids_path)) {
    stop("Usage: Rscript assay_extraction_input_ids.R --batch_input_ids <batch_input_ids.csv> [--out_dir <dir>]")
}

# ---------------------------
# CONFIG (USE SAME OUTPUTS AS INPUT-ID WORKFLOW)
# ---------------------------
out_dir <- if (!is.null(out_dir_override)) out_dir_override else "output_data/union/complete"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

bio_out <- file.path(out_dir, "union_bioassays.csv")
union_out_path <- file.path(out_dir, "union_out.csv")
assay_out <- file.path(out_dir, "union_pubchem_assay_fields.csv")
err_log <- file.path(out_dir, "assay_extraction_errors.csv")

base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid"

if (!file.exists(batch_ids_path)) stop("Batch file does not exist: ", batch_ids_path)

if (!file.exists(bio_out)) {
    printf("[%s] No union_bioassays.csv yet (%s). Nothing to do.", ts(), bio_out)
    quit(save = "no", status = 0)
}


# ---------------------------
# I/O helpers (safe append)
# ---------------------------
append_fwrite <- function(dt, path, label) {
    if (is.null(dt) || !nrow(dt)) {
        return(invisible(NULL))
    }

    file_exists <- file.exists(path)

    if (file_exists) {
        # Lock to existing header to prevent column-misalignment when appending.
        hdr <- names(fread(path, nrows = 0))
        for (h in setdiff(hdr, names(dt))) dt[, (h) := NA]
        extra <- setdiff(names(dt), hdr)
        dt <- dt[, c(hdr, extra), with = FALSE]
        fwrite(dt, path, append = TRUE, col.names = FALSE)
    } else {
        fwrite(dt, path, append = FALSE, col.names = TRUE)
    }

    printf(
        "[%s] %s: wrote %d rows to %s (append=%s)",
        ts(), label, nrow(dt), path, if (file_exists) "TRUE" else "FALSE"
    )
}

log_error <- function(aid, detail) {
    row <- data.table(time = ts(), aid = as.character(aid), detail = as.character(detail))
    append_fwrite(row, err_log, "ERROR_LOG")
}

# ---------------------------
# Load batch input_ids
# ---------------------------
bdt <- fread(batch_ids_path, fill = TRUE)
# The coordinator now passes full rows. We need to identify the ID column.
cn <- names(bdt)
input_id_col <- NULL
for (cand in c("input_id", "inputId", "substance", "name", "drug", "synonym")) {
    if (cand %in% cn) {
        input_id_col <- cand
        break
    }
}
if (is.null(input_id_col)) input_id_col <- cn[1]

batch_ids <- trimws(as.character(bdt[[input_id_col]]))
batch_ids <- batch_ids[nzchar(batch_ids)]
batch_ids <- unique(batch_ids)

printf("[%s] Batch input_ids: %d", ts(), length(batch_ids))
if (!length(batch_ids)) quit(save = "no", status = 0)

# ---------------------------
# Resolve batch CIDs from union_out.csv
# ---------------------------
if (!file.exists(union_out_path)) {
    printf("[%s] No union_out.csv yet (%s). Nothing to do.", ts(), union_out_path)
    quit(save = "no", status = 0)
}

# Load only the necessary columns from union_out
u <- fread(union_out_path, select = c("cid", "mapped_name"), colClasses = "character")
u[, cid := trimws(cid)]

# Create an expanded mapping: CID -> Single JCP ID
u_expanded <- u[, .(input_id = unlist(strsplit(as.character(mapped_name), ";\\s*"))), by = cid]
u_expanded[, input_id := trimws(input_id)]

# Map batch input_ids -> unique CIDs
batch_cids <- unique(u_expanded[input_id %in% batch_ids, cid])
batch_cids <- batch_cids[nzchar(batch_cids) & grepl("^[0-9]+$", batch_cids)]

printf("[%s] Batch CIDs: %d", ts(), length(batch_cids))
if (!length(batch_cids)) {
    printf("[%s] WARNING: No CIDs found in union_out for this batch of input_ids.", ts())
    quit(save = "no", status = 0)
}

# ---------------------------
# AIDs linked to batch CIDs
# ---------------------------
bio <- fread(bio_out, fill = TRUE)
if (!all(c("pubchem_cid", "bioassay_aid") %in% names(bio))) {
    stop("union_bioassays.csv must contain pubchem_cid and bioassay_aid columns.")
}
bio[, pubchem_cid := trimws(as.character(pubchem_cid))]
bio[, bioassay_aid := trimws(as.character(bioassay_aid))]
bio <- bio[pubchem_cid %in% batch_cids]
aids <- unique(bio$bioassay_aid)
aids <- aids[nzchar(aids) & grepl("^[0-9]+$", aids)]

printf("[%s] AIDs linked to batch CIDs: %d", ts(), length(aids))
if (!length(aids)) quit(save = "no", status = 0)

# ---------------------------
# Which AIDs are already done?
# ---------------------------
done_aids <- character(0)
if (file.exists(assay_out)) {
    adt <- tryCatch(fread(assay_out, select = "aid", fill = TRUE), error = function(e) NULL)
    if (!is.null(adt) && nrow(adt)) {
        done_aids <- unique(trimws(as.character(adt$aid)))
        done_aids <- done_aids[nzchar(done_aids)]
    }
}

completed_aids_path <- "/tmp/work/completed_aids.txt"
if (file.exists(completed_aids_path)) {
    lines <- trimws(readLines(completed_aids_path, warn = FALSE))
    lines <- lines[nzchar(lines)]
    if (length(lines) > 0) {
        done_aids <- unique(c(done_aids, lines))
        printf("[%s] Loaded %d completed AIDs globally to prevent re-scraping.", ts(), length(lines))
    }
}
todo_aids <- setdiff(aids, done_aids)

printf("[%s] AID fields — done: %d | todo: %d", ts(), length(done_aids), length(todo_aids))
if (!length(todo_aids)) {
    printf("[%s] DONE assay_extraction_input_ids (nothing to do).", ts())
    quit(save = "no", status = 0)
}

# ---------------------------
# PubChem AID -> fields (manual parsing; no AnnotationGx)
# ---------------------------
.consecutive_api_blocks_assay <- 0L

get_assay_row <- function(aid, max_attempts = 5, base_delay = 1) {
    if (.consecutive_api_blocks_assay >= 10L) return(NULL)

    url <- sprintf("%s/%s/description/JSON", base_url, aid)

    resp <- NULL
    x <- NULL

    for (attempt in seq_len(max_attempts)) {
        resp <- tryCatch(GET(url, timeout(60)), error = function(e) e)

        # Network-layer error
        if (inherits(resp, "error")) {
            .consecutive_api_blocks_assay <<- .consecutive_api_blocks_assay + 1L
            if (attempt < max_attempts) {
                Sys.sleep(base_delay * 2^(attempt - 1))
                next
            } else {
                log_error(aid, paste("network error:", resp$message))
                return(NULL)
            }
        }

        # HTTP error
        if (http_error(resp)) {
            status <- status_code(resp)
            if (status %in% c(429, 500:599, 404)) {
                .consecutive_api_blocks_assay <<- .consecutive_api_blocks_assay + 1L
            } else {
                .consecutive_api_blocks_assay <<- 0L
            }

            if (attempt < max_attempts && status %in% c(429, 500:599)) {
                Sys.sleep(base_delay * 2^(attempt - 1))
                next
            } else {
                log_error(aid, paste("HTTP", status))
                return(NULL)
            }
        }
        
        .consecutive_api_blocks_assay <<- 0L

        # Try to parse JSON
        x <- tryCatch(content(resp, as = "parsed", type = "application/json"), error = function(e) e)
        if (inherits(x, "error")) {
            if (attempt < max_attempts) {
                Sys.sleep(base_delay * 2^(attempt - 1))
                next
            } else {
                log_error(aid, paste("JSON parse error:", x$message))
                return(NULL)
            }
        }

        break
    }

    if (is.null(x) || is.null(x$PC_AssayContainer)) {
        log_error(aid, "no usable JSON (missing PC_AssayContainer)")
        return(NULL)
    }

    container <- x$PC_AssayContainer[[1]]
    assay <- container$assay
    descr <- assay$descr

    # Basic fields
    aid_id <- descr$aid$id %||% aid
    version <- descr$aid$version %||% NA_character_
    assay_name <- descr$name %||% NA_character_

    # Source
    source_name <- descr$aid_source$db$name %||% NA_character_
    source_id <- descr$aid_source$db$source_id$str %||% NA_character_

    # Description/protocol/comment (concatenate arrays)
    description_combined <- if (!is.null(descr$description)) {
        paste(unlist(descr$description), collapse = " | ")
    } else {
        NA_character_
    }

    protocol_combined <- if (!is.null(descr$protocol)) {
        paste(unlist(descr$protocol), collapse = " | ")
    } else {
        NA_character_
    }

    comment_combined <- if (!is.null(descr$comment)) {
        paste(unlist(descr$comment), collapse = " | ")
    } else {
        NA_character_
    }

    # activity_outcome_method (keep only 1 or 2)
    aom <- descr$activity_outcome_method %||% NA_integer_
    if (!(aom %in% c(1, 2))) aom <- NA_integer_

    # Targets
    targets <- descr$target %||% list()

    target_name <- if (length(targets)) {
        paste(
            vapply(targets, function(t) t$name %||% NA_character_, character(1)),
            collapse = " | "
        )
    } else {
        NA_character_
    }

    target_protein_accession <- if (length(targets)) {
        paste(
            vapply(
                targets,
                function(t) t$mol_id$protein_accession %||% NA_character_,
                character(1)
            ),
            collapse = " | "
        )
    } else {
        NA_character_
    }

    data.table(
        aid                      = as.character(aid_id),
        version                  = as.character(version),
        assay_name               = as.character(assay_name),
        source_name              = as.character(source_name),
        source_id                = as.character(source_id),
        description_combined     = as.character(description_combined),
        protocol_combined        = as.character(protocol_combined),
        comment_combined         = as.character(comment_combined),
        activity_outcome_method  = as.integer(aom),
        target_name              = as.character(target_name),
        target_protein_accession = as.character(target_protein_accession)
    )
}

# ---------------------------
# Fetch AID fields in mini-batches
# ---------------------------
mini_size <- 25L
batches <- split(todo_aids, ceiling(seq_along(todo_aids) / mini_size))

for (b in seq_along(batches)) {
    sub_aids <- batches[[b]]
    printf("[%s] AID batch %d/%d | n=%d", ts(), b, length(batches), length(sub_aids))

    out_list <- list()

    for (aid in sub_aids) {
        row <- get_assay_row(aid)
        if (!is.null(row)) out_list[[aid]] <- row
        Sys.sleep(runif(1, 0.75, 1.5)) # throttle dynamically across pods; tune as needed
    }

    if (!length(out_list)) next

    out_dt <- rbindlist(out_list, fill = TRUE)
    if (!nrow(out_dt)) next

    append_fwrite(out_dt, assay_out, "AssayFields")
}

printf("[%s] DONE assay_extraction_input_ids.", ts())
