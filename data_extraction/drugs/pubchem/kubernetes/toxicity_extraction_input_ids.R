#!/usr/bin/env Rscript
# toxicity_extraction_input_ids.R
#
# Extracts Drug Induced Liver Injury (DILI) and Hepatotoxicity data from
# PubChem pug_view for CIDs derived from union_out.csv.
#
# Output (toxicity_output.csv) — one row per CID, wide format:
#   pubchem_cid | dili_severity_grade | dili_annotation | hepatotoxicity_likelihood_score
#
# Field name mapping handles variations across PubChem datasets:
#   dili_severity_grade  ← "Severity Class" | "Severity Grade" | "DILI severity Score" | ...
#   dili_annotation      ← "vDILI-Concern"  | "DILI Annotation" | "DILI concern" | ...
#   hepatotoxicity_*     ← any StringWithMarkup containing "Likelihood score:"
#
# Usage:
#   Rscript toxicity_extraction_input_ids.R \
#       [--union_out     <path/to/union_out.csv>] \
#       [--out_dir       <output_dir>] \
#       [--batch_input_ids <batch.csv>] \
#       [--shard_index   <0-based int>] \
#       [--num_shards    <int>]

library(httr)
library(jsonlite)
library(data.table)

`%||%` <- function(a, b) if (!is.null(a)) a else b
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# ---------------------------------------------------------------------------
# Priority-ordered field name lists (lowercase, matched case-insensitively)
# ---------------------------------------------------------------------------
SEVERITY_NAMES <- c(
    "severity class", # DILIrank 2.0
    "severity grade", # older DILIrank
    "dili severity score", # LTKB benchmark
    "dili_severity_score",
    "dili severity grade"
)

ANNOTATION_NAMES <- c(
    "vdili-concern", # DILIrank 2.0 (v-prefixed)
    "dili-concern",
    "dili annotation", # older DILIrank
    "dili concern", # LTKB benchmark
    "dilist classification", # DILIst dataset
    "dili label",
    "concern"
)

DATASET_NAMES <- c(
    "dataset"
)

# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------
batch_input_ids_csv <- NULL
out_dir_override <- NULL
union_out_override <- NULL
shard_index <- 0L
num_shards <- 1L

args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
    if (args[i] == "--batch_input_ids" && i < length(args)) {
        batch_input_ids_csv <- args[i + 1L]
        i <- i + 2L
    } else if (args[i] == "--out_dir" && i < length(args)) {
        out_dir_override <- args[i + 1L]
        i <- i + 2L
    } else if (args[i] == "--union_out" && i < length(args)) {
        union_out_override <- args[i + 1L]
        i <- i + 2L
    } else if (args[i] == "--shard_index" && i < length(args)) {
        shard_index <- as.integer(args[i + 1L])
        i <- i + 2L
    } else if (args[i] == "--num_shards" && i < length(args)) {
        num_shards <- as.integer(args[i + 1L])
        i <- i + 2L
    } else {
        i <- i + 1L
    }
}

if (num_shards < 1L) stop("--num_shards must be >= 1")
if (shard_index < 0L || shard_index >= num_shards) {
    stop("--shard_index must be in 0..(num_shards-1)")
}

out_dir <- if (!is.null(out_dir_override)) out_dir_override else "output_data/union/complete"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

prop_out <- if (!is.null(union_out_override)) union_out_override else file.path(out_dir, "union_out.csv")
out_file <- file.path(out_dir, "toxicity_output.csv")
err_log <- file.path(out_dir, "toxicity_errors.csv")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
append_csv <- function(df, path) {
    if (is.null(df) || nrow(df) == 0L) {
        return(invisible(NULL))
    }
    if (file.exists(path)) {
        write.table(df, path,
            sep = ",", row.names = FALSE, col.names = FALSE,
            append = TRUE, na = "", qmethod = "double"
        )
    } else {
        write.csv(df, path, row.names = FALSE, na = "")
    }
}

log_error <- function(pubchem_cid, detail) {
    append_csv(data.frame(
        time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        pubchem_cid = as.character(pubchem_cid),
        detail = as.character(detail),
        stringsAsFactors = FALSE
    ), err_log)
}

empty_row <- function(pubchem_cid) {
    data.frame(
        pubchem_cid                     = pubchem_cid,
        reference_number                = NA_character_,
        dili_dataset                    = NA_character_,
        dili_severity_grade             = NA_character_,
        dili_annotation                 = NA_character_,
        dili_source_url                 = NA_character_,
        hepatotoxicity_likelihood_score = NA_character_,
        stringsAsFactors                = FALSE
    )
}

# Match a field name against a priority list (case-insensitive)
matches_any <- function(name, priority_list) {
    tolower(trimws(name)) %in% priority_list
}

# ---------------------------------------------------------------------------
# Resume: which CIDs are already done?
# ---------------------------------------------------------------------------
get_done_cids <- function(path) {
    if (!file.exists(path)) {
        return(character(0))
    }
    dt <- tryCatch(fread(path, fill = TRUE), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) {
        return(character(0))
    }
    col <- if ("pubchem_cid" %in% names(dt)) "pubchem_cid" else if ("cid" %in% names(dt)) "cid" else NULL
    if (is.null(col)) {
        return(character(0))
    }
    done <- unique(trimws(as.character(dt[[col]])))
    done[nzchar(done) & done != "0"]
}

# ---------------------------------------------------------------------------
# Load union_out.csv
# ---------------------------------------------------------------------------
if (!file.exists(prop_out)) stop("union_out.csv not found: ", prop_out)
props <- fread(prop_out, fill = TRUE)
if (!"cid" %in% names(props)) stop("union_out.csv must have a 'cid' column")
props[, cid := trimws(as.character(cid))]
props <- props[nzchar(cid) & cid != "0"]

# ---------------------------------------------------------------------------
# Determine which CIDs to process
# ---------------------------------------------------------------------------
get_batch_cids <- function(props_dt, batch_input_ids_csv) {
    if (is.null(batch_input_ids_csv) || !nzchar(batch_input_ids_csv %||% "")) {
        out <- unique(trimws(as.character(props_dt$cid)))
        return(out[nzchar(out) & out != "0"])
    }
    if (!file.exists(batch_input_ids_csv)) stop("batch_input_ids file not found: ", batch_input_ids_csv)
    b <- read.csv(batch_input_ids_csv, stringsAsFactors = FALSE, check.names = FALSE)
    cn <- names(b)
    id_col <- NULL
    for (cand in c("input_id", "inputId", "substance", "name", "drug", "synonym")) {
        if (cand %in% cn) {
            id_col <- cand
            break
        }
    }
    if (is.null(id_col)) id_col <- cn[1L]
    batch_ids <- unique(trimws(as.character(b[[id_col]])))
    batch_ids <- batch_ids[nzchar(batch_ids)]
    if (!length(batch_ids)) {
        return(character(0))
    }
    u_exp <- props_dt[, .(input_id = unlist(strsplit(as.character(mapped_name), ";\\s*"))), by = cid]
    u_exp[, input_id := trimws(input_id)]
    out <- unique(u_exp[input_id %in% batch_ids, cid])
    out[nzchar(out) & out != "0"]
}

all_cids <- get_batch_cids(props, batch_input_ids_csv)

# Deterministic sharding
all_cids <- sort(all_cids)
if (num_shards > 1L) {
    idx <- which((seq_along(all_cids) - 1L) %% num_shards == shard_index)
    all_cids <- all_cids[idx]
    cat(sprintf("[%s] Shard %d/%d → %d CIDs assigned\n", ts(), shard_index + 1L, num_shards, length(all_cids)))
}

done <- get_done_cids(out_file)
todo <- setdiff(all_cids, done)

cat(sprintf(
    "[%s] Toxicity — total: %d | already done: %d | remaining: %d\n",
    ts(), length(all_cids), length(done), length(todo)
))

# ---------------------------------------------------------------------------
# PubChem API fetch
# ---------------------------------------------------------------------------
.consecutive_blocks <- 0L

# Build a lookup from ReferenceNumber -> list(source_name, url)
build_ref_lookup <- function(x) {
    refs <- x$Record$Reference %||% list()
    lut <- list()
    for (r in refs) {
        num <- as.character(r$ReferenceNumber %||% "NA")
        lut[[num]] <- list(
            source_name = r$SourceName %||% NA_character_,
            url         = r$URL %||% NA_character_
        )
    }
    lut
}

find_sections <- function(sections, heading) {
    results <- list()
    for (sec in sections) {
        if (!is.null(sec$TOCHeading) && sec$TOCHeading == heading) {
            results <- c(results, list(sec))
        }
        if (!is.null(sec$Section)) {
            results <- c(results, find_sections(sec$Section, heading))
        }
    }
    results
}

get_tox <- function(pubchem_cid) {
    url <- sprintf(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON?heading=Toxicity",
        pubchem_cid
    )

    if (.consecutive_blocks >= 10L) {
        cat(sprintf("[%s] MAX consecutive blocks — skipping CID %s\n", ts(), pubchem_cid))
        return(empty_row(pubchem_cid))
    }

    resp <- tryCatch(GET(url, timeout(30)), error = function(e) e)

    is_block <- inherits(resp, "error") ||
        (!inherits(resp, "error") && http_error(resp) && status_code(resp) %in% c(429L, 502L, 503L, 504L))

    if (is_block) {
        .consecutive_blocks <<- .consecutive_blocks + 1L
        if (.consecutive_blocks >= 3L) {
            extra <- .consecutive_blocks * 5L
            cat(sprintf("[%s] Backing off %ds (consecutive blocks: %d)\n", ts(), extra, .consecutive_blocks))
            Sys.sleep(extra)
        }
    } else {
        .consecutive_blocks <<- 0L
    }

    if (inherits(resp, "error")) {
        log_error(pubchem_cid, resp$message)
        return(empty_row(pubchem_cid))
    }
    if (http_error(resp)) {
        log_error(pubchem_cid, paste("HTTP", status_code(resp)))
        return(empty_row(pubchem_cid))
    }

    x <- tryCatch(content(resp, as = "parsed", type = "application/json"), error = function(e) e)
    if (inherits(x, "error") || is.null(x$Record$Section)) {
        log_error(pubchem_cid, "Missing Record/Section")
        return(empty_row(pubchem_cid))
    }

    ref_lookup <- build_ref_lookup(x)
    all_sections <- x$Record$Section %||% list()

    # ---- DILI: scan all Information entries, pick by priority name match ----
    dili_sections <- find_sections(all_sections, "Drug Induced Liver Injury")
    severity_grade <- NA_character_

    # Group Information entries by ReferenceNumber, preserving insertion order
    ref_order <- character(0) # ordered unique ref nums
    ref_groups <- list() # ref_num -> list of (name, value)

    for (sec in dili_sections) {
        if (is.null(sec$Information)) next
        for (info in sec$Information) {
            nm <- info$Name %||% NA_character_
            ref_num <- as.character(info$ReferenceNumber %||% "NA")
            swm <- info$Value$StringWithMarkup
            if (is.na(nm) || is.null(swm) || length(swm) == 0L) next
            val <- paste(
                Filter(nzchar, vapply(swm, function(s) s$String %||% "", character(1))),
                collapse = "; "
            )
            if (!nzchar(val)) next
            if (is.null(ref_groups[[ref_num]])) {
                ref_order <- c(ref_order, ref_num)
                ref_groups[[ref_num]] <- list()
            }
            ref_groups[[ref_num]] <- c(ref_groups[[ref_num]], list(list(name = nm, value = val)))
        }
    }

    # Build one row per reference
    dili_rows <- list()
    for (ref_num in ref_order) {
        entries <- ref_groups[[ref_num]]
        sev_grade <- NA_character_
        ann <- NA_character_
        dataset <- NA_character_
        ref_entry <- ref_lookup[[ref_num]]
        source_url <- if (!is.null(ref_entry)) ref_entry$url %||% NA_character_ else NA_character_

        for (e in entries) {
            if (is.na(sev_grade) && matches_any(e$name, SEVERITY_NAMES)) sev_grade <- e$value
            if (is.na(ann) && matches_any(e$name, ANNOTATION_NAMES)) ann <- e$value
            if (is.na(dataset) && matches_any(e$name, DATASET_NAMES)) dataset <- e$value
        }

        # Only emit a row if the reference contributed at least one meaningful DILI field
        if (!is.na(dataset) || !is.na(sev_grade) || !is.na(ann)) {
            dili_rows[[length(dili_rows) + 1L]] <- data.frame(
                pubchem_cid = pubchem_cid,
                reference_number = ref_num,
                dili_dataset = dataset,
                dili_severity_grade = sev_grade,
                dili_annotation = ann,
                dili_source_url = source_url,
                stringsAsFactors = FALSE
            )
        }
    }

    # ---- Hepatotoxicity: extract "Likelihood score: X" string ----
    hep_sections <- find_sections(all_sections, "Hepatotoxicity")
    hep_score <- NA_character_

    for (sec in hep_sections) {
        if (is.null(sec$Information)) next
        for (info_entry in sec$Information) {
            swm_list <- info_entry$Value$StringWithMarkup
            if (is.null(swm_list)) next
            for (swm in swm_list) {
                if (!is.null(swm$String) && grepl("Likelihood score:", swm$String, fixed = TRUE)) {
                    hep_score <- swm$String
                    break
                }
            }
            if (!is.na(hep_score)) break
        }
        if (!is.na(hep_score)) break
    }

    # Fan hep_score across all DILI rows; fall back to one NA sentinel if no DILI data
    if (length(dili_rows) > 0L) {
        out <- do.call(rbind, dili_rows)
        out$hepatotoxicity_likelihood_score <- hep_score
        out
    } else {
        data.frame(
            pubchem_cid = pubchem_cid,
            reference_number = NA_character_,
            dili_dataset = NA_character_,
            dili_severity_grade = NA_character_,
            dili_annotation = NA_character_,
            dili_source_url = NA_character_,
            hepatotoxicity_likelihood_score = hep_score,
            stringsAsFactors = FALSE
        )
    }
}

# ---------------------------------------------------------------------------
# Main loop — flush every 25 CIDs
# ---------------------------------------------------------------------------
if (!length(todo)) {
    cat(sprintf("[%s] Nothing to do. Exiting.\n", ts()))
    quit(save = "no", status = 0)
}

batch_size <- 25L
batches <- split(todo, ceiling(seq_along(todo) / batch_size))

for (b in seq_along(batches)) {
    batch <- batches[[b]]
    cat(sprintf("[%s] Toxicity batch %d/%d | n=%d\n", ts(), b, length(batches), length(batch)))

    rows <- vector("list", length(batch))

    for (i in seq_along(batch)) {
        cid <- batch[i]
        cat(sprintf("[%s] (%d/%d) CID %s\n", ts(), i, length(batch), cid))

        rows[[i]] <- tryCatch(get_tox(cid), error = function(e) {
            log_error(cid, e$message)
            empty_row(cid)
        })

        Sys.sleep(runif(1, 0.75, 1.5))
    }

    out <- do.call(rbind, rows)
    append_csv(out, out_file)
    cat(sprintf("[%s] Wrote %d rows to %s\n", ts(), nrow(out), out_file))
}

cat(sprintf("[%s] DONE toxicity.\n", ts()))
