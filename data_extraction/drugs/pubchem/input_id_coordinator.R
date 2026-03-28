#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
printf <- function(...) {
    cat(sprintf(...), "\n")
    flush.console()
}

# ---------------------------
# CONFIGURATION
# ---------------------------
# Edit these values to change the default behavior
in_csv <- "output_data/union/complete/DFM-JUMP-CP-Training-Combined.csv" # The master CSV mapping input_ids to PubChem CIDs
out_dir <- "output_data/union/complete"
batch_size <- 50
max_cycles <- 999999

# ---------------------------
# ARGUMENT PARSING (Overrides CONFIG)
# ---------------------------
parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    out <- list()
    key <- NULL
    for (a in args) {
        if (startsWith(a, "--")) {
            key <- sub("^--", "", a)
            out[[key]] <- TRUE
        } else if (!is.null(key)) {
            out[[key]] <- a
            key <- NULL
        }
    }
    out
}

args <- parse_args()
if (!is.null(args$in_csv)) in_csv <- args$in_csv
if (!is.null(args$batch_size)) batch_size <- as.integer(args$batch_size)
if (!is.null(args$max_cycles)) max_cycles <- as.integer(args$max_cycles)
if (!is.null(args$out_dir)) out_dir <- args$out_dir

# ---------------------------
# SETUP
# ---------------------------
if (!file.exists(in_csv)) {
    stop(sprintf("Input CSV missing: %s. Please edit the 'in_csv' variable in the script or provide it via --in_csv.", in_csv))
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

batch_input_ids_path <- file.path(out_dir, "batch_input_ids.csv")
union_out_path <- file.path(out_dir, "union_out.csv")
failed_input_ids_path <- file.path(out_dir, "failed_input_ids.csv")

# These are the scripts you run in the coordinator
pubchem_cmd <- function() {
    sprintf(
        "Rscript pubchem_drug_extraction_from_input_ids.R --in_csv '%s'",
        batch_input_ids_path
    )
}
assay_cmd <- function() {
    sprintf(
        "Rscript assay_extraction_input_ids.R --batch_input_ids '%s'",
        batch_input_ids_path
    )
}
tox_cmd <- function() {
    sprintf(
        "Rscript toxicity_extraction_input_ids.R --batch_input_ids '%s'",
        batch_input_ids_path
    )
}

# Load master mapping
printf("[%s] Loading master input mapping: %s", ts(), in_csv)
master <- fread(in_csv)

cn <- names(master)
# CID column guess
cid_col <- NULL
for (cand in c("cid", "CID", "pubchem_cid", "PubChemCID")) {
    if (cand %in% cn) {
        cid_col <- cand
        break
    }
}
if (is.null(cid_col)) {
    hits <- grep("^cid$", cn, ignore.case = TRUE, value = TRUE)
    if (length(hits) >= 1) cid_col <- hits[1]
}

# input_id column guess
input_id_col <- NULL
for (cand in c("input_id", "inputId", "substance", "name", "drug", "synonym")) {
    if (cand %in% cn) {
        input_id_col <- cand
        break
    }
}
if (is.null(input_id_col)) {
    non_cid <- setdiff(cn, cid_col %||% character(0))
    if (length(non_cid) < 1) stop("Could not infer input_id column in master CSV.")
    input_id_col <- non_cid[1]
}

# InChIKey column guess
inchikey_col <- NULL
for (cand in c("InChiKey", "inchikey", "InChIKey")) {
    if (cand %in% cn) {
        inchikey_col <- cand
        break
    }
}

if (is.null(cid_col) && is.null(inchikey_col)) {
    stop("Master CSV must contain either a 'cid' or an 'InChiKey' column for resolution.")
}

master_ids <- unique(as.character(master[[input_id_col]]))
master_ids <- master_ids[!is.na(master_ids) & nzchar(master_ids)]
printf("[%s] Master unique input_ids: %d", ts(), length(master_ids))

get_done_ids <- function() {
    done <- character()

    # Cross-reference with toxicity_output.csv for true completion
    tox_path <- file.path(out_dir, "toxicity_output.csv")
    finished_cids <- NULL
    if (file.exists(tox_path)) {
        tox <- tryCatch(fread(tox_path, select = "pubchem_cid", fill = TRUE), error = function(e) NULL)
        if (!is.null(tox)) finished_cids <- unique(trimws(as.character(tox$pubchem_cid)))
    }

    # Done via union_out mapped_name
    if (file.exists(union_out_path)) {
        u <- tryCatch(fread(union_out_path, select = c("cid", "mapped_name")), error = function(e) NULL)
        if (!is.null(u) && "mapped_name" %in% names(u)) {
            # If we have finished CIDs, only consider those as 'done'
            if (!is.null(finished_cids)) {
                u <- u[trimws(as.character(cid)) %in% finished_cids]
            }

            names_raw <- u$mapped_name
            names_raw <- names_raw[!is.na(names_raw) & nzchar(names_raw)]
            if (length(names_raw) > 0) {
                parts <- unlist(strsplit(names_raw, ";\\s*"))
                parts <- parts[!is.na(parts) & nzchar(parts)]
                done <- unique(c(done, parts))
            }
        }
    }

    # Done via failed_input_ids (terminal failures)
    if (file.exists(failed_input_ids_path)) {
        f <- tryCatch(fread(failed_input_ids_path), error = function(e) NULL)
        if (!is.null(f) && "input_id" %in% names(f)) {
            fi <- unique(as.character(f$input_id))
            fi <- fi[!is.na(fi) & nzchar(fi)]
            done <- unique(c(done, fi))
        }
    }

    done
}

# Helper: run command and check exit code
run <- function(cmd, stop_on_error = TRUE) {
    printf("[%s] RUN: %s", ts(), cmd)
    code <- system(cmd)
    if (code != 0) {
        msg <- sprintf("Command failed with exit code %s: %s", code, cmd)
        if (stop_on_error) {
            stop(msg)
        } else {
            printf("[%s] WARNING: %s", ts(), msg)
        }
    }
    code
}

# ---------------------------
# INITIALIZATION (Load existing progress once)
# ---------------------------
printf("[%s] Initializing progress... (Reading %s, this may take a few minutes)", ts(), union_out_path)
start_time <- Sys.time()
done_ids <- get_done_ids()
end_time <- Sys.time()
printf(
    "[%s] Initial progress loaded in %.2f seconds. Total done: %d",
    ts(), as.numeric(difftime(end_time, start_time, units = "secs")), length(done_ids)
)

cycle <- 0
repeat {
    cycle <- cycle + 1
    # Optimization: Use the cached 'done_ids' instead of re-reading files every time
    remaining <- setdiff(master_ids, done_ids)

    total_m <- length(master_ids)
    printf("[%s] Done: %d/%d | Remaining: %d", ts(), length(done_ids), total_m, length(remaining))

    if (length(remaining) == 0) {
        printf("[%s] All done. Exiting.", ts())
        break
    }

    if (cycle > max_cycles) {
        stop(sprintf("Reached max_cycles=%d with %d remaining input_ids.", max_cycles, length(remaining)))
    }

    batch_ids <- head(remaining, batch_size)
    printf("[%s] Cycle %d | Taking batch: %d input_ids (batch_size=%d)", ts(), cycle, length(batch_ids), batch_size)

    # Unified batching: Write the full master rows for this batch
    batch_rows <- master[get(input_id_col) %in% batch_ids]
    fwrite(batch_rows, batch_input_ids_path)
    printf("[%s] Wrote batch file: %s", ts(), batch_input_ids_path)

    # Run pipeline steps
    # Core script: if this fails, we STOP because it's the foundation for the others
    run(pubchem_cmd(), stop_on_error = TRUE)

    # Follow-up scripts: we warn but continue,as the gap-detection will pick them up later
    run(assay_cmd(), stop_on_error = FALSE)
    run(tox_cmd(), stop_on_error = FALSE)

    # Optimization: Update cached progress locally so we don't reload the 1GB file
    done_ids <- unique(c(done_ids, batch_ids))

    printf("[%s] Finished cycle %d.", ts(), cycle)
}
