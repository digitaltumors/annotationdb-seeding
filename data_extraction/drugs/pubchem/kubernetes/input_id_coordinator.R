#!/usr/bin/env Rscript
# For a local E2E scraping

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
in_csv <- NULL # REQUIRED: pass via --in_csv
out_dir <- tempfile("annotationdb_run_") # default: OS temp dir
batch_size <- 50
max_cycles <- 999999
gcs_bucket <- NULL
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S") # unique label for this run's batches
merge_py <- NULL

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
if (!is.null(args$gcs_bucket)) gcs_bucket <- args$gcs_bucket
if (!is.null(args$run_id)) run_id <- args$run_id
if (!is.null(args$merge_py)) merge_py <- args$merge_py

if (is.null(in_csv)) stop("--in_csv is required. Provide the master input CSV path.")

# Auto-detect merge_outputs.py next to this script if not provided
if (is.null(merge_py)) {
    # Reliably find the script's own directory under Rscript
    script_args <- commandArgs(trailingOnly = FALSE)
    file_flag <- script_args[startsWith(script_args, "--file=")]
    script_dir <- if (length(file_flag) > 0) {
        dirname(normalizePath(sub("^--file=", "", file_flag[1])))
    } else {
        getwd() # fallback when sourced interactively
    }
    candidate <- file.path(script_dir, "merge_outputs.py")
    if (file.exists(candidate)) merge_py <- candidate
}

# ---------------------------
# SETUP
# ---------------------------
if (!file.exists(in_csv)) {
    stop(sprintf("Input CSV missing: %s. Provide it via --in_csv.", in_csv))
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

batch_input_ids_path <- file.path(out_dir, "batch_input_ids.csv")
union_out_path <- file.path(out_dir, "union_out.csv")
failed_input_ids_path <- file.path(out_dir, "failed_input_ids.csv")

printf("[%s] Working directory: %s", ts(), out_dir)
if (!is.null(gcs_bucket)) printf("[%s] GCS bucket: %s  run_id: %s", ts(), gcs_bucket, run_id)

# ---------------------------------------------------------------------------
# Download GCS master progress files for resume
# Only downloads a file if it doesn't already exist locally.
# ---------------------------------------------------------------------------
if (!is.null(gcs_bucket)) {
    printf("[%s] Downloading existing progress from GCS for resume...", ts())
    for (fname in c("union_out.csv", "toxicity_output.csv", "failed_input_ids.csv")) {
        local_path <- file.path(out_dir, fname)
        if (!file.exists(local_path)) {
            gcs_src <- sprintf("%s/output/%s", gcs_bucket, fname)
            ret <- system(sprintf("gsutil -o 'GSUtil:parallel_process_count=1' -q cp '%s' '%s' 2>/dev/null", gcs_src, local_path))
            if (ret == 0) {
                printf("[%s]   Restored %-35s from GCS.", ts(), fname)
            } else {
                printf("[%s]   %-35s not in GCS yet (fresh start).", ts(), fname)
            }
        } else {
            printf("[%s]   %-35s already exists locally (skipping download).", ts(), fname)
        }
    }
}

# ---------------------------------------------------------------------------
# Pipeline step commands — all write to the shared out_dir
# ---------------------------------------------------------------------------
pubchem_cmd <- function() {
    sprintf(
        "Rscript pubchem_drug_extraction_from_input_ids.R --in_csv '%s' --out_dir '%s'",
        batch_input_ids_path, out_dir
    )
}
assay_cmd <- function() {
    sprintf(
        "Rscript assay_extraction_input_ids.R --batch_input_ids '%s' --out_dir '%s'",
        batch_input_ids_path, out_dir
    )
}
tox_cmd <- function() {
    sprintf(
        "Rscript toxicity_extraction_input_ids.R --batch_input_ids '%s' --out_dir '%s'",
        batch_input_ids_path, out_dir
    )
}

# Upload out_dir CSVs to a named GCS output_batches folder
gcs_upload_batch <- function(cycle_num) {
    if (is.null(gcs_bucket)) {
        return(invisible(NULL))
    }
    batch_name <- sprintf("local_%s_c%04d", run_id, cycle_num)
    gcs_dest <- sprintf("%s/output_batches/%s/", gcs_bucket, batch_name)
    printf("[%s] Uploading cycle %d outputs → %s", ts(), cycle_num, gcs_dest)
    # Use || true so a partial upload (e.g. some files missing) doesn't abort
    ret <- system(sprintf(
        "gsutil -m -q cp '%s'/*.csv '%s' 2>/dev/null || true",
        out_dir, gcs_dest
    ))
    printf("[%s] GCS batch upload done: %s", ts(), batch_name)
    invisible(ret)
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
    printf("[%s] Master CSV does not contain a 'cid' or 'InChiKey' column. Will attempt Name-to-CID resolution via AnnotationGx.", ts())
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

        # Run merge_outputs.py to consolidate all GCS batches into master output/
        if (!is.null(gcs_bucket) && !is.null(merge_py)) {
            printf("[%s] Running merge_outputs.py to consolidate GCS master files...", ts())
            ret <- system(sprintf("python3 '%s'", merge_py))
            if (ret == 0) {
                printf("[%s] Merge complete. Master files updated at %s/output/", ts(), gcs_bucket)
            } else {
                printf("[%s] WARNING: merge_outputs.py exited with code %d — run it manually.", ts(), ret)
            }
        } else if (!is.null(gcs_bucket)) {
            printf("[%s] NOTE: pass --merge_py <path/to/merge_outputs.py> to auto-merge after run.", ts())
        }

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

    # Follow-up scripts: we warn but continue, gap-detection will pick them up later
    run(assay_cmd(), stop_on_error = FALSE)
    run(tox_cmd(), stop_on_error = FALSE)

    # Upload this cycle's outputs to GCS output_batches/
    gcs_upload_batch(cycle)

    # Optimization: Update cached progress locally so we don't reload the 1GB file
    done_ids <- unique(c(done_ids, batch_ids))

    printf("[%s] Finished cycle %d.", ts(), cycle)
}
