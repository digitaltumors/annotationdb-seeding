#!/usr/bin/env Rscript
# chembl_mechanism_annotationgx.R
#
# Extracts mechanism data for all ChEMBL IDs present in union_out.csv
# using AnnotationGx::getChemblMechanism() (hits ChEMBL REST API).
#
# Output: chembl_mechanism_annotationgx.csv in the same output directory.
# Supports incremental resumption via chembl_mechanism_agx_completed_ids.csv.

suppressPackageStartupMessages({
    library(data.table)
    library(AnnotationGx)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
INPUT_CSV   <- "../pubchem/output_data/union/mar-27-2026/union_out.csv"
OUTPUT_DIR  <- "../pubchem/output_data/union/mar-27-2026"
OUTPUT_CSV  <- file.path(OUTPUT_DIR, "chembl_mechanism_annotationgx.csv")
DONE_CSV    <- file.path(OUTPUT_DIR, "chembl_mechanism_agx_completed_ids.csv")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Config ─────────────────────────────────────────────────────────────────────
CHUNK_SIZE  <- 50     # IDs per getChemblMechanism() call (uses __in filter)
SLEEP_SEC   <- 1.5    # polite delay between chunks

# ── Logging ────────────────────────────────────────────────────────────────────
ts   <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logf <- function(fmt, ...) cat(sprintf("[%s] %s\n", ts(), sprintf(fmt, ...)))

# ── Helpers ────────────────────────────────────────────────────────────────────
chunk_vec <- function(x, n) split(x, ceiling(seq_along(x) / n))

extract_chembl_ids <- function(path) {
    if (!file.exists(path) || file.info(path)$size == 0) return(character(0))
    dt <- tryCatch(
        fread(path, fill = TRUE, quote = "\"'", showProgress = FALSE),
        error = function(e) NULL
    )
    if (is.null(dt)) return(character(0))
    vals <- unlist(lapply(dt, as.character), use.names = FALSE)
    vals <- vals[!is.na(vals) & nzchar(vals)]
    unique(unlist(regmatches(vals, gregexpr("CHEMBL\\d+", vals, perl = TRUE))))
}

flatten_variant_sequence <- function(dt) {
    # variant_sequence arrives as a list-column of data.frames or NULLs.
    # Flatten to JSON string so the CSV stays flat.
    if (!("variant_sequence" %in% names(dt))) return(dt)
    dt[, variant_sequence := vapply(variant_sequence, function(x) {
        if (is.null(x) || (is.logical(x) && is.na(x))) return(NA_character_)
        tryCatch(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"),
                 error = function(e) NA_character_)
    }, character(1))]
    dt
}

flatten_mechanism_refs <- function(dt) {
    # mechanism_refs is also a list-column.
    if (!("mechanism_refs" %in% names(dt))) return(dt)
    dt[, mechanism_refs := vapply(mechanism_refs, function(x) {
        if (is.null(x) || (is.logical(x) && is.na(x))) return(NA_character_)
        tryCatch(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"),
                 error = function(e) NA_character_)
    }, character(1))]
    dt
}

# ── Collect IDs ────────────────────────────────────────────────────────────────
logf("Extracting ChEMBL IDs from %s ...", INPUT_CSV)
all_ids <- extract_chembl_ids(INPUT_CSV)
logf("Found %d unique ChEMBL IDs.", length(all_ids))

# Resume support
if (file.exists(DONE_CSV)) {
    done_ids <- as.character(fread(DONE_CSV, header = FALSE)[[1]])
    before   <- length(all_ids)
    all_ids  <- setdiff(all_ids, done_ids)
    logf("Resuming: skipping %d already-completed IDs. %d remaining.", before - length(all_ids), length(all_ids))
}

if (!length(all_ids)) {
    logf("All IDs already processed. Nothing to do.")
    quit(save = "no", status = 0)
}

# ── Chunked extraction ─────────────────────────────────────────────────────────
chunks   <- chunk_vec(all_ids, CHUNK_SIZE)
n_chunks <- length(chunks)
logf("Processing %d IDs in %d chunks of up to %d.", length(all_ids), n_chunks, CHUNK_SIZE)

for (i in seq_along(chunks)) {
    ids <- chunks[[i]]
    logf("[%d/%d] Querying %d IDs ...", i, n_chunks, length(ids))

    result <- tryCatch(
        getChemblMechanism(ids),
        error = function(e) {
            logf("  ERROR in chunk %d: %s", i, conditionMessage(e))
            NULL
        }
    )

    if (!is.null(result) && nrow(result) > 0) {
        setDT(result)
        result <- flatten_variant_sequence(result)
        result <- flatten_mechanism_refs(result)
        is_append <- file.exists(OUTPUT_CSV)
        fwrite(result, OUTPUT_CSV, append = is_append)
        logf("  %s %d rows → %s", ifelse(is_append, "Appended", "Wrote"), nrow(result), basename(OUTPUT_CSV))
    } else {
        logf("  No mechanism rows returned for this chunk.")
    }

    # Record chunk as done
    fwrite(data.table(ids), DONE_CSV, append = file.exists(DONE_CSV), col.names = FALSE)

    if (i < n_chunks) Sys.sleep(SLEEP_SEC)
}

logf("Done! Output: %s", OUTPUT_CSV)
