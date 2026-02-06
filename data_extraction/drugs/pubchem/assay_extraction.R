#!/usr/bin/env Rscript

library(jsonlite)
library(data.table)
library(httr)

`%||%` <- function(a, b) if (!is.null(a)) a else b
printf <- function(...) { cat(sprintf(...), "\n"); flush.console() }
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# ---- Config ----
in_csv  <- "output_data/overlap/complete/overlap_bioassays.csv"
out_csv <- "output_data/overlap/complete/overlap_pubchem_assay_fields.csv"
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid"

# ---- Functions ----

get_assay_row <- function(aid, max_attempts = 3, base_delay = 1) {
    url <- sprintf("%s/%s/description/JSON", base_url, aid)

    resp <- NULL
    x    <- NULL

    # ------------- RETRY LOOP -------------
    for (attempt in seq_len(max_attempts)) {

        resp <- tryCatch(GET(url, timeout(60)), error = function(e) e)

        # Network-layer error
        if (inherits(resp, "error")) {
            if (attempt < max_attempts) {
                printf("[%s] AID %s attempt %d/%d: network error: %s — retrying after %ds",
                       ts(), aid, attempt, max_attempts, resp$message, base_delay * 2^(attempt - 1))
                Sys.sleep(base_delay * 2^(attempt - 1))
                next
            } else {
                printf("[%s] AID %s failed after %d attempts (network error: %s)",
                       ts(), aid, max_attempts, resp$message)
                return(NULL)
            }
        }

        # HTTP error
        if (http_error(resp)) {
            status <- status_code(resp)

            # Retry only for transient-ish codes
            if (attempt < max_attempts && status %in% c(429, 500:599)) {
                printf("[%s] AID %s attempt %d/%d: HTTP %s — retrying after %ds",
                       ts(), aid, attempt, max_attempts, status, base_delay * 2^(attempt - 1))
                Sys.sleep(base_delay * 2^(attempt - 1))
                next
            } else {
                printf("[%s] AID %s failed with HTTP %s after %d attempts",
                       ts(), aid, status, attempt)
                return(NULL)
            }
        }

        # Try to parse JSON
        x <- tryCatch(
            content(resp, as = "parsed", type = "application/json"),
            error = function(e) e
        )

        if (inherits(x, "error")) {
            if (attempt < max_attempts) {
                printf("[%s] AID %s attempt %d/%d: JSON parse error: %s — retrying after %ds",
                       ts(), aid, attempt, max_attempts, x$message, base_delay * 2^(attempt - 1))
                Sys.sleep(base_delay * 2^(attempt - 1))
                next
            } else {
                printf("[%s] AID %s failed after %d attempts (JSON parse error: %s)",
                       ts(), aid, max_attempts, x$message)
                return(NULL)
            }
        }

        # If we got here, response + JSON are OK
        break
    }
    # ------------- END RETRY LOOP -------------

    if (is.null(x) || is.null(x$PC_AssayContainer)) {
        printf("[%s] AID %s: no usable JSON after %d attempts", ts(), aid, max_attempts)
        return(NULL)
    }

    container <- x$PC_AssayContainer[[1]]
    assay     <- container$assay
    descr     <- assay$descr

    # -------------------------
    # Basic fields
    # -------------------------
    aid_id      <- descr$aid$id        %||% aid
    version     <- descr$aid$version   %||% NA_character_
    assay_name  <- descr$name          %||% NA_character_

    # -------------------------
    # Source name + ID
    # -------------------------
    source_name <- descr$aid_source$db$name          %||% NA_character_
    source_id   <- descr$aid_source$db$source_id$str %||% NA_character_

    # -------------------------
    # Description, protocol, comment (concatenated)
    # -------------------------
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

    # -------------------------
    # activity_outcome_method (only keep 1 or 2, else NA)
    # -------------------------
    aom <- descr$activity_outcome_method %||% NA_integer_
    if (!(aom %in% c(1, 2))) {
        aom <- NA_integer_
    }

    # -------------------------
    # Targets (name + protein accession)
    # -------------------------
    targets <- descr$target %||% list()

    target_name <- if (length(targets)) {
        paste(
            vapply(
                targets,
                function(t) t$name %||% NA_character_,
                character(1)
            ),
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

    # -------------------------
    # Build row
    # -------------------------
    data.table(
        aid                      = as.character(aid_id),
        version                  = version,
        assay_name               = assay_name,
        source_name              = source_name,
        source_id                = source_id,
        description_combined     = description_combined,
        protocol_combined        = protocol_combined,
        comment_combined         = comment_combined,
        activity_outcome_method  = aom,
        target_name              = target_name,
        target_protein_accession = target_protein_accession
    )
}

# ---- Main ----

printf("[%s] Loading input: %s", ts(), in_csv)
cid_aid <- fread(in_csv)

if (!all(c("pubchem_cid", "bioassay_aid") %in% names(cid_aid))) {
    stop("Input file must have columns: pubchem_cid,bioassay_aid")
}

aids <- unique(cid_aid$bioassay_aid)
aids <- aids[!is.na(aids)]
total <- length(aids)
printf("[%s] Found %d unique AIDs", ts(), total)

rows <- vector("list", total)

# ---- Loop with progress logging + throttling ----
for (i in seq_along(aids)) {

    aid <- aids[i]

    row <- get_assay_row(aid)

    if (is.null(row)) {
        printf("[%s] FAILED %d / %d   (AID: %s)", ts(), i, total, aid)
    } else {
        rows[[i]] <- row
        printf("[%s] Completed %d / %d   (AID: %s)", ts(), i, total, aid)
    }

    # Every 10 queries, sleep 1 second (your original comment said 3s, code was 1s)
    if (i %% 10 == 0) {
        printf("[%s] Processed %d queries — sleeping 1 second...", ts(), i)
        Sys.sleep(1)
    }
}

# Remove NULL results
rows <- rows[!vapply(rows, is.null, logical(1))]

assays <- rbindlist(rows, fill = TRUE)

printf("[%s] Writing output: %s", ts(), out_csv)
fwrite(assays, out_csv)

printf("[%s] Done. Wrote %d rows.", ts(), nrow(assays))
