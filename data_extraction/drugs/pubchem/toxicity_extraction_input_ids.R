#!/usr/bin/env Rscript

library(httr)
library(jsonlite)
library(data.table)

`%||%` <- function(a, b) if (!is.null(a)) a else b
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

prop_out <- "output_data/union/complete/union_out.csv"
out_file <- "output_data/union/complete/toxicity_output.csv"
err_log <- "output_data/union/complete/toxicity_errors.csv"

batch_input_ids_csv <- NULL
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2 && args[1] == "--batch_input_ids") {
    batch_input_ids_csv <- args[2]
}

append_csv <- function(df, path) {
    if (is.null(df) || nrow(df) == 0) {
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
    row <- data.frame(
        time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        pubchem_cid = as.character(pubchem_cid),
        detail = as.character(detail),
        stringsAsFactors = FALSE
    )
    append_csv(row, err_log)
}

get_done_pubchem_cids <- function(path) {
    if (!file.exists(path)) {
        return(character(0))
    }

    dt <- tryCatch(fread(path, fill = TRUE), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) {
        return(character(0))
    }

    if ("pubchem_cid" %in% names(dt)) {
        done <- unique(trimws(as.character(dt$pubchem_cid)))
    } else if ("cid" %in% names(dt)) {
        done <- unique(trimws(as.character(dt$cid)))
    } else {
        done <- character(0)
    }

    done <- done[nzchar(done)]
    done <- done[done != "0"]
    done
}

# For consistency with your original pipeline requirements
if (!file.exists(prop_out)) stop("union_out.csv missing; run PubChem props first.")
# union_out can contain occasional malformed rows; fill=TRUE prevents fread from stopping early
props <- fread(prop_out, fill = TRUE)
if (!"cid" %in% names(props)) stop("union_out.csv must have 'cid' column")
props[, cid := trimws(as.character(cid))]
props <- props[nzchar(cid) & cid != "0"]

# ---------------------------
# NEW: Determine batch CIDs using input_id_to_cid mapping
# ---------------------------
get_batch_cids_from_input_ids <- function(props_dt, batch_input_ids_csv) {
    # props_dt is the data.table loaded from union_out.csv

    # If no batch file provided, behave like "run everything"
    if (is.null(batch_input_ids_csv) || !nzchar(batch_input_ids_csv)) {
        out <- unique(trimws(as.character(props_dt$cid)))
        out <- out[nzchar(out) & out != "0"]
        return(out)
    }

    if (!file.exists(batch_input_ids_csv)) stop("batch_input_ids file not found: ", batch_input_ids_csv)

    b <- fread(batch_input_ids_csv, fill = TRUE)
    # Handle full row format
    cn <- names(b)
    input_id_col <- NULL
    for (cand in c("input_id", "inputId", "substance", "name", "drug", "synonym")) {
        if (cand %in% cn) {
            input_id_col <- cand
            break
        }
    }
    if (is.null(input_id_col)) input_id_col <- cn[1]

    batch_ids <- unique(trimws(as.character(b[[input_id_col]])))
    batch_ids <- batch_ids[nzchar(batch_ids)]
    if (!length(batch_ids)) {
        return(character(0))
    }

    # Create an expanded mapping: CID -> Single JCP ID from the props_dt
    # We use the already loaded props_dt to avoid re-reading the file
    u_expanded <- props_dt[, .(input_id = unlist(strsplit(mapped_name, ";\\s*"))), by = cid]
    u_expanded[, input_id := trimws(input_id)]

    # map batch input_ids -> cids
    batch_cids <- unique(u_expanded[input_id %in% batch_ids, cid])
    batch_cids <- batch_cids[nzchar(batch_cids) & batch_cids != "0"]

    batch_cids
}

all_cids <- get_batch_cids_from_input_ids(props, batch_input_ids_csv)
total_all <- length(all_cids)
cat(sprintf("[%s] Batch CIDs: %d\n", ts(), total_all))

done <- get_done_pubchem_cids(out_file)
todo <- setdiff(all_cids, done)
total <- length(todo)

cat(sprintf(
    "[%s] Toxicity — total batch CIDs: %d | already done: %d | remaining: %d\n",
    ts(), total_all, length(done), total
))

get_tox <- function(pubchem_cid) {
    url <- sprintf(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON?heading=Toxicity",
        pubchem_cid
    )

    dili_grade <- NA_character_
    dili_ann <- NA_character_
    hep_score <- NA_character_

    resp <- tryCatch(GET(url, timeout(30)), error = function(e) e)
    if (inherits(resp, "error")) {
        log_error(pubchem_cid, resp$message)
        return(data.frame(
            pubchem_cid = pubchem_cid,
            dili_severity_grade = dili_grade,
            dili_annotation = dili_ann,
            hepatotoxicity_likelihood_score = hep_score,
            stringsAsFactors = FALSE
        ))
    }

    if (http_error(resp)) {
        log_error(pubchem_cid, paste("HTTP", status_code(resp)))
        return(data.frame(
            pubchem_cid = pubchem_cid,
            dili_severity_grade = dili_grade,
            dili_annotation = dili_ann,
            hepatotoxicity_likelihood_score = hep_score,
            stringsAsFactors = FALSE
        ))
    }

    x <- tryCatch(
        content(resp, as = "parsed", type = "application/json"),
        error = function(e) e
    )

    if (inherits(x, "error") || is.null(x) || is.null(x$Record) || is.null(x$Record$Section)) {
        log_error(pubchem_cid, "Missing Record/Section")
        return(data.frame(
            pubchem_cid = pubchem_cid,
            dili_severity_grade = dili_grade,
            dili_annotation = dili_ann,
            hepatotoxicity_likelihood_score = hep_score,
            stringsAsFactors = FALSE
        ))
    }

    s1 <- x$Record$Section[[1]]
    tox_info <- s1$Section[[1]] %||% NULL
    subsecs <- tox_info$Section %||% NULL

    if (is.null(subsecs)) {
        return(data.frame(
            pubchem_cid = pubchem_cid,
            dili_severity_grade = dili_grade,
            dili_annotation = dili_ann,
            hepatotoxicity_likelihood_score = hep_score,
            stringsAsFactors = FALSE
        ))
    }

    dili_sec <- NULL
    hep_sec <- NULL

    for (sec in subsecs) {
        if (!is.null(sec$TOCHeading) && sec$TOCHeading == "Drug Induced Liver Injury") {
            dili_sec <- sec
        }
        if (!is.null(sec$TOCHeading) && sec$TOCHeading == "Hepatotoxicity") {
            hep_sec <- sec
        }
    }

    if (!is.null(dili_sec) && !is.null(dili_sec$Information)) {
        for (info in dili_sec$Information) {
            name <- info$Name
            swm <- info$Value$StringWithMarkup
            if (!is.null(name) && !is.null(swm) && length(swm) > 0) {
                val <- swm[[1]]$String
                if (name == "Severity Grade") {
                    dili_grade <- val
                } else if (name == "DILI Annotation") {
                    dili_ann <- val
                }
            }
        }
    }

    if (!is.null(hep_sec) && !is.null(hep_sec$Information)) {
        swm_list <- hep_sec$Information[[1]]$Value$StringWithMarkup
        if (!is.null(swm_list)) {
            for (swm in swm_list) {
                if (!is.null(swm$String) && grepl("Likelihood score:", swm$String, fixed = TRUE)) {
                    hep_score <- swm$String
                    break
                }
            }
        }
    }

    data.frame(
        pubchem_cid = pubchem_cid,
        dili_severity_grade = dili_grade,
        dili_annotation = dili_ann,
        hepatotoxicity_likelihood_score = hep_score,
        stringsAsFactors = FALSE
    )
}

if (!length(todo)) {
    cat(sprintf("[%s] Nothing to do. Exiting.\n", ts()))
    quit(save = "no", status = 0)
}

batch_size <- 25L
sleep_s <- 0.25

batches <- split(todo, ceiling(seq_along(todo) / batch_size))

for (b in seq_along(batches)) {
    batch <- batches[[b]]
    cat(sprintf("[%s] Toxicity batch %d/%d | n=%d\n", ts(), b, length(batches), length(batch)))

    rows <- vector("list", length(batch))

    for (i in seq_along(batch)) {
        cid <- batch[i]
        cat(sprintf("[%s] (%d/%d) CID %s\n", ts(), i, length(batch), cid))

        row <- tryCatch(get_tox(cid), error = function(e) {
            log_error(cid, e$message)
            data.frame(
                pubchem_cid = cid,
                dili_severity_grade = NA_character_,
                dili_annotation = NA_character_,
                hepatotoxicity_likelihood_score = NA_character_,
                stringsAsFactors = FALSE
            )
        })

        rows[[i]] <- row
        Sys.sleep(sleep_s)
    }

    out <- do.call(rbind, rows)
    append_csv(out, out_file)

    cat(sprintf("[%s] Wrote %d rows to %s (append mode)\n", ts(), nrow(out), out_file))
}

cat(sprintf("[%s] DONE toxicity.\n", ts()))
