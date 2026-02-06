#!/usr/bin/env Rscript

library(httr)
library(jsonlite)

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# 1. Read CIDs
cids <- read.csv("output_data/union/cid_unions.csv", stringsAsFactors = FALSE)$cid
total <- length(cids)

cat(sprintf("[%s] Loaded %d CIDs\n", ts(), total))

get_tox <- function(cid) {
    # Proper route: PUG-View toxicity heading for a compound CID
    url <- sprintf(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON?heading=Toxicity",
        cid
    )

    # Defaults (used if anything fails / missing)
    dili_grade <- NA_character_
    dili_ann   <- NA_character_
    hep_score  <- NA_character_

    resp <- tryCatch(GET(url, timeout(30)), error = function(e) NULL)
    if (is.null(resp) || http_error(resp)) {
        return(data.frame(
            cid                              = cid,
            dili_severity_grade              = dili_grade,
            dili_annotation                  = dili_ann,
            hepatotoxicity_likelihood_score  = hep_score,
            stringsAsFactors = FALSE
        ))
    }

    x <- tryCatch(
        content(resp, as = "parsed", type = "application/json"),
        error = function(e) NULL
    )
    if (is.null(x) || is.null(x$Record) || is.null(x$Record$Section)) {
        return(data.frame(
            cid                              = cid,
            dili_severity_grade              = dili_grade,
            dili_annotation                  = dili_ann,
            hepatotoxicity_likelihood_score  = hep_score,
            stringsAsFactors = FALSE
        ))
    }

    # Record$Section[[1]] -> "Toxicity"
    s1 <- x$Record$Section[[1]]
    if (is.null(s1$Section) || length(s1$Section) == 0) {
        return(data.frame(
            cid                              = cid,
            dili_severity_grade              = dili_grade,
            dili_annotation                  = dili_ann,
            hepatotoxicity_likelihood_score  = hep_score,
            stringsAsFactors = FALSE
        ))
    }

    # "Toxicological Information"
    tox_info <- s1$Section[[1]]
    subsecs  <- tox_info$Section
    if (is.null(subsecs)) {
        return(data.frame(
            cid                              = cid,
            dili_severity_grade              = dili_grade,
            dili_annotation                  = dili_ann,
            hepatotoxicity_likelihood_score  = hep_score,
            stringsAsFactors = FALSE
        ))
    }

    dili_sec <- NULL
    hep_sec  <- NULL

    # Find DILI + Hepatotoxicity sections
    for (sec in subsecs) {
        if (!is.null(sec$TOCHeading) && sec$TOCHeading == "Drug Induced Liver Injury") {
            dili_sec <- sec
        }
        if (!is.null(sec$TOCHeading) && sec$TOCHeading == "Hepatotoxicity") {
            hep_sec <- sec
        }
    }

    # DILI section: Severity Grade + DILI Annotation
    if (!is.null(dili_sec) && !is.null(dili_sec$Information)) {
        for (info in dili_sec$Information) {
            name <- info$Name
            swm  <- info$Value$StringWithMarkup
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

    # Hepatotoxicity section: line containing "Likelihood score:"
    if (!is.null(hep_sec) && !is.null(hep_sec$Information)) {
        swm_list <- hep_sec$Information[[1]]$Value$StringWithMarkup
        if (!is.null(swm_list)) {
            for (swm in swm_list) {
                if (!is.null(swm$String) &&
                    grepl("Likelihood score:", swm$String, fixed = TRUE)) {
                    hep_score <- swm$String
                    break
                }
            }
        }
    }

    data.frame(
        cid                              = cid,
        dili_severity_grade              = dili_grade,
        dili_annotation                  = dili_ann,
        hepatotoxicity_likelihood_score  = hep_score,
        stringsAsFactors = FALSE
    )
}

# 2. Apply to all CIDs with simple progress logs
rows <- vector("list", total)

for (i in seq_along(cids)) {
    cid <- cids[i]
    cat(sprintf("[%s] (%d/%d) CID %s\n", ts(), i, total, cid))
    rows[[i]] <- get_tox(cid)
}

out <- do.call(rbind, rows)

# 3. Write CSV (empty cells for missing fields)
out_file <- "output_data/union/complete/toxicity_output.csv"
write.csv(out, out_file, row.names = FALSE, na = "")

cat(sprintf("[%s] Wrote %d rows to %s\n", ts(), nrow(out), out_file))
