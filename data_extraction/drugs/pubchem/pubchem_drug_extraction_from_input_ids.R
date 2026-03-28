#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(AnnotationGx)
    library(httr)
    library(jsonlite)
})

# ---------------------------
# small helpers
# ---------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
printf <- function(...) {
    cat(sprintf(...), "\n")
    flush.console()
}

safe_append_fwrite <- function(dt, path) {
    if (is.null(dt) || !nrow(dt)) {
        return(invisible(NULL))
    }
    exists <- file.exists(path)
    fwrite(dt, path, append = exists, col.names = !exists)
}

write_error_log <- function(path, stage, payload, msg) {
    dt <- data.table(
        ts = ts(),
        stage = stage,
        payload = as.character(payload),
        error = as.character(msg)
    )
    safe_append_fwrite(dt, path)
    printf(
        "[%s] ERROR_LOG: wrote %d rows to %s (append=%s)",
        ts(), nrow(dt), path, if (file.exists(path)) "TRUE" else "FALSE"
    )
}

# Parse CLI args like: --in_csv '...' --batch_input_ids '...'
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

# ---------------------------
# CONFIG (match your workflow paths)
# ---------------------------
args <- parse_args()

in_csv <- args$in_csv %||% stop("Missing --in_csv")
batch_input_ids_path <- args$batch_input_ids %||% NULL

out_dir <- "output_data/union/complete"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Outputs expected by downstream scripts
input_id_to_cid_path <- file.path(out_dir, "input_id_to_cid.csv")
union_out_path <- file.path(out_dir, "union_out.csv")
union_synonyms_path <- file.path(out_dir, "union_synonyms.csv")

# NEW (AID portion like older script)
union_bioassays_path <- file.path(out_dir, "union_bioassays.csv")

# Logging / failure bookkeeping
error_log_path <- file.path(out_dir, "pubchem_extraction_errors_input_ids.csv")
failed_input_ids_path <- file.path(out_dir, "failed_input_ids.csv")

# Allowed AIDs (Homo sapiens)
allowed_aids_path <- "homosapien_aids.txt"

# Hard timeout for AnnotationGx HTTP calls (seconds)
CALL_TIMEOUT <- 40

# ---------------------------
# PubChem properties fields (EXACTLY AS IN YOUR CID->BIOASSAYS SCRIPT)
# ---------------------------
properties <- c(
    "Title", "MolecularFormula", "MolecularWeight", "SMILES", "ConnectivitySMILES",
    "InChI", "InChIKey", "IUPACName", "XLogP", "ExactMass", "MonoisotopicMass", "TPSA",
    "Complexity", "Charge", "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount",
    "HeavyAtomCount", "IsotopeAtomCount", "AtomStereoCount", "DefinedAtomStereoCount",
    "UndefinedAtomStereoCount", "BondStereoCount", "DefinedBondStereoCount",
    "UndefinedBondStereoCount", "CovalentUnitCount", "Volume3D", "XStericQuadrupole3D",
    "YStericQuadrupole3D", "ZStericQuadrupole3D", "FeatureCount3D", "FeatureAcceptorCount3D",
    "FeatureDonorCount3D", "FeatureAnionCount3D", "FeatureCationCount3D", "FeatureRingCount3D",
    "FeatureHydrophobeCount3D", "ConformerModelRMSD3D", "EffectiveRotorCount3D",
    "ConformerCount3D", "Fingerprint2D", "PatentCount", "PatentFamilyCount", "LiteratureCount",
    "AnnotationTypes", "AnnotationTypeCount"
)

# ---------------------------
# rename map for final output (EXACTLY AS IN YOUR CID->BIOASSAYS SCRIPT)
# ---------------------------
rename_map <- c(
    CID = "cid", Title = "title", name = "mapped_name", ChEMBL_ID = "molecule_chembl_id",
    MolecularFormula = "molecular_formula", MolecularWeight = "molecular_weight",
    SMILES = "smiles", ConnectivitySMILES = "connectivity_smiles", InChI = "inchi",
    InChIKey = "inchikey", IUPACName = "iupac_name", XLogP = "xlogp", ExactMass = "exact_mass",
    MonoisotopicMass = "monoisotopic_mass", TPSA = "tpsa", Complexity = "complexity",
    Charge = "charge", HBondDonorCount = "h_bond_donor_count",
    HBondAcceptorCount = "h_bond_acceptor_count", RotatableBondCount = "rotatable_bond_count",
    HeavyAtomCount = "heavy_atom_count", IsotopeAtomCount = "isotope_atom_count",
    AtomStereoCount = "atom_stereo_count", DefinedAtomStereoCount = "defined_atom_stereo_count",
    UndefinedAtomStereoCount = "undefined_atom_stereo_count", BondStereoCount = "bond_stereo_count",
    DefinedBondStereoCount = "defined_bond_stereo_count",
    UndefinedBondStereoCount = "undefined_bond_stereo_count",
    CovalentUnitCount = "covalent_unit_count", Volume3D = "volume_3d",
    XStericQuadrupole3D = "x_steric_quadrupole_3d",
    YStericQuadrupole3D = "y_steric_quadrupole_3d",
    ZStericQuadrupole3D = "z_steric_quadrupole_3d", FeatureCount3D = "feature_count_3d",
    FeatureAcceptorCount3D = "feature_acceptor_count_3d",
    FeatureDonorCount3D = "feature_donor_count_3d",
    FeatureAnionCount3D = "feature_anion_count_3d",
    FeatureCationCount3D = "feature_cation_count_3d",
    FeatureRingCount3D = "feature_ring_count_3d",
    FeatureHydrophobeCount3D = "feature_hydrophobe_count_3d",
    ConformerModelRMSD3D = "conformer_model_rmsd_3d",
    EffectiveRotorCount3D = "effective_rotor_count_3d",
    ConformerCount3D = "conformer_count_3d", Fingerprint2D = "fingerprint_2d",
    PatentCount = "patent_count", PatentFamilyCount = "patent_family_count",
    LiteratureCount = "literature_count", AnnotationTypes = "annotation_types",
    AnnotationTypeCount = "annotation_type_count"
)

# ---------------------------
# union_out.csv schema check (from your header list)
# ---------------------------
UNION_OUT_FIELDS <- c(
    "cid", "molecular_formula", "molecular_weight", "smiles", "connectivity_smiles", "inchi",
    "inchikey", "iupac_name", "xlogp", "exact_mass", "monoisotopic_mass", "tpsa", "complexity",
    "charge", "h_bond_donor_count", "h_bond_acceptor_count", "rotatable_bond_count",
    "heavy_atom_count", "isotope_atom_count", "atom_stereo_count", "defined_atom_stereo_count",
    "undefined_atom_stereo_count", "bond_stereo_count", "defined_bond_stereo_count",
    "undefined_bond_stereo_count", "covalent_unit_count", "volume_3d", "x_steric_quadrupole_3d",
    "y_steric_quadrupole_3d", "z_steric_quadrupole_3d", "feature_count_3d",
    "feature_acceptor_count_3d", "feature_donor_count_3d", "feature_anion_count_3d",
    "feature_cation_count_3d", "feature_ring_count_3d", "feature_hydrophobe_count_3d",
    "conformer_model_rmsd_3d", "effective_rotor_count_3d", "conformer_count_3d",
    "fingerprint_2d", "title", "patent_count", "patent_family_count", "literature_count",
    "annotation_types", "annotation_type_count", "molecule_chembl_id", "chembl_max_phase",
    "drug_like", "mapped_name"
)

# ---------------------------
# canonical CID helper (prevents "123.0" / scientific notation mismatches)
# ---------------------------
canon_cid <- function(x) {
    x <- trimws(as.character(x))
    nx <- suppressWarnings(as.numeric(x))
    nx <- as.integer(nx)
    res <- as.character(nx)
    res[is.na(nx)] <- NA_character_
    res
}

append_fwrite_aligned <- function(dt, path, label) {
    if (!nrow(dt)) {
        return(invisible(NULL))
    }
    if (file.exists(path)) {
        hdr <- tryCatch(names(fread(path, nrows = 0)), error = function(e) NULL)
        if (!is.null(hdr) && length(hdr)) {
            missing <- setdiff(hdr, names(dt))
            if (length(missing)) dt[, (missing) := NA]
            dt <- dt[, hdr, with = FALSE]
        }
        fwrite(dt, path, append = TRUE, col.names = FALSE)
        printf("[%s] %s: appended %d rows to %s", ts(), label, nrow(dt), path)
    } else {
        fwrite(dt, path, append = FALSE, col.names = TRUE)
        printf("[%s] %s: wrote %d rows (new file) to %s", ts(), label, nrow(dt), path)
    }
}

merge_mapped_name <- function(old, add) {
    old <- trimws(as.character(old %||% ""))
    add <- trimws(as.character(add %||% ""))

    old_tokens <- unlist(strsplit(old, ";\\s*"))
    add_tokens <- unlist(strsplit(add, ";\\s*"))

    toks <- unique(trimws(c(old_tokens, add_tokens)))
    toks <- toks[nzchar(toks)]

    if (!length(toks)) {
        return(NA_character_)
    }
    paste(toks, collapse = "; ")
}

# ---------------------------
# Fetch synonyms
# ---------------------------
.syn_handle <- httr::handle("https://pubchem.ncbi.nlm.nih.gov")

fetch_synonyms_batch <- function(cids_chr) {
    if (!length(cids_chr)) {
        return(data.table())
    }
    path <- sprintf("/rest/pug/compound/cid/%s/synonyms/JSON", paste(cids_chr, collapse = ","))
    req_fun <- function() {
        resp <- httr::GET(
            handle = .syn_handle, path = path,
            httr::timeout(30),
            httr::user_agent("AnnotationDB/overlap (synonyms_batch)")
        )
        sc <- httr::status_code(resp)
        if (sc >= 500 || sc == 429) stop("server busy")
        resp
    }
    res <- tryCatch(retry_http(req_fun), error = function(e) NULL)

    out <- data.table()

    if (!is.null(res) && httr::status_code(res) == 200) {
        txt <- httr::content(res, as = "text", encoding = "UTF-8")
        j <- try(jsonlite::fromJSON(txt, simplifyVector = TRUE), silent = TRUE)
        if (!inherits(j, "try-error")) {
            info <- j$InformationList$Information
            if (!is.null(info) && length(info)) {
                out <- rbindlist(lapply(info, function(x) {
                    if (is.null(x$Synonym)) {
                        return(NULL)
                    }
                    data.table(
                        synonym = unique(trimws(x$Synonym)),
                        pubchem_cid = as.character(x$CID),
                        source = "PubChem",
                        version = format(
                            as.POSIXct(Sys.time(), tz = "UTC"),
                            "%Y-%m-%dT%H:%M:%SZ"
                        )
                    )
                }), fill = TRUE)
            }
        }
    }

    missing <- setdiff(as.character(cids_chr), unique(out$pubchem_cid))
    if (length(missing)) {
        txt_list <- lapply(missing, function(cid) {
            url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/TXT", cid)
            syn <- tryCatch(unique(readLines(url, warn = FALSE)),
                error = function(e) character(0)
            )
            syn <- unique(trimws(syn[nzchar(syn)]))
            if (!length(syn)) {
                return(NULL)
            }
            data.table(
                synonym = syn,
                pubchem_cid = cid,
                source = "PubChem",
                version = format(
                    as.POSIXct(Sys.time(), tz = "UTC"),
                    "%Y-%m-%dT%H:%M:%SZ"
                )
            )
        })
        out <- rbindlist(c(list(out), txt_list), fill = TRUE)
    }

    if (nrow(out)) {
        out[, key := paste0(tolower(synonym), "::", pubchem_cid)]
        out <- out[!duplicated(key)][, key := NULL]
    }

    out
}

# ---------------------------
# Load allowed AIDs (like older script)
# ---------------------------
aids_dt <- NULL
if (file.exists(allowed_aids_path)) {
    tmp <- tryCatch(fread(allowed_aids_path), error = function(e) NULL)
    if (!is.null(tmp)) {
        # accept either: headerless list, or a column named 'aid'
        if ("aid" %in% names(tmp)) {
            aids_dt <- unique(data.table(aid = as.character(trimws(tmp$aid))))
        } else {
            aids_dt <- unique(data.table(aid = as.character(trimws(tmp[[1]]))))
        }
        aids_dt <- aids_dt[nzchar(aid)]
        data.table::setkey(aids_dt, aid)
        printf("[%s] Loaded %d unique allowed AIDs from %s", ts(), aids_dt[, .N], allowed_aids_path)
    } else {
        printf("[%s] WARNING: could not read AIDs file '%s'. No AID filtering will be applied.", ts(), allowed_aids_path)
    }
} else {
    printf("[%s] WARNING: AIDs file '%s' not found. No AID filtering will be applied.", ts(), allowed_aids_path)
}

# ---------------------------
# Load input mapping (input_id -> cid)
# ---------------------------
printf("[%s] Loading input: %s", ts(), in_csv)
inp <- fread(in_csv)

cn <- names(inp)

# ---------------------------
# InChIKey -> CID resolution helper
# ---------------------------
inchikey_to_cid <- function(key) {
    key <- trimws(as.character(key %||% ""))
    if (!nzchar(key)) {
        return(character(0))
    }

    url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/cids/JSON", key)

    res <- tryCatch(
        {
            resp <- httr::GET(url, httr::timeout(20), httr::user_agent("AnnotationDB/InChIKey2CID"))
            if (httr::status_code(resp) == 200) {
                txt <- httr::content(resp, as = "text", encoding = "UTF-8")
                jsonlite::fromJSON(txt, simplifyVector = TRUE)$IdentifierList$CID
            } else {
                NULL
            }
        },
        error = function(e) NULL
    )

    as.character(res %||% character(0))
}

# ---------------------------
# Load input mapping (input_id -> cid / inchikey)
# ---------------------------
printf("[%s] Loading input: %s", ts(), in_csv)
inp <- fread(in_csv)

cn <- names(inp)

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
    non_cid_key <- setdiff(cn, c(cid_col, "InChiKey", "inchikey", "InChIKey"))
    if (length(non_cid_key) >= 1) {
        input_id_col <- non_cid_key[1]
    } else {
        input_id_col <- cn[1]
    }
}

# InChIKey column guess
inchikey_col <- NULL
for (cand in c("InChiKey", "inchikey", "InChIKey")) {
    if (cand %in% cn) {
        inchikey_col <- cand
        break
    }
}

# Normalize columns
res_inp <- data.table(input_id = as.character(inp[[input_id_col]]))
if (!is.null(cid_col)) res_inp[, cid := as.character(inp[[cid_col]])] else res_inp[, cid := NA_character_]
if (!is.null(inchikey_col)) res_inp[, inchikey := as.character(inp[[inchikey_col]])] else res_inp[, inchikey := NA_character_]

# ---------------------------
# CID RESOLUTION (by InChIKey)
# ---------------------------
needs_res <- is.na(res_inp$cid) | !nzchar(trimws(as.character(res_inp$cid))) | res_inp$cid == "0"
if (any(needs_res) && !is.null(inchikey_col)) {
    keys_to_res <- unique(res_inp[needs_res & !is.na(inchikey) & nzchar(inchikey), inchikey])
    if (length(keys_to_res) > 0) {
        printf("[%s] Resolving %d unique InChIKey(s) to CIDs...", ts(), length(keys_to_res))

        key_map <- list()
        failed_res_keys <- character()
        for (ik in keys_to_res) {
            cids <- inchikey_to_cid(ik)
            if (length(cids) > 0) {
                key_map[[ik]] <- cids[1] # Take first for simplicity
                printf("[%s]   %s -> %s", ts(), ik, cids[1])
            } else {
                printf("[%s]   %s -> NOT FOUND", ts(), ik)
                failed_res_keys <- c(failed_res_keys, ik)
            }
            Sys.sleep(0.2)
        }

        # Log resolution failures to failed_input_ids.csv so coordinator skips them
        if (length(failed_res_keys) > 0) {
            void_dt <- res_inp[inchikey %in% failed_res_keys, .(
                ts = ts(),
                input_id,
                cid = NA_character_,
                stage = "resolution",
                error = "inchikey_not_found"
            )]
            append_fwrite_aligned(void_dt, failed_input_ids_path, "VOID_ENTRIES")
        }

        # Update res_inp with resolved CIDs
        if (length(key_map) > 0) {
            res_inp[needs_res, cid := vapply(inchikey, function(k) key_map[[k]] %||% NA_character_, character(1))]
        }
    }
}

inp <- res_inp[!is.na(input_id) & nzchar(input_id)]
inp[, cid := canon_cid(cid)]
inp <- inp[!is.na(cid) & nzchar(cid)]
printf("[%s] Final resolved input pairs: %d", ts(), nrow(inp))

# ---------------------------
# Apply batch filter (by input_id)
# ---------------------------
batch_ids <- NULL
# The coordinator now passes the batch file as the main --in_csv
# So 'inp' already contains exactly the rows we want to process.
printf("[%s] Processing input batch: %d rows", ts(), nrow(inp))

# (Removed auxiliary mapping file writing)

run_cids <- unique(inp$cid)
run_cids <- run_cids[!is.na(run_cids) & nzchar(run_cids)]

# ---------------------------
# Load union_out (if exists)
# ---------------------------
existing_union <- NULL
already_done_cids <- character()

if (file.exists(union_out_path)) {
    printf("[%s] Checking existing coverage in %s (CID column only)...", ts(), union_out_path)
    # LOW-MEMORY: Only load the CID column to avoid 1GB RAM hit
    already_done_cids <- tryCatch(
        {
            tmp_cids <- fread(union_out_path, select = "cid", colClasses = "character")
            unique(canon_cid(tmp_cids$cid))
        },
        error = function(e) {
            msg <- sprintf("CRITICAL ERROR: Failed to read CIDs from %s: %s", union_out_path, e$message)
            write_error_log(error_log_path, "union_out_cid_read_fail", union_out_path, msg)
            stop(msg)
        }
    )
    printf("[%s] Found %d existing CIDs in master file.", ts(), length(already_done_cids))
    gc() # Clean up after CID read
}

# ---------------------------
# MAPPED_NAME WORKFLOW ADJUSTMENT
# ---------------------------
cid_to_input <- inp[, .(mapped_name_new = paste(unique(input_id), collapse = "; ")), by = cid]

# Now compute props_todo ONLY for CIDs NOT already present
props_todo <- setdiff(run_cids, already_done_cids)
printf("[%s] Coverage for RUN CIDs — props todo: %d", ts(), length(props_todo))

# ---------------------------
# Bulk -> split -> singleton properties fetch
# ---------------------------
fetch_props_split <- function(cids, max_retries = 3, sleep_base = 1) {
    cids <- unique(as.character(cids))
    cids <- cids[nchar(cids) > 0 & !is.na(cids)]

    try_bulk <- function(cids_bulk) {
        cids_num <- suppressWarnings(as.numeric(cids_bulk))
        cids_num <- cids_num[!is.na(cids_num)]
        if (!length(cids_num)) {
            return(NULL)
        }

        for (i in seq_len(max_retries)) {
            out <- tryCatch(
                {
                    setTimeLimit(elapsed = CALL_TIMEOUT, transient = TRUE)
                    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE), add = TRUE)
                    AnnotationGx::mapCID2Properties(ids = cids_num, properties = properties)
                },
                error = function(e) {
                    write_error_log(
                        error_log_path,
                        "mapCID2Properties_error",
                        paste(cids_bulk, collapse = ","),
                        sprintf("attempt %d/%d: %s", i, max_retries, e$message)
                    )
                    NULL
                }
            )

            if (!is.null(out) && NROW(out) > 0) {
                return(as.data.table(out))
            }
            Sys.sleep(sleep_base * (2^(i - 1)) + runif(1, 0, 0.5))
        }

        NULL
    }

    recurse <- function(chunk) {
        res <- try_bulk(chunk)
        if (!is.null(res) && NROW(res) > 0) {
            return(list(ok = res, failed = character()))
        }
        if (length(chunk) <= 1) {
            return(list(ok = NULL, failed = chunk))
        }

        mid <- floor(length(chunk) / 2)
        left <- chunk[seq_len(mid)]
        right <- chunk[(mid + 1):length(chunk)]

        a <- recurse(left)
        b <- recurse(right)

        ok <- NULL
        if (!is.null(a$ok) && NROW(a$ok) > 0) ok <- a$ok
        if (!is.null(b$ok) && NROW(b$ok) > 0) ok <- if (is.null(ok)) b$ok else rbindlist(list(ok, b$ok), fill = TRUE)

        list(ok = ok, failed = c(a$failed, b$failed))
    }

    recurse(cids)
}

props_dt_all <- NULL
failed_cids <- character()

if (length(props_todo) > 0) {
    chunk_size <- 6
    chunks <- split(props_todo, ceiling(seq_along(props_todo) / chunk_size))
    props_list <- list()

    printf("[%s] PROPS total todo: %d | chunk_size=%d | chunks=%d", ts(), length(props_todo), chunk_size, length(chunks))

    for (i in seq_along(chunks)) {
        chunk <- as.character(chunks[[i]])
        printf("[%s] PROPS chunk %d/%d | CIDs: %d | first=%s", ts(), i, length(chunks), length(chunk), chunk[1])

        res <- fetch_props_split(chunk, max_retries = 3, sleep_base = 2)

        if (!is.null(res$ok) && NROW(res$ok) > 0) {
            props_list[[length(props_list) + 1]] <- res$ok
        }

        if (length(res$failed) > 0) {
            failed_cids <- c(failed_cids, res$failed)
            write_error_log(
                error_log_path,
                "mapCID2Properties_chunk_failed",
                paste(res$failed, collapse = ","),
                sprintf("Chunk %d/%d failed after bulk/split/singleton", i, length(chunks))
            )
        }

        Sys.sleep(0.5 + runif(1, 0, 0.5))
    }

    if (length(props_list) > 0) props_dt_all <- rbindlist(props_list, fill = TRUE)
    failed_cids <- unique(as.character(failed_cids))
}

# ---------------------------
# Normalize + rename to union_out.csv schema
# ---------------------------
normalize_props <- function(dt) {
    if (is.null(dt) || NROW(dt) == 0) {
        return(NULL)
    }
    dt <- as.data.table(dt)

    if (!"CID" %in% names(dt)) {
        for (cand in c("cid", "PubChemCID", "pubchem_cid")) {
            if (cand %in% names(dt)) {
                setnames(dt, cand, "CID")
                break
            }
        }
    }
    if (!"CID" %in% names(dt)) setnames(dt, names(dt)[1], "CID")
    dt[, CID := canon_cid(CID)]
    dt
}

props_dt_all <- normalize_props(props_dt_all)

if (!is.null(props_dt_all) && NROW(props_dt_all) > 0) {
    present <- intersect(names(rename_map), names(props_dt_all))
    setnames(props_dt_all, present, rename_map[present])

    # CRITICAL: Merge user input_id(s) into mapped_name for new entries
    if (nrow(cid_to_input) > 0) {
        if ("mapped_name" %in% names(props_dt_all)) {
            props_dt_all <- merge(props_dt_all, cid_to_input, by = "cid", all.x = TRUE, suffixes = c("", ".new"))
            props_dt_all[, mapped_name := vapply(seq_len(.N), function(i) {
                merge_mapped_name(mapped_name[i], mapped_name_new[i])
            }, character(1))]
            props_dt_all[, mapped_name_new := NULL]
        } else {
            props_dt_all <- merge(props_dt_all, cid_to_input, by = "cid", all.x = TRUE)
            setnames(props_dt_all, "mapped_name_new", "mapped_name")
        }
    }
}

# ---------------------------
# Add ChEMBL ID annotation + chembl_max_phase
# ---------------------------
.fetch_chembl_phase_cache <- new.env(parent = emptyenv())

retry_http <- function(fun, times = 2, base_wait = 0.25) {
    for (i in seq_len(times)) {
        res <- try(fun(), silent = TRUE)
        if (!inherits(res, "try-error")) {
            return(res)
        }
        if (i < times) Sys.sleep(base_wait * (1.5^(i - 1)) + runif(1, 0, 0.1))
    }
    stop("retry failed")
}

fetch_chembl_phase <- function(chembl_id) {
    if (is.na(chembl_id) || !nzchar(chembl_id)) {
        return(NA_integer_)
    }
    if (exists(chembl_id, envir = .fetch_chembl_phase_cache)) {
        return(get(chembl_id, envir = .fetch_chembl_phase_cache))
    }

    url <- sprintf("https://www.ebi.ac.uk/chembl/api/data/molecule/%s?format=json", chembl_id)

    resp <- tryCatch(
        retry_http(
            function() {
                r <- httr::GET(
                    url,
                    httr::timeout(20),
                    httr::user_agent("AnnotationDB/overlap (fetch_chembl_phase)")
                )
                sc <- httr::status_code(r)
                if (sc == 429 || sc >= 500) stop(sprintf("Transient HTTP %s from ChEMBL", sc))
                r
            },
            times = 3,
            base_wait = 0.5
        ),
        error = function(e) {
            printf("[%s]   ChEMBL max_phase fetch failed for %s after retries: %s", ts(), chembl_id, e$message)
            NULL
        }
    )

    if (is.null(resp) || httr::status_code(resp) != 200) {
        assign(chembl_id, NA_integer_, envir = .fetch_chembl_phase_cache)
        return(NA_integer_)
    }

    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    j <- try(jsonlite::fromJSON(txt, simplifyVector = TRUE), silent = TRUE)

    phase <- if (!inherits(j, "try-error") && !is.null(j$max_phase)) suppressWarnings(as.integer(j$max_phase)) else NA_integer_

    assign(chembl_id, phase, envir = .fetch_chembl_phase_cache)
    phase
}

annotate_chembl_id_one <- function(cid) {
    cid_num <- suppressWarnings(as.numeric(cid))
    if (is.na(cid_num)) {
        return(NA_character_)
    }

    out <- tryCatch(
        AnnotationGx::annotatePubchemCompound(cid_num, "ChEMBL ID"),
        error = function(e) NULL
    )
    if (is.null(out)) {
        return(NA_character_)
    }

    x <- as.character(out)
    x <- unique(trimws(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) {
        return(NA_character_)
    }
    paste(x, collapse = ";")
}

if (!is.null(props_dt_all) && NROW(props_dt_all) > 0) {
    if (!"molecule_chembl_id" %in% names(props_dt_all)) {
        props_dt_all[, molecule_chembl_id := NA_character_]
    }

    missing_idx <- is.na(props_dt_all$molecule_chembl_id) | !nzchar(trimws(as.character(props_dt_all$molecule_chembl_id)))
    if (any(missing_idx)) {
        cids_to_annotate <- props_dt_all$cid[missing_idx]
        printf("[%s] Annotating ChEMBL ID for %d CID(s)", ts(), length(cids_to_annotate))

        props_dt_all[missing_idx, molecule_chembl_id := vapply(
            cids_to_annotate,
            annotate_chembl_id_one,
            character(1)
        )]
    }

    props_dt_all[, molecule_chembl_id := vapply(
        as.character(molecule_chembl_id),
        function(x) {
            if (is.na(x) || !nzchar(x)) {
                return(NA_character_)
            }
            strsplit(x, ";")[[1]][1]
        },
        character(1)
    )]

    if (!"chembl_max_phase" %in% names(props_dt_all)) {
        props_dt_all[, chembl_max_phase := NA_integer_]
    }
    props_dt_all[, chembl_max_phase := vapply(
        as.character(molecule_chembl_id),
        fetch_chembl_phase,
        integer(1)
    )]

    # ---------------------------
    # Add drug_like flag
    # ---------------------------
    if (!"drug_like" %in% names(props_dt_all)) {
        props_dt_all[, drug_like := !is.na(annotation_types) &
            grepl("(^|\\|)\\s*Drug and Medication Information\\s*(\\||$)",
                annotation_types,
                ignore.case = TRUE
            )]
    }
}

# ---------------------------
# SYNONYM PORTION
# ---------------------------
if (length(run_cids) > 0) {
    printf("[%s] Synonym fetch: processing %d CID(s) from this run", ts(), length(run_cids))

    # Batch them for efficiency
    syn_mini <- split(run_cids, ceiling(seq_along(run_cids) / 10L))

    local_syn <- list()
    for (i in seq_along(syn_mini)) {
        ids <- syn_mini[[i]]
        printf("[%s]   Synonyms batch %d/%d | CIDs: %d", ts(), i, length(syn_mini), length(ids))

        syn_dt <- fetch_synonyms_batch(ids)
        if (nrow(syn_dt)) {
            local_syn[[length(local_syn) + 1]] <- syn_dt
        }
        Sys.sleep(0.5)
    }

    if (length(local_syn)) {
        syn_all_dt <- rbindlist(local_syn, fill = TRUE)

        # Deduplicate globally (within this batch)
        syn_all_dt[, key := paste0(tolower(synonym), "::", pubchem_cid)]
        syn_all_dt <- syn_all_dt[!duplicated(key)][, key := NULL]

        append_fwrite_aligned(syn_all_dt, union_synonyms_path, "Synonyms")
    }
}

# ---------------------------
# AID PORTION
# ---------------------------
fetch_bioassay_summary_one <- function(cid, tries = 2, base_wait = 0.25) {
    cid_chr <- as.character(cid)
    url <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/JSON", cid_chr)

    for (i in seq_len(tries)) {
        resp <- try(httr::GET(
            url,
            httr::timeout(30),
            httr::user_agent("AnnotationDB/overlap (assaysummary)")
        ), silent = TRUE)

        if (inherits(resp, "try-error")) {
            if (i < tries) {
                Sys.sleep(base_wait * (1.5^(i - 1)) + runif(1, 0, 0.1))
                next
            }
            return(data.table())
        }

        if (httr::status_code(resp) != 200) next

        txt <- httr::content(resp, as = "text", encoding = "UTF-8")
        j <- try(jsonlite::fromJSON(txt, simplifyVector = FALSE), silent = TRUE)
        if (inherits(j, "try-error")) {
            return(data.table())
        }

        tab <- j$Table
        if (is.null(tab)) {
            return(data.table())
        }

        col_names <- try(unlist(tab$Columns$Column, use.names = FALSE), silent = TRUE)
        if (inherits(col_names, "try-error")) {
            return(data.table())
        }

        rows <- tab$Row
        if (is.null(rows) || !length(rows)) {
            return(data.table())
        }

        ix_aid <- which(tolower(col_names) %in% tolower("AID"))[1]
        if (is.na(ix_aid)) ix_aid <- which(tolower(col_names) %in% tolower("Bioassay AID"))[1]
        if (is.na(ix_aid)) {
            return(data.table())
        }

        out <- rbindlist(lapply(rows, function(r) {
            vals <- try(unlist(r$Cell, use.names = FALSE), silent = TRUE)
            if (inherits(vals, "try-error")) {
                return(NULL)
            }
            aid <- if (ix_aid > 0 && ix_aid <= length(vals)) as.character(vals[ix_aid]) else NA_character_
            data.table(pubchem_cid = cid_chr, bioassay_aid = aid)
        }), fill = TRUE)

        out <- out[!is.na(bioassay_aid) & nzchar(bioassay_aid)]
        out <- unique(out)
        return(out)
    }

    data.table()
}

bio_all <- data.table()
if (length(run_cids) > 0) {
    printf("[%s] AID fetch: processing %d CID(s) from this run", ts(), length(run_cids))

    bio_mini <- split(
        seq_along(run_cids),
        ceiling(seq_along(run_cids) / 5L)
    )

    local_bio <- list()
    bi <- 0L

    for (idx in bio_mini) {
        bi <- bi + 1L
        ids <- run_cids[idx]

        printf(
            "[%s]   Bioassays mini-batch %d/%d | CIDs: %d",
            ts(), bi, length(bio_mini), length(ids)
        )

        bdt <- fetch_bioassay_for_cids <- function(cids_chr) {
            if (!length(cids_chr)) {
                return(data.table())
            }
            rbindlist(lapply(cids_chr, fetch_bioassay_summary_one), fill = TRUE)
        }
        bdt <- fetch_bioassay_for_cids(ids)
        if (nrow(bdt)) {
            local_bio[[length(local_bio) + 1]] <- bdt
        }

        Sys.sleep(1)
    }

    if (length(local_bio)) {
        bio_all <- rbindlist(local_bio, fill = TRUE)
    }
}

# AID filtering (Homo sapiens only)
if (nrow(bio_all) && !is.null(aids_dt)) {
    before_n <- nrow(bio_all)
    bio_all <- bio_all[aids_dt, on = .(bioassay_aid = aid), nomatch = 0]
    printf("[%s] Filtered bio by Homo sapiens AIDs: %d -> %d rows", ts(), before_n, nrow(bio_all))
}

# Dedup + keep only CID + AID
if (nrow(bio_all)) {
    bio_all[, key := paste0(pubchem_cid, "::", bioassay_aid)]
    bio_all <- bio_all[!duplicated(key)][, key := NULL]
    bio_all <- bio_all[, .(pubchem_cid, bioassay_aid)]
}

# Append to union_bioassays.csv (aligned)
if (nrow(bio_all)) {
    append_fwrite_aligned(bio_all, union_bioassays_path, "Bioassays")
}

# ---------------------------
# Final union_out.csv update (Full Merge/Update)
# ---------------------------
printf("[%s] Finalizing Props for union_out...", ts())

# Load the master file (only 37MB now, safe to load)
master_union <- if (file.exists(union_out_path)) {
    fread(union_out_path, colClasses = "character")
} else {
    data.table()
}

# Group new mappings by CID for efficient merging
inp_summary <- inp[, .(added_names = paste(input_id, collapse = "; ")), by = cid]
inp_summary[, cid := canon_cid(cid)]

if (nrow(master_union) > 0) {
    # 1. Update existing mappings
    master_union[, cid := canon_cid(cid)]

    # Merge names for CIDs already in the master file
    master_union[inp_summary, on = .(cid), mapped_name := mapply(
        merge_mapped_name, mapped_name, added_names
    )]

    # 2. Add brand new CIDs
    if (!is.null(props_dt_all) && nrow(props_dt_all) > 0) {
        # Ensure schema alignment
        for (m in setdiff(UNION_OUT_FIELDS, names(props_dt_all))) props_dt_all[, (m) := NA]
        props_dt_all <- props_dt_all[, UNION_OUT_FIELDS, with = FALSE]
        props_dt_all[, cid := canon_cid(cid)]

        # Only add CIDs NOT already in the master file
        new_cids_dt <- props_dt_all[!cid %in% master_union$cid]
        if (nrow(new_cids_dt) > 0) {
            master_union <- rbind(master_union, new_cids_dt, fill = TRUE)
            printf("[%s] Added %d brand new CIDs to union_out.", ts(), nrow(new_cids_dt))
        }
    }
} else {
    # First time initialization
    if (!is.null(props_dt_all) && nrow(props_dt_all) > 0) {
        for (m in setdiff(UNION_OUT_FIELDS, names(props_dt_all))) props_dt_all[, (m) := NA]
        master_union <- props_dt_all[, UNION_OUT_FIELDS, with = FALSE]
        master_union[, cid := canon_cid(cid)]
    }
}

# Final Save (Overwrite with updated content)
if (nrow(master_union) > 0) {
    # CRITICAL FIX: Prevent exponential quote growth during load-update-save cycles
    for (col in names(master_union)) {
        if (is.character(master_union[[col]])) {
            # Strip double quotes to prevent CSV escaping explosion, but PRESERVE single quotes for IUPAC names
            master_union[, (col) := gsub('"', "", get(col), fixed = TRUE)]
        }
    }

    fwrite(master_union, union_out_path)
    printf("[%s] Successfully updated %s (%d rows total)", ts(), union_out_path, nrow(master_union))
}

# Final cleanup
rm(master_union, props_dt_all)
gc()

# ---------------------------
# Write failed_input_ids.csv
# ---------------------------
if (length(failed_cids) > 0) {
    failed_dt <- inp[cid %in% failed_cids, .(ts = ts(), input_id, cid, stage = "props", error = "mapCID2Properties_failed")]
    if (nrow(failed_dt) > 0) {
        append_fwrite_aligned(failed_dt, failed_input_ids_path, "FAILED_INPUT_IDS")
    }
}

printf("[%s] DONE pubchem_drug_extraction_from_input_ids.", ts())
