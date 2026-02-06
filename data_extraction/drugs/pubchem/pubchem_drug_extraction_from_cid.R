# From CID to bioassays

library(AnnotationGx)
library(jsonlite)
library(data.table)
library(future.apply)
library(httr)

`%||%` <- function(a, b) if (!is.null(a)) a else b
printf <- function(...) { cat(sprintf(...), "\n"); flush.console() }
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
tic <- function() proc.time()[["elapsed"]]

in_csv   <- "output_data/overlap/cid_overlap.csv"
out_dir  <- "output_data/overlap/complete"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
prop_out <- file.path(out_dir, "overlap_out.csv")
syn_out  <- file.path(out_dir, "overlap_synonyms.csv")
bio_out  <- file.path(out_dir, "overlap_bioassays.csv")

# Allowed AIDs
aids_file <- "homosapien_aids.txt"
aids_dt <- tryCatch(
    fread(aids_file),
    error = function(e) {
        printf("[%s] WARNING: could not read AIDs file '%s': %s",
               ts(), aids_file, e$message)
        NULL
    }
)

if (!is.null(aids_dt) && "aid" %in% names(aids_dt)) {
    aids_dt <- unique(data.table(aid = as.character(trimws(aids_dt$aid))))
    aids_dt <- aids_dt[nzchar(aid)]
    data.table::setkey(aids_dt, aid)
    printf("[%s] Loaded %d unique allowed AIDs from %s",
           ts(), aids_dt[, .N], aids_file)
} else {
    printf("[%s] WARNING: AIDs file missing or malformed. No AID filtering will be applied.",
           ts())
    aids_dt <- NULL
}

printf("[%s] Starting run", ts())
run_start <- tic()

# Read CIDs
printf("[%s] Loading input CSV: %s", ts(), in_csv)
drugs <- read.csv(in_csv, stringsAsFactors = FALSE)

if (!"cid" %in% names(drugs)) {
    stop(sprintf("Input must contain a 'cid' column. Columns found: %s",
                 paste(names(drugs), collapse = ", ")))
}

# Create cid_map from the cid column
cid_vec <- unique(na.omit(trimws(as.character(drugs$cid))))
cid_vec <- cid_vec[nzchar(cid_vec)]

cid_map <- data.frame(
    cid = cid_vec,
    stringsAsFactors = FALSE
)

if ("drug" %in% names(drugs)) {
    # Use the first drug name associated with each CID
    drugs$cid <- trimws(as.character(drugs$cid))
    drugs$drug <- trimws(as.character(drugs$drug))
    name_map <- unique(drugs[!is.na(drugs$cid) & nzchar(drugs$cid),
                              c("cid", "drug")])
    cid_map$name <- name_map$drug[match(cid_map$cid, name_map$cid)]
}

printf("[%s] Loaded %d unique CIDs", ts(), nrow(cid_map))

# PubChem properties fields
properties <- c(
    "Title","MolecularFormula","MolecularWeight","SMILES","ConnectivitySMILES",
    "InChI","InChIKey","IUPACName","XLogP","ExactMass","MonoisotopicMass","TPSA",
    "Complexity","Charge","HBondDonorCount","HBondAcceptorCount","RotatableBondCount",
    "HeavyAtomCount","IsotopeAtomCount","AtomStereoCount","DefinedAtomStereoCount",
    "UndefinedAtomStereoCount","BondStereoCount","DefinedBondStereoCount",
    "UndefinedBondStereoCount","CovalentUnitCount","Volume3D","XStericQuadrupole3D",
    "YStericQuadrupole3D","ZStericQuadrupole3D","FeatureCount3D","FeatureAcceptorCount3D",
    "FeatureDonorCount3D","FeatureAnionCount3D","FeatureCationCount3D","FeatureRingCount3D",
    "FeatureHydrophobeCount3D","ConformerModelRMSD3D","EffectiveRotorCount3D",
    "ConformerCount3D","Fingerprint2D","PatentCount","PatentFamilyCount","LiteratureCount",
    "AnnotationTypes","AnnotationTypeCount"
)

# rename map for final output
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

.syn_handle <- httr::handle("https://pubchem.ncbi.nlm.nih.gov")

retry <- function(fun, times = 2, base_wait = 0.25) {
    for (i in seq_len(times)) {
        res <- try(fun(), silent = TRUE)
        if (!inherits(res, "try-error")) return(res)
        if (i < times)
            Sys.sleep(base_wait * (1.5^(i - 1)) + runif(1,0,0.1))
    }
    stop("retry failed")
}

# Fetch synonyms
fetch_synonyms_batch <- function(cids_chr) {
    if (!length(cids_chr)) return(data.table())
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
    res <- tryCatch(retry(req_fun), error = function(e) NULL)

    out <- data.table()

    if (!is.null(res) && httr::status_code(res) == 200) {
        txt <- httr::content(res, as = "text", encoding = "UTF-8")
        j <- try(jsonlite::fromJSON(txt, simplifyVector = TRUE), silent = TRUE)
        if (!inherits(j, "try-error")) {
            info <- j$InformationList$Information
            if (!is.null(info) && length(info)) {
                out <- rbindlist(lapply(info, function(x) {
                    if (is.null(x$Synonym)) return(NULL)
                    data.table(
                        synonym     = unique(trimws(x$Synonym)),
                        pubchem_cid = as.character(x$CID),
                        source      = "PubChem",
                        version     = format(as.POSIXct(Sys.time(), tz = "UTC"),
                                             "%Y-%m-%dT%H:%M:%SZ")
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
                            error = function(e) character(0))
            syn <- unique(trimws(syn[nzchar(syn)]))
            if (!length(syn)) return(NULL)
            data.table(
                synonym = syn,
                pubchem_cid = cid,
                source = "PubChem",
                version = format(as.POSIXct(Sys.time(), tz = "UTC"),
                                 "%Y-%m-%dT%H:%M:%SZ")
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

# ChEMBL max phase
.fetch_chembl_phase_cache <- new.env(parent = emptyenv())

fetch_chembl_phase <- function(chembl_id) {
    if (is.na(chembl_id) || !nzchar(chembl_id)) return(NA_integer_)

    # return cached if we already looked this one up
    if (exists(chembl_id, envir = .fetch_chembl_phase_cache)) {
        return(get(chembl_id, envir = .fetch_chembl_phase_cache))
    }

    url <- sprintf(
        "https://www.ebi.ac.uk/chembl/api/data/molecule/%s?format=json",
        chembl_id
    )

    # use your retry() helper only around the ChEMBL HTTP call
    resp <- tryCatch(
        retry(
            function() {
                r <- httr::GET(
                    url,
                    httr::timeout(20),
                    httr::user_agent("AnnotationDB/overlap (fetch_chembl_phase)")
                )
                sc <- httr::status_code(r)

                # treat 429 + 5xx as transient and trigger retry()
                if (sc == 429 || sc >= 500) {
                    stop(sprintf("Transient HTTP %s from ChEMBL", sc))
                }

                r
            },
            times     = 3,      # how many attempts
            base_wait = 0.5     # base backoff in seconds
        ),
        error = function(e) {
            printf(
                "[%s]   ChEMBL max_phase fetch failed for %s after retries: %s",
                ts(), chembl_id, e$message
            )
            NULL
        }
    )

    # if still no good after retries, cache NA and bail
    if (is.null(resp) || httr::status_code(resp) != 200) {
        assign(chembl_id, NA_integer_, envir = .fetch_chembl_phase_cache)
        return(NA_integer_)
    }

    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    j   <- try(jsonlite::fromJSON(txt, simplifyVector = TRUE), silent = TRUE)

    phase <- if (!inherits(j, "try-error") && !is.null(j$max_phase)) {
        suppressWarnings(as.integer(j$max_phase))
    } else {
        NA_integer_
    }

    assign(chembl_id, phase, envir = .fetch_chembl_phase_cache)
    phase
}

# Bioassay summary (AIDs only, no lookup)
fetch_bioassay_summary_one <- function(cid, tries = 2, base_wait = 0.25) {
    cid_chr <- as.character(cid)
    url <- sprintf(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/JSON",
        cid_chr
    )

    for (i in seq_len(tries)) {
        resp <- try(httr::GET(
            url,
            httr::timeout(30),
            httr::user_agent("AnnotationDB/overlap (assaysummary)")
        ), silent = TRUE)

        if (inherits(resp, "try-error")) {
            if (i < tries) {
                Sys.sleep(base_wait * (1.5^(i - 1)) + runif(1,0,0.1))
                next
            }
            return(data.table())
        }

        if (httr::status_code(resp) != 200) next

        txt <- httr::content(resp, as = "text", encoding = "UTF-8")
        j <- try(jsonlite::fromJSON(txt, simplifyVector = FALSE), silent = TRUE)
        if (inherits(j, "try-error")) return(data.table())

        tab <- j$Table
        if (is.null(tab)) return(data.table())

        col_names <- try(unlist(tab$Columns$Column, use.names = FALSE),
                         silent = TRUE)
        if (inherits(col_names, "try-error")) return(data.table())

        rows <- tab$Row
        if (is.null(rows) || !length(rows)) return(data.table())

        ix_aid <- which(tolower(col_names) %in% tolower("AID"))[1]
        if (is.na(ix_aid)) ix_aid <- which(tolower(col_names) %in%
                                           tolower("Bioassay AID"))[1]
        if (is.na(ix_aid)) return(data.table())

        out <- rbindlist(lapply(rows, function(r) {
            vals <- try(unlist(r$Cell, use.names = FALSE), silent = TRUE)
            if (inherits(vals, "try-error")) return(NULL)
            aid <- if (ix_aid > 0 && ix_aid <= length(vals))
                as.character(vals[ix_aid]) else NA_character_
            data.table(pubchem_cid = cid_chr, bioassay_aid = aid)
        }), fill = TRUE)

        out <- out[!is.na(bioassay_aid) & nzchar(bioassay_aid)]
        out <- unique(out)

        return(out)
    }

    data.table()
}

fetch_bioassay_for_cids <- function(cids_chr) {
    if (!length(cids_chr)) return(data.table())
    rbindlist(lapply(cids_chr, fetch_bioassay_summary_one), fill = TRUE)
}

# MAIN (single pass, no chunking)
printf("[%s] Processing %d CIDs (single pass, no chunking)",
       ts(), nrow(cid_map))

props_all <- list()
syn_all   <- list()
bio_all   <- list()

chunk <- cid_map
printf("[%s]   Single chunk | CIDs: %d", ts(), nrow(chunk))

#---- PROPERTIES----
prop_rows <- 0L
local_props <- list()

cids_chr_chunk <- as.character(chunk$cid)
mini_chunks <- split(seq_along(cids_chr_chunk),
                     ceiling(seq_along(cids_chr_chunk) / 10L))

mini_i <- 0L
for (mc in mini_chunks) {
    mini_i <- mini_i + 1L
    sub_chr <- cids_chr_chunk[mc]
    sub_num <- suppressWarnings(as.numeric(sub_chr))

    printf("[%s]   Properties mini-batch %d/%d | CIDs: %d",
           ts(), mini_i, length(mini_chunks), length(sub_chr))

    prop <- tryCatch(
        AnnotationGx::mapCID2Properties(ids = sub_num, properties = properties),
        error = function(e) {
            printf("[%s]   mapCID2Properties error: %s", ts(), e$message)
            NULL
        }
    )

    if (!is.null(prop) && nrow(prop)) {
        prop$CID <- as.character(prop$CID)

        # Only set 'name' if we have it in chunk
        if ("name" %in% names(chunk)) {
            prop$name <- chunk$name[match(prop$CID, chunk$cid)]
        }

        chembl <- tryCatch(
            AnnotationGx::annotatePubchemCompound(as.numeric(prop$CID), "ChEMBL ID", nParallel = 2),
            error = function(e) NA_character_
        )

        if (is.atomic(chembl) && length(chembl) %in% c(1, nrow(prop))) {
            prop$ChEMBL_ID <- chembl
        } else {
            prop$ChEMBL_ID <- vapply(prop$CID, function(x)
                tryCatch(
                    AnnotationGx::annotatePubchemCompound(as.numeric(x), "ChEMBL ID"),
                    error = function(e) NA_character_
                ), character(1))
        }

		prop$ChEMBL_ID <- vapply(
			prop$ChEMBL_ID,
			function(x) {
				if (is.na(x) || !nzchar(x)) return(NA_character_)
				strsplit(x, ";")[[1]][1]   # take the first
			},
			character(1)
		)

        prop$chembl_max_phase <- vapply(
            as.character(prop$ChEMBL_ID),
            fetch_chembl_phase,
            integer(1)
        )

        local_props[[length(local_props) + 1]] <- as.data.table(prop)
    }

    Sys.sleep(5)
}

if (length(local_props)) {
    props_dt <- rbindlist(local_props, fill = TRUE)
    present <- intersect(names(rename_map), names(props_dt))
    setnames(props_dt, present, rename_map[present])

    props_dt[, drug_like := !is.na(annotation_types) &
                             grepl("(^|\\|)\\s*Drug and Medication Information\\s*(\\||$)",
                                   annotation_types, ignore.case = TRUE)]

    prop_rows <- nrow(props_dt)
    props_all[[1]] <- props_dt
    printf("[%s]   Properties rows: %d", ts(), prop_rows)
}

#---- SYNONYMS----
syn_rows <- 0L
local_syn <- list()

syn_cids <- unique(chunk$cid)
syn_mini <- split(seq_along(syn_cids),
                  ceiling(seq_along(syn_cids) / 10L))

syn_i <- 0L
for (sc in syn_mini) {
    syn_i <- syn_i + 1L
    ids <- syn_cids[sc]

    printf("[%s]   Synonyms mini-batch %d/%d | CIDs: %d",
           ts(), syn_i, length(syn_mini), length(ids))

    syn_dt <- fetch_synonyms_batch(ids)
    if (nrow(syn_dt)) {
        local_syn[[length(local_syn) + 1]] <- syn_dt
        syn_rows <- syn_rows + nrow(syn_dt)
    }

    Sys.sleep(5)
}

if (length(local_syn)) {
    syn_dt <- rbindlist(local_syn, fill = TRUE)
    syn_dt[, key := paste0(tolower(synonym), "::", pubchem_cid)]
    syn_dt <- syn_dt[!duplicated(key)][, key := NULL]
    syn_all[[1]] <- syn_dt
    printf("[%s]   Synonyms rows: %d", ts(), syn_rows)
}

# BIOASSAYS
bio_rows <- 0L
local_bio <- list()

bio_cids <- unique(chunk$cid)
bio_mini <- split(seq_along(bio_cids),
                  ceiling(seq_along(bio_cids) / 10L))

bio_i <- 0L
for (bi in bio_mini) {
    bio_i <- bio_i + 1L
    ids <- bio_cids[bi]

    printf("[%s]   Bioassays mini-batch %d/%d | CIDs: %d",
           ts(), bio_i, length(bio_mini), length(ids))

    bdt <- fetch_bioassay_for_cids(ids)
    if (nrow(bdt)) {
        local_bio[[length(local_bio) + 1]] <- bdt
        bio_rows <- bio_rows + nrow(bdt)
    }

    Sys.sleep(5)
}

if (length(local_bio)) {
    bio_dt <- rbindlist(local_bio, fill = TRUE)
    bio_all[[1]] <- bio_dt
    printf("[%s]   Bioassays rows: %d", ts(), bio_rows)
}

# Combine all
props_all <- if (length(props_all)) rbindlist(props_all, fill = TRUE) else data.table()
syn_all   <- if (length(syn_all))   rbindlist(syn_all,   fill = TRUE) else data.table()
bio_all   <- if (length(bio_all))   rbindlist(bio_all,   fill = TRUE) else data.table()

printf("[%s] Combined totals — props: %d | syns: %d | bio: %d",
       ts(), nrow(props_all), nrow(syn_all), nrow(bio_all))

# AID filtering (Homo sapiens only)
if (nrow(bio_all) && !is.null(aids_dt)) {
    before_n <- nrow(bio_all)
    bio_all <- bio_all[aids_dt, on = .(bioassay_aid = aid), nomatch = 0]
    printf("[%s] Filtered bio by Homo sapiens AIDs: %d -> %d rows",
           ts(), before_n, nrow(bio_all))
}

# Dedup
if (nrow(props_all)) props_all <- unique(props_all)
if (nrow(syn_all)) {
    syn_all[, key := paste0(tolower(synonym), "::", pubchem_cid)]
    syn_all <- syn_all[!duplicated(key)][, key := NULL]
}
if (nrow(bio_all)) {
    bio_all[, key := paste0(pubchem_cid, "::", bioassay_aid)]
    bio_all <- bio_all[!duplicated(key)][, key := NULL]
}

# Keep only CID + AID in bioassays
if (nrow(bio_all)) {
    bio_all <- bio_all[, .(pubchem_cid, bioassay_aid)]
}

# Append
printf("[%s] Writing outputs (append mode) ...", ts())

append_fwrite <- function(dt, path, label) {
    if (!nrow(dt)) {
        printf("[%s] %s: no rows to write, skipping.", ts(), label)
        return(invisible(NULL))
    }
    file_exists <- file.exists(path)
    fwrite(
        dt,
        path,
        append    = file_exists,
        col.names = !file_exists
    )
    printf("[%s] %s: wrote %d rows to %s (append = %s)",
           ts(), label, nrow(dt), path, if (file_exists) "TRUE" else "FALSE")
}

append_fwrite(props_all, prop_out, "Properties")
append_fwrite(syn_all,   syn_out,  "Synonyms")
append_fwrite(bio_all,   bio_out,  "Bioassays")

printf("[%s] Done in %.1fs", ts(), tic() - run_start)
