#!/usr/bin/env Rscript

library(data.table)

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
printf <- function(...) { cat(sprintf(...), "\n"); flush.console() }

# ---------------------------
# CONFIGURATION
# ---------------------------
# Edit these values to change the default behavior
master_smiles  <- "output_data/union/smiles_union.csv"
out_dir        <- "output_data/union/complete"
batch_n        <- 50L

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
if (!is.null(args$in_csv))     master_smiles <- args$in_csv
if (!is.null(args$batch_size)) batch_n       <- as.integer(args$batch_size)
if (!is.null(args$out_dir))    out_dir       <- args$out_dir

# ---------------------------
# SETUP
# ---------------------------
if (!file.exists(master_smiles)) {
  stop(sprintf("Master SMILES file missing: %s. Please edit the 'master_smiles' variable in the script or provide it via --in_csv.", master_smiles))
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

prop_out <- file.path(out_dir, "union_out.csv")
batch_smiles_csv <- file.path(out_dir, "batch_smiles.csv")

script_smiles  <- "pubchem_drug_extraction_from_SMILES.R"
script_assays  <- "assay_extraction.R"
script_tox     <- "toxicity_extraction.R"

# ---------------------------
# HELPERS
# ---------------------------
run <- function(cmd, stop_on_error = TRUE) {
	printf("[%s] RUN: %s", ts(), cmd)
	code <- system(cmd)
	if (code != 0) {
		msg <- sprintf("Command failed with exit code %s: %s", code, cmd)
		if (stop_on_error) stop(msg)
		else printf("[%s] WARNING: %s", ts(), msg)
	}
	code
}

load_master_smiles <- function(path) {
	if (!file.exists(path)) stop(sprintf("Master SMILES file missing: %s", path))
	inp <- fread(path, select = "SMILES")
	all_smiles <- unique(trimws(as.character(inp$SMILES)))
	all_smiles <- all_smiles[nzchar(all_smiles)]
	all_smiles
}

# Build "seen input SMILES" from union_out.csv:
# - mapped_name tokens (semicolon-separated input SMILES)
# - smiles field (PubChem canonical) as a fallback
# NEW: Cross-check with toxicity and assays to ensure item is TRULY complete.
load_seen_smiles <- function(out_dir) {
	prop_path <- file.path(out_dir, "union_out.csv")
	tox_path  <- file.path(out_dir, "toxicity_output.csv")
	
	if (!file.exists(prop_path)) return(character(0))

	# Load properties (only columns we need)
	props <- tryCatch(fread(prop_path, select = c("cid", "smiles", "mapped_name"), fill = TRUE), error = function(e) NULL)
	if (is.null(props) || !nrow(props)) return(character(0))

	# Identify CIDs that are truly finished (must be in toxicity output if it exists)
	# Assays are harder because not all CIDs have assays, so we primarily use Toxicity as the 'completion' marker
	finished_cids <- NULL
	if (file.exists(tox_path)) {
		tox <- tryCatch(fread(tox_path, select = "pubchem_cid", fill = TRUE), error = function(e) NULL)
		if (!is.null(tox)) finished_cids <- unique(trimws(as.character(tox$pubchem_cid)))
	}

	# If we have a list of finished CIDs, filter the properties
	if (!is.null(finished_cids)) {
		props <- props[trimws(as.character(cid)) %in% finished_cids]
	}

	seen <- character(0)

	if ("mapped_name" %in% names(props)) {
		mn <- trimws(as.character(props$mapped_name))
		mn <- mn[nzchar(mn)]
		if (length(mn)) {
			toks <- unlist(strsplit(mn, ";", fixed = TRUE), use.names = FALSE)
			toks <- trimws(as.character(toks))
			toks <- toks[nzchar(toks)]
			seen <- c(seen, toks)
		}
	}

	if ("smiles" %in% names(props)) {
		s <- trimws(as.character(props$smiles))
		s <- s[nzchar(s)]
		seen <- c(seen, s)
	}

	unique(seen)
}

# ---------------------------
# MAIN LOOP
# ---------------------------
printf("[%s] Loading master SMILES: %s", ts(), master_smiles)
all_smiles <- load_master_smiles(master_smiles)
total_compounds <- length(all_smiles)
printf("[%s] Master unique SMILES: %d", ts(), total_compounds)

cycle <- 0L

repeat {
	cycle <- cycle + 1L

	seen_smiles <- load_seen_smiles(out_dir)

	# "done_n" is just used for progress messaging
	done_n <- length(intersect(all_smiles, seen_smiles))

	todo_smiles <- setdiff(all_smiles, seen_smiles)
	remaining_n <- length(todo_smiles)

	printf("[%s] Cycle %d | Done: %d/%d | Remaining: %d",
		   ts(), cycle, done_n, total_compounds, remaining_n)

	if (!remaining_n) {
		printf("[%s] No remaining SMILES to process. All done.", ts())
		break
	}

	batch_smiles <- head(todo_smiles, batch_n)
	batch_size <- length(batch_smiles)

	start_idx <- done_n + 1L
	end_idx   <- min(done_n + batch_size, total_compounds)

	# Progress message BEFORE starting the batch
	printf("[%s] %d/%d compounds, starting %d-%d (batch size: %d)", ts(), done_n, total_compounds, start_idx, end_idx, batch_size)

	# Write batch file
	fwrite(data.table(SMILES = batch_smiles), batch_smiles_csv)
	printf("[%s] Wrote batch file: %s", ts(), batch_smiles_csv)

	# Run pipeline on this batch
	# We allow the pipeline to continue even if one of the follow-up scripts fails for a batch
	# Because the 'seen' logic will now pick up the missing pieces in the next coordinator run
	run(sprintf("Rscript %s --in_csv %s", script_smiles, batch_smiles_csv))
	run(sprintf("Rscript %s --batch_smiles %s", script_assays, batch_smiles_csv), stop_on_error = FALSE)
	run(sprintf("Rscript %s --batch_smiles %s", script_tox, batch_smiles_csv), stop_on_error = FALSE)

	printf("[%s] Finished cycle %d.\n", ts(), cycle)
}

printf("[%s] DONE coordinator.", ts())
