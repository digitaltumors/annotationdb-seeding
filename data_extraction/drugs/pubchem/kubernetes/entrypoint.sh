#!/bin/bash
set -euo pipefail

# ============================================================================
# GKE Pod Entrypoint
# Usage: entrypoint.sh <batch_index>
#   BATCH_INDEX can come from: $1, $JOB_COMPLETION_INDEX (K8s indexed Job), or env $BATCH_INDEX
# ============================================================================

BUCKET="gs://annotationdb_data_retrieval"

# Determine batch index
BATCH_IDX="${1:-${JOB_COMPLETION_INDEX:-${BATCH_INDEX:-}}}"
if [ -z "$BATCH_IDX" ]; then
    echo "ERROR: No batch index provided. Pass as arg, or set JOB_COMPLETION_INDEX / BATCH_INDEX env var."
    exit 1
fi

# Zero-pad for consistent naming
BATCH_ID=$(printf "batch_%04d" "$BATCH_IDX")
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting $BATCH_ID"

# ---- Local working directory ----
WORK_DIR="/tmp/work"
OUT_DIR="$WORK_DIR/output"
mkdir -p "$WORK_DIR" "$OUT_DIR"

# ---- Download batch input file ----
BATCH_CSV="$WORK_DIR/${BATCH_ID}.csv"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading batch file..."
gsutil -q cp "${BUCKET}/batches/${BATCH_ID}.csv" "$BATCH_CSV"

# ---- Download reference files needed by R scripts ----
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading reference files..."

# homosapien_aids.txt — needed by pubchem_drug_extraction for AID filtering
gsutil -q cp "${BUCKET}/reference/homosapien_aids.txt" "$WORK_DIR/homosapien_aids.txt" 2>/dev/null || \
    echo "WARNING: homosapien_aids.txt not found in GCS, AID filtering will be skipped"

# completed_aids.txt — needed by assay_extraction_input_ids to skip redundant API calls
gsutil -q cp "${BUCKET}/reference/completed_aids.txt" "$WORK_DIR/completed_aids.txt" 2>/dev/null || \
    echo "WARNING: completed_aids.txt not found in GCS, all missing AIDs will be scraped"

# master_union_out.csv — needed by pubchem_drug_extraction to skip redundant properties scraping
gsutil -q cp "${BUCKET}/output/union_out.csv" "$WORK_DIR/master_union_out.csv" 2>/dev/null || \
    echo "WARNING: master_union_out.csv not found in GCS, properties will be scraped"

# ---- Step 1: PubChem drug extraction (properties, synonyms, bioassays) ----
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running pubchem_drug_extraction_from_input_ids.R ..."
cd "$WORK_DIR"
Rscript /app/pubchem_drug_extraction_from_input_ids.R \
    --in_csv "$BATCH_CSV" \
    --out_dir "$OUT_DIR"

# ---- Step 2: Assay field extraction ----
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running assay_extraction_input_ids.R ..."
Rscript /app/assay_extraction_input_ids.R \
    --batch_input_ids "$BATCH_CSV" \
    --out_dir "$OUT_DIR" || \
    echo "WARNING: assay_extraction had errors (non-fatal, will be retried in future runs)"

# ---- Step 3: Toxicity extraction ----
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running toxicity_extraction_input_ids.R ..."
Rscript /app/toxicity_extraction_input_ids.R \
    --batch_input_ids "$BATCH_CSV" \
    --out_dir "$OUT_DIR" || \
    echo "WARNING: toxicity_extraction had errors (non-fatal, will be retried in future runs)"

# ---- Upload results to GCS ----
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Uploading results to ${BUCKET}/output_batches/${BATCH_ID}/ ..."
gsutil -m -q cp "$OUT_DIR"/*.csv "${BUCKET}/output_batches/${BATCH_ID}/" 2>/dev/null || \
    echo "WARNING: No output CSVs to upload (batch may have had no new work)"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] $BATCH_ID DONE."
