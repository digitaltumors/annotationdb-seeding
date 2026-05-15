#!/usr/bin/env bash
# combine_shards.sh
#
# Combines toxicity_output.csv files from all shards into a single file.
# Deduplicates on pubchem_cid + reference_number so re-runs don't produce
# duplicate rows.
#
# Usage:
#   ./combine_shards.sh [--output_dir <dir>] [--num_shards <n>] [--out <file>]
#
# Defaults:
#   --output_dir  ./output
#   --num_shards  4
#   --out         ./output/toxicity_combined.csv

set -euo pipefail

OUTPUT_DIR="./output"
NUM_SHARDS=4
COMBINED="./output/toxicity_combined.csv"

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output_dir)  OUTPUT_DIR="$2";  shift 2 ;;
        --num_shards)  NUM_SHARDS="$2";  shift 2 ;;
        --out)         COMBINED="$2";    shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

echo "=== Combining $NUM_SHARDS shards from: $OUTPUT_DIR ==="
echo "=== Output: $COMBINED ==="

# Collect all shard files that exist
SHARD_FILES=()
MISSING=()
for i in $(seq 0 $((NUM_SHARDS - 1))); do
    f="$OUTPUT_DIR/shard_${i}/toxicity_output.csv"
    if [[ -f "$f" ]]; then
        SHARD_FILES+=("$f")
    else
        MISSING+=("shard_${i}")
    fi
done

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "WARNING: missing shards: ${MISSING[*]}"
fi

if [[ ${#SHARD_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No shard files found in $OUTPUT_DIR"
    exit 1
fi

# Print shard row counts
TOTAL_RAW=0
for f in "${SHARD_FILES[@]}"; do
    n=$(tail -n +2 "$f" | wc -l | tr -d ' ')
    echo "  $(basename $(dirname $f)): $n rows"
    TOTAL_RAW=$((TOTAL_RAW + n))
done
echo "  Total raw rows (before dedup): $TOTAL_RAW"

# Write header from first shard
head -1 "${SHARD_FILES[0]}" > "$COMBINED"

# Append data rows from all shards, then deduplicate on cols 1+2 (pubchem_cid, reference_number)
# awk keeps first occurrence of each (col1, col2) key
for f in "${SHARD_FILES[@]}"; do
    tail -n +2 "$f"
done | awk -F',' '
    {
        # Strip surrounding quotes from key columns
        k1 = $1; gsub(/"/, "", k1)
        k2 = $2; gsub(/"/, "", k2)
        key = k1 SUBSEP k2
        if (!(key in seen)) {
            seen[key] = 1
            print
        }
    }
' >> "$COMBINED"

FINAL=$(tail -n +2 "$COMBINED" | wc -l | tr -d ' ')
echo ""
echo "=== Done ==="
echo "  Output rows (after dedup): $FINAL"
echo "  File: $COMBINED"
