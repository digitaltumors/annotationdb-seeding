#!/usr/bin/env python3
"""
Merge per-batch outputs from GKE pods into final unified CSVs.

Downloads all batch outputs from gs://annotationdb_data_retrieval/output_batches/batch_XXXX/,
concatenates them with the existing outputs, deduplicates, and uploads to output/.

Usage:
    python merge_outputs.py [--dry-run]
"""

import argparse
import csv
import io
import os
import subprocess
import sys
import tempfile


BUCKET = "gs://annotationdb_data_retrieval"
OUTPUT_BATCHES = f"{BUCKET}/output_batches"
OUTPUT_FINAL = f"{BUCKET}/output"

# Files to merge and their dedup keys
MERGE_SPECS = [
    {
        "filename": "union_out.csv",
        "dedup_key": "cid",
        "merge_field": "mapped_name",  # semicolon-separated; merge across batches
    },
    {
        "filename": "union_synonyms.csv",
        "dedup_key": ["synonym", "pubchem_cid"],
    },
    {
        "filename": "union_bioassays.csv",
        "dedup_key": ["pubchem_cid", "bioassay_aid"],
    },
    {
        "filename": "toxicity_output.csv",
        "dedup_key": "pubchem_cid",
    },
    {
        "filename": "union_pubchem_assay_fields.csv",
        "dedup_key": "aid",
    },
]


def gsutil_ls(prefix: str) -> list[str]:
    """List GCS paths under a prefix."""
    try:
        result = subprocess.run(
            ["gsutil", "ls", prefix],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            return []
        return [line.strip() for line in result.stdout.strip().split("\n") if line.strip()]
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return []


def gsutil_cat(gcs_path: str) -> str | None:
    """Read a GCS file's contents."""
    try:
        result = subprocess.run(
            ["gsutil", "-q", "cat", gcs_path],
            capture_output=True, text=True, timeout=300
        )
        return result.stdout if result.returncode == 0 else None
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


def gsutil_cp(local_path: str, gcs_path: str):
    subprocess.run(["gsutil", "-q", "cp", local_path, gcs_path], check=True, timeout=300)


def merge_mapped_names(existing: str, new: str) -> str:
    """Merge semicolon-separated mapped_name fields."""
    tokens = set()
    for val in [existing, new]:
        if val:
            for t in val.split(";"):
                t = t.strip()
                if t:
                    tokens.add(t)
    return "; ".join(sorted(tokens)) if tokens else ""


def merge_file(spec: dict, batch_dirs: list[str], dry_run: bool):
    """Merge a single output file across all batches + existing output."""
    filename = spec["filename"]
    dedup_key = spec["dedup_key"]
    merge_field = spec.get("merge_field")

    print(f"\n--- Merging {filename} ---")

    # Load existing output first
    all_rows: list[dict] = []
    fieldnames: list[str] = []

    existing_text = gsutil_cat(f"{OUTPUT_FINAL}/{filename}")
    if existing_text:
        reader = csv.DictReader(io.StringIO(existing_text))
        fieldnames = list(reader.fieldnames or [])
        all_rows.extend(reader)
        print(f"  Existing: {len(all_rows)} rows")

    # Load per-batch outputs
    batch_count = 0
    for batch_dir in batch_dirs:
        gcs_path = f"{batch_dir}{filename}"
        text = gsutil_cat(gcs_path)
        if text is None:
            continue
        reader = csv.DictReader(io.StringIO(text))
        batch_rows = list(reader)
        if not batch_rows:
            continue

        # Capture fieldnames from first file if we don't have them yet
        if not fieldnames:
            fieldnames = list(reader.fieldnames or [])
        else:
            # Add any new columns from this batch
            for fn in (reader.fieldnames or []):
                if fn not in fieldnames:
                    fieldnames.append(fn)

        all_rows.extend(batch_rows)
        batch_count += 1

    print(f"  Loaded from {batch_count} batches, total: {len(all_rows)} rows")

    if not all_rows:
        print(f"  Nothing to merge for {filename}")
        return

    # Deduplicate
    if isinstance(dedup_key, list):
        key_fn = lambda row: tuple(row.get(k, "").strip() for k in dedup_key)
    else:
        key_fn = lambda row: row.get(dedup_key, "").strip()

    seen: dict = {}  # key -> row index
    deduped: list[dict] = []

    for row in all_rows:
        k = key_fn(row)
        if k in seen:
            # If there's a merge field, merge it
            if merge_field and merge_field in row:
                idx = seen[k]
                old_val = deduped[idx].get(merge_field, "")
                new_val = row.get(merge_field, "")
                deduped[idx][merge_field] = merge_mapped_names(old_val, new_val)
        else:
            seen[k] = len(deduped)
            deduped.append(row)

    print(f"  After dedup: {len(deduped)} rows")

    if dry_run:
        print(f"  [dry-run] Would upload {len(deduped)} rows to {OUTPUT_FINAL}/{filename}")
        return

    # Write to temp file and upload
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(deduped)
        tmp_path = f.name

    gsutil_cp(tmp_path, f"{OUTPUT_FINAL}/{filename}")
    os.unlink(tmp_path)
    print(f"  Uploaded to {OUTPUT_FINAL}/{filename}")


def main():
    parser = argparse.ArgumentParser(description="Merge GKE batch outputs")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    # List batch directories
    print(f"Listing batch outputs from {OUTPUT_BATCHES}/ ...")
    batch_dirs = gsutil_ls(f"{OUTPUT_BATCHES}/")
    batch_dirs = [d for d in batch_dirs if d.endswith("/")]
    print(f"Found {len(batch_dirs)} batch directories")

    if not batch_dirs:
        print("No batch outputs found. Nothing to merge.")
        return

    for spec in MERGE_SPECS:
        merge_file(spec, batch_dirs, args.dry_run)

    # Also merge error logs
    for err_file in ["pubchem_extraction_errors_input_ids.csv",
                     "assay_extraction_errors.csv",
                     "toxicity_errors.csv",
                     "failed_input_ids.csv"]:
        if "toxicity" in err_file:
            dedup = ["time", "pubchem_cid"]
        elif "assay" in err_file:
            dedup = ["time", "aid"]
        elif "failed" in err_file:
            dedup = "input_id"
        else:
            dedup = ["ts", "input_id"]

        merge_file(
            {"filename": err_file, "dedup_key": dedup},
            batch_dirs,
            args.dry_run
        )

    print("\n✅ Merge complete!")


if __name__ == "__main__":
    main()
