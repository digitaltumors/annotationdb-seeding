#!/usr/bin/env python3
"""
GKE Batch Coordinator

Reads the master input CSV, determines which IDs are already done by checking
existing outputs in GCS, splits the remaining work into batch CSVs, uploads
them to GCS, and prints/generates the kubectl command to launch the K8s Job.

Usage:
    python coordinator.py [--batch-size 1000] [--max-batches 5] [--dry-run]
"""

import argparse
import csv
import io
import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


# ============================================================================
# Configuration
# ============================================================================
BUCKET = "gs://annotationdb_data_retrieval"
INPUT_CSV = f"{BUCKET}/input/LINCS_needed.csv"
OUTPUT_PREFIX = f"{BUCKET}/output"
BATCHES_PREFIX = f"{BUCKET}/batches"

# Files to check for completed IDs
UNION_OUT = f"{OUTPUT_PREFIX}/union_out.csv"
FAILED_INPUT_IDS = f"{OUTPUT_PREFIX}/failed_input_ids.csv"
TOXICITY_OUTPUT = f"{OUTPUT_PREFIX}/toxicity_output.csv"
ASSAY_FIELDS_OUTPUT = f"{OUTPUT_PREFIX}/union_pubchem_assay_fields.csv"
COMPLETED_AIDS_REF = f"{BUCKET}/reference/completed_aids.txt"


def gsutil_cat(gcs_path: str) -> str | None:
    """Read a GCS file's contents as a string. Returns None if file doesn't exist."""
    try:
        result = subprocess.run(
            ["gsutil", "-q", "cat", gcs_path],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0:
            return result.stdout
        return None
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


def gsutil_cp(local_path: str, gcs_path: str):
    """Upload a local file to GCS."""
    subprocess.run(
        ["gsutil", "-q", "cp", local_path, gcs_path],
        check=True, timeout=60
    )


def load_master_ids(csv_text: str) -> list[dict]:
    """Parse master CSV and return list of rows as dicts."""
    reader = csv.DictReader(io.StringIO(csv_text))
    rows = list(reader)
    return rows


def detect_columns(fieldnames: list[str]) -> tuple[str, str | None]:
    """Detect input_id and inchikey column names."""
    input_id_col = None
    for cand in ["input_id", "inputId", "substance", "name", "drug", "synonym"]:
        if cand in fieldnames:
            input_id_col = cand
            break
    if input_id_col is None:
        # Use first non-CID, non-InChIKey column
        skip = {"cid", "CID", "pubchem_cid", "InChiKey", "inchikey", "InChIKey", "SMILES"}
        for col in fieldnames:
            if col not in skip:
                input_id_col = col
                break
    if input_id_col is None:
        input_id_col = fieldnames[0]

    inchikey_col = None
    for cand in ["InChiKey", "inchikey", "InChIKey", "InCHIKey"]:
        if cand in fieldnames:
            inchikey_col = cand
            break

    return input_id_col, inchikey_col


def get_done_ids_from_union_out(csv_text: str | None, tox_text: str | None) -> tuple[set[str], set[str], set[str]]:
    """Extract completed input_ids, cids, and inchikeys from union_out.csv,
    cross-referenced with toxicity_output.csv for true completion."""
    done_inputs = set()
    done_cids = set()
    done_inchikeys = set()

    if csv_text is None:
        return done_inputs, done_cids, done_inchikeys

    reader = csv.DictReader(io.StringIO(csv_text))
    
    # Get finished CIDs from toxicity
    finished_cids = None
    if tox_text is not None:
        tox_reader = csv.DictReader(io.StringIO(tox_text))
        if "pubchem_cid" in (tox_reader.fieldnames or []):
            finished_cids = {row["pubchem_cid"].strip() for row in tox_reader
                           if row.get("pubchem_cid", "").strip()}

    for row in reader:
        cid = row.get("cid", "").strip()
        inchikey = row.get("inchikey", "").strip()
        mapped = row.get("mapped_name", "").strip()

        # If we have toxicity data, only count CIDs that are in toxicity as complete
        if finished_cids is not None and cid not in finished_cids:
            continue
            
        if cid:
            try:
                cid_int = int(float(cid))
                done_cids.add(str(cid_int))
            except ValueError:
                done_cids.add(cid)
                
        if inchikey:
            done_inchikeys.add(inchikey)

        for token in mapped.split(";"):
            token = token.strip()
            if token:
                done_inputs.add(token)

    return done_inputs, done_cids, done_inchikeys


def get_done_ids_from_failed(csv_text: str | None) -> set[str]:
    """Extract terminal-failure input_ids from failed_input_ids.csv."""
    done = set()
    if csv_text is None:
        return done

    reader = csv.DictReader(io.StringIO(csv_text))
    if "input_id" not in (reader.fieldnames or []):
        return done

    for row in reader:
        val = row.get("input_id", "").strip()
        if val:
            done.add(val)
    return done


def update_completed_aids_reference(dry_run: bool = False):
    """Download union_pubchem_assay_fields.csv, extract unique AIDs, and upload to reference."""
    print("[coordinator] Extracting completed AIDs from union_pubchem_assay_fields.csv ...")
    with tempfile.TemporaryDirectory() as tmpdir:
        local_csv = os.path.join(tmpdir, "union_pubchem_assay_fields.csv")
        local_txt = os.path.join(tmpdir, "completed_aids.txt")
        
        # Download (ignore if doesn't exist) - using parallel_process_count=1 to fix macOS Python multiprocessing deadlock
        res = subprocess.run([
            "gsutil", 
            "-o", "GSUtil:parallel_process_count=1",
            "cp", ASSAY_FIELDS_OUTPUT, local_csv
        ])
        if res.returncode != 0:
            print("  [coordinator] No existing union_pubchem_assay_fields.csv found. Skipping.")
            return

        aids = set()
        with open(local_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames and "aid" in reader.fieldnames:
                for row in reader:
                    val = row.get("aid", "").strip()
                    if val:
                        aids.add(val)
        
        print(f"  [coordinator] Found {len(aids)} previously processed AIDs.")
        if not aids:
            return

        with open(local_txt, "w", encoding="utf-8") as f:
            for aid in sorted(aids):
                f.write(f"{aid}\n")
        
        if dry_run:
            print(f"  [dry-run] Would upload completed_aids.txt to {COMPLETED_AIDS_REF} ...")
        else:
            print(f"  [coordinator] Uploading completed_aids.txt to {COMPLETED_AIDS_REF} ...")
            subprocess.run(["gsutil", "-q", "cp", local_txt, COMPLETED_AIDS_REF], check=True)


def main():
    parser = argparse.ArgumentParser(description="GKE Batch Coordinator")
    parser.add_argument("--batch-size", type=int, default=1000,
                        help="Number of input IDs per batch (default: 1000)")
    parser.add_argument("--max-batches", type=int, default=5,
                        help="Maximum number of batches to create (default: 5)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be done without uploading")
    parser.add_argument("--job-name", type=str, default="annotationdb-extract",
                        help="K8s Job name")
    args = parser.parse_args()

    # ---- Extract completed AIDs reference file ----
    update_completed_aids_reference(dry_run=args.dry_run)

    # ---- Load master CSV ----
    print(f"[coordinator] Loading master CSV from {INPUT_CSV} ...")
    master_text = gsutil_cat(INPUT_CSV)
    if master_text is None:
        print(f"ERROR: Could not read {INPUT_CSV}", file=sys.stderr)
        sys.exit(1)

    master_rows = load_master_ids(master_text)
    fieldnames = list(master_rows[0].keys()) if master_rows else []
    input_id_col, inchikey_col = detect_columns(fieldnames)
    print(f"[coordinator] Detected input_id column: '{input_id_col}', inchikey column: '{inchikey_col}'")

    # ---- Determine what's already done ----
    print("[coordinator] Checking existing progress in GCS ...")
    union_text = gsutil_cat(UNION_OUT)
    tox_text = gsutil_cat(TOXICITY_OUTPUT)
    failed_text = gsutil_cat(FAILED_INPUT_IDS)

    # Detect cid column for the master CSV
    cid_col = next((c for c in fieldnames if c.lower() in ["cid", "pubchem_cid", "pubchemcid"]), None)

    # --- IN-MEMORY CACHE MERGE OPTIMIZATION ---
    if input_id_col and (inchikey_col or cid_col):
        print("[coordinator] Input file contains InChIKey or CID. Attempting instant local cache merge...")
        inchi_to_names = defaultdict(set)
        cid_to_names = defaultdict(set)

        for row in master_rows:
            val_input = row.get(input_id_col, "").strip()
            if not val_input: continue
                
            if inchikey_col:
                val_inchi = row.get(inchikey_col, "").strip()
                if val_inchi:
                    inchi_to_names[val_inchi].add(val_input)
                
            if cid_col:
                val_cid = row.get(cid_col, "").strip()
                if val_cid:
                    try:
                        vc = str(int(float(val_cid)))
                        cid_to_names[vc].add(val_input)
                    except ValueError:
                        cid_to_names[val_cid].add(val_input)
                    
        if union_text:
            reader = csv.DictReader(io.StringIO(union_text))
            out_fieldnames = reader.fieldnames
            updated_rows = []
            made_changes = False
            
            for cached_row in reader:
                cached_inchi = cached_row.get("inchikey", "").strip()
                cached_cid = cached_row.get("cid", "").strip()
                try: 
                    cached_cid = str(int(float(cached_cid))) if cached_cid else ""
                except ValueError: 
                    pass
                
                matching_names = set()
                if cached_inchi and cached_inchi in inchi_to_names:
                    matching_names.update(inchi_to_names[cached_inchi])
                if cached_cid and cached_cid in cid_to_names:
                    matching_names.update(cid_to_names[cached_cid])
                    
                if matching_names:
                    old_mapping = cached_row.get("mapped_name", "").strip()
                    old_names = {x.strip() for x in old_mapping.split(";")} if old_mapping else set()
                    old_names.discard("")
                    
                    if matching_names - old_names:
                        old_names.update(matching_names)
                        cached_row["mapped_name"] = "; ".join(sorted(old_names))
                        made_changes = True
                
                updated_rows.append(cached_row)
                
            if made_changes and not args.dry_run:
                print(f"[coordinator] Cache hits found! Injecting new input names directly into {UNION_OUT}...")
                output_buf = io.StringIO()
                writer = csv.DictWriter(output_buf, fieldnames=out_fieldnames)
                writer.writeheader()
                writer.writerows(updated_rows)
                union_text = output_buf.getvalue()
                
                with tempfile.NamedTemporaryFile(mode='w', delete=False, newline='', encoding='utf-8') as tmp:
                    tmp.write(union_text)
                    tmp_name = tmp.name
                
                subprocess.run(["gsutil", "-q", "cp", tmp_name, UNION_OUT], check=True)
                os.remove(tmp_name)
                print("[coordinator] Local cache merge complete and safely uploaded!")
            elif made_changes and args.dry_run:
                print("[dry-run] Would merge new names into union_out.csv locally.")
            else:
                print("[coordinator] No new overlapping names to merge into existing cache.")

    done_inputs, done_cids, done_inchikeys = get_done_ids_from_union_out(union_text, tox_text)
    done_from_failed = get_done_ids_from_failed(failed_text)

    # ---- Compute remaining ----
    remaining = []
    already_done_count = 0
    seen_remaining = set()

    for row in master_rows:
        val_input = row.get(input_id_col, "").strip()
        if not val_input:
            continue
            
        val_cid = row.get(cid_col, "").strip() if cid_col else ""
        cid_match = False
        if val_cid:
            try:
                vc = str(int(float(val_cid)))
                if vc in done_cids:
                    cid_match = True
            except ValueError:
                if val_cid in done_cids:
                    cid_match = True
                    
        val_inchi = row.get(inchikey_col, "").strip() if inchikey_col else ""
        inchi_match = val_inchi in done_inchikeys if val_inchi else False

        if val_input in done_inputs or val_input in done_from_failed or cid_match or inchi_match:
            already_done_count += 1
        else:
            if val_input not in seen_remaining:
                remaining.append(val_input)
                seen_remaining.add(val_input)

    print(f"[coordinator] Already done: {already_done_count} "
          f"(union_out inputs: {len(done_inputs)}, cids: {len(done_cids)}, inchikeys: {len(done_inchikeys)}, failed: {len(done_from_failed)})")
          
    print(f"[coordinator] Total unique remaining input IDs: {len(remaining)}")

    if not remaining:
        print("[coordinator] Nothing to do! All IDs are complete.")
        return

    # ---- Split into batches ----
    total_batches = min(
        math.ceil(len(remaining) / args.batch_size),
        args.max_batches
    )
    batches = []
    for i in range(total_batches):
        start = i * args.batch_size
        end = min(start + args.batch_size, len(remaining))
        batch_ids = remaining[start:end]
        if batch_ids:
            batches.append(batch_ids)

    print(f"[coordinator] Created {len(batches)} batches "
          f"(batch_size={args.batch_size}, max_batches={args.max_batches})")

    # ---- Build batch CSV rows from master ----
    # Create a lookup: input_id -> original row
    id_to_row = {}
    for row in master_rows:
        val = row.get(input_id_col, "").strip()
        if val and val not in id_to_row:
            id_to_row[val] = row

    # ---- Upload batch CSVs to GCS ----
    with tempfile.TemporaryDirectory() as tmpdir:
        for i, batch_ids in enumerate(batches):
            batch_name = f"batch_{i:04d}"
            local_path = os.path.join(tmpdir, f"{batch_name}.csv")

            # Write CSV with same columns as master
            with open(local_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for id_ in batch_ids:
                    if id_ in id_to_row:
                        writer.writerow(id_to_row[id_])

            gcs_path = f"{BATCHES_PREFIX}/{batch_name}.csv"
            row_count = len(batch_ids)

            if args.dry_run:
                print(f"  [dry-run] Would upload {batch_name}.csv ({row_count} rows) -> {gcs_path}")
            else:
                gsutil_cp(local_path, gcs_path)
                print(f"  Uploaded {batch_name}.csv ({row_count} rows) -> {gcs_path}")

    # ---- Print kubectl command ----
    print()
    print("=" * 70)
    print(f"Batches ready! To launch the GKE Job:")
    print()
    print(f"  kubectl apply -f job.yaml")
    print()
    print(f"Or to run a single batch locally in Docker:")
    print()
    print(f"  docker run --rm bhklabmattbocc/annotationdb-r:latest 0")
    print("=" * 70)


if __name__ == "__main__":
    main()
