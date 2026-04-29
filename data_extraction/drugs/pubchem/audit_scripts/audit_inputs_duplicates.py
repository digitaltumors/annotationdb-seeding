"""
Audit script specifically for GCS input datasets.
Detects full-row duplicates and duplicate identifiers (cid or inchikey).
"""

import csv
import sys
import tempfile
import subprocess
import os
from collections import Counter
from pathlib import Path

INPUT_BUCKET = "gs://annotationdb_data_retrieval/input"
TMPDIR = Path(tempfile.mkdtemp(prefix="audit_inputs_dups_"))

def read_csv(path):
    rows = []
    csv.field_size_limit(sys.maxsize)
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

def check_duplicates(name, rows):
    issues = []
    
    # 1. Full-row duplicates
    serialized = [tuple(r.values()) for r in rows]
    row_counts = Counter(serialized)
    full_dups = {k: v for k, v in row_counts.items() if v > 1}
    
    if full_dups:
        issues.append(f"  [FULL ROW] {sum(v - 1 for v in full_dups.values())} extra duplicated rows (exactly identical data across all columns)")
    else:
         issues.append("  [FULL ROW] No exactly identical rows found.")

    # 2. Key-based duplicates
    if not rows:
        return issues
        
    cols = list(rows[0].keys())
    id_col = None
    for c in cols:
        if not c: continue
        cl = c.lower()
        if cl in ['cid', 'pubchem_cid', 'pubchemcid', 'pubchem cid', 'inchikey']:
            id_col = c
            break

    if id_col:
        keys = [r.get(id_col, "").strip() for r in rows if r.get(id_col, "").strip()]
        key_counts = Counter(keys)
        key_dups = {k: v for k, v in key_counts.items() if v > 1}
        
        if key_dups:
            issues.append(f"  [ID: {id_col}] {sum(v - 1 for v in key_dups.values())} extra duplicated rows sharing the same identifier")
            # print top 5 examples
            top = list(key_dups.items())[:5]
            for val, count in top:
                issues.append(f"       - '{val}' appears {count} times")
        else:
            issues.append(f"  [ID: {id_col}] No duplicate identifiers found.")
    else:
        issues.append(f"  [ID] No CID or InChIKey column found to check for key duplicates.")
        
    return issues


def main():
    print("=" * 70)
    print("  INPUT DATASETS DUPLICATE AUDIT")
    print("=" * 70)
    print(f"\nDownloading input files from {INPUT_BUCKET} ...")
    cmd = f"gsutil ls {INPUT_BUCKET}/*.csv"
    try:
        output = subprocess.check_output(cmd, shell=True, text=True)
    except subprocess.CalledProcessError:
        print("Failed to list input bucket.")
        return
        
    gcs_files = [line.strip() for line in output.splitlines() if line.strip().endswith('.csv')]
    
    for gcs_f in gcs_files:
        filename = gcs_f.split('/')[-1]
        local_path = TMPDIR / filename
        print(f"\n{'─' * 70}")
        print(f"[{filename}] Downloading ...", flush=True)
        ret = os.system(f'gsutil -q -o "GSUtil:parallel_process_count=1" cp "{gcs_f}" "{local_path}"')
        if ret != 0:
            print(f"  WARNING: Failed to download.")
            continue
            
        print(f"[{filename}] Analyzing ...", flush=True)
        rows = read_csv(local_path)
        print(f"  Total Rows: {len(rows):,}")
        
        issues = check_duplicates(filename, rows)
        for i in issues:
            print(i)
            
    # Cleanup temp directory
    try:
        for f in TMPDIR.iterdir():
            f.unlink()
        TMPDIR.rmdir()
    except Exception:
        pass

if __name__ == "__main__":
    main()
