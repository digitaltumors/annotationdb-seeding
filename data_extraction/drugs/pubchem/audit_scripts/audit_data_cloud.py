"""
Audit script for CSV data files used in DB seeding.
Cloud version: downloads files from GCS before auditing.
Checks for:
  1. Duplicate rows and duplicate key values in each file.
  2. Cross-CSV foreign key validity.
"""

import csv
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

BUCKET = "gs://annotationdb_data_retrieval"
TMPDIR = Path(tempfile.mkdtemp(prefix="audit_cloud_"))

GCS_FILES = {
    "union_out":                  f"{BUCKET}/output/union_out.csv",
    "union_synonyms":             f"{BUCKET}/output/union_synonyms.csv",
    "toxicity_output":            f"{BUCKET}/output/toxicity_output.csv",
    "union_bioassays":            f"{BUCKET}/output/union_bioassays.csv",
    "union_pubchem_assay_fields": f"{BUCKET}/output/union_pubchem_assay_fields.csv",
}

# Primary / natural keys per file
PRIMARY_KEYS = {
    "union_out":                  ["cid"],
    "union_synonyms":             ["synonym", "pubchem_cid", "source"],
    "toxicity_output":            ["pubchem_cid"],
    "union_bioassays":            ["pubchem_cid", "bioassay_aid"],
    "union_pubchem_assay_fields": ["aid"],
}

# Foreign keys: (child_file, child_col) -> (parent_file, parent_col)
FOREIGN_KEYS = [
    ("union_synonyms",           "pubchem_cid",  "union_out",                  "cid"),
    ("toxicity_output",          "pubchem_cid",  "union_out",                  "cid"),
    ("union_bioassays",          "pubchem_cid",  "union_out",                  "cid"),
    ("union_bioassays",          "bioassay_aid", "union_pubchem_assay_fields", "aid"),
]


def download_files():
    local_files = {}
    for name, gcs_path in GCS_FILES.items():
        local_path = TMPDIR / f"{name}.csv"
        print(f"  Downloading {gcs_path} ...", flush=True)
        ret = os.system(f'gsutil -o "GSUtil:parallel_process_count=1" cp "{gcs_path}" "{local_path}"')
        if ret != 0:
            print(f"  WARNING: Failed to download {name}. Skipping.")
            local_files[name] = None
        else:
            local_files[name] = local_path
    return local_files


def read_csv(path):
    rows = []
    csv.field_size_limit(sys.maxsize)
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def audit_duplicates(name, rows, key_cols):
    issues = []

    # Full-row duplicates
    serialized = [tuple(r.values()) for r in rows]
    row_counts = Counter(serialized)
    full_dups = {k: v for k, v in row_counts.items() if v > 1}
    if full_dups:
        issues.append(f"  {sum(v - 1 for v in full_dups.values())} duplicate rows (full row match)")

    # Key-value duplicates (only if all key cols are present)
    available_key_cols = [c for c in key_cols if rows and c in rows[0].keys()]
    if available_key_cols:
        keys = [tuple(r.get(c, "").lower() for c in available_key_cols) for r in rows]
        key_counts = Counter(keys)
        key_dups = {k: v for k, v in key_counts.items() if v > 1}
        if key_dups:
            issues.append(f"  {len(key_dups)} duplicate key combos on ({', '.join(available_key_cols)}): "
                          f"{sum(v - 1 for v in key_dups.values())} extra rows")
            for k, v in list(key_dups.items())[:5]:
                issues.append(f"       - key={dict(zip(available_key_cols, k))} appears {v}x")
    else:
        issues.append(f"  Key columns {key_cols} not found — skipping key duplicate check")

    return issues


def audit_foreign_key(data, child_name, child_col, parent_name, parent_col):
    issues = []
    child_rows = data[child_name]
    parent_rows = data[parent_name]

    if not child_rows:
        return [f"  {child_name} is empty — skipping FK check on {child_col}"]
    if child_col not in child_rows[0]:
        return [f"  Column '{child_col}' not found in {child_name} — skipping FK check"]
    if not parent_rows:
        return [f"  {parent_name} is empty — skipping FK check for {parent_col}"]
    if parent_col not in parent_rows[0]:
        return [f"  Column '{parent_col}' not found in {parent_name} — skipping FK check"]

    parent_vals = set(r[parent_col].strip() for r in parent_rows if r[parent_col].strip())
    child_vals = [r[child_col].strip() for r in child_rows]

    missing = [v for v in child_vals if v and v not in parent_vals]
    if missing:
        missing_unique = sorted(set(missing))
        issues.append(
            f"  WARNING: {len(missing)} broken FK refs in {child_name}.{child_col} -> "
            f"{parent_name}.{parent_col} ({len(missing_unique)} unique missing values)"
        )
        for v in missing_unique[:5]:
            issues.append(f"       - '{v}'")
        if len(missing_unique) > 5:
            issues.append(f"       ... and {len(missing_unique) - 5} more")
    else:
        issues.append(f"  OK: All {len(child_vals):,} refs in {child_name}.{child_col} -> {parent_name}.{parent_col} are valid")

    return issues


def main():
    print("=" * 70)
    print("  DATA FILE AUDIT REPORT (CLOUD)")
    print("=" * 70)
    print(f"\nDownloading files from GCS to {TMPDIR} ...\n")

    local_files = download_files()

    data = {}
    row_counts = {}

    print(f"\n{'═' * 70}")
    print("  DUPLICATE CHECKS")
    print(f"{'═' * 70}")

    for name, local_path in local_files.items():
        print(f"\n{'─' * 70}")
        if local_path is None:
            print(f"[{name}] SKIPPED — download failed")
            data[name] = []
            continue

        rows = read_csv(local_path)
        data[name] = rows
        row_counts[name] = len(rows)
        print(f"[{name}]  rows={len(rows):,}  file={local_path.name}")

        issues = audit_duplicates(name, rows, PRIMARY_KEYS.get(name, []))
        if issues:
            for i in issues:
                print(i)
        else:
            print("  No duplicates found ✓")

    print(f"\n{'═' * 70}")
    print("  FOREIGN KEY VALIDATION")
    print(f"{'═' * 70}")

    for child_name, child_col, parent_name, parent_col in FOREIGN_KEYS:
        print(f"\n  {child_name}.{child_col}  ->  {parent_name}.{parent_col}")
        issues = audit_foreign_key(data, child_name, child_col, parent_name, parent_col)
        for i in issues:
            print(i)

    print(f"\n{'═' * 70}")
    print("  SUMMARY")
    print(f"{'═' * 70}")
    for name, count in row_counts.items():
        print(f"  {name:<35} {count:>10,} rows")
    print(f"{'═' * 70}\n")

    # Cleanup temp files
    for f in TMPDIR.iterdir():
        f.unlink()
    TMPDIR.rmdir()


if __name__ == "__main__":
    main()
