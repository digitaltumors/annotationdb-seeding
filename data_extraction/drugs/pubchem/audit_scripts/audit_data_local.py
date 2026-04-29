"""
Extended Audit script for CSV data files used in DB seeding.
Checks for:
  1. Duplicate rows and duplicate key values in each file.
  2. Cross-CSV foreign key validity.
  3. Overlap between GCS input files and union_out.csv.
"""

import csv
import os
import sys
import tempfile
import subprocess
from collections import Counter
from pathlib import Path

MATPLOTLIB_AVAILABLE = True
try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend for file export
    import matplotlib.pyplot as plt
    import upsetplot
except ImportError:
    MATPLOTLIB_AVAILABLE = False

INPUT_BUCKET = "gs://annotationdb_data_retrieval/input"
TMPDIR = Path(tempfile.mkdtemp(prefix="audit_cloud_extended_"))

# The local folder containing the mar-27-2026 union data
LOCAL_OUTPUT_DIR = Path(__file__).resolve().parent.parent / "output_data" / "union" / "mar-27-2026"

LOCAL_FILES = {
    "union_out":                  LOCAL_OUTPUT_DIR / "union_out.csv",
    "chembl_mechanism":             LOCAL_OUTPUT_DIR / "chembl_mechanism.csv",
    "union_synonyms":             LOCAL_OUTPUT_DIR / "union_synonyms.csv",
    "toxicity_output":            LOCAL_OUTPUT_DIR / "toxicity_output.csv",
    "union_bioassays":            LOCAL_OUTPUT_DIR / "union_bioassays.csv",
    "union_pubchem_assay_fields": LOCAL_OUTPUT_DIR / "union_pubchem_assay_fields.csv",
}

# Primary / natural keys per file
PRIMARY_KEYS = {
    "union_out":                  ["cid"],
    "chembl_mechanism":             ["molecule_chembl_id", "mec_id", "activity_id"],
    "union_synonyms":             ["synonym", "pubchem_cid", "source"],
    "toxicity_output":            ["pubchem_cid"],
    "union_bioassays":            ["pubchem_cid", "bioassay_aid"],
    "union_pubchem_assay_fields": ["aid"],
}

# Foreign keys: (child_file, child_col) -> (parent_file, parent_col)
FOREIGN_KEYS = [
    ("union_synonyms",           "pubchem_cid",  "union_out",                  "cid"),
    ("chembl_mechanism",         "molecule_chembl_id", "union_out",              "chembl_id"),
    ("toxicity_output",          "pubchem_cid",  "union_out",                  "cid"),
    ("union_bioassays",          "pubchem_cid",  "union_out",                  "cid"),
    ("union_bioassays",          "bioassay_aid", "union_pubchem_assay_fields", "aid"),
]


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
    child_rows = data.get(child_name)
    parent_rows = data.get(parent_name)

    if not child_rows:
        return [f"  {child_name} is empty or missing — skipping FK check on {child_col}"]
    if child_col not in child_rows[0]:
        return [f"  Column '{child_col}' not found in {child_name} — skipping FK check"]
    if not parent_rows:
        return [f"  {parent_name} is empty or missing — skipping FK check for {parent_col}"]
    if parent_col not in parent_rows[0]:
        return [f"  Column '{parent_col}' not found in {parent_name} — skipping FK check"]

    parent_vals = set(r.get(parent_col, "").strip() for r in parent_rows if r.get(parent_col, "").strip())
    child_vals = [r.get(child_col, "").strip() for r in child_rows]

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


def download_input_files():
    print(f"\nDownloading input files from {INPUT_BUCKET} ...")
    cmd = f"gsutil ls {INPUT_BUCKET}/*.csv"
    try:
        output = subprocess.check_output(cmd, shell=True, text=True)
    except subprocess.CalledProcessError:
        print("Failed to list input bucket.")
        return []
        
    gcs_files = [line.strip() for line in output.splitlines() if line.strip().endswith('.csv')]
    local_inputs = []
    
    for gcs_f in gcs_files:
        filename = gcs_f.split('/')[-1]
        local_path = TMPDIR / filename
        print(f"  Downloading {filename} ...", flush=True)
        ret = os.system(f'gsutil -q -o "GSUtil:parallel_process_count=1" cp "{gcs_f}" "{local_path}"')
        if ret == 0:
            local_inputs.append(local_path)
        else:
            print(f"  WARNING: Failed to download {filename}. Skipping.")
    return local_inputs


def check_input_overlap(union_rows, input_files):
    union_cids = set()
    inchikey_to_cid = {}
    
    for r in union_rows:
        cid = r.get("cid", "").strip()
        ik = r.get("inchikey", "").strip()
        if cid:
            union_cids.add(cid)
            if ik:
                inchikey_to_cid[ik] = cid
                
    dataset_cids = {} 
    
    print(f"\n{'=' * 70}")
    print("  INPUT DATASETS COVERAGE IN UNION_OUT")
    print(f"{'=' * 70}")
    
    for f in input_files:
        name = f.name
        rows = read_csv(f)
        if not rows:
            print(f"  [{name}] Empty or unreadable file.")
            continue
            
        cids_in_union = set()
        total_input = len(rows)
        
        cols = list(rows[0].keys())
        
        cid_col = None
        ik_col = None
        for c in cols:
            if not c: continue
            cl = c.lower()
            if cl in ['cid', 'pubchem_cid', 'pubchemcid', 'pubchem cid']:
                cid_col = c
            elif cl == 'inchikey':
                ik_col = c
                
        unique_input_ids = set()
        if cid_col:
            for r in rows:
                val = r.get(cid_col, "").strip()
                if val:
                    unique_input_ids.add(val)
                    if val in union_cids:
                        cids_in_union.add(val)
        elif ik_col:
            for r in rows:
                val = r.get(ik_col, "").strip()
                if val:
                    unique_input_ids.add(val)
                    if val in inchikey_to_cid:
                        cids_in_union.add(inchikey_to_cid[val])
        else:
            print(f"  [{name}] No known ID col found for matching (looked for cid/inchikey). Cols: {cols}")
            continue
            
        dataset_cids[name] = cids_in_union
        unique_input_count = len(unique_input_ids)
        pct = len(cids_in_union) / unique_input_count * 100 if unique_input_count > 0 else 0.0
        print(f"  [{name:<35}] {len(cids_in_union):>8,} unique CIDs / {unique_input_count:>8,} unique input rows ({pct:>5.1f}%)")
        
    print(f"\n{'=' * 70}")
    print("  INPUT DATASETS MATRIX OVERLAP (AMONG INCLUDED IN UNION_OUT)")
    print(f"{'=' * 70}")
    
    names = sorted(list(dataset_cids.keys()))
    if not names:
        print("  No datasets could be compared.")
        return

    # Print a matrix header
    print("  " + "            " + " ".join([f"[{n[:10]:>10}]" for n in names]))
    
    matrix_csv_path = LOCAL_OUTPUT_DIR / "dataset_overlap_matrix.csv"
    try:
        with open(matrix_csv_path, mode='w', newline='', encoding='utf-8') as f_csv:
            writer = csv.writer(f_csv)
            writer.writerow(["Dataset"] + names)
            
            for n1 in names:
                row_str = f"[{n1[:10]:<10}] "
                csv_row = [n1]
                s1 = dataset_cids[n1]
                for n2 in names:
                    s2 = dataset_cids[n2]
                    overlap = len(s1.intersection(s2))
                    row_str += f"{overlap:>12,} "
                    csv_row.append(overlap)
                print(row_str)
                writer.writerow(csv_row)
        print(f"\n  Saved matrix overlap to: {matrix_csv_path}")
    except Exception as e:
        print(f"\n  Could not save CSV matrix: {e}")
        for n1 in names:
            row_str = f"[{n1[:10]:<10}] "
            s1 = dataset_cids[n1]
            for n2 in names:
                s2 = dataset_cids[n2]
                overlap = len(s1.intersection(s2))
                row_str += f"{overlap:>12,} "
            print(row_str)

    print("\n  Detailed overlap pairs (>0):")
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            n1 = names[i]
            n2 = names[j]
            s1 = dataset_cids[n1]
            s2 = dataset_cids[n2]
            overlap = len(s1.intersection(s2))
            if overlap > 0:
                print(f"    {n1} ∩ {n2}: {overlap:,} compounds")

    print(f"\n{'=' * 70}")
    print("  GLOBAL OVERLAP (ACROSS ALL DATASETS)")
    print(f"{'=' * 70}")
    
    if len(names) > 1:
        all_sets = list(dataset_cids.values())
        global_intersection = set.intersection(*all_sets)
        global_union = set.union(*all_sets)
        print(f"  Compounds present in ALL {len(names)} datasets simultaneously: {len(global_intersection):,} CIDs")
        print(f"  Total unique compounds across all {len(names)} datasets combined: {len(global_union):,} CIDs")

    return dataset_cids


def generate_upset_plot(dataset_cids: dict, output_dir: Path):
    if not MATPLOTLIB_AVAILABLE:
        print("  Skipping UpSet plot — matplotlib/upsetplot not available.")
        return

    if len(dataset_cids) < 2:
        print("  Skipping UpSet plot — need at least 2 datasets.")
        return

    # Shorten labels for readability: strip suffix like '_needed.csv', '.csv'
    def short_name(n):
        for suffix in ('_needed.csv', '.csv'):
            if n.endswith(suffix):
                return n[: -len(suffix)]
        return n

    named = {short_name(k): v for k, v in dataset_cids.items()}
    all_cids = set.union(*named.values())

    # Build membership boolean Series indexed by MultiIndex of dataset names
    import pandas as pd
    dataset_names = sorted(named.keys())
    memberships = upsetplot.from_memberships(
        memberships=[
            [ds for ds in dataset_names if cid in named[ds]]
            for cid in all_cids
        ]
    )

    fig = plt.figure(figsize=(max(14, len(dataset_names) * 2), 8))
    upset = upsetplot.UpSet(
        memberships,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        totals_plot_elements=3,
        min_subset_size=1,
    )
    upset.plot(fig)
    fig.suptitle(
        "Compound Overlap Across Input Datasets (CIDs in union_out)",
        fontsize=13,
        fontweight="bold",
        y=1.01,
    )
    out_path = output_dir / "dataset_upset_plot.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  UpSet plot saved to: {out_path}")


def main():
    print("=" * 70)
    print("  EXTENDED DATA FILE AUDIT REPORT (LOCAL + CLOUD INPUTS)")
    print("=" * 70)
    
    data = {}
    row_counts = {}

    print(f"\n{'═' * 70}")
    print(f"  LOCAL OUTPUT DATA CHECKS ({LOCAL_OUTPUT_DIR.name})")
    print(f"{'═' * 70}")
    
    if not LOCAL_OUTPUT_DIR.exists():
        print(f"ERROR: Local directory {LOCAL_OUTPUT_DIR} does not exist.")
        sys.exit(1)

    for name, local_path in LOCAL_FILES.items():
        print(f"\n{'─' * 70}")
        if not local_path.exists():
            print(f"[{name}] SKIPPED — file not found at {local_path}")
            data[name] = []
            continue

        print(f"Reading {local_path.name} ...", flush=True)
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
    
    if "union_out" in data and data["union_out"]:
        input_files = download_input_files()
        if input_files:
            dataset_cids = check_input_overlap(data["union_out"], input_files)
            generate_upset_plot(dataset_cids, LOCAL_OUTPUT_DIR)
    else:
        print("Skipping input overlap check because union_out.csv was not loaded.")

    # Cleanup temp files
    try:
        for f in TMPDIR.iterdir():
            f.unlink()
        TMPDIR.rmdir()
    except Exception as e:
        print(f"Warning: Failed to clean up temp dir {TMPDIR}: {e}")


if __name__ == "__main__":
    main()
