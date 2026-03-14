"""
Audit script for CSV data files used in DB seeding.
Checks for:
  1. Duplicate rows and duplicate key values in each file.
  2. Cross-CSV foreign key validity.
"""

import csv
import sys
from collections import Counter
from pathlib import Path

BASE = Path("/Users/mattbocc/uhn/annotationdb-seeding-testing/data_extraction/drugs/pubchem/output_data/union/complete")

FILES = {
    "union_out":                BASE / "union_out.csv",
    "union_synonyms":           BASE / "union_synonyms.csv",
    "toxicity_output":          BASE / "toxicity_output.csv",
    "union_bioassays":          BASE / "union_bioassays.csv",
    "union_pubchem_assay_fields": BASE / "union_pubchem_assay_fields.csv",
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
    ("union_synonyms",           "pubchem_cid",    "union_out",                  "cid"),
    ("toxicity_output",          "pubchem_cid",    "union_out",                  "cid"),
    ("union_bioassays",          "pubchem_cid",    "union_out",                  "cid"),
    ("union_bioassays",          "bioassay_aid",   "union_pubchem_assay_fields", "aid"),
]


def read_csv(path):
    rows = []
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
        issues.append(f"{sum(v - 1 for v in full_dups.values())} duplicate rows (full row match)")

    # Key-value duplicates (only if all key cols are present)
    available_key_cols = [c for c in key_cols if rows and c in rows[0].keys()]
    if available_key_cols:
        keys = [tuple(r.get(c, "").lower() for c in available_key_cols) for r in rows]
        key_counts = Counter(keys)
        key_dups = {k: v for k, v in key_counts.items() if v > 1}
        if key_dups:
            issues.append(f"{len(key_dups)} duplicate key combinations on ({', '.join(available_key_cols)}): {sum(v - 1 for v in key_dups.values())} extra rows")
            for k, v in list(key_dups.items())[:10]:
                issues.append(f"       - key={dict(zip(available_key_cols, k))} appears {v}x")
    else:
        issues.append(f"Key columns {key_cols} not found — skipping key duplicate check")

    return issues


def audit_foreign_key(data, child_name, child_col, parent_name, parent_col):
    issues = []
    child_rows = data[child_name]
    parent_rows = data[parent_name]

    if not child_rows:
        return [f"{child_name} is empty — skipping FK check on {child_col}"]

    if child_col not in child_rows[0]:
        return [f"Column '{child_col}' not found in {child_name} — skipping FK check"]

    if not parent_rows:
        return [f"{parent_name} is empty — skipping FK check for {parent_col}"]

    if parent_col not in parent_rows[0]:
        return [f"Column '{parent_col}' not found in {parent_name} — skipping FK check"]

    parent_vals = set(r[parent_col].strip() for r in parent_rows if r[parent_col].strip())
    child_vals = [r[child_col].strip() for r in child_rows]

    missing = [v for v in child_vals if v and v not in parent_vals]
    if missing:
        missing_unique = sorted(set(missing))
        issues.append(
            f"{len(missing)} broken FK refs in {child_name}.{child_col} "
            f"{parent_name}.{parent_col} "
            f"({len(missing_unique)} unique missing values)"
        )
        for v in missing_unique[:10]:
            issues.append(f"       - '{v}'")
        if len(missing_unique) > 10:
            issues.append(f"       ... and {len(missing_unique) - 10} more")
    else:
        issues.append(f"All {len(child_vals)} refs in {child_name}.{child_col} to{parent_name}.{parent_col} are valid")

    return issues


def main():
    print("=" * 70)
    print("  DATA FILE AUDIT REPORT")
    print("=" * 70)

    data = {}
    row_counts = {}
    for name, path in FILES.items():
        if not path.exists():
            print(f"\n[{name}] FILE NOT FOUND: {path}")
            data[name] = []
            continue
        rows = read_csv(path)
        data[name] = rows
        row_counts[name] = len(rows)
        print(f"\n{'─'*70}")
        print(f"[{name}]  rows={len(rows):,}  file={path.name}")

        issues = audit_duplicates(name, rows, PRIMARY_KEYS.get(name, []))
        if issues:
            for i in issues:
                print(i)
        else:
            print(f"No duplicates found")

    print(f"\n{'═'*70}")
    print("  FOREIGN KEY VALIDATION")
    print(f"{'═'*70}")

    for child_name, child_col, parent_name, parent_col in FOREIGN_KEYS:
        print(f"\n  {child_name}.{child_col}  to {parent_name}.{parent_col}")
        issues = audit_foreign_key(data, child_name, child_col, parent_name, parent_col)
        for i in issues:
            print(i)

    print(f"\n{'═'*70}")
    print("  SUMMARY")
    print(f"{'═'*70}")
    for name, count in row_counts.items():
        print(f"  {name:<35} {count:>8,} rows")
    print(f"{'═'*70}")


if __name__ == "__main__":
    main()
