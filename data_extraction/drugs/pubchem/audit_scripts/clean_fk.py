"""
Cleanup script to purge rows with broken foreign key references. Backs up original files with a .bak extension before overwriting.
Files Processed:
  - toxicity_output.csv (invalid pubchem_cid)
  - union_bioassays.csv (invalid pubchem_cid or bioassay_aid)
  - chembl_mechanism.csv (invalid molecule_chembl_id)
"""

import csv
from pathlib import Path
import shutil

# Base directory where the CSV files live
BASE = Path("/Users/mattbocc/uhn/annotationdb-seeding-testing/data_extraction/drugs/pubchem/output_data/union/mar-27-2026")

# Files we need to read for reference data
UNION_OUT = BASE / "union_out.csv"
UNION_ASSAY_FIELDS = BASE / "union_pubchem_assay_fields.csv"

def load_set(path: Path, column: str) -> set[str]:
    """Load a CSV column into a set of stripped strings."""
    values = set()
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            val = row.get(column, "").strip()
            if val:
                values.add(val)
    return values

# Reference sets
valid_cids = load_set(UNION_OUT, "cid")
valid_aids = load_set(UNION_ASSAY_FIELDS, "aid")

# Helper to clean a file based on a predicate that decides whether to keep a row
def clean_file(src: Path, key_cols: list[str], keep_predicate):
    bak = src.with_suffix(src.suffix + ".bak")
    tmp = src.with_suffix(src.suffix + ".tmp")
    shutil.copy2(src, bak)
    print(f"Backed up {src.name} to {bak.name}")

    total_in = 0
    total_out = 0
    with open(src, newline="", encoding="utf-8") as fin, open(tmp, "w", newline="", encoding="utf-8") as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            total_in += 1
            if keep_predicate(row):
                writer.writerow(row)
                total_out += 1
    tmp.replace(src)
    print(f"{src.name}: {total_in:,} rows in to {total_out:,} rows out ({total_in - total_out:,} removed)")

# 1. toxicity_output.csv – keep rows where pubchem_cid exists in union_out
toxicity_path = BASE / "toxicity_output.csv"
clean_file(
    toxicity_path,
    ["pubchem_cid"],
    lambda r: r.get("pubchem_cid", "").strip() in valid_cids,
)

# 2. union_bioassays.csv – keep rows where both pubchem_cid and bioassay_aid are valid
bioassays_path = BASE / "union_bioassays.csv"
clean_file(
    bioassays_path,
    ["pubchem_cid", "bioassay_aid"],
    lambda r: (
        r.get("pubchem_cid", "").strip() in valid_cids
        and r.get("bioassay_aid", "").strip() in valid_aids
    ),
)

# 3. union_synonyms.csv – keep rows where pubchem_cid exists in union_out
synonyms_path = BASE / "union_synonyms.csv"
clean_file(
    synonyms_path,
    ["pubchem_cid"],
    lambda r: r.get("pubchem_cid", "").strip() in valid_cids,
)

print("\nCleanup complete. Originals saved as .bak files.")
