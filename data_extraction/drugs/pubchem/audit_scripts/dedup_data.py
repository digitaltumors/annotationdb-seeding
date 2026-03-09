"""
Remove row duplicates from DB seeing CSV files (keeping first occurrence). Script also backs up original files before overwriting.
"""

import csv
import shutil
from pathlib import Path

BASE = Path("/Users/mattbocc/uhn/annotationdb-seeding-testing/data_extraction/drugs/pubchem/output_data/union/complete")

JOBS = [
    {
        "file":    BASE / "union_synonyms.csv",
        "key_cols": ["synonym", "pubchem_cid", "source"],
    },
    {
        "file":    BASE / "union_bioassays.csv",
        "key_cols": ["pubchem_cid", "bioassay_aid"],
    },
]

CHUNK = 100_000  # rows to read at a time


def dedup_file(path: Path, key_cols: list[str]):
    bak = path.with_suffix(path.suffix + ".bak")
    tmp = path.with_suffix(path.suffix + ".tmp")

    print(f"\n[{path.name}]  key=({', '.join(key_cols)})")

    # Back up original
    shutil.copy2(path, bak)
    print(f"  Backed up → {bak.name}")

    seen = set()
    total_in = 0
    total_out = 0

    with open(path, newline="", encoding="utf-8") as fin, \
         open(tmp,  "w",  newline="", encoding="utf-8") as fout:

        reader = csv.DictReader(fin)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            total_in += 1
            key = tuple(row.get(c, "").lower() for c in key_cols)
            if key not in seen:
                seen.add(key)
                writer.writerow(row)
                total_out += 1

    # Replace original with deduped version
    tmp.replace(path)
    removed = total_in - total_out
    print(f"  {total_in:,} rows in  →  {total_out:,} rows out  ({removed:,} duplicates removed)")


def main():
    print("=" * 60)
    print("  DEDUPLICATION")
    print("=" * 60)

    for job in JOBS:
        dedup_file(job["file"], job["key_cols"])

    print("\nDone. Originals saved as .bak files.")


if __name__ == "__main__":
    main()
