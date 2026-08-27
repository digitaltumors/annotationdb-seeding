import sqlite3
import pandas as pd
import re
import os

DB_PATH = "/Users/mattbocc/Downloads/chembl_36/chembl_36_sqlite/chembl_36.db"
INPUT_CSV = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/seeding/seeding_data/july_30_2026/union_out.csv"
OUTPUT_DIR = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/seeding/seeding_data/july_30_2026"
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "chembl_drug_indication.csv")

print(f"Connecting to SQLite Database: {DB_PATH}...")
conn = sqlite3.connect(DB_PATH)

# ── 1. Fetch drug indications with molecule and hierarchy info ──────────────
print("Fetching ALL drug indications from drug_indication table...")
ind_query = """
SELECT
    di.drugind_id,
    md.chembl_id                  AS molecule_chembl_id,
    parent_md.chembl_id           AS parent_molecule_chembl_id,
    di.max_phase_for_ind,
    di.mesh_id,
    di.mesh_heading,
    di.efo_id,
    di.efo_term
FROM drug_indication di
JOIN molecule_dictionary md ON di.molregno = md.molregno
LEFT JOIN molecule_hierarchy mh ON mh.molregno = md.molregno
LEFT JOIN molecule_dictionary parent_md ON mh.parent_molregno = parent_md.molregno
"""
ind_df = pd.read_sql(ind_query, conn)
print(f"  drug_indication rows: {len(ind_df):,}  |  distinct mols: {ind_df['molecule_chembl_id'].nunique():,}")

# ── 2. Aggregate indication_refs per drugind_id (excluding ATC) ───────────────
# Multiple ref_ids for the same type are joined with commas (e.g. NCT001,NCT002)
print("Fetching and aggregating indication_refs by ref_type (excluding ATC)...")
refs_query = """
SELECT 
    drugind_id,
    GROUP_CONCAT(CASE WHEN ref_type = 'ClinicalTrials' THEN ref_id END, ',') AS clinical_trials_ref_ids,
    GROUP_CONCAT(CASE WHEN ref_type = 'DailyMed' THEN ref_id END, ',') AS daily_med_ref_ids,
    GROUP_CONCAT(CASE WHEN ref_type = 'EMA' THEN ref_id END, ',') AS ema_ref_ids,
    GROUP_CONCAT(CASE WHEN ref_type = 'FDA' THEN ref_id END, ',') AS fda_ref_ids,
    GROUP_CONCAT(CASE WHEN ref_type = 'USAN' THEN ref_id END, ',') AS usan_ref_ids,
    GROUP_CONCAT(CASE WHEN ref_type = 'INN' THEN ref_id END, ',') AS inn_ref_ids
FROM indication_refs
WHERE ref_type != 'ATC'
GROUP BY drugind_id
"""
refs_df = pd.read_sql(refs_query, conn)
print(f"  Aggregated non-ATC references for {len(refs_df):,} unique drugind_ids")

# Attach reference columns to indications
ind_df = ind_df.merge(refs_df, on="drugind_id", how="left")

# ── 3. Fetch molecule hierarchy for bidirectional parent-child remap ─────────
print("Fetching complete molecule hierarchy for parent-child remap...")
hier_query = """
SELECT
    md.chembl_id AS molecule_chembl_id,
    parent_md.chembl_id AS parent_molecule_chembl_id
FROM molecule_hierarchy mh
JOIN molecule_dictionary md ON mh.molregno = md.molregno
JOIN molecule_dictionary parent_md ON mh.parent_molregno = parent_md.molregno
"""
hier_df = pd.read_sql(hier_query, conn)
conn.close()

ind_df["inferred_from_parent"] = False

# Build parent-child remap table
print("Building parent-child remap table...")
hier_pairs = (
    hier_df[["molecule_chembl_id", "parent_molecule_chembl_id"]]
    .dropna(subset=["parent_molecule_chembl_id"])
    .query("molecule_chembl_id != parent_molecule_chembl_id")
    .drop_duplicates()
)
remap_1 = hier_pairs.rename(columns={"molecule_chembl_id": "from_id", "parent_molecule_chembl_id": "to_id"})
remap_2 = hier_pairs.rename(columns={"parent_molecule_chembl_id": "from_id", "molecule_chembl_id": "to_id"})
remap_df = pd.concat([remap_1, remap_2]).drop_duplicates()

def apply_remap(df, remap_df):
    remapped = df.merge(remap_df, left_on="molecule_chembl_id", right_on="from_id", how="inner")
    remapped["molecule_chembl_id"] = remapped["to_id"]
    # Inferred if the new ID is NOT the parent ID
    remapped["inferred_from_parent"] = remapped["molecule_chembl_id"] != remapped["parent_molecule_chembl_id"]
    remapped.drop(columns=["from_id", "to_id"], inplace=True)
    return remapped

print(f"Remapping drug_indication ({len(ind_df):,} rows) via hierarchy...")
ind_remapped = apply_remap(ind_df, remap_df)
all_ind_df = pd.concat([ind_df, ind_remapped]).drop_duplicates(subset=["molecule_chembl_id", "drugind_id"])
print(f"\nCombined indication universe: {len(all_ind_df):,} rows  |  {all_ind_df['molecule_chembl_id'].nunique():,} distinct mols")

# ── 4. Filter to query compounds from union_out.csv ──────────────────────────
print(f"\nLoading query compounds from {INPUT_CSV}...")
union_df = pd.read_csv(INPUT_CSV, low_memory=False)
input_ids = set()
for col in union_df.columns:
    col_str = union_df[col].astype(str)
    for val in col_str.dropna():
        matches = re.findall(r"CHEMBL\d+", val)
        input_ids.update(matches)

input_ids_df = pd.DataFrame({"molecule_chembl_id": list(input_ids)})
print(f"Extracted {len(input_ids_df):,} unique ChEMBL IDs from union_out.csv.")

final_ind_df = all_ind_df.merge(input_ids_df, on="molecule_chembl_id", how="inner")
print(f"\nAfter filtering to union_out compounds: {len(final_ind_df):,} rows")

# ── 5. Clean and format columns ───────────────────────────────────────────────
# Ensure max_phase_for_ind is formatted as string or numeric
final_ind_df["max_phase_for_ind"] = (
    pd.to_numeric(final_ind_df["max_phase_for_ind"], errors="coerce")
    .map(lambda x: f"{x:g}" if pd.notna(x) else None)
)

COL_ORDER = [
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "drugind_id",
    "max_phase_for_ind",
    "mesh_id",
    "mesh_heading",
    "efo_id",
    "efo_term",
    "clinical_trials_ref_ids",
    "daily_med_ref_ids",
    "ema_ref_ids",
    "fda_ref_ids",
    "usan_ref_ids",
    "inn_ref_ids",
    "inferred_from_parent",
]

col_order = [c for c in COL_ORDER if c in final_ind_df.columns]
final_ind_df = final_ind_df[col_order]
final_ind_df = final_ind_df.drop_duplicates()

print(f"\nFinal rows after deduplication: {len(final_ind_df):,}")

# ── 6. Write output CSV ───────────────────────────────────────────────────────
print(f"\nSaving {len(final_ind_df):,} drug indication rows to {OUTPUT_CSV}...")
os.makedirs(OUTPUT_DIR, exist_ok=True)
final_ind_df.to_csv(OUTPUT_CSV, index=False, quoting=1)
print("DONE!")
