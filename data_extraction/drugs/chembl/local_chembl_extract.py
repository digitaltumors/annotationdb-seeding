import sqlite3
import pandas as pd
import re
import os

DB_PATH = "/Users/mattbocc/Downloads/chembl_36/chembl_36_sqlite/chembl_36.db"
INPUT_CSV = "/Users/mattbocc/uhn/annotationdb-seeding-testing/data_extraction/drugs/pubchem/output_data/union/mar-27-2026/union_out.csv"
OUTPUT_DIR = "/Users/mattbocc/uhn/annotationdb-seeding-testing/data_extraction/drugs/pubchem/output_data/union/mar-27-2026"
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "chembl_mechanism.csv")

print(f"Connecting to SQLite Database: {DB_PATH}...")
conn = sqlite3.connect(DB_PATH)

print("Fetching ALL drug mechanisms instantly from local DB...")
mech_query = """
SELECT 
    md.chembl_id AS molecule_chembl_id,
    dm.action_type,
    dm.mechanism_of_action,
    td.chembl_id AS target_chembl_id,
    td.pref_name AS target_name,
    dm.record_id,
    dm.mec_id
FROM drug_mechanism dm
JOIN molecule_dictionary md ON dm.molregno = md.molregno
LEFT JOIN target_dictionary td ON dm.tid = td.tid
"""
mech_df = pd.read_sql(mech_query, conn)

print("Fetching ALL molecule parent-child hierarchy mappings from local DB...")
hier_query = """
SELECT 
    child.chembl_id AS child_id,
    parent.chembl_id AS parent_id
FROM molecule_hierarchy mh
JOIN molecule_dictionary child ON mh.molregno = child.molregno
JOIN molecule_dictionary parent ON mh.parent_molregno = parent.molregno
WHERE child.chembl_id != parent.chembl_id
"""
hier_df = pd.read_sql(hier_query, conn)
conn.close()

# Re-map using exact bi-directional logic as the R script
print(f"Performing parent-child bidirectional remap for {len(mech_df)} mechanisms...")
remap_1 = hier_df[['child_id', 'parent_id']].rename(columns={'child_id': 'from_id', 'parent_id': 'to_id'})
remap_2 = hier_df[['parent_id', 'child_id']].rename(columns={'parent_id': 'from_id', 'child_id': 'to_id'})
remap_df = pd.concat([remap_1, remap_2]).drop_duplicates()

mech_remapped = mech_df.merge(remap_df, left_on='molecule_chembl_id', right_on='from_id', how='inner')
mech_remapped['molecule_chembl_id'] = mech_remapped['to_id']
mech_remapped.drop(columns=['from_id', 'to_id'], inplace=True)

# Combine and drop duplicates globally
all_mech_df = pd.concat([mech_df, mech_remapped]).drop_duplicates(subset=['molecule_chembl_id', 'mec_id'])
print(f"Total mechanisms across universe mapped: {len(all_mech_df)}")

# Load needed input IDs from CSV
print(f"Loading query compounds from {INPUT_CSV}...")
union_df = pd.read_csv(INPUT_CSV, low_memory=False)
input_ids = set()
for col in union_df.columns:
    col_str = union_df[col].astype(str)
    for val in col_str.dropna():
        matches = re.findall(r'CHEMBL\d+', val)
        input_ids.update(matches)

input_ids_df = pd.DataFrame({'molecule_chembl_id': list(input_ids)})
print(f"Extracted {len(input_ids_df)} unique ChEMBL IDs from union_out.csv.")

# Join mechanism universe with strictly our queried IDs
final_mech_df = all_mech_df.merge(input_ids_df, on='molecule_chembl_id', how='inner')

# Emulate `clean_mechanism_rows` logical filter
print("Filtering rows missing description properties...")
desc_cols = ['action_type', 'mechanism_of_action', 'target_chembl_id', 'target_name']
for col in desc_cols:
    final_mech_df[col] = final_mech_df[col].fillna('').astype(str).str.strip()

keep_mask = final_mech_df[desc_cols].astype(bool).any(axis=1)
final_mech_df = final_mech_df[keep_mask]
final_mech_df = final_mech_df.drop_duplicates()

# Write strictly to output format
print(f"Saving exactly {len(final_mech_df)} matching mechanism rows to {OUTPUT_CSV}...")
final_mech_df.to_csv(OUTPUT_CSV, index=False, quoting=1) # Quote cleanly
print("DONE! Mechanism extraction complete via SQLite in a fraction of a second!")
