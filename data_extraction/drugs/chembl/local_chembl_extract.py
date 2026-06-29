import sqlite3
import json
import pandas as pd
import re
import os

DB_PATH = "/Users/mattbocc/Downloads/chembl_36/chembl_36_sqlite/chembl_36.db"
# INPUT_CSV = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/data_extraction/drugs/pubchem/substance/output/adcs/substance_out.csv"
# OUTPUT_DIR = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/data_extraction/drugs/pubchem/substance/output/adcs"
INPUT_CSV = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/seeding/seeding_data/jun_29_2026/union_out.csv"
OUTPUT_DIR = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/seeding/seeding_data/jun_29_2026"
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "chembl_mechanism.csv")

print(f"Connecting to SQLite Database: {DB_PATH}...")
conn = sqlite3.connect(DB_PATH)

# ── Source 1: drug_mechanism (curated, high-confidence) ──────────────────────
print("Fetching ALL drug mechanisms from drug_mechanism table...")
mech_query = """
SELECT
    md.chembl_id                  AS molecule_chembl_id,
    parent_md.chembl_id           AS parent_molecule_chembl_id,
    dm.action_type,
    dm.mechanism_of_action,
    dm.mechanism_comment,
    dm.binding_site_comment,
    dm.selectivity_comment,
    dm.direct_interaction,
    dm.disease_efficacy,
    dm.molecular_mechanism,
    md.max_phase,
    md.prodrug,
    parent_md.prodrug             AS parent_prodrug,
    dm.record_id,
    dm.mec_id,
    dm.site_id,
    td.chembl_id                  AS target_chembl_id,
    td.pref_name                  AS target_name,
    -- variant_sequence fields (flat, for internal use; JSON column added below)
    vs.accession                  AS variant_sequence_accession,
    vs.isoform                    AS variant_sequence_isoform,
    vs.mutation                   AS variant_sequence_mutation,
    vs.organism                   AS variant_sequence_organism,
    vs.sequence                   AS variant_sequence_sequence,
    vs.tax_id                     AS variant_sequence_tax_id,
    vs.version                    AS variant_sequence_version,
    CAST(NULL AS INTEGER)         AS activity_id,
    'drug_mechanism'              AS source
FROM drug_mechanism dm
JOIN molecule_dictionary md ON dm.molregno = md.molregno
LEFT JOIN molecule_hierarchy mh ON mh.molregno = md.molregno
LEFT JOIN molecule_dictionary parent_md ON mh.parent_molregno = parent_md.molregno
LEFT JOIN target_dictionary td ON dm.tid = td.tid
LEFT JOIN variant_sequences vs ON dm.variant_id = vs.variant_id
"""
mech_df = pd.read_sql(mech_query, conn)
print(f"  drug_mechanism rows: {len(mech_df):,}  |  distinct mols: {mech_df['molecule_chembl_id'].nunique():,}")

# ── mechanism_refs: aggregate per mec_id → JSON array (matches AnnotationGx) ──
print("Fetching mechanism_refs and aggregating to JSON per mec_id...")
refs_query = """
SELECT mec_id, ref_type, ref_id, ref_url
FROM mechanism_refs
ORDER BY mec_id, ref_type
"""
refs_df = pd.read_sql(refs_query, conn)

def refs_to_json(grp):
    return json.dumps(
        grp[['ref_type', 'ref_id', 'ref_url']].to_dict(orient='records'),
        default=str
    )

refs_json = (
    refs_df.groupby('mec_id')
    .apply(refs_to_json)
    .reset_index()
    .rename(columns={0: 'mechanism_refs'})
)
mech_df = mech_df.merge(refs_json, on='mec_id', how='left')
# Rows with no refs get an empty JSON array (consistent with AnnotationGx nulls)
mech_df['mechanism_refs'] = mech_df['mechanism_refs'].fillna('[]')
print(f"  mechanism_refs attached for {refs_json['mec_id'].nunique():,} unique mec_ids")


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

# ── Build variant_sequence JSON column (matches AnnotationGx's list-column) ───
VS_COLS = {
    'accession': 'variant_sequence_accession',
    'isoform':   'variant_sequence_isoform',
    'mutation':  'variant_sequence_mutation',
    'organism':  'variant_sequence_organism',
    'sequence':  'variant_sequence_sequence',
    'tax_id':    'variant_sequence_tax_id',
    'version':   'variant_sequence_version',
}

def build_variant_sequence_json(row):
    """Serialize variant_sequence flat cols to a JSON object, or null if all empty."""
    d = {k: row[v] for k, v in VS_COLS.items()}
    if all(pd.isna(v) or v == '' for v in d.values()):
        return None
    return json.dumps({k: (None if pd.isna(v) else v) for k, v in d.items()})

mech_df['variant_sequence'] = mech_df.apply(build_variant_sequence_json, axis=1)
mech_df['inferred_from_parent'] = False

# ── Bidirectional parent-child remap ─────────────────────────────────────────
print("Building parent-child remap table...")
hier_pairs = (
    hier_df[['molecule_chembl_id', 'parent_molecule_chembl_id']]
    .dropna(subset=['parent_molecule_chembl_id'])
    .query("molecule_chembl_id != parent_molecule_chembl_id")
    .drop_duplicates()
)
remap_1 = hier_pairs.rename(columns={'molecule_chembl_id': 'from_id', 'parent_molecule_chembl_id': 'to_id'})
remap_2 = hier_pairs.rename(columns={'parent_molecule_chembl_id': 'from_id', 'molecule_chembl_id': 'to_id'})
remap_df = pd.concat([remap_1, remap_2]).drop_duplicates()

def apply_remap(df, remap_df):
    remapped = df.merge(remap_df, left_on='molecule_chembl_id', right_on='from_id', how='inner')
    remapped['molecule_chembl_id'] = remapped['to_id']
    # It should only be inferred if the new ID is NOT the parent ID
    remapped['inferred_from_parent'] = remapped['molecule_chembl_id'] != remapped['parent_molecule_chembl_id']
    remapped.drop(columns=['from_id', 'to_id'], inplace=True)
    return remapped

print(f"Remapping drug_mechanism ({len(mech_df):,} rows) via hierarchy...")
mech_remapped = apply_remap(mech_df, remap_df)
all_mech_df = pd.concat([mech_df, mech_remapped]).drop_duplicates(subset=['molecule_chembl_id', 'mec_id'])


# ── Combine both sources ──────────────────────────────────────────────────────
all_combined = all_mech_df
print(f"\nCombined mechanism universe: {len(all_combined):,} rows  |  {all_combined['molecule_chembl_id'].nunique():,} distinct mols")

# ── Load union_out compound IDs ───────────────────────────────────────────────
print(f"\nLoading query compounds from {INPUT_CSV}...")
union_df = pd.read_csv(INPUT_CSV, low_memory=False)
input_ids = set()
for col in union_df.columns:
    col_str = union_df[col].astype(str)
    for val in col_str.dropna():
        matches = re.findall(r'CHEMBL\d+', val)
        input_ids.update(matches)

input_ids_df = pd.DataFrame({'molecule_chembl_id': list(input_ids)})
print(f"Extracted {len(input_ids_df):,} unique ChEMBL IDs from union_out.csv.")

# ── Filter to our compounds ───────────────────────────────────────────────────
final_mech_df = all_combined.merge(input_ids_df, on='molecule_chembl_id', how='inner')
print(f"\nAfter filtering to union_out compounds: {len(final_mech_df):,} rows")

# ── Filter rows missing all description columns ───────────────────────────────
print("Filtering rows missing all description properties...")
desc_cols = ['action_type', 'mechanism_of_action', 'target_chembl_id', 'target_name']
for col in desc_cols:
    final_mech_df[col] = final_mech_df[col].fillna('').astype(str).str.strip()

keep_mask = final_mech_df[desc_cols].astype(bool).any(axis=1)
final_mech_df = final_mech_df[keep_mask]
final_mech_df = final_mech_df.drop_duplicates()

# Ensure max_phase is an integer
final_mech_df['max_phase'] = pd.to_numeric(final_mech_df['max_phase'], errors='coerce').round().astype('Int64')

# ── Reorder columns to match AnnotationGx output order ───────────────────────
ANNOTATIONGX_COL_ORDER = [
    'molecule_chembl_id',
    'parent_molecule_chembl_id',
    'action_type',
    'binding_site_comment',
    'direct_interaction',
    'disease_efficacy',
    'max_phase',
    'prodrug',
    'parent_prodrug',
    'mec_id',
    'mechanism_comment',
    'mechanism_of_action',
    'mechanism_refs',
    'molecular_mechanism',
    'record_id',
    'selectivity_comment',
    'site_id',
    'target_chembl_id',
    'target_name',
    'variant_sequence',
    # flat variant columns kept for convenience
    'variant_sequence_accession',
    'variant_sequence_isoform',
    'variant_sequence_mutation',
    'variant_sequence_organism',
    'variant_sequence_sequence',
    'variant_sequence_tax_id',
    'variant_sequence_version',
    # extra columns from SQLite
    'activity_id',
    'source',
    'inferred_from_parent',
]
# Only keep columns that exist in the dataframe (safety)
col_order = [c for c in ANNOTATIONGX_COL_ORDER if c in final_mech_df.columns]
final_mech_df = final_mech_df[col_order]

print(f"\nFinal rows after filtering: {len(final_mech_df):,}")

# ── Write output ──────────────────────────────────────────────────────────────
print(f"\nSaving {len(final_mech_df):,} mechanism rows to {OUTPUT_CSV}...")
final_mech_df.to_csv(OUTPUT_CSV, index=False, quoting=1)
print("DONE!")
