import sqlite3
import pandas as pd
import numpy as np

INPUT_CSV = "seeding/seeding_data/jun_23_2026/union_out.csv"
DB_PATH = "/Users/mattbocc/Downloads/chembl_36/chembl_36_sqlite/chembl_36.db"

print(f"Loading {INPUT_CSV}...")
df = pd.read_csv(INPUT_CSV, low_memory=False)

# Get ChEMBL IDs that are present but missing max_phase
if 'chembl_max_phase' not in df.columns:
    df['chembl_max_phase'] = np.nan
    
mask = df['molecule_chembl_id'].notna() & df['chembl_max_phase'].isna()
missing_phases = df.loc[mask, 'molecule_chembl_id'].unique().tolist()

if not missing_phases:
    print("No missing max_phases found for the mapped ChEMBL IDs.")
else:
    print(f"Found {len(missing_phases)} ChEMBL IDs missing max_phase. Querying local SQLite DB...")
    
    conn = sqlite3.connect(DB_PATH)
    
    found_mappings = []
    chunk_size = 900
    for i in range(0, len(missing_phases), chunk_size):
        chunk = missing_phases[i:i+chunk_size]
        placeholders = ','.join(['?'] * len(chunk))
        query = f"""
            SELECT chembl_id AS molecule_chembl_id, max_phase AS new_max_phase
            FROM molecule_dictionary
            WHERE chembl_id IN ({placeholders})
        """
        found_mappings.append(pd.read_sql(query, conn, params=chunk))
        
    mapping_df = pd.concat(found_mappings).drop_duplicates(subset=['molecule_chembl_id'])
    conn.close()
    
    print(f"Retrieved max_phase for {len(mapping_df)} ChEMBL IDs.")
    
    # Merge and update
    df = df.merge(mapping_df, on='molecule_chembl_id', how='left')
    df['chembl_max_phase'] = df['chembl_max_phase'].fillna(df['new_max_phase'])
    df = df.drop(columns=['new_max_phase'])
    
    print("Saving updated CSV...")
    df.to_csv(INPUT_CSV, index=False)
    print("Done!")
