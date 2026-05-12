import pandas as pd
import sqlite3
import os

CSV_PATH = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/data_extraction/drugs/pubchem/output_data/union/mar-27-2026/union_out.csv"
DB_PATH = "/Users/mattbocc/Downloads/chembl_36/chembl_36_sqlite/chembl_36.db"

print(f"Loading {CSV_PATH}...")
df = pd.read_csv(CSV_PATH, low_memory=False)

total_rows = len(df)
print(f"Total rows: {total_rows:,}")

# Rows with ChEMBL ID
with_chembl = df[df['molecule_chembl_id'].notna() & (df['molecule_chembl_id'] != '')]
print(f"Rows with ChEMBL ID: {len(with_chembl):,}")

# Missing max_phase
missing_max_phase = with_chembl[with_chembl['chembl_max_phase'].isna()]
print(f"Rows with ChEMBL ID but missing chembl_max_phase: {len(missing_max_phase):,}")

if len(missing_max_phase) > 0:
    print("Checking if any of these missing ones have a max_phase in ChEMBL DB...")
    conn = sqlite3.connect(DB_PATH)
    
    # Get all max phases from ChEMBL
    chembl_ids = missing_max_phase['molecule_chembl_id'].unique().tolist()
    # Batch query for speed
    placeholders = ', '.join(['?'] * len(chembl_ids))
    query = f"SELECT chembl_id, max_phase FROM molecule_dictionary WHERE chembl_id IN ({placeholders}) AND max_phase IS NOT NULL"
    
    found_in_db = pd.read_sql(query, conn, params=chembl_ids)
    conn.close()
    
    print(f"Found {len(found_in_db):,} compounds in DB that have a max_phase but are missing it in CSV.")
    if len(found_in_db) > 0:
        print("\nExamples of missing data:")
        print(found_in_db.head(10))
else:
    print("No ChEMBL IDs are missing max_phase values in the CSV.")
