
import pandas as pd
import numpy as np

file1 = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/data_extraction/drugs/pubchem/output_data/union/mar-27-2026/union_out.csv"
file2 = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/seeding/seeding_data/apr_29_2026/union_out.csv"

# Columns to audit
audit_cols = ['molecular_weight', 'smiles', 'xlogp', 'molecule_chembl_id', 'chembl_max_phase', 'drug_like']

def load_and_normalize(path, cols):
    # Read only CID and target columns
    df = pd.read_csv(path, usecols=['cid'] + cols, low_memory=False)
    
    # Drop rows with non-numeric CID (like the empty row we found)
    df = df[pd.to_numeric(df['cid'], errors='coerce').notnull()]
    df['cid'] = df['cid'].astype(int)
    
    # Normalize drug_like (boolean/string)
    if 'drug_like' in df.columns:
        df['drug_like'] = df['drug_like'].astype(str).str.upper()
    
    # Normalize numeric columns (rounding to handle float precision diffs)
    numeric_cols = ['molecular_weight', 'xlogp', 'chembl_max_phase']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            # Rounding to 4 decimal places to avoid noise from float formatting
            df[col] = df[col].round(4)
            
    return df.set_index('cid')

print("Loading File 1 (Mar 27)...")
df1 = load_and_normalize(file1, audit_cols)

print("Loading File 2 (Apr 29)...")
df2 = load_and_normalize(file2, audit_cols)

print("\nJoining data...")
comparison = df1.join(df2, lsuffix='_mar27', rsuffix='_apr29')

print("\nAuditing differences:")
for col in audit_cols:
    c1 = f"{col}_mar27"
    c2 = f"{col}_apr29"
    
    # Find rows where values differ (ignoring NaN == NaN)
    diff_mask = (comparison[c1] != comparison[c2]) & ~(comparison[c1].isna() & comparison[c2].isna())
    diff_count = diff_mask.sum()
    
    print(f"Field: {col}")
    print(f"  Differences found: {diff_count}")
    
    if diff_count > 0:
        print(f"  Sample differences (first 5):")
        print(comparison[diff_mask][[c1, c2]].head(5))
        print("-" * 30)

