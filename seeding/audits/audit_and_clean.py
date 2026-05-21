import pandas as pd
import os
import shutil

def main():
    data_dir = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/seeding/seeding_data/may_15_2026"
    
    union_out_path = os.path.join(data_dir, "union_out.csv")
    toxicity_path = os.path.join(data_dir, "toxicity_output.csv")
    synonyms_path = os.path.join(data_dir, "union_synonyms.csv")
    bioassays_path = os.path.join(data_dir, "union_bioassays.csv")
    
    print(f"Reading {union_out_path}...")
    df_union = pd.read_csv(union_out_path, dtype=str)
    
    initial_count = len(df_union)
    
    # Identify rows missing a title (NaN or empty string)
    missing_title_mask = df_union['title'].isna() | (df_union['title'].str.strip() == "")
    cids_to_drop = set(df_union.loc[missing_title_mask, 'cid'])
    
    print(f"Found {len(cids_to_drop)} compounds without a title.")
    
    if len(cids_to_drop) == 0:
        print("No compounds to delete. Exiting.")
        return

    # Filter union_out
    df_union_cleaned = df_union[~missing_title_mask]
    final_count = len(df_union_cleaned)
    print(f"union_out.csv: Dropped {initial_count - final_count} rows. Original: {initial_count}, Final: {final_count}")
    df_union_cleaned.to_csv(union_out_path, index=False)
    
    # Clean toxicity (small enough to fit in memory)
    if os.path.exists(toxicity_path):
        print(f"\nCleaning {toxicity_path}...")
        df_tox = pd.read_csv(toxicity_path, dtype=str)
        initial_tox = len(df_tox)
        df_tox_cleaned = df_tox[~df_tox['pubchem_cid'].isin(cids_to_drop)]
        final_tox = len(df_tox_cleaned)
        print(f"toxicity_output.csv: Dropped {initial_tox - final_tox} rows. Original: {initial_tox}, Final: {final_tox}")
        df_tox_cleaned.to_csv(toxicity_path, index=False)
        
    # Clean synonyms using chunking (2M+ rows)
    if os.path.exists(synonyms_path):
        print(f"\nCleaning {synonyms_path} (chunked)...")
        temp_syn = synonyms_path + ".tmp"
        dropped_syn = 0
        total_syn = 0
        first_chunk = True
        for chunk in pd.read_csv(synonyms_path, dtype=str, chunksize=100000):
            initial_chunk_len = len(chunk)
            total_syn += initial_chunk_len
            chunk_cleaned = chunk[~chunk['pubchem_cid'].isin(cids_to_drop)]
            dropped_syn += (initial_chunk_len - len(chunk_cleaned))
            
            chunk_cleaned.to_csv(temp_syn, mode='a', index=False, header=first_chunk)
            first_chunk = False
            
        shutil.move(temp_syn, synonyms_path)
        print(f"union_synonyms.csv: Dropped {dropped_syn} rows. Original: {total_syn}, Final: {total_syn - dropped_syn}")

    # Clean bioassays using chunking (24M+ rows)
    if os.path.exists(bioassays_path):
        print(f"\nCleaning {bioassays_path} (chunked)...")
        temp_bio = bioassays_path + ".tmp"
        dropped_bio = 0
        total_bio = 0
        first_chunk = True
        for chunk in pd.read_csv(bioassays_path, dtype=str, chunksize=500000):
            initial_chunk_len = len(chunk)
            total_bio += initial_chunk_len
            chunk_cleaned = chunk[~chunk['pubchem_cid'].isin(cids_to_drop)]
            dropped_bio += (initial_chunk_len - len(chunk_cleaned))
            
            chunk_cleaned.to_csv(temp_bio, mode='a', index=False, header=first_chunk)
            first_chunk = False
            
        shutil.move(temp_bio, bioassays_path)
        print(f"union_bioassays.csv: Dropped {dropped_bio} rows. Original: {total_bio}, Final: {total_bio - dropped_bio}")

    print("\n✅ Audit and cleanup complete!")

if __name__ == "__main__":
    main()
