import os
import time
import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

INPUT_CSV = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/seeding/seeding_data/jun_23_2026/union_out.csv"

# Configure requests session with retries and exponential backoff
session = requests.Session()
retries = Retry(
    total=5,
    backoff_factor=2,  # 2s, 4s, 8s, 16s, 32s
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET"]
)
adapter = HTTPAdapter(max_retries=retries)
session.mount("https://", adapter)
session.mount("http://", adapter)

def fetch_chembl_ids_bulk(cids: list) -> dict:
    """
    Hits the bulk PubChem PUG REST API to retrieve Registry IDs for a list of CIDs.
    Returns a dictionary mapping CID (int) -> ChEMBL ID (str).
    """
    cids_str = ",".join(map(str, cids))
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids_str}/xrefs/RegistryID/JSON"
    
    results = {}
    try:
        response = session.get(url, timeout=30)
        
        # 404 is normal if none of the requested CIDs have any xrefs at all
        if response.status_code == 404:
            return results
            
        response.raise_for_status()
        data = response.json()
        
        info_list = data.get("InformationList", {}).get("Information", [])
        
        for item in info_list:
            cid = item.get("CID")
            registry_ids = item.get("RegistryID", [])
            for rid in registry_ids:
                if str(rid).startswith("CHEMBL"):
                    results[cid] = str(rid)
                    break
                                    
    except Exception as e:
        print(f"Error fetching chunk of {len(cids)} CIDs: {e}")
        
    return results

def main():
    print(f"Loading {INPUT_CSV}...")
    df = pd.read_csv(INPUT_CSV, low_memory=False)
    
    # Identify rows that have a CID but no ChEMBL ID
    if 'molecule_chembl_id' not in df.columns:
        df['molecule_chembl_id'] = None
        
    missing_mask = (df['molecule_chembl_id'].isna() | (df['molecule_chembl_id'] == '')) & df['cid'].notna()
    cids_to_fetch = df.loc[missing_mask, 'cid'].dropna().astype(int).unique()
    
    print(f"Found {len(cids_to_fetch)} unique CIDs missing a ChEMBL ID.")
    if len(cids_to_fetch) == 0:
        print("Nothing to do!")
        return

    print("Beginning bulk PubChem API extraction (100 CIDs per request)...")
    
    found_count = 0
    checked_count = 0
    chunk_size = 100
    
    for i in range(0, len(cids_to_fetch), chunk_size):
        chunk = cids_to_fetch[i:i+chunk_size]
        chembl_mappings = fetch_chembl_ids_bulk(chunk)
        
        for cid, chembl_id in chembl_mappings.items():
            df.loc[df['cid'] == cid, 'molecule_chembl_id'] = chembl_id
            found_count += 1
        
        checked_count += len(chunk)
        print(f"Progress: {checked_count}/{len(cids_to_fetch)} queried. Found: {found_count} new ChEMBL IDs.")
            
        if (i // chunk_size) > 0 and (i // chunk_size) % 10 == 0:
            print("Checkpointing progress to CSV...")
            df.to_csv(INPUT_CSV, index=False)
            
        # Stutter query: PubChem allows max 5 requests per second. 
        # 0.22s sleep ensures we are safely under the limit.
        time.sleep(0.22)
        
    # Final save
    print(f"\nFinished! Found {found_count} new ChEMBL IDs out of {len(cids_to_fetch)} queried.")
    print(f"Saving final results to {INPUT_CSV}...")
    df.to_csv(INPUT_CSV, index=False)
    print("Done.")

if __name__ == "__main__":
    main()
