import csv
import os
import sys

# Increase CSV limits for massive files
csv.field_size_limit(sys.maxsize)

BUCKET = "gs://annotationdb_data_retrieval"
OUT_CSV = "audit_master_union_out.csv"

# 1. Download union_out.csv
print(f"Downloading master union_out.csv from GCS to local...", flush=True)
os.system(f'gsutil -o "GSUtil:parallel_process_count=1" cp {BUCKET}/output/union_out.csv {OUT_CSV}')

# 2. Extract mappings from union_out
# We will use the row index as the universal 'union_id'
completed_input_ids = {}
completed_inchikeys = {}
completed_cids = {}

union_total_rows = 0
with open(OUT_CSV, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader):
        union_total_rows += 1
        union_id = i
        
        iid = row.get('input_id', '').strip()
        if iid and iid != 'NA':
            completed_input_ids[iid] = union_id
            
        ikey = row.get('inchikey', '').strip()
        if ikey and ikey != 'NA':
            completed_inchikeys[ikey] = union_id

        cid = row.get('cid', '').strip()
        if not cid:
            cid = row.get('CID', '').strip()
        if cid and cid not in ('NA', '0', ''):
            try:
                completed_cids[str(int(float(cid)))] = union_id
            except ValueError:
                pass

print(f"\nLoaded {union_total_rows} total valid rows from union_out.csv\n", flush=True)

# 3. List of input files
inputs = [
    "DFM-JUMP-CP-Training-Combined.csv",
    "GEOM_needed_full.csv",
    "LINCS_needed.csv",
    "SIDER_needed.csv",
    "TOXCAST_needed.csv",
    "oasis_compounds_needed.csv"
]

total_missing_overall = 0
total_inputs_overall = 0

dataset_union_matches = {}

for file_name in inputs:
    local_file = f"audit_{file_name}"
    print(f"--- Checking {file_name} ---", flush=True)
    os.system(f'gsutil -o "GSUtil:parallel_process_count=1" cp {BUCKET}/input/{file_name} {local_file}')
    
    matched_union_ids = set()
    
    with open(local_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        fieldnames_lower = {f.lower(): f for f in fieldnames}
        
        input_id_col = None
        inchikey_col = None
        
        for candidate in ['input_id', 'jcp2022 id', 'broad id', 'oasis_id', 'id']:
            if candidate in fieldnames_lower:
                input_id_col = fieldnames_lower[candidate]
                break
                
        for candidate in ['inchikey', 'inchi_key']:
            if candidate in fieldnames_lower:
                inchikey_col = fieldnames_lower[candidate]
                break
        if not inchikey_col:
            for f_name in fieldnames:
                if f_name.lower() in ('inchikey', 'inchi_key', 'inchi key'):
                    inchikey_col = f_name
                    break

        cid_col = None
        for candidate in ['cid', 'pubchemcid', 'pubchem_cid', 'pubchem cid']:
            if candidate in fieldnames_lower:
                cid_col = fieldnames_lower[candidate]
                break
                
        if not input_id_col and not inchikey_col and not cid_col:
            print(f"WARNING: Could not identify ID column in {file_name}. Columns: {fieldnames}")
            dataset_union_matches[file_name] = set()
            continue

        print(f"  Detected: input_id={input_id_col}, inchikey={inchikey_col}, cid={cid_col}", flush=True)
            
        total_rows = 0
        missing_rows = 0
        
        for row in reader:
            iid = row.get(input_id_col, '').strip() if input_id_col else ''
            ikey = row.get(inchikey_col, '').strip() if inchikey_col else ''
            raw_cid = row.get(cid_col, '').strip() if cid_col else ''
            try:
                cid = str(int(float(raw_cid))) if raw_cid and raw_cid not in ('NA', '0') else ''
            except ValueError:
                cid = ''
            
            if not iid and not ikey and not cid:
                continue
                
            total_rows += 1
            
            union_id = None
            if iid and iid in completed_input_ids:
                union_id = completed_input_ids[iid]
            elif ikey and ikey in completed_inchikeys:
                union_id = completed_inchikeys[ikey]
            elif cid and cid in completed_cids:
                union_id = completed_cids[cid]
                
            if union_id is not None:
                matched_union_ids.add(union_id)
            else:
                missing_rows += 1
                
        found_rows = total_rows - missing_rows
        pct = (found_rows / total_rows * 100) if total_rows > 0 else 100
        
        print(f"  Total valid entries: {total_rows}")
        print(f"  Found in union_out: {found_rows} ({pct:.5f}%)")
        print(f"  Missing: {missing_rows}\n")
        
        total_missing_overall += missing_rows
        total_inputs_overall += total_rows
        dataset_union_matches[file_name] = matched_union_ids

print("===============================================")
pct_overall = ((total_inputs_overall - total_missing_overall) / total_inputs_overall * 100)
print(f"OVERALL MISSING: {total_missing_overall} out of {total_inputs_overall} total compounds.")
print(f"OVERALL COVERAGE: {pct_overall:.5f}%")
print("===============================================\n")

print("===============================================")
print("OVERLAP MATRIX (Number of shared compounds mapped in union_out.csv)")
print("===============================================")

# Format matrix table
short_names = [name.replace('.csv', '').replace('-needed', '').replace('_needed', '').replace('_full', '').replace('-CP-Training-Combined', '') for name in inputs]
header = [""] + short_names
print("\t".join(f"{h:>12}" for h in header))

for i, f1 in enumerate(inputs):
    row_strs = [f"{short_names[i]:>12}"]
    for f2 in inputs:
        s1 = dataset_union_matches.get(f1, set())
        s2 = dataset_union_matches.get(f2, set())
        overlap_size = len(s1 & s2)
        row_strs.append(f"{overlap_size:>12}")
    print("\t".join(row_strs))

print("===============================================\n")

# Cleanup downloaded audit files
os.remove(OUT_CSV)
for f in inputs:
    try:
        os.remove(f"audit_{f}")
    except:
        pass
