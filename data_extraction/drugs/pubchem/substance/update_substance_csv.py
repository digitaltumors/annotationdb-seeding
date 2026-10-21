import csv
import re
import urllib.request
import json
import ssl
import time

input_file = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/data_extraction/drugs/pubchem/substance/output/combined/substance_toxicity copy.csv"
output_file = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/data_extraction/drugs/pubchem/substance/output/combined/substance_toxicity.csv"

def fetch_livertox_url_by_sid(sid):
    ctx = ssl._create_unverified_context()
    # 1. Get RefChem ID
    url1 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/JSON"
    try:
        req = urllib.request.Request(url1, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=10, context=ctx) as resp:
            data = json.loads(resp.read().decode())
            syns = data.get("PC_Substances", [])[0].get("synonyms", [])
            # extract refchem id (first digit-only synonym)
            refchem_id = next((s for s in syns if s.isdigit()), None)
            if not refchem_id:
                return ""
            
            # 2. Get Pug View data
            url2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/refchem/{refchem_id}/JSON/"
            req2 = urllib.request.Request(url2, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req2, timeout=10, context=ctx) as resp2:
                data2 = json.loads(resp2.read().decode())
                record = data2.get("Record", {})
                for ref in record.get("Reference", []):
                    if ref.get("SourceName") == "LiverTox":
                        return ref.get("URL", "")
    except Exception as e:
        pass
    return ""

with open(input_file, 'r', newline='') as f_in:
    reader = csv.DictReader(f_in)
    rows = list(reader)

fieldnames = [
    "sid", "reference_number", "tox_dataset",
    "dili_severity_grade", "dili_annotation", "tox_source_url",
    "hepatotoxicity_likelihood_score", "hepatotoxicity_likelihood_score_reasoning"
]

url_cache = {}
seen_hep_sids = set()
hep_rows_to_write = []

with open(output_file, 'w', newline='') as f_out:
    writer = csv.DictWriter(f_out, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
    writer.writeheader()

    for row in rows:
        sid = row.get("sid", "")
        hep_str = row.get("hepatotoxicity_likelihood_score", "")

        parsed_score = ""
        parsed_reasoning = ""

        if hep_str and "likelihood score:" in hep_str.lower():
            m = re.search(r"Likelihood score:\s*([A-Za-z0-9\*\[\]]+)\s*\((.*?)(?:\)|$)", hep_str, re.IGNORECASE)
            if m:
                parsed_score = m.group(1).strip()
                parsed_reasoning = m.group(2).strip()
            else:
                parsed_score = hep_str.strip()
            
            # Save the Livertox row for later
            if sid not in seen_hep_sids:
                seen_hep_sids.add(sid)
                
                # Fetch URL if we haven't
                if sid not in url_cache:
                    url_cache[sid] = fetch_livertox_url_by_sid(sid)
                    time.sleep(0.3)
                
                hep_rows_to_write.append({
                    "sid": sid,
                    "reference_number": "",
                    "tox_dataset": "Livertox",
                    "dili_severity_grade": "",
                    "dili_annotation": "",
                    "tox_source_url": url_cache[sid],
                    "hepatotoxicity_likelihood_score": parsed_score,
                    "hepatotoxicity_likelihood_score_reasoning": parsed_reasoning
                })

        # Clear the hep fields from the current row
        row["hepatotoxicity_likelihood_score"] = ""
        row["hepatotoxicity_likelihood_score_reasoning"] = ""
        
        # Rename dataset and source_url
        dataset = row.get("dili_dataset", "").strip()
        row["tox_dataset"] = dataset
        row["tox_source_url"] = row.get("dili_source_url", "")

        # If it's an empty row (tox_dataset is empty/NA), we SKIP it.
        if not dataset or dataset.lower() == "na":
            continue
        
        filtered_row = {k: v for k, v in row.items() if k in fieldnames}
        writer.writerow(filtered_row)
        
    for hep_row in hep_rows_to_write:
        writer.writerow(hep_row)

print("Substance CSV rewrite complete!")
