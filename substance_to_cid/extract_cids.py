import json
import csv

input_file = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/substance_to_cid/adcdb_results.json"
output_file = "/Users/mattbocc/uhn/annotationdb-seeding-testing copy/substance_to_cid/extracted_cids.csv"

def extract_cids():
    with open(input_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    extracted = []
    seen = set()

    for item in data:
        # Extract payload
        payload_cid_raw = item.get("payload_pubchem_cid")
        payload_name = item.get("payload_name")
        if payload_cid_raw and payload_name:
            # Handle potential comma-separated values just in case
            for cid in str(payload_cid_raw).split(','):
                cid = cid.strip()
                if cid and cid.lower() not in ['none', 'n/a', 'na', 'null'] and (cid, payload_name) not in seen:
                    extracted.append({"cid": cid, "input_id": payload_name})
                    seen.add((cid, payload_name))

        # Extract linker
        linker_cid_raw = item.get("linker_pubchem_cid")
        linker_name = item.get("linker_name")
        if linker_cid_raw and linker_name:
            for cid in str(linker_cid_raw).split(','):
                cid = cid.strip()
                if cid and cid.lower() not in ['none', 'n/a', 'na', 'null'] and (cid, linker_name) not in seen:
                    extracted.append({"cid": cid, "input_id": linker_name})
                    seen.add((cid, linker_name))

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=["cid", "input_id"])
        writer.writeheader()
        writer.writerows(extracted)

    print(f"Successfully extracted {len(extracted)} unique CID/name pairs to {output_file}")

if __name__ == "__main__":
    extract_cids()
