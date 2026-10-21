import json

with open("adcdb_results.json") as f:
    data = json.load(f)

payloads = []
linkers = []
for d in data:
    pc = str(d.get("payload_pubchem_cid", "")).strip()
    lc = str(d.get("linker_pubchem_cid", "")).strip()
    if pc and pc.lower() not in ["none", "na", "n/a", "null"]:
        payloads.append(pc)
    if lc and lc.lower() not in ["none", "na", "n/a", "null"]:
        linkers.append(lc)

print("Total valid payload CIDs:", len(payloads))
print("Unique valid payload CIDs:", len(set(payloads)))
print("Total valid linker CIDs:", len(linkers))
print("Unique valid linker CIDs:", len(set(linkers)))
print("Total combined unique:", len(set(payloads).union(set(linkers))))
