import pandas as pd

# ── Load union_out as the bridge between CID <-> InChIKey ──────────────────
union_out = pd.read_csv("data/output/cloud_union_out.csv", encoding="utf-8", low_memory=False)
union_out["inchikey"] = union_out["inchikey"].astype(str).str.strip().str.upper()
# Normalize float CIDs like "8197.0" -> "8197"
union_out["cid"] = union_out["cid"].apply(
    lambda x: str(int(float(x))) if str(x).strip() not in ("", "nan", "NA") else ""
)

# Build lookup sets
cid_to_ikey  = dict(zip(union_out["cid"], union_out["inchikey"]))
ikey_to_cid  = dict(zip(union_out["inchikey"], union_out["cid"]))

def to_inchikeys_via_cid(cids):
    """Map a series of CIDs -> inchikeys using union_out."""
    return set(
        cid_to_ikey[c]
        for c in cids.astype(str).str.strip()
        if c in cid_to_ikey
    ) - {"NAN", "NA", ""}

def to_inchikeys_direct(ikeys):
    """Normalise a series of raw InChIKeys."""
    return set(ikeys.astype(str).str.strip().str.upper()) - {"NAN", "NA", ""}

# ── Load each input dataset and resolve to InChIKey sets ──────────────────
oasis = pd.read_csv("data/input/oasis_compounds_needed.csv", encoding="utf-8")
sider = pd.read_csv("data/input/SIDER_needed.csv",           encoding="utf-8")
tox   = pd.read_csv("data/input/toxcast_needed.csv",         encoding="utf-8")
geom  = pd.read_csv("data/input/GEOM_needed_full.csv",       encoding="utf-8")
jump  = pd.read_csv("data/input/DFM-JUMP-CP-Training-Combined.csv", encoding="utf-8")

datasets = {
    "OASIS":    to_inchikeys_direct(oasis["inchikey"]),
    "SIDER":    to_inchikeys_via_cid(sider["cid"]),
    "ToxCast":  to_inchikeys_via_cid(tox["cid"]),
    "GEOM":     to_inchikeys_via_cid(geom["cid"]),
    "DFM-JUMP": to_inchikeys_direct(jump["InChiKey"]),
}

names = list(datasets.keys())

# ── Print sizes ────────────────────────────────────────────────────────────
print("\n=== Dataset sizes (resolved to InChIKey via union_out) ===")
for n, s in datasets.items():
    print(f"  {n:<12} {len(s):>7,} unique InChIKeys")

# ── Pairwise overlap matrix ────────────────────────────────────────────────
print(f"\n=== Pairwise overlap (# shared InChIKeys) ===")
col_w = 10
header = f"{'':12}" + "".join(f"{n:>{col_w}}" for n in names)
print(header)

for a in names:
    row = f"{a:<12}"
    for b in names:
        overlap = len(datasets[a] & datasets[b]) if a != b else "-"
        row += f"{str(overlap):>{col_w}}"
    print(row)

# ── Compounds in ALL 5 datasets ────────────────────────────────────────────
all_five = set.intersection(*datasets.values())
print(f"\n=== Compounds in ALL 5 datasets: {len(all_five):,} ===")

# ── Compounds unique to each dataset (not in any other) ───────────────────
print("\n=== Compounds unique to each dataset (not in any other) ===")
for a in names:
    others = set.union(*(datasets[b] for b in names if b != a))
    unique = datasets[a] - others
    print(f"  {a:<12} {len(unique):>7,} unique compounds")


# ════════════════════════════════════════════════════════════════════════════
# RAW OVERLAP — direct comparison using native identifiers (all input rows,
# regardless of whether they have been scraped into union_out yet)
# ════════════════════════════════════════════════════════════════════════════

def norm_cids(series):
    """Normalise a CID series to plain integer strings."""
    result = set()
    for v in series:
        v = str(v).strip()
        if v in ("", "nan", "NA"):
            continue
        try:
            result.add(str(int(float(v))))
        except ValueError:
            pass
    return result

def norm_ikeys(series):
    return set(series.astype(str).str.strip().str.upper()) - {"NAN", "NA", ""}

# Build a full CID<->InChIKey bridge from BOTH union_out AND all input CIDs
# so we can convert cross-type pairs
all_input_cids = norm_cids(
    pd.concat([sider["cid"], tox["cid"], geom["cid"]])
)
# Map any input CID that union_out knows about
raw_cid_to_ikey = {c: cid_to_ikey[c] for c in all_input_cids if c in cid_to_ikey}

def cids_to_ikeys(cids_set):
    return {raw_cid_to_ikey[c] for c in cids_set if c in raw_cid_to_ikey}

# Raw identifier sets (native units)
raw = {
    "OASIS":    norm_ikeys(oasis["inchikey"]),
    "SIDER":    norm_cids(sider["cid"]),
    "ToxCast":  norm_cids(tox["cid"]),
    "GEOM":     norm_cids(geom["cid"]),
    "DFM-JUMP": norm_ikeys(jump["InChiKey"]),
}

# Identifier type per dataset
id_type = {
    "OASIS": "ikey", "SIDER": "cid", "ToxCast": "cid",
    "GEOM": "cid", "DFM-JUMP": "ikey",
}

def raw_overlap(a, b):
    """Overlap between two datasets on their native identifiers."""
    ta, tb = id_type[a], id_type[b]
    sa, sb = raw[a], raw[b]
    if ta == tb:
        return len(sa & sb)
    # cross-type: convert CID side to InChIKey via union_out bridge
    if ta == "cid":
        sa = cids_to_ikeys(sa)
    else:
        sb = cids_to_ikeys(sb)
    return len(sa & sb)

print("\n\n=== RAW INPUT overlap (total compounds in input files, not filtered to union_out) ===")
print("=== Dataset sizes (raw input identifiers) ===")
for n, s in raw.items():
    print(f"  {n:<12} {len(s):>8,}  ({id_type[n]})")

print(f"\n=== Raw pairwise overlap matrix ===")
col_w = 10
print(f"{'':12}" + "".join(f"{n:>{col_w}}" for n in names))
for a in names:
    row = f"{a:<12}"
    for b in names:
        val = "-" if a == b else raw_overlap(a, b)
        row += f"{str(val):>{col_w}}"
    print(row)

all_ikeys = {
    "OASIS":    raw["OASIS"],
    "SIDER":    cids_to_ikeys(raw["SIDER"]),
    "ToxCast":  cids_to_ikeys(raw["ToxCast"]),
    "GEOM":     cids_to_ikeys(raw["GEOM"]),
    "DFM-JUMP": raw["DFM-JUMP"],
}
all_five_raw = set.intersection(*all_ikeys.values())
print(f"\n=== Compounds in ALL 5 raw datasets: {len(all_five_raw):,} ===")

print("\n=== Unique to each dataset (not in any other, raw) ===")
for a in names:
    others = set.union(*(all_ikeys[b] for b in names if b != a))
    unique = all_ikeys[a] - others
    print(f"  {a:<12} {len(unique):>8,} unique compounds")