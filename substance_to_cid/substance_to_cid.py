"""
Attempt to map mixed PubChem identifiers/synonyms to PubChem CIDs.

Lookup strategy per ID:
        1) Try Substance SID -> CID
        2) Compound name -> CID
        3) Substance name -> CID
        4) Substance name -> SID(s) -> CID(s)
        5) NEW: Substance xref/RegistryID -> SID(s) -> CID(s) (can catch "External IDs" that are not indexed as name synonyms)
"""

import csv
import time
import random
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote
import requests

# Constants and CWD relatives

INPUT_FILE = "inquiry_compounds.csv"
MAPPED_OUT = "mapped_drug_ids_to_cids.csv"
UNMAPPED_OUT = "unmapped_drug_ids.csv"

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

BATCH_SIZE = 50
MAX_RETRIES = 6

MAX_SIDS_PER_ID = 25


# Helpers for reads, sleeps, and chunking data


def _sleep_backoff(attempt: int) -> None:
    base = min(2**attempt, 32)
    time.sleep(base + random.random())


def read_ids_from_tsv(path: str) -> List[str]:
    """
    Reads IDs from the first column of a TSV.
    Deduplicates while preserving order.

    Heuristic: if first row first cell is non-numeric, treat as header.
    """
    ids: List[str] = []
    seen = set()

    with open(path, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if not rows:
        return []

    start_idx = 0
    if rows[0] and rows[0][0]:
        first_cell = rows[0][0].strip()
        if not first_cell.isdigit():
            start_idx = 1

    for row in rows[start_idx:]:
        if not row:
            continue
        raw = (row[0] or "").strip()
        if not raw:
            continue

        parts = raw.replace(";", ",").replace("|", ",").split(",")
        for p in parts:
            p = p.strip()
            if p and p not in seen:
                seen.add(p)
                ids.append(p)

    return ids


def chunked(lst: List[str], n: int) -> List[List[str]]:
    return [lst[i : i + n] for i in range(0, len(lst), n)]


def _dedupe_preserve_order_int(values: List[int]) -> List[int]:
    seen = set()
    out: List[int] = []
    for v in values:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


# HTTP + PubChem parsing


def _get_json(
    session: requests.Session, url: str
) -> Tuple[Optional[dict], Optional[str]]:
    """
    Returns (json, error_reason). error_reason is None on success.
    Retries on common transient errors.
    """
    for attempt in range(MAX_RETRIES):
        try:
            r = session.get(url, timeout=30)

            if r.status_code == 200:
                return r.json(), None

            if r.status_code == 404:
                return None, "404 (not found)"

            if r.status_code == 400:
                return None, "400 (bad request)"

            if r.status_code in (429, 500, 502, 503, 504):
                _sleep_backoff(attempt)
                continue

            return None, f"HTTP {r.status_code}"

        except (requests.Timeout, requests.ConnectionError):
            _sleep_backoff(attempt)
            continue
        except ValueError:
            return None, "Invalid JSON"

    return None, "Retries exhausted"


def _extract_cids(data: dict) -> List[int]:
    info = data.get("InformationList", {}).get("Information", [])
    cids: List[int] = []
    for item in info:
        cid_val = item.get("CID", [])
        if isinstance(cid_val, int):
            cid_val = [cid_val]
        if isinstance(cid_val, list):
            for c in cid_val:
                try:
                    cids.append(int(c))
                except Exception:
                    pass
    return _dedupe_preserve_order_int(cids)


def _extract_sids(data: dict) -> List[int]:
    info = data.get("InformationList", {}).get("Information", [])
    sids: List[int] = []
    for item in info:
        sid_val = item.get("SID", [])
        if isinstance(sid_val, int):
            sid_val = [sid_val]
        if isinstance(sid_val, list):
            for s in sid_val:
                try:
                    sids.append(int(s))
                except Exception:
                    pass
    return _dedupe_preserve_order_int(sids)


# PubChem lookups


def lookup_sid_to_cids(
    session: requests.Session, sid: str
) -> Tuple[Optional[List[int]], str]:
    url = f"{PUBCHEM_BASE}/substance/sid/{sid}/cids/JSON"
    data, err = _get_json(session, url)
    if not data:
        return None, f"sid->cids failed: {err}"
    cids = _extract_cids(data)
    if not cids:
        return None, "sid->cids returned no CIDs"
    return cids, "ok"


def lookup_compound_name_to_cids(
    session: requests.Session, name: str
) -> Tuple[Optional[List[int]], str]:
    enc = quote(name, safe="")
    url = f"{PUBCHEM_BASE}/compound/name/{enc}/cids/JSON"
    data, err = _get_json(session, url)
    if not data:
        return None, f"compound name->cids failed: {err}"
    cids = _extract_cids(data)
    if not cids:
        return None, "compound name->cids returned no CIDs"
    return cids, "ok"


def lookup_substance_name_to_cids(
    session: requests.Session, name: str
) -> Tuple[Optional[List[int]], str]:
    enc = quote(name, safe="")
    url = f"{PUBCHEM_BASE}/substance/name/{enc}/cids/JSON"
    data, err = _get_json(session, url)
    if not data:
        return None, f"substance name->cids failed: {err}"
    cids = _extract_cids(data)
    if not cids:
        return None, "substance name->cids returned no CIDs"
    return cids, "ok"


def lookup_substance_name_to_sids(
    session: requests.Session, name: str
) -> Tuple[Optional[List[int]], str]:
    enc = quote(name, safe="")
    url = f"{PUBCHEM_BASE}/substance/name/{enc}/sids/JSON"
    data, err = _get_json(session, url)
    if not data:
        return None, f"substance name->sids failed: {err}"
    sids = _extract_sids(data)
    if not sids:
        return None, "substance name->sids returned no SIDs"
    return sids, "ok"


def lookup_substance_xref_registryid_to_sids(
    session: requests.Session, registry_id: str
) -> Tuple[Optional[List[int]], str]:
    """
    substance/xref/RegistryID/<registry_id>/sids

    This is the key fallback for PubChem "External IDs" that don't resolve via the name namespace.
    :contentReference[oaicite:3]{index=3}
    """
    rid = quote(registry_id, safe="")
    url = f"{PUBCHEM_BASE}/substance/xref/RegistryID/{rid}/sids/JSON"
    data, err = _get_json(session, url)
    if not data:
        return None, f"xref/RegistryID->sids failed: {err}"
    sids = _extract_sids(data)
    if not sids:
        return None, "xref/RegistryID->sids returned no SIDs"
    return sids, "ok"


# Mapping logic


def map_one_id(
    session: requests.Session, identifier: str
) -> Tuple[Optional[List[int]], str, str]:
    """
    Returns:
      (cids or None, method_used, reason_if_unmapped)
    """
    ident = identifier.strip()

    # 1) If numeric, try as SID first
    sid_attempt_msg = ""
    if ident.isdigit():
        cids, sid_attempt_msg = lookup_sid_to_cids(session, ident)
        if cids:
            return cids, "substance/sid", ""

    # 2) Compound name/synonym lookup
    cids, msg_comp = lookup_compound_name_to_cids(session, ident)
    if cids:
        return cids, "compound/name", ""

    # 3) Substance name lookup (direct CID)
    cids, msg_sub_cids = lookup_substance_name_to_cids(session, ident)
    if cids:
        return cids, "substance/name->cids", ""

    # 4) Substance name -> SIDs -> CIDs
    sids, msg_sub_sids = lookup_substance_name_to_sids(session, ident)
    if sids:
        all_cids: List[int] = []
        for sid in sids[:MAX_SIDS_PER_ID]:
            c, _ = lookup_sid_to_cids(session, str(sid))
            if c:
                all_cids.extend(c)
        uniq_cids = _dedupe_preserve_order_int(all_cids)
        if uniq_cids:
            return uniq_cids, "substance/name->sids->cids", ""

    # 5) NEW: xref/RegistryID -> SIDs -> CIDs
    sids2, msg_xref = lookup_substance_xref_registryid_to_sids(session, ident)
    if sids2:
        all_cids2: List[int] = []
        for sid in sids2[:MAX_SIDS_PER_ID]:
            c, _ = lookup_sid_to_cids(session, str(sid))
            if c:
                all_cids2.extend(c)
        uniq2 = _dedupe_preserve_order_int(all_cids2)
        if uniq2:
            return uniq2, "substance/xref/RegistryID->sids->cids", ""
        return (
            None,
            "substance/xref/RegistryID->sids->cids",
            "Found SID(s) but no CID mapping",
        )

    # Unmapped: combine reasons
    reason_parts = []
    if ident.isdigit() and sid_attempt_msg:
        reason_parts.append(sid_attempt_msg)
    reason_parts.extend([msg_comp, msg_sub_cids, msg_sub_sids, msg_xref])

    reason = "; ".join([p for p in reason_parts if p]) or "No mapping returned"
    return None, "none", reason


# Main coordinator


def main():
    ids = read_ids_from_tsv(INPUT_FILE)
    if not ids:
        raise SystemExit(f"No IDs found in {INPUT_FILE}")

    session = requests.Session()
    session.headers.update({"User-Agent": "matthew.boccalon@uhn.ca"})

    results: Dict[str, dict] = {}

    chunks = chunked(ids, BATCH_SIZE)
    for i, ch in enumerate(chunks, start=1):
        for ident in ch:
            cids, method, reason = map_one_id(session, ident)

            results[ident] = {
                "cids": cids or [],
                "primary": cids[0] if cids else "",
                "status": "mapped" if cids else "unmapped",
                "method": method,
                "reason": reason,
            }

            time.sleep(0.2)

        processed = min(i * BATCH_SIZE, len(ids))
        mapped_n = sum(1 for v in results.values() if v["status"] == "mapped")
        print(f"[{i}/{len(chunks)}] processed={processed}/{len(ids)} mapped={mapped_n}")

    with open(MAPPED_OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            ["input_id", "status", "method", "mapped_cids", "primary_cid", "reason"]
        )
        for ident in ids:
            r = results[ident]
            w.writerow(
                [
                    ident,
                    r["status"],
                    r["method"],
                    ";".join(map(str, r["cids"])) if r["cids"] else "",
                    r["primary"],
                    r["reason"] if r["status"] == "unmapped" else "",
                ]
            )

    with open(UNMAPPED_OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["input_id", "reason"])
        for ident in ids:
            r = results[ident]
            if r["status"] == "unmapped":
                w.writerow([ident, r["reason"]])

    print("\nDone.")
    print(f"  → {MAPPED_OUT}")
    print(f"  → {UNMAPPED_OUT}")


if __name__ == "__main__":
    main()
