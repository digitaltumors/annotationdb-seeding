#!/usr/bin/env python3
"""
extract_refchem_toxicity.py
===========================
Batch-query PubChem for a list of compound names:
    compound name → SIDs → find the SID with a RefChem synonym → RefChem view
    → extract toxicity + synonyms → combined CSVs.

Output files:
    substance_out.csv   – sid, title, molecule_chembl_id, chembl_max_phase
    substance_toxicity.csv  – sid, dili_severity_grade, dili_annotation,
                    hepatotoxicity_likelihood_score, reference_number,
                    dili_dataset, dili_source_url
    substance_synonyms.csv  – synonym, sid

Usage:
    python substance_extraction.py --outdir ./output
    python substance_extraction.py --sid 481101883          # single SID mode
"""

import argparse
import csv
import json
import os
import re
import ssl
import sys
import time
import urllib.request
import urllib.error
import urllib.parse

# ---------------------------------------------------------------------------
# Compound list
# ---------------------------------------------------------------------------
COMPOUNDS = [
    "pertuzumab",
    "trastuzumab",
    "durvalumab",
    "panitumumab",
    "bevacizumab",
    "ipilimumab",
    "nivolumab",
    "trastuzumab deruxtecan",
    "sacituzumab govitecan",
    "daratumumab",
    "elotuzumab",
    "isatuximab",
    "rituximab",
    "obinutuzumab",
    "ofatumumab",
    "brentuximab vedotin",
    "polatuzumab vedotin",
    "ibritumomab tiuxetan",
    "aldesleukin",
    "talimogene laherparepvec",
    "ramucirumab",
    "aflibercept",
    "cetuximab",
    "amivantamab",
    "dinutuximab",
    "denileukin diftitox",
    "ado-trastuzumab emtansine",
    "fam-trastuzumab deruxtecan",
    "trifluridine",
    "dostarlimab",
    "olaratumab",
    "blinatumomab",
    "eculizumab",
]

# ---------------------------------------------------------------------------
# Canonical field-name look-ups (case-insensitive matching)
# ---------------------------------------------------------------------------
SEVERITY_NAMES = {
    "severity class",        # DILIrank 2.0
    "severity grade",        # older DILIrank
    "dili severity score",   # LTKB benchmark
    "dili_severity_score",
    "dili severity grade",
}

ANNOTATION_NAMES = {
    "vdili-concern",           # DILIrank 2.0 (v-prefixed)
    "dili-concern",
    "dili annotation",         # older DILIrank
    "dili concern",            # LTKB benchmark
    "dilist classification",   # DILIst dataset
    "dili label",
    "concern",
}

# ---------------------------------------------------------------------------
# Rate limiter  (PubChem: max 5 requests / second)
# ---------------------------------------------------------------------------

class RateLimiter:
    """Sliding-window rate limiter.

    Ensures no more than *max_requests* are made within any *window* second
    period.  Call :meth:`wait` before every HTTP request.
    """

    def __init__(self, max_requests: int = 5, window: float = 1.0):
        self.max_requests = max_requests
        self.window = window
        self._timestamps: list[float] = []

    def wait(self):
        """Block until it is safe to make the next request."""
        now = time.time()
        # Purge timestamps outside the window
        self._timestamps = [t for t in self._timestamps if now - t < self.window]
        if len(self._timestamps) >= self.max_requests:
            sleep_for = self.window - (now - self._timestamps[0]) + 0.05
            if sleep_for > 0:
                time.sleep(sleep_for)
        self._timestamps.append(time.time())


# Global rate limiter – shared by every call to fetch_json
_limiter = RateLimiter(max_requests=4, window=1.0)   # stay under the 5/s cap

# Seconds to pause between compounds (lets PubChem cool down)
COMPOUND_DELAY = 2.0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Global unverified SSL context to bypass macOS/local certificate errors
try:
    _ssl_context = ssl._create_unverified_context()
except AttributeError:
    _ssl_context = None

def fetch_json(url: str, retries: int = 5, delay: float = 2.0) -> dict | None:
    """Fetch JSON from *url* with retry + rate-limit logic.

    Returns None on HTTP 404 (resource does not exist).
    Automatically retries on 429/503 (throttle) with exponential backoff.
    """
    for attempt in range(retries):
        _limiter.wait()
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30, context=_ssl_context) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            # PubChem throttle: 429 or 503
            if exc.code in (429, 503):
                wait = min(delay * (2 ** attempt), 60)
                print(f"    ⚠ Throttled (HTTP {exc.code}), waiting {wait:.0f}s …")
                time.sleep(wait)
                continue
            if attempt < retries - 1:
                wait = min(delay * (2 ** attempt), 60)
                print(f"    ⚠ HTTP {exc.code} on attempt {attempt + 1}, retrying in {wait:.0f}s …")
                time.sleep(wait)
            else:
                raise
        except (urllib.error.URLError, ConnectionError) as exc:
            if attempt < retries - 1:
                wait = min(delay * (2 ** attempt), 60)
                print(f"    ⚠ Network error on attempt {attempt + 1} ({exc}), retrying in {wait:.0f}s …")
                time.sleep(wait)
            else:
                raise


def extract_refchem_id(synonyms: list[str]) -> str | None:
    """Return the numeric RefChem ID from a synonym like 'RefChem:543'."""
    for syn in synonyms:
        m = re.match(r"(?i)refchem[:\s]+(\d+)", syn)
        if m:
            return m.group(1)
    return None


# ---------------------------------------------------------------------------
# Name → SID resolution
# ---------------------------------------------------------------------------

def resolve_name_to_sids(name: str) -> list[int]:
    """Query PubChem substance-name endpoint; return list of SIDs."""
    encoded = urllib.parse.quote(name, safe="")
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/"
        f"{encoded}/sids/JSON"
    )
    data = fetch_json(url)
    if data is None:
        return []
    return data.get("IdentifierList", {}).get("SID", [])


def find_refchem_sid(sids: list[int]) -> tuple[int | None, str | None, list[str], int | None, list[str]]:
    """
    Iterate SIDs and return the first one whose synonyms contain a RefChem ID.
    Also returns the first valid SID and its synonyms as a fallback.

    Returns (sid, refchem_id, synonyms, fallback_sid, fallback_synonyms)
    """
    fallback_sid = None
    fallback_synonyms = []

    for sid in sids:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/JSON"
        data = fetch_json(url)          # rate-limited internally
        if data is None:
            continue
        substances = data.get("PC_Substances", [])
        if not substances:
            continue
        synonyms = substances[0].get("synonyms", [])

        if fallback_sid is None:
            fallback_sid = sid
            fallback_synonyms = synonyms

        refchem_id = extract_refchem_id(synonyms)
        if refchem_id is not None:
            return sid, refchem_id, synonyms, fallback_sid, fallback_synonyms
    return None, None, [], fallback_sid, fallback_synonyms


# ---------------------------------------------------------------------------
# ChEMBL max_phase lookup
# ---------------------------------------------------------------------------

def fetch_chembl_max_phase(chembl_id: str) -> str | None:
    """Query the ChEMBL API for the molecule's max_phase."""
    if not chembl_id:
        return None
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    data = fetch_json(url)
    if data is None:
        return None
    return str(data.get("max_phase", "")) or None


# ---------------------------------------------------------------------------
# Recursive section walker
# ---------------------------------------------------------------------------

def walk_sections(sections: list[dict]):
    """Yield every section dict in a depth-first traversal."""
    for sec in sections:
        yield sec
        if "Section" in sec:
            yield from walk_sections(sec["Section"])


# ---------------------------------------------------------------------------
# Extraction helpers for the RefChem pug_view response
# ---------------------------------------------------------------------------

def _get_string_values(info_block: dict) -> list[str]:
    """Pull all String values from a single Information entry."""
    val = info_block.get("Value", {})
    strings = []
    for swm in val.get("StringWithMarkup", []):
        s = swm.get("String", "").strip()
        if s:
            strings.append(s)
    return strings


def extract_chembl_id(record: dict) -> str | None:
    """Find the ChEMBL ID from 'Other Identifiers' → 'ChEMBL ID' section."""
    for sec in walk_sections(record.get("Section", [])):
        if sec.get("TOCHeading") == "ChEMBL ID":
            for info in sec.get("Information", []):
                vals = _get_string_values(info)
                if vals:
                    return vals[0]
    return None


def extract_synonyms(record: dict) -> list[str]:
    """Return the list of synonyms from the top-level Synonyms section."""
    for sec in record.get("Section", []):
        if sec.get("TOCHeading") == "Synonyms":
            synonyms = []
            for info in sec.get("Information", []):
                synonyms.extend(_get_string_values(info))
            return synonyms
    return []


def _build_reference_map(record: dict) -> dict[int, dict]:
    """Map ReferenceNumber → reference metadata (source name, URL, etc.)."""
    ref_map = {}
    for ref in record.get("Reference", []):
        rn = ref.get("ReferenceNumber")
        if rn is not None:
            ref_map[rn] = ref
    return ref_map


def extract_dili_data(record: dict) -> list[dict]:
    """
    Walk the record and extract every DILI-related information block.

    For each DILI dataset entry, produce a dict with:
        - reference_number
        - dili_dataset
        - dili_severity_grade
        - dili_annotation
        - dili_source_url
    """
    ref_map = _build_reference_map(record)
    results: list[dict] = []

    for sec in walk_sections(record.get("Section", [])):
        if sec.get("TOCHeading") != "Drug Induced Liver Injury":
            continue

        # Group Information entries by ReferenceNumber (each = one dataset)
        groups: dict[int, list[dict]] = {}
        for info in sec.get("Information", []):
            rn = info.get("ReferenceNumber")
            groups.setdefault(rn, []).append(info)

        for rn, info_list in groups.items():
            row = {
                "reference_number": rn,
                "dili_dataset": None,
                "dili_severity_grade": None,
                "dili_annotation": None,
                "dili_source_url": None,
            }

            # Resolve source URL from the Reference block
            ref_meta = ref_map.get(rn, {})
            row["dili_source_url"] = ref_meta.get("URL")

            for info in info_list:
                name = (info.get("Name") or "").strip()
                name_lower = name.lower()
                vals = _get_string_values(info)
                val = vals[0] if vals else None

                if name_lower == "dataset":
                    row["dili_dataset"] = val
                elif name_lower in SEVERITY_NAMES:
                    row["dili_severity_grade"] = val
                elif name_lower in ANNOTATION_NAMES:
                    row["dili_annotation"] = val

            results.append(row)

    return results


def extract_hepatotoxicity(record: dict) -> str | None:
    """
    Extract the hepatotoxicity likelihood score from the Hepatotoxicity section.
    Matches any StringWithMarkup containing 'Likelihood score:'.
    """
    for sec in walk_sections(record.get("Section", [])):
        if sec.get("TOCHeading") != "Hepatotoxicity":
            continue
        for info in sec.get("Information", []):
            for s in _get_string_values(info):
                if "likelihood score:" in s.lower():
                    return s.strip()
    return None


# ---------------------------------------------------------------------------
# Per-compound pipeline
# ---------------------------------------------------------------------------

def process_compound(compound_name: str, refchem_id: str) -> dict | None:
    """
    Fetch the RefChem pug_view data and extract all fields.

    Returns a dict with keys:
        record_title, chembl_id, synonyms, dili_rows, hepatotox_score
    """
    refchem_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/refchem/"
        f"{refchem_id}/JSON/"
    )
    refchem_data = fetch_json(refchem_url)
    if refchem_data is None:
        print(f"    ⚠ No RefChem view data for ID={refchem_id}")
        return None

    record = refchem_data.get("Record", {})
    record_title = record.get("RecordTitle", compound_name)

    return {
        "record_title": record_title,
        "chembl_id": extract_chembl_id(record),
        "synonyms": extract_synonyms(record),
        "dili_rows": extract_dili_data(record),
        "hepatotox_score": extract_hepatotoxicity(record),
    }


# ---------------------------------------------------------------------------
# Single-SID mode (backwards compatible)
# ---------------------------------------------------------------------------

def run_single_sid(sid: int, outdir: str):
    """Original single-SID mode for backwards compatibility."""
    os.makedirs(outdir, exist_ok=True)

    substance_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/JSON"
    )
    print(f"[1/4] Fetching substance SID={sid} …")
    substance_data = fetch_json(substance_url)
    if substance_data is None:
        print(f"ERROR: No substance data found for SID {sid}.", file=sys.stderr)
        sys.exit(1)

    substances = substance_data.get("PC_Substances", [])
    if not substances:
        print(f"ERROR: No substance data found for SID {sid}.", file=sys.stderr)
        sys.exit(1)

    synonyms_raw = substances[0].get("synonyms", [])
    refchem_id = extract_refchem_id(synonyms_raw)
    if refchem_id is None:
        print(f"ERROR: No RefChem ID found in synonyms for SID {sid}.", file=sys.stderr)
        sys.exit(1)

    print(f"       → RefChem ID: {refchem_id}")

    print(f"[2/4] Fetching RefChem view for ID={refchem_id} …")
    result = process_compound("unknown", refchem_id)
    if result is None:
        sys.exit(1)

    chembl_id = result["chembl_id"]
    print(f"       → Record: {result['record_title']} (#{refchem_id})")

    print(f"[3/4] Fetching ChEMBL max_phase …")
    max_phase = fetch_chembl_max_phase(chembl_id)
    print(f"       → max_phase: {max_phase}")

    print("[4/4] Extracting fields …")
    print(f"       → ChEMBL ID           : {chembl_id}")
    print(f"       → Synonyms count      : {len(result['synonyms'])}")
    print(f"       → DILI dataset entries : {len(result['dili_rows'])}")
    print(f"       → Hepatotox likelihood : {result['hepatotox_score']}")

    # -- substance_out.csv --
    sid_file = os.path.join(outdir, "substance_out.csv")
    with open(sid_file, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["sid", "mapped_name", "title", "molecule_chembl_id", "chembl_max_phase"])
        w.writeheader()
        w.writerow({
            "sid": sid, "mapped_name": "unknown", "title": result["record_title"],
            "molecule_chembl_id": chembl_id, "chembl_max_phase": max_phase,
        })
    print(f"  ✓ substance_out.csv    → {sid_file}")

    # -- substance_toxicity.csv --
    tox_file = os.path.join(outdir, "substance_toxicity.csv")
    tox_fields = [
        "sid", "dili_severity_grade", "dili_annotation",
        "hepatotoxicity_likelihood_score", "reference_number",
        "dili_dataset", "dili_source_url",
    ]
    with open(tox_file, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=tox_fields)
        w.writeheader()
        dili_rows = result["dili_rows"] or [{}]
        for row in dili_rows:
            w.writerow({
                "sid": sid,
                "dili_severity_grade": row.get("dili_severity_grade"),
                "dili_annotation": row.get("dili_annotation"),
                "hepatotoxicity_likelihood_score": result["hepatotox_score"],
                "reference_number": row.get("reference_number"),
                "dili_dataset": row.get("dili_dataset"),
                "dili_source_url": row.get("dili_source_url"),
            })
    print(f"  ✓ substance_toxicity.csv   → {tox_file}")

    # -- substance_synonyms.csv --
    syn_file = os.path.join(outdir, "substance_synonyms.csv")
    with open(syn_file, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["synonym", "sid"])
        w.writeheader()
        for syn in result["synonyms"]:
            w.writerow({"synonym": syn, "sid": sid})
    print(f"  ✓ substance_synonyms.csv   → {syn_file}")

    print("\nDone.")


# ---------------------------------------------------------------------------
# Batch mode
# ---------------------------------------------------------------------------

def run_batch(outdir: str, compounds_data: list = None):
    """Process compounds from compounds_data list, or default COMPOUNDS list.
    compounds_data should be a list of dicts: {"name": str, "sids": list | None}
    """
    os.makedirs(outdir, exist_ok=True)

    if compounds_data is None:
        compounds_data = [{"name": c, "sids": None} for c in COMPOUNDS]

    sid_fields = ["sid", "mapped_name", "title", "molecule_chembl_id", "chembl_max_phase"]
    tox_fields = [
        "sid", "dili_severity_grade", "dili_annotation",
        "hepatotoxicity_likelihood_score", "reference_number",
        "dili_dataset", "dili_source_url",
    ]
    syn_fields = ["synonym", "sid"]
    log_fields = [
        "compound_name", "status", "sids_found", "sid_used", "refchem_id", "error",
    ]

    sid_file = os.path.join(outdir, "substance_out.csv")
    tox_file = os.path.join(outdir, "substance_toxicity.csv")
    syn_file = os.path.join(outdir, "substance_synonyms.csv")
    log_file = os.path.join(outdir, "resolution_log.csv")

    total = len(COMPOUNDS)
    success = 0
    failed = 0

    with (
        open(sid_file, "w", newline="", encoding="utf-8") as f_sid,
        open(tox_file, "w", newline="", encoding="utf-8") as f_tox,
        open(syn_file, "w", newline="", encoding="utf-8") as f_syn,
        open(log_file, "w", newline="", encoding="utf-8") as f_log,
    ):
        w_sid = csv.DictWriter(f_sid, fieldnames=sid_fields)
        w_tox = csv.DictWriter(f_tox, fieldnames=tox_fields)
        w_syn = csv.DictWriter(f_syn, fieldnames=syn_fields)
        w_log = csv.DictWriter(f_log, fieldnames=log_fields)
        w_sid.writeheader()
        w_tox.writeheader()
        w_syn.writeheader()
        w_log.writeheader()

        total = len(compounds_data)
        for idx, item in enumerate(compounds_data, 1):
            compound = item["name"]
            given_sids = item.get("sids")

            print(f"\n{'='*60}")
            print(f"[{idx}/{total}] {compound}")
            print(f"{'='*60}")

            # Cool down between compounds
            if idx > 1:
                time.sleep(COMPOUND_DELAY)

            # ----- Step 1: Name → SIDs or use given SIDs -----
            if given_sids is not None:
                print(f"  → Using provided SIDs ({len(given_sids)} total) …")
                sids = given_sids
            else:
                print(f"  → Resolving name to SIDs …")
                try:
                    sids = resolve_name_to_sids(compound)
                except Exception as exc:
                    print(f"    ✗ Failed to resolve name: {exc}")
                    w_log.writerow({
                        "compound_name": compound, "status": "FAILED",
                        "sids_found": 0, "sid_used": "", "refchem_id": "",
                        "error": f"Name resolution failed: {exc}",
                    })
                    failed += 1
                    continue

            if not sids:
                print(f"    ✗ No SIDs found for '{compound}'")
                w_log.writerow({
                    "compound_name": compound, "status": "NO_SIDS",
                    "sids_found": 0, "sid_used": "", "refchem_id": "",
                    "error": "No substance IDs returned by PubChem",
                })
                failed += 1
                continue

            print(f"    Found {len(sids)} SID(s)")

            # ----- Step 2: Find the SID with a RefChem synonym -----
            print(f"  → Scanning SIDs for RefChem ID …")
            sid, refchem_id, _, fallback_sid, fallback_synonyms = find_refchem_sid(sids)

            if sid is None or refchem_id is None:
                if fallback_sid is not None:
                    print(f"    ✗ No RefChem SID found. Falling back to SID={fallback_sid} for basic info.")
                    chembl_id = next((s for s in fallback_synonyms if s.startswith("CHEMBL")), None)
                    result = {
                        "record_title": compound,
                        "chembl_id": chembl_id,
                        "synonyms": fallback_synonyms,
                        "dili_rows": [],
                        "hepatotox_score": None,
                    }
                    sid = fallback_sid
                    w_log.writerow({
                        "compound_name": compound, "status": "NO_REFCHEM_FALLBACK_USED",
                        "sids_found": len(sids), "sid_used": sid, "refchem_id": "",
                        "error": "None of the SIDs had a RefChem synonym, used first valid SID",
                    })
                else:
                    print(f"    ✗ No RefChem SID found and no fallback available.")
                    w_log.writerow({
                        "compound_name": compound, "status": "NO_REFCHEM_AND_NO_FALLBACK",
                        "sids_found": len(sids), "sid_used": "", "refchem_id": "",
                        "error": "None of the SIDs had a RefChem synonym and could not fetch fallback",
                    })
                    failed += 1
                    continue
            else:
                print(f"    ✓ SID={sid}, RefChem ID={refchem_id}")

                # ----- Step 3: Fetch RefChem view data -----
                print(f"  → Fetching RefChem view …")
                try:
                    result = process_compound(compound, refchem_id)
                except Exception as exc:
                    print(f"    ✗ RefChem view failed: {exc}")
                    w_log.writerow({
                        "compound_name": compound, "status": "REFCHEM_ERROR",
                        "sids_found": len(sids), "sid_used": sid,
                        "refchem_id": refchem_id,
                        "error": f"RefChem view fetch failed: {exc}",
                    })
                    failed += 1
                    continue

                if result is None:
                    w_log.writerow({
                        "compound_name": compound, "status": "REFCHEM_EMPTY",
                        "sids_found": len(sids), "sid_used": sid,
                        "refchem_id": refchem_id,
                        "error": "RefChem view returned no data",
                    })
                    failed += 1
                    continue

            # ----- Step 4: Fetch ChEMBL max_phase -----
            chembl_id = result["chembl_id"]
            print(f"  → Fetching ChEMBL max_phase for {chembl_id} …")
            max_phase = fetch_chembl_max_phase(chembl_id)

            record_title = result["record_title"]
            hepatotox = result["hepatotox_score"]
            dili_rows = result["dili_rows"]
            synonyms = result["synonyms"]

            print(f"    ChEMBL ID           : {chembl_id}")
            print(f"    ChEMBL max_phase    : {max_phase}")
            print(f"    Synonyms            : {len(synonyms)}")
            print(f"    DILI dataset entries : {len(dili_rows)}")
            print(f"    Hepatotox score     : {hepatotox}")

            # ----- Write substance_out row -----
            w_sid.writerow({
                "sid": sid,
                "mapped_name": compound,
                "title": record_title,
                "molecule_chembl_id": chembl_id,
                "chembl_max_phase": max_phase,
            })

            # ----- Write toxicity rows -----
            if dili_rows:
                for row in dili_rows:
                    w_tox.writerow({
                        "sid": sid,
                        "dili_severity_grade": row["dili_severity_grade"],
                        "dili_annotation": row["dili_annotation"],
                        "hepatotoxicity_likelihood_score": hepatotox,
                        "reference_number": row["reference_number"],
                        "dili_dataset": row["dili_dataset"],
                        "dili_source_url": row["dili_source_url"],
                    })
            else:
                w_tox.writerow({
                    "sid": sid,
                    "dili_severity_grade": None,
                    "dili_annotation": None,
                    "hepatotoxicity_likelihood_score": hepatotox,
                    "reference_number": None,
                    "dili_dataset": None,
                    "dili_source_url": None,
                })

            # ----- Write synonym rows -----
            for syn in synonyms:
                w_syn.writerow({"synonym": syn, "sid": sid})

            # ----- Log -----
            w_log.writerow({
                "compound_name": compound, "status": "OK",
                "sids_found": len(sids), "sid_used": sid,
                "refchem_id": refchem_id, "error": "",
            })
            success += 1

            # Flush after each compound so partial results are saved
            f_sid.flush()
            f_tox.flush()
            f_syn.flush()
            f_log.flush()

    print(f"\n{'='*60}")
    print(f"BATCH COMPLETE: {success} succeeded, {failed} failed out of {total}")
    print(f"{'='*60}")
    print(f"  ✓ substance_out.csv       → {sid_file}")
    print(f"  ✓ substance_toxicity.csv      → {tox_file}")
    print(f"  ✓ substance_synonyms.csv      → {syn_file}")
    print(f"  ✓ resolution_log.csv → {log_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract RefChem toxicity and synonym data from PubChem.",
    )
    parser.add_argument(
        "--sid",
        type=int,
        default=None,
        help="Single PubChem Substance ID (SID) to query.  "
             "If omitted, runs in batch mode over the built-in compound list.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="./output",
        help="Directory for output CSV files (default: ./output).",
    )
    parser.add_argument(
        "--json",
        type=str,
        default=None,
        help="Path to adcdb_results.json to read adc_puchem_sid from instead of hardcoded list.",
    )
    args = parser.parse_args()

    if args.json:
        with open(args.json, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        compounds_data = []
        for item in data:
            name = item.get("adc_drug_name") or item.get("adc_name") or "Unknown"
            sid_str = item.get("adc_puchem_sid", "")
            if sid_str and sid_str.strip():
                try:
                    sids = [int(s.strip()) for s in sid_str.split(",") if s.strip()]
                    compounds_data.append({"name": name, "sids": sids})
                except ValueError:
                    print(f"    ⚠ Warning: Could not parse SIDs '{sid_str}' for {name}, falling back to name lookup.")
                    compounds_data.append({"name": name, "sids": None})
            else:
                compounds_data.append({"name": name, "sids": None})
        
        run_batch(args.outdir, compounds_data)
    elif args.sid is not None:
        run_single_sid(args.sid, args.outdir)
    else:
        run_batch(args.outdir)


if __name__ == "__main__":
    main()
