#!/usr/bin/env python3
"""
append_new_compounds.py
-----------------------
Incrementally appends newly scraped compound data to the annotation DB.

Skips any rows whose primary key already exists (INSERT IGNORE semantics),
so it is safe to run multiple times on the same data or after partial failures.

Usage:
    pixi run python seeding/append_new_compounds.py --data_dir <path/to/extraction/output>

    # Example — output produced by input_id_coordinator.R:
    pixi run python seeding/append_new_compounds.py \
        --data_dir data_extraction/drugs/pubchem/kubernetes/output

    # Dry-run (shows row counts without touching the DB):
    pixi run python seeding/append_new_compounds.py \
        --data_dir data_extraction/drugs/pubchem/kubernetes/output \
        --dry-run

Dependency order respected:
    Compounds → CompoundSynonyms → BioAssays → CompoundBioAssays → Toxicity → ChemblMechanism
"""

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path
from urllib.parse import quote_plus

import pandas as pd
from dotenv import load_dotenv
from sqlalchemy import create_engine, text
from sqlalchemy.dialects.mysql import insert as mysql_insert

# ---------------------------------------------------------------------------
# Path resolution — works regardless of cwd
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "seeding"))

from create_tables import (  # noqa: E402  (import after sys.path patch)
    Base,
    BioAssays,
    ChemblMechanism,
    CompoundBioAssays,
    Compounds,
    CompoundSynonyms,
    Toxicity,
)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Append new compound data to annotationdb")
    p.add_argument(
        "--data_dir",
        required=True,
        help="Directory containing the extraction output CSVs "
             "(union_out.csv, union_synonyms.csv, union_bioassays.csv, "
             "union_pubchem_assay_fields.csv, toxicity_output.csv)",
    )
    p.add_argument(
        "--chembl_mech",
        default=None,
        help="Path to chembl_mechanism.csv (optional; skipped if not provided)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print row counts and skip all DB operations",
    )
    p.add_argument(
        "--env",
        default=str(REPO_ROOT / ".env"),
        help="Path to .env file (default: <repo_root>/.env)",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

def make_engine(env_path: str):
    load_dotenv(env_path, override=True)
    password = quote_plus(os.getenv("DATABASE_PASS", ""))
    url = (
        f"mysql+pymysql://{os.getenv('DATABASE_USER')}:{password}"
        f"@{os.getenv('DATABASE_IP')}:{os.getenv('PORT')}/{os.getenv('SELECTED_DB')}"
    )
    return create_engine(url, echo=False)


def align_to_model(df: pd.DataFrame, model) -> pd.DataFrame:
    """Keep only columns that exist in the ORM model."""
    cols = [c.name for c in model.__table__.columns]
    return df[[c for c in cols if c in df.columns]].copy()


def insert_ignore(df: pd.DataFrame, model, engine, label: str, dry_run: bool):
    """
    Insert rows using INSERT IGNORE — rows with duplicate PKs are silently skipped.
    This is safe to run repeatedly on the same data.
    """
    if df.empty:
        print(f"  [{label}] 0 rows — nothing to insert.")
        return

    print(f"  [{label}] {len(df):,} rows to insert (duplicates will be skipped)…")

    if dry_run:
        print(f"  [{label}] [DRY-RUN] Skipped.")
        return

    chunk_size = 500
    inserted = 0
    for start in range(0, len(df), chunk_size):
        chunk = df.iloc[start : start + chunk_size]
        records = chunk.to_dict("records")
        stmt = mysql_insert(model.__table__).values(records).prefix_with("IGNORE")
        with engine.begin() as conn:
            result = conn.execute(stmt)
            inserted += result.rowcount

    print(f"  [{label}] ✓ {inserted:,} new rows inserted ({len(df) - inserted:,} skipped as duplicates).")


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def load_csv(path: Path, label: str) -> pd.DataFrame | None:
    if not path.exists():
        print(f"  [{label}] WARNING: file not found — {path} — skipping.")
        return None
    df = pd.read_csv(path)
    print(f"  [{label}] Loaded {len(df):,} rows from {path.name}")
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    data_dir = Path(args.data_dir)
    dry_run = args.dry_run

    if not data_dir.is_dir():
        print(f"ERROR: --data_dir is not a directory: {data_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"  append_new_compounds.py  [{datetime.now():%Y-%m-%d %H:%M:%S}]")
    print(f"  data_dir : {data_dir}")
    print(f"  dry_run  : {dry_run}")
    print(f"{'='*60}\n")

    # --- Load CSVs ---
    print("── Loading CSVs ────────────────────────────────────────────")
    compounds_df    = load_csv(data_dir / "union_out.csv",                 "compounds")
    synonyms_df     = load_csv(data_dir / "union_synonyms.csv",            "synonyms")
    bioassays_df    = load_csv(data_dir / "union_pubchem_assay_fields.csv","bioassays")
    cb_df           = load_csv(data_dir / "union_bioassays.csv",           "compound_bioassays")
    toxicity_df     = load_csv(data_dir / "toxicity_output.csv",           "toxicity")

    chembl_path = Path(args.chembl_mech) if args.chembl_mech else None
    chembl_df = load_csv(chembl_path, "chembl_mechanism") if chembl_path else None

    if compounds_df is None:
        print("ERROR: union_out.csv is required — cannot continue.", file=sys.stderr)
        sys.exit(1)

    # --- Clean compounds ---
    print("\n── Cleaning & aligning ─────────────────────────────────────")

    # Deduplicate on CID (keep first)
    if "cid" in compounds_df.columns:
        before = len(compounds_df)
        compounds_df = compounds_df.drop_duplicates(subset=["cid"], keep="first")
        print(f"  [compounds] Deduped {before - len(compounds_df):,} duplicate CIDs.")

    # Derive fda_approval from chembl_max_phase
    if "chembl_max_phase" in compounds_df.columns:
        compounds_df["fda_approval"] = compounds_df["chembl_max_phase"] == 4

    compounds_df = align_to_model(compounds_df, Compounds)

    # --- Clean synonyms ---
    if synonyms_df is not None:
        if {"synonym", "pubchem_cid"}.issubset(synonyms_df.columns):
            before = len(synonyms_df)
            synonyms_df = synonyms_df.drop_duplicates(subset=["synonym", "pubchem_cid"], keep="first")
            print(f"  [synonyms] Deduped {before - len(synonyms_df):,} duplicate rows.")
        synonyms_df = align_to_model(synonyms_df, CompoundSynonyms)

    # --- Clean bioassays ---
    if bioassays_df is not None:
        if "aid" in bioassays_df.columns:
            bioassays_df["aid"] = pd.to_numeric(bioassays_df["aid"], errors="coerce")
            bioassays_df = bioassays_df.dropna(subset=["aid"])
            bioassays_df["aid"] = bioassays_df["aid"].astype(int)
            bioassays_df = bioassays_df.drop_duplicates(subset=["aid"], keep="first")
        bioassays_df = align_to_model(bioassays_df, BioAssays)

    # --- Clean compound_bioassays ---
    if cb_df is not None:
        if "bioassay_aid" in cb_df.columns:
            cb_df["bioassay_aid"] = pd.to_numeric(cb_df["bioassay_aid"], errors="coerce")
            cb_df = cb_df.dropna(subset=["bioassay_aid"])
            cb_df["bioassay_aid"] = cb_df["bioassay_aid"].astype(int)

        if "pubchem_cid" in cb_df.columns:
            cb_df["pubchem_cid"] = pd.to_numeric(cb_df["pubchem_cid"], errors="coerce")
            cb_df = cb_df.dropna(subset=["pubchem_cid"])
            cb_df["pubchem_cid"] = cb_df["pubchem_cid"].astype(int)

        # FK safety: only keep rows whose CID exists in compounds being inserted
        valid_cids = set(compounds_df["cid"].dropna().astype(int))
        before = len(cb_df)
        cb_df = cb_df[cb_df["pubchem_cid"].isin(valid_cids)]
        if before != len(cb_df):
            print(f"  [compound_bioassays] Dropped {before - len(cb_df):,} rows with unknown CIDs.")

        if {"bioassay_aid", "pubchem_cid"}.issubset(cb_df.columns):
            cb_df = cb_df.drop_duplicates(subset=["bioassay_aid", "pubchem_cid"], keep="first")

        cb_df = align_to_model(cb_df, CompoundBioAssays)

    # --- Clean toxicity ---
    if toxicity_df is not None:
        if "pubchem_cid" in toxicity_df.columns:
            toxicity_df["pubchem_cid"] = pd.to_numeric(toxicity_df["pubchem_cid"], errors="coerce")
            toxicity_df = toxicity_df.dropna(subset=["pubchem_cid"])
            toxicity_df["pubchem_cid"] = toxicity_df["pubchem_cid"].astype(int)

        valid_cids = set(compounds_df["cid"].dropna().astype(int))
        before = len(toxicity_df)
        toxicity_df = toxicity_df[toxicity_df["pubchem_cid"].isin(valid_cids)]
        if before != len(toxicity_df):
            print(f"  [toxicity] Dropped {before - len(toxicity_df):,} rows with unknown CIDs.")

        toxicity_df = toxicity_df.drop_duplicates(subset=["pubchem_cid"], keep="first")
        toxicity_df = align_to_model(toxicity_df, Toxicity)

    # --- Clean chembl mechanism ---
    if chembl_df is not None:
        valid_chembls = set(compounds_df["molecule_chembl_id"].dropna().astype(str)) \
            if "molecule_chembl_id" in compounds_df.columns else set()
        if "source" in chembl_df.columns:
            chembl_df = chembl_df[chembl_df["source"] == "drug_mechanism"]
        if valid_chembls and "molecule_chembl_id" in chembl_df.columns:
            before = len(chembl_df)
            chembl_df = chembl_df[chembl_df["molecule_chembl_id"].astype(str).isin(valid_chembls)]
            print(f"  [chembl_mechanism] Dropped {before - len(chembl_df):,} orphan rows.")
        chembl_df = align_to_model(chembl_df, ChemblMechanism)

    # --- Connect & ensure tables exist ---
    print("\n── Connecting to DB ─────────────────────────────────────────")
    if dry_run:
        print("  [DRY-RUN] Skipping DB connection.")
        engine = None
    else:
        engine = make_engine(args.env)
        Base.metadata.create_all(engine)
        print("  Connected ✓")

    # --- Insert in FK dependency order ---
    print("\n── Inserting (INSERT IGNORE) ────────────────────────────────")
    insert_ignore(compounds_df,             Compounds,          engine, "compounds",           dry_run)
    insert_ignore(synonyms_df or pd.DataFrame(), CompoundSynonyms, engine, "synonyms",         dry_run)
    insert_ignore(bioassays_df or pd.DataFrame(), BioAssays,    engine, "bioassays",            dry_run)
    insert_ignore(cb_df or pd.DataFrame(), CompoundBioAssays,   engine, "compound_bioassays",   dry_run)
    insert_ignore(toxicity_df or pd.DataFrame(), Toxicity,      engine, "toxicity",             dry_run)
    if chembl_df is not None:
        insert_ignore(chembl_df,            ChemblMechanism,    engine, "chembl_mechanism",     dry_run)

    print(f"\n{'='*60}")
    print(f"  ✅  Done  [{datetime.now():%Y-%m-%d %H:%M:%S}]")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
