from sqlalchemy import create_engine
import os
from datetime import datetime
from dotenv import load_dotenv
from urllib.parse import quote_plus
import pandas as pd

from create_tables import (
    Base,
    Compounds,
    CompoundSynonyms,
    CompoundBioAssays,
    BioAssays,
    Toxicity,
    DIRIL_Toxicity,
    DICT_Rank_Toxicity,
    ChemblMechanism,
    Substances,
    SubstanceSynonyms,
    SubstanceToxicity,
    CellLines,
    CellLineSynonyms,
    CellLineDisease,
    OncoTree,
    AntibodyDrugConjugates,
)

load_dotenv(override=True)

password_cleaned = quote_plus(os.getenv("DATABASE_PASS"))
engine = create_engine(
    f"mysql+pymysql://{os.getenv('DATABASE_USER')}:{password_cleaned}"
    f"@{os.getenv('DATABASE_IP')}:{os.getenv('PORT')}/{os.getenv('SELECTED_DB')}",
    echo=True,
)

Base.metadata.create_all(engine)


# Helper to align a dataframe to a model's columns
def align_to_model(df: pd.DataFrame, model) -> pd.DataFrame:
    cols = [c.name for c in model.__table__.columns]
    return df[[c for c in cols if c in df.columns]]


compounds_df = pd.read_csv(
    "seeding/seeding_data/jun_29_2026/union_out.csv"
)

# Build input_id -> CID lookup from union_out.csv's mapped_name column.
# Used to resolve dictrank_dataset_508.csv input_id values to pubchem_cid.
if {"cid", "mapped_name"}.issubset(compounds_df.columns):
    input_id_to_cid = compounds_df[["cid", "mapped_name"]].dropna().copy()
    input_id_to_cid["input_id"] = input_id_to_cid["mapped_name"].str.split(";")
    input_id_to_cid = input_id_to_cid.explode("input_id")
    input_id_to_cid["input_id"] = input_id_to_cid["input_id"].str.strip()
    input_id_to_cid = input_id_to_cid[input_id_to_cid["input_id"] != ""]
    input_id_to_cid = input_id_to_cid.rename(columns={"cid": "pubchem_cid"})[["input_id", "pubchem_cid"]]
    
    input_id_to_cid["pubchem_cid"] = pd.to_numeric(input_id_to_cid["pubchem_cid"], errors="coerce")
    input_id_to_cid = input_id_to_cid.dropna(subset=["pubchem_cid"])
    input_id_to_cid["pubchem_cid"] = input_id_to_cid["pubchem_cid"].astype(int)
    input_id_to_cid = input_id_to_cid.drop_duplicates(subset=["input_id"])
else:
    print("[WARNING] union_out.csv missing mapped_name or cid; dictrank CID resolution will be empty.")
    input_id_to_cid = pd.DataFrame(columns=["input_id", "pubchem_cid"])


# Load ADCs
adcs_df = pd.read_json("substance_to_cid/adcdb_results.json")
adcs_df = align_to_model(adcs_df, AntibodyDrugConjugates)

# Remove any duplicate entries based on cids
if "cid" in compounds_df.columns:
    before = len(compounds_df)
    compounds_df = compounds_df.drop_duplicates(subset=["cid"], keep="first")
    after = len(compounds_df)
    print(
        f"[Cleanup] Removed {before - after} duplicate cid rows; kept {after} unique entries."
    )

compounds_df["fda_approval"] = compounds_df["chembl_max_phase"] == 4

synonyms_df = pd.read_csv(
    "seeding/seeding_data/jun_29_2026/union_synonyms.csv"
)

# Remove any duplicate synonym entries based on cid/synonym combos
if {"synonym", "pubchem_cid"}.issubset(synonyms_df.columns):
    before = len(synonyms_df)
    synonyms_df = synonyms_df.drop_duplicates(
        subset=["synonym", "pubchem_cid"], keep="first"
    )
    after = len(synonyms_df)
    print(
        f"[Cleanup] Removed {before - after} duplicate synonym/pubchem_cid rows; kept {after} unique entries."
    )

compounds_bioassays_df = pd.read_csv(
    "seeding/seeding_data/jun_29_2026/union_bioassays.csv"
)

bioassays_df = pd.read_csv(
    "seeding/seeding_data/jun_29_2026/union_pubchem_assay_fields.csv"
)

toxicity_df = pd.read_csv(
    "seeding/seeding_data/jun_29_2026/toxicity_output_new.csv"
)

# ----- DIRIL (diril_dataset_508) -----
diril_df = (
    pd.read_csv("seeding/seeding_data/jun_29_2026/diril_dataset_508.csv")
    .rename(columns={"cid": "pubchem_cid", "My Findings  (Toxicity)": "toxicity"})
    .dropna(subset=["pubchem_cid"])
    .drop_duplicates(subset=["pubchem_cid"])
)

# ----- DICT_Rank (dictrank_dataset_508) -----
dict_rank_df = pd.read_csv("seeding/seeding_data/jun_29_2026/dictrank_dataset_508.csv")
dict_rank_df = dict_rank_df.rename(columns=lambda x: x.strip())

label_section_map = {
    "bw": "boxed warning",
    "wp": "warnings and precautions",
    "ar": "adverse reactions",
    "sp": "special conditions",
    "no": "no cardiotoxicity info",
    "overdosage": "cardiotoxicity in overdosage context",
    "withdraw": "withdrawn from market",
    "clinical pharmacology": "clinical pharmacology"
}

if "label_section" in dict_rank_df.columns:
    dict_rank_df["label_section"] = dict_rank_df["label_section"].apply(
        lambda x: label_section_map.get(str(x).strip().lower(), x) if pd.notna(x) else x
    )

dict_rank_df = (
    dict_rank_df
    .merge(input_id_to_cid, on="input_id", how="inner")
    .dropna(subset=["pubchem_cid"])
    .drop_duplicates(subset=["pubchem_cid"])
)

chembl_mech_df = pd.read_csv(
    "seeding/seeding_data/jun_29_2026/chembl_mechanism.csv"
)

substances_df = pd.read_csv(
    "data_extraction/drugs/pubchem/substance/output/combined/substance_out.csv"
)

substance_synonyms_df = pd.read_csv(
    "data_extraction/drugs/pubchem/substance/output/combined/substance_synonyms.csv"
)

substance_toxicity_df = pd.read_csv(
    "data_extraction/drugs/pubchem/substance/output/combined/substance_toxicity.csv"
)

substance_mechanism_df = pd.read_csv(
    "data_extraction/drugs/pubchem/substance/output/combined/chembl_mechanism.csv"
)

chembl_mech_df = pd.concat([chembl_mech_df, substance_mechanism_df], ignore_index=True)
chembl_mech_df = chembl_mech_df.drop_duplicates()
valid_chembls = set(compounds_df["molecule_chembl_id"].dropna().astype(str)).union(
    set(substances_df["molecule_chembl_id"].dropna().astype(str))
)

before = len(chembl_mech_df)
chembl_mech_df = chembl_mech_df[chembl_mech_df["source"] == "drug_mechanism"]  # Filter verified mechanisms only
chembl_mech_df = chembl_mech_df[
    chembl_mech_df["molecule_chembl_id"].astype(str).isin(valid_chembls)
]
print(
    f"[Mechanism] Dropped {before - len(chembl_mech_df)} orphan rows with no matching CHEMBL ID in pubchem_compounds"
)


cell_lines_df = pd.read_csv(
    "data_extraction/cell_lines/cellosaurus/output_data/cell_lines_table_cleaned.csv"
)

cell_lines_synonyms_df = pd.read_csv(
    "data_extraction/cell_lines/cellosaurus/output_data/cell_line_synonyms.csv"
)

cell_lines_disease_df = pd.read_csv(
    "data_extraction/cell_lines/cellosaurus/output_data/cell_line_diseases.csv"
)

oncotree_df = pd.read_csv("data_extraction/oncotree/output_data/oncotree.csv")

# Align columns to ORM models
compounds_df = align_to_model(compounds_df, Compounds)
synonyms_df = align_to_model(synonyms_df, CompoundSynonyms)
compounds_bioassays_df = align_to_model(compounds_bioassays_df, CompoundBioAssays)
bioassays_df = align_to_model(bioassays_df, BioAssays)
chembl_mech_df = align_to_model(chembl_mech_df, ChemblMechanism)
substances_df = align_to_model(substances_df, Substances)
substance_synonyms_df = align_to_model(substance_synonyms_df, SubstanceSynonyms)
substance_toxicity_df = align_to_model(substance_toxicity_df, SubstanceToxicity)
cell_lines_df = align_to_model(cell_lines_df, CellLines)
cell_lines_synonyms_df = align_to_model(cell_lines_synonyms_df, CellLineSynonyms)
cell_lines_disease_df = align_to_model(cell_lines_disease_df, CellLineDisease)
oncotree_df = align_to_model(oncotree_df, OncoTree)
diril_df = align_to_model(diril_df, DIRIL_Toxicity)
dict_rank_df = align_to_model(dict_rank_df, DICT_Rank_Toxicity)

# 1) Deduplicate aids that don't have a corresponding entry
bioassays_df = bioassays_df.drop_duplicates(subset=["aid"], keep="first")

# --- NEW (tiny): normalize aid type early to avoid float/NaN issues ---
if "aid" in bioassays_df.columns:
    bioassays_df["aid"] = pd.to_numeric(bioassays_df["aid"], errors="coerce")
    bioassays_df = bioassays_df.dropna(subset=["aid"])
    bioassays_df["aid"] = bioassays_df["aid"].astype(int)

if {"bioassay_aid", "pubchem_cid"}.issubset(compounds_bioassays_df.columns):
    compounds_bioassays_df = compounds_bioassays_df.drop_duplicates(
        subset=["bioassay_aid", "pubchem_cid"], keep="first"
    )

# 2) Filter mapping rows to only those AIDs that exist in bioassays_df
aids_in_bioassays = set(bioassays_df["aid"].dropna().astype(int))

# Apply the AID filter
compounds_bioassays_df["bioassay_aid"] = pd.to_numeric(
    compounds_bioassays_df["bioassay_aid"], errors="coerce"
)
compounds_bioassays_df = compounds_bioassays_df.dropna(subset=["bioassay_aid"])
compounds_bioassays_df["bioassay_aid"] = compounds_bioassays_df["bioassay_aid"].astype(
    int
)

aid_missing_mask = ~compounds_bioassays_df["bioassay_aid"].isin(aids_in_bioassays)
aid_removed_df = compounds_bioassays_df[aid_missing_mask].copy()
aid_kept_df = compounds_bioassays_df[~aid_missing_mask].copy()

if len(aid_removed_df) > 0:
    aid_removed_df.to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/removed_union_compound_bioassays_rows_missing_aid_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )
    missing_aids = sorted(aid_removed_df["bioassay_aid"].unique().tolist())
    pd.DataFrame({"aid": missing_aids}).to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/missing_aids_from_union_compound_bioassays_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )

compounds_bioassays_df = aid_kept_df

# 3) Filter compound_bioassays rows to only CIDs that exist in compounds_df (FK safety)
cids_in_compounds = set(compounds_df["cid"].dropna().astype(int))

compounds_bioassays_df["pubchem_cid"] = pd.to_numeric(
    compounds_bioassays_df["pubchem_cid"], errors="coerce"
)
compounds_bioassays_df = compounds_bioassays_df.dropna(subset=["pubchem_cid"])
compounds_bioassays_df["pubchem_cid"] = compounds_bioassays_df["pubchem_cid"].astype(
    int
)

cb_missing_cid_mask = ~compounds_bioassays_df["pubchem_cid"].isin(cids_in_compounds)
cb_removed_df = compounds_bioassays_df[cb_missing_cid_mask].copy()
cb_kept_df = compounds_bioassays_df[~cb_missing_cid_mask].copy()

if len(cb_removed_df) > 0:
    cb_removed_df.to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/removed_union_compound_bioassays_rows_missing_cid_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )
    missing_cb_cids = sorted(cb_removed_df["pubchem_cid"].unique().tolist())
    pd.DataFrame({"cid": missing_cb_cids}).to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/missing_cids_from_union_compound_bioassays_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )

compounds_bioassays_df = cb_kept_df

# 4) Filter synonyms rows to only CIDs that exist in compounds_df (FK safety)
synonyms_df["pubchem_cid"] = pd.to_numeric(synonyms_df["pubchem_cid"], errors="coerce")
synonyms_df = synonyms_df.dropna(subset=["pubchem_cid"])
synonyms_df["pubchem_cid"] = synonyms_df["pubchem_cid"].astype(int)

syn_missing_mask = ~synonyms_df["pubchem_cid"].isin(cids_in_compounds)
syn_removed_df = synonyms_df[syn_missing_mask].copy()
syn_kept_df = synonyms_df[~syn_missing_mask].copy()

if len(syn_removed_df) > 0:
    syn_removed_df.to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/removed_synonyms_rows_missing_cid_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )
    missing_syn_cids = sorted(syn_removed_df["pubchem_cid"].unique().tolist())
    pd.DataFrame({"cid": missing_syn_cids}).to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/missing_cids_from_synonyms_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )

synonyms_df = syn_kept_df

# 1) Deduplicate toxicity rows
toxicity_df = toxicity_df.drop_duplicates(subset=["pubchem_cid", "tox_dataset"], keep="first")
substance_toxicity_df = substance_toxicity_df.drop_duplicates(subset=["sid", "tox_dataset"], keep="first")

# 2) Filter toxicity rows to only CIDs that exist in compounds_df
cids_in_compounds = set(compounds_df["cid"].dropna().astype(int))

toxicity_df["pubchem_cid"] = pd.to_numeric(toxicity_df["pubchem_cid"], errors="coerce")
toxicity_df = toxicity_df.dropna(subset=["pubchem_cid", "tox_dataset"])
toxicity_df["pubchem_cid"] = toxicity_df["pubchem_cid"].astype(int)

if "dili_severity_grade" in toxicity_df.columns:
    # Force non-numeric strings like "Not Applicable" into NaN (which become NULL)
    toxicity_df["dili_severity_grade"] = pd.to_numeric(toxicity_df["dili_severity_grade"], errors="coerce")

substance_toxicity_df = substance_toxicity_df.dropna(subset=["tox_dataset"])

tox_missing_mask = ~toxicity_df["pubchem_cid"].isin(cids_in_compounds)
tox_removed_df = toxicity_df[tox_missing_mask].copy()
tox_kept_df = toxicity_df[~tox_missing_mask].copy()

if len(tox_removed_df) > 0:
    tox_removed_df.to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/removed_toxicity_rows_missing_cid_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )
    missing_cids = sorted(tox_removed_df["pubchem_cid"].unique().tolist())
    pd.DataFrame({"cid": missing_cids}).to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/missing_cids_from_toxicity_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )

toxicity_df = tox_kept_df

# FK safety: filter DILIst rows to CIDs that exist in compounds_df
diril_missing_mask = ~diril_df["pubchem_cid"].isin(cids_in_compounds)
if diril_missing_mask.any():
    diril_df[diril_missing_mask].to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/removed_diril_rows_missing_cid_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )
    print(f"[DILIst] Dropped {diril_missing_mask.sum()} rows with no matching CID in pubchem_compounds")
diril_df = diril_df[~diril_missing_mask]
print(f"[DILIst] Seeding {len(diril_df)} rows")

# FK safety: filter DILIrank rows to CIDs that exist in compounds_df
dict_rank_missing_mask = ~dict_rank_df["pubchem_cid"].isin(cids_in_compounds)
if dict_rank_missing_mask.any():
    dict_rank_df[dict_rank_missing_mask].to_csv(
        os.path.join(
            os.getcwd(),
            f"seeding/error_logs/removed_dict_rank_rows_missing_cid_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        ),
        index=False,
    )
    print(f"[DILIrank] Dropped {dict_rank_missing_mask.sum()} rows with no matching CID in pubchem_compounds")
dict_rank_df = dict_rank_df[~dict_rank_missing_mask]
print(f"[DILIrank] Seeding {len(dict_rank_df)} rows")


# Insert into tables named by the ORM models
compounds_df.to_sql(
    name=Compounds.__tablename__, con=engine, if_exists="append", index=False
)
synonyms_df.to_sql(
    name=CompoundSynonyms.__tablename__, con=engine, if_exists="append", index=False
)
bioassays_df.to_sql(
    name=BioAssays.__tablename__, con=engine, if_exists="append", index=False
)
compounds_bioassays_df.to_sql(
    name=CompoundBioAssays.__tablename__, con=engine, if_exists="append", index=False
)
toxicity_df.to_sql(
    name=Toxicity.__tablename__, con=engine, if_exists="append", index=False
)
chembl_mech_df.to_sql(
    name=ChemblMechanism.__tablename__, con=engine, if_exists="append", index=False
)
cell_lines_df.to_sql(
    name=CellLines.__tablename__, con=engine, if_exists="append", index=False
)
cell_lines_synonyms_df.to_sql(
    name=CellLineSynonyms.__tablename__, con=engine, if_exists="append", index=False
)
cell_lines_disease_df.to_sql(
    name=CellLineDisease.__tablename__, con=engine, if_exists="append", index=False
)
oncotree_df.to_sql(
    name=OncoTree.__tablename__, con=engine, if_exists="append", index=False
)
adcs_df.to_sql(
    name=AntibodyDrugConjugates.__tablename__, con=engine, if_exists="append", index=False
)
substances_df.to_sql(
    name=Substances.__tablename__, con=engine, if_exists="append", index=False
)
substance_synonyms_df.to_sql(
    name=SubstanceSynonyms.__tablename__, con=engine, if_exists="append", index=False
)
substance_toxicity_df.to_sql(
    name=SubstanceToxicity.__tablename__, con=engine, if_exists="append", index=False
)
diril_df.to_sql(
    name=DIRIL_Toxicity.__tablename__, con=engine, if_exists="append", index=False
)
dict_rank_df.to_sql(
    name=DICT_Rank_Toxicity.__tablename__, con=engine, if_exists="append", index=False
)
