from sqlalchemy import create_engine
import os
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
    ChemblMechanism,
    CellLines,
    CellLineSynonyms,
    CellLineDisease,
    OncoTree,
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
    "data_extraction/drugs/pubchem/output_data/union/complete/union_out.csv"
)

# Remove any duplicate entries based on cids
if "cid" in compounds_df.columns:
    before = len(compounds_df)
    compounds_df = compounds_df.drop_duplicates(subset=["cid"], keep="first")
    after = len(compounds_df)
    print(
        f"[Cleanup] Removed {before - after} duplicate cid rows; kept {after} unique entries."
    )

synonyms_df = pd.read_csv(
    "data_extraction/drugs/pubchem/output_data/union/complete/union_synonyms.csv"
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
    "data_extraction/drugs/pubchem/output_data/union/complete/union_bioassays.csv"
)

bioassays_df = pd.read_csv(
    "data_extraction/drugs/pubchem/output_data/union/complete/union_pubchem_assay_fields.csv"
)

toxicity_df = pd.read_csv(
    "data_extraction/drugs/pubchem/output_data/union/complete/toxicity_output.csv"
)

chembl_mech_df = pd.read_csv(
    "data_extraction/drugs/pubchem/output_data/union/complete/chembl_mechanism.csv"
)

valid_chembls = set(compounds_df["molecule_chembl_id"].dropna().astype(str))
before = len(chembl_mech_df)
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
cell_lines_df = align_to_model(cell_lines_df, CellLines)
cell_lines_synonyms_df = align_to_model(cell_lines_synonyms_df, CellLineSynonyms)
cell_lines_disease_df = align_to_model(cell_lines_disease_df, CellLineDisease)
oncotree_df = align_to_model(oncotree_df, OncoTree)

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
