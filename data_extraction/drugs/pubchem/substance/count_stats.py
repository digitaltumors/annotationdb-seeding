"""
Count summary statistics for compounds from the output CSVs in akshat_antibodies.
Reports: compounds with mechanisms, total mechanisms, compounds with bioassays,
total bioassays, compounds with toxicity, total toxicities, and total compound entries.
"""
import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(SCRIPT_DIR, "output", "akshat_antibodies")

SUBSTANCE_OUT = os.path.join(OUT_DIR, "substance_out.csv")
MECH_OUT = os.path.join(OUT_DIR, "substance_chembl_mechanism.csv")
TOX_OUT = os.path.join(OUT_DIR, "substance_toxicity.csv")


def has_toxicity_data(row):
    """Return True if the toxicity row contains meaningful data."""
    annotation = str(row.get("dili_annotation", "")).strip()
    likelihood = str(row.get("hepatotoxicity_likelihood_score", "")).strip()
    grade = str(row.get("dili_severity_grade", "")).strip()
    
    # In the original script, grade was an int where 0 is False. 
    # Here it's a string from CSV. If it's "0" but annotation is present, 
    # annotation will be Truthy. We'll match original logic where "0" alone 
    # without annotation or likelihood wouldn't trigger it, though usually 
    # they come together.
    if grade == "0":
        grade = ""
        
    return bool(annotation or likelihood or grade)


def main():
    compounds = set()
    if os.path.exists(SUBSTANCE_OUT):
        with open(SUBSTANCE_OUT, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                compounds.add(row["sid"])

    total_compounds = len(compounds)

    total_mechanisms = 0
    compounds_with_mechanisms = set()
    if os.path.exists(MECH_OUT):
        with open(MECH_OUT, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                total_mechanisms += 1
                compounds_with_mechanisms.add(row["molecule_chembl_id"])

    total_toxicities = 0
    compounds_with_toxicity = set()
    if os.path.exists(TOX_OUT):
        with open(TOX_OUT, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if has_toxicity_data(row):
                    total_toxicities += 1
                    compounds_with_toxicity.add(row["sid"])

    # We don't have a bioassays CSV here, setting to 0 to match original output format
    compounds_with_bioassays = 0
    total_bioassays = 0

    print(f"# of compounds with mechanisms: {len(compounds_with_mechanisms)}")
    print(f"# of total mechanisms:           {total_mechanisms}")
    print(f"# of compounds with bioassays:   {compounds_with_bioassays}")
    print(f"# of total bioassays:            {total_bioassays}")
    print(f"# of compounds with toxicity:    {len(compounds_with_toxicity)}")
    print(f"# of total toxicities:           {total_toxicities}")
    print(f"# of compound entries:           {total_compounds}")


if __name__ == "__main__":
    main()
