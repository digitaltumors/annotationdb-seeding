from sqlalchemy import String, Float, Integer, ForeignKey, Text, Boolean, DateTime, func
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship
from datetime import datetime


class Base(DeclarativeBase):
    pass


class Compounds(Base):
    __tablename__ = "pubchem_compounds"

    cid: Mapped[int] = mapped_column(primary_key=True)
    title: Mapped[str] = mapped_column(Text())
    mapped_name: Mapped[str] = mapped_column(Text())
    molecule_chembl_id: Mapped[str] = mapped_column(
        String(200), unique=True, nullable=True
    )
    molecular_formula: Mapped[str] = mapped_column(String(300))
    molecular_weight: Mapped[str] = mapped_column(String(50))
    smiles: Mapped[str] = mapped_column(String(2000))
    # canonical_smiles: Mapped[str] = mapped_column(String(2000))
    # isomeric_smiles: Mapped[str] = mapped_column(String(2000))
    connectivity_smiles: Mapped[str] = mapped_column(String(2000))
    inchi: Mapped[str] = mapped_column(Text())
    inchikey: Mapped[str] = mapped_column(String(28))
    iupac_name: Mapped[str] = mapped_column(Text())
    xlogp: Mapped[float] = mapped_column(Float)
    exact_mass: Mapped[str] = mapped_column(String(150))
    monoisotopic_mass: Mapped[str] = mapped_column(String(150))
    tpsa: Mapped[float] = mapped_column(Float)
    complexity: Mapped[int] = mapped_column(Integer)
    charge: Mapped[int] = mapped_column(Integer)
    h_bond_donor_count: Mapped[int] = mapped_column(Integer)
    h_bond_acceptor_count: Mapped[int] = mapped_column(Integer)
    rotatable_bond_count: Mapped[int] = mapped_column(Integer)
    heavy_atom_count: Mapped[int] = mapped_column(Integer)
    isotope_atom_count: Mapped[int] = mapped_column(Integer)
    atom_stereo_count: Mapped[int] = mapped_column(Integer)
    defined_atom_stereo_count: Mapped[int] = mapped_column(Integer)
    undefined_atom_stereo_count: Mapped[int] = mapped_column(Integer)
    bond_stereo_count: Mapped[int] = mapped_column(Integer)
    defined_bond_stereo_count: Mapped[int] = mapped_column(Integer)
    undefined_bond_stereo_count: Mapped[int] = mapped_column(Integer)
    covalent_unit_count: Mapped[int] = mapped_column(Integer)
    volume_3d: Mapped[float] = mapped_column(Float)
    x_steric_quadrupole_3d: Mapped[float] = mapped_column(Float)
    y_steric_quadrupole_3d: Mapped[float] = mapped_column(Float)
    z_steric_quadrupole_3d: Mapped[float] = mapped_column(Float)
    feature_count_3d: Mapped[int] = mapped_column(Integer)
    feature_acceptor_count_3d: Mapped[int] = mapped_column(Integer)
    feature_donor_count_3d: Mapped[int] = mapped_column(Integer)
    feature_anion_count_3d: Mapped[int] = mapped_column(Integer)
    feature_cation_count_3d: Mapped[int] = mapped_column(Integer)
    feature_ring_count_3d: Mapped[int] = mapped_column(Integer)
    feature_hydrophobe_count_3d: Mapped[int] = mapped_column(Integer)
    conformer_model_rmsd_3d: Mapped[float] = mapped_column(Float)
    effective_rotor_count_3d: Mapped[int] = mapped_column(Integer)
    conformer_count_3d: Mapped[int] = mapped_column(Integer)
    fingerprint_2d: Mapped[str] = mapped_column(Text())
    patent_count: Mapped[int] = mapped_column(Integer)
    patent_family_count: Mapped[int] = mapped_column(Integer)
    literature_count: Mapped[int] = mapped_column(Integer)
    annotation_types: Mapped[str] = mapped_column(Text())
    annotation_type_count: Mapped[int] = mapped_column(Integer)
    fda_approval: Mapped[bool] = mapped_column(Boolean)
    date_added: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )

    ### ORM layer (fields below don't show up in table but are used in queries later on for convenience)
    synonyms: Mapped[list["CompoundSynonyms"]] = relationship(
        back_populates="compound", cascade="all, delete-orphan"
    )

    mechanisms: Mapped[list["ChemblMechanism"]] = relationship(
        "ChemblMechanism",
        primaryjoin="Compounds.molecule_chembl_id == foreign(ChemblMechanism.molecule_chembl_id)",
        back_populates="compound",
        lazy="selectin",
    )

    compound_bioassays: Mapped[list["CompoundBioAssays"]] = relationship(
        "CompoundBioAssays",
        back_populates="compound",
        cascade="all, delete-orphan",
    )

    bioassays: Mapped[list["BioAssays"]] = relationship(
        "BioAssays",
        secondary="compound_bioassays",
        back_populates="compounds",
        lazy="selectin",
    )

    toxicity: Mapped["Toxicity"] = relationship(
        "Toxicity",
        back_populates="compound",
        uselist=False,
        cascade="all, delete-orphan",
    )


class CompoundSynonyms(Base):
    __tablename__ = "compound_synonyms"

    synonym: Mapped[str] = mapped_column(String(700), primary_key=True)
    pubchem_cid: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("pubchem_compounds.cid", ondelete="CASCADE"),
        primary_key=True,
    )
    source: Mapped[str] = mapped_column(String(50))
    version: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )

    ### ORM layer (fields below don't show up in table but are used in queries later on for convenience)
    compound: Mapped["Compounds"] = relationship(back_populates="synonyms")


class CompoundBioAssays(Base):
    __tablename__ = "compound_bioassays"

    bioassay_aid: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("bioassays.aid"),
        primary_key=True,
    )
    pubchem_cid: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("pubchem_compounds.cid", ondelete="CASCADE"),
        primary_key=True,
    )


class BioAssays(Base):
    __tablename__ = "bioassays"

    aid: Mapped[int] = mapped_column(
        Integer,
        primary_key=True,
    )
    version: Mapped[int] = mapped_column(Integer)
    assay_name: Mapped[str] = mapped_column(String(512))
    source_name: Mapped[str] = mapped_column(String(255))
    source_id: Mapped[str] = mapped_column(String(255))
    description_combined: Mapped[str] = mapped_column(Text)
    protocol_combined: Mapped[str] = mapped_column(Text)
    comment_combined: Mapped[str] = mapped_column(Text)
    activity_outcome_method: Mapped[int] = mapped_column(Integer)
    target_name: Mapped[str] = mapped_column(Text)
    target_protein_accession: Mapped[str] = mapped_column(Text)


class Toxicity(Base):
    __tablename__ = "toxicity"

    pubchem_cid: Mapped[int] = mapped_column(
        Integer,
        ForeignKey("pubchem_compounds.cid", ondelete="CASCADE"),
        primary_key=True,
    )
    dili_severity_grade: Mapped[int] = mapped_column(Integer)
    dili_annotation: Mapped[str] = mapped_column(Text)
    hepatotoxicity_likelihood_score: Mapped[str] = mapped_column(Text)


# https://www.ebi.ac.uk/chembl/api/data/drug/schema?format=json
# class ChemblDrugData(Base):
#     __tablename__ = "chembl_drug_data"

#     molecule_chembl_id: Mapped[str] = mapped_column(
#         String(200), ForeignKey("compounds.molecule_chembl_id", ondelete="CASCADE")
#     )
#     chirality: Mapped[int] = mapped_column(Integer)
#     alogp: Mapped[float] = mapped_column(Float)
#     aromatic_rings: Mapped[int] = mapped_column(Integer)
#     full_molformula: Mapped[str] = mapped_column(String(30))
#     full_mwt: Mapped[float] = mapped_column(Float)
#     hba: Mapped[int] = mapped_column(Integer)
#     hbd: Mapped[int] = mapped_column(Integer)
#     heavy_atoms: Mapped[int] = mapped_column(Integer)
#     mw_freebase: Mapped[float] = mapped_column(Float)
#     np_likeness_score: Mapped[float] = mapped_column(Float)
#     num_ro5_violations: Mapped[int] = mapped_column(Integer)
#     psa: Mapped[float] = mapped_column(Float)
#     qed_weighted: Mapped[float] = mapped_column(Float)
#     ro3_pass: Mapped[str] = mapped_column(String(10))
#     rtb: Mapped[int] = mapped_column(Integer)


# https://www.ebi.ac.uk/chembl/api/data/mechanism/schema?format=json
class ChemblMechanism(Base):
    __tablename__ = "chembl_mechanism"

    molecule_chembl_id: Mapped[str] = mapped_column(
        String(200),
        ForeignKey("pubchem_compounds.molecule_chembl_id", ondelete="CASCADE"),
        primary_key=True,
    )
    parent_molecule_chembl_id: Mapped[str] = mapped_column(String(200))
    action_type: Mapped[str] = mapped_column(String(30))
    binding_site_comment: Mapped[str] = mapped_column(String(30))
    mechanism_of_action: Mapped[str] = mapped_column(Text())
    mechanism_comment: Mapped[str] = mapped_column(Text())
    direct_interaction: Mapped[int] = mapped_column(Integer)
    disease_efficacy: Mapped[int] = mapped_column(Integer)
    max_phase: Mapped[int] = mapped_column(Integer)
    mec_id: Mapped[int] = mapped_column(
        Integer,
        primary_key=True,
    )
    # mechanism_refs: Mapped[] = mapped_column() # best way to represent?
    molecular_mechanism: Mapped[int] = mapped_column(Integer)
    record_id: Mapped[int] = mapped_column(Integer)
    selectivity_comment: Mapped[str] = mapped_column(Text())
    site_id: Mapped[str] = mapped_column(Text())
    target_chembl_id: Mapped[str] = mapped_column(String(30))
    # Fields below are extracted from the variant sequence object (if it exists [like does not])
    variant_sequence_accession: Mapped[str] = mapped_column(String(50))
    variant_sequence_isoform: Mapped[int] = mapped_column(Integer)
    variant_sequence_mutation: Mapped[str] = mapped_column(String(50))
    variant_sequence_organism: Mapped[str] = mapped_column(String(50))
    variant_sequence_sequence: Mapped[str] = mapped_column(Text())
    variant_sequence_tax_id: Mapped[int] = mapped_column(Integer)
    variant_sequence_version: Mapped[int] = mapped_column(Integer)

    compound: Mapped["Compounds"] = relationship(
        "Compounds",
        primaryjoin="foreign(ChemblMechanism.molecule_chembl_id) == Compounds.molecule_chembl_id",
        back_populates="mechanisms",
    )


# https://www.ebi.ac.uk/chembl/api/data/molecule/schema?format=json
# class ChemblMolecule(Base):
#     __tablename__ = "chembl_molecule"

#     molecule_chembl_id: Mapped[str] = mapped_column(
#         String(200), ForeignKey("compounds.molecule_chembl_id", ondelete="CASCADE")
#     )
#     dosed_ingredient: Mapped[bool] = mapped_column(Boolean)
#     first_approval: Mapped[int] = mapped_column(Integer)
#     first_in_class: Mapped[int] = mapped_column(Integer)
#     molecule_type: Mapped[str] = mapped_column(String(50))
#     natural_product: Mapped[int] = mapped_column(Integer)
#     oral: Mapped[bool] = mapped_column(Boolean)
#     orphan: Mapped[int] = mapped_column(Integer)
#     parenteral: Mapped[bool] = mapped_column(Boolean)
#     polymerflag: Mapped[int] = mapped_column(Integer)
#     prodrug: Mapped[int] = mapped_column(Integer)
#     therapeutic_flag: Mapped[bool] = mapped_column(Boolean)
#     topical: Mapped[bool] = mapped_column(Boolean)


# https://www.ebi.ac.uk/chembl/api/data/activity/schema?format=json
# class ChemblActivity(Base):
#     __tablename__ = "chembl_activity"

#     molecule_chembl_id: Mapped[str] = mapped_column(
#         String(200), ForeignKey("compounds.molecule_chembl_id", ondelete="CASCADE")
#     )
#     action_type: Mapped[str] = mapped_column(String(100))
#     activity_comment: Mapped[str] = mapped_column(Text())
#     activity_id: Mapped[int] = mapped_column(Integer)
#     # activity_properties Mapped[] = mapped_column() # Not sure what to do with this
#     assay_chembl_id: Mapped[str] = mapped_column(String(50))
#     assay_description: Mapped[str] = mapped_column(Text())
#     assay_type: Mapped[str] = mapped_column(String(10))
#     assay_variant_accession: Mapped[str] = mapped_column(String(50))
#     assay_variant_mutation: Mapped[str] = mapped_column(String(50))
#     bao_endpoint: Mapped[str] = mapped_column(String(50))
#     bao_format: Mapped[str] = mapped_column(String(50))
#     document_chembl_id: Mapped[str] = mapped_column(String(50))
#     document_journal: Mapped[str] = mapped_column(String(50))
#     document_year: Mapped[int] = mapped_column(Integer)
#     # ligand_efficiency: Mapped[] = mapped_column() # Not sure we need to include this yet
#     pchembl_value: Mapped[float] = mapped_column(Float)
#     qutd_units: Mapped[str] = mapped_column(String(50))
#     standard_flag: Mapped[int] = mapped_column(Integer)
#     standard_relation: Mapped[bool] = mapped_column(Boolean)
#     standard_text_value: Mapped[str] = mapped_column(Text())
#     standard_type: Mapped[str] = mapped_column(String(20))
#     standard_units: Mapped[str] = mapped_column(String(10))
#     standard_upper_value: Mapped[float] = mapped_column(Float)
#     standard_value: Mapped[float] = mapped_column(Float)
#     target_chembl_id: Mapped[str] = mapped_column(String(200))
#     target_organism: Mapped[str] = mapped_column(String(30))
#     target_pref_name: Mapped[str] = mapped_column(String(50))
#     target_tax_id: Mapped[int] = mapped_column(Integer)
#     toid: Mapped[int] = mapped_column(Integer)
#     type: Mapped[str] = mapped_column(String(50))
#     units: Mapped[str] = mapped_column(String(20))
#     uo_units: Mapped[str] = mapped_column(String(50))
#     upper_value: Mapped[float] = mapped_column(Float)
#     value: Mapped[float] = mapped_column(Float)


class CellLines(Base):
    __tablename__ = "cellosaurus_cell_lines"

    accession: Mapped[str] = mapped_column(String(32), primary_key=True)
    cell_line_name: Mapped[str] = mapped_column(String(200))
    category: Mapped[str] = mapped_column(String(100))
    date: Mapped[str] = mapped_column(Text())
    age_at_sampling: Mapped[str] = mapped_column(String(50))
    sex_of_cell: Mapped[str] = mapped_column(String(30))
    hierarchy: Mapped[str] = mapped_column(Text())

    # Comment fields
    cell_type: Mapped[str] = mapped_column(Text())
    derived_from_site: Mapped[str] = mapped_column(Text())
    donor_information: Mapped[str] = mapped_column(Text())
    doubling_time: Mapped[str] = mapped_column(Text())
    genome_ancestry: Mapped[str] = mapped_column(Text())
    hla_typing: Mapped[str] = mapped_column(Text())
    microsatellite_instability: Mapped[str] = mapped_column(Text())
    omics: Mapped[str] = mapped_column(Text())
    part_of: Mapped[str] = mapped_column(Text())
    population: Mapped[str] = mapped_column(Text())
    sequence_variation: Mapped[str] = mapped_column(Text())
    anecdotal: Mapped[str] = mapped_column(Text())
    biotechnology: Mapped[str] = mapped_column(Text())
    discontinued: Mapped[str] = mapped_column(Text())
    group_col: Mapped[str] = mapped_column(Text())
    misspelling: Mapped[str] = mapped_column(Text())
    registration: Mapped[str] = mapped_column(Text())
    virology: Mapped[str] = mapped_column(Text())
    caution: Mapped[str] = mapped_column(Text())
    characteristics: Mapped[str] = mapped_column(Text())
    karyotypic_information: Mapped[str] = mapped_column(Text())
    problematic_cell_line: Mapped[str] = mapped_column(Text())
    transformant: Mapped[str] = mapped_column(Text())
    miscellaneous: Mapped[str] = mapped_column(Text())
    from_col: Mapped[str] = mapped_column(Text())
    genetic_integration: Mapped[str] = mapped_column(Text())
    knockout_cell: Mapped[str] = mapped_column(Text())
    selected_for_resistance_to: Mapped[str] = mapped_column(Text())

    synonyms: Mapped[list["CellLineSynonyms"]] = relationship(
        back_populates="cell_line", cascade="all, delete-orphan"
    )

    diseases: Mapped[list["CellLineDisease"]] = relationship(
        back_populates="cell_line", cascade="all, delete-orphan"
    )


class CellLineSynonyms(Base):
    __tablename__ = "cell_line_synonyms"

    cellosaurus_accession: Mapped[str] = mapped_column(
        String(32),
        ForeignKey("cellosaurus_cell_lines.accession", ondelete="CASCADE"),
        primary_key=True,
    )
    synonym: Mapped[str] = mapped_column(String(100), primary_key=True)
    source: Mapped[str] = mapped_column(String(50), primary_key=True)
    version: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )

    ### ORM layer (fields below don't show up in table but are used in queries later on for convenience)
    cell_line: Mapped["CellLines"] = relationship(back_populates="synonyms")


class CellLineDisease(Base):
    __tablename__ = "cell_line_disease"

    cellosaurus_accession: Mapped[str] = mapped_column(
        String(32),
        ForeignKey("cellosaurus_cell_lines.accession", ondelete="CASCADE"),
        primary_key=True,
    )
    id: Mapped[str] = mapped_column(String(32), primary_key=True)
    source: Mapped[str] = mapped_column(String(4), primary_key=True)
    description: Mapped[str] = mapped_column(Text())

    ### ORM layer (fields below don't show up in table but are used in queries later on for convenience)
    cell_line: Mapped["CellLines"] = relationship(back_populates="diseases")


class OncoTree(Base):
    __tablename__ = "oncotree"

    # Tumor-type fields
    code: Mapped[str] = mapped_column(String(32), primary_key=True)
    name: Mapped[str] = mapped_column(Text())
    main_type: Mapped[str] = mapped_column(String(200))
    color: Mapped[str] = mapped_column(String(64))
    level: Mapped[int] = mapped_column(Integer)
    parent_code: Mapped[str] = mapped_column(String(32), index=True)

    external_references: Mapped[str] = mapped_column(Text(), nullable=True)
    history: Mapped[str] = mapped_column(Text(), nullable=True)
    revocations: Mapped[str] = mapped_column(Text(), nullable=True)
    precursors: Mapped[str] = mapped_column(Text(), nullable=True)
    tissue: Mapped[str] = mapped_column(Text(), nullable=True)

    version_api_identifier: Mapped[str] = mapped_column(String(100), index=True)
    version_release_date: Mapped[str] = mapped_column(String(10))
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now()
    )
