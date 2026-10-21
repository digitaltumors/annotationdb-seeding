SET SESSION net_read_timeout = 600;

-- ============================================================================
-- ANNOTATIONDB AUDIT: annotationdb (new) vs annotationdb_4 (old)
-- Strictly read-only — no data modifications.
-- ============================================================================

-- 1. ROW-COUNT COMPARISON — quick sanity check across all tables

SELECT
    'pubchem_compounds'   AS table_name,
    (SELECT COUNT(*) FROM annotationdb_4.pubchem_compounds)  AS db4_count,
    (SELECT COUNT(*) FROM annotationdb.pubchem_compounds)  AS db_count,
    (SELECT COUNT(*) FROM annotationdb.pubchem_compounds)
      - (SELECT COUNT(*) FROM annotationdb_4.pubchem_compounds) AS diff
UNION ALL SELECT
    'compound_synonyms',
    (SELECT COUNT(*) FROM annotationdb_4.compound_synonyms),
    (SELECT COUNT(*) FROM annotationdb.compound_synonyms),
    (SELECT COUNT(*) FROM annotationdb.compound_synonyms)
      - (SELECT COUNT(*) FROM annotationdb_4.compound_synonyms)
UNION ALL SELECT
    'bioassays',
    (SELECT COUNT(*) FROM annotationdb_4.bioassays),
    (SELECT COUNT(*) FROM annotationdb.bioassays),
    (SELECT COUNT(*) FROM annotationdb.bioassays)
      - (SELECT COUNT(*) FROM annotationdb_4.bioassays)
UNION ALL SELECT
    'compound_bioassays',
    (SELECT COUNT(*) FROM annotationdb_4.compound_bioassays),
    (SELECT COUNT(*) FROM annotationdb.compound_bioassays),
    (SELECT COUNT(*) FROM annotationdb.compound_bioassays)
      - (SELECT COUNT(*) FROM annotationdb_4.compound_bioassays)
UNION ALL SELECT
    'toxicity',
    (SELECT COUNT(*) FROM annotationdb_4.toxicity),
    (SELECT COUNT(*) FROM annotationdb.toxicity),
    (SELECT COUNT(*) FROM annotationdb.toxicity)
      - (SELECT COUNT(*) FROM annotationdb_4.toxicity)
UNION ALL SELECT
    'chembl_mechanism',
    (SELECT COUNT(*) FROM annotationdb_4.chembl_mechanism),
    (SELECT COUNT(*) FROM annotationdb.chembl_mechanism),
    (SELECT COUNT(*) FROM annotationdb.chembl_mechanism)
      - (SELECT COUNT(*) FROM annotationdb_4.chembl_mechanism)
UNION ALL SELECT
    'cellosaurus_cell_lines',
    (SELECT COUNT(*) FROM annotationdb_4.cellosaurus_cell_lines),
    (SELECT COUNT(*) FROM annotationdb.cellosaurus_cell_lines),
    (SELECT COUNT(*) FROM annotationdb.cellosaurus_cell_lines)
      - (SELECT COUNT(*) FROM annotationdb_4.cellosaurus_cell_lines)
UNION ALL SELECT
    'cell_line_synonyms',
    (SELECT COUNT(*) FROM annotationdb_4.cell_line_synonyms),
    (SELECT COUNT(*) FROM annotationdb.cell_line_synonyms),
    (SELECT COUNT(*) FROM annotationdb.cell_line_synonyms)
      - (SELECT COUNT(*) FROM annotationdb_4.cell_line_synonyms)
UNION ALL SELECT
    'cell_line_disease',
    (SELECT COUNT(*) FROM annotationdb_4.cell_line_disease),
    (SELECT COUNT(*) FROM annotationdb.cell_line_disease),
    (SELECT COUNT(*) FROM annotationdb.cell_line_disease)
      - (SELECT COUNT(*) FROM annotationdb_4.cell_line_disease)
UNION ALL SELECT
    'oncotree',
    (SELECT COUNT(*) FROM annotationdb_4.oncotree),
    (SELECT COUNT(*) FROM annotationdb.oncotree),
    (SELECT COUNT(*) FROM annotationdb.oncotree)
      - (SELECT COUNT(*) FROM annotationdb_4.oncotree)
-- New tables (db only — db4 count shown as 0)
UNION ALL SELECT
    'substances',
    0,
    (SELECT COUNT(*) FROM annotationdb.substances),
    (SELECT COUNT(*) FROM annotationdb.substances)
UNION ALL SELECT
    'substance_synonyms',
    0,
    (SELECT COUNT(*) FROM annotationdb.substance_synonyms),
    (SELECT COUNT(*) FROM annotationdb.substance_synonyms)
UNION ALL SELECT
    'substance_toxicity',
    0,
    (SELECT COUNT(*) FROM annotationdb.substance_toxicity),
    (SELECT COUNT(*) FROM annotationdb.substance_toxicity)
UNION ALL SELECT
    'diril_toxicity',
    0,
    (SELECT COUNT(*) FROM annotationdb.diril_toxicity),
    (SELECT COUNT(*) FROM annotationdb.diril_toxicity)
UNION ALL SELECT
    'dict_rank_toxicity',
    0,
    (SELECT COUNT(*) FROM annotationdb.dict_rank_toxicity),
    (SELECT COUNT(*) FROM annotationdb.dict_rank_toxicity)
UNION ALL SELECT
    'antibody_drug_conjugates',
    0,
    (SELECT COUNT(*) FROM annotationdb.antibody_drug_conjugates),
    (SELECT COUNT(*) FROM annotationdb.antibody_drug_conjugates);


-- 2. PRIMARY KEY MISMATCHES — rows that exist in one DB but not the other

-- 2a. pubchem_compounds  (PK: cid)
-- Rows in db4 but NOT in db
SELECT 'pubchem_compounds: in db4 only' AS check_type, a.cid
FROM annotationdb_4.pubchem_compounds a
LEFT JOIN annotationdb.pubchem_compounds b ON a.cid = b.cid
WHERE b.cid IS NULL
LIMIT 100;

-- Rows in db but NOT in db4
SELECT 'pubchem_compounds: in db only' AS check_type, b.cid
FROM annotationdb.pubchem_compounds b
LEFT JOIN annotationdb_4.pubchem_compounds a ON b.cid = a.cid
WHERE a.cid IS NULL
LIMIT 100;

-- 2b. compound_synonyms  (PK: synonym, pubchem_cid)
SELECT 'compound_synonyms: in db4 only' AS check_type, a.synonym, a.pubchem_cid
FROM annotationdb_4.compound_synonyms a
LEFT JOIN annotationdb.compound_synonyms b
       ON a.synonym COLLATE utf8mb4_general_ci = b.synonym COLLATE utf8mb4_general_ci AND a.pubchem_cid = b.pubchem_cid
WHERE b.synonym IS NULL
LIMIT 100;

SELECT 'compound_synonyms: in db only' AS check_type, b.synonym, b.pubchem_cid
FROM annotationdb.compound_synonyms b
LEFT JOIN annotationdb_4.compound_synonyms a
       ON b.synonym COLLATE utf8mb4_general_ci = a.synonym COLLATE utf8mb4_general_ci AND b.pubchem_cid = a.pubchem_cid
WHERE a.synonym IS NULL
LIMIT 100;

-- 2c. bioassays  (PK: aid)
SELECT 'bioassays: in db4 only' AS check_type, a.aid
FROM annotationdb_4.bioassays a
LEFT JOIN annotationdb.bioassays b ON a.aid = b.aid
WHERE b.aid IS NULL
LIMIT 100;

SELECT 'bioassays: in db only' AS check_type, b.aid
FROM annotationdb.bioassays b
LEFT JOIN annotationdb_4.bioassays a ON b.aid = a.aid
WHERE a.aid IS NULL
LIMIT 100;

-- 2d. compound_bioassays  (PK: bioassay_aid, pubchem_cid)
SELECT 'compound_bioassays: in db4 only' AS check_type, a.bioassay_aid, a.pubchem_cid
FROM annotationdb_4.compound_bioassays a
LEFT JOIN annotationdb.compound_bioassays b
       ON a.bioassay_aid = b.bioassay_aid AND a.pubchem_cid = b.pubchem_cid
WHERE b.bioassay_aid IS NULL
LIMIT 100;

SELECT 'compound_bioassays: in db only' AS check_type, b.bioassay_aid, b.pubchem_cid
FROM annotationdb.compound_bioassays b
LEFT JOIN annotationdb_4.compound_bioassays a
       ON b.bioassay_aid = a.bioassay_aid AND b.pubchem_cid = a.pubchem_cid
WHERE a.bioassay_aid IS NULL
LIMIT 100;

-- 2e. toxicity  (PK: pubchem_cid, tox_dataset)
SELECT 'toxicity: in db4 only' AS check_type, a.pubchem_cid, a.tox_dataset
FROM annotationdb_4.toxicity a
LEFT JOIN annotationdb.toxicity b ON a.pubchem_cid = b.pubchem_cid AND a.tox_dataset = b.tox_dataset
WHERE b.pubchem_cid IS NULL
LIMIT 100;

SELECT 'toxicity: in db only' AS check_type, b.pubchem_cid, b.tox_dataset
FROM annotationdb.toxicity b
LEFT JOIN annotationdb_4.toxicity a ON b.pubchem_cid = a.pubchem_cid AND b.tox_dataset = a.tox_dataset
WHERE a.pubchem_cid IS NULL
LIMIT 100;

-- 2f. chembl_mechanism  (PK: molecule_chembl_id, mec_id)
SELECT 'chembl_mechanism: in db4 only' AS check_type, a.molecule_chembl_id, a.mec_id
FROM annotationdb_4.chembl_mechanism a
LEFT JOIN annotationdb.chembl_mechanism b
       ON a.molecule_chembl_id COLLATE utf8mb4_general_ci = b.molecule_chembl_id COLLATE utf8mb4_general_ci AND a.mec_id = b.mec_id
WHERE b.molecule_chembl_id IS NULL
LIMIT 100;

SELECT 'chembl_mechanism: in db only' AS check_type, b.molecule_chembl_id, b.mec_id
FROM annotationdb.chembl_mechanism b
LEFT JOIN annotationdb_4.chembl_mechanism a
       ON b.molecule_chembl_id COLLATE utf8mb4_general_ci = a.molecule_chembl_id COLLATE utf8mb4_general_ci AND b.mec_id = a.mec_id
WHERE a.molecule_chembl_id IS NULL
LIMIT 100;


-- 3. FULL ROW-LEVEL DIFF (CHECKSUM) — find rows that exist in both DBs

-- 3a. pubchem_compounds — rows where PK matches but data differs
--     NOTE: atc_code is new in db; excluded from checksum since db4 won't have it
SELECT 'pubchem_compounds: data mismatch' AS check_type, a.cid
FROM annotationdb_4.pubchem_compounds a
JOIN annotationdb.pubchem_compounds b ON a.cid = b.cid
WHERE MD5(CONCAT_WS('|',
    COALESCE(a.title,''), COALESCE(a.mapped_name,''), COALESCE(a.molecule_chembl_id,''),
    COALESCE(a.molecular_formula,''), COALESCE(a.molecular_weight,''),
    COALESCE(a.smiles,''), COALESCE(a.connectivity_smiles,''),
    COALESCE(a.inchi,''), COALESCE(a.inchikey,''), COALESCE(a.iupac_name,''),
    COALESCE(a.xlogp,0), COALESCE(a.exact_mass,''), COALESCE(a.monoisotopic_mass,''),
    COALESCE(a.tpsa,0), COALESCE(a.complexity,0), COALESCE(a.charge,0),
    COALESCE(a.h_bond_donor_count,0), COALESCE(a.h_bond_acceptor_count,0),
    COALESCE(a.rotatable_bond_count,0), COALESCE(a.heavy_atom_count,0),
    COALESCE(a.isotope_atom_count,0), COALESCE(a.atom_stereo_count,0),
    COALESCE(a.defined_atom_stereo_count,0), COALESCE(a.undefined_atom_stereo_count,0),
    COALESCE(a.bond_stereo_count,0), COALESCE(a.defined_bond_stereo_count,0),
    COALESCE(a.undefined_bond_stereo_count,0), COALESCE(a.covalent_unit_count,0),
    COALESCE(a.volume_3d,0), COALESCE(a.fingerprint_2d,''),
    COALESCE(a.patent_count,0), COALESCE(a.patent_family_count,0),
    COALESCE(a.literature_count,0), COALESCE(a.annotation_types,''),
    COALESCE(a.annotation_type_count,0), COALESCE(a.chembl_max_phase,0),
    COALESCE(a.drug_like,0), COALESCE(a.fda_approval,0)
)) != MD5(CONCAT_WS('|',
    COALESCE(b.title,''), COALESCE(b.mapped_name,''), COALESCE(b.molecule_chembl_id,''),
    COALESCE(b.molecular_formula,''), COALESCE(b.molecular_weight,''),
    COALESCE(b.smiles,''), COALESCE(b.connectivity_smiles,''),
    COALESCE(b.inchi,''), COALESCE(b.inchikey,''), COALESCE(b.iupac_name,''),
    COALESCE(b.xlogp,0), COALESCE(b.exact_mass,''), COALESCE(b.monoisotopic_mass,''),
    COALESCE(b.tpsa,0), COALESCE(b.complexity,0), COALESCE(b.charge,0),
    COALESCE(b.h_bond_donor_count,0), COALESCE(b.h_bond_acceptor_count,0),
    COALESCE(b.rotatable_bond_count,0), COALESCE(b.heavy_atom_count,0),
    COALESCE(b.isotope_atom_count,0), COALESCE(b.atom_stereo_count,0),
    COALESCE(b.defined_atom_stereo_count,0), COALESCE(b.undefined_atom_stereo_count,0),
    COALESCE(b.bond_stereo_count,0), COALESCE(b.defined_bond_stereo_count,0),
    COALESCE(b.undefined_bond_stereo_count,0), COALESCE(b.covalent_unit_count,0),
    COALESCE(b.volume_3d,0), COALESCE(b.fingerprint_2d,''),
    COALESCE(b.patent_count,0), COALESCE(b.patent_family_count,0),
    COALESCE(b.literature_count,0), COALESCE(b.annotation_types,''),
    COALESCE(b.annotation_type_count,0), COALESCE(b.chembl_max_phase,0),
    COALESCE(b.drug_like,0), COALESCE(b.fda_approval,0)
))
LIMIT 100;

-- 3b. compound_synonyms — data mismatch (only non-PK column is `source`)
SELECT 'compound_synonyms: data mismatch' AS check_type, a.synonym, a.pubchem_cid,
       a.source AS db4_source, b.source AS db_source
FROM annotationdb_4.compound_synonyms a
JOIN annotationdb.compound_synonyms b
     ON a.synonym COLLATE utf8mb4_general_ci = b.synonym COLLATE utf8mb4_general_ci AND a.pubchem_cid = b.pubchem_cid
WHERE COALESCE(a.source,'') != COALESCE(b.source,'')
LIMIT 100;

-- 3c. bioassays — data mismatch
SELECT 'bioassays: data mismatch' AS check_type, a.aid
FROM annotationdb_4.bioassays a
JOIN annotationdb.bioassays b ON a.aid = b.aid
WHERE MD5(CONCAT_WS('|',
    COALESCE(a.version,0), COALESCE(a.assay_name,''), COALESCE(a.source_name,''),
    COALESCE(a.source_id,''), COALESCE(a.description_combined,''),
    COALESCE(a.protocol_combined,''), COALESCE(a.comment_combined,''),
    COALESCE(a.activity_outcome_method,0), COALESCE(a.target_name,''),
    COALESCE(a.target_protein_accession,'')
)) != MD5(CONCAT_WS('|',
    COALESCE(b.version,0), COALESCE(b.assay_name,''), COALESCE(b.source_name,''),
    COALESCE(b.source_id,''), COALESCE(b.description_combined,''),
    COALESCE(b.protocol_combined,''), COALESCE(b.comment_combined,''),
    COALESCE(b.activity_outcome_method,0), COALESCE(b.target_name,''),
    COALESCE(b.target_protein_accession,'')
))
LIMIT 100;

-- 3d. toxicity — data mismatch
SELECT 'toxicity: data mismatch' AS check_type, a.pubchem_cid,
       a.reference_number       AS db4_ref, b.reference_number       AS db_ref,
       a.dili_severity_grade    AS db4_sev, b.dili_severity_grade    AS db_sev
FROM annotationdb_4.toxicity a
JOIN annotationdb.toxicity b ON a.pubchem_cid = b.pubchem_cid AND a.tox_dataset = b.tox_dataset
WHERE MD5(CONCAT_WS('|',
    COALESCE(a.reference_number,0), COALESCE(a.dili_dataset,''),
    COALESCE(a.dili_severity_grade,0), COALESCE(a.dili_annotation,''),
    COALESCE(a.tox_source_url,''), COALESCE(a.hepatotoxicity_likelihood_score,'')
)) != MD5(CONCAT_WS('|',
    COALESCE(b.reference_number,0), COALESCE(b.dili_dataset,''),
    COALESCE(b.dili_severity_grade,0), COALESCE(b.dili_annotation,''),
    COALESCE(b.tox_source_url,''), COALESCE(b.hepatotoxicity_likelihood_score,'')
))
LIMIT 100;

-- 3e. chembl_mechanism — data mismatch
SELECT 'chembl_mechanism: data mismatch' AS check_type,
       a.molecule_chembl_id, a.mec_id
FROM annotationdb_4.chembl_mechanism a
JOIN annotationdb.chembl_mechanism b
     ON a.molecule_chembl_id COLLATE utf8mb4_general_ci = b.molecule_chembl_id COLLATE utf8mb4_general_ci AND a.mec_id = b.mec_id
WHERE MD5(CONCAT_WS('|',
    COALESCE(a.parent_molecule_chembl_id,''), COALESCE(a.action_type,''),
    COALESCE(a.mechanism_of_action,''), COALESCE(a.mechanism_comment,''),
    COALESCE(a.direct_interaction,0), COALESCE(a.disease_efficacy,0),
    COALESCE(a.max_phase,0), COALESCE(a.mechanism_refs,''),
    COALESCE(a.molecular_mechanism,0), COALESCE(a.record_id,0),
    COALESCE(a.target_chembl_id,''), COALESCE(a.target_name,''),
    COALESCE(a.prodrug,0), COALESCE(a.parent_prodrug,0), COALESCE(a.source,''),
    COALESCE(a.inferred_from_parent,0)
)) != MD5(CONCAT_WS('|',
    COALESCE(b.parent_molecule_chembl_id,''), COALESCE(b.action_type,''),
    COALESCE(b.mechanism_of_action,''), COALESCE(b.mechanism_comment,''),
    COALESCE(b.direct_interaction,0), COALESCE(b.disease_efficacy,0),
    COALESCE(b.max_phase,0), COALESCE(b.mechanism_refs,''),
    COALESCE(b.molecular_mechanism,0), COALESCE(b.record_id,0),
    COALESCE(b.target_chembl_id,''), COALESCE(b.target_name,''),
    COALESCE(b.prodrug,0), COALESCE(b.parent_prodrug,0), COALESCE(b.source,''),
    COALESCE(b.inferred_from_parent,0)
))
LIMIT 100;


-- 4. NEW COLUMN AUDIT: atc_code in pubchem_compounds (db only)
--    How many compounds got an ATC code assigned?
SELECT 'atc_code: populated' AS check_type,
       COUNT(*) AS compounds_with_atc
FROM annotationdb.pubchem_compounds
WHERE atc_code IS NOT NULL AND atc_code != '';

-- Sample of populated atc_code values
SELECT cid, title, atc_code
FROM annotationdb.pubchem_compounds
WHERE atc_code IS NOT NULL AND atc_code != ''
LIMIT 20;


-- 5. NEW TABLE AUDITS (db only — sanity checks on row counts and data quality)

-- 5a. substances — basic stats
SELECT 'substances: total rows' AS check_type,
       COUNT(*) AS row_count,
       COUNT(DISTINCT sid) AS unique_sids
FROM annotationdb.substances;

-- 5b. substance_synonyms
SELECT 'substance_synonyms: total rows' AS check_type,
       COUNT(*) AS row_count
FROM annotationdb.substance_synonyms;

-- FK orphan check: substance_synonyms referencing nonexistent SIDs
SELECT 'substance_synonyms: orphan SIDs' AS check_type, ss.sid
FROM annotationdb.substance_synonyms ss
LEFT JOIN annotationdb.substances s ON ss.sid = s.sid
WHERE s.sid IS NULL
LIMIT 20;

-- 5c. substance_toxicity
SELECT 'substance_toxicity: total rows' AS check_type,
       COUNT(*) AS row_count
FROM annotationdb.substance_toxicity;

-- FK orphan check
SELECT 'substance_toxicity: orphan SIDs' AS check_type, st.sid
FROM annotationdb.substance_toxicity st
LEFT JOIN annotationdb.substances s ON st.sid = s.sid
WHERE s.sid IS NULL
LIMIT 20;

-- 5d. diril_toxicity
SELECT 'diril_toxicity: total rows' AS check_type,
       COUNT(*) AS row_count
FROM annotationdb.diril_toxicity;

-- FK orphan check
SELECT 'diril_toxicity: orphan CIDs' AS check_type, dt.pubchem_cid
FROM annotationdb.diril_toxicity dt
LEFT JOIN annotationdb.pubchem_compounds pc ON dt.pubchem_cid = pc.cid
WHERE pc.cid IS NULL
LIMIT 20;

-- 5e. dict_rank_toxicity
SELECT 'dict_rank_toxicity: total rows' AS check_type,
       COUNT(*) AS row_count
FROM annotationdb.dict_rank_toxicity;

-- FK orphan check
SELECT 'dict_rank_toxicity: orphan CIDs' AS check_type, dr.pubchem_cid
FROM annotationdb.dict_rank_toxicity dr
LEFT JOIN annotationdb.pubchem_compounds pc ON dr.pubchem_cid = pc.cid
WHERE pc.cid IS NULL
LIMIT 20;

-- 5f. antibody_drug_conjugates
SELECT 'antibody_drug_conjugates: total rows' AS check_type,
       COUNT(*) AS row_count,
       COUNT(DISTINCT adc_id) AS unique_adc_ids
FROM annotationdb.antibody_drug_conjugates;
