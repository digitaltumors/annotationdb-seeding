-- ============================================================================
-- SQL AUDIT: annotationdb (Old) vs annotationdb_5 (New)
-- Strictly read-only — no data modifications.
-- ============================================================================

-- 1. ROW-COUNT COMPARISON
SELECT
    'pubchem_compounds'          AS table_name,
    (SELECT COUNT(*) FROM annotationdb.pubchem_compounds)        AS old_db_count,
    (SELECT COUNT(*) FROM annotationdb_5.pubchem_compounds)      AS new_db_count,
    (SELECT COUNT(*) FROM annotationdb_5.pubchem_compounds)
      - (SELECT COUNT(*) FROM annotationdb.pubchem_compounds)    AS diff
UNION ALL SELECT
    'compound_synonyms',
    (SELECT COUNT(*) FROM annotationdb.compound_synonyms),
    (SELECT COUNT(*) FROM annotationdb_5.compound_synonyms),
    (SELECT COUNT(*) FROM annotationdb_5.compound_synonyms)
      - (SELECT COUNT(*) FROM annotationdb.compound_synonyms)
UNION ALL SELECT
    'bioassays',
    (SELECT COUNT(*) FROM annotationdb.bioassays),
    (SELECT COUNT(*) FROM annotationdb_5.bioassays),
    (SELECT COUNT(*) FROM annotationdb_5.bioassays)
      - (SELECT COUNT(*) FROM annotationdb.bioassays)
UNION ALL SELECT
    'compound_bioassays',
    (SELECT COUNT(*) FROM annotationdb.compound_bioassays),
    (SELECT COUNT(*) FROM annotationdb_5.compound_bioassays),
    (SELECT COUNT(*) FROM annotationdb_5.compound_bioassays)
      - (SELECT COUNT(*) FROM annotationdb.compound_bioassays)
UNION ALL SELECT
    'toxicity',
    (SELECT COUNT(*) FROM annotationdb.toxicity),
    (SELECT COUNT(*) FROM annotationdb_5.toxicity),
    (SELECT COUNT(*) FROM annotationdb_5.toxicity)
      - (SELECT COUNT(*) FROM annotationdb.toxicity)
UNION ALL SELECT
    'chembl_mechanism',
    (SELECT COUNT(*) FROM annotationdb.chembl_mechanism),
    (SELECT COUNT(*) FROM annotationdb_5.chembl_mechanism),
    (SELECT COUNT(*) FROM annotationdb_5.chembl_mechanism)
      - (SELECT COUNT(*) FROM annotationdb.chembl_mechanism)
UNION ALL SELECT
    'chembl_drug_indication',
    0,
    (SELECT COUNT(*) FROM annotationdb_5.chembl_drug_indication),
    (SELECT COUNT(*) FROM annotationdb_5.chembl_drug_indication)
UNION ALL SELECT
    'cellosaurus_cell_lines',
    (SELECT COUNT(*) FROM annotationdb.cellosaurus_cell_lines),
    (SELECT COUNT(*) FROM annotationdb_5.cellosaurus_cell_lines),
    (SELECT COUNT(*) FROM annotationdb_5.cellosaurus_cell_lines)
      - (SELECT COUNT(*) FROM annotationdb.cellosaurus_cell_lines)
UNION ALL SELECT
    'cell_line_synonyms',
    (SELECT COUNT(*) FROM annotationdb.cell_line_synonyms),
    (SELECT COUNT(*) FROM annotationdb_5.cell_line_synonyms),
    (SELECT COUNT(*) FROM annotationdb_5.cell_line_synonyms)
      - (SELECT COUNT(*) FROM annotationdb.cell_line_synonyms)
UNION ALL SELECT
    'cell_line_disease',
    (SELECT COUNT(*) FROM annotationdb.cell_line_disease),
    (SELECT COUNT(*) FROM annotationdb_5.cell_line_disease),
    (SELECT COUNT(*) FROM annotationdb_5.cell_line_disease)
      - (SELECT COUNT(*) FROM annotationdb.cell_line_disease)
UNION ALL SELECT
    'oncotree',
    (SELECT COUNT(*) FROM annotationdb.oncotree),
    (SELECT COUNT(*) FROM annotationdb_5.oncotree),
    (SELECT COUNT(*) FROM annotationdb_5.oncotree)
      - (SELECT COUNT(*) FROM annotationdb.oncotree)
UNION ALL SELECT
    'antibody_drug_conjugates',
    (SELECT COUNT(*) FROM annotationdb.antibody_drug_conjugates),
    (SELECT COUNT(*) FROM annotationdb_5.antibody_drug_conjugates),
    (SELECT COUNT(*) FROM annotationdb_5.antibody_drug_conjugates)
      - (SELECT COUNT(*) FROM annotationdb.antibody_drug_conjugates)
UNION ALL SELECT
    'adc_indications',
    0,
    (SELECT COUNT(*) FROM annotationdb_5.adc_indications),
    (SELECT COUNT(*) FROM annotationdb_5.adc_indications)
UNION ALL SELECT
    'substances',
    (SELECT COUNT(*) FROM annotationdb.substances),
    (SELECT COUNT(*) FROM annotationdb_5.substances),
    (SELECT COUNT(*) FROM annotationdb_5.substances)
      - (SELECT COUNT(*) FROM annotationdb.substances)
UNION ALL SELECT
    'substance_synonyms',
    (SELECT COUNT(*) FROM annotationdb.substance_synonyms),
    (SELECT COUNT(*) FROM annotationdb_5.substance_synonyms),
    (SELECT COUNT(*) FROM annotationdb_5.substance_synonyms)
      - (SELECT COUNT(*) FROM annotationdb.substance_synonyms)
UNION ALL SELECT
    'substance_toxicity',
    (SELECT COUNT(*) FROM annotationdb.substance_toxicity),
    (SELECT COUNT(*) FROM annotationdb_5.substance_toxicity),
    (SELECT COUNT(*) FROM annotationdb_5.substance_toxicity)
      - (SELECT COUNT(*) FROM annotationdb.substance_toxicity)
UNION ALL SELECT
    'diril_toxicity',
    (SELECT COUNT(*) FROM annotationdb.diril_toxicity),
    (SELECT COUNT(*) FROM annotationdb_5.diril_toxicity),
    (SELECT COUNT(*) FROM annotationdb_5.diril_toxicity)
      - (SELECT COUNT(*) FROM annotationdb.diril_toxicity)
UNION ALL SELECT
    'dict_rank_toxicity',
    (SELECT COUNT(*) FROM annotationdb.dict_rank_toxicity),
    (SELECT COUNT(*) FROM annotationdb_5.dict_rank_toxicity),
    (SELECT COUNT(*) FROM annotationdb_5.dict_rank_toxicity)
      - (SELECT COUNT(*) FROM annotationdb.dict_rank_toxicity)
UNION ALL SELECT
    'atc_codes',
    (SELECT COUNT(*) FROM annotationdb.atc_codes),
    (SELECT COUNT(*) FROM annotationdb_5.atc_codes),
    (SELECT COUNT(*) FROM annotationdb_5.atc_codes)
      - (SELECT COUNT(*) FROM annotationdb.atc_codes);
