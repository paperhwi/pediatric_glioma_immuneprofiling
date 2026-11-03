# Revision bundle — final cohort n = 343 (November 2026)

This folder supersedes the n = 349 submission bundle for the revised manuscript. The submission bundle (`../`) is kept unchanged as the record of what was submitted.

## What changed
- **One exclusion rule.** Six specimens that are not pediatric-type diffuse HGG are removed before any analysis, including clustering: five non-HGG histologies and one methylation-class pleomorphic xanthoastrocytoma (`T_excluded_specimens.tsv`).
- **Immune states re-derived.** The locked consensus-clustering procedure was re-run on n = 343. The final label is the modal label over seeds 42–66, the same procedure that reproduces the n = 344 labels exactly.
  - Result: **Immune-high 122 / Immune-intermediate 144 / Immune-low 77**; one tumour changes state.
  - Seed stability is in `T_seed_stability.tsv` and `T_seed_partitions_n343.tsv`.
- **Renamed states.**
  - Lymphocyte-inflamed → Immune-high
  - Myeloid-dominant → Immune-intermediate
  - Immune-desert → Immune-low
- **Overall survival.** OS is measured from the diagnosis of the sequenced HGG. One rule was applied to every record, and four records were re-timed (`T_OS_origin_audit_summary.tsv`). Individual survival times are not redistributed.
- **Primary survival model unchanged from submission.** Immune state + molecular group + age + sex; all sensitivity analyses are reported (`T_survival_sensitivity.tsv`).

## Contents
- `analysis_code/`: the full pipeline (`s00`–`s34`, run in numeric order).
  - `common.py` holds the shared definitions.
  - Paths point to a local OpenPedCan v15 copy (`/home/claude/opc`) and to this repository (`/home/claude/repo`); edit them before running.
  - Slow steps: `clust.py` (seeds 0–99), `s04`, `s09`, `s12`, `s15`, `s16`, `s20a`, `s22`.
- `reproducibility_data/`:
  - per-tumour state assignments (biospecimen IDs only);
  - adjusted differential expression and shared genes;
  - clustering criteria, Figure 2/4/5 statistics, immune axis and loadings;
  - aggregate survival tables;
  - Supplementary Table S25 (software, parameters, seeds, input checksums).
- `figures/`: main Figures 1–7 and the supplementary survival, structure and immune-feature figures (PDF).

## Validate
```bash
python cancers_submission_2026/revision_2026_11_n343/validate_revision_bundle.py
```

## Inputs not redistributed
- OpenPedCan v15 expression and clinical files.
- GSE268577 (Visium) and ScPCA SCPCP000001, available from GEO and the ScPCA portal.
- Restricted survival times.
