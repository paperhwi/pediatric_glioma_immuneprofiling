# Revision bundle for the final cohort n = 343

This directory contains the current November 2026 revision. The original n = 349 September submission in the parent directory remains a historical record.

## Current definitions

Six specimens outside the pediatric-type diffuse high-grade glioma cohort are excluded from the final analyses, including clustering. The retained 46-feature standardized input is restricted to the 343 eligible specimens; this is not a claim that all original features were newly standardized after exclusion.

Immune-state labels were re-derived with the locked consensus procedure and assigned by their mode over seeds 42–66: **Immune-high 122, Immune-intermediate 144, Immune-low 77**. Historical label strings in source-loading or comparison code are provenance identifiers, not alternative current state names.

The primary survival specification remains immune state plus molecular group, age, and sex. One origin-review rule is applied throughout the cohort, with four records re-timed. Aggregate sensitivity results are retained regardless of significance. Individual survival durations are not included in this repository.

## Contents

- `analysis_code/`: the delivered analysis sources, including shared definitions and numeric pipeline stages. Original absolute paths require adaptation and some stages depend on intermediate outputs not bundled here.
- `reproducibility_data/`: final assignments, aggregate results, differential-expression tables, clustering diagnostics, seed records, and the software/parameter workbook.
- `figures/`: delivered main and supplementary PDFs with their original analysis filenames.
- `validate_revision_bundle.py`: checks consistency of the stored public outputs.

## Publication numbering crosswalk

The cleaned manuscript package uses these delivery labels; repository analysis filenames remain stable for provenance.

| Submission label | Repository analysis figure |
|---|---|
| Figure S1 | `FigureS_survival.pdf` |
| Figure S2 | `FigureS_state_structure.pdf` |
| Figure S3 | `FigureS_immune_features.pdf` |
| Figure S4 | ScPCA axis-linkage figure supplied in the manuscript package; no direct n343 rendering source is claimed here |
| Supplementary Table S1 | Detailed n343 analysis workbook supplied with the manuscript |
| Supplementary Table S2 | Software/parameter workbook retained here under its original `Supplementary_Table_S25_software_parameters_seeds_n343.xlsx` filename |
| Supplementary Table S3 | [`additional_analyses_20261103/supplementary/Supplementary_Table_S3_gene_sets_and_sources.xlsx`](additional_analyses_20261103/supplementary/Supplementary_Table_S3_gene_sets_and_sources.xlsx) |
| Supplementary Table S4 | Single-cell per-patient summaries supplied with the manuscript |
| Figure S5 | [`additional_analyses_20261103/figures/FigureS5_standardization_consistency.pdf`](additional_analyses_20261103/figures/FigureS5_standardization_consistency.pdf) |
| Figure S6 | [`additional_analyses_20261103/figures/FigureS6_axis_library_survival.pdf`](additional_analyses_20261103/figures/FigureS6_axis_library_survival.pdf) |
| Supplementary Table S6 | [`additional_analyses_20261103/supplementary/Supplementary_Table_S6_standardization_library_survival.xlsx`](additional_analyses_20261103/supplementary/Supplementary_Table_S6_standardization_library_survival.xlsx) |
| Supplementary Table S5 | Long-observation and survival-origin-correction records supplied with the manuscript; individual survival records are not redistributed here |

## Validation and analytic execution

```bash
python cancers_submission_2026/revision_2026_11_n343/validate_revision_bundle.py
```

This validator checks stored primary and follow-up results and does not fit models. The two additional sensitivity analyses were executed separately; the original primary results and gene membership remain unchanged.

Full analytic execution requires source OpenPedCan v15 expression and clinical files, GSE268577, ScPCA SCPCP000001, and intermediate files identified by the scripts. Do not assume that this folder is a standalone end-to-end rerun environment. Before corrected-OS survival stages, the uniform audit in `s02_os_audit.py` must derive the local, git-ignored `restricted_os_fix.tsv`; subsequent survival scripts must run in a fresh process so `common.py` reads the corrections. Without the required inputs and corrections, stored-output validation is the supported check.

The additional re-standardization and continuous-axis library-control results are supplied with the revised manuscript and available in [`additional_analyses_20261103/`](additional_analyses_20261103/). Re-standardization preserves 334/343 modal labels (97.4%; ARI 0.915), with 116/149/78 sensitivity counts. The library-adjusted axis estimate attenuates to HR 0.861 (95% CI 0.726–1.022), P 0.088. Primary labels remain 122/144/77.

Corrected-OS loading now fails explicitly if the four locally derived corrections are absent. Run the uniform audit and restart the analysis process; released durations must not silently stand in for corrected outcomes.
