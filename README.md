# Pediatric_glioma_immuneprofiling

Analysis and reproducibility files for pediatric glioma immune profiling using brain-tuned
deconvolution and immune-program signatures.

## Publication link 

Publication link will be provided after acceptance

## Description 

We have performed immune deconvolution using signature matrix that we have created using brain specific immune cells. 

Pediatric glioma, which are known to be immune desert in nature, seemed to have distinct immune microenvironment after applying brain-tuned deconvolution 

By classifying the pediatric glioma according to their immune profile, we were able to find relevance with tumor location and other clinical metadata 

## How to use this repository 

- You can download signature matrix
- Check our immune profiling label for OpenPBTA v15 and perform your own Hazard ratio analysis or Kaplan Meier survival analysis
- Validate our calculation for checking the statistical integrity of this article

## Cancers submission (September 2026)

`cancers_submission_2026/` holds the display items of the submission to *Cancers* — Figures 1-7,
Supplementary Figures S1-S16, the legends, the supplementary document and the 25-table
supplementary workbook — together with figure-building code, analysis scripts, clean notebooks,
and a curated non-identifiable reproducibility bundle.

The same folder deposits, under `archived_analyses/`, the analyses that were part of the earlier
revision package but are not submitted: the same-resource pooled sensitivity analyses, the
BRAF/RTK-altered projection, the immune-checkpoint and CAR-T transcript panels, the
differential-expression audit copy, the very small survival strata, and the hematopoietic stem
and progenitor programs.

The public bundle deliberately excludes raw expression matrices, clinical-level datasets, survival
times/status, histology files, RDS/H5AD analysis objects, and other restricted inputs. Biospecimen
identifiers in the deposited matrices are study pseudonyms. Run
`python cancers_submission_2026/analysis_code/validate_public_bundle.py` to verify the locked
cohort size, ecotype counts, matrix dimensions, identifier agreement, and restricted-file guard.
