# Pediatric glioma immune profiling

Reproducibility materials for **Transcriptomic Immune States and Survival in Pediatric Diffuse High-Grade Gliomas**.

## Current revision cohort n = 343

[`cancers_submission_2026/revision_2026_11_n343/`](cancers_submission_2026/revision_2026_11_n343/) is the current November 2026 revision bundle.

- One RNA specimen per participant; six ineligible specimens are excluded from the final cohort.
- Final states: Immune-high 122, Immune-intermediate 144, and Immune-low 77.
- Survival: the submitted molecular-group, age, and sex model is retained; the same survival-origin rule is applied throughout the outcome-available cohort.
- The primary clustering procedure uses the retained standardized 46-feature input restricted to eligible specimens. It must not be described as a new standardization of every input on n = 343.

Validate the delivered outputs from the repository root:

```bash
python cancers_submission_2026/revision_2026_11_n343/validate_revision_bundle.py
```

This checks stored cohort assignments, state counts, exclusions, adjusted-expression counts, and shared genes. It does not rerun the analytic pipeline.

## September submission archive

The original n = 349 materials under [`cancers_submission_2026/`](cancers_submission_2026/) are retained as the historical submission record. Their states, cohort counts, and validation command must not be substituted for the current revision.

## Reproduction inputs

Raw expression, restricted clinical data, individual survival times, and histology records are not redistributed. Public sample-level derived files use biospecimen identifiers. Full reproduction requires the documented OpenPedCan, GEO, and ScPCA inputs and intermediate files. See the [current bundle README](cancers_submission_2026/revision_2026_11_n343/README.md) for scope and execution requirements.

## Executed follow-up analyses

[Within-cohort scaling and RNA-library survival sensitivity analyses](cancers_submission_2026/revision_2026_11_n343/additional_analyses_20261103/) include public inputs for clustering, aggregate Cox/PH outputs, Figures S5/S6 and Supplementary Table S6. Primary labels remain 122/144/77; re-standardized 116/149/78 labels are a sensitivity partition. Library adjustment attenuates the axis association to HR 0.861 (95% CI 0.726–1.022), P = 0.088.


## Final version v2.0

The [v2.0 release](https://github.com/paperhwi/pediatric_glioma_immuneprofiling/releases/tag/v2.0) is the version-fixed final n343 analysis. The [v1.0 tag](https://github.com/paperhwi/pediatric_glioma_immuneprofiling/tree/v1.0) preserves the original September submission snapshot. [Reviewer follow-up](cancers_submission_2026/revision_2026_11_n343/reviewer_followup_20261103/) includes additional categorical Cox and IPW models, the 702-specimen scaling reference, signature provenance, five-year KM displays, and submitted Tables 2–5 retained as Supplementary Table S8. Raw and restricted patient-level clinical or survival records are not newly redistributed.
