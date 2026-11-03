# Reviewer follow-up for the final n343 revision

This addition preserves primary Immune-high / Immune-intermediate / Immune-low assignments of 122 / 144 / 77. The 116 / 149 / 78 partition remains a preprocessing sensitivity.

- `tables/categorical_state_sensitivity_cox.tsv`: both state contrasts, including re-standardized primary adjustment and simultaneous autopsy exclusion plus RNA-library covariate adjustment or strata. Sample/event counts and PH diagnostics are retained.
- `tables/ipw_*.tsv`: the submitted missing-outcome weighting analysis re-estimated on n343. Stabilized 1–99% capped IPW gives Immune-high HR 0.599544 (95% CI 0.387363–0.927948), P=0.021704. Residual imbalance and informative missingness remain possible.
- `supplementary/`: updated gene sources (S3), all new categorical/IPW estimates (S7), and explicit mapping and n343 updates of submitted Tables 2–5 (S8). Tables 2, 4, 5 were already present as S1 source sheets; Table 3 is restored descriptively across all four molecular groups, without testing sparse strata.
- `brain_immune_signatures_current.gmt` and `SIGNATURE_PROVENANCE.md`: neutral cytotoxic name, exact retained memberships, and honest publication-attribution status.
- `tables/Figure2C_CD8_boxplot_diagnostics.tsv` and `tables/Figure5A_source_AUROC.tsv`: lower-bound CD8 concentration and original-target AUROC display sources.

The original 46 features were standardized over the initial 702 RNA-specimen pool using feature means and sample SDs (ddof=1), before restriction to n349 and n343. The retained n343 matrix matches the corresponding original rows within 4.45e-16. It is not a within-n343 z-score matrix; that separate sensitivity matrix is supplied in the previous addition.

Figure 3 PCA uses the retained 343x46 input, PCA(n_components=2), mean centering, whiten=False and svd_solver=auto. UMAP uses the same matrix, n_components=2, n_neighbors=15, min_dist=0.3, Euclidean distance, random_state=42. Embeddings were fitted without state or survival labels.

Figure 7A and subtype KM panels in Figure S1 display 0–5 years. Fits and log-rank tests use full corrected follow-up. Existing five-year administrative-censoring sensitivity results remain a distinct analysis.

## Execution

Use the recorded Python environment (Python 3.12.14, pandas 2.3.3, numpy 2.5.3, scipy 1.18.1, scikit-learn 1.9.1, lifelines 0.30.3; matplotlib for plotting). The private connector checks exact cohort IDs, canonical state labels, molecular group, RNA library and fixed immune-axis values.

```powershell
$env:PGI_RESTRICTED_SURVIVAL_INPUT='C:\authorized-local-data\observations_n343.tsv'
python cancers_submission_2026/revision_2026_11_n343/reviewer_followup_20261103/code/run_categorical_cox.py
python cancers_submission_2026/revision_2026_11_n343/reviewer_followup_20261103/code/run_ipw.py
python cancers_submission_2026/revision_2026_11_n343/reviewer_followup_20261103/code/verify_cox_math.py
python cancers_submission_2026/revision_2026_11_n343/reviewer_followup_20261103/code/regenerate_survival_figures.py
```

The clinical input must remain outside the Git checkout. It uses the schema documented by the existing private connector, plus `immune_state` and `location_class`. No individual event time, censoring status, patient covariate matrix or survival correction value is saved in these public outputs. `generated_figures/` is ignored; stored public figure files are under the revision's `figures/` directory. Original primary sources retain some environment-specific paths and depend on source/intermediate inputs; this is not a claim of a fresh end-to-end rerun.

## Versions and publication files

`v1.0` points to the original September submission snapshot (the same commit as `cancers-final-submission-2026-09-22`). `v2.0` identifies the final n343 revision. Current manuscript composition is preserved locally; full manuscript DOCX and restricted supplements are not uploaded here. References [1]–[54] and the Data Availability statement are supplied separately to the authors.

The IPW propensity dummy coding explicitly retains Immune-low as the reference (historical Immune-desert); this avoids changing the L2 penalty merely because state labels were renamed.


## Changes from the submitted analysis

The submitted n349 cohort is reduced to 343 by six eligibility exclusions applied before clustering and downstream analyses. One RNA specimen per patient includes initial, recurrent/progressive and autopsy specimens. States are re-derived and named Immune-high, Immune-intermediate and Immune-low (122/144/77), rather than carrying forward the historical lymphocyte-inflamed/myeloid-dominant/desert labels. All 343 records receive the same outcome-availability screen and origin-review rules; four durations are re-timed without changing vital status, leaving 245 survival-evaluable patients with 206 deaths. Source correction rules and aggregate action counts are retained; individual corrections are local only.

The retained 702-pool standardized input remains primary. Within-343 re-standardization (116/149/78) is a separate sensitivity, now with the primary Cox adjustment reported. Combined autopsy/library categorical models, the continuous-axis library models, and IPW expose sensitivity to specimen selection and technical preparation. Embedding settings and signature-source limitations are explicit. Figure display changes do not truncate model follow-up. Submitted Tables 2–5 are mapped to S8; current source sheet identifiers remain stable even when prose and display labels use American English.
