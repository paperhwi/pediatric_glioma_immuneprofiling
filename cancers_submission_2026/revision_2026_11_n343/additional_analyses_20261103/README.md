# Within-cohort scaling and RNA-library survival sensitivity analyses

These November 3, 2026 follow-up analyses supplement the primary n = 343 revision. **Primary assignments remain Immune-high 122 / Immune-intermediate 144 / Immune-low 77.** The re-standardized 116 / 149 / 78 assignments are a preprocessing sensitivity partition; primary downstream state results have not been recomputed with those labels.

## Results

- Re-standardizing all 46 features within the 343 eligible participants retained **334/343 modal labels (97.4%; ARI 0.914948)** and changed nine. Median pairwise seed ARI increased from 0.941602 to 0.990156. All-25-seed agreement increased from 311 to 339 participants.
- PAC and consensus silhouette still favor k = 2; the gap first-one-standard-error rule selects k = 3. Assignment and k-diagnostic protocols are distinct, as in the primary analysis.
- The fixed 21-immune-program axis reproduced the primary Cox estimate, HR 0.827260 (95% CI 0.716755–0.954802), P = 0.009537, in 245 patients / 206 deaths.
- Adding RNA-library type attenuated this to **HR 0.861443 (0.725813–1.022417), P = 0.087947**. Library stratification gave HR 0.854557, P = 0.079678. After autopsy exclusion, the library-adjusted estimate was HR 0.942358, P = 0.574971 (189 patients / 150 deaths).

The continuous axis is not established as a survival factor independent of RNA-library preparation. Directional agreement does not remove specimen-timing or technical sensitivity. Exact estimates and PH diagnostics are retained regardless of significance.

## Contents and numbering

- `code/`: executed analysis methods with portable paths and a separate restricted-input connector.
- `inputs/`: derived 46-feature matrices, pseudonymous state/group/library/axis annotations, and public baseline reference outputs. No individual survival durations or status are included.
- `tables/`: assignments, seed comparisons, feature-scaling audit, k/gap diagnostics, aggregate Cox estimates, coefficients, covariance and PH diagnostics.
- `figures/`: publication Figures **S5** (standardization consistency) and **S6** (axis/library survival models), each as SVG/PDF/PNG.
- `supplementary/Supplementary_Table_S3_gene_sets_and_sources.xlsx`: clarified source metadata with gene membership unchanged.
- `supplementary/Supplementary_Table_S6_standardization_library_survival.xlsx`: the final nine-sheet workbook, including parameters and environment.

See [the revision README](../README.md) and [manuscript-alignment notes](MANUSCRIPT_ALIGNMENT.md) for the distinction between the original primary analyses and these sensitivity results.

## Run and validate

From the repository root, validate stored outputs without fitting models:

```bash
python cancers_submission_2026/revision_2026_11_n343/validate_revision_bundle.py
```

Use an isolated environment with `code/requirements.txt`. To rerun clustering from the public derived inputs:

```bash
python cancers_submission_2026/revision_2026_11_n343/additional_analyses_20261103/code/cluster_check.py
```

This runs the paired 25-seed procedures, k2–10 diagnostics and 50-reference gap calculation. Per-seed working files are ignored under `.cache`; consolidated published outputs are in `tables`.

Survival fitting requires a corrected 343-row non-public covariate/outcome TSV. Set `PGI_RESTRICTED_SURVIVAL_INPUT` to that file **outside the checkout**, then run `code/axis_library_survival.py`; `code/verify_cox_math.py` independently checks the Efron likelihood, score and information matrix. The required schema is documented in [code/README.md](code/README.md). Clinical design matrices and event-time residuals are written beside the restricted input, outside the checkout. They are not public outputs.

`code/summarize_and_plot.py` regenerates the two figures from public result tables. This public bundle supports clustering/plot reproduction and aggregate-output verification; full survival refitting additionally requires the authorized source-derived clinical input.
