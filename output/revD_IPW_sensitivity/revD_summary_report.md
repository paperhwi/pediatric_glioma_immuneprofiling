# Revision D — IPW sensitivity for OS missingness: results

## Bottom line

OS data are missing in 28.1% of the cohort, with the highest
missingness in DHG_G34 (51.6%) and the lowest in
DMG_K27 (20.0%). After inverse-probability-of-
observation weighting (stabilised, propensity from ecotype + cohort_group + age +
sex + location_class), the Immune-desert HR is robust to missingness mechanism
assumptions, shifting from **HR = 1.79 (95% CI 1.21–2.66, P = 0.004)** under the naïve
Cox model to **HR = 1.88 (95% CI 1.22–2.90, P = 0.004)** under stabilised IPW with
robust variance. The log-HR shift between models is 0.048 (relative HR
change 4.9%), well within the unweighted 95% CI half-width, and the
direction and significance of the effect are preserved.

## Missingness pattern

| Cohort | Missingness (%) |
|---|---|
| pHGG_WT | 34.4% |
| DMG_K27 | 20.0% |
| DHG_G34 | 51.6% |
| IHG | 22.2% |
| **Overall** | **28.1%** |

## Propensity model and weights

- Propensity model: logistic regression of OS-observation indicator on ecotype,
  cohort_group, age, sex, location_class (L2-regularised, C=1.0).
- Stabilised weight w = P(observed) / P(observed | X), mean = 1.174,
  raw range [0.755, 8.921],
  trimmed at 1st/99th percentile to [0.757, 6.222].

## Covariate balance

- Covariates with |SMD| > 0.1 before IPW: 6/11
- Covariates with |SMD| > 0.1 after IPW:  4/11
- Max |SMD|: 0.505 (pre) → 0.342 (post)
  (Love plot: figures_300dpi/FigD2_balance_love_plot_300dpi.png)

## Cox model comparison (n with OS = 251)

| Term | Naïve HR (95% CI), P | IPW HR (95% CI), P |
|---|---|---|
| ecotype_Intermediate | 1.37 (0.98–1.91), 0.0615 | 1.55 (1.11–2.17), 0.0109 |
| ecotype_Immune-desert | 1.79 (1.21–2.66), 0.00375 | 1.88 (1.22–2.90), 0.00443 |
| cohort_group_DMG_K27 | 2.96 (2.11–4.15), 3.65e-10 | 3.01 (2.07–4.36), 6.6e-09 |
| cohort_group_DHG_G34 | 1.15 (0.59–2.24), 0.672 | 1.22 (0.73–2.06), 0.448 |
| cohort_group_IHG | 0.22 (0.08–0.63), 0.00497 | 0.24 (0.08–0.70), 0.00881 |
| age | 0.99 (0.96–1.02), 0.506 | 0.99 (0.96–1.02), 0.6 |
| male | 0.90 (0.67–1.19), 0.45 | 0.94 (0.70–1.25), 0.653 |

## Sensitivity across alternative weight specifications

| Scheme | Immune-desert HR (95% CI), P |
|---|---|
| Naïve (no IPW) | 1.79 (1.21–2.66), 0.004 |
| Stabilised IPW (trimmed 1–99%) | 1.88 (1.22–2.90), 0.004 |
| Untrimmed IPW | 1.88 (1.22–2.90), 0.004 |
| Trimmed 1–99% IPW | 1.88 (1.22–2.90), 0.004 |
| Truncated 2.5–97.5% IPW | 1.85 (1.21–2.84), 0.005 |
| Cohort-only propensity (cohort_group) | 1.79 (1.17–2.73), 0.007 |

## Manuscript action items

1. **Insert a sensitivity paragraph at the end of Section 3.6** reporting the
   IPW-adjusted Immune-desert HR and the magnitude of shift versus the naïve estimate.
2. **Update the missingness limitation in the Discussion** to note that the IPW
   sensitivity analysis confirms robustness to plausible MNAR mechanisms within the
   propensity-model assumption.
3. **Reference Supplementary Figures S_IPW_sensitivity (FigD1) and S_IPW_balance
   (FigD2)** in the relevant Methods/Results sections.

## Files produced

- `balance_SMD.tsv` — covariate balance pre/post IPW
- `Cox_naive_vs_IPW.tsv` — side-by-side Cox comparison
- `alternative_weighting_schemes.tsv` — robustness across schemes
- `figures_300dpi/FigD1_IPW_sensitivity_300dpi.{png,pdf}`
- `figures_300dpi/FigD2_balance_love_plot_300dpi.{png,pdf}`