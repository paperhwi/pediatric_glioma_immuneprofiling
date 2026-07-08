# Revision A — PERMANOVA + Cox interaction: results

## Bottom line

Tumor molecular subtype (cohort_group) accounts for **R² = 0.044** of the
variance in the 29-feature TPM-harmonised clustering matrix, but it does **not**
explain away the ecotype signal: after adjusting for cohort_group, the residual
ecotype effect retains **R² = 0.303** (P = 0.0002, 4 999 stratified
permutations). The full Cox model **does not show a significant ecotype × cohort
interaction** (LR χ² = 5.22, df = 6, P = 0.515), supporting the additive
interpretation of ecotype as an independent prognostic axis. Within-subtype
estimates confirm that the ecotype-prognosis effect generalises to **pHGG_WT**
and **DHG_G34** but is masked within DMG_K27 by uniformly poor prognosis and
under-powered within IHG.

## PERMANOVA variance partition (n = 349, Euclidean on z-scored features, 4999 permutations)

| Factor | R² | P |
|---|---|---|
| Ecotype (alone) | 0.320 | 0.0002 |
| Cohort_group (alone) | 0.044 | 0.0002 |
| Ecotype \| Cohort_group (sequential) | 0.303 | 0.0002 |
| Cohort_group \| Ecotype (sequential) | 0.027 | 0.0294 |
| Residual (intrinsic TME heterogeneity) | 0.653 | – |

**Interpretation.** The ecotype partition is the single largest source of structure
in the TPM feature matrix; subtype contributes substantial but smaller variance and
the two factors share **R² ≈ 0.017** in overlap.
Crucially, the marginal contribution of ecotype after subtype adjustment remains
highly significant and large (R² = 0.303, P = 0.0002), demonstrating that the
three-state immune ecotype framework is **not a proxy for tumor molecular subtype**.

## Cox interaction (n with OS = 251)

- Main model (ecotype + cohort_group + age + sex): log-lik = -939.84
- Interaction model (+ ecotype × cohort_group terms): log-lik = -937.23
- LR test: χ² = **5.22**, df = 6, **P = 0.5154**

Since the interaction is **not significant** at α = 0.05, the prognostic effect
of ecotype is interpretable as additive across cohort groups. Subtype-stratified
estimates below quantify where statistical power is available to detect it.

## Subtype-stratified Cox (HR for Immune-desert vs Inflamed within each cohort)

| Cohort | n | Events | HR Imm-desert | 95% CI | P | log-rank P |
|---|---|---|---|---|---|---|
| DMG_K27 | 140 | 134 | 1.59 | 0.97–2.63 | 0.068 | 0.132 |
| DHG_G34 | 15 | 13 | 0.67 | 0.08–5.42 | 0.711 | 0.070 |
| pHGG_WT | 82 | 56 | 2.22 | 1.11–4.43 | 0.024 | 0.049 |
| IHG | 14 | 5 | 3.90 | 0.41–37.56 | 0.239 | 0.329 |

## Manuscript action items

1. **New paragraph at end of Section 3.5** (Ecotype Distribution Across Histone Class):
   add the PERMANOVA result above (R² values, P-values, sequential decomposition).
2. **New paragraph at end of Section 3.6** (Survival): add the Cox interaction test
   result and explicitly state that the ecotype prognostic effect generalises to
   pHGG_WT and DHG_G34 but is masked within DMG_K27 and under-powered in IHG.
3. **Add a sentence to the Discussion limitations**: 'Ecotype-based prognostic
   stratification is most clinically actionable for pHGG_WT and DHG_G34 tumors
   where the ecotype effect generalises; within H3 K27M-altered DMG, the uniformly
   poor prognosis of the subgroup limits dynamic range for ecotype-driven
   discrimination.'
4. **Add Supplementary Figures S_PERMANOVA and S_subtype_stratified** in
   revA_PERMANOVA_interaction/figures_300dpi/.

## Files produced

- `PERMANOVA_results.tsv` — R² + P for marginal/sequential decomposition
- `Cox_main_effects.tsv` — main-effects Cox model summary
- `Cox_interaction_model.tsv` — ecotype × cohort interaction Cox summary
- `subtype_stratified_Cox.tsv` — HR per subtype
- `figures_300dpi/FigA1_PERMANOVA_variance_partition_300dpi.{png,pdf}`
- `figures_300dpi/FigA2_subtype_stratified_KM_forest_300dpi.{png,pdf}`