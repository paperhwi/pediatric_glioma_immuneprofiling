# Cancers resubmission — figure & supplementary reorganisation

Source: `E:\Open PBTA\Revision` (REVISION PACKAGE 260916 + Week1 + Week2 + FINAL MANUSCRIPT 260722).
Target: new submission to *Cancers* (MDPI). Manuscript text written by the author; this task
delivers figures, supplementary material, legends and numbering only.

GitHub archive for excluded items: https://github.com/paperhwi/pediatric_glioma_immuneprofiling

---

## MAIN FIGURES — final composition

| New | Panels | Source | Action |
|---|---|---|---|
| **Fig 1** | A study workflow; B analysis-set derivation | old Figure1 + Figure1B_flow_biorender | **redraw both**, one vector file. 1A: `k = 2–10` evaluation, *bootstrap consensus clustering*, "k = 2 most stable; k = 3 retained as an exploratory three-state representation". 1B redrawn in matching style, US spelling. |
| **Fig 2** | A LM22 QC; B cross-method concordance n=349; C cross-method estimates by ecotype | old 2A, 2B, 2D | reuse vector; old **2C → Suppl S1** |
| **Fig 3** | A PAC; B silhouette; C relative ΔAUC; D gap statistic; E minimum cluster size; F PCA/UMAP; G defining-feature heatmap | new A–E (regenerate) + old 3B, 3C | ΔAUC annotation → `largest remaining relative ΔAUC (+29.8%)`. old **3D dropped** (duplicated by Fig 4) |
| **Fig 4** | A–D molecular group, location, within-group, sequential PERMANOVA | old Figure4 | unchanged |
| **Fig 5** | A single-cell signature AUROC; B cell-type source of ecotype-axis genes; C spot-to-centroid assignment, 3 sections; D spatial coherence + Moran's I | **new**, from T_S25/S26/S27/S28 + vector crop of S21 A–C | conservative wording (see below) |
| **Fig 6** | A–C adjusted volcanoes; D overlap (publication style, not Venn); E gradient; F paired effect sizes | regenerate from `adjusted_DEG_all_contrasts.tsv` | old **6B → Suppl S8** |
| **Fig 7** | A overall KM (numbers at risk + censor marks); B multivariable Cox | old 7A, 7C | old **7B → Suppl S10** |

## SUPPLEMENTARY FIGURES — renumbered S1…S15

| New | Content | Old |
|---|---|---|
| S1 | Brain-tuned microglia/MDM signatures vs LM22 macrophage/monocyte axes | main Fig 2C |
| S2 | LM22 quality control in the same-resource pooled set (n = 702) | S4 |
| S3 | ssGSEA TPM-harmonisation concordance | S3 |
| S4 | Re-clustering of the high-confidence LM22 subset (n = 91) | S6 |
| S5 | Clustering-definition and held-out-feature sensitivity | S8 |
| S6 | PERMANOVA variance partition, ecotype vs molecular group | S5 |
| S7 | Bulk immune-programme characterisation (5 panels) | main Fig 5A–E |
| S8 | Pathway enrichment from the adjusted gene lists | main Fig 6B |
| S9 | Unadjusted tie-corrected one-versus-rest volcanoes | S9 |
| S10 | Subtype-stratified survival: KM (adequately sized strata) + within-subtype Cox | main Fig 7B + S10B |
| S11 | Location-stratified Kaplan–Meier | S11 |
| S12 | Scaled Schoenfeld residual diagnostics | S12 |
| S13 | Survival missingness and IPW sensitivity (4 panels) | S13A–D |
| S14 | Single-cell full analysis | S20 |
| S15 | Spatial full analysis | S21 |

**Excluded → GitHub repository:** old S1, S2, S7A–C (pooled n=702 re-clustering/concordance),
S14A–C (BRAF/RTK projection), S15 (checkpoint), S16 (CAR-T), S17 (duplicate DEG audit),
S18 (very small DHG/IHG strata KM), S19 (HSPC/progenitor), S10A (duplicate of new S10).

## REQUIRED PHRASING

- ✗ `predominantly myeloid cellular source` → ✓ *immune-cell-associated, with many genes preferentially expressed in myeloid cells*
- ✗ `every section contains all three ecotypes` → ✓ *spots mapping most closely to each of the three bulk ecotype centroids were represented within each section*
- ✗ `bulk mixing inflates` → ✓ *consistent with an aggregation contribution*
- ✗ `one shared gradient` → ✓ *consistent with an ordered shared immune-presence axis*
- ✗ `independent 1,000-replicate rerun` → ✓ *a separate 1,000-replicate rerun on the same cohort* / *resampling stability analysis*
- ✗ `k = 1 was rejected` → ✓ *k = 3 was selected under the Tibshirani one-standard-error rule; k = 1 was not selected*
- ✗ `last substantial gain` → ✓ *largest remaining relative ΔAUC (+29.8%)*

## DELIVERABLES

1. `Main figures/` — Figure 1–7, PNG 600 dpi + vector PDF
2. `Supplementary figures/` — Figure S1–S15, PNG 600 dpi + vector PDF
3. `Figure_and_Table_legends_Cancers.docx`
4. `Supplementary_Material_Cancers.docx` (figures + table legends, one file)
5. `Supplementary_Tables_Cancers.xlsx` (renumbered, with old→new crosswalk sheet)
6. `GitHub_archive/` — items removed from the submission
