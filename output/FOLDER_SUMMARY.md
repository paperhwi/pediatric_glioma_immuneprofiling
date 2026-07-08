# Folder Understanding — Pediatric HGG Immune Ecotype Project (Open PBTA)

Project: pediatric high-grade glioma (OpenPedCan v15) immune microenvironment
"ecotype" classification + survival analysis. Manuscript target: *Cancers* (MDPI),
currently in **major revision** stage. This folder = analysis output tree, not source repo.

## Cohort
- n = 349 primary-analysis samples (main cohort), plus secondary pools (n=702 full,
  n=322 BRAF/ALT LGG, n=31 BRAF/ALT main).
- Molecular subgroups: DMG_K27 (175), pHGG_WT (125), DHG_G34 (31), IHG (18).
- Confirmed **0 IDH-mutant cases** (Phase 1 audit) — pure IDH-wildtype pediatric cohort.

## Core finding: 3 immune ecotypes
Consensus clustering (k=3, LM22/CIBERSORTx + ssGSEA features) → three TME classes:
- **Lymphocyte-inflamed** (~113/349, 31.8%)
- **Myeloid-dominant / Intermediate** (~157-161/349, ~46%)
- **Immune-desert** (~78-79/349, 22.3%)

Ecotype assignment reproducible across methods (LM22 CIBERSORTx, ssGSEA, quanTIseq,
xCell, EPIC, MCP-counter) and robust to TPM-vs-original normalization (revC: ARI=0.634,
86.8% per-sample agreement).

## Survival signal
- Overall log-rank P = 0.028; median OS: Inflamed 595d > Myeloid-dom 465d > Desert 403d.
- Multivariable Cox (n=251, ecotype+cohort+age+sex): Immune-desert HR ≈ 1.75-1.79
  vs Inflamed (P≈0.004), Lymphocyte-inflamed protective HR≈0.55 (P=0.003).
- Robust to: anatomical location adjustment (HR 1.75→1.77), IPW missingness sensitivity
  (HR 1.79→1.88, missingness 28.1% overall, worst in DHG_G34 51.6%), and PERMANOVA shows
  ecotype ≠ proxy for molecular subtype (residual R²=0.303, P=0.0002 after subtype adjustment).
- No significant ecotype×subtype interaction (Cox LR P=0.515) → effect generalizes best in
  pHGG_WT and DHG_G34; masked in DMG_K27 (uniformly poor prognosis) and underpowered in IHG.

## Location finding (Phase 2, counter-intuitive vs adult glioma)
- Midline tumors enriched for Inflamed ecotype (41.4% vs 25.0% hemispheric, OR=2.12, P=0.010).
- Hemispheric tumors enriched for Immune-desert (32.7% vs 15.7%, P=0.002).
- Location does NOT add prognostic info beyond ecotype (Cox LR P=0.581) — orthogonal axis.
- Ecotype effect significant within Midline (P=0.033) and Hemispheric (P=0.030) strata.

## Other analysis threads present in folder
- **DEG/pathway (Step 8)**: DESeq2 + gprofiler per ecotype (Immune-desert, Lymphocyte-inflamed,
  Myeloid-dominant marker sets, ~200 top markers each).
- **BRAF/ALT projection (Step 10)**: ecotype composition differs sharply by molecular family
  (chi-sq P=3.5e-8); MAPK pathway activity varies strongly across families (P=4.9e-23).
- **Checkpoint/TIS/CAR-T target expression**: checkpoint_expression_by_ecotype,
  cart_target_expression_by_ecotype, Bagaev TIS signature work (Step 4d-4g).
- **V600E-specific sub-analysis** (Step 4e-4g): BRAF V600E LGG replication, bootstrap/permutation
  validation, C1 DEG cluster, CDKN2A/B Fisher test, PROGENy pathway scores.
- **Proliferation control**: checks ecotype differences aren't just proliferation artifacts.
- **Sensitivity/robustness supplements**: SupplS3 (high-confidence subset), SupplS4 (full-pool
  sensitivity), extent-of-resection (EOR) sensitivity.

## Reviewer-response revisions (most recent work, rev A/C/D — likely "major revision" response)
- **revA — PERMANOVA + Cox interaction**: proves ecotype signal independent of molecular
  subtype (see Survival signal section above).
- **revC — TPM/FPKM harmonization**: fixes a Methods documentation error (text said FPKM,
  analysis was actually TPM); independent Python (gseapy) re-implementation confirms
  median Spearman rho=0.958 concordance with original R/GSVA ssGSEA.
- **revD — IPW sensitivity for OS missingness**: confirms Immune-desert HR robust to
  missing-not-at-random assumptions via inverse-probability weighting.
- All three ship regenerated 300dpi figures + manuscript-ready action-item paragraphs.

## Deliverables ready for manuscript insertion
- `Consolidated_Phase1-4_Report.md` / `.docx` / `.pdf` — reviewer-response writeup (IDH audit,
  location analysis, master annotation table, combined figure).
- `phase2-4_location_master_figure/sample_master_annotation.xlsx` (4 sheets) +
  `sample_master.h5ad` — full annotated cohort object.
- `Manuscript_Draft_300dpi_feedback_revised.docx` — main manuscript draft.
- `Pediatric_Glioma_Immune_Ecotype_Talk(_v2).pptx` — conference/lab talk slides.
- `발표대본_소아교종_면역생태형.docx` — Korean presentation script (binary, not parsed here).
- `소아고등급교종.docx` — Korean-titled doc, likely background/notes (binary, not parsed).
- `additional analysis result 300 dpi/` — final Arial 300dpi figure set (PERMANOVA, subtype
  KM/forest, IPW forest/balance, TPM concordance, location KM, molecular×location figure) +
  `Session_Analysis_Figures_Report.docx` cover doc.
- `email_draft_PI_20250625.md` — draft email to PI (not read in detail; flag if content needed).

## File-type map (folder mechanics, not content)
- `Step*.ipynb`, `Suppl*.ipynb`, `rev*.ipynb` — analysis notebooks per pipeline stage.
- `*_summary.json` — machine-readable stats per stage (used above).
- `figs_*/`, `figures_300dpi/` — per-stage and publication-ready plot exports.
- `ccp_main_k*.rds`, `consensus_*_k*.npz` — cached clustering objects (k=2..6 sweep).
- `tpm_for_cibersortx.tsv(.gz)` (180MB) — full TPM input matrix for CIBERSORTx deconvolution.
- `*.log` — run logs (deconv, clustering, survival, pathway, cibersortx submission).

## Not read in detail (flag if you need specifics)
- Binary `.docx` reports (Analysis_logic_review, CIBERSORTx_brain_matrix_feasibility,
  Comprehensive_Step1-8_Report, Survival_Pathway_Report, Step2-5_Report, Table1*.docx,
  Korean-titled docs) — only filenames/dates inspected, not text content.
- Individual notebook cell contents (relied on their `*_summary.json` outputs instead).
