"""
Assemble final 'for submission' folder.
Structure:
  for submission/
    01_manuscript/        — main docx + source md
    02_table/             — Table 1 docx + tsv
    03_figures/           — all 300-DPI PNG + SVG figures (renamed: Fig1-7 + S3/S4)
    04_supplements/       — figure captions, references, analysis plan (docx)
    05_correspondence/    — PI email (docx)
    06_supplementary_data/— key TSV outputs (DEG, enrichment, Cox, BRAF projection, etc.)
    README.txt
"""
from __future__ import annotations
from pathlib import Path
import shutil, subprocess

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"
SUB  = BASE / "for submission"
if SUB.exists():
    print(f"NOTE: {SUB} already exists — overwriting contents")
SUB.mkdir(exist_ok=True)

# Sub-folders
folders = ["01_manuscript", "02_table", "03_figures", "04_supplements",
           "05_correspondence", "06_supplementary_data"]
for f in folders: (SUB / f).mkdir(exist_ok=True)

PANDOC = "pandoc"

def md_to_docx(src: Path, dst: Path):
    subprocess.run([PANDOC, str(src), "-o", str(dst), "--from", "markdown", "--to", "docx",
                    "--reference-doc", str(OUT / "Table1_n349.docx")], check=True)
    print(f"  md→docx: {dst.relative_to(SUB)}")

def copy_file(src: Path, dst: Path):
    shutil.copy2(src, dst)
    print(f"  copy   : {dst.relative_to(SUB)}")

# -----------------------------------------------------------------
# 01_manuscript
# -----------------------------------------------------------------
print("[01_manuscript]")
_full_src = OUT / "Manuscript_full_prose_v2.docx" if (OUT / "Manuscript_full_prose_v2.docx").exists() else OUT / "Manuscript_full_prose.docx"
copy_file(_full_src, SUB / "01_manuscript" / "Manuscript_full_prose.docx")
copy_file(OUT / "Manuscript_rough_draft.docx", SUB / "01_manuscript" / "Manuscript_rough_draft.docx")
copy_file(BASE / "manuscript_full_prose.md", SUB / "01_manuscript" / "manuscript_full_prose.md")
copy_file(BASE / "manuscript_rough_draft.md", SUB / "01_manuscript" / "manuscript_rough_draft.md")

# -----------------------------------------------------------------
# 02_table
# -----------------------------------------------------------------
print("\n[02_table]")
copy_file(OUT / "Table1_n349.docx", SUB / "02_table" / "Table1_n349.docx")
copy_file(OUT / "Table1_n349.tsv",  SUB / "02_table" / "Table1_n349.tsv")

# -----------------------------------------------------------------
# 03_figures — copy all 300dpi PNGs
# -----------------------------------------------------------------
print("\n[03_figures]")
for png in sorted((OUT / "figs_300dpi").glob("*.png")):
    copy_file(png, SUB / "03_figures" / png.name)
for svg in sorted((OUT / "figs_300dpi").glob("*.svg")):
    copy_file(svg, SUB / "03_figures" / svg.name)

# -----------------------------------------------------------------
# 04_supplements — captions + references + analysis_plan + reviewer appeal points
# -----------------------------------------------------------------
print("\n[04_supplements]")
md_to_docx(BASE / "figure_captions.md",        SUB / "04_supplements" / "Figure_captions.docx")
md_to_docx(BASE / "references.md",             SUB / "04_supplements" / "References.docx")
md_to_docx(BASE / "analysis_plan.md",          SUB / "04_supplements" / "Analysis_plan_with_results.docx")
md_to_docx(BASE / "reviewer_appeal_points.md", SUB / "04_supplements" / "Reviewer_appeal_points.docx")
copy_file(BASE / "figure_captions.md",         SUB / "04_supplements" / "figure_captions.md")
copy_file(BASE / "references.md",              SUB / "04_supplements" / "references.md")
copy_file(BASE / "analysis_plan.md",           SUB / "04_supplements" / "analysis_plan.md")
copy_file(BASE / "reviewer_appeal_points.md",  SUB / "04_supplements" / "reviewer_appeal_points.md")

# -----------------------------------------------------------------
# 05_correspondence — PI email
# -----------------------------------------------------------------
print("\n[05_correspondence]")
md_to_docx(BASE / "email_to_PI.md", SUB / "05_correspondence" / "Email_to_PI_SuhmiChung.docx")
copy_file(BASE / "email_to_PI.md",  SUB / "05_correspondence" / "email_to_PI.md")

# -----------------------------------------------------------------
# 06_supplementary_data — key TSV deliverables (referenced as Suppl Tables S1–S6)
# -----------------------------------------------------------------
print("\n[06_supplementary_data]")
supp_files = [
    ("step6_KW_results.tsv",                          "SupplTableS1_KW_per_feature.tsv"),
    ("step6_Dunn_posthoc.tsv",                        "SupplTableS2_Dunn_posthoc.tsv"),
    ("step8_DEG_Lymphocyte_inflamed.tsv",             "SupplTableS3a_DEG_Lymphocyte_inflamed.tsv"),
    ("step8_DEG_Myeloid_dominant.tsv",                "SupplTableS3b_DEG_Myeloid_dominant.tsv"),
    ("step8_DEG_Immune_desert.tsv",                   "SupplTableS3c_DEG_Immune_desert.tsv"),
    ("step8_gprofiler_Lymphocyte_inflamed.tsv",       "SupplTableS4a_gProfiler_Lymphocyte_inflamed.tsv"),
    ("step8_gprofiler_Myeloid_dominant.tsv",          "SupplTableS4b_gProfiler_Myeloid_dominant.tsv"),
    ("step8_gprofiler_Immune_desert.tsv",             "SupplTableS4c_gProfiler_Immune_desert.tsv"),
    ("step9_cox_multivariable.tsv",                   "SupplTableS5_Cox_multivariable.tsv"),
    ("step10_BRAF_ALT_ecotype_assignment.tsv",        "SupplTableS6_BRAF_ALT_ecotype_assignment.tsv"),
    ("step4_LM22_vs_methods_spearman.tsv",            "SupplTableS7_LM22_vs_methods_spearman.tsv"),
    ("step4_microglia_ssGSEA_vs_LM22_macrophage.tsv", "SupplTableS8_Microglia_vs_LM22_macrophage.tsv"),
    ("ecotype_LM22_main_k3_annotated.tsv",            "SupplTableS9_Main_cohort_ecotype_k3.tsv"),
    ("consensus_LM22_main_B1000_PAC.tsv",             "SupplTableS10_PAC_B1000.tsv"),
    ("consensus_LM22_main_B1000_silhouette.tsv",      "SupplTableS11_Silhouette_B1000.tsv"),
    ("consensus_LM22_main_B1000_within_cluster_robustness.tsv", "SupplTableS12_Within_cluster_consensus_B1000.tsv"),
    ("ecotype_LM22_main_assoc_pvalues.tsv",           "SupplTableS13_Ecotype_associations.tsv"),
    ("step7_theme_KW.tsv",                            "SupplTableS14_Theme_composite_KW.tsv"),
    ("step9_KM_logrank_overall.tsv",                  "SupplTableS15_KM_logrank_overall.tsv"),
    ("step9_KM_logrank_by_cohort.tsv",                "SupplTableS16_KM_logrank_by_cohort.tsv"),
    ("step9_median_OS_by_ecotype.tsv",                "SupplTableS17_Median_OS_by_ecotype.tsv"),
    ("supplS3_highconf_vs_Main_k3_crosstab.tsv",      "SupplTableS18_Highconf_subset_crosstab.tsv"),
    ("supplS4_fullpool_vs_Main_k3_crosstab.tsv",      "SupplTableS19_Fullpool_vs_Main_crosstab.tsv"),
    ("supplS4_fullpool_vs_BRAF_projection_crosstab.tsv", "SupplTableS20_Fullpool_vs_BRAF_crosstab.tsv"),
]
for src_name, dst_name in supp_files:
    src = OUT / src_name
    if src.exists():
        copy_file(src, SUB / "06_supplementary_data" / dst_name)
    else:
        print(f"  WARN: {src_name} not found, skipped")

# README
readme_text = """\
# Submission package — Decoding the Tumor Immune Microenvironment of Pediatric Gliomas

Lead author: Suhmi Chung, MD (Asan Medical Center, Department of Neurosurgery)
Dataset: OpenPedCan v15
Date assembled: 2026-06-12
Target journal: Cancers (MDPI), dry-lab paper

## Folder structure

01_manuscript/
  - Manuscript_full_prose.docx   ← MAIN manuscript (Abstract + Intro + Methods + Results + Discussion + Limitations + Conclusions)
  - Manuscript_rough_draft.docx  ← earlier outline version (figure placements + decisive sentences only)
  - .md source files

02_table/
  - Table1_n349.docx — Clinical and molecular characteristics of the n=349 Main cohort
  - Table1_n349.tsv

03_figures/  (all PNG ≥ 300 DPI, manuscript-ready)
  - Fig1_workflow.png + .svg   — Study workflow
  - Fig2A/B/C  — LM22 QC, cross-method consensus, ssGSEA vs LM22 macrophage
  - Fig3A/B/C/D — PAC/silhouette, PCA+UMAP, feature heatmap, cohort stack
  - Fig4A/B/C/D/E — Themes, signature heatmap, cohort-stratified, focus boxplots, Dunn heatmap
  - Fig5A/B — Volcanos, top-15 enrichment
  - Fig6A/B/C — KM overall, KM per cohort grid, Cox forest
  - Fig7A/B/C — BRAF projection, MAPK by family, PCA overlay
  - S3_highconf_subset_reclustering.png — Suppl Fig S3
  - S4_fullpool_sensitivity.png + S4_fullpool_PCA.png — Suppl Fig S4

04_supplements/
  - Figure_captions.docx — Cancers-style full captions for all figures incl. supplementary
  - References.docx — 47-entry reference list (Cancers MDPI numbered style); verify DOIs before submission
  - Analysis_plan_with_results.docx — Step 1–10 산출 수치 반영본
  - Reviewer_appeal_points.docx — 5 killer sentences + preemptive limitations

05_correspondence/
  - Email_to_PI_SuhmiChung.docx — Korean-language status report to Suhmi Chung, MD

06_supplementary_data/ (TSV — Suppl Tables S1–S20)
  - S1 KW per feature, S2 Dunn post-hoc
  - S3a–c Per-ecotype DEG (Wilcoxon, BH-FDR)
  - S4a–c g:Profiler enrichment
  - S5 Multivariable Cox
  - S6 BRAF_ALT ecotype assignment
  - S7 LM22 vs orthogonal methods Spearman
  - S8 Microglia signatures vs LM22 macrophage axes
  - S9 Main cohort k=3 ecotype call
  - S10–S12 Bootstrap robustness (B=1000)
  - S13 Ecotype × clinical association p
  - S14 Theme composite KW
  - S15–S17 Survival summary
  - S18–S20 Sensitivity crosstabs

## Pre-submission checklist (recommended)
[ ] Verify DOIs and PMIDs in References.docx (items marked [verify] need PubMed re-check)
[ ] Confirm OpenPedCan v15 release tag and add appropriate Acknowledgement of Children's Brain Tumor Network
[ ] Decide narrative on G34 reversal (strong vs nuanced) — see PI email §6.1
[ ] Confirm co-author list, funding sources, IRB exemption statement
[ ] Re-render Fig 1 in BioRender if lab preference (matplotlib SVG provided as drop-in)
[ ] Cancers MDPI Word template — paste body into template before submission

## Reproducibility
All figures and supplementary tables can be regenerated end-to-end from:
  scripts/step{1..10}_*.py
  scripts/suppl_S{3,4}_*.py
  scripts/regen_all_figures_300dpi.py
"""
(SUB / "README.txt").write_text(readme_text)
print(f"\n[README] wrote {SUB / 'README.txt'}")

# Final inventory
print("\n=== INVENTORY ===")
for f in sorted(SUB.rglob("*")):
    if f.is_file():
        print(f"  {f.relative_to(BASE)}")
print(f"\nFiles total: {sum(1 for f in SUB.rglob('*') if f.is_file())}")
