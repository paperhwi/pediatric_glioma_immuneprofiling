"""Build the renumbered supplementary table workbook for the Cancers submission.

Only tables needed to reproduce a figure or a reported analysis are kept.
Audit copies, superseded pre-corrected versions and tables that duplicate
Table 1 or a main figure are dropped; the crosswalk sheet records where
every old item went.
"""
import os
import numpy as np
import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

U = "/mnt/user-data/uploads/Open PBTA/Revision"
ALL = f"{U}/FINAL MANUSCRIPT 260722 - 수정본/3. Tables/Source workbooks/All_tables_and_annotation_final.xlsx"
SUP = f"{U}/FINAL MANUSCRIPT 260722 - 수정본/3. Tables/Supplementary/Supplementary_tables.xlsx"
REV = f"{U}/REVISION PACKAGE 260916/3. Tables/Supplementary_tables_S12-S24_revision.xlsx"
TSV = f"{U}/REVISION PACKAGE 260916/3. Tables"
DEG = f"{U}/Week1/_inputs/adjusted_DEG_all_contrasts.tsv"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f"{ROOT}/tables/Supplementary_Tables_Cancers.xlsx"
os.makedirs(f"{ROOT}/tables", exist_ok=True)

GITHUB = "https://github.com/paperhwi/pediatric_glioma_immuneprofiling"


def rev(sheet):
    """Read a sheet of the revision workbook, skipping its title block."""
    raw = pd.read_excel(REV, sheet_name=sheet, header=None)
    blanks = [i for i in range(min(10, len(raw))) if raw.iloc[i].isna().all()]
    hdr = max(blanks) + 1 if blanks else 0
    df = pd.read_excel(REV, sheet_name=sheet, header=hdr)
    return df.dropna(how="all").reset_index(drop=True)


def block(title, df):
    return {"title": title, "df": df}


# ---------------------------------------------------------------- loaders
def t_signature_gene_lists():
    d = pd.read_excel(SUP, sheet_name="S24_signature_gene_lists")
    d = d[["signature", "source_annotation", "n_genes", "genes"]]
    return [block(None, d)]


def t_feature_lists():
    d = pd.read_excel(SUP, sheet_name="item6_feature_lists")
    d = d.rename(columns={"in_46_feature_MAIN": "in_46_feature_primary_matrix",
                          "in_29_feature_TPM": "in_34_feature_TPM_matrix"})
    return [block(None, d)]


def t_location_mapping():
    return [block(None, pd.read_excel(SUP, sheet_name="loc1_classification_mapping"))]


def t_analysis_sets():
    return [block(None, rev("S17"))]


def t_ecotype_assignment():
    return [block(None, pd.read_excel(ALL, sheet_name="S9_Main_cohort_ecotype_k3"))]


def t_ksweep():
    return [block(None, rev("S12"))]


def t_gap():
    return [block(None, rev("S13"))]


def t_criteria():
    return [block(None, rev("S14"))]


def t_feature_sensitivity():
    return [block(None, pd.read_excel(SUP, sheet_name="item5_FINAL_46feature_sensitivi"))]


def t_permanova_blocks():
    return [block(None, pd.read_excel(SUP, sheet_name="item6_FINAL_permanova"))]


def t_crossmethod():
    return [block(None, pd.read_excel(SUP, sheet_name="item4_Figure2B_n349_vs_n702"))]


def t_microglia_lm22():
    return [block(None, pd.read_excel(ALL, sheet_name="S8_Microglia_vs_LM22_macrophage"))]


def t_location_stats():
    return [block("Contingency analyses of ecotype by anatomical location",
                  pd.read_excel(SUP, sheet_name="loc3_restricted_analyses")),
            block("Full statistic set for the ecotype x molecular group and "
                  "ecotype x location analyses",
                  pd.read_excel(SUP, sheet_name="item7_location_statistics")),
            block("Midline-versus-hemispheric odds ratios, one versus rest",
                  pd.read_excel(SUP, sheet_name="item7_midline_vs_hemispheric_OR"))]


def t_kw():
    return [block("Per-signature and per-cell-type Kruskal-Wallis tests",
                  pd.read_excel(ALL, sheet_name="S1_KW_per_feature")),
            block("Theme-composite Kruskal-Wallis tests",
                  pd.read_excel(ALL, sheet_name="S14_Theme_composite_KW"))]


def t_dunn():
    return [block(None, pd.read_excel(ALL, sheet_name="S2_Dunn_posthoc"))]


def t_deg():
    d = pd.read_csv(DEG, sep="\t")
    d = d[["contrast", "gene", "adjusted_log2FC", "t", "p", "q_BH", "df_resid"]]
    d = d.sort_values(["contrast", "q_BH"]).reset_index(drop=True)
    return [block(None, d)]


def t_overlap():
    return [block(None, rev("S15"))]


def t_shared_genes():
    return [block(None, rev("S16"))]


def t_gprofiler():
    parts = []
    for sh in ["S4a_gProfiler_Lymphocyte_inflam", "S4b_gProfiler_Myeloid_dominant",
               "S4c_gProfiler_Immune_desert"]:
        try:
            parts.append(pd.read_excel(ALL, sheet_name=sh))
        except Exception:
            pass
    d = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return [block(None, d)]


def t_cox():
    return [block("Cox proportional-hazards models",
                  pd.read_excel(SUP, sheet_name="item8_cox_models")),
            block("Grambsch-Therneau tests on scaled Schoenfeld residuals",
                  pd.read_excel(SUP, sheet_name="item8_schoenfeld")),
            block("Ecotype hazard ratios with and without location adjustment",
                  pd.read_excel(SUP, sheet_name="item8_model_comparison"))]


def t_missingness():
    return [block("Survival evaluability by baseline characteristic", rev("S18")),
            block("Ecotype and survival evaluability within each molecular group",
                  rev("S19"))]


def t_sc_attribution():
    a = pd.read_csv(f"{TSV}/T_S25_signature_attribution.tsv", sep="\t")
    b = pd.read_csv(f"{TSV}/T_S29_theme_audit.tsv", sep="\t")
    return [block("Signature attribution and AUROC for the named cell type", a),
            block("Effect of correcting the two theme assignments", b)]


def t_sc_gene_source():
    return [block(None, pd.read_csv(f"{TSV}/T_S26_ecotype_gene_source.tsv", sep="\t"))]


def t_spatial():
    a = pd.read_csv(f"{TSV}/T_S27_spatial_autocorrelation.tsv", sep="\t")
    b = pd.read_csv(f"{TSV}/T_S28_spot_ecotype.tsv", sep="\t")
    b = b.rename(columns={"Lymp": "n_Lymphocyte_inflamed", "Myel": "n_Myeloid_dominant",
                          "Immu": "n_Immune_desert",
                          "same_neighbour_fraction": "same_neighbor_fraction",
                          "p": "myeloid_lymphoid_p"})
    return [block("Moran's I for each signature in each section", a),
            block("Spot-level ecotype assignment and spatial coherence", b)]


def t_software():
    a = rev("S24")
    b = pd.read_csv(f"{TSV}/T_S30_week2_methods.tsv", sep="\t")
    return [block("Software versions, parameters and random seeds", a),
            block("Single-cell and spatial analysis parameters", b)]


# ---------------------------------------------------------------- table set
TABLES = [
    ("S1", "Immune signature gene lists",
     "Complete membership of every gene signature used in the study, with its source annotation.",
     "Methods; Figures 2, 3, 5", t_signature_gene_lists),
    ("S2", "Feature composition of the clustering matrices",
     "Every feature of the 46-feature primary clustering matrix (22 CIBERSORTx LM22 fractions plus "
     "24 brain-tuned ssGSEA signatures) and of the 34-feature TPM-harmonized sensitivity matrix "
     "(10 quanTIseq axes plus the 24 ssGSEA signatures), with block membership.",
     "Methods; Figure 3; Figure S5", t_feature_lists),
    ("S3", "Anatomical-location classification mapping",
     "Mapping of the source CNS region and primary site annotations onto the four location classes.",
     "Methods; Figure 4", t_location_mapping),
    ("S4", "Analysis-set reconciliation",
     "Reconciliation of the main cohort (n = 349), the location-annotated set (n = 332), the "
     "sequential PERMANOVA set (n = 258) and the survival-evaluable set (n = 251).",
     "Figure 1B", t_analysis_sets),
    ("S5", "Immune-ecotype assignment of the main cohort",
     "Locked exploratory k = 3 assignment for all 349 biospecimens.",
     "Figures 3, 4, 6, 7", t_ecotype_assignment),
    ("S6", "Cluster-number sweep, k = 2 to 10",
     "Consensus clustering of the 46-feature matrix, 1,000 bootstrap replicates, 80% item "
     "resampling, seed 42.",
     "Figure 3A-C, 3E", t_ksweep),
    ("S7", "Gap statistic, k = 1 to 10",
     "Gap statistic against 50 uniform reference datasets, with the one-standard-error criterion.",
     "Figure 3D", t_gap),
    ("S8", "Cluster-number criteria and the k each selects",
     "Each cluster-number criterion, the k it selects and the criterion family it belongs to.",
     "Figure 3A-E", t_criteria),
    ("S9", "Feature-definition sensitivity of the k = 3 assignment",
     "De novo k = 3 clustering on the primary and on the TPM-harmonized matrices compared with the "
     "locked assignment.",
     "Figure S5", t_feature_sensitivity),
    ("S10", "PERMANOVA on circular and held-out feature blocks",
     "One-way PERMANOVA of ecotype and molecular group against the feature blocks used to define "
     "the clusters and against blocks that were not.",
     "Figure S6", t_permanova_blocks),
    ("S11", "Cross-method concordance, main cohort and pooled set",
     "Spearman correlations between CIBERSORTx LM22 and the peer methods in the n = 349 main cohort "
     "and in the same-resource n = 702 pool.",
     "Figure 2B; Figures S1, S2", t_crossmethod),
    ("S12", "Brain-tuned myeloid signatures versus LM22 macrophage and monocyte axes",
     "Spearman correlations between the brain-tuned microglia and monocyte-derived-macrophage "
     "signatures and the LM22 macrophage and monocyte axes.",
     "Figure S1", t_microglia_lm22),
    ("S13", "Location and molecular-group contingency analyses",
     "Contingency statistics for ecotype by anatomical location and by integrated molecular group, "
     "with Monte-Carlo P values and midline-versus-hemispheric odds ratios.",
     "Figure 4", t_location_stats),
    ("S14", "Kruskal-Wallis tests of the immune features and theme composites",
     "Per-feature and per-theme Kruskal-Wallis tests across the three ecotypes with "
     "Benjamini-Hochberg adjustment.",
     "Figures S7, S8", t_kw),
    ("S15", "Dunn post-hoc comparisons",
     "Pairwise Dunn comparisons for every feature, with within-feature Benjamini-Hochberg "
     "adjustment.",
     "Figure S8", t_dunn),
    ("S16", "Molecular-group-adjusted differential expression, all three contrasts",
     "Gene-wise ordinary least squares on log2(TPM + 1) with ecotype and integrated molecular group, "
     "18,586 annotated genes, n = 349.",
     "Figure 6A-C", t_deg),
    ("S17", "Overlap of the three adjusted contrasts",
     "Partition of the up-regulated gene sets of the two Immune-desert-referenced contrasts.",
     "Figure 6D", t_overlap),
    ("S18", "Genes shared by both Immune-desert-referenced contrasts",
     "The 426 genes up-regulated in both Lymphocyte-inflamed and Myeloid-dominant versus "
     "Immune-desert, with both adjusted effect sizes.",
     "Figure 6D-F", t_shared_genes),
    ("S19", "Pathway enrichment from the adjusted gene lists",
     "g:Profiler enrichment of the molecular-group-adjusted up-regulated gene lists.",
     "Figure S9", t_gprofiler),
    ("S20", "Cox models and proportional-hazards diagnostics",
     "Multivariable Cox models, Grambsch-Therneau tests on scaled Schoenfeld residuals, and the "
     "comparison with and without location adjustment.",
     "Figure 7B; Figures S11, S13", t_cox),
    ("S21", "Survival evaluability and missingness",
     "Survival evaluability by baseline characteristic and ecotype by evaluability within each "
     "molecular group.",
     "Figure S14", t_missingness),
    ("S22", "Single-cell signature attribution",
     "AUROC of each directional signature for the cell type it is named for, and the effect of "
     "correcting the two misassigned theme memberships.",
     "Figure 5A; Figure S15", t_sc_attribution),
    ("S23", "Cell type expressing the ecotype-axis genes",
     "Mean expression and detection rate of each detectable ecotype-axis gene in matched malignant, "
     "myeloid and T cells.",
     "Figure 5B; Figure S15", t_sc_gene_source),
    ("S24", "Spatial autocorrelation and spot-level ecotype assignment",
     "Moran's I for every signature in each section, and the spot-level assignment and spatial "
     "coherence statistics.",
     "Figure 5C-E; Figure S16", t_spatial),
    ("S25", "Software versions, parameters and random seeds",
     "Execution environment, package versions, parameters and seeds for every analysis.",
     "Methods", t_software),
]

# ---------------------------------------------------------------- crosswalk
CROSSWALK = [
    ("Main figure", "Figure 1", "Figure 1A", "workflow redrawn; k = 2-10 evaluation and bootstrap consensus clustering"),
    ("Main figure", "Figure 1B (revision)", "Figure 1B", "redrawn to match Figure 1A"),
    ("Main figure", "Figure 2A", "Figure 2A", "unchanged"),
    ("Main figure", "Figure 2B", "Figure 2B", "unchanged"),
    ("Main figure", "Figure 2C", "Figure S1", "moved to supplementary"),
    ("Main figure", "Figure 2D", "Figure 2C", "panel letter changed"),
    ("Main figure", "Figure 3A (revision)", "Figure 3A-E", "five separate panels; dAUC annotation reworded"),
    ("Main figure", "Figure 3B", "Figure 3F", "unchanged"),
    ("Main figure", "Figure 3C", "Figure 3G", "unchanged"),
    ("Main figure", "Figure 3D", "removed", "duplicates Figure 4B"),
    ("Main figure", "Figure 4", "Figure 4", "unchanged"),
    ("Main figure", "Figure 5A", "Figure S7A", "moved to supplementary"),
    ("Main figure", "Figure 5B", "Figure S8A", "moved to supplementary"),
    ("Main figure", "Figure 5C", "Figure S7B", "moved to supplementary"),
    ("Main figure", "Figure 5D", "Figure S8B", "moved to supplementary"),
    ("Main figure", "Figure 5E", "Figure S8C", "moved to supplementary"),
    ("Main figure", "(new)", "Figure 5A-E", "single-cell and spatial results promoted to a main figure"),
    ("Main figure", "Figure 6A (revision)", "Figure 6A-C", "regenerated; reference group named in each title"),
    ("Main figure", "Figure 6B", "Figure S9", "moved to supplementary"),
    ("Main figure", "Figure 6C (revision)", "Figure 6D-F", "Venn replaced by proportional overlap bars"),
    ("Main figure", "Figure 7A", "Figure 7A", "unchanged; numbers at risk and censoring marks retained"),
    ("Main figure", "Figure 7B", "Figure S11A", "moved to supplementary"),
    ("Main figure", "Figure 7C", "Figure 7B", "panel letter changed"),
    ("Suppl. figure", "Figure S1", "GitHub repository", "pooled n = 702 cross-method sensitivity"),
    ("Suppl. figure", "Figure S2", "GitHub repository", "cross-method consensus heatmap"),
    ("Suppl. figure", "Figure S3", "Figure S3", "unchanged"),
    ("Suppl. figure", "Figure S4", "Figure S2", "renumbered"),
    ("Suppl. figure", "Figure S5", "Figure S6", "renumbered"),
    ("Suppl. figure", "Figure S6", "Figure S4", "renumbered"),
    ("Suppl. figure", "Figure S7A-C", "GitHub repository", "pooled n = 702 re-clustering sensitivity"),
    ("Suppl. figure", "Figure S8", "Figure S5", "renumbered"),
    ("Suppl. figure", "Figure S9", "Figure S10", "renumbered"),
    ("Suppl. figure", "Figure S10A", "removed", "duplicated by Figure S11 and by the removed small-stratum panels"),
    ("Suppl. figure", "Figure S10B", "Figure S11B", "renumbered"),
    ("Suppl. figure", "Figure S11", "Figure S12", "renumbered"),
    ("Suppl. figure", "Figure S12", "Figure S13", "renumbered"),
    ("Suppl. figure", "Figure S13A-D", "Figure S14A-D", "renumbered"),
    ("Suppl. figure", "Figure S14A-C", "GitHub repository", "BRAF/RTK-altered projection"),
    ("Suppl. figure", "Figure S15", "GitHub repository", "immune-checkpoint transcripts"),
    ("Suppl. figure", "Figure S16", "GitHub repository", "CAR-T target transcripts"),
    ("Suppl. figure", "Figure S17", "GitHub repository", "audit copy of the adjusted differential expression"),
    ("Suppl. figure", "Figure S18", "GitHub repository", "very small DHG and IHG survival strata"),
    ("Suppl. figure", "Figure S19 (revision)", "GitHub repository", "hematopoietic stem and progenitor programs"),
    ("Suppl. figure", "Figure S20 (revision)", "Figure S15", "single-cell analysis, full"),
    ("Suppl. figure", "Figure S21 (revision)", "Figure S16", "spatial analysis, full"),
    ("Suppl. table", "Table S1", "removed", "auxiliary molecular alteration flags; descriptive only"),
    ("Suppl. table", "Table S2A", "Table S1", "renumbered"),
    ("Suppl. table", "Table S2B", "merged into Table S2", "duplicated the feature-list table"),
    ("Suppl. table", "Table S3", "Table S3", "renumbered"),
    ("Suppl. table", "Table S4", "merged into Table S13", "location contingency analyses"),
    ("Suppl. table", "Table S5", "Table S11", "renumbered"),
    ("Suppl. table", "Table S6", "Table S9", "renumbered"),
    ("Suppl. table", "Table S7", "Table S10", "renumbered"),
    ("Suppl. table", "Table S8", "Table S2", "renumbered"),
    ("Suppl. table", "Table S9", "Table S20", "renumbered"),
    ("Suppl. table", "Table S10", "merged into Table S13", "location statistics and odds ratios"),
    ("Suppl. table", "Table S11", "removed", "duplicates Table 1 and Figure 4"),
    ("Suppl. table", "Table S12 (revision)", "Table S6", "renumbered"),
    ("Suppl. table", "Table S13 (revision)", "Table S7", "renumbered"),
    ("Suppl. table", "Table S14 (revision)", "Table S8", "renumbered"),
    ("Suppl. table", "Table S15 (revision)", "Table S17", "renumbered"),
    ("Suppl. table", "Table S16 (revision)", "Table S18", "renumbered"),
    ("Suppl. table", "Table S17 (revision)", "Table S4", "renumbered"),
    ("Suppl. table", "Table S18 (revision)", "merged into Table S21", "survival evaluability"),
    ("Suppl. table", "Table S19 (revision)", "merged into Table S21", "survival evaluability"),
    ("Suppl. table", "Tables S20-S23 (revision)", "GitHub repository", "hematopoietic progenitor programs"),
    ("Suppl. table", "Table S24 (revision)", "Table S25", "renumbered"),
    ("Suppl. table", "Table S25 (revision)", "Table S22", "renumbered"),
    ("Suppl. table", "Table S26 (revision)", "Table S23", "renumbered"),
    ("Suppl. table", "Table S27 (revision)", "merged into Table S24", "spatial autocorrelation"),
    ("Suppl. table", "Table S28 (revision)", "merged into Table S24", "spot-level ecotype assignment"),
    ("Suppl. table", "Table S29 (revision)", "merged into Table S22", "theme reassignment sensitivity"),
    ("Suppl. table", "Table S30 (revision)", "merged into Table S25", "single-cell and spatial parameters"),
    ("Suppl. data", "Supplementary data S1, S2, S14", "Tables S14, S15", "Kruskal-Wallis and Dunn results"),
    ("Suppl. data", "Supplementary data S3a-c", "Table S16", "replaced by the molecular-group-adjusted models"),
    ("Suppl. data", "Supplementary data S4a-c", "Table S19", "g:Profiler enrichment"),
    ("Suppl. data", "Supplementary data S9", "Table S5", "ecotype assignment"),
    ("Suppl. data", "Supplementary data S6, S18-S22", "GitHub repository", "BRAF projection and pooled-set crosstabs"),
]

TITLE_F = Font(bold=True, size=12)
DESC_F = Font(italic=True, size=10, color="404040")
HEAD_F = Font(bold=True, size=10)


def write_sheet(wb, num, title, desc, supports, blocks):
    ws = wb.create_sheet(num)
    ws["A1"] = f"Supplementary Table {num}. {title}"
    ws["A1"].font = TITLE_F
    ws["A2"] = desc
    ws["A2"].font = DESC_F
    ws["A3"] = f"Supports: {supports}"
    ws["A3"].font = DESC_F
    r = 5
    for b in blocks:
        if b["title"]:
            ws.cell(r, 1, b["title"]).font = HEAD_F
            r += 1
        df = b["df"]
        for j, c in enumerate(df.columns, start=1):
            ws.cell(r, j, str(c)).font = HEAD_F
        r += 1
        for row in df.itertuples(index=False):
            for j, v in enumerate(row, start=1):
                if isinstance(v, (np.integer,)):
                    v = int(v)
                elif isinstance(v, (np.floating,)):
                    v = None if np.isnan(v) else float(v)
                elif isinstance(v, np.bool_):
                    v = bool(v)
                elif v is not None and not isinstance(v, (int, float, str, bool)):
                    v = str(v)
                ws.cell(r, j, v)
            r += 1
        r += 2
    widths = {}
    for b in blocks:
        for j, c in enumerate(b["df"].columns, start=1):
            widths[j] = min(46, max(widths.get(j, 10), len(str(c)) + 2))
    for j, w in widths.items():
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = "A6"
    return sum(len(b["df"]) for b in blocks)


def main():
    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    readme = wb.create_sheet("README")
    readme["A1"] = "Supplementary Tables — Cancers submission"
    readme["A1"].font = Font(bold=True, size=13)
    readme["A2"] = ("Integrated transcriptomic immune profiling identifies survival-associated "
                    "immune ecotypes in pediatric diffuse high-grade glioma")
    readme["A2"].font = DESC_F
    readme["A3"] = (f"Tables removed from this submission are archived in the study repository: {GITHUB}")
    readme["A3"].font = DESC_F
    for j, h in enumerate(["Table", "Title", "Supports", "Rows", "Description"], start=1):
        readme.cell(5, j, h).font = HEAD_F

    r = 6
    for num, title, desc, supports, loader in TABLES:
        blocks = loader()
        n = write_sheet(wb, num, title, desc, supports, blocks)
        readme.cell(r, 1, num)
        readme.cell(r, 2, title)
        readme.cell(r, 3, supports)
        readme.cell(r, 4, n)
        readme.cell(r, 5, desc)
        r += 1
        print(f"  {num:<5} {title[:58]:<60} {n:>6} rows")
    for col, w in zip("ABCDE", [9, 58, 30, 9, 100]):
        readme.column_dimensions[col].width = w
    for row in readme.iter_rows(min_row=6, max_row=r, min_col=5, max_col=5):
        for c in row:
            c.alignment = Alignment(wrap_text=False)

    cw = wb.create_sheet("Crosswalk")
    cw["A1"] = "Old to new numbering, and items archived to the repository"
    cw["A1"].font = Font(bold=True, size=12)
    cw["A2"] = (f"Items marked 'GitHub repository' are deposited at {GITHUB} and are not part of "
                "the submitted supplementary material.")
    cw["A2"].font = DESC_F
    for j, h in enumerate(["Category", "Previous item", "In this submission", "Note"], start=1):
        cw.cell(4, j, h).font = HEAD_F
    for i, row in enumerate(CROSSWALK, start=5):
        for j, v in enumerate(row, start=1):
            cw.cell(i, j, v)
    for col, w in zip("ABCD", [14, 30, 24, 76]):
        cw.column_dimensions[col].width = w
    cw.freeze_panes = "A5"

    wb._sheets = [wb["README"], wb["Crosswalk"]] + [wb[n] for n, *_ in TABLES]
    wb.save(OUT)
    print(f"\nwrote {OUT}  ({os.path.getsize(OUT)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
