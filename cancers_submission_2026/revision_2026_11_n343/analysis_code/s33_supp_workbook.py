"""Supplementary Tables workbook (n = 343), gene-set source list and Table S25 (software, parameters, seeds, inputs)."""
import sys, glob, hashlib, json, importlib.metadata as md, platform; sys.path.insert(0,"/home/claude/q343/src"); from common import *
O=ROOT+"out/"; D=ROOT+"deliv/"
def rd(p,**k): return pd.read_csv(O+p,sep="\t",**k)
def md5(p):
    h=hashlib.md5()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
# ---- gene sets and sources
B=pd.read_excel("/mnt/user-data/uploads/Pediatric GBM/Revision_round1_Q1-Q5/Item_B_adjusted_Cox/B_gene_sets_and_sources.xlsx",sheet_name=None)
from published_sets import SETS, CHECKPOINT_GENES
pub=pd.DataFrame([dict(programme=k,n_genes=len(v["genes"]),genes=",".join(v["genes"]),source=v["ref"]) for k,v in SETS.items()]+[dict(programme="Checkpoint and evasion genes (gene-level)",n_genes=len(CHECKPOINT_GENES),genes=",".join(CHECKPOINT_GENES),source="curated list of ligands/receptors examined individually")])
SC=pd.read_excel("/mnt/user-data/uploads/Pediatric GBM/Revision_round1_Q1-Q5/Item_C_singlecell/C_results_singlecell.xlsx","Gene_programs")
SI=pd.DataFrame([dict(item="Immune-presence axis",definition="PC1 of the 21 immune ssGSEA signatures (Brain_tuned_ssGSEA sheet, excluding Cell_Cycle_Proliferation, Stemness_Brain_Tumor, MAPK_Activity), z-scored within n = 343",source="this study"),
                 dict(item="ScPCA reference cell classes",definition="lineage-marker scores on ScPCA log-normalised counts (Myeloid C1QA/B/C, TYROBP, AIF1, CSF1R, CD14, FCER1G, CD68, LAPTM5; T CD3D/E, CD2, CD247, TRAC, IL7R; NK NKG7, GNLY, KLRD1, KLRF1, PRF1; B MS4A1, CD79A/B, IGHM, JCHAIN, MZB1, IGKC; DC CLEC10A, CD1C, FCER1A, CLEC9A, LILRA4, IRF8; endothelial PECAM1, CLDN5, VWF, FLT1, CDH5, ESAM; mural PDGFRB, RGS5, ACTA2, NOTCH3, MYH11, KCNJ8; oligodendrocyte MOG, MOBP, MAG, PLP1, MBP, ERMN, OPALIN); microglia (P2RY12, TMEM119, CX3CR1, SALL1, GPR34, P2RY13) vs MDM (CD163, F13A1, MRC1, S100A8, S100A9, VCAN, FCN1, LYZ, TGFBI, CCR2, EREG)",source="ScPCA SCPCP000001 (DeSisto et al. 2024; Hawkins et al. 2024 portal)"),
                 dict(item="Visium spot classes",definition="Pearson correlation of each spot (z of log1p CP10k within section) with the bulk immune-state centroids on the shared Figure 6D genes; nearest centroid",source="GSE268577"),
                 dict(item="Pan-cancer frameworks",definition="Thorsson immune subtypes (ImmuneSubtypeClassifier); Bagaev MFP (29 signatures)",source="Thorsson et al. 2018 Immunity; Bagaev et al. 2021 Cancer Cell")])
flags=pd.DataFrame([dict(item="ssGSEA_T_Cell_Exhaustion",flag="gene list contains TCF7 and CD274 (not exhaustion-specific); interpret as an inhibitory-receptor/checkpoint score"),
                    dict(item="ssGSEA_Cytotoxic_Lymphocyte_Curated",flag="citation (Wang 2017) to be confirmed by the authors; the signature is discordant in external cells (AUROC 0.075, identifies T cells)"),
                    dict(item="Nabbi 2023",flag="HGG/DIPG representation in the 925 pedNST to be confirmed from the full text before citing specifics")])
with pd.ExcelWriter(D+"Gene_sets_and_sources_n343.xlsx") as w:
    pd.DataFrame({"sheet":["Brain_tuned_ssGSEA_24","Curated_programmes_17","Published_programmes","LM22_22","ESTIMATE_immune","SingleCell_programmes","Other_definitions","Citations_to_confirm"],
                  "content":["24 ssGSEA signatures (21 immune + 3 tumour-intrinsic) used for clustering and the immune axis","17 curated programmes (mean z of log2 TPM+1)","Published programmes and 21 checkpoint genes","CIBERSORTx LM22 cell types","ESTIMATE immune signature (validation)","Single-cell programme gene lists (GSE227983, ScPCA)","Axis, reference classes, spots, frameworks","Items needing author confirmation"]}).to_excel(w,sheet_name="README",index=False)
    B["Brain_tuned_ssGSEA_24"].to_excel(w,sheet_name="Brain_tuned_ssGSEA_24",index=False); B["Curated_programs_17"].to_excel(w,sheet_name="Curated_programmes_17",index=False)
    pub.to_excel(w,sheet_name="Published_programmes",index=False); B["LM22_22"].to_excel(w,sheet_name="LM22_22",index=False); B["ESTIMATE_immune"].to_excel(w,sheet_name="ESTIMATE_immune",index=False)
    SC.to_excel(w,sheet_name="SingleCell_programmes",index=False); SI.to_excel(w,sheet_name="Other_definitions",index=False); flags.to_excel(w,sheet_name="Citations_to_confirm",index=False)
# ---- S25 software, parameters, seeds, inputs
pk=["python","pandas","numpy","scipy","scikit-learn","statsmodels","lifelines","scikit-survival","anndata","umap-learn","matplotlib","seaborn","pyreadr","openpyxl"]
sw=pd.DataFrame([dict(component=p,version=(platform.python_version() if p=="python" else md.version(p))) for p in pk])
par=pd.DataFrame([
 ("Cohort","Exclusion","6 non-diffuse-HGG specimens removed before any analysis (T_excluded_specimens)","—"),
 ("Clustering","Consensus clustering (final labels)","k = 3; 500 bootstraps; 80% item resampling; KMeans n_init = 10 (random_state = seed + b); average linkage on 1 − consensus; label = mode over seeds","seeds 42–66"),
 ("Clustering","State naming","highest mean T-cell cytotoxicity = Immune-high, lowest = Immune-low","—"),
 ("Clustering","k sweep (Figure 3)","k = 2–10; 1,000 bootstraps; KMeans n_init = 1; PAC (0.1–0.9), silhouette on 1 − consensus, CDF AUC","seed 42"),
 ("Clustering","Gap statistic","PCA-rotated uniform reference, 50 references, KMeans n_init = 10, one-SE rule","seed 42"),
 ("Clustering","Within-group / library re-clustering","1,000 bootstraps, n_init = 1, k = 2–6; clusters named by centroid correlation (Hungarian)","seed 42"),
 ("Figure 4","Contingency tests","Pearson χ²; Monte-Carlo P, 9,999 label permutations","20260722"),
 ("Figure 4","Sequential PERMANOVA","Euclidean distance, 34-feature TPM matrix; additive model; 4,999 free permutations","20260722"),
 ("Table 1","Monte-Carlo χ²","9,999 permutations; BH across molecular rows","20260722"),
 ("Figure 6","Adjusted DEG","OLS on log2(TPM+1), state + molecular group; genes TPM ≥ 1 in ≥ 5%; BH; |log2FC| > 1","—"),
 ("Figure 5B","Cell-class attribution","ScPCA SCPCP000001 reference (19 donors, 9,680 cells, malignant capped at 300/donor); log2(CPM/10 + 1); 200 expression-matched random draws (20 bins)","42"),
 ("Figure 5C–E","Visium spot classes","nearest bulk-state centroid on shared genes; z within section","—"),
 ("Survival","Primary Cox","immune state + molecular group + age + sex; corrected OS (4 records re-timed by one rule)","—"),
 ("Survival","Incremental value","Harrell C optimism by 500 bootstrap resamples; ΔC bootstrap CI","20261101"),
 ("Survival","k = 2 vs k = 3 optimism","300 bootstraps","7"),
 ("ScPCA linkage","Donor-level tests","Spearman and OLS slope with 10,000 axis permutations; JT trend 2,000 permutations; classifier 5-fold CV","42 / 0"),
],columns=["area","analysis","parameters","seed"])
inputs=[("/home/claude/opc/histologies.tsv","OpenPedCan v15 histologies"),("/home/claude/opc/gene-expression-rsem-tpm-collapsed.rds","OpenPedCan v15 TPM"),("/home/claude/opc/fusion-putative-oncogenic.tsv","OpenPedCan v15 fusions"),
        ("/home/claude/repo/cancers_submission_2026/reproducibility_data/feature_matrices/main_cohort_LM22_ssGSEA_z.tsv","46-feature matrix (repository)"),("/home/claude/repo/cancers_submission_2026/reproducibility_data/ecotype_assignments.tsv","submitted locked labels (repository)")]
inp=pd.DataFrame([dict(file=p.split("/")[-1],description=d,md5=md5(p)) for p,d in inputs])
osfix=pd.DataFrame([dict(participant=k,OS_days_used=v) for k,v in OS_FIX.items()])
with pd.ExcelWriter(D+"Supplementary_Table_S25_software_parameters_seeds_n343.xlsx") as w:
    sw.to_excel(w,sheet_name="Software",index=False); par.to_excel(w,sheet_name="Parameters_and_seeds",index=False); inp.to_excel(w,sheet_name="Input_checksums",index=False)
    rd("T_excluded_specimens.tsv").to_excel(w,sheet_name="Excluded_specimens",index=False); osfix.to_excel(w,sheet_name="OS_corrections",index=False)
    rd("T_seed_partitions_n343.tsv").to_excel(w,sheet_name="Seed_partitions",index=False); rd("T_seed_stability.tsv").to_excel(w,sheet_name="Seed_stability",index=False)
# ---- Supplementary tables workbook
sheets=[
 ("README",None,"Supplementary tables for the n = 343 revision. State names: Immune-high, Immune-intermediate, Immune-low. Proposed numbering is provisional; Citation_map lists every sheet."),
 ("Immune_state_assignments","immune_state_assignments_n343.tsv","Per-tumour final state (modal over seeds 42–66), n = 344 label carried over, submitted locked label, modal frequency"),
 ("Excluded_specimens","T_excluded_specimens.tsv","Six specimens removed before any analysis"),
 ("Seed_stability","T_seed_stability.tsv","Stability of the modal partition"),
 ("Table1","T_Table1_n343.tsv","Clinical and molecular characteristics by immune state"),
 ("Fig2B_cross_method","T_fig2B_cross_method_n343.tsv","LM22 vs other methods, Spearman"),
 ("Fig3_k_sweep","fig3/ksweep_n343.tsv","Consensus criteria k = 2–10"),("Fig3_gap","fig3/gap_n343.tsv","Gap statistic"),
 ("Fig4_statistics","T_fig4_statistics.tsv","Location / molecular group / PERMANOVA"),
 ("Fig4A_counts","T_fig4A_counts.tsv",""),("Fig4B_counts","T_fig4B_counts.tsv",""),("Fig4C_counts","T_fig4C_counts.tsv",""),
 ("Fig5B_designs","fig5b/T_fig5B_designs.tsv","Cell class of highest expression"),("Fig5B_gene_level","fig5b/T_fig5B_gene_level.tsv",""),
 ("Fig5_Visium_by_section","T_visium_by_section.tsv","Spot classes by section"),
 ("Fig6_DEG_all_contrasts","T_adjusted_DEG_n343_all_contrasts.tsv","Molecular-group-adjusted DEG"),
 ("Fig6_shared_genes",glob.glob(O+"T_shared_genes_*.tsv")[0].replace(O,""),"Shared up-regulated genes"),
 ("Fig7_primary_Cox","T_primary_cox.tsv","Pre-specified primary model"),("Fig7_PH_test","T_primary_cox_PH_test.tsv","Schoenfeld tests"),
 ("Survival_sensitivity","T_survival_sensitivity.tsv","All sensitivity analyses"),("Autopsy_by_state","T_autopsy_by_state.tsv","Autopsy specimens by state"),
 ("Followup","T_followup.tsv","Follow-up including ranges"),("OS_origin_audit_all","T_OS_origin_audit_all_n343.tsv","OS-origin check for every tumour"),("OS_origin_audit_summary","T_OS_origin_audit_summary.tsv",""),
 ("Immune_axis_Cox","T_immune_axis_cox.tsv","Continuous immune axis"),("Immune_axis_loadings","immune_axis_loadings.tsv",""),
 ("Within_subtype_Cox","T_within_subtype_state_cox.tsv",""),("Within_subtype_summary","T_within_subtype_survival_summary.tsv",""),
 ("Incremental_value","T_incremental_value.tsv",""),
 ("Reclustering_DMG_sweep","r54/R54_sweep_DMG.tsv",""),("Reclustering_WT_sweep","r54/R54_sweep_pHGG-WT.tsv",""),("Reclustering_DMG_crosstab","r54/R54_crosstab_DMG.tsv",""),("Reclustering_WT_crosstab","r54/R54_crosstab_pHGG-WT.tsv",""),("Reclustering_residualised","r54/R54_crosstab_residualised.tsv",""),("Variance_partition","r54/R54_variance_partition.tsv",""),
 ("k2_vs_k3_survival","r53/R53_survival_k2_k3.tsv",""),("k2_vs_k3_nesting","r53/R53_nesting.tsv",""),("High_vs_Intermediate_features","r53/R53_high_vs_intermediate_features.tsv",""),
 ("Library_survival","r6lib/r6/R6_library_survival.tsv",""),("Library_variance","r6lib/r6/R6_library_variance_partition.tsv",""),
 ("Feature_and_marker_Cox","q1/Q1_cox_markers_and_immune_features.tsv","Markers, 46 features, 17 programmes; base and + axis/state"),
 ("Programmes_by_state","q2/Q2_programmes_by_state.tsv","Curated programmes by state"),
 ("Published_programmes_Cox","T_published_programmes_cox.tsv",""),("Published_by_state_group","T_published_programmes_by_state_group.tsv",""),("Gene_screen","T_gene_screen_shared_genes.tsv",""),
 ("Deconvolution_by_state","T_deconvolution_by_state.tsv",""),("Cross_method_ordering","r52/R52_state_ordering_percentile.tsv",""),("Cross_method_concordance","r52/R52_concordance_344.tsv","n = 343 (file name inherited)"),
 ("Thorsson_by_state","T_Thorsson_subtype_by_state.tsv",""),("MFP_by_state","T_Bagaev_MFP_by_state.tsv",""),
 ("ScPCA_linkage_tests","T_scpca_axis_linkage_donor_tests.tsv",""),("ScPCA_projection","T_scpca_donors_projection.tsv",""),
 ("EXPLORATORY_tumour_cell","T_EXPLORATORY_tumour_cell_programmes_by_group.tsv","Exploratory; deconvolution-derived; not claimed in the main text")]
cm=[]
with pd.ExcelWriter(D+"Supplementary_Tables_n343.xlsx") as w:
    for nm,f,desc in sheets:
        if f is None: continue
        t=rd(f); t.columns=[str(c).replace("immune_PC1","immune_axis").replace("_PC1","_immune_axis").replace("PC1","immune_axis") for c in t.columns]; t=t.replace({"+ continuous immune PC1":"+ continuous immune axis"}); t.to_excel(w,sheet_name=nm[:31],index=False); cm.append(dict(sheet=nm[:31],source=f,rows=len(t),description=desc))
    pd.DataFrame(cm).to_excel(w,sheet_name="Citation_map",index=False)
    x=pd.read_excel(ROOT+"out/T_singlecell_per_patient_summary.xlsx",sheet_name=None)
    x["GSE227983_per_patient"].to_excel(w,sheet_name="SC_GSE227983_per_patient",index=False); x["ScPCA_per_patient"].to_excel(w,sheet_name="SC_ScPCA_per_patient",index=False)
print(len(cm),"sheets")
