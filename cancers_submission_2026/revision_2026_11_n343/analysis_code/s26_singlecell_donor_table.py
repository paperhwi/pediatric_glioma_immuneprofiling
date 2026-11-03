"""Per-patient single-cell summary (GSE227983 LaBelle et al. 2025; ScPCA SCPCP000001): one row per donor and platform with
cell numbers and the key functional measures used in the text. Bulk-cohort-independent values are carried over unchanged;
ScPCA immune axis / projected immune state use the n = 343 reference (s22)."""
import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
C="/mnt/user-data/uploads/Pediatric GBM/Revision_round1_Q1-Q5/Item_C_singlecell/C_results_singlecell.xlsx"
S="/mnt/user-data/uploads/Pediatric GBM/Revision_round1_Q1-Q5/Item_C_singlecell/ScPCA_validation/ScPCA_results.xlsx"
K8=["n_cells","% cytotoxic-high (GZMB+ or PRF1+)","% GZMB+","% PRF1+","% ≥2 inhibitory receptors","% MKI67+","% TCF7+","% TOX+","within-donor rho (cytotoxicity vs exhaustion)"]
K4=["n_cells","% FOXP3+","% CTLA4+","% ≥2 inhibitory receptors","% MKI67+"]
KM=["n_cells","% microglia-like","Immunosuppressive TAM (mean score)","Myeloid checkpoints (mean score)","MHC-II antigen presentation (mean score)","within-donor rho (MHC-II vs suppressive TAM)"]
def block(df,key,cols,pre):
    d=df.set_index(key)[cols].copy(); d.columns=[f"{pre}: {c}" for c in cols]; return d
# GSE227983
c8=pd.read_excel(C,"Donor_CD8"); c4=pd.read_excel(C,"Donor_CD4"); cm=pd.read_excel(C,"Donor_Myeloid")
for d in (c8,c4,cm): d["key"]=d.sampleid+" | "+d.platform
meta=pd.concat([c8,c4,cm])[["key","sampleid","platform","patient","Subtype","Age","Sex","Location","RadioChemoTreatBeforeSurgery","SteroidTreatment","ImmunotherapyTreatment"]].dropna(subset=["patient"]).drop_duplicates("key").set_index("key")
G=meta.join(block(c8,"key",K8,"CD8")).join(block(c4,"key",K4,"CD4")).join(block(cm,"key",KM,"Myeloid"))
G.insert(0,"cohort","GSE227983 (LaBelle et al. 2025)"); G["CD8 assessable (>=10 cells)"]=G["CD8: n_cells"]>=10; G["Myeloid assessable (>=10 cells)"]=G["Myeloid: n_cells"]>=10
# ScPCA
dn=pd.read_excel(S,"Donors").set_index("scpca_sample_id"); s8=pd.read_excel(S,"Donor_CD8"); s4=pd.read_excel(S,"Donor_CD4"); sm=pd.read_excel(S,"Donor_Myeloid")
P=pd.read_csv(ROOT+"out/T_scpca_donors_projection.tsv",sep="\t",index_col=0)
Q=dn[["diagnosis","subdiagnosis","age","tissue_location","location_group","analysis_set","total_cells","immune_%","T_%_of_immune","Myeloid","CD8","CD4"]].rename(columns={"Myeloid":"n_myeloid_cells","CD8":"n_CD8_cells","CD4":"n_CD4_cells"})
Q=Q.join(P[["immune_axis","projected_state","projected_state_prob"]].rename(columns={"immune_axis":"bulk immune axis (n = 343 reference)","projected_state":"projected immune state (n = 343 classifier)","projected_state_prob":"projected state probability"}))
Q=Q.join(block(s8,"sample",[c for c in K8 if c in s8],"CD8")).join(block(s4,"sample",K4,"CD4")).join(block(sm,"sample",KM,"Myeloid"))
Q.insert(0,"cohort","ScPCA SCPCP000001"); Q["CD8 assessable (>=10 cells)"]=Q["CD8: n_cells"]>=10; Q["Myeloid assessable (>=10 cells)"]=Q["Myeloid: n_cells"]>=10
with pd.ExcelWriter(ROOT+"out/T_singlecell_per_patient_summary.xlsx") as w:
    pd.DataFrame({"README":["One row per donor (GSE227983: per donor and platform, SS2 = Smart-seq2 sorted immune cells, TENX = 10x).",
      "Measures: CD8 % cytotoxic-high (GZMB+ or PRF1+), % with >=2 inhibitory receptors, % MKI67+, % TCF7+; CD4 % FOXP3+, % CTLA4+; myeloid % microglia-like, suppressive-TAM and checkpoint programme means (z within compartment).",
      "Donors with < 10 cells in a compartment are shown but were not used for that compartment's tests.",
      "ScPCA: bulk immune axis and projected immune state are recomputed with the n = 343 reference (Immune-high / Immune-intermediate / Immune-low). Normal-margin and PXA-lineage samples are listed with analysis_set and were excluded from the primary set.",
      "GSE227983 sorted immune cells have no survival data and no bulk RNA; ependymomas and checkpoint-treated recurrences were excluded from inference."]}).to_excel(w,sheet_name="README",index=False)
    G.reset_index(drop=True).to_excel(w,sheet_name="GSE227983_per_patient",index=False); Q.reset_index().rename(columns={"index":"scpca_sample_id"}).to_excel(w,sheet_name="ScPCA_per_patient",index=False)
print(G.shape,Q.shape); print(Q[["analysis_set","total_cells","immune_%","projected immune state (n = 343 classifier)"]].to_string())
