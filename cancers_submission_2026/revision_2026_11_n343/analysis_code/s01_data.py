"""Per-tumour data objects for n = 343: TPM, DNA markers, immune-presence axis (PC1 of 21 immune ssGSEA signatures),
curated/published immune programmes, ESTIMATE ImmuneScore."""
import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
import pyreadr
from sklearn.decomposition import PCA
R="/home/claude/repo/cancers_submission_2026/reproducibility_data/"
a=load(); ids=a.index.tolist()
df=next(iter(pyreadr.read_r("/home/claude/opc/gene-expression-rsem-tpm-collapsed.rds").values()))
tpm=df[ids].astype(np.float32); tpm.to_pickle(ROOT+"tpm343.pkl"); print("TPM",tpm.shape)
# DNA markers: per-specimen matrix built in Rev344_05 (genomics.py/markers.py; OpenPedCan v15 consensus MAF + CNV) restricted to n = 343
M=pd.read_excel("/home/claude/w/pg/Revision_round1_Q1-Q5/Rev344_05_tables_Q1-Q5.xlsx",sheet_name="Q1_markers_matrix").set_index("Kids_First_Biospecimen_ID").loc[ids]
M.to_pickle(ROOT+"markers343.pkl"); M.to_csv(ROOT+"out/molecular_markers_n343.tsv",sep="\t")
# immune-presence axis
F=pd.read_csv(R+"feature_matrices/main_cohort_LM22_ssGSEA_z.tsv",sep="\t",index_col=0).loc[ids]
F=F.rename(columns={"ssGSEA_Glioma_Inflammatory_Wang2017":"ssGSEA_Cytotoxic_Lymphocyte_Curated"})
NONIMM=["ssGSEA_Cell_Cycle_Proliferation","ssGSEA_Stemness_Brain_Tumor","ssGSEA_MAPK_Activity"]
SSI=[c for c in F.columns if c.startswith("ssGSEA") and c not in NONIMM]
Xi=(F[SSI]-F[SSI].mean())/F[SSI].std()
pca=PCA(5).fit(Xi.values); pc1=pd.Series(pca.transform(Xi.values)[:,0],index=F.index)
if pc1.corr(F["ssGSEA_T_Cell_Cytotoxicity"])<0: pc1=-pc1
pc1=(pc1-pc1.mean())/pc1.std()
pc1.to_frame("immune_axis").to_csv(ROOT+"out/immune_axis_n343.tsv",sep="\t")
load_=pd.Series(pca.components_[0],index=SSI)*(1 if pca.transform(Xi.values)[:,0].dot(F["ssGSEA_T_Cell_Cytotoxicity"].values)>0 else -1)
pd.DataFrame({"signature":SSI,"PC1_loading":load_.values}).to_csv(ROOT+"out/immune_axis_loadings.tsv",sep="\t",index=False)
pd.DataFrame({"component":[f"PC{i+1}" for i in range(5)],"variance_explained":pca.explained_variance_ratio_}).to_csv(ROOT+"out/immune_axis_scree.tsv",sep="\t",index=False)
print("immune axis PC1 variance explained",round(pca.explained_variance_ratio_[0],4))
# programmes (curated, Rev344 programs.py) + published sets (Rev344_06 published_sets.py)
sys.path.insert(0,"/home/claude/w/c05/code"); import programs as PR
S,meta=PR.compute(tpm); S.to_pickle(ROOT+"programs343.pkl"); S.to_csv(ROOT+"out/immune_programme_scores_n343.tsv",sep="\t"); meta.to_csv(ROOT+"out/immune_programme_definitions.tsv",sep="\t",index=False)
# ESTIMATE ImmuneScore (Yoshihara 2013)
gs=next(iter(pyreadr.read_r(ROOT+"src/ext_estimate_gene_sets.rda").values())); imm=gs["immune_signature"].tolist()
X=np.log2(tpm+1); X=X[(tpm>0.5).sum(1)>=0.05*tpm.shape[1]]
def ssgsea_estimate(X,genes,alpha=0.25):
    out={}
    for s in X.columns:
        v=X[s].sort_values(ascending=False); r=np.arange(len(v),0,-1); inset=v.index.isin(genes); w=np.abs(r)**alpha
        Phit=np.cumsum(np.where(inset,w,0))/w[inset].sum(); Pmiss=np.cumsum(~inset)/(~inset).sum(); out[s]=(Phit-Pmiss).sum()
    return pd.Series(out)
est=ssgsea_estimate(X,imm); est=(est-est.mean())/est.std(); est.to_frame("ESTIMATE_ImmuneScore_z").to_csv(ROOT+"out/ESTIMATE_n343.tsv",sep="\t")
from scipy.stats import spearmanr
print("axis vs ESTIMATE rho",round(spearmanr(pc1,est)[0],3))
