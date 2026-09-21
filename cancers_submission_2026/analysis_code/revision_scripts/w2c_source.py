"""B — are the genes that define the ecotype axis expressed by immune cells or by malignant cells?
Smart-seq2 only, so immune and tumour cells are compared on the same platform."""
import numpy as np, pandas as pd
from scipy import stats
UP="/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs"
PHGG=["H3K27M","H3WT_hemispheric","H3WT_midline","H3G34R/V"]

I=pd.read_pickle("w2_immune_ss2.pkl"); im=pd.read_csv("data/immune_ss2_meta.csv.gz",index_col=0)
T=pd.read_pickle("w2_tumor_ss2.pkl");  tm=pd.read_csv("data/tumor_ss2_meta.csv.gz",index_col=0)
im=im.loc[im.index.intersection(I.columns)]; im=im[im.Subtype.isin(PHGG)]
tm=tm.loc[tm.index.intersection(T.columns)]; tm=tm[tm.Subtype.isin(PHGG)]
print("SS2 immune cells (pHGG):",len(im),"|",im.broad_annot.value_counts().to_dict())
print("SS2 tumour cells (pHGG):",len(tm),"|",tm.Subtype.value_counts().to_dict())
if "CellAnnot" in tm: print("  tumour CellAnnot:",tm.CellAnnot.value_counts().head(6).to_dict())

genes=sorted(set(I.index)&set(T.index))
X=pd.concat([np.log2(I.loc[genes,im.index]/10+1), np.log2(T.loc[genes,tm.index]/10+1)],axis=1)
ct=pd.concat([im.broad_annot, pd.Series("Malignant",index=tm.index)])
ct=ct.loc[X.columns]
print("\ncombined SS2 matrix:",X.shape,"|",ct.value_counts().to_dict())

deg=pd.read_csv(f"{UP}/adjusted_DEG_all_contrasts.tsv",sep="\t")
A="Lymphocyte_inflamed_vs_Immune_desert"
up=deg[(deg.contrast==A)&(deg.q_BH<0.05)&(deg.adjusted_log2FC>1)].nsmallest(300,"q_BH")["gene"]
up=[g for g in up if g in X.index]
print(f"\necotype-defining genes (top 300 of the Lymphocyte-inflamed vs Immune-desert contrast) present in SS2: {len(up)}")

frac=(X.loc[up]>0).T.groupby(ct).mean().T     # detection rate per cell type
mean=X.loc[up].T.groupby(ct).mean().T
print("\n=== mean detection rate of the ecotype-defining genes, per cell type ===")
print(frac.mean().round(3).to_string())
print("\n=== mean log2 expression, per cell type ===")
print(mean.mean().round(3).to_string())

# per gene: which cell type expresses it most
src=mean.idxmax(axis=1)
print("\n=== cell type with the highest mean expression, per ecotype-defining gene ===")
print(src.value_counts().to_string())
print(f"\nimmune-attributed: {int(src.isin(['Myeloid','Tcell']).sum())}/{len(src)} "
      f"({src.isin(['Myeloid','Tcell']).mean()*100:.1f}%)")

# background: random genes matched on overall expression
rng=np.random.default_rng(42)
allmean=X.mean(axis=1); bins=pd.qcut(allmean.rank(method="first"),20,labels=False)
bg=[]
for _ in range(200):
    pick=[rng.choice(allmean.index[bins==bins[g]]) for g in up]
    s=X.loc[pick].T.groupby(ct).mean().T.idxmax(axis=1)
    bg.append(s.isin(["Myeloid","Tcell"]).mean())
bg=np.array(bg)
obs=src.isin(["Myeloid","Tcell"]).mean()
print(f"expression-matched random genes: {bg.mean()*100:.1f}% immune-attributed "
      f"(95% range {np.percentile(bg,2.5)*100:.1f}-{np.percentile(bg,97.5)*100:.1f}%)")
print(f"permutation P = {max((bg>=obs).mean(),1/len(bg)):.4g}")

out=pd.DataFrame({"gene":up,"top_cell_type":src.values})
for c in mean.columns: out[f"mean_{c}"]=mean[c].values.round(3)
for c in frac.columns: out[f"detect_{c}"]=frac[c].values.round(3)
out.to_csv("W2_B_gene_source.tsv",sep="\t",index=False)
pd.DataFrame({"observed_immune_fraction":[obs],"background_mean":[bg.mean()],
              "background_lo":[np.percentile(bg,2.5)],"background_hi":[np.percentile(bg,97.5)],
              "perm_P":[max((bg>=obs).mean(),1/len(bg))]}).to_csv("W2_B_summary.tsv",sep="\t",index=False)
X.to_pickle("w2_ss2_combined.pkl"); ct.to_pickle("w2_ss2_celltype.pkl")
