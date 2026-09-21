"""Gap statistic (Tibshirani 2001) — the only internal index defined at k=1.
Uniform reference over the PCA-rotated bounding box, B=50 reference sets."""
import numpy as np, pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
UP="/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs"
feat=pd.read_csv(f"{UP}/step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv",sep="\t",index_col=0)
main=pd.read_csv(f"{UP}/cohort_main_final.tsv",sep="\t")["Kids_First_Biospecimen_ID"]
X=feat.loc[[s for s in main if s in feat.index]].dropna(axis=0,how="any").values
n,p=X.shape; rng=np.random.default_rng(42)

def Wk(D,k):
    if k==1:
        c=D.mean(0); return float(((D-c)**2).sum())
    km=KMeans(n_clusters=k,n_init=10,random_state=42).fit(D)
    return float(sum(((D[km.labels_==i]-km.cluster_centers_[i])**2).sum() for i in range(k)))

pca=PCA().fit(X); Xp=X@pca.components_.T
lo,hi=Xp.min(0),Xp.max(0)
B=50; ks=range(1,11)
logW=np.array([np.log(Wk(X,k)) for k in ks])
ref=np.zeros((B,len(list(ks))))
for b in range(B):
    Zp=rng.uniform(lo,hi,size=(n,p)); Z=Zp@pca.components_
    ref[b]=[np.log(Wk(Z,k)) for k in ks]
gap=ref.mean(0)-logW
sk=ref.std(0)*np.sqrt(1+1/B)
out=pd.DataFrame({"k":list(ks),"logW":logW.round(4),"gap":gap.round(4),"s_k":sk.round(4)})
# Tibshirani rule: smallest k with gap(k) >= gap(k+1) - s(k+1)
crit=[]
for i in range(len(out)-1):
    crit.append(bool(out.gap[i] >= out.gap[i+1]-out.s_k[i+1]))
crit.append(False)
out["meets_1SE_rule"]=crit
sel=out.loc[out.meets_1SE_rule,"k"]
print(out.to_string(index=False))
print("\nTibshirani 1-SE rule selects k =", int(sel.iloc[0]) if len(sel) else "none in 1..10")
out.to_csv("/tmp/claude-0/w1/T1_gap_statistic.tsv",sep="\t",index=False)
