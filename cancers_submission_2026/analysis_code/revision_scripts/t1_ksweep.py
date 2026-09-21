"""T1: consensus k-sweep k=2..10 on the primary 46-feature matrix.
Method-matched to step5e_bootstrap_1000.py (B=1000, pItem=0.8, KMeans n_init=1, seed=42)."""
import numpy as np, pandas as pd, time, json
from collections import Counter
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score, calinski_harabasz_score, davies_bouldin_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

UP="/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs"
feat=pd.read_csv(f"{UP}/step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv",sep="\t",index_col=0)
main=pd.read_csv(f"{UP}/cohort_main_final.tsv",sep="\t")["Kids_First_Biospecimen_ID"]
main_ids=[s for s in main if s in feat.index]
X_df=feat.loc[main_ids].dropna(axis=0,how="any")
X=X_df.values; n=X.shape[0]
print("matrix:",X.shape,"| features:",X_df.shape[1])

def consensus(k,B=1000,pItem=0.8,seed=42):
    M=np.zeros((n,n)); C=np.zeros((n,n))
    rng=np.random.default_rng(seed)
    for b in range(B):
        idx=rng.choice(n,int(n*pItem),replace=False)
        lbl=KMeans(n_clusters=k,n_init=1,random_state=seed+b).fit_predict(X[idx])
        C[np.ix_(idx,idx)]+=1
        for ci in range(k):
            mem=idx[lbl==ci]
            if len(mem)>=2: M[np.ix_(mem,mem)]+=1
    with np.errstate(invalid="ignore",divide="ignore"):
        cm=np.where(C>0,M/C,0.0)
    np.fill_diagonal(cm,1.0)
    dist=1.0-cm; np.fill_diagonal(dist,0.0)
    Z=linkage(squareform(dist,checks=False),method="average")
    lab=fcluster(Z,t=k,criterion="maxclust")
    return cm,dist,lab

def cdf_auc(cm):
    """Area under the consensus CDF (Monti 2003)."""
    x=np.sort(cm[np.triu_indices(n,1)])
    # AUC = sum (x_i - x_{i-1}) * CDF(x_i)
    cdf=np.arange(1,len(x)+1)/len(x)
    return float(np.sum(np.diff(np.concatenate(([0.0],x)))*cdf))

rows=[]; prev_auc=None; store={}
for k in range(2,11):
    t0=time.time()
    cm,dist,lab=consensus(k)
    up=cm[np.triu_indices(n,1)]
    pac=float(((up>0.1)&(up<0.9)).mean())
    sil=float(silhouette_score(dist,lab,metric="precomputed"))
    same=np.equal.outer(lab,lab)&~np.eye(n,dtype=bool)
    wcr=float(cm[same].mean())
    auc=cdf_auc(cm)
    dauc=float((auc-prev_auc)/prev_auc) if prev_auc else np.nan
    prev_auc=auc
    sizes=sorted(Counter(lab).values(),reverse=True)
    rows.append(dict(k=k,PAC=round(pac,4),silhouette_consensus=round(sil,4),
                     mean_within_cluster_consensus=round(wcr,4),
                     consensus_CDF_AUC=round(auc,4),delta_AUC=round(dauc,4) if dauc==dauc else np.nan,
                     calinski_harabasz=round(calinski_harabasz_score(X,lab),1),
                     davies_bouldin=round(davies_bouldin_score(X,lab),3),
                     min_cluster_size=min(sizes),n_clusters_ge10=sum(s>=10 for s in sizes),
                     sizes=";".join(map(str,sizes))))
    store[k]=lab
    print(f"k={k:2d} PAC={pac:.4f} sil={sil:.4f} sizes={sizes} ({time.time()-t0:.0f}s)")

sweep=pd.DataFrame(rows)

# --- GATE: reproduce locked B1000 values for k=2..6 ---
locked=pd.DataFrame({"k":[2,3,4,5,6],
 "PAC_locked":[0.24264730099133813,0.3841023614267365,0.47584230807232486,0.49412113427526927,0.46390343510193327],
 "sil_locked":[0.9361341280993414,0.7902823685847833,0.5451144227462578,0.5076222778063351,0.4241097394767781]})
chk=sweep.merge(locked,on="k")
chk["PAC_diff"]=(chk.PAC-chk.PAC_locked).abs()
chk["sil_diff"]=(chk.silhouette_consensus-chk.sil_locked).abs()
print("\n=== REPRODUCIBILITY GATE (k=2..6 vs locked B1000) ===")
print(chk[["k","PAC","PAC_locked","PAC_diff","silhouette_consensus","sil_locked","sil_diff"]].to_string(index=False))
ok=bool((chk.PAC_diff<5e-3).all() and (chk.sil_diff<5e-3).all())
print("GATE PASS:",ok)

# --- ARI of de novo k=3 vs locked ecotype labels ---
eco=pd.read_csv(f"{UP}/ecotype_LM22_main_k3_annotated.tsv",sep="\t").set_index("Kids_First_Biospecimen_ID")["ecotype"]
eco=eco.reindex(X_df.index)
ari=adjusted_rand_score(eco,store[3])
print(f"\nde novo k=3 vs locked ecotype: ARI = {ari:.4f}")
print(pd.crosstab(eco,store[3]))

sweep.to_csv("/tmp/claude-0/w1/T1_ksweep_k2_k10.tsv",sep="\t",index=False)
json.dump({"gate_pass":ok,"ari_k3_vs_locked":ari,"n":int(n),"n_features":int(X_df.shape[1])},
          open("/tmp/claude-0/w1/T1_meta.json","w"),indent=1)
print("\n"+sweep.to_string(index=False))
