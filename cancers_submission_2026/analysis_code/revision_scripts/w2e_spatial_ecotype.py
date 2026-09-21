"""C (continued) — assign each Visium spot to the nearest bulk ecotype centroid in the
24-signature space, and test whether the assignments form spatially coherent domains."""
import numpy as np, pandas as pd, pickle
from scipy import sparse, stats
rng=np.random.default_rng(42)
ANN="/mnt/user-data/uploads/Open PBTA/Revision/FINAL MANUSCRIPT 260722 - 수정본/7. Reproducibility data/Sample annotation/sample_master_annotation.tsv"
d=pd.read_csv(ANN,sep="\t").set_index("sample")
sigcols=[c for c in d.columns if c.startswith("ssGSEA_")]
Zb=(d[sigcols]-d[sigcols].mean())/d[sigcols].std()
Zb.columns=[c.replace("ssGSEA_","") for c in Zb.columns]
CENT=Zb.groupby(d["ecotype"]).mean()
ORDER=["Lymphocyte-inflamed","Myeloid-dominant","Immune-desert"]
CENT=CENT.loc[ORDER]
print("bulk ecotype centroids, 24-signature space:",CENT.shape)

spot=pickle.load(open("w2_spot_scores.pkl","rb"))
def neighbours(pos):
    r=pos.array_row.values; c=pos.array_col.values
    key={(a,b):i for i,(a,b) in enumerate(zip(r,c))}
    n=len(r); W=sparse.lil_matrix((n,n))
    for i,(a,b) in enumerate(zip(r,c)):
        for da,db in [(0,-2),(0,2),(-1,-1),(-1,1),(1,-1),(1,1)]:
            j=key.get((a+da,b+db))
            if j is not None: W[i,j]=1
    return W.tocsr()

rows=[];colo=[]
for s,(S,pos) in spot.items():
    cols=[c for c in CENT.columns if c in S.columns]
    Zs=(S[cols]-S[cols].mean())/S[cols].std()
    C=CENT[cols]
    # nearest centroid by Pearson correlation across signatures
    A=Zs.values; B=C.values
    A=(A-A.mean(1,keepdims=True))/A.std(1,keepdims=True)
    B=(B-B.mean(1,keepdims=True))/B.std(1,keepdims=True)
    corr=A@B.T/A.shape[1]
    assign=np.array(ORDER)[corr.argmax(1)]
    W=neighbours(pos)
    # join-count: fraction of neighbouring spot pairs with the same assignment
    src,dst=W.nonzero()
    same=(assign[src]==assign[dst]).mean()
    perm=np.array([(assign[rng.permutation(len(assign))][src]==assign[rng.permutation(len(assign))][dst]).mean()
                   for _ in range(999)])
    # correct null: permute once per replicate
    perm=[]
    for _ in range(999):
        p=rng.permutation(assign); perm.append((p[src]==p[dst]).mean())
    perm=np.array(perm); pval=(np.sum(perm>=same)+1)/1000
    counts=pd.Series(assign).value_counts().reindex(ORDER).fillna(0).astype(int)
    rows.append(dict(sample=s,n_spots=len(assign),
                     **{k.split("-")[0][:4]:int(v) for k,v in counts.items()},
                     n_ecotypes_present=int((counts>0.05*len(assign)).sum()),
                     same_neighbour_fraction=round(same,3),
                     null_mean=round(perm.mean(),3),perm_p=pval))
    my=[c for c in ["Microglia_Klemm2020","MDM_Klemm2020","MgTAM_Antunes2021","MoTAM_Antunes2021",
                    "DAM_KerenShaul2017","M2_Macrophage"] if c in S.columns]
    ly=[c for c in ["T_Cell_Cytotoxicity","NK_Cell_Activity","Chemokine_T_Cell_Recruitment",
                    "T_Cell_Exhaustion","Tregs_Friebel2020"] if c in S.columns]
    r=stats.pearsonr(S[my].mean(axis=1),S[ly].mean(axis=1))
    colo.append(dict(sample=s,n_spots=len(assign),myeloid_lymphoid_r=round(r.statistic,3),p=r.pvalue))
    np.save(f"w2_spot_assign_{s}.npy",assign)
R=pd.DataFrame(rows); Cc=pd.DataFrame(colo)
print("\n=== spot-level ecotype assignment and spatial coherence ===")
print(R.to_string(index=False))
print("\n=== myeloid x lymphoid programme co-localisation per spot ===")
print(Cc.to_string(index=False))
print(f"\nbulk myeloid x lymphoid theme correlation for comparison: "
      f"{np.corrcoef(Zb[[c for c in ['Microglia_Klemm2020','MDM_Klemm2020','MgTAM_Antunes2021','MoTAM_Antunes2021','DAM_KerenShaul2017','M2_Macrophage']]].mean(axis=1), Zb[['T_Cell_Cytotoxicity','NK_Cell_Activity','Chemokine_T_Cell_Recruitment','T_Cell_Exhaustion','Tregs_Friebel2020']].mean(axis=1))[0,1]:.3f}")
R.to_csv("W2_C_spot_ecotype.tsv",sep="\t",index=False); Cc.to_csv("W2_C_colocalisation.tsv",sep="\t",index=False)
