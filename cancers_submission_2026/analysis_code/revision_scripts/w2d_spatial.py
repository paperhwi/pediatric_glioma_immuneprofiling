"""C — do the immune programs form spatially coherent domains in pDMG tissue?
GSE268577, three H3.3K27M pDMG Visium sections. Spots are unsorted, so composition is meaningful."""
import numpy as np, pandas as pd, gzip, glob, io
from scipy import sparse, stats
np.random.seed(42)
UP="/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs"
sigs={}
for line in open(f"{UP}/brain_immune_signatures.gmt"):
    f=line.rstrip("\n").split("\t"); sigs[f[0]]=[g for g in f[2:] if g.strip()]

def load(sample):
    base=glob.glob(f"data/visium/*_pDMG_Sample-{sample}_")[0] if False else \
         glob.glob(f"data/visium/*_pDMG_Sample-{sample}_matrix.mtx.gz")[0].replace("matrix.mtx.gz","")
    feat=pd.read_csv(base+"features.tsv.gz",sep="\t",header=None)
    bc=pd.read_csv(base+"barcodes.tsv.gz",header=None)[0].values
    with gzip.open(base+"matrix.mtx.gz","rt") as f:
        from scipy.io import mmread
        M=mmread(f).tocsc()
    pos=pd.read_csv(base+"tissue_positions.csv.gz").set_index("barcode")
    pos=pos.loc[bc]
    keep=pos.in_tissue.values==1
    M=M[:,keep]; pos=pos[keep]; bc=bc[keep]
    A=np.asarray(M.todense(),dtype=np.float32)          # 17,943 x n_spots, small
    X=pd.DataFrame(A,index=feat[1].values,columns=bc)
    X=X.groupby(level=0).sum()                           # collapse duplicate symbols
    tot=X.sum(axis=0)
    L=np.log2(X/tot*1e4+1).astype(np.float32)
    return L,pos,tot

def neighbours(pos):
    r=pos.array_row.values; c=pos.array_col.values
    n=len(r); W=sparse.lil_matrix((n,n))
    key={(a,b):i for i,(a,b) in enumerate(zip(r,c))}
    for i,(a,b) in enumerate(zip(r,c)):
        for da,db in [(0,-2),(0,2),(-1,-1),(-1,1),(1,-1),(1,1)]:
            j=key.get((a+da,b+db))
            if j is not None: W[i,j]=1
    return W.tocsr()

def morans_I(x,W):
    x=np.asarray(x,dtype=float); z=x-x.mean(); n=len(x)
    S0=W.sum(); den=float((z**2).sum())
    return (n/S0)*(float(z @ (W @ z))/den) if den>0 and S0>0 else np.nan

def morans_I_perm(x,W,B,rng):
    """Vectorised null: B column-permutations at once."""
    x=np.asarray(x,dtype=float); n=len(x); S0=W.sum()
    idx=np.argsort(rng.random((B,n)),axis=1)
    Z=x[idx].T                                   # n x B
    Z=Z-Z.mean(axis=0,keepdims=True)
    num=np.einsum("ij,ij->j",Z,W@Z)
    den=(Z**2).sum(axis=0)
    return (n/S0)*(num/den)

rows=[]; spot_scores={}
for s in ["1","2","3"]:
    L,pos,tot=load(s)
    W=neighbours(pos)
    deg=np.asarray(W.sum(axis=1)).ravel()
    print(f"Sample-{s}: {L.shape[1]} in-tissue spots, {L.shape[0]} genes, "
          f"median {np.median(tot):.0f} counts/spot, mean {deg.mean():.1f} neighbours")
    mean_expr=L.mean(axis=1); bins=pd.qcut(mean_expr.rank(method="first"),25,labels=False)
    rng=np.random.default_rng(42)
    def score(genes,nctrl=50):
        g=[x for x in genes if x in L.index]
        if len(g)<3: return None
        ctrl=[]
        for x in g:
            pool=mean_expr.index[bins==bins[x]]
            ctrl+=list(rng.choice(pool,size=min(nctrl,len(pool)),replace=False))
        return (L.loc[g].mean(axis=0)-L.loc[ctrl].mean(axis=0)).values
    S={k:score(v) for k,v in sigs.items()}
    S={k:v for k,v in S.items() if v is not None}
    Sdf=pd.DataFrame(S,index=pos.index); spot_scores[s]=(Sdf,pos)
    for name,v in S.items():
        I=morans_I(v,W)
        perm=morans_I_perm(v,W,999,rng)
        p=(np.sum(perm>=I)+1)/(len(perm)+1)
        rows.append(dict(sample=f"Sample-{s}",signature=name,n_spots=L.shape[1],
                         morans_I=round(I,3),perm_p=p,
                         perm_mean=round(perm.mean(),4),perm_sd=round(perm.std(),4)))
R=pd.DataFrame(rows)
R["q_BH"]=stats.false_discovery_control(R.perm_p)
R.to_csv("W2_C_morans_I.tsv",sep="\t",index=False)
print("\n=== Moran's I by signature, mean across the three sections ===")
piv=R.pivot(index="signature",columns="sample",values="morans_I")
piv["mean_I"]=piv.mean(axis=1)
sig_all=R.groupby("signature").q_BH.max()
piv["max_q"]=sig_all
print(piv.sort_values("mean_I",ascending=False).round(3).to_string())
print(f"\nsignatures spatially autocorrelated (q<0.05) in all three sections: "
      f"{int((R.groupby('signature').q_BH.max()<0.05).sum())}/{R.signature.nunique()}")
import pickle; pickle.dump(spot_scores,open("w2_spot_scores.pkl","wb"))
