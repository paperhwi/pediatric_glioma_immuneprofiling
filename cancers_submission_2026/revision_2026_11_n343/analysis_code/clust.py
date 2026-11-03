"""Locked n=344 consensus clustering (cluster344.py settings: B=500, pItem=0.8, KMeans n_init=10,
average linkage on 1-consensus, k=3) for a list of seeds; states named by mean T-cell cytotoxicity
(highest = Immune-high, lowest = Immune-low, remaining = Immune-intermediate)."""
import pandas as pd, numpy as np, sys, json
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
R='op/Revision/pediatric_glioma_immuneprofiling/Pediatric Glioma Immune Profiling Source Folder/reproducibility_data/'
feat=pd.read_csv(R+'feature_matrices/main_cohort_LM22_ssGSEA_z.tsv',sep='\t',index_col=0)
def consensus3(X,seed,B=500,pItem=0.8,k=3):
    RNG=np.random.default_rng(seed); n=len(X); M=np.zeros((n,n)); C=np.zeros((n,n))
    for b in range(B):
        idx=RNG.choice(n,int(n*pItem),replace=False)
        lbl=KMeans(n_clusters=k,n_init=10,random_state=seed+b).fit_predict(X[idx])
        for ci in range(k):
            m=idx[lbl==ci]
            if len(m)<2: continue
            M[np.ix_(m,m)]+=1
        C[np.ix_(idx,idx)]+=1
    cm=np.where(C>0,M/np.maximum(C,1),0.); np.fill_diagonal(cm,1); d=1-cm; np.fill_diagonal(d,0)
    return fcluster(linkage(squareform(d,checks=False),'average'),t=k,criterion='maxclust')
def name(Xdf,lbl):
    m=Xdf.assign(c=lbl).groupby('c')['ssGSEA_T_Cell_Cytotoxicity'].mean().sort_values()
    nm={m.index[0]:'Immune-low',m.index[1]:'Immune-intermediate',m.index[2]:'Immune-high'}
    return pd.Series(lbl,index=Xdf.index).map(nm)
if __name__=='__main__':
    idfile,seeds,out=sys.argv[1],[int(s) for s in sys.argv[2].split(',')],sys.argv[3]
    ids=pd.read_csv(idfile,header=None)[0].tolist(); X=feat.loc[ids]
    res={}
    for s in seeds:
        res[s]=name(X,consensus3(X.values,s)); print(s,res[s].value_counts().to_dict(),flush=True)
        pd.DataFrame(res).to_csv(out,sep='\t')
