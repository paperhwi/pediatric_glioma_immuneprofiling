"""Figure 5C-E: Visium (GSE268577, 3 H3.3-K27M pDMG sections) spots assigned to the nearest bulk immune-state centroid on the
shared immune-presence genes (Figure 6), as in Rev344_07 r51_refdecon.py. Spot fractions (InstaPrism) are unchanged.
Step 1 reproduces the n = 344 spot classes (417 genes, n = 344 labels); step 2 uses n = 343 labels and the revised shared genes."""
import sys, glob, gzip; sys.path.insert(0,"/home/claude/q343/src"); from common import *
import anndata as ad, scipy.io as sio
from scipy import stats
V="/home/claude/visium/"
h=ad.read_h5ad("/mnt/user-data/uploads/Pediatric GBM/Revision_round1_Q1-Q5/Round2_R5_2Nov/R5_visium_spots.h5ad")
parts=[]
for i,g in enumerate(["GSM8294434","GSM8294435","GSM8294436"]):
    sec=f"Sample-{i+1}"; pre=glob.glob(V+g+"*_matrix.mtx.gz")[0].replace("matrix.mtx.gz","")
    M=sio.mmread(pre+"matrix.mtx.gz").tocsc(); feat=pd.read_csv(pre+"features.tsv.gz",sep="\t",header=None); bc=pd.read_csv(pre+"barcodes.tsv.gz",header=None)[0]
    df=pd.DataFrame.sparse.from_spmatrix(M,index=feat[1].values,columns=[f"{sec}|{b}" for b in bc]).sparse.to_dense()
    df=df.groupby(level=0).sum(); parts.append(df)
vc=pd.concat(parts,axis=1).fillna(0)[h.obs_names]; print("visium",vc.shape)
lc=np.log1p(vc/vc.sum(0)*1e4)
def classify(tpm,labels,genes):
    zb=np.log2(tpm.loc[[g for g in genes if g in vc.index]]+1); zb=zb.sub(zb.mean(1),axis=0).div(zb.std(1)+1e-9,axis=0)
    cent=pd.DataFrame({e:zb.loc[:,labels.values==e].mean(1) for e in ECO}); cls=pd.Series(index=vc.columns,dtype=object); score=pd.Series(index=vc.columns,dtype=float)
    for s_,idx in h.obs.groupby("section").groups.items():
        x=lc.loc[cent.index,idx]; x=x.sub(x.mean(1),axis=0).div(x.std(1)+1e-9,axis=0)
        r=pd.DataFrame({e:x.corrwith(cent[e]) for e in ECO}); cls.loc[idx]=r.idxmax(1); score.loc[idx]=x.mean(0)
    return cls,score,len(cent)
# step 1: reproduce n = 344
t343=pd.read_pickle(ROOT+"tpm343.pkl")
px=pd.read_csv("/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs/tpm_for_cibersortx.tsv",sep="\t",usecols=["GeneSymbol","BS_2DMBCHNF"]).drop_duplicates("GeneSymbol").set_index("GeneSymbol").BS_2DMBCHNF
df=t343.copy(); df["BS_2DMBCHNF"]=px.reindex(df.index).values
old=pd.read_excel("/home/claude/w/pg/Revision_round1_Q1-Q5/Rev344_05_tables_Q1-Q5.xlsx",sheet_name="ecotype_344").set_index("Kids_First_Biospecimen_ID").ecotype.map({"Lymphocyte-inflamed":HI,"Myeloid-dominant":MID,"Immune-desert":LO})
g417=pd.read_excel("/home/claude/w/pg/Revision_round1_Q1-Q5/Review_gaps_2Nov/R6_gap_analyses_results.xlsx","G1_fig4B_gene_level"); g417=g417.loc[g417.gene_set.str.startswith("417"),"gene"].tolist()
c0,s0,n0=classify(df[old.index],old,g417)
oldcls=h.obs.spot_class.astype(str).map({"LI":HI,"MD":MID,"ID":LO})
rep=dict(genes_used=n0,class_agreement=float((c0==oldcls).mean()),axis_score_r=float(np.corrcoef(s0,h.obs.axis_score)[0,1]))
print("reproduction n=344:",rep)
# step 2: n = 343
a=load(); sh=pd.read_csv(glob.glob(ROOT+"out/T_shared_genes_*.tsv")[0],sep="\t").gene.tolist()
c1,s1,n1=classify(pd.read_pickle(ROOT+"tpm343.pkl"),a.state,sh)
X=np.asarray(h.X); comp=list(h.var_names)
d=pd.DataFrame(X,index=h.obs_names,columns=comp); d["immune_fraction"]=X[:,2:].sum(1); d["section"]=h.obs.section.values; d["spot_class"]=c1; d["axis_score"]=s1
d["n_types_ge10pct"]=h.obs.n_types_ge10pct.values; d["x"]=np.asarray(h.obsm["spatial"])[:,0]; d["y"]=np.asarray(h.obsm["spatial"])[:,1]; d["spot_class_n344"]=oldcls.values
def jt(groups):
    U=0; ns=[len(x) for x in groups]
    for i in range(3):
        for j in range(i+1,3):
            A_,B_=np.asarray(groups[i]),np.asarray(groups[j]); U+=(A_[:,None]<B_[None,:]).sum()+0.5*(A_[:,None]==B_[None,:]).sum()
    N=sum(ns); mu=(N**2-sum(n**2 for n in ns))/4; var=(N**2*(2*N+3)-sum(n**2*(2*n+3) for n in ns))/72
    return 2*stats.norm.sf(abs((U-mu)/np.sqrt(var)))
rows=[]
for s_,g in d.groupby("section"):
    r=dict(section=s_,n_spots=len(g),genes_used=n1,**{f"n_{SHORT[e]}":int((g.spot_class==e).sum()) for e in ORD},
           **{f"median_immune_{SHORT[e]}":g.loc[g.spot_class==e,"immune_fraction"].median() for e in ORD},
           JT_p_immune_low_int_high=jt([g.loc[g.spot_class==e,"immune_fraction"] for e in ECO]) if all((g.spot_class==e).sum()>0 for e in ECO) else np.nan,
           rho_axis_vs_immune=stats.spearmanr(g.axis_score,g.immune_fraction).statistic,pct_spots_mixed=100*(g.n_types_ge10pct>=2).mean(),
           class_agreement_with_n344=float((g.spot_class==g.spot_class_n344).mean()))
    rows.append(r)
VS=pd.DataFrame(rows); VS.to_csv(ROOT+"out/T_visium_by_section.tsv",sep="\t",index=False); d.to_csv(ROOT+"out/T_visium_spots_n1650.tsv",sep="\t")
pd.DataFrame([rep]).to_csv(ROOT+"out/T_visium_reproduction_n344.tsv",sep="\t",index=False)
pd.set_option("display.width",250); print(VS.round(4).to_string())
