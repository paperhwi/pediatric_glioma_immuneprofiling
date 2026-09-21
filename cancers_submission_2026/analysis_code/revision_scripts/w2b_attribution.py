"""A — do the brain-tuned ssGSEA signatures attribute to the cell types they are named for?
10x immune atlas, paediatric high-grade glioma samples only. Per-cell, composition-independent."""
import numpy as np, pandas as pd
from scipy import stats
np.random.seed(42)
UP="/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs"
sigs={}
for line in open(f"{UP}/brain_immune_signatures.gmt"):
    f=line.rstrip("\n").split("\t"); sigs[f[0]]=[g for g in f[2:] if g.strip()]

M=pd.read_pickle("w2_immune_tenx.pkl")
meta=pd.read_csv("data/immune_tenx_meta.csv.gz",index_col=0)
meta=meta.loc[meta.index.intersection(M.columns)]
PHGG=["H3K27M","H3WT_hemispheric"]          # exclude PF-A ependymoma
meta=meta[meta.Subtype.isin(PHGG)]
M=M[meta.index]
lab=meta.detailed_annot.replace({"unclear":np.nan}).dropna()
M=M[lab.index]; meta=meta.loc[lab.index]
print("pHGG immune cells:",M.shape[1],"| samples:",meta.sampleid.nunique())
print(lab.value_counts().to_string())

L=np.log2(M/10.0+1.0)
# verification of the published labels against canonical markers
print("\n=== verification: mean log2 expression by published label ===")
chk=[g for g in ["PTPRC","CD3D","CD2","CD8A","CD4","CSF1R","C1QB","P2RY12","CD68","TYROBP"] if g in L.index]
print(pd.DataFrame({g:L.loc[g].groupby(lab).mean() for g in chk}).round(2).to_string())

# control-gene-matched signature score (Tirosh 2016)
mean_expr=L.mean(axis=1); bins=pd.qcut(mean_expr.rank(method="first"),25,labels=False)
rng=np.random.default_rng(42)
def score(genes,nctrl=50):
    g=[x for x in genes if x in L.index]
    if len(g)<3: return None,len(g)
    ctrl=[]
    for x in g:
        pool=mean_expr.index[bins==bins[x]]
        ctrl+=list(rng.choice(pool,size=min(nctrl,len(pool)),replace=False))
    return L.loc[g].mean(axis=0)-L.loc[ctrl].mean(axis=0),len(g)

rows=[];S={}
for name,genes in sigs.items():
    sc,n=score(genes)
    if sc is None: print(f"  skipped {name}: only {n} genes present"); continue
    S[name]=sc
    grp={k:sc[lab==k].values for k in ["Myeloid","CD4","CD8"]}
    H,p=stats.kruskal(*grp.values())
    best=max(grp,key=lambda k:np.median(grp[k]))
    rows.append(dict(signature=name,n_genes=n,top_cell_type=best,KW_H=round(H,1),p=p,
                     **{f"median_{k}":round(float(np.median(v)),3) for k,v in grp.items()}))
S=pd.DataFrame(S)
res=pd.DataFrame(rows); res["q_BH"]=stats.false_discovery_control(res.p)
res=res.sort_values("signature")
print("\n=== signature attribution (median control-matched score per cell type) ===")
print(res[["signature","n_genes","top_cell_type","median_Myeloid","median_CD4","median_CD8","q_BH"]]
      .to_string(index=False,float_format=lambda x:f"{x:.3g}"))
S.assign(cell_type=lab,sampleid=meta.sampleid).to_pickle("w2_sigscores_10x.pkl")
res.to_csv("W2_A_signature_attribution.tsv",sep="\t",index=False)

# expected-vs-observed table
EXPECT={"Microglia_Core_Homeostatic":"Myeloid","Microglia_Klemm2020":"Myeloid","MDM_Klemm2020":"Myeloid",
 "MgTAM_Antunes2021":"Myeloid","MoTAM_Antunes2021":"Myeloid","DAM_KerenShaul2017":"Myeloid",
 "M1_Macrophage":"Myeloid","M2_Macrophage":"Myeloid","Dendritic_Cell_Activation":"Myeloid",
 "Neutrophil_Activation":"Myeloid","Glioma_Inflammatory_Wang2017":"Myeloid","MHC_Class_II":"Myeloid",
 "T_Cell_Cytotoxicity":"T cell (CD8)","T_Cell_Exhaustion":"T cell","Tregs_Friebel2020":"T cell (CD4)",
 "NK_Cell_Activity":"T cell (CD8)","Chemokine_T_Cell_Recruitment":"either","MHC_Class_I":"any",
 "IFN_Gamma_Response":"any","IFN_Alpha_Response":"any","TGFb_Immunosuppression":"any",
 "MAPK_Activity":"any","Cell_Cycle_Proliferation":"any","Stemness_Brain_Tumor":"none (tumour)"}
res["expected"]=res.signature.map(EXPECT)
def concord(r):
    e=r["expected"]; o=r["top_cell_type"]
    if e in ("any","either","none (tumour)"): return "not directional"
    if e=="Myeloid": return "concordant" if o=="Myeloid" else "DISCORDANT"
    if e.startswith("T cell"):
        if "(" in e: return "concordant" if o==e.split("(")[1][:-1] else ("partly (T cell, wrong subset)" if o in ("CD4","CD8") else "DISCORDANT")
        return "concordant" if o in ("CD4","CD8") else "DISCORDANT"
    return "?"
res["verdict"]=res.apply(concord,axis=1)
res.to_csv("W2_A_signature_attribution.tsv",sep="\t",index=False)
print("\n=== concordance with the name of each signature ===")
print(res.verdict.value_counts().to_string())
print("\ndirectional signatures:")
print(res[res.verdict!="not directional"][["signature","expected","top_cell_type","verdict"]].to_string(index=False))
