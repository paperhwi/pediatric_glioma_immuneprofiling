"""W2.1 — cell-type classification in GSE102130 by canonical markers."""
import numpy as np, pandas as pd
np.random.seed(42)
D=pd.read_csv("data/GSE102130.txt.gz",sep="\t",index_col=0)
cols=pd.Series(D.columns); pat=cols.str.split("-").str[0]
DERIV=["BCH869_DGC100","BCH869_DGC75","BCH869_GF","BCH869_GS","BCH869_PDX","BCH869_SF"]
keep=(~pat.isin(DERIV+["Oligo"])).values
E=D.loc[:,keep].astype(np.float32)
patient=pd.Series(pat[keep].values,index=E.columns,name="patient")
L=np.log2(E/10.0+1.0)                      # Filbin convention
# drop genes detected in <10 cells
L=L.loc[(L>0).sum(axis=1)>=10]
print("expression matrix:",L.shape,"| tumours:",patient.nunique())

MARK={
 "Myeloid":   ["PTPRC","CSF1R","AIF1","C1QA","C1QB","C1QC","CD14","TYROBP","FCER1G","ITGAM",
               "CX3CR1","P2RY12","CD68","LAPTM5","FCGR3A","MS4A7","CD163"],
 "Lymphoid":  ["CD3D","CD3E","CD3G","CD2","IL7R","CCL5","NKG7","GZMA","GZMK","KLRD1","CD8A","TRAC"],
 "Oligodendrocyte":["MBP","PLP1","MOG","MAG","CNP","TF","CLDN11","MOBP","ERMN","UGT8"],
 "Malignant": ["SOX2","OLIG1","OLIG2","DLL3","PDGFRA","CSPG4","ASCL1","BCAN","PTPRZ1","SOX4","CD24"],
}
for k,v in MARK.items(): MARK[k]=[g for g in v if g in L.index]
# control-gene-matched score (Tirosh 2016): mean(signature) - mean(matched control bin)
mean_expr=L.mean(axis=1)
bins=pd.qcut(mean_expr.rank(method="first"),25,labels=False)
rng=np.random.default_rng(42)
def score(genes,nctrl=50):
    g=[x for x in genes if x in L.index]
    ctrl=[]
    for x in g:
        pool=mean_expr.index[bins==bins[x]]
        ctrl+=list(rng.choice(pool,size=min(nctrl,len(pool)),replace=False))
    return L.loc[g].mean(axis=0)-L.loc[ctrl].mean(axis=0)
S=pd.DataFrame({k:score(v) for k,v in MARK.items()})
Z=(S-S.mean())/S.std()

call=Z.idxmax(axis=1); margin=Z.max(axis=1)-Z.apply(lambda r:r.nlargest(2).iloc[1],axis=1)
call[(Z.max(axis=1)<0.5)|(margin<0.25)]="Unassigned"
# immune calls must actually express PTPRC
if "PTPRC" in L.index:
    imm=call.isin(["Myeloid","Lymphoid"])
    bad=imm & (L.loc["PTPRC"]<=0)
    call[bad]="Unassigned"
    print("immune calls without detectable PTPRC, reassigned:",int(bad.sum()))
res=pd.DataFrame({"patient":patient,"cell_type":call}).join(Z.add_prefix("z_"))
print("\n=== cell-type calls ===")
print(res.cell_type.value_counts().to_string())
print("\n=== by tumour ===")
print(pd.crosstab(res.patient,res.cell_type).to_string())

print("\n=== verification: mean log2 expression of canonical markers by call ===")
chk=["PTPRC","CSF1R","C1QB","P2RY12","CD3E","CD2","NKG7","MBP","PLP1","OLIG1","PDGFRA","SOX2"]
chk=[g for g in chk if g in L.index]
v=pd.DataFrame({g:L.loc[g].groupby(call).mean() for g in chk}).round(2)
print(v.to_string())
res.to_csv("w2_cell_calls.tsv",sep="\t")
L.to_pickle("w2_logtpm.pkl"); patient.to_pickle("w2_patient.pkl")
print("\nsaved")
