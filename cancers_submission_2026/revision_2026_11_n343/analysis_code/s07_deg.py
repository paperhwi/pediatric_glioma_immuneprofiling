"""Figure 6 differential expression (molecular-group-adjusted OLS on log2(TPM+1); model and gene filter of
core_pipeline/final_260722_rerun.py::adjusted_deg). Step 1 reproduces the submitted n = 349 table; step 2 fits n = 343."""
import sys, json; sys.path.insert(0,"/home/claude/q343/src"); from common import *
from scipy import stats
I="/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs/"
def filtered_tpm(tpm, ids):
    x=tpm[ids]; x=x.loc[(x>=1).sum(axis=1)>=int(np.ceil(0.05*len(ids)))]; g=pd.Series(x.index.astype(str),index=x.index)
    keep=g.str.match(r"^RPS6K")
    ex=(g.str.match(r"^MT-")|g.str.match(r"^MTRNR")|g.str.match(r"^MT(ATP|CO|CYB|ND|RNR)\d")|(g.str.match(r"^RP[LS]\d")&~keep)|g.str.match(r"^RPLP\d")|g.str.match(r"^RPSA$")
        |g.str.match(r"^RNA(5S|5-8S|18S|28S|45S|5-8SN|45SN)")|g.str.contains(r"rRNA")|g.str.match(r"^RN(U|VU)\d")|g.str.fullmatch(r"U\d+")|g.str.fullmatch(r"7SK|Y_RNA|Vault|Metazoa_SRP|uc_338")
        |g.str.match(r"^RN7S[KL]")|g.str.match(r"^(SNORD|SNORA|SCARNA|SNAR)")|g.str.match(r"^RNY\d")|g.str.match(r"^MIR\d")|g.str.match(r"^MIRLET")|g.str.match(r"^VTRNA")|g.str.fullmatch(r"RPPH1|RMRP")|g.str.match(r"^ENSG\d"))
    return x.loc[~ex]
def fit(tpm,meta,levels,contrasts):
    ids=list(meta.index); x=filtered_tpm(tpm,ids); Y=np.log2(x.T.to_numpy(dtype=np.float64)+1.0)
    eco=pd.Categorical(meta.lab,categories=levels); coh=pd.Categorical(meta.cohort_group,categories=["pHGG_WT","DMG_K27","DHG_G34","IHG"])
    Xdf=pd.get_dummies(pd.DataFrame({"lab":eco,"cohort_group":coh}),drop_first=True).astype(float); Xdf.insert(0,"Intercept",1.0)
    X=Xdf.to_numpy(); inv=np.linalg.pinv(X.T@X); beta=inv@X.T@Y; resid=Y-X@beta; dfr=X.shape[0]-np.linalg.matrix_rank(X); s2=(resid**2).sum(0)/dfr
    rows=[]; summ={"n_samples":len(ids),"n_genes_tested":len(x),"df_resid":int(dfr),"contrasts":{}}
    for name,wts in contrasts.items():
        c=np.zeros(X.shape[1])
        for t,wv in wts.items(): c[Xdf.columns.get_loc(t)]=wv
        eff=c@beta; se=np.sqrt(np.maximum(s2*(c@inv@c),1e-300)); tv=eff/se; p=2*stats.t.sf(np.abs(tv),dfr); q=bh(p)
        rows.append(pd.DataFrame({"gene":x.index,"adjusted_log2FC":eff,"t":tv,"p":p,"q_BH":q,"df_resid":dfr,"contrast":name}))
        summ["contrasts"][name]={"up":int(((q<0.05)&(eff>1)).sum()),"down":int(((q<0.05)&(eff<-1)).sum())}
    return pd.concat(rows,ignore_index=True),summ
old=pd.read_csv(I+"cohort_main_final.tsv",sep="\t",low_memory=False).merge(pd.read_csv(I+"ecotype_LM22_main_k3_annotated.tsv",sep="\t"),on="Kids_First_Biospecimen_ID").set_index("Kids_First_Biospecimen_ID")
old["lab"]=old.ecotype
a=load(); cur=a[["state","cohort_group"]].rename(columns={"state":"lab"})
tpm=pd.read_csv(I+"tpm_for_cibersortx.tsv",sep="\t",usecols=["GeneSymbol",*old.index]).set_index("GeneSymbol"); tpm=tpm.loc[~tpm.index.duplicated(keep="first")]
sub,ss=fit(tpm,old,["Immune-desert","Lymphocyte-inflamed","Myeloid-dominant"],{"Lymphocyte_inflamed_vs_Immune_desert":{"lab_Lymphocyte-inflamed":1.0},"Myeloid_dominant_vs_Immune_desert":{"lab_Myeloid-dominant":1.0},"Lymphocyte_inflamed_vs_Myeloid_dominant":{"lab_Lymphocyte-inflamed":1.0,"lab_Myeloid-dominant":-1.0}})
src=pd.read_csv(I+"adjusted_DEG_all_contrasts.tsv",sep="\t"); m=sub.merge(src,on=["gene","contrast"],suffixes=("_rerun","_sub"))
val={"rows_matched":len(m),"source_rows":len(src),"max_abs_effect_diff":float((m.adjusted_log2FC_rerun-m.adjusted_log2FC_sub).abs().max()),"max_abs_q_diff":float((m.q_BH_rerun-m.q_BH_sub).abs().max())}
assert len(m)==len(src) and val["max_abs_effect_diff"]<1e-8, val
C={"Immune-high vs Immune-low":{"lab_Immune-high":1.0},"Immune-intermediate vs Immune-low":{"lab_Immune-intermediate":1.0},"Immune-high vs Immune-intermediate":{"lab_Immune-high":1.0,"lab_Immune-intermediate":-1.0}}
rev,rs=fit(tpm,cur,[LO,HI,MID],C); rev.to_csv(ROOT+"out/T_adjusted_DEG_n343_all_contrasts.tsv",sep="\t",index=False)
json.dump({"submitted_n349_reproduction":val,"revised_n343":rs},open(ROOT+"out/T_DEG_summary.json","w"),indent=1); print(json.dumps({"val":val,"rev":rs},indent=1))
