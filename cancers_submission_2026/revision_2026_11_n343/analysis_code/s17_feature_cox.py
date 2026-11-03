import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
from lifelines.statistics import proportional_hazard_test
R="/home/claude/repo/cancers_submission_2026/reproducibility_data/"
a=load(); M=pd.read_pickle(ROOT+"markers343.pkl"); P=pd.read_pickle(ROOT+"programs343.pkl")
F=pd.read_csv(R+"feature_matrices/main_cohort_LM22_ssGSEA_z.tsv",sep="\t",index_col=0).loc[a.index]
from sklearn.decomposition import PCA
# Immune-presence axis = PC1 of the 21 immune ssGSEA signatures (excludes the 22 compositional LM22 fractions and the
# 3 tumour-intrinsic signatures: cell cycle, stemness, MAPK). See pc1fix.py for why the 46-feature PC1 explains only ~32%.
NONIMM=["ssGSEA_Cell_Cycle_Proliferation","ssGSEA_Stemness_Brain_Tumor","ssGSEA_MAPK_Activity"]
SSI=[c for c in F.columns if c.startswith("ssGSEA") and c not in NONIMM]
Xi=(F[SSI]-F[SSI].mean())/F[SSI].std()
pca=PCA(5).fit(Xi.values); pc1=pd.Series(pca.transform(Xi.values)[:,0],index=F.index)
if pc1.corr(F["ssGSEA_T_Cell_Cytotoxicity"])<0: pc1=-pc1
pc1=(pc1-pc1.mean())/pc1.std(); a["immune_axis"]=pc1; print("immune PC1 (21 ssGSEA immune signatures) variance explained",pca.explained_variance_ratio_[0].round(3))

ar=load(os_version="released")
DATA={"primary":a[a.os_eval].copy(),"as_released":ar[ar.os_eval].copy()}
def fit_term(s,X,term,eco=False,adj_pc1=False):
    ex=X[[term]].copy()
    if adj_pc1: ex["immune_axis"]=pc1
    d=design(s,extra=ex,eco=eco)
    if d[term].nunique()<2: return None
    try: c=cox(d)
    except Exception: c=cox(d,penalizer=0.01)
    r=c.summary.loc[term]
    try: php=float(proportional_hazard_test(c,d,time_transform="rank").summary.loc[term,"p"])
    except Exception: php=np.nan
    return dict(n=len(d),events=int(d.event.sum()),HR=r["exp(coef)"],lo=r["exp(coef) lower 95%"],hi=r["exp(coef) upper 95%"],p=r["p"],PH_p=php)
res=[]
for ds,s in DATA.items():
    # ---- markers
    mk=[c for c in M.columns if c not in ("has_DNA","has_CNV","TMB_coding","Tumor_purity_RFpurify_ESTIMATE")]
    for c in mk:
        x=M.loc[s.index,c]; 
        if c!="log2_TMB":
            nalt=int((x==1).sum()); ev=int(s.loc[x==1,"event"].sum())
            if nalt<10 or ev<5: continue
        else: nalt=int(x.notna().sum())
        for eco in (False,True):
            r=fit_term(s,M,c,eco=eco)
            if r: res.append(dict(dataset=ds,family="Molecular marker",feature=c,model="+ immune state" if eco else "base",n_altered=nalt,**r))
    # ---- immune features & programs (per SD)
    for fam,X in [("Deconvolution / brain signature (46)",F),("Curated immune program (17)",P)]:
        Xs=(X-X.mean())/X.std()
        for c in X.columns:
            for adj in (False,True):
                r=fit_term(s,Xs,c,adj_pc1=adj)
                if r: res.append(dict(dataset=ds,family=fam,feature=c,model="+ immune axis" if adj else "base",n_altered=np.nan,**r))
T=pd.DataFrame(res)
T["q_BH_within_family"]=T.groupby(["dataset","family","model"]).p.transform(lambda p: bh(p))
T.to_csv(ROOT+"out/q1/Q1_cox_markers_and_immune_features.tsv",sep="\t",index=False)
pd.set_option("display.width",220)
v=T[(T.dataset=="primary")]
print(v[v.family=="Molecular marker"].pivot_table(index=["feature","n_altered"],columns="model",values=["HR","p","q_BH_within_family"]).round(3).to_string())
for fam in ["Deconvolution / brain signature (46)","Curated immune program (17)"]:
    w=v[v.family==fam].pivot_table(index="feature",columns="model",values=["HR","p","q_BH_within_family"]).round(3)
    print(fam); print(w.sort_values(("p","base")).to_string())
