import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
import statsmodels.formula.api as smf, statsmodels.api as sm
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess
a=load(); P=pd.read_pickle(ROOT+"programs343.pkl"); pc1=pd.read_csv(ROOT+"out/immune_axis_n343.tsv",sep="\t",index_col=0).immune_axis
R="/home/claude/repo/cancers_submission_2026/reproducibility_data/"
F=pd.read_csv(R+"feature_matrices/main_cohort_LM22_ssGSEA_z.tsv",sep="\t",index_col=0).loc[a.index]
D=P.copy()
# derived balance indices
D["Effector:suppressor balance"]=P["Effector CD8 / cytotoxic"]-P[["Treg","T-cell exhaustion / inhibitory receptors","Immunosuppressive TAM"]].mean(1)
D["Cytotoxicity not matched by exhaustion (residual)"]=sm.OLS(P["CYT (GZMA, PRF1)"],sm.add_constant(P["T-cell exhaustion / inhibitory receptors"])).fit().resid
D["CD8:Treg (LM22, log ratio)"]=np.nan  # placeholder replaced below
lm=pd.read_csv(R+"feature_matrices/main_cohort_LM22_ssGSEA_z.tsv",sep="\t",index_col=0).loc[a.index]
D["CD8:Treg (LM22, log ratio)"]=lm["LM22_T cells CD8"]-lm["LM22_T cells regulatory (Tregs)"]
D=(D-D.mean())/D.std()
D["state"]=a.state; D["group"]=a.cohort_group; D["pc1"]=pc1
ordc={LO:0,MID:1,HI:2}
rows=[]
progs=[c for c in D.columns if c not in ("state","group","pc1")]
cat=pd.read_csv(ROOT+"out/immune_programme_definitions.tsv",sep="\t").set_index("program").category.to_dict()
for p in progs:
    y=D[p]; g=[y[D.state==e] for e in ECO]
    kw=stats.kruskal(*g).pvalue
    df=pd.DataFrame({"y":y,"eco":pd.Categorical(D.state,ECO),"grp":D.group,"pc1":D.pc1})
    m=smf.ols("y~C(eco)+C(grp)",df).fit()
    li=m.params["C(eco)[T.Immune-high]"]; md=m.params["C(eco)[T.Immune-intermediate]"]
    pli=m.pvalues["C(eco)[T.Immune-high]"]; pmd=m.pvalues["C(eco)[T.Immune-intermediate]"]
    # LI vs MD
    m2=smf.ols("y~C(eco,Treatment('Immune-intermediate'))+C(grp)",df).fit(); plm=m2.pvalues.filter(like="Immune-high").iloc[0]
    rho,prho=stats.spearmanr(df.y,df.pc1)
    # continuum test: does ecotype add beyond a smooth function of PC1?
    base=smf.ols("y~bs(pc1,df=4)+C(grp)",df).fit(); full=smf.ols("y~bs(pc1,df=4)+C(grp)+C(eco)",df).fit()
    Fp=full.compare_f_test(base)[1]
    trend=stats.spearmanr(df.y,D.state.map(ordc))
    rows.append(dict(program=p,category=cat.get(p,"Balance index"),
        mean_low=y[D.state==LO].mean(),mean_intermediate=y[D.state==MID].mean(),mean_high=y[D.state==HI].mean(),
        KW_p=kw,adj_diff_intermediate_vs_low=md,p_intermediate_vs_low=pmd,adj_diff_high_vs_low=li,p_high_vs_low=pli,p_high_vs_intermediate=plm,
        trend_rho_low_int_high=trend.correlation,trend_p=trend.pvalue,
        rho_with_immune_axis=rho,R2_spline_PC1=base.rsquared,delta_R2_state_beyond_axis=full.rsquared-base.rsquared,p_state_beyond_axis=Fp))
T=pd.DataFrame(rows)
for c in ["KW_p","p_intermediate_vs_low","p_high_vs_low","p_high_vs_intermediate","trend_p","p_state_beyond_axis"]: T[c.replace("p","q",1) if c.startswith("p") else c+"_q"]=bh(T[c])
T["fraction_of_low_to_high_change_reached_at_intermediate"]=(T.mean_intermediate-T.mean_low)/(T.mean_high-T.mean_low)
T.to_csv(ROOT+"out/q2/Q2_programmes_by_state.tsv",sep="\t",index=False)
pd.set_option("display.width",250)
print(T[["program","mean_low","mean_intermediate","mean_high","KW_p","p_high_vs_intermediate","rho_with_immune_axis","R2_spline_PC1","delta_R2_state_beyond_axis","p_state_beyond_axis"]].round(3).to_string())
# ---- ordering along the axis: PC1 position where each program's LOESS crosses the cohort mean (0)
def onset(x,y):
    f=lowess(y,x,frac=0.5,return_sorted=True); xs,ys=f[:,0],f[:,1]
    s=np.where((ys[:-1]<0)&(ys[1:]>=0))[0]
    if len(s)==0: return np.nan
    i=s[0]; return xs[i]+(0-ys[i])*(xs[i+1]-xs[i])/(ys[i+1]-ys[i])
rng=np.random.default_rng(20261101); x=D.pc1.values
ON=[]
for p in progs:
    y=D[p].values; o=onset(x,y); bs=[]
    for b in range(300):
        i=rng.integers(0,len(x),len(x)); bs.append(onset(x[i],y[i]))
    bs=np.array(bs); ON.append(dict(program=p,onset_PC1=o,lo=np.nanpercentile(bs,2.5),hi=np.nanpercentile(bs,97.5),boot_defined=np.isfinite(bs).mean()))
O=pd.DataFrame(ON).sort_values("onset_PC1"); O.to_csv(ROOT+"out/q2/Q2_program_onset_along_PC1.tsv",sep="\t",index=False); print(O.round(2).to_string())
# ecotype PC1 ranges
print(D.groupby("state").pc1.describe().round(2))
D.to_pickle("q2_scores.pkl")
