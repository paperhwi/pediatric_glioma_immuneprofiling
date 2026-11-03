"""Published immune programmes (Ayers GEP, Rooney CYT, Wang APS, Garg ICD + DAMP subset) and 21 checkpoint genes (n = 343):
scores as in Rev344_06 supp.py S30; Cox per SD (molecular group + age + sex) without and with the immune axis; state ordering;
group differences given immune axis and state (OLS)."""
import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
import statsmodels.formula.api as smf
from scipy import stats
from published_sets import SETS, CHECKPOINT_GENES
a=load(); s=a[a.os_eval]; pc1=pd.read_csv(ROOT+"out/immune_axis_n343.tsv",sep="\t",index_col=0).immune_axis
tpm=pd.read_pickle(ROOT+"tpm343.pkl").astype(float); x=np.log2(tpm[list(a.index)]+1); z=x.sub(x.mean(1),axis=0).div(x.std(1)+1e-9,axis=0)
cov,sc=[],{}
for k,v in SETS.items():
    g=[q for q in v["genes"] if q in x.index]
    cov.append(dict(programme=k,n_listed=len(v["genes"]),n_present=len(g),genes=",".join(v["genes"]),reference=v["ref"])); sc[k]=x.loc[g].mean() if k.startswith("Cytolytic") else z.loc[g].mean()
PS=pd.DataFrame(sc); PS=(PS-PS.mean())/PS.std()
for g in CHECKPOINT_GENES: PS[f"gene:{g}"]=z.loc[g] if g in z.index else np.nan
rows=[]; eff=[]
for k in PS.columns:
    for adj in (False,True):
        ex=PS[[k]].rename(columns={k:"score"});
        if adj: ex["immune_axis"]=pc1
        r=cox(design(s,extra=ex,eco=False)).summary.loc["score"]
        rows.append(dict(feature=k,model="+ immune axis" if adj else "base",HR_per_SD=r["exp(coef)"],lo=r["exp(coef) lower 95%"],hi=r["exp(coef) upper 95%"],p=r["p"]))
    y=PS[k]; g3=[y[a.state==e] for e in ECO]
    q=pd.DataFrame({"y":y,"ax":pc1,"st":a.state,"g":a.cohort_group}); f1=smf.ols("y ~ ax + C(st, Treatment('Immune-low')) + C(g, Treatment('pHGG_WT'))",q).fit()
    e={"feature":k,**{f"mean_{SHORT[ee]}":y[a.state==ee].mean() for ee in ECO},"spearman_axis":stats.spearmanr(y,pc1).statistic,"KW_p":stats.kruskal(*g3).pvalue}
    for gg in ["DMG_K27","DHG_G34","IHG"]:
        t=f"C(g, Treatment('pHGG_WT'))[T.{gg}]"; e[f"{gg}_vs_pHGG-WT_given_axis_state"]=f1.params[t]; e[f"{gg}_p"]=f1.pvalues[t]
    eff.append(e)
C=pd.DataFrame(rows); C["q_BH"]=C.groupby("model").p.transform(bh); E=pd.DataFrame(eff)
for gg in ["DMG_K27","DHG_G34","IHG"]: E[f"{gg}_q"]=bh(E[f"{gg}_p"])
pd.DataFrame(cov).to_csv(ROOT+"out/T_published_programme_definitions.tsv",sep="\t",index=False); C.to_csv(ROOT+"out/T_published_programmes_cox.tsv",sep="\t",index=False); E.to_csv(ROOT+"out/T_published_programmes_by_state_group.tsv",sep="\t",index=False)
pd.set_option("display.width",250); print(C[~C.feature.str.startswith("gene:")].round(3).to_string()); print(E[~E.feature.str.startswith("gene:")][["feature","mean_Low","mean_Intermediate","mean_High"]].round(2))
print(E.loc[(E[[f"{g}_q" for g in ["DMG_K27","DHG_G34","IHG"]]]<0.1).any(axis=1),["feature"]+[c for c in E if c.endswith("_q") or "given" in c]].round(3).to_string())
