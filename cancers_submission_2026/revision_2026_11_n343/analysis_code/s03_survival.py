"""Survival analyses on n = 343 (corrected OS). Primary model fixed in advance = the submitted model:
immune state + molecular group + age + sex. All sensitivity analyses are reported regardless of significance.
Continuous-axis RNA-library-controlled follow-up fits are in ../additional_analyses_20261103/code/axis_library_survival.py;
aggregate results are kept separately so the original primary T_immune_axis_cox.tsv remains identifiable."""
import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
from lifelines.statistics import multivariate_logrank_test, proportional_hazard_test, logrank_test
from scipy.stats import chi2
a=load(); s=a[a.os_eval].copy()
pc1=pd.read_csv(ROOT+"out/immune_axis_n343.tsv",sep="\t",index_col=0).immune_axis
TERMS={HI:"Immune-high vs Immune-low",MID:"Immune-intermediate vs Immune-low"}
def summ(c,d,label,lr_df=None,note=""):
    r=c.summary; out={"analysis":label,"n":len(d),"deaths":int(d.event.sum())}
    for k,nm in [(HI,"high"),(MID,"intermediate")]:
        out.update({f"{nm}_HR":r.loc[k,"exp(coef)"],f"{nm}_lo":r.loc[k,"exp(coef) lower 95%"],f"{nm}_hi":r.loc[k,"exp(coef) upper 95%"],f"{nm}_p":r.loc[k,"p"]})
    if lr_df is not None: out["logrank_p"]=multivariate_logrank_test(lr_df.time,lr_df.state,lr_df.event).p_value
    out["note"]=note; return out
# ---------- primary
d=design(s); c=cox(d); P0=c.summary
P0.rename(index={**TERMS,"DMG_K27":"DMG vs pHGG-WT","DHG_G34":"DHG vs pHGG-WT","IHG":"IHG vs pHGG-WT","age_years":"Age at diagnosis (per year)","male":"Male vs female"}).to_csv(ROOT+"out/T_primary_cox.tsv",sep="\t")
ph=proportional_hazard_test(c,d,time_transform="rank").summary; ph.to_csv(ROOT+"out/T_primary_cox_PH_test.tsv",sep="\t")
lr_all=multivariate_logrank_test(s.time,s.state,s.event).p_value
d0=design(s,eco=False); c0=cox(d0); lrt_state=chi2.sf(2*(c.log_likelihood_-c0.log_likelihood_),2)
print("PRIMARY n",len(d),"deaths",int(d.event.sum()),"logrank",round(lr_all,4),"LRT state",round(lrt_state,4)); print(P0[["exp(coef)","exp(coef) lower 95%","exp(coef) upper 95%","p"]].round(4))
rows=[summ(c,d,"Primary (pre-specified): corrected OS, molecular group + age + sex",s,"submitted model applied to corrected data")]
# ---------- the four sensitivity analyses requested
na=s[~s.autopsy]; dn=design(na); rows.append(summ(cox(dn),dn,"S1 Excluding autopsy specimens",na))
sp=s.copy(); sp["time"]=s.time-s.sampling_delay_days; sp=sp[sp.time>0]
ds=design(sp); rows.append(summ(cox(ds),ds,"S2 Time origin at specimen collection (OS from sampling)",sp,"autopsy specimens (time 0) cannot contribute"))
dA=design(s,extra=pd.DataFrame({"autopsy_specimen":s.autopsy.astype(float)})); rows.append(summ(cox(dA),dA,"S3 Adjusted for autopsy specimen",s))
lib=pd.get_dummies(s.RNA_library).astype(float); ref=lib.sum().idxmax(); lib=lib.drop(columns=ref); lib.columns=[f"library_{x}" for x in lib.columns]
dL=design(s,extra=lib); rows.append(summ(cox(dL),dL,f"S4 Adjusted for RNA library (reference: {ref})",s))
dLs=design(s,extra=pd.DataFrame({"lib":s.RNA_library.astype("category").cat.codes})); rows.append(summ(cox(dLs,strata=["lib"]),dLs,"S4b Stratified by RNA library",s))
# ---------- further pre-existing sensitivity analyses
sl=s.copy(); sl=sl[sl.sampling_delay_days<sl.time]; de=design(sl); de["entry"]=sl.sampling_delay_days.clip(lower=0).loc[de.index]
rows.append(summ(cox(de,entry_col="entry"),de,"S5 Delayed entry at specimen date (left truncation)",sl))
r0=load(os_version="released"); r0=r0[r0.os_eval]; dr=design(r0); rows.append(summ(cox(dr),dr,"S6 OS as released (no OS-origin corrections)",r0))
cl=load(labels="n344_carried"); cl=cl[cl.os_eval]; dc=design(cl); rows.append(summ(cox(dc),dc,"S7 Immune-state labels carried over from n = 344",cl))
for cap in [5,10]:
    x=s.copy(); o=x.time>cap*365.25; x.loc[o,"event"]=0; x.loc[o,"time"]=cap*365.25; dx=design(x); rows.append(summ(cox(dx),dx,f"S8 Administrative censoring at {cap} y",x))
w=s[s.sampling_delay_days<=365]; dw=design(w); rows.append(summ(cox(dw),dw,"S9 Specimen within 1 y of diagnosis",w))
au=pd.read_csv(ROOT+"out/T_OS_origin_audit_all_n343.tsv",sep="\t"); aud=au.loc[au.category.str.startswith(("R2","R3")),"Kids_First_Participant_ID"]
x=s[~s.Kids_First_Participant_ID.isin(aud)]; dx=design(x); rows.append(summ(cox(dx),dx,"S10 Excluding all 13 records whose OS origin precedes the sequenced-HGG diagnosis",x))
dnoadj=design(s,grp=False,age=False,sex=False); rows.append(summ(cox(dnoadj),dnoadj,"S11 Unadjusted Cox",s))
SENS=pd.DataFrame(rows); SENS.to_csv(ROOT+"out/T_survival_sensitivity.tsv",sep="\t",index=False)
pd.set_option("display.width",250); print(SENS.drop(columns="note").round(3).to_string())
# ---------- autopsy table
def ap(df,lab): 
    t=df.groupby("state").agg(n=("autopsy","size"),autopsy=("autopsy","sum")).reindex(ORD); t["autopsy_%"]=100*t.autopsy/t.n; t.loc["Total"]=[t.n.sum(),t.autopsy.sum(),100*t.autopsy.sum()/t.n.sum()]; t.insert(0,"set",lab); return t
AP=pd.concat([ap(a,"Cohort (n = 343)"),ap(s,"Survival set (n = %d)"%len(s))]); AP.index.name="state"
AP2=s.groupby("state").apply(lambda g: pd.Series({"deaths":int(g.event.sum()),"deaths_non_autopsy":int(g.event[~g.autopsy].sum()),"n_non_autopsy":int((~g.autopsy).sum())})).reindex(ORD)
AP2.loc["Total"]=AP2.sum(); AP=AP.reset_index().merge(AP2.reset_index(),on="state",how="left"); AP.loc[AP.set.str.startswith("Cohort"),["deaths","deaths_non_autopsy","n_non_autopsy"]]=np.nan
from scipy.stats import fisher_exact, chi2_contingency
AP["test"]=""; AP.loc[0,"test"]="autopsy by state, chi-square P = %.3g"%chi2_contingency(pd.crosstab(a.state,a.autopsy))[1]
AP.loc[4,"test"]="autopsy by state, chi-square P = %.3g"%chi2_contingency(pd.crosstab(s.state,s.autopsy))[1]
AP.to_csv(ROOT+"out/T_autopsy_by_state.tsv",sep="\t",index=False); print(AP.round(1).to_string())
# ---------- follow-up
def fu(g):
    rk=KaplanMeierFitter().fit(g.time/365.25,1-g.event); km=KaplanMeierFitter().fit(g.time/365.25,g.event); alive=g[g.event==0]; y=g.time/365.25
    from lifelines.utils import median_survival_times
    ci=median_survival_times(rk.confidence_interval_)
    r={"n":len(g),"deaths":int(g.event.sum()),"censored":int((g.event==0).sum()),"reverse_KM_median_FU_y":rk.median_survival_time_,
       "reverse_KM_FU_95CI":f"{ci.iloc[0,0]:.2f}-{ci.iloc[0,1]:.2f}",
       "survivors_median_FU_y":(alive.time/365.25).median() if len(alive) else np.nan,
       "survivors_FU_IQR_y":f"{(alive.time/365.25).quantile(.25):.2f}-{(alive.time/365.25).quantile(.75):.2f}" if len(alive) else "",
       "survivors_FU_range_y":f"{(alive.time/365.25).min():.2f}-{(alive.time/365.25).max():.2f}" if len(alive) else "",
       "observation_time_median_y":y.median(),"observation_time_range_y":f"{y.min():.2f}-{y.max():.2f}","median_OS_y":km.median_survival_time_,
       "OS_2y":float(km.predict(2)),"OS_5y":float(km.predict(5))}
    for t in [2,5,10]: r[f"at_risk_{t}y"]=int((y>=t).sum())
    return pd.Series(r)
FU=pd.concat({"All":fu(s),**{e:fu(s[s.state==e]) for e in ORD}},axis=1).T; FU.to_csv(ROOT+"out/T_followup.tsv",sep="\t"); print(FU.T.to_string())
# ---------- continuous immune axis
ex=pd.DataFrame({"immune_axis":pc1}); dp=design(s,extra=ex,eco=False); cp=cox(dp); rp=cp.summary.loc["immune_axis"]
d2=design(s,extra=ex); c2=cox(d2); lrt_state_beyond_axis=chi2.sf(2*(c2.log_likelihood_-cp.log_likelihood_),2)
lrt_axis=chi2.sf(2*(cp.log_likelihood_-c0.log_likelihood_),1)
AX=pd.DataFrame([dict(model="Immune axis (per SD) + molecular group + age + sex",n=len(dp),deaths=int(dp.event.sum()),HR=rp["exp(coef)"],lo=rp["exp(coef) lower 95%"],hi=rp["exp(coef) upper 95%"],p=rp["p"],
    LRT_axis_vs_base_p=lrt_axis,AIC_delta=cp.AIC_partial_-c0.AIC_partial_,LRT_state_beyond_axis_p=lrt_state_beyond_axis,PH_p=float(proportional_hazard_test(cp,dp,time_transform="rank").summary.loc["immune_axis","p"]))])
for lab_,sub in [("Excluding autopsy specimens",s[~s.autopsy])]:
    dd=design(sub,extra=ex,eco=False); cc=cox(dd).summary.loc["immune_axis"]
    AX.loc[len(AX)]=dict(model=f"Immune axis, {lab_}",n=len(dd),deaths=int(dd.event.sum()),HR=cc["exp(coef)"],lo=cc["exp(coef) lower 95%"],hi=cc["exp(coef) upper 95%"],p=cc["p"])
AX.to_csv(ROOT+"out/T_immune_axis_cox.tsv",sep="\t",index=False); print(AX.round(4).to_string())
# ---------- within subtype
E=[];SUM=[]
for g in ["pHGG_WT","DMG_K27","DHG_G34","IHG"]:
    se=s[s.cohort_group==g]; row=dict(group=g,n_OS=len(se),deaths=int(se.event.sum()),**{f"n_{SHORT[e]}":int((se.state==e).sum()) for e in ORD})
    if len(se)>=30:
        row["logrank_p"]=multivariate_logrank_test(se.time,se.state,se.event).p_value
        dg=design(se,grp=False); r=cox(dg).summary
        for k in (HI,MID): E.append(dict(group=g,term=TERMS[k],n=len(dg),events=int(dg.event.sum()),HR=r.loc[k,"exp(coef)"],lo=r.loc[k,"exp(coef) lower 95%"],hi=r.loc[k,"exp(coef) upper 95%"],p=r.loc[k,"p"]))
    SUM.append(row)
E=pd.DataFrame(E); E.to_csv(ROOT+"out/T_within_subtype_state_cox.tsv",sep="\t",index=False); pd.DataFrame(SUM).to_csv(ROOT+"out/T_within_subtype_survival_summary.tsv",sep="\t",index=False)
x=s[s.cohort_group.isin(["DMG_K27","pHGG_WT"])]; di=design(x,grp=False); di["DMG"]=(x.cohort_group=="DMG_K27").astype(float).loc[di.index]
ci0=cox(di); di2=di.copy(); di2["HIxDMG"]=di[HI]*di.DMG; di2["MIDxDMG"]=di[MID]*di.DMG; ci1=cox(di2)
pint=chi2.sf(2*(ci1.log_likelihood_-ci0.log_likelihood_),2)
pd.DataFrame([dict(test="Immune state x (DMG vs pHGG-WT) interaction, LRT 2 df",p=pint)]).to_csv(ROOT+"out/T_within_subtype_interaction.tsv",sep="\t",index=False)
print(E.round(3).to_string()); print("interaction",round(pint,3)); print(pd.DataFrame(SUM).round(3).to_string())
