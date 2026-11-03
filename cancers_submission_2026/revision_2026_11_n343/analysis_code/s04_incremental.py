import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
from lifelines.utils import concordance_index
from sksurv.metrics import concordance_index_ipcw, cumulative_dynamic_auc
from sksurv.util import Surv
from scipy.stats import chi2
a=load(); M=pd.read_pickle(ROOT+"markers343.pkl"); pc1=pd.read_csv(ROOT+"out/immune_axis_n343.tsv",sep="\t",index_col=0).immune_axis
def prep(s,tp53=False):
    ex=pd.DataFrame({"immune_PC1":pc1})
    if tp53: ex["TP53"]=M.TP53_alt
    return design(s,extra=ex,eco=True)
SPEC={"M0 dx+age+sex":["DMG_K27","DHG_G34","IHG","age_years","male"],
      "M1 M0+immune state":["DMG_K27","DHG_G34","IHG","age_years","male",HI,MID],
      "M2 M0+immune axis":["DMG_K27","DHG_G34","IHG","age_years","male","immune_PC1"],
      "M3 M0+TP53":["DMG_K27","DHG_G34","IHG","age_years","male","TP53"],
      "M4 M3+immune state":["DMG_K27","DHG_G34","IHG","age_years","male","TP53",HI,MID]}
TAU=5*365.25; TIMES=np.array([1,2,3])*365.25
def metrics(dtr,dte,cols):
    c=cox(dtr[["time","event"]+cols]); lp_te=c.predict_partial_hazard(dte[cols]).values.ravel()
    ytr=Surv.from_arrays(dtr.event.astype(bool),dtr.time); yte=Surv.from_arrays(dte.event.astype(bool),dte.time)
    h=concordance_index(dte.time,-lp_te,dte.event)
    u=concordance_index_ipcw(ytr,yte,lp_te,tau=TAU)[0]
    t=TIMES[(TIMES<dte.time.max())&(TIMES>dte.time.min())]
    auc,_=cumulative_dynamic_auc(ytr,yte,lp_te,t)
    return c,h,u,dict(zip([f"tdAUC_{int(x/365.25)}y" for x in t],auc))
def evaluate(d,label,pairs,B=500,seed=20261101):
    rows=[]; fits={}
    for k in pairs:
        cols=SPEC[k]; c,h,u,auc=metrics(d,d,cols); fits[k]=c
        rows.append(dict(analysis=label,model=k,n=len(d),events=int(d.event.sum()),df=len(cols),logLik=c.log_likelihood_,AIC=c.AIC_partial_,Harrell_C=h,Uno_C_tau5y=u,**auc))
    R=pd.DataFrame(rows).set_index("model")
    # LRT nested
    base,new=pairs
    lr=2*(fits[new].log_likelihood_-fits[base].log_likelihood_); ddf=len(SPEC[new])-len(SPEC[base])
    R.loc[new,"LRT_vs_base_chi2"]=lr; R.loc[new,"LRT_df"]=ddf; R.loc[new,"LRT_p"]=chi2.sf(lr,ddf)
    # bootstrap optimism (Harrell) for each model and for delta
    rng=np.random.default_rng(seed); opt={k:[] for k in pairs}; dlt=[]
    for b in range(B):
        i=rng.integers(0,len(d),len(d)); db=d.iloc[i].reset_index(drop=True)
        cb={}
        for k in pairs:
            cols=SPEC[k]
            cf=safe_cox(db[["time","event"]+cols])
            app=concordance_index(db.time,-lp(cf,db),db.event)
            orig=concordance_index(d.time,-lp(cf,d),d.event)
            opt[k].append(app-orig); cb[k]=(app,orig)
        dlt.append((cb[new][1]-cb[base][1]))
    for k in pairs:
        R.loc[k,"optimism"]=np.mean(opt[k]); R.loc[k,"Harrell_C_optimism_corrected"]=R.loc[k,"Harrell_C"]-np.mean(opt[k])
    dC=R.loc[new,"Harrell_C_optimism_corrected"]-R.loc[base,"Harrell_C_optimism_corrected"]
    R.loc[new,"delta_C_corrected"]=dC
    # bootstrap CI for apparent delta C
    rng=np.random.default_rng(seed+1); app_d=[]
    for b in range(B):
        i=rng.integers(0,len(d),len(d)); db=d.iloc[i].reset_index(drop=True)
        hs=[]
        for k in pairs:
            cols=SPEC[k]; cf=safe_cox(db[["time","event"]+cols]); hs.append(concordance_index(db.time,-lp(cf,db),db.event))
        app_d.append(hs[1]-hs[0])
    R.loc[new,"delta_C_apparent"]=R.loc[new,"Harrell_C"]-R.loc[base,"Harrell_C"]
    R.loc[new,"delta_C_boot_lo"]=np.percentile(app_d,2.5); R.loc[new,"delta_C_boot_hi"]=np.percentile(app_d,97.5)
    R.loc[new,"delta_C_boot_P_le0"]=np.mean(np.array(app_d)<=0)
    return R.reset_index()
out=[]
s=a[a.os_eval]; d=prep(s).reset_index(drop=True)
out.append(evaluate(d,"Primary (n=245)",("M0 dx+age+sex","M1 M0+immune state")))
out.append(evaluate(d,"Primary (n=245)",("M0 dx+age+sex","M2 M0+immune axis")))
dt=prep(s,tp53=True).reset_index(drop=True)
out.append(evaluate(dt,"TP53-extended base (DNA subset)",("M3 M0+TP53","M4 M3+immune state")))
T=pd.concat(out); T.to_csv(ROOT+"out/T_incremental_value.tsv",sep="\t",index=False)
pd.set_option("display.width",260); print(T.round(4).to_string())
