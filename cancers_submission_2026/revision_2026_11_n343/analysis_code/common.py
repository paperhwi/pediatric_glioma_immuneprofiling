"""Shared definitions for the n = 343 revision (Rev343).
Single exclusion rule: every specimen that is not a pediatric-type diffuse HGG (5 non-HGG + PT_8V1YD45K, methylation
class PXA) is removed before any analysis, including clustering. Immune states are re-derived on n = 343 with the locked
procedure (modal label over consensus-clustering seeds 42-66)."""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter, KaplanMeierFitter
from statsmodels.stats.multitest import multipletests
HI,MID,LO="Immune-high","Immune-intermediate","Immune-low"
ECO=[LO,MID,HI]                      # reference first
ORD=[HI,MID,LO]                      # plotting order
SHORT={HI:"High",MID:"Intermediate",LO:"Low"}   # figure legends only
COL={HI:"#3B82F6",MID:"#EF4444",LO:"#9CA3AF"}
GRP=["pHGG_WT","DMG_K27","DHG_G34","IHG"]
GRPNAME={"pHGG_WT":"pHGG, H3/IDH-wildtype","DMG_K27":"DMG, H3 K27-altered","DHG_G34":"DHG, H3 G34-mutant","IHG":"Infant-type hemispheric glioma"}
ROOT="/home/claude/q343/"
# Four OS records are re-timed from the diagnosis of the sequenced HGG. The corrected times are derived from OpenPedCan v15
# histologies by the single rule in s02_os_audit.py; individual survival times are not redistributed here, so the values are read
# from a local, non-public file (participant -> OS days), which s02_os_audit.py regenerates and checks.
import os as _os
OS_FIX=(pd.read_csv(_os.path.join(_os.path.dirname(__file__),"restricted_os_fix.tsv"),sep="\t",index_col=0).iloc[:,0].astype(int).to_dict()
        if _os.path.exists(_os.path.join(_os.path.dirname(__file__),"restricted_os_fix.tsv")) else {})
def design(s,extra=None,eco=True,grp=True,age=True,sex=True,grp_levels=None):
    d=pd.DataFrame({"time":s.time.values,"event":s.event.values},index=s.index)
    if eco:
        d[HI]=(s.state==HI).astype(float); d[MID]=(s.state==MID).astype(float)
    if grp:
        for g in (grp_levels or ["DMG_K27","DHG_G34","IHG"]): d[g]=(s.cohort_group==g).astype(float)
    if age: d["age_years"]=s.age_first_dx_years
    if sex: d["male"]=(s.reported_gender=="Male").astype(float)
    if extra is not None:
        for c in extra.columns: d[c]=extra.loc[s.index,c].values
    return d.dropna()
def cox(d,penalizer=0.0,**kw): return CoxPHFitter(penalizer=penalizer).fit(d,"time","event",**kw)
def safe_cox(d):
    keep=[c for c in d.columns if c in ("time","event","entry") or d[c].nunique()>1]; d=d[keep]
    try: return cox(d)
    except Exception: return cox(d,penalizer=0.01)
def lp(cf,X): return cf.predict_partial_hazard(X[list(cf.params_.index)]).values.ravel()
def bh(p):
    p=np.asarray(p,float); q=np.full_like(p,np.nan); m=~np.isnan(p)
    if m.any(): q[m]=multipletests(p[m],method="fdr_bh")[1]
    return q
def load(os_version="corrected",labels="n343"):
    """os_version: 'corrected' (primary) or 'released' (OpenPedCan v15 as distributed).
    labels: 'n343' (primary) or 'n344_carried' (n = 344 labels carried over, sensitivity)."""
    if os_version=="corrected" and len(OS_FIX)!=4:
        raise RuntimeError("Corrected OS requires the four locally audited corrections. Run s02_os_audit.py first, then start a fresh process; individual durations are not bundled publicly.")
    a=pd.read_pickle(ROOT+"ann343.pkl").copy()
    if labels!="n343": a["state"]=a["state_"+labels]
    a["OS_days_released"]=a.OS_days
    if os_version=="corrected":
        for p,t in OS_FIX.items(): a.loc[a.Kids_First_Participant_ID==p,"OS_days"]=float(t)
    a["os_corrected"]=a.Kids_First_Participant_ID.isin(OS_FIX)&(os_version=="corrected")
    a["age_first_dx_days"]=np.where(a.OS_days.notna()&a.age_last_update_days.notna(),a.age_last_update_days-a.OS_days,a.age_at_diagnosis_days)
    a["age_first_dx_years"]=a.age_first_dx_days/365.25
    a["sampling_delay_days"]=a.age_at_event_days-a["age_first_dx_days"]
    a["time"]=a.OS_days; a["event"]=(a.OS_status=="DECEASED").astype(float)
    a["os_eval"]=a.OS_days.notna()&a.OS_status.notna()&(a.OS_days>0)
    a["autopsy"]=(a.tumor_descriptor=="Deceased")
    return a
def state_ct(s): return s.state.value_counts().reindex(ORD).fillna(0).astype(int)
