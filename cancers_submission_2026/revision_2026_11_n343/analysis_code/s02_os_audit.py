"""OS-origin audit applied with one rule to every tumour in the n = 343 cohort."""
import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
H=pd.read_csv("/home/claude/opc/histologies.tsv",sep="\t",low_memory=False)
a=load(os_version="released"); TOL=30
P=H[H.Kids_First_Participant_ID.isin(a.Kids_First_Participant_ID)&(H.sample_type=="Tumor")].copy()
P["age_at_diagnosis_days"]=pd.to_numeric(P.age_at_diagnosis_days,errors="coerce")
rows=[]
for bs,r in a.iterrows():
    has=pd.notna(r.OS_days) and pd.notna(r.OS_status)
    O=r.age_last_update_days-r.OS_days if has else np.nan; D=r.age_at_diagnosis_days
    pp=P[P.Kids_First_Participant_ID==r.Kids_First_Participant_ID]
    other=pp[(pp.age_at_diagnosis_days<D-TOL)]
    other_dx="; ".join(sorted(set(f"{x.pathology_diagnosis} at {x.age_at_diagnosis_days/365.25:.1f} y" for _,x in other.iterrows())))
    sm=(r.tumor_descriptor=="Second Malignancy")
    gap=D-O if has else np.nan
    efs_match=has and pd.notna(r.EFS_days) and abs(r.EFS_days-gap)<=TOL
    if not has: cat,act="No OS record","not in survival analyses"
    elif r.OS_days<=0: cat,act="OS <= 0 d","not in survival analyses"
    elif pd.isna(r.age_last_update_days): cat,act="R1b origin not checkable (age at last follow-up missing); specimen is the initial tumour and OS >= EFS","accept"
    elif abs(gap)<=TOL: cat,act="R1 OS origin = diagnosis of the sequenced HGG","accept"
    elif gap>TOL and len(other)>0: cat,act="R3a OS origin = a different, earlier tumour documented in the record","re-time from HGG diagnosis"
    elif gap>TOL and (sm or r.EFS_event_type=="Second Malignancy" and r.tumor_descriptor!="Initial CNS Tumor"): cat,act="R3b sequenced HGG recorded as a second malignancy; OS origin = earlier malignancy","re-time from HGG diagnosis"
    elif gap>TOL and r.tumor_descriptor=="Initial CNS Tumor": cat,act="R3c sequenced specimen is the initial tumour but OS origin precedes it (internally inconsistent)","re-time from HGG diagnosis"
    elif gap>TOL and efs_match: cat,act="R2 OS origin = earlier diagnosis of the same HGG (later specimen; interval equals recorded event-free interval)","accept"
    elif gap<-TOL: cat,act="OS origin after diagnosis","review"
    else: cat,act="unclassified","review"
    corrected=r.age_last_update_days-D if act.startswith("re-time") else (r.OS_days if has else np.nan)
    rows.append(dict(Kids_First_Biospecimen_ID=bs,Kids_First_Participant_ID=r.Kids_First_Participant_ID,state=r.state,cohort_group=r.cohort_group,
        tumor_descriptor=r.tumor_descriptor,OS_status=r.OS_status,OS_days_released=r.OS_days,age_dx_sequenced_HGG_days=D,age_at_specimen_days=r.age_at_event_days,
        age_last_update_days=r.age_last_update_days,implied_OS_origin_days=O,origin_precedes_dx_by_days=gap,EFS_days=r.EFS_days,EFS_event_type=r.EFS_event_type,
        EFS_equals_interval=efs_match,earlier_other_diagnoses=other_dx,category=cat,action=act,OS_days_used=corrected,
        specimen_after_last_update=(r.age_at_event_days>r.age_last_update_days) if pd.notna(r.age_last_update_days) else np.nan))
T=pd.DataFrame(rows); T.to_csv(ROOT+"out/T_OS_origin_audit_all_n343.tsv",sep="\t",index=False)
S=T.groupby(["category","action"]).size().rename("n").reset_index(); S.to_csv(ROOT+"out/T_OS_origin_audit_summary.tsv",sep="\t",index=False)
print(S.to_string()); print(T[T.action!="accept"][T.action!="not in survival analyses"][["Kids_First_Participant_ID","category","OS_days_released","OS_days_used"]].to_string())
fix=T[T.action.str.startswith("re-time")].set_index("Kids_First_Participant_ID").OS_days_used.astype(int).to_dict()
assert len(fix)==4, "The final n343 audit must identify four origin corrections."
import os
f=os.path.join(os.path.dirname(os.path.abspath(__file__)),"restricted_os_fix.tsv")
if not OS_FIX: pd.Series(fix,name="OS_days").to_csv(f,sep="\t"); print("wrote",f,"- rerun the pipeline")
else: assert fix==OS_FIX, "rule-derived corrections differ from those used"; print("Corrections derived by the rule == corrections used in all analyses")
