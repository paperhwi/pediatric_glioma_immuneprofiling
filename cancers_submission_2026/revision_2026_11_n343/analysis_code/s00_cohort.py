"""Cohort (n = 343) and immune-state labels. Cohort rule identical to Rev344 cohort.py, plus PT_8V1YD45K removed."""
import pandas as pd, numpy as np, re, sys
sys.path.insert(0,"/home/claude/q343/src"); from common import *
from sklearn.metrics import adjusted_rand_score as ari
R="/home/claude/repo/cancers_submission_2026/reproducibility_data/"
H=pd.read_csv("/home/claude/opc/histologies.tsv",sep="\t",low_memory=False)
F=pd.read_csv("/home/claude/opc/fusion-putative-oncogenic.tsv",sep="\t",low_memory=False,usecols=["Sample","FusionName"])
eco=pd.read_csv(R+"ecotype_assignments.tsv",sep="\t").rename(columns={"ecotype":"ecotype_locked349","cluster":"cluster_locked349"})
REMOVED={"BS_0XEBS3PG":"non-HGG (round 1)","BS_F6V1Y4QS":"non-HGG (round 1)","BS_HJV0J603":"non-HGG (round 1)","BS_S4V9E1XJ":"non-HGG (round 1)","BS_ZZWMD6FA":"non-HGG (round 1)",
         "BS_2DMBCHNF":"PT_8V1YD45K: methylation class pleomorphic xanthoastrocytoma"}
RELAB=["BS_DS438CXE","BS_X23724NJ"]
RTK={"NTRK1","NTRK2","NTRK3","ALK","ROS1","MET"}
F["genes"]=F.FusionName.fillna("").map(lambda s:set(re.split(r"--|/",s)))
rtk_samples=set(F.loc[F.genes.map(lambda g: len(g&RTK)>0),"Sample"])
a=eco.merge(H,on="Kids_First_Biospecimen_ID",how="left")
excl=a[a.Kids_First_Biospecimen_ID.isin(REMOVED)][["Kids_First_Biospecimen_ID","Kids_First_Participant_ID","harmonized_diagnosis","molecular_subtype","ecotype_locked349"]].copy()
excl["reason"]=excl.Kids_First_Biospecimen_ID.map(REMOVED); excl["submitted_label_n349"]=excl.pop("ecotype_locked349").map({"Lymphocyte-inflamed":"Immune-high","Myeloid-dominant":"Immune-intermediate","Immune-desert":"Immune-low"}); excl.to_csv(ROOT+"out/T_excluded_specimens.tsv",sep="\t",index=False)
a=a[~a.Kids_First_Biospecimen_ID.isin(REMOVED)].copy()
hd=a.harmonized_diagnosis.fillna(""); ms=a.molecular_subtype.fillna("")
age_m=a.age_at_diagnosis_days/30.4375
hemi=a.CNS_region.fillna("").eq("Hemispheric")|a.primary_site.fillna("").str.contains("Frontal|Parietal|Temporal|Occipital|Hemispheric",case=False)
rtk=a.Kids_First_Biospecimen_ID.isin(rtk_samples)|a.sample_id.isin(rtk_samples)
ihg=(age_m<18)&hemi&(rtk|ms.str.contains("RTK|IHG",case=False))
g=pd.Series("Other",index=a.index)
g[hd.str.contains("IDH-wildtype and H3-wildtype",case=False)|ms.str.startswith("HGG, H3 wildtype")|ms.str.startswith("HGG, To be")]="pHGG_WT"
g[hd.str.contains("H3 G35",case=False)|ms.str.startswith("DHG, H3 G35")]="DHG_G34"
g[hd.str.contains("H3 K28",case=False)|ms.str.startswith("DMG, H3 K28")]="DMG_K27"
g[ihg]="IHG"
g[a.Kids_First_Biospecimen_ID.isin(RELAB)]="pHGG_WT"
a["cohort_group"]=g.values
assert len(a)==343 and (a.cohort_group=="Other").sum()==0, a.cohort_group.value_counts()
for c in ["OS_days","age_at_diagnosis_days","age_last_update_days","age_at_event_days","EFS_days"]: a[c]=pd.to_numeric(a[c],errors="coerce")
# ---- labels: modal partition over the locked seeds 42-66 (reproduces the n = 344 primary labels 344/344 on n = 344)
S343=pd.concat([pd.read_csv(ROOT+f"out/s343_{x}.tsv",sep="\t",index_col=0) for x in "ab"],axis=1); S343.columns=S343.columns.astype(int)
S344=pd.concat([pd.read_csv(ROOT+f"out/s344_{x}.tsv",sep="\t",index_col=0) for x in "ab"],axis=1); S344.columns=S344.columns.astype(int)
LOCK=list(range(42,67))
lab343=S343[LOCK].mode(axis=1)[0]; lab344=S344[LOCK].mode(axis=1)[0]
old=pd.read_excel("/home/claude/w/pg/Revision_round1_Q1-Q5/Rev344_05_tables_Q1-Q5.xlsx",sheet_name="ecotype_344").set_index("Kids_First_Biospecimen_ID").ecotype
MAP={"Lymphocyte-inflamed":HI,"Myeloid-dominant":MID,"Immune-desert":LO}
assert (lab344==old.map(MAP).loc[lab344.index]).all(), "n=344 labels not reproduced"
a=a.set_index("Kids_First_Biospecimen_ID",drop=False)
a["state"]=lab343.loc[a.index].values
a["state_n344_carried"]=old.map(MAP).loc[a.index].values
a["state_locked349"]=a.ecotype_locked349.map(MAP)
a["modal_frequency_seeds42_66"]=S343[LOCK].eq(lab343,axis=0).mean(1).loc[a.index].values
a.to_pickle(ROOT+"ann343.pkl")
# ---- seed stability table (Table S26 equivalent)
rows=[]
for nm,S,lab in [("n = 343 (primary)",S343,lab343),("n = 344 (previous)",S344,lab344)]:
    A=[ari(S[i],S[j]) for i in LOCK for j in LOCK if i<j]; fr=S[LOCK].eq(lab,axis=0).mean(1)
    rows.append({"cohort":nm,"seeds":"42-66 (25)","modal_partition":"/".join(str(int((lab==e).sum())) for e in ORD),
        "median_pairwise_ARI":np.median(A),"min_ARI":min(A),"max_ARI":max(A),"tumours_same_label_all_25":int((fr==1).sum()),
        "min_modal_frequency":fr.min(),"modal_partition_seeds_0_24":"/".join(str(int((S[list(range(25))].mode(axis=1)[0]==e).sum())) for e in ORD),
        "modal_partition_seeds_0_99":"/".join(str(int((S[list(range(100))].mode(axis=1)[0]==e).sum())) for e in ORD)})
ST=pd.DataFrame(rows); ST.to_csv(ROOT+"out/T_seed_stability.tsv",sep="\t",index=False)
per=pd.DataFrame({"seed":LOCK,**{f"n_{SHORT[e]}":[int((S343[s]==e).sum()) for s in LOCK] for e in ORD},
                  "ARI_vs_modal":[ari(S343[s],lab343) for s in LOCK]}); per.to_csv(ROOT+"out/T_seed_partitions_n343.tsv",sep="\t",index=False)
ct=pd.crosstab(a.state_n344_carried.rename("n = 344 label"),a.state.rename("n = 343 label")).reindex(index=ORD,columns=ORD)
ct.to_csv(ROOT+"out/T_labels_n344_vs_n343.tsv",sep="\t")
mv=a[a.state!=a.state_n344_carried][["Kids_First_Participant_ID","cohort_group","state_n344_carried","state","modal_frequency_seeds42_66"]]
mv.to_csv(ROOT+"out/T_label_changes_n344_to_n343.tsv",sep="\t")
a[["Kids_First_Biospecimen_ID","Kids_First_Participant_ID","cohort_group","state","state_n344_carried","state_locked349","modal_frequency_seeds42_66"]].to_csv(ROOT+"out/immune_state_assignments_n343.tsv",sep="\t",index=False)
print(len(a),a.cohort_group.value_counts().to_dict()); print(state_ct(a).to_dict()); print(ST.T); print(ct); print(mv)
print("ARI n343 vs locked349", ari(a.state,a.state_locked349), " vs n344:",ari(a.state,a.state_n344_carried))
