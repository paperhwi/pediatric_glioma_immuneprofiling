"""AnnData objects for checking: tumours (n = 343) and Visium spots (n = 1,650)."""
import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
import anndata as ad
R="/home/claude/repo/cancers_submission_2026/reproducibility_data/"
a=load(); F=pd.read_csv(R+"feature_matrices/main_cohort_LM22_ssGSEA_z.tsv",sep="\t",index_col=0).loc[a.index]
loc=pd.read_csv(ROOT+"out/T_fig4_per_tumour_location.tsv",sep="\t",index_col=0).location_class
obs=pd.DataFrame({"immune_state":pd.Categorical(a.state,ORD),"immune_state_n344_carried":a.state_n344_carried,"label_submitted_n349":a.state_locked349,"modal_frequency_seeds42_66":a.modal_frequency_seeds42_66,
    "molecular_group":a.cohort_group,"location_class":loc,"specimen_timing":a.tumor_descriptor,"autopsy_specimen":a.autopsy,"RNA_library":a.RNA_library,
    "age_at_diagnosis_years":a.age_first_dx_years,"sex":a.reported_gender,"OS_days_corrected":a.time,"OS_event":a.event,"OS_corrected_record":a.os_corrected,"in_survival_analysis":a.os_eval,
    "immune_axis":pd.read_csv(ROOT+"out/immune_axis_n343.tsv",sep="\t",index_col=0).immune_axis,"ESTIMATE_ImmuneScore_z":pd.read_csv(ROOT+"out/ESTIMATE_n343.tsv",sep="\t",index_col=0).iloc[:,0]},index=a.index)
for c in ["immune_state_n344_carried","label_submitted_n349","molecular_group","location_class","specimen_timing","RNA_library","sex"]: obs[c]=obs[c].astype(str)
T=ad.AnnData(X=F.values.astype(np.float32),obs=obs,var=pd.DataFrame(index=F.columns))
P=pd.read_pickle(ROOT+"programs343.pkl").loc[a.index]; T.obsm["curated_programmes"]=P.values; T.uns["curated_programmes"]=list(P.columns)
old=ad.read_h5ad("/home/claude/w/pg/Revision_round1_Q1-Q5/Round2_R5_2Nov/R5_tumours_n344.h5ad")[a.index]
T.obsm["InstaPrism_fractions"]=np.asarray(old.X); T.uns["InstaPrism_fractions"]=list(old.var_names)
T.uns["description"]="n = 343 pediatric-type diffuse HGG (OpenPedCan v15). X = 46-feature immune matrix (LM22 + 24 ssGSEA, z). immune_state = modal label over consensus seeds 42-66 (Immune-high / Immune-intermediate / Immune-low). OS corrected for 4 records by one rule. Not for public release (contains survival times)."
T.write_h5ad(ROOT+"deliv/Rev343_tumours_n343.h5ad",compression="gzip")
V=ad.read_h5ad("/mnt/user-data/uploads/Pediatric GBM/Revision_round1_Q1-Q5/Round2_R5_2Nov/R5_visium_spots.h5ad")
d=pd.read_csv(ROOT+"out/T_visium_spots_n1650.tsv",sep="\t",index_col=0).loc[V.obs_names]
V.obs["spot_class"]=pd.Categorical(d.spot_class,ORD); V.obs["axis_score"]=d.axis_score.values; V.obs["spot_class_n344"]=d.spot_class_n344.astype(str).values
V.uns["description"]="GSE268577, 3 H3.3-K27M pDMG Visium sections. X = InstaPrism fractions (unchanged). spot_class = nearest n = 343 immune-state centroid on the shared Figure 6D genes (z within section)."
V.write_h5ad(ROOT+"deliv/Rev343_visium_spots_n1650.h5ad",compression="gzip"); print(T); print(V)
