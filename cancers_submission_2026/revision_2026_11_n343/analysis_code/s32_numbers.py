"""One file with every number to be carried into the manuscript, figure legends and response (n = 343)."""
import sys, json, glob; sys.path.insert(0,"/home/claude/q343/src"); from common import *
O=ROOT+"out/"; R=[]
def add(sec,item,val,prev="",src="",note=""): R.append(dict(section=sec,item=item,value_n343=val,previous_n344=prev,source=src,note=note))
f=lambda x,d=2: f"{x:.{d}f}"
def pf(p):
    if p<0.001: m,e=f"{p:.1e}".split("e"); return f"{m} × 10^{int(e)}"
    return f"{p:.3f}" if p<0.1 else f"{p:.2f}"
a=load(); s=a[a.os_eval]
# 1 cohort
add("1 Cohort","Analysis cohort",f"{len(a)} pediatric-type diffuse HGG (one RNA specimen per participant)","344","out/immune_state_assignments_n343.tsv","single rule: all 6 non-diffuse-HGG specimens (5 non-HGG + PXA) removed before any analysis")
add("1 Cohort","Molecular groups","DMG 175 · DHG 31 · pHGG-WT 127 · IHG 10","same, plus 1 PXA","")
td=a.tumor_descriptor.replace({"Primary Tumor":"Initial","Initial CNS Tumor":"Initial","Recurrence":"Rec","Progressive":"Rec","Second Malignancy":"Rec","Deceased":"Autopsy"}).value_counts()
add("1 Cohort","Specimen timing",f"{td['Initial']} initial · {td['Rec']} recurrent/progressive · {td['Autopsy']} autopsy","245/41/58")
st=pd.read_csv(O+"T_seed_stability.tsv",sep="\t").iloc[0]
add("1 Cohort","Immune states (modal over seeds 42–66)","Immune-high 122 · Immune-intermediate 144 · Immune-low 77","122/145/77 (LI/MD/ID)","out/T_seed_stability.tsv")
add("1 Cohort","Labels changed vs n = 344","1 tumour Immune-intermediate → Immune-high (BS_P7W3AYKV, modal frequency 0.64); 342/343 unchanged","","out/T_label_changes_n344_to_n343.tsv","the procedure reproduces the n = 344 labels 344/344")
add("1 Cohort","Seed stability (25 seeds)",f"median pairwise ARI {st.median_pairwise_ARI:.2f} (range {st.min_ARI:.2f}–{st.max_ARI:.2f}); {int(st.tumours_same_label_all_25)}/343 same label in all 25 runs; seeds 0–24 and 0–99 give the same partition","0.87 (0.80–1.00); 313/344","out/T_seed_stability.tsv")
from sklearn.metrics import adjusted_rand_score as ari
add("1 Cohort","ARI final vs submitted locked labels (n = 349 labels, 343 tumours)",f(ari(a.state,a.state_locked349)))
# 2 Fig 2
cm=pd.read_csv(O+"T_fig2B_cross_method_n343.tsv",sep="\t"); js=json.load(open(ROOT+"fig/Figure2/Figure2_summary.json"))
add("2 Figure 2","LM22 deconvolution P < 0.05",f"{js['qc_p_below_005']}/343 ({js['qc_pct_below_005']:.1f}%)","88/344 (25.6%)")
add("2 Figure 2","Cross-method correlations BH q < 0.05",f"{int((cm.q_BH<0.05).sum())}/12; rho range {cm.spearman_rho.min():.2f} to {cm.spearman_rho.max():.2f}","4/12")
# 3 Fig 3
k=pd.read_csv(O+"fig3/ksweep_n343.tsv",sep="\t").set_index("k"); g=pd.read_csv(O+"fig3/gap_n343.tsv",sep="\t"); c3=json.load(open(O+"fig3/consensus_validation_summary.json"))
add("3 Figure 3","PAC k = 2 / k = 3",f"{k.PAC[2]:.3f} / {k.PAC[3]:.3f} (lowest at k = 2)"); add("3 Figure 3","Silhouette k = 2 / k = 3",f"{k.silhouette_consensus[2]:.3f} / {k.silhouette_consensus[3]:.3f} (highest at k = 2)")
add("3 Figure 3","Relative ΔAUC at k = 3 (largest remaining)",f"+{100*k.delta_AUC[3]:.1f}%","+29.8% (submitted)"); add("3 Figure 3","Gap statistic, one-SE rule",f"k = {c3['gap_1SE_selected_k']} selected; k = 1 not selected")
add("3 Figure 3","Single-seed 1,000-replicate k = 3 vs final labels",f"ARI {c3['current_vs_denovo_k3_ARI']:.2f}","0.80")
# 4 Fig 4
F4=pd.read_csv(O+"T_fig4_statistics.tsv",sep="\t")
for _,r in F4.iterrows():
    if r.panel in "ABC": add("4 Figure 4 / Tables 2–4",r.analysis,f"χ² {r.chi2:.2f}, df {int(r.df)}, P_MC {r.p_montecarlo:.4f}, Cramér's V {r.cramers_V:.2f}, n = {int(r.n)}",{"A":"χ² 17.31, P_MC 0.0082, n = 327","B":"χ² 30.70, P_MC 0.0001, n = 343","C":"n = 86, P = 0.41"}[r.panel])
    else: add("4 Figure 4 / Tables 2–4",r.analysis,f"partial R² {100*r.partial_R2:.1f}%, F {r.F:.2f}, P {r.p_perm:.4f}, n = {int(r.n)}")
# 5 Fig 5
D5=pd.read_csv(O+"fig5b/T_fig5B_designs.tsv",sep="\t"); r=D5[D5.gene_set.str.startswith("top 300")&D5.design.str.startswith("D4")].iloc[0]
add("5 Figure 5","5B top 300 Immune-high vs Immune-low genes highest in an immune class (10 ScPCA types)",f"{int(r.n_immune_top)}/300 ({r.pct_immune_top:.1f}%) vs matched random {r.random_mean_pct:.1f}% ({r.random_lo:.1f}–{r.random_hi:.1f}); P = 0.005","293/300 (97.7%) vs 70.2%")
for d in ["D1","D3"]:
    r=D5[D5.gene_set.str.startswith("top 300")&D5.design.str.startswith(d)].iloc[0]; add("5 Figure 5",f"5B design {r.design}",f"{r.pct_immune_top:.1f}% vs random {r.random_mean_pct:.1f}%")
r=D5[D5.gene_set.str.contains("shared")&D5.design.str.startswith("D4")].iloc[0]
add("5 Figure 5","404 shared genes (Fig 6D) highest in an immune class (10 types)",f"{r.pct_immune_top:.1f}% vs random {r.random_mean_pct:.1f}% (not above random)","417 genes: 61–74%, close to random","","the 404 genes are immune-associated plus a mesenchymal/ECM/vascular component")
V=pd.read_csv(O+"T_visium_by_section.tsv",sep="\t")
for _,r in V.iterrows(): add("5 Figure 5",f"Visium {r.section} (n = {r.n_spots} spots)",f"High/Intermediate/Low spots {r.n_High}/{r.n_Intermediate}/{r.n_Low}; median immune fraction {r.median_immune_Low:.3f}/{r.median_immune_Intermediate:.3f}/{r.median_immune_High:.3f} (Low/Int/High); JT P {pf(r.JT_p_immune_low_int_high)}; ρ(axis, immune) {r.rho_axis_vs_immune:.2f}; mixed spots {r.pct_spots_mixed:.0f}%; class agreement with n = 344 {100*r.class_agreement_with_n344:.0f}%")
# 6 Fig 6
S6=json.load(open(O+"T_fig6_summary.json"))
add("6 Figure 6","Genes tested; up/down",f"{S6['genes_tested']}; Immune-high vs Immune-low {S6['counts']['Immune-high vs Immune-low'][0]}/{S6['counts']['Immune-high vs Immune-low'][1]}; Immune-intermediate vs Immune-low {S6['counts']['Immune-intermediate vs Immune-low'][0]}/{S6['counts']['Immune-intermediate vs Immune-low'][1]}; Immune-high vs Immune-intermediate {S6['counts']['Immune-high vs Immune-intermediate'][0]}/{S6['counts']['Immune-high vs Immune-intermediate'][1]}","3681/441/758 up")
add("6 Figure 6","Shared up-regulated genes",f"{S6['shared']} ({S6['pct_intermediate_in_high']:.1f}% of the Immune-intermediate set)","417")
add("6 Figure 6","Gradient",f"orthogonal slope {S6['orthogonal_slope']:.2f}; Pearson r {S6['pearson_r']:.2f}; median ratio on shared genes {S6['median_ratio_shared']:.2f}","slope 1.76")
add("6 Figure 6","Paired Wilcoxon on shared genes (two-sided)",pf(S6['wilcoxon_two_sided_p']),"1.5 × 10^-70")
# 7 Fig 7
P0=pd.read_csv(O+"T_primary_cox.tsv",sep="\t",index_col=0); SE=pd.read_csv(O+"T_survival_sensitivity.tsv",sep="\t")
add("7 Figure 7 (main)","Survival set",f"n = {len(s)}, {int(s.event.sum())} deaths; log-rank P = {SE.logrank_p[0]:.3f}","245/206; 0.058")
for t in ["Immune-high vs Immune-low","Immune-intermediate vs Immune-low","DMG vs pHGG-WT","DHG vs pHGG-WT","IHG vs pHGG-WT","Age at diagnosis (per year)","Male vs female"]:
    r=P0.loc[t]; add("7 Figure 7 (main)",f"Primary model: {t}",f"HR {r['exp(coef)']:.2f} ({r['exp(coef) lower 95%']:.2f}–{r['exp(coef) upper 95%']:.2f}), P = {pf(r.p)}","LI 0.62 (0.42–0.91), P = 0.016" if t.startswith("Immune-high") else "")
add("7 Figure 7 (main)","Proportional hazards (Schoenfeld, rank time)","state terms P = 0.058 (high) and 0.079 (intermediate); no term P < 0.05 except DHG (0.039)","","out/T_primary_cox_PH_test.tsv")
# 8 axis
AX=pd.read_csv(O+"T_immune_axis_cox.tsv",sep="\t").iloc[0]; sc=pd.read_csv(O+"immune_axis_scree.tsv",sep="\t")
est=pd.read_csv(O+"ESTIMATE_n343.tsv",sep="\t",index_col=0).iloc[:,0]; pc1=pd.read_csv(O+"immune_axis_n343.tsv",sep="\t",index_col=0).immune_axis
from scipy.stats import spearmanr
add("8 Continuous immune axis (main)","Definition","PC1 of the 21 immune ssGSEA signatures, z-scored within n = 343",""); add("8 Continuous immune axis (main)","Variance explained",f"{100*sc.variance_explained[0]:.1f}%","65.8%")
add("8 Continuous immune axis (main)","Spearman ρ with ESTIMATE ImmuneScore",f"{spearmanr(pc1,est)[0]:.2f}","0.94")
add("8 Continuous immune axis (main)","HR per SD (molecular group, age, sex)",f"{AX.HR:.2f} ({AX.lo:.2f}–{AX.hi:.2f}), P = {AX.p:.4f}","0.83, P = 0.010")
add("8 Continuous immune axis (main)","Immune state beyond the axis (LRT 2 df)",f"P = {AX.LRT_state_beyond_axis_p:.2f}","0.65 (Wald)")
# 9 sensitivity
for _,r in SE.iterrows(): add("9 Survival sensitivity (supplement)",r.analysis,f"n = {r.n}, {r.deaths} deaths; Immune-high HR {r.high_HR:.2f} ({r.high_lo:.2f}–{r.high_hi:.2f}), P = {pf(r.high_p)}; Immune-intermediate {r.intermediate_HR:.2f}, P = {pf(r.intermediate_p)}")
AP=pd.read_csv(O+"T_autopsy_by_state.tsv",sep="\t")
for _,r in AP.iterrows(): add("10 Autopsy specimens",f"{r.set}: {r.state}",f"{int(r.autopsy)}/{int(r.n)} autopsy ({r['autopsy_%']:.1f}%)"+(f"; deaths {int(r.deaths)}; non-autopsy n = {int(r.n_non_autopsy)}, deaths {int(r.deaths_non_autopsy)}" if pd.notna(r.deaths) else "")+(f"; {r.test}" if isinstance(r.test,str) and r.test else ""))
FU=pd.read_csv(O+"T_followup.tsv",sep="\t",index_col=0)
for i,r in FU.iterrows(): add("11 Follow-up",i,f"reverse-KM median {r.reverse_KM_median_FU_y:.1f} y (95% CI {r.reverse_KM_FU_95CI}); survivors median {r.survivors_median_FU_y:.1f} y, IQR {r.survivors_FU_IQR_y}, range {r.survivors_FU_range_y} y; observation time range {r.observation_time_range_y} y; at risk 2/5/10 y {r.at_risk_2y}/{r.at_risk_5y}/{r.at_risk_10y}")
E=pd.read_csv(O+"T_within_subtype_state_cox.tsv",sep="\t"); ip=pd.read_csv(O+"T_within_subtype_interaction.tsv",sep="\t").p[0]; SU=pd.read_csv(O+"T_within_subtype_survival_summary.tsv",sep="\t").set_index("group")
for _,r in E.iterrows(): add("12 Within molecular subtype",f"{r.group}: {r.term}",f"HR {r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f}), P = {pf(r.p)}; n = {r.n}, {r.events} deaths")
add("12 Within molecular subtype","Log-rank within DMG / pHGG-WT",f"{SU.logrank_p['DMG_K27']:.3f} / {SU.logrank_p['pHGG_WT']:.3f}"); add("12 Within molecular subtype","State × subtype interaction (LRT 2 df)",f"P = {ip:.2f}")
RW=json.load(open(O+"r54/R54_summary.json"))
add("12 Within molecular subtype (main text)","Re-clustering within DMG (n = 175)",f"gap k = {RW['DMG']['gap_k']}; k = 3 agreement {100*RW['DMG']['agreement_k3']:.0f}%, ARI {RW['DMG']['ARI_k3']:.2f}","86%, ARI 0.64")
add("12 Within molecular subtype (main text)","Re-clustering within pHGG-WT (n = 127)",f"gap k = {RW['pHGG-WT']['gap_k']}; k = 3 agreement {100*RW['pHGG-WT']['agreement_k3']:.0f}%, ARI {RW['pHGG-WT']['ARI_k3']:.2f}","96%, ARI 0.88")
add("12 Within molecular subtype (main text)","Molecular-group means removed (n = 343)",f"agreement {100*RW['group_residualised_global']['agreement_k3']:.0f}%, ARI {RW['group_residualised_global']['ARI_k3']:.2f}","86%")
add("12 Within molecular subtype (main text)","Median R² of the 46 features",f"molecular group {RW['variance_partition_median']['R2_group']:.2f}; immune state {RW['variance_partition_median']['R2_state']:.2f}","0.03 / 0.11")
Q=pd.read_csv(O+"T_incremental_value.tsv",sep="\t"); n=Q[(Q.model=="M1 M0+immune state")].iloc[0]; b=Q[(Q.model=="M0 dx+age+sex")].iloc[0]
add("13 Incremental value","Immune state added to molecular group + age + sex",f"LRT P = {n.LRT_p:.3f}; ΔAIC {n.AIC-b.AIC:.1f}; optimism-corrected C {b.Harrell_C_optimism_corrected:.3f} → {n.Harrell_C_optimism_corrected:.3f} (ΔC +{n.delta_C_corrected:.3f}, 95% CI {n.delta_C_boot_lo:.3f}–{n.delta_C_boot_hi:.3f})","0.036; −2.6; 0.606 → 0.635")
K=json.load(open(O+"r53/R53_summary.json"))
add("14 k = 2 vs k = 3","k = 2 clusters",f"{K['k2_sizes']}; Immune-intermediate split 52/92","MD split 91/54"); add("14 k = 2 vs k = 3","LRT k = 3 vs k = 2",f"P = {K['LRT']['k3 vs k2']['p']:.2f}","0.13")
add("14 k = 2 vs k = 3","Immune-high vs Immune-intermediate",f"{K['LI_vs_MD_DEG']['up_in_high']} genes up / {K['LI_vs_MD_DEG']['up_in_intermediate']} down; {K['LI_vs_MD_features_q05']} of 46 features differ (q < 0.05)","857 genes; 27 features")
T1=pd.read_csv(O+"T_Table1_n343.tsv",sep="\t"); add("15 Table 1","DNA alterations by state",f"none differ (minimum BH q = {T1.q_BH_molecular.min():.2f}); TP53 altered 50–52% in every state","all q ≥ 0.54")
add("15 Table 1","Clinical variables with P < 0.05","molecular group, RNA library, extent of resection, specimen timing, anatomical region, age","same")
Q1=pd.read_csv(O+"q1/Q1_cox_markers_and_immune_features.tsv",sep="\t"); v=Q1[Q1.dataset=="primary"]
def nsig(fam,m): x=v[(v.family==fam)&(v.model==m)]; return f"{int((x.q_BH_within_family<0.05).sum())}/{len(x)}: "+", ".join(x.loc[x.q_BH_within_family<0.05,"feature"].str.replace("ssGSEA_","").str.replace("LM22_",""))
add("16 Immune features (supplement)","46 features, base / + immune axis",nsig("Deconvolution / brain signature (46)","base")+" | "+nsig("Deconvolution / brain signature (46)","+ immune axis"),"3 / 3")
add("16 Immune features (supplement)","17 curated programmes, base / + immune axis",nsig("Curated immune program (17)","base")+" | "+nsig("Curated immune program (17)","+ immune axis"),"10/17 / Treg programme")
mk=v[(v.family=="Molecular marker")]; tp=mk[mk.feature=="TP53_alt"].set_index("model"); bc=mk[mk.feature=="BCOR_mut"].set_index("model")
add("16 Immune features (supplement)","DNA markers",f"none at q < 0.05; TP53 HR {tp.HR['base']:.2f} (P = {tp.p['base']:.3f}), {tp.HR['+ immune state']:.2f} with state; BCOR {bc.HR['base']:.2f} (P = {bc.p['base']:.3f})","TP53 1.40 / 1.54; BCOR 2.45")
GS=pd.read_csv(O+"T_gene_screen_shared_genes.tsv",sep="\t"); add("16 Immune features (supplement)","Gene screen of shared genes",f"{int((GS.q_BH<0.05).sum())}/{len(GS)} at q < 0.05, all protective","83/417")
PC=pd.read_csv(O+"T_published_programmes_cox.tsv",sep="\t"); x=PC[PC.feature.str.startswith("Immunogenic")].set_index("model")
add("16 Immune features (supplement)","ICD (Garg 2016) HR per SD",f"{x.HR_per_SD['base']:.2f} (q = {x.q_BH['base']:.3f}); + axis {x.HR_per_SD['+ immune axis']:.2f} (q = {x.q_BH['+ immune axis']:.2f})","0.79 (q 0.008); 0.82 (q 0.25)")
add("16 Immune features (supplement)","Checkpoint genes by group given axis and state (q < 0.05)","DHG: PDCD1LG2, TGFB2, TIGIT up; DMG: VTCN1 up, PDCD1LG2 and TIGIT down; IHG: IDO1 up (others q ≥ 0.08)","same","out/T_published_programmes_by_state_group.tsv")
DC=pd.read_csv(O+"T_deconvolution_by_state.tsv",sep="\t"); r=DC[(DC.method.str.startswith("Insta"))&(DC.component=="Immune (total)")].iloc[0]
add("17 Deconvolution (supplement)","InstaPrism immune fraction by state",f"Low {100*r.median_Low:.1f}% · Intermediate {100*r.median_Intermediate:.1f}% · High {100*r.median_High:.1f}%; JT P = {pf(r.JT_p)}; ρ with axis {r.spearman_immune_axis:.2f}; R² state {r.R2_state:.2f}, group {r.R2_group:.2f}","2.8/8.7/21.4%; ρ 0.78; R² 0.55/0.06")
CMJ=json.load(open(O+"r52/R52_summary.json")); add("17 Deconvolution (supplement)","Cross-method concordance (matched types)",f"{CMJ['concordance_q05']}/{CMJ['concordance_n']} at q < 0.05; total immune Low < Intermediate < High in all five methods (JT P 6 × 10^-17 to 2 × 10^-41)","6/15")
FT=pd.read_csv(O+"T_frameworks_tests.tsv",sep="\t").set_index("framework")
add("18 Pan-cancer frameworks (supplement)","Thorsson",f"300/343 C4; states not separated (χ² P = {FT.p['Thorsson_subtype']:.2f})","301/344; P = 0.08")
add("18 Pan-cancer frameworks (supplement)","Bagaev MFP","Immune-low → D 76/77; Immune-high → IE or IE/F 85/122 (70%); Immune-intermediate D 91/144 (63%)","76/77; 69%; 63%")
LB=json.load(open(O+"r6lib/r6/R6_library_summary.json"))
add("19 RNA library (supplement)","Library by state","exome capture: 52/73 Immune-high; stranded: 56/126 Immune-low (Table 1 P = 0.0001)","51/73; 56/126")
add("19 RNA library (supplement)","Within-library re-clustering ARI",f"stranded {LB['stranded']['ARI_k3']:.2f}, poly-A stranded {LB['poly-A stranded']['ARI_k3']:.2f}, exome capture {LB['exome_capture']['ARI_k3']:.2f}; library means removed {LB['library_residualised_global']['ARI_k3']:.2f}","0.24–0.59")
add("19 RNA library (supplement)","Survival with library","adjusted HR 0.64 (P = 0.042); stratified 0.62 (P = 0.039)","0.64 / 0.62")
CH=pd.read_csv(O+"T_scpca_harmonisation_metrics.tsv",sep="\t",index_col=0).value; LT=pd.read_csv(O+"T_scpca_axis_linkage_donor_tests.tsv",sep="\t")
r=LT[LT.measure=="Immune cells, % of all cells"].iloc[0]
add("20 Single cell (supplement)","ScPCA harmonised classifier",f"CV accuracy {CH['CV_accuracy']:.3f} ({int(CH['CV_correct'])}/343); PC1 {100*CH['PC1_variance']:.1f}%; r with TPM axis {CH['r_with_TPM_axis']:.2f}; projected Low/Intermediate/High 4/7/6 (identical to n = 344)","0.878 (302/344); 4/7/6")
add("20 Single cell (supplement)","Bulk axis vs single-cell immune %",f"ρ = {r.spearman_rho:.2f}, P(perm) = {r.P_perm_rho:.4f}; trend P = {r.JT_trend_P_perm:.3f}","0.78")
for m,c in [("% PRF1+","CD8"),("Immunosuppressive TAM","Myeloid"),("MDM identity","Myeloid"),("% CTLA4+","CD4")]:
    r=LT[(LT.measure==m)&(LT.compartment==c)].iloc[0]; add("20 Single cell (supplement)",f"{c} {m} vs axis",f"ρ = {r.spearman_rho:.2f}, P(perm) = {r.P_perm_rho:.3f}, q = {r.q_P_perm_rho:.2f}; trend P = {r.JT_trend_P_perm:.3f}")
AU=pd.read_csv(O+"T_OS_origin_audit_summary.tsv",sep="\t")
for _,r in AU.iterrows(): add("21 OS-origin audit (all 343)",r.category,f"{r.n} ({r.action})")
EV=pd.read_csv(O+"T_EXPLORATORY_tumour_cell_programmes_by_group.tsv",sep="\t")
add("22 EXPLORATORY tumour-cell programmes (supplement only)","Programme-level group differences given state, age, sex",f"no group effect at q < 0.05 (minimum q = {EV.q_BH.min():.2f})","gene-level CD47/CD200 etc. reported previously; not to be claimed in the main text")
T=pd.DataFrame(R); T.to_csv(ROOT+"deliv/Manuscript_numbers_n343.tsv",sep="\t",index=False)
with open(ROOT+"deliv/Manuscript_numbers_n343.md","w") as fh:
    fh.write("# Numbers for the manuscript, legends and response — final cohort n = 343 (3 Nov 2026)\n\nEvery value is read from the result files in `out/`. `previous_n344` is the value in the 2 Nov package, for tracking changes.\n\n")
    for sec,g in T.groupby("section",sort=False):
        fh.write(f"## {sec}\n\n| Item | Value (n = 343) | Previous (n = 344) |\n|---|---|---|\n")
        for _,r in g.iterrows(): fh.write(f"| {r['item']} | {r.value_n343} | {r.previous_n344} |\n")
        fh.write("\n")
print(len(T)); print(T.head(60).to_string()[:6000])
