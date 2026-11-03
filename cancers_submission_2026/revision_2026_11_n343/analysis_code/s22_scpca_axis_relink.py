"""ScPCA axis linkage (C3) re-projected with the n = 343 reference: harmonised log2-CPM ssGSEA scores of the 21 immune
signatures (unchanged per-sample inputs), PCA and logistic state classifier refitted on n = 343 with the final labels;
donor-level tests recomputed from the stored donor summaries (Spearman/OLS slope with 10,000 permutations; KW; JT trend)."""
import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from scipy import stats
import statsmodels.formula.api as smf
A="/mnt/user-data/uploads/Pediatric GBM/Revision_round1_Q1-Q5/Item_C_singlecell/ScPCA_validation/Axis_linkage/"
a=load(); ESo=pd.read_csv(A+"openpedcan_bulk_harmonised.tsv",sep="\t",index_col=0).loc[a.index]; ESq=pd.read_csv(A+"scpca_bulk_harmonised.tsv",sep="\t",index_col=0)
imm=[c for c in ESo.columns if c not in ("ESTIMATE_immune","immune_axis","ecotype")]; assert len(imm)==21
mu,sd=ESo[imm].mean(),ESo[imm].std(); Zo=(ESo[imm]-mu)/sd; Zq=(ESq[imm]-mu)/sd
pca=PCA(1).fit(Zo.values); sg=np.sign(np.corrcoef(pca.transform(Zo.values)[:,0],ESo["T_Cell_Cytotoxicity"])[0,1])
ao=sg*pca.transform(Zo.values)[:,0]; aq=sg*pca.transform(Zq.values)[:,0]; m0,s0=ao.mean(),ao.std()
axo=pd.Series((ao-m0)/s0,index=ESo.index); axq=pd.Series((aq-m0)/s0,index=ESq.index)
tpm_axis=pd.read_csv(ROOT+"out/immune_axis_n343.tsv",sep="\t",index_col=0).immune_axis
lab=a.state; clf=LogisticRegression(max_iter=5000)
cvp=cross_val_predict(clf,Zo.values,lab,cv=StratifiedKFold(5,shuffle=True,random_state=0)); acc=(cvp==lab.values).mean()
clf.fit(Zo.values,lab); pr=pd.DataFrame(clf.predict_proba(Zq.values),index=ESq.index,columns=clf.classes_)
H=pd.DataFrame({"immune_axis":axq,"projected_state":pr.idxmax(1),"projected_state_prob":pr.max(1)})
old=pd.read_csv(A+"link_donors.tsv",sep="\t",index_col=0)
MAP={"Lymphocyte-inflamed":HI,"Myeloid-dominant":MID,"Immune-desert":LO}
D=old[["location_group","total_cells","Myeloid","CD8","CD4","T_unassigned"]].join(H); D["projected_state_n344"]=old.projected_state.map(MAP)
met=dict(CV_accuracy=acc,CV_correct=int((cvp==lab.values).sum()),PC1_variance=pca.explained_variance_ratio_[0],r_with_TPM_axis=float(np.corrcoef(axo,tpm_axis.loc[axo.index])[0,1]),
         scpca_axis_min=axq.loc[D.index].min(),scpca_axis_max=axq.loc[D.index].max(),**{f"OpenPedCan_median_axis_{SHORT[e]}":axo[lab==e].median() for e in ORD},
         **{f"projected_{SHORT[e]}":int((D.projected_state==e).sum()) for e in ECO},projected_state_agreement_with_n344=float((D.projected_state==D.projected_state_n344).mean()),
         r_scpca_axis_with_n344=float(np.corrcoef(D.immune_axis,old.immune_axis)[0,1]))
print(met)
L=pd.read_csv(A+"link_donor_values_long.tsv",sep="\t"); L["axis"]=L.donor.map(D.immune_axis); L["state"]=L.donor.map(D.projected_state)
rng=np.random.default_rng(42); NPERM=10000
def perm_p(x,y,f):
    obs=f(x,y); null=np.array([f(rng.permutation(x),y) for _ in range(NPERM)]); return obs,(1+np.sum(np.abs(null)>=abs(obs)))/(NPERM+1)
slope=lambda x,y: np.polyfit((x-x.mean())/x.std(),y,1)[0]; rho=lambda x,y: stats.spearmanr(x,y)[0]
def jt(g,y):
    s=0
    for i in range(3):
        for j in range(i+1,3):
            A_=y[g==ECO[i]]; B_=y[g==ECO[j]]; s+=sum((bb>aa)+0.5*(bb==aa) for aa in A_ for bb in B_)
    return s
R=[]
for (comp,meas),d in L.groupby(["compartment","measure"],sort=False):
    x=d.axis.values; y=d.value.values; g=d.state.values
    r,pr_=perm_p(x,y,rho); b,pb=perm_p(x,y,slope); ols=smf.ols("value ~ axis + midline",d).fit()
    kw=stats.kruskal(*[y[g==s] for s in ECO if (g==s).sum()>0]).pvalue if len(set(g))>1 else np.nan
    e=jt(g,y); n=np.array([(g==s).sum() for s in ECO]); mu_=(len(g)**2-np.sum(n**2))/4
    null=np.array([jt(rng.permutation(g),y) for _ in range(2000)]); pjt=(1+np.sum(np.abs(null-mu_)>=abs(e-mu_)))/2001
    R.append(dict(compartment=comp,measure=meas,kind=d.kind.iloc[0],n_donors=len(d),donors_low_intermediate_high="/".join(map(str,n)),spearman_rho=r,P_perm_rho=pr_,
        slope_per_axis_SD=b,P_perm_slope=pb,slope_adj_midline=ols.params["axis"],P_adj_midline=ols.pvalues["axis"],
        **{f"median_{SHORT[s]}":(np.median(y[g==s]) if (g==s).any() else np.nan) for s in ECO},P_KW=kw,JT_trend_P_perm=pjt))
R=pd.DataFrame(R); fam=R.compartment.ne("Composition")
for col in ["P_perm_rho","JT_trend_P_perm"]:
    R.loc[fam,"q_"+col]=bh(R.loc[fam,col]); R.loc[~fam,"q_"+col]=bh(R.loc[~fam,col])
R.to_csv(ROOT+"out/T_scpca_axis_linkage_donor_tests.tsv",sep="\t",index=False); D.to_csv(ROOT+"out/T_scpca_donors_projection.tsv",sep="\t")
pd.DataFrame([met]).T.to_csv(ROOT+"out/T_scpca_harmonisation_metrics.tsv",sep="\t",header=["value"])
pd.set_option("display.width",250); print(R[["compartment","measure","n_donors","donors_low_intermediate_high","spearman_rho","P_perm_rho","q_P_perm_rho","median_Low","median_Intermediate","median_High","JT_trend_P_perm"]].round(3).to_string())
