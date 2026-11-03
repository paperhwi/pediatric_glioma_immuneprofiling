"""Supplementary analyses recomputed from per-tumour outputs that do not depend on state labels:
(a) reference-based deconvolution fractions (InstaPrism balanced markers; MuSiC) by immune state;
(b) pan-cancer frameworks (Thorsson immune subtype; Bagaev MFP) vs immune state;
(c) EXPLORATORY: deconvolution-estimated tumour-cell programme scores by molecular group, adjusted for immune state, age and sex."""
import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import *
import anndata as ad, statsmodels.formula.api as smf
from scipy import stats
a=load(); ids=a.index; pc1=pd.read_csv(ROOT+"out/immune_axis_n343.tsv",sep="\t",index_col=0).immune_axis
T=ad.read_h5ad("/home/claude/w/pg/Revision_round1_Q1-Q5/Round2_R5_2Nov/R5_tumours_n344.h5ad")[ids]
def jt(groups):
    U=0; ns=[len(x) for x in groups]
    for i in range(3):
        for j in range(i+1,3):
            A_,B_=np.asarray(groups[i]),np.asarray(groups[j]); U+=(A_[:,None]<B_[None,:]).sum()+0.5*(A_[:,None]==B_[None,:]).sum()
    N=sum(ns); mu=(N**2-sum(n**2 for n in ns))/4; var=(N**2*(2*N+3)-sum(n**2*(2*n+3) for n in ns))/72; z=(U-mu)/np.sqrt(var)
    return z,2*stats.norm.sf(abs(z))
def cliff(x,y): x,y=np.asarray(x),np.asarray(y); return ((x[:,None]>y[None,:]).mean()-(x[:,None]<y[None,:]).mean())
rows=[]
for meth,M in [("InstaPrism (balanced markers)",pd.DataFrame(np.asarray(T.X),index=ids,columns=T.var_names)),("MuSiC",pd.DataFrame(T.obsm["MuSiC"],index=ids,columns=T.uns["obsm_columns"]["MuSiC"]))]:
    M["Immune (total)"]=M[["Microglia","Macrophage/monocyte","T cell","NK cell","B/plasma"]].sum(1); M["Myeloid"]=M[["Microglia","Macrophage/monocyte"]].sum(1)
    M["MDM share of myeloid"]=M["Macrophage/monocyte"]/M["Myeloid"].replace(0,np.nan)
    for c in M.columns:
        g=[M.loc[a.state==e,c].dropna() for e in ECO]; z,p=jt(g)
        d=pd.DataFrame({"y":M[c],"s":a.state,"g":a.cohort_group}).dropna()
        rows.append(dict(method=meth,component=c,**{f"median_{SHORT[e]}":M.loc[a.state==e,c].median() for e in ECO},JT_z=z,JT_p=p,
            cliff_high_vs_low=cliff(g[2],g[0]),cliff_high_vs_intermediate=cliff(g[2],g[1]),cliff_intermediate_vs_low=cliff(g[1],g[0]),
            spearman_immune_axis=stats.spearmanr(M[c],pc1,nan_policy="omit").statistic,R2_state=smf.ols("y~C(s)",d).fit().rsquared,R2_group=smf.ols("y~C(g)",d).fit().rsquared))
D=pd.DataFrame(rows); D["JT_q_BH"]=bh(D.JT_p); D.to_csv(ROOT+"out/T_deconvolution_by_state.tsv",sep="\t",index=False)
pd.set_option("display.width",250); print(D[D.component.isin(["Immune (total)","Myeloid","T cell"])].round(4).to_string())
# within-group rho with axis for total immune
M=pd.DataFrame(np.asarray(T.X),index=ids,columns=T.var_names); imm=M[["Microglia","Macrophage/monocyte","T cell","NK cell","B/plasma"]].sum(1)
W=pd.DataFrame([dict(group=g,n=int((a.cohort_group==g).sum()),rho_immune_fraction_vs_axis=stats.spearmanr(imm[a.cohort_group==g],pc1[a.cohort_group==g]).statistic) for g in GRP]); W.to_csv(ROOT+"out/T_deconvolution_axis_within_group.tsv",sep="\t",index=False); print(W)
# (b) frameworks
F=[]
for col in ["Thorsson_subtype","Bagaev_MFP"]:
    ct=pd.crosstab(T.obs[col].astype(str),a.state).reindex(columns=ORD).fillna(0).astype(int); ct.to_csv(ROOT+f"out/T_{col}_by_state.tsv",sep="\t")
    chi=stats.chi2_contingency(ct); F.append(dict(framework=col,chi2=chi[0],df=chi[2],p=chi[1])); print(ct)
pd.DataFrame(F).to_csv(ROOT+"out/T_frameworks_tests.tsv",sep="\t",index=False); print(pd.DataFrame(F))
# (c) exploratory evasion programmes
E=ad.read_h5ad("/home/claude/w/pg/Revision_round1_Q1-Q5/Review_gaps_2Nov/R6_tumour_intrinsic_evasion_n344.h5ad")[ids]
X=pd.DataFrame(np.asarray(E.X),index=ids,columns=E.var_names); R=[]
for p_ in X.columns:
    d=pd.DataFrame({"y":X[p_],"g":pd.Categorical(a.cohort_group,["pHGG_WT","DMG_K27","DHG_G34","IHG"]),"s":a.state,"age":a.age_first_dx_years,"sex":a.reported_gender}).dropna()
    f=smf.ols("y~C(g)+C(s)+age+C(sex)",d).fit(); f0=smf.ols("y~C(s)+age+C(sex)",d).fit(); om=f.compare_f_test(f0)[1]
    for g in ["DMG_K27","DHG_G34","IHG"]:
        k=f"C(g)[T.{g}]"; R.append(dict(programme=p_,group=g,vs="pHGG_WT",n=len(d),effect=f.params[k],lo=f.conf_int().loc[k,0],hi=f.conf_int().loc[k,1],p=f.pvalues[k],group_omnibus_p=om))
R=pd.DataFrame(R); R["q_BH"]=bh(R.p); R.to_csv(ROOT+"out/T_EXPLORATORY_tumour_cell_programmes_by_group.tsv",sep="\t",index=False); print(R.round(3).to_string())
