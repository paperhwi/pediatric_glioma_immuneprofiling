"""Supplementary figures (n = 343):
FigureS_state_structure: A seed partitions; B k = 2 vs k = 3 nesting; C agreement of re-clustering within molecular groups,
after removing group / library means, and within RNA library; D reference-based immune fraction by state; E MFP class by state;
F cross-method total-immune percentile by state.
FigureS_immune_features: A immune features and programmes with q < 0.05 (OS HR per SD, base and + immune axis); B published
programmes; C EXPLORATORY tumour-cell programme differences between molecular groups (deconvolution-based)."""
import sys, os, json; sys.path.insert(0,"/home/claude/q343/src"); from common import *
from s08_fig7 import plt, lab, NullLocator, fmtp
import anndata as ad
from matplotlib.gridspec import GridSpec
a=load()
fig=plt.figure(figsize=(8.2,10.6)); gs=GridSpec(3,2,hspace=0.7,wspace=0.45,left=0.1,right=0.97,top=0.96,bottom=0.06)
P=pd.read_csv(ROOT+"out/T_seed_partitions_n343.tsv",sep="\t"); ax=fig.add_subplot(gs[0,0]); bt=np.zeros(len(P))
for e in [LO,MID,HI]: ax.bar(P.seed.astype(str),P[f"n_{SHORT[e]}"],bottom=bt,color=COL[e],label=SHORT[e],width=0.85); bt+=P[f"n_{SHORT[e]}"].values
ax.set_xticks(range(0,25,4)); ax.set_xticklabels(P.seed.iloc[::4]); ax.set_xlabel("Consensus-clustering seed"); ax.set_ylabel("Tumours"); ax.legend(frameon=False,fontsize=6,ncol=3,loc="upper center",bbox_to_anchor=(.5,1.17),title="Immune state",title_fontsize=6)
st=pd.read_csv(ROOT+"out/T_seed_stability.tsv",sep="\t").iloc[0]; ax.set_title(f"Seed partitions (modal {st.modal_partition}; median ARI {st.median_pairwise_ARI:.2f})",fontsize=7,pad=30); lab(ax,"A")
ax=fig.add_subplot(gs[0,1]); N=pd.read_csv(ROOT+"out/r53/R53_nesting.tsv",sep="\t",index_col=0).loc[ORD]; im=ax.imshow(N.values,cmap="Blues")
for i in range(N.shape[0]):
    for j in range(N.shape[1]): ax.text(j,i,N.values[i,j],ha="center",va="center",fontsize=8,color="white" if N.values[i,j]>80 else "black")
ax.set_yticks(range(3)); ax.set_yticklabels(ORD,fontsize=6.5); ax.set_xticks(range(N.shape[1])); ax.set_xticklabels(["k = 2 cluster 1","k = 2 cluster 2"],fontsize=6.5)
J=json.load(open(ROOT+"out/r53/R53_summary.json")); ax.set_title(f"k = 3 states nested in k = 2 clusters\nLRT k = 3 vs k = 2: P = {J['LRT']['k3 vs k2']['p']:.2f}",fontsize=7); lab(ax,"B",-0.45)
ax=fig.add_subplot(gs[1,0]); R=json.load(open(ROOT+"out/r54/R54_summary.json")); L=json.load(open(ROOT+"out/r6lib/r6/R6_library_summary.json"))
items=[("Within DMG (n = 175)",R["DMG"]),("Within pHGG-WT (n = 127)",R["pHGG-WT"]),("Molecular-group means removed",R["group_residualised_global"]),("Library means removed",L["library_residualised_global"]),
       ("Within stranded (n = 126)",L["stranded"]),("Within poly-A stranded (n = 116)",L["poly-A stranded"]),("Within exome capture (n = 73)",L["exome_capture"])]
y=np.arange(len(items))[::-1]; ax.barh(y,[100*v["agreement_k3"] for _,v in items],color=["#1D4ED8"]*3+["#9CA3AF"]*4)
for yi,(_,v) in zip(y,items): ax.text(100*v["agreement_k3"]+1,yi,f"{100*v['agreement_k3']:.0f}% (ARI {v['ARI_k3']:.2f})",va="center",fontsize=6)
ax.set_yticks(y); ax.set_yticklabels([k for k,_ in items],fontsize=6.3); ax.set_xlim(0,130); ax.set_xticks([0,25,50,75,100]); ax.set_xlabel("Agreement with the final states (%)"); ax.set_title("Re-clustering within subsets (k = 3)",fontsize=7); lab(ax,"C",-0.62)
ax=fig.add_subplot(gs[1,1]); T=ad.read_h5ad("/home/claude/w/pg/Revision_round1_Q1-Q5/Round2_R5_2Nov/R5_tumours_n344.h5ad")[a.index]
imm=pd.Series(np.asarray(T.X)[:,2:].sum(1),index=a.index)*100
bp=ax.boxplot([imm[a.state==e] for e in ECO],patch_artist=True,showfliers=False,widths=0.6)
for b_,e in zip(bp["boxes"],ECO): b_.set(facecolor=COL[e],alpha=0.8)
D=pd.read_csv(ROOT+"out/T_deconvolution_by_state.tsv",sep="\t"); r=D[(D.method.str.startswith("Insta"))&(D.component=="Immune (total)")].iloc[0]
ax.set_xticks([1,2,3]); ax.set_xticklabels([SHORT[e] for e in ECO]); ax.set_ylabel("Immune cells, % (InstaPrism, ScPCA reference)"); ax.set_title(f"Reference-based deconvolution\nJT P = {fmtp(r.JT_p)}; ρ with immune axis = {r.spearman_immune_axis:.2f}",fontsize=7); lab(ax,"D")
ax=fig.add_subplot(gs[2,0]); M=pd.read_csv(ROOT+"out/T_Bagaev_MFP_by_state.tsv",sep="\t",index_col=0)[ORD]; Mp=M.div(M.sum(0),axis=1)*100; bt=np.zeros(3); cm={"D":"#9CA3AF","F":"#F59E0B","IE":"#1D4ED8","IE/F":"#7C3AED"}
for c_ in Mp.index: ax.bar([SHORT[e] for e in ORD],Mp.loc[c_],bottom=bt,color=cm.get(c_,"#ccc"),label=c_); bt+=Mp.loc[c_].values
ax.set_ylabel("Tumours (%)"); ax.legend(frameon=False,fontsize=6,title="MFP class",title_fontsize=6,bbox_to_anchor=(1.0,1),loc="upper left"); ax.set_title("Bagaev MFP classes by immune state",fontsize=7); lab(ax,"E")
ax=fig.add_subplot(gs[2,1]); O=pd.read_csv(ROOT+"out/r52/R52_state_ordering_percentile.tsv",sep="\t"); O=O[O.axis.str.startswith("Total immune")]
x=np.arange(len(O)); w=0.26
for i,(e,c) in enumerate([(LO,"median_pct_low"),(MID,"median_pct_intermediate"),(HI,"median_pct_high")]): ax.bar(x+(i-1)*w,O[c]*100,w,color=COL[e],label=SHORT[e])
ax.set_xticks(x); ax.set_xticklabels(O.method.str.replace(" (LM22)","",regex=False),rotation=30,ha="right",fontsize=6.3); ax.set_ylabel("Median within-method percentile"); ax.set_ylim(0,115); ax.legend(frameon=False,fontsize=6,ncol=3,loc="upper center",title="Immune state",title_fontsize=6)
ax.set_title("Total immune estimate, five methods\n(all JT P < 1e-15)",fontsize=7); lab(ax,"F")
os.makedirs(ROOT+"fig/FigureS_state_structure",exist_ok=True)
for ext in ("svg","pdf","png"): fig.savefig(ROOT+f"fig/FigureS_state_structure/FigureS_state_structure.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight",facecolor="white")
# ---------------- immune features
Q=pd.read_csv(ROOT+"out/q1/Q1_cox_markers_and_immune_features.tsv",sep="\t"); Q=Q[(Q.dataset=="primary")&(Q.family!="Molecular marker")]
keep=Q.groupby("feature").q_BH_within_family.min(); keep=keep[keep<0.05].index
W=Q[Q.feature.isin(keep)].pivot_table(index="feature",columns="model",values=["HR","lo","hi","q_BH_within_family"])
W=W.sort_values(("HR","base"))
fig=plt.figure(figsize=(8.2,9.6)); gs=GridSpec(2,2,height_ratios=[1.25,1],hspace=0.45,wspace=1.1,left=0.30,right=0.97,top=0.96,bottom=0.07)
ax=fig.add_subplot(gs[0,:]); y=np.arange(len(W))[::-1]
for off,m,c in [(0.18,"base","#1D4ED8"),(-0.18,"+ immune axis","#F59E0B")]:
    ax.errorbar(W[("HR",m)],y+off,xerr=[W[("HR",m)]-W[("lo",m)],W[("hi",m)]-W[("HR",m)]],fmt="o",color=c,ms=3.5,lw=1,label={"base":"molecular group + age + sex","+ immune axis":"+ continuous immune axis"}[m])
    for yi,(h,qq) in zip(y,zip(W[("HR",m)],W[("q_BH_within_family",m)])):
        if qq<0.05: ax.text(1.75,yi+off,f"q = {qq:.3f}",fontsize=5.5,va="center",color=c)
ax.axvline(1,ls="--",lw=0.6,color="k"); ax.set_xscale("log"); ax.set_xlim(0.5,2.2); ax.set_xticks([0.5,0.75,1,1.5]); ax.set_xticklabels(["0.5","0.75","1","1.5"]); ax.xaxis.set_minor_locator(NullLocator())
ax.set_yticks(y); ax.set_yticklabels([f.replace("ssGSEA_","").replace("LM22_","LM22 ") for f in W.index],fontsize=6.3); ax.set_xlabel("OS hazard ratio per SD (95% CI)"); ax.legend(frameon=False,fontsize=6,loc="upper left",bbox_to_anchor=(0.0,-0.09),ncol=2)
ax.set_title("Immune features and curated programmes with BH q < 0.05 in either model (n = 245, 206 deaths)",fontsize=7.2,loc="left"); lab(ax,"A",-0.5)
ax=fig.add_subplot(gs[1,0]); C=pd.read_csv(ROOT+"out/T_published_programmes_cox.tsv",sep="\t"); C=C[~C.feature.str.startswith("gene:")]
fe=C.feature.unique(); y=np.arange(len(fe))[::-1]
for off,m,c in [(0.15,"base","#1D4ED8"),(-0.15,"+ immune axis","#F59E0B")]:
    d=C[C.model==m].set_index("feature").loc[fe]; ax.errorbar(d.HR_per_SD,y+off,xerr=[d.HR_per_SD-d.lo,d.hi-d.HR_per_SD],fmt="o",color=c,ms=3.5,lw=1)
ax.axvline(1,ls="--",lw=0.6,color="k"); ax.set_yticks(y); ax.set_yticklabels(["T cell-inflamed GEP (Ayers)","Cytolytic activity (Rooney)","MHC-I APS (Wang)","ICD, 34 genes (Garg)","ICD DAMP/ER-stress subset"][:len(fe)],fontsize=6.3); ax.set_xlabel("HR per SD"); ax.set_title("Published programmes",fontsize=7.2); lab(ax,"B",-1.2)
ax=fig.add_subplot(gs[1,1]); E=pd.read_csv(ROOT+"out/T_EXPLORATORY_tumour_cell_programmes_by_group.tsv",sep="\t"); pr=E.programme.unique(); y=np.arange(len(pr))[::-1]
gc={"DMG_K27":"#7C3AED","DHG_G34":"#059669","IHG":"#D97706"}
for k,(g,c) in enumerate(gc.items()):
    d=E[E.group==g].set_index("programme").loc[pr]; ax.errorbar(d.effect,y+(k-1)*0.22,xerr=[d.effect-d.lo,d.hi-d.effect],fmt="o",color=c,ms=3,lw=0.9,label=g.split("_")[0]+" vs pHGG-WT")
ax.axvline(0,ls="--",lw=0.6,color="k"); ax.set_yticks(y); ax.set_yticklabels(pr,fontsize=6); ax.set_xlabel("Difference (SD), adjusted for state, age, sex"); ax.legend(frameon=False,fontsize=5.5,loc="upper left",bbox_to_anchor=(0,-0.16),ncol=3)
ax.set_title("EXPLORATORY: deconvolution-estimated\ntumour-cell programmes by molecular group",fontsize=7); lab(ax,"C",-1.2)
os.makedirs(ROOT+"fig/FigureS_immune_features",exist_ok=True)
for ext in ("svg","pdf","png"): fig.savefig(ROOT+f"fig/FigureS_immune_features/FigureS_immune_features.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print("ok")
