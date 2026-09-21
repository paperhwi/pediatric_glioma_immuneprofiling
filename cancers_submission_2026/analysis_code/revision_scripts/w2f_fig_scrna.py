import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
mpl.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":.8})
BLUE,RED,GREY,GREEN,ORANGE="#2563EB","#B91C1C","#7C8798","#0F766E","#F97316"
A=pd.read_csv("W2_A_signature_attribution.tsv",sep="\t")
U=pd.read_csv("W2_A_auroc.tsv",sep="\t").sort_values("AUROC_all_cells")
B=pd.read_csv("W2_B_gene_source.tsv",sep="\t"); Bs=pd.read_csv("W2_B_summary.tsv",sep="\t")
S=pd.read_pickle("w2_sigscores_10x.pkl"); lab=S.pop("cell_type"); S.pop("sampleid")

fig=plt.figure(figsize=(13.0,7.4),dpi=300)
gs=fig.add_gridspec(2,3,width_ratios=[1.25,1.0,1.0],height_ratios=[1,.85],hspace=.52,wspace=.52)

# A — mean z of each signature per cell type
ax=fig.add_subplot(gs[:,0])
Z=(S-S.mean())/S.std()
M=Z.groupby(lab).mean().T[["Myeloid","CD4","CD8"]]
M=M.loc[M.max(axis=1).sort_values().index]
im=ax.imshow(M.values,cmap="RdBu_r",vmin=-.8,vmax=.8,aspect="auto")
ax.set_yticks(range(len(M))); ax.set_yticklabels([s.replace("_"," ") for s in M.index],fontsize=6.6)
ax.set_xticks(range(3)); ax.set_xticklabels(["Myeloid\n13,920","CD4 T\n2,385","CD8 T\n2,314"],fontsize=7.2)
for i,s in enumerate(M.index):
    if s=="Glioma_Inflammatory_Wang2017":
        ax.add_patch(plt.Rectangle((-.5,i-.5),3,1,fill=False,ec=RED,lw=1.6))
        ax.text(2.62,i,"grouped with\nmyeloid theme",fontsize=6,color=RED,va="center")
ax.set_title("A  Signature score by annotated cell type\n18,619 immune cells, 4 pHGG samples (GSE227983)",
             fontsize=8,loc="left")
cax=ax.inset_axes([0.0,-0.085,1.0,0.022])
cb=fig.colorbar(im,cax=cax,orientation="horizontal")
cb.set_label("mean z of control-matched score",fontsize=6.3); cb.ax.tick_params(labelsize=5.8)

# B — AUROC
ax=fig.add_subplot(gs[0,1])
col=[RED if v<0.5 else (ORANGE if v<0.7 else GREEN) for v in U.AUROC_all_cells]
ax.barh(range(len(U)),U.AUROC_all_cells,color=col,height=.72)
ax.errorbar(U.AUROC_all_cells,range(len(U)),
            xerr=[np.clip(U.AUROC_all_cells-U.AUROC_min_sample,0,None),
                  np.clip(U.AUROC_max_sample-U.AUROC_all_cells,0,None)],
            fmt="none",ecolor="#334155",elinewidth=.7,capsize=1.6)
ax.axvline(.5,color="k",ls="--",lw=.7)
ax.set_yticks(range(len(U))); ax.set_yticklabels([s.replace("_"," ") for s in U.signature],fontsize=6.2)
ax.set_xlabel("AUROC for the named cell type",fontsize=7.2); ax.set_xlim(0,1)
ax.set_title("B  Does each signature identify the cell type it is named for?\nbars, all cells; whiskers, per-sample range",fontsize=8,loc="left")
ax.tick_params(labelsize=6.6); ax.spines[["top","right"]].set_visible(False)
ax.annotate("anti-correlated with\nmyeloid identity",xy=(0.075,0),xytext=(0.34,3.4),fontsize=6,color=RED,
            arrowprops=dict(arrowstyle="->",color=RED,lw=.8))

# C — the discordant signature
ax=fig.add_subplot(gs[1,1])
g="Glioma_Inflammatory_Wang2017"
data=[S.loc[lab==k,g].values for k in ["Myeloid","CD4","CD8"]]
bp=ax.boxplot(data,widths=.6,patch_artist=True,showfliers=False,medianprops=dict(color="k",lw=1))
for p,c in zip(bp["boxes"],[GREY,BLUE,BLUE]): p.set_facecolor(c); p.set_alpha(.5); p.set_edgecolor(c)
ax.set_xticks([1,2,3]); ax.set_xticklabels(["Myeloid","CD4 T","CD8 T"],fontsize=7)
ax.set_ylabel("signature score",fontsize=7)
ax.set_title("C  Glioma_Inflammatory_Wang2017\ngenes: GZMB, PRF1, IFNG, GZMA, NKG7, CD8A, CD8B, CCL5",fontsize=7.6,loc="left",color=RED)
ax.tick_params(labelsize=6.6); ax.spines[["top","right"]].set_visible(False)

# D — source of ecotype-defining genes
ax=fig.add_subplot(gs[0,2])
c=B.top_cell_type.value_counts().reindex(["Myeloid","Tcell","Malignant"]).fillna(0)
bars=ax.bar(range(3),c.values,color=[GREY,BLUE,RED],width=.62)
for i,v in enumerate(c.values): ax.text(i,v+4,int(v),ha="center",fontsize=7.6,weight="bold")
ax.set_xticks(range(3)); ax.set_xticklabels(["Myeloid","T cell","Malignant"],fontsize=7)
ax.set_ylabel("ecotype-defining genes",fontsize=7); ax.set_ylim(0,225)
ax.set_title("D  Which cell type expresses the genes that\ndefine the ecotype axis?  (n = 246, Smart-seq2)",fontsize=8,loc="left")
ax.tick_params(labelsize=6.6); ax.spines[["top","right"]].set_visible(False)

# E — observed vs expression-matched background
ax=fig.add_subplot(gs[1,2])
obs=float(Bs.observed_immune_fraction.iloc[0])*100
lo,hi,mu=float(Bs.background_lo.iloc[0])*100,float(Bs.background_hi.iloc[0])*100,float(Bs.background_mean.iloc[0])*100
ax.barh([0],[mu],xerr=[[mu-lo],[hi-mu]],color=GREY,height=.4,error_kw=dict(elinewidth=1,capsize=3))
ax.barh([1],[obs],color=GREEN,height=.4)
ax.set_yticks([0,1]); ax.set_yticklabels(["expression-matched\nrandom genes","ecotype-defining\ngenes"],fontsize=7)
ax.set_xlim(0,105); ax.set_xlabel("% attributed to immune cells",fontsize=7.2)
ax.text(obs-2,1,f"{obs:.1f}%",ha="right",va="center",fontsize=7.4,color="white",weight="bold")
ax.text(mu+ (hi-mu) +2,0,f"{mu:.1f}%",va="center",fontsize=7)
ax.set_title(f"E  Permutation P = {float(Bs.perm_P.iloc[0]):.3f}",fontsize=8,loc="left")
ax.tick_params(labelsize=6.6); ax.spines[["top","right"]].set_visible(False)

fig.suptitle("Single-cell support for the cellular interpretation of the immune signatures and ecotypes",
             fontsize=9.5,y=.985)
for e in ("png","pdf"): fig.savefig(f"FigureS20_scRNA_support.{e}",dpi=300,bbox_inches="tight",facecolor="white")
print("saved FigureS20")
