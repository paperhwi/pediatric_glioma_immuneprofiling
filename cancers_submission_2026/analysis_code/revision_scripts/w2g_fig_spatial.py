import numpy as np, pandas as pd, pickle, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
mpl.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":.8})
COL={"Lymphocyte-inflamed":"#3B82F6","Myeloid-dominant":"#EF4444","Immune-desert":"#9CA3AF"}
ORDER=list(COL)
spot=pickle.load(open("w2_spot_scores.pkl","rb"))
R=pd.read_csv("W2_C_morans_I.tsv",sep="\t"); E=pd.read_csv("W2_C_spot_ecotype.tsv",sep="\t")
Cc=pd.read_csv("W2_C_colocalisation.tsv",sep="\t")

fig=plt.figure(figsize=(13.0,7.0),dpi=300)
gs=fig.add_gridspec(2,4,height_ratios=[.92,1.0],hspace=.30,wspace=.46)

for i,s in enumerate(["1","2","3"]):
    S,pos=spot[s]; a=np.load(f"w2_spot_assign_{s}.npy",allow_pickle=True)
    ax=fig.add_subplot(gs[0,i])
    ax.scatter(pos.array_col,-pos.array_row,c=[COL[x] for x in a],s=7,lw=0)
    ax.set_aspect("equal"); ax.axis("off")
    row=E[E["sample"].astype(str)==s].iloc[0]
    ax.set_title(f"{'ABC'[i]}  pDMG Sample-{s}   {int(row.n_spots)} spots\n"
                 f"same-neighbour {row.same_neighbour_fraction:.3f} vs null {row.null_mean:.3f}, P = {row.perm_p:.3f}",
                 fontsize=7.6,loc="left")
h=[plt.Line2D([],[],marker="o",ls="",color=COL[k],label=k,ms=5) for k in ORDER]
fig.legend(handles=h,loc="upper right",bbox_to_anchor=(.995,.93),fontsize=7,frameon=False,
           title="nearest bulk ecotype centroid",title_fontsize=7.5)

# D — Moran's I
ax=fig.add_subplot(gs[0,3]); ax.axis("off")
piv=R.pivot(index="signature",columns="sample",values="morans_I")
piv["mean"]=piv.mean(axis=1); piv=piv.sort_values("mean",ascending=True).tail(12)
ax2=fig.add_subplot(gs[1,0:2])
y=np.arange(len(piv))
for j,(c,m) in enumerate(zip(["Sample-1","Sample-2","Sample-3"],["o","s","^"])):
    ax2.scatter(piv[c],y,s=17,marker=m,label=c,alpha=.85,
                color=["#1E3A5F","#F97316","#0F766E"][j],lw=0)
ax2.axvline(0,color="k",lw=.7)
ax2.set_yticks(y); ax2.set_yticklabels([s.replace("_"," ") for s in piv.index],fontsize=6.6)
ax2.set_xlabel("Moran's I (spatial autocorrelation)",fontsize=7.4)
ax2.set_title("D  Immune programmes are spatially structured, not random\ntop 12 of 24 signatures by mean Moran's I; 999-permutation test",
              fontsize=8,loc="left")
ax2.legend(fontsize=6.5,frameon=False,loc="lower right"); ax2.tick_params(labelsize=6.6)
ax2.spines[["top","right"]].set_visible(False)

# E — myeloid vs lymphoid, bulk vs spot
ax=fig.add_subplot(gs[1,2])
ANN="/mnt/user-data/uploads/Open PBTA/Revision/FINAL MANUSCRIPT 260722 - 수정본/7. Reproducibility data/Sample annotation/sample_master_annotation.tsv"
d=pd.read_csv(ANN,sep="\t")
sig=[c for c in d.columns if c.startswith("ssGSEA_")]
Zb=(d[sig]-d[sig].mean())/d[sig].std(); Zb.columns=[c.replace("ssGSEA_","") for c in Zb.columns]
MY=["Microglia_Klemm2020","MDM_Klemm2020","MgTAM_Antunes2021","MoTAM_Antunes2021","DAM_KerenShaul2017","M2_Macrophage"]
LY=["T_Cell_Cytotoxicity","NK_Cell_Activity","Chemokine_T_Cell_Recruitment","T_Cell_Exhaustion","Tregs_Friebel2020"]
rb=np.corrcoef(Zb[MY].mean(axis=1),Zb[LY].mean(axis=1))[0,1]
vals=[rb]+list(Cc.myeloid_lymphoid_r)
labs=["bulk\nn = 349"]+[f"spot\nSample-{s}" for s in Cc["sample"]]
cols=["#1E3A5F","#0F766E","#0F766E","#0F766E"]
ax.bar(range(4),vals,color=cols,width=.62)
for i,v in enumerate(vals): ax.text(i,v+.02,f"{v:.2f}",ha="center",fontsize=7.2)
ax.set_xticks(range(4)); ax.set_xticklabels(labs,fontsize=6.8)
ax.set_ylabel("myeloid x lymphoid correlation",fontsize=7.2); ax.set_ylim(0,1)
ax.set_title("E  Myeloid-lymphoid coupling is weaker\nwithin tissue than in bulk",fontsize=8,loc="left")
ax.tick_params(labelsize=6.6); ax.spines[["top","right"]].set_visible(False)

# F — every tumour contains every ecotype
ax=fig.add_subplot(gs[1,3])
bot=np.zeros(3)
for k in ORDER:
    v=E[k.split("-")[0][:4]].values/E.n_spots.values*100
    ax.bar(range(3),v,bottom=bot,color=COL[k],width=.6,label=k); bot+=v
ax.set_xticks(range(3)); ax.set_xticklabels([f"Sample-{s}" for s in E["sample"]],fontsize=7)
ax.set_ylabel("% of in-tissue spots",fontsize=7.2); ax.set_ylim(0,100)
ax.set_title("F  Every section contains all three\necotypes",fontsize=8,loc="left")
ax.tick_params(labelsize=6.6); ax.spines[["top","right"]].set_visible(False)

fig.suptitle("Spatial organisation of the immune programmes in H3.3 K27M paediatric diffuse midline glioma (GSE268577)",
             fontsize=9.5,y=.985)
for e in ("png","pdf"): fig.savefig(f"FigureS21_spatial.{e}",dpi=300,bbox_inches="tight",facecolor="white")
print("saved FigureS21")
