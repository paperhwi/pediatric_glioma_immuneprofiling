import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.linewidth":0.8,
                     "xtick.major.width":0.8,"ytick.major.width":0.8,"pdf.fonttype":42,"ps.fonttype":42})
sw=pd.read_csv("T1_ksweep_k2_k10.tsv",sep="\t")
gp=pd.read_csv("T1_gap_statistic.tsv",sep="\t")
BLUE,ORANGE,RED,GREY="#2563EB","#F97316","#B91C1C","#64748B"

fig,ax=plt.subplots(1,5,figsize=(13.2,2.75),dpi=300)

def mark(a):
    a.axvline(2,color=GREY,ls=":",lw=0.8); a.axvline(3,color=RED,ls="--",lw=0.9)
    a.set_xlabel("Number of clusters (k)"); a.set_xticks(range(1,11))
    a.spines[["top","right"]].set_visible(False); a.tick_params(labelsize=7)

ax[0].plot(sw.k,sw.PAC,"o-",color=BLUE,ms=3.5,lw=1.2)
ax[0].set_ylabel("PAC (lower is better)"); ax[0].set_title("A  Proportion of ambiguous\nclustering",fontsize=8,loc="left")
mark(ax[0])

ax[1].plot(sw.k,sw.silhouette_consensus,"o-",color=ORANGE,ms=3.5,lw=1.2)
ax[1].set_ylabel("Silhouette (1 - consensus)"); ax[1].set_title("B  Silhouette width",fontsize=8,loc="left")
mark(ax[1])

_d=sw.dropna(subset=["delta_AUC"])
ax[2].bar(_d.k,_d.delta_AUC,color=[RED if k==3 else "#94A3B8" for k in _d.k],width=.6)
ax[2].set_ylabel("Relative $\\Delta$AUC of consensus CDF")
ax[2].set_title("C  Consensus CDF elbow",fontsize=8,loc="left")
mark(ax[2]); ax[2].set_xlim(2.3,10.7); ax[2].set_xticks(range(3,11))
ax[2].annotate("last substantial\ngain (+29.8%)",xy=(3,sw.loc[sw.k==3,"delta_AUC"].iloc[0]),
               xytext=(5.2,0.22),fontsize=6.5,color=RED,
               arrowprops=dict(arrowstyle="->",color=RED,lw=0.8))

ax[3].errorbar(gp.k,gp.gap,yerr=gp.s_k,fmt="o-",color="#0F766E",ms=3.5,lw=1.2,capsize=2,elinewidth=0.8)
ax[3].set_ylabel("Gap statistic"); ax[3].set_title("D  Gap statistic (k = 1-10)",fontsize=8,loc="left")
mark(ax[3])
ax[3].annotate("1-SE rule\nselects k = 3",xy=(3,gp.loc[gp.k==3,"gap"].iloc[0]),xytext=(5.0,1.525),
               fontsize=6.5,color=RED,arrowprops=dict(arrowstyle="->",color=RED,lw=0.8))

ax[4].plot(sw.k,sw.min_cluster_size,"o-",color=RED,ms=3.5,lw=1.2)
ax[4].axhline(10,color=GREY,ls="--",lw=0.8)
ax[4].text(7.6,12,"n = 10",fontsize=6.5,color=GREY)
ax[4].set_yscale("log"); ax[4].set_ylabel("Smallest cluster size (n)")
ax[4].set_title("E  Cluster-size floor",fontsize=8,loc="left")
mark(ax[4])

fig.tight_layout(w_pad=1.6)
for ext in ("png","pdf"):
    fig.savefig(f"/tmp/claude-0/w1/Figure3A_revised.{ext}",dpi=300,bbox_inches="tight")
print("saved")

summary=pd.DataFrame({
 "criterion":["PAC (minimum)","Silhouette on 1-consensus (maximum)","Calinski-Harabasz (maximum)",
              "Davies-Bouldin (minimum)","Consensus CDF delta-AUC elbow","Gap statistic, Tibshirani 1-SE rule (k=1-10)",
              "Cluster-size floor >= 10 samples"],
 "selected_k":[int(sw.loc[sw.PAC.idxmin(),"k"]),int(sw.loc[sw.silhouette_consensus.idxmax(),"k"]),
               int(sw.loc[sw.calinski_harabasz.idxmax(),"k"]),int(sw.loc[sw.davies_bouldin.idxmin(),"k"]),
               3,3,int(sw.loc[sw.min_cluster_size>=10,"k"].max())],
 "family":["separation/stability","separation/stability","separation/stability","separation/stability",
           "number-of-clusters selection","number-of-clusters selection","interpretability constraint"]})
summary.to_csv("/tmp/claude-0/w1/T1_criterion_summary.tsv",sep="\t",index=False)
print(summary.to_string(index=False))
