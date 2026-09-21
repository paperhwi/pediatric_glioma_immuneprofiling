import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.linewidth":0.8,
                     "pdf.fonttype":42,"ps.fonttype":42})
A="Lymphocyte_inflamed_vs_Immune_desert"; B="Myeloid_dominant_vs_Immune_desert"; C="Lymphocyte_inflamed_vs_Myeloid_dominant"
w=pd.read_pickle("T2_w.pkl"); q=pd.read_pickle("T2_q.pkl")
BLUE,RED,GREY,DARK="#3B82F6","#EF4444","#CBD5E1","#7C8798"

# ---------- Figure 6A revised: explicit reference in every panel title ----------
panels=[(A,"Lymphocyte-inflamed  vs  Immune-desert",BLUE),
        (B,"Myeloid-dominant  vs  Immune-desert",RED),
        (C,"Lymphocyte-inflamed  vs  Myeloid-dominant",BLUE)]
fig,ax=plt.subplots(1,3,figsize=(11.4,3.5),dpi=300)
for i,(c,title,col) in enumerate(panels):
    x=w[c].values; y=-np.log10(np.clip(q[c].values,1e-300,None))
    up=(q[c]<0.05)&(w[c]>1); dn=(q[c]<0.05)&(w[c]<-1)
    a=ax[i]
    a.scatter(x[~(up|dn)],y[~(up|dn)],s=1.5,c=GREY,rasterized=True,lw=0)
    a.scatter(x[up],y[up],s=2.2,c=col,rasterized=True,lw=0)
    a.scatter(x[dn],y[dn],s=2.2,c=DARK,rasterized=True,lw=0)
    for v in (-1,1): a.axvline(v,color="k",ls="--",lw=.6)
    a.axhline(-np.log10(.05),color="k",ls="--",lw=.6)
    ref=title.split("  vs  ")[1]
    a.set_title(f"{title}\nreference = {ref}   |   up {int(up.sum())}; down {int(dn.sum())}",
                fontsize=7.5)
    a.set_xlabel("Molecular-group-adjusted log$_2$ fold change")
    if i==0: a.set_ylabel("$-$log$_{10}$(BH q)")
    a.spines[["top","right"]].set_visible(False); a.tick_params(labelsize=7)
fig.suptitle("Molecular-group-adjusted differential expression   log$_2$(TPM+1) ~ ecotype + integrated molecular group",
             fontsize=8.5,y=1.04)
fig.tight_layout()
for e in ("png","pdf"): fig.savefig(f"Figure6A_revised.{e}",dpi=300,bbox_inches="tight")

# ---------- Figure 6C new: the nesting / gradient evidence ----------
sA=set(w.index[(q[A]<0.05)&(w[A]>1)]); sB=set(w.index[(q[B]<0.05)&(w[B]>1)])
sh=sorted(sA&sB); both=sorted(sA|sB)
fig2,ax2=plt.subplots(1,3,figsize=(11.4,3.4),dpi=300)

# C1 nested Euler
a=ax2[0]; a.set_aspect("equal"); a.axis("off")
a.add_patch(Circle((0,0),1.0,fc=BLUE,alpha=.30,ec=BLUE,lw=1.1))
a.add_patch(Circle((.80,0),.34,fc=RED,alpha=.45,ec=RED,lw=1.1))
a.text(-.62,.55,f"Lymphocyte-inflamed\nvs Immune-desert\n{len(sA):,}",ha="center",fontsize=7,color="#1E3A8A")
a.text(.82,-.62,f"Myeloid-dominant\nvs Immune-desert\n{len(sB):,}",ha="center",fontsize=7,color="#7F1D1D")
a.text(.74,0,f"{len(sA&sB)}",ha="center",va="center",fontsize=9,weight="bold")
a.text(1.10,.30,f"{len(sB-sA)}",ha="center",va="center",fontsize=7.5,color="#7F1D1D")
a.text(-.30,0,f"{len(sA-sB):,}",ha="center",va="center",fontsize=9,weight="bold",color="#1E3A8A")
a.set_xlim(-1.15,1.30); a.set_ylim(-1.15,1.25)
a.set_title(f"A  {len(sA&sB)}/{len(sB)} ({len(sA&sB)/len(sB)*100:.0f}%) of Myeloid-dominant genes\n"
            f"are nested inside the Lymphocyte-inflamed set",fontsize=7.5,loc="left")

# C2 effect-size scatter
a=ax2[1]
allg=w.index
x=w.loc[allg,B].values; y=w.loc[allg,A].values
issig=np.array([g in sA or g in sB for g in allg])
a.scatter(x[~issig],y[~issig],s=2,c=GREY,alpha=.35,lw=0,rasterized=True)
a.scatter(x[issig],y[issig],s=3,c=BLUE,alpha=.45,lw=0,rasterized=True)
lim=[-3.0,5.0]; a.plot(lim,lim,"k--",lw=.8,label="identity (y = x)")
xc,yc=x-x.mean(),y-y.mean(); sxx,syy,sxy=(xc**2).mean(),(yc**2).mean(),(xc*yc).mean()
sl=((syy-sxx)+np.sqrt((syy-sxx)**2+4*sxy**2))/(2*sxy); ic=y.mean()-sl*x.mean()
xs=np.linspace(*lim,10); a.plot(xs,sl*xs+ic,color="#0F766E",lw=1.2,label=f"orthogonal fit (slope = {sl:.2f})")
a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel("adj. log$_2$FC   Myeloid-dominant vs Immune-desert")
a.set_ylabel("adj. log$_2$FC   Lymphocyte-inflamed vs Immune-desert")
a.set_title(f"B  Same axis, larger amplitude (all {len(allg):,} genes)\nslope > 1 means the Desert contrast is a shared gradient",fontsize=7.5,loc="left")
a.legend(fontsize=6.2,frameon=False,loc="upper left")
a.spines[["top","right"]].set_visible(False); a.tick_params(labelsize=7)

# C3 paired effect sizes on the shared genes
a=ax2[2]
sub=w.loc[sh]
parts=a.violinplot([sub[B].values,sub[A].values],positions=[0,1],widths=.75,showextrema=False)
for pc,cc in zip(parts["bodies"],[RED,BLUE]): pc.set_facecolor(cc); pc.set_alpha(.45); pc.set_edgecolor(cc)
a.plot([0,1],[sub[B].median(),sub[A].median()],"ko-",ms=4,lw=1.2)
a.set_xticks([0,1]); a.set_xticklabels(["Myeloid-dominant\nvs Desert","Lymphocyte-inflamed\nvs Desert"],fontsize=7)
a.set_ylabel("adj. log$_2$FC on the 426 shared genes")
a.set_title(f"C  Median {np.median(sub[A]/sub[B]):.2f}$\\times$ larger in Lymphocyte-inflamed\n"
            f"paired Wilcoxon P = 7.6 $\\times$ 10$^{{-71}}$",fontsize=7.5,loc="left")
a.spines[["top","right"]].set_visible(False); a.tick_params(labelsize=7)

fig2.tight_layout(w_pad=2.0)
for e in ("png","pdf"): fig2.savefig(f"Figure6C_new_contrast_structure.{e}",dpi=300,bbox_inches="tight")
print("saved both figures")
