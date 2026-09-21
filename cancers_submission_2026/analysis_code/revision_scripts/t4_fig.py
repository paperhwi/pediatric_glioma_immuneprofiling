import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg"); import matplotlib.pyplot as plt
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.linewidth":0.8,
                     "pdf.fonttype":42,"ps.fonttype":42})
s=pd.read_csv("T4_HSPC_ssGSEA_scores.tsv",sep="\t",index_col=0)
kw=pd.read_csv("T4_KW_by_ecotype.tsv",sep="\t").set_index("module")
rs=pd.read_csv("T4_residualised_specificity.tsv",sep="\t").set_index("module")
ORDER=["Lymphocyte-inflamed","Myeloid-dominant","Immune-desert"]
COL={"Lymphocyte-inflamed":"#3B82F6","Myeloid-dominant":"#EF4444","Immune-desert":"#7C8798"}
MODS=["HSC","LMPP","GMP","CLP","EMP (lineage control)","CD14 monocyte (positive control)"]
LAB=["HSC","LMPP","GMP","CLP","EMP\n(lineage control)","CD14 monocyte\n(positive control)"]

fig=plt.figure(figsize=(12.6,5.4),dpi=300)
gs=fig.add_gridspec(2,6,height_ratios=[1,.92],hspace=.62,wspace=.42)
for i,(m,lb) in enumerate(zip(MODS,LAB)):
    a=fig.add_subplot(gs[0,i])
    data=[s.loc[s.ecotype==e,m].values for e in ORDER]
    bp=a.boxplot(data,widths=.6,patch_artist=True,showfliers=False,
                 medianprops=dict(color="k",lw=1.0),whiskerprops=dict(lw=.7),capprops=dict(lw=.7))
    for p,e in zip(bp["boxes"],ORDER): p.set_facecolor(COL[e]); p.set_alpha(.55); p.set_edgecolor(COL[e]); p.set_lw(.8)
    a.set_xticks([1,2,3]); a.set_xticklabels(["LI","MD","ID"],fontsize=7)
    q=kw.loc[m,"q_BH"]; eps=kw.loc[m,"epsilon_sq"]
    a.set_title(f"{lb}\n$\\epsilon^2$ = {eps:.3f}   q = {q:.1e}",fontsize=7)
    if i==0: a.set_ylabel("ssGSEA NES")
    a.spines[["top","right"]].set_visible(False); a.tick_params(labelsize=7)

# bottom-left: effect size before vs after residualising
a=fig.add_subplot(gs[1,0:3])
order2=["HSC","LMPP","CLP","EMP (lineage control)","GMP"]
x=np.arange(len(order2)); w=.36
a.bar(x-w/2,[rs.loc[m,"eps_raw"] for m in order2],w,color="#B91C1C",label="unadjusted")
a.bar(x+w/2,[rs.loc[m,"eps_resid"] for m in order2],w,color="#94A3B8",
      label="residualised on mature\nmonocyte + macrophage content")
a.set_xticks(x); a.set_xticklabels(["HSC","LMPP","CLP","EMP","GMP"],fontsize=7.5)
a.set_ylabel("Kruskal-Wallis $\\epsilon^2$ across ecotypes")
a.set_title("F  Progenitor-like signal is not separable from mature myeloid content\n"
            "(no module survives adjustment; all q $\\geq$ 0.18)",fontsize=7.5,loc="left")
a.legend(fontsize=6.3,frameon=False,loc="upper left")
a.annotate("0.205 $\\rightarrow$ 0.007",xy=(4+w/2,.012),xytext=(2.55,.175),fontsize=6.8,color="#B91C1C",
           arrowprops=dict(arrowstyle="->",color="#B91C1C",lw=.8))
a.spines[["top","right"]].set_visible(False); a.tick_params(labelsize=7)

# bottom-right: GMP vs monocyte scatter
a=fig.add_subplot(gs[1,3:6])
for e in ORDER:
    sub=s[s.ecotype==e]
    a.scatter(sub["CD14 monocyte (positive control)"],sub["GMP"],s=9,c=COL[e],alpha=.65,lw=0,label=e)
r=np.corrcoef(s["CD14 monocyte (positive control)"],s["GMP"])[0,1]
a.set_xlabel("CD14 monocyte NES (mature myeloid content)"); a.set_ylabel("GMP NES")
a.set_title(f"G  GMP tracks mature monocyte content (Pearson r = {r:.2f})\n"
            "ecotypes separate along the shared axis, not across it",fontsize=7.5,loc="left")
a.legend(fontsize=6.3,frameon=False,markerscale=1.3,loc="upper left")
a.spines[["top","right"]].set_visible(False); a.tick_params(labelsize=7)

fig.suptitle("Haematopoietic stem and progenitor programs across immune ecotypes "
             "(Azimuth 2023 human bone-marrow reference; n = 349)",fontsize=9,y=.99)
for e in ("png","pdf"): fig.savefig(f"FigureS_HSPC_programs.{e}",dpi=300,bbox_inches="tight")
print("saved; GMP~monocyte r =",round(r,3))
