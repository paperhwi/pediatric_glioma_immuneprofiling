"""Figure 6: molecular-group-adjusted differential expression between immune states (n = 343)."""
import sys, os, json; sys.path.insert(0,"/home/claude/q343/src"); from common import *
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.linewidth":0.8,"pdf.fonttype":42,"svg.fonttype":"none"})
BLUE,RED,GREY,DARK,NAVY,MAROON,TEAL="#3B82F6","#EF4444","#CBD5E1","#7C8798","#1E3A8A","#7F1D1D","#0F766E"
A,B,C="Immune-high vs Immune-low","Immune-intermediate vs Immune-low","Immune-high vs Immune-intermediate"
d=pd.read_csv(ROOT+"out/T_adjusted_DEG_n343_all_contrasts.tsv",sep="\t")
w=d.pivot(index="gene",columns="contrast",values="adjusted_log2FC")[[A,B,C]].dropna(); q=d.pivot(index="gene",columns="contrast",values="q_BH").loc[w.index,[A,B,C]]; NG=len(w)
def sets(c): return (q[c]<0.05)&(w[c]>1),(q[c]<0.05)&(w[c]<-1)
fig,axs=plt.subplots(2,3,figsize=(11.6,7.4)); axs=axs.ravel()
def lab(ax,t): ax.text(-0.16,1.08,t,transform=ax.transAxes,fontsize=13,weight="bold")
counts={}
for i,(c,ref,col) in enumerate([(A,LO,BLUE),(B,LO,RED),(C,MID,BLUE)]):
    a=axs[i]; up,dn=sets(c); counts[c]=(int(up.sum()),int(dn.sum())); x=w[c].values; y=-np.log10(np.clip(q[c].values,1e-300,None)); oth=~(up|dn).values
    a.scatter(x[oth],y[oth],s=1.5,c=GREY,rasterized=True,lw=0); a.scatter(x[dn.values],y[dn.values],s=2.4,c=DARK,rasterized=True,lw=0); a.scatter(x[up.values],y[up.values],s=2.4,c=col,rasterized=True,lw=0)
    for v in (-1,1): a.axvline(v,color="k",ls="--",lw=.6)
    a.axhline(-np.log10(.05),color="k",ls="--",lw=.6)
    a.set_title(f"{c}\nreference = {ref}   |   up {counts[c][0]:,}; down {counts[c][1]:,}",fontsize=7.8)
    a.set_xlabel("Molecular-group-adjusted log$_2$ fold change",fontsize=7.8); a.set_ylabel("$-$log$_{10}$(BH $q$)",fontsize=7.8); lab(a,"ABC"[i])
upA,_=sets(A); upB,_=sets(B); sA,sB=set(w.index[upA]),set(w.index[upB]); shared=sorted(sA&sB); nS,nA,nB=len(shared),len(sA),len(sB)
a=axs[3]; ypos=[1.0,0.0]; bh_=0.42
a.barh(ypos[0],nS,height=bh_,color=NAVY); a.barh(ypos[0],nA-nS,left=nS,height=bh_,color=BLUE,alpha=.45)
a.barh(ypos[1],nS,height=bh_,color=NAVY); a.barh(ypos[1],nB-nS,left=nS,height=bh_,color=RED,alpha=.75)
a.set_yticks(ypos); a.set_yticklabels(["Immune-high\nvs Immune-low","Immune-intermediate\nvs Immune-low"],fontsize=7.4)
a.set_xlabel("Up-regulated genes (BH $q$ < 0.05, adj. log$_2$FC > 1)",fontsize=7.8); a.set_xlim(0,nA*1.30); a.set_ylim(-0.62,1.62)
a.text(nS/2,ypos[0],f"{nS}",ha="center",va="center",fontsize=7.4,color="white",weight="bold"); a.text(nS+(nA-nS)/2,ypos[0],f"{nA-nS:,} unique",ha="center",va="center",fontsize=7.4,color=NAVY,weight="bold")
a.text(nS/2,ypos[1],f"{nS}",ha="center",va="center",fontsize=7.4,color="white",weight="bold")
a.annotate(f"{nB-nS} unique",xy=(nB,ypos[1]+bh_/2),xytext=(nB+nA*0.09,ypos[1]+0.46),fontsize=7.2,color=MAROON,arrowprops=dict(arrowstyle="->",color=MAROON,lw=0.8))
a.text(nA*1.02,ypos[0],f"total {nA:,}",ha="left",va="center",fontsize=7.2,color=NAVY); a.text(nB+nA*0.02,ypos[1],f"total {nB}",ha="left",va="center",fontsize=7.2,color=MAROON)
a.set_title(f"{nS}/{nB} ({nS/nB*100:.1f}%) of the Immune-intermediate set\nis contained in the Immune-high set",fontsize=7.8)
a.spines[["left"]].set_visible(False); a.tick_params(axis="y",length=0); lab(a,"D")
a=axs[4]; x,y=w[B].values,w[A].values; sig=np.array([g in sA or g in sB for g in w.index])
a.scatter(x[~sig],y[~sig],s=2,c=GREY,alpha=.35,lw=0,rasterized=True); a.scatter(x[sig],y[sig],s=3,c=BLUE,alpha=.45,lw=0,rasterized=True)
lim=[-3.0,5.0]; a.plot(lim,lim,"k--",lw=.8,label="identity ($y = x$)")
xc,yc=x-x.mean(),y-y.mean(); sxx,syy,sxy=(xc**2).mean(),(yc**2).mean(),(xc*yc).mean(); sl=((syy-sxx)+np.sqrt((syy-sxx)**2+4*sxy**2))/(2*sxy); ic=y.mean()-sl*x.mean(); r=stats.pearsonr(x,y).statistic
xs=np.linspace(*lim,10); a.plot(xs,sl*xs+ic,color=TEAL,lw=1.3,label=f"orthogonal fit (slope = {sl:.2f})"); a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel("adj. log$_2$FC   Immune-intermediate vs Immune-low",fontsize=7.8); a.set_ylabel("adj. log$_2$FC   Immune-high vs Immune-low",fontsize=7.8)
a.set_title(f"All {NG:,} genes; Pearson $r$ = {r:.2f}\nslope above identity is consistent with an ordered\nshared immune-presence axis",fontsize=7.8,loc="left"); a.legend(fontsize=6.4,frameon=False,loc="upper left"); lab(a,"E")
a=axs[5]; sub=w.loc[shared]; parts=a.violinplot([sub[B].values,sub[A].values],positions=[0,1],widths=.75,showextrema=False)
for pc,cc in zip(parts["bodies"],[RED,BLUE]): pc.set_facecolor(cc); pc.set_alpha(.45); pc.set_edgecolor(cc)
a.plot([0,1],[sub[B].median(),sub[A].median()],"ko-",ms=4,lw=1.2); a.set_xticks([0,1]); a.set_xticklabels(["Immune-intermediate\nvs Immune-low","Immune-high\nvs Immune-low"],fontsize=7.4)
a.set_ylabel(f"adj. log$_2$FC on the {nS} shared genes",fontsize=7.8); wil=stats.wilcoxon(sub[A].values,sub[B].values); ratio=float(np.median(sub[A]/sub[B]))
m_,e_=f"{wil.pvalue:.1e}".split("e"); a.set_title(f"Median {ratio:.2f}$\\times$ larger in Immune-high\npaired Wilcoxon signed-rank $P$ = {m_} × 10$^{{{int(e_)}}}$",fontsize=7.8,loc="left"); lab(a,"F")
for a in axs: a.spines[["top","right"]].set_visible(False); a.tick_params(labelsize=7)
fig.tight_layout(h_pad=2.5,w_pad=2.0); os.makedirs(ROOT+"fig/Figure6",exist_ok=True)
for ext in ("svg","pdf","png"): fig.savefig(ROOT+f"fig/Figure6/Figure6.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight",facecolor="white")
pd.Series(shared,name="gene").to_csv(ROOT+f"out/T_shared_genes_{nS}.tsv",sep="\t",index=False)
S=dict(genes_tested=NG,up_high_vs_low=nA,up_intermediate_vs_low=nB,shared=nS,unique_high=nA-nS,unique_intermediate=nB-nS,pct_intermediate_in_high=100*nS/nB,
       orthogonal_slope=sl,pearson_r=r,median_ratio_shared=ratio,wilcoxon_two_sided_p=wil.pvalue,counts=counts)
json.dump(S,open(ROOT+"out/T_fig6_summary.json","w"),indent=1,default=float); print(S)
with pd.ExcelWriter(ROOT+"fig/Figure6/Figure6_source_data.xlsx") as wr:
    d.to_excel(wr,sheet_name="adjusted_DEG_all_contrasts",index=False); pd.Series(shared,name="gene").to_frame().to_excel(wr,sheet_name="shared_genes",index=False); pd.DataFrame([S]).astype(str).T.to_excel(wr,sheet_name="summary")
