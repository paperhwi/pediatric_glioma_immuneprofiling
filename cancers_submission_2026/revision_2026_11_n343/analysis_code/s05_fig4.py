"""Figure 4 and Tables 2-4 on n = 343 (immune state x location / molecular group; sequential PERMANOVA).
Statistics functions imported unchanged from Rev344_09 fig4_344.py (which reproduces the submitted n = 349 Figure 4)."""
import sys, os, json; sys.path.insert(0,"/home/claude/q343/src"); from common import *
import importlib.util
spec=importlib.util.spec_from_file_location("f4","/home/claude/w/c09/fig4_344.py")
src=open("/home/claude/w/c09/fig4_344.py").read().replace('sys.path.insert(0, "/home/claude/q")','').replace("from common import load, ECO","ECO=None")
f4={}; exec(compile(src.split("# ---------------- revised figure ----------------")[0].split('if __name__ == "__main__":')[0],"fig4_344","exec"),f4)
location,ct_stats,permanova_seq,Z34,H=f4["location"],f4["ct_stats"],f4["permanova_seq"],f4["Z34"],f4["H"]
a=load(); d=a[["state","cohort_group"]].copy()
d["location_class"]=[location(H.at[b,"CNS_region"],H.at[b,"primary_site"]) for b in d.index]
L=["Midline","Hemispheric","Posterior_fossa","Ambiguous"]; G=["DMG_K27","DHG_G34","pHGG_WT","IHG"]
I={}
dl=d[d.location_class.isin(L)]; I["A"]=pd.crosstab(dl.location_class,dl.state).reindex(index=L,columns=ORD).fillna(0).astype(int)
I["B"]=pd.crosstab(d.cohort_group,d.state).reindex(index=G,columns=ORD).fillna(0).astype(int)
dw=d[(d.cohort_group=="pHGG_WT")&d.location_class.isin(["Midline","Hemispheric"])]; I["C"]=pd.crosstab(dw.location_class,dw.state).reindex(index=["Midline","Hemispheric"],columns=ORD).fillna(0).astype(int)
ids=d[d.location_class.isin(["Midline","Hemispheric","Posterior_fossa"])&d.cohort_group.isin(G)&d.index.isin(Z34.index)].index
R={k:ct_stats(I[k],B=9999) for k in "ABC"}
from scipy.spatial.distance import pdist, squareform
D=squareform(pdist(Z34.loc[ids].values)); mol=d.loc[ids,"cohort_group"].values.astype(str); loc=d.loc[ids,"location_class"].values.astype(str)
R["D_mol_first"]=permanova_seq(D,mol,loc,B=4999); R["D_loc_first"]=permanova_seq(D,loc,mol,B=4999)
rows=[]
for k,nm in [("A","Anatomical location (4 classes) x immune state"),("B","Molecular group x immune state"),("C","pHGG H3/IDH-wildtype: midline vs hemispheric x immune state")]:
    rows.append(dict(panel=k,analysis=nm,**R[k]))
for kk,(f,s_) in [("D_mol_first",("molecular group","location")),("D_loc_first",("location","molecular group"))]:
    t=R[kk]; rows.append(dict(panel="D",analysis=f"PERMANOVA {f} (first)",n=t["n"],partial_R2=t["R2_first"],F=t["F_first"],df=t["df_first"],p_perm=t["P_first"]))
    rows.append(dict(panel="D",analysis=f"PERMANOVA {s_} | {f}",n=t["n"],partial_R2=t["R2_second"],F=t["F_second"],df=t["df_second"],p_perm=t["P_second"]))
T=pd.DataFrame(rows); T.to_csv(ROOT+"out/T_fig4_statistics.tsv",sep="\t",index=False)
for k in "ABC": I[k].to_csv(ROOT+f"out/T_fig4{k}_counts.tsv",sep="\t")
d.to_csv(ROOT+"out/T_fig4_per_tumour_location.tsv",sep="\t")
pd.set_option("display.width",250); print(T.round(4).to_string()); [print(I[k]) for k in "ABC"]
# ---------------- figure
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"svg.fonttype":"none","pdf.fonttype":42,"font.family":"DejaVu Sans"})
fig,axes=plt.subplots(1,4,figsize=(18,5.6),gridspec_kw={"width_ratios":[1.2,1.1,.85,1.05]})
def stacked(ax,ct,s,title,labels):
    pct=ct.div(ct.sum(axis=1),axis=0)*100; bottom=np.zeros(len(ct))
    for e in ORD:
        vals=pct[e].values; ax.bar(np.arange(len(ct)),vals,bottom=bottom,color=COL[e],edgecolor="white",label=SHORT[e])
        for i,(v,b,n) in enumerate(zip(vals,bottom,ct[e])):
            if v>=7: ax.text(i,b+v/2,str(int(n)),ha="center",va="center",fontsize=8,color="white",weight="bold")
        bottom+=vals
    ax.set_xticks(np.arange(len(ct)),[f"{labels.get(x,x)}\n(n={int(ct.loc[x].sum())})" for x in ct.index],fontsize=8)
    ax.set_ylim(0,100); ax.set_ylabel("Tumours (%)"); ax.set_title(title,weight="bold",fontsize=10)
    pm="P_MC<0.0001" if s['p_montecarlo']<1e-4 else f"P_MC={s['p_montecarlo']:.4g}"
    ax.text(.5,-.25,f"χ²={s['chi2']:.2f}, df={s['df']}, {pm}\nCramér's V={s['cramers_V']:.2f}; n={s['n']}",transform=ax.transAxes,ha="center",va="top",fontsize=8)
stacked(axes[0],I["A"],R["A"],"A  Immune state by anatomical location",{"Posterior_fossa":"Posterior\nfossa","Ambiguous":"Multi-compartment/\nambiguous"})
stacked(axes[1],I["B"],R["B"],"B  Immune state by molecular group",{"DMG_K27":"DMG\nH3 K27-alt.","DHG_G34":"DHG\nH3 G34-mut.","pHGG_WT":"pHGG\nH3/IDH-wt","IHG":"IHG"})
stacked(axes[2],I["C"],R["C"],"C  Within H3/IDH-wt pHGG",{})
m,l=R["D_mol_first"],R["D_loc_first"]; vals=[m["R2_first"],m["R2_second"],l["R2_first"],l["R2_second"]]; pv=[m["P_first"],m["P_second"],l["P_first"],l["P_second"]]
ax=axes[3]; ax.bar(np.arange(4),vals,color=["#7C3AED","#059669","#10B981","#A78BFA"])
ax.set_xticks(np.arange(4),["Molecular\nfirst","Location\n| molecular","Location\nfirst","Molecular\n| location"],fontsize=8)
ax.set_ylabel("PERMANOVA partial R²"); ax.set_title(f"D  Sequential PERMANOVA\n(34-feature TPM matrix, n={m['n']})",weight="bold",fontsize=10); ax.set_ylim(0,max(vals)*1.42)
for i,(v,p) in enumerate(zip(vals,pv)): ax.text(i,v+.0015,f"{100*v:.1f}%\nP={p:.4g}",ha="center",fontsize=7.5)
for x in axes: x.spines[["top","right"]].set_visible(False)
h,lab=axes[0].get_legend_handles_labels(); fig.legend(h,lab,loc="upper center",bbox_to_anchor=(.5,1.02),ncol=3,frameon=False,title="Immune state")
fig.tight_layout(rect=[0,.08,1,.95]); os.makedirs(ROOT+"fig/Figure4",exist_ok=True)
for ext in ("svg","pdf","png"): fig.savefig(ROOT+f"fig/Figure4/Figure4.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight",facecolor="white")
