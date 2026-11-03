"""Supplementary survival figure: A sensitivity analyses of the Immune-high vs Immune-low hazard ratio (all reported regardless of
significance); B-C Kaplan-Meier within DMG and pHGG-WT; D immune-state HR within subtype; E incremental prognostic value."""
import sys, os; sys.path.insert(0,"/home/claude/q343/src"); from common import *
from s08_fig7 import km, lab, plt, NullLocator, fmtp
from matplotlib.gridspec import GridSpec
a=load(); s=a[a.os_eval]
S=pd.read_csv(ROOT+"out/T_survival_sensitivity.tsv",sep="\t")
fig=plt.figure(figsize=(7.6,12.4)); gs=GridSpec(3,2,height_ratios=[1.35,1,0.9],hspace=0.85,wspace=0.45,left=0.10,right=0.97,top=0.97,bottom=0.05)
ax=fig.add_subplot(gs[0,:]); p0=ax.get_position(); ax.set_position([0.40,p0.y0,0.33,p0.height]); y=np.arange(len(S))[::-1]
for yi,(_,r) in zip(y,S.iterrows()):
    prim=r.analysis.startswith("Primary"); col="#1D4ED8" if prim else ("#B91C1C" if r.analysis[:3] in ("S1 ","S2 ","S3 ","S4 ","S4b") else "#374151")
    ax.plot([r.high_lo,r.high_hi],[yi,yi],color=col,lw=1.3); ax.plot(r.high_HR,yi,"D" if prim else "s",color=col,ms=4.8)
    ax.text(2.3,yi,f"{r.n} / {r.deaths}    {r.high_HR:.2f} ({r.high_lo:.2f}–{r.high_hi:.2f})    P = {fmtp(r.high_p)}",va="center",fontsize=6)
ax.text(2.3,len(S)-0.35,"n / deaths    HR (95% CI)",fontsize=6,weight="bold")
ax.axvline(1,ls="--",lw=0.6,color="k"); ax.set_xscale("log"); ax.set_xlim(0.25,12); ax.spines["bottom"].set_bounds(0.25,2); ax.set_xticks([0.25,0.5,1,2]); ax.set_xticklabels(["0.25","0.5","1","2"]); ax.xaxis.set_minor_locator(NullLocator())
ax.set_yticks(y); ax.set_yticklabels([x.replace(" (pre-specified): corrected OS, molecular group + age + sex"," (pre-specified)").replace(" (reference: poly-A stranded)","").replace("Excluding all 13 records whose OS origin precedes the sequenced-HGG diagnosis","Excluding 13 records with OS origin before HGG diagnosis") for x in S.analysis],fontsize=6.3); ax.tick_params(axis="y",length=0)
ax.set_xlabel("Hazard ratio, Immune-high vs Immune-low (95% CI)"); ax.set_title("Primary model and sensitivity analyses",fontsize=7.6,loc="left",pad=14); ax.text(0,-0.16,"Blue: pre-specified primary model. Red: sensitivity analyses requested for the response (autopsy, time origin, library).",transform=ax.transAxes,fontsize=5.6,color="#374151"); lab(ax,"A",-1.15)
for j,(g,title,xm,tk) in enumerate([("DMG_K27","DMG, H3 K27-altered",5,[0,1,2,3,4,5]),("pHGG_WT","pHGG, H3/IDH-wildtype",10,[0,2,4,6,8,10])]):
    ax=fig.add_subplot(gs[1,j]); km(ax,s[s.cohort_group==g],xm,tk,title,fs=5.2); lab(ax,"BC"[j])
ax=fig.add_subplot(gs[2,0]); E=pd.read_csv(ROOT+"out/T_within_subtype_state_cox.tsv",sep="\t"); ip=pd.read_csv(ROOT+"out/T_within_subtype_interaction.tsv",sep="\t").p.iloc[0]
ent=[("pHGG_WT","Immune-high vs Immune-low"),("pHGG_WT","Immune-intermediate vs Immune-low"),("DMG_K27","Immune-high vs Immune-low"),("DMG_K27","Immune-intermediate vs Immune-low")]
lg={"pHGG_WT":"pHGG-WT","DMG_K27":"DMG"}; y=np.arange(len(ent))[::-1]
for yi,(g,t) in zip(y,ent):
    rr=E[(E.group==g)&(E.term==t)].iloc[0]; c=COL[HI if t.startswith("Immune-high") else MID]
    ax.plot([rr.lo,rr.hi],[yi,yi],color=c,lw=1.2); ax.plot(rr.HR,yi,"s",color=c,ms=4.5); ax.text(2.6,yi,f"{rr.HR:.2f} ({rr.lo:.2f}–{rr.hi:.2f})\nP = {fmtp(rr.p)}",va="center",fontsize=5.6)
ax.axvline(1,ls="--",lw=0.6,color="k"); ax.set_xscale("log"); ax.set_xlim(0.15,9); ax.spines["bottom"].set_bounds(0.15,2.5); ax.set_xticks([0.25,0.5,1,2]); ax.set_xticklabels(["0.25","0.5","1","2"]); ax.xaxis.set_minor_locator(NullLocator())
ax.set_yticks(y); ax.set_yticklabels([f"{lg[g]}: {t.replace(' vs ',chr(10)+'vs ')}" for g,t in ent],fontsize=5.6); ax.set_xlabel("HR within subtype (age, sex adjusted)")
ax.set_title(f"Immune state within molecular subtype\nstate × subtype interaction P = {ip:.2f}",fontsize=7.4); lab(ax,"D",-0.6)
ax=fig.add_subplot(gs[2,1]); Q=pd.read_csv(ROOT+"out/T_incremental_value.tsv",sep="\t").reset_index(drop=True)
k=Q.index[(Q.analysis.str.startswith("Primary"))&(Q.model=="M1 M0+immune state")][0]; b,n=Q.loc[k-1],Q.loc[k]
mets=[("Harrell_C_optimism_corrected","C-index\n(optimism-\ncorrected)"),("Uno_C_tau5y","Uno C\n(τ = 5 y)"),("tdAUC_1y","AUC\n1 y"),("tdAUC_2y","AUC\n2 y"),("tdAUC_3y","AUC\n3 y")]
x=np.arange(len(mets)); ax.bar(x-0.19,[b[m] for m,_ in mets],0.36,color="#D1D5DB",label="Molecular group + age + sex"); ax.bar(x+0.19,[n[m] for m,_ in mets],0.36,color="#1D4ED8",label="+ immune state")
ax.set_xticks(x); ax.set_xticklabels([l for _,l in mets],fontsize=5.8); ax.set_ylim(0.5,0.9); ax.set_ylabel("Discrimination"); ax.legend(frameon=False,fontsize=5.8,loc="upper left")
ax.text(0.02,0.74,f"LRT χ² = {n.LRT_vs_base_chi2:.1f} (2 df), P = {n.LRT_p:.3f}\nΔAIC = {n.AIC-b.AIC:.1f}\nΔC = +{n.delta_C_corrected:.3f} (95% CI {n.delta_C_boot_lo:.3f}–{n.delta_C_boot_hi:.3f})",transform=ax.transAxes,fontsize=5.8,va="top")
ax.set_title("Incremental prognostic value of immune state",fontsize=7.4); lab(ax,"E")
os.makedirs(ROOT+"fig/FigureS_survival",exist_ok=True)
for ext in ("svg","pdf","png"): fig.savefig(ROOT+f"fig/FigureS_survival/FigureS_survival.{ext}",dpi=400 if ext=="png" else None,bbox_inches="tight",facecolor="white")
print(Q[["analysis","model","LRT_p","Harrell_C_optimism_corrected","delta_C_corrected","delta_C_boot_lo","delta_C_boot_hi"]].round(4))
