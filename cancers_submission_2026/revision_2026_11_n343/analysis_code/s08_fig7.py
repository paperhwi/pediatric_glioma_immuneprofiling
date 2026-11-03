"""Figure 7 (main): A Kaplan-Meier by immune state with a 10-year axis and numbers at risk; B forest plot of the
pre-specified multivariable Cox model. Corrected OS, n = 245 (206 deaths)."""
import sys, os; sys.path.insert(0,"/home/claude/q343/src"); from common import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import NullLocator
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import multivariate_logrank_test
plt.rcParams.update({"font.size":7.2,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.7,"font.family":"DejaVu Sans"})
a=load(); s=a[a.os_eval]
def lab(ax,t,x=-0.13): ax.text(x,1.06,t,transform=ax.transAxes,fontsize=10,fontweight="bold",va="bottom")
def km(ax,df,xmax,ticks,title,fs=6):
    ks=[]
    for e in ORD:
        m=df.state==e; k=KaplanMeierFitter().fit(df.time[m]/365.25,df.event[m],label=f"{SHORT[e]} (n = {m.sum()}, deaths = {int(df.event[m].sum())})")
        k.plot_survival_function(ax=ax,ci_show=False,color=COL[e],lw=1.4,show_censors=True,censor_styles={"ms":3.5,"mew":0.8}); ks.append(k)
    p=multivariate_logrank_test(df.time,df.state,df.event).p_value
    ax.set_xlim(0,xmax); ax.set_xticks(ticks); ax.set_ylim(0,1.03); ax.set_xlabel("Overall survival (years)"); ax.set_ylabel("Survival probability")
    ax.set_title(f"{title}\nlog-rank P = {p:.3f}",fontsize=7.6); ax.legend(frameon=False,fontsize=fs,loc="upper right",title="Immune state",title_fontsize=fs)
    add_at_risk_counts(*ks,ax=ax,rows_to_show=["At risk"],fontsize=5.8,labels=[SHORT[e] for e in ORD],xticks=ticks,ypos=-0.5)
    return p
def fmtp(p): 
    if p>=0.001: return f"{p:.2g}" if p<0.1 else f"{p:.2f}"
    m,e=f"{p:.1e}".split("e"); return f"{m} × 10$^{{{int(e)}}}$"
def forest(ax,rows,r,xlim=(0.05,90),bound=8,txt=9,ticks=(0.1,0.25,0.5,1,2,4)):
    y=np.arange(len(rows))[::-1]
    for yi,(k,l) in zip(y,rows):
        h,lo,hi,p=r.loc[k,["exp(coef)","exp(coef) lower 95%","exp(coef) upper 95%","p"]]
        col=("#B91C1C" if h>1 else "#1D4ED8") if p<0.05 else "#6B7280"
        ax.plot([lo,hi],[yi,yi],color=col,lw=1.2); ax.plot(h,yi,"s",color=col,ms=4.5)
        ax.text(txt,yi,f"{h:.2f} ({lo:.2f}–{hi:.2f})   P = {fmtp(p)}",va="center",fontsize=6)
    ax.axvline(1,ls="--",lw=0.6,color="k"); ax.set_xscale("log"); ax.set_xlim(*xlim); ax.spines["bottom"].set_bounds(xlim[0],bound)
    ax.set_yticks(y); ax.set_yticklabels([l for _,l in rows],fontsize=6.6); ax.set_xticks(ticks); ax.set_xticklabels([f"{t:g}" for t in ticks]); ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(axis="y",length=0)
if __name__=="__main__":
    fig=plt.figure(figsize=(7.4,3.7)); gs=GridSpec(1,2,width_ratios=[1,1.15],wspace=0.75,left=0.09,right=0.97,top=0.86,bottom=0.30)
    ax=fig.add_subplot(gs[0,0]); km(ax,s,10,[0,2,4,6,8,10],"All tumours with overall survival data"); lab(ax,"A")
    ax=fig.add_subplot(gs[0,1]); d=design(s); c=cox(d); r=c.summary
    rows=[(HI,"Immune-high vs Immune-low"),(MID,"Immune-intermediate vs Immune-low"),("DMG_K27","DMG vs pHGG-WT"),("DHG_G34","DHG vs pHGG-WT"),("IHG","IHG vs pHGG-WT"),("age_years","Age at diagnosis (per year)"),("male","Male vs female")]
    forest(ax,rows,r)
    ax.set_xlabel("Adjusted hazard ratio (95% CI)"); ax.set_title(f"Multivariable Cox model (pre-specified)\nn = {len(d)}, {int(d.event.sum())} deaths",fontsize=7.6); lab(ax,"B",-0.6)
    os.makedirs(ROOT+"fig/Figure7",exist_ok=True)
    for ext in ("svg","pdf","png"): fig.savefig(ROOT+f"fig/Figure7/Figure7.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight",facecolor="white")
    # source data
    with pd.ExcelWriter(ROOT+"fig/Figure7/Figure7_source_data.xlsx") as w:
        s[["Kids_First_Participant_ID","state","cohort_group","age_first_dx_years","reported_gender","time","event","os_corrected","autopsy"]].rename(columns={"time":"OS_days_corrected","event":"death"}).to_excel(w,sheet_name="A_per_tumour")
        r.rename(index={HI:"Immune-high vs Immune-low",MID:"Immune-intermediate vs Immune-low"}).to_excel(w,sheet_name="B_cox")
    print(r[["exp(coef)","exp(coef) lower 95%","exp(coef) upper 95%","p"]].round(4))
