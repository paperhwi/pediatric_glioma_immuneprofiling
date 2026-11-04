from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import NullLocator
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.plotting import add_at_risk_counts
from followup_inputs import FOLLOWUP as W, REVISION, observations
O=W/'tables';D=W/'generated_figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.7})

def display_p(v):
 if v >= .001: return f'{v:.3f}'
 a,e=f'{v:.1e}'.split('e'); return a+' × 10'+str(int(e)).translate(str.maketrans('0123456789-','⁰¹²³⁴⁵⁶⁷⁸⁹⁻'))

def savefig(fig,dest):
 dest.parent.mkdir(parents=True,exist_ok=True)
 for ext in ['svg','pdf','png']:fig.savefig(dest.with_suffix('.'+ext),bbox_inches='tight',facecolor='white',dpi=400 if ext=='png' else None)
 plt.close(fig)
# Full-follow-up fits; only x-axis display is restricted.
a=observations();s=a[a.in_survival_analysis].copy();states=['Immune-high','Immune-intermediate','Immune-low'];colors=dict(zip(states,['#3B82F6','#EF4444','#6B7280']))
kmchecks=[]
def km(ax,q,title,fontsize=6.7):
 fits=[]
 for state in states:
  m=q.immune_state==state;k=KaplanMeierFitter().fit(q.loc[m,'OS_days_corrected']/365.25,q.loc[m,'OS_event'],label=f'{state} (n={m.sum()}, deaths={int(q.loc[m,"OS_event"].sum())})')
  k.plot_survival_function(ax=ax,ci_show=False,color=colors[state],lw=1.3,show_censors=True,censor_styles={'ms':3,'mew':.7});fits.append(k)
  kmchecks.append({'plot':title,'state':state,'n':int(m.sum()),'deaths':int(q.loc[m,'OS_event'].sum()),'fit_max_years':float(k.timeline.max()),'display_max_years':5,'at_risk':{str(t):int((q.loc[m,'OS_days_corrected']/365.25>=t).sum()) for t in range(6)}})
 p=multivariate_logrank_test(q.OS_days_corrected,q.immune_state,q.OS_event).p_value
 ax.set(xlim=(0,5),ylim=(0,1.03),xticks=list(range(6)),xlabel='Overall survival (years)',ylabel='Survival probability');ax.set_title(f'{title}\nFull-follow-up log-rank P = {p:.3f}',fontsize=8);ax.legend(frameon=False,fontsize=fontsize,loc='upper right')
 add_at_risk_counts(*fits,ax=ax,rows_to_show=['At risk'],labels=states,xticks=list(range(6)),ypos=-.52,fontsize=6.2)
def letter(ax,t,x=-.15):ax.text(x,1.08,t,transform=ax.transAxes,fontsize=12,fontweight='bold')
terms=pd.read_csv(O/'categorical_cox_all_terms.tsv',sep='\t');primary=terms[terms.model=='Primary original partition'].set_index('term')
fig=plt.figure(figsize=(10.8,4.2));gs=GridSpec(1,2,figure=fig,width_ratios=[1,1.08],left=.18,right=.98,top=.83,bottom=.32,wspace=1.10)
ax=fig.add_subplot(gs[0]);km(ax,s,'All tumors with overall-survival data');letter(ax,'A')
ax=fig.add_subplot(gs[1]);labels={'Immune-high':'Immune-high vs Immune-low','Immune-intermediate':'Immune-intermediate vs Immune-low','DMG_K27':'DMG vs pHGG-WT','DHG_G34':'DHG vs pHGG-WT','IHG':'IHG vs pHGG-WT','age_years':'Age at diagnosis (per year)','male':'Male vs female'}
for y,(term,lbl) in zip(range(6,-1,-1),labels.items()):
 r=primary.loc[term];h,lo,hi,p=r[['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']];ax.plot([lo,hi],[y,y],color='#374151');ax.plot(h,y,'s',color='#2563EB',ms=4);ax.text(10,y,f'{h:.2f} ({lo:.2f}–{hi:.2f}); P={display_p(p)}',fontsize=6.5,va='center')
ax.set_xscale('log');ax.set_xlim(.05,160);ax.set_xticks([.1,.5,1,2,4]);ax.set_xticklabels(['0.1','0.5','1','2','4']);ax.xaxis.set_minor_locator(NullLocator());ax.spines['bottom'].set_bounds(.05,4);ax.axvline(1,color='k',lw=.6,ls='--');ax.set_yticks(range(6,-1,-1),list(labels.values()),fontsize=6.7);ax.tick_params(axis='y',length=0);ax.set_xlabel('Adjusted hazard ratio (95% CI)');ax.set_title('Primary Cox model\nn = 245; 206 deaths',fontsize=8);letter(ax,'B',-.80)
savefig(fig,D/'02_figures/main/Figure7')
def sheet(name):
 return pd.read_csv(O/(name+'.tsv'),sep='\t')
S=pd.read_csv(REVISION/'reproducibility_data/T_survival_sensitivity.tsv',sep='\t')
new=pd.read_csv(O/'categorical_state_sensitivity_cox.tsv',sep='\t');extra=[]
for m,label in [('Restandardized partition primary model','Restandardized partition (116/149/78)'),('Autopsy excluded + RNA library','Autopsy excluded + RNA-library adjustment'),('Autopsy excluded library strata','Autopsy excluded + RNA-library strata')]:
 r=new[(new.model==m)&new.contrast.str.startswith('Immune-high')].iloc[0];mid=new[(new.model==m)&new.contrast.str.startswith('Immune-intermediate')].iloc[0];extra.append({'analysis':label,'n':r.n,'deaths':r.deaths,'high_HR':r.HR,'high_lo':r.CI95_lower,'high_hi':r.CI95_upper,'high_p':r.Wald_P,'intermediate_HR':mid.HR,'intermediate_lo':mid.CI95_lower,'intermediate_hi':mid.CI95_upper,'intermediate_p':mid.Wald_P})
S=S[~S.analysis.isin([r['analysis'] for r in extra])].copy();S=pd.concat([S,pd.DataFrame(extra)],ignore_index=True);S.to_csv(O/'FigureS1_sensitivity_source.tsv',sep='\t',index=False)
fig=plt.figure(figsize=(10,15));gs=GridSpec(3,2,figure=fig,height_ratios=[1.70,1,.96],left=.19,right=.97,top=.96,bottom=.06,hspace=.95,wspace=.65)
ax=fig.add_subplot(gs[0,:]);pos=ax.get_position();ax.set_position([.50,pos.y0,.21,pos.height]);yy=np.arange(len(S))[::-1]
shortlabels=['Primary model','Autopsy excluded','Time origin at specimen collection','Autopsy-status adjustment','RNA-library adjustment','RNA-library strata','Delayed entry at specimen collection','Released survival times','Prior labels (pre-exclusion cohort)','Administrative censoring at 5 years','Administrative censoring at 10 years','Specimen within 1 year of diagnosis','Exclude 13 earlier-origin records','Unadjusted Cox']+[r['analysis'] for r in extra]
assert len(shortlabels)==len(S)
for y,(_,r) in zip(yy,S.iterrows()):
 color='#2563EB' if y==yy[0] else '#374151';ax.plot([r.high_lo,r.high_hi],[y,y],color=color,lw=1);ax.plot(r.high_HR,y,'s',color=color,ms=4);ax.text(2.3,y,f'{int(r.n)}/{int(r.deaths)}   {r.high_HR:.2f} ({r.high_lo:.2f}–{r.high_hi:.2f})   P={display_p(r.high_p)}',fontsize=6.7,va='center')
ax.set_xscale('log');ax.set_xlim(.25,2);ax.set_xticks([.25,.5,1,2]);ax.set_xticklabels(['0.25','0.5','1','2']);ax.xaxis.set_minor_locator(NullLocator());ax.axvline(1,color='k',ls='--',lw=.6);ax.set_yticks(yy,shortlabels,fontsize=7);ax.tick_params(axis='y',length=0);ax.set_xlabel('Immune-high vs Immune-low HR (95% CI)',fontsize=7.5);ax.set_title('Survival sensitivity analyses',fontsize=9,pad=18);ax.text(2.3,len(S)-.15,'n/deaths   HR (95% CI); P',fontsize=7,fontweight='bold');letter(ax,'A',-1.48)
for j,(g,title) in enumerate([('DMG_K27','DMG, H3 K27-altered'),('pHGG_WT','pHGG, H3/IDH-wildtype')]):
 ax=fig.add_subplot(gs[1,j]);km(ax,s[s.molecular_group==g],title,fontsize=5.9);letter(ax,'BC'[j])
ax=fig.add_subplot(gs[2,0]);T=sheet('Within_subtype_Cox');entries=[('pHGG_WT','Immune-high'),('pHGG_WT','Immune-intermediate'),('DMG_K27','Immune-high'),('DMG_K27','Immune-intermediate')]
for y,(g,st) in zip(range(3,-1,-1),entries):
 r=T[(T.group==g)&(T.term==st+' vs Immune-low')].iloc[0];ax.plot([r.lo,r.hi],[y,y],color=colors[st]);ax.plot(r.HR,y,'s',color=colors[st],ms=4);ax.text(2.5,y,f'{r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f})\nP={display_p(r.p)}',fontsize=6,va='center')
ax.set_xscale('log');ax.set_xlim(.15,8);ax.set_xticks([.25,.5,1,2]);ax.set_xticklabels(['0.25','0.5','1','2']);ax.xaxis.set_minor_locator(NullLocator());ax.axvline(1,color='k',ls='--',lw=.6);ax.spines['bottom'].set_bounds(.15,2);ax.set_yticks(range(3,-1,-1),[('pHGG-WT' if g=='pHGG_WT' else 'DMG')+': '+st+'\nvs Immune-low' for g,st in entries],fontsize=6.1);ax.set_xlabel('Age- and sex-adjusted HR (95% CI)');ax.set_title('Within-subtype Cox models\nState × subtype interaction P = 0.78',fontsize=8.5);letter(ax,'D',-.75)
ax=fig.add_subplot(gs[2,1]);Q=sheet('Incremental_value');Q=Q[Q.analysis.str.startswith('Primary')];b,n=Q.iloc[0],Q.iloc[1];mets=[('Harrell_C_optimism_corrected','C-index\noptimism\ncorrected'),('Uno_C_tau5y','Uno C\n5 years'),('tdAUC_1y','AUC\n1 year'),('tdAUC_2y','AUC\n2 years'),('tdAUC_3y','AUC\n3 years')];x=np.arange(5)
ax.bar(x-.19,[b[m] for m,_ in mets],.36,color='#D1D5DB',label='Molecular group + age + sex');ax.bar(x+.19,[n[m] for m,_ in mets],.36,color='#2563EB',label='+ immune state');ax.set_xticks(x,[l for _,l in mets],fontsize=6.2);ax.set_ylim(.5,.9);ax.set_ylabel('Discrimination');ax.set_title('Incremental prognostic information',fontsize=8.5);ax.legend(frameon=False,fontsize=6,loc='upper left');ax.text(.02,.72,f'LRT P={n.LRT_p:.3f}\nΔAIC={n.AIC-b.AIC:.1f}\nΔC={n.delta_C_corrected:.3f} (95% CI {n.delta_C_boot_lo:.3f}–{n.delta_C_boot_hi:.3f})',transform=ax.transAxes,fontsize=6.1,va='top');letter(ax,'E')
savefig(fig,D/'02_figures/supplementary/FigureS1_survival')
# IPW figure restores the submitted missing-outcome sensitivity with transparent diagnostics.
I=pd.read_csv(O/'ipw_cox.tsv',sep='\t');B=pd.read_csv(O/'ipw_balance.tsv',sep='\t');B=B[B.comparison=='Full eligible cohort target'].sort_values('SMD_unweighted');M=pd.read_csv(O/'ipw_missingness.tsv',sep='\t')
fig,axes=plt.subplots(1,3,figsize=(13,4.4),gridspec_kw={'width_ratios':[.75,1,1.3]});fig.subplots_adjust(left=.07,right=.98,bottom=.26,top=.88,wspace=1.1)
ax=axes[0];ax.bar(M.molecular_group,M.missing_percent,color='#64748B');ax.set_ylabel('Missing overall-survival data (%)');ax.set_ylim(0,100);ax.tick_params(axis='x',rotation=35,labelsize=7);ax.set_title('Outcome availability',fontsize=9);letter(ax,'A')
ax=axes[1];yy=np.arange(len(B));ax.scatter(B.SMD_unweighted,yy,s=15,color='#6B7280',label='Unweighted');ax.scatter(B.SMD_IPW_1_99,yy,s=15,color='#2563EB',label='IPW 1–99% capped');ax.set_yticks(yy,B.covariate,fontsize=6);ax.axvline(0,color='k',lw=.5);ax.axvline(.1,color='gray',ls='--',lw=.6);ax.axvline(-.1,color='gray',ls='--',lw=.6);ax.set_xlabel('SMD vs all 343 eligible patients');ax.set_title('Covariate balance',fontsize=9);ax.legend(frameon=False,fontsize=6);letter(ax,'B',-.3)
ax=axes[2];q=I[I.contrast.str.startswith('Immune-high')];yy=np.arange(len(q))[::-1]
for y,(_,r) in zip(yy,q.iterrows()):ax.plot([r.CI95_lower,r.CI95_upper],[y,y],color='#374151');ax.plot(r.HR,y,'s',color='#2563EB',ms=4);ax.text(1.03,y,f'{r.HR:.2f} ({r.CI95_lower:.2f}–{r.CI95_upper:.2f})\nP={display_p(r.Wald_P)}',fontsize=6,va='center')
ax.set_yticks(yy,['Unweighted','IPW 1–99% capped','IPW uncapped','IPW 2.5–97.5% capped','Group-only propensity'],fontsize=6);ax.axvline(1,color='k',ls='--',lw=.6);ax.set_xlim(.25,1);ax.set_xlabel('Immune-high vs Immune-low HR');ax.set_title('Weighted Cox sensitivity',fontsize=9);letter(ax,'C',-.25)
savefig(fig,D/'02_figures/supplementary/FigureS7_IPW_sensitivity')
(O/'KM_display_fit_audit.json').write_text(json.dumps(kmchecks,indent=2),encoding='utf-8');print('Figures exported to',D)
