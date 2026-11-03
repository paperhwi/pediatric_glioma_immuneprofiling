from pathlib import Path
import json,shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

work=Path(__file__).parent
from analysis_paths import OUT as out
C=json.loads((out/'tables/clustering_result.json').read_text(encoding='utf-8'))
S=json.loads((out/'tables/survival_result.json').read_text(encoding='utf-8'))
old=out/'inputs/reference'
qa=out/'.qa';qa.mkdir(exist_ok=True)

# Reproduction of all 25 delivered seed summaries, in addition to modal labels.
seed=pd.read_csv(out/'tables/seed_comparison.tsv',sep='\t').set_index('seed')
original_seed=pd.read_csv(old/'T_seed_partitions_n343.tsv',sep='\t').set_index('seed')
ret=pd.read_csv(out/'tables/assignments_25seeds_retained.tsv',sep='\t',index_col=0)
for k in range(42,67):
    for state,short in [('Immune-high','High'),('Immune-intermediate','Intermediate'),('Immune-low','Low')]:
        assert (ret[str(k)]==state).sum()==original_seed.loc[k,f'n_{short}']
    assert abs(seed.loc[k,'retained_ARI_vs_delivered']-original_seed.loc[k,'ARI_vs_modal'])<1e-12

scale=pd.read_csv(out/'tables/feature_standardization_audit.tsv',sep='\t')
scale['relative_feature_weight_multiplier']=1/scale.retained_SD_n343
scale.to_csv(out/'tables/feature_standardization_audit.tsv',sep='\t',index=False)
knew=pd.read_csv(out/'tables/ksweep_restandardized_n343.tsv',sep='\t');kold=pd.read_csv(old/'ksweep_n343.tsv',sep='\t')
gnew=pd.read_csv(out/'tables/gap_restandardized_n343.tsv',sep='\t');gold=pd.read_csv(old/'gap_n343.tsv',sep='\t')
kd=kold[['k','PAC','silhouette_consensus']].rename(columns={'PAC':'retained_PAC','silhouette_consensus':'retained_silhouette'}).merge(knew[['k','PAC','silhouette_consensus']].rename(columns={'PAC':'new_PAC','silhouette_consensus':'new_silhouette'}),on='k',how='outer')
for frame,pre in [(gold,'retained'),(gnew,'new')]:
    kd=kd.merge(frame[['k','gap','s_k']].rename(columns={'gap':pre+'_gap','s_k':pre+'_gap_SE'}),on='k',how='outer')
kd=kd.sort_values('k');kd.to_csv(out/'tables/k_gap_comparison.tsv',sep='\t',index=False)

spec={'reference_language':['Astra governing figures and tables','OAI-VIS-004 Professional table','OAI-VIS-013/014 confidence intervals','OAI-VIS-010 heatmap'],
      'figures':[{'file':'Axis_survival_library_adjustment','claim':'The inverse point estimate is attenuated and its CI crosses 1 after RNA library control','form':'dot and Wald 95% CI forest','metric':'HR per 1 cohort-wide SD, not within-subset SD','size_mm':[178,112],'order':['base','library covariate','library strata'],'groups':['n245/206','no autopsy n189/150'],'axis':'linear HR with null line 1; no bars','formats':['svg','pdf','png']},
                 {'file':'Clustering_standardization_consistency','claim':'334/343 assignments persist, with 9 changes after within-n343 feature scaling','form':'annotated 3 by 3 transition matrix','metric':'number of participants, row old state and column new state','size_mm':[85,78],'uncertainty':'exact cohort counts; descriptive ARI, no inferential CI','formats':['svg','pdf','png']}]}
(qa/'visual_spec.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')

plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':10,'xtick.labelsize':8,'ytick.labelsize':8.5,'text.color':'#171717','axes.labelcolor':'#171717','axes.edgecolor':'#5E6462','axes.linewidth':.6,'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
teal='#1F6F5F';gray='#777773'

fig=plt.figure(figsize=(178/25.4,112/25.4));ax=fig.add_axes([.205,.26,.565,.55])
ys=[6,5,4,2,1,0];rs=S['models']
for i,(r,y) in enumerate(zip(rs,ys)):
    color=gray if i%3==0 else teal;marker=['o','s','D'][i%3]
    ax.errorbar(r['HR_per_cohort_SD'],y,xerr=[[r['HR_per_cohort_SD']-r['CI95_lower']],[r['CI95_upper']-r['HR_per_cohort_SD']]],fmt=marker,markersize=4.5,color=color,elinewidth=.9,capsize=2,markeredgewidth=.7)
ax.axvline(1,color='#8D8D88',lw=.8,ls=(0,(3,3)),zorder=0)
ax.set_ylim(-.7,6.9);ax.set_xlim(.65,1.20);ax.set_xticks([.7,.8,.9,1.,1.1,1.2]);ax.set_yticks(ys,['Base','+ RNA library','Library strata']*2)
ax.tick_params(axis='y',length=0,pad=8)
ax.set_xlabel('Hazard ratio per 1 cohort SD (95% CI)',labelpad=7)
for side in ['top','right','left']:ax.spines[side].set_visible(False)
fig.text(.04,.945,'Library adjustment attenuates the axis survival association',fontsize=10.5,weight='bold',va='top')
for y,label in [(6.8,'All evaluable patients  n = 245, 206 deaths'),(2.8,'Autopsy specimens excluded  n = 189, 150 deaths')]:
    yy=fig.transFigure.inverted().transform(ax.transData.transform((.65,y)))[1]
    fig.text(.04,yy,label,fontsize=8.5,weight='bold')
fig.text(.79,.825,'HR (95% CI)    P',fontsize=8.0,weight='bold')
for r,y in zip(rs,ys):
    yy=fig.transFigure.inverted().transform(ax.transData.transform((.65,y)))[1]
    fig.text(.79,yy,f"{r['HR_per_cohort_SD']:.2f} ({r['CI95_lower']:.2f}–{r['CI95_upper']:.2f})\nP = {r['Wald_P']:.3f}",fontsize=8,va='center',linespacing=1.05)
fig.text(.04,.055,'All models adjust for molecular group, age and sex. Points show unpenalized Cox estimates.\nIntervals are Wald 95% CIs; P values are two-sided. The SD is fixed across all 343 participants.',fontsize=8,linespacing=1.3)
for ext in ['svg','pdf','png']:fig.savefig(out/'figures'/('FigureS6_axis_library_survival.'+ext),dpi=300)
plt.close(fig)

ct=pd.read_csv(out/'tables/state_transition_counts.tsv',sep='\t',index_col=0)
fig=plt.figure(figsize=(85/25.4,78/25.4));ax=fig.add_axes([.31,.31,.64,.45])
cmap=LinearSegmentedColormap.from_list('white_teal',['#FFFFFF','#DDE9E5',teal])
ax.imshow(ct.values,cmap=cmap,vmin=0,vmax=ct.values.max(),aspect='auto')
for i in range(3):
    for j in range(3):ax.text(j,i,str(ct.iloc[i,j]),ha='center',va='center',fontsize=10,color='white' if ct.iloc[i,j]>100 else '#171717')
ax.set_xticks(range(3),['High','Intermediate','Low']);ax.set_yticks(range(3),['High','Intermediate','Low']);ax.tick_params(length=0,labelsize=8)
ax.set_xlabel('Within-343 scaled labels',labelpad=5,fontsize=8);ax.set_ylabel('Retained-input labels',labelpad=5,fontsize=8)
for s in ax.spines.values():s.set_visible(False)
fig.text(.04,.95,'Within-cohort scaling preserves\nmost assignments',fontsize=10,weight='bold',va='top')
fig.text(.04,.055,'334/343 unchanged (97.4%); ARI = 0.915\nCounts use modal labels over seeds 42–66.',fontsize=8,linespacing=1.3)
for ext in ['svg','pdf','png']:fig.savefig(out/'figures'/('FigureS5_standardization_consistency.'+ext),dpi=300)
plt.close(fig)
for file in (out/'figures').glob('*.png'):
    Image.open(file).convert('L').save(qa/(file.stem+'_grayscale.png'))

def tab(name,frame,note,widths=None):
    return {'name':name,'headers':[str(x) for x in frame.columns],'rows':json.loads(frame.replace([np.inf,-np.inf],np.nan).to_json(orient='values')),'note':note,'widths':widths}
cl=pd.DataFrame(C['clustering_summary']).rename(columns={'n_Immune-high':'High n','n_Immune-intermediate':'Intermediate n','n_Immune-low':'Low n'})
cl=cl[['input','n','High n','Intermediate n','Low n','median_pairwise_ARI','min_pairwise_ARI','max_pairwise_ARI','same_label_all25','agreement_vs_delivered','ARI_vs_delivered']]
cx=pd.DataFrame(S['models'])[['subset','library_handling','n','deaths','HR_per_cohort_SD','CI95_lower','CI95_upper','Wald_P','axis_LRT_P','axis_PH_rank_P']]
tabs=[tab('Clustering',cl,'Locked k=3 assignment protocol: 25 seeds, 500 subsamples/seed, 80% sampling, KMeans n_init=10. New scaling uses mean=0, sample SD=1 across 343. Old labels reproduce exactly.',[22,12,14,16,14,20,20,20,20,22,20]),
      tab('Axis survival',cx,'Corrected OS; molecular group + age + sex in every model. Four RNA-library categories; reference poly-A stranded. HR per SD of the fixed n343 axis. Wald 95% CI and two-sided P.',[27,20,10,10,19,18,18,15,16,18]),
      tab('Changed samples',pd.read_csv(out/'tables/changed_assignments.tsv',sep='\t'),'Nine changes after within-n343 scaling. New labels are sensitivity-analysis outputs and have not replaced the submitted main analysis.',[24,25,27,24,20,24]),
      tab('Per sample',pd.read_csv(out/'tables/per_sample_clustering_comparison.tsv',sep='\t'),'One row per participant. Retained-input modal labels exactly match delivered n343 labels; all 343 records retained.',[24,25,25,22,27,24,20,20,24]),
      tab('Seeds',seed.reset_index(),'Paired old/new seed checks. Descriptive stability measures, not independent patient-level inference.',[12,18,22,24,24,18,20,18]),
      tab('Feature scaling',scale,'Affine rescaling of the delivered float64 feature matrix. ddof=1. Original gene/signature scores are not recomputed. Weight multiplier = 1 / retained SD within n343.',[49,22,22,22,22,28]),
      tab('K diagnostics',kd,'Diagnostic protocol differs from assignment: B=1000, n_init=1, seed42. Gap: 50 rotated-uniform references, n_init=10, first 1SE rule. Best PAC/silhouette k2; gap selects k3 in both inputs.',[10,19,21,19,21,19,20,19,20]),
      tab('Library profile',pd.read_csv(out/'tables/library_descriptives.tsv',sep='\t'),'Complete RNA-library annotation in all 245 survival records. Axis units use the same all-cohort SD; autopsy exclusion is an additional sensitivity.',[27,24,12,12,12,19,19,19,19,19])]
payload={'out':str(out),'work':str(work),'tabs':tabs,'transition':{'headers':['Retained state']+list(ct.columns),'rows':[[i]+list(row) for i,row in ct.iterrows()]}}
(qa/'workbook_payload.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps({'seed_summaries_reproduced':25,'plots':2,'workbook_tabs':len(tabs),'scaling_SD_range':[scale.retained_SD_n343.min(),scale.retained_SD_n343.max()]},indent=2))
