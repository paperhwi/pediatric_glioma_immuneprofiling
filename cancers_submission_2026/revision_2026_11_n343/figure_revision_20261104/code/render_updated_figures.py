from pathlib import Path
import json,re,ast,copy,sys
import numpy as np,pandas as pd,h5py
from lxml import etree as E
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import spearmanr
import cairo,cairosvg
import argparse
parser=argparse.ArgumentParser(description="Regenerate the November 4 n343 figure corrections. Restricted/reference inputs remain local.")
parser.add_argument('--input-config',type=Path,required=True)
parser.add_argument('--baseline-figures',type=Path,required=True,help='figures directory from the immutable v2.0 snapshot')
parser.add_argument('--output-dir',type=Path,required=True)
parser.add_argument('--only',nargs='*',default=[])
args=parser.parse_args()
cfg=json.loads(args.input_config.read_text(encoding='utf-8'))
W=Path(__file__).resolve().parents[1];R=W.parent;D=args.output_dir;D.mkdir(parents=True,exist_ok=True)
S=args.baseline_figures;O=R/'reviewer_followup_20261103/tables';F=D/'figures';T=R/'reproducibility_data';generated=D/'generated';generated.mkdir(exist_ok=True)
BASELINE_NAMES={'FigureS_survival':'FigureS1_survival','FigureS_state_structure':'FigureS2_state_structure','FigureS_immune_features':'FigureS3_immune_features'}
COL={'Immune-high':'#3B82F6','Immune-intermediate':'#EF4444','Immune-low':'#9CA3AF'};ST=list(COL)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':8.5,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.7})
SPEC={'final_width_mm':178,'panel_letters':'14 pt bold at final width; separate from titles','P_display':'P >= 0.001: three decimal places; smaller P: one mantissa decimal and superscript power of ten. P thresholds retain explicit cutoffs.','color':'Preserve current immune-state blue/red/gray; preserve existing evidence encodings.','Figure2':'Same n343 QC, correlations and within-n343 z-score boxes; horizontal B bars and 2 by 2 C layout prevent label collisions.','Figure5':'Same B-E data; corrected complementary annotated T-cell target for curated signature in A; no bar interior values; no unannotated NK performance inferred.','format':'Editable SVG and vector PDF, 300 dpi PNG; source numeric precision retained separately.'}
SPEC['format']='Editable SVG and vector PDF, 400 dpi PNG at native PDF size; source numeric precision retained separately.'
(D/'visual_spec.json').write_text(json.dumps(SPEC,indent=2),encoding='utf-8')
def letter(ax,t):ax.text(-.12,1.08,t,transform=ax.transAxes,size=14,weight='bold',ha='left',va='bottom',gid='panel_letter_'+t)
def save(fig,name):fig.savefig(generated/(name+'.svg'),bbox_inches='tight',facecolor='white');plt.close(fig)
def nice_method(m):return {'MCP':'MCP-counter','MCP_DC':'MCP-counter\nmyeloid DC'}.get(m,m)
def readobs(p):
 with h5py.File(p) as f:
  g=f['obs'];d={}
  for k in g:
   a=g[k]
   if isinstance(a,h5py.Group):
    cat=a['categories'][()];cat=np.array([x.decode() if isinstance(x,bytes) else x for x in cat]);codes=a['codes'][()];v=cat[np.maximum(codes,0)].astype(object);v[codes<0]=None
   else:v=a[()];v=np.array([x.decode() if isinstance(x,bytes) else x for x in v])
   d[k]=v
 return pd.DataFrame(d).set_index('_index')
# Figure 2: the same plot data, with a readable layout.
obs=pd.read_csv(R/'reproducibility_data/immune_state_assignments_n343.tsv',sep='\t',index_col=0).rename(columns={'state':'immune_state'})
qc=pd.read_csv(cfg['cibersort_qc']);qc=qc[qc.Mixture.isin(obs.index)];assert len(qc)==343
long=pd.read_csv(cfg['cross_method_long'],sep='\t');long=long[long['sample'].isin(obs.index)].copy();long['state']=long['sample'].map(obs.immune_state);assert long.state.notna().all()
long['z343']=long.groupby(['axis','method']).z.transform(lambda s:(s-s.mean())/s.std(ddof=0))
rr=pd.read_csv(T/'T_fig2B_cross_method_n343.tsv',sep='\t');peer={'CD8 T cell':['xCell','quanTIseq','MCP','EPIC'],'NK cell':['xCell','quanTIseq','MCP'],'Macrophage':['xCell','MCP','EPIC'],'Monocyte / DC':['xCell','MCP_DC']};mc={'xCell':'#2563EB','quanTIseq':'#059669','MCP':'#D97706','EPIC':'#DC2626','MCP_DC':'#7C3AED'}
fig=plt.figure(figsize=(7,10.5));gs=fig.add_gridspec(3,1,height_ratios=[1.15,3.2,3.5],left=.13,right=.98,top=.96,bottom=.09,hspace=.45)
top=gs[0].subgridspec(1,3,wspace=.40)
for j,field in enumerate(['P-value','Correlation','RMSE']):
 ax=fig.add_subplot(top[j]);ax.hist(qc[field],bins=30,color='#8B5CF6',edgecolor='#453078',lw=.3);ax.set(title=field,xlabel=field)
 if j==0:ax.axvline(.05,color='#EF4444',ls='--',lw=.8);ax.set_ylabel('Specimens');letter(ax,'A')
ax=fig.add_subplot(gs[1]);pos=ax.get_position();ax.set_position([.37,pos.y0,.60,pos.height]);y=np.arange(len(rr))[::-1]
ax.barh(y,rr.spearman_rho,color=[mc[m] for m in rr.peer_method],height=.66)
for yi,(_,r) in zip(y,rr.iterrows()):
 mark='***' if r.q_BH<.001 else '**' if r.q_BH<.01 else '*' if r.q_BH<.05 else 'ns';v=r.spearman_rho;ax.text(v+(.012 if v>=0 else -.012),yi,mark,va='center',ha='left' if v>=0 else 'right',size=7)
labels=[r.axis.replace(' / ','/')+' — '+nice_method(r.peer_method).replace('\n',' ') for _,r in rr.iterrows()]
ax.set_yticks(y,labels,fontsize=7.1);ax.set(xlim=(-.15,.48),xlabel='Spearman ρ versus LM22');ax.set_title('Cross-method concordance (n = 343)',pad=10);ax.axvline(0,color='#555',lw=.6);ax.tick_params(axis='y',length=0);letter(ax,'B')
bot=gs[2].subgridspec(2,2,wspace=.40,hspace=.65)
for j,(axis,methods0) in enumerate(peer.items()):
 ax=fig.add_subplot(bot[j//2,j%2]);methods=['LM22',*methods0];sub=long[long.axis==axis]
 for k,m in enumerate(methods):
  for state,offset in zip(ST,[-.25,0,.25]):
   vals=sub.loc[(sub.method==m)&(sub.state==state),'z343'].to_numpy();assert len(vals)>0
   bp=ax.boxplot(vals,positions=[k+offset],widths=.22,patch_artist=True,showfliers=False)
   bp['boxes'][0].set(facecolor=COL[state],edgecolor='#5B6470',lw=.5)
   for q in bp['whiskers']+bp['caps']+bp['medians']:q.set(color='#5B6470',lw=.6)
 ax.axhline(0,color='#888',lw=.5);ax.set_xticks(range(len(methods)),[nice_method(m) for m in methods],rotation=32,ha='right',fontsize=6.7);ax.tick_params(axis='x',pad=7);ax.set_title(axis,pad=9)
 if j%2==0:ax.set_ylabel('Within-cohort z score')
 if j==0:letter(ax,'C')
fig.legend(handles=[Patch(color=COL[s],label=s) for s in ST],loc='lower center',bbox_to_anchor=(.52,.015),ncol=3,frameon=False,fontsize=7)
save(fig,'Figure2')
# Figure 5: preserve B-E exact final-n343 data and make A wholly vector.
auc=pd.read_csv(W/'tables/Figure5A_AUROC_current_targets.tsv',sep='\t').dropna(subset=['AUROC_all_cells']).sort_values('AUROC_all_cells')
designs=pd.read_csv(T/'T_fig5B_designs.tsv',sep='\t');record=designs.loc[designs.gene_set.str.startswith('top 300 Immune-high')&designs.design.str.startswith('D4 all ten')].iloc[0];counts={}
for part in record.top_counts.split('; '):k,v=part.rsplit(' ',1);counts[k]=int(v)
assert sum(counts.values())==300 and int(record.n_immune_top)==294
spot_path=cfg['visium_h5ad']
spots=readobs(spot_path)
with h5py.File(spot_path) as f:
 xy=f['obsm/spatial'][()];spots['x']=xy[:,0];spots['y']=xy[:,1]
 # X contains the seven nonnegative reference-cell shares in var order;
 # every class except malignant and oligodendrocyte is immune.
 names=[s.decode() for s in f['var/_index'][()]]
 immune=[i for i,n in enumerate(names) if n not in ['Malignant','Oligodendrocyte']]
 assert len(immune)==5
 spots['immune_fraction']=f['X'][()][:,immune].sum(axis=1)
expected=pd.read_csv(T/'T_visium_by_section.tsv',sep='\t').set_index('section')
assert np.allclose(spots.groupby('section').apply(lambda z:spearmanr(z.axis_score,z.immune_fraction).statistic,include_groups=False).values,expected.rho_axis_vs_immune.values)
spots['spot_class']=spots.spot_class.replace({'High':'Immune-high','Intermediate':'Immune-intermediate','Low':'Immune-low','Lymphocyte-inflamed':'Immune-high','Myeloid-dominant':'Immune-intermediate','Immune-desert':'Immune-low'})
assert spots.spot_class.isin(ST).all(),spots.spot_class.unique()
sections=['Sample-1','Sample-2','Sample-3'];assert [int((spots.section==s).sum()) for s in sections]==[577,856,217]
fig=plt.figure(figsize=(7,12));gs=fig.add_gridspec(4,1,height_ratios=[3.2,1.9,1.9,2.0],left=.12,right=.98,top=.97,bottom=.06,hspace=.58)
ax=fig.add_subplot(gs[0]);pos=ax.get_position();ax.set_position([.41,pos.y0,.56,pos.height]);yy=np.arange(len(auc));ac=['#6D28D9' if ('T cell' in r.expected or 'T/NK' in r.expected) else '#0F766E' for _,r in auc.iterrows()]
ax.barh(yy,auc.AUROC_all_cells,color=ac,height=.7)
for yi,(_,r) in zip(yy,auc.iterrows()):ax.plot([r.AUROC_min_sample,r.AUROC_max_sample],[yi,yi],color='#475569',lw=.7)
ax.axvline(.5,color='#374151',ls='--',lw=.7);ax.set_yticks(yy,auc.signature.str.replace('_',' '),fontsize=7.1);ax.tick_params(axis='y',length=0);ax.set(xlim=(0,1.02),xlabel='AUROC for annotated target lineage');ax.set_title('Directional signature discrimination',pad=12);letter(ax,'A')
ax.legend(handles=[Patch(color='#0F766E',label='Myeloid target'),Patch(color='#6D28D9',label='T-cell target')],frameon=False,ncol=2,loc='lower center',bbox_to_anchor=(.5,1.15),fontsize=7)
bg=gs[1].subgridspec(1,2,width_ratios=[1.6,1],wspace=.55)
ax=fig.add_subplot(bg[0]);names=list(counts);yy=np.arange(len(names))[::-1];ax.barh(yy,[counts[n] for n in names],color=['#0F766E' if n not in ['Endothelial','Mural','Oligodendrocyte','Malignant'] else '#94A3B8' for n in names],height=.70);ax.set_yticks(yy,[n.replace('Macrophage/monocyte','Macrophage/\nmonocyte') for n in names],fontsize=7);ax.tick_params(axis='y',length=0,pad=5);ax.set_xlabel('Top-expression genes (n)');ax.set_title('Cell class with highest mean expression',fontsize=8,pad=12);letter(ax,'B')
ax=fig.add_subplot(bg[1]);ax.bar([0,1],[record.pct_immune_top,record.random_mean_pct],color=['#0F766E','#94A3B8']);ax.errorbar(1,record.random_mean_pct,yerr=[[record.random_mean_pct-record.random_lo],[record.random_hi-record.random_mean_pct]],fmt='none',color='#334155',capsize=3);ax.set(ylim=(0,115),ylabel='Immune highest (%)');ax.set_xticks([0,1],['Top 300\ngenes','Matched\nrandom']);ax.text(0,101,f'{record.pct_immune_top:.1f}%',ha='center',size=8,color='#0F766E',weight='bold');ax.text(1,record.random_hi+5,f'{record.random_mean_pct:.1f}%',ha='center',size=8)
cg=gs[2].subgridspec(1,3,wspace=.14);cpositions=[]
for j,s in enumerate(sections):
 ax=fig.add_subplot(cg[j]);q=spots[spots.section==s]
 for state in ST:
  z=q[q.spot_class==state];ax.scatter(z.x,z.y,s=3.1,color=COL[state],alpha=.82,lw=0,label=state)
 ax.set_aspect('equal',adjustable='datalim');ax.invert_yaxis();ax.axis('off');ax.set_title(f'{s} (n = {len(q)})',fontsize=8,pad=8)
 cpositions.append(ax.get_position().y0)
 if j==0:letter(ax,'C');handles,labels=ax.get_legend_handles_labels()
fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.54,min(cpositions)-.012),frameon=False,ncol=3,fontsize=7,title='Nearest immune-state centroid',title_fontsize=7)
dg=gs[3].subgridspec(1,2,wspace=.5);ax=fig.add_subplot(dg[0])
for j,s in enumerate(sections):
 q=spots[spots.section==s]
 for k,state in enumerate(ST[::-1]):
  vals=q.loc[q.spot_class==state,'immune_fraction'].to_numpy();bp=ax.boxplot(vals,positions=[j*4+k],widths=.7,patch_artist=True,showfliers=False);bp['boxes'][0].set(facecolor=COL[state],edgecolor='#475569',lw=.6)
  for item in bp['whiskers']+bp['caps']+bp['medians']:item.set(color='#475569',lw=.7)
ax.set_xticks([1,5,9],sections);ax.set_ylabel('Estimated immune fraction');ax.set_title('Spot immune fraction by centroid',size=8,pad=12);letter(ax,'D')
ax=fig.add_subplot(dg[1])
for j,s in enumerate(sections):
 q=spots[spots.section==s];rho=spearmanr(q.axis_score,q.immune_fraction).statistic;mixed=(q.n_types_ge10pct>=2).mean()*100;ax.scatter(q.axis_score,q.immune_fraction,s=3,alpha=.27,color=['#2563EB','#D97706','#0F766E'][j],label=f'{s} (ρ = {rho:.2f}; mixed {mixed:.0f}%)')
ax.set(xlabel='Spot immune-axis score',ylabel='Estimated immune fraction');ax.set_title('Spatial axis and cell mixtures',size=8,pad=12);ax.legend(frameon=False,fontsize=5.9,loc='upper left');letter(ax,'E');save(fig,'Figure5')
# S7: translate all code-facing factor names in the display, retaining exact data.
GROUP={'DMG_K27':'DMG, H3 K27-altered','pHGG_WT':'pHGG, H3/IDH-wildtype','DHG_G34':'DHG, H3 G34-mutant','IHG':'Infant-type hemispheric glioma'}
def covname(s):
 if s.startswith('molecular_group_'):return GROUP.get(s.removeprefix('molecular_group_'),s)
 if s.startswith('location_'):return {'Other/NA':'Other/not annotated','Midline':'Midline','Hemispheric':'Hemispheric','Posterior_fossa':'Posterior fossa','Posterior fossa':'Posterior fossa','Multi-compartment/ambiguous':'Multiple compartments/ambiguous'}.get(s.removeprefix('location_'),s.removeprefix('location_'))
 return {'age_years':'Age at diagnosis','age':'Age at diagnosis','male':'Male sex','Immune-high':'Immune-high','Immune-intermediate':'Immune-intermediate','immune_state_Immune-high':'Immune-high','immune_state_Immune-intermediate':'Immune-intermediate'}.get(s,s.replace('_',' ')).replace('Posterior_fossa','Posterior fossa')
B=pd.read_csv(O/'ipw_balance.tsv',sep='\t');B=B[B.comparison=='Full eligible cohort target'].sort_values('SMD_unweighted');I=pd.read_csv(O/'ipw_cox.tsv',sep='\t');M=pd.read_csv(O/'ipw_missingness.tsv',sep='\t')
fig=plt.figure(figsize=(7,8.8));gs=fig.add_gridspec(2,2,width_ratios=[1,1.35],height_ratios=[1.35,1],left=.13,right=.98,top=.94,bottom=.09,wspace=1.15,hspace=.85)
ax=fig.add_subplot(gs[0,0]);ax.bar(np.arange(len(M)),M.missing_percent,color='#64748B');ax.set_xticks(np.arange(len(M)),[GROUP.get(g,g).replace(', ','\n') for g in M.molecular_group],rotation=55,ha='right',fontsize=6.7);ax.set(ylabel='Missing overall-survival data (%)',ylim=(0,100));ax.set_title('Outcome availability',pad=12);letter(ax,'A')
ax=fig.add_subplot(gs[0,1]);y=np.arange(len(B));ax.scatter(B.SMD_unweighted,y,s=14,color='#6B7280',label='Unweighted');ax.scatter(B.SMD_IPW_1_99,y,s=14,color='#2563EB',label='IPW 1–99% capped');ax.set_yticks(y,[covname(x) for x in B.covariate],fontsize=6.7);ax.tick_params(axis='y',length=0);ax.axvline(0,color='k',lw=.5)
for x in [-.1,.1]:ax.axvline(x,color='gray',ls='--',lw=.6)
ax.set_xlabel('SMD vs all 343 eligible patients');ax.set_title('Covariate balance',pad=12);ax.legend(frameon=False,fontsize=6,loc='lower center',bbox_to_anchor=(.5,-.26),ncol=2);letter(ax,'B')
ax=fig.add_subplot(gs[1,:]);pos=ax.get_position();ax.set_position([.43,pos.y0,.35,pos.height]);q=I[I.contrast.str.startswith('Immune-high')];y=np.arange(len(q))[::-1]
for yi,(_,r) in zip(y,q.iterrows()):ax.plot([r.CI95_lower,r.CI95_upper],[yi,yi],color='#374151');ax.plot(r.HR,yi,'s',color='#2563EB',ms=4);ax.text(1.04,yi,f'{r.HR:.2f} ({r.CI95_lower:.2f}–{r.CI95_upper:.2f})\nP = {r.Wald_P:.3f}',fontsize=7,va='center')
ax.set_yticks(y,['Unweighted','IPW 1–99% capped','IPW uncapped','IPW 2.5–97.5% capped','Molecular-group-only propensity'],fontsize=7);ax.axvline(1,color='k',ls='--',lw=.6);ax.set(xlim=(.25,1),xlabel='Immune-high vs Immune-low HR');ax.set_title('Weighted Cox sensitivity',pad=16);letter(ax,'C');save(fig,'FigureS7_IPW_sensitivity')
# Editable SVG text utilities from the prior verified builder, without running it.
from svg_text_utils import outlined,tag
env={'outlined':outlined,'tag':tag};log=[]
def plain(n,s):
 for c in list(n):n.remove(c)
 n.text=s
def sty(n,**k):
 d=dict(x.strip().split(':',1) for x in n.get('style','').split(';') if ':' in x);d.update({a.replace('_','-'):str(b) for a,b in k.items()});n.set('style','; '.join(f'{a}:{b}' for a,b in d.items()))
SUP=str.maketrans('0123456789-','⁰¹²³⁴⁵⁶⁷⁸⁹⁻')
def pformat(v):
 if v>=.001:return f'{v:.3f}'
 a,e=f'{v:.1e}'.split('e');return a+' × 10'+str(int(e)).translate(SUP)
preg=re.compile(r'(?P<prefix>\b[Pp](?:_MC|\(perm\))?\s*=\s*)(?P<val>\d*\.?\d+(?:[eE][+−-]?\d+)?)')
def ptext(s):return preg.sub(lambda m:m['prefix']+pformat(float(m['val'].replace('−','-'))),s)
def bbox(axes):
 for ch in axes:
  if ch.get('id','').startswith('patch_'):
   p=ch.find('.//'+tag('path'))
   if p is not None:
    nums=[float(x) for x in re.findall(r'-?\d+(?:\.\d+)?',p.get('d',''))]
    if len(nums)>=8:return min(nums[::2]),min(nums[1::2]),max(nums[::2]),max(nums[1::2])
 return None
for src in S.rglob('*.svg'):
 if args.only and BASELINE_NAMES.get(src.stem,src.stem) not in args.only:continue
 name=BASELINE_NAMES.get(src.stem,src.stem);use=generated/(name+'.svg');tree=E.parse(str(use if use.exists() else src));svg=tree.getroot();vb=list(map(float,svg.get('viewBox').split()));fs=14*(vb[2]+(70 if name=='Figure7' else 0))/504
 for n in list(svg.iter(tag('text'))):
  old=''.join(n.itertext());new=old
  for a,b in [('Tumours','Tumors'),('Tumour','Tumor'),('tumours','tumors'),('tumour','tumor'),('Harmonised','Harmonized'),('harmonised','harmonized'),('programmes','programs'),('Labels carried from n344','Prior labels (pre-exclusion cohort)')]:new=new.replace(a,b)
  # Matplotlib mathtext splits mantissas and exponents into positioned spans.
  # Reconstruct this known source annotation explicitly, retaining its exponent.
  compact=''.join(new.split())
  if name=='FigureS2_state_structure' and compact.startswith('JTP=1.2×10'):
   new='JT P = 1.2 × 10⁻⁴⁷; ρ with immune axis = 0.78'
  else:new=ptext(new)
  new=re.sub(r'(\bP\s*<\s*)(\d+(?:\.\d+)?[eE]-\d+)',lambda m:m[1]+pformat(float(m[2])),new)
  if name=='FigureS3_immune_features' and new=='deconvolution-estimated':new='Deconvolution-estimated'
  if new!=old and '\n' in new:new=' '.join(new.split())
  if new!=old:plain(n,new);log.append([name,old,new])
 # Replace letter prefixes with independent, large bold panel letters.
 if name in ['Figure3','Figure4']:
  seen=set()
  for n in list(svg.iter(tag('text'))):
   txt=''.join(n.itertext());m=re.match(r'^([A-H])\s{2,}(.*)',txt,re.S)
   if not m:continue
   let,rest=m.groups();plain(n,rest);axes=next((a for a in n.iterancestors() if a.get('id','').startswith('axes_')),None);b=bbox(axes) if axes is not None else None
   if let not in seen and b:
    y=b[1]-fs*(.7 if name=='Figure3' else .15)
    t=E.SubElement(svg,tag('text'),x=str(max(2,b[0]-fs*1.2)),y=str(y));t.text=let;t.set('id','panel_letter_'+let);sty(t,font_family="'DejaVu Sans'",font_size=f'{fs}px',font_weight=700,text_anchor='start');seen.add(let)
 # Standardize every existing panel letter, including newly generated ones.
 for n in svg.iter(tag('text')):
  if re.fullmatch('[A-H]',''.join(n.itertext()).strip()) and not any(a.get('id','').startswith('legend_') for a in n.iterancestors()):sty(n,font_size=f'{fs}px',font_weight=700,text_anchor='start')
 if name in ['Figure3','Figure7']:
  extra_top=fs*1.2 if name=='Figure3' else 0;extra_right=70 if name=='Figure7' else 0
  svg.set('viewBox',f'{vb[0]} {-extra_top} {vb[2]+extra_right} {vb[3]+extra_top}');svg.set('width',str(vb[2]+extra_right)+'pt');svg.set('height',str(vb[3]+extra_top)+'pt')
  bg=E.Element(tag('rect'),x='0',y=str(-extra_top),width=str(vb[2]+extra_right),height=str(vb[3]+extra_top),fill='white');svg.insert(0,bg)
 if name=='FigureS4_ScPCA_axis_linkage':
  # Move B legend above the three panels; add a clean header band.
  legend=next((g for g in svg.iter(tag('g')) if g.get('id')=='legend_1'),None)
  if legend is not None:legend.set('transform','translate(-45,-72)')
  svg.set('viewBox',f'{vb[0]} -48 {vb[2]} {vb[3]+48}');svg.set('height',str(vb[3]+48)+'pt')
  bg=E.Element(tag('rect'),x='0',y='-48',width=str(vb[2]),height=str(vb[3]+48),fill='white');svg.insert(0,bg)
  # E header has three lines; move it into the inter-row gap, away from trend-P column.
  g=next(g for g in svg.iter(tag('g')) if g.get('id')=='text_120');g.set('transform','translate(0,-15)')
  for n in g.iter(tag('text')):
   txt=''.join(n.itertext())
   if txt.startswith('Immune-low:'):plain(n,'Immune-low donors: CD8 1; CD4 0; myeloid 2');sty(n,font_size='5.5px')
  # Trend P values are a column of probabilities too, using the same display rule.
  for n in svg.iter(tag('text')):
   txt=''.join(n.itertext()).strip()
   if n.get('x') and float(n.get('x'))>540 and re.fullmatch(r'0\.\d+',txt):plain(n,pformat(float(txt)))
 dest=F/src.name;dest.parent.mkdir(parents=True,exist_ok=True);tree.write(str(dest),encoding='utf-8',xml_declaration=True)
 data=E.tostring(env['outlined'](svg));cairosvg.svg2pdf(bytestring=data,write_to=str(dest.with_suffix('.pdf')));cairosvg.svg2png(bytestring=data,write_to=str(dest.with_suffix('.png')),scale=300/72)
 print('Exported',name,flush=True)
(D/'figure_text_change_log.json').write_text(json.dumps(log,ensure_ascii=False,indent=2),encoding='utf-8')
