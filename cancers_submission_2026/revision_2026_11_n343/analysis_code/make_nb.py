import nbformat as nbf
nb=nbf.v4.new_notebook(); C=[]
md=lambda s: C.append(nbf.v4.new_markdown_cell(s)); code=lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Rev343 — final n = 343 analysis (3 Nov 2026)
Single exclusion rule (6 non-diffuse-HGG specimens removed before any analysis), immune states re-derived with the locked procedure (modal label over consensus seeds 42–66), new state names (Immune-high / Immune-intermediate / Immune-low), one OS-origin rule applied to every patient, pre-specified primary Cox model.

Fast steps are executed here. Slow steps (multi-seed clustering, bootstrap optimism, sweeps, ScPCA reference build, permutation tests) were run as scripts with the same code; with `RUN_SLOW = False` the notebook reads their saved outputs. Section 9 recomputes the key numbers independently (different code paths) and checks them against the outputs.""")
code("""import subprocess, sys, os, json, glob, pandas as pd, numpy as np
os.chdir('/home/claude/q343'); RUN_SLOW=False
def run(s, slow=False):
    if slow and not RUN_SLOW: print(f'[slow step, outputs already present] python3 src/{s}'); return
    r=subprocess.run([sys.executable,f'src/{s}'],capture_output=True,text=True,env={**os.environ,'OMP_NUM_THREADS':'1'}); print(r.stdout[-3000:]); assert r.returncode==0, r.stderr[-3000:]
pd.set_option('display.width',220); pd.set_option('display.max_columns',30)""")
md("## 1. Clustering seeds (slow) and cohort/labels"); code("print('python3 src/clust.py <ids> <seeds> <out>  (seeds 0-99 for n = 343 and n = 344)')\nrun('s00_cohort.py')")
md("## 2. Data objects and OS-origin audit (all 343 tumours)"); code("run('s01_data.py'); run('s02_os_audit.py')")
md("## 3. Survival: primary model, sensitivity analyses, autopsy table, follow-up, immune axis, within subtype"); code("run('s03_survival.py'); run('s04_incremental.py', slow=True)\npd.read_csv('out/T_survival_sensitivity.tsv',sep='\\t').round(3)")
code("pd.read_csv('out/T_autopsy_by_state.tsv',sep='\\t')"); code("pd.read_csv('out/T_followup.tsv',sep='\\t',index_col=0).T")
md("## 4. Figures 1–7 and Table 1"); code("for s in ['s05_fig4.py','s06_table1.py','s07_deg.py','s08_fig7.py','s10_fig6.py','s11_fig2.py']: run(s)\nrun('s12_fig3_sweep.py',slow=True); run('s13_fig3_plot.py'); run('s19_visium.py'); run('s20a_build_scpca_reference.py',slow=True); run('s20b_fig5b_cell_source.py',slow=True); run('s23_fig5.py'); run('s24_fig1.py')")
code("from IPython.display import Image, display\nfor f in ['Figure1','Figure2','Figure3','Figure4','Figure5','Figure6','Figure7']: display(Image(f'fig/{f}/{f}.png',width=900))")
md("## 5. Structure and robustness (supplement)"); code("for s in ['s09_within_group.py','s15_k2_vs_k3.py','s16_library.py']: run(s,slow=True)\nfor s in ['s21_deconv_frameworks_evasion.py','s27_crossmethod.py']: run(s)\nprint(json.dumps(json.load(open('out/r54/R54_summary.json')),indent=0)[:1500])")
md("## 6. Immune features, programmes, gene screen (supplement)"); code("run('s17_feature_cox.py',slow=True); run('s18_programmes_by_state.py',slow=True); run('s28_gene_screen.py'); run('s29_published_programmes.py')")
md("## 7. Single cell (ScPCA linkage with the n = 343 reference; per-patient table)"); code("run('s22_scpca_axis_relink.py',slow=True); run('s26_singlecell_donor_table.py'); run('s31_fig_c3.py')")
md("## 8. Supplementary figures, numbers file, workbooks, h5ad"); code("for s in ['s25_figS_survival.py','s30_figS_structure_features.py','s32_numbers.py','s33_supp_workbook.py','s34_h5ad.py']: run(s)\nfor f in ['FigureS_survival','FigureS_state_structure','FigureS_immune_features']: display(Image(f'fig/{f}/{f}.png',width=800))")
md("## 9. Independent verification\nDifferent code paths (statsmodels PHReg instead of lifelines; scipy χ²; direct counting) on the delivered h5ad.")
code("""import anndata as ad, statsmodels.api as sm
from scipy import stats
T=ad.read_h5ad('deliv/Rev343_tumours_n343.h5ad'); o=T.obs.copy(); chk=[]
def c(name,ok,detail): chk.append(dict(check=name,ok=bool(ok),detail=detail))
c('n tumours',len(o)==343,len(o)); vc=o.immune_state.value_counts(); c('state sizes',(vc['Immune-high'],vc['Immune-intermediate'],vc['Immune-low'])==(122,144,77),vc.to_dict())
c('no excluded specimen present',not o.index.isin(pd.read_csv('out/T_excluded_specimens.tsv',sep='\\t').Kids_First_Biospecimen_ID).any(),'6 excluded')
s=o[o.in_survival_analysis].copy(); X=pd.DataFrame({'high':(s.immune_state=='Immune-high')*1.,'int':(s.immune_state=='Immune-intermediate')*1.,'DMG':(s.molecular_group=='DMG_K27')*1.,'DHG':(s.molecular_group=='DHG_G34')*1.,'IHG':(s.molecular_group=='IHG')*1.,'age':s.age_at_diagnosis_years,'male':(s.sex=='Male')*1.})
m=sm.PHReg(s.OS_days_corrected.values,X.values,status=s.OS_event.values,ties='efron').fit(); hr=np.exp(m.params[0]); ci=np.exp(m.conf_int()[0]); 
ref=pd.read_csv('out/T_primary_cox.tsv',sep='\\t',index_col=0).loc['Immune-high vs Immune-low']
c('primary HR (statsmodels vs lifelines)',abs(hr-ref['exp(coef)'])<1e-4,f'{hr:.4f} ({ci[0]:.3f}-{ci[1]:.3f}) P={m.pvalues[0]:.4f} vs {ref["exp(coef)"]:.4f}')
c('survival set',(len(s),int(s.OS_event.sum()))==(245,206),(len(s),int(s.OS_event.sum())))
na=s[~s.autopsy_specimen]; Xa=X.loc[na.index]; ma=sm.PHReg(na.OS_days_corrected.values,Xa.values,status=na.OS_event.values,ties='efron').fit()
r=pd.read_csv('out/T_survival_sensitivity.tsv',sep='\\t'); c('no-autopsy HR',abs(np.exp(ma.params[0])-r.high_HR[1])<1e-4,f'{np.exp(ma.params[0]):.3f} n={len(na)}')
Xx=X[['DMG','DHG','IHG','age','male']].assign(axis=s.immune_axis); mx=sm.PHReg(s.OS_days_corrected.values,Xx.values,status=s.OS_event.values,ties='efron').fit()
c('immune axis HR',abs(np.exp(mx.params[-1])-0.8273)<5e-4,f'{np.exp(mx.params[-1]):.4f} P={mx.pvalues[-1]:.4f}')
ct=pd.crosstab(o.molecular_group,o.immune_state); chi=stats.chi2_contingency(ct)[0]; c('Fig 4B chi2',abs(chi-30.6138)<1e-3,round(chi,4))
deg=pd.read_csv('out/T_adjusted_DEG_n343_all_contrasts.tsv',sep='\\t'); up=lambda k:set(deg[(deg.contrast==k)&(deg.q_BH<.05)&(deg.adjusted_log2FC>1)].gene)
c('shared genes',len(up('Immune-high vs Immune-low')&up('Immune-intermediate vs Immune-low'))==404,len(up('Immune-high vs Immune-low')&up('Immune-intermediate vs Immune-low')))
au=pd.read_csv('out/T_OS_origin_audit_all_n343.tsv',sep='\\t'); c('OS audit covers every tumour',len(au)==343 and (au.category!='unclassified').all(),au.action.value_counts().to_dict())
c('autopsy among survival set',int(s.autopsy_specimen.sum())==56,int(s.autopsy_specimen.sum()))
V=pd.DataFrame(chk); V.to_csv('deliv/verification_checks.tsv',sep='\\t',index=False); print(V.to_string()); assert V.ok.all()""")
nb["cells"]=C; nbf.write(nb,"deliv/Rev343_final_pipeline.ipynb")
