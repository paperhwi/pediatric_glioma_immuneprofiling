"""Reviewer-requested categorical-state Cox sensitivities; full corrected follow-up."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
from pathlib import Path
import json,warnings
import numpy as np,pandas as pd
from scipy.stats import chi2
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

from followup_inputs import FOLLOWUP as work, REVISION, observations
a=observations()
follow=REVISION/'additional_analyses_20261103'
part=pd.read_csv(follow/'tables/per_sample_clustering_comparison.tsv',sep='\t',index_col=0)
a['restandardized_state']=part.restandardized_modal.reindex(a.index)
assert len(a)==343 and a.restandardized_state.notna().all()
assert a.restandardized_state.value_counts().to_dict()=={'Immune-intermediate':149,'Immune-high':116,'Immune-low':78}
s=a[a.in_survival_analysis].copy()
assert len(s)==245 and s.OS_event.sum()==206
assert len(s[~s.autopsy_specimen])==189 and s.loc[~s.autopsy_specimen,'OS_event'].sum()==150

plan={'unit':'one RNA specimen per patient; 343 eligible patients','outcome':'corrected OS days from relevant HGG diagnosis; original death/censoring status','followup':'full observed follow-up; only KM display limited to 5 years','references':{'state':'Immune-low','molecular_group':'pHGG_WT','sex':'Female','RNA_library':'poly-A stranded'},'covariates':['two categorical immune-state indicators','DMG_K27, DHG_G34, IHG indicators','age at diagnosis, years','male indicator'],'estimation':'unpenalized Cox proportional hazards, Efron ties, Wald 95% CI and two-sided P','additional_models':['same primary Cox with restandardized 116/149/78 partition','original partition, autopsy excluded AND categorical RNA library adjusted','original partition, autopsy excluded AND library-stratified baseline hazards'],'interpretation':'reviewer-requested sensitivity analyses; diagnostic P values unadjusted; no selection by significance'}
(work/'analysis_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')

def design(sub,state_col='immune_state',library='none'):
    d=pd.DataFrame({'time':sub.OS_days_corrected,'event':sub.OS_event},index=sub.index)
    for state in ['Immune-high','Immune-intermediate']:d[state]=(sub[state_col]==state).astype(float)
    for g in ['DMG_K27','DHG_G34','IHG']:d[g]=(sub.molecular_group==g).astype(float)
    d['age_years']=sub.age_at_diagnosis_years;d['male']=(sub.sex=='Male').astype(float)
    if library=='covariate':
        for lib in sorted(a.RNA_library.unique()):
            if lib!='poly-A stranded':d['library_'+lib]=(sub.RNA_library==lib).astype(float)
    if library=='strata':d['RNA_library']=sub.RNA_library
    assert not d.isna().any().any() and (d.time>0).all()
    return d

models=[('Primary original partition',s,'immune_state','none'),('Restandardized partition primary model',s,'restandardized_state','none'),('Original partition + RNA library',s,'immune_state','covariate'),('Original partition library strata',s,'immune_state','strata'),('Autopsy excluded original partition',s[~s.autopsy_specimen],'immune_state','none'),('Autopsy excluded + RNA library',s[~s.autopsy_specimen],'immune_state','covariate'),('Autopsy excluded library strata',s[~s.autopsy_specimen],'immune_state','strata')]
terms=[];rows=[];diagnostics=[];fits={}
for name,sub,state,lib in models:
    d=design(sub,state,lib);kw={'strata':['RNA_library']} if lib=='strata' else {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always');c=CoxPHFitter().fit(d,'time','event',**kw)
    ph=proportional_hazard_test(c,d,time_transform='rank').summary
    tab=c.summary.reset_index().rename(columns={'covariate':'term'});tab.insert(0,'model',name);terms.append(tab)
    x=d.drop(columns=['time','event','RNA_library'],errors='ignore')
    diagnostics.append({'model':name,'n':len(d),'deaths':int(d.event.sum()),'design_rank':int(np.linalg.matrix_rank(np.c_[np.ones(len(x)),x])),'design_columns_with_intercept':x.shape[1]+1,'log_partial_likelihood':float(c.log_likelihood_),'warnings':' | '.join(str(w.message) for w in caught)})
    for term in ['Immune-high','Immune-intermediate']:
        r=c.summary.loc[term]
        rows.append({'model':name,'state_partition':state,'RNA_library_handling':lib,'n':len(d),'deaths':int(d.event.sum()),'contrast':term+' vs Immune-low','HR':r['exp(coef)'],'CI95_lower':r['exp(coef) lower 95%'],'CI95_upper':r['exp(coef) upper 95%'],'Wald_P':r['p'],'PH_rank_P':ph.loc[term,'p'],'n_high':int((sub[state]=='Immune-high').sum()),'n_intermediate':int((sub[state]=='Immune-intermediate').sum()),'n_low':int((sub[state]=='Immune-low').sum())})
    fits[name]=c
table=pd.DataFrame(rows)
public=work/'tables';public.mkdir(exist_ok=True)
table.to_csv(public/'categorical_state_sensitivity_cox.tsv',sep='\t',index=False)
pd.concat(terms).to_csv(public/'categorical_cox_all_terms.tsv',sep='\t',index=False)
pd.DataFrame(diagnostics).to_csv(public/'categorical_cox_diagnostics.tsv',sep='\t',index=False)
summary={'plan':plan,'models':rows,'diagnostics':diagnostics,'baseline_HR':float(table.iloc[0].HR),'input_clinical_file_not_redistributed':True}
(public/'categorical_cox_result.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(table[['model','n','deaths','contrast','HR','CI95_lower','CI95_upper','Wald_P','PH_rank_P']].to_string(index=False))
