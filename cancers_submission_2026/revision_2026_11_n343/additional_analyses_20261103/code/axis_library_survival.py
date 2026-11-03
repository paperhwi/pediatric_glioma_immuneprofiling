"""Continuous immune axis survival association with categorical RNA-library adjustment.
Keep the delivered n343 axis, corrected outcomes, age, molecular group, sex, and per-cohort-SD unit fixed.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,warnings,shutil
import numpy as np,pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import chi2,spearmanr
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from threadpoolctl import threadpool_limits

WORK=Path(__file__).parent
from analysis_paths import OUT
from private_inputs import load_survival_observations, private_output_dir
threadpool_limits(limits=1)

def design(s,library='none',axis=True,reference='poly-A stranded'):
    d=pd.DataFrame({'time':s.OS_days_corrected,'event':s.OS_event},index=s.index)
    for g in ['DMG_K27','DHG_G34','IHG']:d[g]=(s.molecular_group==g).astype(float)
    d['age_years']=s.age_at_diagnosis_years;d['male']=(s.sex=='Male').astype(float)
    if axis:d['immune_axis']=s.immune_axis
    if library=='covariate':
        for L in sorted(s.RNA_library.unique()):
            if L!=reference:d['library_'+L]=(s.RNA_library==L).astype(float)
    if library=='strata':d['RNA_library']=s.RNA_library
    assert not d.isna().any().any() and (d.time>0).all()
    return d

def main():
    a=load_survival_observations()
    F=pd.read_csv(OUT/'inputs/retained_z_n343.tsv',sep='\t',index_col=0)
    nonimm=['ssGSEA_Cell_Cycle_Proliferation','ssGSEA_Stemness_Brain_Tumor','ssGSEA_MAPK_Activity']
    cols=[c for c in F if c.startswith('ssGSEA') and c not in nonimm]
    X=(F[cols]-F[cols].mean())/F[cols].std(ddof=1)
    pca=PCA(5).fit(X.values);pc=pca.transform(X.values)[:,0]
    if np.corrcoef(pc,F.ssGSEA_T_Cell_Cytotoxicity)[0,1]<0:pc=-pc
    pc=(pc-pc.mean())/pc.std(ddof=1)
    axis_check={'features':len(cols),'PC1_variance_explained':float(pca.explained_variance_ratio_[0]),'corr_with_delivered_axis':float(np.corrcoef(pc,a.immune_axis)[0,1]),'max_abs_diff':float(np.max(abs(pc-a.immune_axis.values))),'SD_definition':'sample SD across all 343, not refit in survival subset'}
    assert axis_check['max_abs_diff']<1e-10,axis_check
    a['axis_recomputed_check']=pc
    s=a[a.in_survival_analysis].copy()
    assert len(s)==245 and s.OS_event.sum()==206 and s.sex.isin(['Female','Male']).all()
    assert s.RNA_library.notna().all() and not s.RNA_library.isin(['nan','Unknown','Not Reported']).any()
    ref=s.RNA_library.value_counts().idxmax()
    assert ref=='poly-A stranded'
    corrected=s[s.OS_corrected_record]
    assert len(corrected)==4
    plan={'outcome':'delivered corrected OS days with unchanged death/censoring indicator','primary_predictor':'delivered immune axis per 1 sample SD across n343','adjustments':['molecular group, reference pHGG_WT','age at diagnosis years','male, reference female'],
          'additional_adjustment':f'RNA library categorical, reference {ref}; 4 levels, 3 indicators','sensitivity':'RNA-library-stratified baseline hazards, then repeat the same three models without autopsy specimens','cohort':'survival-evaluable n245, complete on all modeled covariates','method':'unpenalized Cox, Efron ties, Wald 95% CI and two-sided P; nested LRT on identical rows and strata','PH':'scaled Schoenfeld rank-time tests; diagnostic P values unadjusted'}
    (OUT/'code/survival_analysis_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    summaries=[];terms=[];ph=[];warnings_all=[];diagnostics=[]
    libraryrows=[]
    for subset,sub in [('All survival-evaluable',s),('Excluding autopsy',s[~s.autopsy_specimen])]:
        for L,g in sub.groupby('RNA_library'):
            libraryrows.append({'subset':subset,'RNA_library':L,'n':len(g),'deaths':int(g.OS_event.sum()),'censored':int((1-g.OS_event).sum()),'axis_mean':g.immune_axis.mean(),'axis_SD':g.immune_axis.std(),'axis_median':g.immune_axis.median(),'axis_Q1':g.immune_axis.quantile(.25),'axis_Q3':g.immune_axis.quantile(.75)})
        for typ in ['none','covariate','strata']:
            model=f'{subset} / '+{'none':'Molecular group + age + sex','covariate':'Molecular group + age + sex + RNA library','strata':'Molecular group + age + sex; RNA-library strata'}[typ]
            d=design(sub,typ,reference=ref);d0=design(sub,typ,axis=False,reference=ref)
            kwargs={'strata':['RNA_library']} if typ=='strata' else {}
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always');c=CoxPHFitter().fit(d,'time','event',**kwargs);c0=CoxPHFitter().fit(d0,'time','event',**kwargs)
            warning_text=' | '.join(str(x.message) for x in caught)
            warnings_all.extend({'model':model,'warning':str(w.message),'class':w.category.__name__} for w in caught)
            phs=proportional_hazard_test(c,d,time_transform='rank').summary
            r=c.summary.loc['immune_axis'];stat=2*(c.log_likelihood_-c0.log_likelihood_)
            row={'model':model,'subset':subset,'library_handling':typ,'n':len(d),'deaths':int(d.event.sum()),'HR_per_cohort_SD':r['exp(coef)'],'CI95_lower':r['exp(coef) lower 95%'],'CI95_upper':r['exp(coef) upper 95%'],'Wald_P':r['p'],'axis_LRT_chi2':stat,'axis_LRT_df':1,'axis_LRT_P':chi2.sf(stat,1),'axis_PH_rank_P':phs.loc['immune_axis','p'],'log_partial_likelihood':c.log_likelihood_,'partial_AIC':c.AIC_partial_,'n_parameters':len(c.params_),'library_reference':ref if typ=='covariate' else '', 'fit_warnings':warning_text}
            summaries.append(row)
            term=c.summary.reset_index().rename(columns={'covariate':'term'});term.insert(0,'model',model);terms.append(term)
            phf=phs.reset_index().rename(columns={'index':'term'});phf.insert(0,'model',model);ph.append(phf)
            predictors=d.drop(columns=['time','event','RNA_library'],errors='ignore');A=predictors.drop(columns='immune_axis');Az=(A-A.mean())/A.std();Az=np.c_[np.ones(len(Az)),Az.values]
            coef=np.linalg.lstsq(Az,d.immune_axis.values,rcond=None)[0];res=d.immune_axis.values-Az@coef
            r2=1-res.dot(res)/np.sum((d.immune_axis-d.immune_axis.mean())**2)
            diagnostics.append({'model':model,'design_rank':int(np.linalg.matrix_rank(np.c_[np.ones(len(predictors)),predictors.values])),'design_columns_with_intercept':len(predictors.columns)+1,'standardized_condition_number':float(np.linalg.cond(np.c_[np.ones(len(predictors)),((predictors-predictors.mean())/predictors.std()).values])),'axis_R2_on_covariates':r2,'axis_VIF':1/(1-r2),'library_strata_omitted_from_VIF':typ=='strata'})
            if subset=='All survival-evaluable' and typ=='covariate':
                resid=c.compute_residuals(d,'scaled_schoenfeld');resid['time']=d.time.loc[resid.index]
                resid[['immune_axis','time']].to_csv(private_output_dir()/'axis_scaled_Schoenfeld_library_adjusted.tsv',sep='\t')
                d.to_csv(private_output_dir()/'cox_design_axis_library_adjusted_n245.tsv',sep='\t')
                np.savez_compressed(OUT/'tables/axis_library_covariance.npz',covariance=c.variance_matrix_.values,terms=np.array(c.variance_matrix_.index,dtype=str))
    table=pd.DataFrame(summaries)
    original=OUT/'inputs/reference/T_immune_axis_cox.tsv'
    old=pd.read_csv(original,sep='\t').iloc[0]
    base=table.iloc[0]
    for new,prior in [('HR_per_cohort_SD','HR'),('CI95_lower','lo'),('CI95_upper','hi'),('Wald_P','p')]:assert abs(base[new]-old[prior])<1e-10,(new,base[new],old[prior])
    pd.DataFrame(libraryrows).to_csv(OUT/'tables/library_descriptives.tsv',sep='\t',index=False)
    table.to_csv(OUT/'tables/immune_axis_library_cox.tsv',sep='\t',index=False)
    pd.concat(terms).to_csv(OUT/'tables/cox_all_terms.tsv',sep='\t',index=False)
    pd.concat(ph).to_csv(OUT/'tables/PH_rank_tests_all_terms.tsv',sep='\t',index=False)
    pd.DataFrame(diagnostics).to_csv(OUT/'tables/cox_design_diagnostics.tsv',sep='\t',index=False)
    results={'axis_check':axis_check,'baseline_matches_delivered_1e-10':True,'library_reference':ref,'missing_covariates_n245':0,'corrected_OS_records_verified':4,'fit_warnings':warnings_all,'models':summaries}
    (OUT/'tables/survival_result.json').write_text(json.dumps(results,indent=2),encoding='utf-8')

    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
