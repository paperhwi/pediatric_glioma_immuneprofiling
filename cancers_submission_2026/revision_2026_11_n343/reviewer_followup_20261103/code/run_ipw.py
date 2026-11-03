from pathlib import Path
import json,warnings
import numpy as np,pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter
from scipy.stats import chi2_contingency
from openpyxl import load_workbook
from followup_inputs import FOLLOWUP as W, observations
O=W/'tables';O.mkdir(exist_ok=True)
a=observations();y=a.in_survival_analysis.astype(int);s=a[y==1].copy()
assert len(a)==343 and len(s)==245 and s.OS_event.sum()==206
# Original submitted revD propensity specification; no survival outcomes used in propensity fitting.
state_for_propensity=pd.Categorical(a.immune_state,categories=['Immune-low','Immune-high','Immune-intermediate'])
X=pd.get_dummies(pd.DataFrame({'immune_state':state_for_propensity,'molecular_group':a.molecular_group.astype(str),'location':a.location_class.fillna('Unknown').astype(str)},index=a.index),drop_first=True).astype(float)
X['age']=a.age_at_diagnosis_years.fillna(a.age_at_diagnosis_years.median());X['male']=(a.sex=='Male').astype(float)
with warnings.catch_warnings(record=True) as wr:
 warnings.simplefilter('always');prop=LogisticRegression(C=1.0,max_iter=5000,solver='lbfgs').fit(StandardScaler().fit_transform(X),y).predict_proba(StandardScaler().fit_transform(X))[:,1]
w=y.mean()/prop;lo,hi=np.percentile(w,[1,99]);clip=np.clip(w,lo,hi)
Xsimple=pd.get_dummies(a.molecular_group,drop_first=True).astype(float);p2=LogisticRegression(C=1,max_iter=5000,solver='lbfgs').fit(Xsimple,y).predict_proba(Xsimple)[:,1];w2=y.mean()/p2[y==1];w2=np.clip(w2,*np.percentile(w2,[1,99]))
def design(s):
 d=pd.DataFrame({'time':s.OS_days_corrected,'event':s.OS_event},index=s.index)
 for state in ['Immune-high','Immune-intermediate']:d[state]=(s.immune_state==state).astype(float)
 for g in ['DMG_K27','DHG_G34','IHG']:d[g]=(s.molecular_group==g).astype(float)
 d['age_years']=s.age_at_diagnosis_years;d['male']=(s.sex=='Male').astype(float);assert not d.isna().any().any();return d
schemes=[('Unweighted primary',None),('Stabilized IPW 1–99% capped',clip[y==1]),('Stabilized IPW uncapped',w[y==1]),('Stabilized IPW 2.5–97.5% capped',np.clip(w[y==1],*np.percentile(w[y==1],[2.5,97.5]))),('Molecular-group-only propensity 1–99% capped',w2)]
rows=[];terms=[];weightrows=[]
for label,weights in schemes:
 d=design(s);kw={}
 if weights is not None:d['weight']=weights;kw={'weights_col':'weight','robust':True}
 c=CoxPHFitter().fit(d,'time','event',**kw);t=c.summary.reset_index().rename(columns={'covariate':'term'});t.insert(0,'scheme',label);terms.append(t)
 for state in ['Immune-high','Immune-intermediate']:
  r=c.summary.loc[state];rows.append({'scheme':label,'contrast':state+' vs Immune-low','n':len(d),'deaths':int(d.event.sum()),'HR':r['exp(coef)'],'CI95_lower':r['exp(coef) lower 95%'],'CI95_upper':r['exp(coef) upper 95%'],'Wald_P':r['p'],'variance':'sandwich robust; estimated weights treated as fixed' if weights is not None else 'model based'})
 if weights is not None:weightrows.append({'scheme':label,'n_observed':len(weights),'mean':float(np.mean(weights)),'min':float(np.min(weights)),'max':float(np.max(weights)),'effective_sample_size':float(np.sum(weights)**2/np.sum(weights**2))})
pd.DataFrame(rows).to_csv(O/'ipw_cox.tsv',sep='\t',index=False);pd.concat(terms).to_csv(O/'ipw_cox_all_terms.tsv',sep='\t',index=False);pd.DataFrame(weightrows).to_csv(O/'ipw_weight_diagnostics.tsv',sep='\t',index=False)
def smd(obs,target,weights=None):
 m=np.average(obs,axis=0,weights=weights);v=np.average((obs-m)**2,axis=0,weights=weights);mt=target.mean(axis=0);vt=target.var(axis=0);den=np.sqrt((v+vt)/2);return (m-mt)/np.where(den==0,1,den)
bal=[]
for targetname,target in [('Full eligible cohort target',X.values),('Legacy observed-versus-missing comparison',X.values[y==0])]:
 before=smd(X.values[y==1],target);after=smd(X.values[y==1],target,clip[y==1])
 for col,b,c in zip(X.columns,before,after):bal.append({'comparison':targetname,'covariate':col,'SMD_unweighted':b,'SMD_IPW_1_99':c})
pd.DataFrame(bal).to_csv(O/'ipw_balance.tsv',sep='\t',index=False)
missing=a.groupby('molecular_group').in_survival_analysis.agg(['size','sum']).rename(columns={'size':'n_eligible','sum':'n_OS_evaluable'});missing['n_missing']=missing.n_eligible-missing.n_OS_evaluable;missing['missing_percent']=100*missing.n_missing/missing.n_eligible;missing.to_csv(O/'ipw_missingness.tsv',sep='\t')
plan={'source':'Submitted revD_IPW_sensitivity.ipynb, code cells 5, 9 and 11','eligible':343,'observed':245,'missing':98,'propensity_covariates':list(X.columns),'age_missing_imputation':'median across 343, only propensity model; all 245 Cox ages observed','logistic':'L2 C=1.0 lbfgs max_iter=5000; all dummy/continuous columns standardized ddof=0','stabilized_weight':'(245/343)/P(OS observed|X)','cap_1_99':[float(lo),float(hi)],'cap_scope':'1st/99th percentiles over all 343, as submitted protocol; 2.5/97.5 and sparse-model cutoffs over observed patients, as legacy code','Cox':'unpenalized, corrected full follow-up; state + molecular group + age + sex; robust sandwich SE for weighted models','limitation':'Conditional MAR sensitivity only. Weights estimated in the same sample; sandwich SE treats weights as fixed. Residual covariate imbalance and informative missingness may remain; tied event times also limit sandwich approximation. No claim of robustness to MNAR.'}
plan['propensity_state_reference']='Immune-low, preserving the historical Immune-desert reference; dummy order High then Intermediate.'
(O/'ipw_analysis_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')

print(pd.DataFrame(rows).to_string(index=False))
