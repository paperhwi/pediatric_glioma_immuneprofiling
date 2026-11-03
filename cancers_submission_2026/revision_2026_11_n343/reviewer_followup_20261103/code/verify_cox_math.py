from pathlib import Path
import json,numpy as np,pandas as pd
from followup_inputs import FOLLOWUP as W, REVISION, observations
O=W/'tables';a=observations()
part=pd.read_csv(REVISION/'additional_analyses_20261103/tables/per_sample_clustering_comparison.tsv',sep='\t',index_col=0);a['restandardized_state']=part.restandardized_modal.reindex(a.index);s=a[a.in_survival_analysis].copy()
def info(X,t,e,beta,weights):
 eta=X@beta;w=np.exp(eta)*weights;score=np.zeros(len(beta));I=np.zeros((len(beta),len(beta)));ll=0.
 for time in np.unique(t[e]):
  died=(t==time)&e;risk=t>=time;nd=died.sum();v=X[risk];w0=w[risk];x=X[died];wd=w[died];avg=weights[died].sum()/nd
  R0=w0.sum();D0=wd.sum();R1=(v*w0[:,None]).sum(0);D1=(x*wd[:,None]).sum(0);R2=np.einsum('n,ni,nj->ij',w0,v,v);D2=np.einsum('n,ni,nj->ij',wd,x,x)
  ll+=(eta[died]*weights[died]).sum();score+=(x*weights[died,None]).sum(0)
  for j in range(nd):
   f=j/nd;den=R0-f*D0;mean=(R1-f*D1)/den;ll-=avg*np.log(den);score-=avg*mean;I+=avg*((R2-f*D2)/den-np.outer(mean,mean))
 return ll,score,I
def matrix(q,state_col):
 X=pd.DataFrame(index=q.index)
 for st in ['Immune-high','Immune-intermediate']:X[st]=(q[state_col]==st).astype(float)
 for g in ['DMG_K27','DHG_G34','IHG']:X[g]=(q.molecular_group==g).astype(float)
 X['age_years']=q.age_at_diagnosis_years;X['male']=(q.sex=='Male').astype(float)
 for lib in sorted(a.RNA_library.unique()):
  if lib!='poly-A stranded':X['library_'+lib]=(q.RNA_library==lib).astype(float)
 return X
params=pd.read_csv(O/'categorical_cox_all_terms.tsv',sep='\t');diag=pd.read_csv(O/'categorical_cox_diagnostics.tsv',sep='\t').set_index('model');results=[]
for m,state,q,strata in [('Restandardized partition primary model','restandardized_state',s,False),('Autopsy excluded + RNA library','immune_state',s[~s.autopsy_specimen],False),('Autopsy excluded library strata','immune_state',s[~s.autopsy_specimen],True)]:
 r=params[params.model==m].set_index('term');terms=r.index.tolist();beta=r.coef.values;X=matrix(q,state)[terms];t=q.OS_days_corrected.values;e=q.OS_event.values.astype(bool);ll=0.;score=np.zeros(len(beta));I=np.zeros((len(beta),len(beta)))
 for lib in q.RNA_library.unique() if strata else [None]:
  mask=(q.RNA_library==lib).values if strata else np.ones(len(q),bool);v,g,h=info(X.values[mask],t[mask],e[mask],beta,np.ones(mask.sum()));ll+=v;score+=g;I+=h
 se=np.sqrt(np.diag(np.linalg.inv(I)));result={'model':m,'loglik_abs_difference':float(abs(ll-diag.loc[m,'log_partial_likelihood'])),'max_abs_independent_score':float(abs(score).max()),'max_abs_SE_difference':float(abs(se-r['se(coef)'].values).max())};result['passed']=result['loglik_abs_difference']<1e-7 and result['max_abs_independent_score']<1e-3 and result['max_abs_SE_difference']<1e-7;assert result['passed'],result;results.append(result)
(O/'independent_new_cox_verification.json').write_text(json.dumps({'method':'Independent NumPy Efron risk-set likelihood, score and inverse-information standard errors; library strata accumulated separately','models':results,'passed':all(x['passed'] for x in results)},indent=2),encoding='utf-8');print(json.dumps(results,indent=2))
