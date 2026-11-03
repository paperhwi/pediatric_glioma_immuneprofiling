"""Independent Efron partial-likelihood, score and Hessian check of the library-adjusted fit."""
from pathlib import Path
import json,shutil
import numpy as np,pandas as pd
work=Path(__file__).parent
from analysis_paths import OUT as out
from private_inputs import private_output_dir
d=pd.read_csv(private_output_dir()/'cox_design_axis_library_adjusted_n245.tsv',sep='\t',index_col=0)
raw=pd.read_csv(out/'tables/cox_all_terms.tsv',sep='\t')
model='All survival-evaluable / Molecular group + age + sex + RNA library'
r=raw[raw.model==model].set_index('term');terms=[c for c in d if c not in ['time','event']]
X=d[terms].values;beta=r.loc[terms,'coef'].values;t=d.time.values;e=d.event.values.astype(bool);w=np.exp(X@beta)
score=np.zeros(len(terms));info=np.zeros((len(terms),len(terms)));ll=0.
for time in np.unique(t[e]):
    deaths=(t==time)&e;risk=t>=time;nd=int(deaths.sum())
    dx=X[deaths];dw=w[deaths];rx=X[risk];rw=w[risk]
    R0=rw.sum();D0=dw.sum();R1=(rx*rw[:,None]).sum(0);D1=(dx*dw[:,None]).sum(0)
    R2=np.einsum('n,ni,nj->ij',rw,rx,rx);D2=np.einsum('n,ni,nj->ij',dw,dx,dx)
    ll+=(dx@beta).sum();score+=dx.sum(0)
    for j in range(nd):
        fraction=j/nd;den=R0-fraction*D0;mean=(R1-fraction*D1)/den
        ll-=np.log(den);score-=mean;info+=(R2-fraction*D2)/den-np.outer(mean,mean)
reported=pd.read_csv(out/'tables/immune_axis_library_cox.tsv',sep='\t');expected=reported.loc[reported.model==model,'log_partial_likelihood'].iloc[0]
saved=np.load(out/'tables/axis_library_covariance.npz');order=[list(saved['terms']).index(c) for c in terms];cov=saved['covariance'][np.ix_(order,order)]
result={'method':'independently coded NumPy Efron risk-set likelihood, score and observed-information matrix','log_likelihood_abs_difference':float(abs(ll-expected)),
        'max_abs_score':float(abs(score).max()),'covariance_max_abs_difference':float(abs(np.linalg.inv(info)-cov).max()),'passed':bool(abs(ll-expected)<1e-8 and abs(score).max()<1e-3 and abs(np.linalg.inv(info)-cov).max()<1e-8)}
assert result['passed'],result
(out/'tables/independent_cox_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
