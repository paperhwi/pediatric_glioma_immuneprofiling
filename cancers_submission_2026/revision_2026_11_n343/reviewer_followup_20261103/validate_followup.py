from pathlib import Path
import json,csv,hashlib
import numpy as np,pandas as pd
F=Path(__file__).resolve().parent
C=pd.read_csv(F/'tables/categorical_state_sensitivity_cox.tsv',sep='\t')
assert len(C)==14 and set(C.n)=={245,189} and set(C.deaths)=={206,150}
for model,hr in [('Restandardized partition primary model',0.596108636981),('Autopsy excluded + RNA library',0.7585558)]:
 r=C[(C.model==model)&C.contrast.str.startswith('Immune-high')].iloc[0];assert abs(r.HR-hr)<1e-6
I=pd.read_csv(F/'tables/ipw_cox.tsv',sep='\t');assert len(I)==10 and (I.n==245).all() and (I.deaths==206).all()
row=I[(I.scheme=='Stabilized IPW 1–99% capped')&I.contrast.str.startswith('Immune-high')].iloc[0];assert abs(row.HR-.599544)<1e-6
assert json.loads((F/'tables/independent_new_cox_verification.json').read_text())['passed']
K=json.loads((F/'tables/KM_display_fit_audit.json').read_text());assert len(K)==9 and all(r['display_max_years']==5 for r in K) and max(r['fit_max_years'] for r in K)>8
S=pd.read_csv(F/'tables/FigureS1_sensitivity_source.tsv',sep='\t');assert len(S)==17 and S.iloc[-3:].intermediate_HR.notna().all()
with (F/'MANIFEST_SHA256.tsv').open(encoding='utf-8') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  p=F/row['Path'];assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==row['SHA256'],row['Path']
print('OK: additional categorical Cox and IPW models; full-follow-up KM fits displayed to 5 years; exact payload hashes')
