"""Check stored follow-up evidence and its public-data boundary without fitting models."""
from pathlib import Path
import hashlib,json,csv
import pandas as pd
import numpy as np

OUT=Path(__file__).resolve().parent
public=pd.read_csv(OUT/'inputs/public_annotations_n343.tsv',sep='\t',index_col=0)
assert public.shape==(343,4) and public.index.is_unique
assert list(public.columns)==['immune_state','molecular_group','RNA_library','immune_axis']
assert public.immune_state.value_counts().to_dict()=={'Immune-intermediate':144,'Immune-high':122,'Immune-low':77}
F=pd.read_csv(OUT/'inputs/retained_z_n343.tsv',sep='\t',index_col=0)
Z=pd.read_csv(OUT/'inputs/restandardized_z_n343.tsv',sep='\t',index_col=0)
assert F.shape==Z.shape==(343,46) and F.index.equals(public.index) and Z.index.equals(F.index)
assert np.allclose(Z.values,((F-F.mean())/F.std(ddof=1)).values,rtol=0,atol=1e-12)
assert abs(Z.mean()).max()<1e-12 and abs(Z.std()-1).max()<1e-12
a=pd.read_csv(OUT/'tables/per_sample_clustering_comparison.tsv',sep='\t',index_col=0).loc[public.index]
assert a.immune_state.equals(public.immune_state)
assert a.retained_modal.equals(public.immune_state)
assert a.restandardized_modal.value_counts().to_dict()=={'Immune-intermediate':149,'Immune-high':116,'Immune-low':78}
assert int((a.immune_state!=a.restandardized_modal).sum())==9
C=json.loads((OUT/'tables/clustering_result.json').read_text(encoding='utf-8'))
assert C['baseline_modal_reproduced'] and C['changed_samples']==9
assert abs(C['modal_agreement']-334/343)<1e-12 and abs(C['modal_ARI']-0.9149481671210321)<1e-12
for kind in ['retained','restandardized']:
    S=pd.read_csv(OUT/f'tables/assignments_25seeds_{kind}.tsv',sep='\t',index_col=0).loc[public.index]
    assert list(S.columns)==[str(i) for i in range(42,67)]
    assert S.mode(axis=1)[0].equals(a[f'{kind}_modal'])
cox=pd.read_csv(OUT/'tables/immune_axis_library_cox.tsv',sep='\t')
assert len(cox)==6 and cox.n.tolist()==[245,245,245,189,189,189] and cox.deaths.tolist()==[206,206,206,150,150,150]
for i,hr,p in [(0,.827259622265942,.009536676901466312),(1,.8614425806677978,.08794659652038442),(2,.8545567162497281,.07967785463389544),(4,.942358366570694,.5749706758021271)]:
    assert abs(cox.iloc[i].HR_per_cohort_SD-hr)<1e-9 and abs(cox.iloc[i].Wald_P-p)<1e-9
assert cox.iloc[1].CI95_lower<1<cox.iloc[1].CI95_upper
assert json.loads((OUT/'tables/independent_cox_verification.json').read_text(encoding='utf-8'))['passed']

# Exact input allowlist: no individual outcomes, event-time residuals or patient covariates.
expected={'inputs/public_annotations_n343.tsv','inputs/retained_z_n343.tsv','inputs/restandardized_z_n343.tsv',
          'inputs/reference/T_immune_axis_cox.tsv','inputs/reference/T_seed_partitions_n343.tsv','inputs/reference/ksweep_n343.tsv','inputs/reference/gap_n343.tsv'}
assert {p.relative_to(OUT).as_posix() for p in (OUT/'inputs').rglob('*') if p.is_file()}==expected
for forbidden in ['observations_n343.tsv','cox_design_axis_library_adjusted_n245.tsv','axis_scaled_Schoenfeld_library_adjusted.tsv','restricted_os_fix.tsv']:
    assert not list(OUT.rglob(forbidden)),forbidden
manifest=OUT/'MANIFEST_SHA256.tsv'
if manifest.exists():
    with manifest.open(encoding='utf-8') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            file=OUT/row['Path'];assert file.is_file() and hashlib.sha256(file.read_bytes()).hexdigest()==row['SHA256'],row['Path']
print('OK: scaling 343x46; 334/343 labels retained; sensitivity 116/149/78; six axis/library models; public inputs exclude individual outcomes')
