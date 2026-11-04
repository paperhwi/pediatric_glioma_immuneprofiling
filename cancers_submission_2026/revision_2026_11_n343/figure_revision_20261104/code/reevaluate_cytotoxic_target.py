"""Re-evaluate the curated target using the exact binary-complement identity.

Original cell scores are unchanged. The retained historical AUROCs have only
three decimals; the new output intentionally claims no additional precision.
"""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd,h5py

parser=argparse.ArgumentParser()
parser.add_argument('--annotation-h5ad',type=Path,required=True)
parser.add_argument('--output-dir',type=Path,required=True)
args=parser.parse_args()
root=Path(__file__).resolve().parents[2]
with h5py.File(args.annotation_h5ad) as f:
    g=f['obs/detailed_annot']
    categories=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in g['categories'][()]])
    codes=g['codes'][()]
    assert (codes>=0).all()
    labels=categories[codes]
    counts={str(k):int(v) for k,v in zip(*np.unique(labels,return_counts=True))}
assert counts=={'CD4':2385,'CD8':2314,'Myeloid':13920,'unclear':93},counts
old=pd.read_csv(root/'reviewer_followup_20261103/tables/Figure5A_source_AUROC.tsv',sep='\t')
selector=old.signature.isin(['Cytotoxic_Lymphocyte_Curated','Glioma_Inflammatory_Wang2017'])
assert selector.sum()==1
r=old.loc[selector].iloc[0]
values={'signature':'Cytotoxic_Lymphocyte_Curated','positive_class':'annotated CD4 or CD8 T cells','negative_class':'Myeloid','AUROC_all_cells':round(1-r.AUROC_all_cells,3),'AUROC_min_sample':round(1-r.AUROC_max_sample,3),'AUROC_max_sample':round(1-r.AUROC_min_sample,3),'NK_separately_annotated':False,'source_precision_decimals':3}
assert values['AUROC_all_cells']==.925
args.output_dir.mkdir(parents=True,exist_ok=True)
pd.DataFrame([values]).to_csv(args.output_dir/'curated_target_AUROC.tsv',sep='\t',index=False)
(args.output_dir/'annotation_counts.json').write_text(json.dumps(counts,indent=2)+'\n',encoding='utf-8')
print(json.dumps(values,indent=2))
