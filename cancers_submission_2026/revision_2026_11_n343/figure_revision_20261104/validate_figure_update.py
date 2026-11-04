"""Read-only validation of the stored November 4 correction bundle."""
from pathlib import Path
import csv,hashlib,re,json,xml.etree.ElementTree as E
import pandas as pd,numpy as np
W=Path(__file__).resolve().parent;R=W.parent
with (W/'MANIFEST_SHA256.tsv').open(encoding='utf-8',newline='') as f:
    for row in csv.DictReader(f,delimiter='\t'):
        p=W/row['Path'];assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==row['SHA256'],row['Path']
a=pd.read_csv(W/'tables/Figure5A_AUROC_current_targets.tsv',sep='\t')
b=pd.read_csv(R/'reviewer_followup_20261103/tables/Figure5A_source_AUROC.tsv',sep='\t')
sel=a.signature.eq('Cytotoxic_Lymphocyte_Curated');old=b[b.signature.str.contains('Cytotoxic_Lymphocyte_Curated|Glioma_Inflammatory_Wang2017')].iloc[0];r=a[sel].iloc[0]
assert r.AUROC_all_cells==round(1-old.AUROC_all_cells,3)==.925
assert r.AUROC_min_sample==round(1-old.AUROC_max_sample,3)==.927
assert r.AUROC_max_sample==round(1-old.AUROC_min_sample,3)==.976
cols=['AUROC_all_cells','AUROC_min_sample','AUROC_max_sample']
assert np.allclose(a[~sel].set_index('signature')[cols].sort_index(),b[~b.signature.str.contains('Cytotoxic_Lymphocyte_Curated|Glioma_Inflammatory_Wang2017')].set_index('signature')[cols].sort_index(),equal_nan=True)
directional=a[a.AUROC_all_cells.notna()];assert len(directional)==16 and (directional.verdict=='concordant').sum()==16 and np.isclose(directional.AUROC_all_cells.median(),.879)
assert json.loads((W/'tables/Figure5A_target_audit.json').read_text())['positive_cells']==4699
ns={'s':'http://www.w3.org/2000/svg'}
svgs=list((R/'figures').glob('*.svg'));assert len(svgs)==14,len(svgs)
for f in svgs:
    assert f.with_suffix('.png').is_file() and f.with_suffix('.pdf').is_file(),f.name
    svg=E.parse(f);text=' '.join(''.join(t.itertext()) for t in svg.findall('.//s:text',ns))
    assert not re.search(r'Tumours|Harmonised|Labels carried from n344|molecular_group_|location_Other/NA|1\.1e-09',text),f.name
    legend_nodes={t for g in svg.findall('.//s:g',ns) if g.get('id','').startswith('legend_') for t in g.findall('.//s:text',ns)}
    width=float(svg.getroot().get('viewBox').split()[2])
    for t in svg.findall('.//s:text',ns):
        if t not in legend_nodes and re.fullmatch('[A-H]',''.join(t.itertext()).strip()):
            style=t.get('style','');size=float(re.search(r'font-size:\s*([\d.]+)',style)[1])
            assert 'font-weight:700' in style.replace(' ','') and abs(size/width*504-14)<.001,(f.name,style)
assert '0.075' not in (R/'FIGURE_LEGENDS.md').read_text(encoding='utf-8')
print('OK: current curated T-cell-target AUROC 0.925 (specimen range 0.927–0.976); 16/16 current lineage concordance; median 0.879; 14 SVG/PDF/PNG figure sets; correction manifest verified')
