"""Validate the final figure typography without refitting any scientific model."""
from pathlib import Path
from collections import Counter
import argparse,csv,hashlib,json,re,struct
import xml.etree.ElementTree as E
import pandas as pd
import numpy as np
import pymupdf

Q=Path(__file__).resolve().parent;R=Q.parent
parser=argparse.ArgumentParser();parser.add_argument('--figure-dir',type=Path)
args=parser.parse_args();F=args.figure_dir or R/'figures'
NS={'s':'http://www.w3.org/2000/svg'};XL='{http://www.w3.org/1999/xlink}href'
alias={'FigureS1_survival':'FigureS_survival','FigureS2_state_structure':'FigureS_state_structure','FigureS3_immune_features':'FigureS_immune_features'}
def signature(g):
    parent={c:p for p in g.iter() for c in p};out=[]
    for n in g.iter():
        ancestors=[];p=n
        while p in parent:p=parent[p];ancestors.append(p)
        if any(a.tag.endswith('defs') or a.get('id','').startswith('legend_') for a in ancestors):continue
        k=n.tag.split('}')[-1]
        if k=='path':out.append((k,n.get('d')))
        if k=='use':out.append((k,n.get(XL),n.get('x'),n.get('y')))
        if k=='image':out.append((k,n.get(XL),n.get('x'),n.get('y'),n.get('width'),n.get('height')))
    return Counter(out)
audit=json.loads((Q/'layout_audit.json').read_text(encoding='utf-8'));assert len(audit)==14
letters=0;geometry_axes=0
for row in audit:
    name=row['figure'];stem=alias.get(name,name);svgs=list(F.rglob(stem+'.svg'))
    if not svgs:svgs=list(F.rglob(name+'.svg'))
    assert len(svgs)==1,(name,svgs)
    svg=svgs[0];source=next((Q/'inputs/baseline_v2.0.1').rglob(name+'.svg'))
    assert hashlib.sha256(source.read_bytes()).hexdigest()==row['source_sha256']
    tree=E.parse(svg);root=tree.getroot();old=E.parse(source)
    texts=root.findall('.//s:text',NS);alltext=' '.join(''.join(t.itertext()) for t in texts)
    assert not re.search(r'Tumours|Harmonised|Labels carried from n344|molecular_group_|location_Other/NA|1\.1e-09',alltext),name
    for t in texts:assert 'font-family:Arial' in t.get('style',''),(name,t.text)
    panels=root.findall('.//s:g[@data-role="panel-letter"]',NS)
    expected=[p for p in row['panels'] if p['panel']];assert len(panels)==len(expected)
    for g,p in zip(panels,expected):
        t=g.find('s:text',NS);st=t.get('style','').replace(' ','')
        assert t.text==p['panel'] and 'font-size:12.0px' in st and 'font-weight:bold' in st and 'fill:#000000' in st,(name,st)
        xy=list(map(float,re.search(r'translate\(([^)]*)\)',g.get('transform'))[1].split()))
        assert abs(xy[0]-p['x'])<1e-5 and abs(xy[1]-(p['y']+12))<1e-5
        assert p['body_origin'][0]>=p['x']+18 and p['letter_pt']==12 and p['title_pt']==10
        letters+=1
    for p in row['panels']:
        for a in p['axes']:
            s=old.find(f'.//s:g[@id="axes_{a}"]',NS);d=root.find(f'.//s:g[@id="axes_{a}"]',NS)
            assert s is not None and d is not None
            assert signature(s)==signature(d),(name,a,'scientific geometry changed')
            geometry_axes+=1
    # S4's colorbar is relocated independently, but its raster payload remains exact.
    assert Counter(i.get(XL) for i in old.findall('.//s:image',NS))==Counter(i.get(XL) for i in root.findall('.//s:image',NS)),name
    for ext in ['.png','.pdf']:assert svg.with_suffix(ext).stat().st_size>1000
    png=svg.with_suffix('.png').read_bytes();pw,ph=struct.unpack('>II',png[16:24])
    width,height=map(float,root.get('viewBox').split()[2:]);assert pw==round(width/72*400) and abs(ph-height/width*pw)<2
    with pymupdf.open(svg.with_suffix('.pdf')) as pdf:
        assert len(pdf)==1 and abs(pdf[0].rect.width-width)<.01 and abs(pdf[0].rect.height-height)<.01
    if name=='Figure5':assert 'Cytotoxic Lymphocyte Curated*' in alltext and '0.075' not in alltext
    if name=='FigureS5_standardization_consistency':assert 'Within-cohort standardization and' not in alltext
    if name=='FigureS6_axis_library_survival':assert 'Continuous immune-axis survival' not in alltext and 'All models adjusted' not in alltext
assert letters==54,letters
a=pd.read_csv(Q/'tables/Figure5A_final_display_source.tsv',sep='\t')
b=pd.read_csv(R/'reviewer_followup_20261103/tables/Figure5A_source_AUROC.tsv',sep='\t')
r=a[a.signature.eq('Cytotoxic_Lymphocyte_Curated')].iloc[0]
assert r.asterisk and not r.prespecified_lineage_concordant
assert tuple(r[['AUROC_all_cells','AUROC_min_sample','AUROC_max_sample']])==(.925,.927,.976)
assert len(a)==16 and int(a.prespecified_lineage_concordant.sum())==15
orig=b[b.AUROC_all_cells.notna()];assert len(orig)==16 and (orig.verdict=='concordant').sum()==15 and np.isclose(orig.AUROC_all_cells.median(),.868)
base=b.set_index('signature');cur=a.set_index('signature')
for name in cur.index:
    old=base.loc[name];new=cur.loc[name]
    if name=='Cytotoxic_Lymphocyte_Curated':
        assert new.AUROC_all_cells==round(1-old.AUROC_all_cells,3)
        assert new.AUROC_min_sample==round(1-old.AUROC_max_sample,3)
        assert new.AUROC_max_sample==round(1-old.AUROC_min_sample,3)
    else:assert np.allclose(new[['AUROC_all_cells','AUROC_min_sample','AUROC_max_sample']].astype(float),old[['AUROC_all_cells','AUROC_min_sample','AUROC_max_sample']].astype(float))
check=json.loads((Q/'tables/AUROC_pooling_check.json').read_text())
assert check['within_specimen_pairs']+check['cross_specimen_pairs']==check['total_positive_negative_pairs']==4699*13920
assert not check['exact_per_cell_scores_reconstructed']
lo=(.925*check['total_positive_negative_pairs']-.976*check['within_specimen_pairs'])/check['cross_specimen_pairs']
hi=(.925*check['total_positive_negative_pairs']-.927*check['within_specimen_pairs'])/check['cross_specimen_pairs']
assert 0<=lo<=hi<=1 and np.allclose([lo,hi],check['compatible_cross_specimen_AUROC_bounds'])
legend=(R/'FIGURE_LEGENDS.md').read_text(encoding='utf-8');assert '15/16' in legend and 'between-specimen' in legend and '*Cytotoxic' in legend
da=(R/'DATA_AVAILABILITY.md').read_text(encoding='utf-8');assert [u.rstrip('.') for u in re.findall(r'https://github\.com/\S+',da)]==['https://github.com/paperhwi/pediatric_glioma_immuneprofiling/releases/tag/v2.0.2']
if (Q/'MANIFEST_SHA256.tsv').exists():
    with (Q/'MANIFEST_SHA256.tsv').open(encoding='utf-8',newline='') as f:
        for v in csv.DictReader(f,delimiter='\t'):
            path=Q/v['Path'];assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==v['SHA256'],v['Path']
print(f'OK: 14 SVG/PDF/PNG sets; {letters} Arial Bold black 12-pt letters; {geometry_axes} unchanged source axes; original-lineage concordance 15/16; curated T-cell AUROC 0.925; pooled-AUROC compatibility; final manifest')
