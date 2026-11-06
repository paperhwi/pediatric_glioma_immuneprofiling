"""Lay out immutable SVG plot components at their final physical size.

Geometry and text positions are transformed separately, so Arial glyphs remain
undistorted when a dense panel needs more vertical space. No model is refitted.
"""
from pathlib import Path
import argparse,copy,json,re,math,functools,hashlib,os
import numpy as np
from lxml import etree as E
from matplotlib.font_manager import FontProperties,findfont
from matplotlib.textpath import TextPath,TextToPath
from matplotlib.path import Path as MP
from matplotlib.transforms import Affine2D
import cairo,cairosvg

W=Path(__file__).parent
parser=argparse.ArgumentParser();parser.add_argument('--baseline',type=Path);parser.add_argument('--output',type=Path);parser.add_argument('--only',nargs='*');parser.add_argument('--no-audit',action='store_true');args=parser.parse_args()
if args.baseline is None:
 P=json.loads((W/'paths.json').read_text());BASE=Path(P['previous'])/'03_figures';OUT=Path(P['delivery'])/'01_figures'
else:BASE=args.baseline;OUT=args.output
assert OUT is not None
findfont('Arial',fallback_to_default=False)
NS='http://www.w3.org/2000/svg';tag=lambda x:'{'+NS+'}'+x
NUM=r'[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?'
TTP=TextToPath();PAGE=178/25.4*72;MARGIN=8;GUTTER=18;BODY=8.;TITLE=10.;LETTER=12.
SUP='⁰¹²³⁴⁵⁶⁷⁸⁹⁻';PLAIN='0123456789−';UNSUP=str.maketrans(SUP,PLAIN)
SUB='₀₁₂₃₄₅₆₇₈₉₋';UNSUB=str.maketrans(SUB,PLAIN)
audit=[]
def write_output(path,payload):
 temp=path.with_name(path.name+'.codex.tmp');temp.write_bytes(payload);os.replace(temp,path)

def style(s):return {k.strip():v.strip() for k,v in (x.split(':',1) for x in s.split(';') if ':' in x)}
def fsize(s,default=8):return float(re.sub('[^0-9.]','',s.get('font-size',str(default))))
def prop(size,weight='normal',italic='normal'):return FontProperties(family='Arial',size=size,weight=weight,style=italic)
@functools.lru_cache(maxsize=12000)
def metrics(text,size,weight='normal',italic='normal'):
 text=text.replace('⁻','−').translate(UNSUB)  # Conservative width preflight only.
 if not text.strip():return (0.,0.,0.,0.),TTP.get_text_width_height_descent(text,prop(size,weight,italic),False)[0]
 p=TextPath((0,0),text,prop=prop(size,weight,italic),size=size)
 b=p.get_extents();adv=TTP.get_text_width_height_descent(text,prop(size,weight,italic),False)[0]
 return (float(b.x0),float(-b.y1),float(b.x1),float(-b.y0)),adv
def matrix(s):
 out=np.eye(3)
 for op,arg in re.findall(r'(translate|scale|rotate|matrix)\(([^)]*)\)',s):
  v=list(map(float,re.findall(NUM,arg)));q=np.eye(3)
  if op=='translate':q[:2,2]=[v[0],v[1] if len(v)>1 else 0]
  elif op=='scale':q[0,0]=v[0];q[1,1]=v[1] if len(v)>1 else v[0]
  elif op=='rotate':
   a=math.radians(v[0]);r=np.array([[math.cos(a),-math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1.]])
   if len(v)>1:
    t=np.eye(3);t[:2,2]=v[1:3];ti=np.eye(3);ti[:2,2]=-np.array(v[1:3]);q=t@r@ti
   else:q=r
  else:q=np.array([[v[0],v[2],v[4]],[v[1],v[3],v[5]],[0,0,1.]])
  out=out@q
 return out
def fullmatrix(n):
 out=np.eye(3)
 for q in list(n.iterancestors())[::-1]+[n]:out=out@matrix(q.get('transform',''))
 return out
def groupid(n):
 return next((a.get('id') for a in [n]+list(n.iterancestors()) if a.get('id','').startswith('text_')),'')
def inherited(n):
 st={}
 for a in list(n.iterancestors())[::-1]+[n]:st.update(style(a.get('style','')))
 return st
def union(boxes):
 return np.array([min(b[0] for b in boxes),min(b[1] for b in boxes),max(b[2] for b in boxes),max(b[3] for b in boxes)])
def pathbbox(d):
 v=list(map(float,re.findall(NUM,d)));return [min(v[::2]),min(v[1::2]),max(v[::2]),max(v[1::2])]
def axisbox(g):
 boxes=[]
 for ch in g:
  if ch.get('id','').startswith('patch_'):
   p=ch.find('.//'+tag('path'))
   if p is not None:
    box=pathbbox(p.get('d'))
    # Use the plotting rectangle rather than bars extending beyond the
    # clipped data area (or an undefined bar represented at SVG origin).
    if style(p.get('style','')).get('fill') in ['#ffffff','white']:return box
    boxes.append(box)
 if boxes:return union(boxes)
 pts=[(float(n.get('x')),float(n.get('y'))) for n in g.iter() if n.tag==tag('use') and n.get('x') is not None and n.get('y') is not None]
 if pts:
  a=np.array(pts);return [a[:,0].min()-2,a[:,1].min()-2,a[:,0].max()+2,a[:,1].max()+2]
 return None
def richpieces(text,size,weight,italic,anchor='start',x=0,y=0):
 if '\n' in text:
  return [p for i,line in enumerate(text.splitlines()) for p in richpieces(line,size,weight,italic,anchor,x,y+i*size*1.18)]
 parts=[];cursor=0
 for token in re.findall('['+SUP+']+|['+SUB+']+|[^'+SUP+SUB+']+',text):
  supers=all(c in SUP for c in token);subs=all(c in SUB for c in token);s=size*.72 if supers or subs else size;t=token.translate(UNSUP) if supers else token.translate(UNSUB) if subs else token
  if not t.strip():cursor+=metrics(t,s,weight,italic)[1];continue
  parts.append({'text':t,'x':x+cursor,'y':y-size*.40 if supers else y+size*.25 if subs else y,'size':s,'weight':weight,'italic':italic,'anchor':'start','exponent':supers or subs})
  cursor+=metrics(t,s,weight,italic)[1]
 offset=cursor/2 if anchor=='middle' else cursor if anchor=='end' else 0
 for p in parts:p['x']-=offset
 return parts
def textmodel(n,desired=BODY,replace=None):
 st=inherited(n);ownstyle=style(n.get('style',''));x=float(n.get('x','0'));y=float(n.get('y','0'))
 sizes=[fsize(ownstyle)] if 'font-size' in ownstyle else []
 sizes += [fsize(style(c.get('style',''))) for c in n.iter(tag('tspan')) if 'font-size' in style(c.get('style',''))]
 base=max(sizes or [8]);factor=desired/base;m=fullmatrix(n);anchorpos=m@np.array([x,y,1]);angle=math.degrees(math.atan2(m[1,0],m[0,0]));parts=[]
 if replace is not None:
  for i,line in enumerate(replace.split('\n')):parts+=richpieces(line,desired,st.get('font-weight','normal'),st.get('font-style','normal'),st.get('text-anchor','start'),0,i*desired*1.18)
 elif len(n)==0:
  parts=richpieces(n.text or '',desired,st.get('font-weight','normal'),st.get('font-style','normal'),st.get('text-anchor','start'))
 else:
  if n.text and n.text.strip():parts+=richpieces(n.text,desired,st.get('font-weight','normal'),st.get('font-style','normal'),st.get('text-anchor','start'))
  for c in n.iter(tag('tspan')):
   if not (c.text or '').strip():continue
   cs=st|style(c.get('style',''));size=fsize(cs,base)*factor;xs=list(map(float,c.get('x',str(x)).split()));ys=list(map(float,c.get('y',str(y)).split()));txt=''.join(c.itertext())
   if len(xs)>1:
    for i,ch in enumerate(txt):
     xx=(xs[min(i,len(xs)-1)]-x)*factor;yy=(ys[min(i,len(ys)-1)]-y)*factor
     z=richpieces(ch,size,cs.get('font-weight','normal'),cs.get('font-style','normal'),'start',xx,yy);parts+=z
   else:parts+=richpieces(txt,size,cs.get('font-weight','normal'),cs.get('font-style','normal'),cs.get('text-anchor','start'),(xs[0]-x)*factor,(ys[0]-y)*factor)
 for p in parts:p['fill']=st.get('fill','#000000')
 # Every glyph is an undistorted Arial glyph, including superscript exponents.
 isy=any(a.get('id','').startswith('ytick_') for a in n.iterancestors());isx=any(a.get('id','').startswith('xtick_') for a in n.iterancestors());hasaxis=any(a.get('id','').startswith('matplotlib.axis_') for a in n.iterancestors());kind='ylabel' if hasaxis and not isy and not isx and abs(abs(angle)-90)<1 else 'xlabel' if hasaxis and not isy and not isx and abs(angle)<1 else ''
 return {'source':''.join(n.itertext()),'anchor':anchorpos[:2],'angle':angle,'parts':parts,'group':groupid(n),'is_ytick':isy,'is_xtick':isx,'axis_kind':kind,'native_size':base,'tick_id':next((a.get('id') for a in n.iterancestors() if a.get('id','').startswith('ytick_')),''),'axis':next((a.get('id') for a in n.iterancestors() if a.get('id','').startswith('axes_')),'')}
def scaledanchor(t,sx,sy):
 if 'flow_origin' in t:a=t['flow_origin']*[sx,sy]+t['flow_offset']
 else:a=np.array([t['anchor'][0]*sx,t.get('rowcenter',t['anchor'][1])*sy+t.get('rowoffset',0)])
 return a+np.array([t.get('extra_x',0)+t.get('auto_x',0),t.get('extra_y',0)+t.get('auto_y',0)])
def modelbox(t,sx,sy):
 boxes=[];ang=math.radians(t['angle']);rot=np.array([[math.cos(ang),-math.sin(ang)],[math.sin(ang),math.cos(ang)]])
 for p in t['parts']:
  b,adv=metrics(p['text'],p['size'],p['weight'],p['italic']);anchor=p['anchor'];shift=adv/2 if anchor=='middle' else adv if anchor=='end' else 0
  corners=np.array([[b[0],b[1]],[b[2],b[1]],[b[2],b[3]],[b[0],b[3]]])+[p['x']-shift,p['y']];pts=corners@rot.T+scaledanchor(t,sx,sy)
  if np.isfinite(pts).all():boxes.append([pts[:,0].min(),pts[:,1].min(),pts[:,0].max(),pts[:,1].max()])
 return union(boxes) if boxes else [0,0,0,0]
def addmodel(root,t,x=0,y=0,sx=1,sy=1,role='body'):
 gx,gy=scaledanchor(t,sx,sy)+[x,y];g=E.SubElement(root,tag('g'),transform=f'translate({gx:.6f} {gy:.6f}) rotate({t["angle"]:.6f})',**{'data-role':role})
 for p in t['parts']:
  n=E.SubElement(g,tag('text'),x=str(p['x']),y=str(p['y']),style=f'font-family:Arial; font-size:{p["size"]}px; font-weight:{p["weight"]}; font-style:{p["italic"]}; text-anchor:{p["anchor"]}; fill:{p["fill"]}')
  if p['exponent']:n.set('data-role','math-exponent')
  n.text=p['text']
 return g
def simple(root,text,x,y,size=8,weight='normal',anchor='start',role='body',fill='#000000'):
 t={'anchor':np.array([x,y]),'angle':0,'parts':richpieces(text,size,weight,'normal',anchor)}
 for p in t['parts']:p['fill']=fill
 return addmodel(root,t,role=role)
def wrap(text,width,size=10):
 lines=[]
 for para in text.split('\n'):
  cur=''
  for word in para.split():
   q=(cur+' '+word).strip()
   if cur and metrics(q,size)[1]>width:lines.append(cur);cur=word
   else:cur=q
  lines.append(cur)
 return lines
def legenditems(g):
 items=[];pending=[];title=''
 for ch in g:
  ts=list(ch.iter(tag('text')))
  if not ts:
   fills=[style(p.get('style','')).get('fill') for p in ch.iter(tag('path'))]
   if ch.get('id','').startswith('patch_') and any(c in ['#ffffff','white'] for c in fills) and not pending and not items:continue
   pending.append(ch);continue
  label=' '.join((''.join(c.text or '' for c in t.iter(tag('tspan'))) if len(t) else (t.text or '')) for t in ts).replace('\xa0',' ')
  label=' '.join(label.split())
  if label.startswith('identity'):label='Identity (y = x)'
  if not pending and not items:title=label;continue
  colors=[];kind='line'
  for q in pending:
   if q.get('id','').startswith('patch_'):kind='patch'
   for e in q.iter():
    for k,v in style(e.get('style','')).items():
     if k in ['stroke','fill'] and v not in ['none','#ffffff','white','#000000','black']:colors.append(v)
  items.append({'label':label,'color':colors[0] if colors else '#000000','kind':kind});pending=[]
 return title,items
def legendheight(items,width,title=''):
 rows=1;xx=0
 for item in items:
  w=min(width,metrics(item['label'],8)[1]+22)
  if xx and xx+w>width:rows+=1;xx=0
  xx+=w+12
 return (11 if title else 0)+rows*13+3
def drawlegend(root,items,x,y,width,title=''):
 if title:simple(root,title,x,y+8,8);y+=12
 xx=0;yy=y+8
 for it in items:
  w=metrics(it['label'],8)[1]+22
  if xx and xx+w>width:xx=0;yy+=13
  if it['kind']=='patch':E.SubElement(root,tag('rect'),x=str(x+xx),y=str(yy-5.5),width='10',height='5.5',fill=it['color'])
  else:E.SubElement(root,tag('line'),x1=str(x+xx),x2=str(x+xx+11),y1=str(yy-3),y2=str(yy-3),stroke=it['color'],**{'stroke-width':'1.2'})
  simple(root,it['label'],x+xx+16,yy,8);xx+=w+12

# Each row defines the physical panel order; all scientific values come from
# the supplied SVGs. New headers merely relocate existing descriptions.
def spec(letter,axes,title='',sub='',extra=None,drop=None,keep=None,minsy=0,span=1,legend='local',fillheight=False):return dict(letter=letter,axes=axes,title=title,sub=sub,extra=extra or [],drop=drop or [],keep=keep or {},minsy=minsy,span=span,legend=legend,fillheight=fillheight)
LAYOUTS={
 'Figure1':[[spec('A',[1],'Study workflow',drop=['Study workflow'],minsy=.9)],[spec('B',[2],'Analysis-set derivation',drop=['Analysis-set derivation'],minsy=.9)]],
 'Figure2':[[spec('A',[1,2,3],'LM22 deconvolution quality',keep={'P-value':'P-value','Correlation':'Correlation','RMSE':'RMSE'})],[spec('B',[4],'Cross-method concordance (n = 343)',drop=['Cross-method concordance'])],[spec('C',[5,6,7,8],'Immune estimates by method and state',keep={s:s for s in ['CD8 T cell','NK cell','Macrophage','Monocyte / DC']})]],
 'Figure3':[[spec('A',[1],'Ambiguous clustering',drop=['Ambiguous clustering']),spec('B',[2],'Silhouette width',drop=['Silhouette width']),spec('C',[3],'Consensus CDF',drop=['Consensus CDF'])],[spec('D',[4],'Gap statistic',drop=['Gap statistic']),spec('E',[5],'Minimum cluster size',drop=['Minimum cluster size']),None],[spec('F',[6,7],'PCA and UMAP projections',keep={'PCA of the 46-feature matrix':'PCA','UMAP of the 46-feature matrix':'UMAP'})],[spec('G',[8,9],'Top 20 discriminating input features',drop=['Top 20 discriminating input features'])]],
 'Figure4':[[spec('A',[1],'Immune state by anatomical location',drop=['Immune state by anatomical location']),spec('B',[2],'Immune state by molecular group',drop=['Immune state by molecular group'])],[spec('C',[3],'Location within H3/IDH-wildtype pHGG',drop=['Within H3/IDH-wt pHGG']),spec('D',[4],'Sequential PERMANOVA','34-feature TPM matrix; n = 253',drop=['Sequential PERMANOVA'])]],
 'Figure5':[[spec('A',[1],'Directional signature discrimination',drop=['Directional signature discrimination'])],[spec('B',[2,3],'Cell class with highest mean expression',drop=['Cell class with highest mean expression'])],[spec('C',[4,5,6],'Spatial nearest-centroid assignments',keep={f'Sample-{i} (n = {n})':f'Sample-{i} (n = {n})' for i,n in [(1,577),(2,856),(3,217)]},extra=['legend_3'])],[spec('D',[7],'Spot immune fraction by centroid',drop=['Spot immune fraction by centroid']),spec('E',[8],'Spatial axis and cell mixtures',drop=['Spatial axis and cell mixtures'])]],
 'Figure6':[[spec('A',[1],'Immune-high vs Immune-low','Reference: Immune-low; up 3,718; down 228',drop=['Immune-high vs Immune-low']),spec('B',[2],'Immune-intermediate vs Immune-low','Reference: Immune-low; up 431; down 69',drop=['Immune-intermediate vs Immune-low']),spec('C',[3],'Immune-high vs Immune-intermediate','Reference: Immune-intermediate; up 808; down 117',drop=['Immune-high vs Immune-intermediate'])],[spec('D',[4],'Overlap of up-regulated genes','404/431 (93.7%) of the Immune-intermediate set is contained in the Immune-high set',drop=['404/431'])],[spec('E',[5],'Adjusted effect sizes across contrasts','All 18,564 genes; Pearson r = 0.78',drop=['Adjusted effect sizes across contrasts','All'],keep={}),spec('F',[6],'Shared-gene effect sizes','Paired Wilcoxon signed-rank P = 1.6 × 10⁻⁶⁵',drop=['Shared-gene effect sizes','paired'])]],
 'Figure7':[[spec('A',[1,2],'Overall survival','All tumors with survival data; full-follow-up log-rank P = 0.058',drop=['All tumors with overall-survival data'])],[spec('B',[3],'Primary Cox model','n = 245; 206 deaths',drop=['Primary Cox model'])]],
 'FigureS1_survival':[[spec('A',[1],'Survival sensitivity analyses',drop=['Survival sensitivity analyses'])],[spec('B',[2,3],'DMG, H3 K27-altered','Full-follow-up log-rank P = 0.134',drop=['DMG, H3 K27-altered']),spec('C',[4,5],'pHGG, H3/IDH-wildtype','Full-follow-up log-rank P = 0.096',drop=['pHGG, H3/IDH-wildtype'])],[spec('D',[6],'Within-subtype Cox models','State × subtype interaction P = 0.780',drop=['Within-subtype Cox models']),spec('E',[7],'Incremental prognostic information',drop=['Incremental prognostic information'])]],
 'FigureS2_state_structure':[[spec('A',[1],'Seed partitions','Modal counts 122/144/77; median ARI 0.94',drop=['Seed partitions'],legend='global'),spec('B',[2],'k = 3 states within k = 2 clusters','LRT: P = 0.130',drop=['k = 3 states nested'])],[spec('C',[3],'Re-clustering within subsets (k = 3)',drop=['Re-clustering within subsets']),spec('D',[4],'Reference-based deconvolution','JT P = 1.2 × 10⁻⁴⁷; ρ with immune axis = 0.78',drop=['Reference-based deconvolution'])],[spec('E',[5],'Bagaev MFP classes by immune state',drop=['Bagaev MFP classes']),spec('F',[6],'Total immune estimate, five methods','All JT P < 1.0 × 10⁻¹⁵',drop=['Total immune estimate'],legend='global')]],
 'FigureS3_immune_features':[[spec('A',[1],'Immune features and curated programs','BH q < 0.05 in either model; n = 245, 206 deaths',drop=['Immune features and curated programs'])],[spec('B',[2],'Published programs',drop=['Published programs'])],[spec('C',[3],'Deconvolution-estimated tumor-cell programs by molecular group',drop=['Deconvolution-estimated'])]],
 'FigureS4_ScPCA_axis_linkage':[[spec('A',[1],'Harmonized bulk scale','Violin: OpenPedCan (n = 343); dots: ScPCA donors (n = 17)',drop=['Harmonized bulk scale']),spec('B',[2],'Immune content','ρ = 0.78; permutation P = 3.0 × 10⁻⁴\nOrdered-state trend P = 0.003',drop=['Immune content'],legend='global'),spec('C',[3],'CD8 share of T cells','ρ = −0.44; permutation P = 0.073\nOrdered-state trend P = 0.032',drop=['CD8 share of T cells'])],[spec('D',[4],'Donor-level function and bulk immune axis','CD8 n = 11; CD4 n = 8; myeloid n = 14 donors (≥10 cells each)\nFilled: permutation P < 0.05 (none BH q < 0.05); dotted: |ρ| ≈ 0.6',drop=['Single-cell function vs bulk immune axis'],span=2),spec('E',[5,6],'Programs by projected state','Median donor z\nImmune-low donors: CD8 1; CD4 0; myeloid 2',drop=['By projected state'],fillheight=True)],[spec('F',[7],'CD8: PRF1-positive cells','n = 11; ρ = 0.62\nPermutation P = 0.050; q = 0.41',drop=['CD8: % PRF1+']),spec('G',[8],'Myeloid: suppressive TAM','n = 14; ρ = 0.59\nPermutation P = 0.032; q = 0.41',drop=['Myeloid: Immunosuppressive TAM']),spec('H',[9],'Myeloid: MDM identity','n = 14; ρ = 0.54\nPermutation P = 0.047; q = 0.41',drop=['Myeloid: MDM identity'])]],
 'FigureS5_standardization_consistency':[[spec('',[1],extra=['text_19'])]],
 'FigureS6_axis_library_survival':[[spec('',[1],extra=['text_'+str(i) for i in range(15,24)])]],
 'FigureS7_IPW_sensitivity':[[spec('A',[1],'Outcome availability',drop=['Outcome availability']),spec('B',[2],'Covariate balance',drop=['Covariate balance'])],[spec('C',[3],'Weighted Cox sensitivity',drop=['Weighted Cox sensitivity'])]]
}

CAPS={
 'Figure3':{'A':140,'B':140,'C':140,'D':140,'E':140,'G':200},
 'Figure5':{'A':192,'B':120,'D':140,'E':140},
 'Figure6':{'A':110,'B':110,'C':110,'D':65,'E':145,'F':145},
 'Figure7':{'A':180,'B':160},
 'FigureS1_survival':{'A':230,'B':140,'C':140,'D':110,'E':110},
 'FigureS2_state_structure':{p:140 for p in 'ABCDEF'},
 'FigureS3_immune_features':{'A':320,'B':160,'C':230},
 'FigureS7_IPW_sensitivity':{'A':150,'B':150,'C':130}}
# Long two-line subtype contrasts and their estimates need a full-width row.
LAYOUTS['FigureS1_survival'][2:]=[[LAYOUTS['FigureS1_survival'][2][0]],[LAYOUTS['FigureS1_survival'][2][1]]]
# Give the long subset names, percentage/ARI annotations and five ticks a
# full-width row instead of squeezing the actual bars into the left column.
LAYOUTS['FigureS2_state_structure'][1:2]=[[LAYOUTS['FigureS2_state_structure'][1][0]],[LAYOUTS['FigureS2_state_structure'][1][1]]]
for name,rows in LAYOUTS.items():
 for row in rows:
  for sp in row:
   if sp is None:continue
   sp['cap']=CAPS.get(name,{}).get(sp['letter'])
   sp['minrow']=20 if name=='FigureS3_immune_features' and sp['letter']=='A' else 25 if name=='FigureS3_immune_features' and sp['letter']=='B' else 30 if name=='FigureS3_immune_features' and sp['letter']=='C' else 10
   if name=='FigureS1_survival' and sp['letter']=='D':sp['minrow']=25
   if name=='FigureS4_ScPCA_axis_linkage' and sp['letter']=='D':sp['span']=1.7
   if name=='FigureS4_ScPCA_axis_linkage' and sp['letter']=='E':sp['span']=1.3;sp['axes']=[5];sp['colorbar']=6

def prepare_component(tree,sp,name):
 root=tree.getroot();allids={g.get('id'):g for g in root.iter() if g.get('id')};groups=[allids['axes_'+str(i)] for i in sp['axes']]+[allids[k] for k in sp['extra']]
 alltexts=[t for g in groups for t in g.iter(tag('text'))];dropids=set();legends=[]
 for g in groups:
  if g.get('id','').startswith('legend_'):legends.append(g)
  legends.extend(g.xpath('.//s:g[starts-with(@id,"legend_")]',namespaces={'s':NS}))
 legends=list({g.get('id'):g for g in legends}.values())
 for t in alltexts:
  txt=''.join(t.itertext()).strip();compact=''.join(txt.split())
  if any(txt.startswith(p) or compact.startswith(p.replace(' ','')) for p in sp['drop']):dropids.add(groupid(t))
 models=[];grouptexts={}
 for t in alltexts:grouptexts.setdefault(groupid(t),[]).append(t)
 for t in alltexts:
  txt=''.join(t.itertext()).strip()
  if not txt or groupid(t) in dropids or re.fullmatch('[A-H]',txt):continue
  if any(a.get('id','').startswith('legend_') for a in [t]+list(t.iterancestors())):continue
  target=BODY;rep=None
  for old,new in sp['keep'].items():
   if txt==old:target=9;rep=new
  if name=='Figure5' and txt=='Cytotoxic Lymphocyte Curated':rep=txt+'*'
  if name=='Figure4':rep=re.sub(r'\bn\s*=\s*(\d+)',r'n = \1',txt.replace('P_MC','Monte Carlo P'))
  if name=='Figure4' and sp['letter']=='B' and any(a.get('id','').startswith('xtick_') for a in t.iterancestors()):
   if t is not grouptexts[groupid(t)][0]:continue
   first=''.join(grouptexts[groupid(t)][0].itertext()).strip()
   rep={'DMG':'DMG\nH3 K27-alt\n(n = 175)','DHG':'DHG\nH3 G34-mut\n(n = 31)','pHGG':'pHGG\nH3/IDH-wt\n(n = 127)','IHG':'IHG\n(n = 10)'}[first]
   models.append(textmodel(t,target,rep));continue
  if name=='Figure4' and sp['letter']=='A' and ''.join(''.join(q.itertext()).strip() for q in grouptexts[groupid(t)]).startswith('Multi-compartment/'):
   if t is not grouptexts[groupid(t)][0]:continue
   rep='Multiple or\nambiguous\n(n = 73)'
  if name=='Figure5' and sp['letter']=='B' and ''.join(''.join(q.itertext()).strip() for q in grouptexts[groupid(t)])=='Macrophage/monocyte':
   if t is not grouptexts[groupid(t)][0]:continue
   m=textmodel(t,target,'Macrophage/monocyte');ax=next(a for a in t.iterancestors() if a.get('id')=='axes_2');m['anchor'][0]=axisbox(ax)[0];m['extra_x']=-5;m['parts']=richpieces('Macrophage/monocyte',BODY,'normal','normal','end')
   m['anchor'][1]=np.mean([textmodel(q,target)['anchor'][1] for q in grouptexts[groupid(t)]])
   for piece in m['parts']:piece['fill']='#000000'
   models.append(m);continue
  if name=='Figure1':
   rep={'text_33':'Repeat or ineligible specimens','text_34':'6 exclusions: 5 non-HGG; 1 PXA','text_53':'Delayed entry: n = 187','text_54':'Non-autopsy: n = 189'}.get(groupid(t),rep)
   if sp['letter']=='B' and groupid(t) in ['text_47','text_48']:
    if groupid(t)=='text_48':continue
    m=textmodel(t,target,'− 98\nno usable\nOS record');m['anchor']=np.array([922.2,565.9461]);m['parts']=richpieces('− 98\nno usable\nOS record',BODY,'normal','normal','middle')
    mid=(min(q['y']+metrics(q['text'],q['size'])[0][1] for q in m['parts'])+max(q['y']+metrics(q['text'],q['size'])[0][3] for q in m['parts']))/2
    for q in m['parts']:q['y']-=mid;q['fill']='#000000'
    models.append(m);continue
  if name=='FigureS5_standardization_consistency' and re.fullmatch(r'\d+',txt):target=10
  if name=='FigureS6_axis_library_survival' and groupid(t) in ['text_15','text_16','text_17']:target=9
  if name=='FigureS7_IPW_sensitivity' and sp['letter']=='A' and any(a.get('id','').startswith('xtick_') for a in t.iterancestors()):
   if t is not grouptexts[groupid(t)][0]:continue
   labels={'text_1':'DHG\nH3 G34-mut','text_2':'DMG\nH3 K27-alt','text_3':'IHG','text_4':'pHGG\nH3/IDH-wt'}
   m=textmodel(t,target,labels[groupid(t)]);m['angle']=0
   tick=next(a for a in t.iterancestors() if a.get('id','').startswith('xtick_'));u=next(q for q in tick.iter(tag('use')) if q.get('x'))
   m['anchor'][0]=float(u.get('x'));m['parts']=richpieces(labels[groupid(t)],BODY,'normal','normal','middle')
   for piece in m['parts']:piece['fill']='#000000'
   models.append(m);continue
  if name=='Figure6':
   aid=next((a.get('id') for a in t.iterancestors() if a.get('id','').startswith('axes_')),'');axis=next((a.get('id') for a in t.iterancestors() if a.get('id','').startswith('matplotlib.axis_')),'');tick=any(a.get('id','').startswith(('xtick_','ytick_')) for a in t.iterancestors())
   if axis and not tick and len(''.join(txt.split()))>12:
    if t is not grouptexts[groupid(t)][0]:continue
    horizontal=axis==next(a for a in t.iterancestors() if a.get('id')==aid).xpath('./s:g[starts-with(@id,"matplotlib.axis_")]',namespaces={'s':NS})[0].get('id')
    if aid in ['axes_1','axes_2','axes_3']:rep='Adjusted log₂ fold change' if horizontal else '−log₁₀(BH q)'
    if aid=='axes_4' and horizontal:rep='Up-regulated genes (BH q < 0.05; adjusted log₂FC > 1)'
    if aid=='axes_5':rep='Adjusted log₂FC: intermediate vs low' if horizontal else 'Adjusted log₂FC: high vs low'
    if aid=='axes_6' and not horizontal:rep='Adjusted log₂FC (404 shared genes)'
   if aid=='axes_3' and tick and any(a.get('id','').startswith('xtick_') for a in t.iterancestors()):
    try:
     if not float(txt.replace('−','-')).is_integer():continue
    except ValueError:pass
  if name=='FigureS4_ScPCA_axis_linkage':
   ax=next((a.get('id') for a in t.iterancestors() if a.get('id','').startswith('axes_')),'')
   isx=any(a.get('id','').startswith('xtick_') for a in t.iterancestors())
   if ax=='axes_4':
    if txt in ['CD8 T cells','CD4 T cells','Myeloid cells']:continue
    if txt=='Donor-level Spearman ρ with bulk immune axis':rep='Spearman ρ with bulk immune axis'
    if isx:
     try:
      if not any(abs(float(txt.replace('−','-'))-v)<1e-6 for v in [-.8,-.4,0,.4,.8]):continue
     except ValueError:pass
   if ax=='axes_6' and isx and txt not in ['−1.5','0.0','1.5']:continue
   if ax=='axes_5' and isx and ''.join(txt.split()) in ['Immune-low','Immune-intermediate','Immune-high']:
    rep=''.join(txt.split());m=textmodel(t,target,rep);m['angle']=-45;tick=next(a for a in t.iterancestors() if a.get('id','').startswith('xtick_'));u=next(q for q in tick.iter(tag('use')) if q.get('x'));axg=next(a for a in t.iterancestors() if a.get('id')=='axes_5');m['anchor']=np.array([float(u.get('x')),axisbox(axg)[3]]);m['extra_y']=12;m['parts']=richpieces(rep,BODY,'normal','normal','end')
    for piece in m['parts']:piece['fill']='#000000'
    models.append(m);continue
  m=textmodel(t,target,rep)
  if name=='FigureS2_state_structure' and sp['letter']=='C' and re.fullmatch(r'\d+% \(ARI [\d.]+\)',txt):m['extra_x']=5
  models.append(m)
 # Tick labels are anchored to tick positions, not to text widths calculated
 # with a previous font. This also handles multi-line forest labels.
 for t in models:
  if not t['axis']:continue
  axg=allids[t['axis']];ab=axisbox(axg)
  if ab is None:continue
  t['axis_box']=np.asarray(ab)
  if t['is_ytick']:t['right_ticks']=t['anchor'][0]>ab[2]
 # Carry source tick centres directly from the matching XML tick group.
 for n in alltexts:
  tick=next((a for a in n.iterancestors() if a.get('id','').startswith(('xtick_','ytick_'))),None)
  if tick is None:continue
  us=[u for u in tick.iter(tag('use')) if u.get('x') is not None and u.get('y') is not None]
  if not us:continue
  u=us[0];xy=fullmatrix(u)@np.array([float(u.get('x')),float(u.get('y')),1])
  for t in models:
   if t['group']==groupid(n):t['tick_center']=xy[:2]
 if name=='FigureS3_immune_features' and sp['letter']=='A':
  centers=sorted(set(float(t['anchor'][1]) for t in models if t['is_ytick']));pairs={}
  for t in models:
   if re.match(r'^q\s*=',t['source'].strip()):
    center=min(centers,key=lambda y:abs(y-t['anchor'][1]));pairs.setdefault(center,[]).append(t)
  for center,ts in pairs.items():
   for i,t in enumerate(sorted(ts,key=lambda q:q['anchor'][1])):t['rowcenter']=center;t['rowoffset']=(i-(len(ts)-1)/2)*9.5
 # Preserve normal line leading in multi-line labels; font glyphs and their
 # line spacing must not inherit a compressed plot's vertical aspect ratio.
 flowgroups={}
 for t in models:flowgroups.setdefault(t['group'],[]).append(t)
 for ts in flowgroups.values():
  if len(ts)<2 or len(ts)>5 or any('rowcenter' in t for t in ts):continue
  ref=ts[0]
  if not all(abs(t['angle']-ref['angle'])<.01 for t in ts):continue
  for t in ts:
   t['flow_origin']=ref['anchor'].copy();t['flow_offset']=(t['anchor']-ref['anchor'])*(BODY/ref['native_size'])
 if name=='Figure1' and sp['letter']=='B':
  header=next(t for t in models if t['group']=='text_52');point=header['anchor'];candidates=[]
  for ch in groups[0]:
   if not ch.get('id','').startswith('patch_'):continue
   path=ch.find('.//'+tag('path'))
   if path is None:continue
   b=pathbbox(path.get('d'))
   if b[0]<=point[0]<=b[2] and b[1]<=point[1]<=b[3]:candidates.append(b)
  b=min(candidates,key=lambda b:(b[2]-b[0])*(b[3]-b[1]));cx=(b[0]+b[2])/2
  for t in models:
   if t['group'] not in ['text_52','text_53','text_54']:continue
   text={'text_52':'Survival sensitivity','text_53':'Delayed entry: n = 187','text_54':'Non-autopsy: n = 189'}[t['group']]
   weight=t['parts'][0]['weight'];t['anchor'][0]=cx
   if 'flow_origin' in t:t['flow_origin'][0]=cx;t['flow_offset'][0]=0
   t['parts']=richpieces(text,BODY,weight,'normal','middle')
   for piece in t['parts']:piece['fill']='#000000'
 # At-risk tables are label-only axes. The first column has a heading;
 # subsequent columns start with an empty heading line. Preserve one shared
 # row origin even though empty SVG text nodes are intentionally omitted.
 for aid in sp['axes']:
  axg=allids['axes_'+str(aid)];ms=[t for t in models if t['axis']=='axes_'+str(aid) and t['is_xtick']]
  if ms and axisbox(axg) is None:
   top=min(t['anchor'][1] for t in ms)
   for t in ms:
    t.setdefault('flow_origin',t['anchor'].copy());t.setdefault('flow_offset',np.zeros(2));t['flow_origin'][1]=top;t['flow_offset'][1]=(t['anchor'][1]-top)*BODY/t['native_size']
   # Align the label and count columns independently. The source first column
   # combines them with spaces sized for a different font; such strings cannot
   # maintain both a right-aligned label edge and a centered year-zero count.
   km=allids.get('axes_'+str(aid-1));positions=[]
   if km is not None:
    for tick in km.xpath('.//s:g[starts-with(@id,"xtick_")]',namespaces={'s':NS}):
     us=[u for u in tick.iter(tag('use')) if u.get('x') is not None and u.get('y') is not None]
     if us:positions.append(float((fullmatrix(us[0])@np.array([float(us[0].get('x')),float(us[0].get('y')),1]))[0]))
   columns=list(dict.fromkeys(t['group'] for t in ms))
   if len(positions)==len(columns):
    zero_counts=[re.fullmatch(r'(Immune-(?:high|intermediate|low))\s+(\d+)',t['source'].strip()) for t in ms]
    gap=max(metrics(m[2],BODY)[1]/2 for m in zero_counts if m)+8
    for t in ms:
     column=columns.index(t['group']);t['flow_origin'][0]=positions[column];t['flow_offset'][0]=0
     match=re.fullmatch(r'(Immune-(?:high|intermediate|low))\s+(\d+)',t['source'].strip())
     if match:
      number=copy.deepcopy(t);number['source']=match[2];number['parts']=richpieces(match[2],BODY,'normal','normal','middle');number['risk_count_column']=0
      t['source']=match[1];t['parts']=richpieces(match[1],BODY,'normal','normal','end');t['flow_offset'][0]=-gap;t['risk_label']=True
      for q in number['parts']+t['parts']:q['fill']='#000000'
      models.append(number)
     else:
      t['parts']=richpieces(t['source'],BODY,'normal','normal','middle')
      for q in t['parts']:q['fill']='#000000'
 if name=='Figure4' and sp['letter']=='D':
  percentages=[t for t in models if re.fullmatch(r'\d+\.\d+%',t['source'].strip())]
  for p in [t for t in models if t['source'].strip().startswith('P=')]:
   q=min(percentages,key=lambda t:abs(t['anchor'][0]-p['anchor'][0]));origin=np.array([p['anchor'][0],max(p['anchor'][1],q['anchor'][1])]);p['flow_origin']=origin;p['flow_offset']=np.array([0.,0.]);q['flow_origin']=origin;q['flow_offset']=np.array([0.,-10.])
   if p['group']=='text_84':
    p['extra_y']=-20;q['extra_y']=-20
 if name=='Figure4' and sp['letter'] in ['A','B','C']:
  for t in models:
   if t['source'].strip().startswith(('χ²','Cramér')):
    t['is_footer']=True;cx=(t['axis_box'][0]+t['axis_box'][2])/2 if 'axis_box' in t else (axisbox(groups[0])[0]+axisbox(groups[0])[2])/2
    t['anchor'][0]=cx
    if 'flow_origin' in t:t['flow_origin'][0]=cx;t['flow_offset'][0]=0
    txt=re.sub(r'\bn\s*=\s*(\d+)',r'n = \1',t['source'].replace('P_MC','Monte Carlo P'));t['parts']=richpieces(txt,BODY,'normal','normal','middle')
    for piece in t['parts']:piece['fill']='#000000'
 geom=[];boxes=[]
 for source in groups:
  if source.get('id','').startswith('legend_'):continue
  g=copy.deepcopy(source)
  for n in list(g.iter(tag('text'))):n.getparent().remove(n)
  for n in list(g.iter(tag('g'))):
   if n.get('id','').startswith('legend_') and n.getparent() is not None:n.getparent().remove(n)
  for n in list(g.iter(tag('defs'))):n.getparent().remove(n)
  if name=='Figure1' and sp['letter']=='B':
   red=next(n for n in g.iter(tag('g')) if n.get('id')=='patch_34')
   red.set('transform',f'translate(877.224 546.1461) scale({(967.176-877.224)/(967.176-841.752)} {39.6/(579.4893-552.4029)}) translate(-841.752 -552.4029)')
   arrow=next(n for n in g.iter(tag('g')) if n.get('id')=='patch_24');paths=list(arrow.iter(tag('path')));start=856.224;oldtip=846.099844;tip=872.8762
   paths[0].set('transform',f'translate({start} 0) scale({(tip-start)/(oldtip-start)} 1) translate({-start} 0)')
   paths[1].set('transform',f'translate({tip} 0) scale(-1 1) translate({-oldtip} 0)')
  geom.append(g)
  b=axisbox(source)
  if b is not None:boxes.append(b)
 assert boxes,(name,sp['letter'])
 local=[];glob=[]
 for g in legends:
  title,items=legenditems(g)
  if not items:continue
  (glob if sp['legend']=='global' else local).append((title,items))
 return dict(spec=sp,models=models,geom=geom,geom_box=union(boxes),local=local,global_=glob,colorbar=allids.get('axes_'+str(sp.get('colorbar'))) if sp.get('colorbar') else None)

def fit_component(c,width,sy_override=None):
 gb=c['geom_box'];models=c['models']
 def bnd(sx,sy):
  for t in models:t['auto_x']=0;t['auto_y']=0
  for axis in c['spec']['axes']:
   ms=[t for t in models if t['axis']=='axes_'+str(axis)];groups={}
   ticks={}
   for t in ms:
    if t['is_xtick'] or t['is_ytick']:ticks.setdefault(t['group'],[]).append(t)
   for tickmodels in ticks.values():
    tick=tickmodels[0]
    if 'axis_box' not in tick:continue
    box=union([modelbox(t,sx,sy) for t in tickmodels]);ab=tick['axis_box']
    if tick['is_ytick']:
     delta=ab[2]*sx+5-box[0] if tick.get('right_ticks') else ab[0]*sx-5-box[2]
     for t in tickmodels:t['auto_x']=delta
    else:
     for t in tickmodels:t['auto_y']=ab[3]*sy+5-box[1]
     if abs(tick['angle'])<1 and 'tick_center' in tick:
      delta=tick['tick_center'][0]*sx-(box[0]+box[2])/2
      for t in tickmodels:t['auto_x']=delta
   for t in ms:
    if t['axis_kind']:groups.setdefault((t['axis_kind'],t['group']),[]).append(t)
   for (kind,_),labels in groups.items():
    ticks=[t for t in ms if t['is_ytick' if kind=='ylabel' else 'is_xtick']]
    if not ticks:continue
    tb=union([modelbox(t,sx,sy) for t in ticks]);lb=union([modelbox(t,sx,sy) for t in labels])
    right=kind=='ylabel' and labels[0]['anchor'][0]>labels[0]['axis_box'][2]
    delta=(tb[2]+5-lb[0] if right else tb[0]-5-lb[2]) if kind=='ylabel' else tb[3]+5-lb[1]
    for t in labels:t['auto_x' if kind=='ylabel' else 'auto_y']=delta
   footers=[t for t in ms if t.get('is_footer')]
   below=[t for t in ms if t['is_xtick'] or t['axis_kind']=='xlabel']
   if footers and below:
    fb=union([modelbox(t,sx,sy) for t in footers]);bb=union([modelbox(t,sx,sy) for t in below]);delta=bb[3]+9-fb[1]
    for t in footers:t['auto_y']=delta
  return union([[gb[0]*sx,gb[1]*sy,gb[2]*sx,gb[3]*sy]]+[modelbox(t,sx,sy) for t in models])
 lo=0.0001;hi=8.
 if bnd(lo,lo)[2]-bnd(lo,lo)[0]>width:raise RuntimeError(('Text intrinsically too wide',c['spec']['letter'],width))
 for _ in range(45):
  mid=(lo+hi)/2;b=bnd(mid,mid)
  if b[2]-b[0]>width:hi=mid
  else:lo=mid
 sx=lo;sy=max(sx,c['spec']['minsy']);requiredsy=c['spec']['minsy']
 if c['spec'].get('cap'):sy=min(sy,c['spec']['cap']/(gb[3]-gb[1]))
 # Ordinal row labels keep >= 10 pt spacing regardless of source canvas size.
 for ax in c['spec']['axes']:
  tickgroups={}
  for t in models:
   if t['is_ytick'] and t['axis']=='axes_'+str(ax):tickgroups.setdefault(t['tick_id'],[]).append(float(t['anchor'][1]))
  yy=sorted(round(float(np.mean(v)),4) for v in tickgroups.values())
  gaps=np.diff(yy);gaps=gaps[gaps>1]
  if len(gaps):requiredsy=max(requiredsy,c['spec'].get('minrow',10)/float(gaps.min()))
 sy=max(sy,requiredsy)
 if sy_override is not None:sy=max(sy,sy_override)
 b=bnd(sx,sy);legheight=sum(legendheight(items,width,title)+5 for title,items in c['local'])+(35 if c['colorbar'] is not None else 0)
 return dict(sx=sx,sy=sy,bbox=b,bodyheight=b[3]-b[1],legendheight=legheight,total=b[3]-b[1]+legheight)

def outline(root):
 out=copy.deepcopy(root)
 for t in list(out.iter(tag('text'))):
  st=style(t.get('style',''));text=t.text or '';s=fsize(st)
  if not text.strip():t.getparent().remove(t);continue
  p=TextPath((0,0),text,prop=prop(s,st.get('font-weight','normal'),st.get('font-style','normal')),size=s)
  adv=metrics(text,s,st.get('font-weight','normal'),st.get('font-style','normal'))[1];x=float(t.get('x',0));y=float(t.get('y',0));a=st.get('text-anchor','start');x-=adv/2 if a=='middle' else adv if a=='end' else 0
  p=p.transformed(Affine2D().scale(1,-1).translate(x,y));bits=[]
  for v,c in p.iter_segments(curves=True,simplify=False):
   cmd={MP.MOVETO:'M',MP.LINETO:'L',MP.CURVE3:'Q',MP.CURVE4:'C',MP.CLOSEPOLY:'Z'}[c];bits.append(cmd if c==MP.CLOSEPOLY else cmd+' '.join(f'{z:.6f}' for z in v))
  n=E.Element(tag('path'),d=' '.join(bits),fill=st.get('fill','#000000'))
  if t.get('transform'):n.set('transform',t.get('transform'))
  t.getparent().replace(t,n)
 return out

for name,rows in LAYOUTS.items():
 if args.only and name not in args.only:continue
 matches=list(BASE.rglob(name+'.svg'))
 if not matches:
  alias={'FigureS1_survival':'FigureS_survival','FigureS2_state_structure':'FigureS_state_structure','FigureS3_immune_features':'FigureS_immune_features'}.get(name,name);matches=list(BASE.rglob(alias+'.svg'))
 assert len(matches)==1,(name,matches)
 source=matches[0];tree=E.parse(str(source));width=89/25.4*72 if name=='FigureS5_standardization_consistency' else PAGE
 root=E.Element(tag('svg'),nsmap={None:NS,'xlink':'http://www.w3.org/1999/xlink'},version='1.1',width=f'{width}pt')
 defs=E.SubElement(root,tag('defs'));ids=set()
 for d in tree.getroot().iter(tag('defs')):
  for child in d:
   if child.get('id') not in ids:defs.append(copy.deepcopy(child));ids.add(child.get('id'))
 globals_=[];y=MARGIN;panelchecks=[]
 for row in rows:
  slots=sum(s['span'] if s else 1 for s in row);unit=(width-2*MARGIN-GUTTER*(slots-1))/slots;x=MARGIN;entries=[]
  for sp in row:
   if sp is None:x+=unit+GUTTER;continue
   cellwidth=unit*sp['span']+GUTTER*(sp['span']-1);indent=18 if sp['letter'] else 0;bodywidth=cellwidth-indent
   comp=prepare_component(tree,sp,name);fit=fit_component(comp,bodywidth);titlelines=wrap(sp['title'],cellwidth-18,TITLE) if sp['title'] else [];sublines=wrap(sp['sub'],cellwidth-18,BODY) if sp['sub'] else []
   header=max(16,12*len(titlelines)+10*len(sublines))+6 if sp['letter'] or titlelines or sublines else 0
   entries.append(dict(c=comp,f=fit,x=x,w=cellwidth,indent=indent,bw=bodywidth,head=header,title=titlelines,sub=sublines));x+=cellwidth+GUTTER
  rowheight=max(e['head']+e['f']['total'] for e in entries)
  for e in entries:
   c=e['c'];sp=c['spec'];f=e['f'];xx=e['x'];yy=y
   if sp['fillheight']:
    need=rowheight-e['head']-f['legendheight'];ratio=max(1,need/f['bodyheight']);f=fit_component(c,e['bw'],f['sy']*ratio)
   if sp['letter']:simple(root,sp['letter'],xx,yy+12,LETTER,'bold',role='panel-letter')
   ty=yy+11
   for line in e['title']:simple(root,line,xx+18,ty,TITLE,role='panel-title');ty+=12
   for line in e['sub']:simple(root,line,xx+18,ty,BODY,role='panel-subtitle');ty+=10
   bx=xx+e['indent'];by=yy+e['head'];b=f['bbox'];offsetx=bx-b[0];offsety=by-b[1]
   g=E.SubElement(root,tag('g'),transform=f'translate({offsetx:.6f} {offsety:.6f}) scale({f["sx"]:.9f} {f["sy"]:.9f})',**{'data-role':'plot-geometry','data-panel':sp['letter']})
   for n in c['geom']:
    # Scale marker positions with axes while keeping circular/square glyphs
    # undistorted under a change in panel aspect ratio.
    for u in n.iter(tag('use')):
     if u.get('x') is None or u.get('y') is None:continue
     ref=next((q for q in defs if q.get('id')==u.get('{http://www.w3.org/1999/xlink}href','').lstrip('#')),None)
     # Only marker symbols have local coordinates around the origin.
     # Filled distribution polygons use absolute data coordinates and must
     # follow the axes transform exactly, including their widths and heights.
     if ref is None or ref.tag!=tag('path') or max(abs(v) for v in pathbbox(ref.get('d','')))>20:continue
     ux,uy=float(u.get('x')),float(u.get('y'));z=max(.8,min(f['sx'],f['sy']));u.set('transform',f'translate({ux} {uy}) scale({z/f["sx"]} {z/f["sy"]}) translate({-ux} {-uy}) '+u.get('transform',''))
    g.append(n)
   for t in c['models']:addmodel(root,t,offsetx,offsety,f['sx'],f['sy'])
   ly=by+f['bodyheight']+5
   if c['colorbar'] is not None:
    cb=c['colorbar'];bb=axisbox(cb);cg=E.SubElement(root,tag('g'),transform=f'translate({bx} {ly}) scale({e["bw"]/(bb[2]-bb[0])} {6/(bb[3]-bb[1])}) translate({-bb[0]} {-bb[1]})')
    images=list(cb.iter(tag('image')))
    if images:
     for im in images:cg.append(copy.deepcopy(im))
    else:
     for ch in cb:
      if ch.get('id','').startswith(('QuadMesh_','PolyCollection_')):cg.append(copy.deepcopy(ch))
    for frac,lbl in [(0,'−1.5'),(.5,'0'),(1,'1.5')]:simple(root,lbl,bx+e['bw']*frac,ly+17,8,anchor='middle')
    simple(root,'Donor z (median)',bx+e['bw']/2,ly+29,8,anchor='middle');ly+=35
   for title,items in c['local']:drawlegend(root,items,bx,ly,e['bw'],title);ly+=legendheight(items,e['bw'],title)+5
   globals_.extend(c['global_']);panelchecks.append({'panel':sp['letter'],'axes':sp['axes'],'letter_pt':12 if sp['letter'] else None,'title_pt':10,'body_pt':8,'x':xx,'y':yy,'width':e['w'],'body_origin':[bx,by],'geometry_scale':[f['sx'],f['sy']]})
  y+=rowheight+20
 # Figure-level keys and shared state keys are drawn once, outside data areas.
 selected_axes={a for row in rows for sp in row if sp for a in sp['axes']}
 for g in tree.getroot().xpath('.//s:g[starts-with(@id,"legend_")]',namespaces={'s':NS}):
  if not any(a.get('id','').startswith('axes_') for a in g.iterancestors()):
   if name=='Figure5' and g.get('id')=='legend_3':continue
   title,items=legenditems(g)
   if items:globals_.append((title,items))
 seen=set()
 for title,items in globals_:
  key=tuple(it['label'] for it in items)
  if key in seen:continue
  seen.add(key);drawlegend(root,items,MARGIN+18,y,width-2*MARGIN-18,title);y+=legendheight(items,width-2*MARGIN-18,title)+5
 height=y+MARGIN-14;root.set('height',f'{height:.6f}pt');root.set('viewBox',f'0 0 {width:.9f} {height:.9f}')
 bg=E.Element(tag('rect'),x='0',y='0',width=str(width),height=str(height),fill='#ffffff');root.insert(0,bg)
 folder='main' if name in ['Figure'+str(i) for i in range(1,8)] else 'supplementary';dest=OUT/folder/(name+'.svg');dest.parent.mkdir(parents=True,exist_ok=True);write_output(dest,E.tostring(root,encoding='utf-8',xml_declaration=True))
 # Write through Python so native Cairo filename handling cannot break on
 # Unicode delivery-folder names on Windows.
 data=E.tostring(outline(root));write_output(dest.with_suffix('.pdf'),cairosvg.svg2pdf(bytestring=data));write_output(dest.with_suffix('.png'),cairosvg.svg2png(bytestring=data,output_width=round(width/72*400)))
 audit.append({'figure':name,'width_mm':width/72*25.4,'height_mm':height/72*25.4,'panels':panelchecks,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()});print(name,f'{width/72*25.4:.1f} × {height/72*25.4:.1f} mm',flush=True)
if not args.no_audit:
 (OUT.parent/'03_source_and_checks').mkdir(exist_ok=True)
 auditfile=OUT.parent/'03_source_and_checks/Final_figure_layout_audit.json'
 if args.only and auditfile.exists():
  merged={r['figure']:r for r in json.loads(auditfile.read_text(encoding='utf-8'))};merged.update({r['figure']:r for r in audit});audit=[merged[n] for n in LAYOUTS if n in merged]
 auditfile.write_text(json.dumps(audit,indent=2),encoding='utf-8')
