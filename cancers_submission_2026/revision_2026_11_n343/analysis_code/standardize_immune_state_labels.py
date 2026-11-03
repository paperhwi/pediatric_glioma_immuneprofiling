from pathlib import Path
import copy, re, json, math
from lxml import etree as E
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath, TextToPath
from matplotlib.transforms import Affine2D
from matplotlib.path import Path as MPath
from PIL import Image
import cairo
import cairosvg

from argparse import ArgumentParser
parser=ArgumentParser(description="Standardize publication immune-state labels without fitting models")
parser.add_argument("input_root",type=Path,help="Package root containing figures/")
parser.add_argument("output_root",type=Path,help="Destination root; may equal input root")
args=parser.parse_args()
root=args.input_root.resolve()
work=args.output_root.resolve()
work.mkdir(parents=True,exist_ok=True)
ns='http://www.w3.org/2000/svg'; tag=lambda name:'{'+ns+'}'+name
rx=re.compile(r'(?<![-\w])\b(High|Intermediate|Low)\b')
names={'High':'Immune-high','Intermediate':'Immune-intermediate','Low':'Immune-low'}
ttp=TextToPath()

def style(s):
    return dict(x.strip().split(':',1) for x in s.split(';') if ':' in x)

def fp(st):
    family=st.get('font-family','DejaVu Sans').strip().strip("'\"")
    return FontProperties(family=family,size=float(re.sub(r'[^\d.]','',st.get('font-size','8'))),weight=st.get('font-weight','normal').strip(),style=st.get('font-style','normal').strip())

def width(text,st):return ttp.get_text_width_height_descent(text,fp(st),False)[0]

def point(t):
    if 'x' in t.attrib:return float(t.get('x')),float(t.get('y','0'))
    m=re.search(r'translate\(([-\d.]+)[ ,]+([-\d.]+)\)',t.get('transform',''))
    return tuple(map(float,m.groups())) if m else (0,0)

def path_d(p):
    parts=[]
    for v,c in p.iter_segments(curves=True,simplify=False):
        cmd={MPath.MOVETO:'M',MPath.LINETO:'L',MPath.CURVE3:'Q',MPath.CURVE4:'C',MPath.CLOSEPOLY:'Z'}[c]
        parts.append(cmd if c==MPath.CLOSEPOLY else cmd+' '.join(f'{x:.6f}' for x in v))
    return ' '.join(parts)

def outlined(tree):
    out=copy.deepcopy(tree)
    for t in list(out.iter(tag('text'))):
        g=E.Element(tag('g'))
        if t.get('transform'):g.set('transform',t.get('transform'))
        base=style(t.get('style','')); x=float(t.get('x','0')); y=float(t.get('y','0'))
        def draw(text,st,x,y):
            if not text or not text.strip():return
            prop=fp(st); advance=width(text,st); anchor=st.get('text-anchor','start').strip()
            xx=x-(advance/2 if anchor=='middle' else advance if anchor=='end' else 0)
            p=TextPath((0,0),text,prop=prop,size=prop.get_size_in_points(),usetex=False)
            p=p.transformed(Affine2D().scale(1,-1).translate(xx,y))
            e=E.SubElement(g,tag('path'));e.set('d',path_d(p));e.set('fill',st.get('fill','#000000').strip())
            if 'opacity' in st:e.set('opacity',st['opacity'].strip())
        draw(t.text,base,x,y)
        for child in t:
            if child.tag!=tag('tspan'):continue
            st=base|style(child.get('style',''))
            xs=[float(z) for z in child.get('x',str(x)).split()];ys=[float(z) for z in child.get('y',str(y)).split()]
            text=''.join(child.itertext())
            if len(xs)>1:
                for j,ch in enumerate(text):draw(ch,st,xs[min(j,len(xs)-1)],ys[min(j,len(ys)-1)])
            else:draw(text,st,xs[0],ys[0])
        t.getparent().replace(t,g)
    return out

def multiline(t,lines,spacing=None,offset=0):
    st=style(t.get('style','')); size=fp(st).get_size_in_points();spacing=spacing or size*1.12
    x,y=point(t)
    # Text with a translate transform uses local coordinates.
    if 'x' not in t.attrib:x,y=0,0
    t.text=None
    for c in list(t):t.remove(c)
    for i,line in enumerate(lines):
        c=E.SubElement(t,tag('tspan'));c.set('x',str(x));c.set('y',str(y+offset+i*spacing));c.text=line

changes=[];outputs=[]
for source in sorted((root/'figures').rglob('*.svg')):
    tree=E.parse(str(source));svg=tree.getroot(); olddims=list(map(float,svg.get('viewBox').split()))
    nodes=list(svg.iter(tag('text'))); oldtexts={id(t):''.join(t.itertext()) for t in nodes}
    targets=[t for t in nodes if rx.search(oldtexts[id(t)])]
    if not targets:continue
    before_numbers=[]
    for t in targets:
        old=oldtexts[id(t)];new=rx.sub(lambda m:names[m[1]],old)
        if source.stem=='Figure5' and 'Spot immune fraction by nearest-centroid class' in old:
            new='D  Spot immune fraction by nearest immune-state centroid'
        if source.stem=='FigureS4_ScPCA_axis_linkage' and re.match(r'Low<Intermediate<High trend',old):
            new=old.replace('Low<Intermediate<High trend','Ordered immune-state trend')
        st=style(t.get('style','')); in_legend=any(a.get('id','').startswith('legend_') for a in t.iterancestors())
        # Compact tick labels retain the complete canonical name on two lines.
        if old in names and not in_legend:
            suffix=old.lower();multiline(t,['Immune-',suffix],offset=-fp(st).get_size_in_points()*.42)
        elif re.fullmatch(r'\s*(High|Intermediate|Low)\s+\d+\s*',old):
            # Keep the first at-risk count at its original right edge.
            x,y=point(t);right=x+width(old,st);t.set('transform',f'translate({right:.6f} {y})');
            t.set('style',t.get('style','')+'; text-anchor: end');t.text=new.strip()
            if source.stem=='Figure7':
                svg.set('viewBox',f'-28 0 {olddims[2]+28} {olddims[3]}');svg.set('width',f'{olddims[2]+28}pt')
        else:
            t.text=new
        assert re.findall(r'\d+(?:\.\d+)?',old)==re.findall(r'\d+(?:\.\d+)?',new),(source.name,old,new)
        changes.append({'file':source.relative_to(root).as_posix(),'old':old,'new':new})

    # Expand horizontal legend spacing around the original centre, preserving marker ownership.
    for legend in svg.iter(tag('g')):
        if not legend.get('id','').startswith('legend_'):continue
        items=[]; children=list(legend)
        for i,g in enumerate(children):
            texts=list(g.iter(tag('text')))
            if len(texts)==1 and id(texts[0]) in oldtexts and rx.search(oldtexts[id(texts[0])]):
                t=texts[0];x,y=point(t);items.append((i,g,t,x,y))
        if len(items)<2:continue
        if max(a[3] for a in items)-min(a[3] for a in items)<5:
            for j,(i,g,t,x,y) in enumerate(items):
                s=''.join(t.itertext());m=re.fullmatch(r'(Immune-(?:high|intermediate|low)) (\(n.*\))',s)
                if m:
                    size=fp(style(t.get('style',''))).get_size_in_points()
                    multiline(t,list(m.groups()),spacing=size*1.05)
                    dy=j*(size*2.5-(items[1][4]-items[0][4]))
                    for group in children[max(0,i-1):i+1]:group.set('transform',f'translate(0 {dy:.6f}) '+group.get('transform',''))
            continue
        first,last=items[0],items[-1]
        oldleft=first[3];oldright=last[3]+width(oldtexts[id(last[2])],style(last[2].get('style','')))
        starts=[oldleft]
        for a,b in zip(items,items[1:]):
            gap=b[3]-a[3]-width(oldtexts[id(a[2])],style(a[2].get('style','')))
            starts.append(starts[-1]+width(''.join(a[2].itertext()),style(a[2].get('style','')))+gap)
        newright=starts[-1]+width(''.join(last[2].itertext()),style(last[2].get('style','')))
        shift=((oldleft+oldright)-(starts[0]+newright))/2
        for (i,g,t,x,y),start in zip(items,starts):
            dx=start+shift-x
            for group in children[max(0,i-1):i+1]:group.set('transform',f'translate({dx:.6f} 0) '+group.get('transform',''))
        clamp=max(0,10-(starts[0]+shift-20))
        if source.stem=='Figure3':
            center=(starts[0]+newright)/2+shift
            legend.set('transform',f'translate({olddims[2]/2-center:.6f} 54)')
        elif clamp:legend.set('transform',f'translate({clamp:.6f} 0)')

    if source.stem=='FigureS5_standardization_consistency':
        for t in nodes:
            if oldtexts[id(t)]=='Within-343 scaled labels':t.getparent().set('transform','translate(0 8.5)')

    rel=source.relative_to(root);dest=work/rel;dest.parent.mkdir(parents=True,exist_ok=True)
    tree.write(str(dest),encoding='utf-8',xml_declaration=True)
    outlined_svg=outlined(svg)
    raw=E.tostring(outlined_svg,encoding='utf-8')
    pdf_bytes=cairosvg.svg2pdf(bytestring=raw)
    dest.with_suffix('.pdf').write_bytes(pdf_bytes)
    # Preserve the delivered PNG convention, which omits panel letters.
    png_svg=copy.deepcopy(svg)
    for t in png_svg.iter(tag('text')):
        s=''.join(t.itertext())
        if 'font-weight: 700' not in t.get('style',''):continue
        if re.fullmatch(r'[A-I]',s):t.text=''
        elif re.match(r'^[A-I]  \S',s):t.text=s[3:]
    with Image.open(source.with_suffix('.png')) as old:scale=old.height/olddims[3]
    newdims=list(map(float,svg.get('viewBox').split()))
    w,h=round(newdims[2]*scale),round(newdims[3]*scale)
    png_bytes=cairosvg.svg2png(bytestring=E.tostring(outlined(png_svg)),output_width=w,output_height=h)
    import io
    with Image.open(io.BytesIO(png_bytes)) as image:
        white=Image.new('RGB',image.size,'white')
        white.paste(image,mask=image.getchannel('A') if image.mode=='RGBA' else None)
        white.save(dest.with_suffix('.png'),dpi=(300,300))
    # Vector PDF remains vector: no new raster image replaces the panels.
    outputs.append({'file':rel.as_posix(),'viewBox':svg.get('viewBox'),'raster_images':len(list(svg.iter(tag('image')))),'png_size':[w,h]})

(work/'figure_changes.json').write_text(json.dumps(changes,ensure_ascii=False,indent=2),encoding='utf-8')
(work/'figure_outputs.json').write_text(json.dumps(outputs,indent=2),encoding='utf-8')
print(json.dumps({'figures':len(outputs),'label_changes':len(changes),'files':[r['file'] for r in outputs]},indent=2))
