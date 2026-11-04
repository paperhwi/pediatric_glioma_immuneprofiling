from pathlib import Path
import copy, re, json, math
from lxml import etree as E
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath, TextToPath
from matplotlib.transforms import Affine2D
from matplotlib.path import Path as MPath
from PIL import Image
import pymupdf as fitz
import cairo
import cairosvg
ns = 'http://www.w3.org/2000/svg'
tag = lambda name: '{' + ns + '}' + name
ttp = TextToPath()

def style(s):
    return dict((x.strip().split(':', 1) for x in s.split(';') if ':' in x))

def fp(st):
    family = st.get('font-family', 'DejaVu Sans').strip().strip('\'"')
    return FontProperties(family=family, size=float(re.sub('[^\\d.]', '', st.get('font-size', '8'))), weight=st.get('font-weight', 'normal').strip(), style=st.get('font-style', 'normal').strip())

def width(text, st):
    return ttp.get_text_width_height_descent(text, fp(st), False)[0]

def point(t):
    if 'x' in t.attrib:
        return (float(t.get('x')), float(t.get('y', '0')))
    m = re.search('translate\\(([-\\d.]+)[ ,]+([-\\d.]+)\\)', t.get('transform', ''))
    return tuple(map(float, m.groups())) if m else (0, 0)

def path_d(p):
    parts = []
    for v, c in p.iter_segments(curves=True, simplify=False):
        cmd = {MPath.MOVETO: 'M', MPath.LINETO: 'L', MPath.CURVE3: 'Q', MPath.CURVE4: 'C', MPath.CLOSEPOLY: 'Z'}[c]
        parts.append(cmd if c == MPath.CLOSEPOLY else cmd + ' '.join((f'{x:.6f}' for x in v)))
    return ' '.join(parts)

def outlined(tree):
    out = copy.deepcopy(tree)
    for t in list(out.iter(tag('text'))):
        g = E.Element(tag('g'))
        if t.get('transform'):
            g.set('transform', t.get('transform'))
        base = style(t.get('style', ''))
        x = float(t.get('x', '0'))
        y = float(t.get('y', '0'))

        def draw(text, st, x, y):
            if not text or not text.strip():
                return
            prop = fp(st)
            advance = width(text, st)
            anchor = st.get('text-anchor', 'start').strip()
            xx = x - (advance / 2 if anchor == 'middle' else advance if anchor == 'end' else 0)
            p = TextPath((0, 0), text, prop=prop, size=prop.get_size_in_points(), usetex=False)
            p = p.transformed(Affine2D().scale(1, -1).translate(xx, y))
            e = E.SubElement(g, tag('path'))
            e.set('d', path_d(p))
            e.set('fill', st.get('fill', '#000000').strip())
            if 'opacity' in st:
                e.set('opacity', st['opacity'].strip())
        draw(t.text, base, x, y)
        for child in t:
            if child.tag != tag('tspan'):
                continue
            st = base | style(child.get('style', ''))
            xs = [float(z) for z in child.get('x', str(x)).split()]
            ys = [float(z) for z in child.get('y', str(y)).split()]
            text = ''.join(child.itertext())
            if len(xs) > 1:
                for j, ch in enumerate(text):
                    draw(ch, st, xs[min(j, len(xs) - 1)], ys[min(j, len(ys) - 1)])
            else:
                draw(text, st, xs[0], ys[0])
        t.getparent().replace(t, g)
    return out

def multiline(t, lines, spacing=None, offset=0):
    st = style(t.get('style', ''))
    size = fp(st).get_size_in_points()
    spacing = spacing or size * 1.12
    x, y = point(t)
    if 'x' not in t.attrib:
        x, y = (0, 0)
    t.text = None
    for c in list(t):
        t.remove(c)
    for i, line in enumerate(lines):
        c = E.SubElement(t, tag('tspan'))
        c.set('x', str(x))
        c.set('y', str(y + offset + i * spacing))
        c.text = line
