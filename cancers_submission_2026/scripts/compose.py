"""Vector figure composition for the Cancers submission.

Places existing single-panel PDFs (and newly drawn ones) into one multi-panel
vector PDF, adds panel letters and optional panel titles, and renders a
high-resolution PNG.  Nothing is rasterised in the PDF path.
"""
from __future__ import annotations
import os
import numpy as np
import pymupdf

PT = 72.0
LETTER_SIZE = 12.5
TITLE_SIZE = 8.6
LETTER_FONT = "Helvetica-Bold"
TITLE_FONT = "Helvetica"
GAP = 12.0
ROW_GAP = 16.0
MARGIN = 9.0
LETTER_PAD = 5.0


def content_clip(src, pad=5.0, probe_dpi=100, thresh=248, within=None):
    """Tight bounding box of the ink on page 0, in source points."""
    d = pymupdf.open(src)
    page = d[0]
    box = pymupdf.Rect(within) if within else page.rect
    pm = page.get_pixmap(dpi=probe_dpi, clip=box, colorspace=pymupdf.csGRAY)
    a = np.frombuffer(pm.samples, dtype=np.uint8).reshape(pm.height, pm.width)
    ys, xs = np.where(a < thresh)
    d.close()
    if len(xs) == 0:
        return tuple(box)
    sc = PT / probe_dpi
    return (max(box.x0, box.x0 + xs.min() * sc - pad),
            max(box.y0, box.y0 + ys.min() * sc - pad),
            min(box.x1, box.x0 + (xs.max() + 1) * sc + pad),
            min(box.y1, box.y0 + (ys.max() + 1) * sc + pad))


class Panel:
    """One source panel.

    src    : single-page PDF
    label  : panel letter (None for none)
    title  : short title drawn beside the letter (None for none)
    clip   : (x0, y0, x1, y1) crop in source points, before trimming
    trim   : shrink the clip to the ink it contains
    scale  : size multiplier
    """

    def __init__(self, src, label=None, title=None, clip=None, trim=True,
                 scale=1.0):
        self.src = src
        self.label = label
        self.title = title
        self.scale = scale
        if trim:
            clip = content_clip(src, within=clip)
        elif clip is None:
            d = pymupdf.open(src)
            clip = tuple(d[0].rect)
            d.close()
        self.clip = clip
        self.w = (clip[2] - clip[0]) * scale
        self.h = (clip[3] - clip[1]) * scale


def compose(rows, out_pdf, png=None, dpi=600, max_px=12000,
            align="center", row_align="center", gap=GAP, row_gap=ROW_GAP,
            margin=MARGIN):
    """rows: list of list[Panel]."""
    GAP, ROW_GAP, MARGIN = gap, row_gap, margin
    heads = [LETTER_PAD + LETTER_SIZE if any(p.label or p.title for p in r) else 0.0
             for r in rows]
    row_w = [sum(p.w for p in r) + GAP * (len(r) - 1) for r in rows]
    row_h = [max(p.h for p in r) + hd for r, hd in zip(rows, heads)]
    W = max(row_w) + 2 * MARGIN
    H = sum(row_h) + ROW_GAP * (len(rows) - 1) + 2 * MARGIN

    out = pymupdf.open()
    page = out.new_page(width=W, height=H)
    page.draw_rect(page.rect, color=None, fill=(1, 1, 1))

    y = MARGIN
    for r, rw, rh, hd in zip(rows, row_w, row_h, heads):
        x = MARGIN + (max(row_w) - rw) / 2 if align == "center" else MARGIN
        for p in r:
            py = y + hd if row_align == "top" else y + hd + (rh - hd - p.h) / 2
            page.show_pdf_page(pymupdf.Rect(x, py, x + p.w, py + p.h),
                               pymupdf.open(p.src), 0,
                               clip=pymupdf.Rect(p.clip))
            tx = x
            if p.label:
                page.insert_text((tx, y + LETTER_SIZE), p.label,
                                 fontname=LETTER_FONT, fontsize=LETTER_SIZE,
                                 color=(0, 0, 0))
                tx += pymupdf.get_text_length(p.label, LETTER_FONT,
                                              LETTER_SIZE) + 6
            if p.title:
                page.insert_text((tx, y + LETTER_SIZE - 0.6), p.title,
                                 fontname=TITLE_FONT, fontsize=TITLE_SIZE,
                                 color=(0.06, 0.09, 0.16))
            x += p.w + GAP
        y += rh + ROW_GAP

    os.makedirs(os.path.dirname(os.path.abspath(out_pdf)), exist_ok=True)
    out.save(out_pdf, garbage=4, deflate=True)

    png_path = None
    if png:
        use = dpi
        longest = max(W, H) / PT
        if longest * use > max_px:
            use = max_px / longest
        page.get_pixmap(dpi=int(use), colorspace=pymupdf.csRGB).save(png)
        png_path = png
    out.close()
    return out_pdf, png_path


def preview(pdf, out_png, dpi=110):
    d = pymupdf.open(pdf)
    d[0].get_pixmap(dpi=dpi, colorspace=pymupdf.csRGB).save(out_png)
    d.close()
    return out_png


def report(path, png=None):
    d = pymupdf.open(path)
    p = d[0]
    msg = (f"{os.path.basename(path):<36} {p.rect.width/PT:6.2f} x "
           f"{p.rect.height/PT:6.2f} in   pdf {os.path.getsize(path)/1e6:5.2f} MB")
    if png and os.path.exists(png):
        msg += f"   png {os.path.getsize(png)/1e6:5.2f} MB"
    print(msg)
    d.close()
