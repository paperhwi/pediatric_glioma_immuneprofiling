"""Figure 5C — each in-tissue Visium spot coloured by the bulk ecotype centroid
it maps most closely to, for the three H3.3 K27M pediatric diffuse midline
glioma sections.

The scatter artwork is taken, as vector content, from the spatial analysis
figure produced by the spatial notebook (w2g_fig_spatial.py); the titles,
labels and legend are redrawn here so that the panel matches the rest of the
Cancers figure set and uses US spelling.
"""
import os
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pymupdf

mpl.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42,
                     "ps.fonttype": 42})

SRC = ("/mnt/user-data/uploads/Open PBTA/Revision/REVISION PACKAGE 260916/"
       "1. Figure 300dpi/Supplementary figures/FigureS21_spatial.pdf")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "panels")
os.makedirs(OUT, exist_ok=True)

LYM, MYE, DES = "#3B82F6", "#EF4444", "#9CA3AF"
PAD = 3.0
CROPS = {  # ink bounds measured on the source page, in points
    1: (126.9, 103.4, 252.2, 179.0),
    2: (323.8, 77.3, 449.1, 205.1),
    3: (520.4, 101.8, 646.1, 180.6),
}
META = {1: (577, 0.540, 0.358, 0.001),
        2: (856, 0.357, 0.338, 0.019),
        3: (217, 0.562, 0.370, 0.001)}

FW, FH = 11.9, 3.75                     # inches
W, H = FW * 72, FH * 72                 # points
TOP = 53.0                              # points reserved for titles
BOT = 8.0
COL = W / 3.0

crops = {k: (v[0] - PAD, v[1] - PAD, v[2] + PAD, v[3] + PAD)
         for k, v in CROPS.items()}
cw = max(c[2] - c[0] for c in crops.values())
ch = max(c[3] - c[1] for c in crops.values())
scale = min((COL - 26) / cw, (H - TOP - BOT) / ch)

# ---- text layer -----------------------------------------------------------
fig = plt.figure(figsize=(FW, FH), dpi=300)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

for i, k in enumerate((1, 2, 3)):
    n, obs, null, p = META[k]
    cx = COL * i + COL / 2
    ax.text(cx, H - 14, f"pDMG Sample-{k}   ·   {n} in-tissue spots",
            ha="center", va="center", fontsize=8.2, weight="bold", color="#0F172A")
    ax.text(cx, H - 25, f"same-neighbor fraction {obs:.3f} vs permutation null {null:.3f},  "
                        f"$P$ = {p:.3f}",
            ha="center", va="center", fontsize=7.4, color="#334155")

# horizontal legend just under the titles
lx = W / 2 - 32
ly = H - 41
ax.text(lx - 6, ly, "spot assigned to the nearest bulk ecotype centroid:",
        ha="right", va="center", fontsize=7.2, color="#334155")
for dx, col, lab in [(0, LYM, "Lymphocyte-inflamed"), (104, MYE, "Myeloid-dominant"),
                     (196, DES, "Immune-desert")]:
    ax.plot([lx + dx], [ly], marker="o", ms=4.5, color=col, lw=0)
    ax.text(lx + dx + 6, ly, lab, ha="left", va="center", fontsize=7.2, color="#0F172A")

fig.savefig(f"{OUT}/_fig5_C_text.pdf")
plt.close(fig)

# ---- overlay the vector scatters -----------------------------------------
dst = pymupdf.open(f"{OUT}/_fig5_C_text.pdf")
page = dst[0]
src = pymupdf.open(SRC)
band_top, band_bot = TOP, H - BOT
for i, k in enumerate((1, 2, 3)):
    x0, y0, x1, y1 = crops[k]
    w, h = (x1 - x0) * scale, (y1 - y0) * scale
    cx = COL * i + COL / 2
    top = band_top + ((band_bot - band_top) - h) / 2
    page.show_pdf_page(pymupdf.Rect(cx - w / 2, top, cx + w / 2, top + h),
                       src, 0, clip=pymupdf.Rect(x0, y0, x1, y1))
dst.save(f"{OUT}/fig5_C.pdf", garbage=4, deflate=True)
src.close(); dst.close()
os.remove(f"{OUT}/_fig5_C_text.pdf")
print(f"Figure 5C written (scale {scale:.2f})")
