"""Figure 1 for the Cancers submission.

A  Study workflow.                 B  Analysis-set derivation.
Both panels are redrawn here so that they share one visual language, one font
and US spelling.  Output is a single vector PDF plus a 600 dpi PNG.
"""
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

NAVY   = "#1F3864"
FILL   = "#E4ECF7"
FILL2  = "#C3D5EC"
AMB_E, AMB_F = "#B07C12", "#F7E7BE"
EXCL_E, EXCL_F = "#A33A3A", "#FAE4E4"
LYM, MYE, DES = "#3B82F6", "#EF4444", "#9CA3AF"
TXT, GREY = "#0F172A", "#56657F"

FIG_W, FIG_H = 13.4, 9.9
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
axA = fig.add_axes([0.0, 0.487, 1.0, 0.513])
axB = fig.add_axes([0.0, 0.0, 1.0, 0.475])
for ax in (axA, axB):
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

# points per y data-unit, per panel (100 units span the axes height)
UY = {id(axA): FIG_H * 0.513 * 72 / 100.0,
      id(axB): FIG_H * 0.475 * 72 / 100.0}


def box(ax, x, y, w, h, lines, *, fill=FILL, edge=NAVY, fs=8.6, lw=1.1,
        dashed=False, tc=TXT, weights=None, lead=1.38):
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0,rounding_size=1.0",
        linewidth=lw, facecolor=fill, edgecolor=edge,
        linestyle=(0, (3.5, 2.5)) if dashed else "-", zorder=2))
    if isinstance(lines, str):
        lines = [lines]
    weights = weights or ["normal"] * len(lines)
    lh = fs * lead / UY[id(ax)]                 # line height in y-units
    top = y + lh * (len(lines) - 1) / 2
    for i, ln in enumerate(lines):
        ax.text(x, top - i * lh, ln, ha="center", va="center", fontsize=fs,
                color=tc, weight=weights[i], zorder=3)


def arrow(ax, p, q, color=NAVY, lw=1.3, ms=11):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms,
                                 linewidth=lw, color=color, shrinkA=0,
                                 shrinkB=0, zorder=1))


def bus(ax, xs, y_from, y_bus, x_to, y_to, color=NAVY, lw=1.3):
    """drop several sources onto a horizontal bus, then one arrow down."""
    for x in xs:
        ax.plot([x, x], [y_from, y_bus], color=color, lw=lw, zorder=1)
    ax.plot([min(xs), max(xs)], [y_bus, y_bus], color=color, lw=lw, zorder=1)
    arrow(ax, (x_to, y_bus), (x_to, y_to), color=color, lw=lw)


def fan(ax, x_from, y_from, y_bus, targets, y_to, color=NAVY, lw=1.3):
    ax.plot([x_from, x_from], [y_from, y_bus], color=color, lw=lw, zorder=1)
    ax.plot([min(targets), max(targets)], [y_bus, y_bus], color=color, lw=lw, zorder=1)
    for x in targets:
        arrow(ax, (x, y_bus), (x, y_to), color=color, lw=lw)


def merge(ax, sources, y_from, y_bus, x_to, y_to, color=NAVY, lw=1.3):
    for x in sources:
        ax.plot([x, x], [y_from, y_bus], color=color, lw=lw, zorder=1)
    ax.plot([min(sources), max(sources)], [y_bus, y_bus], color=color, lw=lw, zorder=1)
    arrow(ax, (x_to, y_bus), (x_to, y_to), color=color, lw=lw)


# ============================================================== panel A
A = axA
A.text(0.8, 99, "A", fontsize=15, weight="bold", ha="left", va="top", color=TXT)
A.text(50, 99, "Study workflow", fontsize=11.5, weight="bold",
       ha="center", va="top", color=TXT)
A.text(50, 93.6, "OpenPedCan v15  ·  pediatric CNS tumor biospecimens",
       ha="center", va="center", fontsize=8.6, style="italic", color=GREY)

xin = [18, 50, 82]
for cx, txt in zip(xin, [["Clinical and histologic", "annotation"],
                         ["Fusion", "annotation"],
                         ["TPM-normalized", "gene-expression matrix"]]):
    box(A, cx, 85.0, 27, 8.6, txt, fill=FILL2, fs=8.8)

box(A, 50, 70.5, 80, 11.0,
    ["Cohort assembly   n = 349",
     "one independent primary biospecimen per patient",
     "DMG, H3 K27-altered 175 · pHGG, H3/IDH-wildtype 125 · DHG, H3 G34-mutant 31 · infant-type hemispheric 18"],
    fs=8.6, weights=["bold", "normal", "normal"])
bus(A, xin, 80.7, 78.2, 50, 76.0)

box(A, 50, 55.0, 88, 11.0,
    ["Immune-feature generation",
     "CIBERSORTx LM22 (22 immune cell types)  ·  orthogonal deconvolution (xCell, quanTIseq, MCP-counter, EPIC)",
     "24 brain-tuned ssGSEA signatures  →  z-scored 46-feature immune matrix"],
    fs=8.6, weights=["bold", "normal", "normal"])
arrow(A, (50, 65.0), (50, 60.7))

box(A, 50, 36.5, 92, 14.0,
    ["Bootstrap consensus clustering",
     "1,000 bootstrap replicates, 80% item resampling, average linkage on the 1 − consensus distance",
     "k = 2–10 evaluated (PAC, silhouette, consensus-CDF ΔAUC, gap statistic, minimum cluster size)",
     "k = 2 most stable;  k = 3 retained as an exploratory three-state representation"],
    fill=AMB_F, edge=AMB_E, fs=8.6,
    weights=["bold", "normal", "normal", "bold"])
arrow(A, (50, 49.5), (50, 43.8))

eco = [(19, "Lymphocyte-inflamed", "n = 111", LYM),
       (50, "Myeloid-dominant", "n = 160", MYE),
       (81, "Immune-desert", "n = 78", DES)]
for cx, nm, n, col in eco:
    box(A, cx, 20.5, 28, 8.4, [nm, n], fill=col, edge=col, fs=8.8,
        tc="white", weights=["bold", "normal"])
fan(A, 50, 29.5, 26.6, [c for c, *_ in eco], 24.8)

box(A, 50, 6.2, 92, 11.6,
    ["Downstream analyses",
     "integrated molecular group and anatomical location  ·  sequential PERMANOVA  ·  molecular-group-adjusted differential expression",
     "Kaplan–Meier and multivariable Cox  ·  single-cell and spatial support  ·  sensitivity analyses"],
    fs=8.6, weights=["bold", "normal", "normal"])
merge(A, [c for c, *_ in eco], 16.3, 14.0, 50, 12.1)

# ============================================================== panel B
B = axB
B.text(0.8, 99, "B", fontsize=15, weight="bold", ha="left", va="top", color=TXT)
B.text(50, 99, "Analysis-set derivation", fontsize=11.5, weight="bold",
       ha="center", va="top", color=TXT)

box(B, 40, 89.0, 40, 9.0,
    ["OpenPedCan v15", "pediatric CNS tumor biospecimens"],
    fs=8.8, weights=["bold", "normal"])
box(B, 81, 89.0, 34, 8.6,
    ["Excluded: non-primary and repeat biospecimens;",
     "diagnoses outside the four integrated groups"],
    fill=EXCL_F, edge=EXCL_E, dashed=True, fs=7.5)
arrow(B, (60.2, 89.0), (63.7, 89.0), color=EXCL_E)

box(B, 40, 73.5, 54, 11.5,
    ["Main cohort   n = 349",
     "one independent primary biospecimen per patient",
     "DMG_K27 175 · pHGG_WT 125 · DHG_G34 31 · IHG 18"],
    fill=FILL2, fs=8.6, weights=["bold", "normal", "normal"])
arrow(B, (40, 84.5), (40, 79.5))

B.text(50, 63.6,
       "immune deconvolution, signature scoring and consensus clustering used all 349 tumors",
       ha="center", va="center", fontsize=8.2, style="italic", color=GREY)

fan(B, 40, 67.8, 58.5, [33, 72], 52.5)

box(B, 33, 45.5, 32, 13.0,
    ["Anatomical-location analyses", "n = 332",
     "contingency tests and Cramér's V"],
    fs=8.6, weights=["bold", "normal", "normal"])
box(B, 72, 45.5, 32, 13.0,
    ["Survival analyses", "n = 251  (208 deaths)",
     "Kaplan–Meier and multivariable Cox"],
    fs=8.6, weights=["bold", "normal", "normal"])
box(B, 8.0, 45.5, 15, 8.0, ["− 17 no usable", "location annotation"],
    fill=EXCL_F, edge=EXCL_E, dashed=True, fs=7.2)
box(B, 93.0, 45.5, 13, 8.0, ["− 98 no overall-", "survival annotation"],
    fill=EXCL_F, edge=EXCL_E, dashed=True, fs=7.2)
arrow(B, (17.0, 45.5), (15.8, 45.5), color=EXCL_E)
arrow(B, (88.0, 45.5), (86.8, 45.5), color=EXCL_E)

box(B, 33, 18.0, 32, 13.0,
    ["Sequential PERMANOVA", "n = 258",
     "midline, hemispheric, posterior fossa"],
    fs=8.6, weights=["bold", "normal", "normal"])
box(B, 72, 18.0, 32, 13.0,
    ["Missingness sensitivity", "n = 349",
     "inverse-probability-of-observation weighting"],
    fs=8.6, weights=["bold", "normal", "normal"])
box(B, 8.0, 18.0, 15, 8.0, ["− 74 multi-compartment", "or ambiguous"],
    fill=EXCL_F, edge=EXCL_E, dashed=True, fs=7.0)
arrow(B, (33, 39.0), (33, 24.7))
arrow(B, (72, 39.0), (72, 24.7))
arrow(B, (17.0, 18.0), (15.8, 18.0), color=EXCL_E)

B.text(50, 4.0,
       "The two branches are not nested: of the 251 survival-evaluable tumors, 180 fall within the 258-tumor location subset and 3 have no usable location annotation.",
       ha="center", va="center", fontsize=7.8, color=GREY, style="italic")

os.makedirs(f"{OUT}/panels", exist_ok=True)
fig.savefig(f"{OUT}/panels/Figure1_full.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/panels/_preview_fig1.png", dpi=110, bbox_inches="tight")
print("Figure 1 drawn")
