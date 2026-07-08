"""
Figure 1 — Study workflow (publication-quality matplotlib substitute for BioRender).

Layout (top-to-bottom):
  1. Source dataset                : OpenPedCan v15 RNA-seq, n=702 biospecimens
  2. Cohort assembly               : Main (n=349, 4 groups) + BRAF_ALT secondary (n=353)
  3. Dual deconvolution arm        : CIBERSORTx LM22 (+ 4-method validation)  ||  Brain-tuned ssGSEA (24 sigs)
  4. Integration                   : z-scored LM22 (22 cell types) + ssGSEA (24 signatures) feature matrix
  5. Consensus clustering          : k=2-6, B=1000 bootstrap, hierarchical Ward consensus
  6. Three robust ecotypes         : Lymphocyte-inflamed / Myeloid-dominant / Immune-desert
  7. Downstream analyses           : Association (KW + Dunn) · DEG + g:Profiler · KM + Cox · BRAF projection
"""
from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_fig1"; FIG.mkdir(exist_ok=True)

# ----------------------------------------------------------------------
plt.rcParams.update({
    "font.family":    "DejaVu Sans",
    "font.size":      10,
    "axes.linewidth": 0.0,
})

# Color palette (matches downstream panels)
C_DATA      = "#0F172A"      # near-black — source
C_COHORT    = "#1D4ED8"      # blue — cohort
C_LM22      = "#7C3AED"      # purple — LM22 deconv
C_SSGSEA    = "#16A34A"      # green — brain-tuned ssGSEA
C_INTEGRATE = "#0EA5E9"      # cyan — integration
C_CLUSTER   = "#F59E0B"      # amber — clustering
C_LYM       = "#3B82F6"      # ecotypes
C_MYE       = "#EF4444"
C_DES       = "#9CA3AF"
C_DOWN      = "#475569"      # slate — downstream

W, H = 16.5, 11.5     # inches
fig = plt.figure(figsize=(W, H), dpi=200)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_axis_off()

def box(x, y, w, h, text, *, face, edge=None, fontsize=10, bold=False, italic=False,
        text_color="white", lh=1.2, rounding=0.04):
    """Add a rounded box centered at (x, y) with width/height (in coord units) and label."""
    bb = FancyBboxPatch(
        (x - w/2, y - h/2), w, h,
        boxstyle=f"round,pad=0.4,rounding_size={rounding*100}",
        linewidth=1.0,
        facecolor=face,
        edgecolor=edge or face,
    )
    ax.add_patch(bb)
    weight = "bold" if bold else "normal"
    style  = "italic" if italic else "normal"
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color=text_color,
            weight=weight, style=style, linespacing=lh)
    return bb

def arrow(x1, y1, x2, y2, *, color="#334155", lw=1.6, head=14, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                                 arrowstyle=f"-|>,head_width={head*0.05},head_length={head*0.08}",
                                 mutation_scale=head, linewidth=lw, linestyle=ls,
                                 color=color, shrinkA=0, shrinkB=0))

def section_label(x, y, text, color):
    ax.text(x, y, text, ha="left", va="center", fontsize=9, color=color, weight="bold",
            style="italic")

# Title strip
ax.text(50, 97, "Figure 1.   Study workflow: a developmentally integrated immune ecotype framework for pediatric glioma",
        ha="center", va="center", fontsize=13, weight="bold", color="#0F172A")
ax.text(50, 94.0, "OpenPedCan v15 bulk RNA-seq → dual deconvolution → integrated feature matrix → consensus clustering → three robust ecotypes → downstream analyses",
        ha="center", va="center", fontsize=10, color="#475569", style="italic")

# ---------------- Row 1: Source dataset ----------------
box(50, 88, 60, 4.5,
    "OpenPedCan v15  ·  bulk RNA-seq (TPM, RSEM-collapsed HGNC)\n"
    "n = 702 biospecimens",
    face=C_DATA, fontsize=11, bold=True)

# ---------------- Row 2: Cohort assembly ----------------
box(28, 79, 32, 6.5,
    "MAIN cohort (independent-primary-plus)\n"
    "n = 349   |   DMG_K27 175 · DHG_G34 31 · pHGG_WT 125 · IHG 18\n"
    "pHGG_WT age < 21 y;  ambiguous-location flagged",
    face=C_COHORT, fontsize=9.5)

box(72, 79, 32, 6.5,
    "SECONDARY  BRAF / RTK-altered\n"
    "BRAF_ALT_main  n = 31   (HGG / PXA / infant-fusion)\n"
    "BRAF_ALT_LGG   n = 322  (LGG, supplementary projection)",
    face=C_COHORT, fontsize=9.5)

arrow(50, 86, 28, 82.5)
arrow(50, 86, 72, 82.5)

# ---------------- Row 3: Dual deconvolution arms ----------------
section_label(2.5, 71.5, "Dual deconvolution arm", "#0F172A")

# LM22 arm
box(28, 67, 32, 9,
    "CIBERSORTx  LM22  (22 immune cell types)\n"
    "abs mode · 100 perms · QN off · batch off\n"
    "median P = 0.15, median r = 0.09\n"
    "→ peripheral panel ill-fitted to brain TME",
    face=C_LM22, fontsize=9.5)

# Cross-validation arm
box(50, 67, 22, 9,
    "Cross-method validation\n"
    "xCell · quanTIseq\n"
    "MCP-counter · EPIC\n"
    "Spearman ρ on shared cell types",
    face="#A78BFA", fontsize=9.5, text_color="#0F172A")

# ssGSEA arm
box(74, 67, 28, 9,
    "Brain-tuned  ssGSEA  (24 signatures)\n"
    "Microglia (Bohlen 2020 · Klemm 2020 · Antunes 2021 · DAM)\n"
    "MDM · MgTAM · MoTAM\n"
    "MHC I/II · IFN α/γ · MAPK · T-cytotox · M1/M2 · …",
    face=C_SSGSEA, fontsize=9.5)

arrow(28, 76, 28, 71.7)
arrow(50, 76, 50, 71.7)
arrow(72, 76, 74, 71.7)

# ---------------- Row 4: Integration ----------------
box(50, 57.5, 56, 5.0,
    "Integrated feature matrix\n"
    "z-scored  LM22 (22)  +  ssGSEA brain-immune (24)   ·   349 samples × 46 features",
    face=C_INTEGRATE, fontsize=10, bold=True)
arrow(28, 62.5, 50, 60.0)
arrow(50, 62.5, 50, 60.0)
arrow(74, 62.5, 50, 60.0)

# ---------------- Row 5: Consensus clustering ----------------
box(50, 49, 46, 5.2,
    "Consensus clustering  ·  k = 2–6  ·  B = 1000 bootstraps  ·  pItem = 0.8\n"
    "k-means + hierarchical Ward consensus tree  ·  PAC + silhouette  ·  pre-specified k = 3",
    face=C_CLUSTER, fontsize=10, bold=True, text_color="white")
arrow(50, 55, 50, 51.6)

# ---------------- Row 6: Three ecotypes ----------------
ec_y = 38
box(20, ec_y, 24, 7.2,
    "Lymphocyte-inflamed\nn = 111\nT-cytotox · IFN-γ · MHC · Chemokine",
    face=C_LYM, fontsize=10, bold=True)
box(50, ec_y, 24, 7.2,
    "Myeloid-dominant\nn = 160\nMo-TAM · MDM · M1/M2 · Glioma inflammatory",
    face=C_MYE, fontsize=10, bold=True)
box(80, ec_y, 24, 7.2,
    "Immune-desert\nn = 78\nstemness · cell cycle\nlow immune signal",
    face=C_DES, fontsize=10, bold=True, text_color="#0F172A")
arrow(50, 46.4, 20, 41.6)
arrow(50, 46.4, 50, 41.6)
arrow(50, 46.4, 80, 41.6)

# Robustness annotation under cluster row
ax.text(50, 29.8,
        "Robustness: within-cluster mean consensus  k=2 → 0.95  · k=3 → 0.86  · k=4 → 0.71",
        ha="center", va="center", fontsize=9, color="#475569", style="italic")
ax.text(50, 28.0,
        "Sensitivity: P<0.05 subset n=91 → 3 ecotypes retained;   Full pool n=702 de-novo → 91.1 % agreement with Main",
        ha="center", va="center", fontsize=9, color="#475569", style="italic")

# ---------------- Row 7: Downstream analyses ----------------
section_label(2.5, 22.5, "Downstream", "#0F172A")
ds_y = 17.5
dw, dh = 21, 7.5

box(13.5, ds_y, dw, dh,
    "Association testing\nEcotype × cohort_group  p = 2×10⁻⁴\nEcotype × location  p = 6×10⁻⁴\nEcotype × age_dev_group  p = 0.43 (NS)\nKW + Dunn on 105 features (84 q<0.05)",
    face=C_DOWN, fontsize=8.5, lh=1.3)

box(38.5, ds_y, dw, dh,
    "Per-ecotype DEG + enrichment\n"
    "Wilcoxon one-vs-rest, BH q<0.05\n"
    "Lymph 1,456 / Mye 50 / Desert 384 up\n"
    "g:Profiler: GO:BP · REAC · KEGG · WP",
    face=C_DOWN, fontsize=8.5, lh=1.3)

box(63.5, ds_y, dw, dh,
    "Survival (exploratory)\n"
    "KM log-rank p = 0.028 (n=251)\n"
    "Multivariable Cox\n"
    "Lymph-inflamed  HR 0.55  p = 0.003",
    face=C_DOWN, fontsize=8.5, lh=1.3)

box(88.5, ds_y, dw, dh,
    "BRAF_ALT projection\n"
    "Nearest-centroid → Main ecotype space\n"
    "BRAF_main: 61 % Lymph / 0 % Desert\n"
    "MAPK_Activity KW p = 5×10⁻²³",
    face=C_DOWN, fontsize=8.5, lh=1.3)

arrow(20, 34.4, 13.5, 21.5)
arrow(50, 34.4, 38.5, 21.5)
arrow(50, 34.4, 63.5, 21.5)
arrow(80, 34.4, 88.5, 21.5)

# ---------------- Bottom strip: outputs ----------------
ax.text(50, 9.6,
        "Outputs:  ecotype assignment table  ·  per-cell-type abundance / signature score matrix  ·  per-ecotype marker genes + pathway enrichment  ·  multivariable Cox model  ·  BRAF_ALT projection",
        ha="center", va="center", fontsize=8.5, color="#475569")
ax.text(50, 7.4,
        "Target:  developmentally informed immune ecotype framework integrating histone class × anatomical compartment × age",
        ha="center", va="center", fontsize=10, color="#0F172A", weight="bold", style="italic")

# Bottom legend
handles = [
    Line2D([0], [0], marker="s", linestyle="", markerfacecolor=C_LYM,    markeredgecolor=C_LYM,    markersize=10, label="Lymphocyte-inflamed"),
    Line2D([0], [0], marker="s", linestyle="", markerfacecolor=C_MYE,    markeredgecolor=C_MYE,    markersize=10, label="Myeloid-dominant"),
    Line2D([0], [0], marker="s", linestyle="", markerfacecolor=C_DES,    markeredgecolor=C_DES,    markersize=10, label="Immune-desert"),
    Line2D([0], [0], marker="s", linestyle="", markerfacecolor=C_LM22,   markeredgecolor=C_LM22,   markersize=10, label="LM22 (peripheral panel)"),
    Line2D([0], [0], marker="s", linestyle="", markerfacecolor=C_SSGSEA, markeredgecolor=C_SSGSEA, markersize=10, label="Brain-tuned ssGSEA"),
]
leg = ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
                ncol=5, frameon=False, fontsize=9)

png = FIG / "Fig1_workflow.png"
svg = FIG / "Fig1_workflow.svg"
fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(svg, bbox_inches="tight", facecolor="white")
print(f"Wrote: {png}")
print(f"Wrote: {svg}")
plt.close(fig)
