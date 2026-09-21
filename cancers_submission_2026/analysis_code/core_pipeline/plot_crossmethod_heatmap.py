"""Plot Spearman correlation heatmaps for cross-method myeloid/CD8/NK signals."""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

mpl.rcParams['font.family'] = 'DejaVu Sans'

OUT = Path(r"C:\Users\scott\Downloads\Open PBTA\output")

def plot_corr(tsv, title, fname, figsize=(8, 7)):
    df = pd.read_csv(tsv, sep=None, engine="python")
    df.set_index("method", inplace=True)
    fig, ax = plt.subplots(figsize=figsize, dpi=200)
    im = ax.imshow(df.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(df.columns)))
    ax.set_yticks(range(len(df.index)))
    ax.set_xticklabels(df.columns, rotation=40, ha="right", fontsize=9)
    ax.set_yticklabels(df.index, fontsize=9)
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            v = df.values[i, j]
            color = "white" if abs(v) > 0.6 else "black"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", color=color, fontsize=8, weight="bold")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Spearman ρ", fontsize=10)
    ax.set_title(title, fontsize=12, weight="bold", pad=12)
    plt.tight_layout()
    plt.savefig(OUT / fname, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Wrote {fname}")

plot_corr(OUT / "spearman_myeloid_crossmethod.tsv",
          "Myeloid signals: cross-method Spearman correlation (n=702)\nLM22-like methods vs. brain-tuned signatures",
          "fig_corr_myeloid.png", figsize=(9.5, 8))

plot_corr(OUT / "spearman_cd8_crossmethod.tsv",
          "CD8 T-cell / cytolytic signals: cross-method correlation",
          "fig_corr_cd8.png", figsize=(7, 6))

plot_corr(OUT / "spearman_nk_crossmethod.tsv",
          "NK cell signals: cross-method correlation",
          "fig_corr_nk.png", figsize=(7, 6))

# Combined panel for the docx
fig, axes = plt.subplots(1, 1, figsize=(10.5, 8.5), dpi=200)
df = pd.read_csv(OUT / "spearman_myeloid_crossmethod.tsv", sep=None, engine="python")
df.set_index("method", inplace=True)
ax = axes
im = ax.imshow(df.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(df.columns)))
ax.set_yticks(range(len(df.index)))
ax.set_xticklabels(df.columns, rotation=40, ha="right", fontsize=10)
ax.set_yticklabels(df.index, fontsize=10)
for i in range(df.shape[0]):
    for j in range(df.shape[1]):
        v = df.values[i, j]
        color = "white" if abs(v) > 0.6 else "black"
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", color=color, fontsize=9, weight="bold")
cbar = fig.colorbar(im, ax=ax, shrink=0.85)
cbar.set_label("Spearman ρ", fontsize=11)
ax.set_title("Myeloid signal — cross-method Spearman correlation (n=702)\n"
             "Three brain-tuned microglia signatures cluster tightly (r > 0.98);\n"
             "MDM signatures separate (r ≈ 0.40 vs. microglia); xCell macrophage tracks MDM (r = 0.86)",
             fontsize=11, weight="bold", pad=15)

# Add boxes to highlight microglia and MDM clusters
import matplotlib.patches as patches
# Microglia block: Microglia_Klemm, MicrogliaCore_Bohlen, MgTAM_Antunes (rows 4,5,7 in current order)
labels = list(df.index)
mg_idx = [i for i, n in enumerate(labels) if n in ["Microglia_Klemm", "MicrogliaCore_Bohlen", "MgTAM_Antunes"]]
mdm_idx = [i for i, n in enumerate(labels) if n in ["MDM_Klemm", "MoTAM_Antunes"]]
if len(mg_idx) >= 2:
    rect = patches.Rectangle((min(mg_idx)-0.5, min(mg_idx)-0.5),
                             max(mg_idx)-min(mg_idx)+1, max(mg_idx)-min(mg_idx)+1,
                             linewidth=3, edgecolor="#2E5C8A", facecolor="none", zorder=10)
    ax.add_patch(rect)
    ax.annotate("Microglia\ncluster", xy=(max(mg_idx)+0.6, (min(mg_idx)+max(mg_idx))/2),
                fontsize=10, color="#2E5C8A", weight="bold", va="center")
if len(mdm_idx) >= 2:
    rect = patches.Rectangle((min(mdm_idx)-0.5, min(mdm_idx)-0.5),
                             max(mdm_idx)-min(mdm_idx)+1, max(mdm_idx)-min(mdm_idx)+1,
                             linewidth=3, edgecolor="#C0392B", facecolor="none", zorder=10)
    ax.add_patch(rect)
    ax.annotate("MDM\ncluster", xy=(max(mdm_idx)+0.6, (min(mdm_idx)+max(mdm_idx))/2),
                fontsize=10, color="#C0392B", weight="bold", va="center")

plt.tight_layout()
plt.savefig(OUT / "fig_corr_myeloid_annotated.png", dpi=200, bbox_inches="tight", facecolor="white")
print("Wrote fig_corr_myeloid_annotated.png")
