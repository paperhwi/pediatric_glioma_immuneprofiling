"""Regenerate PCA/UMAP plots with fixed ecotype labels and build full Step 5 notebook source."""
from pathlib import Path
import json, warnings
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_step5"; FIG.mkdir(exist_ok=True)

# Reload PCA/UMAP coords saved earlier and patch with fixed ecotype labels
ann3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
ann2 = pd.read_csv(OUT / "ecotype_LM22_main_k2_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")

pca = pd.read_csv(OUT / "ecotype_LM22_main_pca.tsv", sep="\t").set_index("sample")
pca["ecotype_k3"] = ann3.reindex(pca.index)["ecotype"]
pca["ecotype_k2"] = ann2.reindex(pca.index)["ecotype"]
pca.reset_index().to_csv(OUT / "ecotype_LM22_main_pca.tsv", sep="\t", index=False)

try:
    um = pd.read_csv(OUT / "ecotype_LM22_main_umap.tsv", sep="\t").set_index("sample")
    um["ecotype_k3"] = ann3.reindex(um.index)["ecotype"]
    um["ecotype_k2"] = ann2.reindex(um.index)["ecotype"]
    um.reset_index().to_csv(OUT / "ecotype_LM22_main_umap.tsv", sep="\t", index=False)
except FileNotFoundError:
    um = None

ec_colors_k3 = {"Lymphocyte-inflamed": "#3B82F6", "Myeloid-dominant": "#EF4444", "Immune-desert": "#9CA3AF"}
ec_colors_k2 = {"Immune-active": "#16A34A", "Immune-desert": "#9CA3AF"}

# 4-panel grid: PCA k=2, PCA k=3, UMAP k=2, UMAP k=3
nrow = 2 if um is not None else 1
fig, axes = plt.subplots(nrow, 2, figsize=(11, 4.5 * nrow))
if nrow == 1:
    axes = np.array([axes])

def scatter(ax, df, x, y, color_col, palette, title):
    for ec, sub in df.groupby(color_col):
        ax.scatter(sub[x], sub[y], c=palette.get(ec, "#888"), label=ec, s=18, alpha=.75, edgecolor="none")
    ax.set_xlabel(x); ax.set_ylabel(y); ax.set_title(title)
    ax.legend(loc="best", frameon=False, fontsize=8)

scatter(axes[0, 0], pca, "PC1", "PC2", "ecotype_k2", ec_colors_k2, "PCA — k=2 (PAC-optimal)")
scatter(axes[0, 1], pca, "PC1", "PC2", "ecotype_k3", ec_colors_k3, "PCA — k=3 (pre-specified)")
if um is not None:
    scatter(axes[1, 0], um.reset_index(), "UMAP1", "UMAP2", "ecotype_k2", ec_colors_k2, "UMAP — k=2")
    scatter(axes[1, 1], um.reset_index(), "UMAP1", "UMAP2", "ecotype_k3", ec_colors_k3, "UMAP — k=3")
plt.tight_layout()
fig.savefig(FIG / "ecotype_PCA_UMAP.png", dpi=200)
plt.close()

# Stacked bar: ecotype distribution by cohort_group for k=3
ct_k3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_x_cohort_group.tsv", sep="\t").set_index("ecotype")
ct_pct = ct_k3.div(ct_k3.sum(axis=0), axis=1) * 100
fig, ax = plt.subplots(figsize=(6.5, 4.2))
bottom = np.zeros(len(ct_pct.columns))
for ec in ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]:
    if ec in ct_pct.index:
        vals = ct_pct.loc[ec].values
        ax.bar(ct_pct.columns, vals, bottom=bottom, label=ec, color=ec_colors_k3[ec])
        bottom += vals
ax.set_ylabel("% of samples within cohort"); ax.set_ylim(0, 100)
ax.set_title("Ecotype composition by cohort group (k=3)")
ax.legend(loc="upper right", fontsize=8, frameon=False)
plt.tight_layout()
fig.savefig(FIG / "ecotype_k3_by_cohort_group.png", dpi=200)
plt.close()

# Heatmap of k=3 cluster feature means
means = pd.read_csv(OUT / "ecotype_LM22_main_k3_feature_means.tsv", sep="\t", index_col=0)
# Map cluster idx → ecotype name
map3 = {2: "Lymphocyte-inflamed", 3: "Myeloid-dominant", 1: "Immune-desert"}
means.index = [f"{map3.get(int(i), int(i))} (c{int(i)})" for i in means.index]
# Show top 30 features by absolute spread
spread = means.max() - means.min()
keep = spread.sort_values(ascending=False).head(30).index
fig, ax = plt.subplots(figsize=(8, 9))
import seaborn as sns
sns.heatmap(means[keep].T, cmap="RdBu_r", center=0, annot=True, fmt=".2f", ax=ax,
            cbar_kws=dict(label="z-score (mean per cluster)"))
ax.set_title("k=3 ecotype × top-30 discriminating features")
plt.tight_layout()
fig.savefig(FIG / "ecotype_k3_feature_means_heatmap.png", dpi=200)
plt.close()

print("plots regenerated:")
for p in FIG.glob("*.png"):
    print(" ", p.relative_to(BASE))
