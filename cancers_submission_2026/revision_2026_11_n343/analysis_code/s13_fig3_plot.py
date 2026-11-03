"""Plot revised cluster diagnostics and state visualisation from n=344 data."""

from pathlib import Path

import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import load, ROOT as QR, SHORT, COL, HI, MID, LO, canonicalize_figure
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
import umap

E = Path("/mnt/user-data/uploads/Open PBTA")
OUT = Path(QR) / "out/fig3"; FIG = Path(QR) / "fig/Figure3"; FIG.mkdir(parents=True, exist_ok=True)
sw = pd.read_csv(OUT / "ksweep_n343.tsv", sep="\t")
gap = pd.read_csv(OUT / "gap_n343.tsv", sep="\t")
obs = load()
features = pd.read_csv(E / "Revision/Week1/_inputs/step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t", index_col=0)
ids = [s for s in features.index if s in obs.index]
X = features.loc[ids].copy()
assert X.shape == (343, 46) and X.isna().sum().sum() == 0
state = obs.loc[ids, "state"].map(SHORT).to_numpy()
palette = {SHORT[k]: COL[k] for k in (HI, MID, LO)}
mpl.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})

pca = PCA(n_components=2).fit_transform(X)
embedding = umap.UMAP(n_neighbors=15, min_dist=.3, random_state=42).fit_transform(X)
proj = pd.DataFrame({"sample": ids, "state": state, "PC1": pca[:,0], "PC2": pca[:,1],
                     "UMAP1": embedding[:,0], "UMAP2": embedding[:,1]})
proj.to_csv(OUT / "projections_n343.tsv", sep="\t", index=False)
means = X.groupby(pd.Series(state, index=X.index)).mean().loc[["Low", "Intermediate", "High"]]
rank = (means.max(axis=0) - means.min(axis=0)).sort_values(ascending=False)
selected = rank.head(20).index
means[selected].T.to_csv(OUT / "top20_feature_means_n343.tsv", sep="\t")

fig = plt.figure(figsize=(12.8, 9.5), dpi=300)
outer = fig.add_gridspec(3, 1, height_ratios=[1.25, 1.55, 1.6], hspace=.48)
top = outer[0].subgridspec(1, 5, wspace=.42)
colors = ["#2563EB", "#F97316", "#94A3B8", "#0F766E", "#B91C1C"]
specs = [
    ("A", "Ambiguous clustering", "PAC (lower is better)", sw.k, sw.PAC),
    ("B", "Silhouette width", "Silhouette on 1 − consensus", sw.k, sw.silhouette_consensus),
    ("C", "Consensus CDF", "Relative ΔAUC", sw.k, sw.delta_AUC),
    ("D", "Gap statistic", "Gap statistic", gap.k, gap.gap),
    ("E", "Minimum cluster size", "Smallest cluster (n)", sw.k, sw.min_cluster_size),
]
for i, (letter, title, ylab, xx, yy) in enumerate(specs):
    ax = fig.add_subplot(top[i])
    if i == 2:
        ax.bar(xx, yy, color=["#B91C1C" if k == 3 else "#94A3B8" for k in xx], width=.65)
    elif i == 3:
        ax.errorbar(xx, yy, yerr=gap.s_k, fmt="o-", color=colors[i], ms=3, lw=1, capsize=2)
    else:
        ax.plot(xx, yy, "o-", color=colors[i], ms=3, lw=1)
    ax.axvline(3, color="#B91C1C", ls="--", lw=.8)
    if i == 4:
        ax.set_yscale("log")
        ax.axhline(10, color="#64748B", ls=":", lw=.8)
    ax.set_title(f"{letter}  {title}", loc="left", fontsize=8, weight="bold")
    ax.set_xlabel("Number of clusters (k)", fontsize=6.8)
    ax.set_ylabel(ylab, fontsize=6.8)
    ax.set_xticks(range(1, 11))
    ax.tick_params(labelsize=6)
    ax.spines[["top", "right"]].set_visible(False)

mid = outer[1].subgridspec(1, 2, wspace=.18)
for j, (xcol, ycol, title) in enumerate([("PC1", "PC2", "PCA"), ("UMAP1", "UMAP2", "UMAP")]):
    ax = fig.add_subplot(mid[j])
    for s in ["High", "Intermediate", "Low"]:
        sub = proj[proj.state == s]
        ax.scatter(sub[xcol], sub[ycol], s=10, alpha=.8, lw=0, color=palette[s], label=f"{s} (n={len(sub)})")
    ax.set_title(("F  " if j == 0 else "") + title + " of the 46-feature matrix", loc="left", fontsize=9, weight="bold")
    ax.set_xlabel(xcol, fontsize=8); ax.set_ylabel(ycol, fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    if j == 0: ax.legend(frameon=False, fontsize=7, ncol=3, loc="lower left", title="Immune state", title_fontsize=7)

ax = fig.add_subplot(outer[2])
sns.heatmap(means[selected].T, cmap="RdBu_r", center=0, cbar_kws={"label": "Mean input z score", "shrink": .75},
            linewidths=.25, linecolor="white", ax=ax)
ax.set_title("G  Top 20 discriminating input features by immune state", loc="left", fontsize=9, weight="bold")
ax.set_xlabel("Immune state", fontsize=8); ax.set_ylabel("")
ax.set_xticklabels(["Immune-low", "Immune-intermediate", "Immune-high"], rotation=0, fontsize=7)
ax.set_yticks(np.arange(len(selected)) + .5)
ax.set_yticklabels([name.replace("ssGSEA_", "").replace("LM22_", "") for name in selected], fontsize=6.6)
fig.subplots_adjust(left=.18, right=.965, top=.97, bottom=.06)
mpl.rcParams["svg.fonttype"]="none"
for ext in ("pdf","svg","png"): canonicalize_figure(fig); fig.savefig(FIG / f"Figure3.{ext}", dpi=300 if ext=="png" else None, bbox_inches="tight", facecolor="white")
