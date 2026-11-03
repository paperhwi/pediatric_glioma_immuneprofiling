"""Figure 5 (n = 343): A external single-cell AUROC (independent of the bulk cohort; stored AUROC summary, vector rendering); B cell class with the highest
mean expression of the top 300 Immune-high vs Immune-low genes in ScPCA SCPCP000001; C-E Visium spots by nearest bulk-state centroid."""
import json, sys
from pathlib import Path
sys.path.insert(0,"/home/claude/q343/src"); from common import ROOT as QR, SHORT, COL, HI, MID, LO, canonicalize_figure
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr
OUT = Path(QR) / "fig/Figure5"; OUT.mkdir(parents=True, exist_ok=True)
mpl.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
designs = pd.read_csv(QR + "out/fig5b/T_fig5B_designs.tsv", sep="\t")
record = designs.loc[designs.gene_set.str.startswith("top 300 Immune-high") & designs.design.str.startswith("D4 all ten")].iloc[0]
counts = {}
for part in record.top_counts.split("; "):
    name, count = part.rsplit(" ", 1); counts[name] = int(count)
assert sum(counts.values()) == 300
data = pd.read_csv(QR + "out/T_visium_spots_n1650.tsv", sep="\t", index_col=0)
data["spot_class"] = data.spot_class.map(SHORT)
sections = ["Sample-1", "Sample-2", "Sample-3"]
states = [LO, MID, HI]
colors = {SHORT[k]: COL[k] for k in (HI, MID, LO)}
stats = []
for sec in sections:
    sub = data[data.section == sec]
    rho = float(spearmanr(sub.axis_score, sub.immune_fraction).statistic)
    mixed = float((sub.n_types_ge10pct >= 2).mean())
    row = {"section": sec, "n": len(sub), "axis_immune_rho": rho, "mixed_two_or_more_pct": 100*mixed}
    row.update({f"median_immune_{state}": float(sub.loc[sub.spot_class == state, "immune_fraction"].median()) for state in states})
    stats.append(row)
summary = {"gene_set_n": 300, "top_immune_n": int(record.n_immune_top),
           "top_immune_pct": float(record.pct_immune_top),
           "random_matched_mean_pct": float(record.random_mean_pct),
           "random_matched_95pct": [float(record.random_lo), float(record.random_hi)],
           "permutation_p": float(record.perm_p), "spots": stats}
(Path(QR) / "out/T_fig5_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
data[["section", "spot_class", "axis_score", "immune_fraction", "n_types_ge10pct", "x", "y"]].to_csv(
    OUT / "Figure5CDE_source_spots_n1650.tsv", sep="\t")

fig = plt.figure(figsize=(12.5, 9.9), dpi=300)
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.25, 1.4, 1.25], hspace=.55, wspace=.25,
              left=.065, right=.985, top=.97, bottom=.065)

# A: source panel is external-cell signature AUROC, independent of the bulk
# cohort and already checked in the submitted manuscript.
ax = fig.add_subplot(gs[0, 0])
auc = pd.read_csv(Path(__file__).resolve().parents[1] / "reviewer_followup_20261103/tables/Figure5A_source_AUROC.tsv", sep="\t")
auc = auc.sort_values("AUROC_all_cells")
yy = np.arange(len(auc))
bar_colors = ["#B91C1C" if verdict == "DISCORDANT" else "#6D28D9" if "T cell" in target else "#0F766E" for verdict,target in zip(auc.verdict,auc.expected)]
ax.barh(yy,auc.AUROC_all_cells,color=bar_colors,height=.72)
for yi,(_,rr) in zip(yy,auc.iterrows()):
    ax.plot([rr.AUROC_min_sample,rr.AUROC_max_sample],[yi,yi],color="#475569",lw=.7)
display_names = auc.signature.str.replace("_"," ").tolist()
display_names[0] = "Cytotoxic Lymphocyte Curated*"
ax.set_yticks(yy,display_names,fontsize=5.5)
ax.set_xlim(0,1.04);ax.axvline(.5,ls="--",color="#374151",lw=.7)
ax.set_xlabel("AUROC against the originally specified target lineage",fontsize=6)
ax.tick_params(axis="x",labelsize=6);ax.spines[["top","right"]].set_visible(False)
ax.text(-.05,1.05,"A",transform=ax.transAxes,fontsize=15,weight="bold")

# B: all-ten-type ScPCA comparison for the revised n=343 contrast.
bgs = gs[0, 1].subgridspec(1, 2, width_ratios=[1.5, .8], wspace=.32)
ax = fig.add_subplot(bgs[0])
names = list(counts)
vals = [counts[n] for n in names]
cell_colors = ["#0F766E" if n not in ("Endothelial", "Mural", "Oligodendrocyte", "Malignant") else "#94A3B8" for n in names]
ax.bar(range(len(names)), vals, color=cell_colors)
ax.set_xticks(range(len(names)), [n for n in names], rotation=53, ha="right", fontsize=6.3)
ax.set_ylabel("Top-expression genes (n)", fontsize=7)
ax.set_title("B  Cell class with highest mean expression", fontsize=8, loc="left", weight="bold")
ax.tick_params(labelsize=6.5)
ax.spines[["top", "right"]].set_visible(False)
ax = fig.add_subplot(bgs[1])
ax.bar([0,1], [record.pct_immune_top, record.random_mean_pct], color=["#0F766E", "#94A3B8"])
ax.errorbar([1], [record.random_mean_pct], yerr=[[record.random_mean_pct-record.random_lo],
                                             [record.random_hi-record.random_mean_pct]], fmt="none", color="#334155", capsize=3)
ax.set_ylim(0, 105)
ax.set_xticks([0,1], ["Top 300\ngenes", "Matched\nrandom"], fontsize=6.7)
ax.set_ylabel("Immune highest (%)", fontsize=7)
ax.text(0, 99, f"{record.pct_immune_top:.1f}%", ha="center", fontsize=7.2, color="#0F766E", weight="bold")
ax.text(1, record.random_mean_pct+6, f"{record.random_mean_pct:.1f}%", ha="center", fontsize=7)
ax.spines[["top", "right"]].set_visible(False)

# C: spatial classes, each spot coloured by its revised nearest-centroid class.
space = gs[1, :].subgridspec(1, 3, wspace=.13)
for j, sec in enumerate(sections):
    ax = fig.add_subplot(space[j])
    sub = data[data.section == sec]
    for state in [HI, MID, LO]:
        spot = sub[sub.spot_class == state]
        ax.scatter(spot.x, spot.y, s=5, color=colors[state], alpha=.82, lw=0, label=state)
    ax.set_aspect("equal", adjustable="datalim")
    ax.invert_yaxis(); ax.axis("off")
    ax.set_title(("C  " if j == 0 else "") + f"{sec} (n={len(sub)})", fontsize=8.5, weight="bold")
    if j == 0: ax.legend(frameon=False, fontsize=7, ncol=3, loc="lower left", bbox_to_anchor=(0, -.08), title="Nearest immune-state centroid", title_fontsize=7)

# D: distribution of estimated immune fraction by class and section.
ax = fig.add_subplot(gs[2, 0])
for j, sec in enumerate(sections):
    sub = data[data.section == sec]
    for k, state in enumerate(states):
        values = sub.loc[sub.spot_class == state, "immune_fraction"].to_numpy()
        pos = j*4 + k
        bp = ax.boxplot(values, positions=[pos], widths=.7, patch_artist=True, showfliers=False)
        bp["boxes"][0].set(facecolor=colors[state], edgecolor="#475569", linewidth=.6)
        for item in bp["whiskers"] + bp["caps"] + bp["medians"]:
            item.set(color="#475569", linewidth=.7)
ax.set_xticks([1,5,9], sections, fontsize=7.5)
ax.set_ylabel("Estimated immune fraction", fontsize=8)
ax.set_title("D  Spot immune fraction by nearest-centroid class (Low, Intermediate, High)", loc="left", fontsize=8.5, weight="bold")
ax.tick_params(labelsize=7)
ax.spines[["top", "right"]].set_visible(False)

# E: direct relationship between axis score and reference-based immune share.
ax = fig.add_subplot(gs[2, 1])
for j, sec in enumerate(sections):
    sub = data[data.section == sec]
    ax.scatter(sub.axis_score, sub.immune_fraction, s=5, alpha=.27,
               color=["#2563EB", "#D97706", "#0F766E"][j],
               label=f"{sec} (ρ={stats[j]['axis_immune_rho']:.2f}; mixed {stats[j]['mixed_two_or_more_pct']:.0f}%)")
ax.set_xlabel("Spot immune-axis score", fontsize=8)
ax.set_ylabel("Estimated immune fraction", fontsize=8)
ax.set_title("E  Spatial axis and mixed-cell profiles", loc="left", fontsize=8.5, weight="bold")
ax.legend(frameon=False, fontsize=6.6, loc="upper left")
ax.tick_params(labelsize=7)
ax.spines[["top", "right"]].set_visible(False)

for ext in ("pdf","svg","png"): canonicalize_figure(fig); fig.savefig(OUT / f"Figure5.{ext}", dpi=300 if ext=="png" else None, bbox_inches="tight", facecolor="white")
print(json.dumps(summary, indent=2))
