"""Recalculate Figure 2 from original per-specimen outputs on the n=344 set."""

from pathlib import Path
import json

import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import load, HI, MID, LO, SHORT, COL, ROOT as QR
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

E = Path("/mnt/user-data/uploads/Open PBTA")
OUT = Path(QR) / "fig/Figure2"
OUT.mkdir(exist_ok=True)
mpl.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})

obs = load()
ids = set(obs.index)
assert len(ids) == 343
state = obs.state.map(SHORT).to_dict()

qc = pd.read_csv("/mnt/user-data/uploads/Pediatric GBM/Final_n343_3Nov/_inputs/CIBERSORTx_Job15_Results.csv")
qc = qc[qc.Mixture.isin(ids)].copy()
assert len(qc) == 343 and qc.Mixture.is_unique

long = pd.read_csv(E / "1차 과업 수행본 - 김종휘/02_tables_final/SupplTableS22_4method_axis_consensus_long.tsv", sep="\t")
assert long["sample"].nunique() == 349
wide_old = long.pivot(index=["axis", "sample"], columns="method", values="z").reset_index()
wide = wide_old[wide_old["sample"].isin(ids)].copy()
assert wide["sample"].nunique() == 343

peer = {
    "CD8 T cell": ["xCell", "quanTIseq", "MCP", "EPIC"],
    "NK cell": ["xCell", "quanTIseq", "MCP"],
    "Macrophage": ["xCell", "MCP", "EPIC"],
    "Monocyte / DC": ["xCell", "MCP_DC"],
}
old_table = pd.read_csv(E / "1차 과업 수행본 - 김종휘/09_revision2/Table_item4_Figure2B_n349_vs_n702.tsv", sep="\t")
rows = []
max_old_diff = 0.0
for axis, methods in peer.items():
    old = wide_old[wide_old.axis == axis]
    cur = wide[wide.axis == axis]
    for method in methods:
        rho_old, _ = spearmanr(old.LM22, old[method])
        target = old_table.loc[(old_table.axis == axis) & (old_table.peer_method == method), "rho_n349"]
        assert len(target) == 1
        max_old_diff = max(max_old_diff, abs(rho_old - target.iloc[0]))
        rho, p = spearmanr(cur.LM22, cur[method])
        rows.append({"axis": axis, "peer_method": method, "n": len(cur), "spearman_rho": rho, "p": p})
assert max_old_diff < 1e-12, max_old_diff
res = pd.DataFrame(rows)
res["q_BH"] = multipletests(res.p, method="fdr_bh")[1]
res.to_csv(OUT / "T_fig2B_cross_method_n343.tsv", sep="\t", index=False)

long = long[long["sample"].isin(ids)].copy()
long["state"] = long["sample"].map(state)
assert not long.state.isna().any()
# The source S22 table is z-scored on n=349. Recenter and rescale within the
# n=343 cohort for the displayed cross-method boxplots; ranks are unchanged.
long["z343"] = long.groupby(["axis", "method"]).z.transform(lambda s: (s - s.mean()) / s.std(ddof=0))
long.to_csv(OUT / "T_fig2C_axis_estimates_n343_long.tsv", sep="\t", index=False)

palette = {SHORT[k]: COL[k] for k in (HI, MID, LO)}
method_color = {"xCell": "#2563EB", "quanTIseq": "#059669", "MCP": "#D97706", "EPIC": "#DC2626", "MCP_DC": "#7C3AED"}
fig = plt.figure(figsize=(12.4, 10.5), dpi=300)
outer = fig.add_gridspec(3, 1, height_ratios=[1.25, 1.5, 1.5], hspace=.45)

top = outer[0].subgridspec(1, 3, wspace=.24)
for j, field in enumerate(["P-value", "Correlation", "RMSE"]):
    ax = fig.add_subplot(top[j])
    ax.hist(qc[field], bins=30, color="#8B5CF6", edgecolor="#453078", linewidth=.3)
    if field == "P-value":
        ax.axvline(.05, color="#EF4444", ls="--", lw=1)
    ax.set_title(field, fontsize=10)
    ax.set_xlabel(field, fontsize=8)
    if j == 0:
        ax.set_ylabel("Specimens", fontsize=8)
        ax.text(-.12, 1.13, "A", transform=ax.transAxes, fontsize=15, weight="bold")
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)

ax = fig.add_subplot(outer[1])
rr = res.copy()
rr["label"] = rr.axis.str.replace(" / ", "/", regex=False) + "\n" + rr.peer_method
x = np.arange(len(rr))
ax.bar(x, rr.spearman_rho, color=[method_color[m] for m in rr.peer_method], width=.62)
for i, row in rr.iterrows():
    mark = "***" if row.q_BH < .001 else "**" if row.q_BH < .01 else "*" if row.q_BH < .05 else "ns"
    ax.text(i, max(row.spearman_rho, 0) + .012, mark, ha="center", fontsize=7)
ax.axhline(0, color="#555", lw=.8)
ax.set_xticks(x, rr.label, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Spearman ρ versus LM22", fontsize=9)
ax.set_title("Cross-method concordance in the revised cohort (n = 343)", fontsize=10)
ax.text(-.045, 1.05, "B", transform=ax.transAxes, fontsize=15, weight="bold")
ax.spines[["top", "right"]].set_visible(False)

bottom = outer[2].subgridspec(1, 4, wspace=.22)
for j, axis in enumerate(peer):
    ax = fig.add_subplot(bottom[j])
    methods = ["LM22", *peer[axis]]
    sub = long[long.axis == axis]
    for k, method in enumerate(methods):
        for s, offset in [("High", -.25), ("Intermediate", 0), ("Low", .25)]:
            vals = sub.loc[(sub.method == method) & (sub.state == s), "z343"].to_numpy()
            if len(vals):
                bp = ax.boxplot(vals, positions=[k + offset], widths=.22, patch_artist=True, showfliers=False)
                bp["boxes"][0].set(facecolor=palette[s], edgecolor="#5B6470", linewidth=.5)
                for item in bp["whiskers"] + bp["caps"] + bp["medians"]:
                    item.set(color="#5B6470", linewidth=.65)
    ax.axhline(0, color="#777", lw=.5)
    ax.set_xticks(range(len(methods)), methods, rotation=45, ha="right", fontsize=6.5)
    ax.set_title(axis, fontsize=8.5)
    ax.tick_params(labelsize=6.5)
    if j == 0:
        ax.set_ylabel("Within-cohort z score", fontsize=8)
        ax.text(-.18, 1.07, "C", transform=ax.transAxes, fontsize=15, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
fig.subplots_adjust(left=.075, right=.985, top=.925, bottom=.075)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=palette[k], label=k) for k in ("High","Intermediate","Low")], title="Immune state", loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(.5, -.02))
mpl.rcParams["svg.fonttype"] = "none"
for ext in ("pdf","svg","png"): fig.savefig(OUT / f"Figure2.{ext}", dpi=300 if ext=="png" else None, bbox_inches="tight", facecolor="white")
res.to_csv(Path(QR)/"out/T_fig2B_cross_method_n343.tsv", sep="\t", index=False)
summary = {"n": 343, "qc_p_below_005": int((qc["P-value"] < .05).sum()),
           "qc_pct_below_005": float(100 * (qc["P-value"] < .05).mean()),
           "old_correlation_max_abs_diff": float(max_old_diff), "comparisons": len(res)}
(OUT / "Figure2_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
