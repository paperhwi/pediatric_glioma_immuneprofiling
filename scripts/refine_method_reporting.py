"""
Methodology refinement:
  - LM22 QC reported on Main n=349 (primary) + Full 702 + BRAF 353 as sub-pool table (Suppl)
  - New Fig 2D: 4-method consensus on key immune axes × 3 ecotypes
  - Rebuild Fig 2A (LM22 QC) using Main n=349 subset
  - Re-order Suppl Table S1 by family grouping
"""
from __future__ import annotations
from pathlib import Path
import warnings, json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
from scipy import stats
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; DST = OUT / "figs_300dpi"; DST.mkdir(exist_ok=True)
DPI  = 300

PAL_EC = {"Lymphocyte-inflamed":"#3B82F6","Myeloid-dominant":"#EF4444","Immune-desert":"#9CA3AF"}
EC_ORDER = ["Lymphocyte-inflamed","Myeloid-dominant","Immune-desert"]

ec3   = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
braf  = pd.read_csv(OUT / "step10_BRAF_ALT_ecotype_assignment.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
qc    = pd.read_csv(OUT / "CIBERSORTx_Job15_Results.txt", sep="\t").rename(columns={"Mixture":"sample"}).set_index("sample")
merged = pd.read_csv(OUT / "step4_integrated_feature_matrix.tsv", sep="\t", index_col=0)

# ---------------------------------------------------------------------------
# 1) Sub-pool LM22 QC table → SupplTableS21
# ---------------------------------------------------------------------------
print("[1] sub-pool LM22 QC table")
rows = []
for name, ids in [("Full pool (Main + BRAF_ALT)",  qc.index),
                  ("Main cohort",                  ec3.index),
                  ("BRAF_ALT (HGG/PXA + LGG)",     braf.index)]:
    sub = qc.loc[qc.index.intersection(ids)]
    rows.append({
        "sub_pool":              name,
        "n":                     int(len(sub)),
        "median_P":              float(sub["P-value"].median()),
        "median_correlation":    float(sub["Correlation"].median()),
        "median_RMSE":           float(sub["RMSE"].median()),
        "n_P_lt_0.05":           int((sub["P-value"] < 0.05).sum()),
        "pct_P_lt_0.05":         float((sub["P-value"] < 0.05).mean() * 100),
    })
sub_pool_qc = pd.DataFrame(rows)
sub_pool_qc.to_csv(OUT / "supplTableS21_LM22_QC_by_subpool.tsv", sep="\t", index=False)
print(sub_pool_qc.to_string(index=False))

# ---------------------------------------------------------------------------
# 2) Fig 2A — LM22 QC distribution, Main n=349 (rebuild)
# ---------------------------------------------------------------------------
print("\n[2] Fig 2A — Main n=349 LM22 QC")
qc_main = qc.loc[qc.index.intersection(ec3.index)]
fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
for ax, col, lab in zip(axes, ["P-value", "Correlation", "RMSE"], ["P-value", "Correlation", "RMSE"]):
    sns.histplot(qc_main[col], bins=40, ax=ax, color="#7C3AED")
    ax.set_title(f"{lab}  (Main cohort, n={len(qc_main)})")
    ax.set_xlabel(lab)
    if col == "P-value":
        ax.axvline(0.05, color="red", ls="--", alpha=.7, label="0.05"); ax.legend()
plt.tight_layout(); fig.savefig(DST / "Fig2A_LM22_QC.png", dpi=DPI); plt.close()

# Supplementary version on Full pool for completeness
fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
for ax, col, lab in zip(axes, ["P-value", "Correlation", "RMSE"], ["P-value", "Correlation", "RMSE"]):
    sns.histplot(qc[col], bins=40, ax=ax, color="#A78BFA")
    ax.set_title(f"{lab}  (full pool, n={len(qc)})")
    if col == "P-value":
        ax.axvline(0.05, color="red", ls="--", alpha=.7, label="0.05"); ax.legend()
plt.tight_layout(); fig.savefig(DST / "SupplFig_S5_LM22_QC_fullpool.png", dpi=DPI); plt.close()

# ---------------------------------------------------------------------------
# 3) Fig 2D — 4-method consensus on key axes × 3 ecotypes
# ---------------------------------------------------------------------------
print("\n[3] Fig 2D — 4-method × ecotype consensus")
df = merged.loc[merged.index.intersection(ec3.index)].copy()
df["ecotype"] = ec3.loc[df.index, "ecotype"]

# Build a tall table: axis × method → series for each ecotype
axis_map = {
    "CD8 T cell": {
        "LM22":      "LM22_T cells CD8",
        "xCell":     "xCell_T cell CD8+",
        "quanTIseq": "quanTIseq_T cell CD8+",
        "MCP":       "MCP_T cell CD8+",
        "EPIC":      "EPIC_T cell CD8+",
    },
    "NK cell": {
        "xCell":     "xCell_NK cell",
        "quanTIseq": "quanTIseq_NK cell",
        "MCP":       "MCP_NK cell",
        # LM22: combine NK resting + NK activated (constructed below)
    },
    "Macrophage": {
        "xCell":     "xCell_Macrophage",
        "MCP":       "MCP_Macrophage/Monocyte",
        "EPIC":      "EPIC_Macrophage",
        # LM22: M0+M1+M2 sum (constructed below)
    },
    "Monocyte / DC": {
        "xCell":     "xCell_Monocyte",
        "MCP":       "MCP_Monocytic lineage",
        "MCP_DC":    "MCP_Myeloid dendritic cell",
        # LM22: Monocytes + DC resting + DC activated
    },
}

# Add LM22 composite columns
df["LM22_NK_total"]                = df[["LM22_NK cells resting","LM22_NK cells activated"]].sum(axis=1)
df["LM22_Macrophages_total"]       = df[["LM22_Macrophages M0","LM22_Macrophages M1","LM22_Macrophages M2"]].sum(axis=1)
df["LM22_Monocyte_DC_total"]       = df[["LM22_Monocytes","LM22_Dendritic cells resting","LM22_Dendritic cells activated"]].sum(axis=1)
axis_map["NK cell"]["LM22"]           = "LM22_NK_total"
axis_map["Macrophage"]["LM22"]        = "LM22_Macrophages_total"
axis_map["Monocyte / DC"]["LM22"]     = "LM22_Monocyte_DC_total"

# z-score each method column within Main cohort, so they're comparable on same scale
def z(s):
    sd = s.std(ddof=0)
    return (s - s.mean()) / (sd if sd > 0 else 1.0)

axis_dfs = {}
for axis, methods in axis_map.items():
    long = []
    for m, col in methods.items():
        if col in df.columns:
            zs = z(df[col]).dropna()
            for s, val in zs.items():
                ec = df.loc[s, "ecotype"]
                long.append({"axis": axis, "method": m, "sample": s, "z": val, "ecotype": ec})
    axis_dfs[axis] = pd.DataFrame(long)

# Save raw long table as supplementary
all_long = pd.concat([d.assign(axis=k) for k, d in axis_dfs.items()], ignore_index=True)
all_long.to_csv(OUT / "supplTableS22_4method_axis_consensus_long.tsv", sep="\t", index=False)

# Plot 4 axes × 1 row, each axis showing 4-5 methods × 3 ecotypes as grouped boxplot
fig, axes = plt.subplots(1, 4, figsize=(18, 4.6), sharey=False)
for ax, (axis, d) in zip(axes, axis_dfs.items()):
    methods_in_axis = list(axis_map[axis].keys())
    # Ensure consistent method order with LM22 first
    method_order = [m for m in ["LM22", "xCell", "quanTIseq", "MCP", "MCP_DC", "EPIC"] if m in methods_in_axis]
    sns.boxplot(data=d, x="method", y="z", hue="ecotype",
                order=method_order, hue_order=EC_ORDER, palette=PAL_EC,
                ax=ax, showfliers=False, fliersize=0)
    ax.set_title(axis, fontsize=11)
    ax.set_xlabel(""); ax.set_ylabel("z-score (within Main cohort)")
    ax.axhline(0, color="grey", lw=0.6, alpha=.5)
    ax.tick_params(axis="x", rotation=20, labelsize=9)
    if ax is not axes[-1]:
        if ax.get_legend() is not None: ax.get_legend().remove()
    else:
        ax.legend(loc="upper right", fontsize=8, frameon=False)
fig.suptitle("Fig 2D.  Cross-method consensus on key immune axes by ecotype  (Main n=349)", fontsize=12, y=1.02)
plt.tight_layout(); fig.savefig(DST / "Fig2D_4method_by_ecotype.png", dpi=DPI, bbox_inches="tight"); plt.close()

# Also produce a single summary heatmap: rows = method × axis, cols = ecotype, values = mean z
heat_rows = []
for axis, d in axis_dfs.items():
    for m in d["method"].unique():
        for ec in EC_ORDER:
            v = d[(d["method"]==m) & (d["ecotype"]==ec)]["z"].mean()
            heat_rows.append({"axis_method": f"{axis} · {m}", "ecotype": ec, "mean_z": v})
heatdf = pd.DataFrame(heat_rows).pivot(index="axis_method", columns="ecotype", values="mean_z").reindex(columns=EC_ORDER)
# Order rows by axis then method
desired_method_order = ["LM22", "xCell", "quanTIseq", "MCP", "MCP_DC", "EPIC"]
order = []
for ax_name in axis_map:
    for m in desired_method_order:
        key = f"{ax_name} · {m}"
        if key in heatdf.index: order.append(key)
heatdf = heatdf.reindex(order)

fig, ax = plt.subplots(figsize=(6.0, max(4, 0.4*len(heatdf))))
sns.heatmap(heatdf, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
            cbar_kws=dict(label="mean z within ecotype"), ax=ax,
            linewidths=0.4, linecolor="white")
# Draw separators between axis groups
running = 0
axis_sizes = []
for ax_name in axis_map:
    cnt = sum(1 for r in heatdf.index if r.startswith(ax_name + " · "))
    running += cnt; axis_sizes.append(running)
for pos in axis_sizes[:-1]:
    ax.axhline(pos, color="black", lw=1.0)
ax.set_title("Cross-method consensus heatmap by ecotype  (Main n=349)")
ax.set_xlabel(""); ax.set_ylabel("")
plt.tight_layout(); fig.savefig(DST / "Fig2D_4method_heatmap.png", dpi=DPI, bbox_inches="tight"); plt.close()

# ---------------------------------------------------------------------------
# 4) Suppl Table S1 — re-order by family + add cumulative q-rank
# ---------------------------------------------------------------------------
print("\n[4] Suppl Table S1 re-order by family")
S1 = pd.read_csv(OUT / "step6_KW_results.tsv", sep="\t")
family_order = ["LM22", "ssGSEA", "xCell", "quanTIseq", "MCP", "EPIC"]
S1["family"] = pd.Categorical(S1["family"], categories=family_order, ordered=True)
S1 = S1.sort_values(["family", "q_BH"]).reset_index(drop=True)
S1.to_csv(OUT / "supplTableS1_KW_per_feature_by_family.tsv", sep="\t", index=False)

# Family summary
fam = (S1
       .groupby("family", observed=True)
       .agg(n_features=("feature", "count"),
            n_sig_q05=("q_BH", lambda s: int((s < 0.05).sum())),
            min_q=("q_BH", "min"),
            max_epsilon_sq=("epsilon_sq", "max"))
       .reset_index())
fam["pct_sig"] = (fam["n_sig_q05"] / fam["n_features"] * 100).round(1)
fam.to_csv(OUT / "supplTableS1_family_summary.tsv", sep="\t", index=False)
print(fam.to_string(index=False))

# ---------------------------------------------------------------------------
# 5) Summary JSON
# ---------------------------------------------------------------------------
summary = {
    "sub_pool_QC_table_TSV": "supplTableS21_LM22_QC_by_subpool.tsv",
    "Main_cohort_LM22": {
        "n": 349,
        "median_P": float(qc_main["P-value"].median()),
        "median_correlation": float(qc_main["Correlation"].median()),
        "median_RMSE": float(qc_main["RMSE"].median()),
        "n_P_lt_0.05": int((qc_main["P-value"] < 0.05).sum()),
        "pct_P_lt_0.05": float((qc_main["P-value"] < 0.05).mean() * 100),
    },
    "new_figures_300dpi": [
        "Fig2A_LM22_QC.png  (rebuilt with Main n=349)",
        "Fig2D_4method_by_ecotype.png",
        "Fig2D_4method_heatmap.png",
        "SupplFig_S5_LM22_QC_fullpool.png",
    ],
    "updated_suppl_tables": [
        "supplTableS1_KW_per_feature_by_family.tsv",
        "supplTableS1_family_summary.tsv",
        "supplTableS21_LM22_QC_by_subpool.tsv",
        "supplTableS22_4method_axis_consensus_long.tsv",
    ],
}
(OUT / "method_refinement_summary.json").write_text(json.dumps(summary, indent=2))
print("\n" + json.dumps(summary, indent=2))
