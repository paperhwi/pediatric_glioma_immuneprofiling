"""
Integrate CIBERSORTx LM22 results (Job #15) into the user-authored Manuscript_Draft_v1.

The user's manuscript uses:
  - Ecotype labels: Inflamed (n=113) / Intermediate (n=157) / Immune-desert (n=79)
  - Clustering feature set: quanTIseq + 19 immune-related ssGSEA (29 features)
  - Reference ecotype for Cox: Inflamed (HR shown as Immune-desert HR=1.79)
  - LM22 was originally "reserved for planned validation"

We will:
  1) Map our internal k=3 labels (Lymphocyte-inflamed/Myeloid-dominant/Immune-desert,
     111/160/78) onto the user's nomenclature.
  2) Compute LM22 QC and LM22 cell-fraction stratification by ecotype.
  3) Generate two new figures: LM22 QC (Main n=349) and LM22 fractions × ecotype.
  4) Re-render the cross-method correlation heatmap with LM22 added (5 methods).
  5) Modify the docx: update Methods §2.3, insert a new Results §3.2 paragraph after
     the "Brain-Tuned Signatures..." subsection, add Supplementary Figure S? caption.
"""
from __future__ import annotations
from pathlib import Path
import shutil, warnings, json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
from scipy import stats
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_lm22_integration"; FIG.mkdir(exist_ok=True)
DPI  = 300

# ----------------------------------------------------------------------------
# 1) Load data + map labels onto user's nomenclature
# ----------------------------------------------------------------------------
qc = pd.read_csv(OUT / "CIBERSORTx_Job15_Results.txt", sep="\t").rename(columns={"Mixture":"sample"}).set_index("sample")
lm22_cells = qc.drop(columns=["P-value", "Correlation", "RMSE", "Absolute score (sig.score)"])

ec3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
# Map internal labels -> user's manuscript labels
label_map = {
    "Lymphocyte-inflamed": "Inflamed",
    "Myeloid-dominant":    "Intermediate",
    "Immune-desert":       "Immune-desert",
}
ec3["ecotype_manuscript"] = ec3["ecotype"].map(label_map)
print("Ecotype mapping check:")
print(ec3["ecotype_manuscript"].value_counts())

# ----------------------------------------------------------------------------
# 2) LM22 QC on Main n=349
# ----------------------------------------------------------------------------
qc_main = qc.loc[qc.index.intersection(ec3.index)]
print(f"\nLM22 QC on Main n={len(qc_main)}:")
print(f"  median P-value = {qc_main['P-value'].median():.3f}")
print(f"  median Correlation = {qc_main['Correlation'].median():.3f}")
print(f"  median RMSE  = {qc_main['RMSE'].median():.3f}")
print(f"  P<0.05: {(qc_main['P-value']<0.05).sum()} ({(qc_main['P-value']<0.05).mean()*100:.1f}%)")

# Fig: LM22 QC distribution on Main
fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
for ax, col, label in zip(axes, ["P-value", "Correlation", "RMSE"], ["P-value", "Correlation", "RMSE"]):
    sns.histplot(qc_main[col], bins=40, ax=ax, color="#7C3AED")
    ax.set_title(f"LM22 {label} (Main n={len(qc_main)})")
    ax.set_xlabel(label)
    if col == "P-value":
        ax.axvline(0.05, color="red", ls="--", alpha=.7, label="0.05"); ax.legend()
plt.tight_layout(); fig.savefig(FIG / "LM22_QC_distribution_Main.png", dpi=DPI); plt.close()

# ----------------------------------------------------------------------------
# 3) LM22 cell-fraction stratification by ecotype (user's 3-group labels)
# ----------------------------------------------------------------------------
df = lm22_cells.loc[lm22_cells.index.intersection(ec3.index)].copy()
df["ecotype"] = ec3.loc[df.index, "ecotype_manuscript"]
ec_order = ["Inflamed", "Intermediate", "Immune-desert"]
palette  = {"Inflamed":"#3B82F6", "Intermediate":"#10B981", "Immune-desert":"#9CA3AF"}

# Renormalize fractions (LM22 is absolute mode; convert to relative for plotting)
cell_cols = [c for c in df.columns if c != "ecotype"]
df_rel = df[cell_cols].div(df[cell_cols].sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
df_rel["ecotype"] = df["ecotype"]

# Group LM22 cell types into biological super-categories for cleaner display
SUPERCAT = {
    "CD8 T cells":           ["T cells CD8"],
    "CD4 T cells":           ["T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated",
                              "T cells follicular helper"],
    "Tregs / γδ":            ["T cells regulatory (Tregs)", "T cells gamma delta"],
    "B cells":               ["B cells naive", "B cells memory", "Plasma cells"],
    "NK cells":              ["NK cells resting", "NK cells activated"],
    "Monocytes / Macrophages": ["Monocytes","Macrophages M0","Macrophages M1","Macrophages M2"],
    "Dendritic cells":       ["Dendritic cells resting","Dendritic cells activated"],
    "Granulocytes":          ["Mast cells resting","Mast cells activated","Eosinophils","Neutrophils"],
}
supercats = {}
for name, members in SUPERCAT.items():
    cols = [c for c in members if c in df_rel.columns]
    if cols: supercats[name] = df_rel[cols].sum(axis=1)
super_df = pd.DataFrame(supercats); super_df["ecotype"] = df_rel["ecotype"]
super_df.to_csv(OUT / "lm22_supercat_fractions_per_ecotype.tsv", sep="\t")

# KW + Dunn for each super-category
def bh(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p); ranked = p[order]
    q = ranked * n / (np.arange(n) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q,0,1); return out

stats_rows = []
for cat in supercats:
    groups = [super_df[super_df["ecotype"]==e][cat].dropna().values for e in ec_order]
    H, p = stats.kruskal(*groups)
    medians = [np.median(g) for g in groups]
    stats_rows.append({"category": cat, "H": float(H), "p": float(p),
                       "med_Inflamed": medians[0], "med_Intermediate": medians[1], "med_Desert": medians[2]})
sdf = pd.DataFrame(stats_rows)
sdf["q_BH"] = bh(sdf["p"].values)
sdf.to_csv(OUT / "lm22_supercat_KW_by_ecotype.tsv", sep="\t", index=False)
print("\nLM22 super-category KW vs user ecotype:")
print(sdf.to_string(index=False))

# Plot LM22 fractions x ecotype (8 super-cat × 3 boxplots)
fig, axes = plt.subplots(2, 4, figsize=(15, 8))
for ax, cat in zip(axes.flat, supercats):
    sns.boxplot(data=super_df, x="ecotype", y=cat, order=ec_order, palette=palette,
                ax=ax, showfliers=False)
    sns.stripplot(data=super_df, x="ecotype", y=cat, order=ec_order, color="black", size=2, alpha=.35, ax=ax)
    q = float(sdf.loc[sdf["category"]==cat,"q_BH"].values[0])
    star = "***" if q<1e-4 else ("**" if q<1e-2 else ("*" if q<0.05 else "ns"))
    ax.set_title(f"{cat}\nKW q={q:.1e} {star}", fontsize=10)
    ax.set_xlabel(""); ax.set_ylabel("LM22 fraction")
    ax.tick_params(axis="x", rotation=18, labelsize=9)
fig.suptitle("CIBERSORTx LM22 cell-type fractions across pediatric glioma immune ecotypes (Main n=349)", fontsize=12, y=1.02)
plt.tight_layout(); fig.savefig(FIG / "LM22_fractions_by_ecotype.png", dpi=DPI, bbox_inches="tight"); plt.close()

# ----------------------------------------------------------------------------
# 4) Cross-method correlation heatmap WITH LM22 added (5 methods)
# ----------------------------------------------------------------------------
def load_method(path, prefix):
    d = pd.read_csv(path, sep="\t", index_col=0)
    d.index.name = "cell_type"
    return d.T.add_prefix(prefix)

xcell = load_method(OUT / "immunedeconv_xcell.tsv", "xCell_")
qts   = load_method(OUT / "immunedeconv_quantiseq.tsv", "quanTIseq_")
mcp   = load_method(OUT / "immunedeconv_mcp_counter.tsv", "MCP_")
epic  = load_method(OUT / "immunedeconv_epic.tsv", "EPIC_")
lm22_t = df_rel.drop(columns=["ecotype"]).add_prefix("LM22_")

# Build composite cell-axis scores per method
def composite(d, cols):
    cols = [c for c in cols if c in d.columns]
    return d[cols].sum(axis=1) if cols else pd.Series(np.nan, index=d.index)

axes_def = {
    "CD8 T":       {"LM22": ["LM22_T cells CD8"],
                    "xCell": ["xCell_T cell CD8+"],
                    "quanTIseq": ["quanTIseq_T cell CD8+"],
                    "MCP": ["MCP_T cell CD8+"],
                    "EPIC": ["EPIC_T cell CD8+"]},
    "NK":          {"LM22": ["LM22_NK cells resting","LM22_NK cells activated"],
                    "xCell": ["xCell_NK cell"],
                    "quanTIseq": ["quanTIseq_NK cell"],
                    "MCP": ["MCP_NK cell"]},
    "Macrophage":  {"LM22": ["LM22_Macrophages M0","LM22_Macrophages M1","LM22_Macrophages M2"],
                    "xCell": ["xCell_Macrophage"],
                    "MCP": ["MCP_Macrophage/Monocyte"],
                    "EPIC": ["EPIC_Macrophage"]},
    "Monocyte":    {"LM22": ["LM22_Monocytes"],
                    "xCell": ["xCell_Monocyte"],
                    "MCP": ["MCP_Monocytic lineage"]},
}
combined = pd.concat([lm22_t, xcell, qts, mcp, epic], axis=1)
combined = combined.loc[combined.index.intersection(ec3.index)]

# Compute pairwise Spearman ρ between methods within each axis
all_rows = []
for axis, mp in axes_def.items():
    series = {m: composite(combined, cols) for m, cols in mp.items() if any(c in combined.columns for c in cols)}
    series = {m: s for m, s in series.items() if s.notna().any()}
    for m1 in series:
        for m2 in series:
            if m1 >= m2: continue
            x = series[m1]; y = series[m2]; ok = x.notna() & y.notna()
            if ok.sum() < 50: continue
            r, p = stats.spearmanr(x[ok], y[ok])
            all_rows.append({"axis": axis, "method_A": m1, "method_B": m2, "spearman_rho": r, "n": int(ok.sum())})
xmdf = pd.DataFrame(all_rows)
xmdf.to_csv(OUT / "lm22_added_cross_method_spearman.tsv", sep="\t", index=False)
print("\nLM22-included cross-method Spearman:")
print(xmdf.to_string(index=False))

# Pivot per axis and visualize
fig, axes = plt.subplots(1, 4, figsize=(20, 5))
for ax, axis in zip(axes, axes_def.keys()):
    sub = xmdf[xmdf["axis"]==axis]
    if sub.empty: ax.set_visible(False); continue
    methods_all = sorted(set(sub["method_A"]).union(sub["method_B"]))
    mat = pd.DataFrame(np.nan, index=methods_all, columns=methods_all)
    for _, r in sub.iterrows():
        mat.loc[r["method_A"], r["method_B"]] = r["spearman_rho"]
        mat.loc[r["method_B"], r["method_A"]] = r["spearman_rho"]
    np.fill_diagonal(mat.values, 1.0)
    sns.heatmap(mat, annot=True, fmt=".2f", cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                cbar_kws=dict(label="Spearman ρ"), ax=ax, square=True)
    ax.set_title(f"{axis} axis", fontsize=11)
plt.tight_layout(); fig.savefig(FIG / "cross_method_with_LM22_heatmap.png", dpi=DPI, bbox_inches="tight"); plt.close()

# ----------------------------------------------------------------------------
# 5) Summary numbers for manuscript insertion
# ----------------------------------------------------------------------------
summary = {
    "LM22_QC_Main_n349": {
        "n": int(len(qc_main)),
        "median_P": float(qc_main["P-value"].median()),
        "median_correlation": float(qc_main["Correlation"].median()),
        "median_RMSE": float(qc_main["RMSE"].median()),
        "n_P_lt_0.05": int((qc_main["P-value"]<0.05).sum()),
        "pct_P_lt_0.05": float((qc_main["P-value"]<0.05).mean()*100),
    },
    "LM22_supercat_KW_by_ecotype": sdf.to_dict(orient="records"),
    "LM22_supercat_significant_q05": int((sdf["q_BH"]<0.05).sum()),
}
(OUT / "lm22_integration_summary.json").write_text(json.dumps(summary, indent=2))
print("\n=== Summary written ===")
print(json.dumps(summary, indent=2))
