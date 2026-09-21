"""
Step 4 integration: CIBERSORTx LM22 + ssGSEA brain-immune + 4-method immunedeconv
========================================================================
Output goal: tidy sample-level merged feature matrix + cross-method correlation
            + microglia/MDM (ssGSEA) vs. macrophage (LM22) diagnostics.

This module is converted to an ipynb at the end of the run.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

BASE  = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT   = BASE / "output"
FIGS  = OUT / "figs_step4"
FIGS.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load CIBERSORTx LM22 result (sample × feature)
# ---------------------------------------------------------------------------
lm22 = pd.read_csv(OUT / "CIBERSORTx_Job15_Results.txt", sep="\t")
lm22 = lm22.rename(columns={"Mixture": "sample"}).set_index("sample")
qc_cols = ["P-value", "Correlation", "RMSE", "Absolute score (sig.score)"]
cell_cols = [c for c in lm22.columns if c not in qc_cols]
lm22_cells = lm22[cell_cols].copy()    # absolute scores per cell type (sum = absolute score)
lm22_qc    = lm22[qc_cols].copy()

# Renormalize to relative (sum-to-1) per sample for cross-method comparison
lm22_rel = lm22_cells.div(lm22_cells.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
lm22_rel = lm22_rel.add_prefix("LM22_")

print(f"[LM22] n={len(lm22)} samples, k={len(cell_cols)} cell types")
print(f"[LM22] P<0.05 high-confidence: n={(lm22_qc['P-value'] < 0.05).sum()}/{len(lm22_qc)}")

# ---------------------------------------------------------------------------
# 2. Load other deconvolution outputs (feature × sample → transpose)
# ---------------------------------------------------------------------------
def load_long(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.T                                   # rows=sample, cols=feature
    df.index.name = "sample"
    df = df.add_prefix(f"{label}_")
    return df

ssgsea = load_long(OUT / "ssGSEA_brain_immune_scores.tsv", "ssGSEA")
xcell  = load_long(OUT / "immunedeconv_xcell.tsv",       "xCell")
qts    = load_long(OUT / "immunedeconv_quantiseq.tsv",   "quanTIseq")
mcp    = load_long(OUT / "immunedeconv_mcp_counter.tsv", "MCP")
epic   = load_long(OUT / "immunedeconv_epic.tsv",        "EPIC")

print(f"[ssGSEA]    {ssgsea.shape}")
print(f"[xCell]     {xcell.shape}")
print(f"[quanTIseq] {qts.shape}")
print(f"[MCP]       {mcp.shape}")
print(f"[EPIC]      {epic.shape}")

# ---------------------------------------------------------------------------
# 3. Load sample annotation
# ---------------------------------------------------------------------------
manifest = pd.read_csv(OUT / "tpm_manifest.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
manifest.index.name = "sample"
print(f"[manifest]  {manifest.shape}")

# ---------------------------------------------------------------------------
# 4. Merge everything on sample
# ---------------------------------------------------------------------------
merged = (
    manifest
    .join(lm22_rel, how="inner")
    .join(lm22_qc.rename(columns=lambda c: f"LM22_QC_{c.replace(' ','_')}"), how="left")
    .join(ssgsea, how="left")
    .join(xcell,  how="left")
    .join(qts,    how="left")
    .join(mcp,    how="left")
    .join(epic,   how="left")
)
print(f"[merged]    {merged.shape}")
merged.to_csv(OUT / "step4_integrated_feature_matrix.tsv", sep="\t")

# ---------------------------------------------------------------------------
# 5. Cross-method correlation for key cell-type axes
# ---------------------------------------------------------------------------
axes = {
    "CD8 T": {
        "LM22":      "LM22_T cells CD8",
        "xCell":     "xCell_T cell CD8+",
        "quanTIseq": "quanTIseq_T cell CD8+",
        "MCP":       "MCP_T cell CD8+",
        "EPIC":      "EPIC_T cell CD8+",
    },
    "NK": {
        "LM22":      None,                       # combine resting + activated
        "xCell":     "xCell_NK cell",
        "quanTIseq": "quanTIseq_NK cell",
        "MCP":       "MCP_NK cell",
    },
    "Macrophage/Monocyte (LM22) vs Macrophage": {
        "LM22":      None,                       # combine M0+M1+M2 + Monocytes
        "xCell":     "xCell_Macrophage",
        "quanTIseq": "quanTIseq_Macrophage M1",   # M1+M2 sum below
        "MCP":       "MCP_Macrophage/Monocyte",
        "EPIC":      "EPIC_Macrophage",
    },
    "Myeloid composite": {
        "LM22":      None,                       # M0+M1+M2+Monocyte+Neutrophil
        "xCell":     "xCell_Monocyte",
        "MCP":       "MCP_Monocytic lineage",
    },
}

# Build aggregate columns where None
agg = pd.DataFrame(index=merged.index)
agg["CD8 T_LM22"]                     = merged["LM22_T cells CD8"]
agg["NK_LM22"]                        = merged[["LM22_NK cells resting", "LM22_NK cells activated"]].sum(axis=1)
agg["Macrophage_LM22"]                = merged[["LM22_Macrophages M0","LM22_Macrophages M1","LM22_Macrophages M2"]].sum(axis=1)
agg["Macrophage+Monocyte_LM22"]       = agg["Macrophage_LM22"] + merged["LM22_Monocytes"]
agg["Myeloid_LM22"]                   = agg["Macrophage+Monocyte_LM22"] + merged["LM22_Neutrophils"]

# Pair with companions for cross-method correlation
pair_specs = [
    ("CD8 T",         "CD8 T_LM22",                     ["xCell_T cell CD8+", "quanTIseq_T cell CD8+", "MCP_T cell CD8+", "EPIC_T cell CD8+"]),
    ("NK",            "NK_LM22",                        ["xCell_NK cell", "quanTIseq_NK cell", "MCP_NK cell"]),
    ("Macrophage",    "Macrophage_LM22",                ["xCell_Macrophage", "MCP_Macrophage/Monocyte", "EPIC_Macrophage"]),
    ("Myeloid",       "Myeloid_LM22",                   ["xCell_Monocyte", "MCP_Monocytic lineage"]),
]

rows = []
for axis, lm22_col, peers in pair_specs:
    for peer in peers:
        if peer not in merged.columns:
            continue
        x = agg[lm22_col]
        y = merged[peer]
        ok = x.notna() & y.notna()
        if ok.sum() < 50:
            continue
        rho, p = spearmanr(x[ok], y[ok])
        rows.append({"axis": axis, "LM22 column": lm22_col, "peer method": peer, "spearman_rho": rho, "p": p, "n": int(ok.sum())})
cross_corr = pd.DataFrame(rows)
cross_corr.to_csv(OUT / "step4_LM22_vs_methods_spearman.tsv", sep="\t", index=False)
print(cross_corr.to_string(index=False))

# ---------------------------------------------------------------------------
# 6. Plot: heatmap of LM22 vs peers spearman
# ---------------------------------------------------------------------------
pivot = cross_corr.pivot(index="peer method", columns="axis", values="spearman_rho")
fig, ax = plt.subplots(figsize=(6, max(3, 0.45*len(pivot))))
sns.heatmap(pivot, annot=True, fmt=".2f", cmap="RdBu_r", center=0, vmin=-1, vmax=1, ax=ax,
            cbar_kws=dict(label="Spearman ρ (LM22 vs peer)"))
ax.set_title("Cross-method consensus on key immune axes (n=702 samples)")
plt.tight_layout()
fig.savefig(FIGS / "LM22_vs_methods_heatmap.png", dpi=200)
plt.close()

# ---------------------------------------------------------------------------
# 7. Microglia (ssGSEA) vs LM22 macrophage axes — Path B headline diagnostic
# ---------------------------------------------------------------------------
microglia_cols = [c for c in merged.columns if c.startswith("ssGSEA_") and ("Microglia" in c or "TAM" in c or "DAM" in c or "Mg" in c or "MDM" in c or "Mo_TAM" in c)]
print(f"[ssGSEA microglia-related signatures] {len(microglia_cols)} sets")
for c in microglia_cols: print("  -", c)

mg_vs_lm22 = []
for mg in microglia_cols:
    for lm22_axis, col in [("LM22 Macrophages (M0+M1+M2)", "Macrophage_LM22"),
                            ("LM22 Macrophage+Monocyte",    "Macrophage+Monocyte_LM22"),
                            ("LM22 M2",                      "LM22_Macrophages M2"),
                            ("LM22 M1",                      "LM22_Macrophages M1"),
                            ("LM22 Monocytes",               "LM22_Monocytes")]:
        x = merged[mg]
        y = agg[col] if col in agg.columns else merged[col]
        ok = x.notna() & y.notna()
        rho, p = spearmanr(x[ok], y[ok])
        mg_vs_lm22.append({"ssGSEA signature": mg.replace("ssGSEA_",""), "LM22 axis": lm22_axis, "spearman_rho": rho, "p": p, "n": int(ok.sum())})
mg_df = pd.DataFrame(mg_vs_lm22)
mg_df.to_csv(OUT / "step4_microglia_ssGSEA_vs_LM22_macrophage.tsv", sep="\t", index=False)
mg_pivot = mg_df.pivot(index="ssGSEA signature", columns="LM22 axis", values="spearman_rho")

fig, ax = plt.subplots(figsize=(7, max(4, 0.4*len(mg_pivot))))
sns.heatmap(mg_pivot, annot=True, fmt=".2f", cmap="RdBu_r", center=0, vmin=-1, vmax=1, ax=ax,
            cbar_kws=dict(label="Spearman ρ"))
ax.set_title("Brain-tuned ssGSEA vs. LM22 macrophage / monocyte axes (Path B)")
plt.tight_layout()
fig.savefig(FIGS / "ssGSEA_microglia_vs_LM22.png", dpi=200)
plt.close()

# ---------------------------------------------------------------------------
# 8. CIBERSORTx QC distribution
# ---------------------------------------------------------------------------
fig, axes2 = plt.subplots(1, 3, figsize=(12, 3.5))
for ax, col, name in zip(axes2, ["LM22_QC_P-value", "LM22_QC_Correlation", "LM22_QC_RMSE"],
                                ["P-value", "Correlation", "RMSE"]):
    sns.histplot(merged[col].dropna(), bins=40, ax=ax)
    ax.set_title(name)
    if name == "P-value":
        ax.axvline(0.05, color="red", linestyle="--", label="0.05")
        ax.legend()
plt.tight_layout()
fig.savefig(FIGS / "LM22_QC_distribution.png", dpi=200)
plt.close()

# ---------------------------------------------------------------------------
# 9. Ecotype-friendly feature matrix for Step 5 consensus clustering
#    (LM22 relative fractions + selected ssGSEA brain-immune signatures)
# ---------------------------------------------------------------------------
ssGSEA_keep = [c for c in merged.columns if c.startswith("ssGSEA_")]
feature_step5 = pd.concat([lm22_rel, merged[ssGSEA_keep]], axis=1)
# z-score within feature
feat_z = (feature_step5 - feature_step5.mean()) / feature_step5.std().replace(0, np.nan)
feat_z.to_csv(OUT / "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t")
print(f"[Step5 feature matrix] {feat_z.shape}")

# ---------------------------------------------------------------------------
# 10. Save summary JSON
# ---------------------------------------------------------------------------
summary = {
    "n_samples": int(len(merged)),
    "n_LM22_high_confidence_p05": int((lm22_qc['P-value'] < 0.05).sum()),
    "n_LM22_high_confidence_p10": int((lm22_qc['P-value'] < 0.10).sum()),
    "LM22_pvalue_median": float(lm22_qc['P-value'].median()),
    "LM22_correlation_median": float(lm22_qc['Correlation'].median()),
    "LM22_RMSE_median": float(lm22_qc['RMSE'].median()),
    "step4_outputs": [
        "step4_integrated_feature_matrix.tsv",
        "step4_LM22_vs_methods_spearman.tsv",
        "step4_microglia_ssGSEA_vs_LM22_macrophage.tsv",
        "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv",
        "figs_step4/LM22_vs_methods_heatmap.png",
        "figs_step4/ssGSEA_microglia_vs_LM22.png",
        "figs_step4/LM22_QC_distribution.png",
    ],
}
(OUT / "step4_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
