"""
Step 6 — Cell-type fraction × ecotype statistics
================================================

For every feature in the merged Step 4 integrated feature matrix (LM22, ssGSEA,
xCell, quanTIseq, MCP-counter, EPIC), compare the three k=3 ecotypes via:
  1. Kruskal–Wallis omnibus (BH FDR across all features, separately per method family)
  2. Dunn post-hoc (BH FDR within feature) for features with KW q < 0.05
  3. Per-feature effect size (epsilon^2)

Generates:
  - step6_KW_results.tsv (one row per feature: KW H, df, p, q, n)
  - step6_Dunn_posthoc.tsv (pairwise Dunn z, p, q for KW-significant features)
  - step6_effect_sizes.tsv (eta² + epsilon² per feature)
  - figs_step6/* : KW -log10(q) bar; Microglia/MDM, CD8 T, M0/M1/M2 boxplots
  - Step6_KW_Dunn_ecotype.ipynb
"""
from __future__ import annotations
from pathlib import Path
import json, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_step6"; FIG.mkdir(exist_ok=True)

# ----------------------------------------------------------------------------
# 1) Load merged feature matrix + ecotype assignment (k=3, Main cohort)
# ----------------------------------------------------------------------------
merged = pd.read_csv(OUT / "step4_integrated_feature_matrix.tsv", sep="\t", index_col=0)
ec3    = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
common = merged.index.intersection(ec3.index)
df = merged.loc[common].copy()
df["ecotype"] = ec3.loc[common, "ecotype"]
print(f"[Step 6] {len(df)} samples × {df.shape[1]-1} features (incl. annotation)")
print(df["ecotype"].value_counts())

# Drop annotation/QC columns; keep feature columns
non_feat = ["Kids_First_Participant_ID","cohort_group","harmonized_diagnosis","molecular_subtype",
            "location_class","age_years","age_dev_group","reported_gender","has_survival","ecotype"]
qc_prefix = ("LM22_QC_",)
feature_cols = [c for c in df.columns if c not in non_feat and not c.startswith(qc_prefix)]
def family_of(c):
    if c.startswith("LM22_"): return "LM22"
    if c.startswith("ssGSEA_"): return "ssGSEA"
    if c.startswith("xCell_"): return "xCell"
    if c.startswith("quanTIseq_"): return "quanTIseq"
    if c.startswith("MCP_"): return "MCP"
    if c.startswith("EPIC_"): return "EPIC"
    return "Other"

print(f"  features by family: {pd.Series([family_of(c) for c in feature_cols]).value_counts().to_dict()}")

# ----------------------------------------------------------------------------
# 2) Kruskal–Wallis + effect size per feature
# ----------------------------------------------------------------------------
def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p); order = np.argsort(p); ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1); return out

ecotypes = sorted(df["ecotype"].dropna().unique().tolist())
print(f"  ecotypes: {ecotypes}")

kw_rows = []
for c in feature_cols:
    groups = [df.loc[df["ecotype"] == e, c].dropna().values for e in ecotypes]
    sizes  = [len(g) for g in groups]
    if min(sizes) < 5 or sum(sizes) < 50:
        kw_rows.append({"feature": c, "family": family_of(c), "H": np.nan, "df": np.nan,
                        "p": np.nan, "n": sum(sizes), "epsilon_sq": np.nan})
        continue
    H, p = stats.kruskal(*groups)
    N = sum(sizes); k = len(groups)
    # epsilon² effect size (Tomczak/Tomczak 2014)
    eps_sq = (H - k + 1) / (N - k)
    kw_rows.append({"feature": c, "family": family_of(c), "H": float(H), "df": k - 1,
                    "p": float(p), "n": N, "epsilon_sq": float(eps_sq)})
kw = pd.DataFrame(kw_rows)
# BH within each family
kw["q_BH"] = np.nan
for fam, sub in kw.groupby("family"):
    valid = sub["p"].notna()
    if valid.any():
        kw.loc[sub.index[valid], "q_BH"] = bh(sub.loc[valid, "p"].values)
kw = kw.sort_values("q_BH")
kw.to_csv(OUT / "step6_KW_results.tsv", sep="\t", index=False)
print(f"  KW significant (q<0.05): {(kw['q_BH'] < 0.05).sum()} / {len(kw)}")

# ----------------------------------------------------------------------------
# 3) Dunn post-hoc on KW-significant features (BH within feature)
# ----------------------------------------------------------------------------
def dunn(groups, labels):
    """Dunn (1964) post-hoc; returns dict of (i,j) -> (z, p two-sided)."""
    all_vals = np.concatenate(groups)
    ranks    = stats.rankdata(all_vals)
    n_total  = len(all_vals)
    # Tie correction
    _, counts = np.unique(all_vals, return_counts=True)
    T = (counts ** 3 - counts).sum()
    tie_corr = 1 - T / (n_total ** 3 - n_total)
    pos = 0; meanranks = []; ns = []
    for g in groups:
        n = len(g)
        meanranks.append(ranks[pos:pos+n].mean()); ns.append(n); pos += n
    res = {}
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            denom = np.sqrt(((n_total*(n_total+1)/12) * (1/ns[i] + 1/ns[j])) * tie_corr)
            z = (meanranks[i] - meanranks[j]) / denom if denom > 0 else 0
            p = 2 * (1 - stats.norm.cdf(abs(z)))
            res[(labels[i], labels[j])] = (float(z), float(p))
    return res

sig_features = kw.loc[kw["q_BH"] < 0.05, "feature"].tolist()
print(f"  Running Dunn on {len(sig_features)} KW-significant features")
dunn_rows = []
for c in sig_features:
    groups = [df.loc[df["ecotype"] == e, c].dropna().values for e in ecotypes]
    res = dunn(groups, ecotypes)
    ps  = [v[1] for v in res.values()]
    qs  = bh(ps)
    for (pair, (z, p)), q in zip(res.items(), qs):
        dunn_rows.append({
            "feature": c, "family": family_of(c),
            "group_A": pair[0], "group_B": pair[1],
            "z": z, "p": p, "q_BH_per_feature": q,
            "mean_A_z_or_frac": float(df.loc[df["ecotype"] == pair[0], c].dropna().mean()),
            "mean_B_z_or_frac": float(df.loc[df["ecotype"] == pair[1], c].dropna().mean()),
        })
dunn_df = pd.DataFrame(dunn_rows).sort_values(["family", "feature"])
dunn_df.to_csv(OUT / "step6_Dunn_posthoc.tsv", sep="\t", index=False)
print(f"  Dunn pairwise rows: {len(dunn_df)}")

# Wide-format: per feature, columns = pair, value = signed -log10(q) (sign = mean diff)
def signed_log10(z, q):
    sign = np.sign(z) if z != 0 else 0
    return float(sign * -np.log10(max(q, 1e-300)))
dunn_wide = (dunn_df
             .assign(slq=lambda d: d.apply(lambda r: signed_log10(r["z"], r["q_BH_per_feature"]), axis=1),
                     pair=lambda d: d["group_A"] + " vs " + d["group_B"])
             .pivot(index="feature", columns="pair", values="slq"))
dunn_wide.to_csv(OUT / "step6_Dunn_signed_log10q_wide.tsv", sep="\t")

# ----------------------------------------------------------------------------
# 4) Plot: KW -log10(q) bar (top 30 by family)
# ----------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
families_plot = ["LM22", "ssGSEA", "xCell", "quanTIseq", "MCP", "EPIC"]
for ax, fam in zip(axes.flat, families_plot):
    sub = kw[kw["family"] == fam].sort_values("q_BH").head(20)
    sub = sub.dropna(subset=["q_BH"])
    if sub.empty:
        ax.set_visible(False); continue
    nl = -np.log10(np.clip(sub["q_BH"], 1e-300, 1))
    short = sub["feature"].str.replace(f"^{fam}_", "", regex=True)
    ax.barh(short[::-1], nl[::-1], color="#3B82F6")
    ax.axvline(-np.log10(0.05), color="red", ls="--", alpha=.5)
    ax.set_title(f"{fam} — KW q (top 20 by significance)")
    ax.set_xlabel("-log10(q)")
plt.tight_layout()
fig.savefig(FIG / "KW_significance_by_family.png", dpi=200)
plt.close()

# ----------------------------------------------------------------------------
# 5) Focused boxplots: microglia/MDM, CD8 T, Macrophage M0/M1/M2, Neutrophil, MAPK
# ----------------------------------------------------------------------------
focus = [
    ("LM22_T cells CD8",                  "LM22 CD8 T cells (relative)"),
    ("LM22_Macrophages M0",               "LM22 Macrophages M0"),
    ("LM22_Macrophages M1",               "LM22 Macrophages M1"),
    ("LM22_Macrophages M2",               "LM22 Macrophages M2"),
    ("LM22_Monocytes",                    "LM22 Monocytes"),
    ("LM22_Neutrophils",                  "LM22 Neutrophils"),
    ("ssGSEA_Microglia_Klemm2020",        "ssGSEA Microglia (Klemm 2020)"),
    ("ssGSEA_MDM_Klemm2020",              "ssGSEA MDM (Klemm 2020)"),
    ("ssGSEA_MgTAM_Antunes2021",          "ssGSEA Mg-TAM (Antunes 2021)"),
    ("ssGSEA_MoTAM_Antunes2021",          "ssGSEA Mo-TAM (Antunes 2021)"),
    ("ssGSEA_T_Cell_Cytotoxicity",        "ssGSEA T-cell cytotoxicity"),
    ("ssGSEA_IFN_Gamma_Response",         "ssGSEA IFN-γ response"),
]
focus = [(c, label) for c, label in focus if c in df.columns]

ec_order = ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]
palette  = {"Lymphocyte-inflamed": "#3B82F6", "Myeloid-dominant": "#EF4444", "Immune-desert": "#9CA3AF"}

fig, axes = plt.subplots(3, 4, figsize=(15, 11))
for ax, (col, lab) in zip(axes.flat, focus):
    sns.boxplot(data=df, x="ecotype", y=col, order=ec_order, palette=palette, ax=ax, showfliers=False)
    sns.stripplot(data=df, x="ecotype", y=col, order=ec_order, color="black", size=2, alpha=.35, ax=ax)
    q = kw.loc[kw["feature"] == col, "q_BH"]
    qv = float(q.values[0]) if len(q) else np.nan
    star = "***" if qv < 1e-4 else ("**" if qv < 1e-2 else ("*" if qv < 5e-2 else "ns"))
    ax.set_title(f"{lab}\nKW q={qv:.1e} {star}", fontsize=9)
    ax.set_xlabel(""); ax.set_ylabel("z-score / fraction")
    ax.tick_params(axis="x", rotation=15, labelsize=8)
for ax in axes.flat[len(focus):]:
    ax.set_visible(False)
plt.tight_layout()
fig.savefig(FIG / "ecotype_boxplots_focus.png", dpi=200)
plt.close()

# ----------------------------------------------------------------------------
# 6) Dunn pairwise signed -log10(q) heatmap for focus features
# ----------------------------------------------------------------------------
heat_features = [c for c, _ in focus if c in dunn_wide.index]
fig, ax = plt.subplots(figsize=(6.5, max(4, 0.45*len(heat_features))))
sns.heatmap(dunn_wide.loc[heat_features], annot=True, fmt=".1f", cmap="RdBu_r", center=0,
            cbar_kws=dict(label="signed -log10(q) Dunn BH"), ax=ax)
ax.set_title("Dunn post-hoc (BH per feature) — focus signatures")
plt.tight_layout()
fig.savefig(FIG / "Dunn_signed_log10q_focus.png", dpi=200)
plt.close()

# ----------------------------------------------------------------------------
# 7) Summary JSON
# ----------------------------------------------------------------------------
summary = {
    "n_samples": int(len(df)),
    "ecotype_sizes": df["ecotype"].value_counts().to_dict(),
    "features_tested": int(len(feature_cols)),
    "KW_significant_BH05_total":  int((kw["q_BH"] < 0.05).sum()),
    "KW_significant_BH05_by_family": kw[kw["q_BH"] < 0.05]["family"].value_counts().to_dict(),
    "top10_features_by_q": kw.head(10)[["feature","family","H","q_BH","epsilon_sq"]].to_dict(orient="records"),
}
(OUT / "step6_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2, default=str))
