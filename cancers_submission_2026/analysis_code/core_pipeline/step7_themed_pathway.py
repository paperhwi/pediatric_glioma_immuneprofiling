"""
Step 7 — Themed pathway analysis on new LM22+ssGSEA ecotype assignment.

Groups 24 ssGSEA signatures into 5 biological themes and produces:
  - Theme composite score per sample
  - Theme × ecotype boxplots (single multi-panel figure)
  - Theme × ecotype × cohort_group cohort-stratified plots
  - Signature × ecotype heatmap (z-scored within signature)
  - Per-theme KW + Dunn BH summary table

Inputs
------
- output/step4_integrated_feature_matrix.tsv  (sample × feature merged)
- output/ecotype_LM22_main_k3_annotated.tsv   (k=3 ecotype, Main cohort)

Outputs (output/, figs_step7/)
"""
from __future__ import annotations
from pathlib import Path
import json, warnings
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_step7"; FIG.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
merged = pd.read_csv(OUT / "step4_integrated_feature_matrix.tsv", sep="\t", index_col=0)
ec3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
common = merged.index.intersection(ec3.index)
df = merged.loc[common].copy(); df["ecotype"] = ec3.loc[common, "ecotype"]
print(f"Step 7: {len(df)} samples in Main cohort with ecotype assignment")

# ---------------------------------------------------------------------------
themes = {
    "T-cell axis":            ["ssGSEA_T_Cell_Cytotoxicity", "ssGSEA_T_Cell_Exhaustion",
                                "ssGSEA_Tregs_Friebel2020", "ssGSEA_Chemokine_T_Cell_Recruitment",
                                "ssGSEA_NK_Cell_Activity", "ssGSEA_Dendritic_Cell_Activation"],
    "Antigen presentation":   ["ssGSEA_MHC_Class_I", "ssGSEA_MHC_Class_II"],
    "IFN / chemokine":        ["ssGSEA_IFN_Gamma_Response", "ssGSEA_IFN_Alpha_Response"],
    "Myeloid / microglia":    ["ssGSEA_Microglia_Core_Homeostatic", "ssGSEA_Microglia_Klemm2020",
                                "ssGSEA_MDM_Klemm2020", "ssGSEA_MgTAM_Antunes2021",
                                "ssGSEA_MoTAM_Antunes2021", "ssGSEA_DAM_KerenShaul2017",
                                "ssGSEA_M1_Macrophage", "ssGSEA_M2_Macrophage",
                                "ssGSEA_Neutrophil_Activation",
                                "ssGSEA_Glioma_Inflammatory_Wang2017"],
    "Signaling / suppression":["ssGSEA_MAPK_Activity", "ssGSEA_TGFb_Immunosuppression",
                                "ssGSEA_Cell_Cycle_Proliferation", "ssGSEA_Stemness_Brain_Tumor"],
}
themes = {k: [c for c in v if c in df.columns] for k, v in themes.items()}

# Theme composite = mean of z-scored signatures within theme (signatures already z-scored in step4 matrix)
theme_z = pd.DataFrame(index=df.index)
for tname, cols in themes.items():
    theme_z[tname] = df[cols].mean(axis=1)
theme_z["ecotype"] = df["ecotype"]
theme_z.to_csv(OUT / "step7_theme_composite_scores.tsv", sep="\t")

# ---------------------------------------------------------------------------
# KW + Dunn for each theme (composite) and each signature
# ---------------------------------------------------------------------------
def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); order = np.argsort(p)
    ranked = p[order]; q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1); return out

def dunn(groups, labels):
    all_v = np.concatenate(groups); ranks = stats.rankdata(all_v); N = len(all_v)
    _, c = np.unique(all_v, return_counts=True); T = (c**3 - c).sum(); tc = 1 - T/(N**3 - N)
    pos=0; mr=[]; ns=[]
    for g in groups:
        n=len(g); mr.append(ranks[pos:pos+n].mean()); ns.append(n); pos+=n
    res = {}
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            de = np.sqrt(((N*(N+1)/12)*(1/ns[i] + 1/ns[j]))*tc)
            z = (mr[i] - mr[j]) / de if de>0 else 0
            p = 2*(1-stats.norm.cdf(abs(z)))
            res[(labels[i], labels[j])] = (float(z), float(p))
    return res

ecs = sorted(df["ecotype"].dropna().unique())
theme_kw = []
for tname in themes:
    groups = [theme_z.loc[theme_z["ecotype"]==e, tname].dropna().values for e in ecs]
    H, p = stats.kruskal(*groups); N = sum(len(g) for g in groups); k = len(groups)
    eps = (H - k + 1) / (N - k)
    theme_kw.append({"theme": tname, "H": float(H), "p": float(p), "n": N, "epsilon_sq": float(eps)})
theme_kw_df = pd.DataFrame(theme_kw)
theme_kw_df["q_BH"] = bh(theme_kw_df["p"].values)
theme_kw_df.to_csv(OUT / "step7_theme_KW.tsv", sep="\t", index=False)
print(theme_kw_df.to_string(index=False))

# Dunn per theme
theme_dunn = []
for tname in themes:
    groups = [theme_z.loc[theme_z["ecotype"]==e, tname].dropna().values for e in ecs]
    d = dunn(groups, ecs); ps = [v[1] for v in d.values()]; qs = bh(ps)
    for (pair,(z,p)),q in zip(d.items(), qs):
        theme_dunn.append({"theme": tname, "group_A": pair[0], "group_B": pair[1], "z": z, "p": p, "q_BH": q})
pd.DataFrame(theme_dunn).to_csv(OUT / "step7_theme_Dunn.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# Plot 1: Theme composite boxplots (5 themes × 1 row)
# ---------------------------------------------------------------------------
ec_order = ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]
palette  = {"Lymphocyte-inflamed":"#3B82F6","Myeloid-dominant":"#EF4444","Immune-desert":"#9CA3AF"}
fig, axes = plt.subplots(1, len(themes), figsize=(4.2*len(themes), 4.6))
for ax, (tname, cols) in zip(axes, themes.items()):
    sns.boxplot(data=theme_z, x="ecotype", y=tname, order=ec_order, palette=palette, ax=ax, showfliers=False)
    sns.stripplot(data=theme_z, x="ecotype", y=tname, order=ec_order, color="black", size=2, alpha=.3, ax=ax)
    q = float(theme_kw_df.loc[theme_kw_df["theme"]==tname,"q_BH"].values[0])
    eps = float(theme_kw_df.loc[theme_kw_df["theme"]==tname,"epsilon_sq"].values[0])
    star = "***" if q<1e-4 else ("**" if q<1e-2 else ("*" if q<5e-2 else "ns"))
    ax.set_title(f"{tname}\nKW q={q:.1e} {star} · ε²={eps:.2f}", fontsize=10)
    ax.set_xlabel(""); ax.set_ylabel("theme composite z")
    ax.tick_params(axis="x", rotation=18, labelsize=9)
plt.tight_layout(); fig.savefig(FIG / "theme_composite_boxplot.png", dpi=200); plt.close()

# ---------------------------------------------------------------------------
# Plot 2: Signature × ecotype heatmap (mean z per ecotype, signatures grouped by theme)
# ---------------------------------------------------------------------------
sigs_ordered = []
for cols in themes.values(): sigs_ordered += cols
sig_means = df.groupby("ecotype")[sigs_ordered].mean().T
sig_means.index = [c.replace("ssGSEA_","") for c in sig_means.index]
sig_means = sig_means.loc[:, ec_order]
fig, ax = plt.subplots(figsize=(6, 0.45*len(sig_means)+1.0))
sns.heatmap(sig_means, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
            cbar_kws=dict(label="mean z within ecotype"), ax=ax)
# Draw theme separators
yticklabels = list(sig_means.index); pos = 0
for tname, cols in themes.items():
    pos += len(cols)
    if pos < len(yticklabels):
        ax.axhline(pos, color="black", lw=1.0)
ax.set_title("Signature × ecotype  (themes separated by black lines)")
ax.set_ylabel(""); ax.set_xlabel("")
plt.tight_layout(); fig.savefig(FIG / "signature_ecotype_heatmap.png", dpi=200); plt.close()

# ---------------------------------------------------------------------------
# Plot 3: Cohort-stratified per-theme means
# ---------------------------------------------------------------------------
theme_z["cohort_group"] = df["cohort_group"]
fig, axes = plt.subplots(1, len(themes), figsize=(4.2*len(themes), 4.6))
cohort_order = ["DMG_K27", "DHG_G34", "pHGG_WT", "IHG"]
cohort_order = [c for c in cohort_order if c in theme_z["cohort_group"].unique()]
for ax, (tname, _) in zip(axes, themes.items()):
    sns.boxplot(data=theme_z, x="cohort_group", y=tname, hue="ecotype",
                order=cohort_order, hue_order=ec_order, palette=palette, ax=ax, showfliers=False)
    ax.set_title(tname, fontsize=10); ax.set_xlabel(""); ax.set_ylabel("theme z")
    ax.tick_params(axis="x", rotation=15, labelsize=8)
    if ax is not axes[-1]:
        ax.get_legend().remove() if ax.get_legend() else None
    else:
        ax.legend(loc="upper right", fontsize=7, frameon=False)
plt.tight_layout(); fig.savefig(FIG / "theme_by_cohort_stratified.png", dpi=200); plt.close()

# ---------------------------------------------------------------------------
# Summary JSON
# ---------------------------------------------------------------------------
summary = {
    "n_samples": int(len(df)),
    "theme_KW": theme_kw_df.assign(q_BH=lambda d: d["q_BH"].astype(float)).to_dict(orient="records"),
}
(OUT / "step7_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
