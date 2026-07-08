"""
Re-render every manuscript / supplementary figure at dpi=300 from saved data.
Outputs into output/figs_300dpi/ for collection into the 'for submission' folder.

This script reads the SAVED TSV / npz results produced by Step 4-10 / Suppl S3-S4
so we don't have to re-run the heavy computation; we only redraw the panels.
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import warnings
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"
DST  = OUT / "figs_300dpi"; DST.mkdir(exist_ok=True)
DPI  = 300

PAL_EC = {"Lymphocyte-inflamed": "#3B82F6",
          "Myeloid-dominant":    "#EF4444",
          "Immune-desert":       "#9CA3AF"}
EC_ORDER = ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]

# ---------------------------------------------------------------------------
# Re-render Fig 1 via the original script (already uses matplotlib)
# ---------------------------------------------------------------------------
print("[Fig 1] regenerating workflow ...")
import subprocess, sys
fig1_src = BASE / "scripts" / "build_figure1_workflow.py"
fig1_text = fig1_src.read_text()
fig1_text_300 = fig1_text.replace("dpi=240", f"dpi={DPI}")
tmp = Path("/tmp/_fig1_300.py")
tmp.write_text(fig1_text_300.replace(
    'FIG  = OUT / "figs_fig1"',
    f'FIG  = OUT / "figs_300dpi" / "Fig1_workflow_300dpi"'
))
subprocess.run([sys.executable, str(tmp)], check=True, capture_output=True)
# Move into single folder
for f in (DST / "Fig1_workflow_300dpi").glob("*"):
    f.rename(DST / f.name)
try: (DST / "Fig1_workflow_300dpi").rmdir()
except Exception: pass

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
ann_full = pd.read_csv(OUT / "cohort_for_deconvolution.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
ec3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
merged = pd.read_csv(OUT / "step4_integrated_feature_matrix.tsv", sep="\t", index_col=0)
df_main = merged.loc[merged.index.intersection(ec3.index)].copy()
df_main["ecotype"] = ec3.loc[df_main.index, "ecotype"]

# ---------------------------------------------------------------------------
# Fig 2 panels — LM22 QC, cross-method heatmap, ssGSEA microglia vs LM22 macrophage
# ---------------------------------------------------------------------------
print("[Fig 2] LM22 QC / cross-method / ssGSEA × LM22 ...")
lm22 = pd.read_csv(OUT / "CIBERSORTx_Job15_Results.txt", sep="\t")
qc = lm22[["Mixture", "P-value", "Correlation", "RMSE"]].set_index("Mixture")
fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
for ax, col, lab in zip(axes, ["P-value", "Correlation", "RMSE"], ["P-value", "Correlation", "RMSE"]):
    sns.histplot(qc[col], bins=40, ax=ax)
    ax.set_title(lab); ax.set_xlabel(lab)
    if col == "P-value":
        ax.axvline(0.05, color="red", ls="--", alpha=.7, label="0.05"); ax.legend()
plt.tight_layout(); fig.savefig(DST / "Fig2A_LM22_QC.png", dpi=DPI); plt.close()

corr = pd.read_csv(OUT / "step4_LM22_vs_methods_spearman.tsv", sep="\t")
pivot = corr.pivot(index="peer method", columns="axis", values="spearman_rho")
fig, ax = plt.subplots(figsize=(6, max(3, 0.5*len(pivot))))
sns.heatmap(pivot, annot=True, fmt=".2f", cmap="RdBu_r", center=0, vmin=-1, vmax=1,
            cbar_kws=dict(label="Spearman ρ (LM22 vs peer)"), ax=ax)
ax.set_title("Cross-method consensus on key immune axes (n=702)")
plt.tight_layout(); fig.savefig(DST / "Fig2B_LM22_vs_methods_heatmap.png", dpi=DPI); plt.close()

mg = pd.read_csv(OUT / "step4_microglia_ssGSEA_vs_LM22_macrophage.tsv", sep="\t")
mp = mg.pivot(index="ssGSEA signature", columns="LM22 axis", values="spearman_rho")
fig, ax = plt.subplots(figsize=(7, max(4, 0.45*len(mp))))
sns.heatmap(mp, annot=True, fmt=".2f", cmap="RdBu_r", center=0, vmin=-1, vmax=1,
            cbar_kws=dict(label="Spearman ρ"), ax=ax)
ax.set_title("Brain-tuned ssGSEA vs. LM22 macrophage/monocyte axes")
plt.tight_layout(); fig.savefig(DST / "Fig2C_ssGSEA_vs_LM22.png", dpi=DPI); plt.close()

# ---------------------------------------------------------------------------
# Fig 3 panels — PAC/silhouette, PCA/UMAP, feature mean heatmap, cohort stacked
# ---------------------------------------------------------------------------
print("[Fig 3] PAC/silhouette, PCA/UMAP, feature heatmap, cohort stack ...")
pac = pd.read_csv(OUT / "consensus_LM22_main_PAC.tsv", sep="\t")
sil = pd.read_csv(OUT / "consensus_LM22_main_silhouette.tsv", sep="\t")
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
ax[0].plot(pac["k"], pac["PAC"], "o-", color="#1f77b4"); ax[0].axvline(3, color="red", ls="--", alpha=.5, label="k=3")
ax[0].set_xlabel("k"); ax[0].set_ylabel("PAC"); ax[0].set_title("Proportion of Ambiguous Clustering"); ax[0].legend()
ax[1].plot(sil["k"], sil["silhouette"], "o-", color="#ff7f0e"); ax[1].axvline(3, color="red", ls="--", alpha=.5)
ax[1].set_xlabel("k"); ax[1].set_ylabel("Silhouette"); ax[1].set_title("Silhouette on 1−consensus")
plt.tight_layout(); fig.savefig(DST / "Fig3A_PAC_silhouette.png", dpi=DPI); plt.close()

pca_df  = pd.read_csv(OUT / "ecotype_LM22_main_pca.tsv", sep="\t")
umap_df = pd.read_csv(OUT / "ecotype_LM22_main_umap.tsv", sep="\t") if (OUT / "ecotype_LM22_main_umap.tsv").exists() else None
ec_col = "ecotype_k3" if "ecotype_k3" in pca_df.columns else "ecotype"
fig, axes = plt.subplots(1, 2 if umap_df is not None else 1, figsize=(11, 4.2))
ax0 = axes[0] if umap_df is not None else axes
for ec, sub in pca_df.groupby(ec_col):
    ax0.scatter(sub["PC1"], sub["PC2"], c=PAL_EC.get(ec, "#888"), label=ec, s=18, alpha=.75, edgecolor="none")
ax0.set_xlabel("PC1"); ax0.set_ylabel("PC2"); ax0.set_title("PCA — k=3"); ax0.legend(loc="best", frameon=False, fontsize=8)
if umap_df is not None:
    ec_col_u = "ecotype_k3" if "ecotype_k3" in umap_df.columns else "ecotype"
    for ec, sub in umap_df.groupby(ec_col_u):
        axes[1].scatter(sub["UMAP1"], sub["UMAP2"], c=PAL_EC.get(ec, "#888"), label=ec, s=18, alpha=.75, edgecolor="none")
    axes[1].set_xlabel("UMAP1"); axes[1].set_ylabel("UMAP2"); axes[1].set_title("UMAP")
plt.tight_layout(); fig.savefig(DST / "Fig3B_PCA_UMAP.png", dpi=DPI); plt.close()

means = pd.read_csv(OUT / "ecotype_LM22_main_k3_feature_means.tsv", sep="\t", index_col=0)
map_idx = {2: "Lymphocyte-inflamed", 3: "Myeloid-dominant", 1: "Immune-desert"}
means.index = [f"{map_idx.get(int(i), i)} (c{int(i)})" for i in means.index]
spread = means.max() - means.min()
keep = spread.sort_values(ascending=False).head(30).index
fig, ax = plt.subplots(figsize=(8, 9))
sns.heatmap(means[keep].T, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
            cbar_kws=dict(label="z-score (mean per cluster)"), ax=ax)
ax.set_title("k=3 ecotype × top-30 discriminating features")
plt.tight_layout(); fig.savefig(DST / "Fig3C_feature_means_heatmap.png", dpi=DPI); plt.close()

ct = pd.read_csv(OUT / "ecotype_LM22_main_k3_x_cohort_group.tsv", sep="\t").set_index("ecotype")
ct_pct = ct.div(ct.sum(axis=0), axis=1) * 100
fig, ax = plt.subplots(figsize=(6.5, 4.2))
bottom = np.zeros(len(ct_pct.columns))
for ec in EC_ORDER:
    if ec in ct_pct.index:
        ax.bar(ct_pct.columns, ct_pct.loc[ec], bottom=bottom, color=PAL_EC[ec], label=ec); bottom += ct_pct.loc[ec].values
ax.set_ylabel("% within cohort"); ax.set_ylim(0, 100)
ax.set_title("Ecotype composition by cohort (k=3)"); ax.legend(loc="upper right", fontsize=8, frameon=False)
plt.tight_layout(); fig.savefig(DST / "Fig3D_cohort_stacked.png", dpi=DPI); plt.close()

# ---------------------------------------------------------------------------
# Fig 4 panels — theme composites, signature heatmap, cohort-stratified, focus boxplots, Dunn heatmap
# ---------------------------------------------------------------------------
print("[Fig 4] themes + focus boxplots + Dunn ...")
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
themes = {k: [c for c in v if c in df_main.columns] for k, v in themes.items()}
theme_z = pd.DataFrame(index=df_main.index)
for t, cols in themes.items(): theme_z[t] = df_main[cols].mean(axis=1)
theme_z["ecotype"] = df_main["ecotype"]
tkw = pd.read_csv(OUT / "step7_theme_KW.tsv", sep="\t")

fig, axes = plt.subplots(1, len(themes), figsize=(4.2*len(themes), 4.6))
for ax, (tname, _) in zip(axes, themes.items()):
    sns.boxplot(data=theme_z, x="ecotype", y=tname, order=EC_ORDER, palette=PAL_EC, ax=ax, showfliers=False)
    sns.stripplot(data=theme_z, x="ecotype", y=tname, order=EC_ORDER, color="black", size=2, alpha=.3, ax=ax)
    q = float(tkw.loc[tkw["theme"]==tname, "q_BH"].values[0])
    eps = float(tkw.loc[tkw["theme"]==tname, "epsilon_sq"].values[0])
    star = "***" if q<1e-4 else ("**" if q<1e-2 else ("*" if q<5e-2 else "ns"))
    ax.set_title(f"{tname}\nKW q={q:.1e} {star} · ε²={eps:.2f}", fontsize=10)
    ax.set_xlabel(""); ax.set_ylabel("theme composite z")
    ax.tick_params(axis="x", rotation=18, labelsize=9)
plt.tight_layout(); fig.savefig(DST / "Fig4A_theme_composite_boxplot.png", dpi=DPI); plt.close()

sigs = []
for v in themes.values(): sigs += v
sig_means = df_main.groupby("ecotype")[sigs].mean().T.reindex(columns=EC_ORDER)
sig_means.index = [c.replace("ssGSEA_", "") for c in sig_means.index]
fig, ax = plt.subplots(figsize=(6, 0.45*len(sig_means)+1.0))
sns.heatmap(sig_means, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
            cbar_kws=dict(label="mean z within ecotype"), ax=ax)
pos = 0
for tname, cols in themes.items():
    pos += len(cols)
    if pos < len(sig_means): ax.axhline(pos, color="black", lw=1.0)
ax.set_title("Signature × ecotype")
plt.tight_layout(); fig.savefig(DST / "Fig4B_signature_heatmap.png", dpi=DPI); plt.close()

theme_z["cohort_group"] = df_main["cohort_group"]
cohort_order = [c for c in ["DMG_K27", "DHG_G34", "pHGG_WT", "IHG"] if c in theme_z["cohort_group"].unique()]
fig, axes = plt.subplots(1, len(themes), figsize=(4.2*len(themes), 4.6))
for ax, (t, _) in zip(axes, themes.items()):
    sns.boxplot(data=theme_z, x="cohort_group", y=t, hue="ecotype", order=cohort_order,
                hue_order=EC_ORDER, palette=PAL_EC, ax=ax, showfliers=False)
    ax.set_title(t, fontsize=10); ax.set_xlabel(""); ax.set_ylabel("theme z")
    ax.tick_params(axis="x", rotation=15, labelsize=8)
    if ax is not axes[-1] and ax.get_legend() is not None: ax.get_legend().remove()
    else:
        if ax is axes[-1]: ax.legend(loc="upper right", fontsize=7, frameon=False)
plt.tight_layout(); fig.savefig(DST / "Fig4C_theme_by_cohort.png", dpi=DPI); plt.close()

# Focus boxplots
kw = pd.read_csv(OUT / "step6_KW_results.tsv", sep="\t")
focus = [
    ("LM22_T cells CD8",                  "LM22 CD8 T cells"),
    ("LM22_Macrophages M0",               "LM22 M0"),
    ("LM22_Macrophages M1",               "LM22 M1"),
    ("LM22_Macrophages M2",               "LM22 M2"),
    ("LM22_Monocytes",                    "LM22 Monocytes"),
    ("LM22_Neutrophils",                  "LM22 Neutrophils"),
    ("ssGSEA_Microglia_Klemm2020",        "ssGSEA Microglia (Klemm 2020)"),
    ("ssGSEA_MDM_Klemm2020",              "ssGSEA MDM (Klemm 2020)"),
    ("ssGSEA_MgTAM_Antunes2021",          "ssGSEA Mg-TAM (Antunes 2021)"),
    ("ssGSEA_MoTAM_Antunes2021",          "ssGSEA Mo-TAM (Antunes 2021)"),
    ("ssGSEA_T_Cell_Cytotoxicity",        "ssGSEA T-cell cytotoxicity"),
    ("ssGSEA_IFN_Gamma_Response",         "ssGSEA IFN-γ"),
]
focus = [(c, l) for c, l in focus if c in df_main.columns]
fig, axes = plt.subplots(3, 4, figsize=(15, 11))
for ax, (col, lab) in zip(axes.flat, focus):
    sns.boxplot(data=df_main, x="ecotype", y=col, order=EC_ORDER, palette=PAL_EC, ax=ax, showfliers=False)
    sns.stripplot(data=df_main, x="ecotype", y=col, order=EC_ORDER, color="black", size=2, alpha=.35, ax=ax)
    q = kw.loc[kw["feature"]==col, "q_BH"]
    qv = float(q.values[0]) if len(q) else float("nan")
    star = "***" if qv<1e-4 else ("**" if qv<1e-2 else ("*" if qv<5e-2 else "ns"))
    ax.set_title(f"{lab}\nKW q={qv:.1e} {star}", fontsize=9)
    ax.set_xlabel(""); ax.set_ylabel("z / fraction"); ax.tick_params(axis="x", rotation=15, labelsize=8)
for ax in axes.flat[len(focus):]: ax.set_visible(False)
plt.tight_layout(); fig.savefig(DST / "Fig4D_focus_boxplots.png", dpi=DPI); plt.close()

# Dunn heatmap (already-improved version)
wide = pd.read_csv(OUT / "step6_Dunn_signed_log10q_wide.tsv", sep="\t", index_col=0)
preferred = ["Immune-desert vs Lymphocyte-inflamed", "Immune-desert vs Myeloid-dominant", "Lymphocyte-inflamed vs Myeloid-dominant"]
cols = [c for c in preferred if c in wide.columns]
feats = [c for c, _ in focus if c in wide.index]
labels = [l for c, l in focus if c in wide.index]
mat = wide.loc[feats, cols].copy(); mat.index = labels
fig, ax = plt.subplots(figsize=(max(8, 3.6*len(cols)), 1.0 + 0.65*len(labels)))
vmax = max(float(np.nanmax(np.abs(mat.values))), 5)
sns.heatmap(mat, annot=True, fmt=".1f", annot_kws={"size":11}, cmap="RdBu_r", center=0,
            vmin=-vmax, vmax=vmax, linewidths=0.5, linecolor="white",
            cbar_kws=dict(label="signed -log10(q) Dunn BH", shrink=.8), ax=ax)
ax.set_title("Dunn post-hoc — focus signatures", fontsize=12, pad=12)
ax.set_xlabel(""); ax.set_ylabel("")
plt.setp(ax.get_xticklabels(), rotation=18, ha="right", fontsize=10)
plt.setp(ax.get_yticklabels(), rotation=0, fontsize=10)
plt.tight_layout(); fig.savefig(DST / "Fig4E_Dunn_focus.png", dpi=DPI, bbox_inches="tight"); plt.close()

# ---------------------------------------------------------------------------
# Fig 5: volcanos + enrichment
# ---------------------------------------------------------------------------
print("[Fig 5] volcanos + enrichment ...")
ec_keys = [("Lymphocyte-inflamed","Lymphocyte_inflamed"),
           ("Myeloid-dominant","Myeloid_dominant"),
           ("Immune-desert","Immune_desert")]
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
for ax, (ec, fname) in zip(axes, ec_keys):
    d = pd.read_csv(OUT / f"step8_DEG_{fname}.tsv", sep="\t")
    d["nlq"] = -np.log10(np.clip(d["q_BH"], 1e-300, 1))
    sig = (d["q_BH"] < 0.05) & (d["log2FC"].abs() > 1)
    ax.scatter(d.loc[~sig, "log2FC"], d.loc[~sig, "nlq"], s=3, c="lightgray", alpha=.5, edgecolor="none")
    ax.scatter(d.loc[sig, "log2FC"],  d.loc[sig, "nlq"],  s=4, c=PAL_EC[ec], alpha=.7, edgecolor="none")
    top = d[(d["q_BH"]<0.05)&(d["log2FC"]>0)].sort_values("log2FC", ascending=False).head(8)
    for _, r in top.iterrows(): ax.text(r["log2FC"], r["nlq"], r["gene"], fontsize=6.5)
    ax.axhline(-np.log10(0.05), color="black", ls="--", alpha=.3)
    ax.axvline(1,  color="black", ls="--", alpha=.3); ax.axvline(-1, color="black", ls="--", alpha=.3)
    n_up = ((d["q_BH"]<0.05)&(d["log2FC"]>1)).sum(); n_dn = ((d["q_BH"]<0.05)&(d["log2FC"]<-1)).sum()
    ax.set_title(f"{ec}\n↑{n_up}  ↓{n_dn}", fontsize=10)
    ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10(q BH)")
plt.tight_layout(); fig.savefig(DST / "Fig5A_volcanos.png", dpi=DPI); plt.close()

fig, axes = plt.subplots(1, 3, figsize=(22, 5))
for ax, (ec, fname) in zip(axes, ec_keys):
    f = OUT / f"step8_gprofiler_{fname}.tsv"
    if not f.exists(): ax.set_visible(False); continue
    e = pd.read_csv(f, sep="\t")
    e = e[e["source"].isin(["GO:BP","REAC","KEGG","WP"])].sort_values("p_value").head(15)
    if e.empty: ax.set_visible(False); continue
    nl = -np.log10(np.clip(e["p_value"], 1e-300, 1))
    labels = (e["term_name"].str.slice(0,56) + " [" + e["source"] + "]").values
    ax.barh(labels[::-1], nl[::-1], color=PAL_EC[ec])
    ax.set_xlabel("-log10(adjusted p, g:SCS)")
    ax.set_title(f"{ec}", fontsize=10); ax.tick_params(axis="y", labelsize=8)
plt.tight_layout(); fig.savefig(DST / "Fig5B_top_enrichment.png", dpi=DPI); plt.close()

# ---------------------------------------------------------------------------
# Fig 6: KM (overall + per cohort) + Cox forest
# ---------------------------------------------------------------------------
print("[Fig 6] KM + Cox ...")
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
surv = ann_full.join(ec3[["ecotype"]], how="inner").copy()
surv["event"] = (surv["OS_status"].astype(str).str.upper()=="DECEASED").astype(int)
surv["time"]  = pd.to_numeric(surv["OS_days"], errors="coerce")
surv = surv.dropna(subset=["time","event","ecotype"])
lrt = multivariate_logrank_test(surv["time"], surv["ecotype"], surv["event"])

fig, ax = plt.subplots(figsize=(6.5, 4.8))
for ec in EC_ORDER:
    g = surv[surv["ecotype"]==ec]
    KaplanMeierFitter().fit(g["time"], g["event"],
        label=f"{ec} (n={len(g)}, ev={int(g['event'].sum())})").plot_survival_function(
        ax=ax, color=PAL_EC[ec], ci_show=False)
ax.set_title(f"Overall survival by ecotype (log-rank p={float(lrt.p_value):.2e})")
ax.set_xlabel("OS (days)"); ax.set_ylabel("Survival probability"); ax.set_ylim(-0.02, 1.02)
ax.legend(loc="best", frameon=False, fontsize=8)
plt.tight_layout(); fig.savefig(DST / "Fig6A_KM_overall.png", dpi=DPI); plt.close()

cohort_order = ["DMG_K27","DHG_G34","pHGG_WT","IHG"]
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for ax, cg in zip(axes.flat, cohort_order):
    sub = surv[surv["cohort_group"]==cg]
    if len(sub) < 5: ax.set_visible(False); continue
    for ec in EC_ORDER:
        g = sub[sub["ecotype"]==ec]
        if g.empty or g["event"].sum() == 0: continue
        KaplanMeierFitter().fit(g["time"], g["event"],
            label=f"{ec} (n={len(g)}, ev={int(g['event'].sum())})").plot_survival_function(
            ax=ax, color=PAL_EC[ec], ci_show=False)
    if sub["ecotype"].nunique() >= 2:
        lr = multivariate_logrank_test(sub["time"], sub["ecotype"], sub["event"])
        p = float(lr.p_value)
    else:
        p = float("nan")
    ax.set_title(f"{cg}  (n={len(sub)}, log-rank p={p:.2e})", fontsize=10)
    ax.set_xlabel("OS (days)"); ax.set_ylabel("Survival prob"); ax.set_ylim(-0.02, 1.02)
    ax.legend(loc="best", frameon=False, fontsize=7)
plt.tight_layout(); fig.savefig(DST / "Fig6B_KM_by_cohort.png", dpi=DPI); plt.close()

sm = pd.read_csv(OUT / "step9_cox_multivariable.tsv", sep="\t")
sm["lab"] = (sm["covariate"]
             .str.replace("ecotype_","ecotype: ").str.replace("cohort_group_","cohort: ")
             .str.replace("_"," "))
fig, ax = plt.subplots(figsize=(7.5, 0.4*len(sm)+1.5))
for i, r in sm.iterrows():
    color = "red" if r["p"] < 0.05 else "black"
    ax.errorbar(r["exp(coef)"], i,
                xerr=[[r["exp(coef)"]-r["exp(coef) lower 95%"]],
                      [r["exp(coef) upper 95%"]-r["exp(coef)"]]],
                fmt="o", color=color, ecolor=color, capsize=3)
ax.axvline(1, color="grey", ls="--", alpha=.5)
ax.set_yticks(range(len(sm))); ax.set_yticklabels(sm["lab"], fontsize=9)
ax.set_xscale("log"); ax.set_xlabel("Hazard ratio (95% CI)")
ax.set_title("Multivariable Cox  (ref: Immune-desert, pHGG_WT, Female)")
plt.tight_layout(); fig.savefig(DST / "Fig6C_Cox_forest.png", dpi=DPI); plt.close()

# ---------------------------------------------------------------------------
# Fig 7: BRAF projection
# ---------------------------------------------------------------------------
print("[Fig 7] BRAF projection ...")
dist_table = pd.read_csv(OUT / "step10_ecotype_distribution_by_cohort_family.tsv", sep="\t", index_col=0)
dist_pct = dist_table.div(dist_table.sum(axis=0), axis=1)*100
from scipy.stats import chi2_contingency
chi2, p_chi, _, _ = chi2_contingency(dist_table.values)
fig, ax = plt.subplots(figsize=(7, 4.5))
bottom = np.zeros(len(dist_pct.columns))
for ec in EC_ORDER:
    if ec in dist_pct.index:
        vals = dist_pct.loc[ec].values
        ax.bar(dist_pct.columns, vals, bottom=bottom, color=PAL_EC[ec], label=ec); bottom += vals
ax.set_ylabel("% of samples"); ax.set_ylim(0, 100)
ax.set_title(f"Ecotype by cohort family  (chi² p={p_chi:.2e})")
ax.legend(loc="upper right", fontsize=8, frameon=False)
plt.tight_layout(); fig.savefig(DST / "Fig7A_ecotype_by_family.png", dpi=DPI); plt.close()

braf = pd.read_csv(OUT / "step10_BRAF_ALT_ecotype_assignment.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
sets = {
    "Main_n349":     list(ec3.index.intersection(merged.index)),
    "BRAF_ALT_main": list(braf[braf["set"]=="BRAF_ALT_main"].index),
    "BRAF_ALT_LGG":  list(braf[braf["set"]=="BRAF_ALT_LGG"].index),
}
recs = []
for s, ids in sets.items():
    vals = merged.loc[merged.index.intersection(ids), "ssGSEA_MAPK_Activity"].dropna()
    recs.append(pd.DataFrame({"set": s, "MAPK_Activity": vals.values}))
mapk_df = pd.concat(recs, ignore_index=True)
from scipy.stats import kruskal
H, p_kw = kruskal(*[g["MAPK_Activity"].values for _, g in mapk_df.groupby("set")])
fig, ax = plt.subplots(figsize=(6.5, 4.5))
sns.boxplot(data=mapk_df, x="set", y="MAPK_Activity",
            order=["Main_n349","BRAF_ALT_main","BRAF_ALT_LGG"],
            palette={"Main_n349":"#0EA5E9","BRAF_ALT_main":"#F59E0B","BRAF_ALT_LGG":"#10B981"},
            ax=ax, showfliers=False)
sns.stripplot(data=mapk_df, x="set", y="MAPK_Activity",
              order=["Main_n349","BRAF_ALT_main","BRAF_ALT_LGG"],
              color="black", size=2, alpha=.3, ax=ax)
ax.set_xlabel(""); ax.set_ylabel("ssGSEA MAPK_Activity (z)")
ax.set_title(f"MAPK activity across cohort families  (KW p={p_kw:.1e})")
plt.tight_layout(); fig.savefig(DST / "Fig7B_MAPK_by_family.png", dpi=DPI); plt.close()

# PCA overlay
from sklearn.decomposition import PCA
feat = pd.read_csv(OUT / "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t", index_col=0)
main_ids = list(ec3.index.intersection(feat.index))
main_feat = feat.loc[main_ids].dropna(axis=0, how="any")
pca = PCA(n_components=2).fit(main_feat.values)
main_pc = pca.transform(main_feat.values)
ec_arr = ec3.loc[main_feat.index, "ecotype"].values
braf_main_pc = pca.transform(feat.loc[[i for i in braf[braf["set"]=="BRAF_ALT_main"].index if i in feat.index]].values)
braf_lgg_pc  = pca.transform(feat.loc[[i for i in braf[braf["set"]=="BRAF_ALT_LGG"].index  if i in feat.index]].values)
fig, ax = plt.subplots(figsize=(8, 6))
for ec, col in PAL_EC.items():
    m = ec_arr == ec
    ax.scatter(main_pc[m, 0], main_pc[m, 1], c=col, s=16, alpha=.6, label=f"Main · {ec}")
ax.scatter(braf_main_pc[:, 0], braf_main_pc[:, 1], facecolors="none", edgecolors="#F59E0B", s=42, lw=1.5,
           label="BRAF_ALT_main")
ax.scatter(braf_lgg_pc[:, 0],  braf_lgg_pc[:, 1],  facecolors="none", edgecolors="#10B981", s=42, lw=1.5,
           label="BRAF_ALT_LGG")
ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.set_title("Main PCA space — BRAF_ALT overlay")
ax.legend(loc="best", fontsize=8, frameon=False)
plt.tight_layout(); fig.savefig(DST / "Fig7C_BRAF_overlay.png", dpi=DPI); plt.close()

# ---------------------------------------------------------------------------
# Suppl S3 + S4 (single combined heatmaps already exist, regenerate at 300 dpi from npz)
# ---------------------------------------------------------------------------
print("[Suppl] S3 / S4 ...")
import subprocess
for src, dst_name in [
    (BASE / "scripts" / "suppl_S3_highconf_subset.py", "SupplS3"),
    (BASE / "scripts" / "suppl_S4_full_cohort_sensitivity.py", "SupplS4"),
]:
    txt = src.read_text()
    txt2 = txt.replace("dpi=200", f"dpi={DPI}").replace("FIG  = OUT / \"figs_suppl\"",
                                                         f"FIG  = OUT / \"figs_300dpi\"")
    tmp = Path(f"/tmp/_{dst_name}.py")
    tmp.write_text(txt2)
    try:
        subprocess.run([sys.executable, str(tmp)], check=True, capture_output=True, timeout=40)
    except Exception as e:
        print(f"  {dst_name} regen failed: {e}")

print("\n=== Verify DPI of all regenerated figures ===")
from PIL import Image
for p in sorted(DST.glob("*.png")):
    img = Image.open(p)
    print(f"  {p.name}  {img.size}  dpi={img.info.get('dpi')}")
