"""
Step 4c — Extended ssGSEA: Fecci/Sampson review-derived modules
================================================================
A1 Exhaustion / Senescence / Naive-Tcm-loss
A2 MDSC (G-MDSC / M-MDSC)
A3 Hypoxia / HIF
A4 Antigen presentation HLA-I / HLA-II  (+ HLA-I/II ratio)
A6 CSF1R-TAM axis / CCR2 chemokine axis
+   Treg signature, Suppressive cytokine axis

Input  : output/tpm_for_cibersortx.tsv  (GeneSymbol x sample TPM)
         output/ecotype_assignment_k3_annotated.tsv (BSID -> ecotype, cohort_group ...)
         data/fecci_sampson_extended_modules.gmt
Output : output/step4c_ssGSEA_scores.tsv
         output/step4c_KW_by_ecotype.tsv
         output/step4c_Dunn_posthoc.tsv
         output/step4c_KW_by_cohort.tsv
         output/step4c_module_summary.json
         output/figs_step4c/*.png  (300 dpi)
"""
from __future__ import annotations
import json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import gseapy as gp
from scipy import stats
import scikit_posthocs as sp
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT     = Path("/sessions/blissful-gifted-dirac/mnt/Open PBTA")
OUT      = ROOT / "output"
FIGDIR   = OUT / "figs_step4c"
FIGDIR.mkdir(exist_ok=True, parents=True)

TPM_FN   = OUT / "tpm_for_cibersortx.tsv"
ECO_FN   = OUT / "ecotype_LM22_main_k3_annotated.tsv"   # canonical labels
META_FN  = OUT / "ecotype_assignment_k3_annotated.tsv"  # cohort/location metadata
GMT_FN   = ROOT / "data" / "fecci_sampson_extended_modules.gmt"

ECOTYPE_ORDER = ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]
ECOTYPE_COLORS = {
    "Lymphocyte-inflamed": "#2C7BB6",
    "Myeloid-dominant":    "#D7301F",
    "Immune-desert":       "#7F7F7F",
}
COHORT_ORDER = ["DMG_K27", "DHG_G34", "pHGG_WT", "IHG"]
COHORT_COLORS = {
    "DMG_K27": "#1B9E77",
    "DHG_G34": "#D95F02",
    "pHGG_WT": "#7570B3",
    "IHG":     "#E7298A",
}

def load_gmt(fn):
    sets = {}
    with open(fn) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            name, _, *genes = parts
            sets[name] = [g for g in genes if g]
    return sets

def cliffs_delta(a, b):
    a = np.asarray(a); b = np.asarray(b)
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return np.nan
    A = np.tile(a, (n2, 1)).T
    B = np.tile(b, (n1, 1))
    return ((A > B).sum() - (A < B).sum()) / (n1 * n2)

def epsilon_squared(H, n):
    # KW epsilon^2 effect size
    return (H - (1)) / (n - 1) if n > 1 else np.nan
def epsilon_squared_kw(H, k, n):
    # better formula: eps2 = (H - k + 1) / (n - k)
    return max(0.0, (H - k + 1) / (n - k)) if n - k > 0 else np.nan

def main():
    print("[*] Loading TPM matrix ...", flush=True)
    tpm = pd.read_csv(TPM_FN, sep="\t")
    tpm = tpm.rename(columns={tpm.columns[0]: "GeneSymbol"})
    tpm["GeneSymbol"] = tpm["GeneSymbol"].astype(str)
    tpm = tpm.drop_duplicates(subset="GeneSymbol").set_index("GeneSymbol")
    print(f"    matrix shape = {tpm.shape}", flush=True)

    print("[*] Loading ecotype assignment (canonical labels) ...", flush=True)
    eco = pd.read_csv(ECO_FN, sep="\t").set_index("Kids_First_Biospecimen_ID")
    meta = pd.read_csv(META_FN).set_index("Kids_First_Biospecimen_ID")
    eco = eco[["ecotype"]].join(
        meta[["cohort_group", "location_class", "age_dev_group"]], how="left"
    )
    keep = [b for b in tpm.columns if b in eco.index]
    tpm = tpm[keep]
    eco = eco.loc[keep]
    print(f"    n samples with ecotype = {len(keep)}", flush=True)
    print(f"    ecotype labels = {eco.ecotype.value_counts().to_dict()}", flush=True)

    print("[*] Loading gene sets ...", flush=True)
    sets = load_gmt(GMT_FN)
    for n, gs in sets.items():
        hit = sum(g in tpm.index for g in gs)
        print(f"    {n}: {hit}/{len(gs)} genes present", flush=True)

    print("[*] Computing ssGSEA via gseapy ...", flush=True)
    log_tpm = np.log2(tpm.astype(float) + 1.0)
    # gseapy.ssgsea expects gene x sample matrix
    res = gp.ssgsea(data=log_tpm, gene_sets=sets, sample_norm_method="rank",
                    no_plot=True, processes=1, min_size=3, max_size=500,
                    permutation_num=0, outdir=None)
    # gseapy returns res.res2d with Term, Name (=sample), ES, NES
    scores = res.res2d.copy()
    scores["NES"] = pd.to_numeric(scores["NES"], errors="coerce")
    scores = scores.pivot(index="Name", columns="Term", values="NES").astype(float)
    scores.index.name = "Kids_First_Biospecimen_ID"
    # Compute HLA-I/HLA-II ratio (using NES means; use raw NES centered to avoid divide-by-zero)
    if "Antigen_presentation_HLA_I" in scores and "Antigen_presentation_HLA_II" in scores:
        scores["HLA_I_minus_II"] = scores["Antigen_presentation_HLA_I"] - scores["Antigen_presentation_HLA_II"]
    scores.to_csv(OUT / "step4c_ssGSEA_scores.tsv", sep="\t")
    print(f"[+] saved scores -> {OUT/'step4c_ssGSEA_scores.tsv'}", flush=True)

    # ---- merge with metadata
    scores = scores.join(eco[["ecotype", "cohort_group", "location_class", "age_dev_group"]], how="inner")
    scores = scores.dropna(subset=["ecotype"])

    # ---- KW + Dunn per module across ecotypes
    print("[*] Running KW + Dunn (3 ecotype) ...", flush=True)
    META_COLS = {"ecotype", "cohort_group", "location_class", "age_dev_group"}
    modules = [c for c in scores.columns if c not in META_COLS]
    kw_rows, dunn_rows = [], []
    for m in modules:
        groups = [scores.loc[scores.ecotype == e, m].dropna().values for e in ECOTYPE_ORDER]
        H, p = stats.kruskal(*groups)
        eps2 = epsilon_squared_kw(H, k=3, n=sum(len(g) for g in groups))
        kw_rows.append({"module": m, "H": H, "p": p, "eps2": eps2,
                        "median_Lymph": np.median(groups[0]),
                        "median_Mye":   np.median(groups[1]),
                        "median_Desert":np.median(groups[2])})
        # Dunn posthoc
        sub = scores[["ecotype", m]].dropna()
        sub = sub[sub.ecotype.isin(ECOTYPE_ORDER)]
        dunn = sp.posthoc_dunn(sub, val_col=m, group_col="ecotype", p_adjust="fdr_bh")
        for i, g1 in enumerate(ECOTYPE_ORDER):
            for g2 in ECOTYPE_ORDER[i+1:]:
                q = dunn.loc[g1, g2]
                d = cliffs_delta(scores.loc[scores.ecotype == g1, m].dropna(),
                                 scores.loc[scores.ecotype == g2, m].dropna())
                dunn_rows.append({"module": m, "group1": g1, "group2": g2,
                                  "q_BH": q, "cliffs_delta": d})

    kw = pd.DataFrame(kw_rows)
    kw["q_BH"] = multipletests(kw["p"], method="fdr_bh")[1]
    kw = kw.sort_values("eps2", ascending=False).reset_index(drop=True)
    kw.to_csv(OUT / "step4c_KW_by_ecotype.tsv", sep="\t", index=False)
    dunn_df = pd.DataFrame(dunn_rows)
    dunn_df.to_csv(OUT / "step4c_Dunn_posthoc.tsv", sep="\t", index=False)

    # ---- KW per module across cohort_group
    print("[*] Running KW (4 cohort) ...", flush=True)
    cohort_rows = []
    for m in modules:
        groups = [scores.loc[scores.cohort_group == c, m].dropna().values for c in COHORT_ORDER]
        H, p = stats.kruskal(*groups)
        eps2 = epsilon_squared_kw(H, k=4, n=sum(len(g) for g in groups))
        cohort_rows.append({"module": m, "H": H, "p": p, "eps2": eps2,
                            **{f"median_{c}": np.median(g) for c, g in zip(COHORT_ORDER, groups)}})
    cohort_kw = pd.DataFrame(cohort_rows)
    cohort_kw["q_BH"] = multipletests(cohort_kw["p"], method="fdr_bh")[1]
    cohort_kw = cohort_kw.sort_values("eps2", ascending=False).reset_index(drop=True)
    cohort_kw.to_csv(OUT / "step4c_KW_by_cohort.tsv", sep="\t", index=False)

    # ---- Figures (300 dpi)
    plt.rcParams.update({"figure.dpi": 300, "savefig.dpi": 300, "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})

    print("[*] Drawing per-module boxplots (3 ecotype) ...", flush=True)
    sub = scores[scores.ecotype.isin(ECOTYPE_ORDER)].copy()
    n_mod = len(modules)
    ncols = 4
    nrows = int(np.ceil(n_mod / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*3.0, nrows*2.6))
    axes = np.array(axes).reshape(nrows, ncols)
    for idx, m in enumerate(modules):
        ax = axes[idx // ncols, idx % ncols]
        sns.boxplot(data=sub, x="ecotype", y=m, order=ECOTYPE_ORDER,
                    palette=ECOTYPE_COLORS, ax=ax, fliersize=1.5,
                    linewidth=0.7, showcaps=True)
        ax.set_title(m, fontsize=9)
        ax.set_xlabel(""); ax.set_ylabel("ssGSEA NES", fontsize=8)
        ax.tick_params(axis='x', rotation=20, labelsize=7)
        # add p annotation
        kw_row = kw[kw.module == m].iloc[0]
        ax.text(0.97, 0.95, f"q={kw_row.q_BH:.1e}\nε²={kw_row.eps2:.2f}",
                ha="right", va="top", transform=ax.transAxes, fontsize=7,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.8", lw=0.5))
    for j in range(n_mod, nrows*ncols):
        axes[j // ncols, j % ncols].axis("off")
    fig.suptitle("Extended brain-immune ssGSEA modules vs immune ecotype (n=%d)" % len(sub),
                 fontsize=11, y=1.0)
    fig.tight_layout()
    fig.savefig(FIGDIR / "step4c_boxplot_ecotype.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / "step4c_boxplot_ecotype.pdf", bbox_inches="tight")
    plt.close(fig)

    print("[*] Drawing per-module boxplots (4 cohort) ...", flush=True)
    sub2 = scores[scores.cohort_group.isin(COHORT_ORDER)].copy()
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*3.0, nrows*2.6))
    axes = np.array(axes).reshape(nrows, ncols)
    for idx, m in enumerate(modules):
        ax = axes[idx // ncols, idx % ncols]
        sns.boxplot(data=sub2, x="cohort_group", y=m, order=COHORT_ORDER,
                    palette=COHORT_COLORS, ax=ax, fliersize=1.5,
                    linewidth=0.7)
        ax.set_title(m, fontsize=9)
        ax.set_xlabel(""); ax.set_ylabel("ssGSEA NES", fontsize=8)
        ax.tick_params(axis='x', rotation=20, labelsize=7)
        ck = cohort_kw[cohort_kw.module == m].iloc[0]
        ax.text(0.97, 0.95, f"q={ck.q_BH:.1e}\nε²={ck.eps2:.2f}",
                ha="right", va="top", transform=ax.transAxes, fontsize=7,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.8", lw=0.5))
    for j in range(n_mod, nrows*ncols):
        axes[j // ncols, j % ncols].axis("off")
    fig.suptitle("Extended brain-immune ssGSEA modules vs cohort group (n=%d)" % len(sub2),
                 fontsize=11, y=1.0)
    fig.tight_layout()
    fig.savefig(FIGDIR / "step4c_boxplot_cohort.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / "step4c_boxplot_cohort.pdf", bbox_inches="tight")
    plt.close(fig)

    # Heatmap of mean NES per (ecotype, module) and (cohort, module)
    print("[*] Drawing heatmaps ...", flush=True)
    sub_num = sub.copy(); sub_num[modules] = sub_num[modules].apply(pd.to_numeric, errors="coerce")
    sub2_num = sub2.copy(); sub2_num[modules] = sub2_num[modules].apply(pd.to_numeric, errors="coerce")
    eco_mean = sub_num.groupby("ecotype")[modules].mean().T.loc[modules, ECOTYPE_ORDER].astype(float)
    coh_mean = sub2_num.groupby("cohort_group")[modules].mean().T.loc[modules, COHORT_ORDER].astype(float)
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 0.32*len(modules)+1.0))
    sns.heatmap(eco_mean, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                cbar_kws={"label":"mean NES"}, ax=axes[0],
                linewidths=0.3, linecolor="white", annot_kws={"size":7})
    axes[0].set_title("Mean ssGSEA NES — ecotype", fontsize=10)
    axes[0].set_xlabel(""); axes[0].set_ylabel("")
    sns.heatmap(coh_mean, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                cbar_kws={"label":"mean NES"}, ax=axes[1],
                linewidths=0.3, linecolor="white", annot_kws={"size":7})
    axes[1].set_title("Mean ssGSEA NES — cohort", fontsize=10)
    axes[1].set_xlabel(""); axes[1].set_ylabel("")
    fig.tight_layout()
    fig.savefig(FIGDIR / "step4c_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / "step4c_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---- summary JSON
    summary = {
        "n_samples": int(len(scores)),
        "modules": modules,
        "ecotype_counts": sub.ecotype.value_counts().to_dict(),
        "cohort_counts":  sub2.cohort_group.value_counts().to_dict(),
        "top_eco_modules_by_eps2": kw.head(6)[["module","eps2","q_BH"]].to_dict(orient="records"),
        "top_cohort_modules_by_eps2": cohort_kw.head(6)[["module","eps2","q_BH"]].to_dict(orient="records"),
        "HLA_I_minus_II_per_ecotype": sub.groupby("ecotype")["HLA_I_minus_II"].mean().to_dict() if "HLA_I_minus_II" in modules else None,
      }
    with open(OUT / "step4c_module_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)

    print("[*] Top 5 modules by ecotype-discriminating eps^2:")
    print(kw.head(5).to_string(index=False))
    print("\n[*] Top 5 modules by cohort-discriminating eps^2:")
    print(cohort_kw.head(5).to_string(index=False))
    print("\n[*] Mean NES per ecotype (selected):")
    print(eco_mean.round(3))
    print("[+] DONE.")

if __name__ == "__main__":
    main()
