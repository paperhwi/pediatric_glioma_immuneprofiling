"""
Step 4d — Levine 2024 (Nat Commun) integration analysis
========================================================
B1  TIS (Danaher 18-gene tumor inflammation signature)
B3  CD276 (B7-H3) standalone axis — CAR-T candidate identification
B5  Bagaev 2021 pan-cancer microenvironment signatures (24 modules)
B8  7-checkpoint co-expression network per ecotype

Inputs : output/tpm_for_cibersortx.tsv
         output/ecotype_LM22_main_k3_annotated.tsv (canonical ecotype labels)
         output/ecotype_assignment_k3_annotated.tsv (cohort/location metadata)
         data/levine2024_bagaev_modules.gmt
Outputs: output/step4d_TIS_per_sample.tsv
         output/step4d_bagaev_ssGSEA.tsv
         output/step4d_KW_by_ecotype.tsv
         output/step4d_KW_by_cohort.tsv
         output/step4d_checkpoint_correlations.tsv
         output/step4d_summary.json
         output/figs_step4d/*.png|pdf (300 dpi)
"""
from __future__ import annotations
import json
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
from matplotlib.patches import Patch
import seaborn as sns

ROOT     = Path("/sessions/blissful-gifted-dirac/mnt/Open PBTA")
OUT      = ROOT / "output"
FIGDIR   = OUT / "figs_step4d"
FIGDIR.mkdir(exist_ok=True, parents=True)

TPM_FN   = OUT / "tpm_for_cibersortx.tsv"
ECO_FN   = OUT / "ecotype_LM22_main_k3_annotated.tsv"
META_FN  = OUT / "ecotype_assignment_k3_annotated.tsv"
GMT_FN   = ROOT / "data" / "levine2024_bagaev_modules.gmt"
PRIOR_SCORES = OUT / "step4c_ssGSEA_scores.tsv"  # for cross-reference with Step 4c

ECOTYPE_ORDER  = ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]
ECOTYPE_COLORS = {"Lymphocyte-inflamed":"#2C7BB6","Myeloid-dominant":"#D7301F","Immune-desert":"#7F7F7F"}
COHORT_ORDER   = ["DMG_K27","DHG_G34","pHGG_WT","IHG"]
COHORT_COLORS  = {"DMG_K27":"#1B9E77","DHG_G34":"#D95F02","pHGG_WT":"#7570B3","IHG":"#E7298A"}

TIS_GENES = ["CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
             "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
             "PSMB10","STAT1","TIGIT"]
CHECKPOINTS = ["PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT","CD276"]

def load_gmt(fn):
    out = {}
    for line in open(fn):
        n, _, *g = line.rstrip("\n").split("\t")
        out[n] = [x for x in g if x]
    return out

def eps2(H, k, n):
    return max(0.0, (H-k+1)/(n-k)) if n-k>0 else np.nan

def cliffs_delta(a, b):
    a = np.asarray(a); b = np.asarray(b); n1,n2 = len(a),len(b)
    if n1==0 or n2==0: return np.nan
    A = np.tile(a,(n2,1)).T; B = np.tile(b,(n1,1))
    return ((A>B).sum() - (A<B).sum()) / (n1*n2)

def main():
    print("[*] Loading TPM ...", flush=True)
    tpm = pd.read_csv(TPM_FN, sep="\t")
    tpm = tpm.rename(columns={tpm.columns[0]:"GeneSymbol"}).drop_duplicates("GeneSymbol").set_index("GeneSymbol")
    log_tpm = np.log2(tpm.astype(float) + 1.0)

    print("[*] Loading ecotype (canonical) + metadata ...", flush=True)
    eco = pd.read_csv(ECO_FN, sep="\t").set_index("Kids_First_Biospecimen_ID")
    meta = pd.read_csv(META_FN).set_index("Kids_First_Biospecimen_ID")
    eco = eco[["ecotype"]].join(meta[["cohort_group","location_class","age_dev_group"]], how="left")
    samples = [b for b in tpm.columns if b in eco.index]
    log_tpm = log_tpm[samples]
    eco = eco.loc[samples]
    print(f"    n = {len(samples)} samples; ecotype = {eco.ecotype.value_counts().to_dict()}", flush=True)

    # ===================== B1: TIS (Danaher 18-gene) =====================
    print("\n[B1] Computing TIS (Danaher 18-gene)", flush=True)
    tis_present = [g for g in TIS_GENES if g in log_tpm.index]
    print(f"    TIS genes present: {len(tis_present)}/{len(TIS_GENES)} -> {tis_present}", flush=True)
    TIS = log_tpm.loc[tis_present].mean(axis=0)
    TIS.name = "TIS_Danaher18"
    TIS.to_frame().join(eco).to_csv(OUT/"step4d_TIS_per_sample.tsv", sep="\t")

    # KW: TIS by ecotype and by cohort
    tis_eco = [TIS.loc[eco[eco.ecotype==e].index].dropna().values for e in ECOTYPE_ORDER]
    H_e, p_e = stats.kruskal(*tis_eco)
    tis_coh = [TIS.loc[eco[eco.cohort_group==c].index].dropna().values for c in COHORT_ORDER]
    H_c, p_c = stats.kruskal(*tis_coh)
    print(f"    TIS by ecotype: H={H_e:.1f} p={p_e:.2e} ε²={eps2(H_e,3,sum(map(len,tis_eco))):.3f}")
    print(f"    TIS by cohort : H={H_c:.1f} p={p_c:.2e} ε²={eps2(H_c,4,sum(map(len,tis_coh))):.3f}")
    print("    TIS median per ecotype:", {e: float(np.median(v)) for e,v in zip(ECOTYPE_ORDER,tis_eco)})
    print("    TIS median per cohort :", {c: float(np.median(v)) for c,v in zip(COHORT_ORDER,tis_coh)})

    # Correlate TIS with Step 4c modules (cross-reference)
    if PRIOR_SCORES.exists():
        prior = pd.read_csv(PRIOR_SCORES, sep="\t", index_col=0)
        prior_cols = [c for c in prior.columns if c not in {"ecotype","cohort_group","location_class","age_dev_group"}]
        tis_corr = {c: stats.spearmanr(TIS.loc[prior.index.intersection(TIS.index)],
                                        prior.loc[prior.index.intersection(TIS.index), c])[0]
                    for c in prior_cols}
        print("\n    TIS vs Step-4c modules (Spearman):")
        for k,v in sorted(tis_corr.items(), key=lambda x: -abs(x[1])):
            print(f"      {k:35s} ρ = {v:+.3f}")

    # ===================== B3: CD276 standalone axis =====================
    print("\n[B3] CD276 (B7-H3) standalone axis", flush=True)
    cp_expr = log_tpm.loc[[g for g in CHECKPOINTS if g in log_tpm.index]].T  # samples x genes
    cp_expr = cp_expr.join(eco["ecotype"], how="inner")

    # Overall Spearman correlation among checkpoints
    print("    Overall checkpoint correlation matrix (Spearman):")
    overall_corr = cp_expr[CHECKPOINTS].corr(method="spearman")
    print(overall_corr.round(2).to_string())
    overall_corr.to_csv(OUT/"step4d_checkpoint_correlations.tsv", sep="\t")

    # Per-ecotype correlation matrices
    eco_corrs = {}
    for e in ECOTYPE_ORDER:
        sub = cp_expr[cp_expr.ecotype==e][CHECKPOINTS]
        eco_corrs[e] = sub.corr(method="spearman")

    # CD276 vs ecotype KW
    cd276_eco = [cp_expr.loc[cp_expr.ecotype==e,"CD276"].dropna().values for e in ECOTYPE_ORDER]
    H_cd, p_cd = stats.kruskal(*cd276_eco)
    print(f"\n    CD276 by ecotype: H={H_cd:.1f} p={p_cd:.2e} ε²={eps2(H_cd,3,sum(map(len,cd276_eco))):.3f}")
    print("    CD276 median per ecotype:", {e: float(np.median(v)) for e,v in zip(ECOTYPE_ORDER,cd276_eco)})

    # Identify CD276-high / Lymph-low intersection
    cp_expr["TIS"] = TIS.reindex(cp_expr.index)
    cd276_hi_lymph_lo = cp_expr[(cp_expr.CD276 > cp_expr.CD276.quantile(0.75)) &
                                 (cp_expr.TIS  < cp_expr.TIS.quantile(0.25))]
    print(f"    CD276-high (Q4) ∩ TIS-low (Q1)  candidate samples: n = {len(cd276_hi_lymph_lo)}")
    cd276_hi_lymph_lo.to_csv(OUT/"step4d_CD276hi_TISlo_candidates.tsv", sep="\t")

    # ===================== B5: Bagaev 2021 ssGSEA =====================
    print("\n[B5] Bagaev 2021 pan-cancer signature ssGSEA", flush=True)
    sets = load_gmt(GMT_FN)
    for n, gs in sets.items():
        hit = sum(g in log_tpm.index for g in gs)
        print(f"    {n}: {hit}/{len(gs)} genes present")
    res = gp.ssgsea(data=log_tpm, gene_sets=sets, sample_norm_method="rank",
                    no_plot=True, threads=1, min_size=3, max_size=500,
                    permutation_num=0, outdir=None)
    bag = res.res2d.copy()
    bag["NES"] = pd.to_numeric(bag["NES"], errors="coerce")
    bag = bag.pivot(index="Name", columns="Term", values="NES").astype(float)
    bag.index.name = "Kids_First_Biospecimen_ID"
    bag.to_csv(OUT/"step4d_bagaev_ssGSEA.tsv", sep="\t")
    bag_with_meta = bag.join(eco, how="inner")

    modules = list(bag.columns)

    # KW per module by ecotype
    kw_rows = []
    for m in modules:
        g = [bag_with_meta.loc[bag_with_meta.ecotype==e, m].dropna().values for e in ECOTYPE_ORDER]
        H,p = stats.kruskal(*g)
        kw_rows.append({"module":m,"H":H,"p":p,
                        "eps2": eps2(H,3,sum(map(len,g))),
                        **{f"median_{e}": float(np.median(v)) for e,v in zip(ECOTYPE_ORDER,g)}})
    kw = pd.DataFrame(kw_rows)
    kw["q_BH"] = multipletests(kw.p, method="fdr_bh")[1]
    kw = kw.sort_values("eps2", ascending=False).reset_index(drop=True)
    kw.to_csv(OUT/"step4d_KW_by_ecotype.tsv", sep="\t", index=False)
    print("\n    Top 10 Bagaev modules by ecotype-eps2:")
    print(kw[["module","eps2","p","q_BH"]].head(10).to_string(index=False))

    # KW per module by cohort
    ckw_rows = []
    for m in modules:
        g = [bag_with_meta.loc[bag_with_meta.cohort_group==c, m].dropna().values for c in COHORT_ORDER]
        H,p = stats.kruskal(*g)
        ckw_rows.append({"module":m,"H":H,"p":p,
                         "eps2": eps2(H,4,sum(map(len,g))),
                         **{f"median_{c}": float(np.median(v)) for c,v in zip(COHORT_ORDER,g)}})
    ckw = pd.DataFrame(ckw_rows); ckw["q_BH"] = multipletests(ckw.p, method="fdr_bh")[1]
    ckw = ckw.sort_values("eps2", ascending=False).reset_index(drop=True)
    ckw.to_csv(OUT/"step4d_KW_by_cohort.tsv", sep="\t", index=False)

    # ===================== Figures (300 dpi) =====================
    plt.rcParams.update({"figure.dpi":300, "savefig.dpi":300, "font.size":9,
                         "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})

    # Fig 1: TIS by ecotype + by cohort
    print("\n[*] Figure 1: TIS distributions")
    fig, axes = plt.subplots(1,2,figsize=(7.0,3.2))
    df_t = TIS.to_frame("TIS").join(eco)
    sns.boxplot(data=df_t[df_t.ecotype.isin(ECOTYPE_ORDER)], x="ecotype", y="TIS",
                order=ECOTYPE_ORDER, hue="ecotype", palette=ECOTYPE_COLORS, legend=False, ax=axes[0],
                fliersize=2.0, linewidth=0.7)
    axes[0].set_title(f"TIS by ecotype\nKW H={H_e:.1f}, p={p_e:.1e}, ε²={eps2(H_e,3,sum(map(len,tis_eco))):.2f}", fontsize=9)
    axes[0].set_xlabel(""); axes[0].set_ylabel("TIS (mean log2 TPM, 18 genes)")
    axes[0].tick_params(axis="x", rotation=15, labelsize=8)
    sns.boxplot(data=df_t[df_t.cohort_group.isin(COHORT_ORDER)], x="cohort_group", y="TIS",
                order=COHORT_ORDER, hue="cohort_group", palette=COHORT_COLORS, legend=False, ax=axes[1],
                fliersize=2.0, linewidth=0.7)
    axes[1].set_title(f"TIS by cohort_group\nKW H={H_c:.1f}, p={p_c:.1e}, ε²={eps2(H_c,4,sum(map(len,tis_coh))):.2f}", fontsize=9)
    axes[1].set_xlabel(""); axes[1].set_ylabel("")
    axes[1].tick_params(axis="x", rotation=15, labelsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4d_TIS_boxplots.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4d_TIS_boxplots.pdf", bbox_inches="tight")
    plt.close(fig)

    # Fig 2: Checkpoint co-expression — overall + per ecotype
    print("[*] Figure 2: Checkpoint co-expression (overall + per ecotype)")
    fig, axes = plt.subplots(1,4,figsize=(15.0,3.6))
    sns.heatmap(overall_corr, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                vmin=-1, vmax=1, ax=axes[0], cbar=False, linewidths=0.3, linecolor="white",
                annot_kws={"size":7})
    axes[0].set_title("Overall (all samples)", fontsize=10)
    for ax, e in zip(axes[1:], ECOTYPE_ORDER):
        sns.heatmap(eco_corrs[e], annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                    vmin=-1, vmax=1, ax=ax, cbar=False, linewidths=0.3, linecolor="white",
                    annot_kws={"size":7})
        ax.set_title(e, fontsize=10, color=ECOTYPE_COLORS[e])
    fig.suptitle("7-checkpoint Spearman correlation network", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4d_checkpoint_network.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4d_checkpoint_network.pdf", bbox_inches="tight")
    plt.close(fig)

    # Fig 3: CD276 vs TIS scatter (highlights orthogonality)
    print("[*] Figure 3: CD276 vs TIS scatter")
    df_s = pd.DataFrame({"CD276": cp_expr["CD276"], "TIS": TIS.reindex(cp_expr.index),
                         "ecotype": cp_expr["ecotype"]}).dropna()
    rho_all = stats.spearmanr(df_s.CD276, df_s.TIS)
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    for e in ECOTYPE_ORDER:
        sub = df_s[df_s.ecotype==e]
        ax.scatter(sub.TIS, sub.CD276, c=ECOTYPE_COLORS[e], s=14, alpha=0.7,
                   label=f"{e} (n={len(sub)})", edgecolor="none")
    ax.axhline(df_s.CD276.quantile(0.75), color="grey", ls="--", lw=0.6)
    ax.axvline(df_s.TIS.quantile(0.25),   color="grey", ls="--", lw=0.6)
    ax.set_xlabel("TIS (Danaher 18-gene)"); ax.set_ylabel("CD276 (B7-H3) log2 TPM")
    ax.set_title(f"CD276 vs TIS — Spearman ρ={rho_all.correlation:+.2f} (p={rho_all.pvalue:.1e})\n"
                 f"CD276-hi ∩ TIS-lo (CAR-T candidate) n={len(cd276_hi_lymph_lo)}",
                 fontsize=9)
    ax.legend(loc="lower right", fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4d_CD276_vs_TIS.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4d_CD276_vs_TIS.pdf", bbox_inches="tight")
    plt.close(fig)

    # Fig 4: Bagaev mean NES heatmap (ecotype × module)
    print("[*] Figure 4: Bagaev mean-NES heatmap (ecotype × cohort)")
    eco_mean = bag_with_meta[bag_with_meta.ecotype.isin(ECOTYPE_ORDER)] \
                .groupby("ecotype")[modules].mean().T.loc[modules, ECOTYPE_ORDER].astype(float)
    coh_mean = bag_with_meta[bag_with_meta.cohort_group.isin(COHORT_ORDER)] \
                .groupby("cohort_group")[modules].mean().T.loc[modules, COHORT_ORDER].astype(float)
    fig, axes = plt.subplots(1,2,figsize=(9.0, 0.30*len(modules)+1.5))
    sns.heatmap(eco_mean, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                cbar_kws={"label":"mean NES"}, ax=axes[0],
                linewidths=0.3, linecolor="white", annot_kws={"size":7})
    axes[0].set_title("Mean Bagaev NES — ecotype", fontsize=10)
    axes[0].set_xlabel(""); axes[0].set_ylabel("")
    sns.heatmap(coh_mean, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                cbar_kws={"label":"mean NES"}, ax=axes[1],
                linewidths=0.3, linecolor="white", annot_kws={"size":7})
    axes[1].set_title("Mean Bagaev NES — cohort", fontsize=10)
    axes[1].set_xlabel(""); axes[1].set_ylabel("")
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4d_bagaev_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4d_bagaev_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # Fig 5: Anti vs Pro Bagaev coupling (Levine Fig 2b re-creation)
    print("[*] Figure 5: Anti vs Pro Bagaev module coupling")
    ANTI = [m for m in modules if "Anti" in m or "MHC" in m or "_T_cells" in m or "_NK_cells" in m
            or "_B_cells" in m or "M1_signature" in m or "Th1" in m or "Effector" in m
            or "Coactivation" in m or "Checkpoint_inhibitor" in m]
    PRO  = [m for m in modules if "Treg" in m or "TAMs" in m or "Suppression" in m or "Neutrophil" in m
            or "Granulocyte" in m or "Th2" in m or "Macrophage_DC" in m or "Myeloid_cells_traffic" in m
            or "Protumor" in m]
    bag_with_meta["anti_score"] = bag_with_meta[ANTI].mean(axis=1)
    bag_with_meta["pro_score"]  = bag_with_meta[PRO].mean(axis=1)
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    for e in ECOTYPE_ORDER:
        sub = bag_with_meta[bag_with_meta.ecotype==e]
        ax.scatter(sub.anti_score, sub.pro_score, c=ECOTYPE_COLORS[e], s=14, alpha=0.7,
                   label=f"{e} (n={len(sub)})", edgecolor="none")
    rho = stats.spearmanr(bag_with_meta.anti_score, bag_with_meta.pro_score)
    ax.set_xlabel("Anti-tumor composite (mean NES)")
    ax.set_ylabel("Pro-tumor composite (mean NES)")
    ax.set_title(f"Bagaev anti vs pro coupling — ρ={rho.correlation:+.2f} (Levine Fig 2b style)", fontsize=9)
    ax.legend(loc="lower right", fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4d_anti_pro_coupling.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4d_anti_pro_coupling.pdf", bbox_inches="tight")
    plt.close(fig)

    # ===================== summary JSON =====================
    summary = {
        "n_samples": int(len(samples)),
        "TIS_genes_present": tis_present,
        "TIS_ecotype_KW": {"H":float(H_e), "p":float(p_e),
                           "eps2": float(eps2(H_e,3,sum(map(len,tis_eco)))),
                           "medians": {e: float(np.median(v)) for e,v in zip(ECOTYPE_ORDER,tis_eco)}},
        "TIS_cohort_KW":  {"H":float(H_c), "p":float(p_c),
                           "eps2": float(eps2(H_c,4,sum(map(len,tis_coh)))),
                           "medians": {c: float(np.median(v)) for c,v in zip(COHORT_ORDER,tis_coh)}},
        "CD276_ecotype_KW": {"H":float(H_cd), "p":float(p_cd),
                              "medians": {e: float(np.median(v)) for e,v in zip(ECOTYPE_ORDER,cd276_eco)}},
        "CD276hi_TISlo_n": int(len(cd276_hi_lymph_lo)),
        "overall_checkpoint_corr_mean": float(overall_corr.values[np.triu_indices_from(overall_corr,1)].mean()),
        "CD276_vs_others_meanrho_overall": float(overall_corr.loc["CD276", [c for c in CHECKPOINTS if c!="CD276"]].mean()),
        "top_bagaev_ecotype_eps2": kw.head(8)[["module","eps2","q_BH"]].to_dict(orient="records"),
        "top_bagaev_cohort_eps2": ckw.head(8)[["module","eps2","q_BH"]].to_dict(orient="records"),
        "anti_pro_coupling_rho": float(rho.correlation),
    }
    with open(OUT/"step4d_summary.json","w") as fh:
        json.dump(summary, fh, indent=2)
    print("\n[+] DONE — outputs in", OUT, "and figs in", FIGDIR)

if __name__ == "__main__":
    main()
