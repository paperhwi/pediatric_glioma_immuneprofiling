"""
Step 4g — PROGENy proxy + CDKN2A/B crosstab + C1 marker DEG (in single batch)
=============================================================================
B4 (deferred from Levine plan) — PROGENy 14-pathway proxy via HALLMARK-style
                                 gene sets; replicate Levine 2024 Fig 2c.
B-extra — V600E cluster × CDKN2A/B status crosstab (Fisher's exact, OR).
C1-DEG — Wilcoxon rank-sum DEG of C1_hot V600E samples vs C2+C3.

Inputs
------
- output/tpm_for_cibersortx.tsv
- output/ecotype_LM22_main_k3_annotated.tsv  (canonical 349 ecotype)
- output/ecotype_assignment_k3_annotated.tsv (cohort metadata)
- output/step4d_TIS_per_sample.tsv            (TIS scalar)
- output/step4e_v600e_cluster_assignment.tsv  (V600E cluster + CDKN2A/B in molecular_subtype)
- data/progeny_hallmark_proxy.gmt

Outputs (output/step4g_*)
-------
- step4g_progeny_scores.tsv         (n=349 × 14 pathways)
- step4g_progeny_TIS_corr.tsv       (Spearman ρ vs TIS)
- step4g_cdkn_crosstab.tsv          (cluster × CDKN+/- + Fisher's p, OR)
- step4g_C1_DEG.tsv                  (Wilcoxon p, log2FC, BH q for top genes)
- figs_step4g/*.png|pdf
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import gseapy as gp
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path("/sessions/blissful-gifted-dirac/mnt/Open PBTA")
OUT  = ROOT / "output"
FIGDIR = OUT / "figs_step4g"; FIGDIR.mkdir(exist_ok=True, parents=True)

plt.rcParams.update({"figure.dpi":300,"savefig.dpi":300,"font.size":9,
                     "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})

def load_gmt(fn):
    d={}
    for line in open(fn):
        n,_,*g = line.rstrip("\n").split("\t"); d[n]=[x for x in g if x]
    return d

def main():
    print("[*] Loading TPM + ecotype + V600E cluster ...", flush=True)
    tpm = pd.read_csv(OUT/"tpm_for_cibersortx.tsv", sep="\t")
    tpm = tpm.rename(columns={tpm.columns[0]:"GeneSymbol"}).drop_duplicates("GeneSymbol").set_index("GeneSymbol")
    log_tpm = np.log2(tpm.astype(float) + 1.0)

    eco = pd.read_csv(OUT/"ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
    meta = pd.read_csv(OUT/"ecotype_assignment_k3_annotated.tsv").set_index("Kids_First_Biospecimen_ID")
    eco = eco[["ecotype"]].join(meta[["cohort_group","location_class","age_dev_group"]], how="left")
    tis = pd.read_csv(OUT/"step4d_TIS_per_sample.tsv", sep="\t", index_col=0)["TIS_Danaher18"]
    tis.index.name = "Kids_First_Biospecimen_ID"

    ECOTYPE_ORDER = ["Lymphocyte-inflamed","Myeloid-dominant","Immune-desert"]
    ECOTYPE_COLORS = {"Lymphocyte-inflamed":"#2C7BB6","Myeloid-dominant":"#D7301F","Immune-desert":"#7F7F7F"}

    # =============================================================
    # 1. PROGENy proxy ssGSEA on n=349 main cohort
    # =============================================================
    print("\n[1] PROGENy proxy via ssGSEA (HALLMARK-style 14 pathways)", flush=True)
    sets = load_gmt(ROOT/"data"/"progeny_hallmark_proxy.gmt")
    samples_main = [b for b in tpm.columns if b in eco.index]
    log_tpm_main = log_tpm[samples_main]

    res = gp.ssgsea(data=log_tpm_main, gene_sets=sets, sample_norm_method="rank",
                    no_plot=True, threads=1, min_size=3, max_size=500,
                    permutation_num=0, outdir=None)
    pg = res.res2d.copy(); pg["NES"] = pd.to_numeric(pg["NES"], errors="coerce")
    pg = pg.pivot(index="Name", columns="Term", values="NES").astype(float)
    pg.index.name = "Kids_First_Biospecimen_ID"
    pg.to_csv(OUT/"step4g_progeny_scores.tsv", sep="\t")

    # TIS vs PROGENy correlation
    shared = pg.index.intersection(tis.index)
    corrs = {}
    for p in pg.columns:
        r, pv = stats.spearmanr(pg.loc[shared, p], tis.loc[shared])
        corrs[p] = {"rho": float(r), "p": float(pv)}
    corr_df = pd.DataFrame(corrs).T.sort_values("rho", ascending=False)
    corr_df.to_csv(OUT/"step4g_progeny_TIS_corr.tsv", sep="\t")
    print("\n    PROGENy proxy vs TIS Spearman:")
    print(corr_df.round(3).to_string())

    # Figure: PROGENy vs TIS bar (Levine Fig 2c style)
    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    pg_sorted = corr_df.sort_values("rho")
    colors = ["#D7301F" if r<0 else "#2C7BB6" for r in pg_sorted["rho"]]
    ax.barh(range(len(pg_sorted)), pg_sorted["rho"], color=colors, edgecolor="white")
    ax.set_yticks(range(len(pg_sorted)))
    ax.set_yticklabels([p.replace("PROGENy_","") for p in pg_sorted.index], fontsize=8)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_xlabel("Spearman ρ (PROGENy-proxy vs TIS)")
    ax.set_title("PROGENy 14-pathway proxy vs TIS (n=349 main cohort)\n(replicates Levine 2024 Fig 2c)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4g_progeny_TIS_correlation.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4g_progeny_TIS_correlation.pdf", bbox_inches="tight"); plt.close(fig)

    # PROGENy per ecotype heatmap
    pg_eco = pg.join(eco[["ecotype"]], how="inner")
    pg_eco = pg_eco[pg_eco.ecotype.isin(ECOTYPE_ORDER)]
    eco_mean = pg_eco.groupby("ecotype")[pg.columns].mean().T.loc[pg.columns, ECOTYPE_ORDER].astype(float)
    fig, ax = plt.subplots(figsize=(4.4, 0.32*len(pg.columns)+1.5))
    sns.heatmap(eco_mean, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                cbar_kws={"label":"mean NES"}, ax=ax,
                linewidths=0.3, linecolor="white", annot_kws={"size":7})
    ax.set_title("PROGENy-proxy NES per ecotype", fontsize=10)
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_yticklabels([p.replace("PROGENy_","") for p in eco_mean.index], fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4g_progeny_ecotype_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4g_progeny_ecotype_heatmap.pdf", bbox_inches="tight"); plt.close(fig)

    # =============================================================
    # 2. V600E cluster × CDKN2A/B crosstab
    # =============================================================
    print("\n[2] V600E cluster × CDKN2A/B status crosstab", flush=True)
    cl = pd.read_csv(OUT/"step4e_v600e_cluster_assignment.tsv", sep="\t")
    cl["CDKN2AB_lost"] = cl["molecular_subtype"].str.contains("CDKN2A/B", na=False)
    ct = pd.crosstab(cl["cluster"], cl["CDKN2AB_lost"], margins=False)
    ct.columns = ["CDKN_intact","CDKN_lost"]
    print(ct)

    # Fisher's exact: C1 vs C2+C3
    a = ct.loc["C1_hot","CDKN_lost"]; b = ct.loc["C1_hot","CDKN_intact"]
    c = ct.loc[["C2_intermediate","C3_cold"], "CDKN_lost"].sum()
    d = ct.loc[["C2_intermediate","C3_cold"], "CDKN_intact"].sum()
    fisher_or, fisher_p = stats.fisher_exact([[a,b],[c,d]], alternative="greater")
    print(f"\n    Fisher's exact C1 vs C2+C3 CDKN2A/B-lost: OR={fisher_or:.2f}, p={fisher_p:.3f}")

    # Save
    summary_ct = ct.reset_index()
    summary_ct.to_csv(OUT/"step4g_cdkn_crosstab.tsv", sep="\t", index=False)
    with open(OUT/"step4g_cdkn_fisher.json","w") as fh:
        json.dump({"crosstab":ct.to_dict(), "C1_vs_C2C3_OR":float(fisher_or),
                   "C1_vs_C2C3_p":float(fisher_p)}, fh, indent=2)

    # Figure: CDKN crosstab + TIS by CDKN
    fig, axes = plt.subplots(1,2,figsize=(7.0, 3.2))
    ct_pct = ct.div(ct.sum(axis=1), axis=0)*100
    ct_pct.plot(kind="bar", stacked=True, ax=axes[0],
                color=["#92C5DE","#D7301F"], edgecolor="white", width=0.65)
    axes[0].set_title(f"CDKN2A/B-loss frequency by V600E cluster\nFisher C1 vs C2+C3 OR={fisher_or:.2f}, p={fisher_p:.3f}", fontsize=9)
    axes[0].set_xlabel(""); axes[0].set_ylabel("% samples")
    axes[0].legend(loc="upper right", fontsize=8, frameon=False)
    axes[0].tick_params(axis="x", rotation=15, labelsize=8)
    sns.boxplot(data=cl, x="CDKN2AB_lost", y="TIS",
                hue="CDKN2AB_lost", palette={False:"#92C5DE", True:"#D7301F"},
                legend=False, ax=axes[1], fliersize=2, linewidth=0.7)
    axes[1].set_title("TIS by CDKN2A/B status (all V600E)", fontsize=9)
    axes[1].set_xlabel("CDKN2A/B lost"); axes[1].set_ylabel("TIS")
    # MW p
    mw = stats.mannwhitneyu(cl.loc[cl.CDKN2AB_lost,"TIS"], cl.loc[~cl.CDKN2AB_lost,"TIS"], alternative="two-sided")
    axes[1].text(0.97, 0.95, f"MW p={mw.pvalue:.2g}", ha="right", va="top",
                  transform=axes[1].transAxes, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4g_cdkn_analysis.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4g_cdkn_analysis.pdf", bbox_inches="tight"); plt.close(fig)

    # =============================================================
    # 3. C1_hot marker DEG vs C2+C3 (V600E)
    # =============================================================
    print("\n[3] C1_hot marker DEG (V600E)", flush=True)
    v_in = [s for s in cl.Kids_First_Biospecimen_ID if s in log_tpm.columns]
    v_log = log_tpm[v_in]
    cl_v = cl.set_index("Kids_First_Biospecimen_ID").loc[v_in]
    c1 = cl_v.cluster == "C1_hot"
    print(f"    C1_hot n={c1.sum()} vs C2+C3 n={(~c1).sum()}")
    # filter genes with min expression
    expr_filter = (v_log > 1).sum(axis=1) >= 5
    v_log_f = v_log.loc[expr_filter]
    print(f"    genes after filter: {v_log_f.shape[0]}")

    # Wilcoxon rank-sum per gene (vectorised approach)
    print("    running Wilcoxon ranksum on ~{} genes ...".format(v_log_f.shape[0]), flush=True)
    A = v_log_f.loc[:, c1.values].values
    B = v_log_f.loc[:, ~c1.values].values
    pvals = np.empty(v_log_f.shape[0])
    logfc = A.mean(axis=1) - B.mean(axis=1)
    for i in range(v_log_f.shape[0]):
        try:
            pvals[i] = stats.ranksums(A[i], B[i]).pvalue
        except Exception:
            pvals[i] = 1.0
    deg = pd.DataFrame({"gene": v_log_f.index, "log2FC": logfc, "pvalue": pvals})
    deg["q_BH"] = multipletests(deg.pvalue, method="fdr_bh")[1]
    deg = deg.sort_values("pvalue").reset_index(drop=True)
    deg.to_csv(OUT/"step4g_C1_DEG.tsv", sep="\t", index=False)
    print(f"    Top 15 up (log2FC>0) markers of C1_hot:")
    print(deg[(deg.log2FC>1) & (deg.q_BH<0.05)].head(15)[["gene","log2FC","pvalue","q_BH"]].to_string(index=False))
    print(f"    n DEG q<0.05: {(deg.q_BH<0.05).sum()}; up>1 logFC: {((deg.q_BH<0.05) & (deg.log2FC>1)).sum()}; down<-1: {((deg.q_BH<0.05) & (deg.log2FC<-1)).sum()}")

    # Volcano figure
    fig, ax = plt.subplots(figsize=(5.0, 3.8))
    sig_up   = (deg.q_BH<0.05) & (deg.log2FC>1)
    sig_dn   = (deg.q_BH<0.05) & (deg.log2FC<-1)
    ax.scatter(deg.log2FC, -np.log10(deg.pvalue+1e-300), s=4, c="grey", alpha=0.4, edgecolor="none")
    ax.scatter(deg.loc[sig_up,"log2FC"], -np.log10(deg.loc[sig_up,"pvalue"]+1e-300), s=8, c="#D7301F", label=f"up ({sig_up.sum()})", edgecolor="none")
    ax.scatter(deg.loc[sig_dn,"log2FC"], -np.log10(deg.loc[sig_dn,"pvalue"]+1e-300), s=8, c="#2C7BB6", label=f"down ({sig_dn.sum()})", edgecolor="none")
    # label top up
    for _, r in deg[sig_up].nsmallest(15, "pvalue").iterrows():
        ax.text(r.log2FC, -np.log10(r.pvalue+1e-300), r.gene, fontsize=6, ha="left", va="bottom")
    ax.axhline(-np.log10(0.05), color="grey", ls=":", lw=0.6)
    ax.axvline(1, color="grey", ls=":", lw=0.6); ax.axvline(-1, color="grey", ls=":", lw=0.6)
    ax.set_xlabel("log2 fold change (C1_hot vs C2+C3)")
    ax.set_ylabel("-log10(Wilcoxon p)")
    ax.set_title(f"C1_hot V600E marker DEG (n={c1.sum()} vs {(~c1).sum()})", fontsize=9)
    ax.legend(loc="upper left", fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4g_C1_volcano.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4g_C1_volcano.pdf", bbox_inches="tight"); plt.close(fig)

    # ---- summary JSON
    summary = {
        "n_main_cohort": int(len(samples_main)),
        "progeny_pathways": list(pg.columns),
        "progeny_TIS_top3_positive": corr_df.head(3).to_dict(orient="index"),
        "progeny_TIS_top3_negative": corr_df.tail(3).to_dict(orient="index"),
        "cdkn_C1_vs_C2C3_OR": float(fisher_or),
        "cdkn_C1_vs_C2C3_p":  float(fisher_p),
        "cdkn_TIS_MW_p":      float(mw.pvalue),
        "n_C1_hot":           int(c1.sum()),
        "n_C2_C3":            int((~c1).sum()),
        "n_DEG_q05":          int((deg.q_BH<0.05).sum()),
        "n_DEG_up_logFC1":    int(((deg.q_BH<0.05) & (deg.log2FC>1)).sum()),
        "n_DEG_dn_logFC1":    int(((deg.q_BH<0.05) & (deg.log2FC<-1)).sum()),
        "top_C1_up_markers":  deg[(deg.log2FC>1) & (deg.q_BH<0.05)].head(15)["gene"].tolist(),
    }
    with open(OUT/"step4g_summary.json","w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    print("\n[+] DONE — outputs in", OUT)
    print(json.dumps(summary, indent=2, default=float))

if __name__ == "__main__":
    main()
