"""
Step 4e — BRAF V600E LGG 3-cluster replication (Levine 2024 Fig 4 in our cohort)
==============================================================================
Levine 2024 NanoString clustering of BRAF V600E LGG (n=38) yielded 3 immune
clusters with prognostic significance (high TIS = worse PFS). Here we attempt
to replicate this finding using OpenPedCan BRAF_ALT_LGG cohort
(V600E n=68, ~1.8× Levine).

Steps
-----
1. Subset cohort to BRAF V600E vs KIAA1549-BRAF fusion vs other-MAPK/RTK.
2. Compute Bagaev 24 + TIS ssGSEA scores on V600E samples.
3. Hierarchical clustering (Ward.D2 + Euclidean) on scaled (z) feature matrix.
4. Determine k=3 cluster assignment (Levine's pre-specified k).
5. Compare cluster-level TIS, Bagaev anti/pro composite, ecotype-like signature.
6. PFS Kaplan–Meier (EFS_days + Progressive/Recurrence) by cluster:
     (a) 3-cluster log-rank
     (b) Cluster 1 vs (2+3) per Levine
7. Compare V600E vs KIAA-fusion median TIS (Levine: V600E slightly higher).

Inputs
------
- output/tpm_for_cibersortx.tsv
- output/cohort_BRAF_ALT_LGG.tsv
- data/levine2024_bagaev_modules.gmt
- output/Step4d Bagaev results (re-computed here for LGG samples)

Outputs (output/step4e_*)
------
- step4e_v600e_cluster_assignment.tsv
- step4e_v600e_ssGSEA_scores.tsv
- step4e_cluster_TIS_summary.tsv
- step4e_KM_v600e_logrank.tsv
- step4e_summary.json
- figs_step4e/*.png|pdf (300 dpi)
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import gseapy as gp
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT     = Path("/sessions/blissful-gifted-dirac/mnt/Open PBTA")
OUT      = ROOT / "output"
FIGDIR   = OUT / "figs_step4e"; FIGDIR.mkdir(exist_ok=True, parents=True)

TPM_FN  = OUT / "tpm_for_cibersortx.tsv"
LGG_FN  = OUT / "cohort_BRAF_ALT_LGG.tsv"
GMT_FN  = ROOT / "data" / "levine2024_bagaev_modules.gmt"

TIS_GENES = ["CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
             "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
             "PSMB10","STAT1","TIGIT"]

CLUSTER_PALETTE = {"C1_hot": "#D7301F", "C2_intermediate": "#F4A582", "C3_cold": "#92C5DE"}
CLUSTER_ORDER   = ["C1_hot","C2_intermediate","C3_cold"]

def load_gmt(fn):
    out = {}
    for line in open(fn):
        n, _, *g = line.rstrip("\n").split("\t")
        out[n] = [x for x in g if x]
    return out

def classify_braf(s):
    if not isinstance(s, str): return "other"
    s_low = s.lower()
    if "v600e" in s_low: return "V600E"
    if "kiaa1549-braf" in s_low: return "KIAA_fusion"
    if "clcn6" in s_low or "tax1bp1" in s_low: return "other_BRAF_fusion"
    return "other_MAPK_or_RTK"

def main():
    print("[*] Loading cohort + TPM ...", flush=True)
    cohort = pd.read_csv(LGG_FN, sep="\t")
    cohort["braf_class"] = cohort["molecular_subtype"].apply(classify_braf)
    print(cohort.braf_class.value_counts().to_string())
    tpm = pd.read_csv(TPM_FN, sep="\t")
    tpm = tpm.rename(columns={tpm.columns[0]:"GeneSymbol"}).drop_duplicates("GeneSymbol").set_index("GeneSymbol")

    v600e_ids = cohort.loc[cohort.braf_class=="V600E", "Kids_First_Biospecimen_ID"].tolist()
    kiaa_ids  = cohort.loc[cohort.braf_class=="KIAA_fusion", "Kids_First_Biospecimen_ID"].tolist()
    print(f"\n[*] V600E n={len(v600e_ids)} | KIAA-fusion n={len(kiaa_ids)}")
    print(f"    (Levine 2024 had V600E n=38, KIAA n=79 — replication factor ≈{len(v600e_ids)/38:.1f}×)")

    # log2(TPM+1) for V600E and KIAA subsets
    samples_all = [b for b in v600e_ids + kiaa_ids if b in tpm.columns]
    log_tpm = np.log2(tpm[samples_all].astype(float) + 1.0)

    # ===================== ssGSEA on Bagaev 24 + TIS =====================
    print("\n[*] Computing ssGSEA on Bagaev+TIS modules for V600E + KIAA subset ...", flush=True)
    sets = load_gmt(GMT_FN)
    res = gp.ssgsea(data=log_tpm, gene_sets=sets, sample_norm_method="rank",
                    no_plot=True, threads=1, min_size=3, max_size=500,
                    permutation_num=0, outdir=None)
    bag = res.res2d.copy(); bag["NES"] = pd.to_numeric(bag["NES"], errors="coerce")
    bag = bag.pivot(index="Name", columns="Term", values="NES").astype(float)
    bag.index.name = "Kids_First_Biospecimen_ID"

    # Per-sample TIS (mean log2 TPM of 18 genes)
    tis_present = [g for g in TIS_GENES if g in log_tpm.index]
    TIS = log_tpm.loc[tis_present].mean(axis=0)
    TIS.name = "TIS_mean_log2TPM"

    # ===================== Hierarchical clustering of V600E (Ward.D2) =====================
    print("\n[*] Hierarchical clustering of V600E (Ward.D2) ...", flush=True)
    v_in = [s for s in v600e_ids if s in bag.index]
    Xv = bag.loc[v_in].copy()
    Xv_z = pd.DataFrame(StandardScaler().fit_transform(Xv.values),
                        index=Xv.index, columns=Xv.columns)
    Z = linkage(Xv_z.values, method="ward", metric="euclidean")
    cluster_k3 = fcluster(Z, t=3, criterion="maxclust")
    cl_df = pd.DataFrame({"Kids_First_Biospecimen_ID": v_in, "cluster_raw": cluster_k3})
    # rename clusters by median TIS (C1=hot, C2=intermediate, C3=cold)
    cl_df["TIS"] = cl_df["Kids_First_Biospecimen_ID"].map(TIS)
    cl_means = cl_df.groupby("cluster_raw")["TIS"].median().sort_values(ascending=False)
    rename = {cl_means.index[0]:"C1_hot", cl_means.index[1]:"C2_intermediate", cl_means.index[2]:"C3_cold"}
    cl_df["cluster"] = cl_df["cluster_raw"].map(rename)
    cohort_v = cohort.set_index("Kids_First_Biospecimen_ID").loc[v_in].reset_index()
    cl_df = cl_df.merge(cohort_v[["Kids_First_Biospecimen_ID","EFS_days","EFS_event_type",
                                  "molecular_subtype","age_years","reported_gender",
                                  "OS_days","OS_status","tumor_descriptor","CNS_region"]],
                        on="Kids_First_Biospecimen_ID", how="left")
    cl_df.to_csv(OUT/"step4e_v600e_cluster_assignment.tsv", sep="\t", index=False)
    print(cl_df.cluster.value_counts().to_string())

    # KW per Bagaev module across 3 clusters
    print("\n[*] KW per module across 3 V600E clusters ...", flush=True)
    bag_v = bag.loc[v_in].join(cl_df.set_index("Kids_First_Biospecimen_ID")["cluster"])
    kw_rows = []
    for m in bag.columns:
        g = [bag_v.loc[bag_v.cluster==c, m].dropna().values for c in CLUSTER_ORDER]
        H,p = stats.kruskal(*g)
        kw_rows.append({"module":m,"H":H,"p":p,
                        **{f"median_{c}":float(np.median(v)) for c,v in zip(CLUSTER_ORDER,g)}})
    kw = pd.DataFrame(kw_rows).sort_values("H", ascending=False).reset_index(drop=True)
    kw.to_csv(OUT/"step4e_KW_by_v600e_cluster.tsv", sep="\t", index=False)

    bag_v["TIS"] = TIS.reindex(bag_v.index)
    tis_summary = bag_v.groupby("cluster")["TIS"].describe()[["count","mean","50%","std"]]
    tis_summary.to_csv(OUT/"step4e_cluster_TIS_summary.tsv", sep="\t")
    print("\nTIS per cluster:"); print(tis_summary.round(3))

    # TIS V600E vs KIAA
    kiaa_in = [s for s in kiaa_ids if s in bag.index]
    tis_v = TIS.loc[v_in].dropna()
    tis_k = TIS.loc[kiaa_in].dropna()
    mw = stats.mannwhitneyu(tis_v, tis_k, alternative="two-sided")
    print(f"\nTIS V600E (n={len(tis_v)}) median={tis_v.median():.3f} vs KIAA (n={len(tis_k)}) median={tis_k.median():.3f}  Mann-Whitney p={mw.pvalue:.3e}")

    # ===================== Kaplan–Meier (PFS / EFS) =====================
    print("\n[*] Survival analysis (PFS) ...", flush=True)
    surv = cl_df.copy()
    EVENT_TYPES = {"Progressive","Progressive - Metastatic","Recurrence","Recurrence - Metastatic"}
    surv["event"] = surv["EFS_event_type"].isin(EVENT_TYPES).astype(int)
    surv = surv[surv["EFS_days"].notna()].copy()
    surv["EFS_days"] = pd.to_numeric(surv["EFS_days"], errors="coerce")
    surv = surv[surv["EFS_days"].notna()]
    print(f"PFS-analyzable V600E n = {len(surv)} (events = {int(surv.event.sum())})")
    print(surv.groupby("cluster").agg(n=("cluster","size"), n_event=("event","sum"),
                                       median_days=("EFS_days","median")))

    # Overall 3-way logrank
    lr3 = multivariate_logrank_test(surv["EFS_days"], surv["cluster"], surv["event"])
    print(f"3-way log-rank p = {lr3.p_value:.4f}")

    # C1 vs C2+C3
    s1 = surv[surv["cluster"]=="C1_hot"]; s23 = surv[surv["cluster"].isin(["C2_intermediate","C3_cold"])]
    if len(s1) and len(s23):
        lr12 = logrank_test(s1["EFS_days"], s23["EFS_days"], s1["event"], s23["event"])
        print(f"C1 (hot) vs C2+C3 log-rank p = {lr12.p_value:.4f}")
    else:
        lr12 = None

    # Save log-rank table
    pd.DataFrame([{"comparison":"3-cluster", "p_value":lr3.p_value},
                  {"comparison":"C1_hot vs C2+C3", "p_value": getattr(lr12, "p_value", np.nan)}]) \
        .to_csv(OUT/"step4e_KM_v600e_logrank.tsv", sep="\t", index=False)

    # Levine-style TIS-dichotomous PFS (high TIS = worse PFS hypothesis)
    surv["TIS"] = TIS.reindex(surv.set_index("Kids_First_Biospecimen_ID").index).values
    surv_t = surv.dropna(subset=["TIS"]).copy()
    cutoff = surv_t["TIS"].median()
    surv_t["TIS_high"] = surv_t["TIS"] > cutoff
    lr_tis = logrank_test(surv_t.loc[surv_t.TIS_high, "EFS_days"],
                          surv_t.loc[~surv_t.TIS_high, "EFS_days"],
                          surv_t.loc[surv_t.TIS_high, "event"],
                          surv_t.loc[~surv_t.TIS_high, "event"])
    print(f"V600E TIS-high vs TIS-low (median split) log-rank p = {lr_tis.p_value:.4f} (Levine reported p=0.047)")

    # ===================== Figures (300 dpi) =====================
    plt.rcParams.update({"figure.dpi":300,"savefig.dpi":300,"font.size":9,
                         "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})

    # Fig 1: heatmap of Bagaev+TIS NES by V600E sample with cluster colorbar
    print("\n[*] Figure 1: V600E heatmap (clustered) ...", flush=True)
    ordered = cl_df.sort_values(["cluster","TIS"], ascending=[True, False])["Kids_First_Biospecimen_ID"].tolist()
    H = Xv_z.loc[ordered]
    col_colors = [CLUSTER_PALETTE[cl_df.set_index("Kids_First_Biospecimen_ID").loc[s,"cluster"]] for s in ordered]
    fig, ax = plt.subplots(figsize=(11.5, 0.32*len(Xv_z.columns)+1.5))
    sns.heatmap(H.T, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5,
                ax=ax, cbar_kws={"label":"z-score (NES)"},
                xticklabels=False, yticklabels=True,
                linewidths=0)
    # add cluster bar
    for i, c in enumerate(col_colors):
        ax.add_patch(plt.Rectangle((i, -0.7), 1, 0.6, color=c, clip_on=False, ec="none"))
    ax.set_title(f"BRAF V600E LGG hierarchical clustering (n={len(ordered)}) — Bagaev+TIS NES (z)", fontsize=10)
    ax.set_xlabel(""); ax.set_ylabel("")
    legend_elements = [plt.Rectangle((0,0),1,1, fc=CLUSTER_PALETTE[c], label=c) for c in CLUSTER_ORDER]
    ax.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(1.20,1.0), fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4e_v600e_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4e_v600e_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # Fig 2: TIS by V600E cluster + V600E vs KIAA boxplot
    print("[*] Figure 2: TIS by cluster + V600E vs KIAA ...", flush=True)
    fig, axes = plt.subplots(1,2, figsize=(7.0, 3.4))
    df_t = bag_v.reset_index()
    sns.boxplot(data=df_t, x="cluster", y="TIS", order=CLUSTER_ORDER,
                hue="cluster", palette=CLUSTER_PALETTE, legend=False, ax=axes[0],
                fliersize=2, linewidth=0.7)
    axes[0].set_title(f"TIS by V600E cluster (n={len(df_t)})\nLevine Fig 4c style", fontsize=9)
    axes[0].set_xlabel(""); axes[0].set_ylabel("TIS (mean log2 TPM, 18 genes)")
    axes[0].tick_params(axis="x", rotation=15, labelsize=8)

    df_braf = pd.concat([
        pd.DataFrame({"BRAF":"V600E","TIS":tis_v.values}),
        pd.DataFrame({"BRAF":"KIAA1549-BRAF","TIS":tis_k.values})
    ])
    sns.boxplot(data=df_braf, x="BRAF", y="TIS", order=["KIAA1549-BRAF","V600E"],
                palette={"KIAA1549-BRAF":"#92C5DE","V600E":"#D7301F"},
                hue="BRAF", legend=False, ax=axes[1], fliersize=2, linewidth=0.7)
    axes[1].set_title(f"TIS V600E vs KIAA-fusion\nMW p={mw.pvalue:.1e}", fontsize=9)
    axes[1].set_xlabel(""); axes[1].set_ylabel("")
    axes[1].tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4e_TIS_boxplots.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4e_TIS_boxplots.pdf", bbox_inches="tight")
    plt.close(fig)

    # Fig 3: KM curves (3-cluster + TIS-dichotomous)
    print("[*] Figure 3: KM curves ...", flush=True)
    fig, axes = plt.subplots(1,2, figsize=(8.0, 3.6))
    kmf = KaplanMeierFitter()
    for c in CLUSTER_ORDER:
        sub = surv[surv.cluster==c]
        if len(sub):
            kmf.fit(sub["EFS_days"]/365.25, sub["event"], label=f"{c} (n={len(sub)}, ev={int(sub.event.sum())})")
            kmf.plot_survival_function(ax=axes[0], color=CLUSTER_PALETTE[c], ci_show=False, linewidth=1.4)
    axes[0].set_title(f"PFS by V600E cluster\n3-way log-rank p={lr3.p_value:.3f}\nC1 vs C2+C3 p={getattr(lr12,'p_value',np.nan):.3f}", fontsize=8)
    axes[0].set_xlabel("Time (years)"); axes[0].set_ylabel("PFS probability")
    axes[0].set_ylim(0, 1.02); axes[0].legend(fontsize=7, loc="lower left", frameon=False)

    # TIS-dichotomous
    for hi, lbl, col in [(True, f"TIS high (>{cutoff:.2f})", "#D7301F"),
                         (False, f"TIS low (≤{cutoff:.2f})", "#92C5DE")]:
        sub = surv_t[surv_t.TIS_high == hi]
        if len(sub):
            kmf.fit(sub["EFS_days"]/365.25, sub["event"], label=f"{lbl} (n={len(sub)}, ev={int(sub.event.sum())})")
            kmf.plot_survival_function(ax=axes[1], color=col, ci_show=False, linewidth=1.4)
    axes[1].set_title(f"PFS by V600E TIS (median split)\nlog-rank p={lr_tis.p_value:.3f}\n(Levine reported p=0.047)", fontsize=8)
    axes[1].set_xlabel("Time (years)"); axes[1].set_ylabel("PFS probability")
    axes[1].set_ylim(0,1.02); axes[1].legend(fontsize=7, loc="lower left", frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4e_KM_v600e.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4e_KM_v600e.pdf", bbox_inches="tight")
    plt.close(fig)

    # ===================== summary JSON =====================
    summary = {
        "n_V600E": len(v_in),
        "n_KIAA":  len(kiaa_in),
        "levine_replication_factor_V600E": round(len(v_in)/38, 2),
        "cluster_sizes": cl_df.cluster.value_counts().to_dict(),
        "TIS_V600E_median": float(tis_v.median()),
        "TIS_KIAA_median":  float(tis_k.median()),
        "TIS_V600E_vs_KIAA_mw_p": float(mw.pvalue),
        "TIS_per_cluster": tis_summary.to_dict(),
        "logrank_3cluster_p": float(lr3.p_value),
        "logrank_C1_vs_C2plus3_p": float(getattr(lr12, "p_value", np.nan)) if lr12 is not None else None,
        "logrank_TIS_high_vs_low_p_median_split": float(lr_tis.p_value),
        "Levine_reported_TIS_high_vs_low_p": 0.047,
    }
    with open(OUT/"step4e_summary.json","w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    print("\n[+] DONE — outputs in", OUT, "and figs in", FIGDIR)
    print(json.dumps(summary, indent=2, default=float))

if __name__ == "__main__":
    main()
