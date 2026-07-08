"""
Step 4f — Optimal-cutoff bootstrap PFS analysis (Levine 2024 surv_cutpoint analog)
==================================================================================
Goal
----
Test whether the failure of the median-split TIS-PFS analysis in Step 4e
(p=0.16) is a *true* null, or simply reflects suboptimal cut-off choice.

Levine 2024 Methods reports:
    "The optimal cutoff for dichotomizing groups as high/low TIS were determined
     using the surv_cutoff function in the R survminer package and p-values
     were calculated through bootstrapping with 1000 iterations, each sampling
     90% of the dataset."

We replicate this in Python with three p-values for each candidate biomarker:
    1. naïve_optimal_p   — minimum log-rank p across all admissible cut-offs
                            (overoptimistic; equivalent to Levine's surv_cutpoint
                             raw p before bootstrap correction).
    2. bootstrap_p       — Levine method: 1000 × 90% subsamples; report median
                            optimal-cutoff p across resamples (Levine "reported p").
    3. permutation_p     — proper maximally-selected log-rank correction:
                            shuffle event labels, compute max log-rank; corrected
                            p = fraction of permutations with max stat ≥ observed.

Biomarkers tested (BRAF V600E LGG cohort, n=52 PFS-analyzable, 22 events):
    • TIS (Danaher 18-gene)            — main Levine claim
    • Bagaev_Effector_cells            — alternative immune-cytolytic axis
    • Bagaev_Checkpoint_inhibitor      — direct ICI-target axis
    • Bagaev_TAMs                       — myeloid-skewed prognostic candidate

Inputs
------
- output/step4e_v600e_cluster_assignment.tsv   (cohort + TIS + EFS)
- output/step4d_bagaev_ssGSEA.tsv              (NES for 24 modules)
  (covers ecotype-main cohort but not V600E — we recompute V600E)
- output/tpm_for_cibersortx.tsv                (for V600E ssGSEA fallback)

Outputs (output/step4f_*)
-------
- step4f_cutoff_search.tsv
- step4f_bootstrap_pvalues.tsv
- step4f_permutation_null.tsv
- step4f_summary.json
- figs_step4f/*.png|pdf
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import gseapy as gp
from lifelines.statistics import logrank_test as _lr_lifelines
from lifelines import KaplanMeierFitter
from scipy.stats import chi2

# -------- fast vectorised log-rank (two-sample) --------
def fast_logrank(t1, e1, t2, e2):
    """Return (test_statistic, p_value) for the standard 2-sample log-rank test
    (chi-square, 1 df). Vectorised over event times for speed."""
    times = np.concatenate([t1, t2])
    events = np.concatenate([e1, e2])
    groups = np.concatenate([np.zeros_like(t1, dtype=int), np.ones_like(t2, dtype=int)])
    order = np.argsort(times, kind="mergesort")
    times = times[order]; events = events[order]; groups = groups[order]
    # unique event times
    uniq, inv = np.unique(times, return_inverse=True)
    n_t = len(uniq)
    # at-risk counts at each unique time per group (running back from totals)
    n_at_risk = np.zeros((n_t, 2))
    d         = np.zeros((n_t, 2))   # observed events per group per time
    n1 = len(t1); n2 = len(t2)
    # totals at risk before each event time
    cum_total = np.bincount(inv, minlength=n_t)
    n_at_risk_total = n1 + n2 - np.concatenate([[0], np.cumsum(cum_total)[:-1]])
    # group-specific at-risk and event counts
    grp_mask = (groups == 1)
    cum_g1 = np.bincount(inv[grp_mask], minlength=n_t)
    n_at_risk[:,1] = n2 - np.concatenate([[0], np.cumsum(cum_g1)[:-1]])
    n_at_risk[:,0] = n_at_risk_total - n_at_risk[:,1]
    # events
    d[:,1] = np.bincount(inv[grp_mask & (events==1)], minlength=n_t)
    d[:,0] = np.bincount(inv[(~grp_mask) & (events==1)], minlength=n_t) if ((~grp_mask)&(events==1)).any() else 0
    d_total = d.sum(axis=1)
    n_total = n_at_risk_total.astype(float)
    keep = (n_total > 0) & (d_total > 0)
    if not keep.any(): return 0.0, 1.0
    e1_exp = d_total[keep] * n_at_risk[keep,0] / n_total[keep]
    var = (n_at_risk[keep,0] * n_at_risk[keep,1] * d_total[keep] * (n_total[keep]-d_total[keep])) \
          / (n_total[keep]**2 * np.maximum(n_total[keep]-1, 1))
    O1 = d[keep,0].sum(); E1 = e1_exp.sum(); V = var.sum()
    if V <= 0: return 0.0, 1.0
    stat = (O1 - E1)**2 / V
    p = float(chi2.sf(stat, df=1))
    return float(stat), p
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT  = Path("/sessions/blissful-gifted-dirac/mnt/Open PBTA")
OUT   = ROOT / "output"
FIGDIR = OUT / "figs_step4f"; FIGDIR.mkdir(exist_ok=True, parents=True)

CLUSTER_FN = OUT / "step4e_v600e_cluster_assignment.tsv"
TPM_FN     = OUT / "tpm_for_cibersortx.tsv"
GMT_FN     = ROOT / "data" / "levine2024_bagaev_modules.gmt"

EVENT_TYPES = {"Progressive","Progressive - Metastatic","Recurrence","Recurrence - Metastatic"}
TIS_GENES = ["CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
             "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
             "PSMB10","STAT1","TIGIT"]
RNG = np.random.default_rng(20260617)

# ---------------------------------------------------------------
def load_gmt(fn):
    d={}
    for line in open(fn):
        n,_,*g = line.rstrip("\n").split("\t"); d[n]=[x for x in g if x]
    return d

def admissible_cutoffs(scores, min_prop=0.25):
    """Return sorted unique values of `scores` whose use as a cutoff leaves
    at least `min_prop` of samples on each side. Returns midpoints between
    consecutive sorted values to avoid duplicates exactly on cutoff."""
    s = np.sort(np.unique(scores))
    n = len(scores)
    cuts = []
    for c in s:
        n_lo = np.sum(scores <= c)
        if n_lo >= max(2, int(np.ceil(n*min_prop))) and (n - n_lo) >= max(2, int(np.ceil(n*min_prop))):
            cuts.append(c)
    return np.array(cuts)

def best_cutoff(times, events, scores, min_prop=0.25):
    """Sweep cutoffs and return (best_cutoff, best_p, best_stat)."""
    cuts = admissible_cutoffs(scores, min_prop=min_prop)
    best = (np.nan, 1.0, 0.0)
    for c in cuts:
        hi = scores > c
        if hi.sum() < 2 or (~hi).sum() < 2: continue
        stat, p = fast_logrank(times[hi], events[hi], times[~hi], events[~hi])
        if p < best[1]:
            best = (float(c), float(p), float(stat))
    return best

def bootstrap_optimal_p(times, events, scores, n_iter=1000, subsample_frac=0.9, min_prop=0.25):
    """Levine 2024 method: for each of n_iter resamples (subsample without replacement,
    90% of original), find optimal cutoff and corresponding log-rank p."""
    n = len(times)
    n_sub = int(np.round(subsample_frac * n))
    ps = np.empty(n_iter)
    cs = np.empty(n_iter)
    for i in range(n_iter):
        idx = RNG.choice(n, size=n_sub, replace=False)
        t, e, s = times[idx], events[idx], scores[idx]
        c, p, _ = best_cutoff(t, e, s, min_prop=min_prop)
        ps[i] = p; cs[i] = c
    return ps, cs

def permutation_pvalue(times, events, scores, n_perm=1000, min_prop=0.25):
    """Maximally-selected log-rank correction. Permute event labels (Bres–Halmos style)
    and compute max log-rank statistic; corrected p = fraction of permutations with
    max stat >= observed."""
    _, obs_p, obs_stat = best_cutoff(times, events, scores, min_prop=min_prop)
    null_stats = np.empty(n_perm)
    for i in range(n_perm):
        perm_idx = RNG.permutation(len(events))
        _, _, stat = best_cutoff(times, events[perm_idx], scores, min_prop=min_prop)
        null_stats[i] = stat
    pcorr = (1 + np.sum(null_stats >= obs_stat)) / (1 + n_perm)
    return pcorr, obs_stat, obs_p, null_stats

# ---------------------------------------------------------------
def main():
    print("[*] Loading V600E cluster assignment + clinical ...", flush=True)
    cl = pd.read_csv(CLUSTER_FN, sep="\t")
    print(f"    n total V600E in assignment = {len(cl)}")
    cl["event"] = cl["EFS_event_type"].isin(EVENT_TYPES).astype(int)
    cl["EFS_days"] = pd.to_numeric(cl["EFS_days"], errors="coerce")
    cl = cl.dropna(subset=["EFS_days","TIS"]).copy()
    cl = cl[cl["EFS_days"] >= 0]
    print(f"    PFS-analyzable V600E n = {len(cl)} (events={int(cl.event.sum())})")

    # Compute Bagaev modules on V600E samples (cached to avoid recomputation)
    CACHE_FN = OUT/"step4f_v600e_bagaev.tsv"
    if CACHE_FN.exists():
        print(f"\n[*] Loading cached Bagaev V600E ssGSEA from {CACHE_FN}", flush=True)
        bag = pd.read_csv(CACHE_FN, sep="\t", index_col=0)
    else:
        print("\n[*] Recomputing Bagaev ssGSEA on V600E samples ...", flush=True)
        tpm = pd.read_csv(TPM_FN, sep="\t")
        tpm = tpm.rename(columns={tpm.columns[0]:"GeneSymbol"}).drop_duplicates("GeneSymbol").set_index("GeneSymbol")
        keep = [s for s in cl["Kids_First_Biospecimen_ID"] if s in tpm.columns]
        log_tpm = np.log2(tpm[keep].astype(float)+1.0)
        sets = load_gmt(GMT_FN)
        res = gp.ssgsea(data=log_tpm, gene_sets=sets, sample_norm_method="rank",
                        no_plot=True, threads=1, min_size=3, max_size=500,
                        permutation_num=0, outdir=None)
        bag = res.res2d.copy(); bag["NES"]=pd.to_numeric(bag["NES"], errors="coerce")
        bag = bag.pivot(index="Name", columns="Term", values="NES").astype(float)
        bag.index.name = "Kids_First_Biospecimen_ID"
        bag.to_csv(CACHE_FN, sep="\t")
    # Merge
    surv = cl.merge(bag.reset_index(), on="Kids_First_Biospecimen_ID", how="inner")
    print(f"    surv table shape = {surv.shape}")

    BIOMARKERS = ["TIS", "Bagaev_Effector_cells", "Bagaev_Checkpoint_inhibitor", "Bagaev_TAMs"]
    times  = surv["EFS_days"].values.astype(float)
    events = surv["event"].values.astype(int)

    # ---------- Cutoff search per biomarker ----------
    print("\n[*] Best cutoff search per biomarker (BRAF V600E LGG)", flush=True)
    rows = []
    cutoff_curves = {}
    for b in BIOMARKERS:
        scores = surv[b].values.astype(float)
        cuts = admissible_cutoffs(scores, min_prop=0.25)
        pps = []
        for c in cuts:
            hi = scores > c
            stat, p = fast_logrank(times[hi], events[hi], times[~hi], events[~hi])
            pps.append(p)
        cutoff_curves[b] = (cuts, np.array(pps))
        c_best, p_best, stat_best = best_cutoff(times, events, scores, min_prop=0.25)
        # median split for comparison
        med = float(np.median(scores))
        hi_m = scores > med
        med_stat, med_p = fast_logrank(times[hi_m], events[hi_m], times[~hi_m], events[~hi_m])
        rows.append({"biomarker":b, "n":len(scores), "n_events":int(events.sum()),
                     "median_split_cutoff":med, "median_split_p": float(med_p),
                     "best_cutoff":c_best, "naive_best_p":p_best, "best_stat":stat_best,
                     "n_admissible_cutoffs":len(cuts)})
    search = pd.DataFrame(rows)
    search.to_csv(OUT/"step4f_cutoff_search.tsv", sep="\t", index=False)
    print(search.to_string(index=False))

    # ---------- Bootstrap (Levine method) + Permutation null ----------
    print("\n[*] Bootstrap (1000 × 90 %) + Permutation correction (1000 perms)", flush=True)
    boot_rows = []
    null_records = {}
    boot_dists = {}   # cache for Fig 2
    for b in BIOMARKERS:
        scores = surv[b].values.astype(float)
        # Levine bootstrap p
        boot_ps, boot_cs = bootstrap_optimal_p(times, events, scores, n_iter=1000,
                                                subsample_frac=0.9, min_prop=0.25)
        boot_dists[b] = boot_ps
        # Permutation correction
        pcorr, obs_stat, obs_p, null_stats = permutation_pvalue(times, events, scores,
                                                                  n_perm=1000, min_prop=0.25)
        null_records[b] = null_stats
        boot_rows.append({"biomarker":b,
                          "naive_best_p":float(obs_p),
                          "best_stat":float(obs_stat),
                          "bootstrap_median_p":float(np.median(boot_ps)),
                          "bootstrap_25th_p": float(np.percentile(boot_ps,25)),
                          "bootstrap_75th_p": float(np.percentile(boot_ps,75)),
                          "fraction_boot_p_lt_0.05": float((boot_ps<0.05).mean()),
                          "permutation_corrected_p":float(pcorr)})
    boots = pd.DataFrame(boot_rows)
    boots.to_csv(OUT/"step4f_bootstrap_pvalues.tsv", sep="\t", index=False)
    print(boots.to_string(index=False))
    pd.DataFrame(null_records).to_csv(OUT/"step4f_permutation_null.tsv", sep="\t", index=False)

    # ---------- Figures (300 dpi) ----------
    plt.rcParams.update({"figure.dpi":300,"savefig.dpi":300,"font.size":9,
                         "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})

    # Figure 1: cutoff search curves (4 biomarkers)
    print("[*] Figure 1: cutoff search curves", flush=True)
    fig, axes = plt.subplots(2,2,figsize=(8.5,6.0))
    for ax, b in zip(axes.flat, BIOMARKERS):
        cuts, pps = cutoff_curves[b]
        ax.plot(cuts, -np.log10(pps), color="#2C7BB6", lw=1.0)
        best_idx = np.argmin(pps)
        ax.axvline(cuts[best_idx], color="#D7301F", ls="--", lw=0.8, label=f"best cut = {cuts[best_idx]:.2f}")
        ax.axhline(-np.log10(0.05), color="grey", ls=":", lw=0.7)
        ax.set_xlabel(f"{b} cutoff"); ax.set_ylabel("-log10(log-rank p)")
        ax.set_title(f"{b}  best p={pps[best_idx]:.3f}", fontsize=9)
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4f_cutoff_search.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4f_cutoff_search.pdf", bbox_inches="tight"); plt.close(fig)

    # Figure 2: bootstrap p-value distributions
    print("[*] Figure 2: bootstrap p distributions", flush=True)
    fig, axes = plt.subplots(2,2,figsize=(8.5,6.0))
    for ax, b in zip(axes.flat, BIOMARKERS):
        boot_ps = boot_dists[b]
        ax.hist(boot_ps, bins=40, color="#2C7BB6", alpha=0.85, edgecolor="white")
        med = np.median(boot_ps); naive = boots[boots.biomarker==b].naive_best_p.iloc[0]
        ax.axvline(med, color="#D7301F", ls="--", lw=0.9, label=f"median p={med:.3f}")
        ax.axvline(0.05, color="grey", ls=":", lw=0.7, label="α=0.05")
        ax.set_xlabel("Optimal-cutoff log-rank p (each bootstrap)")
        ax.set_ylabel("Count")
        ax.set_title(f"{b}  naïve best p={naive:.3f}", fontsize=9)
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4f_bootstrap_distribution.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4f_bootstrap_distribution.pdf", bbox_inches="tight"); plt.close(fig)

    # Figure 3: KM curves at optimal cutoff for each biomarker
    print("[*] Figure 3: KM at optimal cutoff", flush=True)
    fig, axes = plt.subplots(2,2,figsize=(8.5,6.4))
    kmf = KaplanMeierFitter()
    for ax,b in zip(axes.flat, BIOMARKERS):
        scores = surv[b].values.astype(float)
        c_best = search[search.biomarker==b].best_cutoff.iloc[0]
        hi = scores > c_best
        for grp, mask, col in [(f"high (>{c_best:.2f})", hi, "#D7301F"),
                               (f"low (≤{c_best:.2f})",  ~hi, "#92C5DE")]:
            t_g = times[mask]/365.25; e_g = events[mask]
            if len(t_g):
                kmf.fit(t_g, e_g, label=f"{grp}  n={mask.sum()} ev={int(e_g.sum())}")
                kmf.plot_survival_function(ax=ax, color=col, ci_show=False, linewidth=1.4)
        bp = boots[boots.biomarker==b].iloc[0]
        ax.set_title(f"{b} optimal-cutoff KM\nnaive p={bp.naive_best_p:.3f}  boot median p={bp.bootstrap_median_p:.3f}  perm-corr p={bp.permutation_corrected_p:.3f}",
                     fontsize=8)
        ax.set_xlabel("Years"); ax.set_ylabel("PFS probability"); ax.set_ylim(0,1.02)
        ax.legend(fontsize=7, loc="lower left", frameon=False)
    fig.tight_layout()
    fig.savefig(FIGDIR/"step4f_KM_optimal_cutoff.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR/"step4f_KM_optimal_cutoff.pdf", bbox_inches="tight"); plt.close(fig)

    # ---------- summary JSON ----------
    summary = {
        "cohort": "BRAF V600E LGG (OpenPedCan)",
        "n_analyzable": int(len(surv)),
        "n_events": int(events.sum()),
        "biomarkers_tested": BIOMARKERS,
        "method_note":
            ("naive_best_p = minimum log-rank p over admissible cutoffs (>=25% per side); "
             "bootstrap_median_p = Levine 2024 method (1000x90% subsample, median optimal-cut p); "
             "permutation_corrected_p = max log-rank vs 1000 event-shuffles (proper multiplicity correction)."),
        "rows": (search.merge(boots, on=["biomarker"], how="outer")
                 .to_dict(orient="records")),
        "compare_to_step4e": {
            "median_split_TIS_p_step4e": 0.1614,
            "median_split_TIS_p_now":     float(search.loc[search.biomarker=="TIS","median_split_p"].iloc[0]),
            "Levine_reported_TIS_p":      0.047,
        },
    }
    with open(OUT/"step4f_summary.json","w") as fh:
        json.dump(summary, fh, indent=2)
    print("\n[+] DONE. summary saved to step4f_summary.json")
    print(json.dumps(summary, indent=2, default=float))

if __name__ == "__main__":
    main()
.json")
    print(json.dumps(summary, indent=2, default=float))

if __name__ == "__main__":
    main()
