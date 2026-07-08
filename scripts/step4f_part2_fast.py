"""
Step 4f part 2 — Fast bootstrap + permutation for the optimal-cutoff analysis.
Uses pre-cached output/step4f_v600e_bagaev.tsv.
Vectorised log-rank: precompute sort-by-time once, sweep cutoffs as masks.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2

ROOT  = Path("/sessions/blissful-gifted-dirac/mnt/Open PBTA")
OUT   = ROOT/"output"
EVENT_TYPES = {"Progressive","Progressive - Metastatic","Recurrence","Recurrence - Metastatic"}
BIOMARKERS = ["TIS","Bagaev_Effector_cells","Bagaev_Checkpoint_inhibitor","Bagaev_TAMs"]
N_BOOT = 1000
N_PERM = 1000
RNG = np.random.default_rng(20260617)

# --------- fast log-rank that uses precomputed sort-by-time ---------
def precompute(times, events):
    order = np.argsort(times, kind="mergesort")
    t_s = times[order]; e_s = events[order]
    uniq, inv = np.unique(t_s, return_inverse=True)
    n_t = len(uniq)
    d_total = np.bincount(inv, weights=e_s, minlength=n_t)   # events per time
    cnt     = np.bincount(inv, minlength=n_t)                # samples per time
    n_total = len(times) - np.concatenate([[0], np.cumsum(cnt)[:-1]])  # at risk
    valid   = d_total > 0
    return order, inv, n_total, d_total, valid

def logrank_sweep(order, inv, n_total, d_total, valid, group_g1):
    """Given precomputed survival structure and a boolean (sorted) group-1
    membership vector, compute chi-square log-rank statistic for the
    two-sample test (g1 vs g0)."""
    n_t = len(n_total)
    # group-1 events per time
    d1 = np.bincount(inv, weights=(group_g1).astype(float)*((d_total[inv]>0).astype(float)*np.where(np.zeros_like(group_g1)+1,1,1)), minlength=n_t)
    # The above won't work for events; recompute properly:
    return None  # placeholder

# Simpler: re-implement using full bincount per call (fast enough in numpy)
def fast_logrank_from_sorted(t_s, e_s, g_s):
    """t_s, e_s, g_s are sorted-by-time arrays (g_s is 1 for group-1 member)."""
    n = len(t_s)
    if g_s.sum() < 1 or g_s.sum() >= n: return 0.0, 1.0
    uniq, inv = np.unique(t_s, return_inverse=True)
    n_t = len(uniq)
    n1 = int(g_s.sum())
    cnt_total = np.bincount(inv, minlength=n_t)
    cnt_g1    = np.bincount(inv, weights=g_s, minlength=n_t)
    n_at_risk_total = n - np.concatenate([[0], np.cumsum(cnt_total)[:-1]])
    n_at_risk_g1    = n1 - np.concatenate([[0], np.cumsum(cnt_g1)[:-1]])
    d_total = np.bincount(inv, weights=e_s, minlength=n_t)
    d_g1    = np.bincount(inv, weights=e_s*g_s, minlength=n_t)
    keep = d_total > 0
    n_tot = n_at_risk_total[keep]
    n_1   = n_at_risk_g1[keep]
    d_tot = d_total[keep]
    d_1   = d_g1[keep]
    e1 = d_tot * n_1 / np.maximum(n_tot, 1)
    var = (n_1 * (n_tot-n_1) * d_tot * (n_tot-d_tot)) / (n_tot**2 * np.maximum(n_tot-1, 1))
    O1 = d_1.sum(); E1 = e1.sum(); V = var.sum()
    if V <= 0: return 0.0, 1.0
    stat = (O1 - E1)**2 / V
    return float(stat), float(chi2.sf(stat, df=1))

def admissible_cutoffs(scores, min_prop=0.25):
    s_sorted = np.sort(np.unique(scores))
    n = len(scores)
    min_n = max(2, int(np.ceil(n*min_prop)))
    out = []
    for c in s_sorted:
        n_lo = int(np.sum(scores <= c))
        if n_lo >= min_n and (n - n_lo) >= min_n:
            out.append(c)
    return np.array(out)

def best_cutoff_sorted(t_s, e_s, s_s, min_prop=0.25):
    cuts = admissible_cutoffs(s_s, min_prop=min_prop)
    best = (np.nan, 1.0, 0.0)
    for c in cuts:
        g = (s_s > c).astype(float)
        stat, p = fast_logrank_from_sorted(t_s, e_s, g)
        if p < best[1]:
            best = (float(c), float(p), float(stat))
    return best

def main():
    print("[*] Loading inputs ...", flush=True)
    cl = pd.read_csv(OUT/"step4e_v600e_cluster_assignment.tsv", sep="\t")
    cl["event"] = cl["EFS_event_type"].isin(EVENT_TYPES).astype(int)
    cl["EFS_days"] = pd.to_numeric(cl["EFS_days"], errors="coerce")
    cl = cl.dropna(subset=["EFS_days","TIS"]).copy()
    cl = cl[cl["EFS_days"]>=0]
    bag = pd.read_csv(OUT/"step4f_v600e_bagaev.tsv", sep="\t", index_col=0)
    surv = cl.merge(bag.reset_index(), on="Kids_First_Biospecimen_ID", how="inner")
    print(f"    n V600E PFS-analyzable = {len(surv)} (events={int(surv.event.sum())})")
    print(f"    biomarkers = {BIOMARKERS}")

    rows = []
    boot_dists = {}
    null_records = {}
    for b in BIOMARKERS:
        times  = surv["EFS_days"].values.astype(float)
        events = surv["event"].values.astype(int)
        scores = surv[b].values.astype(float)
        # sort by time once
        order = np.argsort(times, kind="mergesort")
        t_s, e_s, s_s = times[order], events[order], scores[order]
        # observed
        c_best, p_best, stat_best = best_cutoff_sorted(t_s, e_s, s_s, min_prop=0.25)
        med = float(np.median(scores))
        g_med = (s_s > med).astype(float)
        _, med_p = fast_logrank_from_sorted(t_s, e_s, g_med)

        # bootstrap (Levine 2024 method)
        n = len(t_s); n_sub = int(np.round(0.9*n))
        boot_ps = np.empty(N_BOOT)
        for i in range(N_BOOT):
            idx = RNG.choice(n, size=n_sub, replace=False)
            idx.sort()
            tb, eb, sb = t_s[idx], e_s[idx], s_s[idx]
            _, bp, _ = best_cutoff_sorted(tb, eb, sb, min_prop=0.25)
            boot_ps[i] = bp
        boot_dists[b] = boot_ps

        # permutation (event labels permuted)
        null_stats = np.empty(N_PERM)
        for i in range(N_PERM):
            perm = RNG.permutation(n)
            _, _, st = best_cutoff_sorted(t_s, e_s[perm], s_s, min_prop=0.25)
            null_stats[i] = st
        null_records[b] = null_stats
        pcorr = (1 + np.sum(null_stats >= stat_best)) / (1 + N_PERM)

        print(f"    {b:35s} best_cut={c_best:.3f}  naive_p={p_best:.4f}  "
              f"median_split_p={med_p:.4f}  boot_median_p={np.median(boot_ps):.4f}  "
              f"perm_corrected_p={pcorr:.4f}", flush=True)
        rows.append({"biomarker":b, "best_cutoff":c_best, "naive_best_p":p_best,
                     "best_stat":stat_best, "median_split_cutoff":med,
                     "median_split_p":med_p,
                     "bootstrap_median_p":float(np.median(boot_ps)),
                     "bootstrap_25th_p":float(np.percentile(boot_ps,25)),
                     "bootstrap_75th_p":float(np.percentile(boot_ps,75)),
                     "fraction_boot_p_lt_0.05":float((boot_ps<0.05).mean()),
                     "permutation_corrected_p":float(pcorr),
                     "n":int(n), "n_events":int(events.sum())})

    out = pd.DataFrame(rows)
    out.to_csv(OUT/"step4f_bootstrap_pvalues.tsv", sep="\t", index=False)
    pd.DataFrame(boot_dists).to_csv(OUT/"step4f_bootstrap_distributions.tsv", sep="\t", index=False)
    pd.DataFrame(null_records).to_csv(OUT/"step4f_permutation_null.tsv", sep="\t", index=False)

    summary = {
        "cohort":"BRAF V600E LGG (OpenPedCan)",
        "n_analyzable":int(len(surv)),
        "n_events":int(surv.event.sum()),
        "biomarkers":BIOMARKERS,
        "method_note":
            ("naive_best_p = minimum log-rank p across admissible cutoffs (>=25% per side); "
             "bootstrap_median_p = Levine 2024 method (1000x90% subsample, median optimal-cut p); "
             "permutation_corrected_p = max log-rank vs 1000 event-shuffle nulls."),
        "rows": out.to_dict(orient="records"),
        "compare_to_step4e": {
            "median_split_TIS_p_step4e": 0.1614,
            "median_split_TIS_p_now":   float(out.loc[out.biomarker=="TIS","median_split_p"].iloc[0]),
            "Levine_reported_TIS_p":    0.047,
        },
    }
    with open(OUT/"step4f_summary.json","w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    print("\n[+] DONE.")

if __name__ == "__main__":
    main()
