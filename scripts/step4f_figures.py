"""Step 4f figures (300 dpi) — uses cached bootstrap distributions + perm null."""
from pathlib import Path
import numpy as np, pandas as pd
from lifelines import KaplanMeierFitter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/sessions/blissful-gifted-dirac/mnt/Open PBTA")
OUT  = ROOT/"output"
FIGDIR = OUT/"figs_step4f"; FIGDIR.mkdir(exist_ok=True, parents=True)

plt.rcParams.update({"figure.dpi":300,"savefig.dpi":300,"font.size":9,
                     "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})

BIOMARKERS = ["TIS","Bagaev_Effector_cells","Bagaev_Checkpoint_inhibitor","Bagaev_TAMs"]
EVENT_TYPES = {"Progressive","Progressive - Metastatic","Recurrence","Recurrence - Metastatic"}

cl = pd.read_csv(OUT/"step4e_v600e_cluster_assignment.tsv", sep="\t")
cl["event"] = cl["EFS_event_type"].isin(EVENT_TYPES).astype(int)
cl["EFS_days"] = pd.to_numeric(cl["EFS_days"], errors="coerce")
cl = cl.dropna(subset=["EFS_days","TIS"]).copy()
cl = cl[cl["EFS_days"]>=0]
bag = pd.read_csv(OUT/"step4f_v600e_bagaev.tsv", sep="\t", index_col=0)
surv = cl.merge(bag.reset_index(), on="Kids_First_Biospecimen_ID", how="inner")

cutoff = pd.read_csv(OUT/"step4f_cutoff_search.tsv", sep="\t")
boots  = pd.read_csv(OUT/"step4f_bootstrap_pvalues.tsv", sep="\t")
boot_dists = pd.read_csv(OUT/"step4f_bootstrap_distributions.tsv", sep="\t")
null_records = pd.read_csv(OUT/"step4f_permutation_null.tsv", sep="\t")

# Figure 1: Three-p summary forest plot
print("[*] Figure 1: three-p summary plot")
fig, ax = plt.subplots(figsize=(7.0, 3.6))
y = np.arange(len(BIOMARKERS))
width = 0.25
for i,(col,lab,c) in enumerate([("naive_best_p","Naïve best (uncorrected)","#D7301F"),
                                  ("bootstrap_median_p","Bootstrap median\n(Levine 2024 method)","#F4A582"),
                                  ("permutation_corrected_p","Permutation-corrected\n(proper multiplicity)","#2C7BB6")]):
    vals = [boots.loc[boots.biomarker==b, col].iloc[0] for b in BIOMARKERS]
    ax.barh(y+(i-1)*width, [-np.log10(v) for v in vals], height=width, color=c, label=lab, edgecolor="white")
    for j,v in enumerate(vals):
        ax.text(-np.log10(v)+0.05, y[j]+(i-1)*width, f"{v:.3f}", va="center", fontsize=7)
ax.axvline(-np.log10(0.05), color="grey", ls=":", lw=0.7)
ax.set_yticks(y); ax.set_yticklabels(BIOMARKERS); ax.invert_yaxis()
ax.set_xlabel("-log10(p-value)")
ax.set_title("BRAF V600E LGG (OpenPedCan, n=52, 22 events) — TIS-PFS biomarker p-values under 3 correction strategies", fontsize=9)
ax.legend(loc="lower right", fontsize=7, frameon=False, ncol=1)
fig.tight_layout()
fig.savefig(FIGDIR/"step4f_three_pvalue_summary.png", dpi=300, bbox_inches="tight")
fig.savefig(FIGDIR/"step4f_three_pvalue_summary.pdf", bbox_inches="tight"); plt.close(fig)

# Figure 2: bootstrap p distributions
print("[*] Figure 2: bootstrap p distributions")
fig, axes = plt.subplots(2,2,figsize=(8.5,6.0))
for ax, b in zip(axes.flat, BIOMARKERS):
    boot_ps = boot_dists[b].values
    row = boots[boots.biomarker==b].iloc[0]
    ax.hist(boot_ps, bins=40, color="#2C7BB6", alpha=0.85, edgecolor="white")
    ax.axvline(row.bootstrap_median_p, color="#D7301F", ls="--", lw=0.9, label=f"median = {row.bootstrap_median_p:.3f}")
    ax.axvline(0.05, color="grey", ls=":", lw=0.7, label="α = 0.05")
    ax.set_xlabel("Optimal-cutoff log-rank p (each bootstrap)")
    ax.set_ylabel("Count")
    ax.set_title(f"{b}\nnaïve p={row.naive_best_p:.3f}  perm-corr p={row.permutation_corrected_p:.3f}",
                 fontsize=8)
    ax.legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(FIGDIR/"step4f_bootstrap_distributions.png", dpi=300, bbox_inches="tight")
fig.savefig(FIGDIR/"step4f_bootstrap_distributions.pdf", bbox_inches="tight"); plt.close(fig)

# Figure 3: permutation null distribution vs observed stat
print("[*] Figure 3: permutation null vs observed")
fig, axes = plt.subplots(2,2,figsize=(8.5,6.0))
for ax, b in zip(axes.flat, BIOMARKERS):
    null_stats = null_records[b].values
    row = boots[boots.biomarker==b].iloc[0]
    ax.hist(null_stats, bins=40, color="#92C5DE", alpha=0.85, edgecolor="white")
    ax.axvline(row.best_stat, color="#D7301F", ls="--", lw=1.0,
               label=f"observed = {row.best_stat:.2f}")
    ax.set_xlabel("max log-rank χ² (across cutoffs)")
    ax.set_ylabel("Count")
    ax.set_title(f"{b}\nperm-corrected p = {row.permutation_corrected_p:.3f}", fontsize=8)
    ax.legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(FIGDIR/"step4f_permutation_null.png", dpi=300, bbox_inches="tight")
fig.savefig(FIGDIR/"step4f_permutation_null.pdf", bbox_inches="tight"); plt.close(fig)

# Figure 4: KM at optimal cutoff for each biomarker
print("[*] Figure 4: KM at optimal cutoff")
fig, axes = plt.subplots(2,2,figsize=(8.5,6.4))
times  = surv["EFS_days"].values.astype(float)
events = surv["event"].values.astype(int)
kmf = KaplanMeierFitter()
for ax,b in zip(axes.flat, BIOMARKERS):
    row = boots[boots.biomarker==b].iloc[0]
    c_best = row.best_cutoff
    scores = surv[b].values.astype(float)
    hi = scores > c_best
    for grp, mask, col in [(f"high (>{c_best:.2f})", hi, "#D7301F"),
                            (f"low (≤{c_best:.2f})",  ~hi, "#92C5DE")]:
        if mask.sum():
            kmf.fit(times[mask]/365.25, events[mask], label=f"{grp}  n={mask.sum()} ev={int(events[mask].sum())}")
            kmf.plot_survival_function(ax=ax, color=col, ci_show=False, linewidth=1.4)
    ax.set_title(f"{b} optimal-cutoff KM\nnaïve p={row.naive_best_p:.3f}  boot median p={row.bootstrap_median_p:.3f}  perm-corr p={row.permutation_corrected_p:.3f}",
                 fontsize=7.5)
    ax.set_xlabel("Years"); ax.set_ylabel("PFS probability"); ax.set_ylim(0,1.02)
    ax.legend(fontsize=7, loc="lower left", frameon=False)
fig.tight_layout()
fig.savefig(FIGDIR/"step4f_KM_optimal_cutoff.png", dpi=300, bbox_inches="tight")
fig.savefig(FIGDIR/"step4f_KM_optimal_cutoff.pdf", bbox_inches="tight"); plt.close(fig)

print("[+] All figures saved to", FIGDIR)
