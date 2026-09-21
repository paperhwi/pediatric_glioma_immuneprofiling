"""
Step 9 — Exploratory survival on new LM22+ssGSEA ecotype (Main cohort, k=3)

- KM curves overall + per cohort_group with log-rank p
- Multivariable Cox: ecotype + cohort_group + age_years + sex
- Median OS per ecotype (and per cohort × ecotype)

Outputs (output/, figs_step9/)
- step9_KM_logrank_overall.tsv
- step9_KM_logrank_by_cohort.tsv
- step9_median_OS_by_ecotype.tsv
- step9_cox_multivariable.tsv
- figs_step9/KM_overall_by_ecotype.png
- figs_step9/KM_by_cohort_grid.png
- figs_step9/Cox_forest.png
"""
from __future__ import annotations
from pathlib import Path
import json, warnings
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_step9"; FIG.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# 1) Load + assemble survival table
# ---------------------------------------------------------------------------
ann = pd.read_csv(OUT / "cohort_for_deconvolution.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
ec3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
df = ann.join(ec3[["ecotype"]], how="inner")
df["event"] = (df["OS_status"].astype(str).str.upper() == "DECEASED").astype(int)
df["time"]  = pd.to_numeric(df["OS_days"], errors="coerce")
df["age_years"] = pd.to_numeric(df["age_years"], errors="coerce")
df = df.dropna(subset=["time", "event", "ecotype"])
print(f"Step 9: n={len(df)} with OS data and ecotype")
print(df["ecotype"].value_counts())
print(f"events: {df['event'].sum()} / {len(df)}")

ec_order  = ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]
palette   = {"Lymphocyte-inflamed": "#3B82F6", "Myeloid-dominant": "#EF4444", "Immune-desert": "#9CA3AF"}
cohort_order = ["DMG_K27", "DHG_G34", "pHGG_WT", "IHG"]

# ---------------------------------------------------------------------------
# 2) Median OS per ecotype (months); KM overall
# ---------------------------------------------------------------------------
med_rows = []
for ec, g in df.groupby("ecotype"):
    kmf = KaplanMeierFitter().fit(g["time"], g["event"], label=ec)
    med = float(kmf.median_survival_time_) if not np.isnan(kmf.median_survival_time_) else None
    med_rows.append({"ecotype": ec, "n": len(g), "events": int(g["event"].sum()),
                     "median_OS_days": med,
                     "median_OS_months": (med/30.4375) if med is not None else None})
pd.DataFrame(med_rows).to_csv(OUT / "step9_median_OS_by_ecotype.tsv", sep="\t", index=False)
print(pd.DataFrame(med_rows).to_string(index=False))

# Log-rank overall
lrt = multivariate_logrank_test(df["time"], df["ecotype"], df["event"])
overall = pd.DataFrame([{"test": "multivariate_logrank_3groups",
                         "test_statistic": float(lrt.test_statistic),
                         "p_value": float(lrt.p_value),
                         "n": int(len(df)),
                         "events": int(df["event"].sum())}])
overall.to_csv(OUT / "step9_KM_logrank_overall.tsv", sep="\t", index=False)
print(overall.to_string(index=False))

# Pairwise log-rank
pair_rows = []
for i, ei in enumerate(ec_order):
    for ej in ec_order[i+1:]:
        a = df[df["ecotype"] == ei]; b = df[df["ecotype"] == ej]
        if a.empty or b.empty: continue
        r = logrank_test(a["time"], b["time"], a["event"], b["event"])
        pair_rows.append({"A": ei, "B": ej, "stat": float(r.test_statistic),
                          "p": float(r.p_value), "n_A": len(a), "n_B": len(b)})
pd.DataFrame(pair_rows).to_csv(OUT / "step9_KM_logrank_pairwise.tsv", sep="\t", index=False)

# Plot KM overall
fig, ax = plt.subplots(figsize=(6.5, 4.8))
for ec in ec_order:
    g = df[df["ecotype"] == ec]
    if g.empty: continue
    KaplanMeierFitter().fit(g["time"], g["event"],
        label=f"{ec} (n={len(g)}, ev={int(g['event'].sum())})"
    ).plot_survival_function(ax=ax, color=palette[ec], ci_show=False)
ax.set_title(f"Overall survival by ecotype  (log-rank p={lrt.p_value:.2e})")
ax.set_xlabel("OS (days)"); ax.set_ylabel("Survival probability"); ax.set_ylim(-.02, 1.02)
ax.legend(loc="best", frameon=False, fontsize=8)
plt.tight_layout(); fig.savefig(FIG / "KM_overall_by_ecotype.png", dpi=200); plt.close()

# ---------------------------------------------------------------------------
# 3) KM per cohort_group (2x2 grid)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
cohort_rows = []
for ax, cg in zip(axes.flat, cohort_order):
    sub = df[df["cohort_group"] == cg]
    if len(sub) < 5:
        ax.set_visible(False); continue
    for ec in ec_order:
        g = sub[sub["ecotype"] == ec]
        if g.empty or g["event"].sum() == 0:
            continue
        KaplanMeierFitter().fit(g["time"], g["event"],
            label=f"{ec} (n={len(g)}, ev={int(g['event'].sum())})"
        ).plot_survival_function(ax=ax, color=palette[ec], ci_show=False)
    if sub["ecotype"].nunique() >= 2:
        lrt_c = multivariate_logrank_test(sub["time"], sub["ecotype"], sub["event"])
        p_c   = float(lrt_c.p_value)
    else:
        p_c   = float("nan")
    cohort_rows.append({"cohort": cg, "n": len(sub), "events": int(sub["event"].sum()),
                        "logrank_p": p_c})
    ax.set_title(f"{cg}  (n={len(sub)}, log-rank p={p_c:.2e})", fontsize=10)
    ax.set_xlabel("OS (days)"); ax.set_ylabel("Survival prob")
    ax.set_ylim(-.02, 1.02); ax.legend(loc="best", frameon=False, fontsize=7)
pd.DataFrame(cohort_rows).to_csv(OUT / "step9_KM_logrank_by_cohort.tsv", sep="\t", index=False)
plt.tight_layout(); fig.savefig(FIG / "KM_by_cohort_grid.png", dpi=200); plt.close()

# ---------------------------------------------------------------------------
# 4) Multivariable Cox: ecotype + cohort_group + age + sex
# ---------------------------------------------------------------------------
cox = df[["time", "event", "ecotype", "cohort_group", "age_years", "reported_gender"]].copy()
cox = cox.dropna()
# Reference levels: Immune-desert + pHGG_WT + Male
cox["ecotype"]      = pd.Categorical(cox["ecotype"], categories=["Immune-desert", "Lymphocyte-inflamed", "Myeloid-dominant"])
cox["cohort_group"] = pd.Categorical(cox["cohort_group"], categories=[c for c in ["pHGG_WT","DMG_K27","DHG_G34","IHG"] if c in cox["cohort_group"].unique()])
cox["male"]         = (cox["reported_gender"].astype(str).str.lower() == "male").astype(int)
cox = cox.drop(columns=["reported_gender"])
cox_design = pd.get_dummies(cox, columns=["ecotype", "cohort_group"], drop_first=True)
# Convert any bool/object to float for lifelines
for c in cox_design.columns:
    if cox_design[c].dtype == bool:
        cox_design[c] = cox_design[c].astype(int)
cph = CoxPHFitter()
cph.fit(cox_design.astype(float), duration_col="time", event_col="event", show_progress=False)
sm = cph.summary.reset_index()
sm.to_csv(OUT / "step9_cox_multivariable.tsv", sep="\t", index=False)
print("\nCox multivariable summary:")
print(sm[["covariate","exp(coef)","exp(coef) lower 95%","exp(coef) upper 95%","p"]].to_string(index=False))

# Forest plot
sm["lab"] = sm["covariate"].str.replace("ecotype_","ecotype: ").str.replace("cohort_group_","cohort: ").str.replace("_"," ")
fig, ax = plt.subplots(figsize=(7.5, 0.4*len(sm)+1.5))
for i, r in sm.iterrows():
    color = "red" if r["p"] < 0.05 else "black"
    ax.errorbar(r["exp(coef)"], i, xerr=[[r["exp(coef)"]-r["exp(coef) lower 95%"]],
                                          [r["exp(coef) upper 95%"]-r["exp(coef)"]]],
                fmt="o", color=color, ecolor=color, capsize=3)
ax.axvline(1, color="grey", ls="--", alpha=.5)
ax.set_yticks(range(len(sm))); ax.set_yticklabels(sm["lab"], fontsize=9)
ax.set_xscale("log")
ax.set_xlabel("Hazard ratio (95% CI)")
ax.set_title("Multivariable Cox  (reference: Immune-desert, pHGG_WT, Female)")
plt.tight_layout(); fig.savefig(FIG / "Cox_forest.png", dpi=200); plt.close()

# ---------------------------------------------------------------------------
# 5) Summary
# ---------------------------------------------------------------------------
summary = {
    "n_total_with_OS": int(len(df)),
    "n_events": int(df["event"].sum()),
    "logrank_overall_p": float(lrt.p_value),
    "median_OS_days_by_ecotype": {r["ecotype"]: r["median_OS_days"] for r in med_rows},
    "cohort_logrank_p": {r["cohort"]: r["logrank_p"] for r in cohort_rows},
    "cox_n": int(len(cox_design)),
    "cox_top_HRs_p05": sm[sm["p"]<0.05][["covariate","exp(coef)","p"]].to_dict(orient="records"),
}
(OUT / "step9_summary.json").write_text(json.dumps(summary, indent=2))
print("\n=== SUMMARY ===")
print(json.dumps(summary, indent=2))
