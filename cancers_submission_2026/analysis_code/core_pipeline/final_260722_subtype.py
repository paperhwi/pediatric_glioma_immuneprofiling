#!/usr/bin/env python3
"""Final molecular-subtype sensitivity analysis using the 34-feature matrix."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "FINAL MANUSCRIPT 260722"
OUT = FINAL / "8. Final rerun outputs" / "01_tables"
SEED = 20260722
NPERM = 4999


def load_statlib():
    p = FINAL / "2. Scripts" / "Revision analyses" / "statlib_nodep.py"
    spec = importlib.util.spec_from_file_location("statlib_final", p)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def main() -> None:
    sl = load_statlib()
    z = pd.read_csv(OUT / "TPM_34feature_matrix_z.tsv", sep="\t", index_col=0)
    ann = pd.read_csv(OUT / "sample_master_annotation_final.tsv", sep="\t").set_index("sample").loc[z.index]
    rows = []
    for a, b, an, bn in [
        (ann["cohort_group"].values, ann["ecotype"].values, "molecular_group", "ecotype"),
        (ann["ecotype"].values, ann["cohort_group"].values, "ecotype", "molecular_group"),
    ]:
        result = sl.permanova2(z.values, a, b, an, bn, B=NPERM, seed=SEED)
        for _, r in result.iterrows():
            rows.append({"n": len(z), "n_features": 34, "order": f"{an} + {bn}", **r.to_dict()})
    perm = pd.DataFrame(rows)
    perm.to_csv(OUT / "PERMANOVA_34feature_ecotype_vs_subtype.tsv", sep="\t", index=False)
    wanted = [
        ("ecotype + molecular_group", "ecotype", "Ecotype\nfirst", "#2563EB"),
        ("ecotype + molecular_group", "molecular_group | ecotype", "Molecular group\n| ecotype", "#A78BFA"),
        ("molecular_group + ecotype", "molecular_group", "Molecular group\nfirst", "#7C3AED"),
        ("molecular_group + ecotype", "ecotype | molecular_group", "Ecotype\n| molecular group", "#60A5FA"),
    ]
    vals=[]; pvals=[]; labels=[]; colors=[]
    for order, term, label, color in wanted:
        r=perm[(perm["order"]==order)&(perm["term"]==term)].iloc[0]
        vals.append(float(r.partial_R2)); pvals.append(float(r.p)); labels.append(label); colors.append(color)
    fig, ax=plt.subplots(figsize=(8.5,5.2))
    ax.bar(np.arange(4), vals, color=colors)
    ax.set_xticks(np.arange(4), labels)
    ax.set_ylabel("PERMANOVA partial R²")
    ax.set_title("Ecotype and molecular-group variance partition\n34-feature TPM matrix; n=349; 4,999 permutations")
    ax.set_ylim(0,max(vals)*1.25)
    for i,(v,p) in enumerate(zip(vals,pvals)):
        ax.text(i,v+.006,f"R²={v:.3f}\nP={p:.4g}",ha="center",fontsize=8)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    figdir=FINAL/"8. Final rerun outputs"/"02_figures"
    fig.savefig(figdir/"FigureS5.png",dpi=300,bbox_inches="tight")
    fig.savefig(figdir/"FigureS5.pdf",bbox_inches="tight")
    plt.close(fig)

    d = ann.copy()
    d["time"] = pd.to_numeric(d["OS_days"], errors="coerce")
    d["event"] = d["OS_status"].astype(str).str.upper().eq("DECEASED").astype(int)
    d["age"] = pd.to_numeric(d["age_years"], errors="coerce")
    d["male"] = d["reported_gender"].astype(str).str.lower().eq("male").astype(int)
    d = d[d["time"].notna() & (d["time"] > 0) & d["age"].notna()].copy()
    eco = pd.Categorical(d["ecotype"], categories=["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"])
    coh = pd.Categorical(d["cohort_group"], categories=["pHGG_WT", "DMG_K27", "DHG_G34", "IHG"])
    base = pd.get_dummies(pd.DataFrame({"ecotype": eco, "cohort_group": coh}), drop_first=True).astype(float)
    base["age"] = d["age"].values
    base["male"] = d["male"].values
    base["time"] = d["time"].values
    base["event"] = d["event"].values
    full = base.copy()
    eco_cols = [c for c in base if c.startswith("ecotype_")]
    coh_cols = [c for c in base if c.startswith("cohort_group_")]
    for e in eco_cols:
        for c in coh_cols:
            full[f"{e} x {c}"] = base[e] * base[c]
    c0 = CoxPHFitter().fit(base, "time", "event")
    c1 = CoxPHFitter().fit(full, "time", "event")
    lr = 2 * (c1.log_likelihood_ - c0.log_likelihood_)
    df = len(full.columns) - len(base.columns)
    interaction = {
        "n": len(d), "events": int(d["event"].sum()), "chi2": float(lr),
        "df": df, "p": float(stats.chi2.sf(lr, df)),
        "interpretation": "No interaction was detected; this is not evidence that effects generalize across subtypes.",
    }
    (OUT / "ecotype_subtype_interaction_final.json").write_text(json.dumps(interaction, indent=2), encoding="utf-8")

    subgroup = []
    for cohort, g in d.groupby("cohort_group"):
        q = g[g["ecotype"].isin(["Lymphocyte-inflamed", "Immune-desert"])].copy()
        if len(q) < 8 or q["ecotype"].nunique() < 2 or q["event"].sum() < 4:
            continue
        X = pd.DataFrame({
            "Immune_desert_vs_Lymphocyte_inflamed": (q["ecotype"] == "Immune-desert").astype(int),
            "age": q["age"], "male": q["male"], "time": q["time"], "event": q["event"],
        })
        try:
            c = CoxPHFitter().fit(X, "time", "event")
            r = c.summary.loc["Immune_desert_vs_Lymphocyte_inflamed"]
            subgroup.append({
                "cohort_group": cohort, "n": len(q), "events": int(q["event"].sum()),
                "HR": float(np.exp(r["coef"])),
                "CI_low": float(np.exp(r["coef lower 95%"])),
                "CI_high": float(np.exp(r["coef upper 95%"])), "p": float(r["p"]),
            })
        except Exception as exc:
            subgroup.append({"cohort_group": cohort, "n": len(q), "events": int(q["event"].sum()), "error": str(exc)})
    pd.DataFrame(subgroup).to_csv(OUT / "subtype_stratified_Cox_final.tsv", sep="\t", index=False)
    print(json.dumps({"interaction": interaction, "subgroups": subgroup}, indent=2))


if __name__ == "__main__":
    main()
