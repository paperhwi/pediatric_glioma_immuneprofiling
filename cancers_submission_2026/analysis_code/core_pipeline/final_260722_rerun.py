#!/usr/bin/env python3
"""Final 2026-07-22 harmonized rerun for the OpenPBTA manuscript.

Authoritative choices
---------------------
* Fixed ecotype membership: output/ecotype_LM22_main_k3_annotated.tsv
  (111 Lymphocyte-inflamed / 160 Myeloid-dominant / 78 Immune-desert).
* TPM feature space: 10 quanTIseq cell fractions + all 24 brain-tuned ssGSEA
  signatures.  No positional slicing of the signature list is allowed.
* Random seed: 20260722.

The script writes a self-contained audit bundle under
FINAL MANUSCRIPT 260722/8. Final rerun outputs and does not overwrite the
original analysis outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from lifelines import CoxPHFitter
from scipy import stats
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.preprocessing import StandardScaler


SEED = 20260722
ECOTYPES = ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]
COLORS = {
    "Lymphocyte-inflamed": "#3B82F6",
    "Myeloid-dominant": "#EF4444",
    "Immune-desert": "#9CA3AF",
}
COHORTS = ["DMG_K27", "DHG_G34", "pHGG_WT", "IHG"]
LOCATIONS = ["Midline", "Hemispheric", "Posterior_fossa", "Ambiguous"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--consensus-bootstrap", type=int, default=500)
    parser.add_argument("--permutations", type=int, default=4999)
    parser.add_argument("--mc-permutations", type=int, default=99999)
    return parser.parse_args()


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.clip(ranked, 0, 1)
    return out


def safe_name(value: str) -> str:
    return value.replace(" ", "_").replace("-", "_").replace("/", "_")


def save_fig(fig: plt.Figure, base: Path) -> None:
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def load_statlib(root: Path):
    import importlib.util
    path = root / "FINAL MANUSCRIPT 260722" / "2. Scripts" / "Revision analyses" / "statlib_nodep.py"
    spec = importlib.util.spec_from_file_location("statlib_final", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def canonical_annotation(root: Path, tables: Path, backup: Path) -> pd.DataFrame:
    final = root / "FINAL MANUSCRIPT 260722"
    old_path = final / "7. Reproducibility data" / "Sample annotation" / "sample_master_annotation.tsv"
    fixed_path = root / "output" / "ecotype_LM22_main_k3_annotated.tsv"
    old = pd.read_csv(old_path, sep="\t")
    fixed = pd.read_csv(fixed_path, sep="\t").rename(
        columns={"Kids_First_Biospecimen_ID": "sample", "ecotype": "ecotype_final"}
    )
    assert len(fixed) == 349 and fixed["sample"].is_unique
    counts = fixed["ecotype_final"].value_counts().reindex(ECOTYPES).to_dict()
    assert counts == {"Lymphocyte-inflamed": 111, "Myeloid-dominant": 160, "Immune-desert": 78}
    ann = old.drop(columns=["ecotype"], errors="ignore").merge(
        fixed[["sample", "ecotype_final"]], on="sample", how="inner", validate="one_to_one"
    ).rename(columns={"ecotype_final": "ecotype"})
    assert len(ann) == 349 and ann["ecotype"].notna().all()

    hist = pd.read_csv(
        root / "data" / "histologies.tsv", sep="\t", low_memory=False,
        usecols=["Kids_First_Biospecimen_ID", "extent_of_tumor_resection"],
    ).rename(columns={"Kids_First_Biospecimen_ID": "sample"}).drop_duplicates("sample")
    ann = ann.merge(hist, on="sample", how="left", validate="one_to_one")
    ann["OS_days"] = pd.to_numeric(ann["OS_days"], errors="coerce")
    ann["observed"] = (ann["OS_days"].notna() & (ann["OS_days"] > 0)).astype(int)
    ann["assignment_version"] = "final_111_160_78_20260722"
    ann["assignment_source"] = "output/ecotype_LM22_main_k3_annotated.tsv"
    ann.to_csv(tables / "sample_master_annotation_final.tsv", sep="\t", index=False)
    fixed.to_csv(tables / "canonical_ecotype_assignment_111_160_78.tsv", sep="\t", index=False)

    # Preserve the stale packaged annotation for traceability.
    backup.mkdir(parents=True, exist_ok=True)
    if not (backup / old_path.name).exists():
        shutil.copy2(old_path, backup / old_path.name)
    return ann


def build_34_features(root: Path, ann: pd.DataFrame, tables: Path, h5dir: Path,
                      B: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    scores = pd.read_csv(
        root / "output" / "revC_TPM_harmonization" / "ssGSEA_TPM_NES_brain24.tsv", sep="\t"
    ).set_index("Kids_First_Biospecimen_ID")
    assert scores.shape[1] == 24, f"Expected 24 ssGSEA features, found {scores.shape[1]}"
    qt = pd.read_csv(root / "output" / "immunedeconv_quantiseq.tsv", sep="\t").set_index("cell_type").T
    qt.index.name = "Kids_First_Biospecimen_ID"
    qt = qt.apply(pd.to_numeric, errors="coerce")
    qt = qt.drop(columns=["uncharacterized cell"], errors="ignore")
    assert qt.shape[1] == 10, f"Expected 10 quanTIseq features, found {qt.shape[1]}"
    ids = ann["sample"].tolist()
    q = qt.reindex(ids).add_prefix("quanTIseq_")
    s = scores.reindex(ids).add_prefix("ssGSEA_")
    raw = q.join(s)
    assert raw.shape == (349, 34) and raw.notna().all().all()
    z = pd.DataFrame(StandardScaler().fit_transform(raw), index=raw.index, columns=raw.columns)
    raw.to_csv(tables / "TPM_34feature_matrix_raw.tsv", sep="\t")
    z.to_csv(tables / "TPM_34feature_matrix_z.tsv", sep="\t")

    n, k = len(z), 3
    rng = np.random.default_rng(SEED)
    co_assoc = np.zeros((n, n), dtype=np.float32)
    co_count = np.zeros((n, n), dtype=np.float32)
    all_idx = np.arange(n)
    for b in range(B):
        sub = rng.choice(all_idx, size=int(round(0.8 * n)), replace=False)
        lab = KMeans(n_clusters=k, n_init=10, random_state=b).fit_predict(z.values[sub])
        co_count[np.ix_(sub, sub)] += 1
        for cluster in range(k):
            members = sub[lab == cluster]
            co_assoc[np.ix_(members, members)] += 1
    consensus = np.divide(co_assoc, co_count, out=np.zeros_like(co_assoc), where=co_count > 0)
    new_cluster = KMeans(n_clusters=3, n_init=50, random_state=SEED).fit_predict(1 - consensus)
    canonical = ann.set_index("sample").loc[z.index, "ecotype"]
    ct = pd.crosstab(pd.Series(new_cluster, index=z.index, name="cluster_34"), canonical)
    rr, cc = linear_sum_assignment(-ct.values)
    mapping = {int(ct.index[r]): str(ct.columns[c]) for r, c in zip(rr, cc)}
    new_eco = pd.Series(new_cluster, index=z.index).map(mapping)
    assignment = pd.DataFrame({
        "sample": z.index,
        "canonical_ecotype": canonical.values,
        "TPM34_consensus_cluster": new_cluster,
        "TPM34_mapped_ecotype": new_eco.values,
        "agreement": (canonical.values == new_eco.values).astype(int),
    })
    assignment.to_csv(tables / "TPM34_consensus_vs_canonical.tsv", sep="\t", index=False)
    ct2 = pd.crosstab(assignment["TPM34_mapped_ecotype"], assignment["canonical_ecotype"])
    ct2.to_csv(tables / "TPM34_consensus_crosstab.tsv", sep="\t")
    ari = adjusted_rand_score(canonical, new_eco)
    nmi = normalized_mutual_info_score(canonical, new_eco)
    agreement = float((canonical.values == new_eco.values).mean())
    metrics = {
        "n": n, "n_features": 34, "n_quantiseq": 10, "n_ssgsea": 24,
        "consensus_bootstrap": B, "subsample_fraction": 0.8,
        "ARI_vs_canonical": float(ari), "NMI_vs_canonical": float(nmi),
        "per_sample_agreement": agreement, "cluster_mapping": mapping,
        "canonical_counts": canonical.value_counts().reindex(ECOTYPES).to_dict(),
        "TPM34_counts": new_eco.value_counts().reindex(ECOTYPES).to_dict(),
    }
    (tables / "TPM34_summary.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    pca = PCA(n_components=2).fit_transform(z)
    obs = ann.set_index("sample").loc[z.index].copy()
    obs["TPM34_mapped_ecotype"] = new_eco.astype(str)
    obs["TPM34_consensus_cluster"] = new_cluster.astype(str)
    for c in obs.columns:
        if obs[c].dtype == object:
            obs[c] = obs[c].fillna("").astype(str)
    adata = ad.AnnData(X=raw.astype(np.float32).values, obs=obs,
                       var=pd.DataFrame(index=raw.columns.astype(str)))
    adata.layers["z_score"] = z.astype(np.float32).values
    adata.obsm["X_pca"] = pca.astype(np.float32)
    adata.uns["TPM34_metrics"] = {k: v for k, v in metrics.items() if k != "cluster_mapping"}
    adata.write_h5ad(h5dir / "sample_master_final_34feature.h5ad")
    return raw, z, assignment, metrics


def contingency_and_permanova(root: Path, ann: pd.DataFrame, z: pd.DataFrame,
                              tables: Path, B: int, Bmc: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    sl = load_statlib(root)
    d = ann.set_index("sample").loc[z.index].copy()
    d["ecotype"] = pd.Categorical(d["ecotype"], ECOTYPES, ordered=True)
    rows = []

    def one_ct(data: pd.DataFrame, row: str, label: str):
        ct = pd.crosstab(data[row], data["ecotype"]).reindex(columns=ECOTYPES).fillna(0).astype(int)
        chi2, df, pasym, expected = sl.chisq_contingency(ct.values)
        pmc, _, _ = sl.mc_permutation_p(ct.values, B=Bmc, seed=SEED)
        rec = {
            "analysis": label, "n": int(ct.values.sum()), "chi2": chi2, "df": df,
            "p_asymptotic": pasym, "p_montecarlo": pmc, "n_permutations": Bmc,
            "cramers_V": sl.cramers_v(chi2, int(ct.values.sum()), *ct.shape),
            "min_expected": float(expected.min()),
        }
        rows.append(rec)
        return ct, rec

    full = d[d["location_class"].isin(LOCATIONS)]
    ct_loc, loc_stats = one_ct(full, "location_class", "Location (4 categories) x ecotype")
    ct_mol, mol_stats = one_ct(d, "cohort_group", "Molecular group x ecotype")
    wt = d[(d["cohort_group"] == "pHGG_WT") & d["location_class"].isin(["Midline", "Hemispheric"])]
    ct_wt, wt_stats = one_ct(wt, "location_class", "pHGG_WT: Midline vs Hemispheric x ecotype")
    pd.DataFrame(rows).to_csv(tables / "contingency_final.tsv", sep="\t", index=False)
    ct_loc.to_csv(tables / "location_x_ecotype_counts.tsv", sep="\t")
    ct_mol.to_csv(tables / "molecular_x_ecotype_counts.tsv", sep="\t")
    ct_wt.to_csv(tables / "pHGG_location_x_ecotype_counts.tsv", sep="\t")

    m = d["location_class"].isin(["Midline", "Hemispheric", "Posterior_fossa"])
    X = z.loc[m.values].values
    loc = d.loc[m, "location_class"].astype(str).values
    mol = d.loc[m, "cohort_group"].astype(str).values
    perm = []
    for a, b, an, bn in [(mol, loc, "molecular_group", "location"),
                          (loc, mol, "location", "molecular_group")]:
        res = sl.permanova2(X, a, b, an, bn, B=B, seed=SEED)
        for _, r in res.iterrows():
            perm.append({"n": len(X), "n_features": 34, "order": f"{an} + {bn}", **r.to_dict()})
    perm_df = pd.DataFrame(perm)
    perm_df.to_csv(tables / "PERMANOVA_34feature_final.tsv", sep="\t", index=False)
    return pd.DataFrame(rows), perm_df


def survival_design(ann: pd.DataFrame, reference: str) -> pd.DataFrame:
    d = ann.copy()
    d["time"] = pd.to_numeric(d["OS_days"], errors="coerce")
    d["event"] = d["OS_status"].astype(str).str.upper().eq("DECEASED").astype(int)
    d["age"] = pd.to_numeric(d["age_years"], errors="coerce")
    d["male"] = d["reported_gender"].astype(str).str.lower().eq("male").astype(int)
    d = d[d["time"].notna() & (d["time"] > 0) & d["age"].notna()].copy()
    order = [reference] + [e for e in ECOTYPES if e != reference]
    eco = pd.Categorical(d["ecotype"], categories=order)
    coh = pd.Categorical(d["cohort_group"], categories=["pHGG_WT", "DMG_K27", "DHG_G34", "IHG"])
    X = pd.get_dummies(pd.DataFrame({"ecotype": eco, "cohort_group": coh}), drop_first=True).astype(float)
    X["age"] = d["age"].values
    X["male"] = d["male"].values
    X["time"] = d["time"].values
    X["event"] = d["event"].values
    X.index = d.index
    return X


def cox_table(cph: CoxPHFitter) -> pd.DataFrame:
    s = cph.summary
    return pd.DataFrame({
        "term": s.index,
        "coef": s["coef"].values,
        "HR": np.exp(s["coef"].values),
        "CI_low": np.exp(s["coef lower 95%"].values),
        "CI_high": np.exp(s["coef upper 95%"].values),
        "p": s["p"].values,
    })


def main_cox_and_ipw(ann: pd.DataFrame, tables: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    X = survival_design(ann, "Immune-desert")
    cph = CoxPHFitter().fit(X, duration_col="time", event_col="event")
    main = cox_table(cph)
    main.to_csv(tables / "Cox_main_final.tsv", sep="\t", index=False)

    d = ann.copy()
    d["time"] = pd.to_numeric(d["OS_days"], errors="coerce")
    d["event"] = d["OS_status"].astype(str).str.upper().eq("DECEASED").astype(int)
    d["age"] = pd.to_numeric(d["age_years"], errors="coerce")
    d["male"] = d["reported_gender"].astype(str).str.lower().eq("male").astype(int)
    d["observed"] = (d["time"].notna() & (d["time"] > 0)).astype(int)
    prop = pd.get_dummies(d[["ecotype", "cohort_group", "location_class"]].fillna("Unknown"), drop_first=True).astype(float)
    prop["age"] = d["age"].fillna(d["age"].median()).values
    prop["male"] = d["male"].values
    scale = StandardScaler().fit_transform(prop)
    clf = LogisticRegression(penalty="l2", C=1.0, max_iter=5000, random_state=SEED).fit(scale, d["observed"])
    pobs = clf.predict_proba(scale)[:, 1]
    d["propensity_observed"] = pobs
    d["ipw"] = d["observed"].mean() / pobs
    lo, hi = np.percentile(d["ipw"], [1, 99])
    d["ipw_trimmed"] = d["ipw"].clip(lo, hi)

    # Balance between observed and missing before/after weighting.
    A = prop.values[d["observed"].values == 1]
    M = prop.values[d["observed"].values == 0]
    w = d.loc[d["observed"] == 1, "ipw_trimmed"].values
    def smd(a, b, weights=None):
        if weights is None:
            ma, va = a.mean(0), a.var(0)
        else:
            wn = weights / weights.sum()
            ma = a.T @ wn
            va = ((a - ma) ** 2).T @ wn
        mb, vb = b.mean(0), b.var(0)
        den = np.sqrt((va + vb) / 2)
        den[den == 0] = 1e-12
        return (ma - mb) / den
    balance = pd.DataFrame({
        "covariate": prop.columns,
        "SMD_unweighted": smd(A, M),
        "SMD_IPW_weighted": smd(A, M, w),
    })
    balance["abs_SMD_unweighted"] = balance["SMD_unweighted"].abs()
    balance["abs_SMD_IPW_weighted"] = balance["SMD_IPW_weighted"].abs()
    balance.to_csv(tables / "IPW_balance_final.tsv", sep="\t", index=False)

    obs = d[d["observed"] == 1].copy()
    eco = pd.Categorical(obs["ecotype"], categories=["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"])
    coh = pd.Categorical(obs["cohort_group"], categories=["pHGG_WT", "DMG_K27", "DHG_G34", "IHG"])
    Xw = pd.get_dummies(pd.DataFrame({"ecotype": eco, "cohort_group": coh}), drop_first=True).astype(float)
    Xw["age"] = obs["age"].fillna(obs["age"].median()).values
    Xw["male"] = obs["male"].values
    Xw["time"] = obs["time"].values
    Xw["event"] = obs["event"].values
    Xw["weight"] = obs["ipw_trimmed"].values
    naive = CoxPHFitter().fit(Xw.drop(columns="weight"), "time", "event")
    weighted = CoxPHFitter().fit(Xw, "time", "event", weights_col="weight", robust=True)
    ndf, wdf = cox_table(naive), cox_table(weighted)
    compare = ndf.merge(wdf, on="term", suffixes=("_naive", "_IPW"))
    compare.to_csv(tables / "Cox_naive_vs_IPW_final.tsv", sep="\t", index=False)
    summary = {
        "n_total": len(d), "n_observed": int(d["observed"].sum()),
        "missingness_percent": float(100 * (1 - d["observed"].mean())),
        "weight_raw_min": float(d["ipw"].min()), "weight_raw_max": float(d["ipw"].max()),
        "weight_trim_limits": [float(lo), float(hi)],
        "n_abs_SMD_gt_0.1_before": int((balance["abs_SMD_unweighted"] > 0.1).sum()),
        "n_abs_SMD_gt_0.1_after": int((balance["abs_SMD_IPW_weighted"] > 0.1).sum()),
        "max_abs_SMD_after": float(balance["abs_SMD_IPW_weighted"].max()),
        "interpretation": "IPW addresses missingness conditional on measured covariates (MAR-type assumption); it does not establish robustness to MNAR.",
    }
    (tables / "IPW_summary_final.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return main, compare, balance, summary


def eor_analysis(ann: pd.DataFrame, tables: Path) -> tuple[pd.DataFrame, dict]:
    d = ann.copy()
    d["time"] = pd.to_numeric(d["OS_days"], errors="coerce")
    d["event"] = d["OS_status"].astype(str).str.upper().eq("DECEASED").astype(int)
    d["age"] = pd.to_numeric(d["age_years"], errors="coerce")
    d["male"] = d["reported_gender"].astype(str).str.lower().eq("male").astype(int)
    d["EOR"] = d["extent_of_tumor_resection"].replace("Unavailable", np.nan)
    d = d[d["time"].notna() & (d["time"] > 0) & d["age"].notna()].copy()

    def design(x: pd.DataFrame, include_eor: bool) -> pd.DataFrame:
        eco = pd.Categorical(x["ecotype"], categories=["Immune-desert", "Lymphocyte-inflamed", "Myeloid-dominant"])
        coh = pd.Categorical(x["cohort_group"], categories=["pHGG_WT", "DMG_K27", "DHG_G34", "IHG"])
        parts = pd.get_dummies(pd.DataFrame({"ecotype": eco, "cohort_group": coh}), drop_first=True).astype(float)
        if include_eor:
            eor = pd.Categorical(x["EOR"], categories=["Biopsy only", "Partial resection", "Gross/Near total resection"])
            parts = parts.join(pd.get_dummies(eor, prefix="EOR", drop_first=True).astype(float).set_axis(parts.index))
        parts["age"] = x["age"].values
        parts["male"] = x["male"].values
        parts["time"] = x["time"].values
        parts["event"] = x["event"].values
        return parts

    eor = d[d["EOR"].notna()].copy()
    models = {}
    for name, data, inc in [("A_full_no_EOR", d, False), ("B_EOR_subset_no_EOR", eor, False), ("C_EOR_subset_with_EOR", eor, True)]:
        c = CoxPHFitter().fit(design(data, inc), "time", "event")
        t = cox_table(c)
        t["model"] = name
        t["n"] = len(data)
        models[name] = (c, t)
    out = pd.concat([v[1] for v in models.values()], ignore_index=True)
    out.to_csv(tables / "EOR_Cox_final.tsv", sep="\t", index=False)
    cB, cC = models["B_EOR_subset_no_EOR"][0], models["C_EOR_subset_with_EOR"][0]
    lr = 2 * (cC.log_likelihood_ - cB.log_likelihood_)
    summary = {
        "n_full": len(d), "n_EOR_subset": len(eor),
        "events_full": int(d["event"].sum()), "events_EOR_subset": int(eor["event"].sum()),
        "likelihood_ratio_chi2": float(lr), "likelihood_ratio_df": 2,
        "likelihood_ratio_p": float(stats.chi2.sf(lr, 2)),
    }
    (tables / "EOR_summary_final.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return out, summary


def adjusted_deg(root: Path, ann: pd.DataFrame, tables: Path, figures: Path) -> tuple[pd.DataFrame, dict]:
    ids = ann["sample"].tolist()
    use = ["GeneSymbol", *ids]
    # Prefer the complete uncompressed matrix.  The tracked gzip copy in this
    # workspace is truncated (EOF before the gzip end-of-stream marker).
    tpm_path = root / "output" / "tpm_for_cibersortx.tsv"
    if not tpm_path.exists():
        tpm_path = root / "output" / "tpm_for_cibersortx.tsv.gz"
    tpm = pd.read_csv(tpm_path, sep="\t", usecols=use).set_index("GeneSymbol")
    tpm = tpm[~tpm.index.duplicated(keep="first")]
    expr = (tpm >= 1).sum(axis=1) >= int(np.ceil(0.05 * len(ids)))
    tpm = tpm.loc[expr]
    g = pd.Series(tpm.index.astype(str), index=tpm.index)
    keep_kinase = g.str.match(r"^RPS6K")
    exclude = (
        g.str.match(r"^MT-") | g.str.match(r"^MTRNR") | g.str.match(r"^MT(ATP|CO|CYB|ND|RNR)\d") |
        (g.str.match(r"^RP[LS]\d") & ~keep_kinase) | g.str.match(r"^RPLP\d") | g.str.match(r"^RPSA$") |
        g.str.match(r"^RNA(5S|5-8S|18S|28S|45S|5-8SN|45SN)") | g.str.contains(r"rRNA") |
        g.str.match(r"^RN(U|VU)\d") | g.str.fullmatch(r"U\d+") |
        g.str.fullmatch(r"7SK|Y_RNA|Vault|Metazoa_SRP|uc_338") | g.str.match(r"^RN7S[KL]") |
        g.str.match(r"^(SNORD|SNORA|SCARNA|SNAR)") | g.str.match(r"^RNY\d") |
        g.str.match(r"^MIR\d") | g.str.match(r"^MIRLET") | g.str.match(r"^VTRNA") |
        g.str.fullmatch(r"RPPH1|RMRP") | g.str.match(r"^ENSG\d")
    )
    tpm = tpm.loc[~exclude]
    Y = np.log2(tpm[ids].T.values + 1.0)
    meta = ann.set_index("sample").loc[ids]
    eco = pd.Categorical(meta["ecotype"], categories=["Immune-desert", "Lymphocyte-inflamed", "Myeloid-dominant"])
    coh = pd.Categorical(meta["cohort_group"], categories=["pHGG_WT", "DMG_K27", "DHG_G34", "IHG"])
    Xdf = pd.get_dummies(pd.DataFrame({"ecotype": eco, "cohort_group": coh}), drop_first=True).astype(float)
    Xdf.insert(0, "Intercept", 1.0)
    X = Xdf.values
    inv = np.linalg.pinv(X.T @ X)
    beta = inv @ X.T @ Y
    resid = Y - X @ beta
    df_resid = X.shape[0] - np.linalg.matrix_rank(X)
    sigma2 = (resid ** 2).sum(axis=0) / df_resid
    contrasts = {
        "Lymphocyte_inflamed_vs_Immune_desert": {"ecotype_Lymphocyte-inflamed": 1.0},
        "Myeloid_dominant_vs_Immune_desert": {"ecotype_Myeloid-dominant": 1.0},
        "Lymphocyte_inflamed_vs_Myeloid_dominant": {"ecotype_Lymphocyte-inflamed": 1.0, "ecotype_Myeloid-dominant": -1.0},
    }
    all_rows = []
    panels = []
    summary = {"n_samples": len(ids), "n_genes_tested": int(tpm.shape[0]), "model": "log2(TPM+1) ~ ecotype + molecular_group", "contrasts": {}}
    for name, weights in contrasts.items():
        c = np.zeros(X.shape[1])
        for term, value in weights.items():
            c[Xdf.columns.get_loc(term)] = value
        effect = c @ beta
        se = np.sqrt(np.maximum(sigma2 * (c @ inv @ c), 1e-300))
        tval = effect / se
        p = 2 * stats.t.sf(np.abs(tval), df_resid)
        q = bh(p)
        out = pd.DataFrame({"gene": tpm.index, "adjusted_log2FC": effect, "t": tval, "p": p, "q_BH": q, "df_resid": df_resid})
        out = out.sort_values(["q_BH", "adjusted_log2FC"], ascending=[True, False])
        out.to_csv(tables / f"adjusted_DEG_{name}.tsv", sep="\t", index=False)
        out["contrast"] = name
        all_rows.append(out)
        n_up = int(((out["q_BH"] < 0.05) & (out["adjusted_log2FC"] > 1)).sum())
        n_down = int(((out["q_BH"] < 0.05) & (out["adjusted_log2FC"] < -1)).sum())
        summary["contrasts"][name] = {"q05_FCgt1_up": n_up, "q05_FClt_minus1_down": n_down}
        panels.append((name, out, n_up, n_down))
    combined = pd.concat(all_rows, ignore_index=True)
    combined.to_csv(tables / "adjusted_DEG_all_contrasts.tsv", sep="\t", index=False)
    (tables / "adjusted_DEG_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, (name, out, n_up, n_down) in zip(axes, panels):
        y = -np.log10(np.clip(out["q_BH"], 1e-300, 1))
        sig = (out["q_BH"] < 0.05) & (out["adjusted_log2FC"].abs() > 1)
        ax.scatter(out.loc[~sig, "adjusted_log2FC"], y.loc[~sig], s=3, c="0.82", alpha=.45, linewidth=0)
        ax.scatter(out.loc[sig, "adjusted_log2FC"], y.loc[sig], s=4, c="#4C78A8", alpha=.7, linewidth=0)
        top = out[(out["q_BH"] < 0.05) & (out["adjusted_log2FC"] > 0)].nlargest(7, "adjusted_log2FC")
        for idx, r in top.iterrows():
            ax.text(r["adjusted_log2FC"], y.loc[idx], str(r["gene"]), fontsize=6)
        ax.axvline(-1, ls="--", lw=.7, c="0.4"); ax.axvline(1, ls="--", lw=.7, c="0.4")
        ax.axhline(-np.log10(.05), ls="--", lw=.7, c="0.4")
        ax.set_title(name.replace("_", "\n"), fontsize=9)
        ax.set_xlabel("Adjusted log2 fold change")
        ax.set_ylabel("-log10(BH q)")
        ax.text(.02, .98, f"up={n_up}; down={n_down}", transform=ax.transAxes, va="top", fontsize=8)
    fig.suptitle("Molecular-group-adjusted differential expression\nlog2(TPM+1) ~ ecotype + molecular group", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, .91])
    save_fig(fig, figures / "FigureS_adjusted_DEG")
    return combined, summary


def figure4(ann: pd.DataFrame, contingency: pd.DataFrame, perm: pd.DataFrame,
            figures: Path) -> None:
    d = ann.copy()
    d["ecotype"] = pd.Categorical(d["ecotype"], ECOTYPES, ordered=True)
    ct_loc = pd.crosstab(d[d["location_class"].isin(LOCATIONS)]["location_class"], d["ecotype"]).reindex(index=LOCATIONS, columns=ECOTYPES)
    ct_mol = pd.crosstab(d["cohort_group"], d["ecotype"]).reindex(index=COHORTS, columns=ECOTYPES)
    wt = d[(d["cohort_group"] == "pHGG_WT") & d["location_class"].isin(["Midline", "Hemispheric"])]
    ct_wt = pd.crosstab(wt["location_class"], wt["ecotype"]).reindex(index=["Midline", "Hemispheric"], columns=ECOTYPES)
    stat = contingency.set_index("analysis")
    fig, axes = plt.subplots(1, 4, figsize=(18, 5.6), gridspec_kw={"width_ratios": [1.2, 1.1, .85, 1.05]})
    def stacked(ax, ct, title, rowname, labels):
        pct = ct.div(ct.sum(axis=1), axis=0) * 100
        bottom = np.zeros(len(ct))
        for eco in ECOTYPES:
            vals = pct[eco].values
            ax.bar(np.arange(len(ct)), vals, bottom=bottom, color=COLORS[eco], edgecolor="white", label=eco)
            for i, (v, b, n) in enumerate(zip(vals, bottom, ct[eco])):
                if v >= 7: ax.text(i, b + v/2, str(int(n)), ha="center", va="center", fontsize=8, color="white", weight="bold")
            bottom += vals
        ax.set_xticks(np.arange(len(ct)), [f"{labels.get(x, x)}\n(n={int(ct.loc[x].sum())})" for x in ct.index], fontsize=8)
        ax.set_ylim(0, 100); ax.set_ylabel("Samples (%)")
        s = stat.loc[rowname]
        ax.set_title(title, weight="bold", fontsize=10)
        ax.text(.5, -.25, f"χ²={s.chi2:.2f}, df={int(s.df)}, P_MC={s.p_montecarlo:.4g}\nCramér's V={s.cramers_V:.2f}; n={int(s.n)}", transform=ax.transAxes, ha="center", va="top", fontsize=8)
    stacked(axes[0], ct_loc, "A  Ecotype by anatomical location", "Location (4 categories) x ecotype", {"Posterior_fossa": "Posterior\nfossa", "Ambiguous": "Multi-compartment/\nambiguous"})
    stacked(axes[1], ct_mol, "B  Ecotype by molecular group", "Molecular group x ecotype", {"DMG_K27": "DMG\nH3 K27-alt.", "DHG_G34": "DHG\nH3 G34-mut.", "pHGG_WT": "pHGG\nH3/IDH-wt", "IHG": "IHG"})
    stacked(axes[2], ct_wt, "C  Within H3/IDH-wt pHGG", "pHGG_WT: Midline vs Hemispheric x ecotype", {})
    key = perm[perm["term"].isin(["molecular_group", "location | molecular_group", "location", "molecular_group | location"])].copy()
    wanted = [("molecular_group + location", "molecular_group", "Molecular\nfirst", "#7C3AED"),
              ("molecular_group + location", "location | molecular_group", "Location\n| molecular", "#059669"),
              ("location + molecular_group", "location", "Location\nfirst", "#10B981"),
              ("location + molecular_group", "molecular_group | location", "Molecular\n| location", "#A78BFA")]
    vals=[]; labs=[]; cols=[]; pvals=[]
    for order, term, lab, col in wanted:
        r=key[(key["order"]==order)&(key["term"]==term)].iloc[0]
        vals.append(r.partial_R2); labs.append(lab); cols.append(col); pvals.append(r.p)
    axes[3].bar(np.arange(4), vals, color=cols)
    axes[3].set_xticks(np.arange(4), labs, fontsize=8)
    axes[3].set_ylabel("PERMANOVA partial R²")
    axes[3].set_title("D  Sequential PERMANOVA\n(34-feature TPM matrix, n=258)", weight="bold", fontsize=10)
    axes[3].set_ylim(0, max(vals)*1.42)
    for i,(v,p) in enumerate(zip(vals,pvals)): axes[3].text(i,v+.0015,f"{100*v:.1f}%\nP={p:.4g}",ha="center",fontsize=7.5)
    for ax in axes:
        ax.spines[["top","right"]].set_visible(False)
    h,l=axes[0].get_legend_handles_labels()
    fig.legend(h,l,loc="upper center",bbox_to_anchor=(.5,1.02),ncol=3,frameon=False)
    fig.tight_layout(rect=[0,.08,1,.95])
    save_fig(fig, figures / "Figure4")


def figure7c(main_cox: pd.DataFrame, ann: pd.DataFrame, figures: Path) -> None:
    pretty = {
        "ecotype_Lymphocyte-inflamed": "Lymphocyte-inflamed vs Immune-desert",
        "ecotype_Myeloid-dominant": "Myeloid-dominant vs Immune-desert",
        "cohort_group_DMG_K27": "DMG H3 K27-altered vs pHGG (ref)",
        "cohort_group_DHG_G34": "DHG H3 G34-mutant vs pHGG (ref)",
        "cohort_group_IHG": "Infant-type hemispheric vs pHGG (ref)",
        "age": "Age at diagnosis (per year)", "male": "Male vs female",
    }
    r = main_cox.set_index("term").loc[list(pretty)].reset_index()
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    y=np.arange(len(r))[::-1]
    for yy,(_,row) in zip(y,r.iterrows()):
        color="#C62828" if row.p<.05 else "#6B7280"
        ax.plot([row.CI_low,row.CI_high],[yy,yy],c=color,lw=2.3)
        ax.scatter(row.HR,yy,c=color,s=55,zorder=3)
        ax.text(23,yy,f"{row.HR:.2f} ({row.CI_low:.2f}–{row.CI_high:.2f}); P={row.p:.3g}",va="center",fontsize=8.5,family="monospace")
    ax.axvline(1,c="0.35",ls="--",lw=1)
    ax.set_xscale("log"); ax.set_xlim(.05,20)
    ax.set_yticks(y,[pretty[x] for x in r.term],fontsize=9)
    n=int(((pd.to_numeric(ann.OS_days,errors="coerce")>0)&pd.to_numeric(ann.age_years,errors="coerce").notna()).sum())
    events=int(ann.loc[pd.to_numeric(ann.OS_days,errors="coerce")>0,"OS_status"].astype(str).str.upper().eq("DECEASED").sum())
    ax.set_title(f"C  Multivariable Cox proportional-hazards model\nImmune-desert reference; n={n}, {events} events",weight="bold")
    ax.set_xlabel("Adjusted hazard ratio (95% CI), log scale", labelpad=10)
    ax.text(.12,.96,"Lower hazard",transform=ax.transAxes,color="#065F46",ha="center",va="top",fontsize=9)
    ax.text(.75,.96,"Higher hazard",transform=ax.transAxes,color="#7F1D1D",ha="center",va="top",fontsize=9)
    ax.spines[["top","right","left"]].set_visible(False); ax.tick_params(axis="y",length=0)
    fig.subplots_adjust(left=.34,right=.74,bottom=.16,top=.82)
    save_fig(fig, figures / "Figure7C")


def environment_report(root: Path, reports: Path, args: argparse.Namespace) -> None:
    import anndata, lifelines, matplotlib as mpl, sklearn, statsmodels
    versions = {
        "timestamp_local": pd.Timestamp.now().isoformat(), "platform": platform.platform(),
        "python_executable": sys.executable, "python_version": sys.version,
        "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__, "lifelines": lifelines.__version__,
        "matplotlib": mpl.__version__, "anndata": anndata.__version__,
        "statsmodels": statsmodels.__version__, "random_seed": SEED,
        "consensus_bootstrap": args.consensus_bootstrap,
        "permanova_permutations": args.permutations,
        "contingency_mc_permutations": args.mc_permutations,
        "root": str(root),
    }
    (reports / "execution_environment.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")
    lines = ["# Execution environment", "", *[f"- **{k}**: `{v}`" for k,v in versions.items()]]
    (reports / "execution_environment.md").write_text("\n".join(lines)+"\n", encoding="utf-8")


def write_report(metrics: dict, contingency: pd.DataFrame, perm: pd.DataFrame,
                 ipw: dict, eor: dict, deg: dict, reports: Path) -> None:
    loc = contingency.iloc[0]
    wt = contingency.iloc[2]
    molfirst = perm[(perm["order"]=="molecular_group + location") & (perm["term"]=="molecular_group")].iloc[0]
    locadj = perm[(perm["order"]=="molecular_group + location") & (perm["term"]=="location | molecular_group")].iloc[0]
    text = f"""# Final harmonized rerun report (2026-07-22)

## Fixed analysis set

- Canonical membership: **111 Lymphocyte-inflamed / 160 Myeloid-dominant / 78 Immune-desert**.
- TPM feature matrix: **34 features = 10 quanTIseq + 24 ssGSEA**, n=349.
- The previous 29-feature matrix was caused by positional truncation of the ssGSEA list and is superseded.

## 34-feature stability

- Consensus-clustering ARI versus the fixed assignment: **{metrics['ARI_vs_canonical']:.3f}**.
- NMI: **{metrics['NMI_vs_canonical']:.3f}**; exact sample agreement: **{100*metrics['per_sample_agreement']:.1f}%**.
- The rerun is a sensitivity analysis. The pre-specified 111/160/78 membership remains the sole sample annotation used downstream.

## Location and molecular group

- Four-category location association: chi-square={loc.chi2:.2f}, Monte-Carlo P={loc.p_montecarlo:.4g}, Cramer's V={loc.cramers_V:.3f}, n={int(loc.n)}.
- Within H3/IDH-wildtype pHGG: chi-square={wt.chi2:.2f}, Monte-Carlo P={wt.p_montecarlo:.4g}, n={int(wt.n)}.
- Sequential 34-feature PERMANOVA: molecular group first partial R2={molfirst.partial_R2:.4f}, P={molfirst.p:.4g}; location after molecular group partial R2={locadj.partial_R2:.4f}, P={locadj.p:.4g}.

## Missing survival data / IPW

- Observed OS: {ipw['n_observed']}/{ipw['n_total']}; missingness {ipw['missingness_percent']:.1f}%.
- |SMD| > 0.1 before/after weighting: {ipw['n_abs_SMD_gt_0.1_before']} / {ipw['n_abs_SMD_gt_0.1_after']}; maximum post-weighting |SMD|={ipw['max_abs_SMD_after']:.3f}.
- Interpretation is restricted to missingness conditional on measured covariates. **This is not an MNAR analysis and cannot establish robustness to MNAR.**

## Extent of resection

- Full OS model n={eor['n_full']}; complete EOR subset n={eor['n_EOR_subset']}.
- Likelihood-ratio comparison for adding EOR: chi-square={eor['likelihood_ratio_chi2']:.2f}, df=2, P={eor['likelihood_ratio_p']:.4g}.
- EOR results are complete-case sensitivity estimates and should not be presented as causal effects.

## Molecular-group-adjusted DEG

- Model: `{deg['model']}`; n={deg['n_samples']}, genes={deg['n_genes_tested']}.
- Pairwise adjusted contrasts are supplied as full TSV files. The original one-vs-rest DEG and pathway panels should be described as exploratory because they were not adjusted for molecular group.

## Claim calibration

- A non-significant ecotype-by-subtype interaction is reported as “no interaction detected,” not proof of generalizability.
- Subtype-specific estimates, especially DHG and IHG, are underpowered and cannot establish consistency.
- The ecotypes are exploratory transcriptomic states, not a validated clinical classifier or treatment-selection biomarker.
"""
    (reports / "FINAL_RERUN_REPORT.md").write_text(text, encoding="utf-8")


def manifest(out: Path) -> None:
    rows=[]
    for p in sorted(x for x in out.rglob("*") if x.is_file() and x.name != "SHA256_manifest.tsv"):
        h=hashlib.sha256()
        with p.open("rb") as f:
            for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
        rows.append({"file":str(p.relative_to(out)),"bytes":p.stat().st_size,"sha256":h.hexdigest()})
    pd.DataFrame(rows).to_csv(out/"SHA256_manifest.tsv",sep="\t",index=False)


def main() -> None:
    args=parse_args(); root=args.root.resolve()
    final=root/"FINAL MANUSCRIPT 260722"; out=final/"8. Final rerun outputs"
    tables=out/"01_tables"; figures=out/"02_figures"; reports=out/"03_reports"; h5dir=out/"04_reproducibility"; backup=out/"00_original_backup"
    for p in [tables,figures,reports,h5dir,backup]: p.mkdir(parents=True,exist_ok=True)
    ann=canonical_annotation(root,tables,backup)
    raw,z,assignment,metrics=build_34_features(root,ann,tables,h5dir,args.consensus_bootstrap)
    contingency,perm=contingency_and_permanova(root,ann,z,tables,args.permutations,args.mc_permutations)
    main_cox,ipw_table,balance,ipw_summary=main_cox_and_ipw(ann,tables)
    eor_table,eor_summary=eor_analysis(ann,tables)
    deg_table,deg_summary=adjusted_deg(root,ann,tables,figures)
    figure4(ann,contingency,perm,figures)
    figure7c(main_cox,ann,figures)
    environment_report(root,reports,args)
    write_report(metrics,contingency,perm,ipw_summary,eor_summary,deg_summary,reports)
    manifest(out)
    print(json.dumps({"status":"complete","output":str(out),"TPM34":metrics,"IPW":ipw_summary,"EOR":eor_summary,"DEG":deg_summary},indent=2))


if __name__ == "__main__":
    main()
