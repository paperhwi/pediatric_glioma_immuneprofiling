"""Fix Step 5 ecotype annotation (handle lym==mye collision) and produce k=2 + k=3 outputs."""
from pathlib import Path
from collections import Counter
import json, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"

feat = pd.read_csv(OUT / "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t", index_col=0)
main_ids = [s for s in pd.read_csv(OUT / "cohort_main_final.tsv", sep="\t")["Kids_First_Biospecimen_ID"] if s in feat.index]
X_df = feat.loc[main_ids].dropna(axis=0, how="any")

ann_full = pd.read_csv(OUT / "cohort_for_deconvolution.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")

def feature_means(K):
    z = np.load(OUT / f"consensus_LM22_main_k{K}.npz", allow_pickle=True)
    lbl = z["labels"]
    df = X_df.copy(); df["cluster"] = lbl
    means = df.groupby("cluster").mean()
    return lbl, means

def annotate_three(means):
    """Annotate Lymphocyte / Myeloid / Desert using composite scores.
    Robust to two cell-types peaking in same cluster."""
    lym_feat = [c for c in [
        "ssGSEA_T_Cell_Cytotoxicity",
        "ssGSEA_IFN_Gamma_Response",
        "ssGSEA_MHC_Class_I",
        "ssGSEA_Chemokine_T_Cell_Recruitment",
        "LM22_T cells CD8",
    ] if c in means.columns]
    mye_feat = [c for c in [
        "ssGSEA_MoTAM_Antunes2021",
        "ssGSEA_MgTAM_Antunes2021",
        "ssGSEA_MDM_Klemm2020",
        "ssGSEA_Microglia_Klemm2020",
        "ssGSEA_M2_Macrophage",
        "LM22_Macrophages M0", "LM22_Macrophages M1", "LM22_Macrophages M2",
        "LM22_Monocytes",
    ] if c in means.columns]
    # composite z-score-of-z-score per cluster
    lym_score = means[lym_feat].mean(axis=1)
    mye_score = means[mye_feat].mean(axis=1)
    # Immune-desert = lowest combined immune signal
    all_clu = means.index.tolist()
    if len(all_clu) >= 3:
        # Desert: cluster with minimum (lym_score + mye_score)
        desert = (lym_score + mye_score).idxmin()
        remaining = [c for c in all_clu if c != desert]
        # Among remaining, lymph is the one with higher lym_score relative to mye_score
        rel = lym_score[remaining] - mye_score[remaining]
        lym = rel.idxmax()
        mye = [c for c in remaining if c != lym][0]
    elif len(all_clu) == 2:
        # k=2: assign by composite ratio
        rel = lym_score - mye_score
        lym = rel.idxmax()
        mye = rel.idxmin()
        desert = None
    else:
        lym = mye = desert = None
    return lym, mye, desert, lym_score, mye_score

results = {}
for K in (2, 3):
    lbl, means = feature_means(K)
    means.to_csv(OUT / f"ecotype_LM22_main_k{K}_feature_means.tsv", sep="\t")
    lym, mye, desert, lym_score, mye_score = annotate_three(means)
    if K == 3:
        ecotype_map = {int(lym): "Lymphocyte-inflamed", int(mye): "Myeloid-dominant", int(desert): "Immune-desert"}
    else:
        # k=2: call it immune-active vs immune-desert if total score gap is large; otherwise lymph vs myeloid
        all_clu = means.index.tolist()
        c1, c2 = all_clu
        total_diff = abs((lym_score[c1] + mye_score[c1]) - (lym_score[c2] + mye_score[c2]))
        lym_diff   = abs(lym_score[c1] - lym_score[c2])
        mye_diff   = abs(mye_score[c1] - mye_score[c2])
        if total_diff > max(lym_diff, mye_diff):
            # active vs desert axis
            active = (lym_score + mye_score).idxmax()
            desert = (lym_score + mye_score).idxmin()
            ecotype_map = {int(active): "Immune-active", int(desert): "Immune-desert"}
        else:
            ecotype_map = {int(lym): "Lymphocyte-inflamed", int(mye): "Myeloid-dominant"}
    print(f"\n=== k={K} ===")
    print(f"  per-cluster lym_score: {lym_score.to_dict()}")
    print(f"  per-cluster mye_score: {mye_score.to_dict()}")
    print(f"  mapping: {ecotype_map}")

    clus_df = pd.DataFrame({"Kids_First_Biospecimen_ID": X_df.index, "cluster": lbl}).set_index("Kids_First_Biospecimen_ID")
    clus_df["ecotype"] = clus_df["cluster"].map(lambda c: ecotype_map.get(int(c)))
    if clus_df["ecotype"].isna().any():
        print("  WARN: some samples still NaN — falling back to numeric cluster label")
        clus_df["ecotype"] = clus_df["ecotype"].fillna(clus_df["cluster"].astype(str).map(lambda c: f"E{c}"))
    clus_df.to_csv(OUT / f"ecotype_LM22_main_k{K}_annotated.tsv", sep="\t")
    results[K] = {"sizes": dict(Counter(clus_df["ecotype"])), "map": ecotype_map}

    # Association tests
    meta = ann_full.loc[clus_df.index.intersection(ann_full.index)].copy()
    meta["ecotype"] = clus_df["ecotype"]
    for var in ["cohort_group", "location_class", "age_dev_group"]:
        sub = meta[meta["location_class"].isin(["Midline", "Hemispheric", "Posterior_fossa"])] if var == "location_class" else meta
        tab = pd.crosstab(sub["ecotype"], sub[var])
        tab.to_csv(OUT / f"ecotype_LM22_main_k{K}_x_{var}.tsv", sep="\t")

print("\n=== SUMMARY ===")
print(json.dumps(results, indent=2, default=str))
(OUT / "step5_annotation_summary.json").write_text(json.dumps(results, indent=2, default=str))

# Compare with previous (R) ecotype_assignment_k3_annotated.tsv
try:
    prev_raw = pd.read_csv(OUT / "ecotype_assignment_k3_annotated.tsv", sep="\t")
    id_col = next((c for c in prev_raw.columns if "Biospecimen" in c), prev_raw.columns[0])
    prev = prev_raw.set_index(id_col)
    new_k3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
    common = new_k3.index.intersection(prev.index)
    if "ecotype" in prev.columns and len(common) > 0:
        ct = pd.crosstab(new_k3.loc[common, "ecotype"], prev.loc[common, "ecotype"])
        ct.to_csv(OUT / "ecotype_LM22_vs_previous_crosstab.tsv", sep="\t")
        print("\nLM22+ssGSEA (rows) vs previous quanTIseq+ssGSEA (cols):")
        print(ct.to_string())
        agree = (new_k3.loc[common, "ecotype"].values == prev.loc[common, "ecotype"].values).mean()
        print(f"exact label agreement: {agree:.1%}  (n={len(common)})")
    else:
        print(f"prev columns: {prev.columns.tolist()}, skip compare")
except Exception as e:
    print(f"comparison failed: {e}")
