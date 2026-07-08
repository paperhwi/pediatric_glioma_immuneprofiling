"""
Suppl Fig S3 — Re-clustering on the high-confidence (LM22 P<0.05) subset.

Demonstrates that the three-ecotype structure is recapitulated when restricted to
the conservative CIBERSORTx-confidence subset, i.e. the ecotypes are not artifacts
of low-quality deconvolution.

Pipeline:
  - Filter to samples with LM22 P-value < 0.05 within Main cohort
  - Same LM22+ssGSEA z-scored feature matrix
  - Consensus clustering (k=2-6, B=500, n_init=1, KMeans + Ward consensus tree)
  - Annotate k=3 ecotypes; compute agreement vs Main k=3 assignment
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import json, warnings, time
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_suppl"; FIG.mkdir(exist_ok=True)
RNG  = np.random.default_rng(42)

# Inputs
feat  = pd.read_csv(OUT / "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t", index_col=0)
qc    = pd.read_csv(OUT / "CIBERSORTx_Job15_Results.txt", sep="\t").rename(columns={"Mixture":"sample"}).set_index("sample")
main  = pd.read_csv(OUT / "cohort_main_final.tsv", sep="\t")["Kids_First_Biospecimen_ID"].tolist()
ec3   = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
ann   = pd.read_csv(OUT / "cohort_for_deconvolution.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")

# Restrict to Main ∩ P<0.05
hc_mask = qc["P-value"] < 0.05
hc_ids  = set(qc[hc_mask].index)
use_ids = [s for s in main if s in feat.index and s in hc_ids]
X_df = feat.loc[use_ids].dropna(axis=0, how="any")
X    = X_df.values; n = X.shape[0]
print(f"High-confidence (P<0.05) subset within Main: n={n}")
print("cohort_group distribution:")
print(ann.loc[use_ids, "cohort_group"].value_counts())

# Consensus clustering
def consensus(X, k, B=500, pItem=0.8, n_init=1, seed=42):
    n = X.shape[0]; M = np.zeros((n, n)); C = np.zeros((n, n))
    rng = np.random.default_rng(seed)
    for b in range(B):
        idx = rng.choice(n, max(2, int(n * pItem)), replace=False)
        km = KMeans(n_clusters=k, n_init=n_init, random_state=seed+b)
        lbl = km.fit_predict(X[idx])
        ixg = np.ix_(idx, idx); C[ixg] += 1
        for ci in range(k):
            mem = idx[lbl == ci]
            if len(mem) >= 2:
                M[np.ix_(mem, mem)] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        cm = np.where(C > 0, M / C, 0.0)
    np.fill_diagonal(cm, 1.0)
    dist = 1.0 - cm; np.fill_diagonal(dist, 0.0)
    Z = linkage(squareform(dist, checks=False), method="average")
    return cm, fcluster(Z, t=k, criterion="maxclust")

print("\nRunning consensus clustering k=2..6 ...")
pacs, sils, rob = [], [], []
clusters = {}
for k in range(2, 7):
    t0 = time.time()
    cm, lbl = consensus(X, k, B=500)
    upper = cm[np.triu_indices(n, k=1)]
    pac = float(((upper > 0.1) & (upper < 0.9)).mean())
    d = 1.0 - cm; np.fill_diagonal(d, 0.0)
    sil = float(silhouette_score(d, lbl, metric="precomputed"))
    same = np.equal.outer(lbl, lbl)
    rob_v = float(cm[same & ~np.eye(n, dtype=bool)].mean())
    pacs.append({"k": k, "PAC": pac}); sils.append({"k": k, "silhouette": sil})
    rob.append({"k": k, "within_cluster_consensus": rob_v})
    clusters[k] = (cm, lbl)
    print(f"  k={k}: PAC={pac:.3f}, sil={sil:.3f}, robustness={rob_v:.3f}  ({time.time()-t0:.1f}s)")

pd.DataFrame(pacs).to_csv(OUT / "supplS3_highconf_PAC.tsv", sep="\t", index=False)
pd.DataFrame(sils).to_csv(OUT / "supplS3_highconf_silhouette.tsv", sep="\t", index=False)
pd.DataFrame(rob).to_csv(OUT / "supplS3_highconf_robustness.tsv", sep="\t", index=False)

# k=3 annotation reusing Step 5c logic
K = 3
lbl3 = clusters[K][1]
feat_c = X_df.copy(); feat_c["cluster"] = lbl3
means  = feat_c.groupby("cluster").mean()
lym_feat = [c for c in ["ssGSEA_T_Cell_Cytotoxicity","ssGSEA_IFN_Gamma_Response","ssGSEA_MHC_Class_I",
                         "ssGSEA_Chemokine_T_Cell_Recruitment","LM22_T cells CD8"] if c in means.columns]
mye_feat = [c for c in ["ssGSEA_MoTAM_Antunes2021","ssGSEA_MgTAM_Antunes2021","ssGSEA_MDM_Klemm2020",
                         "ssGSEA_Microglia_Klemm2020","ssGSEA_M2_Macrophage","LM22_Macrophages M0",
                         "LM22_Macrophages M1","LM22_Macrophages M2","LM22_Monocytes"] if c in means.columns]
ls = means[lym_feat].mean(axis=1); ms = means[mye_feat].mean(axis=1)
desert = (ls + ms).idxmin()
rest = [c for c in means.index if c != desert]
lym = (ls[rest] - ms[rest]).idxmax()
mye = [c for c in rest if c != lym][0]
mapping = {int(lym):"Lymphocyte-inflamed", int(mye):"Myeloid-dominant", int(desert):"Immune-desert"}
print(f"  mapping (high-conf k=3): {mapping}")

hc_df = pd.DataFrame({"Kids_First_Biospecimen_ID": X_df.index, "cluster": lbl3}).set_index("Kids_First_Biospecimen_ID")
hc_df["ecotype_highconf"] = hc_df["cluster"].map(lambda c: mapping.get(int(c)))
hc_df.to_csv(OUT / "supplS3_highconf_k3_assignment.tsv", sep="\t")
print("\nHigh-conf k=3 sizes:", hc_df["ecotype_highconf"].value_counts().to_dict())

# Agreement vs full Main k=3
both = hc_df.index.intersection(ec3.index)
ct = pd.crosstab(hc_df.loc[both, "ecotype_highconf"], ec3.loc[both, "ecotype"])
ct.to_csv(OUT / "supplS3_highconf_vs_Main_k3_crosstab.tsv", sep="\t")
agree = (hc_df.loc[both, "ecotype_highconf"].values == ec3.loc[both, "ecotype"].values).mean()
print(f"\nAgreement with Main k=3 on shared n={len(both)}: {agree:.1%}")
print(ct.to_string())

# Plot — PAC/silhouette + crosstab heatmap
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
ax = axes[0]
ax.plot([d["k"] for d in pacs], [d["PAC"] for d in pacs], "o-", color="#1f77b4", label="PAC")
ax.axvline(K, color="red", ls="--", alpha=.5); ax.set_xlabel("k"); ax.set_ylabel("PAC")
ax.set_title(f"PAC (high-conf subset, n={n})")
ax = axes[1]
ax.plot([d["k"] for d in sils], [d["silhouette"] for d in sils], "o-", color="#ff7f0e")
ax.axvline(K, color="red", ls="--", alpha=.5); ax.set_xlabel("k"); ax.set_ylabel("silhouette")
ax.set_title("Silhouette")
ax = axes[2]
sns.heatmap(ct, annot=True, fmt="d", cmap="Blues", ax=ax, cbar=False)
ax.set_title(f"High-conf vs Main k=3  (agreement={agree:.1%})")
ax.set_xlabel("Main ecotype"); ax.set_ylabel("High-conf re-cluster ecotype")
plt.tight_layout(); fig.savefig(FIG / "S3_highconf_subset_reclustering.png", dpi=200); plt.close()

# Ecotype × cohort association in high-conf subset (sanity)
hc_df["cohort_group"]   = ann.loc[hc_df.index, "cohort_group"]
hc_df["location_class"] = ann.loc[hc_df.index, "location_class"]
hc_ct = pd.crosstab(hc_df["ecotype_highconf"], hc_df["cohort_group"])
hc_ct.to_csv(OUT / "supplS3_highconf_ecotype_x_cohort.tsv", sep="\t")
print("\nHigh-conf ecotype × cohort_group:")
print(hc_ct.to_string())

summary = {
    "n_highconf_in_Main": int(n),
    "P_thresh": 0.05,
    "ecotype_sizes_highconf": hc_df["ecotype_highconf"].value_counts().to_dict(),
    "agreement_vs_Main_k3": float(agree),
    "PAC_k2_3": [d["PAC"] for d in pacs[:2]],
    "silhouette_k2_3": [d["silhouette"] for d in sils[:2]],
    "robustness_k2_3": [d["within_cluster_consensus"] for d in rob[:2]],
}
(OUT / "supplS3_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
