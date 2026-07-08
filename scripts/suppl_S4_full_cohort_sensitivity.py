"""
Suppl Fig S4 — Independent-primary-plus / Full-cohort sensitivity.

The available LM22+ssGSEA feature matrix covers the 702-biospecimen pool that was
sent to CIBERSORTx — composed of Main (349, independent-primary-plus) + BRAF_ALT
secondary (353). We therefore use this 702-pool de-novo re-clustering as the
practical full-cohort sensitivity check: if the ecotypes discovered on the Main
subset (n=349) are preserved when the BRAF_ALT samples are included in clustering
rather than projected, the partitioning is not an artifact of the cohort filter.

Pipeline:
  - Run consensus clustering (k=2-6, B=500, n_init=1) on the full 702 sample
    feature matrix
  - Annotate k=3 ecotypes by feature means
  - Crosstab vs Main n=349 ecotype assignment on shared samples
  - Crosstab vs Step 10 nearest-centroid projection of BRAF_ALT samples
  - Ecotype × cohort_group on the full 702 pool
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import json, warnings, time
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_suppl"; FIG.mkdir(exist_ok=True)
RNG  = np.random.default_rng(42)

feat = pd.read_csv(OUT / "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t", index_col=0)
X_df = feat.dropna(axis=0, how="any")
X    = X_df.values; n = X.shape[0]
print(f"Full pool feature matrix: n={n}")

ec_main = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
braf    = pd.read_csv(OUT / "step10_BRAF_ALT_ecotype_assignment.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
ann     = pd.read_csv(OUT / "cohort_for_deconvolution.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")

# Add cohort family annotation
def family(s):
    if s in ec_main.index:   return "Main_n349"
    if s in braf.index:
        return braf.loc[s, "set"]
    return "Other"
fam_series = pd.Series([family(s) for s in X_df.index], index=X_df.index, name="cohort_family")
print(fam_series.value_counts())

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

print("\nRunning consensus clustering on full pool k=2..6 ...")
rows_pac, rows_sil, rows_rob = [], [], []
clusters = {}
for k in range(2, 7):
    t0 = time.time()
    cm, lbl = consensus(X, k, B=500)
    upper = cm[np.triu_indices(n, k=1)]
    pac = float(((upper > 0.1) & (upper < 0.9)).mean())
    d = 1.0 - cm; np.fill_diagonal(d, 0.0)
    sil = float(silhouette_score(d, lbl, metric="precomputed"))
    same = np.equal.outer(lbl, lbl)
    rob = float(cm[same & ~np.eye(n, dtype=bool)].mean())
    rows_pac.append({"k": k, "PAC": pac}); rows_sil.append({"k": k, "silhouette": sil})
    rows_rob.append({"k": k, "within_cluster_consensus": rob})
    clusters[k] = (cm, lbl)
    print(f"  k={k}: PAC={pac:.3f}, sil={sil:.3f}, robustness={rob:.3f}  ({time.time()-t0:.1f}s)")

pd.DataFrame(rows_pac).to_csv(OUT / "supplS4_fullpool_PAC.tsv", sep="\t", index=False)
pd.DataFrame(rows_sil).to_csv(OUT / "supplS4_fullpool_silhouette.tsv", sep="\t", index=False)
pd.DataFrame(rows_rob).to_csv(OUT / "supplS4_fullpool_robustness.tsv", sep="\t", index=False)

# k=3 annotation
K = 3
lbl3 = clusters[K][1]
feat_c = X_df.copy(); feat_c["cluster"] = lbl3
means = feat_c.groupby("cluster").mean()
lym_feat = [c for c in ["ssGSEA_T_Cell_Cytotoxicity","ssGSEA_IFN_Gamma_Response","ssGSEA_MHC_Class_I",
                         "ssGSEA_Chemokine_T_Cell_Recruitment","LM22_T cells CD8"] if c in means.columns]
mye_feat = [c for c in ["ssGSEA_MoTAM_Antunes2021","ssGSEA_MgTAM_Antunes2021","ssGSEA_MDM_Klemm2020",
                         "ssGSEA_Microglia_Klemm2020","ssGSEA_M2_Macrophage","LM22_Macrophages M0",
                         "LM22_Macrophages M1","LM22_Macrophages M2","LM22_Monocytes"] if c in means.columns]
ls = means[lym_feat].mean(axis=1); ms = means[mye_feat].mean(axis=1)
desert = (ls + ms).idxmin()
rest = [c for c in means.index if c != desert]
lym = (ls[rest] - ms[rest]).idxmax(); mye = [c for c in rest if c != lym][0]
mapping = {int(lym):"Lymphocyte-inflamed", int(mye):"Myeloid-dominant", int(desert):"Immune-desert"}
print(f"mapping: {mapping}")

full_df = pd.DataFrame({"Kids_First_Biospecimen_ID": X_df.index,
                        "cluster": lbl3,
                        "ecotype_fullpool": [mapping.get(int(c)) for c in lbl3],
                        "cohort_family": fam_series.values}).set_index("Kids_First_Biospecimen_ID")
full_df.to_csv(OUT / "supplS4_fullpool_k3_assignment.tsv", sep="\t")
print("\nFull-pool k=3 sizes:", full_df["ecotype_fullpool"].value_counts().to_dict())

# Agreement vs Main k=3 (on Main samples)
both_main = full_df.index.intersection(ec_main.index)
ct_main = pd.crosstab(full_df.loc[both_main, "ecotype_fullpool"], ec_main.loc[both_main, "ecotype"])
ct_main.to_csv(OUT / "supplS4_fullpool_vs_Main_k3_crosstab.tsv", sep="\t")
agree_main = (full_df.loc[both_main, "ecotype_fullpool"].values == ec_main.loc[both_main, "ecotype"].values).mean()
print(f"\nAgreement on Main samples (n={len(both_main)}): {agree_main:.1%}")
print(ct_main.to_string())

# Agreement vs Step 10 BRAF nearest-centroid
both_braf = full_df.index.intersection(braf.index)
ct_braf = pd.crosstab(full_df.loc[both_braf, "ecotype_fullpool"], braf.loc[both_braf, "ecotype"])
ct_braf.to_csv(OUT / "supplS4_fullpool_vs_BRAF_projection_crosstab.tsv", sep="\t")
agree_braf = (full_df.loc[both_braf, "ecotype_fullpool"].values == braf.loc[both_braf, "ecotype"].values).mean()
print(f"\nAgreement vs BRAF nearest-centroid (n={len(both_braf)}): {agree_braf:.1%}")
print(ct_braf.to_string())

# Ecotype × cohort_family on full pool
full_df["cohort_group"] = ann.loc[full_df.index, "cohort_group"].values
ct_family = pd.crosstab(full_df["ecotype_fullpool"], full_df["cohort_family"])
ct_family_pct = ct_family.div(ct_family.sum(axis=0), axis=1) * 100
ct_family.to_csv(OUT / "supplS4_fullpool_ecotype_x_family.tsv", sep="\t")
print("\nFull-pool ecotype × cohort_family (%):")
print(ct_family_pct.round(1).to_string())

# 3-panel plot: PAC/silhouette ; vs Main heatmap ; vs BRAF heatmap
fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
ax = axes[0]
ax.plot([d["k"] for d in rows_pac], [d["PAC"] for d in rows_pac], "o-", color="#1f77b4", label="PAC")
ax2 = ax.twinx()
ax2.plot([d["k"] for d in rows_sil], [d["silhouette"] for d in rows_sil], "s--", color="#ff7f0e", label="silhouette")
ax.axvline(K, color="red", ls="--", alpha=.5)
ax.set_xlabel("k"); ax.set_ylabel("PAC", color="#1f77b4"); ax2.set_ylabel("silhouette", color="#ff7f0e")
ax.set_title(f"Full pool n={n}  PAC + silhouette")

ax = axes[1]
sns.heatmap(ct_main, annot=True, fmt="d", cmap="Blues", ax=ax, cbar=False)
ax.set_title(f"Full-pool vs Main k=3\n(agreement={agree_main:.1%}, n={len(both_main)})")
ax.set_xlabel("Main ecotype"); ax.set_ylabel("Full-pool re-cluster ecotype")

ax = axes[2]
sns.heatmap(ct_braf, annot=True, fmt="d", cmap="Greens", ax=ax, cbar=False)
ax.set_title(f"Full-pool vs BRAF projection\n(agreement={agree_braf:.1%}, n={len(both_braf)})")
ax.set_xlabel("BRAF (nearest-centroid) ecotype"); ax.set_ylabel("Full-pool re-cluster ecotype")
plt.tight_layout(); fig.savefig(FIG / "S4_fullpool_sensitivity.png", dpi=200); plt.close()

# Optional PCA overlay
pca = PCA(n_components=2).fit_transform(X)
palette = {"Lymphocyte-inflamed":"#3B82F6","Myeloid-dominant":"#EF4444","Immune-desert":"#9CA3AF"}
fig, ax = plt.subplots(figsize=(8, 6))
for ec, col in palette.items():
    mask = np.array(full_df["ecotype_fullpool"].values == ec)
    ax.scatter(pca[mask, 0], pca[mask, 1], c=col, label=ec, s=14, alpha=.7, edgecolor="none")
ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
ax.set_title(f"Full-pool de-novo clustering, k=3 (n={n})")
ax.legend(loc="best", fontsize=8, frameon=False)
plt.tight_layout(); fig.savefig(FIG / "S4_fullpool_PCA.png", dpi=200); plt.close()

summary = {
    "n_full_pool": int(n),
    "ecotype_sizes_fullpool": full_df["ecotype_fullpool"].value_counts().to_dict(),
    "agreement_vs_Main_k3":          float(agree_main),
    "agreement_vs_BRAF_centroid":    float(agree_braf),
    "PAC_k2_3":   [d["PAC"] for d in rows_pac[:2]],
    "silhouette_k2_3": [d["silhouette"] for d in rows_sil[:2]],
    "robustness_k2_3": [d["within_cluster_consensus"] for d in rows_rob[:2]],
}
(OUT / "supplS4_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
