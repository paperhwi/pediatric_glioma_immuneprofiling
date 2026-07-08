"""
Step 5 (re-run): Consensus clustering on Main cohort with LM22 + ssGSEA feature matrix.

Mirrors scripts/consensus_clustering.R but in Python, using:
  - bootstrap k-means consensus (B=500, pItem=0.8) — equivalent to ConsensusClusterPlus
  - PAC (Proportion of Ambiguous Clustering) for optimal-k selection
  - Silhouette from 1 - consensus matrix
  - Ecotype annotation by feature means
  - PCA + UMAP visualization
  - Fisher's exact (Monte Carlo) for ecotype × cohort / location / age_dev

Inputs
------
- output/step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv (sample × feature, z-scored)
- output/cohort_main_final.tsv     (Main cohort biospecimen IDs, n=349)
- output/cohort_for_deconvolution.tsv (annotation table)

Outputs (output/)
-----------------
- consensus_LM22_main_k{k}.npz        (per-k consensus matrix + cluster labels)
- consensus_LM22_main_PAC.tsv
- consensus_LM22_main_silhouette.tsv
- ecotype_LM22_main_k{k}.tsv
- ecotype_LM22_main_k3_annotated.tsv
- ecotype_LM22_main_k3_feature_means.tsv
- ecotype_LM22_main_pca.tsv / umap.tsv
- ecotype_LM22_main_assoc_pvalues.tsv
- ecotype_LM22_main_x_{cohort_group,location,age_dev}.tsv
- figs_step5/*.png
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import fisher_exact
import warnings
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"
FIG  = OUT / "figs_step5"
FIG.mkdir(exist_ok=True)
RNG  = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# 1. Load feature matrix (sample × feature, z-scored) + cohort + annotation
# ---------------------------------------------------------------------------
print("=== [1] Load ===")
feat = pd.read_csv(OUT / "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t", index_col=0)
feat.index.name = "sample"
print(f"  feature matrix: {feat.shape}")

main_ids = pd.read_csv(OUT / "cohort_main_final.tsv", sep="\t")["Kids_First_Biospecimen_ID"].tolist()
main_ids = [s for s in main_ids if s in feat.index]
print(f"  main cohort available in feature matrix: {len(main_ids)}")

ann_full = pd.read_csv(OUT / "cohort_for_deconvolution.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
print(f"  annotation: {ann_full.shape}")

# Restrict feature matrix to main cohort, drop rows with any NA
X_df = feat.loc[main_ids].dropna(axis=0, how="any")
X = X_df.values
print(f"  X for clustering: {X.shape}  (samples × features)")

# ---------------------------------------------------------------------------
# 2. Consensus clustering (k-means bootstrap)
# ---------------------------------------------------------------------------
print("\n=== [2] Consensus clustering (k=2..6, B=500, pItem=0.8) ===")
B = 500
pItem = 0.8
ks = range(2, 7)
n = X.shape[0]

consensus = {}                # k -> (consensus matrix, cluster labels)
for k in ks:
    M = np.zeros((n, n))      # co-cluster count
    C = np.zeros((n, n))      # co-sampling count
    for b in range(B):
        idx = RNG.choice(n, int(n * pItem), replace=False)
        km = KMeans(n_clusters=k, n_init=10, random_state=42 + b)
        lbl = km.fit_predict(X[idx])
        for ci in range(k):
            members = idx[lbl == ci]
            if len(members) < 2:
                continue
            xs, ys = np.meshgrid(members, members)
            M[xs, ys] += 1
        # co-sampling count
        xs, ys = np.meshgrid(idx, idx)
        C[xs, ys] += 1

    with np.errstate(invalid="ignore", divide="ignore"):
        cm = np.where(C > 0, M / C, 0.0)
    np.fill_diagonal(cm, 1.0)

    # Final cluster: hierarchical (ward) on (1 - consensus) — mirror R innerLinkage/finalLinkage
    dist = 1.0 - cm
    np.fill_diagonal(dist, 0.0)
    cond = squareform(dist, checks=False)
    Z = linkage(cond, method="average")  # consensusClusterPlus default consensus tree
    lbl = fcluster(Z, t=k, criterion="maxclust")

    consensus[k] = (cm, lbl)
    np.savez_compressed(OUT / f"consensus_LM22_main_k{k}.npz", consensus=cm, labels=lbl,
                         samples=np.array(X_df.index, dtype=object))
    print(f"  k={k}: cluster sizes = {dict(Counter(lbl))}")

# ---------------------------------------------------------------------------
# 3. PAC + silhouette  → optimal k
# ---------------------------------------------------------------------------
print("\n=== [3] PAC & silhouette ===")
pac_rows, sil_rows = [], []
for k in ks:
    cm, lbl = consensus[k]
    upper = cm[np.triu_indices(n, k=1)]
    pac = float(((upper > 0.1) & (upper < 0.9)).mean())
    dist = 1.0 - cm
    np.fill_diagonal(dist, 0.0)
    sil = float(silhouette_score(dist, lbl, metric="precomputed"))
    pac_rows.append({"k": k, "PAC": pac})
    sil_rows.append({"k": k, "silhouette": sil})

pac_df = pd.DataFrame(pac_rows)
sil_df = pd.DataFrame(sil_rows)
print(pac_df.to_string(index=False))
print(sil_df.to_string(index=False))
pac_df.to_csv(OUT / "consensus_LM22_main_PAC.tsv", sep="\t", index=False)
sil_df.to_csv(OUT / "consensus_LM22_main_silhouette.tsv", sep="\t", index=False)

best_pac_k = int(pac_df.loc[pac_df["PAC"].idxmin(), "k"])
best_sil_k = int(sil_df.loc[sil_df["silhouette"].idxmax(), "k"])
print(f"  best PAC at k={best_pac_k} ; best silhouette at k={best_sil_k}")
print(f"  pre-specified hypothesis: k=3 (Lymph / Myeloid / Desert)")
K = 3

# ---------------------------------------------------------------------------
# 4. PAC / silhouette plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
ax[0].plot(pac_df["k"], pac_df["PAC"], "o-", color="#1f77b4")
ax[0].set_xlabel("k"); ax[0].set_ylabel("PAC (lower = better)"); ax[0].set_title("Proportion of Ambiguous Clustering")
ax[0].axvline(K, color="red", linestyle="--", alpha=0.5, label=f"k={K} (pre-specified)")
ax[0].legend()
ax[1].plot(sil_df["k"], sil_df["silhouette"], "o-", color="#ff7f0e")
ax[1].set_xlabel("k"); ax[1].set_ylabel("Silhouette (higher = better)"); ax[1].set_title("Silhouette on 1−consensus")
ax[1].axvline(K, color="red", linestyle="--", alpha=0.5)
plt.tight_layout()
fig.savefig(FIG / "PAC_silhouette_vs_k.png", dpi=200)
plt.close()

# ---------------------------------------------------------------------------
# 5. Save ecotype assignment for k=2,3,4 (sensitivity)
# ---------------------------------------------------------------------------
print("\n=== [4] Save ecotype assignments k=2,3,4 ===")
for k in (2, 3, 4):
    lbl = consensus[k][1]
    out_df = pd.DataFrame({"Kids_First_Biospecimen_ID": X_df.index, "cluster_k": [f"E{c}" for c in lbl]})
    out_df.to_csv(OUT / f"ecotype_LM22_main_k{k}.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# 6. Annotate k=3 ecotypes by dominant features
# ---------------------------------------------------------------------------
print("\n=== [5] Annotate k=3 ecotypes ===")
lbl3 = consensus[K][1]
clus_df = pd.DataFrame({"Kids_First_Biospecimen_ID": X_df.index, "cluster": lbl3}).set_index("Kids_First_Biospecimen_ID")
feat_clu = X_df.copy()
feat_clu["cluster"] = lbl3
means = feat_clu.groupby("cluster").mean()
means.to_csv(OUT / "ecotype_LM22_main_k3_feature_means.tsv", sep="\t")

def top_cluster(col): return int(means[col].idxmax()) if col in means.columns else None
def bot_cluster(col): return int(means[col].idxmin()) if col in means.columns else None

lym = top_cluster("ssGSEA_T_Cell_Cytotoxicity") or top_cluster("LM22_T cells CD8")
mye = top_cluster("ssGSEA_MoTAM_Antunes2021")  or top_cluster("LM22_Macrophages M2")
all_clu = set(int(c) for c in means.index.tolist())
desert = list(all_clu - {lym, mye})
desert = int(desert[0]) if desert else None

mapping = {lym: "Lymphocyte-inflamed", mye: "Myeloid-dominant"}
if desert is not None:
    mapping[desert] = "Immune-desert"
print(f"  cluster -> ecotype mapping: {mapping}")
clus_df["ecotype"] = clus_df["cluster"].map(mapping)
clus_df.to_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t")

# ---------------------------------------------------------------------------
# 7. PCA + UMAP
# ---------------------------------------------------------------------------
print("\n=== [6] PCA + UMAP ===")
pca = PCA(n_components=2).fit_transform(X)
pca_df = pd.DataFrame({"sample": X_df.index, "PC1": pca[:, 0], "PC2": pca[:, 1],
                       "ecotype": clus_df["ecotype"].values})
pca_df.to_csv(OUT / "ecotype_LM22_main_pca.tsv", sep="\t", index=False)

try:
    import umap
    um = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=42).fit_transform(X)
    umap_df = pd.DataFrame({"sample": X_df.index, "UMAP1": um[:, 0], "UMAP2": um[:, 1],
                            "ecotype": clus_df["ecotype"].values})
    umap_df.to_csv(OUT / "ecotype_LM22_main_umap.tsv", sep="\t", index=False)
except Exception as e:
    print(f"  umap-learn not available ({e}); skipping UMAP")
    umap_df = None

ec_colors = {"Lymphocyte-inflamed": "#3B82F6", "Myeloid-dominant": "#EF4444", "Immune-desert": "#9CA3AF"}
fig, axes = plt.subplots(1, 2 if umap_df is not None else 1, figsize=(11, 4.2))
ax0 = axes[0] if umap_df is not None else axes
for ec, sub in pca_df.groupby("ecotype"):
    ax0.scatter(sub["PC1"], sub["PC2"], c=ec_colors.get(ec, "#888"), label=ec, s=18, alpha=0.7, edgecolor="none")
ax0.set_xlabel("PC1"); ax0.set_ylabel("PC2"); ax0.set_title("PCA")
ax0.legend(loc="best", frameon=False, fontsize=8)
if umap_df is not None:
    ax1 = axes[1]
    for ec, sub in umap_df.groupby("ecotype"):
        ax1.scatter(sub["UMAP1"], sub["UMAP2"], c=ec_colors.get(ec, "#888"), label=ec, s=18, alpha=0.7, edgecolor="none")
    ax1.set_xlabel("UMAP1"); ax1.set_ylabel("UMAP2"); ax1.set_title("UMAP")
plt.tight_layout()
fig.savefig(FIG / "ecotype_PCA_UMAP.png", dpi=200)
plt.close()

# ---------------------------------------------------------------------------
# 8. Association tests (ecotype × cohort / location / age_dev)
# ---------------------------------------------------------------------------
print("\n=== [7] Association tests ===")
meta = ann_full.loc[clus_df.index.intersection(ann_full.index)].copy()
meta["ecotype"] = clus_df["ecotype"]

def fisher_mc(tab, B=10000):
    # simulated Fisher's exact via permutation of group labels
    obs = tab.values
    rows = obs.sum(axis=1); cols = obs.sum(axis=0); N = obs.sum()
    if N == 0:
        return float("nan")
    from scipy.stats import chi2_contingency, fisher_exact
    if obs.shape == (2, 2):
        return float(fisher_exact(obs)[1])
    chi2, _, _, _ = chi2_contingency(obs, correction=False)
    null = []
    flat_groups = []
    for r in range(obs.shape[0]):
        flat_groups.extend([r] * int(rows[r]))
    col_assign = []
    for c in range(obs.shape[1]):
        col_assign.extend([c] * int(cols[c]))
    flat_groups = np.array(flat_groups); col_assign = np.array(col_assign)
    if len(flat_groups) != len(col_assign):
        return float("nan")
    for _ in range(B):
        permuted = RNG.permutation(col_assign)
        ct = np.zeros_like(obs)
        for g, c in zip(flat_groups, permuted):
            ct[g, c] += 1
        chi2_p, _, _, _ = chi2_contingency(ct, correction=False)
        null.append(chi2_p)
    return float((1 + (np.array(null) >= chi2).sum()) / (B + 1))

assoc_rows = []
for var in ["cohort_group", "location_class", "age_dev_group"]:
    if var == "location_class":
        sub = meta[meta["location_class"].isin(["Midline", "Hemispheric", "Posterior_fossa"])]
    else:
        sub = meta
    tab = pd.crosstab(sub["ecotype"], sub[var])
    tab.to_csv(OUT / f"ecotype_LM22_main_x_{var}.tsv", sep="\t")
    print(f"\n  ecotype × {var}:")
    print(tab.to_string())
    p = fisher_mc(tab, B=10000)
    print(f"  Fisher MC p = {p:.2e}")
    assoc_rows.append({"comparison": f"ecotype_x_{var}", "fisher_mc_p": p, "N": int(tab.values.sum())})

pd.DataFrame(assoc_rows).to_csv(OUT / "ecotype_LM22_main_assoc_pvalues.tsv", sep="\t", index=False)

# ---------------------------------------------------------------------------
# 9. Compare with previous quanTIseq-only ecotype assignment (Step 5 baseline)
# ---------------------------------------------------------------------------
print("\n=== [8] Comparison with previous (quanTIseq+ssGSEA) k=3 ecotype assignment ===")
try:
    prev = pd.read_csv(OUT / "ecotype_assignment_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
    common = clus_df.index.intersection(prev.index)
    if len(common) > 0:
        comp = pd.DataFrame({
            "LM22+ssGSEA": clus_df.loc[common, "ecotype"],
            "previous(quanTIseq+ssGSEA)": prev.loc[common, "ecotype"],
        })
        ct = pd.crosstab(comp["LM22+ssGSEA"], comp["previous(quanTIseq+ssGSEA)"])
        ct.to_csv(OUT / "ecotype_LM22_vs_previous_crosstab.tsv", sep="\t")
        print(ct.to_string())
        # agreement
        agree = (comp["LM22+ssGSEA"] == comp["previous(quanTIseq+ssGSEA)"]).mean()
        print(f"  exact label agreement on n={len(common)}: {agree:.1%}")
except Exception as e:
    print(f"  comparison skipped: {e}")

# ---------------------------------------------------------------------------
# 10. Summary JSON
# ---------------------------------------------------------------------------
summary = {
    "feature_matrix":  "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv",
    "n_main_used":     int(X.shape[0]),
    "n_features":      int(X.shape[1]),
    "bootstrap_B":     B,
    "pItem":           pItem,
    "best_PAC_k":      int(best_pac_k),
    "best_silhouette_k": int(best_sil_k),
    "primary_K":       K,
    "ecotype_sizes":   {k: int(v) for k, v in Counter(clus_df["ecotype"].tolist()).items()},
    "outputs": [
        "consensus_LM22_main_PAC.tsv",
        "consensus_LM22_main_silhouette.tsv",
        "ecotype_LM22_main_k3_annotated.tsv",
        "ecotype_LM22_main_k3_feature_means.tsv",
        "ecotype_LM22_main_pca.tsv",
        "ecotype_LM22_main_umap.tsv",
        "ecotype_LM22_main_assoc_pvalues.tsv",
        "ecotype_LM22_vs_previous_crosstab.tsv",
        "figs_step5/PAC_silhouette_vs_k.png",
        "figs_step5/ecotype_PCA_UMAP.png",
    ],
}
(OUT / "step5_summary.json").write_text(json.dumps(summary, indent=2))
print("\n=== SUMMARY ===")
print(json.dumps(summary, indent=2))
