"""Resume Step 5: finish k=6 + PAC/silhouette/annotation/PCA/UMAP/assoc."""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import json, warnings
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency, fisher_exact
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_step5"; FIG.mkdir(exist_ok=True)
RNG  = np.random.default_rng(42)

feat = pd.read_csv(OUT / "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t", index_col=0)
main_ids = [s for s in pd.read_csv(OUT / "cohort_main_final.tsv", sep="\t")["Kids_First_Biospecimen_ID"] if s in feat.index]
X_df = feat.loc[main_ids].dropna(axis=0, how="any"); X = X_df.values
n = X.shape[0]
print(f"X: {X.shape}")

# Run k=6 with vectorized consensus (B=200, n_init=4)
def consensus_one_k(X, k, B=200, pItem=0.8, n_init=4, seed=42):
    n = X.shape[0]
    M = np.zeros((n, n)); C = np.zeros((n, n))
    rng = np.random.default_rng(seed)
    for b in range(B):
        idx = rng.choice(n, int(n * pItem), replace=False)
        km = KMeans(n_clusters=k, n_init=n_init, random_state=seed+b)
        lbl = km.fit_predict(X[idx])
        ixgrid = np.ix_(idx, idx)
        C[ixgrid] += 1
        for ci in range(k):
            mem = idx[lbl == ci]
            if len(mem) >= 2:
                M[np.ix_(mem, mem)] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        cm = np.where(C > 0, M / C, 0.0)
    np.fill_diagonal(cm, 1.0)
    dist = 1.0 - cm; np.fill_diagonal(dist, 0.0)
    Z = linkage(squareform(dist, checks=False), method="average")
    lbl_final = fcluster(Z, t=k, criterion="maxclust")
    return cm, lbl_final

# Load existing k=2-5
consensus = {}
for k in [2, 3, 4, 5]:
    z = np.load(OUT / f"consensus_LM22_main_k{k}.npz", allow_pickle=True)
    consensus[k] = (z["consensus"], z["labels"])
    print(f"loaded k={k}")

# k=6 with faster B=200
cm6, lbl6 = consensus_one_k(X, 6, B=200, n_init=4)
np.savez_compressed(OUT / "consensus_LM22_main_k6.npz", consensus=cm6, labels=lbl6,
                   samples=np.array(X_df.index, dtype=object))
consensus[6] = (cm6, lbl6)
print(f"k=6: {dict(Counter(lbl6))}")

# PAC + silhouette
ks = [2, 3, 4, 5, 6]
pac_rows = []; sil_rows = []
for k in ks:
    cm, lbl = consensus[k]
    upper = cm[np.triu_indices(n, k=1)]
    pac = float(((upper > 0.1) & (upper < 0.9)).mean())
    dist = 1.0 - cm; np.fill_diagonal(dist, 0.0)
    sil = float(silhouette_score(dist, lbl, metric="precomputed"))
    pac_rows.append({"k": k, "PAC": pac}); sil_rows.append({"k": k, "silhouette": sil})
pac_df = pd.DataFrame(pac_rows); sil_df = pd.DataFrame(sil_rows)
print(pac_df.to_string(index=False)); print(sil_df.to_string(index=False))
pac_df.to_csv(OUT / "consensus_LM22_main_PAC.tsv", sep="\t", index=False)
sil_df.to_csv(OUT / "consensus_LM22_main_silhouette.tsv", sep="\t", index=False)
best_pac_k = int(pac_df.loc[pac_df["PAC"].idxmin(), "k"])
best_sil_k = int(sil_df.loc[sil_df["silhouette"].idxmax(), "k"])
K = 3

# Plot PAC/silhouette
fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
ax[0].plot(pac_df["k"], pac_df["PAC"], "o-", color="#1f77b4")
ax[0].axvline(K, color="red", ls="--", alpha=.5, label=f"k={K}"); ax[0].legend()
ax[0].set_xlabel("k"); ax[0].set_ylabel("PAC"); ax[0].set_title("PAC")
ax[1].plot(sil_df["k"], sil_df["silhouette"], "o-", color="#ff7f0e")
ax[1].axvline(K, color="red", ls="--", alpha=.5)
ax[1].set_xlabel("k"); ax[1].set_ylabel("silhouette"); ax[1].set_title("Silhouette")
plt.tight_layout(); fig.savefig(FIG / "PAC_silhouette_vs_k.png", dpi=200); plt.close()

# k=2,3,4 ecotype save
for k in (2, 3, 4):
    pd.DataFrame({"Kids_First_Biospecimen_ID": X_df.index,
                  "cluster_k": [f"E{c}" for c in consensus[k][1]]}).to_csv(
        OUT / f"ecotype_LM22_main_k{k}.tsv", sep="\t", index=False)

# Annotate k=3
lbl3 = consensus[K][1]
clus_df = pd.DataFrame({"Kids_First_Biospecimen_ID": X_df.index, "cluster": lbl3}).set_index("Kids_First_Biospecimen_ID")
feat_clu = X_df.copy(); feat_clu["cluster"] = lbl3
means = feat_clu.groupby("cluster").mean()
means.to_csv(OUT / "ecotype_LM22_main_k3_feature_means.tsv", sep="\t")

def topcol(c): return int(means[c].idxmax()) if c in means.columns else None
lym = topcol("ssGSEA_T_Cell_Cytotoxicity") or topcol("LM22_T cells CD8")
mye = topcol("ssGSEA_MoTAM_Antunes2021")  or topcol("LM22_Macrophages M2")
all_c = set(int(x) for x in means.index.tolist())
des_set = list(all_c - {lym, mye})
des = int(des_set[0]) if des_set else None
mapping = {lym: "Lymphocyte-inflamed", mye: "Myeloid-dominant"}
if des is not None: mapping[des] = "Immune-desert"
print(f"mapping: {mapping}")
clus_df["ecotype"] = clus_df["cluster"].map(mapping)
clus_df.to_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t")

# PCA + UMAP
pca = PCA(n_components=2).fit_transform(X)
pca_df = pd.DataFrame({"sample": X_df.index, "PC1": pca[:, 0], "PC2": pca[:, 1], "ecotype": clus_df["ecotype"].values})
pca_df.to_csv(OUT / "ecotype_LM22_main_pca.tsv", sep="\t", index=False)
try:
    import umap
    um = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=42).fit_transform(X)
    umap_df = pd.DataFrame({"sample": X_df.index, "UMAP1": um[:, 0], "UMAP2": um[:, 1], "ecotype": clus_df["ecotype"].values})
    umap_df.to_csv(OUT / "ecotype_LM22_main_umap.tsv", sep="\t", index=False)
except Exception as e:
    umap_df = None; print(f"UMAP skip: {e}")

ec_colors = {"Lymphocyte-inflamed": "#3B82F6", "Myeloid-dominant": "#EF4444", "Immune-desert": "#9CA3AF"}
fig, axes = plt.subplots(1, 2 if umap_df is not None else 1, figsize=(11, 4.2))
ax0 = axes[0] if umap_df is not None else axes
for ec, sub in pca_df.groupby("ecotype"):
    ax0.scatter(sub["PC1"], sub["PC2"], c=ec_colors.get(ec, "#888"), label=ec, s=18, alpha=.7, edgecolor="none")
ax0.set_xlabel("PC1"); ax0.set_ylabel("PC2"); ax0.set_title("PCA"); ax0.legend(loc="best", frameon=False, fontsize=8)
if umap_df is not None:
    ax1 = axes[1]
    for ec, sub in umap_df.groupby("ecotype"):
        ax1.scatter(sub["UMAP1"], sub["UMAP2"], c=ec_colors.get(ec, "#888"), label=ec, s=18, alpha=.7, edgecolor="none")
    ax1.set_xlabel("UMAP1"); ax1.set_ylabel("UMAP2"); ax1.set_title("UMAP")
plt.tight_layout(); fig.savefig(FIG / "ecotype_PCA_UMAP.png", dpi=200); plt.close()

# Association tests
ann_full = pd.read_csv(OUT / "cohort_for_deconvolution.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
meta = ann_full.loc[clus_df.index.intersection(ann_full.index)].copy()
meta["ecotype"] = clus_df["ecotype"]

def chi2_p(tab, B=10000):
    obs = tab.values
    if obs.size == 0: return float("nan")
    chi2, _, _, _ = chi2_contingency(obs, correction=False)
    rows = obs.sum(axis=1); cols = obs.sum(axis=0)
    if rows.sum() == 0: return float("nan")
    flat_g = np.repeat(np.arange(obs.shape[0]), rows.astype(int))
    flat_c = np.repeat(np.arange(obs.shape[1]), cols.astype(int))
    null = 0
    for _ in range(B):
        pc = RNG.permutation(flat_c)
        ct = np.zeros_like(obs)
        np.add.at(ct, (flat_g, pc), 1)
        c2, _, _, _ = chi2_contingency(ct, correction=False)
        if c2 >= chi2: null += 1
    return (1 + null) / (B + 1)

assoc = []
for var in ["cohort_group", "location_class", "age_dev_group"]:
    sub = meta[meta["location_class"].isin(["Midline", "Hemispheric", "Posterior_fossa"])] if var == "location_class" else meta
    tab = pd.crosstab(sub["ecotype"], sub[var])
    tab.to_csv(OUT / f"ecotype_LM22_main_x_{var}.tsv", sep="\t")
    print(f"\necotype × {var}:\n{tab.to_string()}")
    p = chi2_p(tab, B=5000)
    print(f"  Monte-Carlo p = {p:.2e}")
    assoc.append({"comparison": f"ecotype_x_{var}", "monte_carlo_p": p, "N": int(tab.values.sum())})
pd.DataFrame(assoc).to_csv(OUT / "ecotype_LM22_main_assoc_pvalues.tsv", sep="\t", index=False)

# Compare with previous
try:
    prev = pd.read_csv(OUT / "ecotype_assignment_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
    common = clus_df.index.intersection(prev.index)
    comp = pd.DataFrame({"LM22+ssGSEA": clus_df.loc[common, "ecotype"], "previous(quanTIseq+ssGSEA)": prev.loc[common, "ecotype"]})
    ct = pd.crosstab(comp["LM22+ssGSEA"], comp["previous(quanTIseq+ssGSEA)"])
    ct.to_csv(OUT / "ecotype_LM22_vs_previous_crosstab.tsv", sep="\t")
    print("\nCrosstab vs previous:\n", ct.to_string())
    agree = (comp["LM22+ssGSEA"] == comp["previous(quanTIseq+ssGSEA)"]).mean()
    print(f"exact label agreement (n={len(common)}): {agree:.1%}")
except Exception as e: print(f"compare skipped: {e}")

summary = {
    "n_main_used": int(X.shape[0]), "n_features": int(X.shape[1]),
    "best_PAC_k": best_pac_k, "best_silhouette_k": best_sil_k, "primary_K": K,
    "ecotype_sizes": {k: int(v) for k, v in Counter(clus_df["ecotype"].tolist()).items()},
}
(OUT / "step5_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
