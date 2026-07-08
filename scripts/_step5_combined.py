# Step 5 - Consensus clustering with LM22+ssGSEA feature matrix
# (combined run-record)


# ===== step5b_finish_consensus.py =====
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

# ===== step5c_annotate_fix.py =====
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

# ===== step5d_finalize_plots.py =====
"""Regenerate PCA/UMAP plots with fixed ecotype labels and build full Step 5 notebook source."""
from pathlib import Path
import json, warnings
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_step5"; FIG.mkdir(exist_ok=True)

# Reload PCA/UMAP coords saved earlier and patch with fixed ecotype labels
ann3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
ann2 = pd.read_csv(OUT / "ecotype_LM22_main_k2_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")

pca = pd.read_csv(OUT / "ecotype_LM22_main_pca.tsv", sep="\t").set_index("sample")
pca["ecotype_k3"] = ann3.reindex(pca.index)["ecotype"]
pca["ecotype_k2"] = ann2.reindex(pca.index)["ecotype"]
pca.reset_index().to_csv(OUT / "ecotype_LM22_main_pca.tsv", sep="\t", index=False)

try:
    um = pd.read_csv(OUT / "ecotype_LM22_main_umap.tsv", sep="\t").set_index("sample")
    um["ecotype_k3"] = ann3.reindex(um.index)["ecotype"]
    um["ecotype_k2"] = ann2.reindex(um.index)["ecotype"]
    um.reset_index().to_csv(OUT / "ecotype_LM22_main_umap.tsv", sep="\t", index=False)
except FileNotFoundError:
    um = None

ec_colors_k3 = {"Lymphocyte-inflamed": "#3B82F6", "Myeloid-dominant": "#EF4444", "Immune-desert": "#9CA3AF"}
ec_colors_k2 = {"Immune-active": "#16A34A", "Immune-desert": "#9CA3AF"}

# 4-panel grid: PCA k=2, PCA k=3, UMAP k=2, UMAP k=3
nrow = 2 if um is not None else 1
fig, axes = plt.subplots(nrow, 2, figsize=(11, 4.5 * nrow))
if nrow == 1:
    axes = np.array([axes])

def scatter(ax, df, x, y, color_col, palette, title):
    for ec, sub in df.groupby(color_col):
        ax.scatter(sub[x], sub[y], c=palette.get(ec, "#888"), label=ec, s=18, alpha=.75, edgecolor="none")
    ax.set_xlabel(x); ax.set_ylabel(y); ax.set_title(title)
    ax.legend(loc="best", frameon=False, fontsize=8)

scatter(axes[0, 0], pca, "PC1", "PC2", "ecotype_k2", ec_colors_k2, "PCA — k=2 (PAC-optimal)")
scatter(axes[0, 1], pca, "PC1", "PC2", "ecotype_k3", ec_colors_k3, "PCA — k=3 (pre-specified)")
if um is not None:
    scatter(axes[1, 0], um.reset_index(), "UMAP1", "UMAP2", "ecotype_k2", ec_colors_k2, "UMAP — k=2")
    scatter(axes[1, 1], um.reset_index(), "UMAP1", "UMAP2", "ecotype_k3", ec_colors_k3, "UMAP — k=3")
plt.tight_layout()
fig.savefig(FIG / "ecotype_PCA_UMAP.png", dpi=200)
plt.close()

# Stacked bar: ecotype distribution by cohort_group for k=3
ct_k3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_x_cohort_group.tsv", sep="\t").set_index("ecotype")
ct_pct = ct_k3.div(ct_k3.sum(axis=0), axis=1) * 100
fig, ax = plt.subplots(figsize=(6.5, 4.2))
bottom = np.zeros(len(ct_pct.columns))
for ec in ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]:
    if ec in ct_pct.index:
        vals = ct_pct.loc[ec].values
        ax.bar(ct_pct.columns, vals, bottom=bottom, label=ec, color=ec_colors_k3[ec])
        bottom += vals
ax.set_ylabel("% of samples within cohort"); ax.set_ylim(0, 100)
ax.set_title("Ecotype composition by cohort group (k=3)")
ax.legend(loc="upper right", fontsize=8, frameon=False)
plt.tight_layout()
fig.savefig(FIG / "ecotype_k3_by_cohort_group.png", dpi=200)
plt.close()

# Heatmap of k=3 cluster feature means
means = pd.read_csv(OUT / "ecotype_LM22_main_k3_feature_means.tsv", sep="\t", index_col=0)
# Map cluster idx → ecotype name
map3 = {2: "Lymphocyte-inflamed", 3: "Myeloid-dominant", 1: "Immune-desert"}
means.index = [f"{map3.get(int(i), int(i))} (c{int(i)})" for i in means.index]
# Show top 30 features by absolute spread
spread = means.max() - means.min()
keep = spread.sort_values(ascending=False).head(30).index
fig, ax = plt.subplots(figsize=(8, 9))
import seaborn as sns
sns.heatmap(means[keep].T, cmap="RdBu_r", center=0, annot=True, fmt=".2f", ax=ax,
            cbar_kws=dict(label="z-score (mean per cluster)"))
ax.set_title("k=3 ecotype × top-30 discriminating features")
plt.tight_layout()
fig.savefig(FIG / "ecotype_k3_feature_means_heatmap.png", dpi=200)
plt.close()

print("plots regenerated:")
for p in FIG.glob("*.png"):
    print(" ", p.relative_to(BASE))
