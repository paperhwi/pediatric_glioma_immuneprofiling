"""
Step 10 — Secondary BRAF_ALT cohort: ecotype projection by nearest-centroid

Strategy:
- Compute z-scored feature vector centroid for each k=3 ecotype on Main cohort.
- For each BRAF_ALT_main (HGG/PXA/infant-fusion) and BRAF_ALT_LGG sample, compute
  Euclidean distance to each centroid and assign nearest.
- Compare ecotype distribution: Main vs BRAF_ALT_main vs BRAF_ALT_LGG.
- Test MAPK_Activity and other myeloid axes by cohort family.

Outputs (output/, figs_step10/)
"""
from __future__ import annotations
from pathlib import Path
import json, warnings
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_step10"; FIG.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
feat = pd.read_csv(OUT / "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv",
                   sep="\t", index_col=0)
ec3  = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv",
                   sep="\t").set_index("Kids_First_Biospecimen_ID")

main_ids = ec3.index.intersection(feat.index)
print(f"Main cohort samples with ecotype: {len(main_ids)}")

# Compute z-scored feature centroids per ecotype in Main cohort
centroids = feat.loc[main_ids].join(ec3[["ecotype"]]).groupby("ecotype").mean()
ec_order = ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]
centroids = centroids.reindex(ec_order)
centroids.to_csv(OUT / "step10_ecotype_centroids_z.tsv", sep="\t")

# Load BRAF cohort lists
braf_main = pd.read_csv(OUT / "cohort_BRAF_ALT_main.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
braf_lgg  = pd.read_csv(OUT / "cohort_BRAF_ALT_LGG.tsv",  sep="\t").set_index("Kids_First_Biospecimen_ID")
print(f"BRAF_ALT_main(HGG/PXA/infant-fusion): {len(braf_main)}")
print(f"BRAF_ALT_LGG: {len(braf_lgg)}")

def assign(ids, label):
    use = [i for i in ids if i in feat.index]
    X   = feat.loc[use].dropna(axis=0, how="any")
    dist = np.linalg.norm(X.values[:, None, :] - centroids.values[None, :, :], axis=2)
    assign_idx = np.argmin(dist, axis=1)
    nearest = pd.Series([centroids.index[i] for i in assign_idx], index=X.index, name="ecotype")
    out = pd.DataFrame({"Kids_First_Biospecimen_ID": X.index,
                        "ecotype": nearest.values,
                        "set": label,
                        "dist_Lymph": dist[:, 0], "dist_Mye": dist[:, 1], "dist_Desert": dist[:, 2]})
    return out

asg_main = assign(braf_main.index, "BRAF_ALT_main")
asg_lgg  = assign(braf_lgg.index,  "BRAF_ALT_LGG")
all_assign = pd.concat([asg_main, asg_lgg], ignore_index=True)
all_assign.to_csv(OUT / "step10_BRAF_ALT_ecotype_assignment.tsv", sep="\t", index=False)
print(f"\nBRAF_ALT_main n={len(asg_main)}  BRAF_ALT_LGG n={len(asg_lgg)}")
print("\nBRAF_ALT_main ecotype distribution:")
print(asg_main["ecotype"].value_counts())
print("\nBRAF_ALT_LGG ecotype distribution:")
print(asg_lgg["ecotype"].value_counts())

# Also compute Main reference distribution for comparison
main_dist = ec3.loc[main_ids, "ecotype"].value_counts().reindex(ec_order, fill_value=0)
ml_main = asg_main["ecotype"].value_counts().reindex(ec_order, fill_value=0)
ml_lgg  = asg_lgg["ecotype"].value_counts().reindex(ec_order, fill_value=0)
dist_table = pd.DataFrame({"Main_n349": main_dist, "BRAF_ALT_main": ml_main, "BRAF_ALT_LGG": ml_lgg})
dist_table_pct = dist_table.div(dist_table.sum(axis=0), axis=1) * 100
print("\nEcotype distribution by cohort family (%):")
print(dist_table_pct.round(1).to_string())
dist_table.to_csv(OUT / "step10_ecotype_distribution_by_cohort_family.tsv", sep="\t")

# Chi-square test
chi2, p_chi, _, _ = stats.chi2_contingency(dist_table.values)
print(f"\nChi-square Main vs BRAF_main vs BRAF_LGG: chi2={chi2:.2f}, p={p_chi:.2e}")

# ---------------------------------------------------------------------------
# MAPK_Activity comparison: Main vs BRAF_ALT_main vs BRAF_ALT_LGG
# ---------------------------------------------------------------------------
merged = pd.read_csv(OUT / "step4_integrated_feature_matrix.tsv", sep="\t", index_col=0)
samples_by_set = {
    "Main_n349":      main_ids.tolist(),
    "BRAF_ALT_main":  asg_main["Kids_First_Biospecimen_ID"].tolist(),
    "BRAF_ALT_LGG":   asg_lgg["Kids_First_Biospecimen_ID"].tolist(),
}
mapk = []
for s, ids in samples_by_set.items():
    vals = merged.loc[merged.index.intersection(ids), "ssGSEA_MAPK_Activity"].dropna().values
    mapk.append(pd.DataFrame({"set": s, "MAPK_Activity": vals}))
mapk_df = pd.concat(mapk, ignore_index=True)

H, p_kw = stats.kruskal(*[g["MAPK_Activity"].values for _, g in mapk_df.groupby("set")])
print(f"\nMAPK_Activity KW across cohort families: H={H:.2f}, p={p_kw:.2e}")

# ---------------------------------------------------------------------------
# Plot: ecotype distribution stacked bar
# ---------------------------------------------------------------------------
palette = {"Lymphocyte-inflamed":"#3B82F6","Myeloid-dominant":"#EF4444","Immune-desert":"#9CA3AF"}
fig, ax = plt.subplots(figsize=(7, 4.5))
bottom = np.zeros(len(dist_table_pct.columns))
for ec in ec_order:
    vals = dist_table_pct.loc[ec].values
    ax.bar(dist_table_pct.columns, vals, bottom=bottom, label=ec, color=palette[ec])
    bottom += vals
ax.set_ylabel("% of samples"); ax.set_ylim(0, 100)
ax.set_title(f"Ecotype distribution by cohort family\nChi-square p = {p_chi:.2e}")
ax.legend(loc="upper right", fontsize=8, frameon=False)
plt.tight_layout(); fig.savefig(FIG / "ecotype_distribution_by_cohort_family.png", dpi=200); plt.close()

# ---------------------------------------------------------------------------
# Plot: MAPK_Activity boxplot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 4.5))
sns.boxplot(data=mapk_df, x="set", y="MAPK_Activity",
            order=["Main_n349", "BRAF_ALT_main", "BRAF_ALT_LGG"],
            palette={"Main_n349":"#0EA5E9","BRAF_ALT_main":"#F59E0B","BRAF_ALT_LGG":"#10B981"},
            ax=ax, showfliers=False)
sns.stripplot(data=mapk_df, x="set", y="MAPK_Activity",
              order=["Main_n349", "BRAF_ALT_main", "BRAF_ALT_LGG"],
              color="black", size=2, alpha=.3, ax=ax)
ax.set_xlabel(""); ax.set_ylabel("ssGSEA MAPK_Activity (z-score)")
ax.set_title(f"MAPK activity across cohort families  (KW p={p_kw:.1e})")
plt.tight_layout(); fig.savefig(FIG / "MAPK_by_cohort_family.png", dpi=200); plt.close()

# ---------------------------------------------------------------------------
# Plot: PCA projection of Main + BRAF onto Main PCA space
# ---------------------------------------------------------------------------
from sklearn.decomposition import PCA
common_main_feat = feat.loc[main_ids].dropna(axis=0, how="any")
pca = PCA(n_components=2).fit(common_main_feat.values)
def project(df_sub):
    return pca.transform(df_sub.loc[df_sub.index.intersection(feat.index)].values)

main_pc = pca.transform(common_main_feat.values)
braf_main_feat = feat.loc[asg_main["Kids_First_Biospecimen_ID"]].dropna(axis=0, how="any")
braf_lgg_feat  = feat.loc[asg_lgg["Kids_First_Biospecimen_ID"]].dropna(axis=0, how="any")
braf_main_pc = pca.transform(braf_main_feat.values)
braf_lgg_pc  = pca.transform(braf_lgg_feat.values)

fig, ax = plt.subplots(figsize=(8, 6))
# Color Main by ecotype
ec_arr = ec3.loc[common_main_feat.index, "ecotype"].values
for ec in ec_order:
    mask = ec_arr == ec
    ax.scatter(main_pc[mask, 0], main_pc[mask, 1], c=palette[ec], s=16, alpha=.6, label=f"Main · {ec}")
ax.scatter(braf_main_pc[:, 0], braf_main_pc[:, 1], facecolors="none", edgecolors="#F59E0B",
           s=42, lw=1.5, label="BRAF_ALT_main (HGG/PXA/infant-fusion)")
ax.scatter(braf_lgg_pc[:, 0], braf_lgg_pc[:, 1], facecolors="none", edgecolors="#10B981",
           s=42, lw=1.5, label="BRAF_ALT_LGG")
ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.set_title("Main cohort PCA space — BRAF_ALT overlay")
ax.legend(loc="best", fontsize=8, frameon=False)
plt.tight_layout(); fig.savefig(FIG / "BRAF_ALT_overlay_on_main_PCA.png", dpi=200); plt.close()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
summary = {
    "n_BRAF_ALT_main":  int(len(asg_main)),
    "n_BRAF_ALT_LGG":   int(len(asg_lgg)),
    "ecotype_pct_Main_n349":     dist_table_pct["Main_n349"].round(1).to_dict(),
    "ecotype_pct_BRAF_ALT_main": dist_table_pct["BRAF_ALT_main"].round(1).to_dict(),
    "ecotype_pct_BRAF_ALT_LGG":  dist_table_pct["BRAF_ALT_LGG"].round(1).to_dict(),
    "chi_square_p_Main_vs_BRAF_families": float(p_chi),
    "MAPK_KW_p_across_families":           float(p_kw),
}
(OUT / "step10_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
