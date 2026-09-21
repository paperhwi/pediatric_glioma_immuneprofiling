"""Re-run consensus clustering with B=1000 bootstraps (analysis_plan §1.5 spec).
Designed to be called per-k so each invocation fits within bash timeout.

Usage:
  python step5e_bootstrap_1000.py 2     # run k=2
  python step5e_bootstrap_1000.py 3
  ...
  python step5e_bootstrap_1000.py finalize  # PAC/silhouette + annotation
"""
from __future__ import annotations
import sys, json, time, warnings
from pathlib import Path
from collections import Counter
import numpy as np, pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"

feat = pd.read_csv(OUT / "step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t", index_col=0)
main_ids = [s for s in pd.read_csv(OUT / "cohort_main_final.tsv", sep="\t")["Kids_First_Biospecimen_ID"] if s in feat.index]
X_df = feat.loc[main_ids].dropna(axis=0, how="any")
X = X_df.values; n = X.shape[0]

def run_k(k, B=1000, pItem=0.8, n_init=1, seed=42):
    t0 = time.time()
    M = np.zeros((n, n)); C = np.zeros((n, n))
    rng = np.random.default_rng(seed)
    for b in range(B):
        idx = rng.choice(n, int(n * pItem), replace=False)
        km = KMeans(n_clusters=k, n_init=n_init, random_state=seed + b)
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
    final = fcluster(Z, t=k, criterion="maxclust")
    np.savez_compressed(OUT / f"consensus_LM22_main_B1000_k{k}.npz",
                        consensus=cm, labels=final, samples=np.array(X_df.index, dtype=object))
    print(f"k={k}: sizes={dict(Counter(final))}  ({time.time()-t0:.1f}s)")

def finalize():
    """PAC / silhouette / sample-level robustness across the 1000 bootstraps."""
    pac_rows, sil_rows, rob_rows = [], [], []
    for k in [2, 3, 4, 5, 6]:
        z = np.load(OUT / f"consensus_LM22_main_B1000_k{k}.npz", allow_pickle=True)
        cm = z["consensus"]; lbl = z["labels"]
        upper = cm[np.triu_indices(n, k=1)]
        pac = float(((upper > 0.1) & (upper < 0.9)).mean())
        dist = 1.0 - cm; np.fill_diagonal(dist, 0.0)
        sil = float(silhouette_score(dist, lbl, metric="precomputed"))
        # mean within-cluster consensus = sample-level robustness
        same = np.equal.outer(lbl, lbl)
        same_idx = same & ~np.eye(n, dtype=bool)
        wcr = float(cm[same_idx].mean()) if same_idx.any() else float("nan")
        pac_rows.append({"k": k, "PAC_B1000": pac})
        sil_rows.append({"k": k, "silhouette_B1000": sil})
        rob_rows.append({"k": k, "mean_within_cluster_consensus": wcr})
    pd.DataFrame(pac_rows).to_csv(OUT / "consensus_LM22_main_B1000_PAC.tsv", sep="\t", index=False)
    pd.DataFrame(sil_rows).to_csv(OUT / "consensus_LM22_main_B1000_silhouette.tsv", sep="\t", index=False)
    pd.DataFrame(rob_rows).to_csv(OUT / "consensus_LM22_main_B1000_within_cluster_robustness.tsv", sep="\t", index=False)
    print("PAC:");  print(pd.DataFrame(pac_rows).to_string(index=False))
    print("\nSilhouette:");  print(pd.DataFrame(sil_rows).to_string(index=False))
    print("\nWithin-cluster mean consensus (higher = more robust):");  print(pd.DataFrame(rob_rows).to_string(index=False))

if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "finalize":
        finalize()
    else:
        run_k(int(arg), B=1000, n_init=1)
