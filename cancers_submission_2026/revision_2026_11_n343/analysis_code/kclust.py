"""Consensus clustering and k-selection criteria, method-matched to the locked pipeline
(t1_ksweep.py / step5e_bootstrap_1000.py: KMeans n_init=1, B bootstraps, 80% item resampling, seed 42,
average-linkage on 1 - consensus; gap statistic as in t1_gap.py)."""
import numpy as np, pandas as pd
from collections import Counter
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, adjusted_rand_score
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.optimize import linear_sum_assignment

def consensus(X, k, B=1000, pItem=0.8, seed=42):
    n = X.shape[0]; M = np.zeros((n, n)); C = np.zeros((n, n)); rng = np.random.default_rng(seed)
    for b in range(B):
        idx = rng.choice(n, int(n * pItem), replace=False)
        lbl = KMeans(n_clusters=k, n_init=1, random_state=seed + b).fit_predict(X[idx])
        C[np.ix_(idx, idx)] += 1
        for ci in range(k):
            mem = idx[lbl == ci]
            if len(mem) >= 2: M[np.ix_(mem, mem)] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        cm = np.where(C > 0, M / C, 0.0)
    np.fill_diagonal(cm, 1.0); dist = 1.0 - cm; np.fill_diagonal(dist, 0.0)
    lab = fcluster(linkage(squareform(dist, checks=False), method="average"), t=k, criterion="maxclust")
    return cm, dist, lab

def cdf_auc(cm):
    n = cm.shape[0]; x = np.sort(cm[np.triu_indices(n, 1)]); cdf = np.arange(1, len(x) + 1) / len(x)
    return float(np.sum(np.diff(np.concatenate(([0.0], x))) * cdf))

def sweep(X, ks=range(2, 7), B=1000, seed=42):
    rows, labs, prev = [], {}, None
    for k in ks:
        cm, dist, lab = consensus(X, k, B=B, seed=seed); n = len(lab)
        up = cm[np.triu_indices(n, 1)]; auc = cdf_auc(cm)
        same = np.equal.outer(lab, lab) & ~np.eye(n, dtype=bool); sizes = sorted(Counter(lab).values(), reverse=True)
        rows.append(dict(k=k, PAC=float(((up > 0.1) & (up < 0.9)).mean()), silhouette=float(silhouette_score(dist, lab, metric="precomputed")),
                         within_consensus=float(cm[same].mean()), CDF_AUC=auc, delta_AUC=(auc - prev) / prev if prev else np.nan,
                         min_size=min(sizes), sizes=";".join(map(str, sizes))))
        prev = auc; labs[k] = lab
    return pd.DataFrame(rows), labs

def gap(X, kmax=8, B=50, seed=42):
    n, p = X.shape; rng = np.random.default_rng(seed)
    def Wk(D, k):
        if k == 1: return float(((D - D.mean(0)) ** 2).sum())
        km = KMeans(n_clusters=k, n_init=10, random_state=42).fit(D)
        return float(sum(((D[km.labels_ == i] - km.cluster_centers_[i]) ** 2).sum() for i in range(k)))
    pca = PCA().fit(X); Xp = X @ pca.components_.T; lo, hi = Xp.min(0), Xp.max(0); ks = list(range(1, kmax + 1))
    logW = np.array([np.log(Wk(X, k)) for k in ks]); ref = np.zeros((B, len(ks)))
    for b in range(B):
        Z = rng.uniform(lo, hi, size=(n, p)) @ pca.components_; ref[b] = [np.log(Wk(Z, k)) for k in ks]
    g = ref.mean(0) - logW; sk = ref.std(0) * np.sqrt(1 + 1 / B)
    out = pd.DataFrame({"k": ks, "gap": g, "s_k": sk})
    out["meets_1SE"] = [bool(g[i] >= g[i + 1] - sk[i + 1]) for i in range(len(ks) - 1)] + [False]
    sel = out.k[out.meets_1SE]; return out, (int(sel.iloc[0]) if len(sel) else None)

def match_to_reference(Xdf, lab, ref_centroids):
    """Name clusters by maximal centroid correlation with the global state centroids (Hungarian assignment)."""
    cl = sorted(set(lab)); cen = pd.DataFrame({c: Xdf[lab == c].mean() for c in cl}).T
    R = np.array([[np.corrcoef(cen.loc[c], ref_centroids.loc[s])[0, 1] for s in ref_centroids.index] for c in cl])
    if len(cl) <= len(ref_centroids):
        r, cidx = linear_sum_assignment(-R); name = {cl[i]: ref_centroids.index[j] for i, j in zip(r, cidx)}
    else:
        name = {c: ref_centroids.index[int(np.argmax(R[i]))] + f" ({c})" for i, c in enumerate(cl)}
    return pd.Series(lab, index=Xdf.index).map(name), pd.DataFrame(R, index=[f"cluster {c}" for c in cl], columns=ref_centroids.index)
