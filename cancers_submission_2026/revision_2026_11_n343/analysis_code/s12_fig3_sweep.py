"""Reproduce and rerun the submitted k=2..10 consensus/gap diagnostics."""

import json
import time
from collections import Counter
from pathlib import Path

import sys; sys.path.insert(0,"/home/claude/q343/src"); from common import load, ROOT as QR
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, adjusted_rand_score

E = Path("/mnt/user-data/uploads/Open PBTA")
OUT = Path(QR) / "out/fig3"
OUT.mkdir(exist_ok=True)
feat = pd.read_csv(E / "Revision/Week1/_inputs/step4_clustering_feature_matrix_LM22_plus_ssGSEA_z.tsv", sep="\t", index_col=0)
old_ids = pd.read_csv(E / "Revision/Week1/_inputs/cohort_main_final.tsv", sep="\t")["Kids_First_Biospecimen_ID"].tolist()
current_obs = load()
current_set = set(current_obs.index)
new_ids = [s for s in old_ids if s in current_set]
assert len(old_ids) == 349 and len(new_ids) == 343 and len(current_set) == 343


def consensus(X, k, B=1000, pitem=.8, seed=42):
    n = len(X)
    m = np.zeros((n, n), np.uint16)
    c = np.zeros((n, n), np.uint16)
    rng = np.random.default_rng(seed)
    for b in range(B):
        idx = rng.choice(n, int(n * pitem), replace=False)
        labels = KMeans(n_clusters=k, n_init=1, random_state=seed + b).fit_predict(X[idx])
        c[np.ix_(idx, idx)] += 1
        for ci in range(k):
            members = idx[labels == ci]
            if len(members) >= 2:
                m[np.ix_(members, members)] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        cm = np.where(c > 0, m / c, 0.0)
    np.fill_diagonal(cm, 1.0)
    dist = 1 - cm
    np.fill_diagonal(dist, 0)
    labels = fcluster(linkage(squareform(dist, checks=False), method="average"), t=k, criterion="maxclust")
    return cm, dist, labels


def auc(cm):
    v = np.sort(cm[np.triu_indices(len(cm), 1)])
    cdf = np.arange(1, len(v) + 1) / len(v)
    return float(np.sum(np.diff(np.concatenate(([0.0], v))) * cdf))


def run_k(X, k):
    t = time.time()
    cm, dist, labels = consensus(X, k)
    up = cm[np.triu_indices(len(X), 1)]
    result = {
        "k": k,
        "PAC": float(((up > .1) & (up < .9)).mean()),
        "silhouette_consensus": float(silhouette_score(dist, labels, metric="precomputed")),
        "consensus_CDF_AUC": auc(cm),
        "min_cluster_size": min(Counter(labels).values()),
        "sizes": ";".join(map(str, sorted(Counter(labels).values(), reverse=True))),
    }
    print("k", k, "n", len(X), "PAC", round(result["PAC"], 4), "sil", round(result["silhouette_consensus"], 4), "sec", round(time.time()-t, 1), flush=True)
    return result, labels


old_x = feat.loc[old_ids].to_numpy()
old_locked = {2: (.24264730099133813, .9361341280993414),
              3: (.3841023614267365, .7902823685847833)}
old_check = {}
for k, (pac_target, sil_target) in old_locked.items():
    row, _ = run_k(old_x, k)
    old_check[k] = {"pac_diff": abs(row["PAC"] - pac_target),
                    "sil_diff": abs(row["silhouette_consensus"] - sil_target)}
    if old_check[k]["pac_diff"] > 1e-10 or old_check[k]["sil_diff"] > 1e-10:
        raise RuntimeError(f"Submitted k={k} reproduction failed: {old_check[k]}")

X = feat.loc[new_ids].to_numpy()
rows, labels_by_k = [], {}
for k in range(2, 11):
    row, labels = run_k(X, k)
    rows.append(row)
    labels_by_k[k] = labels
sweep = pd.DataFrame(rows)
sweep["delta_AUC"] = sweep.consensus_CDF_AUC.pct_change()
sweep.to_csv(OUT / "ksweep_n343.tsv", sep="\t", index=False)

current_labels = current_obs.loc[new_ids, "state"].to_numpy()
ari = float(adjusted_rand_score(current_labels, labels_by_k[3]))
ct = pd.crosstab(pd.Series(current_labels, name="current_state"), pd.Series(labels_by_k[3], name="de_novo_k3"))
ct.to_csv(OUT / "current_vs_denovo_k3.tsv", sep="\t")

summary = {"submitted_reproduction": old_check, "n_current": 343, "n_features": X.shape[1],
           "best_PAC_k": int(sweep.loc[sweep.PAC.idxmin(), "k"]),
           "best_silhouette_k": int(sweep.loc[sweep.silhouette_consensus.idxmax(), "k"]),
           "current_vs_denovo_k3_ARI": ari}
(OUT / "consensus_validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2), flush=True)

# Gap statistic: PCA-rotated uniform reference and the original 50-reference
# one-standard-error rule. This is intentionally kept separate from the
# consensus sweep because it uses a distinct KMeans fit.
def wk(data, k):
    if k == 1:
        centre = data.mean(0)
        return float(((data - centre) ** 2).sum())
    model = KMeans(n_clusters=k, n_init=10, random_state=42).fit(data)
    return float(model.inertia_)


start = time.time()
n, p = X.shape
rng = np.random.default_rng(42)
pca = PCA().fit(X)
xp = X @ pca.components_.T
lo, hi = xp.min(0), xp.max(0)
ks = list(range(1, 11))
logw = np.log([wk(X, k) for k in ks])
refs = np.zeros((50, len(ks)))
for b in range(50):
    zp = rng.uniform(lo, hi, size=(n, p))
    z = zp @ pca.components_
    refs[b] = np.log([wk(z, k) for k in ks])
    if (b+1) % 10 == 0:
        print("gap reference", b+1, "of 50", "sec", round(time.time()-start, 1), flush=True)
gap = refs.mean(0) - logw
sk = refs.std(0) * np.sqrt(1 + 1/50)
one_se = [bool(gap[i] >= gap[i+1] - sk[i+1]) for i in range(9)] + [False]
gap_table = pd.DataFrame({"k": ks, "logW": logw, "gap": gap, "s_k": sk, "meets_1SE_rule": one_se})
gap_table.to_csv(OUT / "gap_n343.tsv", sep="\t", index=False)
selected = gap_table.loc[gap_table.meets_1SE_rule, "k"]
summary["gap_1SE_selected_k"] = int(selected.iloc[0]) if len(selected) else None
(OUT / "consensus_validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2), flush=True)
