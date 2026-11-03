"""Reviewer: would similar immune states emerge within individual diagnostic groups?
(1) Consensus clustering separately within DMG (n = 175) and pHGG-WT (n = 127): same 46 features, re-standardised within the
group, same consensus settings as the locked pipeline (B = 1000, 80% resampling, KMeans n_init = 1, seed 42), k = 2-6,
gap statistic (k = 1-8). (2) Global clustering after removing molecular-group means from every feature (n = 343).
Within-group clusters are named by centroid correlation with the global state centroids."""
import sys, json; sys.path.insert(0, "/home/claude/q343/src")
from common import *
from kclust import sweep, gap, consensus, match_to_reference
from sklearn.metrics import adjusted_rand_score
R = "/home/claude/repo/cancers_submission_2026/reproducibility_data/"; O = ROOT+"out/r54/"; import os; os.makedirs(O, exist_ok=True)
P = load(); F = pd.read_csv(R + "feature_matrices/main_cohort_LM22_ssGSEA_z.tsv", sep="\t", index_col=0).loc[P.index]
Fz = (F - F.mean()) / F.std()
cent = Fz.groupby(P.state).mean().loc[ORD]
res = {}; rows = []
for g, nm in [("DMG_K27", "DMG"), ("pHGG_WT", "pHGG-WT")]:
    ids = P.index[P.cohort_group == g]; X = F.loc[ids]; X = (X - X.mean()) / X.std()
    S, labs = sweep(X.values, ks=range(2, 7), B=1000); G, ksel = gap(X.values, kmax=8)
    S["group"] = nm; G["group"] = nm; S.to_csv(O + f"R54_sweep_{nm}.tsv", sep="\t", index=False); G.to_csv(O + f"R54_gap_{nm}.tsv", sep="\t", index=False)
    named, Rm = match_to_reference(X, labs[3], cent)
    glob_ = P.loc[ids, "state"]
    ari = adjusted_rand_score(glob_, named); agree = float((named == glob_).mean())
    ct = pd.crosstab(glob_.rename("global state"), named.rename(f"within-{nm} cluster"))
    # k = 2 within group: does it split ID from the rest, as globally?
    named2, R2 = match_to_reference(X, labs[2], cent)
    res[nm] = dict(n=len(ids), gap_k=ksel, k3_sizes=named.value_counts().to_dict(), ARI_k3=ari, agreement_k3=agree,
                   PAC=S.set_index("k").PAC.round(3).to_dict(), silhouette=S.set_index("k").silhouette.round(3).to_dict(), delta_AUC=S.set_index("k").delta_AUC.round(3).to_dict(),
                   centroid_r_matched=float(np.mean([Rm.values[i, list(cent.index).index(v)] for i, v in enumerate(pd.Series(named).groupby(pd.Series(labs[3], index=ids)).first().values)])))
    Rm.to_csv(O + f"R54_centroid_correlation_{nm}.tsv", sep="\t"); ct.to_csv(O + f"R54_crosstab_{nm}.tsv", sep="\t")
    pd.DataFrame({"global_state": glob_, "within_group_k3": named, "within_group_k2": named2}).to_csv(O + f"R54_labels_{nm}.tsv", sep="\t")
    print(nm, json.dumps(res[nm], default=float)); print(ct); print(Rm.round(2))
# (2) group-residualised global clustering
ids = P.index; X = F.loc[ids]
Xr = X - X.groupby(P.loc[ids, "cohort_group"]).transform("mean"); Xr = (Xr - Xr.mean()) / Xr.std()
cm, dist, lab = consensus(Xr.values, 3, B=1000); named, Rm = match_to_reference(Xr, lab, cent)
ari = adjusted_rand_score(P.loc[ids, "state"], named); agree = float((named == P.loc[ids, "state"]).mean())
res["group_residualised_global"] = dict(n=len(ids), ARI_k3=ari, agreement_k3=agree, sizes=named.value_counts().to_dict())
pd.crosstab(P.loc[ids, "state"].rename("global state"), named.rename("residualised cluster")).to_csv(O + "R54_crosstab_residualised.tsv", sep="\t")
pd.DataFrame({"global_state": P.loc[ids, "state"], "residualised_k3": named, "cohort_group": P.loc[ids, "cohort_group"]}).to_csv(O + "R54_labels_residualised.tsv", sep="\t")
# variance of each feature explained by molecular group vs by immune state
import statsmodels.formula.api as smf
vr = []
for f in F.columns:
    d = pd.DataFrame({"y": F.loc[ids, f], "g": P.loc[ids, "cohort_group"], "s": P.loc[ids, "state"]})
    r_g = smf.ols("y ~ C(g)", d).fit().rsquared; r_s = smf.ols("y ~ C(s)", d).fit().rsquared; r_gs = smf.ols("y ~ C(g) + C(s)", d).fit().rsquared
    vr.append(dict(feature=f, R2_group=r_g, R2_state=r_s, R2_state_given_group=r_gs - r_g, R2_group_given_state=r_gs - r_s))
VR = pd.DataFrame(vr); VR.to_csv(O + "R54_variance_partition.tsv", sep="\t", index=False)
res["variance_partition_median"] = VR[["R2_group", "R2_state", "R2_state_given_group", "R2_group_given_state"]].median().round(3).to_dict()
json.dump(res, open(O + "R54_summary.json", "w"), indent=1, default=float); print(json.dumps(res, indent=1, default=float))
