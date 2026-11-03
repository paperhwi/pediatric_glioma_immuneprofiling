"""Gap G4 (found while building Table 1) — RNA library type is strongly associated with immune state
(exome-capture libraries: 51/73 LI; stranded total RNA: 56/126 ID). Checks whether the states and the survival association are
library artefacts:
(1) consensus clustering within each library with >= 70 tumours (same 46 features re-standardised within library, B = 1000,
    k = 3, clusters named by centroid correlation with the global states) and after removing library means (n = 344);
(2) immune PC1 ordering within each library; (3) Cox model with library as covariate and as strata;
(4) variance of the 46 features explained by library vs state."""
import sys, json, os; sys.path.insert(0, "/home/claude/q343/src")
from common import *
from kclust import sweep, consensus, match_to_reference
from sklearn.metrics import adjusted_rand_score
from scipy import stats
import statsmodels.formula.api as smf
R = "/home/claude/repo/cancers_submission_2026/reproducibility_data/"; O = ROOT+"out/r6lib/r6/"; os.makedirs(O, exist_ok=True)
P = load(); F = pd.read_csv(R + "feature_matrices/main_cohort_LM22_ssGSEA_z.tsv", sep="\t", index_col=0).loc[P.index]
Fz = (F - F.mean()) / F.std(); cent = Fz.groupby(P.state).mean().loc[[HI, MID, LO]]
pc1 = pd.read_csv(ROOT+"out/immune_axis_n343.tsv", sep="\t", index_col=0).immune_axis.loc[P.index]; lib = P.RNA_library
res = {"crosstab": pd.crosstab(lib, P.state).to_dict()}; rows = []
for L in ["stranded", "poly-A stranded", "exome_capture"]:
    ids = P.index[lib == L]; X = F.loc[ids]; X = (X - X.mean()) / X.std()
    S, labs = sweep(X.values, ks=range(2, 5), B=1000)
    named, Rm = match_to_reference(X, labs[3], cent); g = P.loc[ids, "state"]
    res[L] = dict(n=len(ids), ARI_k3=adjusted_rand_score(g, named), agreement_k3=float((named == g).mean()), PAC=S.set_index("k").PAC.round(3).to_dict(),
                  silhouette=S.set_index("k").silhouette.round(3).to_dict(), crosstab=pd.crosstab(g, named).to_dict())
    pd.DataFrame({"global_state": g, "within_library_k3": named}).to_csv(O + f"R6_library_labels_{L.replace(' ', '_')}.tsv", sep="\t")
    rows += [dict(library=L, state=e, n=int((g == e).sum()), median_PC1=float(pc1[ids][g == e].median())) for e in ECO]
    print(L, json.dumps({k: v for k, v in res[L].items() if k != "crosstab"}, default=float))
# library-residualised global clustering
X = F - F.groupby(lib).transform("mean"); X = (X - X.mean()) / X.std()
cm, dist, labr = consensus(X.values, 3, B=1000); named, Rm = match_to_reference(X, labr, cent)
res["library_residualised_global"] = dict(n=len(P), ARI_k3=adjusted_rand_score(P.state, named), agreement_k3=float((named == P.state).mean()),
                                          crosstab=pd.crosstab(P.state, named).to_dict(), library_by_new_state=pd.crosstab(lib, named).to_dict())
pd.DataFrame({"global_state": P.state, "library_residualised_k3": named, "library": lib}).to_csv(O + "R6_library_labels_residualised.tsv", sep="\t")
# PC1 ordering within library (JT)
def jt(groups):
    U = 0; ns = [len(x) for x in groups]
    for i in range(3):
        for j in range(i + 1, 3):
            a, b = np.asarray(groups[i]), np.asarray(groups[j]); U += (a[:, None] < b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum()
    N = sum(ns); mu = (N ** 2 - sum(n ** 2 for n in ns)) / 4; var = (N ** 2 * (2 * N + 3) - sum(n ** 2 * (2 * n + 3) for n in ns)) / 72
    return 2 * stats.norm.sf(abs((U - mu) / np.sqrt(var)))
PCT = pd.DataFrame(rows); res["PC1_JT_within_library"] = {L: jt([pc1[(lib == L) & (P.state == e)] for e in ECO]) for L in lib.unique()}
# survival
s = P[P.os_eval]; base = design(s); c0 = cox(base)
dl = pd.get_dummies(s.RNA_library, prefix="lib", drop_first=True).astype(float); c1 = cox(design(s, extra=dl))
d3 = design(s); d3["lib"] = s.RNA_library.loc[d3.index].values; c2 = CoxPHFitter().fit(d3, "time", "event", strata=["lib"])
SV = pd.DataFrame([dict(model=m, term=t, HR=float(c.summary.loc[t, "exp(coef)"]), lo=float(c.summary.loc[t, "exp(coef) lower 95%"]), hi=float(c.summary.loc[t, "exp(coef) upper 95%"]), p=float(c.summary.loc[t, "p"]))
                   for m, c in [("group + age + sex", c0), ("+ library (covariate)", c1), ("stratified by library", c2)] for t in [HI, MID]])
# variance partition
VR = []
for f in F.columns:
    d = pd.DataFrame({"y": F[f], "l": lib, "s": P.state})
    rl = smf.ols("y ~ C(l)", d).fit().rsquared; rs = smf.ols("y ~ C(s)", d).fit().rsquared; rb = smf.ols("y ~ C(l) + C(s)", d).fit().rsquared
    VR.append(dict(feature=f, R2_library=rl, R2_state=rs, R2_state_given_library=rb - rl, R2_library_given_state=rb - rs))
VR = pd.DataFrame(VR); res["variance_partition_median"] = VR.drop(columns="feature").median().round(3).to_dict()
res["survival"] = SV.round(4).to_dict("records")
PCT.to_csv(O + "R6_library_PC1_by_state.tsv", sep="\t", index=False); SV.to_csv(O + "R6_library_survival.tsv", sep="\t", index=False); VR.to_csv(O + "R6_library_variance_partition.tsv", sep="\t", index=False)
json.dump(res, open(O + "R6_library_summary.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in res.items() if k not in ["crosstab"]}, indent=1, default=float)[:5000])
