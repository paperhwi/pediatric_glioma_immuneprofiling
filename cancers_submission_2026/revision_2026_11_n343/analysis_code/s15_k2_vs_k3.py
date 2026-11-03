"""Reviewer: stability metrics favour k = 2; why retain k = 3?
(1) k = 2 labels with the locked consensus settings (as cluster344.py: B = 500, KMeans n_init = 10, seed 42, n = 344) and nesting in k = 3.
(2) Per-cluster consensus stability at k = 3 (which split is unstable).
(3) What k = 3 adds over k = 2: survival (Cox LRT, AIC, optimism-corrected C), feature separation LI vs MD, adjusted DEG LI vs MD.
(4) Whether k = 3 is better described as positions on the continuous immune axis."""
import sys, json; sys.path.insert(0, "/home/claude/q343/src")
from common import *
from kclust import consensus as consensus_n1
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy import stats
from lifelines.utils import concordance_index
R = "/home/claude/repo/cancers_submission_2026/reproducibility_data/"; O = ROOT+"out/r53/"
P = load(); F = pd.read_csv(R + "feature_matrices/main_cohort_LM22_ssGSEA_z.tsv", sep="\t", index_col=0).loc[P.index]
X = F.values; n = len(X)
def consensus_locked(X, k, seed=42, B=500, pItem=0.8):      # identical to cluster344.consensus
    RNG = np.random.default_rng(seed); n = len(X); M = np.zeros((n, n)); C = np.zeros((n, n))
    for b in range(B):
        idx = RNG.choice(n, int(n * pItem), replace=False)
        lbl = KMeans(n_clusters=k, n_init=10, random_state=seed + b).fit_predict(X[idx])
        for ci in range(k):
            m = idx[lbl == ci]
            if len(m) >= 2: M[np.ix_(m, m)] += 1
        C[np.ix_(idx, idx)] += 1
    cm = np.where(C > 0, M / np.maximum(C, 1), 0.); np.fill_diagonal(cm, 1)
    d = 1 - cm; np.fill_diagonal(d, 0)
    return cm, fcluster(linkage(squareform(d, checks=False), "average"), t=k, criterion="maxclust")
cm2, l2 = consensus_locked(X, 2); cm3, l3 = consensus_locked(X, 3)
k2 = pd.Series(l2, index=P.index); k3 = P.state
nest = pd.crosstab(k3.rename("k = 3 state"), k2.rename("k = 2 cluster")); nest.to_csv(O + "R53_nesting.tsv", sep="\t")
# name k=2 clusters by which k=3 states they contain
lab2 = {c: " + ".join(SH for SH, v in zip(["High", "Intermediate", "Low"], [nest.loc[HI, c], nest.loc[MID, c], nest.loc[LO, c]]) if v >= 0.5 * nest[c].sum() or v >= 20) for c in nest.columns}
P["k2"] = k2.map(lab2)
# per-cluster consensus stability (mean within-cluster consensus) at k = 3 and k = 2
stab = []
for kk, cm, lab in [(3, cm3, l3), (2, cm2, l2)]:
    for c in sorted(set(lab)):
        m = np.where(lab == c)[0]; sub = cm[np.ix_(m, m)]; w = sub[np.triu_indices(len(m), 1)].mean()
        o = np.where(lab != c)[0]; between = cm[np.ix_(m, o)].mean()
        stab.append(dict(k=kk, cluster=int(c), n=len(m), states=pd.Series(k3.values[m]).value_counts().to_dict(), within_consensus=w, between_consensus=between))
# pairwise between-state consensus at k=3 using the published labels
pair = {}
for a in ECO:
    for b in ECO:
        ia = np.where(k3.values == a)[0]; ib = np.where(k3.values == b)[0]
        pair[(a, b)] = cm3[np.ix_(ia, ib)].mean()
PAIR = pd.Series(pair).unstack()
# survival: base vs + k2 vs + k3 vs + continuous PC1 (corrected OS, Fig 7 set)
s = P[P.os_eval].copy(); pc1 = pd.read_csv(ROOT+"out/immune_axis_n343.tsv", sep="\t", index_col=0).immune_axis
def dm(df, mode):
    d = design(df, eco=False)
    if mode in ("k2", "k3+k2"):
        lv = sorted(df.k2.unique()); ref = [l for l in lv if l == "Low"] or [lv[0]]
        for l in lv:
            if l != ref[0]: d["k2_" + l] = (df.k2 == l).astype(float)
    if mode == "k3":
        d[HI] = (df.state == HI).astype(float); d[MID] = (df.state == MID).astype(float)
    if mode == "pc1": d["pc1"] = pc1.loc[df.index].values
    return d
fits = {m: cox(dm(s, m)) for m in ["base", "k2", "k3", "pc1"]}
def lrt(a, b): x = 2 * (fits[b].log_likelihood_ - fits[a].log_likelihood_); dfree = len(fits[b].params_) - len(fits[a].params_); return x, dfree, stats.chi2.sf(x, dfree)
rng = np.random.default_rng(7); Bt = 300
def harrell(c, d): return concordance_index(d.time, -c.predict_partial_hazard(d).values.ravel(), d.event)
opt = {m: [] for m in fits}
for b in range(Bt):
    ii = rng.choice(len(s), len(s), replace=True); sb = s.iloc[ii]
    for m in fits:
        try:
            db = dm(sb, m); keep = [c for c in db.columns if c in ("time", "event") or db[c].nunique() > 1]; cb = cox(db[keep], penalizer=1e-4)
            do = dm(s, m)[keep]; opt[m].append(harrell(cb, db[keep]) - harrell(cb, do))
        except Exception: pass
SV = []
for m in fits:
    d = dm(s, m); app = harrell(fits[m], d)
    SV.append(dict(model={"base": "dx + age + sex", "k2": "+ k = 2 clusters", "k3": "+ k = 3 states", "pc1": "+ continuous immune PC1"}[m], df=len(fits[m].params_), logLik=fits[m].log_likelihood_, AIC=fits[m].AIC_partial_,
                   Harrell_C=app, C_optimism_corrected=app - np.mean(opt[m])))
SV = pd.DataFrame(SV)
L = {"k2 vs base": lrt("base", "k2"), "k3 vs base": lrt("base", "k3"), "pc1 vs base": lrt("base", "pc1")}
# k3 vs k2 nested? k3 refines k2 when each k2 cluster is a union of k3 states -> LRT k3 vs k2 valid (models nested)
nested = bool(((nest > 0).sum(0) >= 1).all() and (nest.gt(0).sum(1) == 1).mean() > 0.9)
L["k3 vs k2"] = lrt("k2", "k3")
c3 = fits["k3"].summary.loc[[HI, MID], ["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
# LI vs MD: which features separate them
fs = []
li, md = P.index[P.state == HI], P.index[P.state == MID]
for f in F.columns:
    u = stats.mannwhitneyu(F.loc[li, f], F.loc[md, f]); dlt = F.loc[li, f].mean() - F.loc[md, f].mean()
    fs.append(dict(feature=f, mean_high=F.loc[li, f].mean(), mean_intermediate=F.loc[md, f].mean(), delta=dlt, p=u.pvalue))
FS = pd.DataFrame(fs); FS["q"] = bh(FS.p.values); FS = FS.sort_values("delta"); FS.to_csv(O + "R53_high_vs_intermediate_features.tsv", sep="\t", index=False)
deg = pd.read_csv(ROOT+"out/T_adjusted_DEG_n343_all_contrasts.tsv", sep="\t"); dc = deg[deg.contrast == "Immune-high vs Immune-intermediate"]
ndeg = dict(up_in_high=int(((dc.q_BH < 0.05) & (dc.adjusted_log2FC > 1)).sum()), up_in_intermediate=int(((dc.q_BH < 0.05) & (dc.adjusted_log2FC < -1)).sum()))
top_up = dc[(dc.q_BH < 0.05) & (dc.adjusted_log2FC > 1)].nlargest(15, "adjusted_log2FC").gene.tolist()
top_dn = dc[(dc.q_BH < 0.05) & (dc.adjusted_log2FC < -1)].nsmallest(15, "adjusted_log2FC").gene.tolist()
# PC1 by state (ordering)
PCs = pc1.groupby(P.state).describe()[["count", "25%", "50%", "75%"]]
# ratio of lymphoid to myeloid signal: separates LI from MD independent of total immune amount?
lym = F[["ssGSEA_T_Cell_Cytotoxicity", "ssGSEA_NK_Cell_Activity"]].mean(1); mye = F[["ssGSEA_MDM_Klemm2020", "ssGSEA_MoTAM_Antunes2021", "ssGSEA_Microglia_Klemm2020"]].mean(1)
bal = (lym - mye); balS = bal.groupby(P.state).median()
# matched on PC1: do LI and MD differ in lymphoid-myeloid balance at equal overall immune presence?
import statsmodels.formula.api as smf
dd = pd.DataFrame({"bal": bal, "pc1": pc1, "high": (P.state == HI).astype(float)}).loc[list(li) + list(md)]
fb = smf.ols("bal ~ pc1 + high", dd).fit()
out = dict(k2_sizes=P.k2.value_counts().to_dict(), nesting=nest.to_dict(), k3_refines_k2=nested,
           stability=stab, pairwise_consensus=PAIR.round(3).to_dict(), survival=SV.round(4).to_dict("records"),
           LRT={k: dict(chi2=v[0], df=v[1], p=v[2]) for k, v in L.items()}, k3_HR=c3.round(4).to_dict("index"),
           LI_vs_MD_features_q05=int((FS.q < 0.05).sum()), LI_vs_MD_DEG=ndeg, top_up_LI=top_up, top_up_MD=top_dn,
           PC1_by_state=PCs.round(3).to_dict("index"), lymphoid_minus_myeloid_median=balS.round(3).to_dict(),
           high_effect_on_balance_given_axis=dict(coef=fb.params["high"], p=fb.pvalues["high"]))
json.dump(out, open(O + "R53_summary.json", "w"), indent=1, default=lambda v: float(v) if hasattr(v, "__float__") else str(v))
SV.to_csv(O + "R53_survival_k2_k3.tsv", sep="\t", index=False); pd.DataFrame(stab).to_csv(O + "R53_cluster_stability.tsv", sep="\t", index=False); PAIR.to_csv(O + "R53_pairwise_consensus.tsv", sep="\t")
P[["k2", "state"]].to_csv(O + "R53_labels_k2_k3.tsv", sep="\t")
print(nest); print(SV.round(3).to_string()); print({k: (round(v[0], 2), v[1], round(v[2], 4)) for k, v in L.items()}); print(pd.DataFrame(stab).to_string()); print(PAIR.round(3)); print(ndeg, top_up[:8], top_dn[:8]); print(balS, fb.params["high"], fb.pvalues["high"])
