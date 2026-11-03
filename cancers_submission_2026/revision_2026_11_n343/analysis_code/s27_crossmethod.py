"""Reviewer: the five deconvolution methods estimate different quantities; how were they harmonised and compared?
Uses the deposited legacy outputs (repository commit 13a4824^): CIBERSORTx Job15 (LM22, absolute mode), immunedeconv
MCP-counter, quanTIseq (tumor = TRUE), EPIC (tumor = TRUE), xCell. Restricted to the n = 344 cohort.
(1) What each method outputs (quantity, unit, which comparisons are valid).
(2) Axis mapping used in Figure 2B/C: subtypes summed within an axis; between-tumour Spearman (scale-free), BH across all pairs.
(3) New scale-free analyses: (a) each method's estimates converted to within-method percentile ranks across tumours;
    (b) ordering across states tested with Jonckheere-Terpstra (ID < MD < LI) and Kruskal-Wallis; Cliff's delta LI vs ID;
    (c) each method's total-immune quantity vs immune PC1 (Spearman)."""
import sys, json; sys.path.insert(0, "/home/claude/q343/src")
from common import *
from scipy import stats
L = "/home/claude/legacy/"; O = ROOT+"out/r52/"
P = load(); ids = list(P.index); pc1 = pd.read_csv(ROOT+"out/immune_axis_n343.tsv", sep="\t", index_col=0).immune_axis
cb = pd.read_csv(L + "CIBERSORTx_Job15_Results.csv").set_index("Mixture"); qc = ["P-value", "Correlation", "RMSE", "Absolute score (sig.score)"]
cells = [c for c in cb.columns if c not in qc]; rel = cb[cells].div(cb[cells].sum(1), axis=0)
D = {m: pd.read_csv(L + f"immunedeconv_{m}.tsv", sep="\t").set_index("cell_type").T for m in ["xcell", "quantiseq", "mcp_counter", "epic"]}
METH = pd.DataFrame([
 dict(method="CIBERSORTx (LM22)", output="Absolute-mode scores re-normalised per tumour to relative fractions (sum to 1 over 22 immune types)", unit="fraction of the immune compartment", reference="LM22 (blood leukocyte profiles)",
      within_tumour_comparison="between immune types: yes (relative)", across_tumours="yes for relative composition; no information on total immune amount", note="no tumour/brain cell types in the reference"),
 dict(method="quanTIseq", output="Absolute fractions of all cells incl. 'uncharacterized'", unit="fraction of all cells", reference="TIL10 (blood/tumour immune)",
      within_tumour_comparison="yes", across_tumours="yes", note="uncharacterized absorbs tumour and brain cells; tumor = TRUE"),
 dict(method="EPIC", output="Absolute fractions of all cells incl. 'uncharacterized'", unit="fraction of all cells (mRNA-corrected)", reference="blood and tumour-infiltrating immune, CAF, endothelial",
      within_tumour_comparison="yes", across_tumours="yes", note="tumor = TRUE"),
 dict(method="MCP-counter", output="Marker-gene abundance scores", unit="arbitrary units, cell-type specific", reference="transcriptomic markers",
      within_tumour_comparison="no (scores of different cell types are on different scales)", across_tumours="yes, within each cell type", note="not a fraction"),
 dict(method="xCell", output="Adjusted enrichment scores (spillover-compensated ssGSEA)", unit="arbitrary units, cell-type specific", reference="64 cell-type signatures",
      within_tumour_comparison="no", across_tumours="yes, within each cell type", note="not a fraction; needs heterogeneous cohort"),
 dict(method="Brain-tuned ssGSEA (this study)", output="Normalised ssGSEA enrichment", unit="arbitrary units, signature specific", reference="24 curated brain/glioma immune signatures",
      within_tumour_comparison="no", across_tumours="yes, within each signature", note="z-scored across tumours before clustering")])
AX = {"CD8 T cell": ["T cells CD8"], "NK cell": ["NK cells resting", "NK cells activated"], "Macrophage": ["Macrophages M0", "Macrophages M1", "Macrophages M2"], "Monocyte / DC": ["Monocytes", "Dendritic cells resting", "Dendritic cells activated"]}
PEER = {"CD8 T cell": [("xCell", "xcell", ["T cell CD8+"]), ("quanTIseq", "quantiseq", ["T cell CD8+"]), ("MCP-counter", "mcp_counter", ["T cell CD8+"]), ("EPIC", "epic", ["T cell CD8+"])],
        "NK cell": [("xCell", "xcell", ["NK cell"]), ("quanTIseq", "quantiseq", ["NK cell"]), ("MCP-counter", "mcp_counter", ["NK cell"]), ("EPIC", "epic", ["NK cell"])],
        "Macrophage": [("xCell", "xcell", ["Macrophage"]), ("quanTIseq", "quantiseq", ["Macrophage M1", "Macrophage M2"]), ("MCP-counter", "mcp_counter", ["Macrophage/Monocyte"]), ("EPIC", "epic", ["Macrophage"])],
        "Monocyte / DC": [("xCell", "xcell", ["Monocyte"]), ("quanTIseq", "quantiseq", ["Monocyte", "Myeloid dendritic cell"]), ("MCP-counter (mDC)", "mcp_counter", ["Myeloid dendritic cell"])]}
MAP = []
for ax, cols in AX.items():
    MAP.append(dict(axis=ax, method="CIBERSORTx (LM22)", features_summed=" + ".join(cols)))
    for lab, m, c in PEER[ax]: MAP.append(dict(axis=ax, method=lab, features_summed=" + ".join(c)))
MAP = pd.DataFrame(MAP)
def ax_val(method, ax, cols=None):
    if method == "CIBERSORTx (LM22)": return rel.loc[ids, AX[ax]].sum(1)
    m = [x for x in PEER[ax] if x[0] == method][0]; return D[m[1]].loc[ids, m[2]].sum(1)
# (2) Figure 2B-type concordance (between-tumour Spearman, scale-free) for n = 344
rows = []
for ax in AX:
    x = ax_val("CIBERSORTx (LM22)", ax)
    for lab, m, c in PEER[ax]:
        r = stats.spearmanr(x, D[m].loc[ids, c].sum(1)); rows.append(dict(axis=ax, peer=lab, rho=r.statistic, p=r.pvalue, n=len(ids)))
CON = pd.DataFrame(rows); CON["q_BH"] = bh(CON.p.values)
# (3) percentile ranks, state ordering, total immune
def jt(groups):   # Jonckheere-Terpstra, normal approximation with tie correction ignored (continuous data)
    U = 0; ns = [len(g) for g in groups]
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            a, b = np.asarray(groups[i]), np.asarray(groups[j]); U += (a[:, None] < b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum()
    N = sum(ns); mu = (N ** 2 - sum(n ** 2 for n in ns)) / 4; var = (N ** 2 * (2 * N + 3) - sum(n ** 2 * (2 * n + 3) for n in ns)) / 72
    z = (U - mu) / np.sqrt(var); return z, 2 * stats.norm.sf(abs(z))
def cliff(a, b): a, b = np.asarray(a), np.asarray(b); return ((a[:, None] > b[None, :]).mean() - (a[:, None] < b[None, :]).mean())
series = {}
for ax in AX:
    series[("CIBERSORTx (LM22)", ax)] = ax_val("CIBERSORTx (LM22)", ax)
    for lab, m, c in PEER[ax]: series[(lab, ax)] = D[m].loc[ids, c].sum(1)
series[("quanTIseq", "Total immune (1 − uncharacterized)")] = 1 - D["quantiseq"].loc[ids, "uncharacterized cell"]
series[("EPIC", "Total immune (immune types summed)")] = D["epic"].loc[ids, ["B cell", "T cell CD4+", "T cell CD8+", "Macrophage", "NK cell"]].sum(1)
series[("xCell", "Total immune (immune score)")] = D["xcell"].loc[ids, "immune score"]
series[("MCP-counter", "Total immune (sum of immune-type z-scores)")] = ((D["mcp_counter"].loc[ids, ["T cell", "NK cell", "B cell", "Monocyte", "Macrophage/Monocyte", "Myeloid dendritic cell", "Neutrophil"]] - D["mcp_counter"].loc[ids].mean()) / D["mcp_counter"].loc[ids].std())[["T cell", "NK cell", "B cell", "Monocyte", "Macrophage/Monocyte", "Myeloid dendritic cell", "Neutrophil"]].sum(1)
series[("CIBERSORTx (LM22)", "Total immune (absolute score)")] = cb.loc[ids, "Absolute score (sig.score)"]
ORD = []
for (meth, ax), v in series.items():
    pr = v.rank(pct=True); g = [pr[P.state == e] for e in ECO]
    z, pj = jt(g); kw = stats.kruskal(*g)
    ORD.append(dict(method=meth, axis=ax, median_pct_low=g[0].median(), median_pct_intermediate=g[1].median(), median_pct_high=g[2].median(),
                    JT_z_low_int_high=z, JT_p=pj, KW_p=kw.pvalue, cliff_high_vs_low=cliff(g[2], g[0]), cliff_high_vs_intermediate=cliff(g[2], g[1]),
                    spearman_with_immune_PC1=stats.spearmanr(v, pc1.loc[ids]).statistic))
ORD = pd.DataFrame(ORD); ORD["JT_q_BH"] = bh(ORD.JT_p.values)
PCT = pd.DataFrame({f"{m} | {a}": v.rank(pct=True) for (m, a), v in series.items()}); PCT["state"] = P.state
for n_, d in [("methods", METH), ("axis_mapping", MAP), ("concordance_344", CON), ("state_ordering_percentile", ORD)]: d.to_csv(O + f"R52_{n_}.tsv", sep="\t", index=False)
PCT.to_csv(O + "R52_percentile_ranks_per_tumour.tsv", sep="\t")
tot = ORD[ORD.axis.str.startswith("Total")]
S = dict(n=len(ids), concordance_q05=int((CON.q_BH < 0.05).sum()), concordance_n=len(CON),
         total_immune_rho_with_PC1={r.method: round(r.spearman_with_immune_PC1, 3) for r in tot.itertuples()},
         total_immune_JT_p={r.method: float(r.JT_p) for r in tot.itertuples()},
         axes_ordered_low_int_high_q05=int(((ORD.JT_q_BH < 0.05) & (ORD.JT_z_low_int_high > 0)).sum()), axes_n=len(ORD))
json.dump(S, open(O + "R52_summary.json", "w"), indent=1, default=float)
pd.set_option("display.width", 250); print(CON.round(3).to_string()); print(ORD[["method", "axis", "median_pct_low", "median_pct_intermediate", "median_pct_high", "JT_z_low_int_high", "JT_q_BH", "spearman_with_immune_PC1"]].round(3).to_string()); print(S)
