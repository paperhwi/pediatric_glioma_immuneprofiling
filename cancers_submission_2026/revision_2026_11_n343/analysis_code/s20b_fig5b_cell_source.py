"""Gap G1 — submitted Fig 4B (revised Fig 5B): "cell type with the highest mean expression of each ecotype-defining gene".
The submitted analysis (repository analysis_code/revision_scripts/w2c_source.py) used Smart-seq2 cells of GSE227983 with only
three classes (sorted Myeloid, sorted T cells, malignant), the top 300 Lymphocyte-inflamed vs Immune-desert adjusted DEGs,
and expression-matched random genes (20 expression bins, 200 draws). With three classes of which two are immune, the random
expectation is necessarily high (72.6%). Here the same definition is applied to the unsorted ScPCA SCPCP000001 10x data, first
with the same three classes, then adding the non-malignant brain cells the reviewer asked about (oligodendrocytes, vascular cells)."""
import sys, json, os, glob; sys.path.insert(0, "/home/claude/q343/src"); from common import ROOT
import numpy as np, pandas as pd
from scipy import io as sio
O = ROOT+"out/fig5b/"; os.makedirs(O, exist_ok=True); REF = "/home/claude/decon/ref/"
genes = pd.read_csv(REF + "genes.tsv", header=None)[0].values; cells = pd.read_csv(REF + "cells.tsv", sep="\t")
X = sio.mmread(REF + "counts_cells_x_genes.mtx").tocsr()
cpm = X.multiply(1e6 / np.asarray(X.sum(1)).ravel()[:, None]).tocsr(); cpm.data = np.log2(cpm.data / 10 + 1)   # log2(CPM/10 + 1), as in the submitted code
ty = cells.cell_type.values
TYPES = sorted(set(ty)); mean_t = pd.DataFrame({t: np.asarray(cpm[ty == t].mean(0)).ravel() for t in TYPES}, index=genes)
det_t = pd.DataFrame({t: np.asarray((cpm[ty == t] > 0).mean(0)).ravel() for t in TYPES}, index=genes)
MY = ["Microglia", "Macrophage/monocyte", "Dendritic cell"]; LY = ["T cell", "NK cell"]
DESIGNS = {
 "D1 submitted design: myeloid, T/NK, malignant": {"Myeloid": MY, "T/NK": LY, "Malignant": ["Malignant"]},
 "D2 + oligodendrocytes": {"Myeloid": MY, "T/NK": LY, "Malignant": ["Malignant"], "Oligodendrocyte": ["Oligodendrocyte"]},
 "D3 + oligodendrocytes + vascular (endothelial, mural)": {"Myeloid": MY, "T/NK": LY, "Malignant": ["Malignant"], "Oligodendrocyte": ["Oligodendrocyte"], "Vascular": ["Endothelial", "Mural"]},
 "D4 all ten reference types (immune = any immune type)": {t: [t] for t in TYPES}}
IMMUNE = set(MY + LY + ["B/plasma", "Myeloid", "T/NK"])
def class_means(design):   # mean over the cells of each class (pooled cells, as the submitted code pooled cells of a class)
    out = {}
    for k, members in design.items():
        idx = np.isin(ty, members); out[k] = np.asarray(cpm[idx].mean(0)).ravel()
    return pd.DataFrame(out, index=genes), np.isin(ty, sum(design.values(), []))
SETS = {}
G0 = pd.read_excel("/home/claude/w/pg/Revision_round1_Q1-Q5/Review_gaps_2Nov/R6_gap_analyses_results.xlsx", "G1_fig4B_gene_level")
SETS["previous top 300 (n = 343, n = 344 labels; reproduction check)"] = G0.loc[G0.gene_set == "revised (n = 343, PXA excluded)", "gene"].tolist()
d = pd.read_csv(ROOT + "out/T_adjusted_DEG_n343_all_contrasts.tsv", sep="\t"); d = d[(d.contrast == "Immune-high vs Immune-low") & (d.q_BH < 0.05) & (d.adjusted_log2FC > 1)].nsmallest(300, "q_BH")
SETS["top 300 Immune-high vs Immune-low (n = 343, final labels)"] = d.gene.tolist()
sh = glob.glob(ROOT + "out/T_shared_genes_*.tsv")[0]; SETS[f"{len(pd.read_csv(sh, sep=chr(9)))} shared immune-presence genes (Figure 6)"] = pd.read_csv(sh, sep="\t").gene.tolist()
rng = np.random.default_rng(42); ROWS = []; GENE = []
for sl, gl in SETS.items():
    for dl, design in DESIGNS.items():
        M, used = class_means(design)
        allmean = np.asarray(cpm[used].mean(0)).ravel(); am = pd.Series(allmean, index=genes)
        expressed = am[am > 0].index
        up = [g for g in dict.fromkeys(gl) if g in expressed]
        src = M.loc[up].idxmax(1); obs = src.isin(IMMUNE).mean()
        bins = pd.qcut(am.loc[expressed].rank(method="first"), 20, labels=False)
        pool = {b: bins.index[bins.values == b].values for b in range(20)}
        bg = []
        for _ in range(200):
            pick = [rng.choice(pool[bins[g]]) for g in up]; bg.append(M.loc[pick].idxmax(1).isin(IMMUNE).mean())
        bg = np.array(bg)
        ROWS.append(dict(gene_set=sl, design=dl, n_classes=len(design), n_genes_detectable=len(up), n_immune_top=int(src.isin(IMMUNE).sum()), pct_immune_top=100 * obs,
                         random_mean_pct=100 * bg.mean(), random_lo=100 * np.percentile(bg, 2.5), random_hi=100 * np.percentile(bg, 97.5),
                         perm_p=max((bg >= obs).mean(), 1 / len(bg)), excess_pct_points=100 * (obs - bg.mean()),
                         top_counts="; ".join(f"{k} {v}" for k, v in src.value_counts().items())))
        if dl.startswith("D3"):
            GENE.append(pd.DataFrame({"gene_set": sl, "gene": up, "top_class_D1": class_means(DESIGNS["D1 submitted design: myeloid, T/NK, malignant"])[0].loc[up].idxmax(1).values,
                                      "top_class_D3": src.values}).join(mean_t.loc[up].add_prefix("mean_").reset_index(drop=True)).join(det_t.loc[up].add_prefix("detect_").reset_index(drop=True)))
R = pd.DataFrame(ROWS); G = pd.concat(GENE)
R.to_csv(O + "T_fig5B_designs.tsv", sep="\t", index=False); G.to_csv(O + "T_fig5B_gene_level.tsv", sep="\t", index=False)
# genes that move from an immune class (D1) to a non-immune class (D3), revised gene set
g = G[G.gene_set.str.startswith("top 300")]
moved = g[g.top_class_D1.isin(IMMUNE) & ~g.top_class_D3.isin(IMMUNE)]
S = dict(designs=R.to_dict("records"), moved_to_non_immune_revised=moved[["gene", "top_class_D1", "top_class_D3"]].to_dict("records"),
         submitted_SS2_reported=dict(pct=98.4, n=246, random_mean=72.6, random_lo=67.5, random_hi=77.2, perm_p=0.005))
json.dump(S, open(O + "T_fig5B_summary.json", "w"), indent=1, default=float)
print(R[["gene_set", "design", "n_genes_detectable", "pct_immune_top", "random_mean_pct", "random_lo", "random_hi", "perm_p", "excess_pct_points"]].round(1).to_string())
print(len(moved), moved.gene.tolist()[:40])
