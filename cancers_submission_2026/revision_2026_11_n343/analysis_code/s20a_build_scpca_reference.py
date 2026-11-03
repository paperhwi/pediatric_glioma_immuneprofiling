"""ScPCA SCPCP000001 single-cell reference for bulk deconvolution (BayesPrism / InstaPrism / MuSiC).
19 pediatric HGG donors (excluded: normal margin SCPCS000011 and the PXA-lineage SCPCS000012, 13, 23, as in the ScPCA validation).
Cell types from lineage marker scores on ScPCA log-normalised counts, cross-checked against the ScPCA consensus/SingleR/CellAssign
labels; inferCNV (where it ran) moves CNV-high 'oligodendrocyte' cells into the malignant compartment; scDblFinder doublets removed.
cell.type (coarse) = Malignant, Microglia, Macrophage/monocyte (MDM), T cell, NK cell, B/plasma, Dendritic cell, Oligodendrocyte,
Endothelial, Mural (pericyte/SMC). cell.state = cell type, except Malignant = Malignant_<donor> (BayesPrism tumour states)."""
import anndata as ad, pandas as pd, numpy as np, scipy.sparse as sp, scipy.io as sio, glob, json, os
D = "/mnt/user-data/uploads/Pediatric GBM/SCPCP000001_extracted/SCPCP000001_single-cell/"
O = "/home/claude/decon/ref/"; os.makedirs(O, exist_ok=True)
EXCL = {"SCPCS000011": "Normal margin", "SCPCS000012": "PXA lineage", "SCPCS000013": "PXA lineage", "SCPCS000023": "PXA lineage"}
MK = {"Myeloid": ["C1QA", "C1QB", "C1QC", "TYROBP", "AIF1", "CSF1R", "CD14", "FCER1G", "CD68", "LAPTM5"],
      "T": ["CD3D", "CD3E", "CD2", "CD247", "TRAC", "IL7R"], "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1"],
      "B": ["MS4A1", "CD79A", "CD79B", "IGHM", "JCHAIN", "MZB1", "IGKC"], "DC": ["CLEC10A", "CD1C", "FCER1A", "CLEC9A", "LILRA4", "IRF8"],
      "Endothelial": ["PECAM1", "CLDN5", "VWF", "FLT1", "CDH5", "ESAM"], "Mural": ["PDGFRB", "RGS5", "ACTA2", "NOTCH3", "MYH11", "KCNJ8"],
      "Oligodendrocyte": ["MOG", "MOBP", "MAG", "PLP1", "MBP", "ERMN", "OPALIN"]}
MG = ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "GPR34", "P2RY13", "SIGLEC8", "CSF3R"][:6]
MDM = ["CD163", "F13A1", "MRC1", "S100A8", "S100A9", "VCAN", "FCN1", "LYZ", "TGFBI", "CCR2", "EREG"]
X_parts, meta, genes = [], [], None
for f in sorted(glob.glob(D + "*/*.h5ad")):
    A = ad.read_h5ad(f); s = A.obs.sample_id.iloc[0]
    if s in EXCL: continue
    sym = A.var.gene_symbol.astype(str).where(A.var.gene_symbol.notna(), A.var_names)
    keep_g = ~sym.duplicated(keep="first").values
    L = A.X.tocsc()[:, keep_g] if sp.issparse(A.X) else sp.csc_matrix(A.X[:, keep_g])
    Cn = A.raw.X.tocsc()[:, keep_g]; sy = sym.values[keep_g]; gi = {g: i for i, g in enumerate(sy)}
    def score(gs):
        idx = [gi[g] for g in gs if g in gi]
        return np.asarray(L[:, idx].mean(1)).ravel() if idx else np.zeros(A.n_obs)
    S = pd.DataFrame({k: score(v) for k, v in MK.items()}, index=A.obs_names)
    best = S.idxmax(1); top = S.max(1); second = S.apply(lambda r: r.nlargest(2).iloc[1], axis=1)
    ptprc = np.asarray(L[:, gi["PTPRC"]].todense()).ravel() > 0 if "PTPRC" in gi else np.zeros(A.n_obs, bool)
    ct = pd.Series("Malignant", index=A.obs_names)
    clear = (top > 0.8) & (top - second > 0.3)
    imm = best.isin(["Myeloid", "T", "NK", "B", "DC"])
    ct[clear & imm & (ptprc | (best == "Myeloid"))] = best[clear & imm & (ptprc | (best == "Myeloid"))]
    ct[clear & best.isin(["Endothelial", "Mural"])] = best[clear & best.isin(["Endothelial", "Mural"])]
    olig = clear & (best == "Oligodendrocyte")
    cnv = A.obs.get("infercnv_total_cnv")
    if cnv is not None and cnv.notna().any():
        thr = cnv[ct.isin(["Myeloid", "T"])].quantile(0.95) if (ct.isin(["Myeloid", "T"])).sum() > 20 else cnv.quantile(0.25)
        olig = olig & ~(cnv > thr).values          # CNV-high "oligodendrocytes" stay malignant
    ct[olig] = "Oligodendrocyte"
    my = ct == "Myeloid"; mgs, mds = score(MG), score(MDM)
    ct[my & (mgs >= mds)] = "Microglia"; ct[my & (mgs < mds)] = "Macrophage/monocyte"
    ct = ct.replace({"T": "T cell", "NK": "NK cell", "B": "B/plasma", "DC": "Dendritic cell"})
    dbl = (A.obs.scDblFinder_class.astype(str) == "doublet").values
    m = pd.DataFrame({"cell": [f"{s}_{c}" for c in A.obs_names], "donor": s, "cell_type": ct.values, "doublet": dbl,
                      "consensus": A.obs.consensus_celltype_annotation.astype(str).values, "singler": A.obs.singler_celltype_annotation.astype(str).values,
                      "cnv": (cnv.values if cnv is not None else np.nan), "umi": np.asarray(Cn.sum(1)).ravel(), "diagnosis": A.obs.diagnosis.astype(str).values,
                      "location": A.obs.tissue_location.astype(str).values})
    keep = ~dbl
    X_parts.append((Cn[keep.nonzero()[0], :], sy)); meta.append(m[keep])
# common genes
common = sorted(set.intersection(*[set(g) for _, g in X_parts]))
Xs = []
for Cn, sy in X_parts:
    gi = {g: i for i, g in enumerate(sy)}; Xs.append(Cn[:, [gi[g] for g in common]].tocsr())
X = sp.vstack(Xs).tocsr(); M = pd.concat(meta, ignore_index=True)
M["cell_state"] = np.where(M.cell_type == "Malignant", "Malignant_" + M.donor, M.cell_type)
# cap malignant cells per donor (300) for tractability; keep all non-malignant cells
rng = np.random.default_rng(42); keep = np.ones(len(M), bool)
for d, idx in M[M.cell_type == "Malignant"].groupby("donor").groups.items():
    idx = np.array(list(idx))
    if len(idx) > 300: keep[rng.choice(idx, len(idx) - 300, replace=False)] = False
M["in_reference"] = keep
tab = pd.crosstab(M.donor, M.cell_type); tab.to_csv(O + "cell_counts_by_donor.tsv", sep="\t")
ct_cons = pd.crosstab(M.cell_type, M.consensus); ct_cons.to_csv(O + "celltype_vs_scpca_consensus.tsv", sep="\t")
Mi = M[M.in_reference].reset_index(drop=True); Xi = X[np.where(M.in_reference)[0], :]
sio.mmwrite(O + "counts_cells_x_genes.mtx", Xi.tocoo()); pd.Series(common).to_csv(O + "genes.tsv", sep="\t", index=False, header=False)
Mi.to_csv(O + "cells.tsv", sep="\t", index=False); M.to_csv(O + "cells_all.tsv", sep="\t", index=False)
# marker means by type (sanity)
lg = np.log1p(Xi.multiply(1e4 / np.maximum(Xi.sum(1), 1)).tocsr()); gi = {g: i for i, g in enumerate(common)}
chk = {}
for t in sorted(Mi.cell_type.unique()):
    r = (Mi.cell_type == t).values
    chk[t] = {g: float(lg[r][:, gi[g]].mean()) for g in ["PTPRC", "P2RY12", "CD163", "CD3E", "NKG7", "MS4A1", "PECAM1", "PDGFRB", "MOG", "SOX2", "OLIG2", "PDGFRA", "EGFR"] if g in gi}
pd.DataFrame(chk).T.round(2).to_csv(O + "marker_means_by_type.tsv", sep="\t")
print(tab.to_string()); print(pd.DataFrame(chk).T.round(2).to_string()); print("reference cells", len(Mi), "genes", len(common))
json.dump(dict(n_cells_all=int(len(M)), n_cells_reference=int(len(Mi)), n_genes=len(common), donors=sorted(Mi.donor.unique().tolist()),
               type_counts=Mi.cell_type.value_counts().to_dict()), open(O + "reference_summary.json", "w"), indent=1)
