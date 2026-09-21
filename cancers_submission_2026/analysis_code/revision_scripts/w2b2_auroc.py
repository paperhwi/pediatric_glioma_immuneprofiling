"""A (continued) — AUROC of each signature for identifying its named cell type,
computed within each sample and then summarised, so it cannot be driven by one patient."""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
S=pd.read_pickle("w2_sigscores_10x.pkl")
lab=S.pop("cell_type"); samp=S.pop("sampleid")
EXPECT={"Microglia_Core_Homeostatic":["Myeloid"],"Microglia_Klemm2020":["Myeloid"],"MDM_Klemm2020":["Myeloid"],
 "MgTAM_Antunes2021":["Myeloid"],"MoTAM_Antunes2021":["Myeloid"],"DAM_KerenShaul2017":["Myeloid"],
 "M1_Macrophage":["Myeloid"],"M2_Macrophage":["Myeloid"],"Dendritic_Cell_Activation":["Myeloid"],
 "Neutrophil_Activation":["Myeloid"],"Glioma_Inflammatory_Wang2017":["Myeloid"],"MHC_Class_II":["Myeloid"],
 "T_Cell_Cytotoxicity":["CD8"],"NK_Cell_Activity":["CD8"],"Tregs_Friebel2020":["CD4"],
 "T_Cell_Exhaustion":["CD4","CD8"]}
rows=[]
for sig,tgt in EXPECT.items():
    y=lab.isin(tgt).astype(int)
    overall=roc_auc_score(y,S[sig])
    per=[]
    for sid,idx in samp.groupby(samp).groups.items():
        yy=y.loc[idx]
        if yy.nunique()==2: per.append(roc_auc_score(yy,S.loc[idx,sig]))
    rows.append(dict(signature=sig,expected_cell_type="/".join(tgt),AUROC_all_cells=round(overall,3),
                     AUROC_min_sample=round(min(per),3),AUROC_max_sample=round(max(per),3),
                     n_samples=len(per)))
r=pd.DataFrame(rows).sort_values("AUROC_all_cells",ascending=False)
print("=== AUROC for identifying the named cell type (18,619 cells, 4 pHGG samples) ===")
print(r.to_string(index=False))
print(f"\nsignatures with AUROC >= 0.70 in every sample: "
      f"{int((r.AUROC_min_sample>=0.70).sum())}/{len(r)}")
print(f"signatures with AUROC >= 0.60 in every sample: {int((r.AUROC_min_sample>=0.60).sum())}/{len(r)}")
print(f"median AUROC across the 16 directional signatures: {r.AUROC_all_cells.median():.3f}")
r.to_csv("W2_A_auroc.tsv",sep="\t",index=False)
