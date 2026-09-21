"""Regenerate Dunn signed -log10(q) heatmap for focus features with readable spacing."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_step6"

wide = pd.read_csv(OUT / "step6_Dunn_signed_log10q_wide.tsv", sep="\t", index_col=0)
print("pairs:", wide.columns.tolist())

focus = [
    ("LM22_T cells CD8",            "LM22  CD8 T cells"),
    ("LM22_Macrophages M0",         "LM22  Macrophage M0"),
    ("LM22_Macrophages M1",         "LM22  Macrophage M1"),
    ("LM22_Macrophages M2",         "LM22  Macrophage M2"),
    ("LM22_Monocytes",              "LM22  Monocytes"),
    ("LM22_Neutrophils",            "LM22  Neutrophils"),
    ("ssGSEA_Microglia_Klemm2020",  "ssGSEA  Microglia (Klemm 2020)"),
    ("ssGSEA_MDM_Klemm2020",        "ssGSEA  MDM (Klemm 2020)"),
    ("ssGSEA_MgTAM_Antunes2021",    "ssGSEA  Mg-TAM (Antunes 2021)"),
    ("ssGSEA_MoTAM_Antunes2021",    "ssGSEA  Mo-TAM (Antunes 2021)"),
    ("ssGSEA_T_Cell_Cytotoxicity",  "ssGSEA  T-cell cytotoxicity"),
    ("ssGSEA_IFN_Gamma_Response",   "ssGSEA  IFN-γ response"),
]
feat  = [c for c, _ in focus if c in wide.index]
label = [l for c, l in focus if c in wide.index]

# Use desired pair order if present
preferred = [
    "Immune-desert vs Lymphocyte-inflamed",
    "Immune-desert vs Myeloid-dominant",
    "Lymphocyte-inflamed vs Myeloid-dominant",
]
cols = [p for p in preferred if p in wide.columns] or list(wide.columns)
mat = wide.loc[feat, cols].copy()
mat.index = label

# Better sizing: 4.5 inches per pair column, 0.65 inch per feature row, big fonts
fig, ax = plt.subplots(figsize=(max(8, 3.6*len(cols)), 1.0 + 0.65*len(label)))
vmax = float(np.nanmax(np.abs(mat.values)))
vmax = max(vmax, 5)
sns.heatmap(mat,
            annot=True, fmt=".1f",
            annot_kws={"size": 11},
            cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
            linewidths=0.5, linecolor="white",
            cbar_kws=dict(label="signed -log10(q)  (Dunn, BH per feature)", shrink=0.8),
            ax=ax)
ax.set_title("Dunn post-hoc — focus signatures\n(positive = first group higher; negative = second group higher)",
             fontsize=12, pad=12)
ax.set_xlabel("")
ax.set_ylabel("")
plt.setp(ax.get_xticklabels(), rotation=18, ha="right", fontsize=10)
plt.setp(ax.get_yticklabels(), rotation=0, fontsize=10)
plt.tight_layout()
fig.savefig(FIG / "Dunn_signed_log10q_focus.png", dpi=220, bbox_inches="tight")
plt.close()
print("saved:", FIG / "Dunn_signed_log10q_focus.png")
