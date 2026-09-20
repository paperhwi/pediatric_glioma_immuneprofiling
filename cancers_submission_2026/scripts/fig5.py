"""Figure 5 — single-cell and spatial support for the cellular interpretation
of the immune signatures and of the ecotype axis.

A  AUROC of each directional signature for the cell type it is named for
B  cell type expressing the genes that define the ecotype axis
C  built separately (fig5_spatial.py): spot-to-centroid assignment, 3 sections
D  spatial coherence of the spot assignments against a permutation null
E  Moran's I for the most spatially structured signatures

Wording is deliberately conservative: the signatures and the ecotype axis are
described as immune-cell-associated, with many genes preferentially expressed
in myeloid cells, and the spatial result is described in terms of spots that
map most closely to each bulk ecotype centroid.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.linewidth": 0.8, "pdf.fonttype": 42,
                     "ps.fonttype": 42})

SRC = "/mnt/user-data/uploads/Open PBTA/Revision/REVISION PACKAGE 260916/3. Tables"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "panels")
os.makedirs(OUT, exist_ok=True)

MYE, LYM, DES = "#EF4444", "#3B82F6", "#9CA3AF"
TEAL, PURPLE, RED, GREY = "#0F766E", "#6D28D9", "#B91C1C", "#64748B"


def pretty(s):
    return s.replace("_", " ")


# ------------------------------------------------------------------ A
t = pd.read_csv(f"{SRC}/T_S25_signature_attribution.tsv", sep="\t")
dr = t.dropna(subset=["AUROC_all_cells"]).copy()
dr["lineage"] = np.where(dr["expected"].str.contains("Myeloid"), "Myeloid",
                         "T / NK cell")
dr = dr.sort_values("AUROC_all_cells", ascending=True)
col = [RED if v == "DISCORDANT" else (TEAL if lg == "Myeloid" else PURPLE)
       for v, lg in zip(dr["verdict"], dr["lineage"])]

fig, a = plt.subplots(figsize=(5.95, 3.70), dpi=300)
ypos = np.arange(len(dr))
a.barh(ypos, dr["AUROC_all_cells"], color=col, height=0.72, edgecolor="none")
a.errorbar(dr["AUROC_all_cells"], ypos,
           xerr=[np.clip(dr["AUROC_all_cells"] - dr["AUROC_min_sample"], 0, None),
                 np.clip(dr["AUROC_max_sample"] - dr["AUROC_all_cells"], 0, None)],
           fmt="none", ecolor="#334155", elinewidth=0.8, capsize=1.8)
a.axvline(0.5, color="k", ls="--", lw=0.7)
a.set_yticks(ypos)
a.set_yticklabels([pretty(s) for s in dr["signature"]], fontsize=7)
a.set_xlim(0, 1.04)
a.set_xlabel("AUROC for the cell type the signature is named for", fontsize=7.8)
med = dr["AUROC_all_cells"].median()
a.set_title("Signature scores resolve the cell types they are named for\n"
            f"18,619 annotated immune cells, 4 pediatric high-grade gliomas (GSE227983); median AUROC {med:.3f}",
            fontsize=7.8, loc="left", pad=6)
a.text(dr["AUROC_all_cells"].iloc[0] + 0.02, 0, f"{dr['AUROC_all_cells'].iloc[0]:.3f}",
       va="center", ha="left", fontsize=6.8, color=RED)
a.legend(handles=[Line2D([], [], marker="s", ls="", ms=6, color=TEAL,
                         label="myeloid signature"),
                  Line2D([], [], marker="s", ls="", ms=6, color=PURPLE,
                         label="T / NK-cell signature"),
                  Line2D([], [], marker="s", ls="", ms=6, color=RED,
                         label="discordant \u2014 identifies T cells")],
         fontsize=6.8, frameon=False, loc="upper center", handletextpad=0.4,
         ncol=3, bbox_to_anchor=(0.5, -0.155), columnspacing=1.4)
a.spines[["top", "right"]].set_visible(False)
a.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig(f"{OUT}/fig5_A.pdf", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ B
g = pd.read_csv(f"{SRC}/T_S26_ecotype_gene_source.tsv", sep="\t")
order = ["Myeloid", "Tcell", "Malignant"]
lab = {"Myeloid": "Myeloid", "Tcell": "T cell", "Malignant": "Malignant"}
cnt = g["top_cell_type"].value_counts().reindex(order).fillna(0).astype(int)
n_tot = int(cnt.sum())
n_imm = int(cnt["Myeloid"] + cnt["Tcell"])
pct = 100 * n_imm / n_tot
BG, BG_LO, BG_HI, PPERM = 72.6, 67.5, 77.2, 0.005

fig, axs = plt.subplots(1, 2, figsize=(5.95, 3.70), dpi=300,
                        gridspec_kw={"width_ratios": [1.45, 1.0]})
a = axs[0]
bars = a.bar([lab[k] for k in order], cnt.values,
             color=[TEAL, PURPLE, "#94A3B8"], width=0.62)
for r, v in zip(bars, cnt.values):
    a.text(r.get_x() + r.get_width() / 2, v + n_tot * 0.02, f"{v}",
           ha="center", va="bottom", fontsize=7.6, weight="bold")
a.set_ylabel("Ecotype-axis genes", fontsize=7.8)
a.set_ylim(0, cnt.max() * 1.20)
a.set_title("Cell type with the highest mean expression\n"
            f"of each of the {n_tot} detectable ecotype-axis genes",
            fontsize=7.8, loc="left", pad=6)
a.spines[["top", "right"]].set_visible(False)
a.tick_params(labelsize=7)

a = axs[1]
a.bar([0], [pct], width=0.5, color="#0F766E")
a.bar([1], [BG], width=0.5, color="#94A3B8")
a.errorbar([1], [BG], yerr=[[BG - BG_LO], [BG_HI - BG]], fmt="none",
           ecolor="#334155", elinewidth=0.9, capsize=3)
a.text(0, pct + 2, f"{pct:.1f}%", ha="center", fontsize=7.6, weight="bold",
       color="#0F766E")
a.text(1, BG_HI + 2, f"{BG:.1f}%", ha="center", fontsize=7.6, color="#334155")
a.set_xticks([0, 1])
a.set_xticklabels(["ecotype-axis\ngenes", "expression-matched\nrandom genes"],
                  fontsize=7)
a.set_ylim(0, 112)
a.set_ylabel("Attributed to immune cells (%)", fontsize=7.8)
a.set_title(f"Permutation $P$ = {PPERM:.3f}", fontsize=7.8, loc="left", pad=6)
a.spines[["top", "right"]].set_visible(False)
a.tick_params(labelsize=7)
fig.text(0.5, -0.045,
         "Immune-cell-associated, with many genes preferentially expressed in myeloid cells",
         ha="center", fontsize=7.2, color="#334155", style="italic")
fig.tight_layout()
fig.savefig(f"{OUT}/fig5_B.pdf", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ D
sp = pd.read_csv(f"{SRC}/T_S28_spot_ecotype.tsv", sep="\t")
sp["label"] = [f"Sample-{int(s)}" for s in sp["sample"]]
fig, a = plt.subplots(figsize=(5.75, 3.25), dpi=300)
x = np.arange(len(sp))
wd = 0.34
a.bar(x - wd / 2, sp["same_neighbour_fraction"], width=wd, color="#0F766E",
      label="observed")
a.bar(x + wd / 2, sp["null_mean"], width=wd, color="#CBD5E1",
      label="permutation null (999 label permutations)")
for i, row in sp.iterrows():
    a.text(i - wd / 2, row["same_neighbour_fraction"] + 0.012,
           f"{row['same_neighbour_fraction']:.3f}", ha="center", fontsize=7,
           weight="bold", color="#0F766E")
    a.text(i + wd / 2, row["null_mean"] + 0.012, f"{row['null_mean']:.3f}",
           ha="center", fontsize=7, color="#475569")
    a.text(i, max(row["same_neighbour_fraction"], row["null_mean"]) + 0.05,
           f"$P$ = {row['perm_p']:.3f}", ha="center", fontsize=7, color="#0F172A")
a.set_xticks(x)
a.set_xticklabels([f"{l}\n{int(n)} in-tissue spots" for l, n in
                   zip(sp["label"], sp["n_spots"])], fontsize=7.2)
a.set_ylabel("Fraction of neighboring spot pairs\nwith the same assignment", fontsize=7.8)
a.set_ylim(0, 0.78)
a.set_title("Spot-to-centroid assignments form spatially coherent domains",
            fontsize=7.8, loc="left", pad=6)
a.legend(fontsize=6.8, frameon=False, loc="upper right")
a.spines[["top", "right"]].set_visible(False)
a.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig(f"{OUT}/fig5_D.pdf", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ E
mo = pd.read_csv(f"{SRC}/T_S27_spatial_autocorrelation.tsv", sep="\t")
mean_i = mo.groupby("signature")["morans_I"].mean().sort_values()
top = mean_i.tail(12).index.tolist()
sub = mo[mo["signature"].isin(top)].copy()
marks = {"Sample-1": ("o", "#1E3A8A"), "Sample-2": ("s", "#F97316"),
         "Sample-3": ("^", "#0F766E")}
fig, a = plt.subplots(figsize=(6.15, 3.25), dpi=300)
ypos = {s: i for i, s in enumerate(top)}
for s, (m, c) in marks.items():
    d = sub[sub["sample"] == s]
    a.scatter(d["morans_I"], [ypos[g] for g in d["signature"]], marker=m, s=26,
              color=c, label=s, lw=0)
a.axvline(0, color="k", lw=0.8)
a.set_yticks(range(len(top)))
a.set_yticklabels([pretty(s) for s in top], fontsize=7)
a.set_xlabel("Moran's I (999-permutation test)", fontsize=7.8)
a.set_title("Immune programs are spatially structured, not randomly distributed\n"
            "top 12 of 24 signatures by mean Moran's I; 10 of 24 reach FDR < 0.05 in all three sections",
            fontsize=7.8, loc="left", pad=6)
a.legend(fontsize=6.8, frameon=False, loc="lower right")
a.spines[["top", "right"]].set_visible(False)
a.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig(f"{OUT}/fig5_E.pdf", bbox_inches="tight")
plt.close(fig)

print(f"Figure 5 panels A,B,D,E written  "
      f"(genes {n_tot}, immune {n_imm} = {pct:.1f}%, median AUROC {med:.3f})")
