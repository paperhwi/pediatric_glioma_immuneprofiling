"""Figure 6 — molecular-group-adjusted differential expression and the
structure of the three contrasts.

A-C  volcano plots, reference group named in every title
D    overlap of the two Immune-desert-referenced contrasts, drawn as
     proportional bars rather than a Venn diagram
E    effect sizes of the two Immune-desert-referenced contrasts against
     each other, with an orthogonal-regression fit
F    paired effect sizes on the shared genes
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.linewidth": 0.8, "pdf.fonttype": 42,
                     "ps.fonttype": 42})

SRC = "/mnt/user-data/uploads/Open PBTA/Revision/Week1/_inputs/adjusted_DEG_all_contrasts.tsv"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "panels")
os.makedirs(OUT, exist_ok=True)

BLUE, RED, GREY, DARK = "#3B82F6", "#EF4444", "#CBD5E1", "#7C8798"
NAVY, MAROON, TEAL = "#1E3A8A", "#7F1D1D", "#0F766E"

A = "Lymphocyte_inflamed_vs_Immune_desert"
B = "Myeloid_dominant_vs_Immune_desert"
C = "Lymphocyte_inflamed_vs_Myeloid_dominant"

d = pd.read_csv(SRC, sep="\t")
w = d.pivot(index="gene", columns="contrast", values="adjusted_log2FC")
q = d.pivot(index="gene", columns="contrast", values="q_BH")
w = w[[A, B, C]].dropna()
q = q.loc[w.index, [A, B, C]]
NG = len(w)


def sets(c):
    up = (q[c] < 0.05) & (w[c] > 1)
    dn = (q[c] < 0.05) & (w[c] < -1)
    return up, dn


# ------------------------------------------------------------ A, B, C
panels = [(A, "Lymphocyte-inflamed  vs  Immune-desert", "Immune-desert", BLUE, "fig6_A"),
          (B, "Myeloid-dominant  vs  Immune-desert", "Immune-desert", RED, "fig6_B"),
          (C, "Lymphocyte-inflamed  vs  Myeloid-dominant", "Myeloid-dominant", BLUE, "fig6_C")]
counts = {}
for c, title, ref, col, name in panels:
    up, dn = sets(c)
    counts[c] = (int(up.sum()), int(dn.sum()))
    fig, a = plt.subplots(figsize=(3.75, 3.45), dpi=300)
    x = w[c].values
    y = -np.log10(np.clip(q[c].values, 1e-300, None))
    oth = ~(up | dn).values
    a.scatter(x[oth], y[oth], s=1.5, c=GREY, rasterized=True, lw=0)
    a.scatter(x[dn.values], y[dn.values], s=2.4, c=DARK, rasterized=True, lw=0)
    a.scatter(x[up.values], y[up.values], s=2.4, c=col, rasterized=True, lw=0)
    for v in (-1, 1):
        a.axvline(v, color="k", ls="--", lw=.6)
    a.axhline(-np.log10(.05), color="k", ls="--", lw=.6)
    a.set_title(f"{title}\nreference = {ref}   |   up {counts[c][0]:,}; down {counts[c][1]:,}",
                fontsize=7.8, pad=6)
    a.set_xlabel("Molecular-group-adjusted log$_2$ fold change", fontsize=7.8)
    a.set_ylabel("$-$log$_{10}$(BH $q$)", fontsize=7.8)
    a.spines[["top", "right"]].set_visible(False)
    a.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)

upA, _ = sets(A)
upB, _ = sets(B)
sA, sB = set(w.index[upA]), set(w.index[upB])
shared = sorted(sA & sB)
nS, nA, nB = len(shared), len(sA), len(sB)

# ------------------------------------------------------------ D  overlap
fig, a = plt.subplots(figsize=(3.75, 3.45), dpi=300)
ypos = [1.0, 0.0]
bh = 0.42
a.barh(ypos[0], nS, height=bh, color=NAVY, edgecolor="none")
a.barh(ypos[0], nA - nS, left=nS, height=bh, color=BLUE, alpha=.45, edgecolor="none")
a.barh(ypos[1], nS, height=bh, color=NAVY, edgecolor="none")
a.barh(ypos[1], nB - nS, left=nS, height=bh, color=RED, alpha=.75, edgecolor="none")
a.set_yticks(ypos)
a.set_yticklabels(["Lymphocyte-inflamed\nvs Immune-desert", "Myeloid-dominant\nvs Immune-desert"],
                  fontsize=7.4)
a.set_xlabel("Up-regulated genes (BH $q$ < 0.05, adj. log$_2$FC > 1)", fontsize=7.8)
a.set_xlim(0, nA * 1.30)
a.set_ylim(-0.62, 1.62)
a.text(nS / 2, ypos[0], f"{nS}", ha="center", va="center", fontsize=7.4,
       color="white", weight="bold")
a.text(nS + (nA - nS) / 2, ypos[0], f"{nA-nS:,} unique", ha="center", va="center",
       fontsize=7.4, color=NAVY, weight="bold")
a.text(nS / 2, ypos[1], f"{nS}", ha="center", va="center", fontsize=7.4,
       color="white", weight="bold")
a.annotate(f"{nB-nS} unique", xy=(nB, ypos[1] + bh / 2), xytext=(nB + nA * 0.09, ypos[1] + 0.46),
           fontsize=7.2, color=MAROON,
           arrowprops=dict(arrowstyle="->", color=MAROON, lw=0.8))
a.text(nA * 1.02, ypos[0], f"total {nA:,}", ha="left", va="center", fontsize=7.2, color=NAVY)
a.text(nB + nA * 0.02, ypos[1], f"total {nB}", ha="left", va="center",
       fontsize=7.2, color=MAROON)
a.set_title(f"{nS}/{nB} ({nS/nB*100:.1f}%) of the Myeloid-dominant set is\n"
            f"contained in the Lymphocyte-inflamed set",
            fontsize=7.8, loc="center", pad=6)
a.spines[["top", "right", "left"]].set_visible(False)
a.tick_params(axis="y", length=0)
a.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig(f"{OUT}/fig6_D.pdf", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------ E  gradient
fig, a = plt.subplots(figsize=(3.75, 3.45), dpi=300)
x, y = w[B].values, w[A].values
sig = np.array([g in sA or g in sB for g in w.index])
a.scatter(x[~sig], y[~sig], s=2, c=GREY, alpha=.35, lw=0, rasterized=True)
a.scatter(x[sig], y[sig], s=3, c=BLUE, alpha=.45, lw=0, rasterized=True)
lim = [-3.0, 5.0]
a.plot(lim, lim, "k--", lw=.8, label="identity ($y = x$)")
xc, yc = x - x.mean(), y - y.mean()
sxx, syy, sxy = (xc ** 2).mean(), (yc ** 2).mean(), (xc * yc).mean()
sl = ((syy - sxx) + np.sqrt((syy - sxx) ** 2 + 4 * sxy ** 2)) / (2 * sxy)
ic = y.mean() - sl * x.mean()
r = stats.pearsonr(x, y).statistic
xs = np.linspace(*lim, 10)
a.plot(xs, sl * xs + ic, color=TEAL, lw=1.3, label=f"orthogonal fit (slope = {sl:.2f})")
a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel("adj. log$_2$FC   Myeloid-dominant vs Immune-desert", fontsize=7.8)
a.set_ylabel("adj. log$_2$FC   Lymphocyte-inflamed vs Immune-desert", fontsize=7.8)
a.set_title(f"All {NG:,} genes; Pearson $r$ = {r:.2f}\n"
            f"slope above identity is consistent with an ordered\nshared immune-presence axis",
            fontsize=7.8, loc="left", pad=6)
a.legend(fontsize=6.4, frameon=False, loc="upper left")
a.spines[["top", "right"]].set_visible(False)
a.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig(f"{OUT}/fig6_E.pdf", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------ F  paired
fig, a = plt.subplots(figsize=(3.75, 3.45), dpi=300)
sub = w.loc[shared]
parts = a.violinplot([sub[B].values, sub[A].values], positions=[0, 1],
                     widths=.75, showextrema=False)
for pc, cc in zip(parts["bodies"], [RED, BLUE]):
    pc.set_facecolor(cc); pc.set_alpha(.45); pc.set_edgecolor(cc)
a.plot([0, 1], [sub[B].median(), sub[A].median()], "ko-", ms=4, lw=1.2)
a.set_xticks([0, 1])
a.set_xticklabels(["Myeloid-dominant\nvs Immune-desert", "Lymphocyte-inflamed\nvs Immune-desert"],
                  fontsize=7.4)
a.set_ylabel(f"adj. log$_2$FC on the {nS} shared genes", fontsize=7.8)
wil = stats.wilcoxon(sub[A].values, sub[B].values)
ratio = float(np.median(sub[A] / sub[B]))
a.set_title(f"Median {ratio:.2f}$\\times$ larger in Lymphocyte-inflamed\n"
            f"paired Wilcoxon signed-rank $P$ = {wil.pvalue:.1e}",
            fontsize=7.8, loc="left", pad=6)
a.spines[["top", "right"]].set_visible(False)
a.tick_params(labelsize=7)
fig.tight_layout()
fig.savefig(f"{OUT}/fig6_F.pdf", bbox_inches="tight")
plt.close(fig)

summary = pd.DataFrame({
    "set_name": ["Shared: LI-vs-Desert AND MD-vs-Desert",
                 "Unique: LI-vs-Desert only", "Unique: MD-vs-Desert only",
                 "Total LI-vs-Desert up", "Total MD-vs-Desert up"],
    "n": [nS, nA - nS, nB - nS, nA, nB]})
summary.to_csv(f"{OUT}/../tables/fig6_overlap_summary.tsv", sep="\t", index=False)

print("counts:", {k: v for k, v in counts.items()})
print(f"genes={NG}  shared={nS}  LIup={nA}  MDup={nB}  slope={sl:.3f}  r={r:.3f}  "
      f"ratio={ratio:.3f}  P={wil.pvalue:.3e}")
