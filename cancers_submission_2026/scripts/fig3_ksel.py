"""Figure 3A-E: selection of the number of clusters, k = 1-10.

Regenerated from the locked sweep tables so that the annotation wording matches
the Cancers submission:
  * consensus-CDF panel  -> "largest remaining relative dAUC (+29.8%)"
  * gap-statistic panel  -> one-standard-error rule selects k = 3
Each panel is written as its own vector PDF; panel letters are added by the
composition step so that they match the rest of the figure.
"""
import os
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.linewidth": 0.8, "xtick.major.width": 0.8,
                     "ytick.major.width": 0.8, "pdf.fonttype": 42,
                     "ps.fonttype": 42})

SRC = "/mnt/user-data/uploads/Open PBTA/Revision/REVISION PACKAGE 260916/3. Tables"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "panels")
os.makedirs(OUT, exist_ok=True)

sw = pd.read_csv(f"{SRC}/T1_ksweep_k2_k10.tsv", sep="\t")
gp = pd.read_csv(f"{SRC}/T1_gap_statistic.tsv", sep="\t")
BLUE, ORANGE, RED, GREY, TEAL = "#2563EB", "#F97316", "#B91C1C", "#64748B", "#0F766E"
W, H = 2.62, 2.80


def new():
    fig, a = plt.subplots(figsize=(W, H), dpi=300)
    return fig, a


def mark(a, title, ylab):
    a.axvline(2, color=GREY, ls=":", lw=0.8)
    a.axvline(3, color=RED, ls="--", lw=0.9)
    a.set_xlabel("Number of clusters (k)", fontsize=8)
    a.set_xticks(range(1, 11))
    a.set_ylabel(ylab, fontsize=8)
    a.set_title(title, fontsize=8.5, loc="left", pad=6)
    a.spines[["top", "right"]].set_visible(False)
    a.tick_params(labelsize=7)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# A ---------------------------------------------------------------- PAC
fig, a = new()
a.plot(sw.k, sw.PAC, "o-", color=BLUE, ms=3.5, lw=1.2)
mark(a, "Proportion of ambiguous clustering", "PAC (lower is better)")
save(fig, "fig3_A_pac")

# B --------------------------------------------------------- silhouette
fig, a = new()
a.plot(sw.k, sw.silhouette_consensus, "o-", color=ORANGE, ms=3.5, lw=1.2)
mark(a, "Silhouette width", "Silhouette on 1 − consensus")
save(fig, "fig3_B_silhouette")

# C --------------------------------------------------- consensus-CDF AUC
fig, a = new()
d = sw.dropna(subset=["delta_AUC"])
a.bar(d.k, d.delta_AUC, color=[RED if k == 3 else "#94A3B8" for k in d.k], width=.6)
mark(a, "Consensus-CDF elbow", "Relative $\\Delta$AUC of consensus CDF")
a.set_xlim(2.3, 10.7); a.set_xticks(range(3, 11))
a.annotate("largest remaining\nrelative $\\Delta$AUC (+29.8%)",
           xy=(3.35, sw.loc[sw.k == 3, "delta_AUC"].iloc[0] * 0.96),
           xytext=(5.0, 0.215), fontsize=6.6, color=RED,
           arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
save(fig, "fig3_C_deltaAUC")

# D ------------------------------------------------------- gap statistic
fig, a = new()
a.errorbar(gp.k, gp.gap, yerr=gp.s_k, fmt="o-", color=TEAL, ms=3.5, lw=1.2,
           capsize=2, elinewidth=0.8)
mark(a, "Gap statistic, k = 1–10", "Gap statistic")
a.annotate("one-standard-error\nrule selects k = 3",
           xy=(3, gp.loc[gp.k == 3, "gap"].iloc[0]), xytext=(4.7, 1.522),
           fontsize=6.6, color=RED,
           arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
save(fig, "fig3_D_gap")

# E ---------------------------------------------------- cluster-size floor
fig, a = new()
a.plot(sw.k, sw.min_cluster_size, "o-", color=RED, ms=3.5, lw=1.2)
a.axhline(10, color=GREY, ls="--", lw=0.8)
a.set_yscale("log")
mark(a, "Minimum cluster size", "Smallest cluster size (n)")
a.text(7.4, 12, "n = 10", fontsize=6.8, color=GREY)
save(fig, "fig3_E_size")

print("Figure 3A-E panels written")
