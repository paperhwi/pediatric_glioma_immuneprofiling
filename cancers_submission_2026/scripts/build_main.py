"""Assemble the seven main figures for the Cancers submission."""
import os
import shutil
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compose
from compose import Panel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = f"{ROOT}/panels"
M = f"{ROOT}/main"
SRC = ("/mnt/user-data/uploads/Open PBTA/Revision/REVISION PACKAGE 260916/"
       "1. Figure 300dpi/Main figures")
os.makedirs(M, exist_ok=True)


def build(name, rows, **kw):
    pdf = f"{M}/{name}.pdf"
    png = f"{M}/{name}.png"
    compose.compose(rows, pdf, png=png, **kw)
    compose.report(pdf, png)


# ---------------------------------------------------------------- Figure 1
shutil.copy(f"{P}/Figure1_full.pdf", f"{M}/Figure1.pdf")
compose.compose([[Panel(f"{P}/Figure1_full.pdf")]], f"{M}/Figure1.pdf",
                png=f"{M}/Figure1.png")
compose.report(f"{M}/Figure1.pdf", f"{M}/Figure1.png")

# ---------------------------------------------------------------- Figure 2
# panel B: the source panel leaves a wide band between the axis labels and its
# legend; recompose the two pieces with a normal gap.
compose.compose([[Panel(f"{SRC}/Figure2B.pdf", clip=(0, 0, 784.5, 352))],
                 [Panel(f"{SRC}/Figure2B.pdf", clip=(0, 398, 784.5, 420))]],
                f"{P}/fig2_B.pdf", row_gap=6, margin=2)

build("Figure2", [
    [Panel(f"{SRC}/Figure2A.pdf", "A",
           "CIBERSORTx LM22 deconvolution quality control (n = 349)")],
    [Panel(f"{P}/fig2_B.pdf", "B")],
    [Panel(f"{SRC}/Figure2D.pdf", "C",
           "Cross-method estimates on key immune axes by ecotype (n = 349)",
           clip=(0, 26, 1288.6, 340.6), scale=0.74)],
])

# ---------------------------------------------------------------- Figure 3
build("Figure3", [
    [Panel(f"{P}/fig3_A_pac.pdf", "A"), Panel(f"{P}/fig3_B_silhouette.pdf", "B"),
     Panel(f"{P}/fig3_C_deltaAUC.pdf", "C"), Panel(f"{P}/fig3_D_gap.pdf", "D"),
     Panel(f"{P}/fig3_E_size.pdf", "E")],
    [Panel(f"{SRC}/Figure3B.pdf", "F",
           "Principal-component and UMAP projections of the 46-feature matrix")],
    [Panel(f"{SRC}/Figure3C.pdf", "G", scale=0.82)],
], row_align="top")

# ---------------------------------------------------------------- Figure 4
build("Figure4", [[Panel(f"{SRC}/Figure4.pdf")]])

# ---------------------------------------------------------------- Figure 5
build("Figure5", [
    [Panel(f"{P}/fig5_A.pdf", "A"), Panel(f"{P}/fig5_B.pdf", "B")],
    [Panel(f"{P}/fig5_C.pdf", "C",
           "Spots mapped to the nearest bulk ecotype centroid, three pDMG sections (GSE268577)")],
    [Panel(f"{P}/fig5_D.pdf", "D"), Panel(f"{P}/fig5_E.pdf", "E")],
], row_align="top")

# ---------------------------------------------------------------- Figure 6
build("Figure6", [
    [Panel(f"{P}/fig6_A.pdf", "A"), Panel(f"{P}/fig6_B.pdf", "B"),
     Panel(f"{P}/fig6_C.pdf", "C")],
    [Panel(f"{P}/fig6_D.pdf", "D"), Panel(f"{P}/fig6_E.pdf", "E"),
     Panel(f"{P}/fig6_F.pdf", "F")],
], row_align="top")

# ---------------------------------------------------------------- Figure 7
build("Figure7", [
    [Panel(f"{SRC}/Figure7A.pdf", "A")],
    [Panel(f"{SRC}/Figure7C.pdf", "B",
           "Multivariable Cox proportional-hazards model, Immune-desert reference (n = 251, 208 events)",
           clip=(0, 34, 652.9, 375.3))],
])
