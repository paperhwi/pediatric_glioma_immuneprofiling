"""Assemble the renumbered supplementary figures S1-S15 for Cancers."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compose
from compose import Panel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = f"{ROOT}/panels"
S = f"{ROOT}/suppl"
B = "/mnt/user-data/uploads/Open PBTA/Revision/REVISION PACKAGE 260916/1. Figure 300dpi"
MF, SF = f"{B}/Main figures", f"{B}/Supplementary figures"
os.makedirs(S, exist_ok=True)


def build(name, rows, **kw):
    pdf, png = f"{S}/{name}.pdf", f"{S}/{name}.png"
    compose.compose(rows, pdf, png=png, **kw)
    compose.report(pdf, png)


# S1  brain-tuned myeloid signatures versus the LM22 macrophage/monocyte axes
build("FigureS1", [[Panel(f"{MF}/Figure2C.pdf")]])

# S2  LM22 quality control in the same-resource pooled set
build("FigureS2", [[Panel(f"{SF}/FigureS4.pdf")]])

# S3  ssGSEA TPM-harmonization concordance
build("FigureS3", [[Panel(f"{SF}/FigureS3.pdf")]])

# S4  re-clustering of the high-confidence LM22 subset
build("FigureS4", [[Panel(f"{SF}/FigureS6.pdf")]])

# S5  clustering-definition and held-out-feature sensitivity
build("FigureS5", [[Panel(f"{P}/suppl_clustering.pdf")]])

# S6  immune-theme composites, overall and molecular-group-stratified
build("FigureS6", [
    [Panel(f"{MF}/Figure5A.pdf", "A")],
    [Panel(f"{MF}/Figure5C.pdf", "B")],
], row_align="top")

# S7  per-signature immune programs and Dunn post-hoc comparisons
build("FigureS7", [
    [Panel(f"{MF}/Figure5B.pdf", "A", scale=0.88),
     Panel(f"{MF}/Figure5D.pdf", "B", scale=0.88)],
    [Panel(f"{MF}/Figure5E.pdf", "C")],
], row_align="top")

# S8  pathway enrichment from the adjusted gene lists
build("FigureS8", [[Panel(f"{MF}/Figure6B.pdf")]])

# S9 unadjusted tie-corrected one-versus-rest volcanoes
build("FigureS9", [[Panel(f"{SF}/FigureS9.pdf")]])

# S10 subtype-stratified survival
build("FigureS10", [
    [Panel(f"{MF}/Figure7B.pdf", "A")],
    [Panel(f"{SF}/FigureS10B.pdf", "B")],
], row_align="top")

# S11 location-stratified Kaplan-Meier
build("FigureS11", [[Panel(f"{SF}/FigureS11.pdf")]])

# S12 proportional-hazards diagnostics
build("FigureS12", [[Panel(f"{SF}/FigureS12.pdf")]])

# S13 survival missingness and inverse-probability weighting
build("FigureS13", [
    [Panel(f"{SF}/FigureS13A.pdf", "A"), Panel(f"{SF}/FigureS13B.pdf", "B")],
    [Panel(f"{SF}/FigureS13C.pdf", "C"), Panel(f"{P}/suppl_ipw_D.pdf", "D")],
], row_align="top")

# S14 single-cell analysis, full
build("FigureS14", [[Panel(f"{SF}/FigureS20_scRNA_support.pdf")]])

# S15 spatial analysis, full
build("FigureS15", [[Panel(f"{P}/suppl_spatial.pdf")]])
