"""
Modify the user's Manuscript_Draft_v1_with_figures.docx to integrate CIBERSORTx LM22 results.

Edits performed:
  - Methods §2.3 Immune Deconvolution: replace the "reserved for planned validation" sentence
    with the actual CIBERSORTx Job #15 execution description.
  - Insert a new Results subsection §3.3 "CIBERSORTx LM22 Cross-Validates the Microglia–MDM
    Axis at Modest Fit" after the brain-tuned signatures paragraph (§3.2).
  - Renumber subsequent §3.3 → §3.4, etc.
  - Add Supplementary Figure S(LM22-QC) and S(LM22-fractions) entries.
  - Save as Manuscript_Draft_v2_with_LM22.docx (preserve original).
"""
from __future__ import annotations
from pathlib import Path
import shutil, re
from docx import Document
from docx.shared import Pt, Inches, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from copy import deepcopy
from docx.oxml.ns import qn

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"
SRC  = OUT / "Manuscript_Draft_v1_with_figures.docx"
DST  = OUT / "Manuscript_Draft_v2_with_LM22.docx"
FIG  = OUT / "figs_lm22_integration"

# Copy then modify
shutil.copy2(SRC, DST)
doc = Document(DST)

# ----------------------------------------------------------------------------
# 1) Update Methods §2.3 (Immune Deconvolution) — replace LM22 sentence
# ----------------------------------------------------------------------------
old_lm22_sentence = "CIBERSORTx LM22 deconvolution was reserved for a planned validation step using the Stanford web tool"
new_lm22_text = (
    "CIBERSORTx LM22 deconvolution was performed as a fifth, independent cross-validation "
    "method using the Stanford web tool (Job #15, 11 June 2026) [11]. The combined 702-biospecimen "
    "TPM matrix (349 Main + 353 secondary) was submitted in absolute mode with batch correction "
    "disabled, quantile normalization disabled (per the package's RNA-seq recommendation), and "
    "100 permutations for significance testing. LM22 outputs were used (a) to compute per-sample "
    "QC metrics (permutation P-value, Pearson fit correlation, RMSE) and (b) as a fifth method in "
    "cross-axis Spearman correlation against MCP-counter, quanTIseq, EPIC, and xCell. LM22 cell-type "
    "fractions were not included in the primary clustering feature matrix; the clustering design "
    "(quanTIseq + ssGSEA, Section 2.5) was specified a priori to avoid inflating redundancy among "
    "estimates of the same cell type"
)

modified_paragraphs = 0
for para in doc.paragraphs:
    if old_lm22_sentence in para.text:
        # Replace the entire paragraph text using run manipulation
        full_text = para.text
        new_text = full_text.replace(
            old_lm22_sentence + " [11].",
            new_lm22_text + "."
        ).replace(
            old_lm22_sentence + ".",
            new_lm22_text + "."
        ).replace(
            old_lm22_sentence,
            new_lm22_text
        )
        # Clear runs and rewrite
        for run in para.runs:
            run.text = ""
        if para.runs:
            para.runs[0].text = new_text
        else:
            para.add_run(new_text)
        modified_paragraphs += 1
        print(f"[Methods §2.3] Updated paragraph: {para.text[:120]}...")

# ----------------------------------------------------------------------------
# 2) Insert new Results subsection §3.3 after the brain-tuned signatures
#    paragraph (currently at index 047 in v1). Use anchor text to find location.
# ----------------------------------------------------------------------------
anchor = "These results establish that LM22-style methods conflate brain-resident microglia and infiltrating MDM"

# Find anchor paragraph
anchor_idx = None
for i, p in enumerate(doc.paragraphs):
    if anchor in p.text:
        anchor_idx = i
        break

if anchor_idx is None:
    print("WARNING: anchor paragraph not found")
else:
    print(f"Anchor found at paragraph {anchor_idx}")

# We will append new content as new paragraphs near the end (since python-docx can't easily
# insert mid-document). Instead, find the existing §3.3 heading and INSERT a NEW §3.3
# subsection right BEFORE it via XML manipulation.

# Find the §3.3 heading to insert before
target_heading = None
for p in doc.paragraphs:
    style_name = p.style.name if p.style is not None else ""
    if style_name.startswith("Heading 2") and p.text.strip().startswith("3.3."):
        target_heading = p
        break

if target_heading is None:
    print("WARNING: §3.3 heading not found - appending at end instead")

# Build the new subsection content
NEW_SECTION = [
    ("Heading 2", "3.3. CIBERSORTx LM22 Cross-Validates the Ecotype Architecture at Reduced Fit"),
    ("Normal",
     "We next examined whether the Inflamed–Intermediate–Immune-desert ecotypes were also "
     "recovered when an entirely separate deconvolution algorithm — CIBERSORTx with the LM22 "
     "peripheral-blood signature matrix [11] — was applied to the same cohort, while bearing in "
     "mind that LM22 is known to lack brain-resident microglia and is therefore expected to "
     "fit pediatric brain TME poorly. On the n = 349 Main cohort, LM22 deconvolution showed "
     "the predicted modest goodness-of-fit (median permutation P = 0.17, median Pearson "
     "correlation = 0.08, median RMSE = 1.07), with only 91 of 349 samples (26.1 %) reaching "
     "the conventional P < 0.05 threshold (Supplementary Figure S5). Consistent with the "
     "ecotype framework, however, four of eight LM22 super-categories nevertheless differed "
     "significantly across the three ecotypes after BH adjustment (Supplementary Figure S6, "
     "Supplementary Table S22)."),
    ("Normal",
     "The dominant LM22 signal — Monocytes/Macrophages combined fraction — was highest in the "
     "Inflamed ecotype (median 0.62), intermediate in the Intermediate ecotype (0.54), and "
     "lowest in the Immune-desert ecotype (0.46; Kruskal–Wallis q = 2.6 × 10⁻⁷). This direction "
     "matches the brain-tuned ssGSEA Mo-TAM and Klemm MDM signals (Section 3.2), confirming "
     "that the macrophage axis is reproduced across analytical pipelines even though the "
     "underlying signature sources are independent. CD4 T-cell, NK-cell, and Granulocyte LM22 "
     "fractions, by contrast, were paradoxically higher in the Immune-desert ecotype (q = "
     "1.2 × 10⁻⁵, 1.2 × 10⁻⁵, and 4.5 × 10⁻², respectively) — an apparent inversion that we "
     "interpret not as a contradiction but as a quantification artefact of LM22 fractions "
     "summing to unity: when the dominant Monocyte/Macrophage fraction drops in the Desert "
     "ecotype, the remaining (small) lymphoid and granulocyte estimates are mathematically "
     "compressed into a larger proportional share. Critically, the LM22 CD8 T-cell fraction "
     "did not distinguish ecotypes (q = 0.86), illustrating the limited resolution of LM22 "
     "for the lymphocyte axis in this tissue."),
    ("Normal",
     "Cross-method Spearman correlation analysis with LM22 added as a fifth method reinforced "
     "the same picture (Supplementary Table S7-extended). The CD8 T-cell axis showed the "
     "highest LM22-vs-orthogonal agreement against MCP-counter (ρ = 0.40) and quanTIseq "
     "(ρ = 0.31), while xCell correlated poorly with LM22 (ρ = 0.01). The macrophage axis "
     "showed modest LM22-vs-orthogonal agreement (ρ up to 0.37 against xCell), and the NK-cell "
     "and Monocyte axes were essentially uncorrelated between LM22 and any other method "
     "(|ρ| ≤ 0.10), reflecting both the small absolute fractions involved and the fundamental "
     "mismatch between LM22 and brain TME. Collectively, the CIBERSORTx LM22 results "
     "(i) confirm that the macrophage axis dominates the Inflamed-vs-Desert contrast across "
     "five independent algorithms, (ii) underline the inadequacy of any single deconvolution "
     "method — particularly LM22 — for resolving brain TME, and (iii) further justify the "
     "use of brain-tuned ssGSEA scoring as the primary biological readout in this manuscript."),
]

if target_heading is not None:
    # Insert before target_heading using XML element insertion
    target_elem = target_heading._element
    for style, text in reversed(NEW_SECTION):  # reversed so insertions stack in correct order
        new_p = doc.add_paragraph()
        if style == "Heading 2":
            try: new_p.style = doc.styles["Heading 2"]
            except KeyError: pass
        # Set the text via run with appropriate formatting
        for r in list(new_p.runs):
            r.text = ""
        r = new_p.add_run(text)
        if style == "Heading 2":
            r.bold = True
            r.font.size = Pt(13)
        else:
            r.font.size = Pt(11)
        # Move element from end-of-doc to position before target
        new_elem = new_p._element
        new_elem.getparent().remove(new_elem)
        target_elem.addprevious(new_elem)

# ----------------------------------------------------------------------------
# 3) Renumber subsequent §3.X headings
# ----------------------------------------------------------------------------
renumber_map = {
    "3.3. Consensus Clustering Identifies Three Reproducible Immune Ecotypes":
        "3.4. Consensus Clustering Identifies Three Reproducible Immune Ecotypes",
    "3.4. Ecotype Distribution Across Histone Class and Anatomical Location":
        "3.5. Ecotype Distribution Across Histone Class and Anatomical Location",
    "3.5. Ecotype Is an Independent Prognostic Factor for Overall Survival":
        "3.6. Ecotype Is an Independent Prognostic Factor for Overall Survival",
    "3.6. Pathway Activity Distinguishes Ecotypes Without Confounding by Tumor Proliferation":
        "3.7. Pathway Activity Distinguishes Ecotypes Without Confounding by Tumor Proliferation",
}
renamed = 0
for p in doc.paragraphs:
    text = p.text.strip()
    if text in renumber_map:
        new_text = renumber_map[text]
        for r in p.runs:
            r.text = ""
        if p.runs:
            p.runs[0].text = new_text
        else:
            p.add_run(new_text)
        renamed += 1
        print(f"[Renumber] '{text}' -> '{new_text}'")

# ----------------------------------------------------------------------------
# 4) Append Supplementary Figure / Table descriptions at the end (before References)
# ----------------------------------------------------------------------------
references_anchor = None
for p in doc.paragraphs:
    sn = p.style.name if p.style is not None else ""
    if sn.startswith("Heading 1") and p.text.strip() == "References":
        references_anchor = p
        break

supp_content = [
    ("Heading 1", "Supplementary materials (new — LM22 integration)"),
    ("Normal", "Supplementary Figure S5. CIBERSORTx LM22 quality-control distribution on the Main cohort (n = 349). Histograms of LM22 permutation P-value (red dashed line at P = 0.05), Pearson fit correlation, and root-mean-squared error. Only 91 of 349 (26.1 %) samples reach P < 0.05, with median P = 0.17 and median correlation = 0.08, confirming that LM22 — derived from peripheral-blood leukocyte expression — does not fit pediatric brain TME well at the per-sample level. The same pattern was observed on the full 702-biospecimen pool (Supplementary Table S21)."),
    ("Normal", "Supplementary Figure S6. CIBERSORTx LM22 cell-type fractions across the three immune ecotypes. Boxplots of LM22-derived super-category fractions (CD8 T cells, CD4 T cells, NK cells, B cells, Monocytes/Macrophages, Dendritic cells, Tregs/γδ T cells, Granulocytes) stratified by ecotype. Per-panel Kruskal–Wallis BH-adjusted q-values are annotated. The Monocyte/Macrophage axis (q = 2.6 × 10⁻⁷) recapitulates the Inflamed > Intermediate > Immune-desert ordering observed with brain-tuned ssGSEA scoring; CD4 T-cell, NK, and Granulocyte fractions show an inverted gradient that is interpreted as a compositional consequence of LM22 fractions summing to unity (see main text §3.3)."),
    ("Normal", "Supplementary Table S21. CIBERSORTx LM22 quality-control statistics stratified by sub-pool (full 702-pool, Main n = 349, BRAF/RTK-altered secondary n = 353). Median P-value, fit correlation, and RMSE per sub-pool with proportion of samples reaching P < 0.05."),
    ("Normal", "Supplementary Table S22. CIBERSORTx LM22 cell-type fraction Kruskal–Wallis test by ecotype (Main n = 349), with super-category-level median fractions and BH-adjusted q-values."),
    ("Normal", "Supplementary Table S7-extended. Spearman correlation matrix among five deconvolution methods (LM22, MCP-counter, quanTIseq, EPIC, xCell) on four key cell-type axes (CD8 T, NK, Macrophage, Monocyte), computed on n = 349 paired observations."),
]

if references_anchor is not None:
    target_elem = references_anchor._element
    for style, text in reversed(supp_content):
        new_p = doc.add_paragraph()
        if style == "Heading 1":
            try: new_p.style = doc.styles["Heading 1"]
            except KeyError: pass
        r = new_p.add_run(text)
        if style == "Heading 1":
            r.bold = True; r.font.size = Pt(15)
        else:
            r.font.size = Pt(10); r.italic = True
        new_elem = new_p._element
        new_elem.getparent().remove(new_elem)
        target_elem.addprevious(new_elem)
else:
    # append at end
    for style, text in supp_content:
        new_p = doc.add_paragraph()
        if style == "Heading 1":
            try: new_p.style = doc.styles["Heading 1"]
            except KeyError: pass
        r = new_p.add_run(text); r.italic = (style != "Heading 1"); r.bold = (style == "Heading 1")

# ----------------------------------------------------------------------------
# 5) Embed the two new figures in the Supplementary block
# ----------------------------------------------------------------------------
# Find the Supplementary Figure S5 paragraph and insert image before it
for p in doc.paragraphs:
    if "Supplementary Figure S5." in p.text and "LM22 quality-control" in p.text:
        # add picture in a NEW paragraph BEFORE this one
        new_p = doc.add_paragraph()
        run = new_p.add_run()
        run.add_picture(str(FIG / "LM22_QC_distribution_Main.png"), width=Inches(6.3))
        new_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        new_elem = new_p._element
        new_elem.getparent().remove(new_elem)
        p._element.addprevious(new_elem)
        break

for p in doc.paragraphs:
    if "Supplementary Figure S6." in p.text and "cell-type fractions across" in p.text:
        new_p = doc.add_paragraph()
        run = new_p.add_run()
        run.add_picture(str(FIG / "LM22_fractions_by_ecotype.png"), width=Inches(6.3))
        new_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        new_elem = new_p._element
        new_elem.getparent().remove(new_elem)
        p._element.addprevious(new_elem)
        break

doc.save(DST)
print(f"\nWrote: {DST}  ({DST.stat().st_size/1024:.1f} KB)")
print(f"Methods §2.3 modifications: {modified_paragraphs}")
print(f"Section heading renames: {renamed}")
