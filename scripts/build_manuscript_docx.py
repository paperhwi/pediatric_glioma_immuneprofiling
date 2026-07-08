"""
Convert manuscript_rough_draft.md → Manuscript_rough_draft.docx
with embedded figure images, Table 1 inline, and clean Word formatting.
"""
from __future__ import annotations
from pathlib import Path
import re
import pandas as pd
from docx import Document
from docx.shared import Pt, Cm, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"

MD_PATH    = BASE / "manuscript_rough_draft.md"
TABLE1_TSV = OUT / "Table1_n349.tsv"
OUTPUT_DOCX= OUT / "Manuscript_rough_draft.docx"

# Map "[**Figure X** | filename · filename]" lines to actual file paths to embed
FIG_MAP = {
    # Fig 1: placeholder (BioRender to be drafted)
    "Fig 1": None,
    "Fig 2": [
        OUT / "figs_step4" / "LM22_QC_distribution.png",
        OUT / "figs_step4" / "LM22_vs_methods_heatmap.png",
        OUT / "figs_step4" / "ssGSEA_microglia_vs_LM22.png",
    ],
    "Fig 3": [
        OUT / "figs_step5" / "PAC_silhouette_vs_k.png",
        OUT / "figs_step5" / "ecotype_PCA_UMAP.png",
        OUT / "figs_step5" / "ecotype_k3_feature_means_heatmap.png",
        OUT / "figs_step5" / "ecotype_k3_by_cohort_group.png",
    ],
    "Fig 4": [
        OUT / "figs_step7" / "theme_composite_boxplot.png",
        OUT / "figs_step7" / "signature_ecotype_heatmap.png",
        OUT / "figs_step7" / "theme_by_cohort_stratified.png",
        OUT / "figs_step6" / "ecotype_boxplots_focus.png",
        OUT / "figs_step6" / "Dunn_signed_log10q_focus.png",
    ],
    "Fig 5": [
        OUT / "figs_step8" / "volcanos_per_ecotype.png",
        OUT / "figs_step8" / "top_enrichment_per_ecotype.png",
    ],
    "Fig 6": [
        OUT / "figs_step9" / "KM_overall_by_ecotype.png",
        OUT / "figs_step9" / "KM_by_cohort_grid.png",
        OUT / "figs_step9" / "Cox_forest.png",
    ],
    "Fig 7": [
        OUT / "figs_step10" / "ecotype_distribution_by_cohort_family.png",
        OUT / "figs_step10" / "MAPK_by_cohort_family.png",
        OUT / "figs_step10" / "BRAF_ALT_overlay_on_main_PCA.png",
    ],
    "Suppl Fig S3": [OUT / "figs_suppl" / "S3_highconf_subset_reclustering.png"],
    "Suppl Fig S4": [
        OUT / "figs_suppl" / "S4_fullpool_sensitivity.png",
        OUT / "figs_suppl" / "S4_fullpool_PCA.png",
    ],
}

text = MD_PATH.read_text(encoding="utf-8")

# Split into ordered blocks: each block is (kind, content)
# kind ∈ {"h1","h2","h3","h4","p","fig","table1","blockquote","listitem","hr"}
def parse_md(md: str):
    blocks = []
    lines = md.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("# "):
            blocks.append(("h1", ln[2:].strip())); i += 1
        elif ln.startswith("## "):
            blocks.append(("h2", ln[3:].strip())); i += 1
        elif ln.startswith("### "):
            blocks.append(("h3", ln[4:].strip())); i += 1
        elif ln.startswith("#### "):
            blocks.append(("h4", ln[5:].strip())); i += 1
        elif ln.strip() == "---":
            blocks.append(("hr", "")); i += 1
        elif ln.startswith("> "):
            buf = [ln[2:]]
            i += 1
            while i < len(lines) and lines[i].startswith("> "):
                buf.append(lines[i][2:]); i += 1
            blocks.append(("blockquote", "\n".join(buf).strip()))
        elif ln.startswith("[**Fig ") or ln.startswith("[**Suppl Fig "):
            # Figure placeholder line
            m = re.match(r"^\[\*\*(Fig \d+|Suppl Fig S\d+)\*\*\s*\|(.*)\]", ln.strip())
            if m:
                blocks.append(("fig", (m.group(1), m.group(2).strip())))
            else:
                blocks.append(("p", ln))
            i += 1
        elif ln.strip().startswith("- "):
            blocks.append(("li", ln.strip()[2:])); i += 1
        elif ln.strip().startswith("|") and ln.strip().endswith("|"):
            # Skip table — we'll handle Table 1 manually
            while i < len(lines) and lines[i].strip().startswith("|"):
                i += 1
        elif ln.strip() == "":
            i += 1
        else:
            blocks.append(("p", ln)); i += 1
    return blocks

blocks = parse_md(text)

# ----- Build docx -----
doc = Document()
section = doc.sections[0]
section.top_margin = Cm(2.2); section.bottom_margin = Cm(2.2)
section.left_margin = Cm(2.0); section.right_margin = Cm(2.0)

# Style defaults
style = doc.styles["Normal"]
style.font.name = "Calibri"; style.font.size = Pt(11)

def add_paragraph_with_inline_md(parent, text, *, bold_all=False, italic_all=False, size=None):
    # parent may be a Document (top-level) or an existing paragraph
    if hasattr(parent, "add_paragraph"):
        p = parent.add_paragraph()
    else:
        p = parent
    # Split by ** and * for bold / italic
    tokens = re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text)
    for tok in tokens:
        if not tok: continue
        run = p.add_run()
        if tok.startswith("**") and tok.endswith("**"):
            run.text = tok[2:-2]; run.bold = True
        elif tok.startswith("*") and tok.endswith("*"):
            run.text = tok[1:-1]; run.italic = True
        else:
            # Inline code/backticks → just strip them for now
            run.text = tok.replace("`", "")
        if bold_all: run.bold = True
        if italic_all: run.italic = True
        if size: run.font.size = size
    return p

# Helpers
def add_heading(text, level):
    sizes = {1: 18, 2: 14, 3: 12, 4: 11}
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = True; r.font.size = Pt(sizes.get(level, 11))
    if level == 1:
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    if level <= 2:
        p.paragraph_format.space_before = Pt(12)
        p.paragraph_format.space_after  = Pt(6)

def add_figure(name, caption):
    paths = FIG_MAP.get(name)
    if not paths:
        # placeholder
        p = doc.add_paragraph()
        r = p.add_run(f"[ {name}: placeholder — to be drafted ]")
        r.italic = True; r.font.color.rgb = RGBColor(0x99, 0x66, 0x00)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        return
    for path in paths:
        if Path(path).exists():
            p = doc.add_paragraph()
            run = p.add_run()
            run.add_picture(str(path), width=Inches(6.3))
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cap = doc.add_paragraph()
    r = cap.add_run(f"{name}.  {caption}")
    r.italic = True; r.font.size = Pt(9)
    cap.alignment = WD_ALIGN_PARAGRAPH.LEFT

# Add Table 1 as a real Word table
def add_table1():
    tbl_df = pd.read_csv(TABLE1_TSV, sep="\t").fillna("")
    title = doc.add_paragraph()
    tr = title.add_run("Table 1.  Clinical and molecular characteristics of the pediatric glioma Main cohort (n = 349).")
    tr.bold = True; tr.font.size = Pt(10)
    tbl = doc.add_table(rows=len(tbl_df)+1, cols=len(tbl_df.columns))
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    # Header
    hdr = tbl.rows[0]
    for i, c in enumerate(tbl_df.columns):
        cell = hdr.cells[i]; cell.text = ""
        p = cell.paragraphs[0]; r = p.add_run(c)
        r.bold = True; r.font.size = Pt(8)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER if i > 0 else WD_ALIGN_PARAGRAPH.LEFT
        cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
    # Data
    for ri, (_, row) in enumerate(tbl_df.iterrows(), start=1):
        is_section_header = all(v == "" for v in row.values[1:])
        for ci, col in enumerate(tbl_df.columns):
            cell = tbl.rows[ri].cells[ci]; cell.text = ""
            p = cell.paragraphs[0]; r = p.add_run(str(row[col]))
            r.font.size = Pt(8)
            if is_section_header and ci == 0:
                r.bold = True; r.italic = True
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if ci > 0 else WD_ALIGN_PARAGRAPH.LEFT
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
    # Borders
    for row in tbl.rows:
        for cell in row.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            tcBorders = OxmlElement("w:tcBorders")
            for side in ("top", "left", "bottom", "right"):
                bd = OxmlElement(f"w:{side}")
                bd.set(qn("w:val"), "single"); bd.set(qn("w:sz"), "4"); bd.set(qn("w:color"), "000000")
                tcBorders.append(bd)
            tcPr.append(tcBorders)

# Render blocks
table1_inserted = False
for kind, content in blocks:
    if kind == "h1":
        add_heading(content, 1)
    elif kind == "h2":
        add_heading(content, 2)
        # After Methods section heading, insert Table 1 once
        if not table1_inserted and content.lower().startswith("2."):
            pass  # we'll insert Table 1 right after Results 3.1 instead
    elif kind == "h3":
        add_heading(content, 3)
        if (not table1_inserted) and "3.1" in content.lower():
            pass
    elif kind == "h4":
        add_heading(content, 4)
    elif kind == "hr":
        # horizontal rule -> page break
        doc.add_paragraph()
    elif kind == "blockquote":
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Cm(0.6)
        p.paragraph_format.right_indent = Cm(0.6)
        add_paragraph_with_inline_md(p, content, italic_all=True, size=Pt(10))
        p.paragraph_format.space_after = Pt(6)
    elif kind == "li":
        p = doc.add_paragraph(style="List Bullet")
        add_paragraph_with_inline_md(p, content)
    elif kind == "fig":
        name, caption = content
        add_figure(name, caption)
    else:
        # Regular paragraph
        if not content.strip():
            continue
        p = doc.add_paragraph()
        add_paragraph_with_inline_md(p, content)
        # Insert Table 1 right after the "3.1." paragraph (the cohort description)
        if (not table1_inserted) and "349-patient" in content:
            doc.add_paragraph()  # spacer
            add_table1()
            table1_inserted = True

# Fallback: if Table 1 wasn't inserted (no matching sentence), append before References
if not table1_inserted:
    doc.add_paragraph()
    add_table1()

doc.save(OUTPUT_DOCX)
print(f"Wrote: {OUTPUT_DOCX}")
print(f"Size: {OUTPUT_DOCX.stat().st_size / 1024:.1f} KB")
