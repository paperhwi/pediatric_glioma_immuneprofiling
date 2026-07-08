"""
Build full manuscript docx by combining:
  - Expanded Abstract / Introduction / Discussion / Limitations / Conclusion (manuscript_full_prose.md)
  - Methods + Results (with embedded figures) from manuscript_rough_draft.md
  - Table 1 inserted after §3.1
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

PROSE_MD = BASE / "manuscript_full_prose.md"
ROUGH_MD = BASE / "manuscript_rough_draft.md"
TABLE1_TSV = OUT / "Table1_n349.tsv"
OUTPUT_DOCX = OUT / "Manuscript_full_prose_v2.docx"

FIG_MAP = {
    "Fig 1": [OUT / "figs_fig1" / "Fig1_workflow.png"],
    "Fig 2": [
        OUT / "figs_300dpi" / "Fig2A_LM22_QC.png",
        OUT / "figs_300dpi" / "Fig2B_LM22_vs_methods_heatmap.png",
        OUT / "figs_300dpi" / "Fig2C_ssGSEA_vs_LM22.png",
        OUT / "figs_300dpi" / "Fig2D_4method_by_ecotype.png",
        OUT / "figs_300dpi" / "Fig2D_4method_heatmap.png",
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

# ---------- Markdown parser (same as build_manuscript_docx.py) ----------
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
            buf = [ln[2:]]; i += 1
            while i < len(lines) and lines[i].startswith("> "):
                buf.append(lines[i][2:]); i += 1
            blocks.append(("blockquote", "\n".join(buf).strip()))
        elif ln.startswith("[**Fig ") or ln.startswith("[**Suppl Fig "):
            m = re.match(r"^\[\*\*(Fig \d+|Suppl Fig S\d+)\*\*\s*\|(.*)\]", ln.strip())
            if m:
                blocks.append(("fig", (m.group(1), m.group(2).strip())))
            else:
                blocks.append(("p", ln))
            i += 1
        elif ln.strip().startswith("- "):
            blocks.append(("li", ln.strip()[2:])); i += 1
        elif ln.strip().startswith("|") and ln.strip().endswith("|"):
            while i < len(lines) and lines[i].strip().startswith("|"):
                i += 1
        elif ln.strip() == "":
            i += 1
        else:
            blocks.append(("p", ln)); i += 1
    return blocks

# ---------- Section extraction ----------
def extract_sections(blocks):
    """Return list of (section_name_lower, [block, ...]) ordered."""
    sections = []
    cur = None
    for b in blocks:
        if b[0] == "h2":
            cur = (b[1], [b]); sections.append(cur)
        elif cur is not None:
            cur[1].append(b)
        else:
            sections.append(("__head__", [b]))
    return sections

prose_blocks = parse_md(PROSE_MD.read_text(encoding="utf-8"))
rough_blocks = parse_md(ROUGH_MD.read_text(encoding="utf-8"))

prose_sections = {s[0].lower(): s[1] for s in extract_sections(prose_blocks) if s[0] != "__head__"}
rough_sections_list = extract_sections(rough_blocks)

# Head/title (everything before first ## in PROSE)
head_blocks = []
for b in prose_blocks:
    if b[0] == "h2":
        break
    head_blocks.append(b)

# Build assembled blocks in correct manuscript order:
#   head → Abstract → 1. Introduction → 2. Methods (rough) → 3. Results (rough) → 4. Discussion (prose) → Limitations (prose) → Conclusions (prose) → Figure/Table inventory (rough)
def find_prose_section(prefix):
    for k, v in prose_sections.items():
        if k.lower().startswith(prefix):
            return v
    return None

def find_rough_section(prefix):
    for name, blocks in rough_sections_list:
        if name.lower().startswith(prefix):
            return blocks
    return None

assembled = []
assembled.extend(head_blocks)
assembled.extend(find_prose_section("abstract") or [])
assembled.extend(find_prose_section("1.") or [])
assembled.extend(find_rough_section("2.") or [])
assembled.extend(find_rough_section("3.") or [])
assembled.extend(find_prose_section("4.") or [])
assembled.extend(find_prose_section("limit") or [])
assembled.extend(find_prose_section("conclus") or [])
# Figure & Table inventory (use the rough doc's tables and headers)
fig_inv = find_rough_section("figure inventory") or []
tab_inv = find_rough_section("table inventory") or []
assembled.extend(fig_inv); assembled.extend(tab_inv)

# Combine consecutive h2 sections of same name (we may have duplicate "## Limitations" headers etc.)
seen_h2 = set()
deduped = []
for b in assembled:
    if b[0] == "h2":
        if b[1].lower() in seen_h2:
            continue
        seen_h2.add(b[1].lower())
    deduped.append(b)
assembled = deduped

# ---------- Word document construction (same as rough builder) ----------
doc = Document()
section = doc.sections[0]
section.top_margin = Cm(2.2); section.bottom_margin = Cm(2.2)
section.left_margin = Cm(2.0); section.right_margin = Cm(2.0)
style = doc.styles["Normal"]
style.font.name = "Calibri"; style.font.size = Pt(11)

def add_paragraph_with_inline_md(parent, text, *, italic_all=False, bold_all=False, size=None):
    if hasattr(parent, "add_paragraph"):
        p = parent.add_paragraph()
    else:
        p = parent
    tokens = re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text)
    for tok in tokens:
        if not tok: continue
        run = p.add_run()
        if tok.startswith("**") and tok.endswith("**"):
            run.text = tok[2:-2]; run.bold = True
        elif tok.startswith("*") and tok.endswith("*"):
            run.text = tok[1:-1]; run.italic = True
        else:
            run.text = tok.replace("`", "")
        if italic_all: run.italic = True
        if bold_all:   run.bold = True
        if size:       run.font.size = size
    return p

def add_heading(text, level):
    sizes = {1: 18, 2: 14, 3: 12, 4: 11}
    p = doc.add_paragraph()
    r = p.add_run(text); r.bold = True; r.font.size = Pt(sizes.get(level, 11))
    if level <= 2:
        p.paragraph_format.space_before = Pt(12)
        p.paragraph_format.space_after  = Pt(6)

def add_figure(name, caption):
    paths = FIG_MAP.get(name)
    if not paths:
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

def add_table1():
    tbl_df = pd.read_csv(TABLE1_TSV, sep="\t").fillna("")
    title = doc.add_paragraph()
    tr = title.add_run("Table 1.  Clinical and molecular characteristics of the pediatric glioma Main cohort (n = 349).")
    tr.bold = True; tr.font.size = Pt(10)
    tbl = doc.add_table(rows=len(tbl_df)+1, cols=len(tbl_df.columns))
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    hdr = tbl.rows[0]
    for i, c in enumerate(tbl_df.columns):
        cell = hdr.cells[i]; cell.text = ""
        p = cell.paragraphs[0]; r = p.add_run(c); r.bold = True; r.font.size = Pt(8)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER if i > 0 else WD_ALIGN_PARAGRAPH.LEFT
        cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
    for ri, (_, row) in enumerate(tbl_df.iterrows(), start=1):
        is_header_row = all(v == "" for v in row.values[1:])
        for ci, col in enumerate(tbl_df.columns):
            cell = tbl.rows[ri].cells[ci]; cell.text = ""
            p = cell.paragraphs[0]; r = p.add_run(str(row[col])); r.font.size = Pt(8)
            if is_header_row and ci == 0:
                r.bold = True; r.italic = True
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if ci > 0 else WD_ALIGN_PARAGRAPH.LEFT
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
    for row in tbl.rows:
        for cell in row.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            tcBorders = OxmlElement("w:tcBorders")
            for side in ("top","left","bottom","right"):
                bd = OxmlElement(f"w:{side}")
                bd.set(qn("w:val"), "single"); bd.set(qn("w:sz"), "4"); bd.set(qn("w:color"), "000000")
                tcBorders.append(bd)
            tcPr.append(tcBorders)

table1_inserted = False
for kind, content in assembled:
    if kind == "h1":
        add_heading(content, 1)
    elif kind == "h2":
        add_heading(content, 2)
    elif kind == "h3":
        add_heading(content, 3)
    elif kind == "h4":
        add_heading(content, 4)
    elif kind == "hr":
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
        if not content.strip(): continue
        p = doc.add_paragraph()
        add_paragraph_with_inline_md(p, content)
        if (not table1_inserted) and "349-patient" in content:
            doc.add_paragraph()
            add_table1()
            table1_inserted = True

if not table1_inserted:
    doc.add_paragraph(); add_table1()

doc.save(OUTPUT_DOCX)
print(f"Wrote: {OUTPUT_DOCX}  ({OUTPUT_DOCX.stat().st_size/1024:.1f} KB)")
