"""Replace text inside an existing vector PDF, keeping everything else intact.

Used to bring a small number of reused panels into line with the wording and
spelling required for the Cancers submission.  The original glyphs are
redacted and the replacement is typeset in the same font (DejaVu Sans, which
is what the source figures use) at the same size and baseline.
"""
import os
import matplotlib
import pymupdf

FONTFILE = os.path.join(os.path.dirname(matplotlib.__file__),
                        "mpl-data", "fonts", "ttf", "DejaVuSans.ttf")
FONTFILE_B = os.path.join(os.path.dirname(matplotlib.__file__),
                          "mpl-data", "fonts", "ttf", "DejaVuSans-Bold.ttf")


def patch(src, dst, replacements, verbose=True):
    """replacements: list of (old_text, new_text). Matching is on whole spans
    that contain old_text; the span text is rewritten with the substitution."""
    doc = pymupdf.open(src)
    page = doc[0]
    todo = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                txt = span["text"]
                new = txt
                for old, rep in replacements:
                    if old in new:
                        new = new.replace(old, rep)
                if new != txt:
                    todo.append((span, new))
    if not todo:
        doc.close()
        raise SystemExit(f"nothing matched in {src}")

    for span, _ in todo:
        page.add_redact_annot(pymupdf.Rect(span["bbox"]))
    page.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE,
                          graphics=pymupdf.PDF_REDACT_LINE_ART_NONE)

    bold = page.insert_font(fontname="DJVB", fontfile=FONTFILE_B)
    reg = page.insert_font(fontname="DJV", fontfile=FONTFILE)
    for span, new in todo:
        fn = "DJVB" if "Bold" in span["font"] else "DJV"
        page.insert_text((span["origin"][0], span["origin"][1]), new,
                         fontname=fn, fontsize=span["size"], color=(0, 0, 0))
        if verbose:
            print(f"   {span['text'][:60]!r} -> {new[:60]!r}")
    doc.save(dst, garbage=4, deflate=True)
    doc.close()
    return dst


if __name__ == "__main__":
    SRC = ("/mnt/user-data/uploads/Open PBTA/Revision/REVISION PACKAGE 260916/"
           "1. Figure 300dpi/Supplementary figures")
    OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "panels")

    print("FigureS21_spatial.pdf")
    patch(f"{SRC}/FigureS21_spatial.pdf", f"{OUT}/suppl_spatial.pdf", [
        ("Spatial organisation of the immune programmes in H3.3 K27M paediatric ",
         "Spatial organization of the immune programs in H3.3 K27M pediatric "),
        ("same-neighbour", "same-neighbor"),
        ("Immune programmes are spatially structured, not random",
         "Immune programs are spatially structured, not randomly distributed"),
        ("Every section contains all three",
         "Spots mapping to each of the three bulk"),
        ("ecotypes", "ecotype centroids were represented in every section"),
        ("Myeloid-lymphoid coupling is weaker",
         "Myeloid–lymphoid coupling is weaker within"),
        ("within tissue than in bulk", "tissue than across bulk tumors"),
    ])

    print("FigureS13D.pdf")
    patch(f"{SRC}/FigureS13D.pdf", f"{OUT}/suppl_ipw_D.pdf", [
        ("Standardised", "Standardized"), ("stabilised", "stabilized"),
    ])
