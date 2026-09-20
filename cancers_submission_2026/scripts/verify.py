# -*- coding: utf-8 -*-
"""Verification pass for the Cancers figure and supplementary package.

  1. every figure file exists, opens, and has one page
  2. the panel letters drawn on each composed figure match the letters its
     legend describes
  3. numbering is contiguous and nothing is referenced that does not exist
  4. builds the citation checklist written into the table workbook
"""
import os
import re
import sys
import glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pymupdf
import legends as L
from build_tables import TABLES

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# panel letters actually drawn or already present inside each figure
DRAWN = {
    "Figure1": "AB", "Figure2": "ABC", "Figure3": "ABCDEFG", "Figure4": "ABCD",
    "Figure5": "ABCDE", "Figure6": "ABCDEF", "Figure7": "AB",
    "FigureS1": "", "FigureS2": "", "FigureS3": "", "FigureS4": "",
    "FigureS5": "", "FigureS6": "", "FigureS7": "AB", "FigureS8": "ABC",
    "FigureS9": "", "FigureS10": "", "FigureS11": "AB", "FigureS12": "",
    "FigureS13": "", "FigureS14": "ABCD", "FigureS15": "ABCDE",
    "FigureS16": "ABCDEF",
}

LETTER = re.compile(r"\(([A-H])(?:–([A-H]))?\)")


def legend_letters(body):
    out = set()
    for a, b in LETTER.findall(body):
        if b:
            out.update(chr(c) for c in range(ord(a), ord(b) + 1))
        else:
            out.add(a)
    return "".join(sorted(out))


def main():
    problems = []

    # ---- 1. files -------------------------------------------------------
    mains = [f"{ROOT}/main/Figure{i}" for i in range(1, 8)]
    supps = [f"{ROOT}/suppl/FigureS{i}" for i in range(1, 17)]
    for stem in mains + supps:
        for ext in (".pdf", ".png"):
            if not os.path.exists(stem + ext):
                problems.append(f"missing file {os.path.basename(stem)}{ext}")
        if os.path.exists(stem + ".pdf"):
            d = pymupdf.open(stem + ".pdf")
            if len(d) != 1:
                problems.append(f"{os.path.basename(stem)}.pdf has {len(d)} pages")
            d.close()
    stray = [f for f in glob.glob(f"{ROOT}/main/*.pdf") + glob.glob(f"{ROOT}/suppl/*.pdf")
             if os.path.splitext(f)[0] not in mains + supps]
    for f in stray:
        problems.append(f"unexpected file in the package: {os.path.basename(f)}")

    # ---- 2. panel letters ----------------------------------------------
    print("panel letters  (legend vs figure)")
    for lab, title, body in L.MAIN + L.SUPPL:
        key = lab.replace(" ", "")
        want = legend_letters(body)
        have = DRAWN.get(key, "?")
        ok = "ok" if want == have else "MISMATCH"
        if want != have:
            problems.append(f"{lab}: legend describes {want or '-'} but figure carries {have or '-'}")
        print(f"  {lab:<12} legend {want or '-':<8} figure {have or '-':<8} {ok}")

    # ---- 3. cross-references -------------------------------------------
    have_fig = {f"Figure S{i}" for i in range(1, 17)} | {f"Figure {i}" for i in range(1, 8)}
    have_tab = {f"Table S{i}" for i in range(1, len(TABLES) + 1)} | \
               {f"Table {i}" for i in range(1, 6)}
    ref = re.compile(r"(Figure|Table)s?\s+(S?\d+(?:\s*(?:[\u2013-]|,|and)\s*S?\d+)*)")
    piece = re.compile(r"S?\d+")

    cited = {}

    def expand(kind, run):
        toks = piece.findall(run)
        seps = re.findall(r"[\u2013-]|,|and", run)
        items, i = [], 0
        while i < len(toks):
            if i < len(seps) and seps[i] in "\u2013-" and i + 1 < len(toks):
                a, b = toks[i], toks[i + 1]
                pre = "S" if a.startswith("S") else ""
                items += [f"{pre}{n}" for n in
                          range(int(a.lstrip("S")), int(b.lstrip("S")) + 1)]
                i += 2
            else:
                items.append(toks[i]); i += 1
        return [f"{kind} {t}" for t in items]

    for lab, title, body in L.MAIN + L.SUPPL + L.MAIN_TABLES:
        for kind, run in ref.findall(body):
            for it in expand(kind, run):
                cited.setdefault(it, []).append(lab)

    for it in sorted(cited):
        pool = have_fig if it.startswith("Figure") else have_tab
        if it not in pool:
            problems.append(f"legend of {cited[it][0]} cites {it}, which does not exist")

    # ---- 4. citation checklist -----------------------------------------
    rows = []
    for i in range(1, 8):
        it = f"Figure {i}"
        rows.append(("Main figure", it, L.MAIN[i - 1][1],
                     "; ".join(cited.get(it, [])) or "-", "required"))
    for i in range(1, 17):
        it = f"Figure S{i}"
        rows.append(("Supplementary figure", it, L.SUPPL[i - 1][1],
                     "; ".join(cited.get(it, [])) or "-", "required"))
    for num, title, desc, supports, _ in TABLES:
        it = f"Table {num}"
        rows.append(("Supplementary table", it, title,
                     "; ".join(cited.get(it, [])) or "-",
                     "required" if "Methods" not in supports else "required (Methods)"))

    uncited = [r[1] for r in rows if r[3] == "-"]
    print(f"\nitems not yet cross-referenced from a legend ({len(uncited)}):")
    print("  " + ", ".join(uncited) if uncited else "  none")

    print("\nproblems:")
    if problems:
        for p in problems:
            print("  !!", p)
    else:
        print("  none")
    return rows, problems


if __name__ == "__main__":
    main()
