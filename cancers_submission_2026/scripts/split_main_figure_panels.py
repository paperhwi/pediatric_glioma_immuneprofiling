"""Split the assembled main figures into one vector PDF and 600 dpi PNG per panel.

The crop boxes are defined on the final 600 dpi assembled PNGs. The same
normalized rectangles are applied to the vector PDFs, so the PDF panels remain
vector artwork rather than raster screenshots.
"""
from __future__ import annotations

import copy
import hashlib
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "1. Main figures"
OUTPUT = ROOT / "Figure separate files"

# Pixel coordinates on the final assembled PNGs: left, top, right, bottom.
CROPS = {
    "Figure1": {
        "A": (0, 0, 8214, 3179),
        "B": (0, 3179, 8214, 6060),
    },
    "Figure2": {
        "A": (0, 0, 8084, 2260),
        "B": (0, 2260, 8084, 5680),
        "C": (0, 5680, 8084, 7793),
    },
    "Figure3": {
        "A": (0, 0, 1650, 1890),
        "B": (1650, 0, 3300, 1890),
        "C": (3300, 0, 4950, 1890),
        "D": (4950, 0, 6605, 1890),
        "E": (6605, 0, 8263, 1890),
        "F": (0, 1890, 8263, 4575),
        "G": (0, 4575, 8263, 10309),
    },
    "Figure4": {
        # The assembled figure has one shared legend above A-C. It is omitted
        # from the stand-alone panels so that no panel contains a legend fragment.
        "A": (0, 300, 3000, 3024),
        "B": (3000, 300, 5800, 3024),
        "C": (5800, 300, 8050, 3024),
        "D": (8050, 300, 10860, 3024),
    },
    "Figure5": {
        "A": (0, 0, 3920, 2410),
        "B": (3920, 0, 8109, 2410),
        "C": (0, 2410, 8109, 5060),
        "D": (0, 5060, 3550, 7267),
        "E": (3600, 5060, 8109, 7267),
    },
    "Figure6": {
        "A": (0, 0, 2360, 2300),
        "B": (2360, 0, 4650, 2300),
        "C": (4650, 0, 7018, 2300),
        "D": (0, 2300, 2360, 4624),
        "E": (2360, 2300, 4650, 4624),
        "F": (4650, 2300, 7018, 4624),
    },
    "Figure7": {
        "A": (0, 0, 6390, 3429),
        "B": (0, 3429, 6390, 6550),
    },
}

def split_pdf(source: Path, output: Path, crop: tuple[int, int, int, int],
              image_size: tuple[int, int]) -> tuple[float, float]:
    reader = PdfReader(str(source))
    if len(reader.pages) != 1:
        raise ValueError(f"Expected one PDF page: {source}")
    page = copy.deepcopy(reader.pages[0])
    page.transfer_rotation_to_content()

    image_width, image_height = image_size
    pdf_width = float(page.mediabox.width)
    pdf_height = float(page.mediabox.height)
    left, top, right, bottom = crop
    x0 = left * pdf_width / image_width
    x1 = right * pdf_width / image_width
    y0 = pdf_height - bottom * pdf_height / image_height
    y1 = pdf_height - top * pdf_height / image_height
    width = x1 - x0
    height = y1 - y0

    page.add_transformation(Transformation().translate(tx=-x0, ty=-y0))
    box = RectangleObject([0, 0, width, height])
    page.mediabox = box
    page.cropbox = box
    page.trimbox = box
    page.bleedbox = box
    page.artbox = box

    writer = PdfWriter()
    writer.add_page(page)
    with output.open("wb") as stream:
        writer.write(stream)
    return width, height


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = [
        "file\tsource_figure\tpanel\tformat\twidth\theight\tsha256"
    ]

    for figure, panels in CROPS.items():
        source_png = SOURCE / f"{figure}.png"
        source_pdf = SOURCE / f"{figure}.pdf"
        with Image.open(source_png) as assembled:
            assembled.load()
            image_size = assembled.size
            for panel, crop in panels.items():
                stem = f"{figure}{panel}"
                png_path = OUTPUT / f"{stem}.png"
                pdf_path = OUTPUT / f"{stem}.pdf"

                panel_image = assembled.crop(crop)
                panel_image.save(
                    png_path,
                    format="PNG",
                    dpi=(600, 600),
                    compress_level=6,
                )
                pdf_width, pdf_height = split_pdf(
                    source_pdf, pdf_path, crop, image_size
                )

                rows.append(
                    f"{png_path.name}\t{figure}\t{panel}\tPNG\t"
                    f"{panel_image.width}\t{panel_image.height}\t{sha256(png_path)}"
                )
                rows.append(
                    f"{pdf_path.name}\t{figure}\t{panel}\tPDF\t"
                    f"{pdf_width:.3f}\t{pdf_height:.3f}\t{sha256(pdf_path)}"
                )

    (OUTPUT / "MANIFEST.tsv").write_text(
        "\n".join(rows) + "\n", encoding="utf-8", newline="\n"
    )
    print(f"Wrote {len(rows) - 1} panel files to {OUTPUT}")


if __name__ == "__main__":
    main()
