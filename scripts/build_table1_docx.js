const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType,
} = require("docx");

const FONT = "Malgun Gothic";

// Read Table 1
const rows1 = fs.readFileSync(
  "C:/Users/scott/Downloads/Open PBTA/output/Table1_cohort_characteristics.tsv", "utf-8"
).trim().split("\n").map(r => r.split("\t"));

const rows2 = fs.readFileSync(
  "C:/Users/scott/Downloads/Open PBTA/output/Table1_supp_BRAF_secondary.tsv", "utf-8"
).trim().split("\n").map(r => r.split("\t"));

const border = { style: BorderStyle.SINGLE, size: 4, color: "888888" };
const borders = { top: border, bottom: border, left: border, right: border };

function buildTable(rows, contentWidth = 9360) {
  const ncol = rows[0].length;
  const colWidth = Math.floor(contentWidth / ncol);
  const cols = new Array(ncol).fill(colWidth);
  cols[ncol - 1] = contentWidth - colWidth * (ncol - 1);

  const trows = rows.map((cells, ri) => new TableRow({
    tableHeader: ri === 0,
    children: cells.map((t, ci) => new TableCell({
      borders,
      width: { size: cols[ci], type: WidthType.DXA },
      shading: ri === 0
        ? { fill: "2E5C8A", type: ShadingType.CLEAR }
        : (ri === rows.length - 1 ? { fill: "F0F4F8", type: ShadingType.CLEAR } : undefined),
      margins: { top: 80, bottom: 80, left: 100, right: 100 },
      children: [new Paragraph({
        alignment: ci === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
        spacing: { after: 0 },
        children: [new TextRun({
          text: t,
          font: FONT,
          size: 18,
          bold: ri === 0 || ri === rows.length - 1,
          color: ri === 0 ? "FFFFFF" : "1A1A1A",
        })],
      })],
    })),
  }));

  return new Table({
    width: { size: contentWidth, type: WidthType.DXA },
    columnWidths: cols,
    rows: trows,
  });
}

// Landscape orientation for wide tables
const doc = new Document({
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840, orientation: "landscape" },
        margin: { top: 1080, right: 1080, bottom: 1080, left: 1080 },
      },
    },
    children: [
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 200 },
        children: [new TextRun({
          text: "Table 1. Clinical and molecular characteristics of pediatric glioma cohort",
          font: FONT, size: 28, bold: true, color: "2E5C8A",
        })],
      }),
      new Paragraph({
        spacing: { after: 200 },
        children: [new TextRun({
          text: "OpenPedCan v15 기반, independent-primary-plus 필터 + pHGG_WT age <21yr cap 적용.",
          font: FONT, size: 20, italics: true, color: "555555",
        })],
      }),
      buildTable(rows1, 13680),
      new Paragraph({ children: [new TextRun({ text: "", font: FONT })] }),
      new Paragraph({
        spacing: { before: 320, after: 160 },
        children: [new TextRun({
          text: "Supplementary — BRAF/RTK-altered secondary cohort",
          font: FONT, size: 24, bold: true, color: "2E5C8A",
        })],
      }),
      buildTable(rows2, 13680),
      new Paragraph({
        spacing: { before: 240 },
        children: [new TextRun({
          text: "Footnote: All counts based on independent-primary-plus filter (one biospecimen per patient). " +
                "Ambiguous = pineal, suprasellar, optic pathway, ventricles. " +
                "IUPAC nomenclature K28/G35 corresponds to historical K27/G34.",
          font: FONT, size: 18, italics: true, color: "666666",
        })],
      }),
    ],
  }],
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync("C:/Users/scott/Downloads/Open PBTA/output/Table1.docx", buf);
  console.log("Table1.docx written");
});
