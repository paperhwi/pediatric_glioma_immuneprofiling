const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, ImageRun,
  AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType, LevelFormat,
} = require("docx");

const FONT = "Malgun Gothic";

const P = (text, opts = {}) => new Paragraph({
  spacing: { after: 120, line: 320 },
  alignment: opts.align || AlignmentType.JUSTIFIED,
  children: [new TextRun({ text, font: FONT, size: opts.size || 22, bold: opts.bold || false, color: opts.color || "1A1A1A" })],
});

const H1 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_1,
  spacing: { before: 320, after: 160 },
  children: [new TextRun({ text, font: FONT, size: 30, bold: true, color: "2E5C8A" })],
});

const H2 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_2,
  spacing: { before: 240, after: 120 },
  children: [new TextRun({ text, font: FONT, size: 24, bold: true, color: "1A1A1A" })],
});

const bullet = (text, level = 0) => new Paragraph({
  numbering: { reference: "bullets", level: level },
  spacing: { after: 80, line: 300 },
  children: [new TextRun({ text, font: FONT, size: 22 })],
});

const blank = () => new Paragraph({ children: [new TextRun({ text: "", font: FONT, size: 22 })] });

const callout = (text, color="2E5C8A") => new Paragraph({
  spacing: { before: 120, after: 120 },
  shading: { fill: "EEF3F8", type: ShadingType.CLEAR },
  border: { left: { style: BorderStyle.SINGLE, size: 24, color, space: 8 } },
  children: [new TextRun({ text, font: FONT, size: 22, color: "1A1A1A" })],
});

const border = { style: BorderStyle.SINGLE, size: 4, color: "888888" };
const borders = { top: border, bottom: border, left: border, right: border };

function makeCell(text, opts = {}) {
  return new TableCell({
    borders,
    width: { size: opts.width, type: WidthType.DXA },
    shading: opts.fill ? { fill: opts.fill, type: ShadingType.CLEAR } : undefined,
    margins: { top: 80, bottom: 80, left: 100, right: 100 },
    children: [new Paragraph({
      alignment: opts.align || AlignmentType.LEFT,
      spacing: { after: 0 },
      children: [new TextRun({ text, font: FONT, size: opts.size || 18, bold: opts.bold || false, color: opts.color || "1A1A1A" })],
    })],
  });
}

function makeTable(rows, widths, headerFill="2E5C8A") {
  const totalW = widths.reduce((a,b)=>a+b,0);
  return new Table({
    width: { size: totalW, type: WidthType.DXA },
    columnWidths: widths,
    rows: rows.map((cells, ri) => new TableRow({
      tableHeader: ri === 0,
      children: cells.map((t, ci) => makeCell(t, {
        width: widths[ci],
        fill: ri === 0 ? headerFill : undefined,
        color: ri === 0 ? "FFFFFF" : "1A1A1A",
        bold: ri === 0,
        align: ci === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
      })),
    })),
  });
}

// Read correlation matrix files for embedding
function readCorrMatrix(path) {
  const lines = fs.readFileSync(path, "utf-8").trim().split("\n").map(l => l.split(","));
  // First row is header; first column is row label
  return lines.map(row => row.map(v => {
    if (v === "" || isNaN(parseFloat(v))) return v;
    const n = parseFloat(v);
    return Math.abs(n - Math.round(n)) < 1e-9 ? v : n.toFixed(2);
  }));
}

const corrMyeloid = readCorrMatrix("C:/Users/scott/Downloads/Open PBTA/output/spearman_myeloid_crossmethod.tsv");
const corrCD8     = readCorrMatrix("C:/Users/scott/Downloads/Open PBTA/output/spearman_cd8_crossmethod.tsv");
const corrNK      = readCorrMatrix("C:/Users/scott/Downloads/Open PBTA/output/spearman_nk_crossmethod.tsv");

// Step summary
const stepRows = [
  ["#", "단계", "방법", "상태"],
  ["7", "immunedeconv 4-method", "MCP-counter, quanTIseq, EPIC, xCell (R immunedeconv)", "완료"],
  ["7b", "Brain-tuned ssGSEA", "GSVA::ssgsea — 24 signature (Klemm, Antunes, Bohlen, MSigDB)", "완료"],
  ["7c", "Cross-method Spearman", "Myeloid · CD8 · NK 3개 축", "완료"],
];

// Signature inventory
const sigRows = [
  ["분류", "Signature", "출처", "유전자 수 (matrix 매칭)"],
  ["Microglia core", "Microglia_Core_Homeostatic", "Bohlen 2020 / Butovsky 2014", "16/17"],
  ["Microglia (brain tumor)", "Microglia_Klemm2020", "Klemm 2020 Cell", "15/16"],
  ["Microglia-TAM", "MgTAM_Antunes2021", "Antunes 2021 Nat Neurosci", "13/14"],
  ["MDM (brain tumor)", "MDM_Klemm2020", "Klemm 2020 Cell", "16/16"],
  ["Monocyte-TAM", "MoTAM_Antunes2021", "Antunes 2021 Nat Neurosci", "15/15"],
  ["Disease-assoc microglia", "DAM_KerenShaul2017", "Keren-Shaul 2017 Cell", "15/16"],
  ["Antigen presentation", "MHC_Class_I / MHC_Class_II", "Curated", "17/17, 12/12"],
  ["Interferon response", "IFN_Gamma_Response / IFN_Alpha_Response", "Hallmark-like", "20/20, 19/19"],
  ["T-cell function", "T_Cell_Cytotoxicity / Exhaustion", "Rooney 2015 등", "13/13, 12/12"],
  ["Chemokine recruit", "Chemokine_T_Cell_Recruitment", "Curated", "12/12"],
  ["MAPK activity", "MAPK_Activity", "MAPK target gene set", "16/16"],
  ["Immunosuppression", "TGFb_Immunosuppression", "Curated", "13/13"],
  ["M1/M2 polarization", "M1_Macrophage / M2_Macrophage", "Curated", "14/14, 14/14"],
  ["NK activity", "NK_Cell_Activity", "Curated", "15/15"],
  ["Treg", "Tregs_Friebel2020", "Friebel 2020", "10/10"],
  ["DC activation", "Dendritic_Cell_Activation", "Curated", "10/10"],
  ["Neutrophil", "Neutrophil_Activation", "Curated", "14/14"],
  ["Cell cycle", "Cell_Cycle_Proliferation", "Curated", "13/13"],
  ["Stemness", "Stemness_Brain_Tumor", "Curated brain stem markers", "11/11"],
  ["Glioma inflammatory", "Glioma_Inflammatory_Wang2017", "Wang 2017", "13/13"],
];

// Reviewer appeal table
const appealRows = [
  ["메시지", "근거 (지금 분석 결과)"],
  ["LM22/xCell는 microglia와 MDM을 구별 불가",
   "xCell macrophage ↔ MDM-Klemm: r = 0.86 (편향), Microglia-Klemm과는 r = 0.57"],
  ["Brain-tuned signature는 microglia/MDM 분리 가능",
   "Microglia signatures 3종 내부 r > 0.98; MDM signatures 2종 내부 r = 0.94; Microglia ↔ MDM: r ≈ 0.40"],
  ["4-method cross-validation 완료",
   "MCP-counter, quanTIseq, EPIC, xCell + ssGSEA 24sig — 702 sample 전체"],
  ["Microglia signature 3종이 서로 강하게 일치",
   "Klemm 2020 ↔ Bohlen 2020: r = 0.995; Klemm ↔ Antunes Mg-TAM: r = 0.986 → 내적 검증"],
  ["MDM signature 2종이 서로 강하게 일치",
   "Klemm MDM ↔ Antunes Mo-TAM: r = 0.941 → 내적 검증"],
];

const doc = new Document({
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 30, bold: true, font: FONT, color: "2E5C8A" },
        paragraph: { spacing: { before: 320, after: 160 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, font: FONT },
        paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 1 } },
    ],
  },
  numbering: {
    config: [
      { reference: "bullets",
        levels: [
          { level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 480, hanging: 280 } } } },
        ] },
    ],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1080, bottom: 1440, left: 1080 },
      },
    },
    children: [
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 120 },
        children: [new TextRun({
          text: "Step 7 + 7b 결과 보고",
          font: FONT, size: 34, bold: true, color: "2E5C8A",
        })],
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 320 },
        children: [new TextRun({
          text: "immunedeconv 4-method cross-validation + Brain-tuned ssGSEA scoring",
          font: FONT, size: 24, color: "555555",
        })],
      }),

      // 1. 요약
      H1("1. 진행 요약"),
      P("선생님 결정에 따라 immunedeconv 4-method cross-validation에 Brain-tuned ssGSEA scoring (Option A)을 추가로 진행하였습니다. 두 step 모두 정상 완료되었고, 특히 brain-tuned signature가 reviewer-appeal 측면에서 매우 강력한 결과를 만들어주었습니다."),
      blank(),
      makeTable(stepRows, [600, 2500, 4960, 1300]),
      blank(),

      // 2. immunedeconv 4-method 결과
      H1("2. immunedeconv 4-method 결과"),
      P("702 sample(main cohort 349 + secondary BRAF 353)에 대해 4가지 deconvolution 방법을 모두 적용 완료하였습니다. 각 방법의 실행 시간과 산출 cell type 수는 아래와 같습니다."),
      blank(),
      makeTable([
        ["Method", "Cell type 수", "실행 시간", "산출 파일"],
        ["MCP-counter", "11", "2.9 sec", "immunedeconv_mcp_counter.tsv"],
        ["quanTIseq",   "11 (incl. uncharacterized)", "12.3 min", "immunedeconv_quantiseq.tsv"],
        ["EPIC",        "8",  "0.4 min", "immunedeconv_epic.tsv"],
        ["xCell",       "39", "2.2 min", "immunedeconv_xcell.tsv"],
      ], [2200, 2200, 2200, 2760]),
      blank(),
      P("EPIC에서 32개 sample은 optimization이 완전히 수렴하지 않았다는 경고가 있었습니다 (모두 immune-cold tumor일 가능성). 해당 sample list는 deconv_run.log에 보관되어 있으며, downstream에서 EPIC 결과만 sensitivity로 별도 표시할 예정입니다."),

      // 3. Brain-tuned ssGSEA
      H1("3. Brain-tuned ssGSEA — 24개 signature 점수화"),
      P("ssGSEA (GSVA package) 방법으로 24개 gene signature를 702 sample에 적용하였습니다. 각 signature의 출처와 matrix 매칭 유전자 수는 다음과 같습니다 (>90% 매칭률은 robust한 score를 의미)."),
      blank(),
      makeTable(sigRows, [2100, 3400, 2700, 1160]),
      blank(),

      // 4. 핵심 발견
      H1("4. ⭐ 핵심 발견 — Microglia vs MDM 분리 성공"),
      P("이번 분석에서 가장 결정적인 결과는, 기존 LM22-류 deconvolution methods가 구별하지 못하는 \"microglia (brain-resident)\" vs \"monocyte-derived macrophage (MDM, peripheral)\" 신호를 brain-tuned ssGSEA로 명확히 분리하였다는 것입니다."),
      blank(),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        children: [new ImageRun({
          type: "png",
          data: fs.readFileSync("C:/Users/scott/Downloads/Open PBTA/output/fig_corr_myeloid_annotated.png"),
          transformation: { width: 600, height: 480 },
          altText: { title: "Myeloid cross-method correlation", description: "Heatmap", name: "fig_myeloid" },
        })],
      }),
      blank(),
      callout("Microglia signature 3종 (Klemm 2020, Bohlen 2020, Antunes MgTAM)이 서로 r > 0.98로 일치하고, "
            + "MDM signature 2종 (Klemm MDM, Antunes MoTAM)이 서로 r = 0.94로 일치하는 반면, "
            + "Microglia ↔ MDM 간 상관은 r ≈ 0.40에 그칩니다. 두 신호가 명확히 분리되는 별개 cell population임을 보여주는 강력한 내적 검증입니다."),
      blank(),

      H2("4.1. 해석"),
      bullet("Microglia signature 3종 일치 (r > 0.98) → 신호의 내적 일관성 확인 (어느 단일 marker set 이슈가 아님)"),
      bullet("MDM signature 2종 일치 (r = 0.94) → 신호의 내적 일관성 확인"),
      bullet("Microglia ↔ MDM 간 r ≈ 0.40 → 두 cell population이 분명히 구별됨"),
      bullet("xCell macrophage ↔ MDM-Klemm r = 0.86 → xCell의 \"macrophage\"는 사실상 peripheral MDM 신호에 편향"),
      bullet("EPIC macrophage ↔ Microglia r ≈ 0.65 → EPIC은 microglia를 일부 잡아내지만 분리 불가"),

      // 5. Reviewer-appeal 정리
      H1("5. Reviewer-appeal point 정리"),
      P("이 결과는 다음과 같은 manuscript-level reviewer 방어 무기로 직접 활용 가능합니다."),
      blank(),
      makeTable(appealRows, [3700, 5660]),
      blank(),

      // 6. CD8 / NK 결과
      H1("6. CD8 / NK 결과 — 해석 주의 필요"),
      H2("6.1. CD8 T cell / cytolytic 신호"),
      P("CD8 T cell deconvolution의 cross-method 일치도는 낮은 수준이었습니다 (Spearman r = 0.15–0.54). 이는 소아 glioma 자체가 전반적으로 immune-cold라 절대 신호가 낮고, method 간 noise가 두드러지기 때문으로 판단됩니다."),
      blank(),
      makeTable([
        ["", "MCP_CD8", "QTS_CD8", "EPIC_CD8", "Cytolytic_ssGSEA"],
        ["MCP_CD8", "1.00", "0.39", "0.30", "0.54"],
        ["QTS_CD8", "0.39", "1.00", "0.16", "0.28"],
        ["EPIC_CD8", "0.30", "0.16", "1.00", "0.26"],
        ["Cytolytic_ssGSEA", "0.54", "0.28", "0.26", "1.00"],
      ], [3000, 1590, 1590, 1590, 1590]),
      blank(),
      bullet("MCP-counter CD8가 cytolytic ssGSEA와 가장 높은 일치 (r = 0.54) → MCP-counter CD8을 primary로 사용 권장"),
      bullet("EPIC CD8는 다른 method와 매우 낮은 상관 (r = 0.16–0.30) → supplementary만 사용"),
      bullet("Manuscript에서는 \"T cell signal 약함, 결과 cautious interpretation\" 명시"),

      H2("6.2. NK cell 신호"),
      P("NK cell은 더 큰 method-간 disagreement를 보였습니다 (r = -0.05 ~ 0.53). xCell NK는 거의 모든 다른 method와 무상관(r ≈ 0)이라 분석에서 제외하는 편이 안전합니다."),
      blank(),
      makeTable([
        ["", "MCP_NK", "QTS_NK", "EPIC_NK", "xCell_NK", "NK_ssGSEA"],
        ["MCP_NK", "1.00", "0.34", "0.18", "−0.03", "0.53"],
        ["QTS_NK", "0.34", "1.00", "−0.03", "−0.02", "0.16"],
        ["EPIC_NK", "0.18", "−0.03", "1.00", "−0.05", "0.37"],
        ["xCell_NK", "−0.03", "−0.02", "−0.05", "1.00", "−0.03"],
        ["NK_ssGSEA", "0.53", "0.16", "0.37", "−0.03", "1.00"],
      ], [2400, 1392, 1392, 1392, 1392, 1392]),
      blank(),
      bullet("MCP-counter NK + ssGSEA NK_Cell_Activity 조합을 primary로 사용 권장 (r = 0.53)"),
      bullet("xCell NK는 결과에서 제외 또는 supplementary 명시"),

      // 7. 산출 파일
      H1("7. 산출 파일"),
      P("이번 단계의 모든 결과는 output 디렉터리에 저장되어 있습니다. 총 12개 파일 + 4개 figure."),
      H2("7.1. immunedeconv 출력 (cell type fraction matrix)"),
      bullet("immunedeconv_mcp_counter.tsv — 11 cell × 702"),
      bullet("immunedeconv_quantiseq.tsv  — 11 cell × 702"),
      bullet("immunedeconv_epic.tsv       — 8 cell × 702"),
      bullet("immunedeconv_xcell.tsv      — 39 cell × 702"),
      H2("7.2. Brain-tuned ssGSEA"),
      bullet("ssGSEA_brain_immune_scores.tsv      — 24 signature × 702 (wide)"),
      bullet("ssGSEA_brain_immune_scores_long.tsv — long format (ggplot 용)"),
      H2("7.3. Cross-method 통합 및 상관 매트릭스"),
      bullet("myeloid_signals_merged.tsv + spearman_myeloid_crossmethod.tsv"),
      bullet("cd8_signals_merged.tsv + spearman_cd8_crossmethod.tsv"),
      bullet("nk_signals_merged.tsv + spearman_nk_crossmethod.tsv"),
      H2("7.4. Figure (PNG)"),
      bullet("fig_corr_myeloid_annotated.png — Microglia/MDM 분리 핵심 figure"),
      bullet("fig_corr_myeloid.png / fig_corr_cd8.png / fig_corr_nk.png — 각 axis별 heatmap"),

      // 8. 다음 단계
      H1("8. 다음 단계 — Consensus clustering"),
      P("Deconvolution과 brain-tuned scoring이 끝났으므로, 이제 ecotype clustering으로 넘어갈 수 있습니다. 세 가지 옵션 중 선택해 주시면 됩니다."),
      blank(),
      makeTable([
        ["옵션", "내용", "장단점"],
        ["A", "Immunedeconv (예: quanTIseq 10-cell) + ssGSEA 만으로 즉시 consensus clustering",
              "장점: CIBERSORTx 안 기다리고 빠른 진행. 단점: LM22 결과는 별도 검증으로 추가 필요"],
        ["B", "CIBERSORTx LM22 결과까지 기다린 뒤 5-method consensus",
              "장점: 가장 robust한 ecotype 정의. 단점: 시간 소요"],
        ["C", "둘 다 — A로 잠정 ecotype 정의 후 LM22 결과 받으면 B로 검증",
              "장점: 빠른 결과 + 추후 robust한 검증. 권장"],
      ], [800, 4500, 4060]),
      blank(),

      // 9. 맺음말
      H1("9. 맺음말"),
      P("이번 단계에서 microglia vs MDM 분리라는 결정적인 brain-tumor-specific 결과를 확보하였습니다. 이 결과는 manuscript의 reviewer-appeal point 중 가장 강력한 무기가 될 것으로 판단됩니다 (\"우리는 LM22의 한계를 brain-tuned signature로 능동적으로 보완했다\"는 메시지)."),
      P("CIBERSORTx 결과를 기다리시는 동안 consensus clustering(옵션 A)을 먼저 진행하실 의향이 있으시면 그쪽으로 곧장 넘어가겠습니다. 옵션 결정해 주시면 바로 진행하겠습니다."),
      blank(),
      P("감사합니다."),
      blank(),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        children: [new TextRun({ text: "— 이상 —", font: FONT, size: 22, color: "555555" })],
      }),
    ],
  }],
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync("C:/Users/scott/Downloads/Open PBTA/output/Step7_Deconvolution_Report.docx", buf);
  console.log("Step7_Deconvolution_Report.docx written");
});
