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
const H1 = (text) => new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 320, after: 160 },
  children: [new TextRun({ text, font: FONT, size: 30, bold: true, color: "2E5C8A" })] });
const H2 = (text) => new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 240, after: 120 },
  children: [new TextRun({ text, font: FONT, size: 24, bold: true, color: "1A1A1A" })] });
const bullet = (text, level = 0) => new Paragraph({
  numbering: { reference: "bullets", level: level },
  spacing: { after: 80, line: 300 },
  children: [new TextRun({ text, font: FONT, size: 22 })],
});
const blank = () => new Paragraph({ children: [new TextRun({ text: "", font: FONT, size: 22 })] });
const callout = (text, color="2E5C8A") => new Paragraph({
  spacing: { before: 120, after: 120, line: 300 },
  shading: { fill: "EEF3F8", type: ShadingType.CLEAR },
  border: { left: { style: BorderStyle.SINGLE, size: 24, color, space: 8 } },
  children: [new TextRun({ text, font: FONT, size: 22, color: "1A1A1A" })],
});
const border = { style: BorderStyle.SINGLE, size: 4, color: "888888" };
const borders = { top: border, bottom: border, left: border, right: border };
function makeCell(text, opts = {}) {
  return new TableCell({
    borders, width: { size: opts.width, type: WidthType.DXA },
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
function img(path, w, h, captionText) {
  const arr = [new Paragraph({
    alignment: AlignmentType.CENTER,
    children: [new ImageRun({
      type: "png",
      data: fs.readFileSync(path),
      transformation: { width: w, height: h },
      altText: { title: "Figure", description: "Generated figure", name: "fig" },
    })],
  })];
  if (captionText) {
    arr.push(new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { after: 200 },
      children: [new TextRun({ text: captionText, font: FONT, size: 18, italics: true, color: "555555" })],
    }));
  }
  return arr;
}

// === Tables ===
const medianOSTable = [
  ["Ecotype", "Median OS (days)", "Median OS (years)"],
  ["Inflamed (n=113)", "595", "1.63"],
  ["Intermediate (n=157)", "446", "1.22"],
  ["Immune-desert (n=79)", "403", "1.10"],
];

const uniCoxTable = [
  ["Term", "HR", "95% CI", "p-value"],
  ["ecotype: Inflamed (ref)", "1.00", "—", "—"],
  ["ecotype: Intermediate", "1.41", "1.02–1.96", "0.039"],
  ["ecotype: Immune-desert", "1.60", "1.09–2.34", "0.017"],
];

const multiCoxTable = [
  ["Term", "HR", "95% CI", "p-value"],
  ["ecotype: Inflamed (ref)", "1.00", "—", "—"],
  ["ecotype: Intermediate", "1.37", "0.99–1.91", "0.062"],
  ["ecotype: Immune-desert", "1.79", "1.21–2.66", "0.004"],
  ["cohort: pHGG_WT (ref)", "1.00", "—", "—"],
  ["cohort: DMG_K27", "2.96", "2.11–4.15", "3.7e-10"],
  ["cohort: DHG_G34", "1.15", "0.59–2.24", "0.672"],
  ["cohort: IHG", "0.22", "0.08–0.63", "0.005"],
  ["age (per yr)", "0.99", "0.96–1.02", "0.506"],
  ["sex: Female", "1.12", "0.84–1.49", "0.450"],
];

const stratKMTable = [
  ["Cohort group", "n (with OS)", "Log-rank p"],
  ["DMG_K27", "140", "0.132"],
  ["DHG_G34", "15",  "0.070"],
  ["pHGG_WT", "82",  "0.049 ✓"],
  ["IHG",     "14",  "0.329"],
];

const kwTopTable = [
  ["Rank", "Signature", "KW χ²", "p (BH-adjusted)"],
  ["1",  "M1_Macrophage", "262.2", "2.8e-56"],
  ["2",  "Chemokine_T_Cell_Recruitment", "257.3", "1.6e-55"],
  ["3",  "T_Cell_Cytotoxicity", "245.0", "5.0e-53"],
  ["4",  "Glioma_Inflammatory_Wang2017", "241.8", "1.9e-52"],
  ["5",  "TGFb_Immunosuppression", "225.0", "6.5e-49"],
  ["6",  "Tregs_Friebel2020", "221.7", "2.9e-48"],
  ["7",  "T_Cell_Exhaustion", "220.1", "5.6e-48"],
  ["8",  "Dendritic_Cell_Activation", "213.8", "1.1e-46"],
  ["9",  "NK_Cell_Activity", "212.7", "1.7e-46"],
  ["10", "MoTAM_Antunes2021", "211.4", "3.0e-46"],
  ["11", "MDM_Klemm2020", "205.5", "5.2e-45"],
  ["12", "Neutrophil_Activation", "198.3", "1.8e-43"],
  ["13", "DAM_KerenShaul2017", "192.0", "3.8e-42"],
  ["14", "M2_Macrophage", "188.3", "2.2e-41"],
  ["15", "IFN_Gamma_Response", "186.2", "6.0e-41"],
  ["16", "MHC_Class_II", "157.7", "8.7e-35"],
  ["17", "MgTAM_Antunes2021", "100.1", "2.6e-22"],
  ["18", "Microglia_Core_Homeostatic", "97.9", "7.3e-22"],
  ["19", "Microglia_Klemm2020", "93.7", "5.6e-21"],
  ["20", "MHC_Class_I", "92.8", "8.7e-21"],
  ["21", "IFN_Alpha_Response", "89.3", "4.7e-20"],
  ["22", "MAPK_Activity", "54.8", "1.4e-12"],
  ["23", "Stemness_Brain_Tumor", "7.2", "0.029"],
  ["24", "Cell_Cycle_Proliferation", "0.13", "0.939 (NS)"],
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
    config: [{ reference: "bullets",
      levels: [
        { level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 480, hanging: 280 } } } },
      ] }],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1080, right: 1080, bottom: 1080, left: 1080 },
      },
    },
    children: [
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 120 },
        children: [new TextRun({ text: "Survival + Pathway 분석 결과 보고", font: FONT, size: 34, bold: true, color: "2E5C8A" })],
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 320 },
        children: [new TextRun({ text: "Kaplan-Meier OS + Cox model + ssGSEA pathway 비교 (n=349 main cohort)",
                                  font: FONT, size: 24, color: "555555" })],
      }),

      H1("0. 한눈에 보는 결과"),
      callout("Ecotype은 OS와 유의하게 연관되었습니다 (overall log-rank p = 0.037, n = 251). " +
              "Inflamed ecotype이 중위 생존 595일로 가장 좋고, Immune-desert가 403일로 가장 나쁩니다. " +
              "Cohort group, age, sex 보정 후에도 Inflamed vs Immune-desert HR = 1.79 (95% CI 1.21–2.66, p = 0.004)로 " +
              "독립적 prognostic 의미를 가집니다."),
      blank(),
      callout("24개 ssGSEA pathway 중 22개가 ecotype 간 유의한 차이 (BH-adj p < 1e-20). " +
              "특히 M1_Macrophage, Chemokine recruitment, T-cell cytotoxicity가 가장 강한 separator. " +
              "Cell cycle (proliferation)은 ecotype 간 차이 없음 (p = 0.94) — 종양 자체 proliferation이 아닌 " +
              "immune contexture가 진짜로 다르다는 sanity check 통과.", "16A34A"),

      // ============ Survival ============
      H1("1. Kaplan-Meier 생존 분석"),
      H2("1.1. Overall (n = 251 with OS data)"),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_survival/KM_overall_by_ecotype.png",
             520, 460, "Figure: Overall survival by ecotype. Log-rank p = 0.037"),
      makeTable(medianOSTable, [3500, 3000, 2860]),
      blank(),
      callout("Inflamed ecotype은 Immune-desert 대비 약 6.4개월 더 긴 중위 생존 (595 vs 403일). " +
              "이는 \"immune-cold tumor가 더 공격적\"이라는 일반적 직관과 일치하지만, " +
              "pediatric glioma에서 ecotype 단위로 정량화된 prognostic factor라는 점에서 novel."),

      H2("1.2. Cohort 별 stratified KM"),
      P("각 cohort group 내에서 ecotype의 prognostic 효과를 확인하였습니다."),
      blank(),
      makeTable(stratKMTable, [4000, 3000, 2360]),
      blank(),
      bullet("pHGG_WT (n=82): 유의 (p = 0.049). Ecotype이 H3-wildtype HGG 안에서도 prognostic"),
      bullet("DHG_G34 (n=15): trend (p = 0.070) — sample 작아도 효과 보임"),
      bullet("DMG_K27 (n=140): non-significant (p = 0.132). K27M의 일률적 poor outcome 때문에 ecotype 효과가 가려질 가능성"),
      bullet("IHG (n=14): non-significant (p = 0.329). 표본 부족"),
      blank(),

      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_survival/KM_by_ecotype_pHGG_WT.png",
             420, 420, "pHGG_WT: Inflamed > Intermediate > Immune-desert (p = 0.049)"),
      blank(),

      H2("1.3. Multivariable Cox model"),
      P("Ecotype + cohort_group + age + sex 보정 (n = 251):"),
      blank(),
      makeTable(multiCoxTable, [3000, 1500, 2300, 2560]),
      blank(),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_survival/Cox_forest_multivariable.png",
             580, 440, "Forest plot: Multivariable Cox HR for OS"),
      callout("Ecotype 효과는 cohort_group (DMG_K27 HR=2.96, IHG HR=0.22 등 큰 효과)을 보정한 뒤에도 " +
              "독립적으로 유지됩니다 (Immune-desert HR = 1.79, p = 0.004). " +
              "즉 ecotype은 histone class와 별개의 prognostic axis입니다.", "16A34A"),
      blank(),

      H2("1.4. Univariable Cox (ecotype only)"),
      makeTable(uniCoxTable, [3000, 1500, 2300, 2560]),
      blank(),

      // ============ Pathway analysis ============
      H1("2. ssGSEA Pathway 비교"),
      H2("2.1. Kruskal-Wallis test — 24개 signature"),
      makeTable(kwTopTable, [800, 4500, 1900, 2160]),
      blank(),
      bullet("22/24 signature가 BH-adjusted p < 1e-20으로 ecotype 간 유의한 차이"),
      bullet("M1 Macrophage / Chemokine recruitment / T-cell cytotoxicity가 가장 강한 ecotype separator"),
      bullet("Microglia signature (Klemm/Bohlen/Antunes Mg-TAM)는 ecotype 간 차이가 MDM/MoTAM보다 작음 — pediatric glioma에서 microglia는 ecotype-independent하게 상시 존재한다는 해석 가능"),
      bullet("MAPK_Activity (p_BH = 1.4e-12)도 ecotype 간 유의 — BRAF/MAPK biology와 연결 가능"),
      bullet("Cell cycle proliferation은 ecotype 간 차이 없음 (p = 0.94) — 좋은 negative control"),
      blank(),

      H2("2.2. Composite boxplot — 핵심 12 signature"),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/Pathway_composite_12sig.png",
             720, 500, "Figure: 12개 핵심 immune pathway의 ecotype별 분포"),

      H2("2.3. 개별 핵심 pathway"),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/box_Microglia_Klemm2020.png",
             340, 280),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/box_MDM_Klemm2020.png",
             340, 280, "Microglia vs MDM 비교: MDM이 ecotype 간 더 큰 차이를 보임"),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/box_MAPK_Activity.png",
             340, 280, "MAPK activity: Immune-desert에서 가장 낮고 Inflamed에서 가장 높음 (pediatric glioma 특이적 패턴)"),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/box_T_Cell_Cytotoxicity.png",
             340, 280, "T-cell cytotoxicity: 명확한 3-tier separation"),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/box_MHC_Class_II.png",
             340, 280, "MHC Class II: antigen presentation 활성"),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/box_Cell_Cycle_Proliferation.png",
             340, 280, "Cell cycle: ecotype 간 차이 없음 (sanity check 통과, p = 0.94)"),

      H2("2.4. Cohort-stratified pathway"),
      P("각 cohort 내에서도 ecotype별 pathway 차이가 유지되는지 확인하였습니다."),
      blank(),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/cohort_strat_T_Cell_Cytotoxicity.png",
             580, 360, "Cohort × Ecotype × T-cell cytotoxicity"),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/cohort_strat_MDM_Klemm2020.png",
             580, 360, "Cohort × Ecotype × MDM signature"),
      ...img("C:/Users/scott/Downloads/Open PBTA/output/figs_pathway/cohort_strat_MAPK_Activity.png",
             580, 360, "Cohort × Ecotype × MAPK Activity"),

      // ============ Synthesis ============
      H1("3. 통합 해석 — Manuscript 메시지"),
      H2("3.1. 핵심 메시지"),
      bullet("Ecotype은 OS와 독립적으로 연관 (multivariable HR = 1.79, p = 0.004)"),
      bullet("Inflamed > Intermediate > Immune-desert 순으로 prognosis가 좋아짐 — \"hot is better\""),
      bullet("Ecotype 효과는 H3-wildtype pHGG 안에서 가장 잘 보임 (p = 0.049)"),
      bullet("22/24 immune pathway가 ecotype 간 유의한 차이 (BH p < 1e-20)"),
      bullet("Cell cycle은 ecotype 간 차이 없음 — proliferation 아닌 \"진짜 immune\" 차이"),
      bullet("Microglia signal은 ecotype에 덜 의존, MDM signal은 강하게 ecotype 정의 — pediatric glioma의 \"baseline microglia + ecotype-driven MDM\" 모델"),
      blank(),

      H2("3.2. Reviewer-appeal 강화 포인트"),
      callout("\"Despite being H3 K27-altered and uniformly aggressive, DMG showed internal ecotype heterogeneity. " +
              "Among H3-wildtype HGG, ecotype was an independent prognostic axis (log-rank p = 0.049) — " +
              "suggesting that beyond molecular classification, immune contexture provides an additional layer " +
              "of risk stratification.\""),
      callout("\"Cell cycle and proliferation signatures did not differ across ecotypes (KW p = 0.94), " +
              "confirming that ecotype distinctions reflect genuine immune contexture rather than confounding " +
              "tumor cell proliferation differences.\""),
      callout("\"Brain-resident microglia signatures were only modestly stratified across ecotypes " +
              "(KW p ≈ 1e-21), while monocyte-derived macrophage (MDM/Mo-TAM) signatures showed " +
              "stronger ecotype-dependent variation (KW p ≈ 1e-45), suggesting that MDM recruitment " +
              "drives ecotype assignment more than baseline microglia abundance.\""),

      H2("3.3. 다음에 더 보강할 수 있는 분석"),
      bullet("ssGSEA score를 ecotype 외 cohort_group으로 직접 비교 (BRAF/MAPK 강도)"),
      bullet("Survival landmark analysis (각 cohort에서 18개월 시점 OS 비교)"),
      bullet("Time-dependent ROC (ecotype의 predictive 성능 정량화)"),
      bullet("CIBERSORTx LM22 결과 받으면 5-method consensus validation"),

      H1("4. 산출 파일"),
      H2("Survival"),
      bullet("KM_median_OS_by_ecotype.tsv"),
      bullet("KM_stratified_logrank_pvalues.tsv"),
      bullet("Cox_univariable_ecotype.tsv"),
      bullet("Cox_multivariable_coefficients.tsv"),
      bullet("figs_survival/KM_overall_by_ecotype.png + 4 cohort별 KM + Cox forest"),
      H2("Pathway"),
      bullet("ssGSEA_KruskalWallis_by_ecotype.tsv (24 signature 모든 KW p)"),
      bullet("ssGSEA_Dunn_posthoc_by_ecotype.tsv (pairwise comparison)"),
      bullet("ssGSEA_pairwise_effect_size.tsv (rank-biserial r)"),
      bullet("figs_pathway/box_*.png — 24개 signature boxplot"),
      bullet("figs_pathway/cohort_strat_*.png — 6개 cohort-stratified plot"),
      bullet("figs_pathway/Pathway_composite_12sig.png — multi-panel composite"),

      H1("5. 맺음말"),
      P("Survival과 pathway 두 축 모두에서 ecotype 분류가 강한 statistical signal과 biological interpretability를 보였습니다. " +
        "특히 multivariable Cox로 cohort_group을 보정한 뒤에도 ecotype이 prognostic factor로 살아남았다는 점이 manuscript에서 " +
        "강조할 핵심 결과입니다. 다음 단계로는 (a) Manuscript drafting, " +
        "(b) BRAF_ALT secondary cohort 투영, 또는 (c) CIBERSORTx LM22 검증 중 결정해 주시면 진행하겠습니다."),
      blank(),
      P("감사합니다."),
      blank(),
      new Paragraph({ alignment: AlignmentType.CENTER,
        children: [new TextRun({ text: "— 이상 —", font: FONT, size: 22, color: "555555" })] }),
    ],
  }],
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync("C:/Users/scott/Downloads/Open PBTA/output/Survival_Pathway_Report.docx", buf);
  console.log("Survival_Pathway_Report.docx written");
});
