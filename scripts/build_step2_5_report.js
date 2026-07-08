const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType,
  LevelFormat,
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

const codeBlock = (text) => new Paragraph({
  spacing: { after: 120, line: 280 },
  shading: { fill: "F2F4F7", type: ShadingType.CLEAR },
  children: [new TextRun({ text, font: "Consolas", size: 20, color: "1A1A1A" })],
});

const border = { style: BorderStyle.SINGLE, size: 4, color: "888888" };
const borders = { top: border, bottom: border, left: border, right: border };

function makeCell(text, opts = {}) {
  return new TableCell({
    borders,
    width: { size: opts.width, type: WidthType.DXA },
    shading: opts.fill ? { fill: opts.fill, type: ShadingType.CLEAR } : undefined,
    margins: { top: 80, bottom: 80, left: 120, right: 120 },
    children: [new Paragraph({
      alignment: opts.align || AlignmentType.LEFT,
      spacing: { after: 0 },
      children: [new TextRun({ text, font: FONT, size: opts.size || 20, bold: opts.bold || false, color: opts.color || "1A1A1A" })],
    })],
  });
}

function makeTable(rows, widths, headerFill = "2E5C8A") {
  const totalW = widths.reduce((a, b) => a + b, 0);
  return new Table({
    width: { size: totalW, type: WidthType.DXA },
    columnWidths: widths,
    rows: rows.map((cells, ri) => new TableRow({
      tableHeader: ri === 0,
      children: cells.map((t, ci) => makeCell(t, {
        width: widths[ci],
        fill: ri === 0 ? headerFill : (ri === rows.length - 1 && rows.length > 2 && cells[0].includes("Total") ? "F0F4F8" : undefined),
        color: ri === 0 ? "FFFFFF" : "1A1A1A",
        bold: ri === 0,
        align: ci === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
      })),
    })),
  });
}

// ============ Tables ============

// Step summary table
const stepRows = [
  ["#", "단계", "산출물", "상태"],
  ["1", "Independent filter + age cap + ambiguous loc 재분류", "cohort_main_final.tsv (n=349)", "완료"],
  ["2", "BRAF_ALT를 HGG/PXA/infant-fusion vs LGG로 분리", "cohort_BRAF_ALT_main.tsv, cohort_BRAF_ALT_LGG.tsv", "완료"],
  ["3", "Table 1 산출", "Table1_cohort_characteristics.tsv, Table1.docx", "완료"],
  ["4", "TPM matrix RDS 다운로드 (285 MB)", "gene-expression-rsem-tpm-collapsed.rds", "완료"],
  ["5", "Cohort sample subset → CIBERSORTx input 포맷 export", "tpm_for_cibersortx.tsv (180 MB)", "완료"],
];

// Cohort flow table
const flowRows = [
  ["필터 단계", "DMG_K27", "DHG_G34", "pHGG_WT", "IHG", "Main 합계"],
  ["1단계: PBTA + RNA-Seq + Tumor (cohort_main.tsv)", "225", "37", "192", "28", "482"],
  ["2단계: + Independent-primary-plus", "175", "31", "135", "18", "359"],
  ["3단계: + pHGG_WT age < 21 yr (최종)", "175", "31", "125", "18", "349"],
];

// Main cohort table (Table 1)
const mainRows = [
  ["변수", "DMG_K27", "DHG_G34", "pHGG_WT", "IHG", "Total"],
  ["n biospecimen", "175", "31", "125", "18", "349"],
  ["n unique patient", "175", "31", "125", "18", "349"],
  ["Age yr (median, IQR)", "7.9 (5.5–11.2)", "15.7 (13.5–17.9)", "9.5 (6.0–13.7)", "0.4 (0.1–0.9)", "8.8 (5.5–12.5)"],
  ["Male", "75 (42.9%)", "16 (51.6%)", "70 (56.0%)", "9 (50.0%)", "170 (48.7%)"],
  ["Female", "97 (55.4%)", "15 (48.4%)", "55 (44.0%)", "9 (50.0%)", "176 (50.4%)"],
  ["Midline", "118 (67.4%)", "0 (0.0%)", "22 (17.6%)", "0 (0.0%)", "140 (40.1%)"],
  ["Hemispheric", "3 (1.7%)", "24 (77.4%)", "62 (49.6%)", "15 (83.3%)", "104 (29.8%)"],
  ["Posterior fossa", "5 (2.9%)", "0 (0.0%)", "9 (7.2%)", "0 (0.0%)", "14 (4.0%)"],
  ["Ambiguous/Other", "49 (28.0%)", "7 (22.6%)", "32 (25.6%)", "3 (16.7%)", "91 (26.1%)"],
  ["BRAF/RTK fusion+", "31 (17.7%)", "4 (12.9%)", "17 (13.6%)", "18 (100.0%)", "70 (20.1%)"],
  ["OS data available", "140 (80.0%)", "15 (48.4%)", "82 (65.6%)", "14 (77.8%)", "251 (71.9%)"],
];

// BRAF secondary table
const brafRows = [
  ["Subgroup", "n biospecimen", "Age yr (median, IQR)", "BRAF/RTK fusion+", "OS 가능"],
  ["BRAF_ALT_main (HGG/PXA/infant-fusion)", "31", "—", "—", "—"],
  ["BRAF_ALT_LGG (LGG only)", "322", "—", "—", "—"],
];

// TPM extraction QC table
const tpmRows = [
  ["항목", "값"],
  ["원본 RDS shape", "60,325 genes × 4,124 samples"],
  ["입력 cohort biospecimen IDs", "702"],
  ["Expression matrix에 매칭된 sample 수", "702 / 702 (100%)"],
  ["NA 값", "0"],
  ["음수 값 (TPM 무결성)", "0"],
  ["중복 gene symbol", "0"],
  ["All-zero gene 제거", "4,917"],
  ["최종 matrix", "55,408 genes × 702 samples"],
  ["출력 파일 크기", "180 MB"],
];

// Deconvolution sample composition
const deconvRows = [
  ["Cohort group", "n samples in TPM matrix", "용도"],
  ["DMG_K27", "175", "Main analysis"],
  ["DHG_G34", "31", "Main analysis"],
  ["pHGG_WT", "125", "Main analysis"],
  ["IHG", "18", "Main analysis"],
  ["BRAF_ALT (main + LGG)", "353", "Secondary / supplementary"],
  ["Total", "702", "—"],
];

// ============ Document ============

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
          { level: 1, format: LevelFormat.BULLET, text: "◦", alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 880, hanging: 280 } } } },
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
      // Title
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 120 },
        children: [new TextRun({
          text: "Step 2–5 결과 보고",
          font: FONT, size: 34, bold: true, color: "2E5C8A",
        })],
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 320 },
        children: [new TextRun({
          text: "Cohort finalization · Table 1 · TPM matrix 추출 완료",
          font: FONT, size: 24, color: "555555",
        })],
      }),

      // ============ 1. 진행 요약 ============
      H1("1. 진행 요약"),
      P("선생님 회신 주신 결정 사항을 반영하여 Step 2부터 Step 5까지 일괄 진행하였습니다. 모든 단계가 정상적으로 마무리되었고, CIBERSORTx 업로드에 사용할 TPM matrix(180 MB)가 준비되었습니다. 다음은 단계별 산출물 요약입니다."),
      blank(),
      makeTable(stepRows, [600, 3700, 3500, 1100]),
      blank(),

      // ============ 2. 필터 적용 흐름 ============
      H1("2. Cohort 필터 적용 흐름"),
      P("Main cohort는 회신 주신 기준에 따라 다음과 같이 축소되었습니다."),
      blank(),
      makeTable(flowRows, [3500, 1280, 1280, 1280, 1280, 1280]),
      blank(),
      bullet("초기 482건 → independent-primary-plus 필터로 환자당 1 sample만 남겨 359건"),
      bullet("pHGG_WT 중 age ≥ 21yr (성인) 10건 제외하여 최종 349건"),
      bullet("최종 main cohort의 n biospecimen = n unique patient = 349 (1:1 매칭)"),

      // ============ 3. Main cohort 특성 (Table 1) ============
      H1("3. Main cohort 특성 (Table 1)"),
      P("선생님 회신에 따라 pineal, suprasellar, optic pathway, ventricles 등 midline-hemispheric 경계가 모호한 sample은 별도 \"Ambiguous/Other\" category로 분리하였고, location association 분석에서는 제외할 예정입니다."),
      blank(),
      makeTable(mainRows, [2500, 1372, 1372, 1372, 1372, 1372]),
      blank(),
      H2("3.1. 생물학적 타당성 확인"),
      bullet("DMG_K27: 67%가 midline 종양, 중위 연령 7.9세 — 알려진 DMG 분포와 일치"),
      bullet("DHG_G34: 77%가 hemispheric, 중위 연령 15.7세 — H3 G34 entity 특성에 부합"),
      bullet("IHG: 100% infant, 83% hemispheric, 100% BRAF/RTK fusion+ — operational definition이 의도대로 작동"),
      bullet("pHGG_WT: midline/hemispheric 혼재 (각 18%/50%) — H3-wildtype HGG의 이질적 분포 반영"),
      bullet("Ambiguous category 비율은 cohort별 17–28%로 일정한 수준 — pineal/suprasellar 등이 일관되게 분리됨"),

      // ============ 4. Secondary BRAF cohort ============
      H1("4. Secondary BRAF cohort 분리"),
      P("BRAF_ALT cohort는 회신 주신 대로 \"HGG/PXA/infant-fusion 중심\"과 \"LGG only\" 두 그룹으로 분리하였습니다. Main interpretation에는 BRAF_ALT_main만 사용하고, LGG는 별도 exploratory/supplementary로 분리 표시할 예정입니다."),
      blank(),
      makeTable(brafRows, [4360, 1500, 1500, 1500, 500]),
      blank(),
      P("BRAF_ALT_main의 31건은 sample 수가 적어 단독 clustering보다는 main cohort와 함께 deconvolution한 뒤 ecotype 분포를 비교하는 방향으로 사용할 예정입니다."),

      // ============ 5. TPM matrix 추출 QC ============
      H1("5. TPM matrix 추출 결과 (CIBERSORTx 입력)"),
      P("OpenPedCan v15의 gene-expression-rsem-tpm-collapsed.rds(285 MB)를 다운로드하여, 위 cohort에 해당하는 sample만 추출하였습니다. QC 결과는 다음과 같습니다."),
      blank(),
      makeTable(tpmRows, [3500, 5860]),
      blank(),
      H2("5.1. Deconvolution 입력 sample 구성"),
      blank(),
      makeTable(deconvRows, [3500, 3500, 2360]),
      blank(),
      bullet("Main cohort 349건 + BRAF_ALT_main 31건 + BRAF_ALT_LGG 322건 = 총 702건을 한 번에 CIBERSORTx 업로드"),
      bullet("Sample matching rate 100% (702/702) — 누락된 sample 없음"),
      bullet("All-zero gene 4,917개 제거 후 55,408개 gene retain"),
      bullet("HGNC gene symbol 기반 (중복 없음), TPM 값 음수 0, NA 0 — CIBERSORTx 입력 요건 충족"),

      // ============ 6. 산출 파일 ============
      H1("6. 산출 파일 위치"),
      P("모든 파일은 다음 경로에 저장되어 있습니다."),
      codeBlock("C:\\Users\\scott\\Downloads\\Open PBTA\\output\\"),
      H2("6.1. Cohort annotation"),
      bullet("cohort_main_final.tsv — Main cohort 349건의 환자/시료 메타정보"),
      bullet("cohort_BRAF_ALT_main.tsv — Secondary HGG/PXA/infant-fusion 31건"),
      bullet("cohort_BRAF_ALT_LGG.tsv — Secondary LGG 322건"),
      bullet("cohort_for_deconvolution.tsv — Deconvolution 입력 sample 전체 목록 (702건)"),
      bullet("tpm_manifest.tsv — TPM matrix에 실제로 포함된 sample 702건의 cohort_group 매핑"),
      H2("6.2. Table 1"),
      bullet("Table1_cohort_characteristics.tsv — Main cohort 4그룹 특성 표 (TSV)"),
      bullet("Table1_supp_BRAF_secondary.tsv — Secondary cohort 표"),
      bullet("Table1.docx — Word 형식 Table 1 (Manuscript용)"),
      H2("6.3. Expression matrix (CIBERSORTx 입력)"),
      bullet("tpm_for_cibersortx.tsv — 55,408 genes × 702 samples TPM matrix (180 MB)"),

      // ============ 7. 다음 단계 ============
      H1("7. 다음 단계 — CIBERSORTx 업로드"),
      P("준비 완료된 tpm_for_cibersortx.tsv를 CIBERSORTx 웹사이트에 업로드하면 됩니다. 절차는 다음과 같습니다."),
      H2("7.1. 사이트 접속 및 회원가입"),
      bullet("URL: https://cibersortx.stanford.edu"),
      bullet("학교 이메일로 회원가입 (Stanford에서 academic license 자동 부여)"),
      H2("7.2. 작업 실행"),
      bullet("메뉴: Menu → \"Impute Cell Fractions\""),
      bullet("Job parameters에서 다음과 같이 설정:"),
      bullet("Mixture file: tpm_for_cibersortx.tsv 업로드", 1),
      bullet("Signature matrix: LM22 (built-in 선택)", 1),
      bullet("Permutations: 100", 1),
      bullet("Quantile normalization: disabled (RNA-seq 권장 설정)", 1),
      bullet("Batch correction: disabled (OpenPedCan 단일 파이프라인이므로 불필요)", 1),
      bullet("Run mode: relative mode (LM22 기본)", 1),
      H2("7.3. 결과 처리"),
      bullet("작업 시간: 분석 큐 상황에 따라 30분 ~ 수시간"),
      bullet("완료 후 결과 CSV 다운로드 → 22 immune cell type × 702 samples 매트릭스 + p-value/correlation 메타데이터"),
      bullet("결과 CSV가 준비되면 알려주시면 됩니다 — 그 다음 immunedeconv cross-validation, consensus clustering으로 바로 이어가겠습니다."),

      // ============ 8. 사전 확인 부탁사항 ============
      H1("8. 진행 전 확인 부탁사항"),
      P("CIBERSORTx 실행 전 다음 2가지만 확인해 주시면 안전하게 진행할 수 있습니다."),
      bullet("CIBERSORTx 계정을 선생님 명의로 만들어서 실행할지, 제 계정으로 실행할지 (논문 reproducibility 측면에서는 선생님 계정 권장)"),
      bullet("180 MB 업로드는 free tier에서 가능합니다만, 혹시 차단되면 sample 단위로 batch 나누어 업로드해도 결과는 동일합니다 (이 경우 제가 batch script 추가 제공)"),

      // ============ 9. 맺음말 ============
      H1("9. 맺음말"),
      P("Step 5까지 모든 준비가 마무리되었습니다. CIBERSORTx 결과만 들어오면 ecotype clustering부터 association testing까지 한 번에 진행할 수 있는 상태입니다. 진행 중 추가 결정이 필요한 부분 있으시면 편하게 말씀 부탁드립니다."),
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
  fs.writeFileSync("C:/Users/scott/Downloads/Open PBTA/output/Step2-5_Report.docx", buf);
  console.log("Step2-5_Report.docx written");
});
