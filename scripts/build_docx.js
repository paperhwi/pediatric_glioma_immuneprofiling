const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  ImageRun, AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType,
  LevelFormat, PageOrientation,
} = require("docx");

// --------- Helpers ---------
const FONT = "Malgun Gothic";

const P = (text, opts = {}) => new Paragraph({
  spacing: { after: 120, line: 320 },
  alignment: opts.align || AlignmentType.JUSTIFIED,
  children: [new TextRun({ text, font: FONT, size: opts.size || 22, bold: opts.bold || false, color: opts.color || "1A1A1A" })],
  ...opts.extra,
});

const PMulti = (runs, opts = {}) => new Paragraph({
  spacing: { after: 120, line: 320 },
  alignment: opts.align || AlignmentType.JUSTIFIED,
  children: runs,
});

const H1 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_1,
  spacing: { before: 320, after: 160 },
  children: [new TextRun({ text, font: FONT, size: 30, bold: true, color: "2E5C8A" })],
});

const H2 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_2,
  spacing: { before: 240, after: 120 },
  children: [new TextRun({ text, font: FONT, size: 26, bold: true, color: "1A1A1A" })],
});

const bullet = (text, level = 0) => new Paragraph({
  numbering: { reference: "bullets", level: level },
  spacing: { after: 80, line: 300 },
  children: [new TextRun({ text, font: FONT, size: 22 })],
});

const blank = () => new Paragraph({ children: [new TextRun({ text: "", font: FONT, size: 22 })] });

const cellBorder = { style: BorderStyle.SINGLE, size: 4, color: "888888" };
const cellBorders = { top: cellBorder, bottom: cellBorder, left: cellBorder, right: cellBorder };

const cell = (text, opts = {}) => new TableCell({
  borders: cellBorders,
  width: { size: opts.width, type: WidthType.DXA },
  shading: opts.fill ? { fill: opts.fill, type: ShadingType.CLEAR } : undefined,
  margins: { top: 100, bottom: 100, left: 140, right: 140 },
  children: [new Paragraph({
    alignment: opts.align || AlignmentType.LEFT,
    spacing: { after: 0 },
    children: [new TextRun({ text, font: FONT, size: opts.size || 20, bold: opts.bold || false, color: opts.color || "1A1A1A" })],
  })],
});

// --------- Cohort breakdown table ---------
const cohortWidths = [1600, 1300, 1500, 1500, 1500, 1960];
const cohortRow = (cells, isHeader = false) => new TableRow({
  tableHeader: isHeader,
  children: cells.map((t, i) => cell(t, {
    width: cohortWidths[i],
    fill: isHeader ? "2E5C8A" : undefined,
    color: isHeader ? "FFFFFF" : "1A1A1A",
    bold: isHeader,
    align: i === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
  })),
});

const cohortTable = new Table({
  width: { size: 9360, type: WidthType.DXA },
  columnWidths: cohortWidths,
  rows: [
    cohortRow(["분류 (Cohort group)", "n (biospecimen)", "Midline %", "Hemispheric %", "Median age (yr)", "Survival 가능 (%)"], true),
    cohortRow(["DMG_K27 (H3 K27/K28-altered DMG)", "225", "91% (205/225)", "2% (5/225)", "8.6", "84% (189/225)"]),
    cohortRow(["DHG_G34 (H3 G34/G35-mutant DHG)", "37", "5% (2/37)", "81% (30/37)", "16.0", "57% (21/37)"]),
    cohortRow(["pHGG_WT (H3/IDH-wildtype HGG)", "192", "36% (70/192)", "52% (100/192)", "10.3", "77% (147/192)"]),
    cohortRow(["IHG (Infant-type hemispheric)", "28", "18% (5/28)", "82% (23/28)", "0.4", "86% (24/28)"]),
    cohortRow(["Main 합계", "482", "—", "—", "—", "79% (381/482)"]),
    cohortRow(["BRAF_ALT (secondary)", "384", "—", "—", "—", "—"]),
  ],
});

// --------- Naming difference table ---------
const namingWidths = [3120, 3120, 3120];
const namingRow = (cells, isHeader = false) => new TableRow({
  tableHeader: isHeader,
  children: cells.map((t, i) => cell(t, {
    width: namingWidths[i],
    fill: isHeader ? "2E5C8A" : (i === 0 ? "F0F4F8" : undefined),
    color: isHeader ? "FFFFFF" : "1A1A1A",
    bold: isHeader || i === 0,
    align: AlignmentType.LEFT,
  })),
});

const namingTable = new Table({
  width: { size: 9360, type: WidthType.DXA },
  columnWidths: namingWidths,
  rows: [
    namingRow(["분류 (초안 표기)", "OpenPedCan v15 표기", "비고"], true),
    namingRow(["H3 K27-altered (DMG)", "H3 K28-altered", "동일 잔기, IUPAC 신 명명법"]),
    namingRow(["H3 G34-mutant (DHG)", "H3 G35-mutant", "동일 잔기, IUPAC 신 명명법"]),
    namingRow(["Infant-type hemispheric glioma", "직접 라벨 없음", "age < 18mo + hemispheric + RTK/BRAF fusion으로 운영 정의"]),
    namingRow(["Pediatric-type HGG, H3/IDH-WT", "High-grade glioma, IDH-wildtype and H3-wildtype", "harmonized_diagnosis로 매칭"]),
    namingRow(["BRAF-altered / RTK fusion", "molecular_subtype + fusion 파일에서 통합", "BRAF V600E, KIAA1549-BRAF, NTRK/ALK/ROS1/MET 등"]),
  ],
});

// --------- Plan timeline table ---------
const planWidths = [900, 2800, 3000, 2660];
const planRow = (cells, isHeader = false) => new TableRow({
  tableHeader: isHeader,
  children: cells.map((t, i) => cell(t, {
    width: planWidths[i],
    fill: isHeader ? "2E5C8A" : undefined,
    color: isHeader ? "FFFFFF" : "1A1A1A",
    bold: isHeader,
    align: i === 0 ? AlignmentType.CENTER : AlignmentType.LEFT,
  })),
});

const planTable = new Table({
  width: { size: 9360, type: WidthType.DXA },
  columnWidths: planWidths,
  rows: [
    planRow(["Step", "분석 단계", "산출물", "검증/확인 지표"], true),
    planRow(["1", "Cohort assembly (완료)", "cohort_main.tsv (n=482) + secondary (n=384)", "diagnosis × location × age 분포 확인"]),
    planRow(["2", "Expression matrix 추출 + HGNC 매핑", "FPKM/TPM matrix (gene × sample)", "결측 0, 유전자/샘플 수 보고"]),
    planRow(["3", "CIBERSORTx (LM22) deconvolution", "22 immune cell × sample 매트릭스 + p-value", "p-value 분포, immune cold 보존 여부 확인"]),
    planRow(["4", "교차 검증 (xCell, quanTIseq)", "주요 cell type 상관계수", "Spearman ρ > 0.5"]),
    planRow(["5", "Consensus clustering (k=2~6)", "PAC plot, silhouette, optimal k", "k=3 가설 검정"]),
    planRow(["6", "Ecotype annotation", "Lymphocyte-inflamed / Myeloid / Desert 라벨", "각 cluster top immune feature heatmap"]),
    planRow(["7", "Association testing", "H3 × ecotype, location × ecotype, age × ecotype", "Chi-square / Fisher, FDR < 0.05"]),
    planRow(["8", "ssGSEA pathway scores", "Antigen presentation, IFN, MAPK, microglia 등", "Boxplot + Kruskal–Wallis"]),
    planRow(["9", "Survival (exploratory)", "KM curve by ecotype (n with OS)", "Log-rank p-value, exploratory 명시"]),
    planRow(["10", "BRAF/RTK supplementary", "Secondary cohort 분포 + MAPK score", "Main 대비 비교"]),
    planRow(["11", "Figure/Table drafting", "Figure 1–6 + Table 1–3 초안", "초안과의 매칭"]),
  ],
});

// --------- Document ---------
const doc = new Document({
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 30, bold: true, font: FONT, color: "2E5C8A" },
        paragraph: { spacing: { before: 320, after: 160 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: FONT },
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
      { reference: "numbers",
        levels: [
          { level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 480, hanging: 280 } } } },
        ] },
    ],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
      },
    },
    children: [
      // Title
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 120 },
        children: [new TextRun({
          text: "Pediatric Glioma Immune Ecotype 분석 진행 보고",
          font: FONT, size: 32, bold: true, color: "2E5C8A",
        })],
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 320 },
        children: [new TextRun({
          text: "OpenPedCan v15 기반 코호트 필터링 결과 및 분석 계획",
          font: FONT, size: 24, color: "555555",
        })],
      }),

      // --- 1. 인사 ---
      H1("1. 인사 말씀"),
      P("선생님, 안녕하십니까."),
      P("지난번 주신 OpenPBTA 관련 분석 요청 잘 받았습니다. 보내주신 초안(Decoding the Tumor Immune Microenvironment of Pediatric Gliomas)과 가이드라인 문서를 함께 검토하였고, 우선 가능한 범위 안에서 데이터셋을 확보하여 코호트 필터링까지 진행하였습니다. 이번 메일에서는 (1) 제가 이해한 초안 요약, (2) 초안에 대한 의견, (3) 분석 가능 여부, (4) 사용 데이터셋(OpenPedCan v15) 및 초안 분류법과의 차이, (5) 현재까지 만든 코호트 결과, (6) 분석 시 유의 사항, (7) 앞으로의 단계별 진행 계획을 정리하여 보고드리고자 합니다."),
      P("내용이 다소 길어졌습니다만, 본격 분석을 시작하기 전에 선생님 의견을 한 번 정리해서 받는 편이 이후 진행에 도움이 될 것 같아 다소 자세히 작성하였습니다. 양해 부탁드립니다."),

      // --- 2. 초안 요약 ---
      H1("2. 제가 이해한 초안 요약"),
      P("보내주신 초안은 OpenPBTA의 bulk RNA-seq 데이터를 활용하여 소아 glioma의 종양면역 미세환경(TME)을 면역세포 deconvolution으로 재구성하고, 이를 \"immune ecotype\"으로 분류하는 dry-lab 연구입니다. 핵심 흐름은 다음과 같이 이해하였습니다."),
      bullet("CIBERSORTx + LM22 signature matrix로 22종 면역/기질 세포 비율을 추정."),
      bullet("Deconvolution p-value < 0.05인 sample만 retain (초안 기준 1,031 → 283, 27.4%)."),
      bullet("Unsupervised clustering을 통해 세 가지 ecotype을 정의: Lymphocyte-inflamed, Myeloid-dominant, Immune-desert."),
      bullet("Ecotype을 (i) histone mutation class (K27, G34, H3-WT), (ii) 종양 위치 (midline vs hemispheric), (iii) developmental context (특히 infant-type)와 연결."),
      bullet("Exploratory survival analysis 및 pathway enrichment (antigen presentation, IFN, MAPK 등)."),
      P("초안의 핵심 메시지는, 소아 glioma가 단순한 \"immune-cold tumor\"가 아니라, histone alteration · anatomical location · developmental context와 결합된 면역 ecotype의 스펙트럼이라는 점으로 이해하였습니다. 가이드라인에서 말씀하신 대로, \"최초의 immune classification\"이 아니라 기존 immune profiling 위에 developmental immune ecotype framework로 reframing하는 방향이 적절하다고 판단됩니다."),

      // --- 3. 초안에 대한 의견 ---
      H1("3. 초안에 대한 의견"),
      P("전체적으로 Cancers(MDPI)에 투고하기에 적합한 규모와 구성으로 보입니다. 특히 Fecci 연구진의 T-cell 배제(immune exclusion) 모델과 developmental immune privilege 개념을 연결한 Discussion 부분이 임상적 통찰을 잘 보여주고 있다고 생각합니다. 다만 실제 분석을 시작하기에 앞서 몇 가지 짚어드려야 할 부분이 있어 정리하여 말씀드리겠습니다."),
      H2("3.1. 강점"),
      bullet("OpenPBTA / OpenPedCan은 공개 데이터셋 중 소아 brain tumor RNA-seq에서 가장 잘 정돈된 자원으로, 재현 가능한 dry-lab 논문의 기반으로 적합합니다."),
      bullet("CIBERSORTx + LM22는 reviewer에게 익숙한 방법으로, methods 측면에서 큰 거부감을 받기 어려운 표준 파이프라인입니다."),
      bullet("\"Developmental immune ecotype framework\"라는 reframing은 기존 pediatric glioma immune classification 논문들과의 차별화에 효과적인 메시지로 보입니다."),
      H2("3.2. 보완이 필요한 부분"),
      bullet("초안의 \"1,031 → 283 (27.4%) p<0.05 통과\" 수치는 ChatGPT 추정치로 보이며, 실제 OpenPedCan v15에서는 PBTA RNA-Seq tumor가 2,327건이고, 그중 가이드라인에 부합하는 4개 진단군의 합계는 482건이었습니다 (자세한 표는 5번 항목 참조)."),
      bullet("LM22 signature는 말초혈액 면역세포 기반이라 microglia가 포함되어 있지 않습니다. 뇌종양 TME에서는 macrophage signal의 상당 부분이 microglia일 가능성이 높아, limitation에 명시하거나 brain-tuned signature (Klemm et al. 2020 등)와 cross-validation이 필요합니다."),
      bullet("\"p<0.05 필터로 immune cold sample이 제거되는\" 문제는 immune-desert ecotype을 정의하려는 연구 목적과 모순될 수 있습니다. 이 부분은 전수 retain + sensitivity analysis로 가는 편이 안전합니다."),
      bullet("초안의 일부 수치, 도표, 일부 method 기술이 실제 분석 산출물이 아닌 가설적 서술이라 명시해 주신 점은 잘 인지하였습니다. 이번 분석에서는 실제 OpenPedCan v15 산출물로 모두 대체하여 정리할 예정입니다."),

      // --- 4. 분석 가능 여부 ---
      H1("4. 분석 가능 여부"),
      P("결론부터 말씀드리면, 가이드라인에서 말씀하신 모든 분석 항목이 OpenPedCan v15 데이터셋만으로 충분히 수행 가능합니다. 각 항목별 가능성과 주의점은 아래와 같습니다."),
      bullet("Sample barcode / annotation 정리: 가능 (이미 완료, cohort_main.tsv)"),
      bullet("Bulk RNA-seq expression matrix 매칭: 가능 (TPM matrix 사용 예정)"),
      bullet("CIBERSORTx LM22 deconvolution: 가능 (단, microglia 부재는 limitation으로 명시)"),
      bullet("3-ecotype clustering: 가능 (consensus clustering으로 k 데이터 기반 결정 권장)"),
      bullet("Histone class · location · age · diagnosis 연관성: 가능"),
      bullet("Survival (exploratory): 가능 (Main cohort 481/482 중 381건에 OS 데이터 존재)"),
      bullet("Pathway enrichment (antigen presentation, IFN, MAPK, microglia signature 등): 가능"),
      bullet("BRAF/RTK secondary cohort 분석: 가능 (n=384 확보)"),

      // --- 5. 사용 데이터셋 및 초안 분류법과의 차이 ---
      H1("5. 사용 데이터셋과 초안 분류법과의 차이"),
      H2("5.1. 데이터셋 — OpenPedCan v15"),
      P("초안에서 언급하신 OpenPBTA는 현재 OpenPedCan 프로젝트로 확장·통합되어 v15 release가 최신입니다 (Children's Brain Tumor Network, AlexsLemonade/d3b-center 공동 관리). 본 분석에서는 다음 세 파일을 우선 다운로드하여 사용하였습니다."),
      bullet("histologies.tsv — 47,895 specimen × 66 컬럼의 통합 임상·분자 annotation"),
      bullet("fusion-putative-oncogenic.tsv — BRAF/RTK fusion 식별용"),
      bullet("independent-specimens.rnaseqpanel.primary-plus.tsv — 환자당 중복 sample 제거용 권장 리스트"),
      P("(추후 expression matrix gene-expression-rsem-tpm-collapsed.rds를 추가 다운로드하여 CIBERSORTx 입력으로 사용 예정.)"),
      H2("5.2. 초안 분류법과 OpenPedCan v15의 차이"),
      P("초안에 사용된 분류명은 WHO 2021 이전의 표기(H3 K27, H3 G34)인 반면, OpenPedCan v15는 IUPAC 신 명명법(H3 K28, H3 G35)을 사용합니다. 동일한 잔기를 지칭하는 것이며, 논문에서는 K27/G34 표기를 유지하되 Methods에 \"K28-altered (formerly K27M-altered)\" 형태로 명시하면 됩니다."),
      blank(),
      namingTable,
      blank(),
      P("Infant-type hemispheric glioma (I-HGG)는 WHO 2021에 새로 도입된 entity로, OpenPedCan에 직접 라벨이 부여되어 있지 않습니다. 이에 따라 다음과 같이 운영 정의를 적용하였습니다."),
      bullet("진단 시 연령 < 18개월"),
      bullet("Hemispheric 위치 (Frontal/Parietal/Temporal/Occipital 또는 CNS_region = Hemispheric)"),
      bullet("NTRK1/2/3, ALK, ROS1, MET, BRAF 중 하나의 fusion 존재 (fusion-putative-oncogenic.tsv)"),
      P("위 기준으로 28건이 분류되었으며, 중위 연령 0.4세, 23/28건이 hemispheric이라는 점에서 entity 정의에 부합합니다."),

      // --- 6. 코호트 필터링 결과 ---
      H1("6. 코호트 필터링 결과"),
      P("위 기준에 따라 PBTA + RNA-Seq + Tumor sample을 필터링하고, 4개 main 진단군과 1개 secondary 진단군으로 분류한 결과는 아래 표와 같습니다."),
      blank(),
      cohortTable,
      blank(),
      P("생물학적으로 기대되는 분포 — DMG_K27의 91%가 midline, DHG_G34의 81%가 hemispheric이며 청소년 연령대, IHG가 100% infant — 가 잘 재현되어 코호트 필터링이 적절히 작동함을 확인하였습니다."),
      P("결과 파일은 모두 output 디렉터리에 저장하였으며, 환자 ID · biospecimen ID · sample/aliquot ID · pathology/harmonized/molecular subtype · primary_site · CNS_region · location_class · age (days/years/months/developmental group) · sex · OS_days · OS_status · BRAF/IHG fusion 플래그 · independent specimen 플래그를 모두 포함합니다. 필요하시면 Table 1 형태로 정리하여 별도 송부드리겠습니다."),

      // --- 7. 분석 시 유의 사항 ---
      H1("7. 분석 시 유의 사항"),
      P("실제 분석을 진행함에 있어 다음 다섯 가지를 먼저 합의해 두는 편이 안전합니다."),
      H2("7.1. CIBERSORTx p-value 필터링"),
      P("초안에서는 p < 0.05를 cutoff로 사용하나, 이 경우 immune cold tumor가 사전에 제거되어 immune-desert ecotype 정의 자체가 어려워질 수 있습니다. → 전수 retain 후 p-value를 메타데이터로 보유하고, sensitivity analysis로 \"high-confidence subset\"을 supplementary에 제시하는 방향을 권장드립니다."),
      H2("7.2. LM22의 한계와 brain-tuned reference"),
      P("LM22는 microglia를 포함하지 않습니다. 뇌종양 TME에서 macrophage signal의 상당 부분이 microglia일 가능성이 있으므로, 다음 중 한 가지 방향을 권장드립니다."),
      bullet("(보수적) LM22만 사용하고 limitation에 명시 — 가장 빠른 경로"),
      bullet("(권장) immunedeconv 패키지로 xCell · quanTIseq · MCP-counter를 cross-validation에 추가"),
      bullet("(공격적) Klemm 2020 또는 BrainImmuneAtlas scRNA-seq 기반 custom signature 생성"),
      H2("7.3. 표본 수와 clustering 안정성"),
      P("Main cohort 482건은 통계적으로 충분하나, 개별 subtype (DHG_G34 37, IHG 28)은 cluster 안정성에 영향을 줄 수 있습니다. → consensus clustering의 cumulative distribution function(CDF) 및 proportion of ambiguous clustering(PAC) plot으로 k를 결정하고, bootstrap resampling으로 ecotype assignment의 robustness를 제시할 예정입니다."),
      H2("7.4. Independent specimen 적용 여부"),
      P("환자당 다중 sample이 있는 경우(주로 longitudinal/recurrence) 통계 검정의 독립성이 깨질 수 있습니다. → independent-primary-plus 권장 리스트 적용 시 main cohort가 약 359건으로 줄어듭니다. 주분석은 \"independent\" 기준, supplementary는 전수 기준으로 두 번 제시하는 방향을 제안드립니다."),
      H2("7.5. Survival 해석"),
      P("Main cohort 482건 중 381건(79%)에 OS 데이터가 있어 KM curve 자체는 가능합니다만, ecotype × histone subgroup으로 stratify하면 검정력이 낮아질 수 있습니다. 가이드라인에서 말씀하신 대로 exploratory 수준으로 두고, 본문 figure 6 (KM) + supplementary Cox 정도로 정리하는 편이 안전합니다."),

      // --- 8. 분석 워크플로우 ---
      H1("8. 분석 워크플로우"),
      P("앞서 정리한 단계를 도식으로 표현하면 다음과 같습니다."),
      blank(),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        children: [new ImageRun({
          type: "png",
          data: fs.readFileSync(path.resolve("C:/Users/scott/Downloads/Open PBTA/output/workflow.png")),
          transformation: { width: 580, height: 395 },
          altText: { title: "Analysis workflow", description: "Workflow diagram", name: "workflow" },
        })],
      }),
      blank(),

      // --- 9. 단계별 진행 계획 ---
      H1("9. 단계별 진행 계획"),
      P("위 워크플로우를 단계별로 분해하여, 각 단계의 산출물과 검증 지표를 정리하면 아래와 같습니다."),
      blank(),
      planTable,
      blank(),
      P("이 중 Step 1(cohort assembly)은 이번 메일과 함께 보고드린 부분으로 완료된 상태입니다. 선생님 confirmation을 받는 대로 Step 2(expression matrix 추출) → Step 3(CIBERSORTx 실행)으로 넘어가려 합니다. CIBERSORTx는 Stanford web tool에서 무료로 실행 가능하나 회원 가입이 필요하므로, 제 계정으로 실행하겠습니다."),

      // --- 10. 사전 confirmation 요청 ---
      H1("10. 사전 확인 부탁드리는 사항"),
      P("본격 분석으로 넘어가기 전에 다음 3가지를 한 번 결정해 주시면 이후 진행이 매끄러울 것 같습니다."),
      bullet("Deconvolution 도구를 CIBERSORTx(LM22)만으로 갈지, immunedeconv로 xCell·quanTIseq를 cross-validation에 함께 둘지"),
      bullet("Independent-primary-plus 필터를 main analysis에 적용할지(보수적), 또는 전수 기준을 main으로 두고 independent를 supplementary로 둘지"),
      bullet("BRAF_ALT(n=384) supplementary cohort의 범위 — LGG 포함 vs HGG/PXA만 포함"),

      // --- 11. 맺음말 ---
      H1("11. 맺음말"),
      P("이상이 현재까지의 진행 사항입니다. 선생님께서 임상적·생물학적 관점에서 추가로 보완하거나 강조할 부분이 있다면 말씀해 주시면 반영하여 진행하겠습니다. 특히 infant-type hemispheric glioma의 운영 정의(연령 cutoff 18개월), pHGG_WT에서 21세 이상 outlier 처리, midline/hemispheric 경계가 모호한 경우(예: pineal region, suprasellar)의 처리 방향에 대해서도 의견 주시면 감사드리겠습니다."),
      P("긴 글 읽어주셔서 감사드립니다. 회신 받는 대로 다음 단계 진행 후 결과 다시 정리하여 보고드리겠습니다."),
      blank(),
      P("감사합니다."),
      blank(),
      blank(),
      new Paragraph({
        spacing: { after: 60 },
        children: [new TextRun({ text: "— 이상 —", font: FONT, size: 22, color: "555555" })],
        alignment: AlignmentType.CENTER,
      }),
    ],
  }],
});

const outPath = "C:/Users/scott/Downloads/Open PBTA/output/Pediatric_Glioma_Immune_Ecotype_Report.docx";
Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(outPath, buf);
  console.log("DOCX written:", outPath);
});
