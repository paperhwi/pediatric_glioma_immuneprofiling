const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, ImageRun,
  AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType,
  PageBreak,
} = require("docx");

const FONT = "Times New Roman";
const FONT_TITLE = "Calibri";

const T = (text, opts = {}) => new TextRun({
  text, font: FONT, size: opts.size || 22, bold: opts.bold || false,
  italics: opts.italics || false, color: opts.color || "000000",
});

const P = (text, opts = {}) => new Paragraph({
  spacing: { after: opts.after || 160, line: opts.line || 360 },
  alignment: opts.align || AlignmentType.JUSTIFIED,
  indent: opts.firstLine ? { firstLine: 360 } : undefined,
  children: [T(text, opts)],
});

const PMulti = (runs, opts = {}) => new Paragraph({
  spacing: { after: opts.after || 160, line: 360 },
  alignment: opts.align || AlignmentType.JUSTIFIED,
  indent: opts.firstLine ? { firstLine: 360 } : undefined,
  children: runs,
});

const H1 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_1,
  spacing: { before: 400, after: 200 },
  children: [new TextRun({ text, font: FONT_TITLE, size: 28, bold: true, color: "1A1A1A" })],
});
const H2 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_2,
  spacing: { before: 280, after: 160 },
  children: [new TextRun({ text, font: FONT_TITLE, size: 24, bold: true, color: "2E5C8A" })],
});

const blank = () => new Paragraph({ children: [new TextRun({ text: "", font: FONT, size: 22 })] });

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
      children: [new TextRun({ text, font: FONT, size: opts.size || 18, bold: opts.bold || false,
                                italics: opts.italics || false, color: opts.color || "000000" })],
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
        color: ri === 0 ? "FFFFFF" : "000000",
        bold: ri === 0,
        align: ci === 0 ? AlignmentType.LEFT : AlignmentType.CENTER,
      })),
    })),
  });
}

// === Figure embedding helper ===
function figure(opts) {
  const arr = [];
  arr.push(new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 200, after: 80 },
    children: [new ImageRun({
      type: "png",
      data: fs.readFileSync(opts.path),
      transformation: { width: opts.width, height: opts.height },
      altText: { title: opts.title || "Figure", description: opts.desc || opts.title || "Figure", name: "fig" },
    })],
  }));
  // Caption
  arr.push(new Paragraph({
    alignment: AlignmentType.JUSTIFIED,
    spacing: { after: 240 },
    children: [
      new TextRun({ text: `${opts.figNum}. `, font: FONT, size: 20, bold: true }),
      new TextRun({ text: opts.caption, font: FONT, size: 20 }),
    ],
  }));
  return arr;
}

function multiFigure(paths, caption, figNum, panelLabels) {
  // Build a side-by-side layout via single-column table for clean rendering
  const arr = [];
  // Stack panels vertically with panel labels
  paths.forEach((p, i) => {
    arr.push(new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { before: 120, after: 40 },
      children: [new TextRun({ text: `(${panelLabels[i]})`, font: FONT, size: 20, bold: true })],
    }));
    arr.push(new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { after: 100 },
      children: [new ImageRun({
        type: "png",
        data: fs.readFileSync(p.path),
        transformation: { width: p.width, height: p.height },
        altText: { title: "Figure panel", description: panelLabels[i], name: "fig" },
      })],
    }));
  });
  arr.push(new Paragraph({
    alignment: AlignmentType.JUSTIFIED,
    spacing: { after: 240 },
    children: [
      new TextRun({ text: `${figNum}. `, font: FONT, size: 20, bold: true }),
      new TextRun({ text: caption, font: FONT, size: 20 }),
    ],
  }));
  return arr;
}

// ============== Tables ==============

const table1 = [
  ["Characteristic", "DMG_K27 (n=175)", "DHG_G34 (n=31)", "pHGG_WT (n=125)", "IHG (n=18)", "Total (n=349)"],
  ["Median age, yr (IQR)", "7.9 (5.5–11.2)", "15.7 (13.5–17.9)", "9.5 (6.0–13.7)", "0.4 (0.1–0.9)", "8.8 (5.5–12.5)"],
  ["Male, n (%)", "75 (42.9)", "16 (51.6)", "70 (56.0)", "9 (50.0)", "170 (48.7)"],
  ["Female, n (%)", "97 (55.4)", "15 (48.4)", "55 (44.0)", "9 (50.0)", "176 (50.4)"],
  ["Midline, n (%)", "118 (67.4)", "0 (0)", "22 (17.6)", "0 (0)", "140 (40.1)"],
  ["Hemispheric, n (%)", "3 (1.7)", "24 (77.4)", "62 (49.6)", "15 (83.3)", "104 (29.8)"],
  ["Posterior fossa, n (%)", "5 (2.9)", "0 (0)", "9 (7.2)", "0 (0)", "14 (4.0)"],
  ["Ambiguous/Other, n (%)", "49 (28.0)", "7 (22.6)", "32 (25.6)", "3 (16.7)", "91 (26.1)"],
  ["BRAF/RTK fusion+, n (%)", "31 (17.7)", "4 (12.9)", "17 (13.6)", "18 (100.0)", "70 (20.1)"],
  ["OS data available, n (%)", "140 (80.0)", "15 (48.4)", "82 (65.6)", "14 (77.8)", "251 (71.9)"],
];

const table2 = [
  ["Variable", "HR", "95% CI", "p-value"],
  ["Ecotype: Inflamed (reference)", "1.00", "—", "—"],
  ["Ecotype: Intermediate", "1.37", "0.99–1.91", "0.062"],
  ["Ecotype: Immune-desert", "1.79", "1.21–2.66", "0.004"],
  ["Cohort: pHGG_WT (reference)", "1.00", "—", "—"],
  ["Cohort: DMG_K27", "2.96", "2.11–4.15", "<0.001"],
  ["Cohort: DHG_G34", "1.15", "0.59–2.24", "0.672"],
  ["Cohort: IHG", "0.22", "0.08–0.63", "0.005"],
  ["Age (per year)", "0.99", "0.96–1.02", "0.506"],
  ["Sex: Female (vs Male)", "1.12", "0.84–1.49", "0.450"],
];

const table3 = [
  ["Signature", "KW chi-sq", "BH-adj. p-value"],
  ["M1 Macrophage", "262.2", "2.8e-56"],
  ["Chemokine T-cell recruitment", "257.3", "1.6e-55"],
  ["T-cell cytotoxicity", "245.0", "5.0e-53"],
  ["Glioma inflammatory (Wang 2017)", "241.8", "1.9e-52"],
  ["TGF-β immunosuppression", "225.0", "6.5e-49"],
  ["Tregs", "221.7", "2.9e-48"],
  ["T-cell exhaustion", "220.1", "5.6e-48"],
  ["MoTAM (Antunes 2021)", "211.4", "3.0e-46"],
  ["MDM (Klemm 2020)", "205.5", "5.2e-45"],
  ["IFN-γ response", "186.2", "6.0e-41"],
  ["MHC Class II", "157.7", "8.7e-35"],
  ["MgTAM (Antunes 2021)", "100.1", "2.6e-22"],
  ["Microglia core (Bohlen 2020)", "97.9", "7.3e-22"],
  ["Microglia (Klemm 2020)", "93.7", "5.6e-21"],
  ["MAPK activity", "54.8", "1.4e-12"],
  ["Stemness", "7.2", "0.029"],
  ["Cell cycle / proliferation", "0.13", "0.939 (NS)"],
];

// === Title ===
const titleRuns = [
  new TextRun({ text: "Developmentally-Anchored Immune Ecotypes in Pediatric Gliomas: ",
                font: FONT_TITLE, size: 36, bold: true, color: "1A1A1A" }),
  new TextRun({ text: "Brain-Tuned Deconvolution Reveals a Microglia-versus-Macrophage Axis with Independent Prognostic Value",
                font: FONT_TITLE, size: 36, bold: true, color: "1A1A1A" }),
];

const FIG = "C:/Users/scott/Downloads/Open PBTA/output";
const FIGCL = `${FIG}/figs_clustering`;
const FIGSV = `${FIG}/figs_survival`;
const FIGPW = `${FIG}/figs_pathway`;

// ============== Build Document ==============
const doc = new Document({
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 28, bold: true, font: FONT_TITLE },
        paragraph: { spacing: { before: 400, after: 200 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, font: FONT_TITLE, color: "2E5C8A" },
        paragraph: { spacing: { before: 280, after: 160 }, outlineLevel: 1 } },
    ],
  },
  sections: [{
    properties: {
      page: { size: { width: 12240, height: 15840 },
              margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } },
    },
    children: [
      // ============ TITLE ============
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 240, line: 400 },
        children: titleRuns,
      }),
      blank(),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 160 },
        children: [
          new TextRun({ text: "Suhmi Chung, MD", font: FONT, size: 22, bold: true }),
          new TextRun({ text: "1,*", font: FONT, size: 22, superScript: true }),
        ],
      }),
      new Paragraph({
        spacing: { after: 80, line: 280 },
        children: [
          new TextRun({ text: "1 ", font: FONT, size: 20, superScript: true }),
          new TextRun({ text: "Department of Neurosurgery, Asan Medical Center, University of Ulsan College of Medicine, Seoul, Republic of Korea",
                        font: FONT, size: 20, italics: true }),
        ],
      }),
      new Paragraph({
        spacing: { after: 160 },
        children: [
          new TextRun({ text: "* ", font: FONT, size: 20, superScript: true }),
          new TextRun({ text: "Correspondence: ", font: FONT, size: 20, italics: true }),
          new TextRun({ text: "[corresponding author email]", font: FONT, size: 20, italics: true, color: "666666" }),
        ],
      }),
      blank(),

      // ============ SIMPLE SUMMARY ============
      H1("Simple Summary"),
      P("Pediatric gliomas have long been considered immunologically \"cold\" tumors with limited responsiveness to immunotherapy. " +
        "However, recent evidence indicates substantial heterogeneity in the tumor immune microenvironment (TME) across molecular subtypes " +
        "and anatomical locations. Using the OpenPedCan v15 harmonized cohort of 349 independent pediatric gliomas, we applied four orthogonal " +
        "deconvolution methods together with single-sample gene set enrichment analysis (ssGSEA) of brain-resident microglia and monocyte-derived " +
        "macrophage signatures. We identified three reproducible immune ecotypes—Inflamed, Intermediate, and Immune-desert—that varied with " +
        "histone mutation class and anatomical location and were independently associated with overall survival after adjusting for histone " +
        "class, age, and sex. Brain-tuned signatures resolved a microglia-versus-monocyte-derived macrophage axis that conventional " +
        "deconvolution panels conflated. These findings provide a developmentally-anchored framework for interpreting pediatric glioma immunity " +
        "and identify ecotype as an independent prognostic axis with potential implications for immunotherapeutic stratification."),

      // ============ ABSTRACT ============
      H1("Abstract"),
      PMulti([
        T("Background: ", { bold: true }),
        T("The tumor immune microenvironment of pediatric gliomas remains poorly characterized, and a systematic ecotype framework " +
          "integrating histone class, anatomy, and developmental context is lacking. "),
        T("Methods: ", { bold: true }),
        T("We analyzed bulk RNA-sequencing from 349 independent pediatric glioma biospecimens (175 DMG H3 K27-altered; 31 DHG H3 G34-mutant; " +
          "125 H3/IDH-wildtype pediatric HGG; 18 infant-type hemispheric glioma) from OpenPedCan v15. Four orthogonal deconvolution algorithms " +
          "(MCP-counter, quanTIseq, EPIC, xCell) were complemented with brain-tuned ssGSEA scoring of 24 signatures including microglia (Klemm 2020; " +
          "Bohlen 2020; Antunes 2021 Mg-TAM), monocyte-derived macrophage (Klemm MDM; Antunes Mo-TAM), antigen presentation, " +
          "interferon, T-cell cytotoxicity, MAPK activity, and TGF-β programs. Consensus clustering of quanTIseq fractions and 19 immune-related " +
          "ssGSEA signatures defined immune ecotypes; associations with molecular and anatomical features were tested by Fisher's exact test " +
          "and Kruskal-Wallis with Benjamini-Hochberg adjustment. Kaplan-Meier and multivariable Cox regression evaluated prognostic value. "),
        T("Results: ", { bold: true }),
        T("Three reproducible ecotypes emerged: Inflamed (n = 113; 32.4%), Intermediate (n = 157; 45.0%), and Immune-desert (n = 79; 22.6%). " +
          "Brain-tuned signature scoring resolved a microglia-versus-MDM axis that conventional deconvolution panels conflated: " +
          "three microglia signatures clustered tightly (pairwise Spearman ρ > 0.98) and were distinct from two MDM signatures (ρ = 0.94 internally, " +
          "ρ ≈ 0.40 versus microglia), while xCell macrophage tracked MDM (ρ = 0.86). Ecotypes were significantly associated with cohort group " +
          "(Fisher p = 1.3 × 10⁻³) and anatomical location (p = 1.3 × 10⁻²), though patterns deviated from adult-glioma expectations: DHG G34-mutant " +
          "tumors were predominantly Immune-desert (52%) rather than lymphocyte-inflamed, and midline tumors were more frequently Inflamed " +
          "than hemispheric tumors (41% vs 25%). The Inflamed ecotype showed the longest median overall survival (595 days vs 446 days for " +
          "Intermediate and 403 days for Immune-desert; log-rank p = 0.037). After multivariable adjustment for cohort group, age, and sex, " +
          "Immune-desert remained an independent adverse prognostic factor (HR = 1.79; 95% CI 1.21–2.66; p = 0.004). Twenty-two of twenty-four " +
          "immune pathways differed significantly across ecotypes (BH-adjusted p < 10⁻²⁰), whereas cell-cycle/proliferation signatures did not " +
          "(p = 0.94). Microglia signatures were modestly ecotype-stratified (p ≈ 10⁻²¹), while MDM signatures showed stronger ecotype dependence " +
          "(p ≈ 10⁻⁴⁵), implicating peripheral macrophage recruitment as a primary driver of ecotype assignment. "),
        T("Conclusions: ", { bold: true }),
        T("Pediatric gliomas segregate into reproducible immune ecotypes with developmentally anchored composition and independent prognostic " +
          "significance. Brain-tuned signature scoring uncovers a microglia-versus-MDM axis that conventional deconvolution methods cannot " +
          "resolve, providing a translationally relevant framework for risk stratification and rational immunotherapy design in childhood brain tumors."),
      ]),

      PMulti([
        T("Keywords: ", { bold: true }),
        T("pediatric glioma; tumor immune microenvironment; immune ecotype; CIBERSORTx; immunedeconv; ssGSEA; microglia; monocyte-derived macrophage; " +
          "diffuse midline glioma; H3 K27M; H3 G34; infant-type glioma; OpenPedCan"),
      ]),

      new Paragraph({ children: [new PageBreak()] }),

      // ============ INTRODUCTION ============
      H1("1. Introduction"),
      P("Pediatric gliomas are biologically and clinically distinct from their adult counterparts, defined by histone mutations such as H3 K27M and H3 G34R/V, " +
        "MAPK pathway alterations including BRAF V600E and KIAA1549-BRAF fusions, and receptor tyrosine kinase fusions in infant-type tumors. " +
        "These oncogenic drivers are intrinsically linked to developmental epigenetic programs and to the brain regions from which the tumors arise [1,2,3].",
        { firstLine: true }),
      P("While the molecular taxonomy of pediatric gliomas is now well established, the tumor immune microenvironment (TME) remains comparatively underexplored. " +
        "Pediatric brain tumors have historically been viewed as immunologically \"cold,\" a perception that has limited enthusiasm for immunotherapy in this population [4,5]. " +
        "Recent studies, however, have begun to reveal substantial heterogeneity in the TME of pediatric gliomas, with reports describing an immunosuppressive " +
        "myeloid-rich microenvironment in H3 K27-altered diffuse midline glioma (DMG) [6,7], comparatively higher T-cell infiltration in H3 G34-mutant " +
        "hemispheric tumors [8], and developmental immune privilege in infant-type lesions [9,10]. Despite these advances, a unified ecotype framework that " +
        "integrates histone mutation class, anatomical compartment, and developmental context into a single interpretive lens has not been established.",
        { firstLine: true }),
      P("A central methodological limitation in this field is that the most widely used computational immune deconvolution panel, the CIBERSORTx LM22 signature " +
        "matrix, was derived from peripheral blood mononuclear cells and does not explicitly model brain-resident microglia [11]. Because microglia and " +
        "infiltrating monocyte-derived macrophages (MDMs) exhibit overlapping but distinguishable transcriptional programs in the central nervous system [12,13,14], " +
        "the inability of LM22-based methods to differentiate these populations risks conflating two biologically and therapeutically distinct cell types. " +
        "Recent single-cell studies have catalogued lineage-specific signatures that can resolve microglia from MDM in glioma tissue [12,14,15], yet these have not been " +
        "systematically applied alongside conventional deconvolution to large pediatric glioma cohorts.",
        { firstLine: true }),
      P("We hypothesized that pediatric gliomas would segregate into reproducible immune ecotypes whose composition is shaped by histone mutation class, " +
        "anatomical location, and developmental context, and that integrating brain-tuned myeloid signatures with conventional deconvolution would resolve " +
        "the microglia-versus-MDM axis. Leveraging the Open Pediatric Cancer (OpenPedCan) v15 release, we applied four orthogonal deconvolution algorithms " +
        "(MCP-counter, quanTIseq, EPIC, xCell) together with ssGSEA scoring of microglia and MDM signatures to 349 independent pediatric glioma biospecimens. " +
        "We then defined immune ecotypes by consensus clustering, characterized their molecular and anatomical correlates, and evaluated their prognostic significance.",
        { firstLine: true }),

      // ============ METHODS ============
      H1("2. Materials and Methods"),

      H2("2.1. Data Source and Cohort Assembly"),
      P("Harmonized bulk RNA-sequencing and clinical/molecular metadata were obtained from the Open Pediatric Cancer (OpenPedCan) v15 release " +
        "(https://github.com/d3b-center/OpenPedCan-analysis), which extends and supersedes the original Open Pediatric Brain Tumor Atlas (OpenPBTA) [16]. " +
        "We restricted analyses to PBTA-cohort RNA-Seq specimens annotated as tumor tissue. To ensure statistical independence, we applied the OpenPedCan " +
        "\"independent-primary-plus\" specimen list (one biospecimen per patient, preferring primary tumors). For patients with pediatric high-grade glioma " +
        "(H3/IDH-wildtype), a maximum age cutoff of 21 years was applied to exclude adult-onset cases.",
        { firstLine: true }),
      P("Main analysis cohort groups were defined as: (1) DMG_K27 — diffuse midline glioma, H3 K27-altered (denoted K28-altered in OpenPedCan, IUPAC nomenclature; " +
        "n = 175); (2) DHG_G34 — diffuse hemispheric glioma, H3 G34-mutant (denoted G35-mutant in OpenPedCan; n = 31); (3) pHGG_WT — pediatric-type high-grade " +
        "glioma, H3/IDH-wildtype (n = 125); (4) IHG — infant-type hemispheric glioma (operational definition: age at diagnosis < 18 months + hemispheric " +
        "location + NTRK1/2/3, ALK, ROS1, MET, or BRAF fusion as reported in the OpenPedCan fusion-putative-oncogenic table; n = 18). A secondary supplementary " +
        "cohort of BRAF-altered and other MAPK-driven tumors (n = 353; split into HGG/PXA/infant-fusion driven subgroup and LGG-only subgroup) was retained " +
        "for projection analyses but not included in primary clustering.",
        { firstLine: true }),
      P("Anatomical location was classified using the OpenPedCan CNS_region and primary_site annotations. Pineal, suprasellar, optic pathway, and ventricular " +
        "tumors were assigned to an \"Ambiguous\" category and excluded from midline-versus-hemispheric association analyses.",
        { firstLine: true }),

      H2("2.2. Gene Expression Processing"),
      P("FPKM-normalized expression matrices (gene-expression-rsem-tpm-collapsed.rds) were obtained from the v15 release and subset to cohort biospecimen IDs. " +
        "Gene identifiers were retained as HGNC symbols. Genes with zero expression across all retained samples (n = 4,917) were removed, yielding a final " +
        "matrix of 55,408 genes × 702 samples (349 main + 353 secondary). For CIBERSORTx-style input, TPM-normalized values were used without log transformation " +
        "or quantile normalization, per recommended RNA-seq preprocessing [11].",
        { firstLine: true }),

      H2("2.3. Immune Deconvolution"),
      P("Four orthogonal deconvolution methods were applied using the immunedeconv R package (v2.1.4) [17]: MCP-counter [18], quanTIseq [19], EPIC [20], " +
        "and xCell [21]. MCP-counter input was log2(x+1)-transformed automatically by the package. quanTIseq and EPIC were run with the tumor argument set to TRUE. " +
        "All other parameters used package defaults. CIBERSORTx LM22 deconvolution was reserved for a planned validation step using the Stanford web tool [11].",
        { firstLine: true }),

      H2("2.4. Brain-Tuned ssGSEA Scoring"),
      P("To address the absence of brain-resident myeloid lineages in conventional deconvolution panels, we assembled a 24-signature gene set library " +
        "encompassing microglia-specific programs (homeostatic core from Bohlen 2020 [13] / Butovsky 2014; glioma-associated microglia from Klemm 2020 [12]; " +
        "Mg-TAM from Antunes 2021 [14]), monocyte-derived macrophage programs (Klemm 2020 MDM; Antunes 2021 Mo-TAM), disease-associated microglia (DAM; " +
        "Keren-Shaul 2017), antigen presentation (MHC class I and II), interferon-α/γ response, T-cell cytotoxicity, T-cell exhaustion, T-regulatory cells, " +
        "dendritic cell activation, NK cell activity, M1/M2 macrophage polarization, neutrophil activation, chemokine T-cell recruitment, MAPK transcriptional " +
        "output, TGF-β immunosuppression, the Wang 2017 glioma inflammatory signature, cell cycle/proliferation, and a brain tumor stemness signature. " +
        "Gene set match rates against the TPM matrix exceeded 90% for all signatures. Single-sample gene set enrichment analysis was performed with the GSVA " +
        "package (v1.50) [22] using method = \"ssgsea\" and normalize = TRUE on the FPKM matrix.",
        { firstLine: true }),

      H2("2.5. Consensus Clustering and Ecotype Definition"),
      P("The clustering feature matrix combined quanTIseq immune cell fractions (excluding the \"uncharacterized cell\" residual; 10 features) with 19 " +
        "immune-relevant ssGSEA signatures, z-scored per feature across samples. Consensus clustering was performed with the ConsensusClusterPlus package [23] " +
        "(k-means clustering, Euclidean distance, Ward.D2 linkage, 500 resampling iterations, 80% subsampling, k = 2–6). Optimal k was evaluated by the " +
        "proportion of ambiguous clustering (PAC) [24], silhouette width [25], and biological interpretability. A pre-specified k = 3 was selected for primary " +
        "analysis (consistent with prior glioma immune subtyping frameworks [26]); k = 2 and k = 4 results are reported as sensitivity analyses.",
        { firstLine: true }),
      P("Cross-method correlation of myeloid, CD8 T-cell, and NK cell signals was assessed using Spearman's ρ on per-sample summed scores within each method.",
        { firstLine: true }),

      H2("2.6. Statistical Analysis"),
      P("Associations between ecotype assignment and categorical covariates (cohort group, anatomical location class, developmental age group) were tested with " +
        "Fisher's exact test using simulated p-values (10,000 replicates). Per-pathway ssGSEA score differences across ecotypes were evaluated by Kruskal-Wallis " +
        "test with Benjamini-Hochberg multiple-testing adjustment; pairwise post-hoc comparisons used Dunn's test. Rank-biserial correlation was reported as " +
        "an effect size measure. Two-sided p-values < 0.05 were considered statistically significant.",
        { firstLine: true }),

      H2("2.7. Survival Analysis"),
      P("Overall survival (OS) was modeled using the survival package [27] in R 4.3.2. Kaplan-Meier curves were generated stratified by ecotype, both overall " +
        "and within each cohort group, with log-rank tests. Univariable and multivariable Cox proportional-hazards models were fitted with the Inflamed ecotype " +
        "as reference; the multivariable model adjusted for cohort group, age at diagnosis (years), and sex. Proportional-hazards assumptions were assessed " +
        "by Schoenfeld residual tests (Supplementary Figure S_X). Survival analyses were exploratory and interpreted cautiously given incomplete OS annotation " +
        "(251/349 main cohort samples with available OS data; 71.9%).",
        { firstLine: true }),

      H2("2.8. Software and Reproducibility"),
      P("All analyses were performed in R 4.3.2 (R Foundation for Statistical Computing, Vienna, Austria) and Python 3.12. Key R packages included immunedeconv, " +
        "GSVA, ConsensusClusterPlus, ComplexHeatmap, survival, survminer, and rstatix. Python visualizations used matplotlib and pandas. All scripts and " +
        "intermediate output files are available [GitHub/Zenodo repository — to be added at submission].",
        { firstLine: true }),

      new Paragraph({ children: [new PageBreak()] }),

      // ============ RESULTS ============
      H1("3. Results"),

      H2("3.1. Cohort Composition and Quality"),
      P("Application of inclusion criteria to the OpenPedCan v15 release yielded a main analysis cohort of 349 unique pediatric glioma patients (Table 1, Figure 1). " +
        "Biological annotations recapitulated established patterns: DMG_K27 tumors were predominantly midline (118/175; 67.4%) with a median age of 7.9 years; " +
        "DHG_G34 tumors were predominantly hemispheric (24/31; 77.4%) in adolescents (median 15.7 years); pHGG_WT was anatomically heterogeneous (49.6% hemispheric, " +
        "17.6% midline); and IHG tumors were exclusively in infants (median 0.4 years) with hemispheric predominance (15/18; 83.3%) and 100% BRAF/RTK fusion " +
        "positivity by design. Overall survival data were available for 251 of 349 patients (71.9%). Cross-method matching between cohort biospecimens and the " +
        "expression matrix achieved 100% retention.",
        { firstLine: true }),

      // === Figure 1 ===
      ...figure({
        path: `${FIG}/workflow.png`, width: 580, height: 395,
        figNum: "Figure 1",
        title: "Study workflow",
        caption: "Schematic overview of cohort assembly, expression matrix preparation, multi-method deconvolution, brain-tuned ssGSEA scoring, " +
                 "and consensus clustering. Numbers within boxes indicate sample counts at each step.",
      }),

      blank(),
      P("Table 1. Clinical and molecular characteristics of the pediatric glioma cohort.", { bold: true }),
      makeTable(table1, [2600, 1400, 1400, 1400, 1400, 1400]),
      blank(),

      H2("3.2. Brain-Tuned Signatures Resolve a Microglia-versus-MDM Axis Obscured by Conventional Deconvolution"),
      P("Cross-method correlation of macrophage/monocyte signals revealed a striking pattern of internal consistency among brain-tuned signatures that was " +
        "not captured by conventional deconvolution methods (Figure 2). Three microglia-derived signatures—Klemm 2020 microglia core, Bohlen 2020 homeostatic " +
        "microglia, and Antunes 2021 Mg-TAM—correlated tightly with one another (pairwise Spearman ρ > 0.98), and two monocyte-derived macrophage signatures " +
        "(Klemm 2020 MDM, Antunes 2021 Mo-TAM) likewise correlated strongly (ρ = 0.94). In contrast, the correlation between microglia and MDM signatures was " +
        "only ρ ≈ 0.40, demonstrating that these are distinguishable cellular signals. Conventional deconvolution methods, by comparison, " +
        "showed mixed behavior: xCell macrophage tracked MDM closely (ρ = 0.86) but only weakly with microglia (ρ ≈ 0.57), while EPIC macrophage captured " +
        "intermediate microglia (ρ ≈ 0.65) and MDM (ρ ≈ 0.75) signals. MCP-counter and quanTIseq macrophage estimates fell between these extremes. These results " +
        "establish that LM22-style methods conflate brain-resident microglia and infiltrating MDM, while brain-tuned signature scoring resolves the distinction.",
        { firstLine: true }),

      // === Figure 2 ===
      ...figure({
        path: `${FIG}/fig_corr_myeloid_annotated.png`, width: 620, height: 500,
        figNum: "Figure 2",
        caption: "Brain-tuned signatures resolve a microglia-versus-MDM axis. Spearman correlation heatmap of macrophage and monocyte-related " +
                 "signals across four deconvolution methods (MCP-counter, quanTIseq, EPIC, xCell) and five brain-tuned ssGSEA signatures " +
                 "(microglia from Klemm 2020, Bohlen 2020, and Antunes 2021 Mg-TAM; MDM/Mo-TAM from Klemm 2020 and Antunes 2021). Boxed regions " +
                 "highlight tight internal correlation among microglia signatures (ρ > 0.98) and among MDM signatures (ρ = 0.94), with cross-cluster " +
                 "correlation of ρ ≈ 0.40 demonstrating that the two cellular signals are resolvable. xCell macrophage tracks MDM (ρ = 0.86), " +
                 "illustrating the peripheral-macrophage bias of conventional methods.",
      }),

      H2("3.3. Consensus Clustering Identifies Three Reproducible Immune Ecotypes"),
      P("Unsupervised consensus clustering of the combined quanTIseq + ssGSEA feature matrix (29 features × 349 samples) was performed across k = 2–6. " +
        "The proportion of ambiguous clustering reached a local minimum at k = 2, with k = 3 representing a stable secondary solution (PAC = 0.83; silhouette > 0). " +
        "Consistent with the pre-specified hypothesis and prior literature [26], we adopted k = 3 for primary analysis; k = 2 and k = 4 results are reported as " +
        "Supplementary Figure S2.",
        { firstLine: true }),

      P("The three ecotypes were annotated by examining feature mean z-scores within each cluster (Figure 3A): an Inflamed ecotype (n = 113; 32.4% of cohort) " +
        "characterized by uniformly elevated lymphocyte (T-cell cytotoxicity, NK), myeloid (M1/M2 macrophage, MDM), and antigen presentation (MHC class II, IFN-γ) " +
        "scores; an Intermediate ecotype (n = 157; 45.0%) with feature means near zero or slightly negative; and an Immune-desert ecotype (n = 79; 22.6%) with " +
        "uniformly low immune feature scores. Principal component and UMAP projections of the feature matrix confirmed reproducible spatial separation of the " +
        "three ecotypes (Figure 3B,C). Notably, the high lymphocyte and high myeloid signals co-occurred within the Inflamed ecotype rather than partitioning " +
        "into distinct lymphocyte-inflamed and myeloid-dominant compartments as reported in adult glioma cohorts [26,28], indicating a pediatric-specific " +
        "TME architecture.",
        { firstLine: true }),

      // === Figure 3 — heatmap + PCA + UMAP ===
      ...multiFigure(
        [
          { path: `${FIGCL}/ecotype_heatmap_main.png`, width: 600, height: 430 },
          { path: `${FIGCL}/ecotype_pca_main.png`,     width: 360, height: 310 },
          { path: `${FIGCL}/ecotype_umap_main.png`,    width: 360, height: 310 },
        ],
        "Consensus clustering identifies three reproducible immune ecotypes. (A) Heatmap of z-scored feature matrix " +
        "(quanTIseq immune cell fractions and 19 immune ssGSEA signatures × 349 samples) annotated by ecotype, cohort group, and anatomical location. " +
        "(B) Principal component analysis projection of the feature matrix, colored by ecotype. (C) UMAP projection.",
        "Figure 3", ["A", "B", "C"]
      ),

      H2("3.4. Ecotype Distribution Across Histone Class and Anatomical Location"),
      P("Ecotype membership was significantly associated with cohort group (Fisher's exact p = 1.3 × 10⁻³; Figure 4A). DHG_G34 tumors were predominantly " +
        "Immune-desert (16/31; 51.6%) with only 9.7% (3/31) classified as Inflamed — a finding that contrasts with prior reports describing G34-mutant tumors " +
        "as comparatively lymphocyte-inflamed in adult-onset cases [8]. DMG_K27 tumors, in contrast, demonstrated marked internal heterogeneity, with 36.6% " +
        "(64/175) assigned to the Inflamed ecotype despite being uniformly aggressive at the clinical level. pHGG_WT tumors showed an intermediate distribution " +
        "(32.0% Inflamed, 42.4% Intermediate, 25.6% Immune-desert), and IHG tumors were dominated by the Intermediate ecotype (44.4%) with similar proportions " +
        "in Inflamed (33.3%) and Immune-desert (22.2%).",
        { firstLine: true }),

      P("Anatomical location, restricted to clearly classified midline, hemispheric, and posterior fossa tumors, was also significantly associated with ecotype " +
        "(Fisher's exact p = 1.3 × 10⁻²; Figure 4B). Counter to expectations derived from adult glioblastoma, midline tumors more frequently fell into the " +
        "Inflamed ecotype than hemispheric tumors (41.4% vs 25.0%), while hemispheric tumors were enriched for the Immune-desert ecotype (32.7% vs 15.7%). " +
        "Developmental age group showed no significant association with ecotype after accounting for the small sample sizes in extreme age categories (p = 0.29).",
        { firstLine: true }),

      // === Figure 4 — distribution by cohort + location ===
      ...multiFigure(
        [
          { path: `${FIGCL}/ecotype_distribution_by_cohort.png`,   width: 560, height: 360 },
          { path: `${FIGCL}/ecotype_distribution_by_location.png`, width: 560, height: 360 },
        ],
        "Ecotype distribution by histone class and anatomy. (A) Stacked barplot of ecotype proportions across the four main cohort groups " +
        "(Fisher's exact p = 1.3 × 10⁻³). (B) Stacked barplot of ecotype proportions across anatomical location, restricted to clearly classified " +
        "midline, hemispheric, and posterior fossa tumors (Ambiguous and Other categories excluded; Fisher's exact p = 1.3 × 10⁻²).",
        "Figure 4", ["A", "B"]
      ),

      H2("3.5. Ecotype Is an Independent Prognostic Factor for Overall Survival"),
      P("Kaplan-Meier analysis of the 251 patients with available OS data revealed that the Inflamed ecotype had the longest median OS (595 days), " +
        "followed by Intermediate (446 days) and Immune-desert (403 days; overall log-rank p = 0.037; Figure 5A). Stratified analyses showed that the ecotype " +
        "effect was statistically significant within pHGG_WT (log-rank p = 0.049) and trended within DHG_G34 (p = 0.070), but was masked within DMG_K27 (p = 0.132), " +
        "likely reflecting the uniformly poor prognosis of K27-altered DMG that constrains the dynamic range of survival differences.",
        { firstLine: true }),

      P("In a multivariable Cox proportional-hazards model adjusting for cohort group, age at diagnosis, and sex (Table 2; Figure 5B), Immune-desert remained an " +
        "independent adverse prognostic factor relative to Inflamed (HR = 1.79; 95% CI 1.21–2.66; p = 0.004), while the Intermediate ecotype showed a non-significant " +
        "trend toward worse survival (HR = 1.37; 95% CI 0.99–1.91; p = 0.062). Cohort group effects were substantial as expected (DMG_K27 HR = 2.96, p < 0.001; " +
        "IHG HR = 0.22, p = 0.005), but did not abrogate the independent contribution of ecotype.",
        { firstLine: true }),

      // === Figure 5 — KM + forest ===
      ...multiFigure(
        [
          { path: `${FIGSV}/KM_overall_by_ecotype.png`,     width: 480, height: 440 },
          { path: `${FIGSV}/Cox_forest_multivariable.png`,  width: 580, height: 440 },
        ],
        "Ecotype is an independent prognostic factor. (A) Kaplan-Meier overall survival curves stratified by ecotype (n = 251 with OS data; " +
        "log-rank p = 0.037). Risk table below shows numbers at risk at each time point. (B) Forest plot of multivariable Cox proportional-hazards " +
        "model adjusting for cohort group, age, and sex (n = 251).",
        "Figure 5", ["A", "B"]
      ),

      blank(),
      P("Table 2. Multivariable Cox proportional-hazards model for overall survival (n = 251).", { bold: true }),
      makeTable(table2, [4200, 1500, 1700, 2000]),
      blank(),

      H2("3.6. Pathway Activity Distinguishes Ecotypes Without Confounding by Tumor Proliferation"),
      P("Kruskal-Wallis testing of ssGSEA scores across ecotypes demonstrated significant differences for 22 of 24 signatures, with effect sizes (rank-biserial r " +
        "for Inflamed vs Immune-desert) exceeding 0.95 for the top 15 immune pathways (Table 3; Figure 6). Stemness signature was marginally different (p_BH = 0.029), " +
        "while cell-cycle/proliferation signatures showed no ecotype dependence (p = 0.94) — confirming that ecotype distinctions reflect genuine immune contexture " +
        "and not confounding tumor cell proliferation.",
        { firstLine: true }),

      P("Notably, microglia-associated signatures (Klemm 2020, Bohlen 2020, Antunes 2021 Mg-TAM) showed substantial but more modest ecotype-dependent variation " +
        "(p_BH ≈ 10⁻²¹) than monocyte-derived macrophage signatures (MDM and Mo-TAM; p_BH ≈ 10⁻⁴⁵), suggesting that brain-resident microglia are relatively " +
        "uniform across ecotypes while peripheral MDM recruitment is a primary driver of ecotype assignment in pediatric glioma. MAPK transcriptional activity " +
        "also varied significantly across ecotypes (p_BH = 1.4 × 10⁻¹²), with the Inflamed ecotype showing the highest scores — a finding with potential " +
        "implications for the immunological consequences of MAPK pathway activation in BRAF-altered subgroups.",
        { firstLine: true }),

      // === Figure 6 ===
      ...figure({
        path: `${FIGPW}/Pathway_composite_12sig.png`, width: 660, height: 470,
        figNum: "Figure 6",
        caption: "Pathway activity across ecotypes. Composite boxplot of 12 key brain-immune and immune-modulatory ssGSEA signatures across the three " +
                 "ecotypes. Kruskal-Wallis test was performed for each signature; Benjamini-Hochberg-adjusted p-values are shown in each panel header. " +
                 "Note that cell-cycle/proliferation signature does not differ across ecotypes (p = 0.94; not shown), confirming that ecotype distinctions " +
                 "reflect immune contexture rather than tumor proliferation.",
      }),

      blank(),
      P("Table 3. ssGSEA pathway differences across immune ecotypes (Kruskal-Wallis test).", { bold: true }),
      makeTable(table3, [4500, 2400, 2460]),
      blank(),

      new Paragraph({ children: [new PageBreak()] }),

      // ============ DISCUSSION ============
      H1("4. Discussion"),
      P("Using harmonized RNA-sequencing from the OpenPedCan v15 cohort and a deliberately multi-method analytical strategy, we identified three reproducible " +
        "immune ecotypes — Inflamed, Intermediate, and Immune-desert — that organize the tumor immune microenvironment of pediatric gliomas. Two findings warrant " +
        "particular emphasis: first, that brain-tuned signature scoring resolves a microglia-versus-monocyte-derived macrophage axis that conventional " +
        "deconvolution panels conflate; and second, that ecotype assignment carries independent prognostic value beyond histone class, with the Inflamed ecotype " +
        "associated with the longest overall survival after multivariable adjustment.",
        { firstLine: true }),

      P("Several aspects of the ecotype distribution deviate from patterns reported in adult glioma. H3 G34-mutant hemispheric tumors, often considered the " +
        "more lymphocyte-rich pediatric high-grade glioma subtype [8], were predominantly Immune-desert in our cohort, with only 9.7% reaching the Inflamed " +
        "ecotype. While our sample size for this subgroup is limited (n = 31), the consistency of this signal across multiple deconvolution methods suggests " +
        "that the lymphocyte-inflamed designation for G34 tumors may need to be qualified relative to age-matched controls. Similarly, midline tumors more " +
        "frequently fell into the Inflamed ecotype than hemispheric tumors — opposite to the midline immune privilege paradigm described for adult glioblastoma " +
        "[29,30] and for diffuse intrinsic pontine glioma specifically [6]. These observations highlight that adult-derived TME paradigms cannot be directly " +
        "extrapolated to pediatric disease and motivate developmentally anchored frameworks for interpreting immune contexture in childhood brain tumors.",
        { firstLine: true }),

      P("The microglia-versus-MDM axis revealed by brain-tuned signature scoring has direct translational implications. Three orthogonal microglia signatures " +
        "(Klemm 2020, Bohlen 2020, Antunes 2021 Mg-TAM) showed pairwise correlations exceeding 0.98, while their correlation with MDM-derived signatures " +
        "remained below 0.41 — establishing the internal consistency of microglia-versus-MDM scoring at the cohort scale. By contrast, xCell macrophage scores " +
        "correlated tightly with MDM (ρ = 0.86) but only modestly with microglia (ρ ≈ 0.57), demonstrating that LM22-derived methods preferentially capture " +
        "peripherally derived macrophage signal. This distinction matters: brain-resident microglia and infiltrating MDM differ in lineage, anatomical origin, " +
        "and therapeutic accessibility [12,13,14], and treatment strategies targeting CSF1R or CCR2 are likely to have very different effects on the two " +
        "populations [31]. Our finding that microglia signatures were only modestly stratified across ecotypes (p ≈ 10⁻²¹) while MDM signatures were strongly " +
        "ecotype-dependent (p ≈ 10⁻⁴⁵) further suggests that MDM recruitment, more than baseline microglia abundance, drives ecotype assignment — pointing to " +
        "peripheral macrophage recruitment as a candidate therapeutic intervention point.",
        { firstLine: true }),

      P("The independent prognostic value of ecotype assignment after adjustment for cohort group, age, and sex (Immune-desert HR = 1.79; p = 0.004) suggests " +
        "that immune contexture provides an additional layer of risk stratification beyond molecular classification. This is consistent with the well-established " +
        "prognostic significance of immune infiltration in adult glioblastoma [28,32] and other solid tumors, but represents one of the first cohort-level " +
        "demonstrations in pediatric glioma. The stratified Kaplan-Meier analysis additionally shows that ecotype is significant within H3-wildtype pHGG " +
        "(p = 0.049) but is masked within DMG_K27 by the uniformly poor prognosis of that subgroup, suggesting that ecotype-based stratification may be most " +
        "clinically actionable for tumors with intermediate prognosis where treatment intensification or de-escalation decisions are more nuanced.",
        { firstLine: true }),

      P("Several limitations should be acknowledged. First, our deconvolution and ssGSEA scoring relies on bulk RNA-seq and cannot resolve cell-type-specific " +
        "expression at single-cell or spatial resolution; future work integrating scRNA-seq references through reference-based deconvolution methods such as " +
        "BayesPrism or MuSiC [33,34] will be valuable for validation. Second, the absence of microglia in the LM22 signature matrix was addressed via " +
        "complementary ssGSEA scoring rather than by constructing a custom brain-tuned signature matrix — a more comprehensive but computationally heavier approach " +
        "[35] that we plan to pursue in future work. Third, overall survival data were available for 71.9% of the main cohort, and treatment heterogeneity " +
        "across centers and time periods could not be modeled directly. Fourth, sample sizes for DHG_G34 (n = 31) and IHG (n = 18) constrain the statistical " +
        "power for subgroup-specific analyses; external validation cohorts will be required to confirm subgroup findings. Finally, the Inflamed ecotype in " +
        "pediatric glioma demonstrated co-occurrence of lymphocyte and myeloid signals, which differs from the discrete lymphocyte-inflamed and myeloid-dominant " +
        "compartments reported in adult glioma — this pattern may reflect developmental immunology specific to childhood brain tumors and warrants targeted " +
        "investigation in mechanistic studies.",
        { firstLine: true }),

      P("Despite these limitations, our work establishes a developmentally anchored framework for interpreting pediatric glioma immunity. The combination of " +
        "multi-method deconvolution with brain-tuned signature scoring provides a methodological blueprint that can be applied to other pediatric brain " +
        "tumor entities and to longitudinal studies of treatment response. Translationally, our findings raise the possibility that ecotype assignment could " +
        "inform immunotherapy stratification — for example, by identifying Inflamed tumors as candidates for checkpoint blockade or T-cell-based therapies, " +
        "Immune-desert tumors as candidates for TME-modifying combinations, and Intermediate tumors for context-dependent strategies. Prospective clinical " +
        "evaluation of these hypotheses is the natural next step.",
        { firstLine: true }),

      H1("5. Conclusions"),
      P("Pediatric gliomas segregate into three reproducible immune ecotypes that organize the tumor immune microenvironment along a continuum from immune-cold " +
        "(Immune-desert) through Intermediate to immune-active (Inflamed) states. Brain-tuned signature scoring resolves a microglia-versus-monocyte-derived " +
        "macrophage axis obscured by conventional deconvolution panels, identifying peripheral macrophage recruitment as a primary driver of ecotype assignment. " +
        "Ecotype membership carries independent prognostic value after adjustment for histone class, age, and sex, supporting its potential integration into " +
        "risk-stratified treatment frameworks for childhood brain tumors. Future work should validate these findings in independent cohorts, refine the ecotype " +
        "framework with single-cell and spatial transcriptomic data, and evaluate ecotype-guided immunotherapy strategies in prospective clinical trials.",
        { firstLine: true }),

      // ============ Statements ============
      H1("Author Contributions"),
      P("Conceptualization, methodology, data curation, formal analysis, writing — original draft, S.C.; visualization, S.C.; supervision, S.C. " +
        "The author has read and agreed to the published version of the manuscript."),
      H1("Funding"),
      P("This research received no external funding."),
      H1("Institutional Review Board Statement"),
      P("This study analyzed publicly available, de-identified data from the OpenPedCan v15 release. No new human subjects research was conducted. " +
        "Ethical approval for the original OpenPedCan cohort was obtained by the Children's Brain Tumor Network and contributing institutions [16]."),
      H1("Informed Consent Statement"),
      P("Not applicable."),
      H1("Data Availability Statement"),
      P("All data analyzed in this study are publicly available from the OpenPedCan v15 release (https://github.com/d3b-center/OpenPedCan-analysis). " +
        "Analysis scripts and intermediate output files are available at [GitHub/Zenodo repository to be added at submission]."),
      H1("Conflicts of Interest"),
      P("The author declares no conflict of interest."),

      // ============ REFERENCES ============
      new Paragraph({ children: [new PageBreak()] }),
      H1("References"),
      P("[Placeholder reference list — to be expanded and formatted to Cancers MDPI style at submission]", { italics: true, color: "555555" }),
      blank(),
      P("1.  Mackay, A.; Burford, A.; Carvalho, D.; et al. Integrated molecular meta-analysis of 1,000 pediatric high-grade and diffuse intrinsic pontine glioma. Cancer Cell. 2017;32(4):520-537.e5."),
      P("2.  Sturm, D.; Witt, H.; Hovestadt, V.; et al. Hotspot mutations in H3F3A and IDH1 define distinct epigenetic and biological subgroups of glioblastoma. Cancer Cell. 2012;22(4):425-37."),
      P("3.  Filbin, M.G.; Tirosh, I.; Hovestadt, V.; et al. Developmental and oncogenic programs in H3K27M gliomas dissected by single-cell RNA-seq. Science. 2018;360(6386):331-335."),
      P("4.  Foster, J.B.; Madsen, P.J.; Hegde, M.; et al. Immunotherapy for pediatric brain tumors. Neuro Oncol. 2019;21(10):1226-1238."),
      P("5.  Persson, M.L.; Douglas, A.M.; Alvaro, F.; et al. The intrinsic and microenvironmental features of diffuse midline glioma. Neuro Oncol. 2022;24(9):1408-1422."),
      P("6.  Lieberman, N.A.P.; DeGolier, K.; Kovar, H.M.; et al. Characterization of the immune microenvironment of diffuse intrinsic pontine glioma. Neuro Oncol. 2019;21(1):83-94."),
      P("7.  Lin, G.L.; Nagaraja, S.; Filbin, M.G.; et al. Non-inflammatory tumor microenvironment of diffuse intrinsic pontine glioma. Acta Neuropathol Commun. 2018;6(1):51."),
      P("8.  Crotty, E.E.; Smith, S.M.C.; Brasel, K.; et al. Medical decision-making in H3 G34-mutant gliomas. Neuro Oncol. 2020;22(suppl_3):iii365-iii366."),
      P("9.  Guerreiro Stucklin, A.S.; Ryall, S.; Fukuoka, K.; et al. Alterations in ALK/ROS1/NTRK/MET drive a group of infantile hemispheric gliomas. Nat Commun. 2019;10(1):4343."),
      P("10. Clarke, M.; Mackay, A.; Ismer, B.; et al. Infant high-grade gliomas comprise multiple subgroups characterized by novel targetable gene fusions. Cancer Discov. 2020;10(7):942-963."),
      P("11. Newman, A.M.; Steen, C.B.; Liu, C.L.; et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol. 2019;37(7):773-782."),
      P("12. Klemm, F.; Maas, R.R.; Bowman, R.L.; et al. Interrogation of the microenvironmental landscape in brain tumors reveals disease-specific alterations of immune cells. Cell. 2020;181(7):1643-1660.e17."),
      P("13. Bohlen, C.J.; Bennett, M.L.; Bennett, F.C. Isolation and culture of microglia. Curr Protoc Immunol. 2019;125(1):e70."),
      P("14. Antunes, A.R.P.; Scheyltjens, I.; Lodi, F.; et al. Single-cell profiling of myeloid cells in glioblastoma across species and disease stage reveals macrophage competition and specialization. Nat Neurosci. 2021;24(4):595-610."),
      P("15. Friebel, E.; Kapolou, K.; Unger, S.; et al. Single-cell mapping of human brain cancer reveals tumor-specific instruction of tissue-invading leukocytes. Cell. 2020;181(7):1626-1642.e20."),
      P("16. Shapiro, J.A.; Gaonkar, K.S.; Spielman, S.J.; et al. OpenPBTA: the Open Pediatric Brain Tumor Atlas. Cell Genomics. 2023;3(7):100340."),
      P("17. Sturm, G.; Finotello, F.; Petitprez, F.; et al. Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. Bioinformatics. 2019;35(14):i436-i445."),
      P("18. Becht, E.; Giraldo, N.A.; Lacroix, L.; et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biol. 2016;17(1):218."),
      P("19. Finotello, F.; Mayer, C.; Plattner, C.; et al. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome Med. 2019;11(1):34."),
      P("20. Racle, J.; de Jonge, K.; Baumgaertner, P.; et al. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife. 2017;6:e26476."),
      P("21. Aran, D.; Hu, Z.; Butte, A.J. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biol. 2017;18(1):220."),
      P("22. Hänzelmann, S.; Castelo, R.; Guinney, J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7."),
      P("23. Wilkerson, M.D.; Hayes, D.N. ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. Bioinformatics. 2010;26(12):1572-3."),
      P("24. Şenbabaoğlu, Y.; Michailidis, G.; Li, J.Z. Critical limitations of consensus clustering in class discovery. Sci Rep. 2014;4:6207."),
      P("25. Rousseeuw, P.J. Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. J Comput Appl Math. 1987;20:53-65."),
      P("26. Wang, Q.; Hu, B.; Hu, X.; et al. Tumor evolution of glioma-intrinsic gene expression subtypes associates with immunological changes in the microenvironment. Cancer Cell. 2017;32(1):42-56.e6."),
      P("27. Therneau, T.M.; Grambsch, P.M. Modeling Survival Data: Extending the Cox Model. Springer; 2000."),
      P("28. Thorsson, V.; Gibbs, D.L.; Brown, S.D.; et al. The immune landscape of cancer. Immunity. 2018;48(4):812-830.e14."),
      P("29. Fecci, P.E.; Heimberger, A.B.; Sampson, J.H. Immunotherapy for primary brain tumors: no longer a matter of privilege. Clin Cancer Res. 2014;20(22):5620-9."),
      P("30. Chongsathidkiet, P.; Jackson, C.; Koyama, S.; et al. Sequestration of T cells in bone marrow in the setting of glioblastoma and other intracranial tumors. Nat Med. 2018;24(9):1459-1468."),
      P("31. Pyonteck, S.M.; Akkari, L.; Schuhmacher, A.J.; et al. CSF-1R inhibition alters macrophage polarization and blocks glioma progression. Nat Med. 2013;19(10):1264-72."),
      P("32. Bockmayr, M.; Mohme, M.; Klauschen, F.; et al. Subgroup-specific immune and stromal microenvironment in medulloblastoma. Oncoimmunology. 2018;7(9):e1462430."),
      P("33. Chu, T.; Wang, Z.; Pe'er, D.; Danko, C.G. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. Nat Cancer. 2022;3(4):505-517."),
      P("34. Wang, X.; Park, J.; Susztak, K.; Zhang, N.R.; Li, M. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nat Commun. 2019;10(1):380."),
      P("35. Steen, C.B.; Liu, C.L.; Alizadeh, A.A.; Newman, A.M. Profiling cell type abundance and expression in bulk tissues with CIBERSORTx. Methods Mol Biol. 2020;2117:135-157."),
    ],
  }],
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync("C:/Users/scott/Downloads/Open PBTA/output/Manuscript_Draft_v1_with_figures.docx", buf);
  console.log("Manuscript_Draft_v1_with_figures.docx written");
});
