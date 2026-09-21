import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
A=lambda **k: Font(name="Arial",**k)
HDR=PatternFill("solid",fgColor="E8EDF5"); WARN=PatternFill("solid",fgColor="FBE9E9")
OKF=PatternFill("solid",fgColor="EEF4EA"); NEWF=PatternFill("solid",fgColor="FFF8E7")
THIN=Side(style="thin",color="BFBFBF")

MAIN=[  # label, status, description, where it is/should be cited in the body
 ("Figure 1A","revised legend","Study workflow (unchanged artwork)","Results, Cohort overview, p.12-13 — currently cited as 'Figure 1'; change to 'Figure 1A'"),
 ("Figure 1B","NEW","Analysis-set derivation; branching because the location and survival subsets are not nested","Results, Cohort overview, p.13 — add citation; also cite in Methods, Cohort assembly"),
 ("Figure 2A","unchanged","CIBERSORTx LM22 quality control, n=349","Results, deconvolution paragraph, p.13 — add panel-level citation"),
 ("Figure 2B","unchanged","Cross-method concordance, n=349","Results, deconvolution paragraph, p.13 — add panel-level citation"),
 ("Figure 2C","unchanged","Brain-tuned myeloid signatures vs LM22 axes","Results, deconvolution paragraph, p.13 — add panel-level citation"),
 ("Figure 2D","unchanged","Cross-method estimates by ecotype","Results, deconvolution paragraph, p.13 — add panel-level citation"),
 ("Figure 3A","REPLACED","Cluster-number selection, k = 1-10, five panels incl. gap statistic and cluster-size floor","Results, Consensus clustering, p.13-14 — cite explicitly in the rewritten paragraph (INS-4)"),
 ("Figure 3B","unchanged","PCA and UMAP by ecotype","Results, Consensus clustering, p.14 — add panel-level citation"),
 ("Figure 3C","unchanged","Top discriminating features","Results, Ecotypes represented coordinated immune programs — add panel-level citation"),
 ("Figure 3D","unchanged","Ecotype distribution by integrated molecular group","Results, Ecotype distribution paragraph — add panel-level citation"),
 ("Figure 4","unchanged","Location and molecular-group statistics; 34-feature PERMANOVA","Results, Location paragraph, p.15 — already cited"),
 ("Figure 5A","unchanged","Five theme composites","Results, theme-level paragraph, p.17 — add panel-level citation"),
 ("Figure 5B","unchanged","Per-signature z-score heatmap","Results, theme-level paragraph — add panel-level citation"),
 ("Figure 5C","unchanged","Molecular-group-stratified theme scores","Results, theme-level paragraph — add panel-level citation"),
 ("Figure 5D","unchanged","Twelve focus-signature boxplots","Results, theme-level paragraph — add panel-level citation"),
 ("Figure 5E","unchanged","Dunn post-hoc heatmap","Results, theme-level paragraph — add panel-level citation"),
 ("Figure 6A","REPLACED","Adjusted differential expression; reference group named in every panel title","Results, DEG paragraph, p.17 — cite as Figure 6A (INS-5)"),
 ("Figure 6B","unchanged","Enrichment from adjusted gene lists","Results, DEG paragraph, p.17 — add panel-level citation"),
 ("Figure 6C","NEW","Contrast structure: nesting, gradient slope, paired effect sizes","Results, DEG paragraph, p.17 — cite at the end of the inserted sentences (INS-5)"),
 ("Figure 7A","legend amended","Overall Kaplan-Meier; legend now states the curves are unadjusted","Results, survival paragraph, p.17-18 — add panel-level citation"),
 ("Figure 7B","unchanged","Adequately sized DMG and pHGG strata","Results, survival paragraph, p.18 — add panel-level citation"),
 ("Figure 7C","unchanged","Multivariable Cox forest plot","Results, survival paragraph, p.18 — add panel-level citation"),
]
SUPPF=[ # label, status, description, cited?, where to cite
 ("Figure S1","unchanged","Same-resource pooled cross-method sensitivity, n=702","NOT CITED","Results, deconvolution paragraph, p.13"),
 ("Figure S2","unchanged","Cross-method consensus heatmap","NOT CITED","Results, deconvolution paragraph, p.13"),
 ("Figure S3","unchanged","ssGSEA TPM-harmonisation concordance","NOT CITED","Methods, gene-expression processing, p.8; or Results, deconvolution paragraph"),
 ("Figure S4","unchanged","LM22 quality control in the pooled set","NOT CITED","Results, deconvolution paragraph, p.13"),
 ("Figure S5","unchanged","Corrected 34-feature ecotype / molecular-group PERMANOVA","NOT CITED","Results, Location paragraph, p.16"),
 ("Figure S6","unchanged","High-confidence LM22 re-clustering","NOT CITED","Results, Consensus clustering, p.14"),
 ("Figure S7A","unchanged","Pooled n=702 re-clustering sensitivity","NOT CITED","Results, Sensitivity analyses, p.18"),
 ("Figure S7B","unchanged","Pooled n=702 PCA","NOT CITED","Results, Sensitivity analyses, p.18"),
 ("Figure S7C","unchanged","Main cohort vs pooled cross-method comparison","NOT CITED","Results, Sensitivity analyses, p.18"),
 ("Figure S8","unchanged","46-feature and 34-feature sensitivity incl. confusion matrix","cited p.16","already cited — keep"),
 ("Figure S9","unchanged","Tie-corrected Mann-Whitney one-versus-rest volcanoes","NOT CITED","Results, DEG paragraph, p.17 — this citation was present in the July draft and lost in the submitted version"),
 ("Figure S10A","unchanged","Subtype-stratified Kaplan-Meier","NOT CITED","Results, Sensitivity and subgroup analyses, p.18"),
 ("Figure S10B","unchanged","Subtype-stratified Cox estimates","NOT CITED","Results, Sensitivity and subgroup analyses, p.18"),
 ("Figure S11","unchanged","Location-stratified survival","NOT CITED","Results, Sensitivity and subgroup analyses, p.18"),
 ("Figure S12","unchanged","Scaled Schoenfeld residual diagnostics","NOT CITED","Results, survival paragraph where the PH assumption is reported, p.18"),
 ("Figure S13A","unchanged","Overall-survival missingness by molecular group","NOT CITED","Results, new selection-bias paragraph (INS-6), p.18"),
 ("Figure S13B","unchanged","Trimmed stabilised IPW weights","NOT CITED","Results, Sensitivity analyses, IPW sentence, p.18"),
 ("Figure S13C","unchanged","Hazard-ratio estimates across weighting schemes","NOT CITED","Results, Sensitivity analyses, IPW sentence, p.18"),
 ("Figure S13D","unchanged","Covariate balance before and after IPW","NOT CITED","Results, Sensitivity analyses, IPW sentence, p.18"),
 ("Figure S14A","unchanged","Same-resource BRAF/RTK cohort-family projection","NOT CITED","Results, Sensitivity analyses, p.18"),
 ("Figure S14B","unchanged","Same-resource MAPK activity sensitivity","NOT CITED","Results, Sensitivity analyses, p.18"),
 ("Figure S14C","unchanged","Same-resource BRAF/RTK PCA overlay","NOT CITED","Results, Sensitivity analyses, p.18"),
 ("Figure S15","unchanged","Immune-checkpoint-related transcript expression","NOT CITED","Results, Ecotypes represented coordinated immune programs"),
 ("Figure S16","unchanged","CAR-T target-related transcripts","NOT CITED","Results, Ecotypes represented coordinated immune programs"),
 ("Figure S17","unchanged","Audit copy of adjusted DEG promoted to main Figure 6","NOT CITED","Results, DEG paragraph, p.17 — or delete, since Figure 6A now supersedes it"),
 ("Figure S18","unchanged","Small DHG and IHG survival strata; descriptive","cited p.19","already cited — keep"),
 ("Figure S19","NEW","Haematopoietic stem and progenitor programs across ecotypes","NEW","Discussion, new progenitor paragraph (INS-7)"),
]
SUPPT=[
 ("Table S1","unchanged","Auxiliary molecular alteration flags","NOT CITED","Methods, cohort assembly, p.7"),
 ("Table S2A","unchanged","Immune signature gene lists","NOT CITED","Methods, CNS-relevant immune signatures, p.8"),
 ("Table S2B","unchanged","Signature membership by feature matrix","NOT CITED","Methods, consensus clustering, p.9"),
 ("Table S3","unchanged","Anatomical location classification mapping","NOT CITED","Methods, clinical and anatomical variables, p.8"),
 ("Table S4","unchanged","Location contingency analyses","NOT CITED","Results, Location paragraph, p.15-16"),
 ("Table S5","unchanged","Cross-method concordance, n=349 vs n=702","NOT CITED","Results, deconvolution paragraph, p.13"),
 ("Table S6","unchanged","k=3 sensitivity on the 46-feature matrix","NOT CITED","Results, Consensus clustering, p.14"),
 ("Table S7","unchanged","PERMANOVA, circular vs held-out feature blocks","cited p.16","already cited — keep"),
 ("Table S8","unchanged","Complete feature lists of both matrices","NOT CITED","Methods, consensus clustering, p.9"),
 ("Table S9","unchanged","Cox models and proportional-hazards diagnostics","NOT CITED","Results, survival paragraph, p.18"),
 ("Table S10","unchanged","Location statistics and midline-vs-hemispheric odds ratios","NOT CITED","Results, Location paragraph, p.16"),
 ("Table S11","unchanged","Table 1 ecotype rows and contingency counts","NOT CITED","Results, Cohort overview, p.13"),
 ("Table S12","NEW","Cluster-number sweep, k = 2 to 10","NEW","Results, Consensus clustering (INS-4)"),
 ("Table S13","NEW","Gap statistic, k = 1 to 10","NEW","Methods, cluster-number selection (INS-3); Results (INS-4)"),
 ("Table S14","NEW","Cluster-number criteria and the k each selects","NEW","Results, Consensus clustering (INS-4)"),
 ("Table S15","NEW","Overlap of the three adjusted DEG contrasts","NEW","Results, DEG paragraph (INS-5)"),
 ("Table S16","NEW","Genes shared by both Immune-desert-referenced contrasts","NEW","Results, DEG paragraph (INS-5)"),
 ("Table S17","NEW","Analysis-set reconciliation","NEW","Figure 1B legend; Results, Cohort overview"),
 ("Table S18","NEW","Survival-evaluability by baseline characteristic","NEW","Results, new selection-bias paragraph (INS-6)"),
 ("Table S19","NEW","Ecotype and survival-evaluability within molecular group","NEW","Results, new selection-bias paragraph (INS-6)"),
 ("Table S20","NEW","Haematopoietic progenitor module gene lists","NEW","Discussion, progenitor paragraph (INS-7); Figure S19 legend"),
 ("Table S21","NEW","Haematopoietic progenitor ssGSEA scores","NEW","Figure S19 legend"),
 ("Table S22","NEW","Progenitor programs by ecotype, raw and residualised","NEW","Discussion, progenitor paragraph (INS-7)"),
 ("Table S23","NEW","Dunn post-hoc comparisons for the progenitor modules","NEW","Discussion, progenitor paragraph (INS-7)"),
 ("Table S24","NEW","Software versions, parameters and random seeds","NEW","Data and code availability (INS-13)"),
]

wb=Workbook(); ws=wb.active; ws.title="README"
ws["A1"]="Figure and table crosswalk, and citation audit"; ws["A1"].font=A(size=14,bold=True,color="1E3A5F")
ws["A2"]="Cancer Informatics · CIX-26-0214 · 16 September 2026"; ws["A2"].font=A(size=10,color="595959")
notes=[
 "",
 "The Editor asked that all figures and tables, including supplementary items and their subsections, be numbered,",
 "labelled and cited chronologically in the main document. Auditing the submitted PDF against the main text shows",
 "that this is the substantive issue, not the numbering:",
 "",
 "  · Main figures are cited only at figure level (Figure 1 to Figure 7). No panel is cited individually,",
 "    although every figure except Figure 4 has panels.",
 "  · 16 of the 18 supplementary figures (24 of the 26 lettered panels) are never cited in the main text.",
 "    Only Figure S8 (p.16) and",
 "    Figure S18 (p.19) are cited. The Figure S9 citation was present in the July draft and was lost.",
 "  · 11 of the 12 supplementary table items are never cited in the main text. Only Table S7 (p.16) is cited.",
 "",
 "Numbering itself is sound and nothing needs renumbering. Supplementary figures run S1 to S18 and supplementary",
 "tables run S1 to S11 in the submitted version, so this revision adds Figure S19 and Tables S12 to S24 on the end.",
 "",
 "One file-naming mismatch to fix: the file currently named FigureS19.png is the manuscript's Supplementary",
 "Figure S18. Rename it before resubmission so that the new progenitor figure can take the name FigureS19.",
 "",
 "The three sheets that follow give, for every display item, its status in this revision and the exact place in the",
 "main text where a citation must be inserted.",
]
for i,t in enumerate(notes,start=3):
    c=ws.cell(row=i,column=1,value=t); c.font=A(size=10,bold=t.startswith("  ·")==False and t.endswith(":"))
ws.column_dimensions["A"].width=120

def sheet(name,cols,rows,widths,statuscol):
    s=wb.create_sheet(name)
    for j,h in enumerate(cols,start=1):
        c=s.cell(row=1,column=j,value=h); c.font=A(size=10,bold=True); c.fill=HDR
        c.border=Border(bottom=THIN); s.column_dimensions[get_column_letter(j)].width=widths[j-1]
    for i,r in enumerate(rows,start=2):
        for j,v in enumerate(r,start=1):
            c=s.cell(row=i,column=j,value=v); c.font=A(size=10)
            c.alignment=Alignment(vertical="top",wrap_text=(j>=3))
        st=str(r[statuscol])
        fill=WARN if "NOT CITED" in st else (NEWF if st=="NEW" else OKF if "cited" in st else None)
        if fill:
            s.cell(row=i,column=statuscol+1).fill=fill
        if str(r[1]) in ("NEW","REPLACED"): s.cell(row=i,column=2).fill=NEWF
    s.freeze_panes="A2"
    return s

sheet("Main_figures",["Label","Status in revision","Description","Citation action"],MAIN,[12,17,52,60],1)
sheet("Supplementary_figures",["Label","Status","Description","Cited in main text?","Citation action"],SUPPF,[13,13,48,15,58],3)
sheet("Supplementary_tables",["Label","Status","Description","Cited in main text?","Citation action"],SUPPT,[12,11,48,15,58],3)

s=wb.create_sheet("Summary")
rows=[("Main figures",len(MAIN),sum(1 for r in MAIN if r[1] in ("NEW","REPLACED")),"0 of 22 panels cited individually"),
      ("Supplementary figures",len(SUPPF),sum(1 for r in SUPPF if r[1]=="NEW"),f"{sum(1 for r in SUPPF if r[3]=='NOT CITED')} of {len(SUPPF)-1} existing items never cited"),
      ("Supplementary tables",len(SUPPT),sum(1 for r in SUPPT if r[1]=="NEW"),f"{sum(1 for r in SUPPT if r[3]=='NOT CITED')} of 12 existing items never cited")]
for j,h in enumerate(["Category","Items after revision","New or replaced","Citation defect in the submitted version"],start=1):
    c=s.cell(row=1,column=j,value=h); c.font=A(size=10,bold=True); c.fill=HDR
    s.column_dimensions[get_column_letter(j)].width=[26,20,18,52][j-1]
for i,r in enumerate(rows,start=2):
    for j,v in enumerate(r,start=1):
        c=s.cell(row=i,column=j,value=v); c.font=A(size=10); c.alignment=Alignment(wrap_text=(j==4),vertical="top")
        if j==4: c.fill=WARN
wb.save("/mnt/user-data/outputs/w1/01_tables/Figure_table_crosswalk_and_citation_audit.xlsx")
print("saved. main",len(MAIN),"suppF",len(SUPPF),"suppT",len(SUPPT))
print("uncited supp figures:",sum(1 for r in SUPPF if r[3]=="NOT CITED"))
print("uncited supp tables :",sum(1 for r in SUPPT if r[3]=="NOT CITED"))
