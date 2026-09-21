CANCERS SUBMISSION PACKAGE — 21 September 2026
Integrated transcriptomic immune profiling identifies survival-associated immune ecotypes
in pediatric diffuse high-grade glioma

This package contains the display items only. The manuscript text (Abstract, Introduction,
Results, Discussion, Methods) is written separately by the authors.

--------------------------------------------------------------------------------
CONTENTS
--------------------------------------------------------------------------------
1. Main figures/            Figure1-Figure7, each as one assembled file
                            600 dpi PNG + vector PDF (Type-42 embedded fonts)

2. Supplementary figures/   FigureS1-FigureS15, same two formats

3. Legends/                 Figure_and_Table_Legends_Cancers.docx
                            All main and supplementary figure legends, the five main
                            table legends, the supplementary table legends, the
                            abbreviation list.

4. Supplementary material/  Supplementary_Material_Cancers.docx
                              one file: contents, all 15 supplementary figures with
                              their legends, and the supplementary table legends
                            Supplementary_Tables_Cancers.xlsx
                              README     index of all 25 tables
                              S1 - S25   one worksheet per table

5. Scripts/                 Every script that generates or assembles the figures and
                            tables in this package.

6. Figure separate files/   Figure1A-Figure7B, 29 main-figure panels supplied separately
                            as 600 dpi PNG and vector PDF; MANIFEST.tsv records dimensions
                            and SHA-256 hashes.

--------------------------------------------------------------------------------
WHAT CHANGED RELATIVE TO THE PREVIOUS PACKAGE
--------------------------------------------------------------------------------
Main figures
  Figure 1   workflow redrawn: k = 2-10 evaluation, bootstrap consensus clustering, and
             "k = 2 most stable; k = 3 retained as an exploratory three-state
             representation"; the analysis-set derivation is now panel B, redrawn in the
             same style and with US spelling.
  Figure 2   LM22 limitation and the core cross-method benchmarking only. The
             ssGSEA-versus-LM22 macrophage/monocyte heatmap is Supplementary Figure S1.
  Figure 3   k-selection in five panels (PAC, silhouette, relative dAUC, gap statistic,
             minimum cluster size) plus PCA/UMAP and the defining-feature heatmap. The
             cohort-composition panel is dropped because Figure 4B carries it.
  Figure 4   unchanged.
  Figure 5   new: single-cell AUROC, cell-type source of the ecotype-axis genes, the
             three spatial sections, spatial coherence and Moran's I. The previous bulk
             immune-program panels are Supplementary Figures S6 and S7.
  Figure 6   molecular-group-adjusted volcanoes plus the overlap and gradient analysis;
             the Venn diagram is replaced by proportional overlap bars. Pathway
             enrichment is Supplementary Figure S8.
  Figure 7   overall Kaplan-Meier with numbers at risk and censoring marks, plus the
             multivariable Cox model. Subtype-stratified survival is Supplementary
             Figure S10.

Supplementary figures
  Renumbered from S1 with no gaps.

Supplementary tables
  Renumbered from S1 with no gaps, 25 tables. Audit copies, superseded pre-corrected
  versions and tables that duplicate Table 1 or a main figure were removed; the tables
  needed to reproduce the figures were added where they were previously only in a
  separate data workbook.

--------------------------------------------------------------------------------
TWO NUMBERS TO NOTE
--------------------------------------------------------------------------------
  * Figure 6F reports the paired Wilcoxon signed-rank test as two-sided,
    P = 1.5 × 10⁻⁷⁰. All other statistics are unchanged.
  * Figure 3 states that the k = 2-6 values reproduce the locked pipeline output to
    within 5 x 10^-5, and describes this as a separate 1,000-replicate rerun on the same
    cohort, that is, a resampling stability analysis, not an independent replication.

--------------------------------------------------------------------------------
STILL REQUIRED FROM THE AUTHORS
--------------------------------------------------------------------------------
  1. A main-text citation for every display item. This was the Editor's point in the previous
     round and it is not fixed by renumbering alone.
  2. A Zenodo DOI for the public repository, if required by the journal.
  3. Declarations: author contributions, acknowledgements, funding, conflict of interest.
  4. Study duration: the OpenPedCan release date and the analysis period.
  5. English-editing statement, if a service was used.
