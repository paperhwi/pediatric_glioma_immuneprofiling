import nbformat as nbf
from nbclient import NotebookClient
src=lambda f: open(f).read()
def build(fn,title,intro,steps):
    nb=nbf.v4.new_notebook(); C=[nbf.v4.new_markdown_cell(f"# {title}\n\n{intro}")]
    for k,b in steps: C.append(nbf.v4.new_markdown_cell(b) if k=="md" else nbf.v4.new_code_cell(b))
    nb["cells"]=C
    nb.metadata={"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"},"language_info":{"name":"python"}}
    NotebookClient(nb,timeout=2400,kernel_name="python3",resources={"metadata":{"path":"/tmp/claude-0/w2"}}).execute()
    nbf.write(nb,f"/mnt/user-data/outputs/{fn}"); print("built",fn)

build("W2_single_cell_support.ipynb",
 "W2 — Single-cell support for the cellular interpretation  [Reviewer 1 item 5]",
 """Reviewer 1: *\"Please add single-cell spatial or tissue level immune ecotypes validation if possible.\"*

**What this is not.** The obvious design — pseudobulk each public tumour, assign an ecotype, and
compare with the single-cell composition — is not valid here. Every public paediatric high-grade
glioma single-cell dataset large enough to use is CD45-sorted or otherwise enriched, so the
cell-type proportions are a function of the sorting gate rather than of the tissue, and with a
handful of tumours the comparison would also be powerless. We do not make that claim.

**What is testable.** Two questions that do not depend on composition being unbiased, both of them
addressing the manuscript's own second limitation — that deconvolution and signature scores are
computational estimates that cannot establish cellular identity without single-cell data:

1. Do the 24 brain-tuned signatures score in the cell types they are named for?
2. Are the genes that define the ecotype axis expressed by immune cells or by malignant cells?

Data: GSE227983, immune and tumour cell landscape in paediatric high-grade glioma, with the
original authors' cell annotations.""",
 [("md","## 1. Signature attribution across annotated immune cell types (10x, 18,619 cells)"),
  ("code",src("w2b_attribution.py")),
  ("md","""15 of the 16 directional signatures score highest in the cell type they are named for.
The exception is `Glioma_Inflammatory_Wang2017`.

## 2. Effect size: AUROC per signature, and per sample"""),
  ("code",src("w2b2_auroc.py")),
  ("md","""Median AUROC 0.868. Twelve signatures exceed 0.85. Three are weak or wrong:
`M1_Macrophage` (0.668), `Dendritic_Cell_Activation` (0.640) and `Glioma_Inflammatory_Wang2017`
(0.075 — anti-correlated with myeloid identity, i.e. it identifies T cells).

## 3. Why that signature fails, and whether it matters

`Glioma_Inflammatory_Wang2017` contains GZMB, PRF1, IFNG, GZMA, NKG7, CD8A, CD8B and CCL5 — a
cytotoxic lymphocyte programme — yet it is grouped in the **myeloid/microglia** theme in the
published analysis. `Dendritic_Cell_Activation`, a myeloid programme, is grouped in the
**T-cell axis**. Both assignments are wrong, in opposite directions. The question is whether
correcting them changes anything."""),
  ("code",src("w2_theme_audit.py")),
  ("md","""Correcting both assignments moves the myeloid theme effect size from eps-sq 0.755 to 0.732 and the
T-cell axis from 0.663 to 0.678. The ranking of the five themes is unchanged, the ecotype means
are unchanged to two decimal places, and the myeloid-lymphoid theme correlation moves from 0.856
to 0.842. The assignments should be corrected for accuracy, but no conclusion depends on them.

## 4. Cellular source of the ecotype-defining genes (Smart-seq2, same platform for both)"""),
  ("code",src("w2c_source.py")),
  ("md","""242 of 246 (98.4%) of the genes defining the Lymphocyte-inflamed versus Immune-desert axis are
expressed most highly by immune cells — 195 by myeloid cells, 47 by T cells — and only 4 by
malignant cells. Against expression-matched random genes drawn from the same dataset (72.6%,
95% range 67.5-77.2%), permutation P = 0.005.

The ecotype axis is therefore immune-derived rather than a tumour-intrinsic programme carrying an
immune label, and its dominant cellular source is myeloid, which is consistent with the
myeloid/microglia theme showing the largest effect in the bulk cohort.

## 5. Figure"""),
  ("code",src("w2f_fig_scrna.py")),
  ("md","""## What this does and does not establish

It supports the **interpretation** of the signatures and of the ecotype axis at cellular
resolution. It does **not** validate the ecotypes themselves: four tumours cannot establish that a
three-state representation generalises, and the sorted design means nothing here speaks to
prevalence. That distinction is stated in the response letter and in the manuscript.""")])

build("W2_spatial_support.ipynb",
 "W2 — Spatial organisation of the immune programmes  [Reviewer 1 item 5]",
 """Visium spots are not sorted, so unlike the single-cell data above, composition is meaningful
here. Three H3.3 K27M paediatric diffuse midline glioma sections (GSE268577).

Two questions: do the immune programmes form spatially coherent domains rather than noise, and
does the spot-level picture match what the bulk ecotypes imply?""",
 [("md","## 1. Spot-level signature scores and Moran's I"),
  ("code",src("w2d_spatial.py")),
  ("md","""Ten of the 24 signatures are spatially autocorrelated at q < 0.05 in all three sections. The
strongest are myeloid and antigen-presentation programmes (MDM_Klemm2020 mean I = 0.455,
MHC_Class_II 0.410, MHC_Class_I 0.397, MoTAM_Antunes2021 0.395). Sample-2 is consistently the
weakest section across every analysis here.

## 2. Spot-level ecotype assignment and spatial coherence"""),
  ("code",src("w2e_spatial_ecotype.py")),
  ("md","""Two results, and the second is the more interesting one.

**The assignments are spatially coherent.** Neighbouring spots share an ecotype assignment more
often than chance in all three sections (0.540, 0.357 and 0.562 against nulls of 0.358, 0.338 and
0.370; P = 0.001, 0.019, 0.001). The programmes are organised in tissue, not scattered.

**Every section contains all three ecotypes.** A single tumour is not one ecotype; the bulk label
is a spatial average over regions. This is a limitation of the bulk framework that the spatial
data make concrete, and it should be stated rather than buried.

**The myeloid-lymphoid coupling is an aggregation effect, in part.** In bulk the two theme
composites correlate at r = 0.80; per spot the correlation is 0.36, 0.05 and 0.48. Bulk mixing
inflates the apparent coupling of myeloid and lymphoid programmes.

## 3. Figure"""),
  ("code",src("w2g_fig_spatial.py")),
  ("md","""## Framing

Three sections from three tumours. This shows that the immune programmes are spatially structured
and that ecotype assignment varies within a tumour. It is not external validation and does not
speak to the survival association.""")])
