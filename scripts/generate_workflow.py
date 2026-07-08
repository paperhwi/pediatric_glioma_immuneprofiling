"""Generate a simple workflow PNG using matplotlib."""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import matplotlib

matplotlib.rcParams['font.family'] = 'DejaVu Sans'

fig, ax = plt.subplots(figsize=(11, 7.5), dpi=150)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Color palette
c_data   = "#2E5C8A"
c_qc     = "#5B8DBE"
c_decon  = "#7BAFD4"
c_clust  = "#E8A87C"
c_assoc  = "#C38D9E"
c_text_w = "white"
c_text_b = "#1A1A1A"

def box(ax, x, y, w, h, text, color, text_color="white", fontsize=10):
    rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                          facecolor=color, edgecolor="#333", linewidth=1.2)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fontsize, color=text_color, weight="bold")

def arrow(ax, x1, y1, x2, y2):
    a = FancyArrowPatch((x1, y1), (x2, y2),
                        arrowstyle="->,head_width=4,head_length=6",
                        color="#333", linewidth=1.6, mutation_scale=10)
    ax.add_patch(a)

# Row 1: Data sources
box(ax, 5, 82, 28, 12,
    "OpenPedCan v15\nhistologies.tsv\n(n=47,895 specimens)", c_data, fontsize=9)
box(ax, 38, 82, 24, 12,
    "fusion-putative-\noncogenic.tsv", c_data, fontsize=9)
box(ax, 67, 82, 28, 12,
    "gene-expression-rsem\n-tpm-collapsed.rds", c_data, fontsize=9)

# Row 2: Cohort filtering
box(ax, 15, 64, 70, 11,
    "Step 1.  Cohort filtering   (PBTA + RNA-Seq + Tumor, pediatric glioma subtypes)",
    c_qc, fontsize=10)

# Row 3: Cohort groups
box(ax, 4, 47, 17, 11,
    "DMG_K27\n(K28-altered)\nn = 225", c_qc, fontsize=8)
box(ax, 23, 47, 17, 11,
    "DHG_G34\n(G35-mutant)\nn = 37", c_qc, fontsize=8)
box(ax, 42, 47, 17, 11,
    "pHGG_WT\n(H3/IDH-WT)\nn = 192", c_qc, fontsize=8)
box(ax, 61, 47, 17, 11,
    "IHG\n(infant + fusion)\nn = 28", c_qc, fontsize=8)
box(ax, 80, 47, 17, 11,
    "BRAF_ALT\n(secondary)\nn = 384", "#9DBED0", fontsize=8)

# Row 4: Expression matrix + CIBERSORTx
box(ax, 15, 31, 70, 10,
    "Step 2.  TPM matrix extraction + HGNC mapping  →  CIBERSORTx (LM22)",
    c_decon, fontsize=10)

# Row 5: Clustering + association
box(ax, 4, 16, 44, 10,
    "Step 3. Consensus / k-means\nclustering  →  3 ecotypes",
    c_clust, "white", fontsize=10)
box(ax, 52, 16, 44, 10,
    "Step 4. Association tests\n(H3, location, age, OS)",
    c_assoc, "white", fontsize=10)

# Row 6: Output
box(ax, 25, 2, 50, 8,
    "Step 5.  Pathway enrichment (ssGSEA) + visualization",
    "#5C7A5C", fontsize=10)

# Arrows
arrow(ax, 19, 82, 35, 75)
arrow(ax, 50, 82, 50, 75)
arrow(ax, 81, 82, 65, 75)

arrow(ax, 12, 64, 12, 58)
arrow(ax, 31, 64, 31, 58)
arrow(ax, 50, 64, 50, 58)
arrow(ax, 69, 64, 69, 58)
arrow(ax, 88, 64, 88, 58)

arrow(ax, 50, 47, 50, 41)

arrow(ax, 50, 31, 26, 26)
arrow(ax, 50, 31, 74, 26)

arrow(ax, 26, 16, 50, 10)
arrow(ax, 74, 16, 50, 10)

# Title
ax.text(50, 97, "Analysis Workflow — Pediatric Glioma Immune Ecotypes (OpenPedCan v15)",
        ha="center", va="center", fontsize=13, weight="bold", color="#1A1A1A")

plt.tight_layout()
plt.savefig(r"C:\Users\scott\Downloads\Open PBTA\output\workflow.png",
            dpi=200, bbox_inches='tight', facecolor='white')
print("workflow.png written")
