"""
Step 8 — Per-ecotype differential expression + g:Profiler functional enrichment

Per ecotype (one-vs-rest), Wilcoxon rank-sum on log2(TPM+1) across genes,
log2FC computed on means, BH FDR per comparison.

Top 200 markers per ecotype → g:Profiler REST API (g:GOSt) over GO:BP, GO:MF,
GO:CC, KEGG, REAC, WP for functional enrichment.

Inputs
------
- output/tpm_for_cibersortx.tsv (gene × sample, TPM, 55408 × 702)
- output/ecotype_LM22_main_k3_annotated.tsv (Main cohort, n=349)

Outputs (output/, figs_step8/)
- step8_DEG_{ecotype}.tsv (gene, log2FC, p, q, n_eco, n_rest, mean_eco, mean_rest)
- step8_top_markers.tsv (top 200 markers per ecotype combined)
- step8_gprofiler_{ecotype}.tsv (g:Profiler enrichment table)
- figs_step8/volcano_{ecotype}.png
- figs_step8/top_enrichment_{ecotype}.png
"""
from __future__ import annotations
from pathlib import Path
import json, time, warnings
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import requests
warnings.filterwarnings("ignore")

BASE = Path("/sessions/zealous-vibrant-rubin/mnt/Open PBTA")
OUT  = BASE / "output"; FIG = OUT / "figs_step8"; FIG.mkdir(exist_ok=True)

# ----------------------------------------------------------------------------
# 1) Load TPM and ecotype
# ----------------------------------------------------------------------------
print("=== [1] Load TPM (this takes ~30s) ===")
t0 = time.time()
tpm = pd.read_csv(OUT / "tpm_for_cibersortx.tsv", sep="\t", index_col=0)
print(f"  raw TPM: {tpm.shape}  ({time.time()-t0:.1f}s)")

ec3 = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t").set_index("Kids_First_Biospecimen_ID")
common = list(set(tpm.columns) & set(ec3.index))
print(f"  Main cohort with ecotype: {len(common)}")
tpm = tpm[common]; ec3 = ec3.loc[common]

# ----------------------------------------------------------------------------
# 2) Gene filtering: expressed in ≥5% samples at TPM>1
# ----------------------------------------------------------------------------
print("\n=== [2] Filter genes ===")
expr_mask = (tpm > 1).sum(axis=1) >= (0.05 * tpm.shape[1])
tpm = tpm[expr_mask].copy()
print(f"  after filter: {tpm.shape[0]} genes retained")
log_tpm = np.log2(tpm + 1.0)
print(f"  log2(TPM+1) ready ({time.time()-t0:.1f}s elapsed)")

# ----------------------------------------------------------------------------
# 3) Wilcoxon one-vs-rest per ecotype (vectorized)
# ----------------------------------------------------------------------------
print("\n=== [3] Wilcoxon DE (one-vs-rest per ecotype) ===")
def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); order = np.argsort(p)
    ranked = p[order]; q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1); return out

ecotypes = sorted(ec3["ecotype"].dropna().unique().tolist())
deg_tables = {}
for ec in ecotypes:
    in_ec  = ec3.index[ec3["ecotype"] == ec].tolist()
    out_ec = ec3.index[ec3["ecotype"] != ec].tolist()
    print(f"  {ec}: n_in={len(in_ec)}, n_rest={len(out_ec)}")
    A = log_tpm[in_ec].values
    B = log_tpm[out_ec].values
    # Mann-Whitney U via scipy.stats.ranksums (z-based normal approximation, ties handled poorly,
    # but on continuous-ish log-TPM is fine and fast). Vectorize over genes.
    t1 = time.time()
    n1, n2 = A.shape[1], B.shape[1]
    n = n1 + n2
    # Vectorized ranks across rows
    combined = np.concatenate([A, B], axis=1)
    ranks = np.apply_along_axis(stats.rankdata, 1, combined)
    R1 = ranks[:, :n1].sum(axis=1)
    U1 = R1 - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    sigma = np.sqrt(n1 * n2 * (n + 1) / 12)
    z = (U1 - mu) / sigma
    p = 2 * (1 - stats.norm.cdf(np.abs(z)))
    log2fc = A.mean(axis=1) - B.mean(axis=1)
    q = bh(p)
    deg = pd.DataFrame({
        "gene": log_tpm.index,
        "log2FC": log2fc,
        "z": z,
        "p": p,
        "q_BH": q,
        "mean_eco_log2tpm": A.mean(axis=1),
        "mean_rest_log2tpm": B.mean(axis=1),
        "n_eco": n1,
        "n_rest": n2,
    }).sort_values("q_BH")
    deg.to_csv(OUT / f"step8_DEG_{ec.replace(' ', '_').replace('-', '_')}.tsv",
               sep="\t", index=False)
    deg_tables[ec] = deg
    print(f"    Wilcoxon done ({time.time()-t1:.1f}s) — "
          f"q<0.05: {(deg['q_BH']<0.05).sum()}; |log2FC|>1 & q<0.05: "
          f"{((deg['q_BH']<0.05)&(deg['log2FC'].abs()>1)).sum()}")

# ----------------------------------------------------------------------------
# 4) Top 200 up-regulated markers per ecotype (q<0.05, log2FC>0, sorted by log2FC)
# ----------------------------------------------------------------------------
print("\n=== [4] Top markers per ecotype ===")
top_markers = {}
combined_rows = []
for ec, deg in deg_tables.items():
    up = deg[(deg["q_BH"] < 0.05) & (deg["log2FC"] > 0)].sort_values("log2FC", ascending=False)
    markers = up.head(200)["gene"].tolist()
    top_markers[ec] = markers
    print(f"  {ec}: {len(markers)} up-markers (top 200 of {len(up)})")
    for g in markers:
        row = up[up["gene"] == g].iloc[0]
        combined_rows.append({"ecotype": ec, "gene": g, "log2FC": float(row["log2FC"]),
                              "q_BH": float(row["q_BH"])})
pd.DataFrame(combined_rows).to_csv(OUT / "step8_top_markers.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------
# 5) Volcano plots per ecotype
# ----------------------------------------------------------------------------
print("\n=== [5] Volcano plots ===")
fig, axes = plt.subplots(1, len(ecotypes), figsize=(5.2*len(ecotypes), 4.5))
if len(ecotypes) == 1: axes = [axes]
ec_colors = {"Lymphocyte-inflamed":"#3B82F6","Myeloid-dominant":"#EF4444","Immune-desert":"#9CA3AF"}
for ax, ec in zip(axes, ecotypes):
    d = deg_tables[ec].copy()
    d["nl_q"] = -np.log10(np.clip(d["q_BH"], 1e-300, 1))
    sig = (d["q_BH"] < 0.05) & (d["log2FC"].abs() > 1)
    ax.scatter(d.loc[~sig, "log2FC"], d.loc[~sig, "nl_q"], s=3, c="lightgray", alpha=.5, edgecolor="none")
    ax.scatter(d.loc[sig, "log2FC"], d.loc[sig, "nl_q"], s=4, c=ec_colors.get(ec, "#666"), alpha=.7, edgecolor="none")
    # label top 8 up markers
    top = d[(d["q_BH"]<0.05)&(d["log2FC"]>0)].sort_values("log2FC", ascending=False).head(8)
    for _, r in top.iterrows():
        ax.text(r["log2FC"], r["nl_q"], r["gene"], fontsize=6.5, ha="left")
    ax.axhline(-np.log10(0.05), color="black", ls="--", alpha=.3)
    ax.axvline(1, color="black", ls="--", alpha=.3)
    ax.axvline(-1, color="black", ls="--", alpha=.3)
    ax.set_xlabel("log2 fold change (ecotype vs rest)")
    ax.set_ylabel("-log10(q BH)")
    n_up = ((d["q_BH"]<0.05)&(d["log2FC"]>1)).sum()
    n_dn = ((d["q_BH"]<0.05)&(d["log2FC"]<-1)).sum()
    ax.set_title(f"{ec}\n↑{n_up}  ↓{n_dn} sig genes (|FC|>2, q<0.05)", fontsize=10)
plt.tight_layout()
fig.savefig(FIG / "volcanos_per_ecotype.png", dpi=200); plt.close()

# ----------------------------------------------------------------------------
# 6) g:Profiler enrichment per ecotype
# ----------------------------------------------------------------------------
print("\n=== [6] g:Profiler enrichment ===")
GP_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
all_enrich = []
for ec, markers in top_markers.items():
    if len(markers) < 10:
        print(f"  {ec}: too few markers, skip"); continue
    payload = {
        "organism": "hsapiens",
        "query": markers,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "ordered": True,
        "no_iea": False,
        "domain_scope": "annotated",
        "no_evidences": True,
    }
    try:
        r = requests.post(GP_URL, json=payload, timeout=30)
        r.raise_for_status()
        results = r.json().get("result", [])
        print(f"  {ec}: {len(results)} terms returned")
        rows = []
        for term in results[:300]:
            rows.append({
                "ecotype": ec,
                "source": term.get("source"),
                "term_id": term.get("native"),
                "term_name": term.get("name"),
                "p_value": term.get("p_value"),
                "term_size": term.get("term_size"),
                "query_size": term.get("query_size"),
                "intersection_size": term.get("intersection_size"),
                "precision": term.get("precision"),
                "recall": term.get("recall"),
            })
        edf = pd.DataFrame(rows)
        edf.to_csv(OUT / f"step8_gprofiler_{ec.replace(' ', '_').replace('-', '_')}.tsv",
                   sep="\t", index=False)
        all_enrich.append(edf)
    except Exception as e:
        print(f"  {ec}: g:Profiler error — {e}")
if all_enrich:
    pd.concat(all_enrich, ignore_index=True).to_csv(OUT / "step8_gprofiler_all.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------
# 7) Top 15 enrichment bar plot per ecotype
# ----------------------------------------------------------------------------
print("\n=== [7] Enrichment bar plots ===")
fig, axes = plt.subplots(1, len(ecotypes), figsize=(7.5*len(ecotypes), 5))
if len(ecotypes) == 1: axes = [axes]
for ax, ec in zip(axes, ecotypes):
    f = OUT / f"step8_gprofiler_{ec.replace(' ', '_').replace('-', '_')}.tsv"
    if not f.exists():
        ax.set_visible(False); continue
    e = pd.read_csv(f, sep="\t")
    # Keep biological-relevance sources, prefer GO:BP + REAC
    e = e[e["source"].isin(["GO:BP", "REAC", "KEGG", "WP"])]
    e = e.sort_values("p_value").head(15)
    if e.empty:
        ax.set_visible(False); continue
    nl = -np.log10(np.clip(e["p_value"], 1e-300, 1))
    labels = (e["term_name"].str.slice(0, 56) + " [" + e["source"] + "]").values
    ax.barh(labels[::-1], nl[::-1], color=ec_colors.get(ec, "#666"))
    ax.set_xlabel("-log10(adjusted p, g:SCS)")
    ax.set_title(f"{ec}  (top 15 BP/REAC/KEGG/WP)", fontsize=10)
    ax.tick_params(axis="y", labelsize=8)
plt.tight_layout()
fig.savefig(FIG / "top_enrichment_per_ecotype.png", dpi=200); plt.close()

# ----------------------------------------------------------------------------
# 8) Summary
# ----------------------------------------------------------------------------
summary = {
    "n_samples": int(len(common)),
    "n_genes_after_filter": int(tpm.shape[0]),
    "DEG_q05_per_ecotype": {ec: int((deg_tables[ec]["q_BH"]<0.05).sum()) for ec in ecotypes},
    "up_genes_q05_log2FC_gt1": {ec: int(((deg_tables[ec]["q_BH"]<0.05)&(deg_tables[ec]["log2FC"]>1)).sum()) for ec in ecotypes},
    "top_marker_counts": {ec: len(top_markers[ec]) for ec in ecotypes},
}
(OUT / "step8_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
print(f"\nTotal elapsed: {time.time()-t0:.1f}s")
