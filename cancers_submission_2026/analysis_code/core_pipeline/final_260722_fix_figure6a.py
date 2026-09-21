#!/usr/bin/env python3
"""Recompute and plot Figure 6A with tie-corrected Mann-Whitney p-values."""
from __future__ import annotations

from pathlib import Path
import json
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(r"E:\Open PBTA")
OUT = ROOT / "output"
FINAL = ROOT / "FINAL MANUSCRIPT 260722"
BUNDLE = FINAL / "8. Final rerun outputs"
ECOTYPES = ["Lymphocyte-inflamed", "Myeloid-dominant", "Immune-desert"]
FILE_KEYS = {
    "Lymphocyte-inflamed": "Lymphocyte_inflamed",
    "Myeloid-dominant": "Myeloid_dominant",
    "Immune-desert": "Immune_desert",
}
COLORS = {
    "Lymphocyte-inflamed": "#3B82F6",
    "Myeloid-dominant": "#EF4444",
    "Immune-desert": "#7C8798",
}


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def load_filtered_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    ann = pd.read_csv(OUT / "ecotype_LM22_main_k3_annotated.tsv", sep="\t")
    ann = ann.set_index("Kids_First_Biospecimen_ID")
    if ann["ecotype"].value_counts().to_dict() != {
        "Myeloid-dominant": 160,
        "Lymphocyte-inflamed": 111,
        "Immune-desert": 78,
    }:
        raise RuntimeError("Ecotype counts do not match the locked 111/160/78 assignment")

    use = ["GeneSymbol", *ann.index.tolist()]
    tpm = pd.read_csv(OUT / "tpm_for_cibersortx.tsv", sep="\t", usecols=use)
    tpm = tpm.set_index("GeneSymbol")
    if not tpm.index.is_unique:
        raise RuntimeError("TPM gene index is not unique")

    expr = (tpm > 1).sum(axis=1) >= 0.05 * tpm.shape[1]
    tpm = tpm.loc[expr]
    gene = pd.Series(tpm.index.astype(str), index=tpm.index)
    keep_kinase = gene.str.match(r"^RPS6K")
    exclude = (
        gene.str.match(r"^MT-") | gene.str.match(r"^MTRNR") |
        gene.str.match(r"^MT(ATP|CO|CYB|ND|RNR)\d") |
        (gene.str.match(r"^RP[LS]\d") & ~keep_kinase) |
        gene.str.match(r"^RPLP\d") | gene.str.match(r"^RPSA$") |
        gene.str.match(r"^RNA(5S|5-8S|18S|28S|45S|5-8SN|45SN)") |
        gene.str.contains(r"rRNA") | gene.str.match(r"^RN(U|VU)\d") |
        gene.str.fullmatch(r"U\d+") |
        gene.str.fullmatch(r"7SK|Y_RNA|Vault|Metazoa_SRP|uc_338") |
        gene.str.match(r"^RN7S[KL]") |
        gene.str.match(r"^(SNORD|SNORA|SCARNA|SNAR)") |
        gene.str.match(r"^RNY\d") | gene.str.match(r"^MIR\d") |
        gene.str.match(r"^MIRLET") | gene.str.match(r"^VTRNA") |
        gene.str.fullmatch(r"RPPH1|RMRP") | gene.str.match(r"^ENSG\d")
    )
    tpm = tpm.loc[~exclude]
    if tpm.shape != (18566, 349):
        raise RuntimeError(f"Unexpected filtered TPM shape: {tpm.shape}")
    return np.log2(tpm + 1), ann


def recompute(log_tpm: pd.DataFrame, ann: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], dict]:
    tables: dict[str, pd.DataFrame] = {}
    summary = {
        "method": "two-sided tie-corrected asymptotic Mann-Whitney U with continuity correction",
        "multiple_testing": "Benjamini-Hochberg within each one-versus-rest contrast",
        "n_samples": 349,
        "n_genes": int(log_tpm.shape[0]),
        "comparisons": {},
    }
    for ecotype in ECOTYPES:
        inside = ann.index[ann["ecotype"] == ecotype]
        outside = ann.index[ann["ecotype"] != ecotype]
        a = log_tpm[inside].to_numpy()
        b = log_tpm[outside].to_numpy()
        result = stats.mannwhitneyu(
            a, b, axis=1, alternative="two-sided", method="asymptotic",
            use_continuity=True,
        )
        u = np.asarray(result.statistic, dtype=float)
        p = np.asarray(result.pvalue, dtype=float)
        if not (np.isfinite(p).all() and (p > 0).all()):
            raise RuntimeError(f"Non-finite or zero p-value in {ecotype}")
        q = bh(p)
        mu = len(inside) * len(outside) / 2
        z = np.sign(u - mu) * stats.norm.isf(np.clip(p / 2, np.finfo(float).tiny, 0.5))
        log2fc = a.mean(axis=1) - b.mean(axis=1)
        table = pd.DataFrame({
            "gene": log_tpm.index,
            "log2FC": log2fc,
            "z": z,
            "p": p,
            "q_BH": q,
            "mean_eco_log2tpm": a.mean(axis=1),
            "mean_rest_log2tpm": b.mean(axis=1),
            "n_eco": len(inside),
            "n_rest": len(outside),
        }).sort_values(["q_BH", "log2FC"], ascending=[True, False])
        up = int(((table["q_BH"] < 0.05) & (table["log2FC"] > 1)).sum())
        down = int(((table["q_BH"] < 0.05) & (table["log2FC"] < -1)).sum())
        max_nlq = float((-np.log10(table["q_BH"])).max())
        summary["comparisons"][ecotype] = {
            "n_ecotype": len(inside),
            "n_rest": len(outside),
            "up_q05_log2FC_gt1": up,
            "down_q05_log2FC_lt_minus1": down,
            "min_p": float(table["p"].min()),
            "min_q_BH": float(table["q_BH"].min()),
            "max_minus_log10_q": max_nlq,
            "zero_p": int((table["p"] == 0).sum()),
            "zero_q": int((table["q_BH"] == 0).sum()),
        }
        tables[ecotype] = table
    return tables, summary


def save_tables(tables: dict[str, pd.DataFrame], summary: dict) -> None:
    bundle_tables = BUNDLE / "01_tables"
    submission = ROOT / "for submission" / "06_supplementary_data"
    bundle_tables.mkdir(parents=True, exist_ok=True)
    submission.mkdir(parents=True, exist_ok=True)
    prefixes = {
        "Lymphocyte-inflamed": "SupplTableS3a_DEG_Lymphocyte_inflamed.tsv",
        "Myeloid-dominant": "SupplTableS3b_DEG_Myeloid_dominant.tsv",
        "Immune-desert": "SupplTableS3c_DEG_Immune_desert.tsv",
    }
    for ecotype, table in tables.items():
        key = FILE_KEYS[ecotype]
        table.to_csv(OUT / f"step8_DEG_{key}.tsv", sep="\t", index=False)
        table.to_csv(submission / prefixes[ecotype], sep="\t", index=False)
        table.to_csv(bundle_tables / f"unadjusted_DEG_tiecorrected_{key}.tsv", sep="\t", index=False)
    (bundle_tables / "Figure6A_tiecorrected_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    step8_summary = {
        "n_samples": 349,
        "n_genes_after_expr_filter": 23726,
        "n_technical_excluded": 5160,
        "n_genes_for_DEG": 18566,
        "test": summary["method"],
        "multiple_testing": summary["multiple_testing"],
        "DEG_q05_per_ecotype": {
            ec: int((tables[ec]["q_BH"] < 0.05).sum()) for ec in ECOTYPES
        },
        "up_genes_q05_log2FC_gt1": {
            ec: summary["comparisons"][ec]["up_q05_log2FC_gt1"] for ec in ECOTYPES
        },
        "down_genes_q05_log2FC_ltm1": {
            ec: summary["comparisons"][ec]["down_q05_log2FC_lt_minus1"] for ec in ECOTYPES
        },
        "max_minus_log10_q": {
            ec: summary["comparisons"][ec]["max_minus_log10_q"] for ec in ECOTYPES
        },
    }
    (OUT / "step8_summary.json").write_text(json.dumps(step8_summary, indent=2), encoding="utf-8")


def plot(tables: dict[str, pd.DataFrame], summary: dict) -> tuple[Path, Path]:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.2))
    for ax, ecotype in zip(axes, ECOTYPES):
        table = tables[ecotype].copy()
        y = -np.log10(table["q_BH"])
        sig = (table["q_BH"] < 0.05) & (table["log2FC"].abs() > 1)
        ax.scatter(table.loc[~sig, "log2FC"], y.loc[~sig], s=4, c="#D3D6DB", alpha=0.45, linewidth=0)
        ax.scatter(table.loc[sig, "log2FC"], y.loc[sig], s=6, c=COLORS[ecotype], alpha=0.72, linewidth=0)
        ax.axhline(-np.log10(0.05), color="#666666", ls="--", lw=0.8)
        ax.axvline(-1, color="#777777", ls="--", lw=0.8)
        ax.axvline(1, color="#777777", ls="--", lw=0.8)
        ax.set_xlabel("Mean log2(TPM+1) difference: ecotype vs rest")
        ax.set_ylabel("-log10(BH q)")

        comp = summary["comparisons"][ecotype]
        ymax = comp["max_minus_log10_q"]
        ax.set_ylim(-0.02 * ymax, ymax * 1.08)
        xmin, xmax = table["log2FC"].quantile([0.0005, 0.9995])
        pad = 0.11 * (xmax - xmin)
        ax.set_xlim(xmin - pad, xmax + 2.2 * pad)
        ax.set_title(
            f"{ecotype}\nup {comp['up_q05_log2FC_gt1']}  down {comp['down_q05_log2FC_lt_minus1']}"
            f"  |  max -log10(q)={ymax:.2f}",
            fontsize=10.5,
        )

        top = table[(table["q_BH"] < 0.05) & (table["log2FC"] > 0)].nlargest(8, "log2FC")
        label_y = np.linspace(ymax * 0.90, ymax * 0.38, len(top))
        text_x = ax.get_xlim()[1] - 0.04 * (ax.get_xlim()[1] - ax.get_xlim()[0])
        for (_, row), ty in zip(top.iterrows(), label_y):
            py = float(-np.log10(row["q_BH"]))
            ax.plot(row["log2FC"], py, marker="o", ms=3.2, mfc="white", mec="#30343B", mew=0.7, zorder=5)
            ax.annotate(
                row["gene"], xy=(row["log2FC"], py), xytext=(text_x, ty),
                ha="right", va="center", fontsize=7.2,
                arrowprops={"arrowstyle": "-", "color": "#70757A", "lw": 0.65},
            )
        ax.grid(axis="y", color="#E8EAED", lw=0.5, alpha=0.7)

    fig.suptitle(
        "Tie-corrected Mann-Whitney differential expression (one-versus-rest, n=349)",
        fontsize=12.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=2.2)
    figure_dir = BUNDLE / "02_figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    png = figure_dir / "Figure6A.png"
    pdf = figure_dir / "Figure6A.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    destinations = [
        FINAL / "1. Figure 300dpi" / "Main figures" / "Figure6A.png",
        OUT / "figs_step8" / "volcanos_per_ecotype.png",
        OUT / "figs_300dpi" / "Fig5A_volcanos.png",
        ROOT / "for submission" / "03_figures" / "Fig5A_volcanos.png",
    ]
    for destination in destinations:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(png, destination)
    shutil.copy2(pdf, FINAL / "1. Figure 300dpi" / "Main figures" / "Figure6A.pdf")
    return png, pdf


def main() -> None:
    log_tpm, ann = load_filtered_expression()
    tables, summary = recompute(log_tpm, ann)
    save_tables(tables, summary)
    png, pdf = plot(tables, summary)
    print(json.dumps({"status": "complete", "png": str(png), "pdf": str(pdf), **summary}, indent=2))


if __name__ == "__main__":
    main()
