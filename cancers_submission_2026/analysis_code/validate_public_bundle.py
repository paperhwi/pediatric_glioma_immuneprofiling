#!/usr/bin/env python3
"""Validate the curated public reproducibility bundle using only the Python standard library."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reproducibility_data"
EXPECTED_COUNTS = {
    "Lymphocyte-inflamed": 111,
    "Myeloid-dominant": 160,
    "Immune-desert": 78,
}
FORBIDDEN_SUFFIXES = {".rds", ".rdata", ".h5", ".h5ad"}
FORBIDDEN_NAME_PARTS = {"tpm_for_cibersortx", "cohort_main_final", "histologies"}
FORBIDDEN_SAMPLE_COLUMNS = {
    "age",
    "age_at_diagnosis",
    "diagnosis",
    "gender",
    "sex",
    "os_days",
    "os_status",
    "overall_survival",
    "vital_status",
}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AssertionError(f"Missing header: {path}")
        return reader.fieldnames, list(reader)


def assert_sample_file(path: Path, id_column: str, expected_columns: int) -> list[dict[str, str]]:
    header, rows = read_tsv(path)
    assert len(header) == expected_columns, (path, len(header), expected_columns)
    assert len(rows) == 349, (path, len(rows))
    ids = [row[id_column] for row in rows]
    assert len(set(ids)) == 349, f"Duplicate identifiers in {path}"
    forbidden = {column.lower() for column in header} & FORBIDDEN_SAMPLE_COLUMNS
    assert not forbidden, f"Restricted columns in {path}: {sorted(forbidden)}"
    return rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        lower_name = path.name.lower()
        assert path.suffix.lower() not in FORBIDDEN_SUFFIXES, f"Restricted file type: {path}"
        assert not any(part in lower_name for part in FORBIDDEN_NAME_PARTS), f"Restricted file: {path}"

    assignments = assert_sample_file(DATA / "ecotype_assignments.tsv", "Kids_First_Biospecimen_ID", 3)
    assignment_by_id = {row["Kids_First_Biospecimen_ID"]: row["ecotype"] for row in assignments}
    assert Counter(assignment_by_id.values()) == Counter(EXPECTED_COUNTS)

    main_matrix = assert_sample_file(
        DATA / "feature_matrices" / "main_cohort_LM22_ssGSEA_z.tsv", "sample", 47
    )
    tpm_matrix = assert_sample_file(
        DATA / "feature_matrices" / "TPM_harmonization_quanTIseq_ssGSEA_z.tsv",
        "Kids_First_Biospecimen_ID",
        30,
    )
    assert {row["sample"] for row in main_matrix} == set(assignment_by_id)
    assert {row["Kids_First_Biospecimen_ID"] for row in tpm_matrix} == set(assignment_by_id)

    coordinate_sets = (
        ("pca_coordinates.tsv", ("PC1", "PC2")),
        ("umap_coordinates.tsv", ("UMAP1", "UMAP2")),
    )
    for name, axes in coordinate_sets:
        rows = assert_sample_file(DATA / "embeddings" / name, "sample", 4)
        assert {row["sample"] for row in rows} == set(assignment_by_id)
        for row in rows:
            assert row["ecotype"] == assignment_by_id[row["sample"]]
            for axis in axes:
                float(row[axis])

    for name in ("revA_summary.json", "revD_summary.json"):
        with (DATA / "survival_robustness" / name).open(encoding="utf-8") as handle:
            json.load(handle)

    manifest_path = DATA / "SHA256SUMS.tsv"
    header, manifest = read_tsv(manifest_path)
    assert header == ["sha256", "path"]
    expected_files = sorted(path for path in DATA.rglob("*") if path.is_file() and path != manifest_path)
    relative_paths = [path.relative_to(DATA).as_posix() for path in expected_files]
    assert [row["path"] for row in manifest] == relative_paths
    for row, path in zip(manifest, expected_files, strict=True):
        assert row["sha256"] == sha256(path), f"Checksum mismatch: {path}"

    print("Public bundle validation passed")
    print(f"  cohort: {len(assignments)}")
    print(f"  ecotypes: {dict(Counter(assignment_by_id.values()))}")
    print("  main feature matrix: 349 x 46 features")
    print("  TPM sensitivity matrix: 349 x 29 features")
    print(f"  checksummed files: {len(expected_files)}")


if __name__ == "__main__":
    main()
