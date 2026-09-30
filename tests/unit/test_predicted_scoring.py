"""Generic model harmonization and scoring (predicted expression, splicing, protein, methylation)."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from efgpp.data.genotype.formats import resolve_fileset
from efgpp.data.predicted.harmonize import BuildMismatchError, apply_coverage, harmonize
from efgpp.data.predicted.models import (
    WEIGHT_SCHEMA,
    predictdb_tissues,
    predictdb_weights,
    standardize,
)
from efgpp.data.predicted.scorer import bed_sites, score_command, score_python

from ._genotypes import write_bed, write_predictdb

SAMPLES = ["P001", "P002", "P003"]


def _geno(tmp_path: Path) -> object:
    # bim: (chrom, id, pos, A1, A2); PLINK counts A1. A1 counts per participant.
    return resolve_fileset(write_bed(tmp_path / "g", SAMPLES, [
        ("1", "rs1", 100, "G", "A"),   # REF A, ALT G
        ("1", "rs2", 200, "T", "C"),   # REF C, ALT T
        ("1", "rs3", 300, "T", "A"),   # A/T: strand-ambiguous
        ("1", "rs4", 400, "C", "T"),   # REF T, ALT C
    ], [[0, 1, 2], [2, 1, 0], [1, 1, 1], [0, 2, 1]]), "bed")


def _weights(rows: list[tuple[str, str, int | None, str, str, float]], **const: object) -> pl.DataFrame:
    df = pl.DataFrame(rows, schema=["feature_id", "chromosome", "position", "effect_allele", "other_allele", "weight"],
                      orient="row")
    return standardize(df, **{"model_id": "m1", "genome_build": "GRCh38", **const})


def _scores(fs: object, weights: pl.DataFrame, **kw: object) -> tuple[pl.DataFrame, pl.DataFrame]:
    h = harmonize(weights, bed_sites(fs), genotype_build="GRCh38", model_build="GRCh38", **kw)  # type: ignore[arg-type]
    return score_python(fs, h.matched).sort("IID"), h.qc  # type: ignore[arg-type]


def test_expression_model_exact_values(tmp_path: Path) -> None:
    """GeneA: rs1 weight 0.5 (effect G), rs2 weight -0.2 (effect T)."""
    fs = _geno(tmp_path)
    s, qc = _scores(fs, _weights([("GeneA", "1", 100, "G", "A", 0.5), ("GeneA", "1", 200, "T", "C", -0.2)]))
    # P001: 0*0.5 + 2*-0.2 = -0.4; P002: 0.5 - 0.2 = 0.3; P003: 1.0 + 0 = 1.0
    assert [round(x, 10) for x in s.get_column("GeneA").to_list()] == [-0.4, 0.3, 1.0]
    assert qc.select("n_exact_matches", "coverage_fraction").row(0) == (2, 1.0)


def test_generic_scorer_swap_flip_ambiguous_missing(tmp_path: Path) -> None:
    fs = _geno(tmp_path)
    base = [("ProteinA", "1", 100, "G", "A", 0.25), ("ProteinA", "1", 200, "T", "C", -0.5)]
    s, qc = _scores(fs, _weights(base))
    assert [round(x, 10) for x in s.get_column("ProteinA").to_list()] == [-1.0, -0.25, 0.5]
    # allele swap: effect = REF A (count of A = 2 - ALT count)
    s, qc = _scores(fs, _weights([("P", "1", 100, "A", "G", 0.25)]))
    assert s.get_column("P").to_list() == [0.5, 0.25, 0.0] and qc.get_column("n_allele_swaps").item() == 1
    # strand flip: model on the minus strand (C/T for REF A / ALT G -> complement of G is C)
    s, qc = _scores(fs, _weights([("P", "1", 100, "C", "T", 0.25)]))
    assert s.get_column("P").to_list() == [0.0, 0.25, 0.5] and qc.get_column("n_strand_flips").item() == 1
    # A/T SNP: excluded (not strand-resolvable) unless alleles are known forward-strand
    _s, qc = _scores(fs, _weights([("P", "1", 300, "T", "A", 1.0), ("P", "1", 100, "G", "A", 0.25)]))
    assert qc.select("n_ambiguous_excluded", "n_matched_variants").row(0) == (1, 1)
    s, qc = _scores(fs, _weights([("P", "1", 300, "T", "A", 1.0)]), forward_strand=True)
    assert s.get_column("P").to_list() == [1.0, 1.0, 1.0] and qc.get_column("n_exact_matches").item() == 1
    # missing SNP and allele mismatch
    _s, qc = _scores(fs, _weights([("P", "1", 100, "G", "A", 0.25), ("P", "1", 999, "G", "A", 1.0),
                                   ("P", "1", 400, "A", "C", 1.0)]))  # T/C site: A/C matches no strand
    assert qc.select("n_missing_variants", "n_allele_mismatch", "n_matched_variants").row(0) == (1, 1, 1)


def test_low_coverage_is_flagged_not_predicted(tmp_path: Path) -> None:
    fs = _geno(tmp_path)
    rows = [("F", "1", 100, "G", "A", 0.1), ("F", "1", 200, "T", "C", 0.1), ("F", "1", 400, "C", "T", 0.1),
            ("F", "1", 300, "T", "A", 0.1)] + [("F", "1", 1000 + i, "G", "A", 0.1) for i in range(6)]
    h = harmonize(_weights(rows), bed_sites(fs), genotype_build="GRCh38", model_build="GRCh38", forward_strand=True)
    qc = apply_coverage(h.qc, 0.50, allow_low=False)
    row = qc.row(0, named=True)
    assert (row["n_model_variants"], row["n_matched_variants"], row["coverage_fraction"]) == (10, 4, 0.4)
    assert (row["status"], row["predict"]) == ("LOW_COVERAGE", False)
    assert apply_coverage(h.qc, 0.50, allow_low=True).get_column("predict").item()


def test_build_mismatch_never_silent(tmp_path: Path) -> None:
    fs = _geno(tmp_path)
    w = _weights([("F", "1", 100, "G", "A", 0.1)])
    with pytest.raises(BuildMismatchError, match="GRCh37 != genotype build GRCh38|model build GRCh38 != genotype build GRCh37"):
        harmonize(w, bed_sites(fs), genotype_build="GRCh37", model_build="GRCh38")
    with pytest.raises(BuildMismatchError, match="unknown"):
        harmonize(w, bed_sites(fs), genotype_build="GRCh38", model_build=None)


def test_predictdb_expression_and_splicing_models(tmp_path: Path) -> None:
    """Tiny PredictDB SQLite databases: a gene model and an intron (sQTL) model, scored exactly."""
    fs = _geno(tmp_path)
    eqtl = write_predictdb(tmp_path / "eqtl" / "mashr_Whole_Blood.db", [
        ("ENSG0001.1", "rs1", "chr1_100_A_G_b38", "A", "G", 0.5), ("ENSG0001.1", "rs2", "chr1_200_C_T_b38", "C", "T", -0.2)])
    sqtl = write_predictdb(tmp_path / "sqtl" / "mashr_Whole_Blood.db", [
        ("intron_1_150_250", "rs2", "chr1_200_C_T_b38", "C", "T", 0.3),
        ("intron_1_150_250", "rs3", "chr1_300_A_T_b38", "A", "T", 0.1)])
    assert list(predictdb_tissues(tmp_path / "eqtl")) == ["Whole_Blood"]
    w, feats = predictdb_weights(eqtl, modality="expression", tissue="Whole_Blood", dataset="gtex_v8_mashr_eqtl",
                                 genome_build="GRCh38")
    assert set(w.columns) == set(WEIGHT_SCHEMA) and feats.get_column("n_variants").item() == 2
    h = harmonize(w, bed_sites(fs), genotype_build="GRCh38", model_build="GRCh38", forward_strand=True)
    got = score_python(fs, h.matched).sort("IID").get_column("ENSG0001.1").to_list()
    assert [round(x, 10) for x in got] == [-0.4, 0.3, 1.0]
    w2, _ = predictdb_weights(sqtl, modality="splicing", tissue="Whole_Blood", dataset="gtex_v8_mashr_sqtl",
                              genome_build="GRCh38")
    h2 = harmonize(w2, bed_sites(fs), genotype_build="GRCh38", model_build="GRCh38", forward_strand=True)
    got2 = score_python(fs, h2.matched).sort("IID").get_column("intron_1_150_250").to_list()
    # rs2 T count x 0.3 + rs3 T count x 0.1 (A/T kept: PredictDB ids are forward strand)
    assert [round(x, 10) for x in got2] == [0.7, 0.4, 0.1]


def test_methylation_weighted_sum(tmp_path: Path) -> None:
    fs = _geno(tmp_path)
    w = _weights([("cg0001", "1", 100, "G", "A", 0.02), ("cg0001", "1", 400, "C", "T", -0.01)],
                 model_id="mimosa:whole_blood_v2")
    s, _ = _scores(fs, w)
    # P001: 0 + 0*-0.01 ; P002: 0.02 + 2*-0.01 ; P003: 0.04 + 1*-0.01
    assert [round(x, 10) for x in s.get_column("cg0001").to_list()] == [0.0, 0.0, 0.03]


def test_plink2_score_command(tmp_path: Path) -> None:
    fs = resolve_fileset(tmp_path / "keyed", "pgen")
    cmd = score_command(fs, tmp_path / "s.tsv", tmp_path / "e.txt", 3, tmp_path / "out", 2)
    assert cmd[cmd.index("--score") + 1:cmd.index("--score") + 7] == [str(tmp_path / "s.tsv"), "1", "2", "header-read",
                                                                     "no-mean-imputation", "cols=+scoresums"]
    assert cmd[cmd.index("--score-col-nums") + 1] == "3-5" and "--extract" in cmd
