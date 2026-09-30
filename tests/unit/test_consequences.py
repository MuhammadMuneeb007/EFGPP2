"""Participant consequence counts and gene burden from carriers + annotations."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from efgpp.config.data import ParticipantVariantsConfig
from efgpp.data.genotype.build import FastaIndex
from efgpp.data.genotype.carriers import extract_carriers_bed
from efgpp.data.genotype.consequences import aggregate, most_severe, pick_vep
from efgpp.data.genotype.formats import read_samples, resolve_fileset

from ._genotypes import write_bed, write_fasta


def _vep_row(key: str, gene: str, tx: str, cons: str, *, pick: str | None = None, canonical: str | None = None,
             symbol: str | None = None) -> dict[str, str | None]:
    return {"Uploaded_variation": key, "Gene": gene, "Feature": tx, "Feature_type": "Transcript",
            "Consequence": cons, "IMPACT": "MODERATE", "SYMBOL": symbol or gene, "PICK": pick, "CANONICAL": canonical,
            "HGVSc": None, "HGVSp": None}


def _cohort(tmp_path: Path) -> tuple[Path, list[str]]:
    """P001: var1 het, var2 hom; P002: var1 hom; P003: var3 het (the spec's demonstration cohort)."""
    fasta = FastaIndex(write_fasta(tmp_path / "ref.fa", {"1": {100: "A", 200: "C", 300: "G"}}))
    fs = resolve_fileset(write_bed(tmp_path / "g", ["P001", "P002", "P003"],
                                   [("1", "var1", 100, "G", "A"), ("1", "var2", 200, "T", "C"),
                                    ("1", "var3", 300, "A", "G")],
                                   [[1, 2, 0], [2, 0, 0], [0, 0, 1]]), "bed")
    samples = read_samples(fs).with_columns(pl.col("IID").alias("native_id"), pl.col("IID").alias("participant_id"))
    extract_carriers_bed(fs, samples, fasta, tmp_path / "carriers", source_id="GENO001", genome_build="GRCh38")
    return tmp_path / "carriers" / "participant_variants", ["P001", "P002", "P003"]


def test_pick_uses_vep_pick_then_canonical_then_severity() -> None:
    key = "1:100:A:G"
    rows = [_vep_row(key, "G1", f"T{i}", "intron_variant") for i in range(4)]
    rows.append(_vep_row(key, "G1", "T9", "missense_variant", pick="1"))
    vep = pl.DataFrame(rows)
    picked = pick_vep(vep)
    assert picked.height == 1 and picked.row(0, named=True)["transcript_id"] == "T9"
    no_pick = pick_vep(pl.DataFrame([_vep_row(key, "G1", "T1", "intron_variant"),
                                     _vep_row(key, "G1", "T2", "stop_gained,splice_region_variant")]))
    assert no_pick.row(0, named=True)["most_severe_consequence"] == "stop_gained"
    assert no_pick.row(0, named=True)["consequence"] == "stop_gained&splice_region_variant"
    assert most_severe("intron_variant&missense_variant") == "missense_variant"


def test_consequence_matrix_and_no_double_counting(tmp_path: Path) -> None:
    root, participants = _cohort(tmp_path)
    vep = pl.DataFrame(
        # var1 overlaps five transcripts: it must still count once per carrier
        [_vep_row("1:100:A:G", "ENSG1", f"ENST{i}", "missense_variant", pick="1" if i == 0 else None) for i in range(5)]
        + [_vep_row("1:200:C:T", "ENSG2", "ENST9", "stop_gained", pick="1"),
           _vep_row("1:300:G:A", "ENSG3", "ENST8", "splice_donor_variant", pick="1")])
    ann = pl.DataFrame({"variant_key": ["1:100:A:G", "1:200:C:T", "1:300:G:A"]}).join(
        pick_vep(vep), on="variant_key", how="left").with_columns(
        pl.lit(None, dtype=pl.Utf8).alias("clinvar_significance"), pl.Series("gnomad_af", [0.2, 0.001, None]),
        pl.Series("alphamissense_score", [0.9, None, None]), pl.lit(None, dtype=pl.Utf8).alias("alphamissense_class"),
        pl.Series("spliceai_max_score", [None, None, 0.85]))
    out = aggregate(root, ann, tmp_path / "burden", participants, ParticipantVariantsConfig())
    counts = pl.read_parquet(out["counts"]).sort("participant_id")
    got = counts.select("participant_id", "n_missense_variant", "n_stop_gained", "n_splice_donor_variant").rows()
    assert got == [("P001", 1, 1, 0), ("P002", 1, 0, 0), ("P003", 0, 0, 1)]
    assert counts.get_column("n_variants_total").to_list() == [2, 1, 1]  # not 5 per missense carrier
    assert counts.select("alt_burden_missense_variant", "alt_burden_stop_gained").rows() == [(1.0, 2.0), (2.0, 0),
                                                                                          (0, 0)]
    assert counts.get_column("n_lof").to_list() == [1, 0, 1]
    assert counts.get_column("n_damaging_missense").to_list() == [1, 1, 0]  # AlphaMissense 0.9 >= 0.564
    assert counts.get_column("n_spliceai_ge_0_8").to_list() == [0, 0, 1]
    assert counts.get_column("n_rare").to_list() == [1, 0, 0]
    for col in ("n_frameshift_variant", "n_intergenic_variant", "alt_burden_3_prime_UTR_variant"):
        assert counts.get_column(col).to_list() == [0, 0, 0]  # always present
    annotated = pl.concat([pl.read_parquet(f) for f in sorted(out["annotated"].rglob("*.parquet"))])
    assert annotated.filter(pl.col("variant_key") == "1:100:A:G").sort("participant_id").select(
        "participant_id", "alt_count", "most_severe_consequence").rows() == [
        ("P001", 1, "missense_variant"), ("P002", 2, "missense_variant")]
    gene = pl.read_parquet(out["gene_burden"])
    assert gene.filter(pl.col("gene_id") == "ENSG1").sort("participant_id").select(
        "participant_id", "n_variants", "n_alt_alleles", "n_missense", "max_alphamissense").rows() == [
        ("P001", 1, 1.0, 1, 0.9), ("P002", 1, 2.0, 1, 0.9)]
    assert gene.filter(pl.col("gene_id") == "ENSG2").select("n_lof").item() == 1
