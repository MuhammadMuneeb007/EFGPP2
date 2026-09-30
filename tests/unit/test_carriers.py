"""Participant carrier genotypes: counted allele verified against the reference FASTA."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from efgpp.data.genotype.build import FastaIndex
from efgpp.data.genotype.carriers import (
    decode_bed,
    extract_carriers_bed,
    extract_carriers_vcf,
    orient,
)
from efgpp.data.genotype.formats import GenotypeFileset, read_samples, resolve_fileset

from ._genotypes import write_bed, write_fasta


def _samples(fs: GenotypeFileset) -> pl.DataFrame:
    s = read_samples(fs)
    return s.with_columns(pl.col("IID").alias("native_id"), pl.col("IID").alias("participant_id"))


def _carriers(out: Path) -> pl.DataFrame:
    files = sorted((out / "participant_variants").rglob("*.parquet"))
    return pl.concat([pl.read_parquet(f) for f in files]).sort(["variant_key", "participant_id"])


def test_orientation_rules() -> None:
    assert orient("G", "A", "A").status == "match"
    o = orient("A", "G", "A")
    assert (o.status, o.ref, o.alt, o.a1_is_alt) == ("swapped", "A", "G", False)
    f = orient("C", "T", "A")  # T/C reported on the other strand = A/G
    assert (f.status, f.ref, f.alt) == ("strand_flipped", "A", "G")
    assert orient("A", "T", "A").strand_ambiguous
    assert orient("A", "T", "G").status == "mismatch"  # A/T is never strand-flipped
    assert orient("G", "C", "A").status == "mismatch"
    assert orient("0", "A", "A").status == "monomorphic"
    assert orient("G", "A", None).status == "unverified"


def test_bed_decoding() -> None:
    raw = np.array([[0b11_10_01_00]], dtype=np.uint8)  # samples: hom A1, missing, het, hom A2
    assert decode_bed(raw, 4).tolist() == [[2, -1, 1, 0]]


def test_ref_a_alt_g_counts_in_both_plink_orientations(tmp_path: Path) -> None:
    """REF=A ALT=G: A/A -> 0, A/G -> 1, G/G -> 2, whichever allele PLINK calls A1."""
    fasta = FastaIndex(write_fasta(tmp_path / "ref.fa", {"chr1": {100: "A", 200: "A"}}))
    samples = ["P001", "P002", "P003"]
    # var1: A1=G (ALT), A2=A (REF): the A1 count is the ALT count.
    # var2: A1=A (REF), A2=G: PLINK counts the REF allele, so ALT count = 2 - A1 count.
    fs = resolve_fileset(write_bed(tmp_path / "g", samples,
                                   [("1", "var1", 100, "G", "A"), ("1", "var2", 200, "A", "G")],
                                   [[0, 1, 2], [2, 1, 0]]), "bed")
    res = extract_carriers_bed(fs, _samples(fs), fasta, tmp_path / "out", source_id="GENO001",
                               genome_build="GRCh38", carrier_only=False)
    c = _carriers(tmp_path / "out")
    for key in ("1:100:A:G", "1:200:A:G"):
        got = c.filter(pl.col("variant_key") == key).select("participant_id", "hardcall_alt_count", "GT").rows()
        assert got == [("P001", 0, "0/0"), ("P002", 1, "0/1"), ("P003", 2, "1/1")], key
    assert res.ref_check == {"match": 1, "swapped": 1}
    sites = pl.read_parquet(res.sites_path)
    assert sites.select("variant_key", "a1", "a2", "ref_check").rows() == [
        ("1:100:A:G", "G", "A", "match"), ("1:200:A:G", "A", "G", "swapped")]


def test_carrier_only_excludes_non_carriers(tmp_path: Path) -> None:
    fasta = FastaIndex(write_fasta(tmp_path / "ref.fa", {"1": {100: "A"}}))
    fs = resolve_fileset(write_bed(tmp_path / "g", ["P001", "P002", "P003"], [("1", "rs1", 100, "G", "A")],
                                   [[0, 1, 2]]), "bed")
    res = extract_carriers_bed(fs, _samples(fs), fasta, tmp_path / "out", source_id="GENO001", genome_build="GRCh38")
    c = _carriers(tmp_path / "out")
    assert c.select("participant_id", "hardcall_alt_count", "dosage", "heterozygous", "homozygous_alt").rows() == [
        ("P002", 1, 1.0, True, False), ("P003", 2, 2.0, False, True)]
    assert c.get_column("rsid").to_list() == ["rs1", "rs1"]
    assert set(c.columns) >= {"variant_key", "ploidy", "GT", "native_sample_id", "genome_build", "source_id"}
    assert res.n_rows == 2
    assert pl.read_parquet(res.sites_path).select("n_carriers", "n_hom_alt", "alt_allele_count").row(0) == (2, 1, 3)


def test_mismatch_and_missing(tmp_path: Path) -> None:
    fasta = FastaIndex(write_fasta(tmp_path / "ref.fa", {"1": {100: "A", 150: "T"}}))
    fs = resolve_fileset(write_bed(tmp_path / "g", ["P1", "P2"],
                                   [("1", "v1", 100, "G", "A"), ("1", "v2", 150, "G", "C")],
                                   [[-1, 2], [2, 2]]), "bed")
    res = extract_carriers_bed(fs, _samples(fs), fasta, tmp_path / "out", source_id="G", genome_build="GRCh38",
                               carrier_only=False)
    c = _carriers(tmp_path / "out")
    assert c.get_column("variant_key").unique().to_list() == ["1:100:A:G"]  # v2 matches neither allele
    assert c.filter(pl.col("participant_id") == "P1").select("GT", "hardcall_alt_count").row(0) == ("./.", None)
    assert res.ref_check == {"match": 1, "mismatch": 1}


def test_haploid_and_male_x(tmp_path: Path) -> None:
    fasta = FastaIndex(write_fasta(tmp_path / "ref.fa", {"X": {100: "A"}, "MT": {50: "A"}}))
    fs = resolve_fileset(write_bed(tmp_path / "g", ["M1", "F1"], [("23", "x1", 100, "G", "A"), ("26", "m1", 50, "G", "A")],
                                   [[2, 1], [2, 0]], sex=[1, 2]), "bed")
    extract_carriers_bed(fs, _samples(fs), fasta, tmp_path / "out", source_id="G", genome_build="GRCh38")
    c = _carriers(tmp_path / "out")
    assert c.select("variant_key", "participant_id", "hardcall_alt_count", "ploidy", "GT").rows() == [
        ("MT:50:A:G", "M1", 1, 1, "1"), ("X:100:A:G", "F1", 1, 2, "0/1"), ("X:100:A:G", "M1", 1, 1, "1")]


def test_vcf_dosage_kept_separately(tmp_path: Path) -> None:
    fasta = FastaIndex(write_fasta(tmp_path / "ref.fa", {"1": {100: "A"}}))
    vcf = tmp_path / "d.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP001\tP002\tP003\n"
        "1\t100\trs1\tA\tG\t.\t.\tR2=0.91\tGT:DS\t0/0:0.02\t0/1:0.93\t1/1:1.97\n"
        "1\t120\trs2\tT\tG\t.\t.\t.\tGT:DS\t0/1:1\t0/0:0\t0/0:0\n", encoding="utf-8")
    samples = pl.DataFrame({"IID": ["P001", "P002", "P003"], "native_id": ["P001", "P002", "P003"],
                            "participant_id": ["P001", "P002", "P003"], "SEX": [None, None, None]},
                           schema_overrides={"SEX": pl.Utf8})
    res = extract_carriers_vcf(vcf, samples, fasta, tmp_path / "out", source_id="G", genome_build="GRCh38")
    c = _carriers(tmp_path / "out")
    got = c.select("participant_id", "hardcall_alt_count", "dosage").rows()
    # P001 is 0/0 but has dosage 0.02 > 0: kept, hard call and dosage stay separate (never rounded).
    assert [(p, h, round(d, 2)) for p, h, d in got] == [("P001", 0, 0.02), ("P002", 1, 0.93), ("P003", 2, 1.97)]
    assert round(c.get_column("info_r2")[0], 2) == 0.91
    assert res.dosage and res.ref_check == {"match": 1, "mismatch": 1}  # rs2: REF T but the FASTA has C
