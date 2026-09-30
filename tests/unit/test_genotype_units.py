from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import polars as pl

from efgpp.constants import GenomeBuild
from efgpp.data.genotype.build import (
    CHROM_LENGTHS,
    FastaIndex,
    evidence_from_vcf_header,
    infer_build,
)
from efgpp.data.genotype.formats import (
    bgen_samples,
    read_samples,
    read_variants,
    resolve_fileset,
    vcf_samples,
)
from efgpp.data.genotype.liftover import Lifter, lift_table
from efgpp.data.genotype.loader import read_plink_table
from efgpp.data.genotype.relatedness import classify_pairs, greedy_unrelated_removal
from efgpp.data.simulate import read_bed_dosage, simulate_genotypes, write_bed


def test_bed_roundtrip(tmp_path: Path) -> None:
    dosage = np.array([[0, 1, 2], [2, -1, 0], [1, 1, 1], [0, 0, 2], [2, 2, 2]], dtype=np.int8)
    write_bed(tmp_path / "t", dosage, [f"S{i}" for i in range(5)], ["1", "1", "2"], [100, 200, 300],
              ["A", "C", "G"], ["G", "T", "A"])
    assert (read_bed_dosage(tmp_path / "t", 5, 3) == dosage).all()
    fs = resolve_fileset(tmp_path / "t.bed")
    assert fs.format == "bed" and fs.complete
    assert read_samples(fs).get_column("IID").to_list() == [f"S{i}" for i in range(5)]
    v = read_variants(fs)
    assert v.get_column("reference").to_list() == ["G", "T", "A"]  # A2 is REF in .bim


def test_pgen_psam_and_vcf_samples(tmp_path: Path) -> None:
    (tmp_path / "c.pgen").write_bytes(b"")
    (tmp_path / "c.pvar").write_text("##fileformat\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t10\trs1\tA\tG\n")
    (tmp_path / "c.psam").write_text("#FID\tIID\tSEX\nF1\tI1\t1\nF2\tI2\t2\n")
    fs = resolve_fileset(tmp_path / "c", "auto")
    assert fs.format == "pgen"
    assert read_samples(fs).get_column("IID").to_list() == ["I1", "I2"]
    assert read_variants(fs).row(0) == ("1", 10, "rs1", "A", "G")
    vcf = tmp_path / "x.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n")
    assert vcf_samples(vcf) == ["A", "B"]


def test_bgen_embedded_samples(tmp_path: Path) -> None:
    import struct

    ids = [b"S1", b"S22"]
    block = b"".join(struct.pack("<H", len(i)) + i for i in ids)
    header_len = 20
    flags = (1 << 31) | 0b100 | 1
    body = struct.pack("<IIII", header_len, 7, len(ids), 0) + struct.pack("<I", flags)
    sample_block = struct.pack("<II", 8 + len(block), len(ids)) + block
    path = tmp_path / "x.bgen"
    path.write_bytes(struct.pack("<I", header_len + len(sample_block)) + body + sample_block)
    assert bgen_samples(path) == ["S1", "S22"]


def test_build_inference_rules(tmp_path: Path) -> None:
    ev = evidence_from_vcf_header(["##contig=<ID=chr1,length=248956422>", "##reference=GRCh38.fa"])
    assert {e.supports for e in ev} == {GenomeBuild.GRCH38}
    assert infer_build(declared="auto", vcf_meta=["##contig=<ID=1,length=249250621>"]).build == GenomeBuild.GRCH37
    # Declared build alone is accepted; nothing at all requires the user.
    assert infer_build(declared="hg19").build == GenomeBuild.GRCH37
    assert infer_build(declared="auto").requires_user
    # Coordinate bounds alone are too weak to decide.
    beyond = pl.DataFrame({"chromosome": ["1"], "position": [249_000_000], "variant_id": ["v"],
                           "reference": ["A"], "alternate": ["G"]})
    weak = infer_build(declared="auto", variants=beyond)
    assert weak.requires_user and any(e.supports == GenomeBuild.GRCH37 for e in weak.evidence)
    # Declared build contradicted by the data -> stop and ask.
    assert infer_build(declared="GRCh38", vcf_meta=["##contig=<ID=1,length=249250621>"]).requires_user


def test_fasta_concordance(tmp_path: Path) -> None:
    rng = np.random.default_rng(0)
    seq37 = "".join(rng.choice(list("ACGT"), 5000))
    seq38 = "".join(rng.choice(list("ACGT"), 5000))
    for name, seq in (("b37.fa", seq37), ("b38.fa", seq38)):
        (tmp_path / name).write_text(">chr1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    assert FastaIndex(tmp_path / "b38.fa").base("1", 61) == seq38[60]
    pos = list(range(1, 5000, 17))
    variants = pl.DataFrame({"chromosome": ["1"] * len(pos), "position": pos, "variant_id": [str(p) for p in pos],
                             "reference": [seq38[p - 1] for p in pos],
                             "alternate": ["A" if seq38[p - 1] != "A" else "C" for p in pos]})
    result = infer_build(declared="auto", variants=variants,
                         fastas={GenomeBuild.GRCH37: tmp_path / "b37.fa", GenomeBuild.GRCH38: tmp_path / "b38.fa"})
    assert result.build == GenomeBuild.GRCH38 and result.confident


def test_chrom_lengths_distinguish_builds() -> None:
    assert CHROM_LENGTHS[GenomeBuild.GRCH37]["1"] != CHROM_LENGTHS[GenomeBuild.GRCH38]["1"]


def test_plink_report_parsing_and_relatedness(tmp_path: Path) -> None:
    (tmp_path / "x.kin0").write_text(
        "#FID1\tIID1\tFID2\tIID2\tNSNP\tHETHET\tIBS0\tKINSHIP\n"
        "A\tA\tB\tB\t100\t0.1\t0\t0.49\nA\tA\tC\tC\t100\t0.1\t0\t0.20\nD\tD\tE\tE\t100\t0.1\t0\t0.05\n")
    kin = classify_pairs(read_plink_table(tmp_path / "x.kin0"))
    assert kin.get_column("relationship").to_list() == ["duplicate_or_mz_twin", "first_degree", "third_degree"]
    removed = greedy_unrelated_removal(kin.filter(pl.col("KINSHIP") > 0.0884), {"A": 0.01, "B": 0.0, "C": 0.0})
    assert removed == {"A"}  # A is in two pairs: removing it resolves both
    (tmp_path / "x.smiss").write_text("#IID  MISSING_CT OBS_CT F_MISS\nS1 1 100 0.01\n")
    assert read_plink_table(tmp_path / "x.smiss").get_column("F_MISS").to_list() == [0.01]


def test_chain_liftover(tmp_path: Path) -> None:
    chain = tmp_path / "t.over.chain"
    chain.write_text(
        "chain 1000 chr1 1000 + 0 300 chr1 1000 + 10 310 1\n100 50 50\n150\n\n"
        "chain 500 chr2 1000 + 0 100 chr3 1000 - 0 100 2\n100\n")
    lifter = Lifter(chain)  # pyliftover over a local chain file
    assert lifter.lift("1", 1) == ("1", 11, "+")
    assert lifter.lift("chr1", 160) == ("1", 170, "+")
    assert lifter.lift("1", 120) is None  # in a gap
    v = pl.DataFrame({"variant_id": ["a", "b", "c", "d", "e"], "chromosome": ["1", "1", "2", "5", "1"],
                      "position": [1, 120, 5, 5, 0]})
    out = lift_table(v, lifter)
    assert out.get_column("status").to_list() == ["mapped", "unmapped", "chromosome_changed", "unmapped", "unplaced"]


def test_simulated_cohort_is_deterministic(tmp_path: Path) -> None:
    a = simulate_genotypes(20, 30, 2, 5, tmp_path / "a")
    b = simulate_genotypes(20, 30, 2, 5, tmp_path / "b")
    assert (a.dosage == b.dosage).all() and a.variant_ids == b.variant_ids
    assert (tmp_path / "a.bed").read_bytes() == (tmp_path / "b.bed").read_bytes()
