"""Tiny synthetic genotype files for tests (PLINK .bed/.bim/.fam, FASTA, PredictDB models)."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np

# A1 count -> 2-bit .bed code (00 hom A1, 10 het, 11 hom A2, 01 missing)
_CODE = {2: 0b00, 1: 0b10, 0: 0b11, -1: 0b01}


def write_bed(prefix: Path, samples: list[str], variants: list[tuple[str, str, int, str, str]],
              a1_counts: list[list[int]], sex: list[int] | None = None) -> Path:
    """variants: (chrom, id, pos, A1, A2); a1_counts: one row per variant, one value per sample."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    sexes = sex or [0] * len(samples)
    with open(prefix.with_suffix(".fam"), "w", encoding="utf-8", newline="\n") as fh:
        for s, x in zip(samples, sexes, strict=True):
            fh.write(f"{s} {s} 0 0 {x} -9\n")
    with open(prefix.with_suffix(".bim"), "w", encoding="utf-8", newline="\n") as fh:
        for chrom, vid, pos, a1, a2 in variants:
            fh.write(f"{chrom}\t{vid}\t0\t{pos}\t{a1}\t{a2}\n")
    bpv = (len(samples) + 3) // 4
    with open(prefix.with_suffix(".bed"), "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        for row in a1_counts:
            buf = np.zeros(bpv, dtype=np.uint8)
            for j, c in enumerate(row):
                buf[j // 4] |= _CODE[c] << (2 * (j % 4))
            fh.write(buf.tobytes())
    return prefix


def write_fasta(path: Path, contigs: dict[str, dict[int, str]], length: int = 400, width: int = 60) -> Path:
    """FASTA with 'C' everywhere except the given {1-based position: base} per contig."""
    from efgpp.data.genotype.build import build_fai

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for name, bases in contigs.items():
            seq = ["C"] * length
            for pos, b in bases.items():
                seq[pos - 1] = b
            s = "".join(seq)
            fh.write(f">{name}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")
    build_fai(path)
    return path


def write_predictdb(path: Path, weights: list[tuple[str, str, str, str, str, float]],
                    extra: list[tuple[str, str, int, float]] | None = None) -> Path:
    """PredictDB SQLite as in GTEx v8 MASHR: weights (gene, rsid, varID, ref_allele, eff_allele, weight) and
    extra (gene, genename, n.snps.in.model, pred.perf.R2)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE weights (gene TEXT, rsid TEXT, varID TEXT, ref_allele TEXT, eff_allele TEXT, weight DOUBLE)")
    con.executemany("INSERT INTO weights VALUES (?, ?, ?, ?, ?, ?)", weights)
    con.execute('CREATE TABLE extra (gene TEXT, genename TEXT, "n.snps.in.model" INTEGER, "pred.perf.R2" DOUBLE)')
    genes = extra or [(g, g, sum(1 for w in weights if w[0] == g), 0.1) for g in dict.fromkeys(w[0] for w in weights)]
    con.executemany("INSERT INTO extra VALUES (?, ?, ?, ?)", genes)
    con.commit()
    con.close()
    return path
