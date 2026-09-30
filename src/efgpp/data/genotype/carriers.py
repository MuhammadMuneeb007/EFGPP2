"""Participant carrier genotypes: which ALT alleles each participant actually carries.

The cohort variant table says which variants exist somewhere in the cohort; this table says
participant X carries variant Y with ALT dosage Z. ALT is always the non-reference allele of
the project's reference FASTA (GRCh38). The counted allele is verified, never assumed:

  PLINK .bed / hard-call .pgen   A1/A2 are compared with the FASTA base. A2 = REF (what PLINK
                                 provisionally assumes) keeps the A1 count; A1 = REF means the ALT
                                 count is 2 - A1 count; strand flips are complemented; A/T and C/G
                                 sites are flagged as strand-ambiguous; sites matching neither
                                 allele are excluded and reported. Without a FASTA the PLINK
                                 convention is used and every site is labelled `unverified`.
  VCF / BCF / BGEN (dosages)     exported by PLINK 2 with REF taken from the FASTA
                                 (--ref-from-fa), normalized with `bcftools norm` (left-align,
                                 REF check, split multiallelics); GT and DS are read per
                                 participant and dosage is never rounded.

Canonical variant key: chromosome:position:REF:ALT (rsID is kept as metadata only).

Storage (never a dense participant x variant matrix):
  data/derived/participant_variants/<GENOTYPE_ID>/
      participant_variants/chromosome=<c>/part-*.parquet   carrier rows (ALT dosage > 0 by default)
      variants.parquet                                      one row per variant: key, REF check, counts
"""

from __future__ import annotations

import gzip
import shutil
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.formats import (
    GenotypeFileset,
    read_samples,
    read_variants,
    resolve_fileset,
)
from efgpp.data.io import write_parquet
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths

COMPLEMENT = str.maketrans("ACGT", "TGCA")
NUCLEOTIDES = set("ACGT")
AMBIGUOUS_PAIRS = {frozenset("AT"), frozenset("CG")}
HAPLOID_CHROMS = {"Y", "MT"}
# PLINK numeric chromosome codes -> canonical names (XY = pseudo-autosomal part of X).
CHROM_ALIASES = {"23": "X", "24": "Y", "25": "X", "XY": "X", "26": "MT", "M": "MT", "0": "0"}
# .bed 2-bit codes -> number of A1 alleles (01 = missing).
BED_A1_COUNT = np.array([2, -1, 1, 0], dtype=np.int8)
GT_DICTIONARY = ["0/0", "0/1", "1/1", "./.", "0", "1", "."]

CARRIER_COLUMNS = [
    "participant_id", "native_sample_id", "variant_key", "variant_id", "rsid", "chromosome", "position",
    "ref", "alt", "GT", "hardcall_alt_count", "dosage", "ploidy", "heterozygous", "homozygous_alt",
    "info_r2", "genome_build", "source_id",
]


def canonical_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return CHROM_ALIASES.get(c.upper(), c.upper() if c.isalpha() else c)


def variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{canonical_chrom(chrom)}:{int(pos)}:{ref}:{alt}"


def complement(allele: str) -> str:
    return allele.translate(COMPLEMENT)


# ------------------------------------------------------------------ allele orientation
@dataclass(frozen=True)
class Orientation:
    ref: str | None
    alt: str | None
    a1_is_alt: bool  # True: ALT count = A1 count; False: ALT count = 2 - A1 count
    status: str  # match | swapped | strand_flipped | strand_flipped_swapped | unverified | mismatch | ...
    strand_ambiguous: bool = False

    @property
    def keep(self) -> bool:
        return self.status in ("match", "swapped", "strand_flipped", "strand_flipped_swapped", "unverified")


def orient(a1: str, a2: str, reference: str | None) -> Orientation:
    """Decide REF/ALT of a PLINK A1/A2 site from the reference sequence starting at its position.

    `reference` is the FASTA sequence at the site (at least as long as the longer allele), or
    None when no FASTA is available (then the PLINK convention A2 = REF is used, unverified).
    """
    a1, a2 = a1.upper(), a2.upper()
    if a1 in ("0", ".", "") or a2 in ("0", ".", ""):
        return Orientation(None, None, True, "monomorphic")
    ambiguous = frozenset((a1, a2)) in AMBIGUOUS_PAIRS
    if reference is None:
        return Orientation(a2, a1, True, "unverified", ambiguous)
    if not (set(a1) <= NUCLEOTIDES and set(a2) <= NUCLEOTIDES):
        return Orientation(None, None, True, "non_nucleotide")
    ref_seq = reference.upper()

    def at_ref(allele: str) -> bool:
        return len(ref_seq) >= len(allele) and ref_seq[: len(allele)] == allele

    if at_ref(a2):
        return Orientation(a2, a1, True, "match", ambiguous)
    if at_ref(a1):
        return Orientation(a1, a2, False, "swapped", ambiguous)
    if not ambiguous:  # A/T and C/G sites cannot be told apart from a strand flip
        c1, c2 = complement(a1), complement(a2)
        if at_ref(c2):
            return Orientation(c2, c1, True, "strand_flipped")
        if at_ref(c1):
            return Orientation(c1, c2, False, "strand_flipped_swapped")
    return Orientation(None, None, True, "mismatch", ambiguous)


def orient_sites(variants: pl.DataFrame, fasta: Any | None) -> pl.DataFrame:
    """Add ref, alt, a1_is_alt, ref_check, strand_ambiguous, keep to a .bim-style table
    (columns chromosome, position, variant_id, reference = A2, alternate = A1)."""
    chrom = [canonical_chrom(c) for c in variants.get_column("chromosome").to_list()]
    pos = variants.get_column("position").to_list()
    a1 = [str(x or "") for x in variants.get_column("alternate").to_list()]
    a2 = [str(x or "") for x in variants.get_column("reference").to_list()]
    refs: list[str | None] = [None] * len(pos)
    if fasta is not None:
        by_chrom: dict[str, list[int]] = {}
        for i, c in enumerate(chrom):
            by_chrom.setdefault(c, []).append(i)
        for c, idx in by_chrom.items():
            seqs = fasta.bases_many(c, [pos[i] for i in idx], [max(len(a1[i]), len(a2[i]), 1) for i in idx])
            for i, s in zip(idx, seqs, strict=True):
                refs[i] = s if s is not None else ""  # "" = position outside the FASTA -> mismatch
    o = [orient(x, y, r) for x, y, r in zip(a1, a2, refs, strict=True)]
    return variants.with_columns(
        pl.Series("chromosome", chrom, dtype=pl.Utf8),
        pl.Series("a1", a1, dtype=pl.Utf8),
        pl.Series("a2", a2, dtype=pl.Utf8),
        pl.Series("ref", [x.ref for x in o], dtype=pl.Utf8),
        pl.Series("alt", [x.alt for x in o], dtype=pl.Utf8),
        pl.Series("a1_is_alt", [x.a1_is_alt for x in o], dtype=pl.Boolean),
        pl.Series("ref_check", [x.status for x in o], dtype=pl.Utf8),
        pl.Series("strand_ambiguous", [x.strand_ambiguous for x in o], dtype=pl.Boolean),
        pl.Series("keep", [x.keep for x in o], dtype=pl.Boolean),
    ).with_columns(
        pl.when(pl.col("keep"))
        .then(pl.concat_str([pl.col("chromosome"), pl.col("position").cast(pl.Utf8), pl.col("ref"), pl.col("alt")],
                            separator=":"))
        .otherwise(None).alias("variant_key"),
        pl.when(pl.col("variant_id").str.starts_with("rs")).then(pl.col("variant_id")).otherwise(None).alias("rsid"),
    )


# ------------------------------------------------------------------ .bed decoding
def decode_bed(raw: np.ndarray, n_samples: int) -> np.ndarray:
    """(variants x bytes) uint8 block of a SNP-major .bed -> (variants x samples) A1 counts, -1 = missing."""
    codes = np.stack([(raw >> s) & 3 for s in (0, 2, 4, 6)], axis=-1).reshape(raw.shape[0], -1)[:, :n_samples]
    return BED_A1_COUNT[codes]


def iter_bed(bed: Path, n_samples: int, n_variants: int, chunk: int) -> Iterator[tuple[int, np.ndarray]]:
    """Yield (first variant index, A1 counts) blocks of a SNP-major .bed file."""
    bpv = (n_samples + 3) // 4
    with open(bed, "rb") as fh:
        magic = fh.read(3)
        if magic != b"\x6c\x1b\x01":
            raise ValueError(f"{bed}: not a SNP-major PLINK .bed file")
        start = 0
        while start < n_variants:
            k = min(chunk, n_variants - start)
            buf = fh.read(bpv * k)
            if len(buf) != bpv * k:
                raise ValueError(f"{bed}: truncated (expected {n_variants} variants x {n_samples} samples)")
            yield start, decode_bed(np.frombuffer(buf, dtype=np.uint8).reshape(k, bpv), n_samples)
            start += k


def ploidy_matrix(chrom: np.ndarray, male: np.ndarray, male_haploid_sites: np.ndarray) -> np.ndarray:
    """(sites x samples) ploidy: 1 on chrY/MT, 1 for males on non-PAR chrX, else 2."""
    out = np.full((len(chrom), len(male)), 2, dtype=np.int8)
    out[np.isin(chrom, list(HAPLOID_CHROMS))] = 1
    if male_haploid_sites.any() and male.any():
        out[np.ix_(male_haploid_sites, male)] = 1
    return out


def alt_counts(a1_counts: np.ndarray, a1_is_alt: np.ndarray, ploidy: np.ndarray) -> np.ndarray:
    """A1 counts -> ALT counts (-1 missing). Haploid calls (ploidy 1; stored as homozygous in .bed)
    count 0 or 1; a heterozygous haploid call is an error and becomes missing."""
    out = np.where(a1_is_alt[:, None], a1_counts, 2 - a1_counts).astype(np.int8)
    out[a1_counts < 0] = -1
    hap = ploidy == 1
    out[hap & (out == 1)] = -1
    out[hap & (out == 2)] = 1
    return out


# ------------------------------------------------------------------ writing
class PartitionWriter:
    """Parquet files partitioned by chromosome (hive layout: chromosome=<c>/part-<n>.parquet)."""

    def __init__(self, root: Path, partition: bool, compression: str = "zstd") -> None:
        self.root = root
        self.partition = partition
        self.compression = compression
        self.writers: dict[str, pq.ParquetWriter] = {}
        self.rows = 0
        if root.exists():
            shutil.rmtree(root)
        root.mkdir(parents=True)

    def write(self, chrom: str, table: pa.Table) -> None:
        if table.num_rows == 0:
            return
        key = chrom if self.partition else "all"
        w = self.writers.get(key)
        if w is None:
            d = self.root / f"chromosome={key}" if self.partition else self.root
            d.mkdir(parents=True, exist_ok=True)
            w = pq.ParquetWriter(d / "part-0.parquet", table.schema, compression=self.compression)
            self.writers[key] = w
        w.write_table(table)
        self.rows += table.num_rows

    def close(self) -> None:
        for w in self.writers.values():
            w.close()


def _dict(indices: np.ndarray, values: list[Any] | np.ndarray, typ: pa.DataType = pa.string()) -> pa.DictionaryArray:
    """Dictionary-encoded column (no per-row Python strings); None values become nulls."""
    # Deduplicated dictionary (first-seen order) so identical input always gives identical bytes.
    uniq: dict[Any, int] = {}
    dict_vals: list[Any] = []
    remap = np.zeros(len(values), dtype=np.int32)
    none = np.zeros(len(values), dtype=bool)
    for i, v in enumerate(values):
        if v is None:
            none[i] = True
            continue
        j = uniq.get(v)
        if j is None:
            j = uniq[v] = len(dict_vals)
            dict_vals.append(v)
        remap[i] = j
    idx = remap[indices]
    mask = none[indices]
    return pa.DictionaryArray.from_arrays(pa.array(idx, mask=mask) if mask.any() else pa.array(idx),
                                          pa.array(dict_vals, type=typ))


def _const(n: int, value: str | None) -> pa.Array:
    return pa.DictionaryArray.from_arrays(pa.array(np.zeros(n, dtype=np.int32)), pa.array([value], type=pa.string()))


def carrier_block(counts: np.ndarray, sites: dict[str, list[Any]], ploidy: np.ndarray, native: list[str],
                  participants: list[str | None], *, carrier_only: bool, genome_build: str | None,
                  source_id: str, dosage: np.ndarray | None = None, gt: list[list[str]] | None = None,
                  info_r2: np.ndarray | None = None) -> pa.Table:
    """Carrier rows of one block of sites (counts and ploidy: sites x samples; ALT count -1 = missing)."""
    mask = counts > 0 if carrier_only else np.ones_like(counts, dtype=bool)
    if dosage is not None and carrier_only:
        mask = mask | (np.nan_to_num(dosage, nan=0.0) > 0)
    vi, si = np.nonzero(mask)
    c = counts[vi, si]
    p = ploidy[vi, si] if ploidy.ndim == 2 else ploidy[vi]
    missing = c < 0
    n = len(vi)
    if gt is None:
        gt_idx = np.where(p == 2, np.where(missing, 3, c), np.where(missing, 6, 4 + np.clip(c, 0, 1)))
        gt_arr: pa.Array = _dict(gt_idx, GT_DICTIONARY)
    else:
        gt_arr = pa.array([gt[v][s] for v, s in zip(vi.tolist(), si.tolist(), strict=True)], type=pa.string())
    hard = pa.array(np.where(missing, 0, c).astype(np.int8), mask=missing)
    if dosage is not None:
        d = dosage[vi, si].astype(np.float32)
        dose = pa.array(d, mask=np.isnan(d))
    else:
        dose = pa.array(np.where(missing, 0, c).astype(np.float32), mask=missing)
    r2 = pa.array(info_r2[vi].astype(np.float32), mask=np.isnan(info_r2[vi])) if info_r2 is not None \
        else pa.nulls(n, pa.float32())
    cols = {
        "participant_id": _dict(si, participants),
        "native_sample_id": _dict(si, native),
        "variant_key": _dict(vi, sites["variant_key"]),
        "variant_id": _dict(vi, sites["variant_id"]),
        "rsid": _dict(vi, sites["rsid"]),
        "chromosome": _dict(vi, sites["chromosome"]),
        "position": pa.array(np.asarray(sites["position"], dtype=np.int64)[vi]),
        "ref": _dict(vi, sites["ref"]),
        "alt": _dict(vi, sites["alt"]),
        "GT": gt_arr,
        "hardcall_alt_count": hard,
        "dosage": dose,
        "ploidy": pa.array(p.astype(np.int8)),
        "heterozygous": pa.array((c == 1) & (p == 2), mask=missing),
        "homozygous_alt": pa.array(c == p, mask=missing),
        "info_r2": r2,
        "genome_build": _const(n, genome_build),
        "source_id": _const(n, source_id),
    }
    return pa.table(cols)


# ------------------------------------------------------------------ genotype resolution
@dataclass
class CarrierResult:
    out_dir: Path
    carrier_root: Path
    sites_path: Path
    n_participants: int
    n_rows: int
    n_sites: int
    ref_check: dict[str, int] = field(default_factory=dict)
    mode: str = "bed"
    dosage: bool = False
    participant_ids: list[str] = field(default_factory=list)


def resolve_genotype_artifact(project: Project, source_id: str, use_qc: bool = True) -> Artifact:
    """QC-passed genotype when present (and wanted), else the copy in the target build, else the source."""
    target = project.config.defaults.target_build
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        art = store.latest(source_id=source_id, artifact_type="qc_genotype") if use_qc else None
        if art is None:
            lifted = store.latest(source_id=source_id, artifact_type="liftover_genotype")
            src = store.latest(source_id=source_id, artifact_type="source")
            art = lifted if lifted is not None and lifted.genome_build == target else src
    if art is None:
        raise RuntimeError(f"{source_id} is not registered; run `efgpp data prepare` first")
    return art


def genotype_samples(project: Project, source_id: str, fs: GenotypeFileset) -> pl.DataFrame:
    """Samples of a genotype fileset with native ids, participant ids and sex (1 = male)."""
    from efgpp.data.aliases import plink_native_ids

    _, src = project.data.get_source(source_id)
    samples = read_samples(fs)
    sid = src.sample_id  # type: ignore[attr-defined]
    if fs.format in ("bed", "pgen"):
        samples = samples.with_columns(plink_native_ids(samples, sid.mode, sid.column, sid.separator))
    elif "native_id" not in samples.columns:
        samples = samples.with_columns(pl.col("IID").alias("native_id"))
    with Registry.open(project) as reg:
        aliases = reg.frame("SELECT native_id, participant_id FROM sample_aliases WHERE source_id = ?", [source_id])
    if "SEX" not in samples.columns:
        samples = samples.with_columns(pl.lit(None, dtype=pl.Utf8).alias("SEX"))
    # Row order must stay the .fam/.psam order: genotype columns are matched to samples by position.
    out = samples.join(aliases.unique("native_id", keep="first"), on="native_id", how="left", maintain_order="left")
    if out.get_column("native_id").to_list() != samples.get_column("native_id").to_list():
        raise RuntimeError("sample order changed while resolving participant ids")
    return out


def _source_format(project: Project, source_id: str) -> str:
    _, src = project.data.get_source(source_id)
    fmt = getattr(src, "format", "auto")
    if fmt in (None, "auto"):
        with Registry.open(project) as reg:
            art = ArtifactStore(reg).latest(source_id=source_id, artifact_type="source")
        fmt = art.format if art else "bed"
    return str(fmt)


# ------------------------------------------------------------------ main entry points
def extract_carriers_bed(fs: GenotypeFileset, samples: pl.DataFrame, fasta: Any | None, out_dir: Path, *,
                         source_id: str, genome_build: str | None, carrier_only: bool = True,
                         partition: bool = True, chunk_cells: int = 20_000_000) -> CarrierResult:
    """Carrier rows from a PLINK .bed fileset (hard calls), with REF/ALT verified against `fasta`."""
    if fs.format != "bed":
        raise ValueError("extract_carriers_bed needs a .bed fileset")
    variants = orient_sites(read_variants(fs), fasta)
    n_samples = samples.height
    male = samples.get_column("SEX").cast(pl.Utf8).fill_null("").to_numpy() == "1"
    chrom = variants.get_column("chromosome").to_numpy()
    # non-PAR chrX (PLINK 23 / X); PLINK 25 / XY is the pseudo-autosomal region and stays diploid
    raw_chrom = read_variants(fs).get_column("chromosome").str.to_uppercase().str.replace("^CHR", "").to_numpy()
    male_haploid_all = np.isin(raw_chrom, ["23", "X"])
    keep = variants.get_column("keep").to_numpy()
    a1_is_alt = variants.get_column("a1_is_alt").to_numpy()
    cols = {c: variants.get_column(c).to_list() for c in ("variant_key", "variant_id", "rsid", "chromosome",
                                                          "position", "ref", "alt")}
    native = samples.get_column("native_id").to_list()
    participants = samples.get_column("participant_id").to_list()
    carrier_root = out_dir / "participant_variants"
    writer = PartitionWriter(carrier_root, partition)
    stats = {k: np.zeros(variants.height, dtype=np.int64) for k in ("n_called", "n_carriers", "n_hom_alt",
                                                                   "alt_allele_count")}
    chunk = max(1, chunk_cells // max(n_samples, 1))
    bed = fs.prefix.with_name(fs.prefix.name + ".bed")
    try:
        for start, a1 in iter_bed(bed, n_samples, variants.height, chunk):
            stop = start + a1.shape[0]
            sl = slice(start, stop)
            ploidy = ploidy_matrix(chrom[sl], male, male_haploid_all[sl])
            counts = alt_counts(a1, a1_is_alt[sl], ploidy)
            k = keep[sl]
            called = counts >= 0
            stats["n_called"][sl] = called.sum(axis=1)
            stats["n_carriers"][sl] = (counts > 0).sum(axis=1)
            stats["n_hom_alt"][sl] = (called & (counts == ploidy)).sum(axis=1)
            stats["alt_allele_count"][sl] = np.where(called, counts, 0).sum(axis=1)
            if not k.any():
                continue
            idx = np.nonzero(k)[0]
            sub = {c: [v[start + i] for i in idx] for c, v in cols.items()}
            block_chrom = np.asarray(sub["chromosome"])
            for c in dict.fromkeys(sub["chromosome"]):
                sel = np.nonzero(block_chrom == c)[0]
                sites = {name: [vals[i] for i in sel] for name, vals in sub.items()}
                table = carrier_block(counts[idx[sel]], sites, ploidy[idx[sel]], native, participants,
                                      carrier_only=carrier_only, genome_build=genome_build, source_id=source_id)
                writer.write(str(c), table)
    finally:
        writer.close()
    sites_table = variants.select(
        "variant_key", "variant_id", "rsid", "chromosome", "position",
        pl.col("ref").alias("reference"), pl.col("alt").alias("alternate"), "a1", "a2", "ref_check",
        "strand_ambiguous", pl.col("keep").alias("included"),
    ).with_columns([pl.Series(k, v) for k, v in stats.items()])
    sites_path = write_parquet(sites_table, out_dir / "variants.parquet")
    checks = {k: int(v) for k, v in variants.get_column("ref_check").value_counts().iter_rows()}
    return CarrierResult(out_dir, carrier_root, sites_path, n_samples, writer.rows, int(keep.sum()), checks, "bed",
                         participant_ids=[p for p in participants if p])


def _parse_info_r2(info: str) -> float:
    for kv in info.split(";"):
        k, _, v = kv.partition("=")
        if k in ("R2", "DR2", "INFO", "IMPINFO") and v:
            try:
                return float(v.split(",")[0])
            except ValueError:
                return float("nan")
    return float("nan")


def _gt_alt_count(gt: str) -> tuple[int, int]:
    """(ALT count, ploidy) of a biallelic VCF GT; ALT count -1 when missing."""
    alleles = gt.replace("|", "/").split("/")
    if any(a in (".", "") for a in alleles):
        return -1, len(alleles)
    return sum(1 for a in alleles if a != "0"), len(alleles)


def extract_carriers_vcf(vcf: Path, samples: pl.DataFrame, fasta: Any | None, out_dir: Path, *,
                         source_id: str, genome_build: str | None, carrier_only: bool = True,
                         partition: bool = True, batch_sites: int = 20000) -> CarrierResult:
    """Carrier rows from a (normalized, biallelic) VCF with GT and optional DS (dosage)."""
    opener = gzip.open if vcf.name.endswith((".gz", ".bgz")) else open
    by_native = {r["native_id"]: r for r in samples.to_dicts()}
    by_iid = {r["IID"]: r for r in samples.to_dicts()} if "IID" in samples.columns else {}
    carrier_root = out_dir / "participant_variants"
    writer = PartitionWriter(carrier_root, partition)
    site_rows: list[dict[str, Any]] = []
    checks: dict[str, int] = {}
    any_dosage = False
    native: list[str] = []
    participants: list[str | None] = []

    def flush(buf: list[tuple[dict[str, Any], list[int], list[float], list[str], list[int]]]) -> None:
        if not buf:
            return
        by_chrom: dict[str, list[int]] = {}
        for i, (site, *_rest) in enumerate(buf):
            by_chrom.setdefault(site["chromosome"], []).append(i)
        for c, idx in by_chrom.items():
            counts = np.array([buf[i][1] for i in idx], dtype=np.int8)
            dosage = np.array([buf[i][2] for i in idx], dtype=np.float32)
            gts = [buf[i][3] for i in idx]
            ploidy = np.array([buf[i][4] for i in idx], dtype=np.int8)
            sites = {k: [buf[i][0][k] for i in idx] for k in ("variant_key", "variant_id", "rsid", "chromosome",
                                                              "position", "ref", "alt")}
            r2 = np.array([buf[i][0]["info_r2"] for i in idx], dtype=np.float32)
            writer.write(c, carrier_block(counts, sites, ploidy, native, participants, carrier_only=carrier_only,
                                          genome_build=genome_build, source_id=source_id,
                                          dosage=dosage if any_dosage else None, gt=gts, info_r2=r2))
        buf.clear()

    buffer: list[tuple[dict[str, Any], list[int], list[float], list[str], list[int]]] = []
    try:
        with opener(vcf, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header_samples = line.rstrip("\n").split("\t")[9:]
                    for s in header_samples:
                        row = by_native.get(s) or by_iid.get(s) or by_native.get(s.split("_", 1)[-1])
                        native.append(row["native_id"] if row else s)
                        participants.append(row["participant_id"] if row else None)
                    continue
                f = line.rstrip("\n").split("\t")
                chrom, pos, vid, ref, alt = canonical_chrom(f[0]), int(f[1]), f[2], f[3].upper(), f[4].upper()
                if "," in alt:
                    checks["multiallelic_excluded"] = checks.get("multiallelic_excluded", 0) + 1
                    continue
                if alt in (".", "*", ""):
                    checks["monomorphic"] = checks.get("monomorphic", 0) + 1
                    continue
                status = "unverified"
                if fasta is not None:
                    seq = fasta.bases_many(chrom, [pos], [len(ref)])[0]
                    status = "match" if seq == ref else "mismatch"
                checks[status] = checks.get(status, 0) + 1
                if status == "mismatch":
                    continue
                fmt = f[8].split(":")
                gi = fmt.index("GT") if "GT" in fmt else None
                di = fmt.index("DS") if "DS" in fmt else None
                counts, dos, gts, ploidies = [], [], [], []
                for field_ in f[9:]:
                    parts = field_.split(":")
                    gt = parts[gi] if gi is not None and gi < len(parts) else "./."
                    c, ploidy = _gt_alt_count(gt)
                    counts.append(c)
                    ploidies.append(ploidy)
                    gts.append(gt)
                    if di is not None and di < len(parts) and parts[di] not in (".", ""):
                        dos.append(float(parts[di]))
                        any_dosage = True
                    else:
                        dos.append(float("nan") if c < 0 else float(c))
                key = variant_key(chrom, pos, ref, alt)
                site = {"variant_key": key, "variant_id": vid if vid != "." else key,
                        "rsid": vid if vid.startswith("rs") else None, "chromosome": chrom, "position": pos,
                        "ref": ref, "alt": alt, "info_r2": _parse_info_r2(f[7])}
                arr, parr = np.array(counts), np.array(ploidies)
                site_rows.append({**{k: site[k] for k in ("variant_key", "variant_id", "rsid", "chromosome",
                                                          "position")},
                                  "reference": ref, "alternate": alt, "ref_check": status, "included": True,
                                  "n_called": int((arr >= 0).sum()), "n_carriers": int((arr > 0).sum()),
                                  "n_hom_alt": int(((arr >= 0) & (arr == parr)).sum()),
                                  "alt_allele_count": int(np.where(arr >= 0, arr, 0).sum())})
                buffer.append((site, counts, dos, gts, ploidies))
                if len(buffer) >= batch_sites:
                    flush(buffer)
            flush(buffer)
    finally:
        writer.close()
    schema = {"variant_key": pl.Utf8, "variant_id": pl.Utf8, "rsid": pl.Utf8, "chromosome": pl.Utf8,
              "position": pl.Int64, "reference": pl.Utf8, "alternate": pl.Utf8, "ref_check": pl.Utf8,
              "included": pl.Boolean, "n_called": pl.Int64, "n_carriers": pl.Int64, "n_hom_alt": pl.Int64,
              "alt_allele_count": pl.Int64}
    sites_path = write_parquet(pl.DataFrame(site_rows, schema=schema), out_dir / "variants.parquet")
    return CarrierResult(out_dir, carrier_root, sites_path, len(native), writer.rows, len(site_rows), checks,
                         "vcf", any_dosage, [p for p in participants if p])


def run_participant_variants(project: Project, source_id: str, *, step_id: str | None = None,
                             threads: int = 1) -> dict[str, Any]:
    """Build and register the participant carrier table of one genotype source."""
    from efgpp.data.genotype.build import FastaIndex
    from efgpp.data.provenance import run_tool
    from efgpp.data.references.genome import installed_fasta
    from efgpp.setup.tools import available

    cfg = project.data.participant_variants
    art = resolve_genotype_artifact(project, source_id, cfg.use_qc_genotype)
    fs = resolve_fileset(Path(art.path), art.format)
    fasta_path = installed_fasta(project)
    fasta = FastaIndex(fasta_path) if fasta_path else None
    work = project.work_root / "participant_variants" / source_id
    work.mkdir(parents=True, exist_ok=True)
    out_dir = project.artifact_dir(Origin.DERIVED, Modality.PARTICIPANT_VARIANTS, source_id)
    dosage_source = _source_format(project, source_id) in ("bgen", "vcf", "bcf")
    tools = [art.tool or ""]
    if dosage_source:
        # Dosages: PLINK 2 exports a VCF with DS (REF from the FASTA), bcftools normalizes it.
        base = work / "export"
        args = [*fs.plink_input_args(), "--export", "vcf", "bgz", "vcf-dosage=DS-force", "--out", str(base),
                "--threads", str(threads)]
        if fasta_path:
            args += ["--fa", str(fasta_path), "--ref-from-fa", "force"]
        run_tool(project, "plink2", args, step_id=step_id, inputs=[art.artifact_id])  # type: ignore[list-item]
        vcf = base.with_suffix(".vcf.gz")
        if cfg.normalize and fasta_path and available(project, "bcftools"):
            norm = work / "normalized.vcf.gz"
            run_tool(project, "bcftools", ["norm", "-f", str(fasta_path), "-m", "-any", "-c", "x", "-Oz",
                                           "-o", str(norm), str(vcf)], step_id=step_id)
            vcf = norm
            tools.append("bcftools")
        samples = genotype_samples(project, source_id, fs)
        result = extract_carriers_vcf(vcf, samples, fasta, out_dir, source_id=source_id,
                                      genome_build=art.genome_build, carrier_only=cfg.carrier_only,
                                      partition=cfg.partition_by_chromosome)
    else:
        bed_fs = fs
        if fs.format != "bed":
            prefix = work / "hardcalls"
            run_tool(project, "plink2", [*fs.plink_input_args(), "--make-bed", "--out", str(prefix),
                                         "--threads", str(threads)], step_id=step_id,
                     inputs=[art.artifact_id])  # type: ignore[list-item]
            bed_fs = resolve_fileset(prefix, "bed")
        samples = genotype_samples(project, source_id, bed_fs)
        result = extract_carriers_bed(bed_fs, samples, fasta, out_dir, source_id=source_id,
                                      genome_build=art.genome_build, carrier_only=cfg.carrier_only,
                                      partition=cfg.partition_by_chromosome)
    return register_carriers(project, source_id, art, result, fasta_path, step_id)


def register_carriers(project: Project, source_id: str, genotype: Artifact, result: CarrierResult,
                      fasta: Path | None, step_id: str | None = None) -> dict[str, Any]:
    cfg = project.data.participant_variants
    files = [result.carrier_root, result.sites_path]
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        a = store.register_replacing(Artifact(
            artifact_name=f"{source_id}_participant_variants", artifact_type="participant_variants",
            modality=Modality.PARTICIPANT_VARIANTS, origin=Origin.DERIVED, status=ArtifactStatus.READY,
            path=str(result.out_dir), format="parquet_dataset",
            size=sum(f.stat().st_size for d in files for f in ([d] if d.is_file() else d.rglob("*")) if f.is_file()),
            checksum=str(checksum_paths(files, full_limit_bytes=_limit(project))), participant_count=result.n_participants,
            feature_count=result.n_sites, genome_build=genotype.genome_build, temporal_type=TemporalType.STATIC,
            source_id=source_id, parent_artifact_ids=[genotype.artifact_id],  # type: ignore[list-item]
            tool="efgpp", resource_versions={"reference_fasta": str(fasta) if fasta else "none"},
            metadata={
                "members": [str(f) for f in files], "carrier_rows": result.n_rows,
                "ref_check": result.ref_check, "mode": result.mode, "dosage": result.dosage,
                "carrier_only": cfg.carrier_only, "variant_key": "chromosome:position:REF:ALT",
                "counted_allele": "ALT relative to the reference FASTA" if fasta else
                                  "ALT = PLINK A1 (unverified: no reference FASTA installed)",
                "note": "phenotype-independent; derived from the genotype only",
            },
        ))
        write_assays(reg, a, f"{source_id}_{Modality.PARTICIPANT_VARIANTS.value}", Modality.PARTICIPANT_VARIANTS,
                     result.participant_ids)
    return {"artifact": a.artifact_id, "rows": result.n_rows, "sites": result.n_sites, "ref_check": result.ref_check}


def write_assays(reg: Registry, art: Artifact, assay_source: str, modality: Modality, participant_ids: list[str],
                 origin: Origin = Origin.DERIVED, tissue: str | None = None) -> None:
    """Record which participants an artifact covers (feeds the availability matrix)."""
    ids = sorted(set(participant_ids))
    reg.execute("DELETE FROM assays WHERE source_id = ?", [assay_source])
    if ids:
        reg.insert_frame("assays", pl.DataFrame({
            "assay_id": [f"{art.artifact_id}:{i}" for i in range(len(ids))],
            "artifact_id": [art.artifact_id] * len(ids),
            "source_id": [assay_source] * len(ids),
            "modality": [modality.value] * len(ids),
            "origin": [origin.value] * len(ids),
            "participant_id": ids,
            "tissue": [tissue] * len(ids),
        }))


def carrier_dataset(project: Project, source_id: str) -> tuple[Artifact, Path, Path] | None:
    """(artifact, carrier dataset root, variants.parquet) of the latest participant-variant table."""
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).latest(source_id=source_id, artifact_type="participant_variants")
    if art is None:
        return None
    root = Path(art.path)
    return art, root / "participant_variants", root / "variants.parquet"


def _limit(project: Project) -> int:
    """Same full-checksum size limit as snapshot verification (larger files: labelled sampled hash)."""
    return int(project.config.storage.full_checksum_limit_gb * 1024**3)
