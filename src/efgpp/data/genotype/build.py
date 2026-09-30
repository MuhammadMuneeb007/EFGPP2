"""Genome build inference (GRCh37 vs GRCh38) from several independent lines of evidence.

Evidence, strongest first:
  1. reference FASTA concordance - REF alleles of sampled variants vs each build's FASTA
  2. user metadata              - an explicit genome_build in data.yaml
  3. VCF header                 - ##reference / ##contig lengths / ##assembly
  4. coordinate bounds          - positions beyond a chromosome's length in one build
  5. marker panel               - known rsID positions (optional resource file)

If the evidence is weak or contradictory the build stays "unknown" and validation asks
the user to set `genome_build` explicitly. Nothing is ever lifted over silently.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from efgpp.constants import GenomeBuild

# Primary-assembly chromosome lengths (UCSC hg19 / hg38 chrom.sizes).
CHROM_LENGTHS: dict[GenomeBuild, dict[str, int]] = {
    GenomeBuild.GRCH37: {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260,
        "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747,
        "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392,
        "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
    },
    GenomeBuild.GRCH38: {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259,
        "6": 170805979, "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422,
        "11": 135086622, "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189,
        "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
    },
}

BUILDS = (GenomeBuild.GRCH37, GenomeBuild.GRCH38)


@dataclass
class Evidence:
    method: str
    supports: GenomeBuild | None
    weight: float  # 0..1
    detail: str


@dataclass
class BuildInference:
    build: GenomeBuild
    confidence: float
    evidence: list[Evidence] = field(default_factory=list)
    requires_user: bool = False

    @property
    def confident(self) -> bool:
        return self.build in BUILDS and not self.requires_user

    def to_dict(self) -> dict[str, object]:
        return {
            "build": self.build.value,
            "confidence": round(self.confidence, 3),
            "requires_user": self.requires_user,
            "evidence": [
                {"method": e.method, "supports": e.supports.value if e.supports else None,
                 "weight": e.weight, "detail": e.detail}
                for e in self.evidence
            ],
        }


def normalize_chrom(c: str) -> str:
    c = re.sub(r"^chr", "", str(c), flags=re.IGNORECASE)
    return {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}.get(c.upper(), c.upper())


def evidence_from_vcf_header(meta_lines: list[str]) -> list[Evidence]:
    out: list[Evidence] = []
    for line in meta_lines:
        low = line.lower()
        if low.startswith(("##reference=", "##assembly=")):
            text = low
            if any(k in text for k in ("grch38", "hg38", "hs38")):
                out.append(Evidence("vcf_header_reference", GenomeBuild.GRCH38, 0.6, line[:120]))
            elif any(k in text for k in ("grch37", "hg19", "hs37", "b37", "human_g1k_v37")):
                out.append(Evidence("vcf_header_reference", GenomeBuild.GRCH37, 0.6, line[:120]))
        m = re.match(r"##contig=<ID=([^,>]+).*length=(\d+)", line)
        if m:
            chrom, length = normalize_chrom(m.group(1)), int(m.group(2))
            for build in BUILDS:
                if CHROM_LENGTHS[build].get(chrom) == length:
                    out.append(Evidence("vcf_header_contig_length", build, 0.9, f"chr{chrom} length {length}"))
                    break
    # Collapse many contig lines into one piece of evidence per build.
    collapsed: dict[tuple[str, GenomeBuild | None], Evidence] = {}
    for e in out:
        key = (e.method, e.supports)
        if key in collapsed:
            collapsed[key].detail = f"{collapsed[key].detail}; {e.detail}"[:240]
        else:
            collapsed[key] = e
    return list(collapsed.values())


def evidence_from_positions(variants: pl.DataFrame) -> list[Evidence]:
    """Positions past a chromosome's end rule a build out (weak, but never wrong)."""
    if variants.height == 0:
        return []
    maxpos = (
        variants.with_columns(pl.col("chromosome").map_elements(normalize_chrom, return_dtype=pl.Utf8))
        .group_by("chromosome").agg(pl.col("position").max())
    )
    excluded: dict[GenomeBuild, list[str]] = {b: [] for b in BUILDS}
    for chrom, pos in maxpos.iter_rows():
        for build in BUILDS:
            length = CHROM_LENGTHS[build].get(chrom)
            if length is not None and pos > length:
                excluded[build].append(f"chr{chrom}:{pos}")
    out = []
    for build in BUILDS:
        other = GenomeBuild.GRCH38 if build == GenomeBuild.GRCH37 else GenomeBuild.GRCH37
        if excluded[build] and not excluded[other]:
            out.append(Evidence("coordinate_bounds", other, 0.7,
                                f"positions beyond {build.value} chromosome ends: {excluded[build][:3]}"))
    return out


class FastaIndex:
    """Random access to an uncompressed FASTA through its .fai index (samtools format)."""

    def __init__(self, fasta: Path) -> None:
        self.fasta = fasta
        fai = Path(str(fasta) + ".fai")
        if not fai.exists():
            build_fai(fasta)
        self.index: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            name, length, offset, bases, width = line.split("\t")[:5]
            self.index[normalize_chrom(name)] = (int(length), int(offset), int(bases), int(width))

    def base(self, chrom: str, pos: int, n: int = 1) -> str | None:
        entry = self.index.get(normalize_chrom(chrom))
        if entry is None:
            return None
        length, offset, bases, width = entry
        if pos < 1 or pos + n - 1 > length:
            return None
        out = []
        with open(self.fasta, "rb") as fh:
            for p in range(pos - 1, pos - 1 + n):
                fh.seek(offset + (p // bases) * width + p % bases)
                out.append(fh.read(1).decode())
        return "".join(out).upper()


def build_fai(fasta: Path) -> Path:
    """Create a samtools-compatible .fai for an uncompressed FASTA."""
    if fasta.name.endswith(".gz"):
        raise ValueError("FASTA must be uncompressed (or bgzipped with a samtools .fai)")
    entries = []
    with open(fasta, "rb") as fh:
        name = None
        length = offset = bases = width = 0
        pos = 0
        for raw in fh:
            if raw.startswith(b">"):
                if name is not None:
                    entries.append((name, length, offset, bases, width))
                name = raw[1:].split()[0].decode()
                length = bases = width = 0
                offset = pos + len(raw)
            else:
                seq = raw.rstrip(b"\r\n")
                if bases == 0:
                    bases, width = len(seq), len(raw)
                length += len(seq)
            pos += len(raw)
        if name is not None:
            entries.append((name, length, offset, bases, width))
    fai = Path(str(fasta) + ".fai")
    fai.write_text("".join(f"{n}\t{size}\t{o}\t{b}\t{w}\n" for n, size, o, b, w in entries), encoding="utf-8")
    return fai


def evidence_from_fasta(variants: pl.DataFrame, fastas: dict[GenomeBuild, Path], n: int = 2000) -> list[Evidence]:
    snvs = variants.filter(
        (pl.col("reference").str.len_chars() == 1) & (pl.col("alternate").str.len_chars() == 1)
        & pl.col("reference").is_in(["A", "C", "G", "T"])
    )
    if snvs.height == 0:
        return []
    sample = snvs.sample(n=min(n, snvs.height), seed=0)
    rates: dict[GenomeBuild, float] = {}
    for build, fasta in fastas.items():
        idx = FastaIndex(fasta)
        hit = tested = 0
        for chrom, pos, ref, alt in sample.select("chromosome", "position", "reference", "alternate").iter_rows():
            b = idx.base(chrom, pos)
            if b is None or b == "N":
                continue
            tested += 1
            # PLINK .bim files may carry REF/ALT swapped, so either allele matching counts.
            hit += b in (ref, alt)
        if tested:
            rates[build] = hit / tested
    if not rates:
        return []
    best = max(rates, key=lambda b: rates[b])
    others = [r for b, r in rates.items() if b != best]
    margin = rates[best] - (max(others) if others else 0.5)
    detail = ", ".join(f"{b.value}: {r:.1%}" for b, r in rates.items())
    if rates[best] >= 0.95 and margin >= 0.1:
        return [Evidence("reference_fasta", best, 1.0, f"allele concordance {detail}")]
    return [Evidence("reference_fasta", None, 0.0, f"inconclusive allele concordance {detail}")]


def evidence_from_marker_panel(variants: pl.DataFrame, panel: Path) -> list[Evidence]:
    """Panel TSV columns: rsid, chromosome, pos_grch37, pos_grch38."""
    ref = pl.read_csv(panel, separator="\t", schema_overrides={"chromosome": pl.Utf8})
    joined = variants.join(ref, left_on="variant_id", right_on="rsid", how="inner")
    if joined.height < 20:
        return []
    m37 = int((joined.get_column("position") == joined.get_column("pos_grch37")).sum())
    m38 = int((joined.get_column("position") == joined.get_column("pos_grch38")).sum())
    total = joined.height
    if max(m37, m38) / total >= 0.9:
        build = GenomeBuild.GRCH37 if m37 > m38 else GenomeBuild.GRCH38
        return [Evidence("marker_panel", build, 0.95, f"{max(m37, m38)}/{total} rsIDs at {build.value} positions")]
    return [Evidence("marker_panel", None, 0.0, f"37: {m37}, 38: {m38} of {total}")]


def infer_build(
    *,
    declared: str | None,
    variants: pl.DataFrame | None = None,
    vcf_meta: list[str] | None = None,
    fastas: dict[GenomeBuild, Path] | None = None,
    marker_panel: Path | None = None,
) -> BuildInference:
    evidence: list[Evidence] = []
    decl = GenomeBuild.normalize(declared) if declared else GenomeBuild.AUTO
    if decl in BUILDS:
        evidence.append(Evidence("user_metadata", decl, 0.8, f"declared {declared}"))
    if vcf_meta:
        evidence += evidence_from_vcf_header(vcf_meta)
    if variants is not None:
        evidence += evidence_from_positions(variants)
        if fastas:
            evidence += evidence_from_fasta(variants, fastas)
        if marker_panel and marker_panel.exists():
            evidence += evidence_from_marker_panel(variants, marker_panel)

    score = dict.fromkeys(BUILDS, 0.0)
    for e in evidence:
        if e.supports in score:
            score[e.supports] += e.weight  # type: ignore[index]
    total = sum(score.values())
    if total == 0:
        return BuildInference(GenomeBuild.UNKNOWN, 0.0, evidence, requires_user=True)
    best = max(score, key=lambda b: score[b])
    confidence = score[best] / total
    contradiction = all(v > 0 for v in score.values())
    strong = score[best] >= 0.8
    if decl in BUILDS and decl != best and score[best] >= 1.0:
        # Data contradict the declared build: stop rather than trust either silently.
        return BuildInference(GenomeBuild.UNKNOWN, confidence, evidence, requires_user=True)
    requires_user = not strong or (contradiction and confidence < 0.8)
    return BuildInference(best if not requires_user else GenomeBuild.UNKNOWN, confidence, evidence, requires_user)

