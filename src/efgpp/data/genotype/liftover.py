"""Genome build checks and coordinate liftover with pyliftover.

Policy: every EFGPP project has one target build (project.yaml `defaults.target_build`,
GRCh38). Genotypes, reference panels and GWAS summary statistics in another build are lifted
to it, always into new artifacts; originals are never modified.

Chain files (UCSC hg19ToHg38 / hg38ToHg19) are downloaded once into the project's visible
`resources/liftover/` folder and passed explicitly to pyliftover (and to GWASLab), so nothing
is downloaded into the home directory.

Build check ("reference bases"): for a sample of SNVs, the reference base at each position is
fetched for both GRCh37 and GRCh38 (Ensembl REST, batch requests) and compared with the
variant's alleles; the positions are also lifted with pyliftover and checked against the other
build's reference. Data in build B match B's reference for ~100% of SNVs (either allele, so
PLINK's A1/A2 order does not matter) and the other build for ~40-60%.
"""

from __future__ import annotations

import os
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import httpx
import polars as pl

from efgpp.constants import GenomeBuild

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath"
CHAINS = {
    (GenomeBuild.GRCH37, GenomeBuild.GRCH38): f"{UCSC}/hg19/liftOver/hg19ToHg38.over.chain.gz",
    (GenomeBuild.GRCH38, GenomeBuild.GRCH37): f"{UCSC}/hg38/liftOver/hg38ToHg19.over.chain.gz",
}
ENSEMBL = {GenomeBuild.GRCH38: "https://rest.ensembl.org", GenomeBuild.GRCH37: "https://grch37.rest.ensembl.org"}
ENSEMBL_BATCH = 50  # regions per POST request (Ensembl limit)

# fetcher(build, [(chrom, pos1), ...]) -> {(chrom, pos1): base}
BaseFetcher = Callable[[GenomeBuild, list[tuple[str, int]]], dict[tuple[str, int], str]]


def offline() -> bool:
    return os.environ.get("EFGPP_OFFLINE", "").lower() in ("1", "true", "yes")


def chain_file(resource_root: Path, source: GenomeBuild, target: GenomeBuild) -> Path:
    """Local UCSC chain file under <project>/resources/liftover/, downloaded on first use."""
    from efgpp.resources.downloader import download

    url = CHAINS.get((source, target))
    if url is None:
        raise ValueError(f"no chain file for {source.value} -> {target.value}")
    path = resource_root / "liftover" / url.rsplit("/", 1)[-1]
    if not path.exists():
        if offline():
            raise RuntimeError(f"chain file {path.name} missing and EFGPP_OFFLINE is set")
        download(url, path)
    return path


def _ucsc_chrom(chrom: str) -> str:
    """'1' / 'chr1' / '23' / 'MT' -> UCSC names used by the chain files ('chr1', 'chrX', 'chrM')."""
    c = str(chrom)
    for prefix in ("chr", "CHR", "Chr"):
        c = c.removeprefix(prefix)
    c = {"23": "X", "24": "Y", "25": "X", "26": "M", "MT": "M"}.get(c, c)
    return f"chr{c}"


def _plain_chrom(chrom: str) -> str:
    c = chrom.removeprefix("chr")
    return "MT" if c == "M" else c


@dataclass
class Lifter:
    """pyliftover over a local chain file; 1-based positions in and out."""

    chain: Path

    def __post_init__(self) -> None:
        from pyliftover import LiftOver

        self._lo = LiftOver(str(self.chain))

    def lift(self, chrom: str, pos: int) -> tuple[str, int, str] | None:
        hits = self._lo.convert_coordinate(_ucsc_chrom(chrom), int(pos) - 1)
        if not hits:
            return None
        new_chrom, new_pos0, strand, _score = hits[0]
        return _plain_chrom(new_chrom), int(new_pos0) + 1, strand


def lift_table(variants: pl.DataFrame, lifter: Lifter) -> pl.DataFrame:
    """variant_id, chromosome, position -> new coordinates and a status per variant."""
    rows: list[tuple[str, str | None, int | None, str]] = []
    for vid, chrom, pos in variants.select("variant_id", "chromosome", "position").iter_rows():
        if pos is None or pos < 1:
            rows.append((vid, None, None, "unplaced"))
            continue
        hit = lifter.lift(chrom, pos)
        if hit is None:
            rows.append((vid, None, None, "unmapped"))
        elif hit[0] != _plain_chrom(_ucsc_chrom(chrom)):
            rows.append((vid, hit[0], hit[1], "chromosome_changed"))
        elif hit[2] == "-":
            rows.append((vid, hit[0], hit[1], "reverse_strand"))
        else:
            rows.append((vid, hit[0], hit[1], "mapped"))
    return pl.DataFrame(rows, schema={"variant_id": pl.Utf8, "new_chromosome": pl.Utf8,
                                      "new_position": pl.Int64, "status": pl.Utf8}, orient="row")


def ensembl_bases(build: GenomeBuild, sites: list[tuple[str, int]],
                  client: httpx.Client | None = None) -> dict[tuple[str, int], str]:
    """Reference base at each (chrom, pos1) from Ensembl REST (batched)."""
    own = client is None
    client = client or httpx.Client(timeout=60)
    out: dict[tuple[str, int], str] = {}
    try:
        for i in range(0, len(sites), ENSEMBL_BATCH):
            chunk = sites[i:i + ENSEMBL_BATCH]
            regions = [f"{_plain_chrom(_ucsc_chrom(c))}:{p}..{p}:1" for c, p in chunk]
            r = client.post(f"{ENSEMBL[build]}/sequence/region/human", json={"regions": regions},
                            headers={"Content-Type": "application/json", "Accept": "application/json"})
            r.raise_for_status()
            by_query = {item["query"]: item.get("seq", "").upper() for item in r.json()}
            for (c, p), q in zip(chunk, regions, strict=True):
                if by_query.get(q):
                    out[(c, p)] = by_query[q]
    finally:
        if own:
            client.close()
    return out


@dataclass
class ReferenceCheck:
    rate37: float
    rate38: float
    lifted_rate: float | None  # pyliftover-lifted positions vs the other build's reference
    tested: int
    lifted_direction: str | None


def check_reference_bases(variants: pl.DataFrame, resource_root: Path, *, n: int = 80,
                          fetch: BaseFetcher = ensembl_bases, lifter_factory: Callable[[Path], Lifter] = Lifter
                          ) -> ReferenceCheck | None:
    """Compare alleles with both references, and confirm with pyliftover. None if impossible."""
    snvs = variants.filter(
        pl.col("chromosome").cast(pl.Utf8).str.replace("^chr", "").is_in([str(i) for i in range(1, 23)])
        & (pl.col("position") > 0)
        & pl.col("reference").is_in(["A", "C", "G", "T"]) & pl.col("alternate").is_in(["A", "C", "G", "T"])
    )
    if snvs.height < 20:
        return None
    sample = snvs.sample(n=min(n, snvs.height), seed=11)
    sites = [(str(c), int(p)) for c, p in sample.select("chromosome", "position").iter_rows()]
    alleles = {(str(c), int(p)): {r, a} for c, p, r, a in
               sample.select("chromosome", "position", "reference", "alternate").iter_rows()}
    b37 = fetch(GenomeBuild.GRCH37, sites)
    b38 = fetch(GenomeBuild.GRCH38, sites)
    tested = [s for s in sites if s in b37 and s in b38]
    if len(tested) < 20:
        return None
    rate37 = sum(b37[s] in alleles[s] for s in tested) / len(tested)
    rate38 = sum(b38[s] in alleles[s] for s in tested) / len(tested)
    lifted_rate, direction = None, None
    likely = GenomeBuild.GRCH37 if rate37 >= rate38 else GenomeBuild.GRCH38
    other = GenomeBuild.GRCH38 if likely == GenomeBuild.GRCH37 else GenomeBuild.GRCH37
    try:
        lifter = lifter_factory(chain_file(resource_root, likely, other))
        lifted = {}
        for s in tested:
            hit = lifter.lift(*s)
            if hit and hit[2] == "+":
                lifted[s] = (hit[0], hit[1])
        if lifted:
            other_bases = fetch(other, list(lifted.values()))
            ok = [other_bases[t] in alleles[s] for s, t in lifted.items() if t in other_bases]
            if ok:
                lifted_rate = sum(ok) / len(ok)
                direction = f"{likely.value}->{other.value}"
    except (RuntimeError, ValueError, OSError, httpx.HTTPError):
        pass
    return ReferenceCheck(rate37, rate38, lifted_rate, len(tested), direction)

