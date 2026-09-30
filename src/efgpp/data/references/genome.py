"""Reference genome FASTA (UCSC analysis sets) and liftOver chain files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from efgpp.data.genotype.build import build_fai
from efgpp.data.references.base import FetchedResource, ReferenceProvider
from efgpp.resources.downloader import download, gunzip

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath"
FASTA_URLS = {"GRCh37": f"{UCSC}/hg19/bigZips/hg19.fa.gz", "GRCh38": f"{UCSC}/hg38/bigZips/hg38.fa.gz"}
CHAIN_URLS = {
    "GRCh37->GRCh38": f"{UCSC}/hg19/liftOver/hg19ToHg38.over.chain.gz",
    "GRCh38->GRCh37": f"{UCSC}/hg38/liftOver/hg38ToHg19.over.chain.gz",
}


def installed_fasta(project: Any) -> Path | None:
    """The project's reference FASTA: resources.genome.fasta, else the installed genome resource."""
    cfg = project.resources.genome
    if cfg.fasta:
        p = project.resolve(cfg.fasta)
        return p if p.exists() else None
    from efgpp.data.registry import Registry

    with Registry.open(project) as reg:
        row = reg.one("SELECT local_path, genome_build FROM resources WHERE name = 'genome' "
                      "ORDER BY download_date DESC LIMIT 1")
    if row and row["local_path"] and row["genome_build"] in (None, project.config.defaults.target_build):
        p = project.resolve(row["local_path"])
        if p.is_file():
            return p
    return None


class GenomeProvider(ReferenceProvider):
    name = "genome"
    license = "UCSC Genome Browser data: free for academic, non-profit and commercial use"
    homepage = "https://hgdownload.soe.ucsc.edu/"

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "builds": list(FASTA_URLS), "source": UCSC,
                "provides": "uncompressed FASTA + .fai (build checks, VEP --fasta) and chain files"}

    def version(self) -> str | None:
        return "UCSC"

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        build = build or self.project.resources.genome.build
        if build not in FASTA_URLS:
            raise ValueError(f"unsupported build {build!r}")
        gz = staging / Path(FASTA_URLS[build]).name
        digest = download(FASTA_URLS[build], gz, client=self.client)
        fa = gunzip(gz, staging / gz.name.removesuffix(".gz"))
        gz.unlink()
        fai = build_fai(fa)
        files = [fa, fai]
        for key, url in CHAIN_URLS.items():
            if key.startswith(build):
                chain = staging / Path(url).name
                download(url, chain, client=self.client)
                files.append(chain)
        return FetchedResource(self.name, f"{build}-ucsc", files, fa, FASTA_URLS[build], self.license,
                               genome_build=build, metadata={"download_sha256": digest})

    def validate(self, path: Path) -> list[str]:
        problems = super().validate(path)
        if not problems and not Path(str(path) + ".fai").exists():
            problems.append("FASTA index (.fai) missing")
        return problems
