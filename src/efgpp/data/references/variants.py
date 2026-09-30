"""Variant-level reference providers: ClinVar, gnomAD, AlphaMissense, AlphaGenome."""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Any

from efgpp.data.references.base import FetchedResource, ReferenceProvider
from efgpp.resources.downloader import download

CLINVAR = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
ALPHAMISSENSE = "https://storage.googleapis.com/dm_alphamissense"


def vcf_file_date(path: Path) -> str | None:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            m = re.match(r"##fileDate=(\d{8})", line)
            if m:
                return m.group(1)
    return None


class ClinVarProvider(ReferenceProvider):
    name = "clinvar"
    license = "Public domain (NCBI); cite ClinVar"
    homepage = "https://www.ncbi.nlm.nih.gov/clinvar/"

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "source": CLINVAR, "provides": "clinical significance per variant (VCF)"}

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        build = build or self.project.resources.genome.build
        url = f"{CLINVAR}/vcf_{build}/clinvar.vcf.gz"
        dest = staging / "clinvar.vcf.gz"
        digest = download(url, dest, client=self.client)
        tbi = staging / "clinvar.vcf.gz.tbi"
        download(url + ".tbi", tbi, client=self.client)
        version = vcf_file_date(dest) or "unknown"
        return FetchedResource(self.name, f"{version}_{build}", [dest, tbi], dest, url, self.license,
                               genome_build=build, release_date=version, metadata={"download_sha256": digest})

    def version(self) -> str | None:
        r = self.client.head(f"{CLINVAR}/vcf_GRCh38/clinvar.vcf.gz")
        return r.headers.get("last-modified")


class GnomADProvider(ReferenceProvider):
    name = "gnomad"
    license = "ODbL 1.0 / CC0 (see gnomAD terms)"
    homepage = "https://gnomad.broadinstitute.org/downloads"

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "source": self.homepage,
                "provides": "population allele frequencies",
                "note": "releases are terabyte-scale; set resources.yaml gnomad.path to an existing sites VCF"}

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        cfg = self.project.resources.gnomad
        if not cfg.url:
            raise RuntimeError("gnomAD is too large to fetch wholesale; set gnomad.path (local VCF) or gnomad.url "
                               "(one sites VCF, e.g. a single chromosome) in resources.yaml")
        dest = staging / Path(cfg.url).name
        digest = download(cfg.url, dest, client=self.client)
        return FetchedResource(self.name, cfg.version or "custom", [dest], dest, cfg.url, self.license,
                               genome_build=build, metadata={"download_sha256": digest})


class AlphaMissenseProvider(ReferenceProvider):
    name = "alphamissense"
    license = "CC BY 4.0 (predictions); see AlphaMissense terms"
    homepage = "https://github.com/google-deepmind/alphamissense"

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "source": ALPHAMISSENSE,
                "provides": "precomputed missense pathogenicity scores (reference score source, not trained here)"}

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        build = build or self.project.resources.genome.build
        code = {"GRCh37": "hg19", "GRCh38": "hg38"}[build]
        url = self.project.resources.alphamissense.url or f"{ALPHAMISSENSE}/AlphaMissense_{code}.tsv.gz"
        dest = staging / Path(url).name
        digest = download(url, dest, client=self.client)
        version = self.project.resources.alphamissense.version or "2023"
        return FetchedResource(self.name, f"{version}_{build}", [dest], dest, url, self.license,
                               genome_build=build, metadata={"download_sha256": digest})


class AlphaGenomeProvider(ReferenceProvider):
    name = "alphagenome"
    license = "AlphaGenome API terms of use"
    homepage = "https://github.com/google-deepmind/alphagenome"
    bulk = False

    def describe(self) -> dict[str, Any]:
        cfg = self.project.resources.alphagenome
        from efgpp.data.annotation.alphagenome import api_key

        try:
            import alphagenome  # noqa: F401

            client = True
        except ImportError:
            client = False
        return {"name": self.name, "mode": cfg.mode, "atlas_path": cfg.atlas_path, "client_installed": client,
                "api_key_configured": bool(api_key(self.project)),
                "note": "atlas mode needs no API key; API mode is for targeted variant sets"}

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        cfg = self.project.resources.alphagenome
        if not cfg.url:
            raise RuntimeError("set alphagenome.url to a precomputed atlas table, or alphagenome.atlas_path to a local copy")
        dest = staging / Path(cfg.url).name
        digest = download(cfg.url, dest, client=self.client)
        return FetchedResource(self.name, cfg.version or "atlas", [dest], dest, cfg.url, self.license,
                               genome_build=build, metadata={"download_sha256": digest})
