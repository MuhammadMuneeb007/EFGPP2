"""resources.yaml: external knowledge that does not belong to any participant."""

from __future__ import annotations

from typing import Literal

from efgpp.config.models import OpenModel, StrictModel


class ResourceEntry(OpenModel):
    """Common switches for every reference resource. Providers may add their own keys."""

    enabled: bool = False
    version: str | None = None
    # Use an existing local copy instead of downloading.
    path: str | None = None
    url: str | None = None


class GenomeResource(ResourceEntry):
    enabled: bool = True
    # Must match project.yaml defaults.target_build (GRCh38); resources are installed for it.
    build: Literal["GRCh37", "GRCh38"] = "GRCh38"
    auto_download: bool = True
    fasta: str | None = None


class VEPResource(ResourceEntry):
    cache: str = "auto"  # auto | path to an existing cache directory
    release: int | None = None
    species: str = "homo_sapiens"
    extra_args: list[str] = []


class OpenCravatResource(ResourceEntry):
    annotators: list[str] = []


class AlphaGenomeResource(ResourceEntry):
    mode: Literal["atlas", "api"] = "atlas"
    api_key_env: str = "ALPHAGENOME_API_KEY"
    # Precomputed atlas table (TSV/Parquet keyed by chrom/pos/ref/alt).
    atlas_path: str | None = None
    max_api_variants: int = 1000


class SpliceAIResource(ResourceEntry):
    # Maximum distance between the variant and gained/lost splice site (SpliceAI -D).
    distance: int = 50


class PGSCatalogResource(ResourceEntry):
    score_ids: list[str] = []
    genome_build: Literal["GRCh37", "GRCh38"] | None = None


class LiftOverResource(ResourceEntry):
    chain_files: dict[str, str] = {}  # e.g. {"GRCh37->GRCh38": "resources/.../hg19ToHg38.over.chain.gz"}


class ResourcesConfig(StrictModel):
    genome: GenomeResource = GenomeResource()
    vep: VEPResource = VEPResource()
    opencravat: OpenCravatResource = OpenCravatResource()
    clinvar: ResourceEntry = ResourceEntry()
    gnomad: ResourceEntry = ResourceEntry()
    dbsnp: ResourceEntry = ResourceEntry()
    alphamissense: ResourceEntry = ResourceEntry()
    alphagenome: AlphaGenomeResource = AlphaGenomeResource()
    spliceai: SpliceAIResource = SpliceAIResource()
    gtex: ResourceEntry = ResourceEntry()
    predictdb: ResourceEntry = ResourceEntry()
    omicspred: ResourceEntry = ResourceEntry()
    gwas_catalog: ResourceEntry = ResourceEntry()
    pgs_catalog: PGSCatalogResource = PGSCatalogResource()
    opentargets: ResourceEntry = ResourceEntry()
    ld_reference: ResourceEntry = ResourceEntry()
    liftover: LiftOverResource = LiftOverResource()

    def enabled_names(self) -> list[str]:
        return [name for name in type(self).model_fields if getattr(self, name).enabled]
