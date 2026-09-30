"""project.yaml: controls EFGPP behaviour for one cohort project."""

from __future__ import annotations

from typing import Literal

from pydantic import Field

from efgpp.config.models import StrictModel
from efgpp.constants import StorageMode


class ProjectInfo(StrictModel):
    name: str = "efgpp_project"
    root: str = "."


class SpeciesConfig(StrictModel):
    name: str = "homo_sapiens"


class DefaultsConfig(StrictModel):
    genome_build: str = "auto"


class StorageConfig(StrictModel):
    generated_data_root: str = "data"
    resource_root: str = "resources"
    work_root: str = "work"
    report_root: str = "reports"
    tabular_format: Literal["parquet"] = "parquet"
    compression: Literal["zstd", "snappy", "gzip", "lz4", "uncompressed"] = "zstd"
    omics_format: Literal["zarr", "parquet"] = "zarr"
    external_source_policy: StorageMode = StorageMode.AUTO
    copy_threshold_gb: float = Field(1.0, ge=0)
    # Files larger than this get a sampled SHA-256 (clearly labelled) instead of a full hash.
    full_checksum_limit_gb: float = Field(20.0, ge=0)


class RegistryConfig(StrictModel):
    database: str = "registry/efgpp.duckdb"


class ContainersConfig(StrictModel):
    enabled: bool = False
    linux_hpc: Literal["apptainer", "singularity", "docker", "none"] = "apptainer"
    desktop: Literal["docker", "podman", "none"] = "docker"


class HPCConfig(StrictModel):
    scheduler: Literal["auto", "slurm", "pbs", "lsf", "sge", "none"] = "auto"
    partition: str | None = None
    account: str | None = None
    default_runtime_min: int = 240
    default_mem_mb: int = 8000
    max_jobs: int = 50


class ExecutionConfig(StrictModel):
    # "snakemake" is preferred; "auto" falls back to the built-in executor when
    # Snakemake is not installed (e.g. native Windows).
    engine: Literal["snakemake", "builtin", "auto"] = "auto"
    local_cores: int = Field(8, ge=1)
    # pixi (default) | micromamba | mamba | conda | system (create nothing; use execution.tools/PATH)
    environment_manager: Literal["pixi", "micromamba", "mamba", "conda", "system"] = "pixi"
    # Used when pixi is unavailable or fails; falls through to any other conda-family tool found.
    fallback_environment_manager: Literal["micromamba", "mamba", "conda", "system"] = "mamba"
    containers: ContainersConfig = ContainersConfig()
    hpc: HPCConfig = HPCConfig()
    # Explicit executable overrides, e.g. {"plink2": "/opt/plink2/plink2"}.
    tools: dict[str, str] = {}


class PrivacyConfig(StrictModel):
    preserve_source_ids: bool = True
    allow_absolute_dates: bool = True


class GenotypeQCConfig(StrictModel):
    """Phenotype-independent genotype QC thresholds (PLINK2 semantics)."""

    sample_missingness: float = Field(0.02, ge=0, le=1)  # --mind
    variant_missingness: float = Field(0.02, ge=0, le=1)  # --geno
    maf: float = Field(0.01, ge=0, le=0.5)
    hwe_p: float = Field(1e-6, ge=0, le=1)
    heterozygosity_sd: float = Field(3.0, gt=0)
    kinship_cutoff: float = Field(0.0884, ge=0, le=0.5)  # 2nd-degree relatives
    duplicate_kinship: float = Field(0.354, ge=0, le=0.5)  # duplicates / MZ twins
    sex_check: bool = True
    # Related pairs are always reported; removing one member is an analysis decision.
    remove_related: bool = False
    remove_duplicates: bool = True
    ld_window: int = 200
    ld_step: int = 50
    ld_r2: float = 0.2
    # Fraction of samples failing QC above which the artifact is QC_WARN / QC_FAIL.
    warn_sample_fail_fraction: float = 0.05
    fail_sample_fail_fraction: float = 0.5


class PCAConfig(StrictModel):
    qc_engine: Literal["plink2", "flashpca2"] = "plink2"
    optional_engines: list[Literal["flashpca2"]] = []
    n_pcs: int = Field(10, ge=1, le=100)
    # PLINK2 randomized PCA; "auto" enables it above 5,000 samples.
    approx: Literal["auto", "always", "never"] = "auto"


class RelatednessConfig(StrictModel):
    enabled: bool = True


class ROHConfig(StrictModel):
    # ROH requires PLINK 1.9 (--homozyg is not implemented in PLINK 2).
    enabled: bool = False


class GenotypeProcessingConfig(StrictModel):
    qc: GenotypeQCConfig = GenotypeQCConfig()
    pca: PCAConfig = PCAConfig()
    relatedness: RelatednessConfig = RelatednessConfig()
    roh: ROHConfig = ROHConfig()


class ReportingConfig(StrictModel):
    plotly_js: Literal["inline", "cdn"] = "inline"
    multiqc: bool = True


class ProjectConfig(StrictModel):
    schema_version: int = 1
    project: ProjectInfo = ProjectInfo()
    species: SpeciesConfig = SpeciesConfig()
    defaults: DefaultsConfig = DefaultsConfig()
    storage: StorageConfig = StorageConfig()
    registry: RegistryConfig = RegistryConfig()
    execution: ExecutionConfig = ExecutionConfig()
    privacy: PrivacyConfig = PrivacyConfig()
    genotype: GenotypeProcessingConfig = GenotypeProcessingConfig()
    reporting: ReportingConfig = ReportingConfig()
