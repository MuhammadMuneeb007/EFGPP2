"""data.yaml: every participant-level source in the cohort.

Phenotypes are plain configuration objects; no phenotype is known to the code.
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import Field, field_validator, model_validator

from efgpp.config.models import StrictModel
from efgpp.constants import (
    DEFAULT_MISSING_TOKENS,
    Modality,
    Origin,
    PhenotypeType,
    StorageMode,
    TemporalType,
)

ID_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$"


class OntologyTerm(StrictModel):
    id: str
    label: str | None = None


class TimelineSpec(StrictModel):
    """Columns in a source table that place each row on the participant timeline."""

    event_column: str | None = None
    time_column: str | None = None  # date or datetime of the event/collection
    study_day_column: str | None = None
    visit_name_column: str | None = None
    age_column: str | None = None


class SampleIdSpec(StrictModel):
    column: Literal["IID", "FID"] = "IID"
    # iid: use IID (or FID); fid_iid: join FID and IID with `separator`.
    mode: Literal["iid", "fid_iid"] = "iid"
    separator: str = "_"


class SourceBase(StrictModel):
    id: str = Field(pattern=ID_PATTERN)
    path: str
    format: str = "auto"
    mode: StorageMode = StorageMode.AUTO
    origin: Origin = Origin.OBSERVED
    description: str | None = None

    @field_validator("origin")
    @classmethod
    def _participant_level_origin(cls, v: Origin) -> Origin:
        if v == Origin.REFERENCE:
            raise ValueError("participant-level sources cannot have origin 'reference'")
        return v


class GenotypeSource(SourceBase):
    genome_build: str = "auto"
    sample_id: SampleIdSpec = SampleIdSpec()


class PhenotypeSource(SourceBase):
    name: str
    participant_id_column: str = "participant_id"
    value_column: str
    type: PhenotypeType
    # multiclass: allowed classes; ordinal: classes in increasing order.
    levels: list[str] | None = None
    # Binary coding. When omitted {0,1} is used, or PLINK {1,2} is detected with a warning.
    case_values: list[str] | None = None
    control_values: list[str] | None = None
    missing_values: list[str] = list(DEFAULT_MISSING_TOKENS)
    units: str | None = None
    ontology_term: OntologyTerm | None = None
    timeline: TimelineSpec | None = None

    @model_validator(mode="after")
    def _coding_pairs(self) -> PhenotypeSource:
        if (self.case_values is None) != (self.control_values is None):
            raise ValueError("case_values and control_values must be given together")
        if self.type == PhenotypeType.ORDINAL and not self.levels:
            raise ValueError("ordinal phenotypes need `levels` in increasing order")
        return self


class CovariateSource(SourceBase):
    participant_id_column: str = "participant_id"
    variables: list[str] = Field(min_length=1)
    categorical: list[str] = []
    missing_values: list[str] = list(DEFAULT_MISSING_TOKENS)
    timeline: TimelineSpec | None = None


class OmicsSource(SourceBase):
    participant_id_column: str = "participant_id"
    tissue: str | None = None
    feature_id_system: str | None = None  # e.g. ensembl_gene, illumina_probe, uniprot, hmdb
    measurement_type: str | None = None  # e.g. counts, tpm, beta, m_value, npx, intensity
    normalization: str | None = None  # e.g. raw, log2_tpm, quantile
    platform: str | None = None
    units: str | None = None
    genome_build: str | None = None
    batch_column: str | None = None
    biospecimen_column: str | None = None
    # Tabular layouts: rows are samples (default) or rows are features.
    orientation: Literal["samples_by_features", "features_by_samples"] = "samples_by_features"
    feature_id_column: str | None = None  # required for features_by_samples
    timeline: TimelineSpec | None = None
    temporal_type: TemporalType | None = None

    @model_validator(mode="after")
    def _orientation(self) -> OmicsSource:
        if self.orientation == "features_by_samples" and not self.feature_id_column:
            raise ValueError("features_by_samples tables need `feature_id_column`")
        return self


class ClinicalSource(SourceBase):
    participant_id_column: str = "participant_id"
    layout: Literal["wide", "long"] = "wide"
    variables: list[str] | None = None  # wide layout; None = every non-id column
    variable_column: str | None = None  # long layout
    value_column: str | None = None
    unit_column: str | None = None
    missing_values: list[str] = list(DEFAULT_MISSING_TOKENS)
    timeline: TimelineSpec | None = None

    @model_validator(mode="after")
    def _long_layout(self) -> ClinicalSource:
        if self.layout == "long" and not (self.variable_column and self.value_column):
            raise ValueError("long clinical tables need `variable_column` and `value_column`")
        return self


class ObservedSources(StrictModel):
    genotype: list[GenotypeSource] = []
    phenotypes: list[PhenotypeSource] = []
    covariates: list[CovariateSource] = []
    expression: list[OmicsSource] = []
    methylation: list[OmicsSource] = []
    proteomics: list[OmicsSource] = []
    metabolomics: list[OmicsSource] = []
    clinical: list[ClinicalSource] = []

    @field_validator("phenotypes", mode="before")
    @classmethod
    def _phenotypes_from_mapping(cls, value: Any) -> Any:
        """Accept the compact `{PH001: {name, column, type}}` form as well as a list."""
        if isinstance(value, dict):
            items = []
            for key, spec in value.items():
                spec = dict(spec or {})
                spec.setdefault("id", key)
                if "column" in spec and "value_column" not in spec:
                    spec["value_column"] = spec.pop("column")
                items.append(spec)
            return items
        return value

    def omics(self) -> dict[Modality, list[OmicsSource]]:
        return {
            Modality.EXPRESSION: self.expression,
            Modality.METHYLATION: self.methylation,
            Modality.PROTEOMICS: self.proteomics,
            Modality.METABOLOMICS: self.metabolomics,
        }


class ModelProviderConfig(StrictModel):
    type: Literal["predictdb", "custom"] = "predictdb"
    models_dir: str | None = None
    # PredictDB naming, e.g. prefix "mashr_" -> mashr_Whole_Blood.db
    model_prefix: str = "mashr_"
    model_suffix: str = ".db"
    # Explicit per-tissue model databases override models_dir.
    model_paths: dict[str, str] = {}
    version: str | None = None


class PredictedModalityConfig(StrictModel):
    enabled: bool = False
    engine: Literal["metaxcan"] = "metaxcan"
    genotype_artifact: str | None = None  # a genotype source id, e.g. GENO001
    tissues: list[str] = []
    model_provider: ModelProviderConfig = ModelProviderConfig()
    # Use the QC-passed genotype when one exists.
    use_qc_genotype: bool = True
    # Genome build of the prediction models; genotypes must match (EFGPP never lifts silently).
    model_genome_build: Literal["GRCh37", "GRCh38"] = "GRCh38"
    # MetaXcan --on_the_fly_mapping METADATA pattern (GTEx v8 PredictDB variant ids).
    variant_id_pattern: str | None = "chr{}_{}_{}_{}_b38"
    extra_args: list[str] = []


class PredictedConfig(StrictModel):
    expression: PredictedModalityConfig = PredictedModalityConfig()
    proteomics: PredictedModalityConfig = PredictedModalityConfig()
    metabolomics: PredictedModalityConfig = PredictedModalityConfig()

    def items(self) -> list[tuple[Modality, PredictedModalityConfig]]:
        return [
            (Modality.EXPRESSION, self.expression),
            (Modality.PROTEOMICS, self.proteomics),
            (Modality.METABOLOMICS, self.metabolomics),
        ]


class AliasFileSpec(StrictModel):
    """Maps native identifiers used by some sources to canonical participant IDs."""

    path: str
    participant_id_column: str = "participant_id"
    alias_column: str
    # Source ids this mapping applies to; ["*"] means every source.
    sources: list[str] = ["*"]


class ParticipantsSpec(StrictModel):
    canonical_id: str = "participant_id"
    # Optional master participant list (id + static attributes).
    path: str | None = None
    id_column: str = "participant_id"
    aliases: list[AliasFileSpec] = []


class EventsSpec(StrictModel):
    path: str
    participant_id_column: str = "participant_id"
    event_id_column: str = "event_id"
    visit_name_column: str | None = None
    date_column: str | None = None
    study_day_column: str | None = None
    age_column: str | None = None


class BiospecimensSpec(StrictModel):
    path: str
    biospecimen_id_column: str = "biospecimen_id"
    participant_id_column: str = "participant_id"
    event_id_column: str | None = None
    tissue_column: str | None = None
    material_column: str | None = None
    collection_date_column: str | None = None
    processing_method_column: str | None = None
    storage_condition_column: str | None = None


class SimGenotypeConfig(StrictModel):
    enabled: bool = False
    engine: Literal["builtin", "stdpopsim"] = "builtin"
    n_participants: int = Field(500, ge=2)
    n_variants: int = Field(2000, ge=1)
    n_chromosomes: int = Field(2, ge=1, le=22)
    seed: int = 1


class SimPhenotypeConfig(StrictModel):
    enabled: bool = False
    count: int = Field(1, ge=1)
    type: PhenotypeType = PhenotypeType.BINARY
    heritability: float = Field(0.5, ge=0, le=1)
    n_causal: int = Field(20, ge=1)
    prevalence: float = Field(0.2, gt=0, lt=1)
    seed: int = 2


class SimExpressionConfig(StrictModel):
    enabled: bool = False
    n_features: int = Field(50, ge=1)
    seed: int = 3


class SimCovariatesConfig(StrictModel):
    enabled: bool = False
    seed: int = 4


class SimulationConfig(StrictModel):
    genotype: SimGenotypeConfig = SimGenotypeConfig()
    phenotype: SimPhenotypeConfig = SimPhenotypeConfig()
    expression: SimExpressionConfig = SimExpressionConfig()
    covariates: SimCovariatesConfig = SimCovariatesConfig()


class DataConfig(StrictModel):
    participants: ParticipantsSpec = ParticipantsSpec()
    observed: ObservedSources = ObservedSources()
    predicted: PredictedConfig = PredictedConfig()
    simulation: SimulationConfig = SimulationConfig()
    events: EventsSpec | None = None
    biospecimens: BiospecimensSpec | None = None

    @model_validator(mode="after")
    def _unique_ids(self) -> DataConfig:
        seen: dict[str, str] = {}
        for modality, source in self.iter_sources():
            if source.id in seen:
                raise ValueError(
                    f"source id {source.id!r} is used by both {seen[source.id]} and {modality}"
                )
            seen[source.id] = str(modality)
        names = [p.name for p in self.observed.phenotypes]
        dupes = sorted({n for n in names if names.count(n) > 1})
        if dupes:
            raise ValueError(f"duplicate phenotype names: {dupes}")
        return self

    def iter_sources(self) -> list[tuple[Modality, SourceBase]]:
        o = self.observed
        out: list[tuple[Modality, SourceBase]] = []
        out += [(Modality.GENOTYPE, s) for s in o.genotype]
        out += [(Modality.PHENOTYPE, s) for s in o.phenotypes]
        out += [(Modality.COVARIATES, s) for s in o.covariates]
        for modality, sources in o.omics().items():
            out += [(modality, s) for s in sources]
        out += [(Modality.CLINICAL, s) for s in o.clinical]
        return out

    def get_source(self, source_id: str) -> tuple[Modality, SourceBase]:
        for modality, source in self.iter_sources():
            if source.id == source_id:
                return modality, source
        raise KeyError(f"unknown source id {source_id!r}")

    def phenotype(self, key: str) -> PhenotypeSource:
        """Look a phenotype up by id or by name."""
        for p in self.observed.phenotypes:
            if key in (p.id, p.name):
                return p
        raise KeyError(f"unknown phenotype {key!r}")

    def next_id(self, prefix: str) -> str:
        used = {s.id for _, s in self.iter_sources()}
        n = 1
        while f"{prefix}{n:03d}" in used:
            n += 1
        return f"{prefix}{n:03d}"
