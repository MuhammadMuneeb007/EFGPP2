"""Controlled vocabularies shared by every EFGPP data component."""

from __future__ import annotations

from enum import StrEnum


class Origin(StrEnum):
    """Where an artifact's values come from. Every artifact has exactly one origin."""

    OBSERVED = "observed"  # actually measured or supplied
    DERIVED = "derived"  # calculated from participant data, not claimed as measured
    PREDICTED = "predicted"  # estimated with an external predictive model
    SIMULATED = "simulated"  # synthetic participant-level data
    REFERENCE = "reference"  # knowledge that does not belong to one participant


class ArtifactStatus(StrEnum):
    REGISTERED = "REGISTERED"
    VALIDATED = "VALIDATED"
    STANDARDIZED = "STANDARDIZED"
    QC_PASS = "QC_PASS"
    QC_WARN = "QC_WARN"
    QC_FAIL = "QC_FAIL"
    READY = "READY"
    SUPERSEDED = "SUPERSEDED"


# Statuses whose artifacts may be handed to the Representation layer.
USABLE_STATUSES = frozenset(
    {
        ArtifactStatus.VALIDATED,
        ArtifactStatus.STANDARDIZED,
        ArtifactStatus.QC_PASS,
        ArtifactStatus.QC_WARN,
        ArtifactStatus.READY,
    }
)


class StorageMode(StrEnum):
    COPY = "copy"
    LINK = "link"
    REFERENCE = "reference"
    AUTO = "auto"


class Modality(StrEnum):
    # participant-level, observed/predicted/simulated
    GENOTYPE = "genotype"
    PHENOTYPE = "phenotype"
    COVARIATES = "covariates"
    EXPRESSION = "expression"
    METHYLATION = "methylation"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    SPLICING = "splicing"
    CLINICAL = "clinical"
    # phenotype-independent derivations
    GENOTYPE_QC = "genotype_qc"
    ANCESTRY = "ancestry"
    QC_PCA = "qc_pca"
    KINSHIP = "kinship"
    ROH = "roh"
    HLA = "hla"
    GENOTYPE_IMPUTATION = "genotype_imputation"
    VARIANTS = "variants"
    VARIANT_ANNOTATIONS = "variant_annotations"
    # participant-specific variants and their aggregation (DERIVED from the genotype)
    PARTICIPANT_VARIANTS = "participant_variants"
    CONSEQUENCE_BURDEN = "consequence_burden"
    GENE_BURDEN = "gene_burden"
    # simulation ground truth
    TRUTH = "truth"


OMICS_MODALITIES = (
    Modality.EXPRESSION,
    Modality.METHYLATION,
    Modality.PROTEOMICS,
    Modality.METABOLOMICS,
    Modality.SPLICING,
)

# Genetically predicted molecular modalities (origin PREDICTED; never "observed").
GENETICALLY_PREDICTED_MODALITIES = (
    Modality.EXPRESSION,
    Modality.SPLICING,
    Modality.PROTEOMICS,
    Modality.METABOLOMICS,
    Modality.METHYLATION,
)

# Modalities that can serve as participant-level predictors.
PREDICTIVE_MODALITIES = (
    Modality.GENOTYPE,
    Modality.COVARIATES,
    Modality.CLINICAL,
    *OMICS_MODALITIES,
)


class PhenotypeType(StrEnum):
    BINARY = "binary"
    CONTINUOUS = "continuous"
    MULTICLASS = "multiclass"
    ORDINAL = "ordinal"
    # architected for, not yet implemented
    SURVIVAL = "survival"
    LONGITUDINAL = "longitudinal"
    TIME_TO_EVENT = "time_to_event"
    REPEATED_MEASURES = "repeated_measures"


SUPPORTED_PHENOTYPE_TYPES = frozenset(
    {
        PhenotypeType.BINARY,
        PhenotypeType.CONTINUOUS,
        PhenotypeType.MULTICLASS,
        PhenotypeType.ORDINAL,
    }
)


class TemporalType(StrEnum):
    STATIC = "static"  # e.g. germline genotype
    EVENT = "event"  # tied to a visit / collection event
    GENETICALLY_PREDICTED_STATIC = "genetically_predicted_static"
    UNKNOWN = "unknown"


class GenomeBuild(StrEnum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"
    AUTO = "auto"
    UNKNOWN = "unknown"

    @classmethod
    def normalize(cls, value: str | None) -> GenomeBuild:
        if value is None:
            return cls.UNKNOWN
        v = str(value).strip().lower()
        aliases = {
            "grch37": cls.GRCH37,
            "hg19": cls.GRCH37,
            "b37": cls.GRCH37,
            "37": cls.GRCH37,
            "grch38": cls.GRCH38,
            "hg38": cls.GRCH38,
            "b38": cls.GRCH38,
            "38": cls.GRCH38,
            "auto": cls.AUTO,
        }
        return aliases.get(v, cls.UNKNOWN)


GENOTYPE_FORMATS = ("pgen", "bed", "bgen", "vcf", "bcf")
TABULAR_FORMATS = ("csv", "tsv", "txt", "parquet")
OMICS_FORMATS = (*TABULAR_FORMATS, "h5ad", "zarr")

# Default tokens treated as missing in user-supplied tables.
DEFAULT_MISSING_TOKENS = ("", "NA", "N/A", "NaN", "nan", "NULL", "null", "None", ".", "-9")

# Validation check outcomes.
PASS = "pass"
FAIL = "fail"
WARN = "warn"
INFO = "info"
