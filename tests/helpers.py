"""Helpers that add sources to a test project exactly as `efgpp data add` would."""

from __future__ import annotations

from pathlib import Path

from efgpp.config.data import (
    CovariateSource,
    GenotypeSource,
    OmicsSource,
    PhenotypeSource,
    TimelineSpec,
)
from efgpp.constants import PhenotypeType, StorageMode
from efgpp.project import Project


def add_genotype(project: Project, cohort: Path, build: str = "GRCh38", sid: str = "GENO001") -> Project:
    project.data.observed.genotype.append(GenotypeSource(
        id=sid, path=str(cohort / "genotype" / "cohort"), format="bed", genome_build=build, mode=StorageMode.REFERENCE))
    project.save_data_config()
    return Project.load(project.root)


PHENOS = {
    "trait_a": (PhenotypeType.BINARY, None),
    "trait_b": (PhenotypeType.CONTINUOUS, None),
    "trait_c": (PhenotypeType.MULTICLASS, None),
    "trait_d": (PhenotypeType.ORDINAL, ["low", "mid", "high"]),
    "trait_e": (PhenotypeType.CONTINUOUS, None),
}


def add_phenotypes(project: Project, cohort: Path, names: list[str]) -> Project:
    for i, name in enumerate(names, start=1):
        ptype, levels = PHENOS[name]
        project.data.observed.phenotypes.append(PhenotypeSource(
            id=f"PH{i:03d}", name=name, path=str(cohort / "phenotypes.csv"), participant_id_column="IID",
            value_column=name, type=ptype, levels=levels))
    project.save_data_config()
    return Project.load(project.root)


def add_covariates(project: Project, cohort: Path) -> Project:
    project.data.observed.covariates.append(CovariateSource(
        id="COV001", path=str(cohort / "covariates.csv"), participant_id_column="IID",
        variables=["age", "sex", "bmi"], categorical=["sex"]))
    project.save_data_config()
    return Project.load(project.root)


def add_expression(project: Project, cohort: Path, longitudinal: bool = False) -> Project:
    if longitudinal:
        src = OmicsSource(id="RNA002", path=str(cohort / "expression_longitudinal.csv"), tissue="whole_blood",
                          biospecimen_column="biospecimen", feature_id_system="gene_symbol",
                          measurement_type="log2_tpm", normalization="log2",
                          timeline=TimelineSpec(event_column="visit_id", time_column="collection_date"))
    else:
        src = OmicsSource(id="RNA001", path=str(cohort / "expression.parquet"), tissue="whole_blood",
                          feature_id_system="ensembl_gene", measurement_type="log2_tpm", normalization="log2")
    project.data.observed.expression.append(src)
    project.save_data_config()
    return Project.load(project.root)
