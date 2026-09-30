from __future__ import annotations

import pytest
import yaml
from pydantic import ValidationError

from efgpp.config import DataConfig, ProjectConfig, ResourcesConfig
from efgpp.constants import Origin, PhenotypeType, StorageMode

SPEC_DATA_YAML = """
participants:
  canonical_id: participant_id
observed:
  genotype:
    - id: GENO001
      path: /external/path/cohort
      format: pgen
      mode: reference
      genome_build: auto
      sample_id:
        column: IID
  phenotypes:
    - id: PH001
      name: trait_a
      path: inputs/phenotypes.parquet
      mode: copy
      participant_id_column: participant_id
      value_column: trait_a
      type: binary
    - id: PH002
      name: trait_b
      path: inputs/phenotypes.parquet
      participant_id_column: participant_id
      value_column: trait_b
      type: continuous
  covariates:
    - id: COV001
      path: inputs/covariates.parquet
      participant_id_column: participant_id
      variables: [age, sex, bmi]
  expression:
    - id: RNA001
      path: inputs/expression.zarr
      origin: observed
      tissue: whole_blood
      participant_id_column: participant_id
      timeline:
        event_column: visit_id
        time_column: collection_date
  methylation: []
  proteomics: []
  metabolomics: []
"""


def test_spec_data_yaml_parses() -> None:
    cfg = DataConfig.model_validate(yaml.safe_load(SPEC_DATA_YAML))
    assert [s.id for _, s in cfg.iter_sources()] == ["GENO001", "PH001", "PH002", "COV001", "RNA001"]
    assert cfg.observed.phenotypes[0].type == PhenotypeType.BINARY
    assert cfg.observed.phenotypes[0].mode == StorageMode.COPY
    assert cfg.observed.expression[0].timeline.event_column == "visit_id"  # type: ignore[union-attr]
    assert cfg.next_id("PH") == "PH003"


def test_compact_phenotype_mapping_form() -> None:
    cfg = DataConfig.model_validate({"observed": {"phenotypes": {
        "PH001": {"name": "trait_a", "column": "trait_a", "type": "binary", "path": "p.csv"},
        "PH002": {"name": "trait_b", "column": "trait_b", "type": "continuous", "path": "p.csv"},
    }}})
    assert [p.id for p in cfg.observed.phenotypes] == ["PH001", "PH002"]
    assert cfg.observed.phenotypes[1].value_column == "trait_b"


@pytest.mark.parametrize("bad", [
    {"observed": {"phenotypes": [{"id": "X", "name": "a", "path": "p", "value_column": "a", "type": "ordinal"}]}},
    {"observed": {"genotype": [{"id": "G", "path": "g", "origin": "reference"}]}},
    {"observed": {"genotype": [{"id": "G", "path": "g"}], "covariates": [{"id": "G", "path": "c", "variables": ["a"]}]}},
    {"observed": {"phenotypes": [{"id": "A", "name": "same", "path": "p", "value_column": "a", "type": "binary"},
                                 {"id": "B", "name": "same", "path": "p", "value_column": "b", "type": "binary"}]}},
    {"observed": {"genotype": [{"id": "G", "path": "g", "typo_field": 1}]}},
])
def test_invalid_configs_are_rejected(bad: dict) -> None:
    with pytest.raises(ValidationError):
        DataConfig.model_validate(bad)


def test_project_and_resources_defaults_roundtrip() -> None:
    p = ProjectConfig()
    assert p.storage.copy_threshold_gb == 1.0
    assert p.execution.environment_manager == "pixi"
    assert ProjectConfig.model_validate(p.to_yaml_dict()) == p
    r = ResourcesConfig.model_validate({"vep": {"enabled": True}, "alphagenome": {"enabled": True, "mode": "atlas"}})
    assert "vep" in r.enabled_names() and "alphagenome" in r.enabled_names()


def test_future_phenotype_types_are_declared_but_origins_constrained() -> None:
    cfg = DataConfig.model_validate({"observed": {"phenotypes": [
        {"id": "S1", "name": "s", "path": "p", "value_column": "t", "type": "survival"}]}})
    assert cfg.observed.phenotypes[0].type == PhenotypeType.SURVIVAL
    assert Origin("predicted") == Origin.PREDICTED
