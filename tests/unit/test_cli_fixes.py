"""Regression tests for problems met on the first real cohort (legacy migraine data)."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import yaml

from efgpp.constants import GenomeBuild
from efgpp.data.genotype.build import infer_build
from efgpp.project import Project
from efgpp.workflow.executor import choose_engine


def test_positions_past_several_grch38_ends_mean_grch37() -> None:
    # Real case: chr1:249222527, chr2:243041411, chr4:190904950 exceed GRCh38 but not GRCh37.
    v = pl.DataFrame({"chromosome": ["1", "2", "4"], "position": [249_222_527, 243_041_411, 190_904_950],
                      "variant_id": ["a", "b", "c"], "reference": ["A"] * 3, "alternate": ["G"] * 3})
    result = infer_build(declared="auto", variants=v)
    assert result.build == GenomeBuild.GRCH37 and result.confident


def test_local_runs_use_the_builtin_executor(project: Project) -> None:
    assert choose_engine(project, None) == "builtin"
    assert choose_engine(project, "builtin") == "builtin"


def test_covariates_default_to_all_columns_and_reject_unknown(project: Project, tmp_path: Path, run_cli,
                                                               monkeypatch) -> None:  # type: ignore[no-untyped-def]
    cov = tmp_path / "migraine.cov"
    cov.write_text("FID IID Sex Age PC1\nA A 1 50 0.1\nB B 2 60 0.2\n", encoding="utf-8")
    monkeypatch.chdir(project.root)
    bad = run_cli("data", "add", "covariates", "--path", str(cov), "--id-column", "IID",
                  "--columns", "COV1", "COV2", expect=2)
    assert "columns not in the file: COV1, COV2" in bad.output
    run_cli("data", "add", "covariates", "--path", str(cov), "--id-column", "IID", "--categorical", "Sex")
    data = yaml.safe_load((project.root / "data.yaml").read_text(encoding="utf-8"))
    assert data["observed"]["covariates"][0]["variables"] == ["Sex", "Age", "PC1"]
    run_cli("data", "remove", "COV001")
    data = yaml.safe_load((project.root / "data.yaml").read_text(encoding="utf-8"))
    assert data["observed"].get("covariates", []) == []
