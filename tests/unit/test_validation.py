"""Validation collects every problem (section 55 of the specification)."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from efgpp.config.data import PhenotypeSource
from efgpp.constants import FAIL, PASS, WARN, PhenotypeType
from efgpp.data.adapters import AdapterContext, CovariateAdapter, PhenotypeAdapter
from efgpp.data.adapters.phenotype import TYPE_HANDLERS
from efgpp.data.registry import Registry
from efgpp.project import Project


def _pheno(project: Project, path: Path, **kw) -> PhenotypeSource:  # type: ignore[no-untyped-def]
    return PhenotypeSource(id="PH001", name="trait_a", path=str(path), participant_id_column="IID",
                           value_column="trait_a", **kw)


def test_all_errors_are_collected(project: Project, tmp_path: Path) -> None:
    ids = [f"P{i}" for i in range(100)] + [f"P{i}" for i in range(22)]  # 22 duplicated IDs
    values = ["0", "1"] * 50 + ["7", "8", "x", "9"] + ["NA"] * 18  # 4 invalid, 18 missing (+1 below)
    pl.DataFrame({"IID": ids, "trait_a": values}).write_csv(tmp_path / "p.csv")
    with Registry.open(project) as reg:
        rep = PhenotypeAdapter(AdapterContext(project, reg), _pheno(project, tmp_path / "p.csv", type=PhenotypeType.BINARY)).validate()
    by = {c.message: c for c in rep.checks}
    assert by["participant ID column exists (IID)"].status == PASS
    assert by["phenotype type recognized (binary)"].status == PASS
    assert by["22 duplicated participant IDs"].status == FAIL
    assert by["4 invalid binary values"].status == FAIL
    assert by["18 missing values"].status == WARN
    assert not rep.passed and rep.n_fail == 2


def test_missing_columns_and_unsupported_type(project: Project, tmp_path: Path) -> None:
    pl.DataFrame({"ID": ["a"], "v": ["1"]}).write_csv(tmp_path / "p.csv")
    with Registry.open(project) as reg:
        ctx = AdapterContext(project, reg)
        rep = PhenotypeAdapter(ctx, _pheno(project, tmp_path / "p.csv", type=PhenotypeType.BINARY)).validate()
        assert sum(c.status == FAIL for c in rep.checks) == 2  # both columns missing, both reported
        rep2 = PhenotypeAdapter(ctx, _pheno(project, tmp_path / "p.csv", type=PhenotypeType.SURVIVAL)).validate()
        assert any("not supported in v0.1" in c.message for c in rep2.checks)


def test_type_handlers() -> None:
    src = PhenotypeSource(id="X", name="n", path="p", value_column="v", type=PhenotypeType.BINARY)
    vals = pl.Series(["1", "2", "2", "-9"])
    miss = vals == "-9"
    coded = TYPE_HANDLERS[PhenotypeType.BINARY].code(vals, miss, src)
    assert coded.numeric.to_list() == [0.0, 1.0, 1.0, None]  # PLINK 1/2 coding detected
    assert "PLINK coding" in coded.notes[0]
    ordinal = PhenotypeSource(id="X", name="n", path="p", value_column="v", type=PhenotypeType.ORDINAL,
                              levels=["low", "mid", "high"])
    c = TYPE_HANDLERS[PhenotypeType.ORDINAL].code(pl.Series(["high", "low", "bad"]), pl.Series([False] * 3), ordinal)
    assert c.numeric.to_list() == [2.0, 0.0, None] and c.invalid.to_list() == [False, False, True]
    custom = PhenotypeSource(id="X", name="n", path="p", value_column="v", type=PhenotypeType.BINARY,
                             case_values=["yes"], control_values=["no"])
    c = TYPE_HANDLERS[PhenotypeType.BINARY].code(pl.Series(["yes", "no", "maybe"]), pl.Series([False] * 3), custom)
    assert c.numeric.to_list() == [1.0, 0.0, None]


def test_covariates_never_imputed(project: Project, tmp_path: Path) -> None:
    from efgpp.config.data import CovariateSource

    pl.DataFrame({"IID": ["a", "b", "c"], "age": ["40", "NA", "60"], "sex": ["F", "M", "F"]}).write_csv(tmp_path / "c.csv")
    src = CovariateSource(id="COV001", path=str(tmp_path / "c.csv"), participant_id_column="IID", variables=["age", "sex"])
    project.data.observed.covariates.append(src)
    project.save_data_config()
    project = Project.load(project.root)
    with Registry.open(project) as reg:
        ad = CovariateAdapter(AdapterContext(project, reg), project.data.observed.covariates[0])
        ad.register()
        rep = ad.validate()
        assert any("treated as categorical" in c.message for c in rep.checks)
        art = ad.standardize()
    out = pl.read_parquet(art.path)
    assert out.get_column("age").to_list() == [40.0, None, 60.0]  # missing stays missing, no scaling
