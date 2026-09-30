"""Shared fixtures: tiny synthetic cohorts and projects. No multi-GB resources are used."""

from __future__ import annotations

import os
import shutil
from collections.abc import Callable
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from efgpp.data.simulate import simulate_genotypes, simulate_phenotype
from efgpp.project import Project, init_project

N = 120  # participants

# Tests never call Ensembl / UCSC for build checks or chain files.
os.environ["EFGPP_OFFLINE"] = "1"
M = 400  # variants


@pytest.fixture
def project(tmp_path: Path) -> Project:
    p = init_project(tmp_path / "proj", name="test_cohort")
    # Keep test paths short and portable.
    p.config.execution.local_cores = 2
    p.save_project_config()
    return Project.load(p.root)


@pytest.fixture(scope="session")
def cohort_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """An external 'data' folder with genotype + phenotype/covariate/omics tables."""
    root = tmp_path_factory.mktemp("external_cohort")
    cohort = simulate_genotypes(N, M, 2, seed=7, prefix=root / "genotype" / "cohort")
    rng = np.random.default_rng(3)
    ids = cohort.iids
    binary, _ = simulate_phenotype(cohort, _pt("binary"), 0.5, 10, 0.3, seed=1)
    cont, _ = simulate_phenotype(cohort, _pt("continuous"), 0.5, 10, 0.3, seed=2)
    pheno = pl.DataFrame({
        "IID": ids,
        "trait_a": [str(int(v)) for v in binary],
        "trait_b": [f"{v:.5f}" for v in cont],
        "trait_c": rng.choice(["x", "y", "z"], len(ids)).tolist(),
        "trait_d": rng.choice(["low", "mid", "high"], len(ids)).tolist(),
        "trait_e": [f"{v:.3f}" for v in rng.normal(10, 2, len(ids))],
    })
    pheno.write_csv(root / "phenotypes.csv")
    pl.DataFrame({"IID": ids, "age": rng.normal(50, 8, len(ids)).round(1),
                  "sex": rng.choice(["F", "M"], len(ids)).tolist(),
                  "bmi": rng.normal(26, 4, len(ids)).round(1)}).write_csv(root / "covariates.csv")
    # Measured expression for a subset of participants (samples x features).
    sub = ids[: N // 2]
    expr = pl.DataFrame({"participant_id": sub}).hstack(
        pl.DataFrame(rng.normal(5, 1, (len(sub), 8)), schema=[f"ENSG{i:011d}" for i in range(8)], orient="row"))
    expr.write_parquet(root / "expression.parquet")
    # Longitudinal expression: three visits with attrition.
    rows = []
    for i, pid in enumerate(ids[:40]):
        for day, visit in ((0, "baseline"), (180, "month_6"), (365, "month_12")):
            if visit == "month_12" and i % 4 == 0:
                continue
            if visit == "month_6" and i % 10 == 0:
                continue
            rows.append({"participant_id": pid, "visit_id": visit, "collection_date": f"2024-{1 + day // 31:02d}-01",
                         "biospecimen": f"B{pid}_{visit}", **{f"G{j}": float(rng.normal()) for j in range(5)}})
    pl.DataFrame(rows).write_csv(root / "expression_longitudinal.csv")
    return root


def _pt(name: str):  # type: ignore[no-untyped-def]
    from efgpp.constants import PhenotypeType

    return PhenotypeType(name)


def _find_plink2(cache: Path) -> Path | None:
    for candidate in (os.environ.get("EFGPP_TEST_PLINK2"), shutil.which("plink2")):
        if candidate and Path(candidate).exists():
            return Path(candidate)
    exe = cache / ("plink2.exe" if os.name == "nt" else "plink2")
    if exe.exists():
        return exe
    if os.environ.get("EFGPP_TEST_OFFLINE"):
        return None
    try:  # the same installer `efgpp setup data` uses
        from efgpp.setup.manager import install_plink2_binary
        from efgpp.setup.platform import detect

        tmp = init_project(cache / "bootstrap", name="bootstrap")
        install_plink2_binary(tmp, detect())
        found = next(tmp.bin_dir.glob("plink2*"))
        shutil.copy2(found, exe)
        exe.chmod(0o755)
        return exe
    except Exception:  # noqa: BLE001 - offline machines simply skip PLINK-dependent tests
        return None


@pytest.fixture(scope="session")
def plink2_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    cache = Path(os.environ.get("EFGPP_TEST_CACHE", tmp_path_factory.getbasetemp() / "plink_cache"))
    cache.mkdir(parents=True, exist_ok=True)
    exe = _find_plink2(cache)
    if exe is None:
        pytest.skip("plink2 not available (set EFGPP_TEST_PLINK2 or allow the download)")
    return exe


@pytest.fixture
def with_plink2(project: Project, plink2_path: Path) -> Project:
    project.config.execution.tools = {"plink2": str(plink2_path)}
    project.save_project_config()
    return Project.load(project.root)


@pytest.fixture
def run_cli() -> Callable[..., object]:
    from typer.testing import CliRunner

    from efgpp.cli.main import app

    runner = CliRunner()

    def invoke(*args: str, expect: int | None = 0):  # type: ignore[no-untyped-def]
        result = runner.invoke(app, list(args), catch_exceptions=False)
        if expect is not None:
            assert result.exit_code == expect, result.output
        return result

    return invoke
