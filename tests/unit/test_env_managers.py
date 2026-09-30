from __future__ import annotations

from pathlib import Path

import pytest

from efgpp.project import Project
from efgpp.setup import micromamba


@pytest.fixture
def fake_path(monkeypatch: pytest.MonkeyPatch):  # type: ignore[no-untyped-def]
    available: set[str] = set()
    monkeypatch.setattr(micromamba.shutil, "which", lambda name: f"/usr/bin/{name}" if name in available else None)
    return available


def test_fastest_available_conda_tool_is_chosen(project: Project, fake_path: set[str]) -> None:
    assert micromamba.find_conda_tool(project) is None
    fake_path.update({"conda"})
    assert micromamba.find_conda_tool(project)[0] == "conda"  # type: ignore[index]
    fake_path.update({"mamba"})
    assert micromamba.find_conda_tool(project)[0] == "mamba"  # type: ignore[index]
    fake_path.update({"micromamba"})
    assert micromamba.find_conda_tool(project)[0] == "micromamba"  # type: ignore[index]
    # An explicit preference wins when that tool exists.
    assert micromamba.find_conda_tool(project, "conda")[0] == "conda"  # type: ignore[index]


def test_create_commands(tmp_path: Path) -> None:
    spec, prefix = tmp_path / "genetics.yaml", tmp_path / "envs" / "genetics"
    assert micromamba.create_command("micromamba", Path("mm"), prefix, spec)[1:3] == ["create", "-y"]
    for kind in ("mamba", "conda"):
        cmd = micromamba.create_command(kind, Path(kind), prefix, spec)
        assert cmd[1:4] == ["env", "create", "-y"] and cmd[-2:] == ["-f", str(spec)] and str(prefix) in cmd


def test_default_fallback_is_mamba() -> None:
    from efgpp.config import ProjectConfig

    cfg = ProjectConfig()
    assert cfg.execution.environment_manager == "pixi"
    assert cfg.execution.fallback_environment_manager == "mamba"
    assert ProjectConfig.model_validate({"execution": {"environment_manager": "mamba"}}).execution.environment_manager == "mamba"
