"""Selection logic of the direct installers (no downloads)."""

from __future__ import annotations

from pathlib import Path

import pytest

from efgpp.project import Project
from efgpp.setup import installers
from efgpp.setup.platform import PlatformInfo

PLINK2_LINKS = [f"https://s3.amazonaws.com/plink2-assets/alpha7/plink2_{p}_20260929.zip" for p in
                ("linux_amd_avx2", "linux_avx2", "linux_x86_64", "linux_arm64", "mac", "mac_arm64", "mac_avx2",
                 "win64", "win_avx2")]
PLINK19_LINKS = ["https://s3.amazonaws.com/plink1-assets/dev/plink_linux_x86_64.zip",
                 "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip",
                 "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20260927.zip",
                 "https://s3.amazonaws.com/plink1-assets/plink_mac_20260927.zip"]


def _info(os_: str = "linux", arch: str = "x86_64", avx2: bool = True) -> PlatformInfo:
    return PlatformInfo(os=os_, arch=arch, python="3.12", cpu_count=8, memory_gb=64, avx2=avx2, wsl=False,
                        hpc_scheduler="slurm")


def _name(url: str | None) -> str:
    assert url is not None
    return url.rsplit("/", 1)[-1]


def test_plink2_build_selection() -> None:
    assert _name(installers.plink2_asset(_info(avx2=False), PLINK2_LINKS)) == "plink2_linux_x86_64_20260929.zip"
    assert _name(installers.plink2_asset(_info("linux", "aarch64"), PLINK2_LINKS)) == "plink2_linux_arm64_20260929.zip"
    assert _name(installers.plink2_asset(_info("darwin", "aarch64"), PLINK2_LINKS)) == "plink2_mac_arm64_20260929.zip"
    assert _name(installers.plink2_asset(_info("windows"), PLINK2_LINKS)) == "plink2_win_avx2_20260929.zip"
    linux_avx2 = _name(installers.plink2_asset(_info(), PLINK2_LINKS))
    assert linux_avx2 in ("plink2_linux_avx2_20260929.zip", "plink2_linux_amd_avx2_20260929.zip")


def test_plink19_picks_newest_stable_build() -> None:
    assert _name(installers.plink19_asset(_info(), PLINK19_LINKS)) == "plink_linux_x86_64_20260927.zip"
    assert installers.plink19_asset(_info("linux", "aarch64"), PLINK19_LINKS) is None


def test_unsupported_platforms_fail_clearly(tmp_path: Path) -> None:
    with pytest.raises(installers.InstallError, match="no binary"):
        installers.install_flashpca(tmp_path, _info("linux", "aarch64"), lambda _m: None)
    with pytest.raises(installers.InstallError, match="WSL2"):
        installers.install_vep_container(tmp_path, _info("windows"), lambda _m: None)


def test_everything_installs_into_the_working_directory(project: Project, monkeypatch: pytest.MonkeyPatch,
                                                         tmp_path: Path) -> None:
    monkeypatch.delenv("EFGPP_TOOLS_HOME", raising=False)
    assert installers.install_root(project, shared=False) == project.root / "software"
    monkeypatch.chdir(tmp_path)
    assert installers.install_root(None, shared=False) == tmp_path / "software"  # never the home directory
    with pytest.raises(installers.InstallError, match="EFGPP_TOOLS_HOME"):
        installers.install_root(project, shared=True)
    monkeypatch.setenv("EFGPP_TOOLS_HOME", str(tmp_path / "tools"))
    assert installers.install_root(project, shared=True) == tmp_path / "tools"
    project.bin_dir.mkdir(parents=True, exist_ok=True)
    assert installers.path_exports(project) == [project.bin_dir]
    assert project.bin_dir == project.root / "software" / "bin" and project.logs_dir == project.root / "logs"


def test_each_conda_tool_gets_its_own_environment() -> None:
    envs = {tool: spec[0] for tool, spec in installers.CONDA_TOOLS.items()}
    assert envs["vep"] == "vep" and "ensembl-vep" in installers.CONDA_TOOLS["vep"][1]
    assert envs["bcftools"] == envs["tabix"] == envs["bgzip"] == "bcftools"
    assert len({envs[t] for t in ("plink2", "plink", "bcftools", "snakemake", "multiqc", "oc", "vep", "predixcan")}) == 8
    assert "efgpp" not in envs.values()


def test_conda_first_then_official_source(project: Project, monkeypatch: pytest.MonkeyPatch) -> None:
    import efgpp.setup.tools as tools_mod

    monkeypatch.setattr(tools_mod, "available", lambda _p, _t: False)
    monkeypatch.setattr(installers, "detect", lambda: _info())
    calls: list[str] = []

    def conda_ok(tool, root, info, say):  # type: ignore[no-untyped-def]
        calls.append(f"conda:{tool}")
        return installers.InstallResult(tool, "conda", root / "envs" / tool)

    monkeypatch.setattr(installers, "install_conda_tool", conda_ok)
    rows = installers.install_tools(project, ["bcftools", "tabix", "vep"])
    assert calls == ["conda:bcftools", "conda:vep"]  # tabix comes with the bcftools environment
    assert rows[1] == ("tabix", "installed", "with the bcftools environment")

    def conda_fails(tool, root, info, say):  # type: ignore[no-untyped-def]
        raise installers.InstallError("solver failed")

    def official(root, info, say):  # type: ignore[no-untyped-def]
        return installers.InstallResult("plink2", "official binary", root / "bin" / "plink2")

    monkeypatch.setattr(installers, "install_conda_tool", conda_fails)
    monkeypatch.setitem(installers.INSTALLERS, "plink2", official)
    assert installers.install_tools(project, ["plink2"])[0][2].startswith("official binary")


def test_every_component_tool_has_an_installer() -> None:
    for tools in installers.COMPONENT_TOOLS.values():
        assert set(tools) <= set(installers.INSTALLERS)


def test_already_available_tools_are_skipped(project: Project, monkeypatch: pytest.MonkeyPatch) -> None:
    import efgpp.setup.tools as tools_mod

    monkeypatch.setattr(tools_mod, "available", lambda _p, _t: True)
    rows = installers.install_tools(project, ["plink2", "snakemake"])
    assert rows == [("plink2", "ok", "already available"), ("snakemake", "ok", "already available")]
