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


def test_install_roots_and_path_line(project: Project, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(installers, "shared_root", lambda: tmp_path / "shared")
    assert installers.install_root(project, shared=False) == project.path(".efgpp")
    assert installers.install_root(project, shared=True) == tmp_path / "shared"
    assert installers.install_root(None, shared=False) == tmp_path / "shared"
    project.bin_dir.mkdir(parents=True, exist_ok=True)
    assert installers.path_exports(project) == [project.bin_dir]


def test_every_component_tool_has_an_installer() -> None:
    for tools in installers.COMPONENT_TOOLS.values():
        assert set(tools) <= set(installers.INSTALLERS)


def test_already_available_tools_are_skipped(project: Project, monkeypatch: pytest.MonkeyPatch) -> None:
    import efgpp.setup.tools as tools_mod

    monkeypatch.setattr(tools_mod, "available", lambda _p, _t: True)
    rows = installers.install_tools(project, ["plink2", "snakemake"])
    assert rows == [("plink2", "ok", "already available"), ("snakemake", "ok", "already available")]
