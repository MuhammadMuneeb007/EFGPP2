"""Pixi (preferred) bootstrap and per-purpose environments from conda-forge + bioconda."""

from __future__ import annotations

import platform as _platform
import shutil
import subprocess
import tarfile
import zipfile
from pathlib import Path

from efgpp.project import Project
from efgpp.resources.downloader import download
from efgpp.setup.platform import PlatformInfo
from efgpp.workflow.templates import pixi_manifest

RELEASE = "https://github.com/prefix-dev/pixi/releases/latest/download"
ASSETS = {
    ("linux", "x86_64"): "pixi-x86_64-unknown-linux-musl.tar.gz",
    ("linux", "aarch64"): "pixi-aarch64-unknown-linux-musl.tar.gz",
    ("darwin", "x86_64"): "pixi-x86_64-apple-darwin.tar.gz",
    ("darwin", "aarch64"): "pixi-aarch64-apple-darwin.tar.gz",
    ("windows", "x86_64"): "pixi-x86_64-pc-windows-msvc.zip",
}


def conda_platform(info: PlatformInfo) -> str:
    os_part = {"linux": "linux", "darwin": "osx", "windows": "win"}[info.os]
    arch = {"x86_64": "64", "aarch64": "aarch64" if info.os == "linux" else "arm64"}.get(info.arch, "64")
    return f"{os_part}-{arch}"


def find_pixi(project: Project) -> Path | None:
    exe = "pixi.exe" if _platform.system() == "Windows" else "pixi"
    local = project.bin_dir / exe
    if local.exists():
        return local
    found = shutil.which("pixi")
    return Path(found) if found else None


def bootstrap(project: Project, info: PlatformInfo) -> Path:
    """Download the pixi binary into .efgpp/bin (no shell installer scripts are piped)."""
    existing = find_pixi(project)
    if existing:
        return existing
    asset = ASSETS.get((info.os, info.arch))
    if asset is None:
        raise RuntimeError(f"no pixi build for {info.os}/{info.arch}")
    archive = project.path(".efgpp", "downloads", asset)
    download(f"{RELEASE}/{asset}", archive)
    project.bin_dir.mkdir(parents=True, exist_ok=True)
    if asset.endswith(".zip"):
        with zipfile.ZipFile(archive) as z:
            z.extractall(project.bin_dir)
    else:
        with tarfile.open(archive) as t:
            t.extractall(project.bin_dir, filter="data")
    exe = project.bin_dir / ("pixi.exe" if info.os == "windows" else "pixi")
    exe.chmod(0o755)
    return exe


def create_env(project: Project, pixi: Path, name: str, info: PlatformInfo) -> subprocess.CompletedProcess[str]:
    env_dir = project.envs_dir / name
    env_dir.mkdir(parents=True, exist_ok=True)
    manifest = env_dir / "pixi.toml"
    manifest.write_text(pixi_manifest(name, conda_platform(info)), encoding="utf-8")
    return subprocess.run([str(pixi), "install", "--manifest-path", str(manifest)],
                          capture_output=True, text=True)
