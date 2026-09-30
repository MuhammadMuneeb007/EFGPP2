"""Installing the data-layer tools into the project, one conda environment per tool.

Everything lands in the project directory you work in, in plain sight:

    <project>/software/bin/            commands (links to each tool's environment / wrappers)
    <project>/software/envs/<tool>/    one conda environment per tool (never the efgpp env)
    <project>/software/opt/<tool>/     sources, unpacked downloads, container images
    <project>/software/downloads/      downloaded archives

(`--shared` installs into $EFGPP_TOOLS_HOME instead, only when you set it.)

For each tool EFGPP first creates its own conda environment from conda-forge + Bioconda
(mamba > micromamba > conda; micromamba is downloaded when none is installed). Only if that
fails does it fall back to the official source:
    official binary   PLINK 2, PLINK 1.9, FlashPCA2 (FlashPCA2 is not in conda at all)
    Python venv       MultiQC, OpenCRAVAT, GWASLab                         (uv, else venv+pip)
    source build      htslib (tabix, bgzip) + bcftools                     (needs gcc, make, zlib)
    source + venv     MetaXcan / PrediXcan                                 (GitHub archive + Python 3.11)
    container         Ensembl VEP                                          (Apptainer/Singularity/Docker)
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import zipfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import httpx

from efgpp.project import Project
from efgpp.resources.downloader import download
from efgpp.setup.platform import PlatformInfo, detect
from efgpp.setup.tools import shared_root

PLINK2_PAGE = "https://www.cog-genomics.org/plink/2.0/"
PLINK19_PAGE = "https://www.cog-genomics.org/plink/1.9/"
FLASHPCA_RELEASE = "https://github.com/gabraham/flashpca/releases/download/v2.0"
SAMTOOLS_API = "https://api.github.com/repos/samtools/{repo}/releases/latest"
METAXCAN_REPO = "hakyimlab/MetaXcan"
VEP_IMAGE = "docker://ensemblorg/ensembl-vep:latest"

PYTHON_TOOLS: dict[str, tuple[str, list[str], list[str]]] = {
    # tool: (environment name, pip requirements, executables exposed in bin/)
    "multiqc": ("reporting", ["multiqc"], ["multiqc"]),
    "oc": ("opencravat", ["open-cravat"], ["oc"]),
    # fallback when Bioconda's spliceai cannot be solved: pip, with TensorFlow inside this environment only
    "spliceai": ("spliceai", ["spliceai", "tensorflow>=2.10,<2.16"], ["spliceai"]),
}
METAXCAN_REQUIREMENTS = ["numpy<2", "scipy<1.14", "pandas>=2,<2.3", "sqlalchemy", "patsy", "statsmodels", "h5py",
                         "pyliftover", "cyvcf2", "bgen-reader"]

Progress = Callable[[str], None]


class InstallError(RuntimeError):
    pass


@dataclass
class InstallResult:
    tool: str
    method: str
    location: Path


# tool: (environment, conda packages, pip packages, commands exposed in software/bin)
CONDA_TOOLS: dict[str, tuple[str, list[str], list[str], list[str]]] = {
    "plink2": ("plink2", ["plink2"], [], ["plink2"]),
    "plink": ("plink", ["plink"], [], ["plink"]),
    "bcftools": ("bcftools", ["bcftools", "htslib"], [], ["bcftools", "tabix", "bgzip"]),
    "tabix": ("bcftools", ["bcftools", "htslib"], [], ["bcftools", "tabix", "bgzip"]),
    "bgzip": ("bcftools", ["bcftools", "htslib"], [], ["bcftools", "tabix", "bgzip"]),
    "multiqc": ("multiqc", ["multiqc"], [], ["multiqc"]),
    "oc": ("opencravat", ["open-cravat"], [], ["oc"]),
    "vep": ("vep", ["ensembl-vep", "perl", "htslib"], [], ["vep", "vep_install"]),
    # MetaXcan (Python 3): source pinned to a resolved commit, see metaxcan_source()
    "predixcan": ("predixcan", ["python=3.11", "numpy<2", "scipy", "pandas<2.3", "sqlalchemy", "patsy", "statsmodels",
                                "h5py", "cyvcf2", "pip"], ["bgen-reader", "pyliftover"], []),
    # TensorFlow < 2.16 keeps Keras 2, which loads SpliceAI's bundled .h5 models
    "spliceai": ("spliceai", ["python=3.10", "spliceai", "tensorflow>=2.10,<2.16", "pip"], [], ["spliceai"]),
    # GWASLab pins pysam/matplotlib/pandas versions: it must never share the efgpp environment.
    # GWASLab >= 4 lifts over with the separate `sumstats-liftover` package (not a declared dependency)
    "gwaslab": ("gwaslab", ["python=3.11", "gwaslab", "pyliftover", "pyarrow", "pip"], ["sumstats-liftover"],
                ["python:gwaslab-python"]),
}


def install_root(project: Project | None, shared: bool) -> Path:
    """<project>/software by default (./software outside a project); $EFGPP_TOOLS_HOME with --shared."""
    if shared:
        root = shared_root()
        if root is None:
            raise InstallError("--shared needs EFGPP_TOOLS_HOME (e.g. export EFGPP_TOOLS_HOME=$PWD/software); "
                               "without --shared everything installs into the project")
        return root
    if project is not None:
        return project.software_dir
    return Path.cwd() / "software"


def _exe(name: str, info: PlatformInfo) -> str:
    return f"{name}.exe" if info.is_windows else name


def _make_executable(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _expose(target: Path, bin_dir: Path, name: str | None = None) -> Path:
    """Make `target` callable from bin/ (symlink on POSIX, copy on Windows)."""
    bin_dir.mkdir(parents=True, exist_ok=True)
    link = bin_dir / (name or target.name)
    if link.exists() or link.is_symlink():
        link.unlink()
    if os.name == "nt":
        shutil.copy2(target, link)
    else:
        link.symlink_to(target)
    return link


def _page_links(url: str, pattern: str) -> list[str]:
    html = httpx.get(url, follow_redirects=True, timeout=60).text
    return re.findall(pattern, html)


def _extract_zip_member(archive: Path, names: tuple[str, ...], dest: Path) -> Path:
    with zipfile.ZipFile(archive) as z:
        for member in z.namelist():
            if member.rsplit("/", 1)[-1] in names:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(z.read(member))
                _make_executable(dest)
                return dest
    raise InstallError(f"{archive.name} does not contain {names}")


# --------------------------------------------------------------------- binaries
def plink2_asset(info: PlatformInfo, links: list[str]) -> str | None:
    """Choose the official PLINK 2 build for this OS / CPU from the download-page links."""
    if info.os == "windows":
        prefs = ["plink2_win_avx2_", "plink2_win64_"] if info.avx2 else ["plink2_win64_"]
    elif info.os == "darwin":
        prefs = ["plink2_mac_arm64_"] if info.arch == "aarch64" else (
            ["plink2_mac_avx2_", "plink2_mac_"] if info.avx2 else ["plink2_mac_"])
    elif info.arch == "aarch64":
        prefs = ["plink2_linux_arm64_", "plink2_linux_aarch64_"]
    else:
        vendor_amd = False
        try:
            with open("/proc/cpuinfo", encoding="utf-8") as fh:
                vendor_amd = "AuthenticAMD" in fh.read(4096)
        except OSError:
            pass
        prefs = (["plink2_linux_amd_avx2_"] if vendor_amd and info.avx2 else []) + \
                (["plink2_linux_avx2_"] if info.avx2 else []) + ["plink2_linux_x86_64_"]
    for pref in prefs:
        for link in links:
            if link.rsplit("/", 1)[-1].startswith(pref):
                return link
    return None


def install_plink2(root: Path, info: PlatformInfo, say: Progress) -> InstallResult:
    links = _page_links(PLINK2_PAGE, r'href="(https://s3\.amazonaws\.com/plink2-assets/[^"]+\.zip)"')
    url = plink2_asset(info, links)
    if url is None:
        raise InstallError(f"no official PLINK 2 build for {info.os}/{info.arch}")
    say(f"downloading {url}")
    archive = root / "downloads" / url.rsplit("/", 1)[-1]
    download(url, archive)
    exe = _extract_zip_member(archive, ("plink2", "plink2.exe"), root / "bin" / _exe("plink2", info))
    return InstallResult("plink2", f"official binary {url.rsplit('/', 1)[-1]}", exe)


def plink19_asset(info: PlatformInfo, links: list[str]) -> str | None:
    key = {"windows": "plink_win64_", "darwin": "plink_mac_"}.get(info.os, "plink_linux_x86_64_")
    if info.os == "linux" and info.arch != "x86_64":
        return None
    dated = [lk for lk in links if lk.rsplit("/", 1)[-1].startswith(key) and "/dev/" not in lk]
    return sorted(dated)[-1] if dated else None


def install_plink19(root: Path, info: PlatformInfo, say: Progress) -> InstallResult:
    links = _page_links(PLINK19_PAGE, r'href="(https://s3\.amazonaws\.com/plink1-assets/[^"]+\.zip)"')
    url = plink19_asset(info, links)
    if url is None:
        raise InstallError(f"no official PLINK 1.9 build for {info.os}/{info.arch}")
    say(f"downloading {url}")
    archive = root / "downloads" / url.rsplit("/", 1)[-1]
    download(url, archive)
    exe = _extract_zip_member(archive, ("plink", "plink.exe"), root / "bin" / _exe("plink", info))
    return InstallResult("plink", f"official binary {url.rsplit('/', 1)[-1]}", exe)


def install_flashpca(root: Path, info: PlatformInfo, say: Progress) -> InstallResult:
    if info.os == "linux" and info.arch == "x86_64":
        url, gz = f"{FLASHPCA_RELEASE}/flashpca_x86-64.gz", True
    elif info.os == "darwin" and info.arch == "x86_64":
        url, gz = f"{FLASHPCA_RELEASE}/flashpca-mac-intel", False
    else:
        raise InstallError(f"FlashPCA2 publishes no binary for {info.os}/{info.arch} (optional; PLINK 2 PCA is the default)")
    say(f"downloading {url}")
    raw = root / "downloads" / url.rsplit("/", 1)[-1]
    download(url, raw)
    exe = root / "bin" / "flashpca"
    exe.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(raw, "rb") as fin, open(exe, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copy2(raw, exe)
    _make_executable(exe)
    return InstallResult("flashpca2", "official static binary v2.0", exe)


# ---------------------------------------------------------------- source build
def _latest_samtools(repo: str) -> tuple[str, str]:
    r = httpx.get(SAMTOOLS_API.format(repo=repo), follow_redirects=True, timeout=60)
    r.raise_for_status()
    data = r.json()
    url = next(a["browser_download_url"] for a in data["assets"] if a["name"].endswith(".tar.bz2"))
    return data["tag_name"], url


def _run(cmd: list[str], cwd: Path, log: Path) -> None:
    with open(log, "a", encoding="utf-8") as fh:
        fh.write(f"$ {' '.join(cmd)}\n")
        fh.flush()
        proc = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        tail = "".join(log.read_text(encoding="utf-8", errors="replace").splitlines(True)[-15:])
        raise InstallError(f"`{' '.join(cmd)}` failed in {cwd} (log: {log})\n{tail}")


def install_htslib_bcftools(root: Path, info: PlatformInfo, say: Progress) -> InstallResult:
    """Build htslib (tabix, bgzip) and bcftools from the official release tarballs."""
    if info.is_windows:
        raise InstallError("building bcftools needs a Unix toolchain; use WSL2 or Docker on Windows")
    if shutil.which("make") is None or not (shutil.which("cc") or shutil.which("gcc")):
        raise InstallError("building bcftools from source needs a C compiler and make (e.g. `module load gcc`)")
    prefix = root / "opt" / "samtools"
    build = root / "opt" / "build"
    build.mkdir(parents=True, exist_ok=True)
    log = build / "build.log"
    log.write_text("", encoding="utf-8")
    jobs = str(max(1, min(8, info.cpu_count)))
    sources = {}
    for repo in ("htslib", "bcftools"):
        version, url = _latest_samtools(repo)
        say(f"downloading {repo} {version}")
        archive = root / "downloads" / url.rsplit("/", 1)[-1]
        download(url, archive)
        with tarfile.open(archive) as tar:
            tar.extractall(build, filter="data")
        sources[repo] = build / f"{repo}-{version}"
    # Minimal dependencies: only zlib is required for VCF/BCF work.
    say("building htslib (tabix, bgzip)")
    _run(["./configure", f"--prefix={prefix}", "--disable-bz2", "--disable-lzma", "--disable-libcurl"],
         sources["htslib"], log)
    _run(["make", f"-j{jobs}"], sources["htslib"], log)
    _run(["make", "install"], sources["htslib"], log)
    say("building bcftools")
    _run(["./configure", f"--prefix={prefix}", f"--with-htslib={sources['htslib']}"], sources["bcftools"], log)
    _run(["make", f"-j{jobs}"], sources["bcftools"], log)
    _run(["make", "install"], sources["bcftools"], log)
    for tool in ("bcftools", "tabix", "bgzip"):
        _expose(prefix / "bin" / tool, root / "bin")
    return InstallResult("bcftools", "built from official source (htslib + bcftools)", prefix / "bin")


# ------------------------------------------------------------ Python tools
def _uv() -> str | None:
    try:
        from uv import find_uv_bin

        return str(find_uv_bin())
    except (ImportError, FileNotFoundError):
        return shutil.which("uv")


def create_venv(prefix: Path, requirements: list[str], say: Progress, python: str | None = None) -> Path:
    """Isolated Python environment at `prefix` with `requirements`; returns its interpreter."""
    uv = _uv()
    py = prefix / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    log = prefix.parent / f"{prefix.name}.install.log"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("", encoding="utf-8")
    if uv:
        say(f"creating {prefix.name} environment with uv")
        _run([uv, "venv", "--allow-existing", *(["--python", python] if python else []), str(prefix)], prefix.parent, log)
        _run([uv, "pip", "install", "--python", str(py), *requirements], prefix.parent, log)
    else:
        say(f"creating {prefix.name} environment with venv + pip")
        _run([sys.executable, "-m", "venv", str(prefix)], prefix.parent, log)
        _run([str(py), "-m", "pip", "install", "--upgrade", "pip"], prefix.parent, log)
        _run([str(py), "-m", "pip", "install", *requirements], prefix.parent, log)
    return py


def _venv_bin(prefix: Path) -> Path:
    return prefix / ("Scripts" if os.name == "nt" else "bin")


def install_python_tool(tool: str, root: Path, info: PlatformInfo, say: Progress) -> InstallResult:
    env, requirements, executables = PYTHON_TOOLS[tool]
    prefix = root / "envs" / env
    if (prefix / "conda-meta").exists():
        prefix = root / "envs" / f"{env}-venv"  # never write into a conda environment
    create_venv(prefix, requirements, say)
    for exe in executables:
        target = _venv_bin(prefix) / _exe(exe, info)
        if not target.exists():
            raise InstallError(f"{exe} missing after installing {requirements}")
        _expose(target, root / "bin")
    return InstallResult(tool, f"Python environment {prefix.name} ({', '.join(requirements)})", prefix)


def install_gwaslab_venv(root: Path, info: PlatformInfo, say: Progress) -> InstallResult:
    """Fallback: GWASLab in its own Python 3.11 virtual environment."""
    prefix = root / "envs" / "gwaslab-venv"
    py = create_venv(prefix, ["gwaslab", "pyliftover", "pyarrow"], say, python="3.11")
    _expose(py, root / "bin", "gwaslab-python")
    return InstallResult("gwaslab", "Python environment gwaslab-venv (gwaslab, pyliftover)", prefix)


def install_metaxcan(root: Path, info: PlatformInfo, say: Progress) -> InstallResult:
    """Fallback: MetaXcan source (GitHub archive, no git needed) + its own Python 3.11 venv."""
    prefix = root / "envs" / "predixcan-venv"
    create_venv(prefix, METAXCAN_REQUIREMENTS, say, python="3.11")
    py = prefix / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    wrapper = metaxcan_source(root, py, info, say)
    return InstallResult("predixcan", "MetaXcan source + Python 3.11 environment", wrapper)


def metaxcan_source(root: Path, py: Path, info: PlatformInfo, say: Progress) -> Path:
    """Download MetaXcan at a resolved commit (pinned in opt/MetaXcan/EFGPP_SOURCE.json and the lock
    file) and write the `predixcan` wrapper that runs Predict.py with `py`."""
    from efgpp.setup.toolkit_installer import resolve_commit, write_source_pin

    commit = resolve_commit(METAXCAN_REPO, "master")
    if commit is None:
        raise InstallError("could not resolve the MetaXcan commit (GitHub API unreachable); nothing installed")
    say(f"downloading MetaXcan @ {commit[:10]}")
    archive = root / "downloads" / f"MetaXcan-{commit[:12]}.zip"
    archive_sha = download(f"https://github.com/{METAXCAN_REPO}/archive/{commit}.zip", archive)
    src = root / "opt" / "MetaXcan"
    shutil.rmtree(src, ignore_errors=True)
    staging = root / "opt" / ".MetaXcan.tmp"
    shutil.rmtree(staging, ignore_errors=True)
    with zipfile.ZipFile(archive) as z:
        z.extractall(staging)
    next(staging.iterdir()).rename(src)
    shutil.rmtree(staging, ignore_errors=True)
    write_source_pin(src, METAXCAN_REPO, "master", commit, archive_sha)
    script = src / "software" / "Predict.py"
    if info.is_windows:
        wrapper = root / "bin" / "predixcan.bat"
        wrapper.parent.mkdir(parents=True, exist_ok=True)
        wrapper.write_text(f'@"{py}" "{script}" %*\r\n', encoding="utf-8")
    else:
        wrapper = root / "bin" / "predixcan"
        wrapper.parent.mkdir(parents=True, exist_ok=True)
        wrapper.write_text(f'#!/bin/sh\nexec "{py}" "{script}" "$@"\n', encoding="utf-8")
        _make_executable(wrapper)
    return wrapper


# ------------------------------------------------------------------- containers
def install_vep_container(root: Path, info: PlatformInfo, say: Progress) -> InstallResult:
    """Official Ensembl VEP image behind `vep` / `vep_install` wrapper scripts."""
    if info.is_windows:
        raise InstallError("run VEP through WSL2 on Windows")
    runtime = next((r for r in ("apptainer", "singularity") if shutil.which(r)), None)
    bin_dir = root / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    project_root = root.parent if root.name == "software" else Path.cwd()
    if runtime:
        sif = root / "opt" / "containers" / "ensembl-vep.sif"
        sif.parent.mkdir(parents=True, exist_ok=True)
        if not sif.exists():
            say(f"pulling {VEP_IMAGE} with {runtime} (large image; one time)")
            log = sif.with_suffix(".log")
            log.write_text("", encoding="utf-8")
            _run([runtime, "pull", str(sif), VEP_IMAGE], sif.parent, log)
        run = f'exec {runtime} exec --bind "{project_root}" --bind "$PWD" "{sif}"'
        method = f"{runtime} image {sif.name}"
    elif shutil.which("docker"):
        say(f"pulling {VEP_IMAGE} with docker")
        subprocess.run(["docker", "pull", VEP_IMAGE.removeprefix("docker://")], check=True)
        run = (f'exec docker run --rm -u "$(id -u):$(id -g)" -v "{project_root}:{project_root}" '
               f'-v "$PWD:$PWD" -w "$PWD" {VEP_IMAGE.removeprefix("docker://")}')
        method = "docker image ensemblorg/ensembl-vep"
    else:
        raise InstallError("VEP could not be installed in its own conda environment and no container runtime "
                           "(Apptainer/Singularity/Docker) is available. Do NOT install ensembl-vep into the efgpp "
                           "environment; see software/logs/vep.install.log")
    for name, command in (("vep", "vep"), ("vep_install", "INSTALL.pl")):
        wrapper = bin_dir / name
        wrapper.write_text(f'#!/bin/sh\n{run} {command} "$@"\n', encoding="utf-8")
        _make_executable(wrapper)
    return InstallResult("vep", method, bin_dir / "vep")


# ----------------------------------------------------------------------- registry
INSTALLERS: dict[str, Callable[[Path, PlatformInfo, Progress], InstallResult]] = {
    "plink2": install_plink2,
    "plink": install_plink19,
    "flashpca2": install_flashpca,
    "bcftools": install_htslib_bcftools,
    "tabix": install_htslib_bcftools,
    "bgzip": install_htslib_bcftools,
    "multiqc": lambda r, i, s: install_python_tool("multiqc", r, i, s),
    "oc": lambda r, i, s: install_python_tool("oc", r, i, s),
    "spliceai": lambda r, i, s: install_python_tool("spliceai", r, i, s),
    "predixcan": install_metaxcan,
    "gwaslab": lambda r, i, s: install_gwaslab_venv(r, i, s),
    "vep": install_vep_container,
}

# What `efgpp setup data --components ...` needs from each component.
COMPONENT_TOOLS = {
    "genetics": ["plink2", "plink", "bcftools", "tabix", "bgzip", "flashpca2"],
    "annotation": ["vep", "oc"],
    "metaxcan": ["predixcan"],
    "reporting": ["multiqc"],
    "gwas": ["gwaslab"],
    "spliceai": ["spliceai"],
}


def install_tools(project: Project | None, tools: list[str], *, shared: bool = False, force: bool = False,
                  progress: Progress | None = None) -> list[tuple[str, str, str]]:
    """Install each tool that is missing. Returns (tool, status, detail) rows."""
    from efgpp.setup.tools import available

    say = progress or (lambda _m: None)
    info = detect()
    root = install_root(project, shared)
    rows: list[tuple[str, str, str]] = []
    done: set[Callable[..., InstallResult]] = set()
    done_envs: set[str] = set()
    for tool in dict.fromkeys(tools):
        if tool not in INSTALLERS:
            rows.append((tool, "failed", f"no installer; choose from {', '.join(sorted(INSTALLERS))}"))
            continue
        if not force and available(project, tool):
            rows.append((tool, "ok", "already available"))
            continue
        if tool in CONDA_TOOLS and not info.is_windows:
            env = CONDA_TOOLS[tool][0]
            if env in done_envs:
                rows.append((tool, "installed", f"with the {env} environment"))
                continue
            try:
                say(f"installing {tool} into its own conda environment software/envs/{env}")
                res = install_conda_tool(tool, root, info, say)
                done_envs.add(env)
                rows.append((tool, "installed", f"{res.method} -> {res.location}"))
                continue
            except Exception as exc:  # noqa: BLE001 - fall back to the official download
                say(f"conda install of {tool} failed ({str(exc).splitlines()[0]}); using the official source")
        fn = INSTALLERS[tool]
        if fn in done:  # e.g. tabix/bgzip come with the bcftools build
            rows.append((tool, "installed", "with bcftools"))
            continue
        try:
            say(f"installing {tool}")
            res = fn(root, info, say)
            done.add(fn)
            rows.append((tool, "installed", f"{res.method} -> {res.location}"))
        except Exception as exc:  # noqa: BLE001 - report every tool, keep going
            rows.append((tool, "failed", str(exc)))
    return rows


def install_conda_tool(tool: str, root: Path, info: PlatformInfo, say: Progress) -> InstallResult:
    """Create the tool's own conda environment (software/envs/<env>) and expose its commands."""
    from efgpp.setup.toolkit_installer import conda_tool, ensure_conda_env, env_bin, env_prefix

    env, packages, pip, commands = CONDA_TOOLS[tool]
    log = root / "logs" / f"{env}.install.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("", encoding="utf-8")
    manager = conda_tool(None, root, say)
    ensure_conda_env(manager, env_prefix(root, env), packages, log, say)
    py = env_bin(root, env) / "python"
    for pkg in pip:
        _run([str(py), "-m", "pip", "install", "--no-input", pkg], root, log)
    for command in commands:
        source_name, _, exposed_name = command.partition(":")
        target = env_bin(root, env) / source_name
        if not target.exists():
            raise InstallError(f"{source_name} missing from the {env} environment (log: {log})")
        _expose(target, root / "bin", exposed_name or None)
    location = env_prefix(root, env)
    if tool == "predixcan":
        location = metaxcan_source(root, py, info, say)
    return InstallResult(tool, f"conda environment software/envs/{env} via {manager[0]}", location)


def path_exports(project: Project | None) -> list[Path]:
    """Directories to put on PATH to call EFGPP-installed tools directly from a shell."""
    dirs = [project.bin_dir, *(d / "bin" for d in project.legacy_software_dirs)] if project else [Path.cwd() / "software" / "bin"]
    shared = shared_root()
    if shared is not None:
        dirs.append(shared / "bin")
    return [d for d in dirs if d.exists()]
