"""`efgpp setup data` / `efgpp setup annotation`: detect, install, test and lock.

Scientific tools go into small isolated environments (Pixi preferred, Micromamba as
fallback) under .efgpp/envs/<purpose>. On native Windows, where Bioconda has no builds,
PLINK 2 is installed from its official Windows binary and Linux-only components are
reported as requiring WSL2 or Docker.
"""

from __future__ import annotations

import importlib.util
import re
import zipfile
from collections.abc import Callable
from dataclasses import dataclass, field

import httpx

from efgpp.data.registry import Registry, utcnow
from efgpp.project import Project
from efgpp.resources.downloader import download
from efgpp.setup import micromamba, pixi
from efgpp.setup.platform import PlatformInfo, detect
from efgpp.setup.tools import ToolNotFoundError, detect_version, resolve

PLINK2_PAGE = "https://www.cog-genomics.org/plink/2.0/"
CORE_MODULES = ("pydantic", "typer", "rich", "yaml", "polars", "pyarrow", "duckdb", "pandera", "numpy",
                "scipy", "plotly", "jinja2", "httpx", "filelock")


@dataclass
class SetupStep:
    name: str
    status: str  # ok | installed | skipped | failed | warning
    detail: str = ""


@dataclass
class SetupReport:
    platform: PlatformInfo
    steps: list[SetupStep] = field(default_factory=list)

    def add(self, name: str, status: str, detail: str = "") -> None:
        self.steps.append(SetupStep(name, status, detail))

    @property
    def ok(self) -> bool:
        return not any(s.status == "failed" for s in self.steps)


def plink2_asset(info: PlatformInfo, links: list[str]) -> str | None:
    """Choose the official PLINK 2 build for this OS / CPU from the download-page links."""
    if info.os == "windows":
        prefs = ["plink2_win_avx2_", "plink2_win64_"] if info.avx2 else ["plink2_win64_"]
    elif info.os == "darwin":
        prefs = ["plink2_mac_arm64_"] if info.arch == "aarch64" else (["plink2_mac_avx2_", "plink2_mac_"] if info.avx2 else ["plink2_mac_"])
    else:
        if info.arch == "aarch64":
            prefs = ["plink2_linux_arm64_", "plink2_linux_aarch64_"]
        else:
            vendor_amd = False
            try:
                with open("/proc/cpuinfo", encoding="utf-8") as fh:
                    vendor_amd = "AuthenticAMD" in fh.read(4096)
            except OSError:
                pass
            prefs = (["plink2_linux_amd_avx2_"] if vendor_amd else []) + (["plink2_linux_avx2_"] if info.avx2 else []) + ["plink2_linux_x86_64_"]
    for pref in prefs:
        for link in links:
            if link.rsplit("/", 1)[-1].startswith(pref):
                return link
    return None


def install_plink2_binary(project: Project, info: PlatformInfo) -> str:
    html = httpx.get(PLINK2_PAGE, follow_redirects=True, timeout=60).text
    links = re.findall(r'href="(https://s3\.amazonaws\.com/plink2-assets/[^"]+\.zip)"', html)
    url = plink2_asset(info, links)
    if url is None:
        raise RuntimeError(f"no PLINK 2 binary for {info.os}/{info.arch} on {PLINK2_PAGE}")
    archive = project.path(".efgpp", "downloads", url.rsplit("/", 1)[-1])
    download(url, archive)
    project.bin_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for member in z.namelist():
            if member.rsplit("/", 1)[-1] in ("plink2", "plink2.exe"):
                target = project.bin_dir / member.rsplit("/", 1)[-1]
                target.write_bytes(z.read(member))
                target.chmod(0o755)
                return url
    raise RuntimeError(f"{archive.name} contains no plink2 executable")


def _record_tool(project: Project, name: str) -> str | None:
    try:
        tool = resolve(project, name)
    except ToolNotFoundError:
        return None
    version = detect_version(tool)
    with Registry.open(project) as reg:
        reg.upsert("software", {"name": name, "environment": tool.env, "version": version,
                                "path": str(tool.path), "detected_at": utcnow()})
    return version or "unknown version"


def _create_env(project: Project, name: str, info: PlatformInfo, report: SetupReport, manager: tuple[str, object] | None) -> bool:
    if manager is None:
        report.add(f"environment {name}", "skipped", "no environment manager available")
        return False
    kind, exe = manager
    proc = pixi.create_env(project, exe, name, info) if kind == "pixi" else micromamba.create_env(project, exe, name)  # type: ignore[arg-type]
    if proc.returncode == 0:
        report.add(f"environment {name}", "installed", f"{kind}: .efgpp/envs/{name}")
        return True
    report.add(f"environment {name}", "failed", (proc.stderr or proc.stdout)[-400:])
    return False


def setup_data(project: Project, *, components: set[str] | None = None, dry_run: bool = False,
               progress: Callable[[str], None] | None = None) -> SetupReport:
    say = progress or (lambda _m: None)
    want = components or {"genetics", "core", "reporting"}
    if project.data.predicted.expression.enabled or "metaxcan" in want:
        want.add("metaxcan")
    info = detect()
    report = SetupReport(info)
    # 1-3: platform
    report.add("platform", "ok", f"{info.os}/{info.arch}, {info.cpu_count} CPUs, {info.memory_gb} GB, AVX2={info.avx2}")
    report.add("hpc", "ok", info.hpc_scheduler or "no scheduler detected (local execution)")
    # 4: python dependencies
    missing = [m for m in CORE_MODULES if importlib.util.find_spec(m) is None]
    report.add("python dependencies", "failed" if missing else "ok",
               f"missing: {missing} (run `uv sync`)" if missing else f"{len(CORE_MODULES)} core modules importable")
    if importlib.util.find_spec("anndata") is None:
        report.add("omics storage", "warning", "anndata not installed; omics stored as Parquet (pip install 'efgpp[omics]')")
    if dry_run:
        report.add("installation", "skipped", f"dry run; would install: {sorted(want)}")
        return report

    # 5: environment manager
    manager: tuple[str, object] | None = None
    if info.bioconda_supported and project.config.execution.environment_manager != "system":
        try:
            say("bootstrapping pixi")
            manager = ("pixi", pixi.bootstrap(project, info))
            report.add("pixi", "ok", str(manager[1]))
        except Exception as exc:  # noqa: BLE001
            mm = micromamba.find_micromamba(project)
            manager = ("micromamba", mm) if mm else None
            report.add("pixi", "warning", f"{exc}; fallback: {'micromamba' if mm else 'none'}")
    elif info.is_windows:
        report.add("bioconda", "warning", "Bioconda has no Windows builds; installing native binaries where available, "
                                         "use WSL2 or Docker for VEP, bcftools and MetaXcan")

    # 6-13: environments and tools
    for env in ("core", "genetics", "annotation", "metaxcan", "reporting"):
        if env in want and manager is not None:
            say(f"creating environment {env}")
            _create_env(project, env, info, report, manager)
    if "genetics" in want and _record_tool(project, "plink2") is None:
        try:
            say("installing PLINK 2 binary")
            url = install_plink2_binary(project, info)
            report.add("plink2", "installed", url)
        except Exception as exc:  # noqa: BLE001
            report.add("plink2", "failed", str(exc))
    if "metaxcan" in want and manager is not None:
        _install_metaxcan(project, report)

    # 14-16: test executables and save versions
    for tool in ("plink2", "plink", "bcftools", "tabix", "flashpca2", "snakemake", "multiqc",
                 *(("vep", "oc") if "annotation" in want else ()), *(("predixcan",) if "metaxcan" in want else ())):
        version = _record_tool(project, tool)
        if version:
            report.add(f"test {tool}", "ok", version)
        elif tool in ("plink2",) and "genetics" in want:
            report.add(f"test {tool}", "failed", "not runnable")
        else:
            report.add(f"test {tool}", "warning", "not available")
    # 17: lock file
    from efgpp.setup.lock import write_lock

    report.add("lock file", "ok", str(write_lock(project)))
    return report


def _install_metaxcan(project: Project, report: SetupReport) -> None:
    """Clone MetaXcan into the metaxcan environment and point the predixcan tool at it."""
    import subprocess

    target = project.envs_dir / "metaxcan" / "MetaXcan"
    if not target.exists():
        proc = subprocess.run(["git", "clone", "--depth", "1", "https://github.com/hakyimlab/MetaXcan", str(target)],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            report.add("metaxcan", "failed", proc.stderr[-300:])
            return
    script = target / "software" / "Predict.py"
    if script.exists():
        project.config.execution.tools = {**project.config.execution.tools,
                                          "predixcan": project.relative(script)}
        project.save_project_config()
        report.add("metaxcan", "installed", project.relative(script))
    else:
        report.add("metaxcan", "failed", "Predict.py not found in the MetaXcan checkout")


def setup_annotation(project: Project, *, build: str | None = None, progress: Callable[[str], None] | None = None) -> SetupReport:
    """Install VEP (+ matching cache/FASTA) and OpenCRAVAT, validate, record versions."""
    from efgpp.resources.manager import install

    report = setup_data(project, components={"annotation"}, progress=progress)
    build = build or project.resources.genome.build
    for name in ("genome", "vep"):
        try:
            res = install(project, name, build=build)
            report.add(f"resource {name}", res.status, f"{res.version} -> {res.path}")
        except Exception as exc:  # noqa: BLE001
            report.add(f"resource {name}", "failed", str(exc))
    from efgpp.setup.lock import write_lock

    write_lock(project)
    return report
