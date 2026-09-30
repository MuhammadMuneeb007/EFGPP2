"""`efgpp setup data` / `efgpp setup annotation`: detect, install, test and lock.

Scientific tools go into small isolated environments (Pixi preferred; mamba, micromamba or
conda as fallback) under software/envs/<purpose>. Any tool still missing afterwards is
downloaded from its official source by efgpp.setup.installers (binaries, source builds,
Python environments, containers) into the project or the shared per-user folder.
"""

from __future__ import annotations

import importlib.util
from collections.abc import Callable
from dataclasses import dataclass, field

from efgpp.data.registry import Registry, utcnow
from efgpp.project import Project
from efgpp.setup import micromamba, pixi
from efgpp.setup.installers import COMPONENT_TOOLS, install_plink2, install_tools
from efgpp.setup.platform import PlatformInfo, detect
from efgpp.setup.tools import ToolNotFoundError, detect_version, resolve

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


def install_plink2_binary(project: Project, info: PlatformInfo) -> str:
    """Official PLINK 2 binary into <project>/software/bin (kept for callers of the old API)."""
    return install_plink2(project.software_dir, info, lambda _m: None).method


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
    proc = pixi.create_env(project, exe, name, info) if kind == "pixi" else micromamba.create_env(project, exe, name, kind)  # type: ignore[arg-type]
    if proc.returncode == 0:
        report.add(f"environment {name}", "installed", f"{kind}: software/envs/{name}")
        return True
    if kind == "pixi":  # retry with the fastest conda-family tool before giving up
        fallback = _conda_family(project)
        if fallback is not None:
            report.add(f"environment {name}", "warning", f"pixi failed; retrying with {fallback[0]}")
            return _create_env(project, name, info, report, fallback)
    report.add(f"environment {name}", "failed", (proc.stderr or proc.stdout)[-400:])
    return False


def _conda_family(project: Project) -> tuple[str, object] | None:
    pref = project.config.execution.fallback_environment_manager
    if pref == "system":
        return None
    return micromamba.find_conda_tool(project, pref)


def setup_data(project: Project, *, components: set[str] | None = None, dry_run: bool = False,
               shared: bool = False, progress: Callable[[str], None] | None = None) -> SetupReport:
    say = progress or (lambda _m: None)
    want = components or {"genetics", "reporting", "gwas"}
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
    chosen = project.config.execution.environment_manager
    if info.bioconda_supported and chosen == "pixi":
        try:
            say("bootstrapping pixi")
            manager = ("pixi", pixi.bootstrap(project, info))
            report.add("pixi", "ok", str(manager[1]))
        except Exception as exc:  # noqa: BLE001
            manager = _conda_family(project)
            report.add("pixi", "warning", f"{exc}; fallback: {manager[0] if manager else 'none'}")
    elif info.bioconda_supported and chosen != "system":
        manager = micromamba.find_conda_tool(project, chosen)
        if manager is None or manager[0] != chosen:
            report.add(chosen, "warning", f"{chosen} not found; using {manager[0] if manager else 'none'}")
        else:
            report.add(chosen, "ok", str(manager[1]))
    elif info.is_windows:
        report.add("bioconda", "warning", "Bioconda has no Windows builds; installing native binaries where available, "
                                         "use WSL2 or Docker for VEP, bcftools and MetaXcan")

    # 6-13: environments and tools
    for env in ("genetics", "annotation", "metaxcan", "reporting"):
        if env in want and manager is not None:
            say(f"creating environment {env}")
            _create_env(project, env, info, report, manager)
    # Anything still missing: download it from its official source into the project
    # (or the shared per-user folder), so setup never depends on conda alone.
    missing = [tool for comp in sorted(want) for tool in COMPONENT_TOOLS.get(comp, [])]
    for tool, status, detail in install_tools(project, missing, shared=shared, progress=say):
        if status != "ok":
            # Only PLINK 2 is essential; other tools are reported and the plan marks their steps blocked.
            report.add(tool, "warning" if status == "failed" and tool != "plink2" else status, detail)

    # 14-16: test executables and save versions
    for tool in ("plink2", "plink", "bcftools", "tabix", "flashpca2", "multiqc", "gwaslab",
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
