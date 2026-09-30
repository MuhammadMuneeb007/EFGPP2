"""Locating scientific executables inside EFGPP's isolated environments.

Search order for a tool: explicit override in project.yaml (execution.tools) ->
software/bin -> the tool's own environment (software/envs/<env>) -> legacy .efgpp/ ->
$EFGPP_TOOLS_HOME (only when set) -> PATH.
Each tool belongs to exactly one environment; there is no monolithic environment.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from efgpp.project import Project


def shared_root() -> Path | None:
    """Optional install root shared by several projects (`--shared`): the directory named by
    EFGPP_TOOLS_HOME. There is no default: without it, everything installs into the project."""
    override = os.environ.get("EFGPP_TOOLS_HOME")
    return Path(override).expanduser() if override else None


@dataclass(frozen=True)
class ToolSpec:
    name: str
    env: str  # core | genetics | annotation | metaxcan | reporting | system
    group: str  # doctor grouping
    executables: tuple[str, ...]
    version_args: tuple[str, ...] | None = ("--version",)
    conda_package: str | None = None
    optional: bool = False
    description: str = ""


TOOLS: dict[str, ToolSpec] = {
    t.name: t
    for t in (
        ToolSpec("snakemake", "core", "CORE", ("snakemake",), conda_package="snakemake-minimal", optional=True,
                 description="optional: the built-in executor runs plans locally without it"),
        ToolSpec("plink2", "genetics", "GENETICS", ("plink2",), conda_package="plink2",
                 description="genotype QC, PCA, relatedness, format conversion"),
        ToolSpec("plink", "genetics", "GENETICS", ("plink", "plink1.9"), conda_package="plink",
                 optional=True, description="PLINK 1.9, only for runs of homozygosity"),
        ToolSpec("bcftools", "genetics", "GENETICS", ("bcftools",), conda_package="bcftools"),
        ToolSpec("tabix", "genetics", "GENETICS", ("tabix",), conda_package="htslib"),
        ToolSpec("bgzip", "genetics", "GENETICS", ("bgzip",), conda_package="htslib"),
        ToolSpec("flashpca2", "genetics", "GENETICS", ("flashpca", "flashpca2"),
                 conda_package="flashpca", optional=True),
        ToolSpec("vep", "annotation", "ANNOTATION", ("vep",), ("--help",),
                 conda_package="ensembl-vep"),
        ToolSpec("oc", "annotation", "ANNOTATION", ("oc",), ("version",),
                 conda_package="open-cravat", optional=True),
        ToolSpec("predixcan", "metaxcan", "PREDICTED OMICS", ("predixcan", "Predict.py"), None, optional=True,
                 description="MetaXcan/PrediXcan Predict.py"),
        ToolSpec("multiqc", "reporting", "REPORTING", ("multiqc",), conda_package="multiqc",
                 optional=True),
        ToolSpec("gwaslab", "gwaslab", "GWAS", ("gwaslab-python",),
                 ("-c", "import importlib.metadata as m; print(m.version('gwaslab'))"), conda_package="gwaslab",
                 optional=True, description="GWASLab (own environment) for GWAS summary statistics"),
        ToolSpec("pixi", "system", "ENVIRONMENTS", ("pixi",), optional=True),
        ToolSpec("micromamba", "system", "ENVIRONMENTS", ("micromamba",), optional=True),
        ToolSpec("sbatch", "system", "HPC", ("sbatch",), ("--version",), optional=True),
        ToolSpec("apptainer", "system", "HPC", ("apptainer", "singularity"), optional=True),
        ToolSpec("docker", "system", "HPC", ("docker",), optional=True),
    )
}

ENV_NAMES = ("core", "genetics", "annotation", "metaxcan", "reporting")


@dataclass
class ResolvedTool:
    name: str
    path: Path
    env: str
    env_prefix: Path | None = None
    launcher: list[str] = field(default_factory=list)  # e.g. [python] for .py scripts

    def command(self, *args: str) -> list[str]:
        return [*self.launcher, str(self.path), *args]

    def environ(self) -> dict[str, str]:
        """Process environment with the tool's own environment activated."""
        env = dict(os.environ)
        if self.env_prefix is not None:
            bins = [str(d) for d in env_bin_dirs(self.env_prefix) if d.exists()]
            env["PATH"] = os.pathsep.join([*bins, env.get("PATH", "")])
            env["CONDA_PREFIX"] = str(self.env_prefix)
        return env


class ToolNotFoundError(RuntimeError):
    pass


def env_prefixes(project: Project, env: str) -> list[Path]:
    out = []
    for envs in (project.envs_dir, *(d / "envs" for d in project.legacy_software_dirs)):
        base = envs / env
        out += [base / ".pixi" / "envs" / "default", base]
    return out


def env_bin_dirs(prefix: Path) -> list[Path]:
    if sys.platform == "win32":
        return [prefix / "Library" / "bin", prefix / "Scripts", prefix / "bin", prefix]
    return [prefix / "bin"]


def _candidate_names(exe: str) -> list[str]:
    if sys.platform == "win32" and not exe.endswith(".py"):
        return [exe + ".exe", exe + ".bat", exe + ".cmd", exe]
    return [exe]


def resolve(project: Project | None, name: str) -> ResolvedTool:
    spec = TOOLS.get(name) or ToolSpec(name, "system", "OTHER", (name,))
    if project is not None:
        override = project.config.execution.tools.get(name)
        if override:
            p = project.resolve(override)
            if not p.exists():
                raise ToolNotFoundError(f"{name}: configured path {p} does not exist")
            return _finish(ResolvedTool(name, p, spec.env if p.suffix == ".py" else "override"), project)
        search: list[tuple[Path, str, Path | None]] = [(project.bin_dir, "project-bin", None)]
        search += [(d / "bin", "project-bin", None) for d in project.legacy_software_dirs]
        for env_prefix in env_prefixes(project, spec.env):
            search += [(d, spec.env, env_prefix) for d in env_bin_dirs(env_prefix)]
        for directory, env, prefix in search:
            for exe in spec.executables:
                for cand in _candidate_names(exe):
                    p = directory / cand
                    if p.is_file():
                        return _finish(ResolvedTool(name, p, env, prefix), project)
    shared_home = shared_root()
    if shared_home is not None:
        shared = shared_home / "bin"
        for exe in spec.executables:
            for cand in _candidate_names(exe):
                if (shared / cand).is_file():
                    return _finish(ResolvedTool(name, shared / cand, "shared-bin"), project)
    for exe in spec.executables:
        found = shutil.which(exe)
        if found:
            return _finish(ResolvedTool(name, Path(found), "PATH"), project)
    raise ToolNotFoundError(
        f"{name} not found (looked in software/bin, software/envs and PATH); "
        f"run `efgpp setup data` or set execution.tools.{name} in project.yaml"
    )


def _finish(tool: ResolvedTool, project: Project | None) -> ResolvedTool:
    """Python scripts run with the interpreter of their own environment when it exists."""
    if tool.path.suffix == ".py":
        tool.launcher = [sys.executable]
        spec = TOOLS.get(tool.name)
        if project is not None and spec is not None:
            for prefix in env_prefixes(project, spec.env):
                for d in env_bin_dirs(prefix):
                    for exe in ("python.exe", "python"):
                        if (d / exe).is_file():
                            tool.launcher = [str(d / exe)]
                            tool.env_prefix = prefix
                            return tool
    return tool


def available(project: Project | None, name: str) -> bool:
    try:
        resolve(project, name)
    except ToolNotFoundError:
        return False
    return True


_VERSION_RE = re.compile(r"v?(\d+(?:\.\d+)+[\w.\-]*(?:\s+\([^)]*\))?)")


def detect_version(tool: ResolvedTool, timeout: int = 30) -> str | None:
    spec = TOOLS.get(tool.name)
    args = spec.version_args if spec else ("--version",)
    if args is None:
        return None
    try:
        proc = subprocess.run(
            tool.command(*args), capture_output=True, text=True, timeout=timeout, env=tool.environ()
        )
    except (OSError, subprocess.SubprocessError):
        return None
    text = (proc.stdout or "") + "\n" + (proc.stderr or "")
    if tool.name == "vep":
        m = re.search(r"ensembl-vep\s*:\s*([\w.]+)", text)
        return m.group(1) if m else None
    for line in text.splitlines():
        m = _VERSION_RE.search(line)
        if m:
            return m.group(1).strip()
    return None
