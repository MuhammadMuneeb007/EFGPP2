"""Conda-family fallback (micromamba, mamba, conda) for environments in workflow/envs/*.yaml.

Pixi is preferred; when it cannot be used, the fastest conda-family tool available is used:
micromamba, then mamba, then conda (conda >= 23.10 also uses the libmamba solver).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from efgpp.project import Project

ORDER = ("micromamba", "mamba", "conda")


def find_conda_tool(project: Project, preferred: str | None = None) -> tuple[str, Path] | None:
    """(kind, executable) of the preferred or fastest available conda-family tool."""
    order = [preferred, *ORDER] if preferred in ORDER else list(ORDER)
    for kind in dict.fromkeys(order):
        for name in (kind, f"{kind}.exe", f"{kind}.bat"):
            local = project.bin_dir / name
            if local.exists():
                return kind, local
        found = shutil.which(kind)  # type: ignore[arg-type]
        if found:
            return kind, Path(found)  # type: ignore[return-value]
    return None


def find_micromamba(project: Project) -> Path | None:
    hit = find_conda_tool(project, "micromamba")
    return hit[1] if hit and hit[0] == "micromamba" else None


def create_command(kind: str, exe: Path, prefix: Path, spec: Path) -> list[str]:
    if kind == "micromamba":
        return [str(exe), "create", "-y", "-p", str(prefix), "-f", str(spec)]
    # mamba 2.x and conda share the `env create` interface
    return [str(exe), "env", "create", "-y", "-p", str(prefix), "-f", str(spec)]


def create_env(project: Project, exe: Path, name: str, kind: str = "micromamba") -> subprocess.CompletedProcess[str]:
    spec = project.path("workflow", "envs", f"{name}.yaml")
    prefix = project.envs_dir / name
    return subprocess.run(create_command(kind, exe, prefix, spec), capture_output=True, text=True)
