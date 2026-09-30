"""Micromamba fallback for environments defined in workflow/envs/*.yaml."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from efgpp.project import Project


def find_micromamba(project: Project) -> Path | None:
    for name in ("micromamba", "micromamba.exe"):
        local = project.bin_dir / name
        if local.exists():
            return local
    found = shutil.which("micromamba")
    return Path(found) if found else None


def create_env(project: Project, micromamba: Path, name: str) -> subprocess.CompletedProcess[str]:
    spec = project.path("workflow", "envs", f"{name}.yaml")
    prefix = project.envs_dir / name
    return subprocess.run(
        [str(micromamba), "create", "-y", "-p", str(prefix), "-f", str(spec)],
        capture_output=True, text=True,
    )
