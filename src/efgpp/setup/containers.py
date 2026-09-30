"""Container runtimes (Apptainer on Linux/HPC, Docker on desktop/CI).

Containers are optional: they are used only for components that cannot run natively
(e.g. Linux-only annotation tools on Windows) and only when enabled in project.yaml.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass

from efgpp.project import Project
from efgpp.setup.platform import PlatformInfo

# Official images for Linux-only components.
IMAGES = {
    "vep": "docker://ensemblorg/ensembl-vep:latest",
    "oc": "docker://karchinlab/opencravat:latest",
}


@dataclass
class ContainerRuntime:
    name: str
    path: str


def runtime(project: Project, info: PlatformInfo) -> ContainerRuntime | None:
    cfg = project.config.execution.containers
    preferred = cfg.linux_hpc if info.os == "linux" and not info.wsl else cfg.desktop
    for name in (preferred, "apptainer", "singularity", "docker", "podman"):
        if name and name != "none" and (p := shutil.which(name)):
            return ContainerRuntime(name, p)
    return None


def image_digest(rt: ContainerRuntime, image: str) -> str | None:
    """Resolved digest of a pulled image (recorded in efgpp.lock.yaml)."""
    if rt.name not in ("docker", "podman"):
        return None
    ref = image.removeprefix("docker://")
    proc = subprocess.run([rt.path, "image", "inspect", "--format", "{{index .RepoDigests 0}}", ref],
                          capture_output=True, text=True)
    return proc.stdout.strip() or None if proc.returncode == 0 else None
