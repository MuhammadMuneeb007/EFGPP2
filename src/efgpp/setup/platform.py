"""Operating system, architecture, CPU features and HPC environment detection."""

from __future__ import annotations

import os
import platform
import shutil
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import psutil


@dataclass
class PlatformInfo:
    os: str  # linux | darwin | windows
    arch: str  # x86_64 | aarch64
    python: str
    cpu_count: int
    memory_gb: float
    avx2: bool
    wsl: bool
    hpc_scheduler: str | None
    container_runtimes: list[str] = field(default_factory=list)

    @property
    def is_windows(self) -> bool:
        return self.os == "windows"

    @property
    def bioconda_supported(self) -> bool:
        # Bioconda publishes linux-64, linux-aarch64, osx-64 and osx-arm64 packages only.
        return self.os in ("linux", "darwin")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _arch() -> str:
    m = platform.machine().lower()
    if m in ("amd64", "x86_64", "x64"):
        return "x86_64"
    if m in ("arm64", "aarch64"):
        return "aarch64"
    return m


def _has_avx2() -> bool:
    try:
        text = Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="ignore")
        return " avx2 " in text.replace("\n", " ")
    except OSError:
        pass
    if sys.platform == "darwin":
        try:
            import subprocess

            out = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.leaf7_features"], capture_output=True, text=True, timeout=5
            ).stdout
            return "AVX2" in out
        except (OSError, subprocess.SubprocessError):
            return False
    if sys.platform == "win32":
        try:
            import ctypes

            PF_AVX2_INSTRUCTIONS_AVAILABLE = 40
            return bool(ctypes.windll.kernel32.IsProcessorFeaturePresent(PF_AVX2_INSTRUCTIONS_AVAILABLE))  # type: ignore[attr-defined]
        except (AttributeError, OSError):
            return False
    return False


def detect_scheduler() -> str | None:
    env = os.environ
    if "SLURM_JOB_ID" in env or "SLURM_CLUSTER_NAME" in env or shutil.which("sbatch"):
        return "slurm"
    if "PBS_JOBID" in env or shutil.which("qsub") and shutil.which("pbsnodes"):
        return "pbs"
    if "LSB_JOBID" in env or shutil.which("bsub"):
        return "lsf"
    if "SGE_ROOT" in env:
        return "sge"
    return None


def detect() -> PlatformInfo:
    os_name = {"win32": "windows", "darwin": "darwin"}.get(sys.platform, "linux")
    wsl = os_name == "linux" and "microsoft" in platform.release().lower()
    runtimes = [r for r in ("apptainer", "singularity", "docker", "podman") if shutil.which(r)]
    return PlatformInfo(
        os=os_name,
        arch=_arch(),
        python=platform.python_version(),
        cpu_count=os.cpu_count() or 1,
        memory_gb=round(psutil.virtual_memory().total / 1024**3, 1),
        avx2=_has_avx2(),
        wsl=wsl,
        hpc_scheduler=detect_scheduler(),
        container_runtimes=runtimes,
    )
