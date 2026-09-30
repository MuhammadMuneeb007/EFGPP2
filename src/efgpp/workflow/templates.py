"""Scientific environment definitions written at `efgpp init` (workflow/envs/*.yaml).

One small environment per purpose - never one enormous environment.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import yaml

if TYPE_CHECKING:
    from efgpp.project import Project

CHANNELS = ["conda-forge", "bioconda"]

# Conda/Bioconda package specifications per environment.
ENVIRONMENTS: dict[str, dict[str, str]] = {
    "genetics": {"plink2": "*", "plink": "*", "bcftools": "*", "htslib": "*"},
    "annotation": {"ensembl-vep": "*", "htslib": "*", "open-cravat": "*"},
    "metaxcan": {"python": "3.10.*", "numpy": "*", "scipy": "*", "pandas": "*", "cyvcf2": "*",
                 "bgen-reader": "*", "h5py": "*", "git": "*"},
    "reporting": {"multiqc": "*"},
}

# Packages without a build for a platform are dropped from that platform's manifest.
UNAVAILABLE = {
    "osx-arm64": {"flashpca", "plink"},
    "linux-aarch64": {"flashpca"},
}

def conda_env_yaml(name: str) -> dict[str, object]:
    deps = [f"{pkg}{'' if ver == '*' else ver}" for pkg, ver in ENVIRONMENTS[name].items()]
    return {"name": f"efgpp-{name}", "channels": CHANNELS, "dependencies": deps}


def pixi_manifest(name: str, platform: str) -> str:
    skip = UNAVAILABLE.get(platform, set())
    lines = [
        "[workspace]",
        f'name = "efgpp-{name}"',
        f"channels = {CHANNELS!r}".replace("'", '"'),
        f'platforms = ["{platform}"]',
        "",
        "[dependencies]",
    ]
    lines += [f'{pkg} = "{ver}"' for pkg, ver in ENVIRONMENTS[name].items() if pkg not in skip]
    return "\n".join(lines) + "\n"


def write_workflow_templates(project: Project) -> None:
    wf = project.path("workflow")
    (wf / "envs").mkdir(parents=True, exist_ok=True)
    for name in ENVIRONMENTS:
        (wf / "envs" / f"{name}.yaml").write_text(
            yaml.safe_dump(conda_env_yaml(name), sort_keys=False), encoding="utf-8")
