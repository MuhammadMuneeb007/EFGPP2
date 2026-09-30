"""`efgpp doctor`: what is installed, what is missing, and what that means."""

from __future__ import annotations

import importlib.util
from dataclasses import dataclass

from efgpp.data.references.variants import AlphaGenomeProvider
from efgpp.project import Project
from efgpp.setup.platform import detect
from efgpp.setup.tools import TOOLS, ToolNotFoundError, detect_version, resolve


@dataclass
class DoctorLine:
    group: str
    name: str
    status: str  # ok | missing | warn | info
    detail: str = ""


def run_doctor(project: Project | None) -> list[DoctorLine]:
    lines: list[DoctorLine] = []
    info = detect()
    lines.append(DoctorLine("SYSTEM", "platform", "info", f"{info.os}/{info.arch}, Python {info.python}, AVX2={info.avx2}"))
    for mod, label in (("polars", "Polars"), ("pyarrow", "PyArrow"), ("duckdb", "DuckDB"), ("pandera", "Pandera"),
                       ("pydantic", "Pydantic"), ("anndata", "AnnData"), ("zarr", "Zarr")):
        ok = importlib.util.find_spec(mod) is not None
        lines.append(DoctorLine("CORE", label, "ok" if ok else ("warn" if mod in ("anndata", "zarr") else "missing")))
    for name, spec in TOOLS.items():
        if spec.group == "ENVIRONMENTS" and name == "micromamba":
            continue
        try:
            tool = resolve(project, name)
            version = detect_version(tool)
            lines.append(DoctorLine(spec.group, name, "ok", f"{version or ''} [{tool.env}]".strip()))
        except ToolNotFoundError:
            lines.append(DoctorLine(spec.group, name, "warn" if spec.optional else "missing", spec.description))
    if project is not None:
        from efgpp.data.registry import Registry

        with Registry.open(project) as reg:
            res = {r["name"]: r for r in reg.rows("SELECT name, version, local_path FROM resources")}
        for name in ("genome", "vep", "clinvar", "alphamissense", "predictdb", "predictdb-gtex-v8-expression",
                     "predictdb-gtex-v8-splicing", "omicspred", "mimosa", "hibag"):
            cfg = getattr(project.resources, name.replace("-", "_"), None)
            if name in res:
                lines.append(DoctorLine("RESOURCES", name, "ok", str(res[name]["version"])))
            elif cfg is not None and cfg.enabled:
                lines.append(DoctorLine("RESOURCES", name, "missing", f"efgpp resources install {name}"))
        ag = AlphaGenomeProvider(project).describe()
        if project.resources.alphagenome.enabled:
            if not ag["api_key_configured"]:
                lines.append(DoctorLine("RESOURCES", "AlphaGenome API key", "warn",
                                        "not configured; atlas mode remains available"))
            else:
                lines.append(DoctorLine("RESOURCES", "AlphaGenome API key", "ok"))
    lines.append(DoctorLine("HPC", "scheduler", "ok" if info.hpc_scheduler else "info", info.hpc_scheduler or "none detected"))
    if info.is_windows:
        lines.append(DoctorLine("SYSTEM", "Windows", "info",
                                "core + PLINK 2 run natively; VEP/bcftools/MetaXcan need WSL2 or Docker"))
    order = list(dict.fromkeys(line.group for line in lines))
    return sorted(lines, key=lambda line: order.index(line.group))  # stable: one block per group
