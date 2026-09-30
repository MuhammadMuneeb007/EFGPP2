"""Data modules: one readable, editable file per kind of data.

    columns.py     shared rules: ID columns, missing values, value types
    genotype.py    genotype format and sample-ID inference
    phenotype.py   which columns are phenotypes, and their type/coding
    covariates.py  covariate roles (sex, age, PCs, batch) and categorical vs numeric
    gwas.py        GWAS column recognition for GWASLab

To adapt a module to your data, copy it into the project and edit the copy:

    efgpp modules export covariates      # -> <project>/modules/covariates.py

EFGPP then uses <project>/modules/<name>.py instead of the built-in version (`efgpp modules
list` shows which one is active). Keep the function names and signatures.
"""

from __future__ import annotations

import importlib
import importlib.util
import shutil
import sys
from pathlib import Path
from types import ModuleType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from efgpp.project import Project

MODULES = ("columns", "genotype", "phenotype", "covariates", "gwas")
PACKAGE_DIR = Path(__file__).parent


def project_module_path(project: Project | None, name: str) -> Path | None:
    if project is None:
        return None
    path = project.path("modules", f"{name}.py")
    return path if path.exists() else None


def load(project: Project | None, name: str) -> ModuleType:
    """The project's copy of a module when present, else the built-in one."""
    if name not in MODULES:
        raise KeyError(f"unknown module {name!r}; choose from {', '.join(MODULES)}")
    override = project_module_path(project, name)
    if override is None:
        return importlib.import_module(f"efgpp.modules.{name}")
    spec = importlib.util.spec_from_file_location(f"efgpp_project_modules.{name}", override)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {override}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # dataclasses in the module need it registered
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(spec.name, None)
        raise
    return module


def export(project: Project, names: list[str], force: bool = False) -> list[Path]:
    """Copy built-in modules into <project>/modules/ for editing."""
    out = []
    target_dir = project.path("modules")
    target_dir.mkdir(parents=True, exist_ok=True)
    for name in names or MODULES:
        if name not in MODULES:
            raise KeyError(f"unknown module {name!r}; choose from {', '.join(MODULES)}")
        dest = target_dir / f"{name}.py"
        if dest.exists() and not force:
            raise FileExistsError(f"{dest} exists (use --force to overwrite your edits)")
        shutil.copy2(PACKAGE_DIR / f"{name}.py", dest)
        out.append(dest)
    return out
