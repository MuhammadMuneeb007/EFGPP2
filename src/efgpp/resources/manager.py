"""Installing, listing and checking reference resources.

Resources are installed into versioned directories (resources/<name>/<version>/) and
never overwritten silently: a version referenced by a frozen snapshot is immutable, and
`update --check` only reports newer upstream versions.
"""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from efgpp.data.references import PROVIDERS, FetchedResource
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths
from efgpp.resources.manifest import ResourceManifest


class ResourceProtectedError(RuntimeError):
    pass


@dataclass
class InstallResult:
    name: str
    version: str
    path: Path
    status: str  # installed | already-installed | replaced


def installed(project: Project) -> list[dict[str, Any]]:
    with Registry.open(project) as reg:
        return reg.rows("SELECT * FROM resources ORDER BY name, download_date")


def list_resources(project: Project) -> list[dict[str, Any]]:
    inst: dict[str, dict[str, Any]] = {}
    for r in installed(project):
        inst[r["name"]] = r
    rows = []
    names = sorted(set(type(project.resources).model_fields) | set(PROVIDERS))
    for name in names:
        cfg = getattr(project.resources, name, None)
        entry = inst.get(name)
        r = entry or {}
        rows.append({
            "name": name,
            "enabled": bool(cfg and cfg.enabled),
            "installed_version": r["version"] if r else None,
            "genome_build": r["genome_build"] if r else None,
            "local_path": r["local_path"] if r else (cfg.path if cfg is not None and getattr(cfg, "path", None) else None),
            "provider": "yes" if name in PROVIDERS or name == "vep" else "config-only",
        })
    return rows


def _register(project: Project, fetched: FetchedResource, dest: Path, local: Path) -> ResourceManifest:
    checksum = checksum_paths([dest / f.name for f in fetched.files if (dest / f.name).exists()] or [local])
    manifest = ResourceManifest(
        resource_id=f"RES_{fetched.name}_{fetched.version}", name=fetched.name, version=fetched.version,
        release_date=fetched.release_date, genome_build=fetched.genome_build, source=fetched.source,
        license=fetched.license, checksum=str(checksum), download_date=datetime.now(UTC).isoformat(timespec="seconds"),
        local_path=str(local), files=[f.name for f in fetched.files], metadata=fetched.metadata,
    )
    manifest.write(dest)
    with Registry.open(project) as reg:
        reg.upsert("resources", {
            "resource_id": manifest.resource_id, "name": manifest.name, "version": manifest.version,
            "release_date": manifest.release_date, "genome_build": manifest.genome_build, "source": manifest.source,
            "license": manifest.license, "checksum": manifest.checksum, "download_date": manifest.download_date,
            "local_path": manifest.local_path, "metadata": manifest.metadata,
        })
    return manifest


def install(project: Project, name: str, *, build: str | None = None, force: bool = False,
            **options: Any) -> InstallResult:
    """Install a reference resource. `options` are provider-specific (e.g. dataset=, url=, progress=)."""
    from efgpp.data.snapshots import protected_resource_paths

    target_build = project.config.defaults.target_build
    if build is not None and build != target_build:
        raise ValueError(f"this project's target build is {target_build}; resources are installed for "
                         f"{target_build} only (requested {build}). Change defaults.target_build in project.yaml "
                         f"to use {build}.")
    build = build or target_build
    if name == "vep":
        return install_vep_cache(project, build=build)
    if name not in PROVIDERS:
        raise KeyError(f"no provider for {name!r}; available: {', '.join(sorted(PROVIDERS))}")
    provider = PROVIDERS[name](project)
    provider.options = options  # type: ignore[attr-defined]
    staging = project.path(".efgpp", "downloads", f"{name}-{datetime.now(UTC):%Y%m%d%H%M%S}")
    staging.mkdir(parents=True, exist_ok=True)
    try:
        fetched = provider.fetch(build=build, staging=staging)
        dest = provider.cache(re.sub(r"[^\w.\-+]", "_", fetched.version))
        protected = protected_resource_paths(project)
        status = "installed"
        if (dest / "manifest.yaml").exists():
            existing = ResourceManifest.read(dest)
            if existing.local_path in protected:
                raise ResourceProtectedError(
                    f"{name} {fetched.version} is used by a frozen snapshot and cannot be overwritten")
            if not force:
                return InstallResult(name, fetched.version, Path(existing.local_path), "already-installed")
            status = "replaced"
        for f in fetched.files:
            target = dest / f.relative_to(staging) if f.is_relative_to(staging) else dest / f.name
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.is_dir() and not target.is_symlink():
                shutil.rmtree(target)
            elif target.exists():
                target.unlink()
            shutil.move(str(f), target)
        local = dest / fetched.local_path.relative_to(staging) if fetched.local_path.is_relative_to(staging) else dest
        _register(project, fetched, dest, local)
        if name == "genome":
            _configure_genome(project, fetched, dest)
        return InstallResult(name, fetched.version, local, status)
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _configure_genome(project: Project, fetched: FetchedResource, dest: Path) -> None:
    changed = False
    for f in fetched.files:
        if f.name.endswith(".over.chain.gz"):
            key = "GRCh37->GRCh38" if f.name.startswith("hg19To") else "GRCh38->GRCh37"
            project.resources.liftover.chain_files[key] = project.relative(dest / f.name)
            changed = True
    if changed:
        project.resources.liftover.chain_files = dict(project.resources.liftover.chain_files)
        project.save_resources_config()


def install_vep_cache(project: Project, *, build: str | None = None) -> InstallResult:
    """Install the VEP cache (and FASTA) matching the requested build with vep_install."""
    from efgpp.data.provenance import run_tool
    from efgpp.setup.tools import resolve

    build = build or project.resources.genome.build
    cfg = project.resources.vep
    cache_dir = project.resolve(cfg.cache) if cfg.cache != "auto" else project.resource_root / "vep"
    cache_dir.mkdir(parents=True, exist_ok=True)
    vep = resolve(project, "vep")
    installer = next((d / n for d in [vep.path.parent] for n in ("vep_install", "vep_install.pl", "INSTALL.pl")
                      if (d / n).exists()), None)
    if installer is None:
        raise RuntimeError("vep_install not found next to the vep executable")
    args = ["-a", "cf", "-s", cfg.species, "-y", build, "-c", str(cache_dir), "--NO_UPDATE", "--CONVERT"]
    if cfg.release:
        args += ["-r", str(cfg.release)]
    run_tool(project, "vep_install", args, tool=type(vep)(name="vep_install", path=installer, env=vep.env,
                                                           env_prefix=vep.env_prefix))
    releases = sorted(cache_dir.glob(f"{cfg.species}/*_{build}"))
    if not releases:
        raise RuntimeError(f"VEP cache for {build} not found under {cache_dir} after installation")
    release = releases[-1].name.split("_")[0]
    fetched = FetchedResource("vep", f"{release}_{build}", [], releases[-1], "Ensembl VEP cache",
                              "Ensembl terms of use", genome_build=build)
    _register(project, fetched, releases[-1], releases[-1])
    return InstallResult("vep", fetched.version, releases[-1], "installed")


def check_updates(project: Project) -> list[dict[str, Any]]:
    """Compare installed versions with upstream. Never updates anything."""
    rows = []
    latest: dict[str, dict[str, Any]] = {}
    for r in installed(project):
        latest[r["name"]] = r
    for name, r in latest.items():
        cls = PROVIDERS.get(name)
        upstream = None
        if cls is not None:
            try:
                upstream = cls(project).version()
            except Exception as exc:  # noqa: BLE001 - network problems are reported, not fatal
                upstream = f"unavailable ({exc.__class__.__name__})"
        rows.append({"name": name, "installed": r["version"], "upstream": upstream,
                     "action": "none (EFGPP never updates resources automatically)"})
    return rows

