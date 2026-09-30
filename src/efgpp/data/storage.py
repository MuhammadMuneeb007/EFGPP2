"""Source storage policies: copy, link, reference, auto.

Sources are never modified. `reference` leaves files where they are and records their
absolute path, size, checksum and modification time; `copy`/`link` place them under
data/<origin>/<modality>/<source_id>/.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from efgpp.constants import GENOTYPE_FORMATS, Modality, Origin, StorageMode
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths, iter_files, total_size

GB = 1024**3


@dataclass
class PlacedSource:
    mode: StorageMode
    record_path: Path  # what the artifact's `path` points at (file, directory or PLINK prefix)
    members: list[Path]
    size: int
    checksum: str
    mtime: float
    notes: list[str] = field(default_factory=list)


def resolve_mode(
    requested: StorageMode,
    *,
    fmt: str,
    size: int,
    project: Project,
) -> StorageMode:
    if requested != StorageMode.AUTO:
        return requested
    policy = project.config.storage.external_source_policy
    if policy != StorageMode.AUTO:
        return policy
    # Large genotype cohorts are referenced in place; small tables are copied.
    if fmt in GENOTYPE_FORMATS:
        return StorageMode.REFERENCE
    if size < project.config.storage.copy_threshold_gb * GB:
        return StorageMode.COPY
    return StorageMode.REFERENCE


def _place_one(src: Path, dest: Path, mode: StorageMode) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        if dest.is_dir() and not dest.is_symlink():
            shutil.rmtree(dest)
        else:
            dest.unlink()
    if mode == StorageMode.COPY:
        if src.is_dir():
            shutil.copytree(src, dest)
        else:
            shutil.copy2(src, dest)
    else:
        os.symlink(src.resolve(), dest, target_is_directory=src.is_dir())


def place_source(
    project: Project,
    *,
    modality: Modality,
    origin: Origin,
    source_id: str,
    record_path: Path,
    members: list[Path],
    requested: StorageMode,
    fmt: str,
) -> PlacedSource:
    """Apply the storage policy and fingerprint the (placed) source."""
    record_path = record_path.resolve()
    members = [m.resolve() for m in members]
    size = total_size(members)
    mode = resolve_mode(requested, fmt=fmt, size=size, project=project)
    notes: list[str] = []

    if mode in (StorageMode.COPY, StorageMode.LINK):
        dest_dir = project.artifact_dir(origin, modality, source_id)
        placed = []
        try:
            for m in members:
                _place_one(m, dest_dir / m.name, mode)
                placed.append(dest_dir / m.name)
        except OSError as exc:
            if mode != StorageMode.LINK:
                raise
            # Symlinks need extra privileges on Windows; fall back without copying data.
            notes.append(f"symbolic link failed ({exc.__class__.__name__}); stored by reference")
            mode = StorageMode.REFERENCE
        else:
            members = placed
            record_path = dest_dir / record_path.name
    limit = int(project.config.storage.full_checksum_limit_gb * GB)
    checksum = checksum_paths(members, full_limit_bytes=limit)
    mtime = max((f.stat().st_mtime for m in members for f in iter_files(m)), default=0.0)
    return PlacedSource(mode, record_path, members, size, str(checksum), mtime, notes)
