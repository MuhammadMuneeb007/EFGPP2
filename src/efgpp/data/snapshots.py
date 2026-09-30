"""Immutable DataSnapshots: the contract handed to the Representation layer.

`efgpp data freeze --name data_v1` writes

    snapshots/data_v1.yaml       artifacts, hashes, software, resources, builds, QC config,
                                 predicted modalities, timeline, availability summary
    snapshots/data_v1/           frozen registry tables (participants, aliases, events,
                                 biospecimens, assays, phenotypes, availability)

and the Representation layer asks the snapshot for artifacts instead of searching folders.
"""

from __future__ import annotations

import json
import os
import shutil
import stat
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import polars as pl
import yaml

from efgpp import __version__
from efgpp.constants import USABLE_STATUSES, ArtifactStatus, Modality, Origin
from efgpp.data import availability as availability_mod
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.registry import Registry, utcnow
from efgpp.data.timeline import timeline_summary
from efgpp.project import LOCK_FILE, Project
from efgpp.resources.checksums import Checksum, checksum_paths, text_sha256

FROZEN_TABLES = ("participants", "sample_aliases", "events", "biospecimens", "assays",
                 "phenotype_definitions", "phenotype_observations")


class SnapshotExistsError(FileExistsError):
    pass


def _artifact_files(art: Artifact) -> list[Path]:
    members = art.metadata.get("members")
    if members:
        return [Path(m) for m in members]
    p = Path(art.path)
    if art.format in ("pgen", "bed"):
        from efgpp.data.genotype.formats import resolve_fileset

        return resolve_fileset(p, art.format).members
    return [p]


def verify_artifact(art: Artifact, full_limit_bytes: int) -> str | None:
    """Return a problem description, or None when the artifact is intact."""
    files = _artifact_files(art)
    missing = [str(f) for f in files if not f.exists()]
    if missing:
        return f"missing files: {missing[:3]}"
    if not art.checksum:
        return None
    expected = Checksum.parse(art.checksum)
    actual = checksum_paths(files, full_limit_bytes=full_limit_bytes)
    if actual.value != expected.value:
        return f"checksum changed ({expected.algorithm})"
    return None


def freeze(project: Project, name: str, *, include_failed: bool = False, verify: bool = True) -> Path:
    path = project.path("snapshots", f"{name}.yaml")
    if path.exists():
        raise SnapshotExistsError(f"snapshot {name!r} already exists; snapshots are immutable - choose a new name")
    statuses = set(USABLE_STATUSES) | ({ArtifactStatus.QC_FAIL} if include_failed else set())
    limit = int(project.config.storage.full_checksum_limit_gb * 1024**3)
    with Registry.open(project) as reg:
        if reg.scalar("SELECT count(*) FROM snapshots WHERE name = ?", [name]):
            raise SnapshotExistsError(f"snapshot {name!r} is already registered")
        av = availability_mod.compute(project, reg)
        availability_mod.write(project, reg, av)
        store = ArtifactStore(reg)
        artifacts = store.find(statuses=statuses)
        problems = {a.artifact_id: verify_artifact(a, limit) for a in artifacts} if verify else {}
        bad = {k: v for k, v in problems.items() if v}
        if bad:
            raise RuntimeError(f"cannot freeze: artifacts changed on disk since registration: {bad}")
        snapshot_id = reg.next_id("SNAP", width=4)
        frozen_dir = project.path("snapshots", name)
        if frozen_dir.exists():
            raise SnapshotExistsError(f"snapshots/{name}/ already exists; choose a new name")
        frozen_dir.mkdir(parents=True)
        for table in FROZEN_TABLES:
            reg.execute(f"COPY (SELECT * FROM {table}) TO '{(frozen_dir / f'{table}.parquet').as_posix()}' (FORMAT PARQUET)")
        av.matrix.write_parquet(frozen_dir / "availability.parquet")
        av.columns.write_parquet(frozen_dir / "availability_columns.parquet")
        software = reg.rows("SELECT name, environment, version, path FROM software ORDER BY name")
        resources = reg.rows("SELECT resource_id, name, version, genome_build, checksum, local_path FROM resources ORDER BY name")
        builds = sorted({a.genome_build for a in artifacts if a.genome_build})
        predicted = [a for a in artifacts if a.origin == Origin.PREDICTED]
        content: dict[str, Any] = {
            "snapshot_id": snapshot_id,
            "name": name,
            "created_at": utcnow().isoformat(),
            "efgpp_version": __version__,
            "project": project.config.project.name,
            "participant_registry_version": reg.table_version(["participants", "sample_aliases"]),
            "participants": av.n_participants,
            "genome_builds": builds,
            "artifacts": [_artifact_entry(a) for a in artifacts],
            "software": software,
            "resources": resources,
            "qc_configuration": project.config.genotype.model_dump(mode="json"),
            "predicted_modalities": [
                {"artifact_id": a.artifact_id, "modality": str(a.modality), "tissue": a.tissue,
                 "resource_versions": a.resource_versions} for a in predicted
            ],
            "timeline": timeline_summary(reg),
            "availability": {
                "columns": av.columns.to_dicts(),
                "matrix_sha256": str(checksum_paths([frozen_dir / "availability.parquet"]).value),
            },
            "configuration_checksums": project.config_checksums(),
            "lock_file_sha256": text_sha256(project.path(LOCK_FILE).read_text(encoding="utf-8"))
            if project.path(LOCK_FILE).exists() else None,
            "frozen_tables": {t: f"{name}/{t}.parquet" for t in FROZEN_TABLES},
        }
        text = yaml.safe_dump(content, sort_keys=False, allow_unicode=True, default_flow_style=False)
        reg.execute(
            "INSERT INTO snapshots (snapshot_id, name, created_at, path, checksum, content) VALUES (?, ?, ?, ?, ?, ?)",
            [snapshot_id, name, utcnow(), project.relative(path), text_sha256(text),
             json.dumps({"artifacts": [a.artifact_id for a in artifacts]})],
        )
        try:
            path.write_text(text, encoding="utf-8")
            _make_read_only(path, frozen_dir)
        except Exception:
            # All-or-nothing: never leave a half-written snapshot behind.
            reg.execute("DELETE FROM snapshots WHERE snapshot_id = ?", [snapshot_id])
            delete_snapshot_dir(frozen_dir)
            if path.exists():
                os.chmod(path, stat.S_IWRITE | stat.S_IREAD)
                path.unlink()
            raise
    return path


def _artifact_entry(a: Artifact) -> dict[str, Any]:
    return {
        "artifact_id": a.artifact_id, "artifact_name": a.artifact_name, "artifact_type": a.artifact_type,
        "modality": str(a.modality), "origin": a.origin.value, "status": a.status.value, "path": a.path,
        "format": a.format, "size": a.size, "checksum": a.checksum, "source_id": a.source_id,
        "participant_count": a.participant_count, "feature_count": a.feature_count,
        "genome_build": a.genome_build, "tissue": a.tissue,
        "temporal_type": a.temporal_type.value if a.temporal_type else None,
        "parent_artifact_ids": a.parent_artifact_ids, "tool": a.tool, "tool_version": a.tool_version,
        "resource_versions": a.resource_versions, "storage_mode": a.storage_mode,
        "files": [str(f) for f in _artifact_files(a)],
    }


def _make_read_only(*paths: Path) -> None:
    for p in paths:
        targets = [p] if p.is_file() else [f for f in p.rglob("*") if f.is_file()]
        for f in targets:
            os.chmod(f, stat.S_IREAD | stat.S_IRGRP | stat.S_IROTH)


def list_snapshots(project: Project) -> list[dict[str, Any]]:
    with Registry.open(project) as reg:
        return reg.rows("SELECT snapshot_id, name, created_at, path, checksum FROM snapshots ORDER BY created_at")


def protected_resource_paths(project: Project) -> set[str]:
    """Local paths of resources referenced by any frozen snapshot (never overwritten)."""
    out: set[str] = set()
    for f in project.path("snapshots").glob("*.yaml"):
        data = yaml.safe_load(f.read_text(encoding="utf-8")) or {}
        out.update(r["local_path"] for r in data.get("resources", []) if r.get("local_path"))
    return out


@dataclass
class Selection:
    participants: list[str]
    artifacts: dict[str, list[dict[str, Any]]] = field(default_factory=dict)
    columns: list[str] = field(default_factory=list)


@dataclass
class DataSnapshot:
    """Read-only view of a frozen snapshot, used by the Representation layer."""

    project: Project
    name: str
    content: dict[str, Any]

    @classmethod
    def load(cls, project: Project, name: str) -> DataSnapshot:
        path = project.path("snapshots", f"{name}.yaml")
        if not path.exists():
            raise FileNotFoundError(f"no snapshot named {name!r}")
        return cls(project, name, yaml.safe_load(path.read_text(encoding="utf-8")))

    @property
    def snapshot_id(self) -> str:
        return str(self.content["snapshot_id"])

    @property
    def created_at(self) -> datetime:
        return datetime.fromisoformat(self.content["created_at"])

    def table(self, name: str) -> pl.DataFrame:
        return pl.read_parquet(self.project.path("snapshots", self.name, f"{name}.parquet"))

    def availability(self) -> pl.DataFrame:
        return self.table("availability")

    def artifacts(self, *, modality: str | Modality | None = None, origin: str | Origin | None = None,
                  source_id: str | None = None, tissue: str | None = None,
                  artifact_type: str | None = None) -> list[dict[str, Any]]:
        out = []
        for a in self.content["artifacts"]:
            if modality is not None and a["modality"] != str(modality):
                continue
            if origin is not None and a["origin"] != str(origin):
                continue
            if source_id is not None and a["source_id"] != source_id:
                continue
            if tissue is not None and a["tissue"] != tissue:
                continue
            if artifact_type is not None and a["artifact_type"] != artifact_type:
                continue
            out.append(a)
        return out

    def artifact(self, artifact_id: str) -> dict[str, Any]:
        for a in self.content["artifacts"]:
            if a["artifact_id"] == artifact_id:
                return a
        raise KeyError(artifact_id)

    def columns(self) -> list[dict[str, Any]]:
        return list(self.content["availability"]["columns"])

    def phenotype_column(self, phenotype: str) -> str:
        """Availability column of a phenotype given its id or name."""
        defs = self.table("phenotype_definitions")
        hit = defs.filter((pl.col("phenotype_id") == phenotype) | (pl.col("name") == phenotype))
        if hit.height == 0:
            raise KeyError(f"phenotype {phenotype!r} is not in snapshot {self.name}")
        return str(hit.get_column("phenotype_id").item())

    def select(self, phenotype: str | None = None, *, require: list[str] | None = None,
               include: list[str] | None = None) -> Selection:
        """Participants having `phenotype` and every `require`d availability column, plus
        the exact artifacts for the requested modalities.

        `require` / `include` accept availability column names (e.g. "GENO001_QC",
        "RNA001_OBSERVED") or modality names (e.g. "genotype_qc", "expression").
        """
        av = self.availability()
        cols_meta = self.columns()

        def expand(key: str) -> list[str]:
            direct = [c["column"] for c in cols_meta if c["column"] == key]
            return direct or [c["column"] for c in cols_meta if c["modality"] == key]

        needed: list[str] = []
        if phenotype:
            needed.append(self.phenotype_column(phenotype))
        for key in require or []:
            cols = expand(key)
            if not cols:
                raise KeyError(f"{key!r} is not available in snapshot {self.name}")
            needed += cols
        mask = pl.all_horizontal([pl.col(c) for c in needed]) if needed else pl.lit(True)
        participants = av.filter(mask).get_column("participant_id").to_list()
        artifacts: dict[str, list[dict[str, Any]]] = {}
        for key in dict.fromkeys([*(require or []), *(include or [])]):
            cols = expand(key)
            ids = {c["artifact_id"] for c in cols_meta if c["column"] in cols}
            by_modality = self.artifacts(modality=key)
            artifacts[key] = [a for a in self.content["artifacts"] if a["artifact_id"] in ids] or by_modality
        if phenotype:
            pid = self.phenotype_column(phenotype)
            artifacts["phenotype"] = self.artifacts(source_id=pid, modality=Modality.PHENOTYPE, artifact_type="standardized")
        return Selection(participants, artifacts, needed)

    def verify(self) -> dict[str, str]:
        limit = int(self.project.config.storage.full_checksum_limit_gb * 1024**3)
        problems = {}
        for entry in self.content["artifacts"]:
            art = Artifact.model_validate({**entry, "metadata": {"members": entry.get("files")}})
            issue = verify_artifact(art, limit)
            if issue:
                problems[entry["artifact_id"]] = issue
        return problems


def delete_snapshot_dir(path: Path) -> None:
    if not path.exists():
        return
    for f in path.rglob("*"):
        os.chmod(f, stat.S_IWRITE | stat.S_IREAD)
    shutil.rmtree(path)
