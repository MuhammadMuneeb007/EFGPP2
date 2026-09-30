"""Artifacts: every registered or generated file set, with full lineage.

An artifact is never deleted because a newer one exists; it is marked SUPERSEDED and
keeps pointing at its successor.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

from pydantic import BaseModel, Field

from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.registry import Registry, to_json, utcnow


class Artifact(BaseModel):
    artifact_id: str | None = None
    artifact_name: str
    artifact_type: str  # e.g. source, standardized, qc_table, qc_genotype, pca, annotation
    modality: Modality | str
    origin: Origin
    status: ArtifactStatus = ArtifactStatus.REGISTERED
    path: str
    format: str
    size: int | None = None
    checksum: str | None = None
    participant_count: int | None = None
    feature_count: int | None = None
    genome_build: str | None = None
    tissue: str | None = None
    biospecimen: str | None = None
    timepoint: str | None = None
    temporal_type: TemporalType | None = None
    source_id: str | None = None
    storage_mode: str | None = None
    parent_artifact_ids: list[str] = Field(default_factory=list)
    tool: str | None = None
    tool_version: str | None = None
    resource_versions: dict[str, str] = Field(default_factory=dict)
    configuration_hash: str | None = None
    created_at: datetime | None = None
    updated_at: datetime | None = None
    command: str | None = None
    environment: dict[str, Any] = Field(default_factory=dict)
    metadata: dict[str, Any] = Field(default_factory=dict)
    superseded_by: str | None = None

    def to_record(self) -> dict[str, Any]:
        rec = self.model_dump(mode="json", exclude={"parent_artifact_ids"})
        rec["created_at"] = self.created_at
        rec["updated_at"] = self.updated_at
        return rec


_COLUMNS = [f for f in Artifact.model_fields if f != "parent_artifact_ids"]


class ArtifactStore:
    def __init__(self, registry: Registry) -> None:
        self.reg = registry

    def register(self, artifact: Artifact) -> Artifact:
        now = utcnow()
        if artifact.artifact_id is None:
            artifact.artifact_id = self.reg.next_id("ART")
        artifact.created_at = artifact.created_at or now
        artifact.updated_at = now
        self.reg.upsert("artifacts", artifact.to_record())
        self.reg.execute("DELETE FROM artifact_parents WHERE artifact_id = ?", [artifact.artifact_id])
        for parent in dict.fromkeys(artifact.parent_artifact_ids):
            self.reg.execute(
                "INSERT OR REPLACE INTO artifact_parents VALUES (?, ?)",
                [artifact.artifact_id, parent],
            )
        return artifact

    def register_replacing(self, artifact: Artifact) -> Artifact:
        """Register `artifact` and supersede live artifacts with the same source/type/name."""
        olds = [
            a for a in self.find(source_id=artifact.source_id, artifact_type=artifact.artifact_type)
            if a.artifact_name == artifact.artifact_name
        ]
        self.register(artifact)
        for old in olds:
            if old.artifact_id != artifact.artifact_id:
                self.supersede(old.artifact_id, artifact.artifact_id)  # type: ignore[arg-type]
        return artifact

    def get(self, artifact_id: str) -> Artifact:
        row = self.reg.one("SELECT * FROM artifacts WHERE artifact_id = ?", [artifact_id])
        if row is None:
            raise KeyError(f"unknown artifact {artifact_id!r}")
        return self._hydrate(row)

    def find(
        self,
        *,
        source_id: str | None = None,
        modality: str | None = None,
        origin: str | None = None,
        artifact_type: str | None = None,
        statuses: set[ArtifactStatus] | None = None,
        include_superseded: bool = False,
        tissue: str | None = None,
    ) -> list[Artifact]:
        where, params = [], []
        for col, val in (
            ("source_id", source_id),
            ("modality", modality),
            ("origin", origin),
            ("artifact_type", artifact_type),
            ("tissue", tissue),
        ):
            if val is not None:
                where.append(f"{col} = ?")
                params.append(str(val))
        if not include_superseded:
            where.append("status <> ?")
            params.append(ArtifactStatus.SUPERSEDED.value)
        if statuses:
            where.append(f"status IN ({', '.join('?' for _ in statuses)})")
            params.extend(s.value for s in statuses)
        sql = "SELECT * FROM artifacts"
        if where:
            sql += " WHERE " + " AND ".join(where)
        sql += " ORDER BY artifact_id"
        return [self._hydrate(r) for r in self.reg.rows(sql, params)]

    def latest(self, **filters: Any) -> Artifact | None:
        found = self.find(**filters)
        return found[-1] if found else None

    def set_status(self, artifact_id: str, status: ArtifactStatus, **metadata: Any) -> None:
        art = self.get(artifact_id)
        art.status = status
        if metadata:
            art.metadata.update(metadata)
        art.updated_at = utcnow()
        self.reg.execute(
            "UPDATE artifacts SET status = ?, metadata = ?, updated_at = ? WHERE artifact_id = ?",
            [status.value, to_json(art.metadata), art.updated_at, artifact_id],
        )

    def update(self, artifact: Artifact) -> None:
        artifact.updated_at = utcnow()
        self.reg.upsert("artifacts", artifact.to_record())

    def supersede(self, old_id: str, new_id: str) -> None:
        self.reg.execute(
            "UPDATE artifacts SET status = ?, superseded_by = ?, updated_at = ? WHERE artifact_id = ?",
            [ArtifactStatus.SUPERSEDED.value, new_id, utcnow(), old_id],
        )

    def parents(self, artifact_id: str) -> list[str]:
        rows = self.reg.rows(
            "SELECT parent_artifact_id FROM artifact_parents WHERE artifact_id = ? ORDER BY 1",
            [artifact_id],
        )
        return [r["parent_artifact_id"] for r in rows]

    def lineage(self, artifact_id: str) -> list[str]:
        """All ancestors, nearest first."""
        seen: list[str] = []
        frontier = [artifact_id]
        while frontier:
            nxt = []
            for a in frontier:
                for p in self.parents(a):
                    if p not in seen:
                        seen.append(p)
                        nxt.append(p)
            frontier = nxt
        return seen

    def _hydrate(self, row: dict[str, Any]) -> Artifact:
        data = {k: row.get(k) for k in _COLUMNS}
        for key in ("resource_versions", "environment", "metadata"):
            data[key] = data.get(key) or {}
        data["parent_artifact_ids"] = self.parents(row["artifact_id"])
        return Artifact.model_validate(data)
