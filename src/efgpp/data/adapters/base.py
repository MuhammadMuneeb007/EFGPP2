"""The DataAdapter interface every modality implements.

Lifecycle per source:

    register -> inspect -> validate -> standardize -> qc -> derive -> summarize/report

Adapters only perform phenotype-independent work. Nothing here may look at a target
phenotype to select features, fit scalers or impute values.
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar

import polars as pl

from efgpp.config.data import SourceBase
from efgpp.constants import ArtifactStatus, Modality, TemporalType
from efgpp.data.aliases import AliasResolver, Resolution
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.io import detect_table_format, read_table
from efgpp.data.registry import Registry, to_json, utcnow
from efgpp.data.storage import place_source
from efgpp.data.validation import ValidationReport
from efgpp.project import Project


@dataclass
class AdapterContext:
    project: Project
    registry: Registry
    resolver: AliasResolver = field(init=False)

    def __post_init__(self) -> None:
        self.resolver = AliasResolver(self.project)

    @property
    def artifacts(self) -> ArtifactStore:
        return ArtifactStore(self.registry)


@dataclass
class QCResult:
    status: ArtifactStatus
    metrics: dict[str, Any] = field(default_factory=dict)
    messages: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)


class DataAdapter(ABC):
    """Base class. Subclasses set `modality` and implement the abstract methods."""

    modality: ClassVar[Modality]
    temporal_type: ClassVar[TemporalType] = TemporalType.UNKNOWN

    def __init__(self, ctx: AdapterContext, source: SourceBase) -> None:
        self.ctx = ctx
        self.source = source

    # ------------------------------------------------------------------ helpers
    @property
    def project(self) -> Project:
        return self.ctx.project

    @property
    def reg(self) -> Registry:
        return self.ctx.registry

    @property
    def store(self) -> ArtifactStore:
        return self.ctx.artifacts

    @property
    def source_path(self) -> Path:
        return self.project.resolve(self.source.path)

    def subject(self) -> str:
        return f"{self.modality.value.upper()} {self.source.id}"

    def source_artifact(self) -> Artifact | None:
        return self.store.latest(source_id=self.source.id, artifact_type="source")

    def latest(self, artifact_type: str, **filters: Any) -> Artifact | None:
        return self.store.latest(source_id=self.source.id, artifact_type=artifact_type, **filters)

    def register_replacing(self, artifact: Artifact) -> Artifact:
        """Register `artifact`, superseding any live artifact of the same kind."""
        previous = self.store.find(
            source_id=artifact.source_id,
            artifact_type=artifact.artifact_type,
            modality=str(artifact.modality),
        )
        previous = [p for p in previous if p.artifact_name == artifact.artifact_name]
        self.store.register(artifact)
        for p in previous:
            if p.artifact_id != artifact.artifact_id:
                self.store.supersede(p.artifact_id, artifact.artifact_id)  # type: ignore[arg-type]
        return artifact

    def configuration_hash(self) -> str:
        return self.source.checksum()

    # ---------------------------------------------------------------- interface
    @classmethod
    def detect(cls, path: Path) -> bool:
        """Whether this adapter can read `path`."""
        return path.exists()

    def members(self) -> tuple[Path, list[Path], str]:
        """(record path, member files, format) of the raw source."""
        path = self.source_path
        return path, [path], detect_table_format(path, self.source.format)

    def register(self) -> Artifact:
        record_path, members, fmt = self.members()
        missing = [m for m in members if not m.exists()]
        if missing:
            raise FileNotFoundError(f"{self.source.id}: missing {', '.join(map(str, missing))}")
        placed = place_source(
            self.project, modality=self.modality, origin=self.source.origin,
            source_id=self.source.id, record_path=record_path, members=members,
            requested=self.source.mode, fmt=fmt,
        )
        current = self.source_artifact()
        if current and current.checksum == placed.checksum and current.path == str(placed.record_path):
            return current  # unchanged: registration is idempotent
        art = Artifact(
            artifact_name=f"{self.source.id}_source",
            artifact_type="source",
            modality=self.modality,
            origin=self.source.origin,
            path=str(placed.record_path),
            format=fmt,
            size=placed.size,
            checksum=placed.checksum,
            source_id=self.source.id,
            storage_mode=placed.mode.value,
            temporal_type=self.temporal_type,
            configuration_hash=self.configuration_hash(),
            tool="efgpp",
            metadata={
                "original_path": str(self.source_path),
                "members": [str(m) for m in placed.members],
                "source_mtime": placed.mtime,
                "notes": placed.notes,
                "source_config": self.source.to_yaml_dict(),
            },
        )
        return self.register_replacing(art)

    @abstractmethod
    def inspect(self) -> dict[str, Any]: ...

    @abstractmethod
    def validate(self) -> ValidationReport: ...

    @abstractmethod
    def standardize(self) -> Artifact | None: ...

    def qc(self) -> QCResult | None:
        return None

    def derive(self) -> list[Artifact]:
        return []

    def summarize(self) -> dict[str, Any]:
        art = self.latest("standardized") or self.source_artifact()
        return {
            "source_id": self.source.id,
            "modality": self.modality.value,
            "origin": self.source.origin.value,
            "artifact_id": art.artifact_id if art else None,
            "status": art.status.value if art else None,
            "participants": art.participant_count if art else None,
            "features": art.feature_count if art else None,
        }

    def report(self) -> dict[str, Any]:
        return self.summarize()

    # ------------------------------------------------------------ shared steps
    def resolve_ids(self, native_ids: pl.Series) -> Resolution:
        from efgpp.data.participants import register_mapping

        res = self.ctx.resolver.resolve(self.source.id, native_ids)
        register_mapping(self.reg, self.source.id, res.mapping)
        return res

    def record_assays(
        self,
        artifact: Artifact,
        participant_ids: pl.Series,
        *,
        event_ids: pl.Series | None = None,
        biospecimen_ids: pl.Series | None = None,
        times: pl.Series | None = None,
        sample_keys: pl.Series | None = None,
    ) -> int:
        """One assay row per participant (x event) measured by this artifact."""
        n = participant_ids.len()

        def col(s: pl.Series | None) -> pl.Series:
            return s.cast(pl.Utf8) if s is not None else pl.Series([None] * n, dtype=pl.Utf8)

        df = pl.DataFrame({
            "artifact_id": [artifact.artifact_id] * n,
            "source_id": [self.source.id] * n,
            "modality": [str(artifact.modality)] * n,
            "origin": [artifact.origin.value] * n,
            "participant_id": participant_ids.cast(pl.Utf8),
            "event_id": col(event_ids),
            "biospecimen_id": col(biospecimen_ids),
            "tissue": [artifact.tissue] * n,
            "collection_time": col(times),
            "sample_key": col(sample_keys),
        }, schema_overrides={"tissue": pl.Utf8, "artifact_id": pl.Utf8})
        df = df.filter(pl.col("participant_id").is_not_null())
        df = df.with_columns(
            (pl.lit(f"{artifact.artifact_id}:") + pl.int_range(pl.len()).cast(pl.Utf8)).alias("assay_id")
        )
        self.reg.execute("DELETE FROM assays WHERE source_id = ?", [self.source.id])
        return self.reg.insert_frame("assays", df)

    def save_qc(self, artifact: Artifact, result: QCResult, config: dict[str, Any] | None = None) -> str:
        qc_run_id = self.reg.next_id("QC")
        summary_path = self.project.qc_dir(self.modality_qc_dir()) / f"{self.source.id}_qc_summary.json"
        summary_path.write_text(
            json.dumps({"status": result.status.value, "metrics": result.metrics,
                        "messages": result.messages}, indent=2, default=str),
            encoding="utf-8",
        )
        result.outputs.append(self.project.relative(summary_path))
        self.reg.upsert("qc_runs", {
            "qc_run_id": qc_run_id, "artifact_id": artifact.artifact_id,
            "modality": self.modality.value, "run_id": None, "status": result.status.value,
            "created_at": utcnow(), "config": config or {}, "summary": result.metrics,
        })
        rows = []
        for key, value in _flatten(result.metrics).items():
            numeric = value if isinstance(value, int | float) and not isinstance(value, bool) else None
            rows.append({
                "qc_run_id": qc_run_id, "artifact_id": artifact.artifact_id, "scope": "dataset",
                "metric": key, "value": float(numeric) if numeric is not None else None,
                "text_value": None if numeric is not None else to_json(value), "threshold": None,
                "status": None,
            })
        if rows:
            self.reg.insert_frame("qc_metrics", pl.DataFrame(rows, schema_overrides={
                "value": pl.Float64, "text_value": pl.Utf8, "threshold": pl.Utf8, "status": pl.Utf8}),
                replace=False)
        self.store.set_status(artifact.artifact_id, result.status, qc_run_id=qc_run_id)  # type: ignore[arg-type]
        return qc_run_id

    def modality_qc_dir(self) -> str:
        return self.modality.value


def _flatten(d: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten(v, key + "."))
        else:
            out[key] = v
    return out


class TabularAdapter(DataAdapter):
    """Shared behaviour for table-shaped sources (phenotypes, covariates, clinical)."""

    def read(self) -> pl.DataFrame:
        art = self.source_artifact()
        path = Path(art.path) if art else self.source_path
        return read_table(path, self.source.format)

    def inspect(self) -> dict[str, Any]:
        path = self.source_path
        info: dict[str, Any] = {
            "source_id": self.source.id, "modality": self.modality.value,
            "origin": self.source.origin.value, "path": str(path), "exists": path.exists(),
        }
        if not path.exists():
            return info
        df = read_table(path, self.source.format)
        id_col = getattr(self.source, "participant_id_column", None)
        info.update({
            "format": detect_table_format(path, self.source.format),
            "rows": df.height,
            "columns": df.width,
            "participants": df.get_column(id_col).n_unique() if id_col in df.columns else None,
            "size": path.stat().st_size if path.is_file() else None,
        })
        return info
