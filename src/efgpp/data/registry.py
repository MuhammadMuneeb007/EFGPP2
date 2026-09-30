"""DuckDB-backed registry of participants, artifacts, provenance and snapshots.

The registry is the single source of truth that the Representation layer queries. Every
write happens under a file lock because DuckDB permits a single writing process.
Selected tables are mirrored to Parquet under registry/ for tool-agnostic access.
"""

from __future__ import annotations

import json
import threading
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import duckdb
import polars as pl
from filelock import FileLock

if TYPE_CHECKING:
    from efgpp.project import Project

SCHEMA_VERSION = 1

_ACTIVE = threading.local()  # per-thread open registries, for re-entrant opens

DDL = """
CREATE TABLE IF NOT EXISTS meta (key VARCHAR PRIMARY KEY, value VARCHAR);
CREATE TABLE IF NOT EXISTS counters (name VARCHAR PRIMARY KEY, value BIGINT);

CREATE TABLE IF NOT EXISTS participants (
    participant_id VARCHAR PRIMARY KEY,
    first_source VARCHAR,
    created_at TIMESTAMP,
    attributes VARCHAR
);
CREATE TABLE IF NOT EXISTS sample_aliases (
    source_id VARCHAR,
    native_id VARCHAR,
    participant_id VARCHAR,
    PRIMARY KEY (source_id, native_id)
);
CREATE TABLE IF NOT EXISTS events (
    participant_id VARCHAR,
    event_id VARCHAR,
    visit_name VARCHAR,
    event_date VARCHAR,
    study_day DOUBLE,
    age_at_event DOUBLE,
    source_id VARCHAR,
    PRIMARY KEY (participant_id, event_id)
);
CREATE TABLE IF NOT EXISTS biospecimens (
    biospecimen_id VARCHAR PRIMARY KEY,
    participant_id VARCHAR,
    event_id VARCHAR,
    tissue VARCHAR,
    material VARCHAR,
    collection_date VARCHAR,
    processing_method VARCHAR,
    storage_condition VARCHAR,
    source_id VARCHAR
);
CREATE TABLE IF NOT EXISTS assays (
    assay_id VARCHAR PRIMARY KEY,
    artifact_id VARCHAR,
    source_id VARCHAR,
    modality VARCHAR,
    origin VARCHAR,
    participant_id VARCHAR,
    event_id VARCHAR,
    biospecimen_id VARCHAR,
    tissue VARCHAR,
    collection_time VARCHAR,
    sample_key VARCHAR
);
CREATE TABLE IF NOT EXISTS phenotype_definitions (
    phenotype_id VARCHAR PRIMARY KEY,
    name VARCHAR,
    type VARCHAR,
    source_id VARCHAR,
    value_column VARCHAR,
    units VARCHAR,
    levels VARCHAR,
    ontology_term VARCHAR,
    definition VARCHAR,
    artifact_id VARCHAR
);
CREATE TABLE IF NOT EXISTS phenotype_observations (
    phenotype_id VARCHAR,
    participant_id VARCHAR,
    event_id VARCHAR,
    value_numeric DOUBLE,
    value_text VARCHAR,
    is_missing BOOLEAN
);
CREATE TABLE IF NOT EXISTS artifacts (
    artifact_id VARCHAR PRIMARY KEY,
    artifact_name VARCHAR,
    artifact_type VARCHAR,
    modality VARCHAR,
    origin VARCHAR,
    status VARCHAR,
    path VARCHAR,
    format VARCHAR,
    size BIGINT,
    checksum VARCHAR,
    participant_count BIGINT,
    feature_count BIGINT,
    genome_build VARCHAR,
    tissue VARCHAR,
    biospecimen VARCHAR,
    timepoint VARCHAR,
    temporal_type VARCHAR,
    source_id VARCHAR,
    storage_mode VARCHAR,
    tool VARCHAR,
    tool_version VARCHAR,
    resource_versions VARCHAR,
    configuration_hash VARCHAR,
    created_at TIMESTAMP,
    updated_at TIMESTAMP,
    command VARCHAR,
    environment VARCHAR,
    metadata VARCHAR,
    superseded_by VARCHAR
);
CREATE TABLE IF NOT EXISTS artifact_parents (
    artifact_id VARCHAR,
    parent_artifact_id VARCHAR,
    PRIMARY KEY (artifact_id, parent_artifact_id)
);
CREATE TABLE IF NOT EXISTS resources (
    resource_id VARCHAR PRIMARY KEY,
    name VARCHAR,
    version VARCHAR,
    release_date VARCHAR,
    genome_build VARCHAR,
    source VARCHAR,
    license VARCHAR,
    checksum VARCHAR,
    download_date VARCHAR,
    local_path VARCHAR,
    metadata VARCHAR
);
CREATE TABLE IF NOT EXISTS software (
    name VARCHAR,
    environment VARCHAR,
    version VARCHAR,
    path VARCHAR,
    detected_at TIMESTAMP,
    PRIMARY KEY (name, environment)
);
CREATE TABLE IF NOT EXISTS tool_runs (
    run_id VARCHAR PRIMARY KEY,
    step_id VARCHAR,
    tool VARCHAR,
    tool_version VARCHAR,
    command VARCHAR,
    inputs VARCHAR,
    outputs VARCHAR,
    started_at TIMESTAMP,
    completed_at TIMESTAMP,
    exit_code INTEGER,
    status VARCHAR,
    environment VARCHAR,
    resources VARCHAR,
    log_path VARCHAR
);
CREATE TABLE IF NOT EXISTS qc_runs (
    qc_run_id VARCHAR PRIMARY KEY,
    artifact_id VARCHAR,
    modality VARCHAR,
    run_id VARCHAR,
    status VARCHAR,
    created_at TIMESTAMP,
    config VARCHAR,
    summary VARCHAR
);
CREATE TABLE IF NOT EXISTS qc_metrics (
    qc_run_id VARCHAR,
    artifact_id VARCHAR,
    scope VARCHAR,
    metric VARCHAR,
    value DOUBLE,
    text_value VARCHAR,
    threshold VARCHAR,
    status VARCHAR
);
CREATE TABLE IF NOT EXISTS modality_availability (
    participant_id VARCHAR,
    column_name VARCHAR,
    origin VARCHAR,
    modality VARCHAR,
    artifact_id VARCHAR,
    available BOOLEAN
);
CREATE TABLE IF NOT EXISTS molecular_models (
    model_id VARCHAR,
    provider VARCHAR,
    provider_dataset_id VARCHAR,
    modality VARCHAR,
    feature_id VARCHAR,
    feature_name VARCHAR,
    tissue VARCHAR,
    platform VARCHAR,
    training_cohort VARCHAR,
    training_ancestry VARCHAR,
    genome_build VARCHAR,
    validation_r2 DOUBLE,
    n_variants BIGINT,
    coverage_fraction DOUBLE,
    status VARCHAR,
    resource_version VARCHAR,
    resource_sha256 VARCHAR,
    local_path VARCHAR
);
CREATE TABLE IF NOT EXISTS snapshots (
    snapshot_id VARCHAR PRIMARY KEY,
    name VARCHAR UNIQUE,
    created_at TIMESTAMP,
    path VARCHAR,
    checksum VARCHAR,
    content VARCHAR
);
"""

# registry table -> exported parquet file name (section 10 of the specification)
PARQUET_EXPORTS = {
    "participants": "participants.parquet",
    "sample_aliases": "aliases.parquet",
    "events": "events.parquet",
    "biospecimens": "biospecimens.parquet",
    "assays": "assays.parquet",
    "phenotype_definitions": "phenotypes.parquet",
    "molecular_models": "molecular_models.parquet",
}

JSON_COLUMNS = {
    "resource_versions",
    "environment",
    "metadata",
    "attributes",
    "levels",
    "ontology_term",
    "definition",
    "inputs",
    "outputs",
    "resources",
    "config",
    "summary",
    "content",
}


def utcnow() -> datetime:
    return datetime.now(UTC).replace(tzinfo=None)


def to_json(value: Any) -> str | None:
    if value is None:
        return None
    return json.dumps(value, sort_keys=True, default=str)


def from_json(value: Any) -> Any:
    if value is None or not isinstance(value, str):
        return value
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return value


class Registry:
    """Thin, explicit wrapper over a DuckDB connection."""

    def __init__(self, con: duckdb.DuckDBPyConnection, project: Project) -> None:
        self.con = con
        self.project = project

    # --------------------------------------------------------------- lifecycle
    @classmethod
    @contextmanager
    def open(cls, project: Project, read_only: bool = False) -> Iterator[Registry]:
        """Open the registry. Nested opens in the same thread reuse the open connection,
        so code that records provenance never deadlocks against its caller's lock."""
        db_path: Path = project.registry_path
        key = str(db_path.resolve())
        active: dict[str, list[Any]] = _ACTIVE.__dict__.setdefault("open", {})
        if key in active:
            active[key][1] += 1
            try:
                yield active[key][0]
            finally:
                active[key][1] -= 1
            return
        db_path.parent.mkdir(parents=True, exist_ok=True)
        project.path(".efgpp", "locks").mkdir(parents=True, exist_ok=True)
        lock = FileLock(str(project.path(".efgpp", "locks", "registry.lock")), timeout=3600)
        with lock:
            con = duckdb.connect(str(db_path))
            try:
                reg = cls(con, project)
                reg._migrate()
                active[key] = [reg, 1]
                yield reg  # DuckDB auto-commits each statement
            finally:
                active.pop(key, None)
                con.close()

    def _migrate(self) -> None:
        self.con.execute(DDL)
        current = self.get_meta("schema_version")
        if current is None:
            self.set_meta("schema_version", str(SCHEMA_VERSION))

    # -------------------------------------------------------------- primitives
    def execute(self, sql: str, params: Sequence[Any] | None = None) -> duckdb.DuckDBPyConnection:
        return self.con.execute(sql, params or [])

    def rows(self, sql: str, params: Sequence[Any] | None = None) -> list[dict[str, Any]]:
        cur = self.con.execute(sql, params or [])
        cols = [d[0] for d in cur.description]
        out = []
        for row in cur.fetchall():
            rec = dict(zip(cols, row, strict=True))
            for k in JSON_COLUMNS & rec.keys():
                rec[k] = from_json(rec[k])
            out.append(rec)
        return out

    def one(self, sql: str, params: Sequence[Any] | None = None) -> dict[str, Any] | None:
        r = self.rows(sql, params)
        return r[0] if r else None

    def scalar(self, sql: str, params: Sequence[Any] | None = None) -> Any:
        row = self.con.execute(sql, params or []).fetchone()
        return row[0] if row else None

    def frame(self, sql: str, params: Sequence[Any] | None = None) -> pl.DataFrame:
        return self.con.execute(sql, params or []).pl()

    def insert_frame(self, table: str, df: pl.DataFrame, replace: bool = True) -> int:
        """Bulk insert a polars frame whose columns are a subset of `table`'s columns."""
        if df.height == 0:
            return 0
        self.con.register("_efgpp_tmp", df.to_arrow())
        try:
            cols = ", ".join(f'"{c}"' for c in df.columns)
            verb = "INSERT OR REPLACE" if replace else "INSERT"
            self.con.execute(f"{verb} INTO {table} ({cols}) SELECT {cols} FROM _efgpp_tmp")
        finally:
            self.con.unregister("_efgpp_tmp")
        return df.height

    def upsert(self, table: str, record: dict[str, Any]) -> None:
        rec = {k: (to_json(v) if k in JSON_COLUMNS else v) for k, v in record.items()}
        cols = ", ".join(f'"{c}"' for c in rec)
        marks = ", ".join("?" for _ in rec)
        self.con.execute(f"INSERT OR REPLACE INTO {table} ({cols}) VALUES ({marks})", list(rec.values()))

    # ------------------------------------------------------------- meta / ids
    def get_meta(self, key: str) -> str | None:
        return self.scalar("SELECT value FROM meta WHERE key = ?", [key])

    def set_meta(self, key: str, value: str) -> None:
        self.con.execute("INSERT OR REPLACE INTO meta VALUES (?, ?)", [key, value])

    def next_id(self, prefix: str, width: int = 6) -> str:
        current = self.scalar("SELECT value FROM counters WHERE name = ?", [prefix]) or 0
        nxt = int(current) + 1
        self.con.execute("INSERT OR REPLACE INTO counters VALUES (?, ?)", [prefix, nxt])
        return f"{prefix}{nxt:0{width}d}"

    # ----------------------------------------------------------------- exports
    def export_parquet(self) -> list[Path]:
        written = []
        compression = self.project.config.storage.compression.upper()
        for table, filename in PARQUET_EXPORTS.items():
            out = self.project.registry_path.parent / filename
            self.con.execute(
                f"COPY (SELECT * FROM {table}) TO '{out.as_posix()}' "
                f"(FORMAT PARQUET, COMPRESSION {compression})"
            )
            written.append(out)
        return written

    def table_version(self, tables: Sequence[str]) -> str:
        """Content hash of registry tables; used as the participant-registry version."""
        import hashlib

        h = hashlib.sha256()
        for t in tables:
            df = self.frame(f"SELECT * FROM {t}")
            if df.height:
                df = df.sort(df.columns)
            h.update(t.encode())
            h.update(df.write_csv().encode())
        return h.hexdigest()[:16]
