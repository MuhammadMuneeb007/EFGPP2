"""Clinical measurements (laboratory values, vitals, ...), wide or long, optionally timed.

Standardized to a long table: participant_id, event_id, variable, value_numeric,
value_text, unit.
"""

from __future__ import annotations

import polars as pl

from efgpp.config.data import ClinicalSource
from efgpp.constants import ArtifactStatus, Modality, TemporalType
from efgpp.data.adapters.base import QCResult, TabularAdapter
from efgpp.data.artifacts import Artifact
from efgpp.data.io import as_text, missing_mask, write_parquet
from efgpp.data.timeline import events_from_spec, register_events
from efgpp.data.validation import ValidationReport, plural
from efgpp.resources.checksums import checksum_paths


class ClinicalAdapter(TabularAdapter):
    modality = Modality.CLINICAL
    temporal_type = TemporalType.EVENT
    source: ClinicalSource

    def _meta_columns(self) -> list[str]:
        s = self.source
        cols = [s.participant_id_column]
        if s.timeline:
            cols += [c for c in (s.timeline.event_column, s.timeline.time_column, s.timeline.study_day_column,
                                 s.timeline.visit_name_column, s.timeline.age_column) if c]
        return cols

    def _variables(self, df: pl.DataFrame) -> list[str]:
        s = self.source
        if s.variables:
            return s.variables
        return [c for c in df.columns if c not in self._meta_columns()]

    def validate(self) -> ValidationReport:
        s = self.source
        rep = ValidationReport(self.subject(), s.id)
        if not self.source_path.exists():
            rep.fail(f"file not found: {self.source_path}")
            return rep
        df = self.read()
        needed = [s.participant_id_column]
        if s.layout == "long":
            needed += [s.variable_column, s.value_column]  # type: ignore[list-item]
        if s.timeline and s.timeline.event_column:
            needed.append(s.timeline.event_column)
        absent = [c for c in needed if c not in df.columns]
        if absent:
            rep.fail(f"required columns not found: {absent}")
            return rep
        rep.ok("required columns present")
        if s.layout == "wide":
            missing_vars = [v for v in self._variables(df) if v not in df.columns]
            if missing_vars:
                rep.fail(f"{plural(len(missing_vars), 'variable')} not found", len(missing_vars), missing_vars)
            key_cols = [s.participant_id_column] + ([s.timeline.event_column] if s.timeline and s.timeline.event_column else [])
            n_dup = df.select(key_cols).is_duplicated().sum()
            if n_dup:
                rep.fail(f"{plural(int(n_dup), 'duplicated participant/event row')}", int(n_dup))
            else:
                rep.ok("one row per participant" + (" and event" if len(key_cols) > 1 else ""))
        return rep

    def _long(self, df: pl.DataFrame, pids: pl.Series) -> pl.DataFrame:
        s = self.source
        ev = s.timeline.event_column if s.timeline else None
        base = pl.DataFrame({
            "participant_id": pids,
            "event_id": as_text(df, ev) if ev else pl.Series([None] * df.height, dtype=pl.Utf8),
        })
        if s.layout == "long":
            text = as_text(df, s.value_column)  # type: ignore[arg-type]
            miss = missing_mask(text, s.missing_values)
            long = base.with_columns(
                as_text(df, s.variable_column).alias("variable"),  # type: ignore[arg-type]
                pl.Series([None if m else t for t, m in zip(text.to_list(), miss.to_list(), strict=True)], dtype=pl.Utf8).alias("value_text"),
                (as_text(df, s.unit_column) if s.unit_column else pl.Series([None] * df.height, dtype=pl.Utf8)).alias("unit"),
            )
        else:
            variables = self._variables(df)
            wide = base.with_columns([as_text(df, v).alias(v) for v in variables])
            long = wide.unpivot(index=["participant_id", "event_id"], on=variables,
                                variable_name="variable", value_name="value_text")
            long = long.with_columns(
                pl.when(pl.col("value_text").is_in(s.missing_values)).then(None).otherwise(pl.col("value_text")).alias("value_text"),
                pl.lit(None, dtype=pl.Utf8).alias("unit"),
            )
        return long.with_columns(pl.col("value_text").cast(pl.Float64, strict=False).alias("value_numeric")).select(
            "participant_id", "event_id", "variable", "value_numeric", "value_text", "unit"
        ).filter(pl.col("participant_id").is_not_null())

    def standardize(self) -> Artifact:
        s = self.source
        df = self.read()
        ids = as_text(df, s.participant_id_column)
        res = self.resolve_ids(ids)
        pids = ids.replace_strict(dict(res.mapping.iter_rows()), default=None)
        long = self._long(df, pids)
        events = events_from_spec(df, pids, s.timeline, s.id)
        if events is not None:
            register_events(self.reg, events.filter(pl.col("participant_id").is_not_null()))
        out = write_parquet(long, self.project.artifact_dir(s.origin, Modality.CLINICAL, s.id) / "clinical_long.parquet",
                            self.project.config.storage.compression)
        src = self.source_artifact()
        observed = long.filter(pl.col("value_text").is_not_null())
        art = self.register_replacing(Artifact(
            artifact_name=f"{s.id}_clinical", artifact_type="standardized", modality=Modality.CLINICAL,
            origin=s.origin, status=ArtifactStatus.STANDARDIZED, path=str(out), format="parquet",
            size=out.stat().st_size, checksum=str(checksum_paths([out])),
            participant_count=observed.get_column("participant_id").n_unique(),
            feature_count=long.get_column("variable").n_unique(), source_id=s.id,
            temporal_type=TemporalType.EVENT if s.timeline and s.timeline.event_column else TemporalType.STATIC,
            parent_artifact_ids=[src.artifact_id] if src and src.artifact_id else [], tool="efgpp",
            configuration_hash=self.configuration_hash(),
        ))
        keys = observed.select("participant_id", "event_id").unique(maintain_order=True)
        self.record_assays(art, keys.get_column("participant_id"), event_ids=keys.get_column("event_id"))
        return art

    def qc(self) -> QCResult | None:
        art = self.latest("standardized")
        if art is None:
            return None
        long = pl.read_parquet(art.path)
        per_var = (
            long.group_by("variable")
            .agg(
                pl.len().alias("n"),
                pl.col("value_text").null_count().alias("n_missing"),
                pl.col("value_numeric").mean().alias("mean"),
                pl.col("value_numeric").min().alias("min"),
                pl.col("value_numeric").max().alias("max"),
            )
            .sort("variable")
        )
        metrics = {"variables": {r["variable"]: {k: v for k, v in r.items() if k != "variable"} for r in per_var.to_dicts()}}
        result = QCResult(ArtifactStatus.QC_PASS, metrics)
        self.save_qc(art, result)
        self.store.set_status(art.artifact_id, ArtifactStatus.READY, qc_status="QC_PASS")  # type: ignore[arg-type]
        return result
