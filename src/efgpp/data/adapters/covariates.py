"""Covariates: typed and profiled, never imputed or scaled (that would leak information)."""

from __future__ import annotations

from typing import Any

import polars as pl

from efgpp.config.data import CovariateSource
from efgpp.constants import ArtifactStatus, Modality, TemporalType
from efgpp.data.adapters.base import QCResult, TabularAdapter
from efgpp.data.artifacts import Artifact
from efgpp.data.io import as_text, missing_mask, write_parquet
from efgpp.data.schemas.covariates import covariate_schema
from efgpp.data.timeline import events_from_spec, register_events
from efgpp.data.validation import ValidationReport, plural
from efgpp.resources.checksums import checksum_paths

MISSING_WARN = 0.05


def type_columns(
    df: pl.DataFrame, variables: list[str], categorical: list[str], missing_tokens: list[str]
) -> tuple[pl.DataFrame, dict[str, str]]:
    """Cast each variable to Float64 (numeric) or Utf8 (categorical). Returns (frame, kinds)."""
    kinds: dict[str, str] = {}
    cols = []
    for v in variables:
        text = as_text(df, v)
        text = pl.Series(v, [None if m else x for x, m in zip(text.to_list(), missing_mask(text, missing_tokens).to_list(), strict=True)], dtype=pl.Utf8)
        if v in categorical:
            kinds[v] = "categorical"
            cols.append(text)
            continue
        num = text.cast(pl.Float64, strict=False)
        if int((text.is_not_null() & num.is_null()).sum()) == 0:
            kinds[v] = "numeric"
            cols.append(num)
        else:
            kinds[v] = "categorical"
            cols.append(text)
    return pl.DataFrame(cols), kinds


class CovariateAdapter(TabularAdapter):
    modality = Modality.COVARIATES
    temporal_type = TemporalType.STATIC
    source: CovariateSource

    def _event_column(self) -> str | None:
        return self.source.timeline.event_column if self.source.timeline else None

    def validate(self) -> ValidationReport:
        s = self.source
        rep = ValidationReport(self.subject(), s.id)
        if not self.source_path.exists():
            rep.fail(f"file not found: {self.source_path}")
            return rep
        try:
            df = self.read()
        except Exception as exc:  # noqa: BLE001
            rep.fail(f"cannot read table: {exc}")
            return rep
        if s.participant_id_column not in df.columns:
            rep.fail(f"participant ID column {s.participant_id_column!r} not found")
            return rep
        rep.ok(f"participant ID column exists ({s.participant_id_column})")
        absent = [v for v in s.variables if v not in df.columns]
        if absent:
            rep.fail(f"{plural(len(absent), 'covariate column')} not found", len(absent), absent)
        present = [v for v in s.variables if v in df.columns]
        if present:
            rep.ok(f"{plural(len(present), 'covariate column')} found")
        ev = self._event_column()
        ids = as_text(df, s.participant_id_column)
        keys = pl.DataFrame({"id": ids, **({"ev": as_text(df, ev)} if ev and ev in df.columns else {})})
        dup = keys.filter(keys.is_duplicated())
        if dup.height:
            n = dup.get_column("id").n_unique()
            rep.fail(f"{plural(n, 'duplicated participant ID')}", n, dup.get_column("id").unique().to_list())
        else:
            rep.ok("participant IDs are unique")
        if present:
            typed, kinds = type_columns(df, present, s.categorical, s.missing_values)
            for v in present:
                frac = typed.get_column(v).null_count() / max(df.height, 1)
                if frac > MISSING_WARN:
                    rep.warn(f"{v}: {frac:.1%} missing")
                if kinds[v] == "categorical" and v not in s.categorical:
                    rep.warn(f"{v}: non-numeric values; treated as categorical (list it under `categorical` to silence)")
        return rep

    def standardize(self) -> Artifact:
        s = self.source
        df = self.read()
        ids = as_text(df, s.participant_id_column)
        res = self.resolve_ids(ids)
        pids = ids.replace_strict(dict(res.mapping.iter_rows()), default=None)
        typed, kinds = type_columns(df, s.variables, s.categorical, s.missing_values)
        ev = self._event_column()
        base = {"participant_id": pids}
        if ev:
            base["event_id"] = as_text(df, ev)
        out_df = pl.DataFrame(base).hstack(typed).filter(
            pl.col("participant_id").is_not_null()
        )
        out_df = covariate_schema(bool(ev)).validate(out_df)
        events = events_from_spec(df, pids, s.timeline, s.id)
        if events is not None:
            register_events(self.reg, events.filter(pl.col("participant_id").is_not_null()))

        out = write_parquet(
            out_df,
            self.project.artifact_dir(s.origin, Modality.COVARIATES, s.id) / "covariates.parquet",
            self.project.config.storage.compression,
        )
        src = self.source_artifact()
        art = self.register_replacing(Artifact(
            artifact_name=f"{s.id}_covariates", artifact_type="standardized",
            modality=Modality.COVARIATES, origin=s.origin, status=ArtifactStatus.STANDARDIZED,
            path=str(out), format="parquet", size=out.stat().st_size,
            checksum=str(checksum_paths([out])),
            participant_count=out_df.get_column("participant_id").n_unique(),
            feature_count=len(s.variables), source_id=s.id,
            temporal_type=TemporalType.EVENT if ev else TemporalType.STATIC,
            parent_artifact_ids=[src.artifact_id] if src and src.artifact_id else [],
            tool="efgpp", configuration_hash=self.configuration_hash(),
            metadata={"variables": kinds},
        ))
        has_any = out_df.select(pl.any_horizontal([pl.col(v).is_not_null() for v in s.variables])).to_series()
        present = out_df.filter(has_any)
        self.record_assays(art, present.get_column("participant_id"),
                           event_ids=present.get_column("event_id") if ev else None)
        return art

    def qc(self) -> QCResult | None:
        art = self.latest("standardized")
        if art is None:
            return None
        df = pl.read_parquet(art.path)
        kinds: dict[str, str] = art.metadata.get("variables", {})
        per_var: dict[str, Any] = {}
        warn = []
        for v, kind in kinds.items():
            col = df.get_column(v)
            miss = col.null_count() / max(df.height, 1)
            entry: dict[str, Any] = {"kind": kind, "missing_fraction": round(miss, 4)}
            if kind == "numeric" and col.drop_nulls().len():
                entry.update({"mean": col.mean(), "sd": col.std(), "min": col.min(), "max": col.max()})
                if not col.std():
                    warn.append(f"{v} is constant")
            elif kind == "categorical":
                entry["levels"] = {str(k): int(n) for k, n in col.drop_nulls().value_counts(sort=True).iter_rows()}
            if miss > MISSING_WARN:
                warn.append(f"{v}: {miss:.1%} missing")
            per_var[v] = entry
        status = ArtifactStatus.QC_WARN if warn else ArtifactStatus.QC_PASS
        result = QCResult(status, {"n_rows": df.height, "variables": per_var}, warn)
        self.save_qc(art, result, {"missing_warn": MISSING_WARN})
        self.store.set_status(art.artifact_id, ArtifactStatus.READY, qc_status=status.value)  # type: ignore[arg-type]
        return result
