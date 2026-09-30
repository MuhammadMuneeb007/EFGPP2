"""Phenotypes: arbitrary user-defined traits of any supported type.

Nothing in this module knows any phenotype's name or meaning. Type-specific behaviour
lives in `TypeHandler`s registered in `TYPE_HANDLERS`; future types (survival,
longitudinal, time-to-event, repeated measures) plug in the same way.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, ClassVar

import polars as pl
import yaml

from efgpp.config.data import PhenotypeSource
from efgpp.constants import (
    SUPPORTED_PHENOTYPE_TYPES,
    ArtifactStatus,
    Modality,
    PhenotypeType,
    TemporalType,
)
from efgpp.data.adapters.base import QCResult, TabularAdapter
from efgpp.data.artifacts import Artifact
from efgpp.data.io import as_text, missing_mask, write_parquet
from efgpp.data.schemas.phenotype import PHENOTYPE_OBSERVATION_SCHEMA
from efgpp.data.timeline import events_from_spec, register_events
from efgpp.data.validation import ValidationReport, plural
from efgpp.resources.checksums import checksum_paths

MIN_CLASS_COUNT = 10


@dataclass
class Coded:
    numeric: pl.Series  # Float64, null where missing/invalid
    text: pl.Series  # Utf8 canonical label, null where missing/invalid
    invalid: pl.Series  # Boolean
    notes: list[str]


class TypeHandler:
    type: ClassVar[PhenotypeType]

    def code(self, values: pl.Series, missing: pl.Series, spec: PhenotypeSource) -> Coded:
        raise NotImplementedError

    def summarize(self, numeric: pl.Series, text: pl.Series) -> dict[str, Any]:
        raise NotImplementedError

    def qc_status(self, summary: dict[str, Any]) -> tuple[ArtifactStatus, list[str]]:
        return ArtifactStatus.QC_PASS, []


class BinaryHandler(TypeHandler):
    type = PhenotypeType.BINARY

    def coding(self, observed: set[str], spec: PhenotypeSource) -> tuple[set[str], set[str], list[str]]:
        if spec.case_values is not None and spec.control_values is not None:
            return set(spec.case_values), set(spec.control_values), []
        numeric = set()
        for v in observed:
            try:
                numeric.add(float(v))
            except ValueError:
                numeric.add(math.nan)
        if numeric <= {1.0, 2.0} and 2.0 in numeric:
            return (
                {v for v in observed if float(v) == 2.0},
                {v for v in observed if float(v) == 1.0},
                ["PLINK coding detected (1 = control, 2 = case); set case_values/control_values to override"],
            )
        return (
            {v for v in observed if _as_float(v) == 1.0},
            {v for v in observed if _as_float(v) == 0.0},
            [],
        )

    def code(self, values: pl.Series, missing: pl.Series, spec: PhenotypeSource) -> Coded:
        observed = set(values.filter(~missing).unique().to_list())
        cases, controls, notes = self.coding(observed, spec)
        numeric = [
            None if m else (1.0 if v in cases else 0.0 if v in controls else None)
            for v, m in zip(values.to_list(), missing.to_list(), strict=True)
        ]
        num = pl.Series(numeric, dtype=pl.Float64)
        invalid = (~missing) & num.is_null()
        text = num.replace_strict({1.0: "case", 0.0: "control"}, default=None, return_dtype=pl.Utf8)
        return Coded(num, text, invalid, notes)

    def summarize(self, numeric: pl.Series, text: pl.Series) -> dict[str, Any]:
        cases = int((numeric == 1.0).sum())
        controls = int((numeric == 0.0).sum())
        total = cases + controls
        return {"cases": cases, "controls": controls,
                "case_fraction": round(cases / total, 4) if total else None}

    def qc_status(self, summary: dict[str, Any]) -> tuple[ArtifactStatus, list[str]]:
        if summary["cases"] == 0 or summary["controls"] == 0:
            return ArtifactStatus.QC_FAIL, ["only one class present"]
        if min(summary["cases"], summary["controls"]) < MIN_CLASS_COUNT:
            return ArtifactStatus.QC_WARN, [f"fewer than {MIN_CLASS_COUNT} observations in a class"]
        return ArtifactStatus.QC_PASS, []


class ContinuousHandler(TypeHandler):
    type = PhenotypeType.CONTINUOUS

    def code(self, values: pl.Series, missing: pl.Series, spec: PhenotypeSource) -> Coded:
        num = values.cast(pl.Float64, strict=False)
        num = pl.Series([None if m else v for v, m in zip(num.to_list(), missing.to_list(), strict=True)], dtype=pl.Float64)
        num = num.fill_nan(None)
        invalid = (~missing) & num.is_null()
        return Coded(num, pl.Series([None] * values.len(), dtype=pl.Utf8), invalid, [])

    def summarize(self, numeric: pl.Series, text: pl.Series) -> dict[str, Any]:
        x = numeric.drop_nulls()
        if x.len() == 0:
            return {"n": 0}
        mean, sd = x.mean(), x.std()
        out = {
            "n": x.len(), "mean": mean, "sd": sd, "min": x.min(), "max": x.max(),
            "median": x.median(), "q01": x.quantile(0.01), "q99": x.quantile(0.99),
            "skewness": x.skew(),
        }
        if sd:
            out["n_beyond_5sd"] = int(((x - mean).abs() > 5 * sd).sum())  # type: ignore[operator]
        return out

    def qc_status(self, summary: dict[str, Any]) -> tuple[ArtifactStatus, list[str]]:
        if summary.get("n", 0) == 0:
            return ArtifactStatus.QC_FAIL, ["no numeric observations"]
        if not summary.get("sd"):
            return ArtifactStatus.QC_FAIL, ["phenotype is constant"]
        if summary.get("n_beyond_5sd", 0) > 0:
            return ArtifactStatus.QC_WARN, [f"{summary['n_beyond_5sd']} values beyond 5 SD (reported, not removed)"]
        return ArtifactStatus.QC_PASS, []


class CategoricalHandler(TypeHandler):
    type = PhenotypeType.MULTICLASS

    def levels(self, values: pl.Series, missing: pl.Series, spec: PhenotypeSource) -> list[str]:
        if spec.levels:
            return [str(v) for v in spec.levels]
        return sorted(values.filter(~missing).unique().to_list())

    def code(self, values: pl.Series, missing: pl.Series, spec: PhenotypeSource) -> Coded:
        levels = self.levels(values, missing, spec)
        index = {lv: float(i) for i, lv in enumerate(levels)}
        text = pl.Series([None if m or v not in index else v for v, m in zip(values.to_list(), missing.to_list(), strict=True)], dtype=pl.Utf8)
        num = text.replace_strict(index, default=None, return_dtype=pl.Float64)
        invalid = (~missing) & text.is_null()
        return Coded(num, text, invalid, [f"levels: {levels}"])

    def summarize(self, numeric: pl.Series, text: pl.Series) -> dict[str, Any]:
        counts = text.drop_nulls().value_counts(sort=True)
        return {"classes": {r[0]: int(r[1]) for r in counts.iter_rows()}}

    def qc_status(self, summary: dict[str, Any]) -> tuple[ArtifactStatus, list[str]]:
        classes = summary["classes"]
        if len(classes) < 2:
            return ArtifactStatus.QC_FAIL, ["fewer than two classes present"]
        small = [k for k, v in classes.items() if v < MIN_CLASS_COUNT]
        if small:
            return ArtifactStatus.QC_WARN, [f"classes with < {MIN_CLASS_COUNT} observations: {small}"]
        return ArtifactStatus.QC_PASS, []


class OrdinalHandler(CategoricalHandler):
    type = PhenotypeType.ORDINAL  # numeric code = rank in the configured level order


TYPE_HANDLERS: dict[PhenotypeType, TypeHandler] = {
    h.type: h for h in (BinaryHandler(), ContinuousHandler(), CategoricalHandler(), OrdinalHandler())
}


def _as_float(v: str) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


class PhenotypeAdapter(TabularAdapter):
    modality = Modality.PHENOTYPE
    temporal_type = TemporalType.STATIC
    source: PhenotypeSource

    def subject(self) -> str:
        return f"PHENOTYPE {self.source.id} ({self.source.name}, {self.source.type.value})"

    def modality_qc_dir(self) -> str:
        return "phenotype"

    def _event_column(self) -> str | None:
        return self.source.timeline.event_column if self.source.timeline else None

    def validate(self) -> ValidationReport:
        s = self.source
        rep = ValidationReport(self.subject(), s.id)
        path = self.source_path
        if not path.exists():
            rep.fail(f"file not found: {path}")
            return rep
        try:
            df = self.read()
        except Exception as exc:  # noqa: BLE001 - report every unreadable input
            rep.fail(f"cannot read table: {exc}")
            return rep
        cols_ok = True
        for label, col in (("participant ID", s.participant_id_column), ("value", s.value_column)):
            if col in df.columns:
                rep.ok(f"{label} column exists ({col})")
            else:
                rep.fail(f"{label} column {col!r} not found; columns: {df.columns[:10]}")
                cols_ok = False
        ev = self._event_column()
        if ev and ev not in df.columns:
            rep.fail(f"event column {ev!r} not found")
            cols_ok = False
        if s.type in SUPPORTED_PHENOTYPE_TYPES:
            rep.ok(f"phenotype type recognized ({s.type.value})")
        else:
            rep.fail(f"phenotype type {s.type.value!r} is planned but not supported in v0.1")
            return rep
        if not cols_ok:
            return rep

        ids = as_text(df, s.participant_id_column)
        n_missing_ids = int((ids.is_null() | (ids == "")).sum())
        if n_missing_ids:
            rep.fail(f"{plural(n_missing_ids, 'row')} without a participant ID", n_missing_ids)
        keys = df.select(ids.alias("id"), *( [as_text(df, ev).alias("ev")] if ev else []))
        dup = keys.filter(keys.is_duplicated())
        if dup.height:
            n = dup.get_column("id").n_unique()
            rep.fail(f"{plural(n, 'duplicated participant ID')}", n, dup.get_column("id").unique().to_list())
        else:
            rep.ok("participant IDs are unique" + (" per event" if ev else ""))

        values = as_text(df, s.value_column)
        missing = missing_mask(values, s.missing_values)
        coded = TYPE_HANDLERS[s.type].code(values, missing, s)
        n_invalid = int(coded.invalid.sum())
        if n_invalid:
            bad = values.filter(coded.invalid).unique().to_list()
            rep.fail(f"{plural(n_invalid, f'invalid {s.type.value} value')}", n_invalid, bad)
        else:
            rep.ok(f"all non-missing values are valid {s.type.value} values")
        for note in coded.notes:
            rep.info(note)
        n_miss = int(missing.sum())
        if n_miss:
            rep.warn(f"{plural(n_miss, 'missing value')}", n_miss)
        return rep

    def standardize(self) -> Artifact:
        s = self.source
        df = self.read()
        ids = as_text(df, s.participant_id_column)
        res = self.resolve_ids(ids)
        lookup = dict(res.mapping.iter_rows())
        pids = ids.replace_strict(lookup, default=None)
        values = as_text(df, s.value_column)
        missing = missing_mask(values, s.missing_values)
        coded = TYPE_HANDLERS[s.type].code(values, missing, s)
        ev = self._event_column()
        obs = pl.DataFrame({
            "participant_id": pids,
            "event_id": as_text(df, ev) if ev else pl.Series([None] * df.height, dtype=pl.Utf8),
            "value_numeric": coded.numeric,
            "value_text": coded.text,
            "is_missing": coded.numeric.is_null() & coded.text.is_null(),
        }).filter(pl.col("participant_id").is_not_null())
        obs = obs.unique(subset=["participant_id", "event_id"], keep="first", maintain_order=True)
        obs = PHENOTYPE_OBSERVATION_SCHEMA.validate(obs)

        events = events_from_spec(df, pids, s.timeline, s.id)
        if events is not None:
            register_events(self.reg, events.filter(pl.col("participant_id").is_not_null()))

        pdir = self.project.phenotype_dir(s.id)
        out = write_parquet(obs, pdir / "observations.parquet", self.project.config.storage.compression)
        # The phenotype as in the source file: original column name and values (for feature engineering).
        cols = {"participant_id": pids, s.value_column: pl.Series(
            [None if m else v for v, m in zip(values.to_list(), missing.to_list(), strict=True)], dtype=pl.Utf8)}
        if ev:
            cols[ev] = as_text(df, ev)
        original = pl.DataFrame(cols).filter(pl.col("participant_id").is_not_null()).unique(
            ["participant_id", *([ev] if ev else [])], keep="first", maintain_order=True)
        if coded.numeric.is_not_null().any() and s.type == PhenotypeType.CONTINUOUS:
            original = original.with_columns(pl.col(s.value_column).cast(pl.Float64, strict=False))
        write_parquet(original, pdir / "phenotype.parquet", self.project.config.storage.compression)
        levels = coded.notes[0].removeprefix("levels: ") if s.type in (PhenotypeType.MULTICLASS, PhenotypeType.ORDINAL) else None
        definition = {
            "phenotype_id": s.id, "name": s.name, "type": s.type.value,
            "source_path": s.path, "value_column": s.value_column, "units": s.units,
            "levels": s.levels, "case_values": s.case_values, "control_values": s.control_values,
            "ontology_term": s.ontology_term.model_dump() if s.ontology_term else None,
            "coding_notes": coded.notes,
        }
        (pdir / "definition.yaml").write_text(yaml.safe_dump(definition, sort_keys=False), encoding="utf-8")

        src = self.source_artifact()
        n_obs = int((~obs.get_column("is_missing")).sum())
        art = self.register_replacing(Artifact(
            artifact_name=f"{s.id}_observations", artifact_type="standardized",
            modality=Modality.PHENOTYPE, origin=s.origin, status=ArtifactStatus.STANDARDIZED,
            path=str(out), format="parquet", size=out.stat().st_size,
            checksum=str(checksum_paths([out])), participant_count=n_obs, feature_count=1,
            source_id=s.id, temporal_type=TemporalType.EVENT if ev else TemporalType.STATIC,
            parent_artifact_ids=[src.artifact_id] if src and src.artifact_id else [],
            tool="efgpp", configuration_hash=self.configuration_hash(),
            metadata={"phenotype_name": s.name, "phenotype_type": s.type.value, "levels": levels},
        ))
        self.reg.upsert("phenotype_definitions", {
            "phenotype_id": s.id, "name": s.name, "type": s.type.value, "source_id": s.id,
            "value_column": s.value_column, "units": s.units, "levels": s.levels,
            "ontology_term": definition["ontology_term"], "definition": definition,
            "artifact_id": art.artifact_id,
        })
        self.reg.execute("DELETE FROM phenotype_observations WHERE phenotype_id = ?", [s.id])
        self.reg.insert_frame(
            "phenotype_observations",
            obs.with_columns(pl.lit(s.id).alias("phenotype_id")),
            replace=False,
        )
        observed = obs.filter(~pl.col("is_missing"))
        self.record_assays(art, observed.get_column("participant_id"), event_ids=observed.get_column("event_id"))
        return art

    def qc(self) -> QCResult | None:
        art = self.latest("standardized")
        if art is None:
            return None
        obs = pl.read_parquet(art.path)
        handler = TYPE_HANDLERS[self.source.type]
        summary = handler.summarize(obs.get_column("value_numeric"), obs.get_column("value_text"))
        n_missing = int(obs.get_column("is_missing").sum())
        metrics = {
            "type": self.source.type.value, "n_rows": obs.height, "n_missing": n_missing,
            "missing_fraction": round(n_missing / obs.height, 4) if obs.height else None,
            **summary,
        }
        status, messages = handler.qc_status(summary)
        if status == ArtifactStatus.QC_PASS and metrics["missing_fraction"] and metrics["missing_fraction"] > 0.2:
            status, messages = ArtifactStatus.QC_WARN, [*messages, "more than 20% missing"]
        result = QCResult(status, metrics, messages)
        self.save_qc(art, result, {"min_class_count": MIN_CLASS_COUNT})
        if status in (ArtifactStatus.QC_PASS, ArtifactStatus.QC_WARN):
            self.store.set_status(art.artifact_id, ArtifactStatus.READY, qc_status=status.value)  # type: ignore[arg-type]
        return result

    def summarize(self) -> dict[str, Any]:
        out = super().summarize()
        out.update({"name": self.source.name, "type": self.source.type.value})
        row = self.reg.one(
            "SELECT summary FROM qc_runs WHERE artifact_id = ? ORDER BY created_at DESC LIMIT 1",
            [out["artifact_id"]],
        )
        if row:
            out["qc"] = row["summary"]
        return out

