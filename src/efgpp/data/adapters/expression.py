"""Processed omics matrices (expression, methylation, proteomics, metabolomics).

Accepted inputs: CSV / TSV / Parquet (samples x features or features x samples), AnnData
(.h5ad) and AnnData Zarr stores. Standardized output: AnnData Zarr (or Parquet when
`storage.omics_format: parquet` or anndata is not installed) whose `obs` carries the
canonical participant_id and timeline keys and whose `uns["efgpp"]` carries metadata.

The Data layer profiles these matrices; it never normalises, scales, imputes or filters
them (that would be fitted on all participants and leak into prediction).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar

import numpy as np
import polars as pl

from efgpp.config.data import OmicsSource
from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.adapters.base import DataAdapter, QCResult
from efgpp.data.artifacts import Artifact
from efgpp.data.biospecimens import register_inline_biospecimens
from efgpp.data.io import as_text, detect_table_format, read_table, write_parquet
from efgpp.data.timeline import events_frame, register_events
from efgpp.data.validation import ValidationReport, plural
from efgpp.resources.checksums import checksum_paths


def anndata_available() -> bool:
    try:
        import anndata  # noqa: F401
    except ImportError:
        return False
    return True


@dataclass
class OmicsMatrix:
    X: np.ndarray  # samples x features, float32 (NaN = missing)
    obs: pl.DataFrame  # native_id, event_id, time, biospecimen_id, batch, sample_key
    features: list[str]
    non_numeric: int = 0


class OmicsAdapter(DataAdapter):
    """Common implementation; subclasses declare modality-specific metadata expectations."""

    temporal_type = TemporalType.EVENT
    source: OmicsSource
    recommended_metadata: ClassVar[tuple[str, ...]] = (
        "tissue", "feature_id_system", "measurement_type", "normalization",
    )

    def members(self) -> tuple[Path, list[Path], str]:
        p = self.source_path
        return p, [p], detect_table_format(p, self.source.format)

    def _fmt(self, path: Path) -> str:
        return detect_table_format(path, self.source.format)

    def _current_path(self) -> Path:
        art = self.source_artifact()
        return Path(art.path) if art else self.source_path

    # ------------------------------------------------------------------ reading
    def read_matrix(self, path: Path | None = None) -> OmicsMatrix:
        path = path or self._current_path()
        fmt = self._fmt(path)
        s = self.source
        tl = s.timeline
        if fmt in ("h5ad", "zarr"):
            import anndata as ad

            a: Any = ad.read_h5ad(path) if fmt == "h5ad" else ad.read_zarr(path)
            obs = pl.from_pandas(a.obs.reset_index(names="_obs_index"))
            X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
            id_col = s.participant_id_column if s.participant_id_column in obs.columns else "_obs_index"
            meta = self._obs_meta(obs, id_col, sample_key=obs.get_column("_obs_index").cast(pl.Utf8))
            return OmicsMatrix(X.astype(np.float32), meta, [str(v) for v in a.var_names])

        df = read_table(path, fmt)
        if s.orientation == "features_by_samples":
            fid = s.feature_id_column
            sample_cols = [c for c in df.columns if c != fid]
            features = as_text(df, fid).to_list()  # type: ignore[arg-type]
            vals = df.select(sample_cols).cast(pl.Float64, strict=False)
            non_num = sum(int((df.get_column(c).is_not_null() & vals.get_column(c).is_null()).sum()) for c in sample_cols) if df.schema[sample_cols[0]] == pl.Utf8 else 0
            X = vals.to_numpy().T.astype(np.float32)
            obs = pl.DataFrame({"native_id": sample_cols}, schema={"native_id": pl.Utf8})
            meta = self._obs_meta(obs, "native_id", sample_key=obs.get_column("native_id"))
            return OmicsMatrix(X, meta, features, non_num)

        meta_cols = [c for c in (
            s.participant_id_column, s.batch_column, s.biospecimen_column,
            *( [tl.event_column, tl.time_column, tl.study_day_column, tl.visit_name_column, tl.age_column] if tl else [] ),
        ) if c]
        feature_cols = [c for c in df.columns if c not in meta_cols]
        vals = df.select(feature_cols).cast(pl.Float64, strict=False)
        non_num = 0
        for c in feature_cols:
            if df.schema[c] == pl.Utf8:
                raw = df.get_column(c).str.strip_chars()
                non_num += int((raw.is_not_null() & (raw != "") & ~raw.is_in(["NA", "NaN", "nan"]) & vals.get_column(c).is_null()).sum())
        X = vals.to_numpy().astype(np.float32) if feature_cols else np.zeros((df.height, 0), np.float32)
        meta = self._obs_meta(df, s.participant_id_column, sample_key=None)
        return OmicsMatrix(X, meta, feature_cols, non_num)

    def _obs_meta(self, obs: pl.DataFrame, id_col: str, sample_key: pl.Series | None) -> pl.DataFrame:
        s = self.source
        tl = s.timeline
        n = obs.height

        def col(name: str | None) -> pl.Series:
            if name and name in obs.columns:
                return as_text(obs, name)
            return pl.Series([None] * n, dtype=pl.Utf8)

        return pl.DataFrame({
            "native_id": as_text(obs, id_col) if id_col in obs.columns else pl.Series([None] * n, dtype=pl.Utf8),
            "event_id": col(tl.event_column if tl else None),
            "time": col(tl.time_column if tl else None),
            "study_day": col(tl.study_day_column if tl else None),
            "biospecimen_id": col(s.biospecimen_column),
            "batch": col(s.batch_column),
            "sample_key": sample_key if sample_key is not None else pl.Series([None] * n, dtype=pl.Utf8),
        })

    # --------------------------------------------------------------- lifecycle
    def inspect(self) -> dict[str, Any]:
        path = self.source_path
        info: dict[str, Any] = {"source_id": self.source.id, "modality": self.modality.value,
                                "origin": self.source.origin.value, "path": str(path), "exists": path.exists()}
        if path.exists():
            m = self.read_matrix(path)
            info.update({
                "format": self._fmt(path), "samples": m.X.shape[0], "features": m.X.shape[1],
                "participants": m.obs.get_column("native_id").n_unique(),
                "events": m.obs.get_column("event_id").drop_nulls().n_unique(),
                "tissue": self.source.tissue,
            })
        return info

    def validate(self) -> ValidationReport:
        s = self.source
        rep = ValidationReport(self.subject(), s.id)
        path = self.source_path
        if not path.exists():
            rep.fail(f"not found: {path}")
            return rep
        fmt = self._fmt(path)
        if fmt in ("h5ad", "zarr") and not anndata_available():
            rep.fail(f"{fmt} input requires the 'anndata' package (pip install 'efgpp[omics]')")
            return rep
        if fmt not in ("h5ad", "zarr") and s.orientation == "samples_by_features":
            cols = read_table(path, fmt).columns
            if s.participant_id_column not in cols:
                rep.fail(f"participant ID column {s.participant_id_column!r} not found")
                return rep
            if s.timeline and s.timeline.event_column and s.timeline.event_column not in cols:
                rep.fail(f"event column {s.timeline.event_column!r} not found")
                return rep
        try:
            m = self.read_matrix(path)
        except Exception as exc:  # noqa: BLE001
            rep.fail(f"cannot read matrix: {exc}")
            return rep
        n_samples, n_features = m.X.shape
        rep.ok(f"{plural(n_samples, 'sample')} x {plural(n_features, 'feature')}")
        if n_features == 0:
            rep.fail("no feature columns")
        if m.non_numeric:
            rep.fail(f"{plural(m.non_numeric, 'non-numeric value')} in feature columns", m.non_numeric)
        if m.obs.get_column("native_id").null_count():
            rep.fail(f"{plural(m.obs.get_column('native_id').null_count(), 'sample')} without participant ID")
        keys = m.obs.select("native_id", "event_id")
        n_dup = int(keys.is_duplicated().sum())
        if n_dup:
            hint = "" if s.timeline else " (longitudinal data? set --event-column)"
            rep.fail(f"{plural(n_dup, 'duplicated participant/event sample')}{hint}", n_dup)
        else:
            rep.ok("one sample per participant" + (" and event" if s.timeline else ""))
        dup_features = n_features - len(set(m.features))
        if dup_features:
            rep.fail(f"{plural(dup_features, 'duplicated feature identifier')}", dup_features)
        if n_samples and n_features:
            all_missing = int(np.isnan(m.X).all(axis=0).sum())
            if all_missing:
                rep.warn(f"{plural(all_missing, 'feature')} entirely missing", all_missing)
        for attr in self.recommended_metadata:
            if getattr(s, attr) in (None, ""):
                rep.warn(f"metadata `{attr}` not set")
        rep.extend(self.modality_checks(m))
        return rep

    def modality_checks(self, m: OmicsMatrix) -> ValidationReport:
        rep = ValidationReport(self.subject())
        mt = (self.source.measurement_type or "").lower()
        if mt in ("counts", "raw_counts") and m.X.size and np.nanmin(m.X) < 0:
            rep.fail("negative values in a count matrix")
        return rep

    def standardize(self) -> Artifact:
        s = self.source
        m = self.read_matrix()
        res = self.resolve_ids(m.obs.get_column("native_id"))
        pids = m.obs.get_column("native_id").replace_strict(dict(res.mapping.iter_rows()), default=None)
        # Biospecimen IDs that are already registered can carry participant/event information.
        obs = m.obs.with_columns(pids.alias("participant_id"))
        known = self.reg.frame("SELECT biospecimen_id, participant_id AS b_pid, event_id AS b_ev FROM biospecimens")
        if known.height:
            key = "biospecimen_id" if obs.get_column("biospecimen_id").drop_nulls().len() else None
            if key is None and s.orientation == "features_by_samples":
                obs = obs.with_columns(pl.col("native_id").alias("biospecimen_id"))
                key = "biospecimen_id"
            if key:
                obs = obs.join(known, on="biospecimen_id", how="left").with_columns(
                    pl.coalesce("b_pid", "participant_id").alias("participant_id"),
                    pl.coalesce("event_id", "b_ev").alias("event_id"),
                ).drop("b_pid", "b_ev")
        tl = s.timeline
        if tl and tl.event_column:
            register_events(self.reg, events_frame(
                obs, obs.get_column("participant_id"), event_column="event_id", source_id=s.id,
                date_column="time" if tl.time_column else None,
                study_day_column="study_day" if tl.study_day_column else None,
            ).filter(pl.col("participant_id").is_not_null()))
        if obs.get_column("biospecimen_id").drop_nulls().len():
            register_inline_biospecimens(
                self.reg, biospecimen_ids=obs.get_column("biospecimen_id"),
                participant_ids=obs.get_column("participant_id"), event_ids=obs.get_column("event_id"),
                tissue=s.tissue, source_id=s.id,
            )
        out = self._write(m, obs)
        src = self.source_artifact()
        temporal = s.temporal_type or (TemporalType.EVENT if tl and tl.event_column else TemporalType.STATIC)
        if s.origin == Origin.PREDICTED:
            temporal = TemporalType.GENETICALLY_PREDICTED_STATIC
        art = self.register_replacing(Artifact(
            artifact_name=f"{s.id}_{self.modality.value}", artifact_type="standardized",
            modality=self.modality, origin=s.origin, status=ArtifactStatus.STANDARDIZED,
            path=str(out), format="zarr" if out.suffix == ".zarr" else "parquet",
            checksum=str(checksum_paths([out])), size=sum(f.stat().st_size for f in out.rglob("*") if f.is_file()) if out.is_dir() else out.stat().st_size,
            participant_count=obs.get_column("participant_id").drop_nulls().n_unique(),
            feature_count=len(m.features), genome_build=s.genome_build, tissue=s.tissue,
            temporal_type=temporal, source_id=s.id,
            parent_artifact_ids=[src.artifact_id] if src and src.artifact_id else [],
            tool="efgpp", configuration_hash=self.configuration_hash(),
            metadata={**self._metadata(), "notes": [n for n in [getattr(self, "_storage_note", None)] if n]},
        ))
        self.record_assays(
            art, obs.get_column("participant_id"), event_ids=obs.get_column("event_id"),
            biospecimen_ids=obs.get_column("biospecimen_id"), times=obs.get_column("time"),
            sample_keys=obs.get_column("sample_key"),
        )  # full obs (before empty columns were dropped for storage)
        return art

    def _metadata(self) -> dict[str, Any]:
        s = self.source
        keys = ("tissue", "feature_id_system", "measurement_type", "normalization", "platform",
                "units", "genome_build", "orientation")
        return {k: getattr(s, k) for k in keys} | {"modality": self.modality.value, "origin": s.origin.value}

    def _write(self, m: OmicsMatrix, obs: pl.DataFrame) -> Path:
        s = self.source
        folder = self.project.artifact_dir(s.origin, self.modality, s.id)
        # Timeline/biospecimen columns that are entirely empty carry no information.
        obs = obs.select([c for c in obs.columns if c in ("participant_id", "native_id") or obs.get_column(c).null_count() < obs.height])
        if self.project.config.storage.omics_format == "zarr" and anndata_available():
            import shutil

            import anndata as ad
            import pandas as pd

            obs_pd = obs.to_pandas().astype("string")
            obs_pd.index = [f"s{i}" for i in range(obs.height)]
            a = ad.AnnData(X=m.X, obs=obs_pd, var=pd.DataFrame(index=pd.Index(m.features, dtype="string")))
            a.uns["efgpp"] = {k: v for k, v in self._metadata().items() if v is not None}
            out = folder / f"{s.id}.zarr"
            shutil.rmtree(out, ignore_errors=True)
            try:
                a.write_zarr(out)
                return out
            except OSError as exc:
                # e.g. Windows MAX_PATH limits inside deep project paths: keep the data, use Parquet.
                shutil.rmtree(out, ignore_errors=True)
                self._storage_note = f"Zarr write failed ({exc.__class__.__name__}); stored as Parquet"
        frame = obs.hstack(pl.DataFrame(m.X, schema=m.features, orient="row"))
        return write_parquet(frame, folder / f"{s.id}.parquet", self.project.config.storage.compression)

    def load_standardized(self) -> OmicsMatrix:
        art = self.latest("standardized")
        if art is None:
            raise RuntimeError(f"{self.source.id} has not been standardized")
        return load_omics_artifact(Path(art.path))

    def qc(self) -> QCResult | None:
        art = self.latest("standardized")
        if art is None:
            return None
        m = load_omics_artifact(Path(art.path))
        X = m.X
        n_s, n_f = X.shape
        miss = np.isnan(X)
        feat_miss = miss.mean(axis=0) if n_s else np.zeros(n_f)
        samp_miss = miss.mean(axis=1) if n_f else np.zeros(n_s)
        finite = X[~miss]
        metrics: dict[str, Any] = {
            "samples": n_s, "features": n_f,
            "participants": m.obs.get_column("participant_id").n_unique() if "participant_id" in m.obs.columns else None,
            "missing_fraction": float(miss.mean()) if X.size else 0.0,
            "features_all_missing": int((feat_miss == 1).sum()),
            "features_over_20pct_missing": int((feat_miss > 0.2).sum()),
            "samples_over_20pct_missing": int((samp_miss > 0.2).sum()),
            "value_min": float(finite.min()) if finite.size else None,
            "value_median": float(np.median(finite)) if finite.size else None,
            "value_max": float(finite.max()) if finite.size else None,
            "zero_fraction": float((finite == 0).mean()) if finite.size else None,
            "negative_fraction": float((finite < 0).mean()) if finite.size else None,
        }
        msgs = []
        if metrics["value_max"] is not None and metrics["value_max"] > 1000 and "log" not in (self.source.normalization or ""):
            msgs.append("large dynamic range; values look untransformed (informational, no action taken)")
        status = ArtifactStatus.QC_PASS
        if metrics["samples_over_20pct_missing"] or metrics["features_all_missing"]:
            status = ArtifactStatus.QC_WARN
            msgs.append("samples or features with substantial missingness")
        if n_s == 0 or n_f == 0:
            status = ArtifactStatus.QC_FAIL
        result = QCResult(status, metrics, msgs)
        self.save_qc(art, result)
        if status != ArtifactStatus.QC_FAIL:
            self.store.set_status(art.artifact_id, ArtifactStatus.READY, qc_status=status.value)  # type: ignore[arg-type]
        return result


def load_omics_artifact(path: Path) -> OmicsMatrix:
    """Load a standardized omics artifact (Zarr or Parquet) written by `OmicsAdapter`."""
    if path.suffix == ".zarr":
        import anndata as ad

        a: Any = ad.read_zarr(path)
        obs = pl.from_pandas(a.obs.reset_index(drop=True).astype("string"))
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        return OmicsMatrix(X.astype(np.float32), obs, [str(v) for v in a.var_names])
    df = pl.read_parquet(path)
    meta = ["native_id", "participant_id", "event_id", "time", "study_day", "biospecimen_id", "batch", "sample_key"]
    obs_cols = [c for c in meta if c in df.columns]
    features = [c for c in df.columns if c not in obs_cols]
    return OmicsMatrix(df.select(features).to_numpy().astype(np.float32), df.select(obs_cols), features)


class ExpressionAdapter(OmicsAdapter):
    modality = Modality.EXPRESSION
