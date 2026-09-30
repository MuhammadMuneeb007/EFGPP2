from __future__ import annotations

import pandera.polars as pa
import polars as pl


def covariate_schema(event_column: bool) -> pa.DataFrameSchema:
    unique = ["participant_id", "event_id"] if event_column else ["participant_id"]
    cols = {"participant_id": pa.Column(pl.Utf8, nullable=False)}
    if event_column:
        cols["event_id"] = pa.Column(pl.Utf8, nullable=True)
    return pa.DataFrameSchema(cols, unique=unique, strict=False, name="covariates")
