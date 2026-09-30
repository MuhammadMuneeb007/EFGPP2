from __future__ import annotations

import pandera.polars as pa
import polars as pl

EVENT_SCHEMA = pa.DataFrameSchema(
    {
        "participant_id": pa.Column(pl.Utf8, nullable=False),
        "event_id": pa.Column(pl.Utf8, nullable=False),
        "visit_name": pa.Column(pl.Utf8, nullable=True),
        "event_date": pa.Column(pl.Utf8, nullable=True),
        "study_day": pa.Column(pl.Float64, nullable=True),
        "age_at_event": pa.Column(pl.Float64, nullable=True),
        "source_id": pa.Column(pl.Utf8, nullable=True),
    },
    unique=["participant_id", "event_id"],
    strict=True,
    name="events",
)
