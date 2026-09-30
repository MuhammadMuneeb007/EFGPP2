from __future__ import annotations

import pandera.polars as pa
import polars as pl

# Standardized phenotype observations: one row per participant (per event, if timed).
PHENOTYPE_OBSERVATION_SCHEMA = pa.DataFrameSchema(
    {
        "participant_id": pa.Column(pl.Utf8, nullable=False),
        "event_id": pa.Column(pl.Utf8, nullable=True),
        "value_numeric": pa.Column(pl.Float64, nullable=True),
        "value_text": pa.Column(pl.Utf8, nullable=True),
        "is_missing": pa.Column(pl.Boolean, nullable=False),
    },
    unique=["participant_id", "event_id"],
    strict=True,
    name="phenotype_observations",
)
