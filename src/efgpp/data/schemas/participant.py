from __future__ import annotations

import pandera.polars as pa
import polars as pl

PARTICIPANT_SCHEMA = pa.DataFrameSchema(
    {
        "participant_id": pa.Column(pl.Utf8, nullable=False, unique=True),
        "first_source": pa.Column(pl.Utf8, nullable=True),
    },
    strict=False,
    name="participants",
)

ALIAS_SCHEMA = pa.DataFrameSchema(
    {
        "source_id": pa.Column(pl.Utf8, nullable=False),
        "native_id": pa.Column(pl.Utf8, nullable=False, unique=True),
        "participant_id": pa.Column(pl.Utf8, nullable=False),
    },
    strict=True,
    name="sample_aliases",
)
