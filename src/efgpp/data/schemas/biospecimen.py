from __future__ import annotations

import pandera.polars as pa
import polars as pl

BIOSPECIMEN_SCHEMA = pa.DataFrameSchema(
    {
        "biospecimen_id": pa.Column(pl.Utf8, nullable=False, unique=True),
        "participant_id": pa.Column(pl.Utf8, nullable=False),
        "event_id": pa.Column(pl.Utf8, nullable=True),
        "tissue": pa.Column(pl.Utf8, nullable=True),
        "material": pa.Column(pl.Utf8, nullable=True),
        "collection_date": pa.Column(pl.Utf8, nullable=True),
        "processing_method": pa.Column(pl.Utf8, nullable=True),
        "storage_condition": pa.Column(pl.Utf8, nullable=True),
        "source_id": pa.Column(pl.Utf8, nullable=True),
    },
    strict=True,
    name="biospecimens",
)
