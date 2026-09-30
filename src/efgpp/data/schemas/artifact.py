from __future__ import annotations

import pandera.polars as pa
import polars as pl

from efgpp.constants import ArtifactStatus, Origin

ARTIFACT_TABLE_SCHEMA = pa.DataFrameSchema(
    {
        "artifact_id": pa.Column(pl.Utf8, nullable=False, unique=True),
        "origin": pa.Column(pl.Utf8, checks=pa.Check.isin([o.value for o in Origin])),
        "status": pa.Column(pl.Utf8, checks=pa.Check.isin([s.value for s in ArtifactStatus])),
        "path": pa.Column(pl.Utf8, nullable=False),
    },
    strict=False,
    name="artifacts",
)

# Standardized variant table (section 27): human-readable keys plus optional VRS ids.
VARIANT_SCHEMA = pa.DataFrameSchema(
    {
        "genome_build": pa.Column(pl.Utf8, nullable=True),
        "chromosome": pa.Column(pl.Utf8, nullable=False),
        "position": pa.Column(pl.Int64, nullable=False, checks=pa.Check.ge(1)),
        "reference": pa.Column(pl.Utf8, nullable=True),
        "alternate": pa.Column(pl.Utf8, nullable=True),
        "rsid": pa.Column(pl.Utf8, nullable=True),
        "vrs_id": pa.Column(pl.Utf8, nullable=True),
    },
    strict=False,
    name="variants",
)
