"""Biospecimens: participant -> event -> biospecimen -> assay."""

from __future__ import annotations

import polars as pl

from efgpp.data.aliases import AliasResolver
from efgpp.data.io import as_text, read_table
from efgpp.data.registry import Registry
from efgpp.data.schemas.biospecimen import BIOSPECIMEN_SCHEMA
from efgpp.project import Project

COLUMNS = list(BIOSPECIMEN_SCHEMA.columns)


def _col(df: pl.DataFrame, name: str | None) -> pl.Series:
    if not name:
        return pl.Series([None] * df.height, dtype=pl.Utf8)
    return as_text(df, name)


def load_biospecimens_file(project: Project, reg: Registry) -> int:
    spec = project.data.biospecimens
    if spec is None:
        return 0
    df = read_table(project.resolve(spec.path))
    res = AliasResolver(project).resolve("biospecimens", as_text(df, spec.participant_id_column))
    lookup = dict(res.mapping.iter_rows())
    out = pl.DataFrame(
        {
            "biospecimen_id": as_text(df, spec.biospecimen_id_column),
            "participant_id": as_text(df, spec.participant_id_column).replace_strict(lookup, default=None),
            "event_id": _col(df, spec.event_id_column),
            "tissue": _col(df, spec.tissue_column),
            "material": _col(df, spec.material_column),
            "collection_date": _col(df, spec.collection_date_column),
            "processing_method": _col(df, spec.processing_method_column),
            "storage_condition": _col(df, spec.storage_condition_column),
            "source_id": pl.Series(["biospecimens"] * df.height, dtype=pl.Utf8),
        }
    )
    out = BIOSPECIMEN_SCHEMA.validate(out)
    reg.execute("DELETE FROM biospecimens WHERE source_id = 'biospecimens'")
    return reg.insert_frame("biospecimens", out)


def register_inline_biospecimens(
    reg: Registry,
    *,
    biospecimen_ids: pl.Series,
    participant_ids: pl.Series,
    event_ids: pl.Series | None,
    tissue: str | None,
    source_id: str,
) -> int:
    """Create biospecimen rows referenced by an assay table but absent from the registry."""
    n = biospecimen_ids.len()
    df = pl.DataFrame(
        {
            "biospecimen_id": biospecimen_ids.cast(pl.Utf8),
            "participant_id": participant_ids.cast(pl.Utf8),
            "event_id": event_ids.cast(pl.Utf8) if event_ids is not None else pl.Series([None] * n, dtype=pl.Utf8),
            "tissue": pl.Series([tissue] * n, dtype=pl.Utf8),
            "material": pl.Series([None] * n, dtype=pl.Utf8),
            "collection_date": pl.Series([None] * n, dtype=pl.Utf8),
            "processing_method": pl.Series([None] * n, dtype=pl.Utf8),
            "storage_condition": pl.Series([None] * n, dtype=pl.Utf8),
            "source_id": pl.Series([source_id] * n, dtype=pl.Utf8),
        }
    ).filter(pl.col("biospecimen_id").is_not_null()).unique("biospecimen_id", keep="first")
    existing = reg.frame("SELECT biospecimen_id FROM biospecimens")
    if existing.height:
        df = df.join(existing, on="biospecimen_id", how="anti")
    return reg.insert_frame("biospecimens", BIOSPECIMEN_SCHEMA.validate(df), replace=False)
