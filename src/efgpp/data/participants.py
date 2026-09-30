"""Canonical participant registry."""

from __future__ import annotations

import json

import polars as pl

from efgpp.data.aliases import AliasResolver
from efgpp.data.io import as_text, read_table
from efgpp.data.registry import Registry, utcnow
from efgpp.data.schemas.participant import ALIAS_SCHEMA, PARTICIPANT_SCHEMA


def register_mapping(reg: Registry, source_id: str, mapping: pl.DataFrame) -> dict[str, int]:
    """Record a source's native->canonical mapping and create unseen participants.

    Re-registering a source replaces its aliases (sources are re-read after edits).
    """
    mapping = ALIAS_SCHEMA.validate(
        mapping.with_columns(pl.lit(source_id).alias("source_id")).select(
            "source_id", "native_id", "participant_id"
        )
    )
    reg.execute("DELETE FROM sample_aliases WHERE source_id = ?", [source_id])
    reg.insert_frame("sample_aliases", mapping)

    existing = set(reg.frame("SELECT participant_id FROM participants").get_column("participant_id").to_list())
    new_ids = [p for p in mapping.get_column("participant_id").unique(maintain_order=True).to_list() if p not in existing]
    if new_ids:
        now = utcnow()
        new = pl.DataFrame(
            {
                "participant_id": new_ids,
                "first_source": [source_id] * len(new_ids),
                "created_at": [now] * len(new_ids),
                "attributes": [None] * len(new_ids),
            },
            schema={"participant_id": pl.Utf8, "first_source": pl.Utf8, "created_at": pl.Datetime("us"), "attributes": pl.Utf8},
        )
        reg.insert_frame("participants", PARTICIPANT_SCHEMA.validate(new))
    return {"participants": mapping.get_column("participant_id").n_unique(), "new": len(new_ids)}


def load_master_list(project, reg: Registry) -> int:  # type: ignore[no-untyped-def]
    """Register the optional master participant list with its static attributes."""
    spec = project.data.participants
    if not spec.path:
        return 0
    df = read_table(project.resolve(spec.path))
    ids = as_text(df, spec.id_column)
    resolver = AliasResolver(project)
    res = resolver.resolve("participants", ids)
    register_mapping(reg, "participants", res.mapping)
    attrs = df.drop(spec.id_column)
    for pid, row in zip(ids.to_list(), attrs.iter_rows(named=True), strict=True):
        canonical = res.mapping.filter(pl.col("native_id") == pid).get_column("participant_id").item()
        reg.execute(
            "UPDATE participants SET attributes = ? WHERE participant_id = ?",
            [json.dumps(row, default=str), canonical],
        )
    return len(ids)


def participant_ids(reg: Registry) -> list[str]:
    return reg.frame("SELECT participant_id FROM participants ORDER BY participant_id").get_column("participant_id").to_list()


def aliases_for(reg: Registry, source_id: str) -> pl.DataFrame:
    return reg.frame(
        "SELECT native_id, participant_id FROM sample_aliases WHERE source_id = ?", [source_id]
    )
