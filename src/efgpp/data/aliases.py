"""Mapping native sample identifiers to canonical participant IDs.

Datasets are never joined by row position. Each source's native IDs are resolved to a
canonical `participant_id` through (in order): alias files configured for that source,
alias files configured for every source ("*"), and finally the identity mapping.
"""

from __future__ import annotations

from dataclasses import dataclass

import polars as pl

from efgpp.data.io import as_text, read_table
from efgpp.project import Project


@dataclass
class Resolution:
    mapping: pl.DataFrame  # native_id, participant_id (one row per distinct native id)
    n_mapped_by_alias: int
    n_identity: int
    alias_applied: bool


class AliasResolver:
    def __init__(self, project: Project) -> None:
        self.project = project
        self._maps: list[tuple[list[str], dict[str, str]]] = []
        for spec in project.data.participants.aliases:
            df = read_table(project.resolve(spec.path))
            natives = as_text(df, spec.alias_column).to_list()
            canon = as_text(df, spec.participant_id_column).to_list()
            self._maps.append((spec.sources, dict(zip(natives, canon, strict=True))))

    def _applicable(self, source_id: str) -> list[dict[str, str]]:
        specific = [m for srcs, m in self._maps if source_id in srcs]
        wildcard = [m for srcs, m in self._maps if "*" in srcs]
        return specific + wildcard

    def resolve(self, source_id: str, native_ids: pl.Series) -> Resolution:
        natives = native_ids.cast(pl.Utf8).str.strip_chars().drop_nulls().unique(maintain_order=True)
        maps = self._applicable(source_id)
        resolved, by_alias = [], 0
        for nid in natives.to_list():
            pid = next((m[nid] for m in maps if nid in m), None)
            if pid is not None:
                by_alias += 1
            resolved.append(pid if pid is not None else nid)
        mapping = pl.DataFrame(
            {"native_id": natives.to_list(), "participant_id": resolved},
            schema={"native_id": pl.Utf8, "participant_id": pl.Utf8},
        )
        return Resolution(mapping, by_alias, len(resolved) - by_alias, bool(maps))


def plink_native_ids(samples: pl.DataFrame, mode: str, column: str, separator: str) -> pl.Series:
    """Native ID for each PLINK sample row according to the configured ID mode."""
    if mode == "fid_iid":
        fid = samples.get_column("FID").fill_null("0")
        return (fid + separator + samples.get_column("IID")).alias("native_id")
    return samples.get_column(column).alias("native_id")
