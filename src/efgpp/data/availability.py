"""Which participants have which modalities (participant x modality availability)."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Any

import polars as pl

from efgpp.constants import OMICS_MODALITIES, ArtifactStatus, Modality, Origin
from efgpp.data.io import write_parquet
from efgpp.data.registry import Registry
from efgpp.project import Project

ORIGIN_ORDER = [o.value for o in (Origin.OBSERVED, Origin.PREDICTED, Origin.DERIVED, Origin.SIMULATED)]
OMICS_VALUES = {m.value for m in OMICS_MODALITIES}


@dataclass
class Availability:
    matrix: pl.DataFrame  # participant_id + one boolean column per modality source
    columns: pl.DataFrame  # column, source_id, modality, origin, tissue, artifact_id, participants
    participant_level_free: list[dict[str, Any]] = field(default_factory=list)  # e.g. annotations

    @property
    def n_participants(self) -> int:
        return self.matrix.height

    def column_names(self) -> list[str]:
        return self.columns.get_column("column").to_list()

    def count(self, *cols: str) -> int:
        if not cols:
            return self.matrix.height
        return int(self.matrix.select(pl.all_horizontal([pl.col(c) for c in cols])).to_series().sum())

    def participants_with(self, *cols: str) -> list[str]:
        if not cols:
            return self.matrix.get_column("participant_id").to_list()
        return self.matrix.filter(pl.all_horizontal([pl.col(c) for c in cols])).get_column("participant_id").to_list()

    def patterns(self) -> pl.DataFrame:
        """Every observed combination of modalities with its participant count (UpSet data)."""
        cols = self.column_names()
        if not cols:
            return pl.DataFrame()
        return self.matrix.group_by(cols).len("participants").sort("participants", descending=True)

    def pairwise(self) -> pl.DataFrame:
        cols = self.column_names()
        rows = [{"a": a, "b": b, "participants": self.count(a, b)} for a, b in combinations(cols, 2)]
        return pl.DataFrame(rows) if rows else pl.DataFrame(schema={"a": pl.Utf8, "b": pl.Utf8, "participants": pl.Int64})

    def key_intersections(self) -> list[dict[str, Any]]:
        """Phenotype-centred intersections: each phenotype with genotype, covariates and omics."""
        meta = {r["column"]: r for r in self.columns.to_dicts()}
        phen = [c for c, m in meta.items() if m["modality"] == Modality.PHENOTYPE.value]
        geno = [c for c, m in meta.items() if m["modality"] == Modality.GENOTYPE.value]
        cov = [c for c, m in meta.items() if m["modality"] == Modality.COVARIATES.value]
        omics = [c for c, m in meta.items() if m["modality"] in OMICS_VALUES]
        combos: list[tuple[str, ...]] = []
        for g in geno:
            combos.append((g,))
            combos += [(g, o) for o in omics]
        for p in phen:
            for g in geno:
                combos.append((g, p))
                if cov:
                    combos.append((g, p, *cov))
            combos += [(p, o) for o in omics]
        combos += list(combinations(omics, 2))
        seen, out = set(), []
        for c in combos:
            if c not in seen:
                seen.add(c)
                out.append({"modalities": " + ".join(c), "participants": self.count(*c)})
        return out


def _column_name(source_id: str, modality: str, origin: str) -> str:
    if modality in OMICS_VALUES and origin != Origin.PREDICTED.value:
        return f"{source_id}_{origin.upper()}"
    return source_id


def compute(project: Project, reg: Registry) -> Availability:
    live = [s.value for s in ArtifactStatus if s not in (ArtifactStatus.SUPERSEDED, ArtifactStatus.QC_FAIL)]
    marks = ", ".join("?" for _ in live)
    groups = reg.frame(
        f"""
        SELECT a.source_id, a.modality, a.origin, max(a.tissue) AS tissue, max(a.artifact_id) AS artifact_id,
               count(DISTINCT a.participant_id) AS participants
        FROM assays a JOIN artifacts t ON t.artifact_id = a.artifact_id
        WHERE t.status IN ({marks}) AND a.participant_id IS NOT NULL
        GROUP BY a.source_id, a.modality, a.origin
        """,
        live,
    )
    order = {o: i for i, o in enumerate(ORIGIN_ORDER)}
    rows = sorted(groups.to_dicts(), key=lambda r: (order.get(r["origin"], 9), r["source_id"]))
    for r in rows:
        r["column"] = _column_name(r["source_id"], r["modality"], r["origin"])
    columns = pl.DataFrame(rows, schema={"source_id": pl.Utf8, "modality": pl.Utf8, "origin": pl.Utf8,
                                         "tissue": pl.Utf8, "artifact_id": pl.Utf8, "participants": pl.Int64,
                                         "column": pl.Utf8})
    matrix = reg.frame("SELECT participant_id FROM participants ORDER BY participant_id")
    for r in rows:
        ids = reg.frame(
            f"SELECT DISTINCT a.participant_id FROM assays a JOIN artifacts t ON t.artifact_id = a.artifact_id "
            f"WHERE a.source_id = ? AND a.origin = ? AND t.status IN ({marks}) AND a.participant_id IS NOT NULL",
            [r["source_id"], r["origin"], *live],
        ).get_column("participant_id")
        matrix = matrix.with_columns(pl.col("participant_id").is_in(ids.to_list()).alias(r["column"]))
    free = reg.rows(
        "SELECT artifact_id, artifact_name, modality, origin, feature_count FROM artifacts "
        "WHERE modality = ? AND status <> 'SUPERSEDED' ORDER BY artifact_id",
        [Modality.VARIANT_ANNOTATIONS.value],
    )
    return Availability(matrix, columns.select("column", "source_id", "modality", "origin", "tissue", "artifact_id", "participants"), free)


def write(project: Project, reg: Registry, av: Availability) -> None:
    out = project.registry_path.parent
    write_parquet(av.matrix, out / "availability.parquet")
    write_parquet(av.columns, out / "availability_columns.parquet")
    reg.execute("DELETE FROM modality_availability")
    if av.columns.height:
        long = av.matrix.unpivot(index="participant_id", variable_name="column_name", value_name="available")
        long = long.join(av.columns.select(pl.col("column").alias("column_name"), "origin", "modality", "artifact_id"),
                         on="column_name", how="left")
        reg.insert_frame("modality_availability", long.select(
            "participant_id", "column_name", "origin", "modality", "artifact_id", "available"), replace=False)


def build(project: Project) -> Availability:
    with Registry.open(project) as reg:
        av = compute(project, reg)
        write(project, reg, av)
        reg.export_parquet()
    return av


def summary_rows(av: Availability) -> list[tuple[str, list[tuple[str, int]]]]:
    """Counts grouped by origin, for the section-63 style console summary."""
    out = []
    for origin in ORIGIN_ORDER:
        sub = av.columns.filter(pl.col("origin") == origin)
        if sub.height:
            out.append((origin.upper(), [(r["column"], r["participants"]) for r in sub.to_dicts()]))
    return out
