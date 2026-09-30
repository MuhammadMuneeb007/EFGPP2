"""Participant timelines: events (visits) and the measurements attached to them."""

from __future__ import annotations

from typing import Any

import polars as pl

from efgpp.config.data import TimelineSpec
from efgpp.data.aliases import AliasResolver
from efgpp.data.io import as_text, read_table
from efgpp.data.registry import Registry
from efgpp.data.schemas.timeline import EVENT_SCHEMA
from efgpp.project import Project

STATIC_EVENT = None  # static measurements (e.g. germline genotype) carry no event


def _num(df: pl.DataFrame, col: str | None) -> pl.Series:
    if not col:
        return pl.Series([None] * df.height, dtype=pl.Float64)
    return as_text(df, col).cast(pl.Float64, strict=False)


def _text(df: pl.DataFrame, col: str | None) -> pl.Series:
    if not col:
        return pl.Series([None] * df.height, dtype=pl.Utf8)
    return as_text(df, col)


def events_frame(
    df: pl.DataFrame,
    participant_ids: pl.Series,
    *,
    event_column: str,
    source_id: str,
    visit_name_column: str | None = None,
    date_column: str | None = None,
    study_day_column: str | None = None,
    age_column: str | None = None,
) -> pl.DataFrame:
    """Build event rows from any table that has an event column."""
    out = pl.DataFrame(
        {
            "participant_id": participant_ids.cast(pl.Utf8),
            "event_id": as_text(df, event_column),
            "visit_name": _text(df, visit_name_column),
            "event_date": _text(df, date_column),
            "study_day": _num(df, study_day_column),
            "age_at_event": _num(df, age_column),
            "source_id": pl.Series([source_id] * df.height, dtype=pl.Utf8),
        }
    ).filter(pl.col("event_id").is_not_null() & (pl.col("event_id") != ""))
    # A visit can be listed on several rows (e.g. many assays); keep the first description.
    return out.unique(subset=["participant_id", "event_id"], keep="first", maintain_order=True)


def events_from_spec(
    df: pl.DataFrame, participant_ids: pl.Series, spec: TimelineSpec | None, source_id: str
) -> pl.DataFrame | None:
    if spec is None or not spec.event_column:
        return None
    return events_frame(
        df,
        participant_ids,
        event_column=spec.event_column,
        source_id=source_id,
        visit_name_column=spec.visit_name_column,
        date_column=spec.time_column,
        study_day_column=spec.study_day_column,
        age_column=spec.age_column,
    )


def register_events(reg: Registry, events: pl.DataFrame) -> int:
    """Insert events, never overwriting richer metadata registered earlier."""
    if events.height == 0:
        return 0
    events = EVENT_SCHEMA.validate(events)
    existing = reg.frame("SELECT participant_id, event_id FROM events")
    if existing.height:
        events = events.join(existing, on=["participant_id", "event_id"], how="anti")
    return reg.insert_frame("events", events, replace=False)


def load_events_file(project: Project, reg: Registry) -> int:
    spec = project.data.events
    if spec is None:
        return 0
    df = read_table(project.resolve(spec.path))
    res = AliasResolver(project).resolve("events", as_text(df, spec.participant_id_column))
    lookup = dict(res.mapping.iter_rows())
    pids = as_text(df, spec.participant_id_column).replace_strict(lookup, default=None)
    ev = events_frame(
        df, pids, event_column=spec.event_id_column, source_id="events",
        visit_name_column=spec.visit_name_column, date_column=spec.date_column,
        study_day_column=spec.study_day_column, age_column=spec.age_column,
    )
    # The events file is authoritative: replace generic rows derived from assay tables.
    reg.execute("DELETE FROM events WHERE source_id = 'events'")
    keys = ev.select("participant_id", "event_id")
    reg.con.register("_ev_keys", keys.to_arrow())
    reg.execute(
        "DELETE FROM events USING _ev_keys k "
        "WHERE events.participant_id = k.participant_id AND events.event_id = k.event_id"
    )
    reg.con.unregister("_ev_keys")
    return reg.insert_frame("events", EVENT_SCHEMA.validate(ev), replace=True)


def _event_order(reg: Registry) -> list[str]:
    """Order events by median study day / date when known, else by first appearance."""
    df = reg.frame(
        """
        SELECT event_id,
               median(study_day) AS day,
               min(TRY_CAST(event_date AS DATE)) AS first_date,
               min(rowid) AS first_seen
        FROM events GROUP BY event_id
        """
    )
    if df.height == 0:
        return []
    return df.sort(["day", "first_date", "first_seen"], nulls_last=True).get_column("event_id").to_list()


def timeline_summary(reg: Registry) -> dict[str, Any]:
    """Counts per event, modalities per event, attrition, repeated measures and time gaps."""
    order = _event_order(reg)
    per_event = reg.frame(
        "SELECT event_id, count(DISTINCT participant_id) AS participants FROM events GROUP BY event_id"
    )
    modalities = reg.frame(
        """
        SELECT event_id, source_id, modality, origin, count(DISTINCT participant_id) AS participants
        FROM assays WHERE event_id IS NOT NULL
        GROUP BY event_id, source_id, modality, origin
        ORDER BY event_id, modality
        """
    )
    specimens = reg.frame(
        "SELECT event_id, count(*) AS biospecimens FROM biospecimens WHERE event_id IS NOT NULL GROUP BY event_id"
    )
    repeated = reg.frame(
        """
        SELECT source_id, modality, count(*) AS participants_with_repeats
        FROM (SELECT source_id, modality, participant_id, count(DISTINCT event_id) AS n
              FROM assays WHERE event_id IS NOT NULL GROUP BY 1, 2, 3)
        WHERE n > 1 GROUP BY 1, 2
        """
    )
    # Attrition relative to the first event.
    attrition: list[dict[str, Any]] = []
    if order:
        base = set(
            reg.frame("SELECT DISTINCT participant_id FROM events WHERE event_id = ?", [order[0]])
            .get_column("participant_id").to_list()
        )
        for ev in order:
            present = set(
                reg.frame("SELECT DISTINCT participant_id FROM events WHERE event_id = ?", [ev])
                .get_column("participant_id").to_list()
            )
            retained = len(base & present)
            attrition.append({
                "event_id": ev,
                "participants": len(present),
                "retained_from_first": retained,
                "retention": (retained / len(base)) if base else None,
            })
    gaps = reg.frame(
        """
        WITH d AS (
            SELECT participant_id, event_id,
                   coalesce(study_day, date_diff('day', DATE '1970-01-01', TRY_CAST(event_date AS DATE))) AS t
            FROM events
        ), g AS (
            SELECT participant_id, t - lag(t) OVER (PARTITION BY participant_id ORDER BY t) AS gap
            FROM d WHERE t IS NOT NULL
        )
        SELECT count(gap) AS n_gaps, median(gap) AS median_gap_days, min(gap) AS min_gap_days,
               max(gap) AS max_gap_days
        FROM g WHERE gap IS NOT NULL
        """
    )
    return {
        "event_order": order,
        "participants_per_event": per_event.to_dicts(),
        "modalities_per_event": modalities.to_dicts(),
        "biospecimens_per_event": specimens.to_dicts(),
        "attrition": attrition,
        "repeated_measurements": repeated.to_dicts(),
        "time_gaps": gaps.to_dicts()[0] if gaps.height else {},
    }
