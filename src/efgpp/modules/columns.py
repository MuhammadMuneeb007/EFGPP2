"""Shared column rules used by every data module: ID columns and value types.

Edit this file (or `efgpp modules export columns` and edit the project copy) to teach EFGPP
the column conventions of your cohort. Column names are only *read*; EFGPP never renames them.
"""

from __future__ import annotations

import re

import polars as pl

# ----------------------------------------------------------------------------------------------
# Participant ID columns, most specific first. Matching ignores case and a leading '#'.
# ----------------------------------------------------------------------------------------------
ID_COLUMNS = ["participant_id", "iid", "eid", "sample_id", "sampleid", "sample", "subject_id", "subject",
              "individual_id", "id", "patient_id"]
FAMILY_ID_COLUMNS = ["fid", "family_id", "famid"]

# Values that mean "missing" in user tables (PLINK uses -9).
MISSING_TOKENS = ["", "NA", "N/A", "NaN", "nan", "NULL", "null", "None", ".", "-9"]

# A numeric column with at most this many distinct values may be categorical (e.g. codes 1/2/3).
MAX_CODE_LEVELS = 10


def normalize(name: str) -> str:
    return name.strip().lstrip("#").lower()


def find_id_column(columns: list[str], preferred: str | None = None) -> str | None:
    """The participant ID column: `preferred` if present, else the first known ID name."""
    if preferred and preferred in columns:
        return preferred
    lookup = {normalize(c): c for c in columns}
    for name in ID_COLUMNS:
        if name in lookup:
            return lookup[name]
    return None


def family_id_columns(columns: list[str]) -> list[str]:
    return [c for c in columns if normalize(c) in FAMILY_ID_COLUMNS]


def non_missing(values: pl.Series) -> pl.Series:
    s = values.cast(pl.Utf8).str.strip_chars()
    return s.filter(s.is_not_null() & ~s.is_in(MISSING_TOKENS))


def is_numeric(values: pl.Series) -> bool:
    s = non_missing(values)
    return s.len() > 0 and s.cast(pl.Float64, strict=False).null_count() == 0


def distinct(values: pl.Series) -> list[str]:
    return sorted(non_missing(values).unique().to_list())


def looks_like(name: str, patterns: list[str]) -> bool:
    """True if the column name matches any regex in `patterns` (case-insensitive)."""
    n = normalize(name)
    return any(re.fullmatch(p, n) for p in patterns)
