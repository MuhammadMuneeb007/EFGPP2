"""COVARIATES module: which columns are covariates, what role each plays, categorical or not.

Used by `efgpp data add covariates --path <file>` (all columns except IDs are taken unless
--columns is given). Roles are inferred from column names and values:

    sex        sex, gender, is_male, msex, ...     (values M/F, male/female, 1/2, 0/1)
    age        age, age_at_*, *_age, yob, birth_year
    pc         PC1, pc_2, genetic_pc3, ...
    batch      batch, array, chip, centre/center, site, plate, assessment_centre
    other      anything else

Categorical vs numeric:
    text values                                   -> categorical
    sex / batch roles                              -> categorical (codes stay as they are)
    numeric with <= MAX_CODE_LEVELS integer codes  -> categorical
    other numeric                                  -> numeric
Values and column names are never changed; no imputation, no scaling (that would leak).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import polars as pl

from efgpp.modules.columns import (
    MAX_CODE_LEVELS,
    distinct,
    family_id_columns,
    find_id_column,
    is_numeric,
    looks_like,
)

ROLE_PATTERNS: dict[str, list[str]] = {
    "sex": [r"sex", r"gender", r"is_?male", r"is_?female", r"msex", r"sex_?.*", r".*_sex", r"genetic_?sex"],
    "age": [r"age", r"age_.*", r".*_age", r"yob", r"birth_?year", r"year_?of_?birth", r"age\d*"],
    "pc": [r"pc_?\d+", r"genetic_?pc_?\d+", r"gpc_?\d+", r"pca_?\d+", r"c\d+"],
    "batch": [r"batch", r"array", r"chip", r"centre", r"center", r"site", r"plate", r"assessment_?cent(re|er)",
              r"genotyping_?batch", r"cohort", r"study"],
}
SEX_VALUES = [{"m", "f"}, {"male", "female"}, {"1", "2"}, {"0", "1"}]


def covariate_role(name: str, values: pl.Series | None = None) -> str:
    for role, patterns in ROLE_PATTERNS.items():
        if looks_like(name, patterns):
            if role == "sex" and values is not None:
                vals = {v.lower() for v in distinct(values)}
                if vals and not any(vals <= allowed for allowed in SEX_VALUES):
                    continue  # named like sex but values are not sex codes
            if role == "pc" and values is not None and not is_numeric(values):
                continue
            return role
    return "other"


def is_categorical(name: str, values: pl.Series, role: str) -> bool:
    if role in ("sex", "batch"):
        return True
    if not is_numeric(values):
        return True
    levels = distinct(values)
    integer_codes = all(float(v).is_integer() for v in levels)
    return role not in ("age", "pc") and integer_codes and len(levels) <= MAX_CODE_LEVELS


@dataclass
class CovariateGuess:
    id_column: str | None
    variables: list[str] = field(default_factory=list)  # original names, file order
    categorical: list[str] = field(default_factory=list)
    roles: dict[str, str] = field(default_factory=dict)
    skipped: dict[str, str] = field(default_factory=dict)


def infer_covariates(df: pl.DataFrame, id_column: str | None = None,
                     columns: list[str] | None = None) -> CovariateGuess:
    idc = find_id_column(df.columns, id_column)
    guess = CovariateGuess(idc)
    fids = set(family_id_columns(df.columns))
    for col in columns or df.columns:
        if col == idc or col in fids:
            continue
        if col not in df.columns:
            guess.skipped[col] = "not in the file"
            continue
        values = df.get_column(col)
        if not distinct(values):
            guess.skipped[col] = "no values"
            continue
        role = covariate_role(col, values)
        guess.variables.append(col)
        guess.roles[col] = role
        if is_categorical(col, values, role):
            guess.categorical.append(col)
    return guess
