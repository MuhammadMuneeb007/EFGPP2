"""PHENOTYPE module: which columns of a table are phenotypes, and of which type.

Used by `efgpp phenotype add --path <file>` when --value-column / --type are not given.
Nothing here knows any trait: a phenotype is whatever column is not an ID or a covariate.

Rules (edit to suit your data):
  * ID columns (IID, participant_id, eid, ...) and family IDs (FID) are not phenotypes.
  * Columns that look like covariates (sex, age, PCs, batch, ...) are skipped when the file
    mixes phenotypes and covariates.
  * Type of each phenotype column:
      binary      exactly 2 distinct values (0/1, 1/2 PLINK, case/control, yes/no, true/false)
      continuous  numeric with more than MAX_CODE_LEVELS distinct values
      multiclass  text, or numeric codes with 3..MAX_CODE_LEVELS distinct values
    (ordinal cannot be inferred from values: pass --type ordinal --levels a,b,c)
  * Binary coding: PLINK 1/2 -> 2 = case; 0/1 -> 1 = case; text -> the value named like
    case/yes/true/affected is the case.
The original column name is kept everywhere (definition, standardized files).
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
    non_missing,
    normalize,
)
from efgpp.modules.covariates import covariate_role

CASE_WORDS = {"case", "cases", "yes", "y", "true", "t", "affected", "positive", "1"}
CONTROL_WORDS = {"control", "controls", "no", "n", "false", "f", "unaffected", "negative", "0"}


@dataclass
class PhenotypeColumn:
    column: str  # original column name, kept as is
    type: str  # binary | continuous | multiclass
    levels: list[str] | None = None
    case_values: list[str] | None = None
    control_values: list[str] | None = None
    n_observed: int = 0
    note: str = ""


@dataclass
class PhenotypeGuess:
    id_column: str | None
    phenotypes: list[PhenotypeColumn] = field(default_factory=list)
    skipped: dict[str, str] = field(default_factory=dict)  # column -> reason


def infer_type(values: pl.Series) -> PhenotypeColumn:
    """Type and coding of one phenotype column."""
    levels = distinct(values)
    n = non_missing(values).len()
    name = values.name
    if len(levels) == 2:
        low = {v.lower() for v in levels}
        if is_numeric(values):
            nums = sorted(float(v) for v in levels)
            if nums == [1.0, 2.0]:
                return PhenotypeColumn(name, "binary", None, [v for v in levels if float(v) == 2.0],
                                       [v for v in levels if float(v) == 1.0], n, "PLINK coding: 2 = case, 1 = control")
            if nums == [0.0, 1.0]:
                return PhenotypeColumn(name, "binary", None, [v for v in levels if float(v) == 1.0],
                                       [v for v in levels if float(v) == 0.0], n, "1 = case, 0 = control")
        case = [v for v in levels if v.lower() in CASE_WORDS]
        control = [v for v in levels if v.lower() in CONTROL_WORDS]
        if len(case) == 1 and len(control) == 1 and low:
            return PhenotypeColumn(name, "binary", None, case, control, n, f"{case[0]} = case")
        # Two arbitrary values: binary, the second in sort order is called the case.
        return PhenotypeColumn(name, "binary", None, [levels[1]], [levels[0]], n,
                               f"two values; {levels[1]!r} treated as case (override with --case-values)")
    if is_numeric(values) and len(levels) > MAX_CODE_LEVELS:
        return PhenotypeColumn(name, "continuous", None, None, None, n)
    if len(levels) >= 3:
        return PhenotypeColumn(name, "multiclass", levels, None, None, n)
    return PhenotypeColumn(name, "continuous" if is_numeric(values) else "multiclass", levels or None, None, None, n,
                           "fewer than two observed values")


def infer_phenotypes(df: pl.DataFrame, id_column: str | None = None,
                     value_columns: list[str] | None = None) -> PhenotypeGuess:
    """Find the ID column and every phenotype column of a table."""
    idc = find_id_column(df.columns, id_column)
    guess = PhenotypeGuess(idc)
    fids = set(family_id_columns(df.columns))
    candidates = value_columns or [c for c in df.columns if c != idc and c not in fids]
    for col in candidates:
        if col not in df.columns:
            guess.skipped[col] = "not in the file"
            continue
        role = covariate_role(col, df.get_column(col))
        if value_columns is None and role not in ("other",):
            guess.skipped[col] = f"looks like a covariate ({role})"
            continue
        info = infer_type(df.get_column(col))
        if info.n_observed == 0:
            guess.skipped[col] = "no values"
            continue
        guess.phenotypes.append(info)
    return guess


def default_name(column: str, file_stem: str, single: bool) -> str:
    """Phenotype name: the file name for single-phenotype files (<trait>.height -> <trait>),
    otherwise the column name."""
    generic = {"height", "pheno", "phenotype", "value", "trait", "status", "y", "outcome"}
    if single and normalize(column) in generic:
        return file_stem
    return column
