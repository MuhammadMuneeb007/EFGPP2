"""Relatedness from KING-robust kinship coefficients."""

from __future__ import annotations

import polars as pl

# KING kinship boundaries (Manichaikul et al. 2010).
DEGREES = [
    (0.354, "duplicate_or_mz_twin"),
    (0.177, "first_degree"),
    (0.0884, "second_degree"),
    (0.0442, "third_degree"),
]


def classify_pairs(kin0: pl.DataFrame) -> pl.DataFrame:
    """Normalize a PLINK2 .kin0 table and label the relationship degree of each pair."""
    if kin0.height == 0:
        return pl.DataFrame(schema={"IID1": pl.Utf8, "IID2": pl.Utf8, "KINSHIP": pl.Float64, "relationship": pl.Utf8})
    expr = pl.lit("unrelated")
    for threshold, label in reversed(DEGREES):
        expr = pl.when(pl.col("KINSHIP") > threshold).then(pl.lit(label)).otherwise(expr)
    cols = [c for c in ("FID1", "IID1", "FID2", "IID2", "NSNP", "HETHET", "IBS0", "KINSHIP") if c in kin0.columns]
    return kin0.select(cols).with_columns(expr.alias("relationship"))


def greedy_unrelated_removal(pairs: pl.DataFrame, missingness: dict[str, float | None]) -> set[str]:
    """Smallest-effort removal set leaving no pair: repeatedly drop the sample involved in
    the most remaining pairs, breaking ties by higher missingness."""
    edges = {(a, b) for a, b in pairs.select("IID1", "IID2").iter_rows()}
    removed: set[str] = set()
    while edges:
        degree: dict[str, int] = {}
        for a, b in edges:
            degree[a] = degree.get(a, 0) + 1
            degree[b] = degree.get(b, 0) + 1
        victim = max(degree, key=lambda s: (degree[s], missingness.get(s) or 0.0, s))
        removed.add(victim)
        edges = {(a, b) for a, b in edges if victim not in (a, b)}
    return removed
