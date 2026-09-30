"""Reading user-supplied tables without guessing types behind the user's back.

Identifier columns must never lose leading zeros, so text formats are read as strings
and typed explicitly afterwards.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import polars as pl

TEXT_SUFFIXES = (".csv", ".tsv", ".txt", ".tab", ".cov", ".pheno", ".phen", ".height")


def detect_table_format(path: Path, declared: str = "auto") -> str:
    if declared not in ("auto", None, ""):
        return declared
    name = path.name.lower()
    for ext in (".gz", ".bgz", ".zst"):
        if name.endswith(ext):
            name = name[: -len(ext)]
    if path.is_dir() and name.endswith(".zarr"):
        return "zarr"
    if name.endswith(".parquet") or name.endswith(".pq"):
        return "parquet"
    if name.endswith(".h5ad"):
        return "h5ad"
    if name.endswith(".zarr"):
        return "zarr"
    if name.endswith(".csv"):
        return "csv"
    if name.endswith(".tsv") or name.endswith(".tab"):
        return "tsv"
    return "txt"  # sniffed delimiter


def _open_text(path: Path):  # type: ignore[no-untyped-def]
    if path.name.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def sniff_separator(path: Path) -> str | None:
    """Return ',' or '\\t', or None for runs of whitespace (PLINK-style tables)."""
    with _open_text(path) as fh:
        header = fh.readline()
    if "\t" in header:
        return "\t"
    if "," in header:
        return ","
    return None


def read_table(path: Path, fmt: str = "auto", columns: list[str] | None = None) -> pl.DataFrame:
    """Read a table; text formats come back with every column as Utf8."""
    fmt = detect_table_format(path, fmt)
    if fmt == "parquet":
        return pl.read_parquet(path, columns=columns)
    if fmt in ("csv", "tsv", "txt"):
        sep = {"csv": ",", "tsv": "\t"}.get(fmt) or sniff_separator(path)
        if sep is None:
            with _open_text(path) as fh:
                rows = [re.split(r"\s+", ln.strip()) for ln in fh if ln.strip()]
            header, body = rows[0], rows[1:]
            df = pl.DataFrame(body, schema=[(h, pl.Utf8) for h in header], orient="row")
        else:
            df = pl.read_csv(path, separator=sep, infer_schema=False)
        return df.select(columns) if columns else df
    raise ValueError(f"{path}: {fmt!r} is not a tabular format")


def table_columns(path: Path, fmt: str = "auto") -> list[str]:
    fmt = detect_table_format(path, fmt)
    if fmt == "parquet":
        return list(pl.read_parquet_schema(path).keys())
    return read_table(path, fmt).columns


def as_text(df: pl.DataFrame, column: str) -> pl.Series:
    """Column as stripped text, independent of the source dtype."""
    s = df.get_column(column)
    if s.dtype != pl.Utf8:
        s = s.cast(pl.Utf8)
    return s.str.strip_chars()


def missing_mask(values: pl.Series, tokens: list[str]) -> pl.Series:
    return values.is_null() | values.is_in(tokens)


def write_parquet(df: pl.DataFrame, path: Path, compression: str = "zstd") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.write_parquet(tmp, compression=compression)  # type: ignore[arg-type]
    tmp.replace(path)
    return path
