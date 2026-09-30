"""Loading genotype metadata and parsing PLINK report files into typed tables.

PLINK outputs are parsed into Parquet; nothing downstream reads PLINK logs.
"""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl

from efgpp.data.genotype.formats import GenotypeFileset, read_samples, resolve_fileset

# Column dtypes for common PLINK 2 / 1.9 report columns.
_FLOAT_COLS = {
    "F_MISS", "ALT_FREQS", "MAF", "P", "F", "O(HOM)", "E(HOM)", "O(HET_A1)", "E(HET_A1)",
    "KINSHIP", "HETHET", "IBS0", "SNPSEX", "PEDSEX", "KB", "DENSITY", "PHOM", "PHET",
    "SNPSEX_F",
}
_INT_COLS = {"MISSING_CT", "OBS_CT", "NSNP", "POS", "BP", "HOM_A1_CT", "HET_A1_CT", "TWO_AX_CT",
             "NSNP_HET", "POS1", "POS2", "NSEG"}


def read_plink_table(path: Path) -> pl.DataFrame:
    """Parse a PLINK report (tab- or space-delimited; leading '#' on the header)."""
    with open(path, encoding="utf-8") as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.startswith("##")]
    if not lines:
        return pl.DataFrame()
    header = re.split(r"\s+", lines[0].strip().lstrip("#"))
    rows = [re.split(r"\s+", ln.strip()) for ln in lines[1:]]
    df = pl.DataFrame(rows, schema=[(h, pl.Utf8) for h in header], orient="row")
    casts = []
    for c in df.columns:
        if c in _FLOAT_COLS or c.startswith("PC") and c[2:].isdigit():
            casts.append(pl.col(c).replace({"NA": None, "nan": None}).cast(pl.Float64, strict=False))
        elif c in _INT_COLS:
            casts.append(pl.col(c).cast(pl.Int64, strict=False))
    return df.with_columns(casts) if casts else df


def read_eigenvec(path: Path) -> pl.DataFrame:
    return read_plink_table(path)


def read_eigenval(path: Path) -> list[float]:
    return [float(x) for x in path.read_text(encoding="utf-8").split()]


def load_fileset(path: str | Path, fmt: str = "auto") -> GenotypeFileset:
    return resolve_fileset(Path(path), fmt)


def sample_frame(fs: GenotypeFileset) -> pl.DataFrame:
    return read_samples(fs)


def write_keep_file(path: Path, samples: pl.DataFrame) -> Path:
    """Write a PLINK2 --keep/--remove file from a frame with FID (nullable) and IID."""
    path.parent.mkdir(parents=True, exist_ok=True)
    has_fid = samples.get_column("FID").null_count() < samples.height
    with open(path, "w", encoding="utf-8") as fh:
        if has_fid:
            fh.write("#FID\tIID\n")
            for fid, iid in samples.select("FID", "IID").iter_rows():
                fh.write(f"{fid}\t{iid}\n")
        else:
            fh.write("#IID\n")
            for (iid,) in samples.select("IID").iter_rows():
                fh.write(f"{iid}\n")
    return path


def write_id_list(path: Path, ids: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{i}\n" for i in ids), encoding="utf-8")
    return path
