"""Standardized genetic prediction models (PredictDB, OmicsPred, MIMOSA, custom weights).

Every provider's raw files are converted once, at installation, into the same weight schema;
prediction never reads provider-specific formats directly (PredictDB SQLite is also read by
MetaXcan itself when that engine is chosen). Missing metadata stays null - nothing is invented.
"""

from __future__ import annotations

import re
import sqlite3
from pathlib import Path
from typing import Any

import polars as pl

WEIGHT_SCHEMA: dict[str, Any] = {
    "provider": pl.Utf8, "provider_dataset_id": pl.Utf8, "model_id": pl.Utf8, "feature_id": pl.Utf8,
    "feature_name": pl.Utf8, "modality": pl.Utf8, "chromosome": pl.Utf8, "position": pl.Int64,
    "variant_id": pl.Utf8, "ref": pl.Utf8, "alt": pl.Utf8, "effect_allele": pl.Utf8, "other_allele": pl.Utf8,
    "weight": pl.Float64, "genome_build": pl.Utf8, "tissue": pl.Utf8, "sample_type": pl.Utf8, "platform": pl.Utf8,
    "training_cohort": pl.Utf8, "training_ancestry": pl.Utf8, "training_n": pl.Int64,
    "validation_cohort": pl.Utf8, "validation_ancestry": pl.Utf8, "validation_r2": pl.Float64,
    "source_publication": pl.Utf8, "model_version": pl.Utf8,
}
# One row per model/feature (what is predicted and how well it validated).
FEATURE_SCHEMA: dict[str, Any] = {
    "model_id": pl.Utf8, "feature_id": pl.Utf8, "feature_name": pl.Utf8, "gene_id": pl.Utf8, "gene_name": pl.Utf8,
    "n_variants": pl.Int64, "validation_r2": pl.Float64, "validation_cohort": pl.Utf8,
    "below_threshold": pl.Boolean, "selected_method": pl.Utf8,
}

_VARID = re.compile(r"^(?:chr)?([0-9XYMT]+)_(\d+)_([ACGTN]+)_([ACGTN]+)(?:_b\d+)?$", re.IGNORECASE)


def standardize(df: pl.DataFrame, schema: dict[str, Any] | None = None, **constants: Any) -> pl.DataFrame:
    """Cast to the schema; add missing columns as null (or the given constant)."""
    schema = schema or WEIGHT_SCHEMA
    cols = []
    for name, dtype in schema.items():
        if name in df.columns:
            cols.append(pl.col(name).cast(dtype, strict=False))
        else:
            cols.append(pl.lit(constants.get(name), dtype=dtype).alias(name))
    out = df.select(cols)
    if "effect_allele" in out.columns:
        out = out.with_columns(pl.col("effect_allele").str.to_uppercase(), pl.col("other_allele").str.to_uppercase(),
                               pl.col("chromosome").str.replace(r"^(?i)chr", ""))
    return out


# ------------------------------------------------------------------ PredictDB (SQLite)
def predictdb_tissues(models_dir: Path, prefix: str = "mashr_", suffix: str = ".db") -> dict[str, Path]:
    """Tissues available in an extracted PredictDB archive (enumerated from the files, not hard-coded)."""
    out = {}
    for db in sorted(models_dir.rglob(f"*{suffix}")):
        name = db.name[: -len(suffix)]
        out[name[len(prefix):] if prefix and name.startswith(prefix) else name] = db
    return out


def _sqlite_tables(con: sqlite3.Connection) -> set[str]:
    return {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}


def predictdb_weights(db: Path, *, modality: str, tissue: str | None, dataset: str, genome_build: str | None,
                      provider: str = "predictdb") -> tuple[pl.DataFrame, pl.DataFrame]:
    """(standard weights, feature table) of one PredictDB model database."""
    con = sqlite3.connect(f"file:{db.as_posix()}?mode=ro", uri=True)
    try:
        tables = _sqlite_tables(con)
        if "weights" not in tables:
            raise ValueError(f"{db}: not a PredictDB model (no `weights` table)")
        cols = {r[1] for r in con.execute("PRAGMA table_info(weights)")}
        ref_col = next((c for c in ("ref_allele", "ref_vAllele") if c in cols), None)  # GTEx v8 MASHR: ref_allele
        if ref_col is None or not {"gene", "rsid", "varID", "eff_allele", "weight"} <= cols:
            raise ValueError(f"{db}: unexpected PredictDB weights columns {sorted(cols)}")
        w = pl.read_database(f"SELECT gene, rsid, varID, {ref_col} AS ref_vAllele, eff_allele, weight FROM weights",
                             con)
        extra = pl.read_database("SELECT * FROM extra", con) if "extra" in tables else pl.DataFrame({"gene": []})
    finally:
        con.close()
    parsed = [(_VARID.match(v or "") or None) for v in w.get_column("varID").to_list()]
    if genome_build is None:  # PredictDB variant ids state their build: chr1_123_A_G_b38
        suffixes = {m.group(1) for v in w.get_column("varID").drop_nulls().to_list()[:5000]
                    if (m := re.search(r"_b(3[78])$", v))}
        genome_build = {"38": "GRCh38", "37": "GRCh37"}[suffixes.pop()] if len(suffixes) == 1 else None
    w = w.with_columns(
        pl.Series("chromosome", [m.group(1) if m else None for m in parsed], dtype=pl.Utf8),
        pl.Series("position", [int(m.group(2)) if m else None for m in parsed], dtype=pl.Int64),
        pl.Series("ref", [m.group(3).upper() if m else None for m in parsed], dtype=pl.Utf8),
        pl.Series("alt", [m.group(4).upper() if m else None for m in parsed], dtype=pl.Utf8),
    )
    r2_col = next((c for c in ("pred.perf.R2", "pred_perf_r2", "R2") if c in extra.columns), None)
    name_col = next((c for c in ("genename", "gene_name") if c in extra.columns), None)
    n_col = next((c for c in ("n.snps.in.model", "n_snps_in_model") if c in extra.columns), None)
    feats = extra.select(
        pl.col("gene").cast(pl.Utf8).alias("feature_id"),
        (pl.col(name_col).cast(pl.Utf8) if name_col else pl.lit(None, dtype=pl.Utf8)).alias("feature_name"),
        (pl.col(n_col).cast(pl.Int64, strict=False) if n_col else pl.lit(None, dtype=pl.Int64)).alias("n_variants"),
        (pl.col(r2_col).cast(pl.Float64, strict=False) if r2_col else pl.lit(None, dtype=pl.Float64))
        .alias("validation_r2"),
    ) if extra.height else pl.DataFrame(schema={"feature_id": pl.Utf8, "feature_name": pl.Utf8,
                                                 "n_variants": pl.Int64, "validation_r2": pl.Float64})
    model_id = f"{provider}:{dataset}:{tissue}" if tissue else f"{provider}:{dataset}"
    weights = standardize(
        w.join(feats.select("feature_id", "feature_name", "validation_r2"), left_on="gene", right_on="feature_id",
               how="left")
        .select(pl.col("gene").alias("feature_id"), "feature_name", "chromosome", "position",
                pl.col("rsid").alias("variant_id"), "ref", "alt", pl.col("eff_allele").alias("effect_allele"),
                pl.col("ref_vAllele").alias("other_allele"), "weight", "validation_r2"),
        provider=provider, provider_dataset_id=dataset, model_id=model_id, modality=modality,
        genome_build=genome_build, tissue=tissue, model_version=db.name,
    )
    features = standardize(feats.with_columns(pl.lit(model_id).alias("model_id")), FEATURE_SCHEMA)
    return weights, features


# ------------------------------------------------------------------ PGS-style scoring files (OmicsPred)
def read_scoring_file(path: Path) -> tuple[dict[str, str], pl.DataFrame]:
    """PGS-Catalog-style scoring file (`#key=value` header, then a table) -> (header, weights).

    Harmonized files (hm_chr/hm_pos, `#HmPOS_build=`) are used in the harmonized build."""
    import gzip

    header: dict[str, str] = {}
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        skip = 0
        for line in fh:
            if not line.startswith("#"):
                break
            skip += 1
            if "=" in line:
                k, v = line.lstrip("#").strip().split("=", 1)
                header[k.strip()] = v.strip()
    df = pl.read_csv(path, separator="\t", skip_rows=skip, infer_schema=False, quote_char=None)
    cols = {c.lower(): c for c in df.columns}
    hm_build = header.get("HmPOS_build")
    use_hm = "hm_pos" in cols and hm_build is not None
    build = hm_build if use_hm else header.get("genome_build")
    chrom = cols.get("hm_chr") if use_hm else cols.get("chr_name")
    pos = cols.get("hm_pos") if use_hm else cols.get("chr_position")
    rsid = cols.get("hm_rsid") or cols.get("rsid")
    out = df.select(
        (pl.col(chrom) if chrom else pl.lit(None, dtype=pl.Utf8)).alias("chromosome"),
        (pl.col(pos).cast(pl.Int64, strict=False) if pos else pl.lit(None, dtype=pl.Int64)).alias("position"),
        (pl.col(rsid) if rsid else pl.lit(None, dtype=pl.Utf8)).alias("variant_id"),
        pl.col(cols["effect_allele"]).alias("effect_allele"),
        (pl.col(cols["other_allele"]) if "other_allele" in cols else pl.lit(None, dtype=pl.Utf8)).alias("other_allele"),
        pl.col(cols["effect_weight"]).cast(pl.Float64, strict=False).alias("weight"),
    )
    header["_genome_build"] = build or "unknown"
    return header, out
