"""AlphaGenome regulatory annotation (optional).

Two modes (resources.yaml: alphagenome.mode):
  atlas  precomputed variant scores (a TSV/Parquet keyed by chromosome/position/ref/alt);
         preferred for large variant sets, no API key needed.
  api    live model queries for a targeted variant set (at most `max_api_variants`),
         requires the `alphagenome` package and an API key in the configured environment
         variable. The key is never required for a basic EFGPP installation.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import polars as pl

from efgpp.data.annotation import KEY, cohort_variants, register_annotation
from efgpp.project import Project


def api_key(project: Project) -> str | None:
    return os.environ.get(project.resources.alphagenome.api_key_env)


def _atlas(project: Project, variants: pl.DataFrame, path: Path) -> pl.DataFrame:
    atlas = pl.read_parquet(path) if path.suffix == ".parquet" else pl.read_csv(
        path, separator="\t", infer_schema_length=10000, schema_overrides={"chromosome": pl.Utf8})
    atlas = atlas.with_columns(pl.col("chromosome").cast(pl.Utf8).str.replace(r"^chr", ""),
                               pl.col("position").cast(pl.Int64))
    return variants.select(KEY).join(atlas, on=KEY, how="inner")


def _api(project: Project, variants: pl.DataFrame, build: str | None) -> pl.DataFrame:  # pragma: no cover - live API
    key = api_key(project)
    if not key:
        raise RuntimeError(f"AlphaGenome API mode needs ${project.resources.alphagenome.api_key_env}")
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client, variant_scorers
    except ImportError as exc:
        raise RuntimeError("install the optional client: pip install 'efgpp[alphagenome]'") from exc
    limit = project.resources.alphagenome.max_api_variants
    if variants.height > limit:
        raise RuntimeError(f"{variants.height} variants exceed max_api_variants={limit}; use atlas mode")
    model = dna_client.create(key)
    scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
    rows: list[dict[str, Any]] = []
    width = dna_client.SEQUENCE_LENGTH_1MB
    for chrom, pos, ref, alt in variants.select(KEY).iter_rows():
        v = genome.Variant(chromosome=f"chr{chrom}", position=pos, reference_bases=ref, alternate_bases=alt)
        interval = v.reference_interval.resize(width)
        scores = model.score_variant(interval=interval, variant=v, variant_scorers=scorers)
        tidy = variant_scorers.tidy_scores([scores])
        for rec in tidy.to_dict("records"):
            rows.append({"chromosome": chrom, "position": pos, "reference": ref, "alternate": alt,
                         **{k: (str(v) if not isinstance(v, int | float) else v) for k, v in rec.items()}})
    return pl.DataFrame(rows)


def run_alphagenome(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    cfg = project.resources.alphagenome
    variants, vt, build = cohort_variants(project, source_id)
    if cfg.mode == "atlas":
        if not cfg.atlas_path:
            raise RuntimeError("AlphaGenome atlas mode needs resources.yaml: alphagenome.atlas_path")
        table = _atlas(project, variants, project.resolve(cfg.atlas_path))
    else:
        table = _api(project, variants, build)
    return register_annotation(
        project, source_id=source_id, name="alphagenome", table=table, parent=vt, tool=f"alphagenome_{cfg.mode}",
        tool_version=None, resource_versions={"alphagenome": cfg.version or cfg.mode}, genome_build=build,
        metadata={"mode": cfg.mode, "matched_variants": table.height},
    )
