"""MIMOSA whole-blood DNA-methylation prediction models -> standardized EFGPP weights.

MIMOSA (Melton et al. 2023; Zenodo 8400313, MIMOSA-Models.zip) provides, for each CpG, five
penalized-regression models (ElNet, MNet, SCAD, MCP, LASSO) trained with GoDMC mQTL summary
data and evaluated on Framingham Heart Study test data. `mimosa_export.R` dumps every method;
here the model with the highest test R2 among valid ones is selected (ties: ElNet, MNet, SCAD,
MCP, LASSO). Models at or below `minimum_model_r2` (upstream MIMOSA uses 0.005) are kept but
marked below_threshold and are not scored by default.

The prediction is a genetically predicted methylation SCORE (sum of dosage x weight), not a
beta value: MIMOSA supplies no intercept/calibration for reconstructing beta values.
Context: whole blood only - never brain, liver or tumour methylation.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from efgpp.data.predicted.models import FEATURE_SCHEMA, standardize

METHODS = ["ElNet", "MNet", "SCAD", "MCP", "LASSO"]
EXPORT_SCRIPT = Path(__file__).with_name("mimosa_export.R")
CONTEXT = {
    "modality": "methylation", "tissue": "whole_blood", "sample_type": "whole blood",
    "training_resource": "GoDMC mQTL summary statistics; Framingham Heart Study test/validation data",
    "value": "genetically_predicted_methylation_score (not a beta value)",
    "not_for": "brain, liver or tumour methylation",
}


def read_export(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t", infer_schema=False, null_values=["NA", ""]).with_columns(
        pl.col("method_index").cast(pl.Int64), pl.col("valid").str.to_uppercase() == "TRUE",
        pl.col("test_r2").cast(pl.Float64, strict=False), pl.col("lambda").cast(pl.Float64, strict=False),
        pl.col("position").cast(pl.Float64, strict=False).cast(pl.Int64, strict=False),
        pl.col("weight").cast(pl.Float64, strict=False))


def select_models(export: pl.DataFrame, minimum_r2: float = 0.005) -> tuple[pl.DataFrame, pl.DataFrame]:
    """(weights of the selected method per CpG, one model row per CpG)."""
    order = {m: i for i, m in enumerate(METHODS)}
    per_method = export.group_by("cpg_id", "method").agg(
        pl.col("valid").first(), pl.col("test_r2").first(), pl.col("weight").drop_nulls().len().alias("n_variants"))
    candidates = per_method.filter(pl.col("valid") & pl.col("test_r2").is_not_null() & (pl.col("n_variants") > 0))
    candidates = candidates.with_columns(
        pl.col("method").replace_strict(order, default=len(order), return_dtype=pl.Int64).alias("_order"))
    best = candidates.sort(["cpg_id", "test_r2", "_order"], descending=[False, True, False]).unique(
        "cpg_id", keep="first", maintain_order=True)
    all_cpgs = export.select("cpg_id").unique()
    models = all_cpgs.join(best.select("cpg_id", pl.col("method").alias("selected_method"),
                                       pl.col("test_r2").alias("validation_r2"), "n_variants"),
                           on="cpg_id", how="left").with_columns(
        (pl.col("validation_r2") <= minimum_r2).fill_null(True).alias("below_threshold"),
        pl.lit("GoDMC mQTL / Framingham Heart Study test").alias("training_resource"),
        pl.lit("whole_blood").alias("tissue"),
    ).sort("cpg_id")
    chosen = export.join(best.select("cpg_id", "method"), on=["cpg_id", "method"], how="semi").filter(
        pl.col("weight").is_not_null())
    weights = chosen.select(
        pl.col("cpg_id").alias("feature_id"), pl.col("cpg_id").alias("feature_name"),
        pl.col("chromosome").str.replace(r"^(?i)chr", ""), "position", pl.col("snp").alias("variant_id"),
        pl.col("a1").alias("effect_allele"), pl.col("a2").alias("other_allele"), "weight",
        pl.col("test_r2").alias("validation_r2"), pl.col("method").alias("model_method"))
    return weights, models


def standardized(export: pl.DataFrame, minimum_r2: float, *, genome_build: str, version: str
                 ) -> tuple[pl.DataFrame, pl.DataFrame]:
    weights, models = select_models(export, minimum_r2)
    model_id = f"mimosa:{version}"
    w = standardize(weights, provider="mimosa", provider_dataset_id=version, model_id=model_id,
                    modality="methylation", genome_build=genome_build, tissue="whole_blood",
                    sample_type="whole blood", training_cohort="GoDMC (mQTL summary statistics)",
                    validation_cohort="Framingham Heart Study", model_version=version,
                    source_publication="Melton et al. 2023, MIMOSA (doi:10.5281/zenodo.8400313)")
    feats = standardize(models.with_columns(pl.lit(model_id).alias("model_id"),
                                            pl.col("cpg_id").alias("feature_id"),
                                            pl.col("cpg_id").alias("feature_name"),
                                            pl.lit("Framingham Heart Study").alias("validation_cohort")),
                        FEATURE_SCHEMA)
    return w.with_columns(weights.get_column("model_method").alias("model_method")) if weights.height else w, feats
