"""Genetically predicted molecular data: planning and running one prediction unit.

A unit is one tissue (PredictDB expression/splicing) or one dataset (OmicsPred / PredictDB
protein, MIMOSA) of one modality; each unit is one restartable plan step
(`predict.<modality>.<unit>`). Engines:

  metaxcan                 PrediXcan (MetaXcan Predict.py, individual-level) on PredictDB SQLite
  plink_score / generic_weights / mimosa
                           standardized weights -> harmonize -> `plink2 --score`

Every unit writes, under data/predicted/<modality>/<provider>/<unit>/:
  predicted_<modality>.parquet   participant x feature (wide) or partitioned long table
  model_qc.parquet               per feature: variant coverage, allele matching, status
  feature_metadata.parquet       what each feature is and how it validated
  manifest.yaml                  model resource, versions, checksums, configuration, limitations

Values are genetically predicted components - never observed data - and no phenotype is read.
Nothing runs without its model resource: a missing model is reported, never replaced.
"""

from __future__ import annotations

import json
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import polars as pl
import yaml

from efgpp.config.data import PredictedModalityConfig
from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.io import write_parquet
from efgpp.data.registry import Registry, utcnow
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths

VALUE_NAMES = {
    Modality.EXPRESSION: "genetically_predicted_expression",
    Modality.SPLICING: "genetically_predicted_splicing",
    Modality.PROTEOMICS: "genetically_predicted_protein_score",
    Modality.METABOLOMICS: "genetically_predicted_metabolite_score",
    Modality.METHYLATION: "genetically_predicted_methylation_score",
}
LIMITATIONS = {
    Modality.EXPRESSION: "PrediXcan predicts the genetically regulated component of expression; it is not RNA-seq.",
    Modality.SPLICING: "sQTL models give genetically predicted splicing traits, not splice-junction measurements.",
    Modality.PROTEOMICS: "Protein genetic scores are the genetically predictable component, not measured protein "
                         "concentration.",
    Modality.METABOLOMICS: "Metabolite genetic scores are the genetically predictable component, not measured "
                           "metabolomics.",
    Modality.METHYLATION: "MIMOSA estimates the genetically predictable whole-blood methylation component, not a "
                          "measured beta value; not valid for other tissues.",
}
GTEX_RESOURCES = {Modality.EXPRESSION: "predictdb-gtex-v8-expression", Modality.SPLICING: "predictdb-gtex-v8-splicing"}
APPROX_SIZE = {"predictdb-gtex-v8-expression": "~262 MB", "predictdb-gtex-v8-splicing": "~669 MB",
               "mimosa": "~3.5 GB archive"}
WIDE_LIMIT = 20_000_000  # participant x feature cells above which the long, partitioned layout is used
ANCESTRY_CODES = {"european": "EUR", "eur": "EUR", "nfe": "EUR", "white": "EUR", "african": "AFR", "afr": "AFR",
                  "african american or afro-caribbean": "AFR", "afa": "AFR", "east asian": "EAS", "eas": "EAS",
                  "chn": "EAS", "south asian": "SAS", "sas": "SAS", "hispanic or latin american": "AMR",
                  "his": "AMR", "amr": "AMR", "admixed american": "AMR"}


@dataclass
class PredictionUnit:
    modality: Modality
    key: str  # tissue or dataset
    engine: str
    provider: str
    dataset: str
    tissue: str | None = None
    model_path: Path | None = None  # PredictDB .db, or standardized weights (file or directory)
    resource: dict[str, Any] | None = None
    reason: str | None = None
    required: str | None = None  # command that would make this unit available

    @property
    def step_id(self) -> str:
        return f"predict.{self.modality.value}.{self.key}"

    @property
    def output_subdir(self) -> str:
        """gtex_v8/<tissue> | <provider>/<dataset> | mimosa/whole_blood (under data/predicted/<modality>/)."""
        if self.provider == "predictdb" and self.dataset.startswith("gtex_v8"):
            return f"gtex_v8/{self.key}"
        if self.provider == "mimosa":
            return f"mimosa/{self.tissue or 'whole_blood'}"
        return f"{self.provider}/{self.key}"


def installed_resource(project: Project, name: str, version: str | None = None) -> dict[str, Any] | None:
    q = "SELECT * FROM resources WHERE name = ?" + (" AND version = ?" if version else "")
    with Registry.open(project) as reg:
        rows = reg.rows(q + " ORDER BY download_date DESC LIMIT 1", [name, version] if version else [name])
    if not rows:
        return None
    row = rows[0]
    if isinstance(row.get("metadata"), str):
        row["metadata"] = json.loads(row["metadata"] or "{}")
    path = row.get("local_path")
    return row if path and project.resolve(path).exists() else None


def plan_units(project: Project, modality: Modality, cfg: PredictedModalityConfig) -> list[PredictionUnit]:
    from efgpp.data.predicted.models import predictdb_tissues

    engine = cfg.engine
    if modality in GTEX_RESOURCES and cfg.provider == "predictdb":
        name = GTEX_RESOURCES[modality]
        res = installed_resource(project, name)
        mp = cfg.model_provider
        models_dir = project.resolve(mp.models_dir) if mp.models_dir else (
            project.resolve(res["local_path"]) if res else None)
        if models_dir is None and modality == Modality.EXPRESSION:  # older `resources install predictdb`
            res = installed_resource(project, "predictdb")
            models_dir = project.resolve(res["local_path"]) if res else None
        available = predictdb_tissues(models_dir, mp.model_prefix, mp.model_suffix) if models_dir else {}
        wanted = sorted(available) if "all" in cfg.tissues else list(cfg.tissues)
        dataset = cfg.dataset or ("gtex_v8_mashr_eqtl" if modality == Modality.EXPRESSION else "gtex_v8_mashr_sqtl")
        if not wanted:
            return [PredictionUnit(modality, "none", engine, cfg.provider, dataset,
                                   reason=f"no tissues selected (predicted.{modality.value}.tissues, or 'all')",
                                   required=f"efgpp data predict enable {modality.value} --tissue Whole_Blood")]
        units = []
        for t in wanted:
            path = project.resolve(mp.model_paths[t]) if t in mp.model_paths else available.get(t)
            u = PredictionUnit(modality, t, engine, cfg.provider, dataset, tissue=t, model_path=path, resource=res)
            if path is None or not path.exists():
                if not available:
                    u.reason = f"GTEx v8 {modality.value} models not installed"
                    u.required = f"efgpp resources install {name}  ({APPROX_SIZE[name]})"
                else:
                    u.reason = f"tissue {t!r} not in the installed archive ({len(available)} tissues)"
            units.append(u)
        return units
    if modality in (Modality.PROTEOMICS, Modality.METABOLOMICS, Modality.EXPRESSION, Modality.SPLICING):
        if not cfg.datasets:
            return [PredictionUnit(modality, "none", engine, cfg.provider, "", reason="no dataset selected",
                                   required=f"efgpp resources omicspred list --modality {modality.value}")]
        units = []
        for d in cfg.datasets:
            name = "omicspred" if cfg.provider == "omicspred" else "predictdb-protein"
            res = installed_resource(project, name, d if name == "omicspred" else d.removesuffix(".db"))
            path = project.resolve(res["local_path"]) if res else None
            u = PredictionUnit(modality, d.removesuffix(".db"), engine, cfg.provider, d, model_path=path, resource=res)
            if res is None:
                u.reason = f"{cfg.provider} dataset {d} not installed"
                u.required = f"efgpp resources install {name} --dataset {d}"
            units.append(u)
        return units
    if modality == Modality.METHYLATION:
        res = installed_resource(project, "mimosa")
        u = PredictionUnit(modality, cfg.dataset or "whole_blood_v2", engine, "mimosa", cfg.dataset or "whole_blood_v2",
                           tissue="whole_blood", model_path=project.resolve(res["local_path"]) if res else None,
                           resource=res)
        if res is None:
            u.reason, u.required = "MIMOSA models not installed", f"efgpp resources install mimosa  ({APPROX_SIZE['mimosa']})"
        return [u]
    return [PredictionUnit(modality, "none", engine, cfg.provider, cfg.dataset or "",
                           reason=f"provider {cfg.provider} is not supported for {modality.value}")]


def find_unit(project: Project, modality: Modality, key: str) -> tuple[PredictedModalityConfig, PredictionUnit]:
    cfg: PredictedModalityConfig = getattr(project.data.predicted, modality.value)
    for u in plan_units(project, modality, cfg):
        if u.key == key:
            return cfg, u
    raise KeyError(f"no {modality.value} prediction unit {key!r}")


def load_weights(unit: PredictionUnit, declared_build: str | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Standardized (weights, features) of a unit's model resource. `declared_build` (the
    configured model_genome_build) is used only when neither the resource nor the model says."""
    from efgpp.data.predicted.models import FEATURE_SCHEMA, predictdb_weights, standardize

    path = unit.model_path
    if path is None:
        raise RuntimeError(f"{unit.step_id}: {unit.reason}")
    if path.suffix == ".db":
        build = (unit.resource or {}).get("genome_build") or ((unit.resource or {}).get("metadata") or {}).get(
            "genome_build")
        w, f = predictdb_weights(path, modality=unit.modality.value, tissue=unit.tissue, dataset=unit.dataset,
                                 genome_build=build, provider="predictdb")
        if declared_build and w.get_column("genome_build").null_count() == w.height:
            w = w.with_columns(pl.lit(declared_build).alias("genome_build"))
        return w, f
    wdir = path / "weights" if (path / "weights").is_dir() else path
    files = sorted(wdir.rglob("*.parquet")) if wdir.is_dir() else [wdir]
    weights = pl.concat([pl.read_parquet(f) for f in files if f.name != "features.parquet"
                         and f.name != "models.parquet"], how="vertical_relaxed")
    for name in ("features.parquet", "models.parquet"):
        if (path / name).exists():
            return weights, pl.read_parquet(path / name)
    feats = weights.group_by("model_id", "feature_id").agg(pl.col("feature_name").first(), pl.len().alias("n_variants"),
                                                           pl.col("validation_r2").first())
    return weights, standardize(feats, FEATURE_SCHEMA)


def participant_ancestries(project: Project, source_id: str) -> set[str]:
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).latest(source_id=source_id, artifact_type="ancestry")
    if art is None or not Path(art.path).exists():
        return set()
    df = pl.read_parquet(art.path)
    if "ancestry" not in df.columns:
        return set()
    return {ANCESTRY_CODES.get(str(a).lower(), str(a).upper()) for a in df.get_column("ancestry").drop_nulls().unique()}


def ancestry_match_status(model_ancestry: str | None, participants: set[str]) -> str:
    """matched | mixed | partial | mismatched | unknown (reported, never used to reject a model)."""
    if not model_ancestry or not participants:
        return "unknown"
    words = [w.strip() for w in model_ancestry.replace(";", ",").split(",") if w.strip()]
    model = {ANCESTRY_CODES.get(w.lower(), w.upper()) for w in words}
    if "MOSTLY EUROPEAN (GTEX V8 DONORS)" in model:
        model = {"EUR"}
    overlap = model & participants
    if not overlap:
        return "mismatched"
    if participants <= model:
        return "matched" if len(model) == 1 else "mixed"
    return "partial"


def _genotype(project: Project, cfg: PredictedModalityConfig) -> tuple[str, Artifact]:
    from efgpp.data.genotype.carriers import resolve_genotype_artifact

    gid = cfg.genotype_artifact or (project.data.observed.genotype[0].id if project.data.observed.genotype else None)
    if gid is None:
        raise RuntimeError("genetically predicted data needs a genotype source")
    return gid, resolve_genotype_artifact(project, gid, cfg.use_qc_genotype)


def score_unit(project: Project, unit: PredictionUnit, cfg: PredictedModalityConfig, genotype: Artifact, *,
               step_id: str | None, threads: int) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, dict[str, Any]]:
    """Generic engines: harmonize the unit's weights to the genotype and score them.
    Returns (scores by IID, model QC, feature metadata, software info)."""
    from efgpp.data.genotype.formats import resolve_fileset
    from efgpp.data.predicted.harmonize import apply_coverage, harmonize
    from efgpp.data.predicted.scorer import (
        bed_sites,
        keyed_genotype,
        keyed_sites,
        score_plink2,
        score_python,
    )
    from efgpp.setup.tools import available

    weights, features = load_weights(unit, cfg.model_genome_build if unit.provider == "predictdb" else None)
    builds = weights.get_column("genome_build").drop_nulls().unique().to_list()
    model_build = builds[0] if len(builds) == 1 else None
    fs = resolve_fileset(Path(genotype.path), genotype.format)
    work = project.work_root / "predicted" / unit.modality.value / unit.key
    use_plink = available(project, "plink2")
    if use_plink:
        fs = keyed_genotype(project, fs, project.work_root / "scoring" / str(genotype.artifact_id), threads, step_id)
        sites = keyed_sites(fs)
    elif fs.format == "bed":
        sites = bed_sites(fs)
    else:
        raise RuntimeError("PLINK 2 is required to score this genotype format (efgpp setup tools plink2)")
    h = harmonize(weights, sites, genotype_build=genotype.genome_build, model_build=model_build,
                  forward_strand=unit.provider == "predictdb")
    qc = apply_coverage(h.qc, cfg.minimum_variant_coverage, cfg.allow_low_coverage, cfg.minimum_model_r2)
    below = features.select("model_id", "feature_id", "below_threshold") if "below_threshold" in features.columns \
        else None
    if below is not None:
        qc = qc.drop("below_threshold").join(below, on=["model_id", "feature_id"], how="left", maintain_order="left").with_columns(
            pl.col("below_threshold").fill_null(False))
        qc = qc.with_columns((pl.col("predict") & ~pl.col("below_threshold")).alias("predict"))
    todo = h.matched.join(qc.filter(pl.col("predict")).select("model_id", "feature_id"), on=["model_id", "feature_id"],
                          how="semi")
    if todo.height == 0:
        scores = None
    elif use_plink:
        scores = score_plink2(project, fs, todo, work, threads=threads, step_id=step_id)
    else:
        scores = score_python(fs, todo)
    write_parquet(h.excluded, work / "excluded_variants.parquet")
    return (scores if scores is not None else pl.DataFrame({"IID": []}, schema={"IID": pl.Utf8})), qc, features, {
        "backend": "plink2 --score" if use_plink else "python (numpy, .bed)", "model_build": model_build}


def metaxcan_unit(project: Project, unit: PredictionUnit, cfg: PredictedModalityConfig, genotype: Artifact, *,
                  step_id: str | None, threads: int) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, dict[str, Any]]:
    """PrediXcan (MetaXcan Predict.py) on a PredictDB database; coverage QC via the same harmonizer."""
    from efgpp.data.genotype.formats import read_variants, resolve_fileset
    from efgpp.data.predicted.harmonize import apply_coverage, harmonize
    from efgpp.data.predicted.metaxcan import ensure_vcf, parse_prediction
    from efgpp.data.provenance import run_tool

    weights, features = load_weights(unit)
    fs = resolve_fileset(Path(genotype.path), genotype.format)
    v = read_variants(fs)
    sites = v.select("variant_id", "chromosome", "position", pl.col("reference").alias("ref"),
                     pl.col("alternate").alias("alt"))
    h = harmonize(weights, sites, genotype_build=genotype.genome_build, model_build=cfg.model_genome_build,
                  forward_strand=True)
    qc = apply_coverage(h.qc, cfg.minimum_variant_coverage, cfg.allow_low_coverage, cfg.minimum_model_r2)
    vcf = ensure_vcf(project, genotype, step_id, threads)
    work = project.work_root / "predixcan" / unit.modality.value / unit.key
    work.mkdir(parents=True, exist_ok=True)
    pred, summary = work / "predict.txt", work / "summary.txt"
    args = ["--model_db_path", str(unit.model_path), "--vcf_genotypes", vcf.path, "--vcf_mode", "genotyped",
            "--prediction_output", str(pred), "--prediction_summary_output", str(summary), "--throw"]
    if cfg.variant_id_pattern:
        args += ["--on_the_fly_mapping", "METADATA", cfg.variant_id_pattern]
    rec = run_tool(project, "predixcan", [*args, *cfg.extra_args], step_id=step_id,
                   inputs=[vcf.artifact_id])  # type: ignore[list-item]
    table = parse_prediction(pred)
    keep = set(qc.filter(pl.col("predict")).get_column("feature_id").to_list())
    feats = [c for c in table.columns if c not in ("FID", "IID")]
    table = table.with_columns([pl.lit(None, dtype=pl.Float64).alias(c) for c in feats if c not in keep])
    write_parquet(h.excluded, work / "excluded_variants.parquet")
    return table, qc, features, {"backend": "MetaXcan Predict.py", "tool_version": rec.tool_version,
                                 "command": " ".join(rec.command), "summary": str(summary),
                                 "model_build": cfg.model_genome_build}


def run_unit(project: Project, modality: Modality, key: str, *, step_id: str | None = None,
             threads: int = 1) -> dict[str, Any]:
    from efgpp.data.genotype.carriers import genotype_samples, write_assays
    from efgpp.data.genotype.formats import resolve_fileset

    cfg, unit = find_unit(project, modality, key)
    if unit.reason:
        raise RuntimeError(f"{unit.step_id}: {unit.reason}" + (f" - {unit.required}" if unit.required else ""))
    gid, genotype = _genotype(project, cfg)
    runner = metaxcan_unit if unit.engine == "metaxcan" else score_unit
    scores, qc, features, software = runner(project, unit, cfg, genotype, step_id=step_id, threads=threads)
    samples = genotype_samples(project, gid, resolve_fileset(Path(genotype.path), genotype.format))
    by_iid = samples.select("IID", "native_id", "participant_id").unique("IID")
    scores = scores.drop([c for c in ("FID",) if c in scores.columns]).join(by_iid, on="IID", how="left", maintain_order="left")
    feature_cols = [c for c in scores.columns if c not in ("IID", "native_id", "participant_id")]
    table = scores.select(pl.col("native_id").fill_null(pl.col("IID")), "participant_id", *feature_cols)
    anc_status = ancestry_match_status(
        (unit.resource or {}).get("metadata", {}).get("training_ancestry") if unit.resource else None,
        participant_ancestries(project, gid))
    qc = qc.with_columns(pl.lit(anc_status).alias("ancestry_match_status"))
    out_dir = project.data_root / Origin.PREDICTED.value / modality.value / unit.output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    value_name = VALUE_NAMES[modality]
    n_cells = table.height * max(len(feature_cols), 1)
    if n_cells <= WIDE_LIMIT:
        pred = write_parquet(table, out_dir / f"predicted_{modality.value}.parquet")
        layout = "wide"
    else:
        long = table.unpivot(index=["native_id", "participant_id"], variable_name="feature_id",
                             value_name="predicted_value").with_columns(
            pl.lit(qc.get_column("model_id")[0] if qc.height else None).alias("model_id"))
        pred = out_dir / f"predicted_{modality.value}"
        long.with_columns(pl.col("feature_id").hash().mod(64).alias("bucket")).write_parquet(
            pred, partition_by="bucket")
        layout = "long (participant_id, feature_id, predicted_value, model_id), partitioned"
    qc_path = write_parquet(qc, out_dir / "model_qc.parquet")
    meta_path = write_parquet(features, out_dir / "feature_metadata.parquet")
    if modality == Modality.SPLICING and cfg.feature_mapping:
        mapping = pl.read_csv(project.resolve(cfg.feature_mapping), separator="\t", infer_schema=False)
        write_parquet(mapping, out_dir / "feature_mapping.parquet")
    res = unit.resource or {}
    res_meta = res.get("metadata") or {}
    ok = qc.filter(pl.col("status") == "OK")
    manifest = {
        "value": value_name, "origin": "predicted", "modality": modality.value, "engine": unit.engine,
        "provider": unit.provider, "dataset": unit.dataset, "tissue": unit.tissue or res_meta.get("tissue"),
        "platform": res_meta.get("platform"), "training_cohort": res_meta.get("training_cohort"),
        "training_ancestry": res_meta.get("training_ancestry"), "ancestry_match_status": anc_status,
        "model_resource": {"name": res.get("name"), "version": res.get("version"), "checksum": res.get("checksum"),
                           "local_path": res.get("local_path"), "genome_build": res.get("genome_build"),
                           "model_file": str(unit.model_path)},
        "genotype": {"source_id": gid, "artifact_id": genotype.artifact_id, "genome_build": genotype.genome_build},
        "software": software, "configuration": cfg.model_dump(mode="json"), "layout": layout,
        "participants": table.height, "features_predicted": int(ok.height if unit.engine != "metaxcan" else
                                                                qc.filter(pl.col("predict")).height),
        "features_in_model": qc.height,
        "low_coverage_features": int((qc.get_column("status") == "LOW_COVERAGE").sum()) if qc.height else 0,
        "median_coverage": statistics.median(qc.get_column("coverage_fraction").to_list()) if qc.height else None,
        "median_validation_r2": (statistics.median(r2) if (r2 := qc.get_column("validation_r2").drop_nulls().to_list())
                                 else None),
        "created_at": utcnow().isoformat(timespec="seconds") + "Z",
        "limitation": LIMITATIONS[modality],
        "phenotype_used": False,
    }
    manifest_path = out_dir / "manifest.yaml"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    files = [pred, qc_path, meta_path, manifest_path]
    source_id = f"PRED_{modality.value.upper()}_{unit.output_subdir.replace('/', '_')}"
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=source_id, artifact_type="predicted", modality=modality, origin=Origin.PREDICTED,
            status=ArtifactStatus.READY, path=str(pred), format="parquet" if layout == "wide" else "parquet_dataset",
            size=sum(f.stat().st_size for d in files for f in ([d] if d.is_file() else d.rglob("*")) if f.is_file()),
            checksum=str(checksum_paths(files, full_limit_bytes=_limit(project))), participant_count=table.get_column("participant_id").drop_nulls()
            .n_unique(), feature_count=len(feature_cols), genome_build=genotype.genome_build,
            tissue=unit.tissue or res_meta.get("tissue"), temporal_type=TemporalType.GENETICALLY_PREDICTED_STATIC,
            source_id=source_id, parent_artifact_ids=[genotype.artifact_id],  # type: ignore[list-item]
            tool=unit.engine, tool_version=software.get("tool_version"),
            resource_versions={str(res.get("name") or unit.provider): str(res.get("version") or unit.dataset),
                               "model_checksum": str(res.get("checksum") or checksum_paths([unit.model_path])
                                                     if unit.model_path else "")},
            metadata={"members": [str(f) for f in files], "value": value_name, "manifest": str(manifest_path),
                      "note": "genetically predicted; not measured", "layout": layout,
                      "provider": unit.provider, "dataset": unit.dataset},
        ))
        write_assays(reg, art, source_id, modality, table.get_column("participant_id").drop_nulls().to_list(),
                     origin=Origin.PREDICTED, tissue=unit.tissue)
        _catalog_models(reg, unit, features, qc, res)
    return {"artifact": art.artifact_id, "features": len(feature_cols), "participants": table.height,
            "low_coverage": manifest["low_coverage_features"]}


def _catalog_models(reg: Registry, unit: PredictionUnit, features: pl.DataFrame, qc: pl.DataFrame,
                    res: dict[str, Any]) -> None:
    """Persistent model catalogue (registry table molecular_models)."""
    meta = res.get("metadata") or {}
    cat = features.select("model_id", "feature_id", "feature_name", "validation_r2", "n_variants").join(
        qc.select("model_id", "feature_id", "coverage_fraction", "status"), on=["model_id", "feature_id"], how="left", maintain_order="left")
    cat = cat.with_columns(
        pl.lit(unit.provider).alias("provider"), pl.lit(unit.dataset).alias("provider_dataset_id"),
        pl.lit(unit.modality.value).alias("modality"), pl.lit(unit.tissue or meta.get("tissue")).alias("tissue"),
        pl.lit(meta.get("platform")).alias("platform"), pl.lit(meta.get("training_cohort")).alias("training_cohort"),
        pl.lit(meta.get("training_ancestry")).alias("training_ancestry"),
        pl.lit(res.get("genome_build") or meta.get("genome_build")).alias("genome_build"),
        pl.lit(res.get("version")).alias("resource_version"), pl.lit(res.get("checksum")).alias("resource_sha256"),
        pl.lit(str(unit.model_path)).alias("local_path"))
    reg.execute("DELETE FROM molecular_models WHERE provider = ? AND provider_dataset_id = ? AND modality = ? "
                "AND coalesce(tissue, '') = coalesce(?, '')",
                [unit.provider, unit.dataset, unit.modality.value, unit.tissue or meta.get("tissue")])
    if cat.height:
        reg.insert_frame("molecular_models", cat.select(
            "model_id", "provider", "provider_dataset_id", "modality", "feature_id", "feature_name", "tissue",
            "platform", "training_cohort", "training_ancestry", "genome_build", "validation_r2", "n_variants",
            "coverage_fraction", "status", "resource_version", "resource_sha256", "local_path"), replace=False)


# ------------------------------------------------------------------ plan (no prediction)
@dataclass
class PlanRow:
    item: str
    status: str
    detail: str = ""
    required: str = ""
    extra: dict[str, Any] = field(default_factory=dict)


def prediction_plan(project: Project) -> list[PlanRow]:
    """What `efgpp data predict all` would do and what is missing - nothing is computed."""
    from efgpp.data.genotype.carriers import carrier_dataset, resolve_genotype_artifact
    from efgpp.data.references.genome import installed_fasta
    from efgpp.setup.tools import available

    rows: list[PlanRow] = []
    genos = project.data.observed.genotype
    if not genos:
        rows.append(PlanRow("Genotype", "MISSING", required="efgpp data add genotype --path <prefix>"))
    for g in genos:
        try:
            art = resolve_genotype_artifact(project, g.id, True)
            kind = "QC-passed" if art.artifact_type == "qc_genotype" else art.artifact_type
            rows.append(PlanRow(f"Genotype {g.id}", "READY", f"{kind}, {art.genome_build}, "
                                f"{art.participant_count or '?'} participants, {art.feature_count or '?'} variants"))
        except RuntimeError as exc:
            rows.append(PlanRow(f"Genotype {g.id}", "NOT READY", str(exc), "efgpp data prepare"))
        pv = project.data.participant_variants
        if pv.enabled:
            rows.append(PlanRow(f"Participant variants {g.id}", "READY" if carrier_dataset(project, g.id) else "PLANNED",
                                "carrier table" if carrier_dataset(project, g.id) else "built by `efgpp data prepare`"))
    rows.append(PlanRow("Reference FASTA", "INSTALLED" if installed_fasta(project) else "NOT INSTALLED",
                        required="" if installed_fasta(project) else "efgpp resources install genome"))
    for tool, cmd in (("plink2", "efgpp setup tools plink2"), ("vep", "efgpp setup tools vep"),
                      ("spliceai", "efgpp setup toolkit spliceai"), ("predixcan", "efgpp setup toolkit metaxcan")):
        ok = available(project, tool)
        rows.append(PlanRow(f"Software {tool}", "INSTALLED" if ok else "NOT INSTALLED", required="" if ok else cmd))
    for modality, cfg in project.data.predicted.items():
        if not cfg.enabled:
            rows.append(PlanRow(f"Predicted {modality.value}", "DISABLED",
                                required=f"efgpp data predict enable {modality.value} ..."))
            continue
        for u in plan_units(project, modality, cfg):
            label = f"Predicted {modality.value} [{u.key}]"
            if u.reason:
                rows.append(PlanRow(label, "NO DATASET SELECTED" if "no dataset" in u.reason else "NOT INSTALLED",
                                    u.reason, u.required or ""))
            else:
                rows.append(PlanRow(label, "READY", f"{u.engine}; {u.model_path}"))
    return rows


def _limit(project: Project) -> int:
    """Same full-checksum size limit as snapshot verification (larger files: labelled sampled hash)."""
    return int(project.config.storage.full_checksum_limit_gb * 1024**3)
