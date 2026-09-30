"""Genetically predicted molecular traits with MetaXcan / PrediXcan (Predict.py).

Outputs live under data/predicted/<modality>/<tissue>/ with origin PREDICTED and temporal
type genetically_predicted_static. They are never registered as measured data.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import polars as pl

from efgpp.config.data import PredictedModalityConfig
from efgpp.constants import ArtifactStatus, GenomeBuild, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.conversion import convert
from efgpp.data.io import write_parquet
from efgpp.data.provenance import run_tool
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths


def predicted_source_id(modality: Modality, tissue: str) -> str:
    return f"PRED_{modality.value.upper()}_{tissue}"


def model_path(project: Project, cfg: PredictedModalityConfig, tissue: str) -> Path:
    mp = cfg.model_provider
    if tissue in mp.model_paths:
        return project.resolve(mp.model_paths[tissue])
    if mp.models_dir:
        base = project.resolve(mp.models_dir)
    else:
        # Models installed with `efgpp resources install predictdb` (versioned folder), else resources/predictdb.
        from efgpp.data.annotation import resource_version

        _, installed = resource_version(project, "predictdb")
        installed = project.resolve(str(installed)) if installed else None
        base = installed if installed and installed.is_dir() else project.resource_root / "predictdb"
    return base / f"{mp.model_prefix}{tissue}{mp.model_suffix}"


def model_metadata(db: Path) -> dict[str, str]:
    """Read the `construction`/`extra` tables PredictDB models ship with, when present."""
    out: dict[str, str] = {}
    try:
        con = sqlite3.connect(f"file:{db.as_posix()}?mode=ro", uri=True)
        try:
            tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if "extra" in tables:
                out["n_genes"] = str(con.execute("SELECT count(*) FROM extra").fetchone()[0])
            if "weights" in tables:
                out["n_weights"] = str(con.execute("SELECT count(*) FROM weights").fetchone()[0])
        finally:
            con.close()
    except sqlite3.Error:
        pass
    return out


def genotype_for_prediction(project: Project, cfg: PredictedModalityConfig) -> Artifact:
    gid = cfg.genotype_artifact or (project.data.observed.genotype[0].id if project.data.observed.genotype else None)
    if gid is None:
        raise RuntimeError("predicted modalities need a genotype source")
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        art = store.latest(source_id=gid, artifact_type="qc_genotype") if cfg.use_qc_genotype else None
        art = art or store.latest(source_id=gid, artifact_type="source")
    if art is None:
        raise RuntimeError(f"{gid} is not registered")
    return art


def ensure_vcf(project: Project, genotype: Artifact, step_id: str | None, threads: int) -> Artifact:
    with Registry.open(project) as reg:
        existing = ArtifactStore(reg).find(source_id=genotype.source_id, artifact_type="converted_genotype")
    for a in existing:
        if a.format == "vcf" and genotype.artifact_id in a.parent_artifact_ids and Path(a.path).exists():
            return a
    new_id = convert(project, genotype.artifact_id, "vcf", step_id=step_id, threads=threads)  # type: ignore[arg-type]
    with Registry.open(project) as reg:
        return ArtifactStore(reg).get(new_id)


def parse_prediction(path: Path) -> pl.DataFrame:
    df = pl.read_csv(path, separator="\t", infer_schema_length=0)
    genes = [c for c in df.columns if c not in ("FID", "IID")]
    return df.with_columns([pl.col(g).cast(pl.Float64, strict=False) for g in genes])


def run_prediction(project: Project, modality: Modality, tissue: str, *, step_id: str | None = None,
                   threads: int = 1) -> str:
    cfg: PredictedModalityConfig = getattr(project.data.predicted, modality.value)
    genotype = genotype_for_prediction(project, cfg)
    geno_build = GenomeBuild.normalize(genotype.genome_build)
    if geno_build != GenomeBuild.normalize(cfg.model_genome_build):
        raise RuntimeError(
            f"genotype build {genotype.genome_build} != model build {cfg.model_genome_build}; "
            "lift the genotype over explicitly (efgpp data liftover) - EFGPP never converts silently"
        )
    db = model_path(project, cfg, tissue)
    if not db.exists():
        raise RuntimeError(f"PredictDB model not found for {tissue}: {db}")
    vcf = ensure_vcf(project, genotype, step_id, threads)
    work = project.work_root / "predixcan" / tissue
    work.mkdir(parents=True, exist_ok=True)
    pred, summary = work / f"{tissue}_predict.txt", work / f"{tissue}_summary.txt"
    args = ["--model_db_path", str(db), "--vcf_genotypes", vcf.path, "--vcf_mode", "genotyped",
            "--prediction_output", str(pred), "--prediction_summary_output", str(summary), "--throw"]
    if cfg.variant_id_pattern:
        args += ["--on_the_fly_mapping", "METADATA", cfg.variant_id_pattern]
    rec = run_tool(project, "predixcan", [*args, *cfg.extra_args], step_id=step_id,
                   inputs=[vcf.artifact_id])  # type: ignore[list-item]

    table = parse_prediction(pred)
    with Registry.open(project) as reg:
        mapping = dict(reg.frame("SELECT native_id, participant_id FROM sample_aliases WHERE source_id = ?",
                                 [genotype.source_id]).iter_rows())
    table = table.with_columns(pl.col("IID").replace_strict(mapping, default=None).alias("participant_id"))
    genes = [c for c in table.columns if c not in ("FID", "IID", "participant_id")]
    source_id = predicted_source_id(modality, tissue)
    out_dir = project.artifact_dir(Origin.PREDICTED, modality, tissue)
    out = write_parquet(
        table.select(pl.col("IID").alias("native_id"), "participant_id", *genes),
        out_dir / f"predicted_{modality.value}_{tissue}.parquet",
    )
    summary_out = write_parquet(pl.read_csv(summary, separator="\t", infer_schema_length=0),
                                out_dir / f"predicted_{modality.value}_{tissue}_model_summary.parquet")
    version = cfg.model_provider.version or db.name
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        art = store.register_replacing(Artifact(
            artifact_name=f"{source_id}", artifact_type="predicted", modality=modality, origin=Origin.PREDICTED,
            status=ArtifactStatus.READY, path=str(out), format="parquet",
            size=out.stat().st_size + summary_out.stat().st_size,
            checksum=str(checksum_paths([out, summary_out])),
            participant_count=table.get_column("participant_id").drop_nulls().n_unique(),
            feature_count=len(genes), genome_build=genotype.genome_build, tissue=tissue,
            temporal_type=TemporalType.GENETICALLY_PREDICTED_STATIC, source_id=source_id,
            parent_artifact_ids=[genotype.artifact_id, vcf.artifact_id],  # type: ignore[list-item]
            tool="metaxcan", tool_version=rec.tool_version,
            resource_versions={"predictdb": version, "model_sha256": str(checksum_paths([db]))},
            command=" ".join(rec.command),
            metadata={"model": str(db), "model_info": model_metadata(db), "summary": str(summary_out),
                      "members": [str(out), str(summary_out)],
                      "note": "genetically predicted; not measured expression"},
        ))
        present = table.filter(pl.col("participant_id").is_not_null())
        reg.execute("DELETE FROM assays WHERE source_id = ?", [source_id])
        reg.insert_frame("assays", pl.DataFrame({
            "assay_id": [f"{art.artifact_id}:{i}" for i in range(present.height)],
            "artifact_id": [art.artifact_id] * present.height,
            "source_id": [source_id] * present.height,
            "modality": [modality.value] * present.height,
            "origin": [Origin.PREDICTED.value] * present.height,
            "participant_id": present.get_column("participant_id"),
            "tissue": [tissue] * present.height,
        }))
    return art.artifact_id  # type: ignore[return-value]
