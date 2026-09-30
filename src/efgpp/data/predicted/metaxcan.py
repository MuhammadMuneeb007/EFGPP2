"""MetaXcan / PrediXcan (Predict.py) helpers: VCF input, output parsing, model metadata.

Individual-level PrediXcan only (never S-PrediXcan, which needs GWAS summary statistics).
Units are run by `efgpp.data.predicted.engine` (engine: metaxcan).
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import polars as pl

from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.conversion import convert
from efgpp.data.registry import Registry
from efgpp.project import Project


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
