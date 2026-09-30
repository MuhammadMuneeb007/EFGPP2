"""Runs of homozygosity with PLINK 1.9 (--homozyg is not implemented in PLINK 2)."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.loader import read_plink_table
from efgpp.data.io import write_parquet
from efgpp.data.provenance import run_tool
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths


def run_roh(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    with Registry.open(project) as reg:
        qc = ArtifactStore(reg).latest(source_id=source_id, artifact_type="qc_genotype")
    if qc is None:
        raise RuntimeError(f"{source_id}: run genotype QC before ROH")
    out_dir = project.artifact_dir(Origin.DERIVED, Modality.ROH, source_id)
    work = project.work_root / "roh" / source_id
    work.mkdir(parents=True, exist_ok=True)
    bed = work / f"{source_id}_qc_bed"
    run_tool(project, "plink2", ["--pfile", qc.path, "--threads", str(threads), "--make-bed", "--out", str(bed)],
             step_id=step_id, inputs=[qc.artifact_id])  # type: ignore[list-item]
    rec = run_tool(project, "plink", ["--bfile", str(bed), "--homozyg", "--out", str(work / "roh")],
                   step_id=step_id, inputs=[qc.artifact_id])  # type: ignore[list-item]
    segments = read_plink_table(work / "roh.hom")
    per_sample = read_plink_table(work / "roh.hom.indiv")
    seg_path = write_parquet(segments, out_dir / f"{source_id}_roh_segments.parquet")
    ind_path = write_parquet(per_sample, out_dir / f"{source_id}_roh_per_sample.parquet")
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=f"{source_id}_roh", artifact_type="roh", modality=Modality.ROH, origin=Origin.DERIVED,
            status=ArtifactStatus.READY, path=str(ind_path), format="parquet",
            size=ind_path.stat().st_size + seg_path.stat().st_size,
            checksum=str(checksum_paths([ind_path, seg_path])), participant_count=per_sample.height,
            genome_build=qc.genome_build, temporal_type=TemporalType.STATIC, source_id=source_id,
            parent_artifact_ids=[qc.artifact_id], tool="plink", tool_version=rec.tool_version,  # type: ignore[list-item]
            command=" ".join(rec.command), metadata={"segments": str(seg_path), "members": [str(ind_path), str(seg_path)]},
        ))
    return art.artifact_id  # type: ignore[return-value]


def summarize_roh(path: Path) -> dict[str, float]:
    df = pl.read_parquet(path)
    if df.height == 0 or "KB" not in df.columns:
        return {}
    kb = df.get_column("KB").cast(pl.Float64)
    return {"samples": df.height, "mean_kb": float(kb.mean() or 0), "max_kb": float(kb.max() or 0)}  # type: ignore[arg-type]
