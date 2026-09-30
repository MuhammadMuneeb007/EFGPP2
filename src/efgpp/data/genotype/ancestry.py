"""Genetic ancestry estimation by projection onto a labelled reference panel.

Requires `resources.yaml: ld_reference` with
    path:         PLINK2 prefix (pgen/pvar/psam) of the reference panel
    populations:  TSV with columns IID and population
    genome_build: optional; detected otherwise. Panels not in the target build (GRCh38) are
                  lifted once with pyliftover into resources/ld_reference/<name>_GRCh38.
Both reference and cohort samples are projected with the *same* PLINK2 --score command so
their coordinates are directly comparable; each participant is assigned the population
of the nearest reference centroid, with the distance ratio as a confidence measure.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import polars as pl

from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.formats import read_variants, resolve_fileset
from efgpp.data.genotype.loader import read_plink_table, write_id_list
from efgpp.data.io import write_parquet
from efgpp.data.provenance import run_tool
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths


def ancestry_configured(project: Project) -> bool:
    ref = project.resources.ld_reference
    extra = ref.model_extra or {}
    return bool(ref.enabled and ref.path and extra.get("populations"))


def score_columns(allele_file: Path) -> tuple[int, int, list[int]]:
    """1-based (ID column, allele column, PC columns) of a PLINK2 .eigenvec.allele file."""
    with open(allele_file, encoding="utf-8") as fh:
        header = fh.readline().lstrip("#").split()
    id_col = header.index("ID") + 1
    allele_col = header.index("A1") + 1
    pcs = [i + 1 for i, h in enumerate(header) if h.startswith("PC")]
    return id_col, allele_col, pcs


def nearest_centroid(target: np.ndarray, ref: np.ndarray, labels: list[str]) -> tuple[list[str], np.ndarray]:
    pops = sorted(set(labels))
    lab = np.array(labels)
    centroids = np.stack([ref[lab == p].mean(axis=0) for p in pops])
    d = np.linalg.norm(target[:, None, :] - centroids[None, :, :], axis=2)
    order = np.argsort(d, axis=1)
    best = d[np.arange(len(d)), order[:, 0]]
    second = d[np.arange(len(d)), order[:, 1]] if len(pops) > 1 else np.full(len(d), np.inf)
    confidence = 1 - best / np.where(second > 0, second, np.inf)
    return [pops[i] for i in order[:, 0]], confidence


def run_ancestry(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    ref_cfg = project.resources.ld_reference
    pops_path = project.resolve((ref_cfg.model_extra or {})["populations"])
    ref_prefix = project.resolve(ref_cfg.path)  # type: ignore[arg-type]
    n_pcs = project.config.genotype.pca.n_pcs
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        qc = store.latest(source_id=source_id, artifact_type="qc_genotype")
        mapping = dict(reg.frame("SELECT native_id, participant_id FROM sample_aliases WHERE source_id = ?",
                                 [source_id]).iter_rows())
    if qc is None:
        raise RuntimeError(f"{source_id}: run genotype QC before ancestry")
    work = project.work_root / "ancestry" / source_id
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)

    cohort_ids = set(read_variants(resolve_fileset(Path(qc.path), "pgen")).get_column("variant_id").to_list())
    # The reference panel must be in the project's target build (lifted once if it is not).
    from efgpp.data.genotype.harmonization import ensure_target_build

    ref_fs = ensure_target_build(
        project, resolve_fileset(ref_prefix, "auto"), (ref_cfg.model_extra or {}).get("genome_build"),
        project.resource_root / "ld_reference" / ref_prefix.name, step_id=step_id, threads=threads)
    overlap = [v for v in read_variants(ref_fs).get_column("variant_id").to_list() if v in cohort_ids]
    if len(overlap) < 100:
        raise RuntimeError(f"only {len(overlap)} variants shared with the reference panel (IDs must match)")
    shared = write_id_list(work / "shared.txt", overlap)

    base = ["--threads", str(threads)]
    run_tool(project, "plink2", [*ref_fs.plink_input_args(), *base, "--extract", str(shared), "--maf", "0.01",
                                 "--freq", "counts", "--pca", str(n_pcs), "allele-wts", "--out", str(work / "ref")],
             step_id=step_id)
    allele_file = work / "ref.eigenvec.allele"
    id_col, a1_col, pc_cols = score_columns(allele_file)
    score = ["--read-freq", str(work / "ref.acount"), "--score", str(allele_file), str(id_col), str(a1_col),
             "header-read", "no-mean-imputation", "variance-standardize",
             "--score-col-nums", f"{pc_cols[0]}-{pc_cols[-1]}"]
    run_tool(project, "plink2", [*ref_fs.plink_input_args(), *base, *score, "--out", str(work / "ref_proj")], step_id=step_id)
    rec = run_tool(project, "plink2", ["--pfile", qc.path, *base, *score, "--out", str(work / "cohort_proj")],
                   step_id=step_id, inputs=[qc.artifact_id])  # type: ignore[list-item]

    def load(p: Path) -> pl.DataFrame:
        df = read_plink_table(p)
        pcs = [c for c in df.columns if c.endswith("_AVG") or c.endswith("_SUM")][: len(pc_cols)]
        return df.select("IID", *[pl.col(c).cast(pl.Float64).alias(f"PC{i + 1}") for i, c in enumerate(pcs)])

    ref_proj = load(work / "ref_proj.sscore").join(
        pl.read_csv(pops_path, separator="\t", schema_overrides={"IID": pl.Utf8}).select("IID", "population"), on="IID")
    coh = load(work / "cohort_proj.sscore")
    pc_names = [c for c in coh.columns if c.startswith("PC")]
    labels, conf = nearest_centroid(coh.select(pc_names).to_numpy(), ref_proj.select(pc_names).to_numpy(),
                                    ref_proj.get_column("population").to_list())
    result = coh.with_columns(
        pl.col("IID").replace_strict(mapping, default=None).alias("participant_id"),
        pl.Series("ancestry", labels), pl.Series("ancestry_confidence", conf),
    )
    out_dir = project.artifact_dir(Origin.DERIVED, Modality.ANCESTRY, source_id)
    out = write_parquet(result, out_dir / f"{source_id}_ancestry.parquet")
    write_parquet(ref_proj, out_dir / f"{source_id}_reference_projection.parquet")
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=f"{source_id}_ancestry", artifact_type="ancestry", modality=Modality.ANCESTRY,
            origin=Origin.DERIVED, status=ArtifactStatus.READY, path=str(out), format="parquet",
            size=out.stat().st_size, checksum=str(checksum_paths([out])), participant_count=result.height,
            feature_count=len(pc_names), genome_build=qc.genome_build, temporal_type=TemporalType.STATIC,
            source_id=source_id, parent_artifact_ids=[qc.artifact_id],  # type: ignore[list-item]
            tool="plink2", tool_version=rec.tool_version,
            resource_versions={"ld_reference": ref_cfg.version or "unversioned"},
            metadata={"method": "nearest reference centroid in projected PC space", "shared_variants": len(overlap)},
        ))
    return art.artifact_id  # type: ignore[return-value]
