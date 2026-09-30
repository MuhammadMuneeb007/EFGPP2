"""QC principal components (population-structure description).

QC PCA != modelling PCA. These components describe the whole QC-passed cohort for QC and
reporting. PCs used as model covariates must be refit inside training folds by the
Representation layer.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.loader import read_eigenval, read_plink_table
from efgpp.data.io import write_parquet
from efgpp.data.provenance import run_tool
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths
from efgpp.setup.tools import available


def _ids(*arts: Artifact | None) -> list[str]:
    return [a.artifact_id for a in arts if a is not None and a.artifact_id]


def _qc_inputs(project: Project, source_id: str) -> tuple[Artifact, Artifact | None, dict[str, str]]:
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        qc = store.latest(source_id=source_id, artifact_type="qc_genotype")
        prune = store.latest(source_id=source_id, artifact_type="ld_prune")
        mapping = dict(reg.frame("SELECT native_id, participant_id FROM sample_aliases WHERE source_id = ?",
                                 [source_id]).iter_rows())
    if qc is None:
        raise RuntimeError(f"{source_id}: run genotype QC before PCA")
    return qc, prune, mapping


def _attach_participants(df: pl.DataFrame, project: Project, source_id: str, mapping: dict[str, str]) -> pl.DataFrame:
    _, src = project.data.get_source(source_id)
    sid = src.sample_id  # type: ignore[attr-defined]
    if sid.mode == "fid_iid" and "FID" in df.columns:
        native = pl.col("FID") + sid.separator + pl.col("IID")
    else:
        native = pl.col(sid.column if sid.column in df.columns else "IID")
    return df.with_columns(native.replace_strict(mapping, default=None).alias("participant_id"))


def run_qc_pca(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> list[str]:
    cfg = project.config.genotype.pca
    qc, prune, mapping = _qc_inputs(project, source_id)
    out_dir = project.artifact_dir(Origin.DERIVED, Modality.QC_PCA, source_id)
    n_samples = qc.participant_count or 0
    n_variants = prune.feature_count if prune and prune.feature_count else (qc.feature_count or 0)
    n_pcs = max(1, min(cfg.n_pcs, n_samples - 1, n_variants - 1))
    registered: list[str] = []

    engines = [cfg.qc_engine, *[e for e in cfg.optional_engines if e != cfg.qc_engine]]
    for engine in engines:
        if engine == "flashpca2" and not available(project, "flashpca2"):
            continue
        prefix = out_dir / f"{source_id}_pca_{engine}"
        if engine == "plink2":
            approx = cfg.approx == "always" or (cfg.approx == "auto" and n_samples > 5000)
            args = ["--pfile", qc.path, "--threads", str(threads), "--pca", str(n_pcs)]
            if approx:
                args.append("approx")
            if prune is not None:
                args += ["--extract", prune.path]
            rec = run_tool(project, "plink2", [*args, "--out", str(prefix)], step_id=step_id,
                           inputs=_ids(qc, prune))
            vec = read_plink_table(Path(str(prefix) + ".eigenvec"))
            vals = read_eigenval(Path(str(prefix) + ".eigenval"))
        else:
            bed_prefix = out_dir / f"{source_id}_pruned_bed"
            args = ["--pfile", qc.path, "--threads", str(threads), "--make-bed", "--out", str(bed_prefix)]
            if prune is not None:
                args[4:4] = ["--extract", prune.path]
            run_tool(project, "plink2", args, step_id=step_id, inputs=[qc.artifact_id])  # type: ignore[list-item]
            pcs_txt, val_txt = Path(str(prefix) + ".pcs.txt"), Path(str(prefix) + ".eigenval.txt")
            rec = run_tool(project, "flashpca2", ["--bfile", str(bed_prefix), "--ndim", str(n_pcs),
                                                   "--outpc", str(pcs_txt), "--outval", str(val_txt),
                                                   "--numthreads", str(threads)], step_id=step_id)
            vec = read_plink_table(pcs_txt)
            vals = [float(x) for x in val_txt.read_text(encoding="utf-8").split()]
        vec = _attach_participants(vec, project, source_id, mapping)
        total = sum(vals) or 1.0
        pcs_path = write_parquet(vec, Path(str(prefix) + ".parquet"))
        eig_path = write_parquet(
            pl.DataFrame({"pc": [f"PC{i + 1}" for i in range(len(vals))], "eigenvalue": vals,
                          "fraction_of_listed_variance": [v / total for v in vals]}),
            Path(str(prefix) + "_eigenvalues.parquet"),
        )
        with Registry.open(project) as reg:
            store = ArtifactStore(reg)
            olds = [a for a in store.find(source_id=source_id, artifact_type="qc_pca")
                    if a.artifact_name == f"{source_id}_qc_pca_{engine}"]
            art = store.register(Artifact(
                artifact_name=f"{source_id}_qc_pca_{engine}", artifact_type="qc_pca", modality=Modality.QC_PCA,
                origin=Origin.DERIVED, status=ArtifactStatus.READY, path=str(pcs_path), format="parquet",
                size=pcs_path.stat().st_size, checksum=str(checksum_paths([pcs_path, eig_path])),
                participant_count=vec.height, feature_count=n_pcs, genome_build=qc.genome_build,
                temporal_type=TemporalType.STATIC, source_id=source_id,
                parent_artifact_ids=_ids(qc, prune),
                tool=engine if engine != "plink2" else "plink2", tool_version=rec.tool_version,
                configuration_hash=cfg.checksum(), command=" ".join(rec.command),
                metadata={"engine": engine, "eigenvalues": str(eig_path), "purpose": "qc",
                          "members": [str(pcs_path), str(eig_path)],
                          "note": "QC PCA; refit PCs inside training folds for modelling"},
            ))
            for o in olds:
                store.supersede(o.artifact_id, art.artifact_id)  # type: ignore[arg-type]
            registered.append(art.artifact_id)  # type: ignore[arg-type]
            if engine == cfg.qc_engine:
                present = vec.filter(pl.col("participant_id").is_not_null())
                reg.execute("DELETE FROM assays WHERE source_id = ?", [f"{source_id}_PCA"])
                reg.insert_frame("assays", pl.DataFrame({
                    "assay_id": [f"{art.artifact_id}:{i}" for i in range(present.height)],
                    "artifact_id": [art.artifact_id] * present.height,
                    "source_id": [f"{source_id}_PCA"] * present.height,
                    "modality": [Modality.QC_PCA.value] * present.height,
                    "origin": [Origin.DERIVED.value] * present.height,
                    "participant_id": present.get_column("participant_id"),
                }))
    return registered
