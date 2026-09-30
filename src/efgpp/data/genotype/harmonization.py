"""Genome harmonization: bring genotypes (and reference panels) to the project's target build.

Coordinates are converted with pyliftover (UCSC chain files kept in resources/liftover/) and
PLINK 2 writes the lifted fileset:

    plink2 <input> --exclude <not mappable> --update-map <new positions> --sort-vars --make-pgen

Liftover always creates a *new* artifact; the original coordinates are never replaced.
Variants that do not map, map to another chromosome, or map to the reverse strand (which
would need allele complementing) are excluded and listed in a report, never guessed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import polars as pl

from efgpp.constants import ArtifactStatus, GenomeBuild, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.formats import (
    GenotypeFileset,
    count_variants,
    read_variants,
    resolve_fileset,
)
from efgpp.data.genotype.liftover import Lifter, chain_file, lift_table
from efgpp.data.genotype.loader import write_id_list
from efgpp.data.io import write_parquet
from efgpp.data.provenance import RunRecord, run_tool
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths


def target_build(project: Project) -> GenomeBuild:
    return GenomeBuild.normalize(project.config.defaults.target_build)


def liftover_fileset(project: Project, fs: GenotypeFileset, source: GenomeBuild, target: GenomeBuild,
                     out_prefix: Path, *, step_id: str | None = None, threads: int = 1,
                     inputs: list[str] | None = None) -> tuple[GenotypeFileset, dict[str, int], Path, RunRecord]:
    """Lift any PLINK-readable fileset with pyliftover + PLINK 2. Returns (new fileset, counts,
    report path, run record)."""
    lifter = Lifter(chain_file(project.resource_root, source, target))
    variants = read_variants(fs)
    mapping = lift_table(variants, lifter)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    report = write_parquet(mapping, out_prefix.with_name(out_prefix.name + "_liftover_report.parquet"))
    mapped = mapping.filter(pl.col("status") == "mapped")
    work = project.work_root / "liftover" / out_prefix.name
    work.mkdir(parents=True, exist_ok=True)
    exclude = write_id_list(work / "exclude.txt", mapping.filter(pl.col("status") != "mapped").get_column("variant_id").to_list())
    update = work / "update_map.txt"
    update.write_text("".join(f"{v}\t{p}\n" for v, p in mapped.select("variant_id", "new_position").iter_rows()),
                      encoding="utf-8")
    rec = run_tool(project, "plink2", [*fs.plink_input_args(), "--threads", str(threads), "--exclude", str(exclude),
                                       "--update-map", str(update), "--sort-vars", "--make-pgen", "--out", str(out_prefix)],
                   step_id=step_id, inputs=inputs or [])
    counts = {s: int(n) for s, n in mapping.group_by("status").len().iter_rows()}
    return resolve_fileset(out_prefix, "pgen"), counts, report, rec


def _register_lifted(project: Project, src: Artifact, lifted: GenotypeFileset, target: GenomeBuild,
                     counts: dict[str, int], report: Path, rec: RunRecord, chain: Path) -> Artifact:
    with Registry.open(project) as reg:
        return ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=f"{src.source_id}_{target.value}", artifact_type="liftover_genotype",
            modality=src.modality, origin=Origin.DERIVED, status=ArtifactStatus.READY, path=str(lifted.prefix),
            format="pgen", size=sum(m.stat().st_size for m in lifted.members),
            checksum=str(checksum_paths(lifted.members)), participant_count=src.participant_count,
            feature_count=count_variants(lifted), genome_build=target.value, temporal_type=TemporalType.STATIC,
            source_id=src.source_id, parent_artifact_ids=[src.artifact_id],  # type: ignore[list-item]
            tool="pyliftover+plink2", tool_version=rec.tool_version,
            resource_versions={"chain": chain.name}, command=" ".join(rec.command),
            metadata={"from_build": src.genome_build, "liftover_counts": counts, "report": str(report),
                      "members": [str(m) for m in lifted.members]},
        ))


def run_liftover(project: Project, artifact_id: str, target_build_name: str, *, step_id: str | None = None,
                 threads: int = 1) -> str:
    """Lift one registered genotype artifact to `target_build_name` (new artifact)."""
    with Registry.open(project) as reg:
        src = ArtifactStore(reg).get(artifact_id)
    source = GenomeBuild.normalize(src.genome_build)
    target = GenomeBuild.normalize(target_build_name)
    if source not in (GenomeBuild.GRCH37, GenomeBuild.GRCH38):
        raise RuntimeError(f"{artifact_id}: source genome build is {src.genome_build!r}; set it before lifting over")
    if source == target:
        raise RuntimeError(f"{artifact_id} is already {target.value}")
    fs = resolve_fileset(Path(src.path), src.format)
    out_dir = project.artifact_dir(Origin.DERIVED, Modality.GENOTYPE_QC, src.source_id or "lifted", "liftover")
    prefix = out_dir / f"{src.source_id or src.artifact_name}_{target.value}"
    lifted, counts, report, rec = liftover_fileset(project, fs, source, target, prefix, step_id=step_id,
                                                   threads=threads, inputs=[artifact_id])
    art = _register_lifted(project, src, lifted, target, counts, report, rec, chain_file(project.resource_root, source, target))
    return art.artifact_id  # type: ignore[return-value]


def harmonize_genotype(project: Project, source_id: str, *, step_id: str | None = None,
                       threads: int = 1) -> dict[str, Any]:
    """Plan step: make sure the genotype is in the target build (lift it when it is not)."""
    from efgpp.data.genotype.qc import variant_table
    from efgpp.data.schemas.artifact import VARIANT_SCHEMA

    target = target_build(project)
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        src = store.latest(source_id=source_id, artifact_type="source")
        existing = store.latest(source_id=source_id, artifact_type="liftover_genotype")
    if src is None:
        raise RuntimeError(f"{source_id} is not registered")
    build = GenomeBuild.normalize(src.genome_build)
    if build == target:
        return {"build": build.value, "action": f"already {target.value}; nothing to lift"}
    if build not in (GenomeBuild.GRCH37, GenomeBuild.GRCH38):
        raise RuntimeError(f"{source_id}: genome build unknown - set genome_build in data.yaml")
    if existing and existing.genome_build == target.value and src.artifact_id in existing.parent_artifact_ids \
            and Path(existing.path + ".pgen").exists():
        return {"build": build.value, "action": f"already lifted to {target.value}", "artifact": existing.artifact_id}
    new_id = run_liftover(project, src.artifact_id, target.value, step_id=step_id, threads=threads)  # type: ignore[arg-type]
    # The variant table always describes the genotype used downstream (target build).
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        lifted = store.get(new_id)
        from efgpp.data.genotype.qc import placed_variants

        variants, _ = placed_variants(variant_table(resolve_fileset(Path(lifted.path), "pgen"), target.value))
        variants = VARIANT_SCHEMA.validate(variants)
        out = write_parquet(variants, project.artifact_dir(Origin.DERIVED, Modality.VARIANTS, source_id)
                            / f"{source_id}_variants.parquet", project.config.storage.compression)
        store.register_replacing(Artifact(
            artifact_name=f"{source_id}_variants", artifact_type="variant_table", modality=Modality.VARIANTS,
            origin=Origin.DERIVED, status=ArtifactStatus.READY, path=str(out), format="parquet",
            size=out.stat().st_size, checksum=str(checksum_paths([out])), feature_count=variants.height,
            genome_build=target.value, source_id=source_id, parent_artifact_ids=[new_id], tool="efgpp",
            temporal_type=TemporalType.STATIC,
        ))
    return {"build": build.value, "action": f"lifted to {target.value}", "artifact": new_id,
            "counts": lifted.metadata.get("liftover_counts")}


def detect_fileset_build(project: Project, fs: GenotypeFileset, declared: str | None = None) -> GenomeBuild:
    """Build of any fileset (e.g. a reference panel): declared value, else the same evidence
    used for genotypes (reference bases + pyliftover, coordinate bounds, FASTA)."""
    from efgpp.data.genotype.build import infer_build
    from efgpp.data.genotype.liftover import check_reference_bases, offline
    from efgpp.data.genotype.qc import variant_table

    if declared and GenomeBuild.normalize(declared) in (GenomeBuild.GRCH37, GenomeBuild.GRCH38):
        return GenomeBuild.normalize(declared)
    variants = variant_table(fs, None, limit=200_000)
    check = None
    if not offline() and project.config.defaults.online_build_check:
        try:
            check = check_reference_bases(variants, project.resource_root)
        except Exception:  # noqa: BLE001
            check = None
    result = infer_build(declared="auto", variants=variants, reference_check=check)
    if not result.confident:
        raise RuntimeError(f"cannot determine the genome build of {fs.prefix}; declare it "
                           "(resources.yaml: genome_build: GRCh37|GRCh38)")
    return result.build


def ensure_target_build(project: Project, fs: GenotypeFileset, declared: str | None, cache_prefix: Path, *,
                        step_id: str | None = None, threads: int = 1) -> GenotypeFileset:
    """Return `fs` if it is already in the target build, else a lifted copy (created once)."""
    target = target_build(project)
    build = detect_fileset_build(project, fs, declared)
    if build == target:
        return fs
    lifted_prefix = cache_prefix.with_name(f"{cache_prefix.name}_{target.value}")
    if lifted_prefix.with_name(lifted_prefix.name + ".pgen").exists():
        return resolve_fileset(lifted_prefix, "pgen")
    lifted, _counts, _report, _rec = liftover_fileset(project, fs, build, target, lifted_prefix,
                                                     step_id=step_id, threads=threads)
    return lifted
