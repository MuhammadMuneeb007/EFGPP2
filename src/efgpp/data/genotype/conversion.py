"""Genotype format conversion with PLINK 2. Outputs are new DERIVED artifacts; the
original fileset is never replaced."""

from __future__ import annotations

from pathlib import Path

from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.formats import count_variants, resolve_fileset
from efgpp.data.provenance import run_tool
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths

EXPORTS = {
    "pgen": ["--make-pgen"],
    "bed": ["--make-bed"],
    "vcf": ["--export", "vcf", "bgz", "id-paste=iid"],
    "bgen": ["--export", "bgen-1.2", "bits=8"],
}


def convert(
    project: Project,
    artifact_id: str,
    target: str,
    *,
    step_id: str | None = None,
    threads: int = 1,
    out_dir: Path | None = None,
) -> str:
    if target not in EXPORTS:
        raise ValueError(f"unsupported target format {target!r}; choose from {sorted(EXPORTS)}")
    with Registry.open(project) as reg:
        src = ArtifactStore(reg).get(artifact_id)
    fs = resolve_fileset(Path(src.path), src.format)
    out_dir = out_dir or project.artifact_dir(Origin.DERIVED, Modality.GENOTYPE_QC, src.source_id or "converted")
    prefix = out_dir / f"{src.artifact_name}.{target}"
    rec = run_tool(project, "plink2", [*fs.plink_input_args(), "--threads", str(threads), *EXPORTS[target],
                                       "--out", str(prefix)], step_id=step_id, inputs=[artifact_id])
    produced = resolve_fileset(prefix if target in ("pgen", "bed") else Path(str(prefix) + (".vcf.gz" if target == "vcf" else ".bgen")), target)
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=f"{src.artifact_name}_{target}", artifact_type="converted_genotype",
            modality=src.modality, origin=src.origin if src.origin != Origin.OBSERVED else Origin.DERIVED,
            status=ArtifactStatus.READY, path=str(produced.prefix if target in ("pgen", "bed") else produced.main),
            format=target, size=sum(m.stat().st_size for m in produced.members),
            checksum=str(checksum_paths(produced.members)), participant_count=src.participant_count,
            feature_count=count_variants(produced) or src.feature_count, genome_build=src.genome_build,
            temporal_type=TemporalType.STATIC, source_id=src.source_id, parent_artifact_ids=[artifact_id],
            tool="plink2", tool_version=rec.tool_version, command=" ".join(rec.command),
        ))
    return art.artifact_id  # type: ignore[return-value]
