"""GWAS summary statistics with GWASLab (reference knowledge, not participant data).

Each GWAS in data.yaml (`gwas:`) is processed by the `gwas.<id>` plan step:

    GWASLab (own conda environment software/envs/gwaslab, via `gwaslab-python`)
      Sumstats(...)  -> basic_check(remove, remove_dup)  -> infer_build()  -> liftover(to GRCh38,
      chain_path = resources/liftover/hg19ToHg38.over.chain.gz)  -> resources/gwas/<id>/<id>.GRCh38.parquet

The output is registered as a REFERENCE artifact (modality gwas_sumstats) in the target build.
Selecting variants or thresholds for a target phenotype stays in the Representation layer.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from efgpp.config.data import GwasSource
from efgpp.constants import ArtifactStatus, GenomeBuild, Origin
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.liftover import chain_file
from efgpp.data.provenance import run_tool
from efgpp.data.registry import Registry
from efgpp.data.validation import ValidationReport
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths

GWASLAB_BUILD = {GenomeBuild.GRCH37: "19", GenomeBuild.GRCH38: "38"}
JOB_SCRIPT = Path(__file__).with_name("gwas_gwaslab_job.py")
MODALITY = "gwas_sumstats"


def get_gwas(project: Project, gwas_id: str) -> GwasSource:
    for g in project.data.gwas:
        if g.id == gwas_id:
            return g
    raise KeyError(f"unknown GWAS {gwas_id!r}")


def validate_gwas(project: Project, g: GwasSource) -> ValidationReport:
    import gzip

    rep = ValidationReport(f"GWAS {g.id} ({g.trait})", g.id)
    path = project.resolve(g.path)
    if not path.exists():
        rep.fail(f"file not found: {path}")
        return rep
    rep.ok(f"file exists ({path.stat().st_size / 1e6:,.1f} MB)")
    opener = gzip.open if path.name.endswith((".gz", ".bgz")) else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[operator]
            header = fh.readline().replace(",", " ").split()
    except (OSError, EOFError) as exc:
        rep.fail(f"cannot read the header: {exc}")
        return rep
    rep.info(f"columns: {' '.join(header[:20])}")
    missing = [c for c in g.columns.values() if c not in header]
    if missing:
        rep.fail(f"mapped columns not in the file: {missing}")
    elif g.columns:
        rep.ok(f"{len(g.columns)} mapped columns present")
    target = project.config.defaults.target_build
    rep.info(f"GWASLab: basic_check, infer_build (build: {g.build}), liftover to {target}")
    return rep


def run_gwas(project: Project, gwas_id: str, *, step_id: str | None = None, threads: int = 1) -> dict[str, Any]:
    g = get_gwas(project, gwas_id)
    target = GenomeBuild.normalize(project.config.defaults.target_build)
    declared = GenomeBuild.normalize(g.build) if g.build != "auto" else None
    out_dir = project.resource_root / "gwas" / g.id
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / f"{g.id}.{target.value}.parquet"
    report_path = out_dir / f"{g.id}.gwaslab_report.json"
    chains = {}
    for source in (GenomeBuild.GRCH37, GenomeBuild.GRCH38):
        if source != target:
            chains[f"{GWASLAB_BUILD[source]}->{GWASLAB_BUILD[target]}"] = str(chain_file(project.resource_root, source, target))
    job = {
        "input": str(project.resolve(g.path)), "fmt": g.fmt,
        "build": GWASLAB_BUILD.get(declared, "99") if declared else "99",
        "columns": g.columns, "n": g.n, "target": GWASLAB_BUILD[target], "chains": chains,
        "output": str(output), "report": str(report_path),
    }
    job_path = out_dir / f"{g.id}.job.json"
    job_path.write_text(json.dumps(job, indent=2), encoding="utf-8")
    rec = run_tool(project, "gwaslab", [str(JOB_SCRIPT), str(job_path)], step_id=step_id)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    qc_copy = project.qc_dir("gwas") / f"{g.id}_gwaslab.json"
    qc_copy.write_text(json.dumps(report, indent=2), encoding="utf-8")
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=f"{g.id}_sumstats", artifact_type="gwas_sumstats", modality=MODALITY,
            origin=Origin.REFERENCE, status=ArtifactStatus.READY, path=str(output), format="parquet",
            size=output.stat().st_size, checksum=str(checksum_paths([output])), feature_count=report.get("rows_final"),
            genome_build=target.value, source_id=g.id, tool="gwaslab", tool_version=report.get("gwaslab_version"),
            resource_versions={"chain": ", ".join(Path(c).name for c in chains.values())} if report.get("lifted") else {},
            command=" ".join(rec.command),
            metadata={"trait": g.trait, "input": job["input"], "gwaslab_report": report, "report": str(report_path)},
        ))
    return {"artifact": art.artifact_id, "build_detected": report.get("build_detected"),
            "lifted": report.get("lifted"), "rows": report.get("rows_final")}
