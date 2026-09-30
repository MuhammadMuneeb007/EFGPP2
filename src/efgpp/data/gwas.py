"""GWAS summary statistics with GWASLab (reference knowledge, not participant data).

Each GWAS in data.yaml (`gwas:`) is processed by the `gwas.<id>` plan step:

    GWASLab (own conda environment software/envs/gwaslab, via `gwaslab-python`)
      Sumstats(...)  -> basic_check(remove, remove_dup)  -> infer_build()  -> liftover(to GRCh38,
      chain_path = resources/liftover/hg19ToHg38.over.chain.gz)
      -> resources/gwas/<id>/<id>.GRCh38.parquet          (original column names)
         resources/gwas/<id>/<id>.GRCh38.gwaslab.parquet  (GWASLab standard names)

Study details (ancestry, cases, controls, sample size, the phenotypes the GWAS is for) are kept
with the artifact; `efgpp data gwas list --phenotype <name>` finds the GWAS for a phenotype.

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


def _input_signature(path: Path, g: GwasSource) -> str:
    st = path.stat()
    return f"{st.st_size}:{int(st.st_mtime)}:{g.checksum()}"


def run_gwas(project: Project, gwas_id: str, *, step_id: str | None = None, threads: int = 1,
             force: bool = False) -> dict[str, Any]:
    g = get_gwas(project, gwas_id)
    target = GenomeBuild.normalize(project.config.defaults.target_build)
    declared = GenomeBuild.normalize(g.build) if g.build != "auto" else None
    input_path = project.resolve(g.path)
    signature = _input_signature(input_path, g)
    out_dir = project.resource_root / "gwas" / g.id
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / f"{g.id}.{target.value}.parquet"
    output_gwaslab = out_dir / f"{g.id}.{target.value}.gwaslab.parquet"
    report_path = out_dir / f"{g.id}.gwaslab_report.json"

    with Registry.open(project) as reg:
        current = ArtifactStore(reg).latest(source_id=g.id, artifact_type="gwas_sumstats")
    if (not force and current is not None and current.metadata.get("input_signature") == signature
            and output.exists()):
        return {"artifact": current.artifact_id, "path": str(output), "rows": current.feature_count,
                "build_detected": (current.metadata.get("gwaslab_report") or {}).get("build_detected"),
                "lifted": (current.metadata.get("gwaslab_report") or {}).get("lifted"), "status": "up to date"}

    chains = {}
    for source in (GenomeBuild.GRCH37, GenomeBuild.GRCH38):
        if source != target:
            chains[f"{GWASLAB_BUILD[source]}->{GWASLAB_BUILD[target]}"] = str(chain_file(project.resource_root, source, target))
    header = _header(input_path)
    job = {
        "input": str(input_path), "fmt": g.fmt,
        "build": GWASLAB_BUILD.get(declared, "99") if declared else "99",
        "columns": g.columns, "other": [c for c in header if c not in g.columns.values()],
        "n": g.n, "ncase": g.n_cases, "ncontrol": g.n_controls,
        "target": GWASLAB_BUILD[target], "chains": chains,
        "output": str(output), "output_gwaslab": str(output_gwaslab), "report": str(report_path),
    }
    job_path = out_dir / f"{g.id}.job.json"
    job_path.write_text(json.dumps(job, indent=2), encoding="utf-8")
    rec = run_tool(project, "gwaslab", [str(JOB_SCRIPT), str(job_path)], step_id=step_id)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    qc_copy = project.qc_dir("gwas") / f"{g.id}_gwaslab.json"
    qc_copy.write_text(json.dumps(report, indent=2), encoding="utf-8")
    files = [output, *([output_gwaslab] if output_gwaslab.exists() else [])]
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=f"{g.id}_sumstats", artifact_type="gwas_sumstats", modality=MODALITY,
            origin=Origin.REFERENCE, status=ArtifactStatus.READY, path=str(output), format="parquet",
            size=sum(f.stat().st_size for f in files), checksum=str(checksum_paths(files)),
            feature_count=report.get("rows_final"), genome_build=target.value, source_id=g.id, tool="gwaslab",
            tool_version=report.get("gwaslab_version"),
            resource_versions={"chain": ", ".join(Path(c).name for c in chains.values())} if report.get("lifted") else {},
            command=" ".join(rec.command),
            metadata={"trait": g.trait, "phenotypes": g.phenotypes or [g.trait], "ancestry": g.ancestry,
                      "n_cases": g.n_cases, "n_controls": g.n_controls, "n": g.n, "study": g.study,
                      "input": job["input"], "input_signature": signature, "gwaslab_report": report,
                      "report": str(report_path), "gwaslab_standard_file": str(output_gwaslab),
                      "members": [str(f) for f in files]},
        ))
    return {"artifact": art.artifact_id, "path": str(output), "build_detected": report.get("build_detected"),
            "lifted": report.get("lifted"), "rows": report.get("rows_final")}


def _header(path: Path) -> list[str]:
    import gzip

    opener = gzip.open if path.name.endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[operator]
        return fh.readline().replace(",", " ").split()


def gwas_for_phenotype(project: Project, phenotype: str) -> list[GwasSource]:
    """GWAS configured for a project phenotype (by name or id)."""
    names = {phenotype}
    for p in project.data.observed.phenotypes:
        if phenotype in (p.id, p.name):
            names |= {p.id, p.name}
    return [g for g in project.data.gwas if names & set(g.phenotypes or [g.trait])]
