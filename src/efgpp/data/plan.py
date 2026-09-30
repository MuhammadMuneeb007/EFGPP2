"""The data plan: which steps exist, which can run now, and why others cannot.

Steps are pure descriptions. `run_step` executes one of them and is what both the
built-in executor and every Snakemake rule call (`efgpp data step <id>`).
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from efgpp.constants import ArtifactStatus, Modality, Origin
from efgpp.project import Project
from efgpp.setup.tools import available


@dataclass
class Step:
    id: str
    kind: str
    group: str  # used for SLURM script grouping
    description: str
    needs: list[str] = field(default_factory=list)
    tools: list[str] = field(default_factory=list)
    env: str = "core"
    threads: int = 1
    mem_mb: int = 4000
    runtime_min: int = 60
    source_id: str | None = None
    params: dict[str, Any] = field(default_factory=dict)
    enabled: bool = True
    reason: str | None = None
    always_run: bool = False

    @property
    def marker(self) -> str:
        return f"work/steps/{self.id}.done"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _tool_gate(project: Project, step: Step) -> Step:
    for tool in step.tools:
        if not available(project, tool):
            step.enabled = False
            step.reason = f"{tool} not installed (run `efgpp setup data`)"
            break
    return step


def _resource_installed(project: Project, name: str) -> bool:
    from efgpp.data.registry import Registry

    with Registry.open(project) as reg:
        row = reg.one("SELECT local_path FROM resources WHERE name = ? ORDER BY download_date DESC LIMIT 1", [name])
    return bool(row and row["local_path"] and Path(row["local_path"]).exists())


def build_plan(project: Project) -> list[Step]:
    data = project.data
    cores = project.config.execution.local_cores
    steps: list[Step] = []
    steps.append(Step("register", "register", "prepare",
                      "register sources, participants, events and biospecimens (applies storage policy)",
                      always_run=True))
    steps.append(Step("validate", "validate", "prepare", "validate every source (collects all errors)",
                      needs=["register"], always_run=True))

    genotype_qc_steps: dict[str, str] = {}
    for modality, source in data.iter_sources():
        std = Step(f"standardize.{source.id}", "standardize", "prepare",
                   f"resolve participant IDs, timeline and standardize {modality.value} {source.id}",
                   needs=["validate"], source_id=source.id)
        steps.append(std)
        if modality == Modality.GENOTYPE:
            qc = _tool_gate(project, Step(
                f"genotype_qc.{source.id}", "genotype_qc", "genotype_qc",
                "sample/variant missingness, MAF, HWE, heterozygosity, LD pruning, KING kinship, sex check",
                needs=[std.id], tools=["plink2"], env="genetics", threads=cores, mem_mb=16000, runtime_min=240,
                source_id=source.id))
            steps.append(qc)
            genotype_qc_steps[source.id] = qc.id
            steps.append(_gated_after(qc, _tool_gate(project, Step(
                f"pca.{source.id}", "pca", "pca", "QC principal components (not modelling PCs)",
                needs=[qc.id], tools=["plink2"], env="genetics", threads=cores, mem_mb=16000,
                runtime_min=120, source_id=source.id))))
            if project.config.genotype.roh.enabled:
                steps.append(_gated_after(qc, _tool_gate(project, Step(
                    f"roh.{source.id}", "roh", "relatedness", "runs of homozygosity (PLINK 1.9)",
                    needs=[qc.id], tools=["plink2", "plink"], env="genetics", source_id=source.id))))
            from efgpp.data.genotype.ancestry import ancestry_configured

            if ancestry_configured(project):
                steps.append(_gated_after(qc, _tool_gate(project, Step(
                    f"ancestry.{source.id}", "ancestry", "pca", "ancestry by reference-panel projection",
                    needs=[qc.id], tools=["plink2"], env="genetics", threads=cores, source_id=source.id))))
            steps += _annotation_steps(project, source.id, qc if qc.enabled else std)
        else:
            steps.append(Step(f"qc.{source.id}", "qc", "prepare", f"profile and QC {modality.value} {source.id}",
                              needs=[std.id], source_id=source.id))

    for modality, cfg in data.predicted.items():
        if not cfg.enabled:
            continue
        gid = cfg.genotype_artifact or (data.observed.genotype[0].id if data.observed.genotype else None)
        for tissue in cfg.tissues:
            step = Step(f"predict.{modality.value}.{tissue}", "predict", "expression_prediction",
                        f"PrediXcan predicted {modality.value} ({tissue})",
                        needs=[genotype_qc_steps.get(gid, f"standardize.{gid}")] if gid else ["validate"],
                        tools=["plink2", "predixcan"], env="metaxcan", threads=1, mem_mb=8000,
                        runtime_min=240, source_id=gid, params={"modality": modality.value, "tissue": tissue})
            if gid is None:
                step.enabled, step.reason = False, "no genotype source for prediction"
            else:
                _tool_gate(project, step)
                from efgpp.data.predicted.metaxcan import model_path

                if step.enabled and not model_path(project, cfg, tissue).exists():
                    step.enabled, step.reason = False, f"PredictDB model for {tissue} not found"
                if gid in genotype_qc_steps:
                    qc_step = next(s for s in steps if s.id == genotype_qc_steps[gid])
                    _gated_after(qc_step, step)
            steps.append(step)

    enabled_ids = [s.id for s in steps if s.enabled and s.id not in ("register", "validate")]
    steps.append(Step("availability", "availability", "report", "participant x modality availability",
                      needs=["validate", *enabled_ids], always_run=True))
    steps.append(Step("report", "report", "report", "HTML data report", needs=["availability"], always_run=True))
    return steps


def _gated_after(upstream: Step, step: Step) -> Step:
    if not upstream.enabled and step.enabled:
        step.enabled, step.reason = False, f"needs {upstream.id} ({upstream.reason})"
    return step


def _annotation_steps(project: Project, source_id: str, upstream: Step) -> list[Step]:
    r = project.resources
    specs = [
        ("vep", r.vep.enabled, ["vep"], "annotation", "Ensembl VEP consequences", None),
        ("opencravat", r.opencravat.enabled, ["oc"], "annotation", "OpenCRAVAT annotation", None),
        ("clinvar", r.clinvar.enabled, [], "core", "ClinVar significance", "clinvar"),
        ("alphamissense", r.alphamissense.enabled, [], "core", "AlphaMissense scores", "alphamissense"),
        ("gnomad", r.gnomad.enabled, [], "core", "gnomAD allele frequencies", "gnomad"),
        ("dbsnp", r.dbsnp.enabled, [], "core", "dbSNP rsIDs", "dbsnp"),
        ("alphagenome", r.alphagenome.enabled, [], "core", "AlphaGenome regulatory scores", None),
    ]
    out = []
    for name, enabled, tools, env, desc, resource in specs:
        if not enabled:
            continue
        step = _tool_gate(project, Step(
            f"annotate.{name}.{source_id}", f"annotate_{name}", "annotation", desc, needs=[upstream.id],
            tools=tools, env=env, threads=project.config.execution.local_cores if name in ("vep", "opencravat") else 1,
            mem_mb=8000, runtime_min=240, source_id=source_id, params={"annotation": name}))
        cfg = getattr(r, name)
        if step.enabled and resource and not (cfg.path or _resource_installed(project, resource)):
            step.enabled, step.reason = False, f"{resource} not installed (efgpp resources install {resource})"
        if step.enabled and name == "vep" and r.vep.cache == "auto" and not any((project.resource_root / "vep").iterdir()):
            step.enabled, step.reason = False, "VEP cache missing (efgpp resources install vep)"
        if step.enabled and name == "alphagenome":
            from efgpp.data.annotation.alphagenome import api_key

            if r.alphagenome.mode == "atlas" and not r.alphagenome.atlas_path:
                step.enabled, step.reason = False, "AlphaGenome atlas_path not configured"
            elif r.alphagenome.mode == "api" and not api_key(project):
                step.enabled, step.reason = False, f"${r.alphagenome.api_key_env} not set (atlas mode needs no key)"
        out.append(_gated_after(upstream, step))
    return out


def select(steps: list[Step], kinds: set[str] | None = None, ids: set[str] | None = None) -> list[Step]:
    """Steps of the requested kinds plus everything they transitively need."""
    by_id = {s.id: s for s in steps}
    wanted = {s.id for s in steps if (kinds is None or s.kind in kinds or any(s.kind.startswith(k) for k in kinds))
              and (ids is None or s.id in ids)}
    frontier = list(wanted)
    while frontier:
        sid = frontier.pop()
        for n in by_id[sid].needs:
            if n in by_id and n not in wanted:
                wanted.add(n)
                frontier.append(n)
    return [s for s in steps if s.id in wanted]


def fingerprint(project: Project, step: Step) -> str:
    """Identity of a step's inputs: its definition, configuration and source checksums."""
    from efgpp.data.registry import Registry

    h = hashlib.sha256(json.dumps(step.to_dict(), sort_keys=True, default=str).encode())
    h.update(json.dumps(project.config_checksums(), sort_keys=True).encode())
    if step.source_id:
        with Registry.open(project) as reg:
            rows = reg.rows("SELECT checksum FROM artifacts WHERE source_id = ? AND artifact_type = 'source' "
                            "AND status <> 'SUPERSEDED'", [step.source_id])
        h.update("".join(sorted(str(r["checksum"]) for r in rows)).encode())
    return h.hexdigest()


def run_step(project: Project, step: Step, threads: int | None = None) -> dict[str, Any]:
    """Execute one step. Raises on failure."""
    from efgpp.data.manager import DataManager

    threads = threads or step.threads
    mgr = DataManager(project)
    kind = step.kind
    result: dict[str, Any] = {"step": step.id}
    if kind == "register":
        result["registered"] = [r.__dict__ for r in mgr.register_all()]
        _register_truth(project)
    elif kind == "validate":
        reports = mgr.validate()
        result["failed"] = [r.source_id for r in reports if not r.passed and r.source_id]
    elif kind == "standardize":
        art = mgr.standardize(step.source_id)  # type: ignore[arg-type]
        result["artifact"] = art.artifact_id if art else None
    elif kind == "qc":
        qc = mgr.qc(step.source_id)  # type: ignore[arg-type]
        result["status"] = qc.status.value if qc else None
    elif kind == "genotype_qc":
        from efgpp.data.genotype.qc import register_qc_outputs, run_genotype_qc

        out = run_genotype_qc(project, step.source_id, step_id=step.id, threads=threads)  # type: ignore[arg-type]
        result.update(register_qc_outputs(project, step.source_id, out, step.id))  # type: ignore[arg-type]
        result["status"] = out.status.value
    elif kind == "pca":
        from efgpp.data.genotype.pca import run_qc_pca

        result["artifacts"] = run_qc_pca(project, step.source_id, step_id=step.id, threads=threads)  # type: ignore[arg-type]
    elif kind == "roh":
        from efgpp.data.genotype.roh import run_roh

        result["artifact"] = run_roh(project, step.source_id, step_id=step.id, threads=threads)  # type: ignore[arg-type]
    elif kind == "ancestry":
        from efgpp.data.genotype.ancestry import run_ancestry

        result["artifact"] = run_ancestry(project, step.source_id, step_id=step.id, threads=threads)  # type: ignore[arg-type]
    elif kind.startswith("annotate_"):
        result["artifact"] = _annotate(project, step, threads)
    elif kind == "predict":
        from efgpp.data.predicted.metaxcan import run_prediction

        result["artifact"] = run_prediction(project, Modality(step.params["modality"]), step.params["tissue"],
                                            step_id=step.id, threads=threads)
    elif kind == "availability":
        from efgpp.data import availability

        av = availability.build(project)
        result["participants"] = av.n_participants
    elif kind == "report":
        from efgpp.reporting.data_report import build_report

        result["report"] = str(build_report(project))
    else:
        raise ValueError(f"unknown step kind {kind!r}")
    return result


def _annotate(project: Project, step: Step, threads: int) -> str:
    name = step.params["annotation"]
    sid = step.source_id
    if name == "vep":
        from efgpp.data.annotation.vep import run_vep as fn
    elif name == "opencravat":
        from efgpp.data.annotation.opencravat import run_opencravat as fn
    elif name == "clinvar":
        from efgpp.data.annotation.clinvar import run_clinvar as fn
    elif name == "alphamissense":
        from efgpp.data.annotation.alphamissense import run_alphamissense as fn
    elif name == "gnomad":
        from efgpp.data.annotation.gnomad import run_gnomad as fn
    elif name == "dbsnp":
        from efgpp.data.annotation.gnomad import run_dbsnp as fn
    elif name == "alphagenome":
        from efgpp.data.annotation.alphagenome import run_alphagenome as fn
    else:
        raise ValueError(name)
    return fn(project, sid, step_id=step.id, threads=threads)  # type: ignore[arg-type]


def _register_truth(project: Project) -> None:
    truth = project.data_root / Origin.SIMULATED.value / "truth" / "truth.yaml"
    if not truth.exists():
        return
    from efgpp.data.artifacts import Artifact, ArtifactStore
    from efgpp.data.registry import Registry
    from efgpp.resources.checksums import checksum_paths

    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        checksum = str(checksum_paths([truth]))
        current = store.latest(source_id="SIMULATION", artifact_type="simulation_truth")
        if current and current.checksum == checksum:
            return
        store.register_replacing(Artifact(
            artifact_name="simulation_truth", artifact_type="simulation_truth", modality=Modality.TRUTH,
            origin=Origin.SIMULATED, status=ArtifactStatus.READY, path=str(truth), format="yaml",
            size=truth.stat().st_size, checksum=checksum, source_id="SIMULATION", tool="efgpp",
            metadata={"note": "ground truth; never use as model input"},
        ))


def write_plan(project: Project, steps: list[Step]) -> Path:
    out = project.path("workflow", "plan.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {"project": project.config.project.name, "root": str(project.root),
               "steps": [s.to_dict() for s in steps]}
    tmp = out.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(tmp, out)
    return out


def load_plan(project: Project) -> list[Step]:
    data = json.loads(project.path("workflow", "plan.json").read_text(encoding="utf-8"))
    return [Step(**s) for s in data["steps"]]
