"""Critical acceptance tests A-J (section 85 of the specification)."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import polars as pl
import pytest

from efgpp.constants import ArtifactStatus, Modality, Origin
from efgpp.data import availability
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.plan import build_plan
from efgpp.data.registry import Registry
from efgpp.data.snapshots import DataSnapshot, freeze
from efgpp.project import Project
from efgpp.workflow.executor import prepare
from tests.helpers import add_covariates, add_expression, add_genotype, add_phenotypes

pytestmark = pytest.mark.integration


def _run(project: Project, **kw) -> dict[str, str]:  # type: ignore[no-untyped-def]
    outcomes = prepare(project, engine="builtin", **kw)
    return {o.step_id: o.status for o in outcomes}


def _store(project: Project):  # type: ignore[no-untyped-def]
    return Registry.open(project)


def test_A_genotype_binary_phenotype(project: Project, cohort_dir: Path) -> None:
    p = add_phenotypes(add_genotype(project, cohort_dir), cohort_dir, ["trait_a"])
    status = _run(p)
    assert status["standardize.GENO001"] == "done" and status["qc.PH001"] == "done"
    av = availability.build(p)
    assert av.count("GENO001", "PH001") == 120
    with Registry.open(p) as reg:
        art = ArtifactStore(reg).latest(source_id="PH001", artifact_type="standardized")
        q = reg.one("SELECT summary FROM qc_runs WHERE artifact_id = ?", [art.artifact_id])  # type: ignore[union-attr]
    assert art.status == ArtifactStatus.READY and q["summary"]["cases"] + q["summary"]["controls"] == 120  # type: ignore[index, union-attr]


def test_B_genotype_quantitative_phenotype(project: Project, cohort_dir: Path) -> None:
    p = add_phenotypes(add_genotype(project, cohort_dir), cohort_dir, ["trait_b"])
    _run(p)
    obs = pl.read_parquet(p.path("phenotypes", "PH001", "observations.parquet"))
    assert obs.height == 120 and obs.get_column("value_numeric").std() > 0
    assert p.path("phenotypes", "PH001", "definition.yaml").exists()


def test_C_genotype_five_phenotypes(project: Project, cohort_dir: Path) -> None:
    p = add_phenotypes(add_genotype(project, cohort_dir), cohort_dir,
                       ["trait_a", "trait_b", "trait_c", "trait_d", "trait_e"])
    status = _run(p)
    assert all(status[f"qc.PH00{i}"] == "done" for i in range(1, 6))
    with Registry.open(p) as reg:
        defs = reg.rows("SELECT phenotype_id, type FROM phenotype_definitions ORDER BY 1")
    assert [d["type"] for d in defs] == ["binary", "continuous", "multiclass", "ordinal", "continuous"]
    ordinal = pl.read_parquet(p.path("phenotypes", "PH004", "observations.parquet"))
    assert set(ordinal.get_column("value_numeric").drop_nulls().to_list()) <= {0.0, 1.0, 2.0}


def test_D_phenotype_and_expression_without_genotype(project: Project, cohort_dir: Path) -> None:
    p = add_expression(add_phenotypes(project, cohort_dir, ["trait_a"]), cohort_dir)
    status = _run(p)
    assert "genotype_qc.GENO001" not in status and status["qc.RNA001"] == "done"
    av = availability.build(p)
    assert av.count("PH001", "RNA001_OBSERVED") == 60


def test_E_genotype_phenotype_measured_expression(project: Project, cohort_dir: Path) -> None:
    p = add_expression(add_phenotypes(add_genotype(project, cohort_dir), cohort_dir, ["trait_a"]), cohort_dir)
    _run(p)
    with Registry.open(p) as reg:
        rna = ArtifactStore(reg).latest(source_id="RNA001", artifact_type="standardized")
    assert rna.origin == Origin.OBSERVED and rna.feature_count == 8 and rna.participant_count == 60  # type: ignore[union-attr]
    assert Path(rna.path).is_relative_to(p.data_root / "observed" / "expression")  # type: ignore[union-attr]
    inter = {i["modalities"]: i["participants"] for i in availability.build(p).key_intersections()}
    assert inter["GENO001 + RNA001_OBSERVED"] == 60


def test_F_genotype_phenotype_predicted_expression(project: Project, cohort_dir: Path, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    from efgpp.data.predicted import metaxcan
    from efgpp.data.provenance import RunRecord

    p = add_phenotypes(add_genotype(project, cohort_dir), cohort_dir, ["trait_a"])
    p.data.predicted.expression.enabled = True
    p.data.predicted.expression.tissues = ["Whole_Blood"]
    p.save_data_config()
    p = Project.load(p.root)
    model = p.resource_root / "predictdb" / "mashr_Whole_Blood.db"
    con = sqlite3.connect(model)
    con.execute("CREATE TABLE extra (gene TEXT)")
    con.commit()
    con.close()
    _run(p)
    step = next(s for s in build_plan(p) if s.id == "predict.expression.Whole_Blood")
    assert not step.enabled and "not installed" in (step.reason or "")  # tools absent: reported, not run

    samples = pl.read_csv(cohort_dir / "genotype" / "cohort.fam", separator="\t", has_header=False).get_column("column_2").cast(pl.Utf8)

    def fake_convert(project: Project, artifact_id: str, target: str, **_kw):  # type: ignore[no-untyped-def]
        with Registry.open(project) as reg:
            return ArtifactStore(reg).register(Artifact(
                artifact_name="geno_vcf", artifact_type="converted_genotype", modality=Modality.GENOTYPE,
                origin=Origin.DERIVED, path=str(project.work_root / "g.vcf.gz"), format="vcf", source_id="GENO001",
                parent_artifact_ids=[artifact_id], status=ArtifactStatus.READY)).artifact_id

    def fake_run_tool(project: Project, tool: str, args: list[str], **_kw) -> RunRecord:  # type: ignore[no-untyped-def]
        pred = Path(args[args.index("--prediction_output") + 1])
        pl.DataFrame({"FID": samples, "IID": samples, "ENSG01": [0.1] * len(samples),
                      "ENSG02": [0.2] * len(samples)}).write_csv(pred, separator="\t")
        Path(args[args.index("--prediction_summary_output") + 1]).write_text("gene\tn_snps_used\nENSG01\t3\n")
        return RunRecord("RUN999999", None, tool, "test", [tool, *args], exit_code=0, status="completed")

    monkeypatch.setattr(metaxcan, "convert", fake_convert)
    monkeypatch.setattr(metaxcan, "run_tool", fake_run_tool)
    art_id = metaxcan.run_prediction(p, Modality.EXPRESSION, "Whole_Blood")
    with Registry.open(p) as reg:
        art = ArtifactStore(reg).get(art_id)
    assert art.origin == Origin.PREDICTED and art.tissue == "Whole_Blood"
    assert art.temporal_type is not None and art.temporal_type.value == "genetically_predicted_static"
    assert Path(art.path).is_relative_to(p.data_root / "predicted" / "expression")
    av = availability.build(p)
    col = av.columns.filter(pl.col("column") == "PRED_EXPRESSION_Whole_Blood").to_dicts()[0]
    assert col["origin"] == "predicted" and col["participants"] == 120
    assert not any(c.endswith("_OBSERVED") for c in av.column_names())  # never reported as measured


def test_G_longitudinal_expression(project: Project, cohort_dir: Path) -> None:
    from efgpp.data.timeline import timeline_summary

    p = add_expression(project, cohort_dir, longitudinal=True)
    status = _run(p)
    assert status["qc.RNA002"] == "done"
    with Registry.open(p) as reg:
        tl = timeline_summary(reg)
        n_specimens = reg.scalar("SELECT count(*) FROM biospecimens")
        n_assays = reg.scalar("SELECT count(*) FROM assays WHERE source_id = 'RNA002'")
    assert tl["event_order"] == ["baseline", "month_6", "month_12"]
    per_event = {r["event_id"]: r["participants"] for r in tl["participants_per_event"]}
    assert per_event == {"baseline": 40, "month_6": 36, "month_12": 30}
    assert tl["attrition"][-1]["retained_from_first"] == 30
    assert tl["repeated_measurements"][0]["participants_with_repeats"] == 38  # 2 have baseline only
    assert n_specimens == n_assays == 106


def test_H_external_genotype_referenced_not_copied(project: Project, cohort_dir: Path) -> None:
    bed = cohort_dir / "genotype" / "cohort.bed"
    mtime = bed.stat().st_mtime
    p = add_genotype(project, cohort_dir)
    _run(p)
    with Registry.open(p) as reg:
        src = ArtifactStore(reg).latest(source_id="GENO001", artifact_type="source")
    assert src.storage_mode == "reference"  # type: ignore[union-attr]
    assert Path(src.path) == (cohort_dir / "genotype" / "cohort").resolve()  # type: ignore[union-attr]
    assert src.checksum and src.checksum.startswith("sha256-manifest:")  # type: ignore[union-attr]
    assert not any((p.data_root / "observed" / "genotype").iterdir())
    assert bed.stat().st_mtime == mtime


@pytest.mark.plink2
def test_I_local_execution_end_to_end(with_plink2: Project, cohort_dir: Path) -> None:
    p = add_covariates(add_phenotypes(add_genotype(with_plink2, cohort_dir), cohort_dir, ["trait_a", "trait_b"]), cohort_dir)
    status = _run(p)
    for sid in ("genotype_qc.GENO001", "pca.GENO001", "report", "availability"):
        assert status[sid] == "done", status
    with Registry.open(p) as reg:
        store = ArtifactStore(reg)
        qc = store.latest(source_id="GENO001", artifact_type="qc_genotype")
        pca = store.latest(source_id="GENO001", artifact_type="qc_pca")
        runs = reg.scalar("SELECT count(*) FROM tool_runs WHERE tool = 'plink2' AND exit_code = 0")
    assert qc is not None and qc.origin == Origin.DERIVED and qc.format == "pgen"
    assert qc.metadata["qc_metrics"]["samples_total"] == 120
    assert pca is not None and "PC1" in pl.read_parquet(pca.path).columns
    assert runs >= 4 and any(p.logs_dir.glob("RUN*.jsonl")) and any(p.logs_dir.glob("RUN*.log"))
    html = (p.report_root / "data" / "index.html").read_text(encoding="utf-8")
    assert "9. Genotype QC" in html and "23. Provenance" in html
    # Re-running is idempotent: nothing changed, so derived steps are up to date.
    again = _run(p)
    assert again["genotype_qc.GENO001"] == "up-to-date" and again["pca.GENO001"] == "up-to-date"
    path = freeze(p, "data_v1")
    snap = DataSnapshot.load(p, "data_v1")
    sel = snap.select("trait_a", require=["genotype_qc"], include=["covariates", "qc_pca"])
    assert len(sel.participants) == qc.participant_count
    assert sel.artifacts["genotype_qc"][0]["artifact_id"] == qc.artifact_id
    assert snap.verify() == {} and path.exists()


def test_J_slurm_plan_generation(project: Project, cohort_dir: Path) -> None:
    from efgpp.workflow.slurm import export_slurm
    from efgpp.workflow.snakemake import write_snakefile

    p = add_phenotypes(add_genotype(project, cohort_dir), cohort_dir, ["trait_a"])
    p.config.execution.tools = {"plink2": str(Path(__file__))}  # pretend plink2 exists: planning only
    p.config.execution.hpc.partition = "compute"
    p.save_project_config()
    p = Project.load(p.root)
    files = {f.name: f for f in export_slurm(p)}
    assert {"prepare.sbatch", "genotype_qc.sbatch", "pca.sbatch", "report.sbatch", "submit_all.sh"} <= set(files)
    gq = files["genotype_qc.sbatch"].read_text()
    assert "#SBATCH --partition=compute" in gq and "data step genotype_qc.GENO001" in gq
    submit = files["submit_all.sh"].read_text()
    assert "--dependency=afterok:${JOB_PREPARE}" in submit and "--dependency=afterok:" in submit.split("JOB_PCA=")[1]
    snakefile = write_snakefile(p, build_plan(p))
    rules = (p.path("workflow", "rules", "genotype_qc.smk")).read_text()
    assert "rule genotype_qc_GENO001" in rules and "mem_mb=16000" in rules
    assert 'ancient("work/steps/validate.done")' in p.path("workflow", "rules", "validate.smk").read_text()
    assert "include: \"rules/pca.smk\"" in snakefile.read_text()
