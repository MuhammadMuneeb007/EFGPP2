"""Migration, CLI, report, snapshots, resources, phenopackets and the phenotype-agnostic rule."""

from __future__ import annotations

import json
import re
from pathlib import Path

import httpx
import polars as pl
import pytest
import yaml

from efgpp.constants import StorageMode
from efgpp.data import migrate
from efgpp.data.references.base import FetchedResource
from efgpp.data.registry import Registry
from efgpp.data.simulate import simulate_genotypes
from efgpp.data.snapshots import DataSnapshot, SnapshotExistsError, freeze
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths
from efgpp.workflow.executor import prepare
from tests.helpers import add_covariates, add_genotype, add_phenotypes

pytestmark = pytest.mark.integration
SRC = Path(__file__).resolve().parents[2] / "src" / "efgpp"


def test_no_phenotype_or_disease_is_hard_coded() -> None:
    banned = re.compile(r"\b(migraine|depression|asthma|cancer|diabetes|schizophrenia|alzheimer)\b", re.IGNORECASE)
    offenders = [f"{f.relative_to(SRC)}:{i}" for f in SRC.rglob("*.py")
                 for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1) if banned.search(line)]
    assert offenders == []


def _legacy_project(root: Path) -> Path:
    """Mimic the legacy EFGPP layout: <trait>/<trait>/{.bed,.bim,.fam,.cov,.height,.gz}."""
    for sub in ("traitx/traitx", "traitx/traitx_5"):
        d = root / sub
        cohort = simulate_genotypes(30, 50, 1, 11, d / "traitx")
        pl.DataFrame({"FID": cohort.iids, "IID": cohort.iids, "COV1": ["1.0"] * 30}).write_csv(d / "traitx.cov", separator=" ")
        pl.DataFrame({"IID": cohort.iids, "FID": cohort.iids, "Height": ["1", "2"] * 15}).write_csv(d / "traitx.height", separator=" ")
    import gzip

    with gzip.open(root / "traitx" / "traitx.gz", "wt") as fh:
        fh.write("CHR\tBP\tSNP\tA1\tA2\tN\tSE\tP\tOR\n1\t100\trs1\tA\tG\t1000\t0.1\t0.5\t1.01\n")
    (root / "Annotations.tsv").write_text("SNP\tscore\nrs1\t0.5\n")
    (root / "traitx" / "traitx" / "prs.profile").write_text("FID IID SCORE\n")
    (root / "notes.docx").write_bytes(b"x")
    return root


def test_migration_dry_run_and_apply(project: Project, tmp_path: Path) -> None:
    legacy = _legacy_project(tmp_path / "legacy")
    before = {f: checksum_paths([f]) for f in legacy.rglob("*") if f.is_file()}
    items = migrate.scan(legacy)
    kinds = [i.kind for i in items]
    assert kinds.count("genotype") == 2 and kinds.count("phenotype") == 2 and kinds.count("covariates") == 2
    assert {"gwas", "annotations", "prs", "unknown"} <= set(kinds)
    names = sorted(i.details["name"] for i in items if i.kind == "phenotype")
    assert names == ["traitx", "traitx_5"] and all(i.details["type"] == "binary" for i in items if i.kind == "phenotype")
    plan = migrate.write_plan(project, legacy, "legacy-efgpp", StorageMode.REFERENCE, items)
    assert yaml.safe_load(plan.read_text())["summary"]["genotype"] == 2
    done = migrate.apply(project, legacy, "legacy-efgpp", StorageMode.REFERENCE)
    project = Project.load(project.root)
    assert len(project.data.observed.genotype) == 2 and len(project.data.observed.phenotypes) == 2
    assert done["gwas"] and done["prs"]
    outcomes = {o.step_id: o.status for o in prepare(project, engine="builtin")}
    assert outcomes["qc.PH001"] == "done" and outcomes["qc.COV001"] == "done"
    after = {f: checksum_paths([f]) for f in legacy.rglob("*") if f.is_file()}
    assert before == after  # the source is never modified


def test_cli_end_to_end(tmp_path: Path, cohort_dir: Path, run_cli, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    root = tmp_path / "cli_proj"
    run_cli("init", str(root))
    monkeypatch.chdir(root)
    run_cli("data", "add", "genotype", "--path", str(cohort_dir / "genotype" / "cohort"), "--format", "bed",
            "--build", "GRCh38", "--mode", "reference")
    run_cli("phenotype", "add", "--name", "trait_a", "--path", str(cohort_dir / "phenotypes.csv"),
            "--id-column", "IID", "--value-column", "trait_a", "--type", "binary")
    run_cli("phenotype", "add", "--name", "trait_b", "--path", str(cohort_dir / "phenotypes.csv"),
            "--id-column", "IID", "--value-column", "trait_b", "--type", "continuous")
    run_cli("data", "add", "covariates", "--path", str(cohort_dir / "covariates.csv"), "--id-column", "IID",
            "--columns", "age", "sex", "bmi")
    data = yaml.safe_load((root / "data.yaml").read_text())
    assert data["observed"]["covariates"][0]["variables"] == ["age", "sex", "bmi"]
    assert [p["id"] for p in data["observed"]["phenotypes"]] == ["PH001", "PH002"]
    # Small phenotype tables are copied into the project (auto policy).
    assert any((root / "data" / "observed" / "phenotype" / "PH001").iterdir())
    out = run_cli("data", "inspect").output
    assert "GENO001" in out and "PH002" in out
    assert "treated as categorical" in run_cli("data", "validate").output  # a warning, not an error
    run_cli("data", "plan")
    run_cli("data", "prepare", "--engine", "builtin", expect=None)
    out = run_cli("data", "availability").output
    assert "OBSERVED" in out and "PH001" in out
    run_cli("data", "report")
    assert (root / "reports" / "data" / "index.html").exists()
    run_cli("data", "freeze", "--name", "cohort_data_v1")
    run_cli("data", "freeze", "--name", "cohort_data_v1", expect=1)
    run_cli("export", "slurm")
    assert (root / "hpc" / "submit_all.sh").exists()
    run_cli("phenotype", "export", "--format", "phenopackets")
    packets = list((root / "reports" / "data" / "phenopackets").glob("*.json"))
    assert len(packets) == 120 and "phenotypicFeatures" in json.loads(packets[0].read_text())
    assert "EFGPP DATA SYSTEM" in run_cli("doctor").output
    assert (root / "efgpp.lock.yaml").exists()


def test_snapshot_immutable_and_detects_changes(project: Project, cohort_dir: Path, tmp_path: Path) -> None:
    import shutil

    pheno = tmp_path / "ph.csv"
    shutil.copy(cohort_dir / "phenotypes.csv", pheno)
    p = add_phenotypes(project, cohort_dir, ["trait_a"])
    p.data.observed.phenotypes[0].path = str(pheno)
    p.data.observed.phenotypes[0].mode = StorageMode.REFERENCE
    p.save_data_config()
    p = Project.load(p.root)
    prepare(p, engine="builtin")
    freeze(p, "v1")
    with pytest.raises(SnapshotExistsError):
        freeze(p, "v1")
    snap = DataSnapshot.load(p, "v1")
    assert snap.verify() == {}
    pheno.write_text(pheno.read_text() + "EXTRA,1,1,x,low,1\n")
    assert snap.verify()  # the referenced source changed on disk
    with pytest.raises(RuntimeError, match="changed on disk"):
        freeze(p, "v2")


def test_resource_install_versioned_and_protected(project: Project, monkeypatch, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    from efgpp.data.references import PROVIDERS
    from efgpp.resources import manager

    class FakeClinVar(PROVIDERS["clinvar"]):  # type: ignore[misc, valid-type]
        def fetch(self, *, build=None, staging):  # type: ignore[no-untyped-def]
            f = staging / "clinvar.vcf.gz"
            f.write_bytes(b"fake")
            return FetchedResource("clinvar", "20260901_GRCh38", [f], f, "test", "public", genome_build="GRCh38")

    monkeypatch.setitem(PROVIDERS, "clinvar", FakeClinVar)
    res = manager.install(project, "clinvar")
    assert res.status == "installed" and res.path == project.resource_root / "clinvar" / "20260901_GRCh38" / "clinvar.vcf.gz"
    assert manager.install(project, "clinvar").status == "already-installed"
    assert (res.path.parent / "manifest.yaml").exists()
    freeze(project, "with_resource")
    with pytest.raises(manager.ResourceProtectedError):
        manager.install(project, "clinvar", force=True)
    rows = {r["name"]: r for r in manager.list_resources(project)}
    assert rows["clinvar"]["installed_version"] == "20260901_GRCh38"


def test_query_providers_with_mocked_http(project: Project) -> None:
    from efgpp.data.references.opentargets import OpenTargetsProvider
    from efgpp.data.references.pgs_catalog import PGSCatalogProvider

    def handler(request: httpx.Request) -> httpx.Response:
        if "pgscatalog" in request.url.host:
            return httpx.Response(200, json={"id": "PGS000001", "ftp_scoring_file": "https://x/PGS000001.txt.gz",
                                             "ftp_harmonized_scoring_files": {"GRCh38": {"positions": "https://x/h38.txt.gz"}}})
        return httpx.Response(200, json={"data": {"meta": {"apiVersion": {"x": 1}, "dataVersion": {"year": "26", "month": "6"}}}})

    client = httpx.Client(transport=httpx.MockTransport(handler))
    pgs = PGSCatalogProvider(project, client)
    meta = pgs.query("PGS000001")
    assert pgs.scoring_file_url(meta, "GRCh38") == "https://x/h38.txt.gz"
    assert pgs.scoring_file_url(meta, "GRCh37") == "https://x/PGS000001.txt.gz"
    assert OpenTargetsProvider(project, client).version() == "26.06"


def test_report_has_all_sections(project: Project, cohort_dir: Path) -> None:
    from efgpp.reporting.data_report import build_report

    p = add_covariates(add_phenotypes(add_genotype(project, cohort_dir), cohort_dir, ["trait_a", "trait_b"]), cohort_dir)
    prepare(p, engine="builtin")
    html = build_report(p).read_text(encoding="utf-8")
    for n in range(1, 24):
        assert f">{n}. " in html, n
    assert "Modality intersections (UpSet)" in html and 'data-theme="dark"' in html


def test_registry_parquet_exports(project: Project, cohort_dir: Path) -> None:
    p = add_phenotypes(project, cohort_dir, ["trait_a"])
    prepare(p, engine="builtin")
    reg_dir = p.registry_path.parent
    for name in ("participants", "aliases", "events", "biospecimens", "assays", "phenotypes", "availability"):
        assert (reg_dir / f"{name}.parquet").exists(), name
    with Registry.open(p) as reg:
        assert reg.scalar("SELECT count(*) FROM participants") == 120
