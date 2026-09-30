"""Target build GRCh38 everywhere: pyliftover build check, genotype harmonization, GWASLab step."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import polars as pl
import pytest

from efgpp.constants import GenomeBuild
from efgpp.data.genotype.build import evidence_from_reference_check, infer_build
from efgpp.data.genotype.liftover import ReferenceCheck, check_reference_bases
from efgpp.project import Project


def test_reference_check_decides_the_build() -> None:
    strong37 = ReferenceCheck(rate37=1.0, rate38=0.5, lifted_rate=1.0, tested=50, lifted_direction="GRCh37->GRCh38")
    assert evidence_from_reference_check(strong37)[0].supports == GenomeBuild.GRCH37
    result = infer_build(declared="auto", reference_check=strong37)
    assert result.build == GenomeBuild.GRCH37 and result.confident
    # A declared build contradicted by the reference bases stops and asks the user.
    assert infer_build(declared="GRCh38", reference_check=strong37).requires_user
    unclear = ReferenceCheck(rate37=0.7, rate38=0.6, lifted_rate=None, tested=50, lifted_direction=None)
    assert evidence_from_reference_check(unclear)[0].supports is None
    # pyliftover disagreeing with the reference invalidates the call.
    bad_lift = ReferenceCheck(rate37=1.0, rate38=0.5, lifted_rate=0.4, tested=50, lifted_direction="GRCh37->GRCh38")
    assert evidence_from_reference_check(bad_lift)[0].supports is None


def test_check_reference_bases_with_fake_reference(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import efgpp.data.genotype.liftover as lo

    sites = [(str(1 + i % 20), 1_000_000 + 1000 * i) for i in range(40)]
    ref37 = dict.fromkeys(sites, "A")
    ref38 = {s: ("A" if i % 2 else "C") for i, s in enumerate(sites)}
    lifted38 = {(c, p + 5): "A" for c, p in sites}
    variants = pl.DataFrame({"chromosome": [c for c, _ in sites], "position": [p for _, p in sites],
                             "variant_id": [f"v{i}" for i in range(40)], "reference": ["A"] * 40, "alternate": ["G"] * 40})

    def fetch(build, wanted):  # type: ignore[no-untyped-def]
        table = ref37 if build == GenomeBuild.GRCH37 else {**ref38, **lifted38}
        return {s: table[s] for s in wanted if s in table}

    class FakeLifter:
        def __init__(self, _path):  # type: ignore[no-untyped-def]
            pass

        def lift(self, chrom, pos):  # type: ignore[no-untyped-def]
            return chrom, pos + 5, "+"

    monkeypatch.setattr(lo, "chain_file", lambda root, s, t: tmp_path / "chain")
    check = check_reference_bases(variants, tmp_path, fetch=fetch, lifter_factory=FakeLifter)
    assert check is not None and check.rate37 == 1.0 and check.rate38 == 0.5 and check.lifted_rate == 1.0
    assert infer_build(declared="auto", reference_check=check).build == GenomeBuild.GRCH37


def test_target_build_is_grch38_and_resources_follow(project: Project) -> None:
    from efgpp.resources.manager import install

    assert project.config.defaults.target_build == "GRCh38"
    with pytest.raises(ValueError, match="target build is GRCh38"):
        install(project, "clinvar", build="GRCh37")


def _toy_chain(path: Path, offset: int) -> Path:
    """GRCh37 -> GRCh38 chain shifting chr1 and chr2 by `offset` (enough for the simulated cohort)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(f"chain 1000 chr{c} 200000000 + 0 100000000 chr{c} 200000000 + {offset} {100000000 + offset} {c}\n"
                   f"100000000\n\n" for c in (1, 2))
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


@pytest.mark.plink2
def test_grch37_genotype_is_lifted_before_qc(with_plink2: Project, cohort_dir: Path) -> None:
    from efgpp.data.artifacts import ArtifactStore
    from efgpp.data.registry import Registry
    from efgpp.workflow.executor import prepare
    from tests.helpers import add_genotype, add_phenotypes

    p = add_phenotypes(add_genotype(with_plink2, cohort_dir, build="GRCh37"), cohort_dir, ["trait_a"])
    _toy_chain(p.resource_root / "liftover" / "hg19ToHg38.over.chain.gz", offset=1000)
    status = {o.step_id: o.status for o in prepare(p, engine="builtin")}
    assert status["harmonize.GENO001"] == "done" and status["genotype_qc.GENO001"] == "done", status
    with Registry.open(p) as reg:
        store = ArtifactStore(reg)
        lifted = store.latest(source_id="GENO001", artifact_type="liftover_genotype")
        qc = store.latest(source_id="GENO001", artifact_type="qc_genotype")
        vt = store.latest(source_id="GENO001", artifact_type="variant_table")
        src = store.latest(source_id="GENO001", artifact_type="source")
    assert lifted is not None and lifted.genome_build == "GRCh38" and lifted.tool == "pyliftover+plink2"
    assert qc is not None and qc.genome_build == "GRCh38" and lifted.artifact_id in qc.parent_artifact_ids
    assert vt is not None and vt.genome_build == "GRCh38"
    original = pl.read_csv(cohort_dir / "genotype" / "cohort.bim", separator="\t", has_header=False)
    new = pl.read_parquet(vt.path)
    assert new.get_column("position").min() == original.get_column("column_4").min() + 1000
    assert src is not None and src.genome_build == "GRCh37"  # the original is untouched


def test_gwas_step_writes_grch38_reference_artifact(project: Project, tmp_path: Path, monkeypatch, run_cli) -> None:  # type: ignore[no-untyped-def]
    from efgpp.data import gwas as gwas_mod
    from efgpp.data.artifacts import ArtifactStore
    from efgpp.data.provenance import RunRecord
    from efgpp.data.registry import Registry

    sumstats = tmp_path / "migraine.gz"
    with gzip.open(sumstats, "wt") as fh:
        fh.write("CHR\tBP\tSNP\tA1\tA2\tN\tSE\tP\tOR\n1\t101592213\trs1\tT\tC\t480359\t0.0153\t0.34\t1.01\n")
    monkeypatch.chdir(project.root)
    run_cli("data", "add", "gwas", "--path", str(sumstats), "--trait", "migraine", "--col", "chrom=CHR",
            "--col", "pos=BP", "--col", "snpid=SNP", "--col", "ea=A1", "--col", "nea=A2", "--col", "p=P")
    p = Project.load(project.root)
    assert p.data.gwas[0].id == "GWAS001" and p.data.gwas[0].columns["pos"] == "BP"
    _toy_chain(p.resource_root / "liftover" / "hg19ToHg38.over.chain.gz", offset=1000)

    def fake_run_tool(project, tool, args, **_kw):  # type: ignore[no-untyped-def]
        job = json.loads(Path(args[1]).read_text())
        assert tool == "gwaslab" and job["target"] == "38" and job["chains"]["19->38"].endswith("hg19ToHg38.over.chain.gz")
        pl.DataFrame({"CHR": [1], "POS": [101593213], "EA": ["T"], "NEA": ["C"], "P": [0.34]}).write_parquet(job["output"])
        Path(job["report"]).write_text(json.dumps({"gwaslab_version": "4.2.3", "build_detected": "19",
                                                   "lifted": "19->38", "rows_final": 1}))
        return RunRecord("RUN1", None, tool, "4.2.3", [tool, *args], exit_code=0, status="completed")

    monkeypatch.setattr(gwas_mod, "run_tool", fake_run_tool)
    out = gwas_mod.run_gwas(p, "GWAS001")
    assert out["lifted"] == "19->38"
    with Registry.open(p) as reg:
        art = ArtifactStore(reg).latest(source_id="GWAS001", artifact_type="gwas_sumstats")
    assert art is not None and art.origin.value == "reference" and art.genome_build == "GRCh38"
    assert Path(art.path).is_relative_to(p.resource_root / "gwas")
    assert "GWAS001" in run_cli("data", "validate", expect=None).output
