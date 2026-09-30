"""Genotype-derived molecular data: carriers, consequences, predictions - phenotype-independent and frozen."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from efgpp.constants import Modality
from efgpp.data.genotype.carriers import run_participant_variants
from efgpp.data.genotype.consequences import run_consequences
from efgpp.data.predicted import engine
from efgpp.data.registry import Registry
from efgpp.data.snapshots import DataSnapshot, freeze
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths
from efgpp.workflow.executor import prepare
from tests.helpers import add_genotype, add_phenotypes
from tests.unit._genotypes import write_predictdb

pytestmark = pytest.mark.integration


def _outputs(p: Project) -> dict[str, str]:
    """Checksums of every genotype-derived output (logs/manifests with timestamps excluded)."""
    roots = {
        "participant_variants": p.data_root / "derived" / "participant_variants" / "GENO001",
        "consequences": p.data_root / "derived" / "consequence_burden" / "GENO001",
        "predicted": p.data_root / "predicted" / "expression" / "gtex_v8" / "Whole_Blood",
    }
    out = {}
    for name, root in roots.items():
        files = sorted(f for f in root.rglob("*.parquet"))
        assert files, name
        out[name] = checksum_paths(files).value
    return out


def _derive(p: Project) -> dict[str, str]:
    run_participant_variants(p, "GENO001")
    run_consequences(p, "GENO001")
    engine.run_unit(p, Modality.EXPRESSION, "Whole_Blood")
    return _outputs(p)


def _setup(project: Project, cohort_dir: Path, phenotype_file: Path) -> Project:
    p = add_phenotypes(add_genotype(project, cohort_dir), cohort_dir, ["trait_a"])
    p.data.observed.phenotypes[0].path = str(phenotype_file)
    bim = pl.read_csv(cohort_dir / "genotype" / "cohort.bim", separator="\t", has_header=False,
                      new_columns=["chrom", "id", "cm", "pos", "a1", "a2"], infer_schema=False).head(4)
    write_predictdb(p.resource_root / "models" / "mashr_Whole_Blood.db",
                    [("ENSG01", r["id"], f"chr{r['chrom']}_{r['pos']}_{r['a2']}_{r['a1']}_b38", r["a2"], r["a1"], w)
                     for r, w in zip(bim.to_dicts(), (0.5, -0.2, 0.1, 0.3), strict=True)])
    cfg = p.data.predicted.expression
    cfg.enabled, cfg.tissues, cfg.engine = True, ["Whole_Blood"], "plink_score"
    cfg.model_provider.models_dir = str(p.resource_root / "models")
    p.data.participant_variants.enabled = True
    p.save_data_config()
    p = Project.load(p.root)
    prepare(p)
    return p


def test_phenotype_independence_and_snapshot(project: Project, cohort_dir: Path, tmp_path: Path) -> None:
    pheno = pl.read_csv(cohort_dir / "phenotypes.csv", infer_schema=False)
    original = tmp_path / "pheno_original.csv"
    pheno.write_csv(original)
    p = _setup(project, cohort_dir, original)
    before = _derive(p)

    # Replace every phenotype value (cases <-> controls) and derive again from scratch.
    flipped = pheno.with_columns(pl.when(pl.col("trait_a") == "1").then(pl.lit("0")).otherwise(pl.lit("1"))
                                 .alias("trait_a"))
    original.write_text(flipped.write_csv(), encoding="utf-8")
    prepare(p, force=True)
    after = _derive(Project.load(p.root))
    assert before == after  # carriers, consequence counts, gene burden and predictions are identical

    counts = pl.read_parquet(p.data_root / "derived" / "consequence_burden" / "GENO001" /
                             "participant_consequence_counts.parquet")
    assert counts.height == 120 and counts.get_column("n_variants_total").min() > 0  # type: ignore[operator]
    from efgpp.reporting.data_report import build_report

    html = build_report(p).read_text(encoding="utf-8")
    assert "Genotype-Derived Molecular Representations" in html and "Genetically predicted data" in html
    snap = freeze(p, "derived_v1")
    content = DataSnapshot.load(p, "derived_v1").content
    assert any(e["artifact_type"] == "gene_burden" for e in content["artifacts"])
    with Registry.open(p) as reg:
        kinds = {r["artifact_type"] for r in reg.rows("SELECT artifact_type FROM artifacts WHERE status <> 'SUPERSEDED'")}
    assert {"participant_variants", "participant_annotated_variants", "consequence_counts", "gene_burden",
            "predicted"} <= kinds
    assert snap.exists() and (p.path("snapshots", "derived_v1") / "molecular_models.parquet").exists()
    text = snap.read_text(encoding="utf-8")
    assert "consequence_burden" in text and "gene_burden" in text and "predicted" in text


def test_array_with_position_zero_probes_standardizes(project: Project, tmp_path: Path) -> None:
    """Arrays carry unplaced control probes at position 0; they must not stop standardization."""
    from efgpp.config.data import GenotypeSource
    from efgpp.constants import StorageMode
    from tests.unit._genotypes import write_bed

    samples = [f"S{i}" for i in range(20)]
    variants = [("0", "probe1", 0, "A", "G"), ("0", "probe2", 0, "C", "T")] + \
               [("1", f"rs{i}", 1000 + i * 10, "G", "A") for i in range(30)]
    counts = [[(i + j) % 3 for j in range(20)] for i in range(len(variants))]
    prefix = write_bed(tmp_path / "arr" / "array", samples, variants, counts)
    project.data.observed.genotype.append(GenotypeSource(id="GENO001", path=str(prefix), format="bed",
                                                         genome_build="GRCh38", mode=StorageMode.REFERENCE))
    project.save_data_config()
    status = {o.step_id: o.status for o in prepare(Project.load(project.root))}
    assert status["standardize.GENO001"] == "done"
    with Registry.open(project) as reg:
        n = reg.scalar("SELECT feature_count FROM artifacts WHERE artifact_type = 'variant_table' "
                       "AND status <> 'SUPERSEDED'")
    assert n == 30  # the two position-0 probes are left out of the variant table only
