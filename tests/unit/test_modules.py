"""Data modules: column inference for phenotypes, covariates, genotypes and GWAS files."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest
import yaml

from efgpp.modules import export, load
from efgpp.modules.covariates import infer_covariates
from efgpp.modules.gwas import infer_gwas_columns, missing_essentials
from efgpp.modules.phenotype import default_name, infer_phenotypes, infer_type
from efgpp.project import Project


def _s(name: str, values: list[str]) -> pl.Series:
    return pl.Series(name, values, dtype=pl.Utf8)


def test_phenotype_types_and_coding() -> None:
    plink = infer_type(_s("Height", ["1", "2", "2", "-9", "1"]))
    assert plink.type == "binary" and plink.case_values == ["2"] and plink.control_values == ["1"]
    words = infer_type(_s("status", ["case", "control", "case"]))
    assert words.type == "binary" and words.case_values == ["case"]
    cont = infer_type(_s("bmi", [f"{20 + i * 0.37:.2f}" for i in range(40)]))
    assert cont.type == "continuous"
    multi = infer_type(_s("subtype", ["a", "b", "c", "a"]))
    assert multi.type == "multiclass" and multi.levels == ["a", "b", "c"]


def test_phenotype_file_inference_keeps_names() -> None:
    # Legacy EFGPP layout: IID FID Height, 1/2 coding.
    df = pl.DataFrame({"IID": ["a", "b", "c"], "FID": ["a", "b", "c"], "Height": ["1", "2", "1"]})
    guess = infer_phenotypes(df)
    assert guess.id_column == "IID" and [p.column for p in guess.phenotypes] == ["Height"]
    assert default_name("Height", "migraine", single=True) == "migraine"
    # A file mixing phenotypes and covariates: covariate-looking columns are skipped.
    df2 = pl.DataFrame({"eid": ["1", "2", "3"], "sex": ["M", "F", "M"], "age": ["50", "61", "47"],
                        "migraine": ["0", "1", "0"], "depression": ["1", "0", "0"]})
    g2 = infer_phenotypes(df2)
    assert g2.id_column == "eid" and [p.column for p in g2.phenotypes] == ["migraine", "depression"]
    assert set(g2.skipped) == {"sex", "age"}


def test_covariate_roles_and_types() -> None:
    df = pl.DataFrame({
        "FID": ["a", "b", "c"], "IID": ["a", "b", "c"], "Sex": ["1", "2", "1"], "Age": ["50", "61", "47"],
        "PC1": ["0.1", "-0.2", "0.05"], "genotyping_batch": ["3", "7", "3"], "smoking": ["never", "past", "never"],
        "bmi": ["21.5", "30.1", "25.0"],
    })
    g = infer_covariates(df)
    assert g.id_column == "IID" and g.variables == ["Sex", "Age", "PC1", "genotyping_batch", "smoking", "bmi"]
    assert g.roles == {"Sex": "sex", "Age": "age", "PC1": "pc", "genotyping_batch": "batch", "smoking": "other",
                       "bmi": "other"}
    assert set(g.categorical) == {"Sex", "genotyping_batch", "smoking"}
    # "sex" name but continuous values is not treated as sex.
    assert infer_covariates(pl.DataFrame({"IID": ["a", "b"], "sex_hormone": ["1.25", "3.5"]})).roles["sex_hormone"] == "other"


def test_gwas_column_recognition() -> None:
    # Legacy EFGPP GWAS header
    header = ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "OR", "INFO", "MAF"]
    m = infer_gwas_columns(header)
    assert m == {"snpid": "SNP", "chrom": "CHR", "pos": "BP", "ea": "A1", "nea": "A2", "OR": "OR", "se": "SE",
                 "p": "P", "n": "N", "info": "INFO"}
    assert missing_essentials(m) == []
    ssf = infer_gwas_columns(["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
                              "standard_error", "p_value", "rsid", "effect_allele_frequency"])
    assert ssf["pos"] == "base_pair_location" and ssf["eaf"] == "effect_allele_frequency" and ssf["rsid"] == "rsid"
    assert missing_essentials({"chrom": "c", "pos": "p"}) == ["effect and other allele", "beta, OR or z", "p-value"]


def test_project_copy_of_a_module_is_used(project: Project) -> None:
    assert load(project, "covariates").__name__ == "efgpp.modules.covariates"
    path = export(project, ["covariates"])[0]
    path.write_text(path.read_text(encoding="utf-8").replace('"age": [r"age"', '"age": [r"years", r"age"'),
                    encoding="utf-8")
    custom = load(project, "covariates")
    assert custom.covariate_role("years") == "age"  # the project's edit is active
    with pytest.raises(FileExistsError):
        export(project, ["covariates"])


def test_cli_infers_everything_for_legacy_files(project: Project, tmp_path: Path, run_cli, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    pheno = tmp_path / "migraine.height"
    pheno.write_text("IID FID Height\nS1 S1 1\nS2 S2 2\nS3 S3 2\n", encoding="utf-8")
    cov = tmp_path / "migraine.cov"
    cov.write_text("FID IID Sex Age PC1 PC2\nS1 S1 1 50 0.1 0.2\nS2 S2 2 61 -0.1 0.0\nS3 S3 1 47 0.3 -0.2\n",
                   encoding="utf-8")
    monkeypatch.chdir(project.root)
    out = run_cli("phenotype", "add", "--path", str(pheno)).output
    assert "PLINK coding" in out
    run_cli("data", "add", "covariates", "--path", str(cov))
    data = yaml.safe_load((project.root / "data.yaml").read_text(encoding="utf-8"))
    ph = data["observed"]["phenotypes"][0]
    assert (ph["name"], ph["value_column"], ph["type"], ph["participant_id_column"]) == ("migraine", "Height", "binary", "IID")
    assert ph["case_values"] == ["2"]
    cv = data["observed"]["covariates"][0]
    assert cv["variables"] == ["Sex", "Age", "PC1", "PC2"] and cv["categorical"] == ["Sex"]
    assert cv["roles"]["Sex"] == "sex" and cv["participant_id_column"] == "IID"
    from efgpp.workflow.executor import prepare

    prepare(Project.load(project.root))
    original = pl.read_parquet(project.root / "phenotypes" / "PH001" / "phenotype.parquet")
    assert original.columns == ["participant_id", "Height"]  # original column name kept
    covs = pl.read_parquet(next((project.root / "data" / "observed" / "covariates" / "COV001").glob("covariates.parquet")))
    assert {"Sex", "Age", "PC1", "PC2"} <= set(covs.columns)


def test_gwas_add_with_study_details(project: Project, tmp_path: Path, run_cli, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    import gzip

    f = tmp_path / "depression.gz"
    with gzip.open(f, "wt") as fh:
        fh.write("CHR\tBP\tSNP\tA1\tA2\tSE\tP\tOR\tINFO\n1\t100\trs1\tA\tG\t0.1\t0.5\t1.01\t0.9\n")
    monkeypatch.chdir(project.root)
    out = run_cli("data", "add", "gwas", "--path", str(f), "--trait", "depression", "--phenotype", "migraine",
                  "--ancestry", "European", "--n-cases", "1000", "--n-controls", "4000").output
    assert "European" in out and "1,000 cases" in out
    g = Project.load(project.root).data.gwas[0]
    assert (g.trait, g.phenotypes, g.ancestry, g.n_cases, g.n_controls) == ("depression", ["migraine"], "European", 1000, 4000)
    assert g.columns["pos"] == "BP" and g.columns["OR"] == "OR"
    from efgpp.data.gwas import gwas_for_phenotype

    p = Project.load(project.root)
    assert [x.id for x in gwas_for_phenotype(p, "migraine")] == ["GWAS001"]
    assert gwas_for_phenotype(p, "some_other_trait") == []
    assert "not run" in run_cli("data", "gwas", "list", "--phenotype", "migraine").output


def test_resources_enable(project: Project, run_cli, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    monkeypatch.chdir(project.root)
    assert "enabled: vep, clinvar" in run_cli("resources", "enable", "vep", "clinvar").output
    r = Project.load(project.root).resources
    assert r.vep.enabled and r.clinvar.enabled and not r.gnomad.enabled
    run_cli("resources", "enable", "clinvar", "--off")
    assert not Project.load(project.root).resources.clinvar.enabled
    assert run_cli("resources", "enable", "nosuch", expect=1).exit_code == 1


def test_data_predict_switch(project: Project, run_cli, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    monkeypatch.chdir(project.root)
    run_cli("data", "predict", "enable", "expression", "--tissue", "Whole_Blood", "--tissue", "Brain_Cortex")
    cfg = Project.load(project.root).data.predicted.expression
    assert cfg.enabled and cfg.tissues == ["Whole_Blood", "Brain_Cortex"] and cfg.engine == "metaxcan"
    run_cli("data", "predict", "enable", "proteomics", "--dataset", "OPD000001")
    prot = Project.load(project.root).data.predicted.proteomics
    assert prot.enabled and prot.datasets == ["OPD000001"] and prot.provider == "omicspred"
    # a methylation engine that makes no sense is refused
    assert run_cli("data", "predict", "enable", "methylation", "--engine", "metaxcan", expect=2).exit_code == 2
    assert "plan" in run_cli("data", "predict", "plan").output
    from efgpp.data.predicted.engine import prediction_plan

    rows = {r.item: r for r in prediction_plan(Project.load(project.root))}
    assert rows["Predicted expression [Whole_Blood]"].status == "NOT INSTALLED"
    assert "predictdb-gtex-v8-expression" in rows["Predicted expression [Whole_Blood]"].required
    assert rows["Predicted proteomics [OPD000001]"].required == "efgpp resources install omicspred --dataset OPD000001"
    run_cli("data", "predict", "enable", "expression", "--off")
    assert not Project.load(project.root).data.predicted.expression.enabled
    run_cli("data", "variants", "enable", "--spliceai")
    pv = Project.load(project.root).data.participant_variants
    assert pv.enabled and pv.carrier_only and pv.annotations.spliceai
