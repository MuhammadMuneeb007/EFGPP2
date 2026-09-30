"""MIMOSA model selection and conversion (R export -> standardized weights)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import polars as pl
import pytest

from efgpp.data.predicted.mimosa import EXPORT_SCRIPT, read_export, select_models, standardized

HEADER = "cpg_id\tmethod\tmethod_index\tvalid\ttest_r2\tlambda\tsnp\tchromosome\tposition\ta1\ta2\tweight\n"


def _export(tmp_path: Path) -> Path:
    rows = []
    # cg01: ElNet 0.12, MNet 0.09, SCAD invalid, MCP 0.14, LASSO 0.07 -> MCP selected
    for i, (method, valid, r2) in enumerate([("ElNet", "TRUE", 0.12), ("MNet", "TRUE", 0.09), ("SCAD", "FALSE", "NA"),
                                              ("MCP", "TRUE", 0.14), ("LASSO", "TRUE", 0.07)], start=1):
        w1, w2 = (0.5, -0.25) if method == "MCP" else (0.1, 0.1)
        rows.append(f"cg01\t{method}\t{i}\t{valid}\t{r2}\t0.01\trs1\t1\t100\tG\tA\t{w1}")
        rows.append(f"cg01\t{method}\t{i}\t{valid}\t{r2}\t0.01\trs2\t1\t200\tT\tC\t{w2}")
    # cg02: tie between ElNet and LASSO -> ElNet (deterministic order); cg03: weak model
    rows += ["cg02\tElNet\t1\tTRUE\t0.2\t0.1\trs1\t1\t100\tG\tA\t1.0",
             "cg02\tLASSO\t5\tTRUE\t0.2\t0.1\trs1\t1\t100\tG\tA\t2.0",
             "cg03\tElNet\t1\tTRUE\t0.004\t0.1\trs1\t1\t100\tG\tA\t1.0"]
    p = tmp_path / "export.tsv"
    p.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")
    return p


def test_highest_test_r2_valid_model_selected(tmp_path: Path) -> None:
    weights, models = select_models(read_export(_export(tmp_path)), 0.005)
    m = {r["cpg_id"]: r for r in models.to_dicts()}
    assert (m["cg01"]["selected_method"], m["cg01"]["validation_r2"]) == ("MCP", 0.14)
    assert m["cg02"]["selected_method"] == "ElNet"
    assert m["cg03"]["below_threshold"] and not m["cg01"]["below_threshold"]
    cg01 = weights.filter(pl.col("feature_id") == "cg01").sort("position")
    assert cg01.select("variant_id", "effect_allele", "other_allele", "weight", "model_method").rows() == [
        ("rs1", "G", "A", 0.5, "MCP"), ("rs2", "T", "C", -0.25, "MCP")]
    w, feats = standardized(read_export(_export(tmp_path)), 0.005, genome_build="GRCh38", version="v2")
    assert set(w.get_column("genome_build").unique()) == {"GRCh38"}
    assert w.get_column("model_id").unique().to_list() == ["mimosa:v2"]
    assert feats.filter(pl.col("feature_id") == "cg03").get_column("below_threshold").item()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="R not installed")
def test_r_export_of_tiny_rds(tmp_path: Path) -> None:
    models = tmp_path / "models"
    models.mkdir()
    make = f"""
    snps <- data.frame(SNP = c("rs1", "rs2"), SNPChr = c(1, 1), SNPPos = c(100, 200), a1 = c("G", "T"),
                       a2 = c("A", "C"), CpG = "cg01")
    m <- function(ok, w, r2) list(ok, snps, w, r2, 0.01)
    saveRDS(list(m(TRUE, c(0.1, 0.1), 0.12), m(TRUE, c(0.1, 0.1), 0.09), m(FALSE, c(0, 0), numeric(0)),
                 m(TRUE, c(0.5, -0.25), 0.14), m(TRUE, c(0.1, 0.1), 0.07)), "{(models / 'cg01.rds').as_posix()}")
    """
    subprocess.run(["Rscript", "-e", make], check=True)
    out = tmp_path / "export.tsv"
    subprocess.run(["Rscript", str(EXPORT_SCRIPT), str(models), str(out)], check=True)
    weights, models_t = select_models(read_export(out), 0.005)
    assert models_t.row(0, named=True)["selected_method"] == "MCP"
    assert weights.sort("position").get_column("weight").to_list() == [0.5, -0.25]
