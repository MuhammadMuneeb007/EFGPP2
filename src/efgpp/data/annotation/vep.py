"""Ensembl VEP: canonical molecular-consequence annotation (offline cache mode)."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from efgpp.data.annotation import (
    cohort_variants,
    register_annotation,
    resource_version,
    write_sites_vcf,
)
from efgpp.data.provenance import run_tool
from efgpp.project import Project

VEP_FLAGS = ["--offline", "--cache", "--tab", "--symbol", "--canonical", "--biotype",
             "--variant_class", "--numbers", "--flag_pick", "--no_stats", "--force_overwrite"]


def vep_command(project: Project, sites: Path, out: Path, build: str, threads: int) -> list[str]:
    cfg = project.resources.vep
    cache = cfg.cache if cfg.cache != "auto" else str(project.resource_root / "vep")
    args = ["--input_file", str(sites), "--format", "vcf", "--output_file", str(out), *VEP_FLAGS,
            "--dir_cache", cache, "--assembly", build, "--species", cfg.species]
    if cfg.release:
        args += ["--cache_version", str(cfg.release)]
    if threads > 1:
        args += ["--fork", str(threads)]
    fasta = project.resources.genome.fasta
    if fasta:
        args += ["--fasta", str(project.resolve(fasta)), "--hgvs"]
    return [*args, *cfg.extra_args]


def parse_vep_tab(path: Path) -> pl.DataFrame:
    """Parse VEP --tab output; an `Extra` column (key=value;...) is expanded into columns."""
    with open(path, encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh if not ln.startswith("##")]
    header = lines[0].lstrip("#").split("\t")
    rows = [ln.split("\t") for ln in lines[1:] if ln]
    df = pl.DataFrame(rows, schema=[(h, pl.Utf8) for h in header], orient="row")
    if "Extra" in df.columns:
        extra = [dict(kv.split("=", 1) for kv in e.split(";") if "=" in kv) if e not in (None, "-") else {}
                 for e in df.get_column("Extra").to_list()]
        keys = sorted({k for d in extra for k in d})
        df = df.drop("Extra").with_columns([pl.Series(k, [d.get(k) for d in extra], dtype=pl.Utf8) for k in keys])
    return df.with_columns([pl.when(pl.col(c) == "-").then(None).otherwise(pl.col(c)).alias(c) for c in df.columns])


def consequence_counts(df: pl.DataFrame) -> dict[str, int]:
    if "Consequence" not in df.columns:
        return {}
    picked = df.filter(pl.col("PICK") == "1") if "PICK" in df.columns else df
    cons = picked.get_column("Consequence").str.split(",").explode()
    return {k: int(v) for k, v in cons.value_counts(sort=True).iter_rows()}


def run_vep(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    variants, vt, build = cohort_variants(project, source_id)
    if not build:
        raise RuntimeError(f"{source_id}: genome build unknown; VEP needs --assembly")
    work = project.work_root / "vep" / source_id
    sites = write_sites_vcf(variants, work / "sites.vcf", build)
    out = work / "vep.tsv"
    rec = run_tool(project, "vep", vep_command(project, sites, out, build, threads), step_id=step_id,
                   inputs=[vt.artifact_id])  # type: ignore[list-item]
    table = parse_vep_tab(out)
    cache_version, _ = resource_version(project, "vep")
    versions = {"vep": rec.tool_version or "unknown", "vep_cache": cache_version or str(project.resources.vep.release or "unknown")}
    return register_annotation(
        project, source_id=source_id, name="vep", table=table, parent=vt, tool="vep",
        tool_version=rec.tool_version, resource_versions=versions, genome_build=build,
        metadata={"consequences": consequence_counts(table), "command": " ".join(rec.command)},
    )
