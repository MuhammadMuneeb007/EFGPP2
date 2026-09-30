"""OpenCRAVAT: modular annotation aggregation (used where modules reduce custom parsers)."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from efgpp.data.annotation import cohort_variants, register_annotation, write_sites_vcf
from efgpp.data.provenance import run_tool
from efgpp.project import Project

GENOME_CODES = {"GRCh37": "hg19", "GRCh38": "hg38"}


def parse_cravat_tsv(path: Path) -> pl.DataFrame:
    """Parse an OpenCRAVAT TSV report: '#' comment lines, then a header row."""
    with open(path, encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() and not ln.startswith("#")]
    if not lines:
        return pl.DataFrame()
    header = lines[0].split("\t")
    header = [h if header.count(h) == 1 else f"{h}_{i}" for i, h in enumerate(header)]
    rows = [ln.split("\t") + [""] * (len(header) - len(ln.split("\t"))) for ln in lines[1:]]
    return pl.DataFrame([r[: len(header)] for r in rows], schema=[(h, pl.Utf8) for h in header], orient="row")


def run_opencravat(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    variants, vt, build = cohort_variants(project, source_id)
    if build not in GENOME_CODES:
        raise RuntimeError(f"{source_id}: OpenCRAVAT needs a known genome build (got {build!r})")
    work = project.work_root / "opencravat" / source_id
    sites = write_sites_vcf(variants, work / "sites.vcf", build)
    args = ["run", str(sites), "-l", GENOME_CODES[build], "-t", "tsv", "-d", str(work), "--mp", str(threads)]
    annotators = project.resources.opencravat.annotators
    if annotators:
        args += ["-a", *annotators]
    rec = run_tool(project, "oc", args, step_id=step_id, inputs=[vt.artifact_id])  # type: ignore[list-item]
    reports = sorted(work.glob("*.variant.tsv")) or sorted(p for p in work.glob("*.tsv") if p.name != "sites.vcf")
    if not reports:
        raise RuntimeError(f"OpenCRAVAT produced no TSV report in {work}")
    table = parse_cravat_tsv(reports[0])
    return register_annotation(
        project, source_id=source_id, name="opencravat", table=table, parent=vt, tool="opencravat",
        tool_version=rec.tool_version, resource_versions={"opencravat_annotators": ",".join(annotators) or "default"},
        genome_build=build, metadata={"report": str(reports[0])},
    )
