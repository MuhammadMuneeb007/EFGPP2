"""gnomAD population allele frequencies from a (user-provided or installed) sites VCF.

gnomAD releases are very large; point `resources.yaml: gnomad.path` at an existing
sites VCF (optionally restricted to your variants) rather than downloading everything.
"""

from __future__ import annotations

from pathlib import Path

from efgpp.data.annotation import (
    cohort_variants,
    info_field,
    lookup_table,
    register_annotation,
    resource_version,
)
from efgpp.project import Project

DEFAULT_FIELDS = ("AF", "AF_afr", "AF_amr", "AF_eas", "AF_nfe", "AF_sas", "nhomalt")


def run_gnomad(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    variants, vt, build = cohort_variants(project, source_id)
    cfg = project.resources.gnomad
    version, path = resource_version(project, "gnomad")
    if cfg.path:
        path, version = project.resolve(cfg.path), cfg.version or version
    if path is None or not Path(path).exists():
        raise RuntimeError("gnomAD sites VCF not configured (resources.yaml: gnomad.path)")
    fields = (cfg.model_extra or {}).get("fields", list(DEFAULT_FIELDS))
    table = lookup_table(Path(path), variants, select_sql=", ".join(info_field(f) for f in fields))
    return register_annotation(
        project, source_id=source_id, name="gnomad", table=table, parent=vt, tool="efgpp", tool_version=None,
        resource_versions={"gnomad": version or "unknown"}, genome_build=build,
        metadata={"matched_variants": table.height, "fields": fields},
    )


def run_dbsnp(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    """Attach dbSNP rsIDs by exact chrom/pos/ref/alt match."""
    variants, vt, build = cohort_variants(project, source_id)
    cfg = project.resources.dbsnp
    version, path = resource_version(project, "dbsnp")
    if cfg.path:
        path, version = project.resolve(cfg.path), cfg.version or version
    if path is None or not Path(path).exists():
        raise RuntimeError("dbSNP VCF not configured (resources.yaml: dbsnp.path)")
    table = lookup_table(Path(path), variants, select_sql="r.ID AS rsid")
    return register_annotation(
        project, source_id=source_id, name="dbsnp", table=table, parent=vt, tool="efgpp", tool_version=None,
        resource_versions={"dbsnp": version or "unknown"}, genome_build=build,
    )
