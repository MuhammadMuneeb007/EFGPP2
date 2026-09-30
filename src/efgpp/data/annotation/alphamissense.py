"""AlphaMissense as a reference score source (precomputed scores; nothing is trained)."""

from __future__ import annotations

from efgpp.data.annotation import (
    cohort_variants,
    lookup_table,
    register_annotation,
    resource_version,
)
from efgpp.project import Project


def run_alphamissense(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    variants, vt, build = cohort_variants(project, source_id)
    version, path = resource_version(project, "alphamissense")
    cfg = project.resources.alphamissense
    if cfg.path:
        path, version = project.resolve(cfg.path), cfg.version or version
    if path is None or not path.exists():
        raise RuntimeError("AlphaMissense scores not installed; run `efgpp resources install alphamissense`")
    table = lookup_table(
        path, variants,
        select_sql=("r.am_pathogenicity AS alphamissense_score, r.am_class AS alphamissense_class, "
                    "r.uniprot_id, r.transcript_id, r.protein_variant"),
    )
    table = table.with_columns(table.get_column("alphamissense_score").cast(float, strict=False))
    return register_annotation(
        project, source_id=source_id, name="alphamissense", table=table, parent=vt, tool="efgpp",
        tool_version=None, resource_versions={"alphamissense": version or "unknown"}, genome_build=build,
        metadata={"matched_variants": table.height},
    )
