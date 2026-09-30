"""ClinVar clinical significance lookup from the official ClinVar VCF."""

from __future__ import annotations

from efgpp.data.annotation import (
    cohort_variants,
    info_field,
    lookup_table,
    register_annotation,
    resource_version,
)
from efgpp.project import Project


def run_clinvar(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    variants, vt, build = cohort_variants(project, source_id)
    version, path = resource_version(project, "clinvar")
    if path is None or not path.exists():
        raise RuntimeError("ClinVar is not installed; run `efgpp resources install clinvar`")
    table = lookup_table(
        path, variants,
        select_sql=", ".join([
            "r.ID AS clinvar_variation_id",
            info_field("CLNSIG"), info_field("CLNREVSTAT"), info_field("CLNDN"), info_field("GENEINFO"),
        ]),
    )
    return register_annotation(
        project, source_id=source_id, name="clinvar", table=table, parent=vt, tool="efgpp",
        tool_version=None, resource_versions={"clinvar": version or "unknown"}, genome_build=build,
        metadata={"matched_variants": table.height},
    )
