"""`efgpp phenotype ...` commands. Phenotypes are configuration; add as many as needed."""

from __future__ import annotations

import typer

from efgpp.cli.common import console, load_project, table, user_path
from efgpp.config.data import OntologyTerm, PhenotypeSource, TimelineSpec
from efgpp.constants import Origin, PhenotypeType, StorageMode

app = typer.Typer(help="Define and inspect phenotypes (any number, any supported type).", no_args_is_help=True)


def _split(v: str | None) -> list[str] | None:
    return [x.strip() for x in v.split(",") if x.strip()] if v else None


@app.command("add")
def add(
    name: str = typer.Option(..., "--name", help="phenotype name (free text; never interpreted)"),
    path: str = typer.Option(..., "--path"),
    id_column: str = typer.Option("participant_id", "--id-column"),
    value_column: str = typer.Option(..., "--value-column"),
    type: PhenotypeType = typer.Option(..., "--type"),
    levels: str | None = typer.Option(None, "--levels", help="comma-separated classes (ordered for ordinal)"),
    case_values: str | None = typer.Option(None, "--case-values", help="binary: comma-separated case codes"),
    control_values: str | None = typer.Option(None, "--control-values", help="binary: comma-separated control codes"),
    missing_values: str | None = typer.Option(None, "--missing-values"),
    units: str | None = typer.Option(None, "--units"),
    ontology_id: str | None = typer.Option(None, "--ontology-id", help="e.g. an HPO/EFO/MONDO term id"),
    ontology_label: str | None = typer.Option(None, "--ontology-label"),
    event_column: str | None = typer.Option(None, "--event-column"),
    mode: StorageMode = typer.Option(StorageMode.AUTO, "--mode"),
    origin: Origin = typer.Option(Origin.OBSERVED, "--origin"),
    phenotype_id: str | None = typer.Option(None, "--id"),
) -> None:
    """Add a phenotype column from a table."""
    project = load_project()
    pid = phenotype_id or project.data.next_id("PH")
    spec = dict(
        id=pid, name=name, path=user_path(project, path), participant_id_column=id_column, value_column=value_column,
        type=type, levels=_split(levels), case_values=_split(case_values), control_values=_split(control_values),
        units=units, mode=mode, origin=origin,
        ontology_term=OntologyTerm(id=ontology_id, label=ontology_label) if ontology_id else None,
        timeline=TimelineSpec(event_column=event_column) if event_column else None,
    )
    if missing_values is not None:
        spec["missing_values"] = _split(missing_values) or []
    try:
        project.data.observed.phenotypes.append(PhenotypeSource(**spec))  # type: ignore[arg-type]
        project.save_data_config()
    except ValueError as exc:
        console.print(f"[red]{exc}[/]")
        raise typer.Exit(2) from exc
    from efgpp.cli.data import _register

    _register(project, pid)


@app.command("list")
def list_() -> None:
    """Configured phenotypes."""
    project = load_project()
    rows = [[p.id, p.name, p.type.value, p.origin.value, p.value_column, p.path] for p in project.data.observed.phenotypes]
    console.print(table("Phenotypes", ["id", "name", "type", "origin", "column", "path"], rows))


@app.command("remove")
def remove(phenotype: str = typer.Argument(..., help="id or name")) -> None:
    """Remove a phenotype from data.yaml (registered artifacts are kept for provenance)."""
    project = load_project()
    p = project.data.phenotype(phenotype)
    project.data.observed.phenotypes = [x for x in project.data.observed.phenotypes if x.id != p.id]
    project.save_data_config()
    console.print(f"[green]✓[/] removed {p.id} ({p.name}) from data.yaml")


@app.command("export")
def export(format: str = typer.Option("phenopackets", "--format")) -> None:
    """Export standardized phenotypes (GA4GH Phenopackets v2 JSON)."""
    if format != "phenopackets":
        console.print("[red]only --format phenopackets is supported[/]")
        raise typer.Exit(2)
    from efgpp.data.phenopackets import export as do_export

    project = load_project()
    out, problems = do_export(project)
    for p in problems:
        console.print(f"[yellow]![/] {p}")
    console.print(f"[green]✓[/] phenopackets written to {project.relative(out)}")
