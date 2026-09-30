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
    path: str = typer.Option(..., "--path"),
    name: str | None = typer.Option(None, "--name",
                                    help="phenotype name (default: the file name for one-column files, else the column)"),
    id_column: str | None = typer.Option(None, "--id-column", help="default: inferred (IID, participant_id, eid, ...)"),
    value_column: str | None = typer.Option(None, "--value-column", help="default: every phenotype column in the file"),
    type: PhenotypeType | None = typer.Option(None, "--type", help="default: inferred (binary/continuous/multiclass)"),
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
    """Add phenotype(s) from a table. The ID column, the phenotype column(s), their type and the
    case/control coding are inferred by modules/phenotype.py unless given."""
    from pathlib import Path

    from efgpp.data.io import read_table
    from efgpp.modules import load

    project = load_project()
    file_path = project.resolve(user_path(project, path))
    df = read_table(file_path)
    module = load(project, "phenotype")
    guess = module.infer_phenotypes(df, id_column, [value_column] if value_column else None)
    if guess.id_column is None:
        console.print(f"[red]no participant ID column found; pass --id-column (columns: {', '.join(df.columns)})[/]")
        raise typer.Exit(2)
    if not guess.phenotypes:
        console.print(f"[red]no phenotype column found[/] ({guess.skipped}); pass --value-column")
        raise typer.Exit(2)
    stem = Path(file_path.name).name.split(".")[0]
    single = len(guess.phenotypes) == 1
    rows, added = [], []
    for info in guess.phenotypes:
        ptype = type or PhenotypeType(info.type)
        pname = name if (name and single) else module.default_name(info.column, stem, single)
        pid = phenotype_id if (phenotype_id and single) else project.data.next_id("PH")
        spec = dict(
            id=pid, name=pname, path=user_path(project, path), participant_id_column=guess.id_column,
            value_column=info.column, type=ptype,
            levels=_split(levels) or (info.levels if ptype in (PhenotypeType.MULTICLASS, PhenotypeType.ORDINAL) else None),
            case_values=_split(case_values) or (info.case_values if ptype == PhenotypeType.BINARY else None),
            control_values=_split(control_values) or (info.control_values if ptype == PhenotypeType.BINARY else None),
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
            console.print(f"[red]{pname}: {exc}[/]")
            raise typer.Exit(2) from exc
        rows.append([pid, pname, info.column, ptype.value, f"{info.n_observed:,}", info.note])
        added.append(pid)
    console.print(table(f"Phenotypes (ID column: {guess.id_column})", ["id", "name", "column", "type", "values", "coding"], rows))
    for col, why in guess.skipped.items():
        console.print(f"[dim]skipped {col}: {why}[/]")
    from efgpp.cli.data import _register

    for pid in added:
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
