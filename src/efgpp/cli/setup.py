"""`efgpp setup ...` and `efgpp doctor`."""

from __future__ import annotations

import typer

from efgpp.cli.common import console, load_project, state
from efgpp.project import Project, ProjectNotFoundError

app = typer.Typer(help="Install and verify scientific software in isolated environments.", no_args_is_help=True)

SYMBOL = {"ok": "[green]✓[/]", "installed": "[green]✓[/]", "already-installed": "[green]✓[/]", "replaced": "[green]✓[/]",
          "skipped": "[dim]-[/]", "warning": "[yellow]![/]", "failed": "[red]✗[/]"}


def _print(report) -> None:  # type: ignore[no-untyped-def]
    for s in report.steps:
        console.print(f"{SYMBOL.get(s.status, '?')} {s.name:<28} {s.detail}")
    if not report.ok:
        raise typer.Exit(1)


@app.command("data")
def data(dry_run: bool = typer.Option(False, "--dry-run"),
         components: str = typer.Option("core,genetics,reporting", "--components",
                                        help="comma-separated: core,genetics,annotation,metaxcan,reporting")) -> None:
    """Detect platform, bootstrap Pixi/Micromamba, install PLINK2, bcftools, PCA backends, Snakemake, reporting; lock."""
    from efgpp.setup.manager import setup_data

    project = load_project()
    report = setup_data(project, components={c.strip() for c in components.split(",") if c.strip()},
                        dry_run=dry_run, progress=lambda m: console.print(f"[dim]{m}...[/]"))
    _print(report)


@app.command("annotation")
def annotation(build: str | None = typer.Option(None, "--build")) -> None:
    """Install VEP + matching cache/FASTA and OpenCRAVAT; validate; record releases."""
    from efgpp.setup.manager import setup_annotation

    project = load_project()
    _print(setup_annotation(project, build=build, progress=lambda m: console.print(f"[dim]{m}...[/]")))


def doctor() -> None:
    """Show which components are installed and what is missing."""
    from efgpp.setup.doctor import run_doctor

    try:
        project: Project | None = Project.load(state.project_root)
    except ProjectNotFoundError:
        project = None
    console.print("\n[bold]EFGPP DATA SYSTEM[/]")
    group = None
    sym = {"ok": "[green]✓[/]", "missing": "[red]✗[/]", "warn": "[yellow]![/]", "info": "[blue]i[/]"}
    for line in run_doctor(project):
        if line.group != group:
            group = line.group
            console.print(f"\n[bold]{group}[/]")
        console.print(f"{sym[line.status]} {line.name:<22} [dim]{line.detail}[/]")
    console.print()
