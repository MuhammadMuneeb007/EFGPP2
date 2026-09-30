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
                                        help="comma-separated: core,genetics,annotation,metaxcan,reporting"),
         shared: bool = typer.Option(False, "--shared", help="download missing tools once for all your projects")) -> None:
    """Detect platform, build tool environments (Pixi/mamba), download anything still missing; lock."""
    from efgpp.setup.manager import setup_data

    project = load_project()
    report = setup_data(project, components={c.strip() for c in components.split(",") if c.strip()},
                        dry_run=dry_run, shared=shared, progress=lambda m: console.print(f"[dim]{m}...[/]"))
    _print(report)


@app.command("tools")
def tools(names: list[str] = typer.Argument(None, help="tools to install (default: all missing)"),
          shared: bool = typer.Option(False, "--shared", help="install into the per-user folder shared by all projects"),
          force: bool = typer.Option(False, "--force", help="reinstall even if the tool is already available")) -> None:
    """Download tools from their official sources: plink2 plink flashpca2 bcftools tabix bgzip
    snakemake multiqc oc predixcan vep."""
    from efgpp.setup.installers import INSTALLERS, install_root, install_tools

    try:
        project: Project | None = Project.load(state.project_root)
    except ProjectNotFoundError:
        project = None
        shared = True  # outside a project, tools go to the shared folder
    wanted = list(names) if names else list(INSTALLERS)
    root = install_root(project, shared)
    console.print(f"[dim]installing into {root}[/]")
    failed = 0
    for tool, status, detail in install_tools(project, wanted, shared=shared, force=force,
                                              progress=lambda m: console.print(f"[dim]{m}...[/]")):
        console.print(f"{SYMBOL.get(status, '?')} {tool:<12} {detail}")
        failed += status == "failed"
    console.print('\nto call them from your shell: [bold]eval "$(efgpp setup path)"[/]')
    if failed and names:
        raise typer.Exit(1)


@app.command("path")
def path() -> None:
    """Print an `export PATH=...` line for EFGPP-installed tools (use: eval "$(efgpp setup path)")."""
    import os

    from efgpp.setup.installers import path_exports

    try:
        project: Project | None = Project.load(state.project_root)
    except ProjectNotFoundError:
        project = None
    dirs = [str(d) for d in path_exports(project)]
    if os.name == "nt":
        print(f'$env:PATH = "{os.pathsep.join(dirs)}{os.pathsep}$env:PATH"')
    else:
        print(f'export PATH="{os.pathsep.join(dirs)}{os.pathsep}$PATH"')


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
