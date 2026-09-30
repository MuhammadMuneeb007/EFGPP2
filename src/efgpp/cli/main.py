"""efgpp command-line entry point."""

from __future__ import annotations

from pathlib import Path

import typer

from efgpp import __version__
from efgpp.cli import data, phenotype, resources, setup
from efgpp.cli.common import console, ensure_utf8, load_project, state

app = typer.Typer(
    name="efgpp",
    help="EFGPP - Exploratory Framework for Genotype-Phenotype Prediction (Data layer).",
    no_args_is_help=True,
    pretty_exceptions_show_locals=False,
)
app.add_typer(data.app, name="data")
app.add_typer(phenotype.app, name="phenotype")
app.add_typer(setup.app, name="setup")
app.add_typer(resources.app, name="resources")
export_app = typer.Typer(help="Export execution artefacts.", no_args_is_help=True)
app.add_typer(export_app, name="export")


def _version(value: bool) -> None:
    if value:
        console.print(f"efgpp {__version__}")
        raise typer.Exit()


@app.callback()
def main_callback(
    project: Path | None = typer.Option(None, "--project", "-C", help="project root (default: search upwards from cwd)"),
    version: bool = typer.Option(False, "--version", callback=_version, is_eager=True),
) -> None:
    state.project_root = project


@app.command()
def init(path: Path = typer.Argument(Path(), help="project directory"),
         name: str | None = typer.Option(None, "--name"),
         force: bool = typer.Option(False, "--force")) -> None:
    """Create an EFGPP project (directory layout + project.yaml, data.yaml, resources.yaml)."""
    from efgpp.project import init_project

    try:
        project = init_project(path, name=name, force=force)
    except FileExistsError as exc:
        console.print(f"[red]{exc}[/]")
        raise typer.Exit(1) from exc
    console.print(f"[green]✓[/] EFGPP project [bold]{project.config.project.name}[/] initialised at {project.root}")
    console.print("next: [bold]efgpp setup data[/], then [bold]efgpp data add genotype ...[/] / "
                  "[bold]efgpp phenotype add ...[/]")


app.command("doctor")(setup.doctor)


@export_app.command("slurm")
def export_slurm(out: Path | None = typer.Option(None, "--out")) -> None:
    """Write hpc/*.sbatch scripts and submit_all.sh for the current data plan."""
    from efgpp.workflow.slurm import export_slurm as do_export

    project = load_project()
    for p in do_export(project, out):
        console.print(f"[green]✓[/] {project.relative(p)}")


def main() -> None:
    ensure_utf8()
    app()


if __name__ == "__main__":
    main()
