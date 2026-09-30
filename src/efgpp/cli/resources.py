"""`efgpp resources ...` commands."""

from __future__ import annotations

import typer

from efgpp.cli.common import console, load_project, table

app = typer.Typer(help="Reference knowledge: list, install, check for updates.", no_args_is_help=True)


@app.command("list")
def list_() -> None:
    """Configured and installed reference resources."""
    from efgpp.resources.manager import list_resources

    project = load_project()
    rows = [[r["name"], "yes" if r["enabled"] else "", r["installed_version"], r["genome_build"], r["provider"],
             r["local_path"]] for r in list_resources(project)]
    console.print(table("Resources", ["name", "enabled", "installed", "build", "provider", "path"], rows))


@app.command("install")
def install(name: str = typer.Argument(...), build: str | None = typer.Option(None, "--build"),
            force: bool = typer.Option(False, "--force", help="replace an installed copy of the same version")) -> None:
    """Download, verify and register a resource (versioned; frozen resources are never overwritten)."""
    from efgpp.resources.manager import install as do_install

    project = load_project()
    try:
        res = do_install(project, name, build=build, force=force)
    except Exception as exc:  # noqa: BLE001
        console.print(f"[red]✗ {name}: {exc}[/]")
        raise typer.Exit(1) from exc
    console.print(f"[green]✓[/] {name} {res.version} {res.status} -> {res.path}")


@app.command("update")
def update(check: bool = typer.Option(True, "--check/--no-check")) -> None:
    """Report newer upstream versions (never updates automatically)."""
    from efgpp.resources.manager import check_updates

    project = load_project()
    rows = [[r["name"], r["installed"], r["upstream"], r["action"]] for r in check_updates(project)]
    console.print(table("Resource updates", ["name", "installed", "upstream", "action"], rows))
