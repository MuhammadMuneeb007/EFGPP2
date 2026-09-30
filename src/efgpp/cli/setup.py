"""`efgpp setup ...` and `efgpp doctor`."""

from __future__ import annotations

import typer

from efgpp.cli.common import console, load_project, state
from efgpp.project import Project, ProjectNotFoundError

app = typer.Typer(help="Install and verify scientific software in isolated environments.", no_args_is_help=True)

SYMBOL = {"ok": "[green]✓[/]", "installed": "[green]✓[/]", "already-installed": "[green]✓[/]", "replaced": "[green]✓[/]",
          "skipped": "[dim]-[/]", "warning": "[yellow]![/]", "failed": "[red]✗[/]", "missing": "[red]✗[/]",
          "manual": "[yellow]![/]"}


def _optional_project() -> Project | None:
    try:
        return Project.load(state.project_root)
    except ProjectNotFoundError:
        return None


def _print(report) -> None:  # type: ignore[no-untyped-def]
    for s in report.steps:
        console.print(f"{SYMBOL.get(s.status, '?')} {s.name:<28} {s.detail}")
    if not report.ok:
        raise typer.Exit(1)


@app.command("data")
def data(dry_run: bool = typer.Option(False, "--dry-run"),
         components: str = typer.Option("genetics,reporting,gwas", "--components",
                                        help="comma-separated: genetics,annotation,metaxcan,reporting,gwas"),
         shared: bool = typer.Option(False, "--shared", help="install into $EFGPP_TOOLS_HOME instead of ./software")) -> None:
    """Detect platform, build tool environments (Pixi/mamba), download anything still missing; lock."""
    from efgpp.setup.manager import setup_data

    project = load_project()
    report = setup_data(project, components={c.strip() for c in components.split(",") if c.strip()},
                        dry_run=dry_run, shared=shared, progress=lambda m: console.print(f"[dim]{m}...[/]"))
    _print(report)


@app.command("tools")
def tools(names: list[str] = typer.Argument(None, help="tools to install (default: all missing)"),
          shared: bool = typer.Option(False, "--shared", help="install into $EFGPP_TOOLS_HOME instead of ./software"),
          force: bool = typer.Option(False, "--force", help="reinstall even if the tool is already available")) -> None:
    """Download tools from their official sources: plink2 plink flashpca2 bcftools tabix bgzip
    multiqc oc predixcan vep gwaslab."""
    from efgpp.setup.installers import INSTALLERS, install_root, install_tools

    project = _optional_project()  # outside a project, tools go to ./software
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


@app.command("toolkit")
def toolkit(names: list[str] = typer.Argument(None, help="perl r prs prs-python prs-py27 simulation | all"),
            shared: bool = typer.Option(False, "--shared", help="install into $EFGPP_TOOLS_HOME instead of ./software"),
            list_: bool = typer.Option(False, "--list", help="show what each toolkit contains")) -> None:
    """Install complete toolkits: Perl, R + PRS R packages (LDpred-2, lassosum, ...), PRS binaries and
    repositories (PRSTools), Python 2.7/3.10 method environments, simulation (simuPOP, msprime)."""
    from efgpp.setup.installers import InstallError, install_root
    from efgpp.setup.toolkit_installer import install_toolkits
    from efgpp.setup.toolkits import ALIASES, TOOLKITS

    if list_ or not names:
        for tk in TOOLKITS.values():
            console.print(f"\n[bold]{tk.name}[/] - {tk.description}")
            if tk.conda:
                console.print(f"  conda ({tk.env}): {', '.join(tk.conda)}")
            if tk.pip:
                console.print(f"  pip: {', '.join(tk.pip)}")
            if tk.r_cran or tk.r_github or tk.r_bioc:
                console.print(f"  R: {', '.join([*tk.r_cran, *tk.r_bioc, *tk.r_github])}")
            if tk.downloads:
                console.print(f"  downloads: {', '.join(d.name for d in tk.downloads)}")
            if tk.repos:
                console.print(f"  repositories: {', '.join(r.github for r in tk.repos)}")
        console.print(f"\naliases: {', '.join(f'{k} = {v}' for k, v in ALIASES.items())}")
        return
    project = _optional_project()
    console.print(f"[dim]installing into {install_root(project, shared)}[/]")
    try:
        results = install_toolkits(project, list(names), shared=shared,
                                   progress=lambda m: console.print(f"[dim]{m}...[/]"))
    except (InstallError, KeyError) as exc:
        console.print(f"[red]{exc}[/]")
        raise typer.Exit(1) from exc
    failed = False
    for res in results:
        console.print(f"\n[bold]{res.toolkit}[/]")
        for item in res.items:
            console.print(f"  {SYMBOL.get(item.status, '?')} {item.kind:<9}{item.name:<24}{item.detail}")
        failed |= not res.ok
    console.print('\ncheck everything with [bold]efgpp setup check[/]; use tools from a shell with '
                  '[bold]eval "$(efgpp setup path)"[/]')
    if failed:
        raise typer.Exit(1)


@app.command("check")
def check(shared: bool = typer.Option(False, "--shared", help="check the $EFGPP_TOOLS_HOME install only"),
          toolkits: str = typer.Option("all", "--toolkits", help="comma-separated toolkits to check, or 'none'")) -> None:
    """Check every tool, toolkit, R package and repository: installed or missing."""
    from efgpp.setup.installers import COMPONENT_TOOLS
    from efgpp.setup.toolkit_installer import check_toolkit
    from efgpp.setup.toolkits import resolve_names
    from efgpp.setup.tools import ToolNotFoundError, resolve

    project = _optional_project()
    missing = 0
    console.print("\n[bold]DATA-LAYER TOOLS[/] (efgpp setup tools)")
    for component, tools in COMPONENT_TOOLS.items():
        for tool in tools:
            try:
                where = str(resolve(project, tool).path)
                console.print(f"  {SYMBOL['ok']} {tool:<12}{component:<12}[dim]{where}[/]")
            except ToolNotFoundError:
                missing += 1
                console.print(f"  {SYMBOL['missing']} {tool:<12}{component:<12}[dim]efgpp setup tools {tool}[/]")
    names = [] if toolkits == "none" else resolve_names([t for t in toolkits.split(",") if t])
    for name in names:
        res = check_toolkit(project, name, shared=shared)
        n_ok = sum(i.status in ("ok", "manual") for i in res.items)
        console.print(f"\n[bold]TOOLKIT {name}[/] ({n_ok}/{len(res.items)} present)")
        for item in res.items:
            if item.status in ("missing", "failed"):
                missing += 1
            console.print(f"  {SYMBOL.get(item.status, '?')} {item.kind:<9}{item.name:<24}[dim]{item.detail}[/]")
    console.print()
    if missing:
        console.print(f"[yellow]{missing} item(s) missing[/] - install with [bold]efgpp setup tools[/] and "
                      "[bold]efgpp setup toolkit all[/]")
        raise typer.Exit(1)
    console.print("[green]everything is installed[/]")


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
