"""`efgpp resources ...` commands."""

from __future__ import annotations

from typing import Any

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


@app.command("enable")
def enable(names: list[str] = typer.Argument(..., help="resource names, e.g. vep clinvar alphamissense"),
           off: bool = typer.Option(False, "--off", help="disable instead")) -> None:
    """Switch resources on (or off with --off) in resources.yaml; enabled annotations run in `data prepare`."""
    project = load_project()
    known = list(type(project.resources).model_fields)
    unknown = [n for n in names if n not in known]
    if unknown:
        console.print(f"[red]✗ unknown resource(s): {', '.join(unknown)}[/] (known: {', '.join(known)})")
        raise typer.Exit(1)
    for n in names:
        getattr(project.resources, n).enabled = not off
    project.save_resources_config()
    console.print(f"[green]✓[/] {'disabled' if off else 'enabled'}: {', '.join(names)}")


@app.command("install")
def install(
    name: str = typer.Argument(..., help="e.g. genome, vep, clinvar, predictdb-gtex-v8-expression, omicspred, mimosa"),
    build: str | None = typer.Option(None, "--build"),
    force: bool = typer.Option(False, "--force", help="replace an installed copy of the same version"),
    dataset: str | None = typer.Option(None, "--dataset", help="omicspred: OPD id; predictdb-protein: .db file"),
    url: str | None = typer.Option(None, "--url", help="hibag: classifier URL"),
    path: str | None = typer.Option(None, "--path", help="hibag: local classifier file"),
    ancestry: str | None = typer.Option(None, "--ancestry", help="hibag: training population"),
    platform: str | None = typer.Option(None, "--platform", help="hibag: genotyping platform"),
    model_build: str | None = typer.Option(None, "--model-build", help="hibag / predictdb-protein: model genome build"),
    loci: str | None = typer.Option(None, "--loci", help="hibag: HLA loci in the classifier, e.g. A,B,C,DRB1"),
) -> None:
    """Download, verify and register a resource (versioned; frozen resources are never overwritten)."""
    from efgpp.resources.manager import install as do_install

    project = load_project()
    options: dict[str, Any] = {k: v for k, v in {"dataset": dataset, "url": url, "path": path, "ancestry": ancestry,
                                 "platform": platform, "model_build": model_build, "loci": loci}.items() if v}
    options["progress"] = lambda m: console.print(f"[dim]{m}[/]")
    try:
        res = do_install(project, name, build=build, force=force, **options)
    except Exception as exc:  # noqa: BLE001
        console.print(f"[red]✗ {name}: {exc}[/]")
        raise typer.Exit(1) from exc
    console.print(f"[green]✓[/] {name} {res.version} {res.status} -> {res.path}")


omicspred_app = typer.Typer(help="OmicsPred genetic scores (proteomics, metabolomics, transcriptomics).",
                            no_args_is_help=True)
app.add_typer(omicspred_app, name="omicspred")


@omicspred_app.command("refresh")
def omicspred_refresh() -> None:
    """Fetch the current OmicsPred dataset catalogue from its REST API (cached with retrieval date)."""
    import httpx

    from efgpp.data.references.molecular import omicspred_catalog_path
    from efgpp.data.references.molecular import omicspred_refresh as refresh

    project = load_project()
    with httpx.Client(follow_redirects=True, timeout=120) as client:
        try:
            cat = refresh(project, client)
        except Exception as exc:  # noqa: BLE001
            console.print(f"[red]✗ OmicsPred: {exc}[/] (cached catalogue kept)")
            raise typer.Exit(1) from exc
    console.print(f"[green]✓[/] {cat['count']} datasets, retrieved {cat['retrieved_at']} -> "
                  f"{omicspred_catalog_path(project)}")


@omicspred_app.command("list")
def omicspred_list(
    modality: str | None = typer.Option(None, "--modality", help="proteomics | metabolomics | transcriptomics"),
    cohort: str | None = typer.Option(None, "--cohort"),
    platform: str | None = typer.Option(None, "--platform"),
    ancestry: str | None = typer.Option(None, "--ancestry"),
) -> None:
    """Datasets in the cached OmicsPred catalogue (install one with --dataset <id>)."""
    import httpx

    from efgpp.data.references.molecular import omicspred_catalog, omicspred_rows

    project = load_project()
    with httpx.Client(follow_redirects=True, timeout=120) as client:
        cat = omicspred_catalog(project, client)
    rows = omicspred_rows(cat, modality=modality, cohort=cohort, platform=platform, ancestry=ancestry)
    console.print(table(f"OmicsPred datasets (catalogue {cat['retrieved_at']})",
                        ["id", "name", "type", "platform", "tissue", "scores", "training cohort", "ancestry", "n"],
                        [[r["id"], r["name"], r["omics_type"], r["platform"], r["tissue"], r["scores"],
                          r["training_cohort"], r["training_ancestry"], r["training_n"]] for r in rows]))
    console.print("install: efgpp resources install omicspred --dataset <id>")


@app.command("update")
def update(check: bool = typer.Option(True, "--check/--no-check")) -> None:
    """Report newer upstream versions (never updates automatically)."""
    from efgpp.resources.manager import check_updates

    project = load_project()
    rows = [[r["name"], r["installed"], r["upstream"], r["action"]] for r in check_updates(project)]
    console.print(table("Resource updates", ["name", "installed", "upstream", "action"], rows))
