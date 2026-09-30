"""`efgpp data predict ...`, `efgpp data variants ...`, `efgpp data hla ...`: genotype-derived molecular data."""

from __future__ import annotations

import typer

from efgpp.cli.common import console, load_project, table
from efgpp.constants import Modality

predict_app = typer.Typer(help="Genetically predicted expression, splicing, proteins, metabolites and methylation.",
                          no_args_is_help=True)
variants_app = typer.Typer(help="Participant carrier variants, consequence counts and gene burden.",
                           no_args_is_help=True)
hla_app = typer.Typer(help="Optional HLA imputation with HIBAG.", no_args_is_help=True)

PREDICTED = [Modality.EXPRESSION, Modality.SPLICING, Modality.PROTEOMICS, Modality.METABOLOMICS, Modality.METHYLATION]


def _execute(ids: set[str] | None, kinds: set[str] | None, cores: int | None, force: bool) -> None:
    from efgpp.setup.lock import write_lock
    from efgpp.workflow.executor import prepare

    project = load_project()
    outcomes = prepare(project, kinds=kinds, ids=ids, cores=cores, force=force, console=console)
    write_lock(project)
    rows = [[o.step_id, o.status, f"{o.seconds:.1f}s" if o.seconds else "", o.detail or ""] for o in outcomes]
    console.print(table("Summary", ["step", "status", "time", "detail"], rows))
    if any(o.status == "failed" for o in outcomes):
        raise typer.Exit(1)


# ------------------------------------------------------------------ predict
@predict_app.command("plan")
def predict_plan() -> None:
    """What would be predicted and what is missing (nothing is computed)."""
    from efgpp.data.predicted.engine import prediction_plan

    rows = prediction_plan(load_project())
    colour = {"READY": "green", "INSTALLED": "green", "PLANNED": "cyan", "DISABLED": "dim"}
    console.print(table("Genotype-derived molecular data: plan", ["item", "status", "detail", "required"],
                        [[r.item, f"[{colour.get(r.status, 'yellow')}]{r.status}[/]", r.detail, r.required]
                         for r in rows]))


@predict_app.command("enable")
def predict_enable(
    modality: str = typer.Argument(..., help="expression | splicing | proteomics | metabolomics | methylation"),
    tissue: list[str] = typer.Option([], "--tissue", help="PredictDB tissue (repeatable), or 'all'"),
    dataset: list[str] = typer.Option([], "--dataset", help="OmicsPred OPD id / PredictDB protein file (repeatable)"),
    engine: str | None = typer.Option(None, "--engine", help="metaxcan | plink_score | generic_weights | mimosa"),
    provider: str | None = typer.Option(None, "--provider", help="predictdb | omicspred | mimosa"),
    genotype: str | None = typer.Option(None, "--genotype", help="genotype source id (default: the first)"),
    min_coverage: float | None = typer.Option(None, "--min-coverage", help="minimum model variant coverage (0-1)"),
    off: bool = typer.Option(False, "--off", help="switch this modality off"),
) -> None:
    """Switch on genetically predicted data for one modality (runs in `data prepare` / `data predict`)."""
    project = load_project()
    if modality not in {m.value for m in PREDICTED}:
        console.print(f"[red]✗ unknown modality {modality!r}[/] (choose from {', '.join(m.value for m in PREDICTED)})")
        raise typer.Exit(2)
    cfg = getattr(project.data.predicted, modality)
    cfg.enabled = not off
    if not off:
        cfg.tissues = list(dict.fromkeys([*cfg.tissues, *tissue]))
        cfg.datasets = list(dict.fromkeys([*cfg.datasets, *dataset]))
        for attr, value in (("engine", engine), ("provider", provider), ("genotype_artifact", genotype),
                            ("minimum_variant_coverage", min_coverage)):
            if value is not None:
                setattr(cfg, attr, value)
    try:
        project.data.predicted.model_validate(project.data.predicted.model_dump())
    except Exception as exc:  # noqa: BLE001
        console.print(f"[red]✗ {exc}[/]")
        raise typer.Exit(2) from exc
    project.save_data_config()
    what = ", ".join(cfg.tissues or cfg.datasets) or "(no tissue/dataset yet)"
    console.print(f"[green]✓[/] predicted {modality} {'off' if off else f'on: {what} ({cfg.provider}, {cfg.engine})'}")
    console.print("  next: efgpp data predict plan")


def _predict_command(modality: Modality):  # type: ignore[no-untyped-def]
    def command(cores: int = typer.Option(None, "--cores"),
                force: bool = typer.Option(False, "--force", help="recompute up-to-date units")) -> None:
        _execute(None, {f"predict_{modality.value}"}, cores, force)

    command.__doc__ = f"Run every enabled genetically predicted {modality.value} unit (restartable)."
    return command


for _m in PREDICTED:
    predict_app.command(_m.value)(_predict_command(_m))


@predict_app.command("all")
def predict_all(cores: int = typer.Option(None, "--cores"), force: bool = typer.Option(False, "--force")) -> None:
    """Every enabled modality; a failing modality never touches finished outputs. Updates availability and report."""
    _execute(None, {"predict", "availability", "report"}, cores, force)


# ------------------------------------------------------------------ participant variants
@variants_app.command("enable")
def variants_enable(
    spliceai: bool = typer.Option(False, "--spliceai", help="also run SpliceAI (needs `efgpp setup toolkit spliceai`)"),
    all_genotypes: bool = typer.Option(False, "--all-genotypes", help="store non-carrier rows too (small cohorts)"),
    genotype: str | None = typer.Option(None, "--genotype", help="only this genotype source"),
    off: bool = typer.Option(False, "--off"),
) -> None:
    """Build participant carrier variants, annotate them (VEP, ClinVar, gnomAD, AlphaMissense) and count them."""
    project = load_project()
    cfg = project.data.participant_variants
    cfg.enabled = not off
    if not off:
        cfg.carrier_only = not all_genotypes
        cfg.annotations.spliceai = spliceai or cfg.annotations.spliceai
        if genotype:
            cfg.genotype_artifact = genotype
    project.save_data_config()
    console.print(f"[green]✓[/] participant variants {'off' if off else 'on'}")
    if not off:
        console.print("  resources used when installed: genome (FASTA, REF check), vep, clinvar, alphamissense, "
                      "gnomad (path), spliceai\n  next: efgpp data variants run")


@variants_app.command("run")
def variants_run(cores: int = typer.Option(None, "--cores"), force: bool = typer.Option(False, "--force")) -> None:
    """Carrier table -> annotations -> consequence counts and gene burden."""
    _execute(None, {"participant_variants", "annotate_", "consequences"}, cores, force)


# ------------------------------------------------------------------ HLA
@hla_app.command("enable")
def hla_enable(classifier: str = typer.Option(..., "--classifier", help="installed HIBAG classifier version"),
               loci: str | None = typer.Option(None, "--loci", help="comma-separated, e.g. A,B,C,DRB1"),
               off: bool = typer.Option(False, "--off")) -> None:
    """Switch on HLA imputation with an installed HIBAG classifier."""
    project = load_project()
    project.data.hla.enabled = not off
    project.data.hla.classifier = classifier
    if loci:
        project.data.hla.loci = [x.strip() for x in loci.split(",") if x.strip()]
    project.save_data_config()
    console.print(f"[green]✓[/] HLA imputation {'off' if off else f'on with {classifier}'}")


@hla_app.command("run")
def hla_run(cores: int = typer.Option(None, "--cores"), force: bool = typer.Option(False, "--force")) -> None:
    """Run HIBAG on every genotype source (classifier build must match; see `efgpp data liftover`)."""
    _execute(None, {"hla"}, cores, force)
