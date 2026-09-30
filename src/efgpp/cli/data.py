"""`efgpp data ...` commands."""

from __future__ import annotations

from pathlib import Path

import typer

from efgpp.cli.common import console, load_project, table, user_path
from efgpp.config.data import (
    ClinicalSource,
    GenotypeSource,
    OmicsSource,
    SampleIdSpec,
    TimelineSpec,
)
from efgpp.config.data import CovariateSource as CovSource
from efgpp.constants import Modality, Origin, StorageMode
from efgpp.project import Project

app = typer.Typer(help="Register, validate, prepare and freeze participant-level data.", no_args_is_help=True)
add_app = typer.Typer(help="Add a participant-level source to data.yaml and register it.", no_args_is_help=True)
app.add_typer(add_app, name="add")
gwas_app = typer.Typer(help="GWAS summary statistics: list and (re)run GWASLab.", no_args_is_help=True)
app.add_typer(gwas_app, name="gwas")


@gwas_app.command("list")
def gwas_list(phenotype: str | None = typer.Option(None, "--phenotype", help="only GWAS for this phenotype")) -> None:
    """GWAS in this project, optionally only those for one phenotype (--phenotype <name or id>)."""
    from efgpp.data.artifacts import ArtifactStore
    from efgpp.data.registry import Registry

    project = load_project()
    with Registry.open(project) as reg:
        done = {a.source_id: a for a in ArtifactStore(reg).find(artifact_type="gwas_sumstats")}
    rows = []
    for g in project.data.gwas:
        targets = g.phenotypes or [g.trait]
        if phenotype and phenotype not in targets and phenotype != g.trait:
            continue
        art = done.get(g.id)
        rows.append([g.id, g.trait, ", ".join(targets), g.ancestry or "", g.n_cases or "", g.n_controls or "",
                     art.genome_build if art else "not run", f"{art.feature_count:,}" if art and art.feature_count else ""])
    console.print(table("GWAS" + (f" for {phenotype}" if phenotype else ""),
                        ["id", "trait", "for phenotype", "ancestry", "cases", "controls", "build", "variants"], rows))


@gwas_app.command("run")
def gwas_run(gwas_id: str = typer.Argument(...), force: bool = typer.Option(False, "--force")) -> None:
    """(Re)process one GWAS with GWASLab now."""
    from efgpp.data.gwas import run_gwas

    project = load_project()
    out = run_gwas(project, gwas_id, force=force)
    console.print(f"[green]✓[/] {gwas_id}: {out}")


def _register(project, source_id: str) -> None:  # type: ignore[no-untyped-def]
    from efgpp.data.manager import DataManager
    from efgpp.data.registry import Registry

    with Registry.open(project) as reg:
        art = DataManager(project).adapter(reg, source_id).register()
    notes = "; ".join(art.metadata.get("notes", []))
    console.print(f"[green]✓[/] {source_id} registered as {art.artifact_id} "
                  f"(storage: {art.storage_mode}{'; ' + notes if notes else ''})")


def _timeline(event_column: str | None, time_column: str | None) -> TimelineSpec | None:
    if not (event_column or time_column):
        return None
    return TimelineSpec(event_column=event_column, time_column=time_column)


@add_app.command("genotype")
def add_genotype(
    path: str = typer.Option(..., "--path", help="PLINK prefix, or BGEN/VCF/BCF file"),
    format: str = typer.Option("auto", "--format", help="pgen | bed | bgen | vcf | bcf | auto"),
    build: str = typer.Option("auto", "--build", help="GRCh37 | GRCh38 | auto (detected, then lifted to GRCh38)"),
    mode: StorageMode = typer.Option(StorageMode.REFERENCE, "--mode"),
    source_id: str | None = typer.Option(None, "--id"),
    origin: Origin = typer.Option(Origin.OBSERVED, "--origin"),
    sample_id_mode: str = typer.Option("auto", "--sample-id-mode", help="iid | fid_iid | auto"),
) -> None:
    """Add a genotype dataset (referenced in place; format and sample IDs inferred by modules/genotype.py)."""
    from efgpp.modules import load

    project = load_project()
    guess = load(project, "genotype").infer_genotype(project.resolve(user_path(project, path)))
    fmt = guess.format if format == "auto" else format
    id_mode = guess.sample_id_mode if sample_id_mode == "auto" else sample_id_mode
    console.print(f"[dim]genotype: {fmt}, {guess.samples:,} samples, sample IDs: {id_mode}"
                  f"{' (' + guess.note + ')' if guess.note else ''}[/]")
    sid = source_id or project.data.next_id("GENO")
    project.data.observed.genotype.append(GenotypeSource(
        id=sid, path=user_path(project, path), format=fmt, genome_build=build, mode=mode, origin=origin,
        sample_id=SampleIdSpec(mode=id_mode)))  # type: ignore[arg-type]
    project.save_data_config()
    _register(project, sid)


@add_app.command("covariates", context_settings={"allow_extra_args": True, "ignore_unknown_options": False})
def add_covariates(
    ctx: typer.Context,
    path: str = typer.Option(..., "--path"),
    id_column: str | None = typer.Option(None, "--id-column", help="default: inferred (IID, participant_id, eid, ...)"),
    columns: str | None = typer.Option(None, "--columns",
                                       help="covariate columns (space or comma separated); default: all except IDs"),
    categorical: str = typer.Option("", "--categorical", help="extra categorical columns (others are inferred)"),
    event_column: str | None = typer.Option(None, "--event-column"),
    mode: StorageMode = typer.Option(StorageMode.AUTO, "--mode"),
    source_id: str | None = typer.Option(None, "--id"),
) -> None:
    """Add a covariate table. ID column, covariates, roles (sex, age, PCs, batch) and categorical
    columns are inferred by modules/covariates.py; column names are kept as they are."""
    from efgpp.data.io import read_table
    from efgpp.modules import load

    project = load_project()
    df = read_table(project.resolve(user_path(project, path)))
    wanted = [c for part in [columns or "", *ctx.args] for c in part.replace(",", " ").split()] or None
    guess = load(project, "covariates").infer_covariates(df, id_column, wanted)
    if guess.id_column is None:
        console.print(f"[red]no participant ID column found; pass --id-column (columns: {', '.join(df.columns)})[/]")
        raise typer.Exit(2)
    missing = [c for c, why in guess.skipped.items() if why == "not in the file"]
    if missing:
        console.print(f"[red]columns not in the file: {', '.join(missing)}[/]; available: {', '.join(df.columns)}")
        raise typer.Exit(2)
    extra_cat = [c for c in categorical.split(",") if c]
    cats = list(dict.fromkeys([*guess.categorical, *extra_cat]))
    console.print(table(f"Covariates (ID column: {guess.id_column})", ["column", "role", "type"],
                        [[c, guess.roles[c], "categorical" if c in cats else "numeric"] for c in guess.variables]))
    sid = source_id or project.data.next_id("COV")
    project.data.observed.covariates.append(CovSource(
        id=sid, path=user_path(project, path), participant_id_column=guess.id_column, variables=guess.variables,
        categorical=cats, roles=guess.roles, mode=mode, timeline=_timeline(event_column, None)))
    project.save_data_config()
    _register(project, sid)


def _add_omics(modality: Modality, prefix: str):  # type: ignore[no-untyped-def]
    def command(
        path: str = typer.Option(..., "--path"),
        origin: Origin = typer.Option(Origin.OBSERVED, "--origin"),
        id_column: str = typer.Option("participant_id", "--id-column"),
        tissue: str | None = typer.Option(None, "--tissue"),
        event_column: str | None = typer.Option(None, "--event-column"),
        time_column: str | None = typer.Option(None, "--time-column"),
        biospecimen_column: str | None = typer.Option(None, "--biospecimen-column"),
        batch_column: str | None = typer.Option(None, "--batch-column"),
        feature_id_system: str | None = typer.Option(None, "--feature-id-system"),
        measurement_type: str | None = typer.Option(None, "--measurement-type"),
        normalization: str | None = typer.Option(None, "--normalization"),
        platform: str | None = typer.Option(None, "--platform"),
        units: str | None = typer.Option(None, "--units"),
        genome_build: str | None = typer.Option(None, "--build"),
        orientation: str = typer.Option("samples_by_features", "--orientation"),
        feature_id_column: str | None = typer.Option(None, "--feature-id-column"),
        mode: StorageMode = typer.Option(StorageMode.AUTO, "--mode"),
        source_id: str | None = typer.Option(None, "--id"),
    ) -> None:
        project = load_project()
        if origin == Origin.PREDICTED:
            console.print("[yellow]note:[/] predicted modalities are usually generated by `efgpp data prepare` "
                          "(predicted: section of data.yaml); registering an external predicted matrix.")
        sid = source_id or project.data.next_id(prefix)
        src = OmicsSource(
            id=sid, path=user_path(project, path), origin=origin, participant_id_column=id_column, tissue=tissue,
            timeline=_timeline(event_column, time_column), biospecimen_column=biospecimen_column,
            batch_column=batch_column, feature_id_system=feature_id_system, measurement_type=measurement_type,
            normalization=normalization, platform=platform, units=units, genome_build=genome_build,
            orientation=orientation, feature_id_column=feature_id_column, mode=mode)  # type: ignore[arg-type]
        getattr(project.data.observed, modality.value).append(src)
        project.save_data_config()
        _register(project, sid)

    command.__doc__ = f"Add a processed {modality.value} matrix (CSV/TSV/Parquet/AnnData/Zarr)."
    return command


for _m, _p in ((Modality.EXPRESSION, "RNA"), (Modality.METHYLATION, "METH"), (Modality.PROTEOMICS, "PROT"),
               (Modality.METABOLOMICS, "METAB")):
    add_app.command(_m.value)(_add_omics(_m, _p))


@add_app.command("gwas", context_settings={"allow_extra_args": True})
def add_gwas(
    ctx: typer.Context,
    path: str = typer.Option(..., "--path", help="GWAS summary statistics file (any text/gz format)"),
    trait: str | None = typer.Option(None, "--trait", help="trait this GWAS studied (default: file name)"),
    phenotype: list[str] = typer.Option([], "--phenotype",
                                        help="project phenotype(s) this GWAS is for (repeatable); default: --trait"),
    ancestry: str | None = typer.Option(None, "--ancestry", help="e.g. European, EUR, East Asian, multi-ancestry"),
    n_cases: int | None = typer.Option(None, "--n-cases"),
    n_controls: int | None = typer.Option(None, "--n-controls"),
    n: int | None = typer.Option(None, "--n", help="total sample size when the file has no N column"),
    study: str | None = typer.Option(None, "--study", help="consortium / publication / GWAS Catalog accession"),
    build: str = typer.Option("auto", "--build", help="GRCh37 | GRCh38 | auto (GWASLab infer_build)"),
    fmt: str = typer.Option("auto", "--fmt", help="GWASLab format (ssf, gwascatalog, plink2, regenie, ...) or auto"),
    col: list[str] = typer.Option([], "--col", help="override one mapping: GWASLab keyword=column, e.g. pos=BP"),
    run: bool = typer.Option(True, "--run/--no-run", help="process with GWASLab now (default) or at `data prepare`"),
    source_id: str | None = typer.Option(None, "--id"),
) -> None:
    """Add GWAS summary statistics. Columns are recognised by modules/gwas.py, GWASLab checks the
    file, infers the build and lifts it to GRCh38; the original column names are kept."""
    import gzip

    from efgpp.config.data import GwasSource
    from efgpp.modules import load

    project = load_project()
    file_path = project.resolve(user_path(project, path))
    opener = gzip.open if file_path.name.endswith((".gz", ".bgz")) else open
    with opener(file_path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[operator]
        header = fh.readline().replace(",", " ").split()
    gwas_module = load(project, "gwas")
    columns = gwas_module.infer_gwas_columns(header)
    for item in [*col, *ctx.args]:
        key, _, value = item.partition("=")
        if not value:
            console.print(f"[red]--col expects keyword=column, got {item!r}[/]")
            raise typer.Exit(2)
        columns[key.strip()] = value.strip()
    unknown = [v for v in columns.values() if v not in header]
    if unknown:
        console.print(f"[red]columns not in the file: {unknown}[/]; header: {' '.join(header)}")
        raise typer.Exit(2)
    missing = gwas_module.missing_essentials(columns)
    if missing:
        console.print(f"[red]could not find: {', '.join(missing)}[/] in {' '.join(header)}; add --col keyword=column")
        raise typer.Exit(2)
    trait = trait or file_path.name.split(".")[0]
    targets = list(phenotype)
    if not targets and any(trait in (p.id, p.name) for p in project.data.observed.phenotypes):
        targets = [trait]
    unmapped = [c for c in header if c not in columns.values()]
    console.print(table("GWAS columns (GWASLab keyword -> column)", ["keyword", "column"],
                        [[k, v] for k, v in columns.items()] + [["(kept as is)", c] for c in unmapped]))
    gid = source_id or project.data.next_id("GWAS")
    project.data.gwas.append(GwasSource(
        id=gid, path=user_path(project, path), trait=trait, phenotypes=targets, ancestry=ancestry, n_cases=n_cases,
        n_controls=n_controls, n=n, study=study, build=build, fmt=fmt, columns=columns))
    project.save_data_config()
    details = ", ".join(x for x in (ancestry, f"{n_cases:,} cases" if n_cases else "",
                                    f"{n_controls:,} controls" if n_controls else "") if x)
    console.print(f"[green]✓[/] {gid} added: {trait}{' (' + details + ')' if details else ''}"
                  f"{'; for ' + ', '.join(targets) if targets else ''}")
    if not run:
        console.print(f"[dim]processed by `efgpp data prepare` (step gwas.{gid})[/]")
        return
    from efgpp.data.gwas import run_gwas
    from efgpp.setup.tools import available

    if not available(project, "gwaslab"):
        console.print("[yellow]GWASLab is not installed: run `efgpp setup tools gwaslab`, then "
                      f"`efgpp data gwas run {gid}` (or `efgpp data prepare`)[/]")
        return
    console.print("[dim]running GWASLab (basic_check, infer_build, liftover to "
                  f"{project.config.defaults.target_build})...[/]")
    out = run_gwas(Project.load(project.root), gid)
    console.print(f"[green]✓[/] {gid}: {out.get('rows'):,} variants, build {out.get('build_detected')}"
                  f"{', lifted ' + out['lifted'] if out.get('lifted') else ''} -> {out.get('path')}")


@add_app.command("clinical")
def add_clinical(
    path: str = typer.Option(..., "--path"),
    id_column: str = typer.Option("participant_id", "--id-column"),
    layout: str = typer.Option("wide", "--layout", help="wide | long"),
    variable_column: str | None = typer.Option(None, "--variable-column"),
    value_column: str | None = typer.Option(None, "--value-column"),
    unit_column: str | None = typer.Option(None, "--unit-column"),
    event_column: str | None = typer.Option(None, "--event-column"),
    time_column: str | None = typer.Option(None, "--time-column"),
    mode: StorageMode = typer.Option(StorageMode.AUTO, "--mode"),
    source_id: str | None = typer.Option(None, "--id"),
) -> None:
    """Add clinical/laboratory measurements (wide or long layout)."""
    project = load_project()
    sid = source_id or project.data.next_id("CLIN")
    project.data.observed.clinical.append(ClinicalSource(
        id=sid, path=user_path(project, path), participant_id_column=id_column, layout=layout,  # type: ignore[arg-type]
        variable_column=variable_column, value_column=value_column, unit_column=unit_column,
        timeline=_timeline(event_column, time_column), mode=mode))
    project.save_data_config()
    _register(project, sid)


@app.command()
def remove(source_id: str = typer.Argument(..., help="source id, e.g. COV001")) -> None:
    """Remove a source from data.yaml (registered artifacts are kept for provenance)."""
    project = load_project()
    if any(g.id == source_id for g in project.data.gwas):
        project.data.gwas = [g for g in project.data.gwas if g.id != source_id]
        project.save_data_config()
        console.print(f"[green]✓[/] removed {source_id} (gwas) from data.yaml")
        return
    try:
        modality, _ = project.data.get_source(source_id)
    except KeyError as exc:
        console.print(f"[red]{exc}[/]")
        raise typer.Exit(2) from exc
    obs = project.data.observed
    for name in type(obs).model_fields:
        items = getattr(obs, name)
        setattr(obs, name, [s for s in items if s.id != source_id])
    project.save_data_config()
    console.print(f"[green]✓[/] removed {source_id} ({modality.value}) from data.yaml")


@app.command()
def inspect() -> None:
    """What exists: every configured source with format, size, participants and features."""
    from efgpp.data.manager import DataManager

    project = load_project()
    rows = []
    for i in DataManager(project).inspect():
        rows.append([i.get("source_id"), i.get("modality"), i.get("origin"), i.get("format"),
                     "yes" if i.get("exists") else "[red]no[/]", i.get("samples") or i.get("participants") or i.get("rows"),
                     i.get("variants") or i.get("features") or i.get("columns"),
                     i.get("genome_build") or i.get("tissue") or ""])
    console.print(table("EFGPP DATA INSPECT", ["source", "modality", "origin", "format", "exists", "participants",
                                               "features", "build/tissue"], rows))


@app.command()
def validate() -> None:
    """Validate every source and cross-source consistency (collects all errors)."""
    from efgpp.data.manager import DataManager

    project = load_project()
    DataManager(project).register_all()
    reports = DataManager(project).validate()
    failed = 0
    for rep in reports:
        console.print()
        for line in rep.render_lines():
            console.print(line)
        failed += 0 if rep.passed else 1
    console.print()
    if failed:
        console.print(f"[red]{failed} validation report(s) with errors[/]")
        raise typer.Exit(1)
    console.print("[green]all sources valid[/]")


@app.command()
def plan() -> None:
    """What can be generated now, what cannot, and why. Writes workflow/plan.json."""
    from efgpp.data.plan import build_plan, write_plan

    project = load_project()
    steps = build_plan(project)
    write_plan(project, steps)
    rows = [["[green]run[/]" if s.enabled else "[yellow]blocked[/]", s.id, s.description if s.enabled else s.reason,
             ",".join(s.tools) or "-", s.env] for s in steps]
    console.print(table("EFGPP DATA PLAN", ["", "step", "what / why not", "tools", "env"], rows))
    console.print(f"[dim]{project.relative(project.path('workflow', 'plan.json'))} written[/]")


def _run(kinds: set[str] | None, cores: int | None, force: bool, dry_run: bool) -> None:
    from efgpp.setup.lock import write_lock
    from efgpp.workflow.executor import prepare

    project = load_project()
    outcomes = prepare(project, kinds=kinds, cores=cores, force=force, dry_run=dry_run, console=console)
    write_lock(project)
    rows = [[o.step_id, o.status, f"{o.seconds:.1f}s" if o.seconds else "", o.detail or ""] for o in outcomes]
    console.print(table("Summary", ["step", "status", "time", "detail"], rows))
    if any(o.status == "failed" for o in outcomes):
        raise typer.Exit(1)


@app.command()
def prepare(cores: int = typer.Option(None, "--cores"),
            force: bool = typer.Option(False, "--force", help="rerun up-to-date steps"),
            dry_run: bool = typer.Option(False, "--dry-run")) -> None:
    """Run the full data lifecycle (register -> ... -> availability -> report)."""
    _run(None, cores, force, dry_run)


@app.command()
def qc(cores: int = typer.Option(None, "--cores"), force: bool = typer.Option(False, "--force")) -> None:
    """Run QC steps only (and what they depend on)."""
    _run({"qc", "genotype_qc"}, cores, force, False)


@app.command()
def derive(cores: int = typer.Option(None, "--cores"), force: bool = typer.Option(False, "--force")) -> None:
    """Run phenotype-independent derivations, annotation, GWAS and prediction steps."""
    _run({"pca", "roh", "ancestry", "annotate_", "predict", "gwas"}, cores, force, False)


@app.command()
def step(step_id: str = typer.Argument(...), threads: int = typer.Option(1, "--threads"),
         marker: str | None = typer.Option(None, "--marker")) -> None:
    """Run a single plan step (used by the SLURM scripts from `efgpp export slurm`)."""
    from efgpp.data.plan import build_plan, run_step
    from efgpp.workflow.executor import write_marker

    project = load_project()
    steps = {s.id: s for s in build_plan(project)}
    if step_id not in steps:
        console.print(f"[red]unknown step {step_id}[/]")
        raise typer.Exit(2)
    s = steps[step_id]
    result = run_step(project, s, threads=threads)
    if marker:
        write_marker(project, s)
    console.print(result)


@app.command()
def availability() -> None:
    """Which participants have which modalities; intersections."""
    from efgpp.data import availability as av_mod

    project = load_project()
    av = av_mod.build(project)
    console.print(f"\n[bold]Participants: {av.n_participants:,}[/]\n")
    for origin, entries in av_mod.summary_rows(av):
        console.print(f"[bold]{origin}[/]")
        for col, n in entries:
            console.print(f"  {col:<34}{n:>10,}")
        console.print()
    for a in av.participant_level_free:
        console.print(f"  {a['artifact_name']:<34}{'available':>10}")
    inter = av.key_intersections()
    if inter:
        console.print(table("Key intersections", ["modalities", "participants"],
                            [[i["modalities"], i["participants"]] for i in inter]))
    console.print("[dim]registry/availability.parquet written[/]")


@app.command()
def report() -> None:
    """Generate reports/data/index.html."""
    from efgpp.reporting.data_report import build_report

    project = load_project()
    out = build_report(project)
    console.print(f"[green]✓[/] {out}")


@app.command()
def freeze(name: str = typer.Option(..., "--name"),
           include_failed: bool = typer.Option(False, "--include-failed")) -> None:
    """Create an immutable DataSnapshot for the Representation layer."""
    from efgpp.data.snapshots import freeze as do_freeze
    from efgpp.setup.lock import write_lock

    project = load_project()
    write_lock(project)
    try:
        path = do_freeze(project, name, include_failed=include_failed)
    except (FileExistsError, RuntimeError) as exc:
        console.print(f"[red]{exc}[/]")
        raise typer.Exit(1) from exc
    console.print(f"[green]✓[/] snapshot {name} -> {project.relative(path)}")


@app.command()
def snapshots() -> None:
    """List frozen snapshots."""
    from efgpp.data.snapshots import list_snapshots

    project = load_project()
    rows = [[s["snapshot_id"], s["name"], s["created_at"], s["path"]] for s in list_snapshots(project)]
    console.print(table("Snapshots", ["id", "name", "created", "path"], rows))


@app.command()
def verify(name: str = typer.Argument(..., help="snapshot name")) -> None:
    """Re-check every file hash recorded in a snapshot."""
    from efgpp.data.snapshots import DataSnapshot

    project = load_project()
    problems = DataSnapshot.load(project, name).verify()
    if problems:
        for k, v in problems.items():
            console.print(f"[red]✗[/] {k}: {v}")
        raise typer.Exit(1)
    console.print(f"[green]✓[/] snapshot {name} intact")


@app.command()
def migrate(source: Path = typer.Option(..., "--source", exists=True, file_okay=False),
            profile: str = typer.Option("legacy-efgpp", "--profile"),
            mode: StorageMode = typer.Option(StorageMode.REFERENCE, "--mode"),
            dry_run: bool = typer.Option(False, "--dry-run"), apply_: bool = typer.Option(False, "--apply")) -> None:
    """Migrate an old project (never modifies the source)."""
    from efgpp.data import migrate as mig

    project = load_project()
    if dry_run == apply_:
        console.print("[red]choose exactly one of --dry-run or --apply[/]")
        raise typer.Exit(2)
    if dry_run:
        items = mig.scan(source, profile)
        path = mig.write_plan(project, source, profile, mode, items)
        counts: dict[str, int] = {}
        for it in items:
            counts[it.kind] = counts.get(it.kind, 0) + 1
        console.print(table("Migration plan", ["kind", "files"], [[k, v] for k, v in sorted(counts.items())]))
        console.print(f"[green]✓[/] {project.relative(path)} written; review it, then run with --apply")
        return
    done = mig.apply(project, source, profile, mode)
    console.print(table("Migrated", ["kind", "registered"], [[k, ", ".join(v)] for k, v in done.items()]))


@app.command()
def liftover(artifact_id: str = typer.Argument(...), to: str = typer.Option(..., "--to", help="GRCh37 | GRCh38"),
             threads: int = typer.Option(1, "--threads")) -> None:
    """Lift a genotype artifact to another build (creates a new artifact)."""
    from efgpp.data.genotype.harmonization import run_liftover

    project = load_project()
    new = run_liftover(project, artifact_id, to, threads=threads)
    console.print(f"[green]✓[/] {artifact_id} -> {new}")


@app.command()
def simulate() -> None:
    """Generate the simulated modalities enabled under `simulation:` in data.yaml."""
    from efgpp.data.simulate import run_simulation

    project = load_project()
    added = run_simulation(project)
    console.print(f"[green]✓[/] simulated sources: {', '.join(added) or 'none (enable simulation.genotype)'}")
