"""Running a data plan: through Snakemake when available, otherwise the built-in executor.

Both paths share the same step markers (work/steps/<id>.done, containing the step's input
fingerprint), so a project can move between them freely.
"""

from __future__ import annotations

import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rich.console import Console

from efgpp.data.plan import Step, build_plan, fingerprint, run_step, select, write_plan
from efgpp.project import Project
from efgpp.setup.tools import available


@dataclass
class StepOutcome:
    step_id: str
    status: str  # done | skipped | up-to-date | failed | blocked | disabled
    seconds: float = 0.0
    detail: str | None = None
    result: dict[str, Any] = field(default_factory=dict)


def marker_path(project: Project, step: Step) -> Path:
    return project.root / step.marker


def write_marker(project: Project, step: Step) -> None:
    p = marker_path(project, step)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(fingerprint(project, step), encoding="utf-8")


def is_current(project: Project, step: Step) -> bool:
    p = marker_path(project, step)
    return p.exists() and p.read_text(encoding="utf-8").strip() == fingerprint(project, step)


def choose_engine(project: Project, requested: str | None) -> str:
    engine = requested or project.config.execution.engine
    if engine == "auto":
        return "snakemake" if available(project, "snakemake") else "builtin"
    return engine


def run_builtin(project: Project, steps: list[Step], *, cores: int, force: bool = False,
                console: Console | None = None) -> list[StepOutcome]:
    console = console or Console()
    outcomes: dict[str, StepOutcome] = {}
    rerun: set[str] = set()
    for step in steps:  # plan order is topological
        if not step.enabled:
            outcomes[step.id] = StepOutcome(step.id, "disabled", detail=step.reason)
            continue
        failed_needs = [n for n in step.needs if n in outcomes and outcomes[n].status in ("failed", "blocked")]
        if failed_needs and step.kind in ("availability", "report"):
            console.print(f"[yellow]![/] {step.id}: running despite failed {', '.join(failed_needs)}")
        elif failed_needs:
            outcomes[step.id] = StepOutcome(step.id, "blocked", detail=f"upstream failed: {', '.join(failed_needs)}")
            console.print(f"[yellow]-[/] {step.id}: blocked ({', '.join(failed_needs)} failed)")
            continue
        upstream_rerun = any(n in rerun for n in step.needs)
        if not (force or step.always_run or upstream_rerun) and is_current(project, step):
            outcomes[step.id] = StepOutcome(step.id, "up-to-date")
            console.print(f"[dim]= {step.id}: up to date[/]")
            continue
        console.print(f"[cyan]>[/] {step.id}: {step.description}")
        t0 = time.monotonic()
        try:
            result = run_step(project, step, threads=min(step.threads, cores))
        except Exception as exc:  # noqa: BLE001 - keep going with independent steps
            outcomes[step.id] = StepOutcome(step.id, "failed", time.monotonic() - t0, f"{exc}",
                                            {"traceback": traceback.format_exc()})
            console.print(f"[red]x {step.id} failed:[/] {exc}")
            continue
        write_marker(project, step)
        if not step.always_run:
            rerun.add(step.id)
        outcomes[step.id] = StepOutcome(step.id, "done", time.monotonic() - t0, result=result)
        if step.kind == "validate" and result.get("failed"):
            console.print(f"[yellow]! validation failed for {', '.join(result['failed'])}; "
                          "their downstream steps will fail until fixed[/]")
    return list(outcomes.values())


def prepare(project: Project, *, kinds: set[str] | None = None, cores: int | None = None,
            executor: str | None = None, engine: str | None = None, force: bool = False,
            dry_run: bool = False, console: Console | None = None) -> list[StepOutcome]:
    """Plan and execute. `kinds` restricts to some step kinds (plus their dependencies)."""
    console = console or Console()
    cores = cores or project.config.execution.local_cores
    ensure_simulation(project, console)
    steps = build_plan(project)
    write_plan(project, steps)
    chosen = select(steps, kinds) if kinds else steps
    selected_engine = choose_engine(project, "snakemake" if executor and executor != "local" else engine)
    if selected_engine == "snakemake":
        from efgpp.workflow.snakemake import run_snakemake, write_snakefile

        write_snakefile(project, steps)
        for s in chosen:  # invalidate stale or always-run markers so Snakemake reruns them
            if s.enabled and (force or s.always_run or not is_current(project, s)):
                marker_path(project, s).unlink(missing_ok=True)
        code = run_snakemake(project, cores=cores, executor=executor, dry_run=dry_run)
        status = "done" if code == 0 else "failed"
        return [StepOutcome("snakemake", status, detail=f"exit code {code}")]
    if dry_run:
        return [StepOutcome(s.id, "planned" if s.enabled else "disabled", detail=s.reason) for s in chosen]
    return run_builtin(project, chosen, cores=cores, force=force, console=console)


def ensure_simulation(project: Project, console: Console) -> None:
    """Generate enabled simulated modalities once, before planning (they add sources)."""
    from efgpp.constants import Origin

    sim = project.data.simulation
    have = {s.id for _, s in project.data.iter_sources() if s.origin == Origin.SIMULATED}
    if sim.genotype.enabled and "SIMGENO" not in have:
        from efgpp.data.simulate import run_simulation

        added = run_simulation(project)
        console.print(f"[cyan]>[/] simulate: generated {', '.join(added)}")
