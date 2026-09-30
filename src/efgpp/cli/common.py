"""Shared CLI state and console helpers."""

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from efgpp.project import Project, ProjectNotFoundError

console = Console(highlight=False)


class State:
    project_root: Path | None = None


state = State()


def load_project() -> Project:
    try:
        return Project.load(state.project_root)
    except ProjectNotFoundError as exc:
        console.print(f"[red]{exc}[/]")
        raise typer.Exit(2) from exc
    except Exception as exc:
        console.print(f"[red]configuration error:[/] {exc}")
        raise typer.Exit(2) from exc


def user_path(project: Project, path: str) -> str:
    """Store user-supplied paths relative to the project root when they are inside it."""
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = (Path.cwd() / p).resolve()
    return project.relative(p)


def table(title: str, columns: list[str], rows: Sequence[Sequence[object]]) -> Table:
    t = Table(title=title, title_justify="left", header_style="bold")
    for c in columns:
        t.add_column(c)
    for r in rows:
        t.add_row(*["" if v is None else f"{v:,}" if isinstance(v, int) and not isinstance(v, bool) else str(v) for v in r])
    return t


def ensure_utf8() -> None:
    """Windows consoles default to a legacy code page; ✓/✗ need UTF-8."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")  # type: ignore[union-attr]
        except (AttributeError, ValueError):
            pass
