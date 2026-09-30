"""Validation results that collect every problem instead of stopping at the first."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from efgpp.constants import FAIL, INFO, PASS, WARN

SYMBOLS = {PASS: "[green]✓[/]", FAIL: "[red]✗[/]", WARN: "[yellow]![/]", INFO: "[blue]i[/]"}


@dataclass
class Check:
    status: str
    message: str
    count: int | None = None
    examples: list[str] = field(default_factory=list)


@dataclass
class ValidationReport:
    subject: str  # e.g. "PHENOTYPE PH001 (trait_a)"
    source_id: str | None = None
    checks: list[Check] = field(default_factory=list)

    def ok(self, message: str) -> None:
        self.checks.append(Check(PASS, message))

    def fail(self, message: str, count: int | None = None, examples: list[Any] | None = None) -> None:
        self.checks.append(Check(FAIL, message, count, [str(e) for e in (examples or [])[:5]]))

    def warn(self, message: str, count: int | None = None, examples: list[Any] | None = None) -> None:
        self.checks.append(Check(WARN, message, count, [str(e) for e in (examples or [])[:5]]))

    def info(self, message: str) -> None:
        self.checks.append(Check(INFO, message))

    def extend(self, other: ValidationReport) -> None:
        self.checks.extend(other.checks)

    @property
    def passed(self) -> bool:
        return not any(c.status == FAIL for c in self.checks)

    @property
    def n_fail(self) -> int:
        return sum(c.status == FAIL for c in self.checks)

    @property
    def n_warn(self) -> int:
        return sum(c.status == WARN for c in self.checks)

    def render_lines(self) -> list[str]:
        lines = [f"[bold]{self.subject}[/]", ""]
        for c in self.checks:
            text = c.message
            if c.examples:
                text += f"  [dim](e.g. {', '.join(c.examples)})[/]"
            lines.append(f"{SYMBOLS[c.status]} {text}")
        return lines

    def to_dict(self) -> dict[str, Any]:
        return {"subject": self.subject, "source_id": self.source_id, "passed": self.passed,
                "checks": [asdict(c) for c in self.checks]}

    def save(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
        return path


def plural(n: int, word: str) -> str:
    return f"{n:,} {word}" + ("" if n == 1 else "s")
