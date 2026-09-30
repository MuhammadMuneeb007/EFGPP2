"""Workflow execution: EFGPP decides WHAT must happen; Snakemake (or the built-in
executor where Snakemake is unavailable) decides WHEN, in which order, in parallel,
with restarts and on which executor (local or SLURM)."""

from __future__ import annotations

from efgpp.workflow.templates import write_workflow_templates

__all__ = ["write_workflow_templates"]
