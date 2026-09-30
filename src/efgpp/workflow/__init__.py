"""Workflow execution: EFGPP plans WHAT must happen and runs it with its built-in executor;
`efgpp export slurm` writes explicit SLURM scripts for clusters."""

from __future__ import annotations

from efgpp.workflow.templates import write_workflow_templates

__all__ = ["write_workflow_templates"]
