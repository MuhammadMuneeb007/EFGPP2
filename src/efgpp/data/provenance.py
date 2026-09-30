"""Provenance: every external command produces a run record and structured logs.

For each run EFGPP writes
    logs/RUN000123.jsonl   structured events (one JSON object per line)
    logs/RUN000123.log     human-readable combined stdout/stderr
    logs/RUN000123.yaml    the run record
and a row in the registry's tool_runs table.

The registry lock is *not* held while a tool runs, so parallel SLURM jobs do not
serialise on it.
"""

from __future__ import annotations

import json
import os
import platform
import shlex
import subprocess
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import psutil

from efgpp import __version__
from efgpp.config.models import write_yaml
from efgpp.data.registry import Registry, utcnow
from efgpp.project import Project
from efgpp.setup.tools import ResolvedTool, detect_version, resolve


class ToolRunError(RuntimeError):
    def __init__(self, record: RunRecord, tail: str) -> None:
        super().__init__(
            f"{record.tool} exited with code {record.exit_code} (run {record.run_id}); "
            f"see {record.log_path}\n{tail}"
        )
        self.record = record


@dataclass
class RunRecord:
    run_id: str
    step_id: str | None
    tool: str
    tool_version: str | None
    command: list[str]
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)
    started_at: datetime | None = None
    completed_at: datetime | None = None
    exit_code: int | None = None
    status: str = "running"
    environment: dict[str, Any] = field(default_factory=dict)
    resources: dict[str, Any] = field(default_factory=dict)
    log_path: str | None = None

    def to_yaml(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "step_id": self.step_id,
            "tool": {"name": self.tool, "version": self.tool_version},
            "command": shlex.join(self.command),
            "inputs": self.inputs,
            "outputs": self.outputs,
            "started": self.started_at.isoformat() if self.started_at else None,
            "completed": self.completed_at.isoformat() if self.completed_at else None,
            "exit_code": self.exit_code,
            "status": self.status,
            "environment": self.environment,
            "resources": self.resources,
        }


def _software_version(project: Project, tool: ResolvedTool) -> str | None:
    """Version from the registry cache, detecting (and caching) it on first use."""
    with Registry.open(project) as reg:
        v = reg.scalar(
            "SELECT version FROM software WHERE name = ? AND environment = ?", [tool.name, tool.env]
        )
        if v:
            return str(v)
    version = detect_version(tool)
    with Registry.open(project) as reg:
        reg.upsert("software", {
            "name": tool.name, "environment": tool.env, "version": version,
            "path": str(tool.path), "detected_at": utcnow(),
        })
    return version


class _PeakMemory(threading.Thread):
    def __init__(self, pid: int) -> None:
        super().__init__(daemon=True)
        self.pid = pid
        self.peak = 0
        self.cpu_seconds = 0.0
        self._halt = threading.Event()

    def run(self) -> None:
        try:
            proc = psutil.Process(self.pid)
        except psutil.Error:
            return
        while not self._halt.is_set():
            try:
                procs = [proc, *proc.children(recursive=True)]
                self.peak = max(self.peak, sum(p.memory_info().rss for p in procs))
                t = proc.cpu_times()
                self.cpu_seconds = max(self.cpu_seconds, t.user + t.system)
            except psutil.Error:
                break
            self._halt.wait(0.25)

    def stop(self) -> None:
        self._halt.set()


def run_tool(
    project: Project,
    tool_name: str,
    args: list[str],
    *,
    step_id: str | None = None,
    inputs: list[str] | None = None,
    outputs: list[str] | None = None,
    cwd: Path | None = None,
    check: bool = True,
    tool: ResolvedTool | None = None,
) -> RunRecord:
    """Run an external tool with full provenance capture."""
    tool = tool or resolve(project, tool_name)
    version = _software_version(project, tool)
    with Registry.open(project) as reg:
        run_id = reg.next_id("RUN")
    logs = project.logs_dir
    logs.mkdir(parents=True, exist_ok=True)
    log_path, jsonl_path = logs / f"{run_id}.log", logs / f"{run_id}.jsonl"
    command = tool.command(*args)
    record = RunRecord(
        run_id=run_id, step_id=step_id, tool=tool.name, tool_version=version, command=command,
        inputs=inputs or [], outputs=outputs or [], log_path=project.relative(log_path),
        environment={
            "efgpp_version": __version__, "env": tool.env, "executable": str(tool.path),
            "platform": platform.platform(), "hostname": platform.node(),
        },
    )

    def event(kind: str, **payload: Any) -> None:
        with open(jsonl_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps({"time": utcnow().isoformat(), "run_id": run_id, "event": kind, **payload}, default=str) + "\n")

    record.started_at = utcnow()
    event("start", tool=tool.name, version=version, command=command)
    t0 = time.monotonic()
    with open(log_path, "w", encoding="utf-8") as log:
        log.write(f"# {run_id} {tool.name} {version or ''}\n# $ {shlex.join(command)}\n")
        log.flush()
        proc = subprocess.Popen(
            command, stdout=log, stderr=subprocess.STDOUT, cwd=cwd, env=tool.environ()
        )
        monitor = _PeakMemory(proc.pid)
        monitor.start()
        exit_code = proc.wait()
        monitor.stop()
    record.completed_at = utcnow()
    record.exit_code = exit_code
    record.status = "completed" if exit_code == 0 else "failed"
    record.resources = {
        "wall_seconds": round(time.monotonic() - t0, 3),
        "cpu_seconds": round(monitor.cpu_seconds, 3),
        "peak_memory_mb": round(monitor.peak / 1024**2, 1),
        "cpu_count": os.cpu_count(),
    }
    event("end", exit_code=exit_code, resources=record.resources)
    _persist(project, record)
    if check and exit_code != 0:
        tail = "".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines(True)[-20:])
        raise ToolRunError(record, tail)
    return record


def record_python_step(
    project: Project,
    step_id: str,
    *,
    description: str,
    inputs: list[str] | None = None,
    outputs: list[str] | None = None,
    started_at: datetime | None = None,
) -> RunRecord:
    """Record an EFGPP-internal (pure Python) processing step."""
    with Registry.open(project) as reg:
        run_id = reg.next_id("RUN")
    record = RunRecord(
        run_id=run_id, step_id=step_id, tool="efgpp", tool_version=__version__,
        command=["efgpp", "internal", description], inputs=inputs or [], outputs=outputs or [],
        started_at=started_at or utcnow(), completed_at=utcnow(), exit_code=0, status="completed",
        environment={"efgpp_version": __version__, "env": "core", "platform": platform.platform()},
    )
    _persist(project, record)
    return record


def _persist(project: Project, record: RunRecord) -> None:
    write_yaml(project.logs_dir / f"{record.run_id}.yaml", record.to_yaml())
    with Registry.open(project) as reg:
        reg.upsert("tool_runs", {
            "run_id": record.run_id, "step_id": record.step_id, "tool": record.tool,
            "tool_version": record.tool_version, "command": shlex.join(record.command),
            "inputs": record.inputs, "outputs": record.outputs,
            "started_at": record.started_at, "completed_at": record.completed_at,
            "exit_code": record.exit_code, "status": record.status,
            "environment": record.environment, "resources": record.resources,
            "log_path": record.log_path,
        })
