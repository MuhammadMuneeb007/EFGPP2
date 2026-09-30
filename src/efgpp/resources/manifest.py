"""Per-resource manifest (resources/<name>/<version>/manifest.yaml)."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class ResourceManifest:
    resource_id: str
    name: str
    version: str
    release_date: str | None
    genome_build: str | None
    source: str
    license: str
    checksum: str
    download_date: str
    local_path: str
    files: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def write(self, directory: Path) -> Path:
        path = directory / "manifest.yaml"
        path.write_text(yaml.safe_dump(asdict(self), sort_keys=False), encoding="utf-8")
        return path

    @classmethod
    def read(cls, directory: Path) -> ResourceManifest:
        return cls(**yaml.safe_load((directory / "manifest.yaml").read_text(encoding="utf-8")))
