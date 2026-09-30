"""ReferenceProvider: the interface for knowledge that belongs to no participant."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar

import httpx

from efgpp.project import Project


@dataclass
class FetchedResource:
    name: str
    version: str
    files: list[Path]
    local_path: Path  # main file or directory
    source: str
    license: str
    genome_build: str | None = None
    release_date: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


class ReferenceProvider(ABC):
    """Every provider implements fetch / version / validate / cache / query / describe."""

    name: ClassVar[str]
    license: ClassVar[str] = "see source"
    homepage: ClassVar[str] = ""
    # Whether `fetch` downloads bulk data (True) or the provider is query-only (False).
    bulk: ClassVar[bool] = True

    def __init__(self, project: Project, client: httpx.Client | None = None) -> None:
        self.project = project
        self._client = client

    @property
    def client(self) -> httpx.Client:
        if self._client is None:
            self._client = httpx.Client(follow_redirects=True, timeout=60.0,
                                        headers={"User-Agent": "efgpp-data/0.1"})
        return self._client

    @property
    def config(self) -> Any:
        return getattr(self.project.resources, self.name, None)

    def cache(self, version: str | None = None) -> Path:
        """Directory where (a version of) this resource is stored."""
        d = self.project.resource_root / self.name
        if version:
            d = d / version
        d.mkdir(parents=True, exist_ok=True)
        return d

    @abstractmethod
    def describe(self) -> dict[str, Any]: ...

    def version(self) -> str | None:
        """Latest available upstream version (may contact the network)."""
        return None

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        raise NotImplementedError(f"{self.name} is query-only; nothing to download")

    def validate(self, path: Path) -> list[str]:
        """Problems with an installed copy (empty list = valid)."""
        return [] if path.exists() else [f"{path} does not exist"]

    def query(self, *args: Any, **kwargs: Any) -> Any:
        raise NotImplementedError(f"{self.name} does not support queries")
