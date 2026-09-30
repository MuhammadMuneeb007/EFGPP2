"""GTEx Portal API v2 (tissue expression and eQTL reference knowledge)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from efgpp.data.references.base import FetchedResource, ReferenceProvider
from efgpp.resources.downloader import download

API = "https://gtexportal.org/api/v2"


class GTExProvider(ReferenceProvider):
    name = "gtex"
    license = "GTEx Portal terms (open-access data)"
    homepage = "https://gtexportal.org/"

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "api": API, "provides": "tissue expression, eQTL/sQTL reference data"}

    def query(self, endpoint: str, **params: Any) -> dict[str, Any]:
        r = self.client.get(f"{API}/{endpoint.lstrip('/')}", params=params)
        r.raise_for_status()
        return dict(r.json())

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        cfg = self.project.resources.gtex
        urls = (cfg.model_extra or {}).get("files", [])
        if not urls:
            raise RuntimeError("list GTEx download URLs under resources.yaml: gtex.files")
        files = []
        for url in urls:
            dest = staging / Path(url).name
            download(url, dest, client=self.client)
            files.append(dest)
        return FetchedResource(self.name, cfg.version or "v8", files, staging, API, self.license, genome_build="GRCh38")
