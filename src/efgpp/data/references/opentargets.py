"""Open Targets Platform: GraphQL for targeted queries, bulk downloads for systematic use.

Prefer the downloadable datasets for bulk analyses instead of thousands of one-entity
GraphQL calls.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from efgpp.data.references.base import FetchedResource, ReferenceProvider
from efgpp.resources.downloader import download

GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
FTP = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"

VERSION_QUERY = "query { meta { apiVersion { x y z } dataVersion { year month } } }"


class OpenTargetsProvider(ReferenceProvider):
    name = "opentargets"
    license = "CC0 1.0"
    homepage = "https://platform.opentargets.org/"

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "graphql": GRAPHQL, "downloads": FTP,
                "entities": ["variant", "target", "disease", "credible set", "evidence"]}

    def query(self, query: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
        r = self.client.post(GRAPHQL, json={"query": query, "variables": variables or {}})
        r.raise_for_status()
        payload = r.json()
        if payload.get("errors"):
            raise RuntimeError(f"Open Targets GraphQL error: {payload['errors']}")
        return dict(payload["data"])

    def version(self) -> str | None:
        dv = self.query(VERSION_QUERY)["meta"]["dataVersion"]
        return f"{dv['year']}.{int(dv['month']):02d}"

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        """Download selected dataset files: resources.yaml opentargets.datasets = [urls or paths
        relative to <FTP>/<version>/output/]."""
        cfg = self.project.resources.opentargets
        version = cfg.version or self.version() or "latest"
        datasets = (cfg.model_extra or {}).get("datasets", [])
        if not datasets:
            raise RuntimeError("list dataset files under resources.yaml: opentargets.datasets")
        files = []
        for d in datasets:
            url = d if d.startswith("http") else f"{FTP}/{version}/output/{d}"
            dest = staging / Path(url).name
            download(url, dest, client=self.client)
            files.append(dest)
        return FetchedResource(self.name, version, files, staging, FTP, self.license)
