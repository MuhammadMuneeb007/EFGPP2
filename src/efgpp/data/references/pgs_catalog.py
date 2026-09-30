"""PGS Catalog: registers published score definitions as reference resources.

Scores are only downloaded here. Applying them (e.g. with pgsc_calc) and any
target-specific PRS work belong to the Representation layer.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from efgpp.data.references.base import FetchedResource, ReferenceProvider
from efgpp.resources.downloader import download

API = "https://www.pgscatalog.org/rest"


class PGSCatalogProvider(ReferenceProvider):
    name = "pgs_catalog"
    license = "PGS Catalog terms (per-score licences apply)"
    homepage = "https://www.pgscatalog.org/"

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "api": API, "score_ids": self.project.resources.pgs_catalog.score_ids}

    def query(self, pgs_id: str) -> dict[str, Any]:
        r = self.client.get(f"{API}/score/{pgs_id}")
        r.raise_for_status()
        return dict(r.json())

    def scoring_file_url(self, meta: dict[str, Any], build: str | None) -> str:
        if build:
            harmonized = meta.get("ftp_harmonized_scoring_files") or {}
            entry = harmonized.get(build) or {}
            if entry.get("positions"):
                return str(entry["positions"])
        return str(meta["ftp_scoring_file"])

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        cfg = self.project.resources.pgs_catalog
        build = build or cfg.genome_build or self.project.resources.genome.build
        if not cfg.score_ids:
            raise RuntimeError("list score ids under resources.yaml: pgs_catalog.score_ids")
        files, meta_all = [], {}
        for pgs_id in cfg.score_ids:
            meta = self.query(pgs_id)
            url = self.scoring_file_url(meta, build)
            dest = staging / Path(url).name
            download(url, dest, client=self.client)
            files.append(dest)
            meta_all[pgs_id] = {"trait": meta.get("trait_reported"), "variants": meta.get("variants_number"),
                                "url": url, "license": meta.get("license")}
        version = "+".join(cfg.score_ids)
        return FetchedResource(self.name, version, files, staging, API, self.license, genome_build=build,
                               metadata={"scores": meta_all})
