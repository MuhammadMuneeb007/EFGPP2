"""PredictDB prediction models for MetaXcan/PrediXcan."""

from __future__ import annotations

import tarfile
from pathlib import Path
from typing import Any

from efgpp.data.references.base import FetchedResource, ReferenceProvider
from efgpp.resources.downloader import download

# GTEx v8 MASHR expression models (Barbeira et al. 2021), distributed via Zenodo.
DEFAULT_URL = "https://zenodo.org/record/3518299/files/mashr_eqtl.tar?download=1"


class PredictDBProvider(ReferenceProvider):
    name = "predictdb"
    license = "CC BY 4.0 (see predictdb.org)"
    homepage = "https://predictdb.org/"

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "default_models": "GTEx v8 MASHR eQTL (GRCh38, chr_pos_ref_alt_b38 ids)",
                "url": self.project.resources.predictdb.url or DEFAULT_URL}

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        cfg = self.project.resources.predictdb
        url = cfg.url or DEFAULT_URL
        archive = staging / "models.tar"
        digest = download(url, archive, client=self.client, timeout=600)
        with tarfile.open(archive) as tar:
            tar.extractall(staging, filter="data")
        archive.unlink()
        dbs = sorted(staging.rglob("*.db"))
        if not dbs:
            raise RuntimeError("downloaded archive contains no .db models")
        models_dir = dbs[0].parent
        return FetchedResource(self.name, cfg.version or "gtex_v8_mashr", dbs, models_dir, url, self.license,
                               genome_build="GRCh38", metadata={"download_sha256": digest, "models": len(dbs)})

    def validate(self, path: Path) -> list[str]:
        problems = super().validate(path)
        if not problems and not list(Path(path).glob("*.db")):
            problems.append("no .db model files")
        return problems
