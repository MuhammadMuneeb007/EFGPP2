"""NHGRI-EBI GWAS Catalog (REST API v2) for metadata / association discovery.

REST API v1 is retired. The API does not serve full summary statistics; those come from
the official summary-statistics FTP, organised in blocks of 1,000 accessions.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from efgpp.data.references.base import FetchedResource, ReferenceProvider
from efgpp.resources.downloader import download

API = "https://www.ebi.ac.uk/gwas/rest/api/v2"
SUMSTATS_FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"


def sumstats_directory(accession: str) -> str:
    """FTP directory of a GCST accession, e.g. GCST90000123 ->
    .../GCST90000001-GCST90001000/GCST90000123."""
    m = re.fullmatch(r"GCST(\d+)", accession)
    if not m:
        raise ValueError(f"not a GWAS Catalog study accession: {accession!r}")
    digits = m.group(1)
    n = int(digits)
    lo = ((n - 1) // 1000) * 1000 + 1
    width = len(digits)
    return f"{SUMSTATS_FTP}/GCST{lo:0{width}d}-GCST{lo + 999:0{width}d}/{accession}"


class GWASCatalogProvider(ReferenceProvider):
    name = "gwas_catalog"
    license = "EMBL-EBI terms of use"
    homepage = "https://www.ebi.ac.uk/gwas/docs/api"
    bulk = False

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "api": API, "summary_statistics": SUMSTATS_FTP,
                "note": "REST v2 for metadata; summary statistics via FTP downloads"}

    def query(self, endpoint: str, **params: Any) -> dict[str, Any]:
        """GET <API>/<endpoint> with query parameters, e.g. query('studies', efo_trait='...')."""
        r = self.client.get(f"{API}/{endpoint.lstrip('/')}", params=params)
        r.raise_for_status()
        return dict(r.json())

    def list_sumstats(self, accession: str) -> list[str]:
        r = self.client.get(sumstats_directory(accession) + "/")
        r.raise_for_status()
        return sorted(set(re.findall(r'href="([^"?/][^"]*)"', r.text)))

    def fetch_sumstats(self, accession: str, staging: Path, filename: str | None = None) -> FetchedResource:
        files = self.list_sumstats(accession)
        candidates = [f for f in files if filename in (None, f) and f.endswith((".tsv.gz", ".tsv", ".txt.gz"))]
        if not candidates:
            raise RuntimeError(f"no summary-statistics file found for {accession}")
        url = f"{sumstats_directory(accession)}/{candidates[0]}"
        dest = staging / candidates[0]
        digest = download(url, dest, client=self.client)
        return FetchedResource(self.name, accession, [dest], dest, url, self.license,
                               metadata={"download_sha256": digest, "kind": "summary_statistics"})
