"""Zenodo records: resolve a file through the current metadata API, then download and verify it.

No download link is hard-coded: the record is queried (redirects to newer record ids are
followed), the exact file is located, its published checksum is read and verified after
download (together with a SHA-256 of our own), and only then moved into place.
"""

from __future__ import annotations

import shutil
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import httpx

from efgpp.resources.downloader import DownloadError, download, sha256_md5_file

API = "https://zenodo.org/api/records/{record}"
GB = 1024**3


@dataclass
class ZenodoFile:
    record: str
    resolved_record: str
    key: str
    size: int
    checksum: str  # "md5:<hex>"
    url: str
    title: str | None = None
    doi: str | None = None
    license: str | None = None
    version: str | None = None

    @property
    def md5(self) -> str | None:
        algo, _, value = self.checksum.partition(":")
        return value.lower() if algo == "md5" and value else None

    def manifest(self) -> dict[str, Any]:
        return {"provider_record": f"zenodo:{self.record}", "resolved_record": self.resolved_record,
                "source_filename": self.key, "source_checksum": self.checksum, "source_size": self.size,
                "source_url": self.url, "doi": self.doi, "record_title": self.title, "record_license": self.license,
                "record_version": self.version}


def record_files(record: str, client: httpx.Client) -> list[ZenodoFile]:
    r = client.get(API.format(record=record), follow_redirects=True)
    if r.status_code != 200:
        raise DownloadError(f"Zenodo record {record}: HTTP {r.status_code}")
    d = r.json()
    meta = d.get("metadata") or {}
    lic = meta.get("license")
    out = []
    for f in d.get("files") or []:
        out.append(ZenodoFile(
            record=str(record), resolved_record=str(d.get("id")), key=f["key"], size=int(f.get("size") or 0),
            checksum=str(f.get("checksum") or ""), url=f["links"]["self"], title=meta.get("title"),
            doi=d.get("doi") or meta.get("doi"), license=lic.get("id") if isinstance(lic, dict) else lic,
            version=meta.get("version")))
    return out


def record_file(record: str, filename: str, client: httpx.Client) -> ZenodoFile:
    files = record_files(record, client)
    for f in files:
        if f.key == filename:
            return f
    raise DownloadError(f"{filename} is not in Zenodo record {record} (files: {', '.join(f.key for f in files)})")


def check_space(directory: Path, need: int, what: str, say: Callable[[str], None] | None = None) -> None:
    """Refuse to start a download that cannot fit; report large downloads before they start."""
    directory.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(directory).free
    if say and need >= GB:
        say(f"{what}: {need / GB:.1f} GB needed in {directory} ({free / GB:.1f} GB free)")
    if free < need:
        raise DownloadError(f"not enough disk space for {what}: need {need / GB:.1f} GB, "
                            f"{free / GB:.1f} GB free in {directory}")


def fetch_file(record: str, filename: str, dest_dir: Path, *, client: httpx.Client,
               expected_md5: str | None = None, extract_factor: float = 1.0,
               say: Callable[[str], None] | None = None) -> tuple[Path, str, ZenodoFile]:
    """Download one file of a record into dest_dir (reusing a verified copy). Returns (path, sha256, meta).

    `expected_md5` is the checksum EFGPP knows as published; if Zenodo now reports another
    one the download is refused (the upstream file changed and must be reviewed)."""
    zf = record_file(record, filename, client)
    if expected_md5 and zf.md5 and zf.md5 != expected_md5.lower():
        raise DownloadError(f"Zenodo {record}/{filename}: published MD5 is now {zf.md5}, expected "
                            f"{expected_md5}; the upstream file changed - review before using it")
    md5 = zf.md5 or expected_md5
    dest = dest_dir / filename
    if dest.exists():
        sha, got = sha256_md5_file(dest)
        if md5 is None or got == md5:
            return dest, sha, zf
        dest.unlink()
    check_space(dest_dir, int(zf.size * (1 + extract_factor)), f"{filename} ({zf.size / GB:.2f} GB download)", say)
    if say:
        say(f"downloading {filename} from Zenodo record {zf.resolved_record} ({zf.size / 1024**2:,.0f} MB)")
    sha = download(zf.url, dest, client=client, timeout=600, expected_md5=md5)
    return dest, sha, zf
