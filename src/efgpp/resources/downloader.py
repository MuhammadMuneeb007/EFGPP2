"""Resumable, retried, checksummed downloads into .efgpp/downloads/."""

from __future__ import annotations

import gzip
import hashlib
import shutil
from pathlib import Path

import httpx
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

USER_AGENT = "efgpp-data/0.1 (+https://github.com/)"


class DownloadError(RuntimeError):
    pass


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, max=60),
       retry=retry_if_exception_type((httpx.TransportError, httpx.RemoteProtocolError)), reraise=True)
def download(url: str, dest: Path, *, expected_sha256: str | None = None, client: httpx.Client | None = None,
             timeout: float = 120.0, expected_md5: str | None = None) -> str:
    """Download `url` to `dest` (resuming a partial .part file). Returns the SHA-256.

    `expected_md5` (e.g. a Zenodo published checksum) is verified in the same pass; on any
    mismatch the partial file is deleted and nothing is moved into place."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    headers = {"User-Agent": USER_AGENT}
    offset = part.stat().st_size if part.exists() else 0
    if offset:
        headers["Range"] = f"bytes={offset}-"
    own = client is None
    client = client or httpx.Client(follow_redirects=True, timeout=timeout)
    try:
        with client.stream("GET", url, headers=headers) as r:
            if r.status_code == 416:  # already complete
                pass
            elif r.status_code not in (200, 206):
                raise DownloadError(f"GET {url} -> HTTP {r.status_code}")
            else:
                mode = "ab" if r.status_code == 206 else "wb"
                with open(part, mode) as fh:
                    for chunk in r.iter_bytes(1024 * 1024):
                        fh.write(chunk)
    finally:
        if own:
            client.close()
    digest, md5 = sha256_md5_file(part)
    if expected_sha256 and digest != expected_sha256:
        part.unlink(missing_ok=True)
        raise DownloadError(f"checksum mismatch for {url}: {digest} != {expected_sha256}")
    if expected_md5 and md5 != expected_md5.lower():
        part.unlink(missing_ok=True)
        raise DownloadError(f"MD5 mismatch for {url}: {md5} != published {expected_md5}")
    part.replace(dest)
    return digest


def sha256_md5_file(path: Path) -> tuple[str, str]:
    """SHA-256 and MD5 of a file in one read."""
    h, m = hashlib.sha256(), hashlib.md5()  # noqa: S324 - MD5 only to verify published checksums
    with open(path, "rb") as fh:
        while chunk := fh.read(8 * 1024 * 1024):
            h.update(chunk)
            m.update(chunk)
    return h.hexdigest(), m.hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(8 * 1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def gunzip(src: Path, dest: Path) -> Path:
    with gzip.open(src, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout, length=8 * 1024 * 1024)
    return dest
