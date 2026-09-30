"""Checksums for artifact identity (SHA-256) and fast local scanning (optional BLAKE3).

Filesets (PLINK prefixes, VCF + index) and directories (Zarr stores) are hashed as a
manifest of their member hashes, so the result is independent of traversal order.
Files above the configured size limit get a *sampled* SHA-256 over size + head + middle
+ tail blocks; its algorithm label ("sha256-sampled") makes that explicit.
"""

from __future__ import annotations

import hashlib
import os
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

CHUNK = 8 * 1024 * 1024
SAMPLE_BLOCK = 64 * 1024 * 1024


@dataclass(frozen=True)
class Checksum:
    algorithm: str
    value: str

    def __str__(self) -> str:
        return f"{self.algorithm}:{self.value}"

    @classmethod
    def parse(cls, text: str) -> Checksum:
        algorithm, _, value = text.partition(":")
        return cls(algorithm, value)


def _new_hasher(algorithm: str):  # type: ignore[no-untyped-def]
    if algorithm == "blake3":
        try:
            import blake3
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise RuntimeError("BLAKE3 requested but the 'blake3' package is not installed") from exc
        return blake3.blake3()
    return hashlib.sha256()


def hash_file(path: Path, algorithm: str = "sha256") -> str:
    h = _new_hasher(algorithm)
    with open(path, "rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
    return str(h.hexdigest())


def sampled_sha256(path: Path, block: int = SAMPLE_BLOCK) -> str:
    size = path.stat().st_size
    h = hashlib.sha256()
    h.update(f"size={size}\n".encode())
    with open(path, "rb") as fh:
        for offset in sorted({0, max(0, size // 2 - block // 2), max(0, size - block)}):
            fh.seek(offset)
            h.update(fh.read(block))
    return h.hexdigest()


def iter_files(path: Path) -> Iterable[Path]:
    if path.is_file():
        yield path
        return
    for root, dirs, files in os.walk(path):
        dirs.sort()
        for name in sorted(files):
            yield Path(root) / name


def checksum_paths(
    paths: Iterable[Path],
    *,
    full_limit_bytes: int | None = None,
    algorithm: str = "sha256",
) -> Checksum:
    """Checksum one file, one directory, or a fileset (several files) as a single identity."""
    members: list[tuple[str, str]] = []
    sampled = False
    path_list = [Path(p) for p in paths]
    for base in path_list:
        for f in iter_files(base):
            rel = f.name if base.is_file() else f"{base.name}/{f.relative_to(base).as_posix()}"
            if (
                algorithm == "sha256"
                and full_limit_bytes is not None
                and f.stat().st_size > full_limit_bytes
            ):
                members.append((rel, sampled_sha256(f)))
                sampled = True
            else:
                members.append((rel, hash_file(f, algorithm)))
    label = f"{algorithm}-sampled" if sampled else algorithm
    if len(path_list) == 1 and path_list[0].is_file():
        return Checksum(label, members[0][1])
    manifest = "".join(f"{name}\t{digest}\n" for name, digest in sorted(members))
    return Checksum(f"{label}-manifest", hashlib.sha256(manifest.encode()).hexdigest())


def total_size(paths: Iterable[Path]) -> int:
    return sum(f.stat().st_size for p in paths for f in iter_files(Path(p)))


def text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()
