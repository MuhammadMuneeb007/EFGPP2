from __future__ import annotations

import hashlib
from pathlib import Path

from hypothesis import given, settings
from hypothesis import strategies as st

from efgpp.constants import Modality, Origin, StorageMode
from efgpp.data.storage import place_source, resolve_mode
from efgpp.project import Project
from efgpp.resources.checksums import Checksum, checksum_paths


def test_single_file_is_plain_sha256(tmp_path: Path) -> None:
    f = tmp_path / "a.txt"
    f.write_bytes(b"hello")
    c = checksum_paths([f])
    assert c.algorithm == "sha256" and c.value == hashlib.sha256(b"hello").hexdigest()
    assert Checksum.parse(str(c)) == c


def test_fileset_manifest_is_order_independent(tmp_path: Path) -> None:
    a, b = tmp_path / "x.bed", tmp_path / "x.bim"
    a.write_bytes(b"1")
    b.write_bytes(b"2")
    assert checksum_paths([a, b]) == checksum_paths([b, a])
    assert checksum_paths([a, b]).algorithm == "sha256-manifest"


def test_sampled_checksum_is_labelled(tmp_path: Path) -> None:
    f = tmp_path / "big.bin"
    f.write_bytes(b"0" * 5000)
    c = checksum_paths([f], full_limit_bytes=100)
    assert c.algorithm == "sha256-sampled"


@settings(max_examples=25, deadline=None)
@given(st.binary(min_size=0, max_size=2048))
def test_checksum_detects_any_change(tmp_path_factory, data: bytes) -> None:  # type: ignore[no-untyped-def]
    d = tmp_path_factory.mktemp("h")
    f = d / "f"
    f.write_bytes(data)
    before = checksum_paths([f])
    f.write_bytes(data + b"!")
    assert checksum_paths([f]) != before


def test_auto_mode_policy(project: Project) -> None:
    gb = 1024**3
    assert resolve_mode(StorageMode.AUTO, fmt="bed", size=10, project=project) == StorageMode.REFERENCE
    assert resolve_mode(StorageMode.AUTO, fmt="csv", size=10, project=project) == StorageMode.COPY
    assert resolve_mode(StorageMode.AUTO, fmt="csv", size=2 * gb, project=project) == StorageMode.REFERENCE
    project.config.storage.copy_threshold_gb = 5
    assert resolve_mode(StorageMode.AUTO, fmt="csv", size=2 * gb, project=project) == StorageMode.COPY
    assert resolve_mode(StorageMode.COPY, fmt="bed", size=10, project=project) == StorageMode.COPY


def test_copy_reference_and_link(project: Project, tmp_path: Path) -> None:
    src = tmp_path / "ext" / "pheno.csv"
    src.parent.mkdir()
    src.write_text("IID,x\n1,2\n")
    copied = place_source(project, modality=Modality.PHENOTYPE, origin=Origin.OBSERVED, source_id="PH001",
                          record_path=src, members=[src], requested=StorageMode.COPY, fmt="csv")
    assert copied.mode == StorageMode.COPY and copied.record_path.is_relative_to(project.root)
    assert copied.record_path.read_text() == src.read_text()
    ref = place_source(project, modality=Modality.PHENOTYPE, origin=Origin.OBSERVED, source_id="PH002",
                       record_path=src, members=[src], requested=StorageMode.REFERENCE, fmt="csv")
    assert ref.mode == StorageMode.REFERENCE and ref.record_path == src.resolve()
    linked = place_source(project, modality=Modality.PHENOTYPE, origin=Origin.OBSERVED, source_id="PH003",
                          record_path=src, members=[src], requested=StorageMode.LINK, fmt="csv")
    # Symlinks may be unavailable (e.g. Windows without developer mode): then reference, never a copy.
    assert linked.mode in (StorageMode.LINK, StorageMode.REFERENCE)
    assert copied.checksum == ref.checksum == linked.checksum
