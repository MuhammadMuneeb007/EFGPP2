from __future__ import annotations

import polars as pl

from efgpp.constants import ArtifactStatus, Modality, Origin
from efgpp.data.aliases import AliasResolver, plink_native_ids
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.participants import register_mapping
from efgpp.data.registry import Registry
from efgpp.project import Project


def _art(name: str, **kw) -> Artifact:  # type: ignore[no-untyped-def]
    return Artifact(artifact_name=name, artifact_type=kw.pop("artifact_type", "source"), modality=Modality.PHENOTYPE,
                    origin=Origin.OBSERVED, path=f"/x/{name}", format="csv", source_id=kw.pop("source_id", "PH001"), **kw)


def test_register_supersede_and_lineage(project: Project) -> None:
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        a = store.register(_art("a"))
        b = store.register_replacing(_art("a", checksum="new"))
        c = store.register(_art("c", artifact_type="standardized", parent_artifact_ids=[b.artifact_id]))
        assert a.artifact_id == "ART000001" and b.artifact_id == "ART000002"
        assert store.get(a.artifact_id).status == ArtifactStatus.SUPERSEDED  # type: ignore[arg-type]
        assert store.get(a.artifact_id).superseded_by == b.artifact_id  # type: ignore[arg-type]
        assert store.lineage(c.artifact_id) == [b.artifact_id]  # type: ignore[arg-type]
        # Superseded artifacts are kept, never deleted.
        assert len(store.find(include_superseded=True)) == 3
        assert [x.artifact_id for x in store.find()] == [b.artifact_id, c.artifact_id]


def test_registry_open_is_reentrant(project: Project) -> None:
    with Registry.open(project) as outer, Registry.open(project) as inner:
        assert outer is inner
        inner.set_meta("k", "v")
    with Registry.open(project) as reg:
        assert reg.get_meta("k") == "v"


def test_alias_mapping_and_participant_creation(project: Project, tmp_path) -> None:  # type: ignore[no-untyped-def]
    alias = tmp_path / "aliases.csv"
    pl.DataFrame({"participant_id": ["P1", "P2"], "lab_id": ["L-9", "L-8"]}).write_csv(alias)
    project.data.participants.aliases = [{"path": str(alias), "alias_column": "lab_id", "sources": ["RNA001"]}]  # type: ignore[list-item]
    project.save_data_config()
    project = Project.load(project.root)
    resolver = AliasResolver(project)
    res = resolver.resolve("RNA001", pl.Series(["L-9", "L-8", "P3"]))
    assert dict(res.mapping.iter_rows()) == {"L-9": "P1", "L-8": "P2", "P3": "P3"}
    assert res.n_mapped_by_alias == 2 and res.n_identity == 1
    # Sources without an applicable alias file use identity.
    assert dict(resolver.resolve("PH001", pl.Series(["L-9"])).mapping.iter_rows()) == {"L-9": "L-9"}
    with Registry.open(project) as reg:
        stats = register_mapping(reg, "RNA001", res.mapping)
        assert stats == {"participants": 3, "new": 3}
        assert register_mapping(reg, "RNA001", res.mapping)["new"] == 0
        assert reg.scalar("SELECT count(*) FROM sample_aliases") == 3


def test_plink_native_ids_modes() -> None:
    s = pl.DataFrame({"FID": ["F1", None], "IID": ["I1", "I2"]}, schema={"FID": pl.Utf8, "IID": pl.Utf8})
    assert plink_native_ids(s, "iid", "IID", "_").to_list() == ["I1", "I2"]
    assert plink_native_ids(s, "fid_iid", "IID", "_").to_list() == ["F1_I1", "0_I2"]
