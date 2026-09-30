"""Germline genotype sources (PGEN, BED, BGEN, VCF, BCF), referenced in place by default.

The adapter reads only sidecar files; matrix-level work (QC, PCA, kinship) is delegated
to PLINK 2 by the workflow steps in efgpp.data.genotype.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import polars as pl

from efgpp.config.data import GenotypeSource
from efgpp.constants import ArtifactStatus, GenomeBuild, Modality, Origin, TemporalType
from efgpp.data.adapters.base import DataAdapter
from efgpp.data.aliases import plink_native_ids
from efgpp.data.artifacts import Artifact
from efgpp.data.genotype.build import BuildInference, infer_build
from efgpp.data.genotype.formats import (
    GenotypeFileset,
    count_variants,
    detect_format,
    read_samples,
    resolve_fileset,
    vcf_header,
)
from efgpp.data.genotype.qc import variant_table
from efgpp.data.io import write_parquet
from efgpp.data.schemas.artifact import VARIANT_SCHEMA
from efgpp.data.validation import ValidationReport, plural
from efgpp.resources.checksums import checksum_paths

BUILD_SAMPLE_VARIANTS = 200_000


class GenotypeAdapter(DataAdapter):
    modality = Modality.GENOTYPE
    temporal_type = TemporalType.STATIC
    source: GenotypeSource

    @classmethod
    def detect(cls, path: Path) -> bool:
        return detect_format(path) is not None

    def fileset(self, registered: bool = True) -> GenotypeFileset:
        art = self.source_artifact() if registered else None
        if art is not None:
            return resolve_fileset(Path(art.path), art.format)
        return resolve_fileset(self.source_path, self.source.format)

    def members(self) -> tuple[Path, list[Path], str]:
        fs = resolve_fileset(self.source_path, self.source.format)
        return fs.prefix, fs.members + fs.missing, fs.format

    def native_ids(self, samples: pl.DataFrame) -> pl.Series:
        sid = self.source.sample_id
        return plink_native_ids(samples, sid.mode, sid.column, sid.separator)

    def _fastas(self) -> dict[GenomeBuild, Path]:
        """Reference FASTAs installed by the resource manager, keyed by build."""
        out: dict[GenomeBuild, Path] = {}
        g = self.project.resources.genome
        if g.fasta:
            p = self.project.resolve(g.fasta)
            if p.exists():
                out[GenomeBuild.normalize(g.build)] = p
        rows = self.reg.rows("SELECT genome_build, local_path FROM resources WHERE name = 'genome'")
        for r in rows:
            b = GenomeBuild.normalize(r["genome_build"])
            p = Path(r["local_path"])
            if b not in out and p.exists() and p.suffix in (".fa", ".fasta"):
                out[b] = p
        return out

    def infer_build(self, fs: GenotypeFileset | None = None) -> BuildInference:
        fs = fs or self.fileset()
        declared = self.source.genome_build
        if declared in ("auto", None):
            declared = self.project.config.defaults.genome_build
        variants = None
        meta = None
        try:
            if fs.format == "vcf":
                meta = vcf_header(fs.main)[0]
            if fs.format in ("pgen", "bed", "vcf"):
                variants = variant_table(fs, None, limit=BUILD_SAMPLE_VARIANTS)
        except (NotImplementedError, ValueError, OSError):
            variants = None
        panel = self.project.resource_root / "genomes" / "build_marker_panel.tsv"
        return infer_build(declared=declared, variants=variants, vcf_meta=meta, fastas=self._fastas(),
                           marker_panel=panel if panel.exists() else None)

    # ------------------------------------------------------------------ lifecycle
    def inspect(self) -> dict[str, Any]:
        info: dict[str, Any] = {"source_id": self.source.id, "modality": "genotype",
                                "origin": self.source.origin.value, "path": str(self.source_path)}
        try:
            fs = self.fileset(registered=False)
        except FileNotFoundError as exc:
            info.update(exists=False, error=str(exc))
            return info
        info.update(exists=fs.complete, format=fs.format, files=[m.name for m in fs.members],
                    size=sum(m.stat().st_size for m in fs.members))
        try:
            info["samples"] = read_samples(fs).height
        except NotImplementedError as exc:
            info["samples"] = None
            info["note"] = str(exc)
        info["variants"] = count_variants(fs)
        b = self.infer_build(fs)
        info["genome_build"] = b.build.value
        info["build_confidence"] = round(b.confidence, 2)
        return info

    def validate(self) -> ValidationReport:
        rep = ValidationReport(self.subject(), self.source.id)
        try:
            fs = self.fileset()
        except FileNotFoundError as exc:
            rep.fail(str(exc))
            return rep
        if fs.missing:
            rep.fail(f"incomplete {fs.format} fileset; missing {', '.join(m.name for m in fs.missing)}")
            return rep
        rep.ok(f"{fs.format.upper()} fileset complete ({', '.join(m.name for m in fs.members)})")
        try:
            samples = read_samples(fs)
        except NotImplementedError as exc:
            rep.warn(str(exc))
            samples = None
        if samples is not None:
            ids = self.native_ids(samples)
            rep.ok(f"{plural(samples.height, 'sample')}")
            dup = ids.filter(ids.is_duplicated()).unique().to_list()
            if dup:
                rep.fail(f"{plural(len(dup), 'duplicated sample ID')}", len(dup), dup)
            else:
                rep.ok("sample IDs are unique")
            if ids.null_count() or (ids == "").sum():
                rep.fail("samples without an ID")
        try:
            v = variant_table(fs, None, limit=BUILD_SAMPLE_VARIANTS)
            n = count_variants(fs)
            rep.ok(f"{plural(n, 'variant')}" if n is not None else f"variant file readable ({v.height:,} checked)")
            bad_pos = int((v.get_column("position").fill_null(0) < 1).sum())
            if bad_pos:
                rep.fail(f"{plural(bad_pos, 'variant')} with invalid positions")
            missing_ids = int(v.get_column("variant_id").is_in([".", ""]).sum())
            dup_ids = int(v.filter(~v.get_column("variant_id").is_in([".", ""])).get_column("variant_id").is_duplicated().sum())
            if missing_ids:
                rep.warn(f"{plural(missing_ids, 'variant')} without an ID (QC assigns chr:pos:ref:alt IDs)")
            if dup_ids:
                rep.warn(f"{plural(dup_ids, 'variant')} with duplicated IDs (QC removes duplicates)")
        except NotImplementedError as exc:
            rep.info(str(exc))
        build = self.infer_build(fs)
        if build.confident:
            methods = ", ".join(sorted({e.method for e in build.evidence if e.supports == build.build}))
            rep.ok(f"genome build {build.build.value} (confidence {build.confidence:.0%}; {methods})")
        else:
            details = "; ".join(f"{e.method}: {e.detail}" for e in build.evidence) or "no evidence"
            rep.fail(
                "genome build could not be determined confidently - set `genome_build` for "
                f"{self.source.id} in data.yaml ({details})"
            )
        return rep

    def standardize(self) -> Artifact | None:
        """Resolve participant IDs, record the variant table and the genome build.

        The genotype matrix itself is not converted here (large cohorts stay in place).
        """
        s = self.source
        src = self.source_artifact()
        if src is None:
            raise RuntimeError(f"{s.id} is not registered")
        fs = self.fileset()
        samples = read_samples(fs)
        native = self.native_ids(samples)
        res = self.resolve_ids(native)
        pids = native.replace_strict(dict(res.mapping.iter_rows()), default=None)
        build = self.infer_build(fs)
        build_value = build.build.value if build.confident else None
        src.genome_build = build_value
        src.participant_count = samples.height
        src.feature_count = count_variants(fs)
        src.status = ArtifactStatus.STANDARDIZED
        src.metadata.update({"build_inference": build.to_dict(), "sample_id": s.sample_id.model_dump()})
        self.store.update(src)
        self.record_assays(src, pids, sample_keys=native)

        try:
            variants = variant_table(fs, build_value)
        except NotImplementedError:
            return src
        variants = VARIANT_SCHEMA.validate(variants)
        out_dir = self.project.artifact_dir(Origin.DERIVED, Modality.VARIANTS, s.id)
        out = write_parquet(variants, out_dir / f"{s.id}_variants.parquet", self.project.config.storage.compression)
        self.register_replacing(Artifact(
            artifact_name=f"{s.id}_variants", artifact_type="variant_table", modality=Modality.VARIANTS,
            origin=Origin.DERIVED, status=ArtifactStatus.READY, path=str(out), format="parquet",
            size=out.stat().st_size, checksum=str(checksum_paths([out])), feature_count=variants.height,
            genome_build=build_value, source_id=s.id, parent_artifact_ids=[src.artifact_id],  # type: ignore[list-item]
            tool="efgpp", temporal_type=TemporalType.STATIC,
        ))
        return src

    def summarize(self) -> dict[str, Any]:
        out = super().summarize()
        src = self.source_artifact()
        if src:
            out.update({"format": src.format, "storage_mode": src.storage_mode, "genome_build": src.genome_build,
                        "participants": src.participant_count, "variants": src.feature_count, "path": src.path})
        return out
