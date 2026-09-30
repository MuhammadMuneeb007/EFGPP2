"""DataManager: orchestrates the lifecycle of every configured source.

    REGISTER -> INSPECT -> VALIDATE -> PARTICIPANT ID RESOLUTION -> TIMELINE/BIOSPECIMEN
    RESOLUTION -> FORMAT STANDARDIZATION -> QC -> DERIVATION -> ANNOTATION -> PREDICTION
    -> AVAILABILITY -> REPORT -> FREEZE
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from efgpp.constants import PREDICTIVE_MODALITIES, ArtifactStatus, Modality
from efgpp.data.adapters import AdapterContext, DataAdapter, adapter_for
from efgpp.data.adapters.base import QCResult
from efgpp.data.artifacts import Artifact
from efgpp.data.biospecimens import load_biospecimens_file
from efgpp.data.io import as_text, read_table
from efgpp.data.participants import load_master_list
from efgpp.data.registry import Registry
from efgpp.data.timeline import load_events_file
from efgpp.data.validation import ValidationReport, plural
from efgpp.project import Project


@dataclass
class Registration:
    source_id: str
    modality: str
    artifact_id: str
    storage_mode: str
    notes: list[str]


class DataManager:
    def __init__(self, project: Project) -> None:
        self.project = project

    def _adapters(self, reg: Registry, source_ids: list[str] | None = None) -> list[DataAdapter]:
        ctx = AdapterContext(self.project, reg)
        return [
            adapter_for(ctx, modality, source)
            for modality, source in self.project.data.iter_sources()
            if source_ids is None or source.id in source_ids
        ]

    def adapter(self, reg: Registry, source_id: str) -> DataAdapter:
        modality, source = self.project.data.get_source(source_id)
        return adapter_for(AdapterContext(self.project, reg), modality, source)

    # ----------------------------------------------------------------- register
    def register_all(self) -> list[Registration]:
        out = []
        with Registry.open(self.project) as reg:
            load_master_list(self.project, reg)
            for ad in self._adapters(reg):
                art = ad.register()
                out.append(Registration(ad.source.id, ad.modality.value, art.artifact_id or "",
                                        art.storage_mode or "", art.metadata.get("notes", [])))
            load_events_file(self.project, reg)
            load_biospecimens_file(self.project, reg)
            reg.export_parquet()
        return out

    # ------------------------------------------------------------------ inspect
    def inspect(self) -> list[dict[str, Any]]:
        with Registry.open(self.project) as reg:
            return [ad.inspect() for ad in self._adapters(reg)]

    # ----------------------------------------------------------------- validate
    def validate(self, source_ids: list[str] | None = None) -> list[ValidationReport]:
        reports = []
        with Registry.open(self.project) as reg:
            for ad in self._adapters(reg, source_ids):
                rep = ad.validate()
                reports.append(rep)
                rep.save(self.project.qc_dir(ad.modality_qc_dir()) / f"{ad.source.id}_validation.json")
                art = ad.source_artifact()
                if art is not None and art.status in (ArtifactStatus.REGISTERED, ArtifactStatus.VALIDATED):
                    ad.store.set_status(art.artifact_id,  # type: ignore[arg-type]
                                        ArtifactStatus.VALIDATED if rep.passed else ArtifactStatus.REGISTERED,
                                        validation_passed=rep.passed, validation_failures=rep.n_fail)
                elif art is not None:
                    ad.store.set_status(art.artifact_id, art.status, validation_passed=rep.passed,  # type: ignore[arg-type]
                                        validation_failures=rep.n_fail)
            if source_ids is None:
                from efgpp.data.gwas import validate_gwas

                for g in self.project.data.gwas:
                    rep = validate_gwas(self.project, g)
                    rep.save(self.project.qc_dir("gwas") / f"{g.id}_validation.json")
                    reports.append(rep)
                reports.append(self.cross_source_report(reg))
        return reports

    def cross_source_report(self, reg: Registry) -> ValidationReport:
        rep = ValidationReport("PROJECT (cross-source checks)")
        data = self.project.data
        phenos = data.observed.phenotypes
        predictive = [s for m, s in data.iter_sources() if m in PREDICTIVE_MODALITIES]
        if not phenos:
            rep.warn("no phenotype configured (prediction will need at least one)")
        else:
            rep.ok(f"{plural(len(phenos), 'phenotype')} configured")
        if not predictive:
            rep.warn("no participant-level predictive modality configured")
        else:
            rep.ok(f"{plural(len(predictive), 'predictive source')} configured")
        if not data.observed.genotype:
            rep.info("no genotype configured (genotype-centric analyses will be unavailable)")

        ctx = AdapterContext(self.project, reg)
        id_sets: dict[str, set[str]] = {}
        for modality, source in data.iter_sources():
            try:
                ad = adapter_for(ctx, modality, source)
                if modality == Modality.GENOTYPE:
                    from efgpp.data.genotype.formats import read_samples

                    natives = ad.native_ids(read_samples(ad.fileset()))  # type: ignore[attr-defined]
                elif hasattr(ad, "read_matrix"):
                    natives = ad.read_matrix().obs.get_column("native_id")  # type: ignore[attr-defined]
                else:
                    natives = as_text(read_table(ad.source_path, source.format), source.participant_id_column)  # type: ignore[attr-defined]
                res = ctx.resolver.resolve(source.id, natives)
                id_sets[source.id] = set(res.mapping.get_column("participant_id").to_list())
            except Exception:  # noqa: BLE001 - per-source problems are reported by that source
                continue
        geno_ids = set().union(*(id_sets.get(g.id, set()) for g in data.observed.genotype)) if data.observed.genotype else None
        for p in phenos:
            ids = id_sets.get(p.id)
            if ids is None or geno_ids is None:
                continue
            missing = len(ids - geno_ids)
            if missing:
                rep.warn(f"{p.id} ({p.name}): {plural(missing, 'participant')} have no genotype", missing)
            else:
                rep.ok(f"{p.id} ({p.name}): every participant has genotype data")
        builds = {g.id: g.genome_build for g in data.observed.genotype if g.genome_build not in ("auto", None)}
        if len(set(builds.values())) > 1:
            rep.warn(f"genotype sources declare different genome builds: {builds}")
        return rep

    # ------------------------------------------------------------- per source
    def standardize(self, source_id: str) -> Artifact | None:
        with Registry.open(self.project) as reg:
            ad = self.adapter(reg, source_id)
            art = ad.source_artifact()
            if art is None:
                art = ad.register()
            if art.metadata.get("validation_passed") is False:
                raise RuntimeError(f"{source_id} failed validation; fix it and run `efgpp data validate`")
            if art.metadata.get("validation_passed") is None:
                rep = ad.validate()
                if not rep.passed:
                    raise RuntimeError(f"{source_id} failed validation:\n" + "\n".join(
                        c.message for c in rep.checks if c.status == "fail"))
            result = ad.standardize()
            reg.export_parquet()
            return result

    def qc(self, source_id: str) -> QCResult | None:
        with Registry.open(self.project) as reg:
            return self.adapter(reg, source_id).qc()

    def summaries(self) -> list[dict[str, Any]]:
        with Registry.open(self.project) as reg:
            return [ad.report() for ad in self._adapters(reg)]
