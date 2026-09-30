"""reports/data/index.html - the data-quality report (23 sections).

Everything is read from the registry and the registered artifacts; the report never
scans directories to guess what exists.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import yaml
from jinja2 import Environment, FileSystemLoader, select_autoescape

from efgpp import __version__
from efgpp.constants import OMICS_MODALITIES, Modality
from efgpp.data import availability as availability_mod
from efgpp.data.artifacts import ArtifactStore
from efgpp.data.plan import build_plan
from efgpp.data.registry import Registry, utcnow
from efgpp.data.timeline import timeline_summary
from efgpp.project import Project
from efgpp.reporting import charts
from efgpp.setup.tools import available, resolve

TEMPLATES = Path(__file__).parent / "templates"


@dataclass
class Table:
    caption: str
    columns: list[str]
    rows: list[list[str]]


@dataclass
class Section:
    id: str
    title: str
    intro: str = ""
    charts: list[str] = field(default_factory=list)
    tables: list[Table] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    empty: str | None = None


def _table(caption: str, records: list[dict[str, Any]], columns: list[str] | None = None, limit: int = 500) -> Table:
    columns = columns or (list(records[0].keys()) if records else [])
    rows = [[charts.fmt(r.get(c)) for c in columns] for r in records[:limit]]
    if len(records) > limit:
        caption += f" (first {limit} of {len(records):,})"
    return Table(caption, columns, rows)


class ReportBuilder:
    def __init__(self, project: Project, reg: Registry) -> None:
        self.project = project
        self.reg = reg
        self.store = ArtifactStore(reg)
        self._n = 0

    def chart(self, fig: Any) -> str:
        self._n += 1
        return charts.to_html(fig, f"chart{self._n}")

    def qc_summary(self, artifact_id: str | None) -> dict[str, Any] | None:
        if not artifact_id:
            return None
        row = self.reg.one("SELECT summary, status FROM qc_runs WHERE artifact_id = ? ORDER BY created_at DESC LIMIT 1",
                           [artifact_id])
        return {"status": row["status"], **(row["summary"] or {})} if row else None

    # ------------------------------------------------------------------ sections
    def project_section(self) -> Section:
        p = self.project
        info = [
            {"field": "project", "value": p.config.project.name},
            {"field": "root", "value": str(p.root)},
            {"field": "species", "value": p.config.species.name},
            {"field": "default genome build", "value": p.config.defaults.genome_build},
            {"field": "EFGPP version", "value": __version__},
            {"field": "report generated", "value": utcnow().isoformat(timespec="seconds") + "Z"},
            *({"field": f"sha256 {k}", "value": (v or "")[:16]} for k, v in p.config_checksums().items()),
        ]
        sources = [{"source": s.id, "modality": m.value, "origin": s.origin.value, "path": s.path}
                   for m, s in p.data.iter_sources()]
        return Section("project", "1. Project", tables=[_table("Project", info), _table("Configured sources", sources)])

    def participants_section(self, av: availability_mod.Availability) -> Section:
        by_source = self.reg.rows("SELECT first_source AS source, count(*) AS participants FROM participants GROUP BY 1 ORDER BY 2 DESC")
        per_source = self.reg.rows("SELECT source_id, count(DISTINCT participant_id) AS participants FROM sample_aliases GROUP BY 1 ORDER BY 1")
        s = Section("participants", "2. Participants", intro=f"{av.n_participants:,} canonical participants.")
        if per_source:
            s.charts.append(self.chart(charts.bar([r["source_id"] for r in per_source], [r["participants"] for r in per_source],
                                                  "Participants resolved per source", horizontal=True)))
        s.tables += [_table("Participants per source", per_source), _table("First seen in", by_source)]
        return s

    def phenotypes_section(self) -> Section:
        s = Section("phenotypes", "3. Phenotypes")
        defs = self.reg.rows("SELECT * FROM phenotype_definitions ORDER BY phenotype_id")
        rows = []
        for d in defs:
            q = self.qc_summary(d["artifact_id"]) or {}
            obs = self.reg.frame("SELECT value_numeric, value_text FROM phenotype_observations WHERE phenotype_id = ? AND NOT is_missing",
                                 [d["phenotype_id"]])
            rows.append({"id": d["phenotype_id"], "name": d["name"], "type": d["type"], "observed": obs.height,
                         "missing": q.get("n_missing"), "summary": _pheno_summary(q), "qc": q.get("status")})
            if d["type"] == "continuous":
                s.charts.append(self.chart(charts.histogram(obs.get_column("value_numeric").to_list(),
                                                            f"{d['phenotype_id']} · {d['name']}", x=d.get("units") or "value")))
            else:
                counts = obs.get_column("value_text").value_counts(sort=True)
                s.charts.append(self.chart(charts.bar([str(x) for x in counts.get_column("value_text").to_list()],
                                                      counts.get_column("count").to_list(),
                                                      f"{d['phenotype_id']} · {d['name']}", y="participants")))
        s.tables.append(_table("Phenotype registry", rows))
        if not defs:
            s.empty = "No phenotype has been standardized yet."
        return s

    def timeline_section(self, tl: dict[str, Any]) -> Section:
        s = Section("timeline", "4. Timeline")
        if not tl["participants_per_event"]:
            s.empty = "No timed events registered (all measurements are static)."
            return s
        order = tl["event_order"]
        per = {r["event_id"]: r["participants"] for r in tl["participants_per_event"]}
        s.charts.append(self.chart(charts.bar(order, [per.get(e, 0) for e in order], "Participants per event", y="participants")))
        s.tables += [_table("Modalities per event", tl["modalities_per_event"]), _table("Attrition from first event", tl["attrition"]),
                     _table("Repeated measurements", tl["repeated_measurements"]), _table("Time gaps between events (days)", [tl["time_gaps"]])]
        return s

    def biospecimens_section(self) -> Section:
        s = Section("biospecimens", "5. Biospecimens")
        rows = self.reg.rows("SELECT coalesce(tissue, '?') AS tissue, coalesce(material, '?') AS material, coalesce(event_id, '(static)') AS event_id, "
                             "count(*) AS biospecimens, count(DISTINCT participant_id) AS participants FROM biospecimens GROUP BY 1, 2, 3 ORDER BY 1, 3")
        if not rows:
            s.empty = "No biospecimens registered."
        s.tables.append(_table("Biospecimens by tissue, material and event", rows))
        return s

    def availability_section(self, av: availability_mod.Availability) -> Section:
        s = Section("availability", "6. Modality availability",
                    intro="Participants per registered modality, grouped by origin. Predicted and derived data are never counted as observed.")
        cols = av.columns.to_dicts()
        if not cols:
            s.empty = "No participant-level artifacts yet."
            return s
        s.charts.append(self.chart(charts.bar([f"{c['column']} ({c['origin']})" for c in cols], [c["participants"] for c in cols],
                                              "Participants per modality", horizontal=True)))
        pats = av.patterns()
        names = av.column_names()
        if pats.height:
            patterns = [(tuple(n for n in names if row[n]), row["participants"]) for row in pats.to_dicts()]
            s.charts.append(self.chart(charts.upset(patterns, names, "Modality intersections (UpSet)")))
        if 1 < len(names) <= 30:
            z = [[av.count(a, b) if a != b else av.count(a) for b in names] for a in names]
            s.charts.append(self.chart(charts.heatmap(z, names, "Pairwise overlap (participants)")))
        s.tables += [_table("Availability columns", cols, ["column", "origin", "modality", "tissue", "participants", "artifact_id"]),
                     _table("Key intersections", av.key_intersections())]
        return s

    def covariates_section(self) -> Section:
        s = Section("covariates", "7. Covariates", intro="Covariates are typed and profiled only; nothing is imputed or scaled.")
        for art in self.store.find(modality=Modality.COVARIATES.value, artifact_type="standardized"):
            q = self.qc_summary(art.artifact_id) or {}
            variables = q.get("variables", {})
            rows = [{"variable": k, **v} for k, v in variables.items()]
            for r in rows:
                if isinstance(r.get("levels"), dict):
                    r["levels"] = ", ".join(f"{k}: {v}" for k, v in list(r["levels"].items())[:6])
            if rows:
                s.charts.append(self.chart(charts.bar([r["variable"] for r in rows], [100 * r["missing_fraction"] for r in rows],
                                                      f"{art.source_id}: missing values per variable", y="% missing", hover_unit="%")))
            s.tables.append(_table(f"{art.source_id} covariates", rows))
        if not s.tables:
            s.empty = "No covariates registered."
        return s

    def genotype_section(self) -> Section:
        s = Section("genotype", "8. Genotype")
        rows = []
        for art in self.store.find(modality=Modality.GENOTYPE.value, artifact_type="source"):
            bi = art.metadata.get("build_inference", {})
            rows.append({"source": art.source_id, "origin": art.origin.value, "format": art.format,
                         "storage": art.storage_mode, "participants": art.participant_count, "variants": art.feature_count,
                         "genome build": art.genome_build or "unknown", "build confidence": bi.get("confidence"),
                         "status": art.status.value, "path": art.path})
            if bi.get("evidence"):
                s.tables.append(_table(f"{art.source_id}: genome-build evidence", bi["evidence"]))
        s.tables.insert(0, _table("Genotype sources", rows))
        if not rows:
            s.empty = "No genotype registered (phenotype/omics-only project)."
        return s

    def genotype_qc_section(self) -> tuple[Section, Section, Section]:
        qc = Section("genotype-qc", "9. Genotype QC", intro="Phenotype-independent QC with PLINK 2; HWE uses all retained samples.")
        pop = Section("population", "10. Population structure", intro="QC principal components describe the cohort; refit PCs inside training folds for modelling.")
        rel = Section("relatedness", "11. Relatedness")
        thresholds = self.project.config.genotype.qc
        for art in self.store.find(artifact_type="qc_table"):
            df = pl.read_parquet(art.path)
            if art.artifact_name.endswith("sample_qc"):
                qc.charts.append(self.chart(charts.histogram(df.get_column("sample_missingness").to_list(), f"{art.source_id}: sample missingness",
                                                             x="fraction missing", threshold=thresholds.sample_missingness)))
                if df.get_column("het_f").drop_nulls().len():
                    qc.charts.append(self.chart(charts.histogram(df.get_column("het_f").to_list(), f"{art.source_id}: heterozygosity F", x="F")))
                reasons = df.filter(~pl.col("pass")).get_column("fail_reasons").str.split(";").explode().value_counts(sort=True)
                qc.tables.append(_table(f"{art.source_id}: sample failures by reason", reasons.rename({"fail_reasons": "reason"}).to_dicts()))
            else:
                qc.charts.append(self.chart(charts.histogram(df.get_column("MAF").to_list(), f"{art.source_id}: minor allele frequency",
                                                             x="MAF", threshold=thresholds.maf)))
                p = df.get_column("hwe_p").drop_nulls()
                if p.len():
                    qc.charts.append(self.chart(charts.histogram((-np.log10(np.clip(p.to_numpy(), 1e-300, 1))).tolist(),
                                                                 f"{art.source_id}: HWE -log10(p)", x="-log10 p",
                                                                 threshold=-np.log10(thresholds.hwe_p), log_y=True)))
                qc.charts.append(self.chart(charts.histogram(df.get_column("variant_missingness").to_list(),
                                                             f"{art.source_id}: variant missingness", x="fraction missing",
                                                             threshold=thresholds.variant_missingness)))
        for art in self.store.find(artifact_type="qc_genotype", include_superseded=False):
            m = art.metadata.get("qc_metrics", {})
            qc.tables.insert(0, _table(f"{art.source_id}: QC summary ({art.metadata.get('qc_status')})",
                                       [{"metric": k, "value": v} for k, v in m.items()]))
            qc.notes += [f"{art.source_id}: {msg}" for msg in art.metadata.get("messages", [])]
            if m:
                qc.charts.insert(0, self.chart(charts.status_bar(["samples", "variants"], [m.get("samples_pass", 0), m.get("variants_pass", 0)],
                                                                 [m.get("samples_fail", 0), m.get("variants_total", 0) - m.get("variants_pass", 0)],
                                                                 f"{art.source_id}: QC outcome")))
        if not qc.tables:
            qc.empty = "Genotype QC has not run (requires PLINK 2: `efgpp setup data`)."
        ancestry = {a.source_id: pl.read_parquet(a.path) for a in self.store.find(artifact_type="ancestry")}
        for art in self.store.find(artifact_type="qc_pca"):
            df = pl.read_parquet(art.path)
            if {"PC1", "PC2"} <= set(df.columns):
                groups = None
                anc = ancestry.get(art.source_id)
                if anc is not None:
                    df = df.join(anc.select("participant_id", "ancestry"), on="participant_id", how="left")
                    groups = df.get_column("ancestry").fill_null("unassigned").to_list()
                pop.charts.append(self.chart(charts.scatter(df.get_column("PC1").to_list(), df.get_column("PC2").to_list(),
                                                            f"{art.artifact_name}: PC1 vs PC2", groups=groups, xlabel="PC1", ylabel="PC2",
                                                            hover=df.get_column("participant_id").fill_null("?").to_list())))
            eig = art.metadata.get("eigenvalues")
            if eig and Path(eig).exists():
                e = pl.read_parquet(eig)
                pop.charts.append(self.chart(charts.bar(e.get_column("pc").to_list(), e.get_column("eigenvalue").to_list(),
                                                        f"{art.artifact_name}: eigenvalues", y="eigenvalue")))
        for a in ancestry.values():
            counts = a.get_column("ancestry").value_counts(sort=True)
            pop.tables.append(_table("Assigned ancestry", counts.to_dicts()))
        if not pop.charts:
            pop.empty = "No QC PCA available."
        for art in self.store.find(artifact_type="kinship"):
            df = pl.read_parquet(art.path)
            if df.height:
                rel.charts.append(self.chart(charts.histogram(df.get_column("KINSHIP").to_list(), f"{art.source_id}: KING kinship (reported pairs)",
                                                              x="kinship", threshold=thresholds.kinship_cutoff)))
                rel.tables.append(_table(f"{art.source_id}: pairs by relationship", df.get_column("relationship").value_counts(sort=True).to_dicts()))
        for art in self.store.find(artifact_type="roh"):
            rel.tables.append(_table(f"{art.source_id}: runs of homozygosity per sample", pl.read_parquet(art.path).head(50).to_dicts()))
        if not rel.tables:
            rel.empty = "No related pairs reported (or relatedness not computed)."
        return qc, pop, rel

    def omics_sections(self) -> list[Section]:
        out = []
        for i, modality in enumerate(OMICS_MODALITIES):
            s = Section(modality.value, f"{12 + i}. {modality.value.capitalize()}")
            rows = []
            for art in self.store.find(modality=modality.value, artifact_type="standardized"):
                q = self.qc_summary(art.artifact_id) or {}
                rows.append({"source": art.source_id, "origin": art.origin.value, "tissue": art.tissue, "format": art.format,
                             "participants": art.participant_count, "features": art.feature_count, "qc": q.get("status"),
                             "missing": q.get("missing_fraction"), "median": q.get("value_median"),
                             "normalization": art.metadata.get("normalization"), "measurement": art.metadata.get("measurement_type")})
            if rows:
                s.tables.append(_table(f"{modality.value} sources", rows))
            else:
                s.empty = f"No measured {modality.value} registered."
            out.append(s)
        return out

    def predicted_section(self) -> Section:
        s = Section("genotype_derived", "16. Genotype-Derived Molecular Representations",
                    intro="Derived from the genotype only (no phenotype is read): participant carrier variants, "
                          "consequence and gene burdens (origin: derived) and genetically predicted expression, "
                          "splicing, proteins, metabolites and methylation (origin: predicted). Predicted values "
                          "are genetically predicted components, not measurements.")
        for a in self.store.find(artifact_type="participant_variants"):
            md = a.metadata
            s.tables.append(_table(f"{a.source_id}: participant carrier table", [{
                "participants": a.participant_count, "variants": a.feature_count, "carrier rows": md.get("carrier_rows"),
                "counted allele": md.get("counted_allele"), "REF check": json.dumps(md.get("ref_check")),
                "dosage": md.get("dosage")}]))
        for a in self.store.find(artifact_type="consequence_counts"):
            counts = pl.read_parquet(a.path)
            for col, label in (("n_variants_total", "variants carried"), ("n_missense_variant", "missense"),
                               ("n_lof", "loss of function"), ("n_splice", "splice"),
                               ("n_damaging_missense", "damaging missense (AlphaMissense)")):
                if col in counts.columns and counts.height:
                    s.charts.append(self.chart(charts.histogram(counts.get_column(col).to_list(),
                                                                f"{a.source_id}: {label} per participant",
                                                                x="variants per participant")))
            s.notes.append(f"{a.source_id}: consequence counts are EFGPP-derived (VEP PICK / canonical / most severe "
                           "consequence, one per variant); not a published burden model.")
        rows = []
        for a in self.store.find(origin="predicted"):
            man_path = a.metadata.get("manifest")
            man = yaml.safe_load(Path(man_path).read_text(encoding="utf-8")) if man_path and Path(man_path).exists() else {}
            rows.append({"modality": str(a.modality), "provider": man.get("provider"), "dataset": man.get("dataset"),
                         "model version": (man.get("model_resource") or {}).get("version"), "tissue": a.tissue,
                         "platform": man.get("platform"), "training ancestry": man.get("training_ancestry"),
                         "ancestry match": man.get("ancestry_match_status"), "participants": a.participant_count,
                         "features": a.feature_count, "median coverage": man.get("median_coverage"),
                         "median validation R2": man.get("median_validation_r2"),
                         "low coverage": man.get("low_coverage_features"), "engine": man.get("engine")})
        if rows:
            s.tables.append(_table("Genetically predicted data", rows))
        from efgpp.data.predicted.compare import compare_observed_predicted

        qc_rows = []
        for modality in (Modality.EXPRESSION, Modality.SPLICING, Modality.PROTEOMICS, Modality.METABOLOMICS,
                         Modality.METHYLATION):
            try:
                qc_rows += [{"modality": modality.value, **r} for r in compare_observed_predicted(self.project, modality)]
            except Exception as exc:  # noqa: BLE001 - QC must never break the report
                s.notes.append(f"observed vs predicted {modality.value}: not compared ({exc})")
        if qc_rows:
            s.tables.append(_table("QC: observed vs genetically predicted (per-feature Pearson r)", qc_rows))
        if not s.tables and not s.charts:
            s.empty = ("Nothing derived yet: `efgpp data variants enable` / `efgpp data predict enable ...`, "
                       "then `efgpp data predict plan`.")
        return s

    def annotation_section(self) -> Section:
        s = Section("annotations", "17. Variant annotations")
        rows = []
        for a in self.store.find(modality=Modality.VARIANT_ANNOTATIONS.value):
            rows.append({"artifact": a.artifact_id, "name": a.artifact_name, "rows": a.feature_count, "tool": a.tool,
                         "versions": json.dumps(a.resource_versions)})
            cons = a.metadata.get("consequences")
            if cons:
                top = list(cons.items())[:15]
                s.charts.append(self.chart(charts.bar([k for k, _ in top], [v for _, v in top], f"{a.source_id}: VEP consequences (picked)",
                                                      horizontal=True)))
        for a in self.store.find(artifact_type="variant_table"):
            rows.append({"artifact": a.artifact_id, "name": a.artifact_name, "rows": a.feature_count, "tool": a.tool, "versions": ""})
        s.tables.append(_table("Variant tables and annotations", rows))
        if not rows:
            s.empty = "No variant annotations."
        return s

    def reference_section(self) -> tuple[Section, Section]:
        res = self.reg.rows("SELECT resource_id, name, version, genome_build, source, license, download_date, local_path FROM resources ORDER BY name")
        ref = Section("reference", "18. Reference knowledge")
        enabled = [{"resource": n, "enabled": "yes"} for n in self.project.resources.enabled_names()]
        ref.tables += [_table("Enabled in resources.yaml", enabled), _table("Installed resources", res)]
        gwas = [{"id": a.source_id, "trait": a.metadata.get("trait"), "variants": a.feature_count,
                 "build (input)": (a.metadata.get("gwaslab_report") or {}).get("build_detected"),
                 "lifted": (a.metadata.get("gwaslab_report") or {}).get("lifted"), "build": a.genome_build,
                 "gwaslab": a.tool_version, "path": a.path}
                for a in self.store.find(artifact_type="gwas_sumstats")]
        if gwas:
            ref.tables.append(_table("GWAS summary statistics (GWASLab)", gwas))
        ver = Section("resource-versions", "22. Resource versions")
        used = [{"artifact": a.artifact_id, "name": a.artifact_name, "resources": json.dumps(a.resource_versions)}
                for a in self.store.find() if a.resource_versions]
        ver.tables += [_table("Installed resource versions", res, ["name", "version", "genome_build", "download_date"]),
                       _table("Resource versions used by artifacts", used)]
        return ref, ver

    def missing_section(self) -> Section:
        s = Section("missing", "19. Missing data")
        rows = []
        for d in self.reg.rows("SELECT phenotype_id, count(*) AS n, sum(CAST(is_missing AS INTEGER)) AS missing FROM phenotype_observations GROUP BY 1"):
            rows.append({"source": d["phenotype_id"], "kind": "phenotype", "missing fraction": d["missing"] / d["n"] if d["n"] else None})
        for a in self.store.find(artifact_type="standardized"):
            q = self.qc_summary(a.artifact_id) or {}
            if "missing_fraction" in q and str(a.modality) != "phenotype":
                rows.append({"source": a.source_id, "kind": str(a.modality), "missing fraction": q["missing_fraction"]})
            for var, v in (q.get("variables") or {}).items():
                if isinstance(v, dict) and "missing_fraction" in v:
                    rows.append({"source": f"{a.source_id}.{var}", "kind": "covariate", "missing fraction": v["missing_fraction"]})
        if rows:
            s.charts.append(self.chart(charts.bar([r["source"] for r in rows], [100 * (r["missing fraction"] or 0) for r in rows],
                                                  "Missing values by source", horizontal=True, hover_unit="%")))
        s.tables.append(_table("Missingness", rows))
        return s

    def warnings_section(self) -> Section:
        s = Section("warnings", "20. Warnings")
        rows = []
        for f in sorted(self.project.path("qc").rglob("*_validation.json")):
            data = json.loads(f.read_text(encoding="utf-8"))
            for c in data["checks"]:
                if c["status"] in ("fail", "warn"):
                    rows.append({"level": "✗ error" if c["status"] == "fail" else "! warning", "source": data.get("source_id") or data["subject"],
                                 "message": c["message"]})
        for step in build_plan(self.project):
            if not step.enabled:
                rows.append({"level": "i not run", "source": step.id, "message": step.reason})
        s.tables.append(_table("Validation problems and steps that cannot run yet", rows))
        if not rows:
            s.empty = "No warnings."
        return s

    def software_section(self) -> Section:
        s = Section("software", "21. Software versions")
        rows = self.reg.rows("SELECT name, environment, version, path FROM software ORDER BY name")
        import importlib.metadata as md

        py = []
        for pkg in ("polars", "pyarrow", "duckdb", "pandera", "pydantic", "anndata", "zarr", "plotly"):
            try:
                py.append({"package": pkg, "version": md.version(pkg)})
            except md.PackageNotFoundError:
                continue
        s.tables += [_table("Scientific executables", rows), _table("Python packages", py)]
        return s

    def provenance_section(self) -> Section:
        s = Section("provenance", "23. Provenance")
        runs = self.reg.rows("SELECT run_id, step_id, tool, tool_version, exit_code, status, started_at, completed_at, log_path "
                             "FROM tool_runs ORDER BY run_id DESC LIMIT 200")
        arts = self.reg.rows("SELECT a.artifact_id, a.artifact_name, a.origin, a.status, a.tool, "
                             "string_agg(p.parent_artifact_id, ', ') AS parents FROM artifacts a "
                             "LEFT JOIN artifact_parents p ON p.artifact_id = a.artifact_id GROUP BY ALL ORDER BY a.artifact_id")
        s.tables += [_table("Artifacts and lineage", arts), _table("Tool runs (latest 200)", runs)]
        return s


def _pheno_summary(q: dict[str, Any]) -> str:
    if "cases" in q:
        return f"{q['cases']:,} cases / {q['controls']:,} controls"
    if "mean" in q:
        return f"mean {q['mean']:.4g}, sd {q['sd']:.4g}"
    if "classes" in q:
        return ", ".join(f"{k}: {v}" for k, v in q["classes"].items())
    return ""


def _run_multiqc(project: Project) -> Path | None:
    if not (project.config.reporting.multiqc and available(project, "multiqc")):
        return None
    out = project.report_root / "data" / "multiqc"
    tool = resolve(project, "multiqc")
    proc = subprocess.run(tool.command(str(project.path("qc")), "-o", str(out), "-f", "--quiet"),
                          capture_output=True, env=tool.environ())
    report = out / "multiqc_report.html"
    return report if proc.returncode == 0 and report.exists() else None


def build_report(project: Project) -> Path:
    with Registry.open(project) as reg:
        av = availability_mod.compute(project, reg)
        b = ReportBuilder(project, reg)
        qc, pop, rel = b.genotype_qc_section()
        ref, ver = b.reference_section()
        sections = [
            b.project_section(), b.participants_section(av), b.phenotypes_section(),
            b.timeline_section(timeline_summary(reg)), b.biospecimens_section(), b.availability_section(av),
            b.covariates_section(), b.genotype_section(), qc, pop, rel, *b.omics_sections(),
            b.predicted_section(), b.annotation_section(), ref, b.missing_section(), b.warnings_section(),
            b.software_section(), ver, b.provenance_section(),
        ]
    multiqc = _run_multiqc(project)
    if project.config.reporting.plotly_js == "inline":
        from plotly.offline import get_plotlyjs

        plotly_tag = f"<script>{get_plotlyjs()}</script>"
    else:
        plotly_tag = '<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>'
    env = Environment(loader=FileSystemLoader(TEMPLATES), autoescape=select_autoescape(["html"]))
    html = env.get_template("data_report.html.j2").render(
        project=project.config.project.name, version=__version__, generated=utcnow().isoformat(timespec="seconds"),
        sections=sections, participants=av.n_participants, plotly_tag=plotly_tag,
        dark_script=charts.dark_mode_script(),
        multiqc=multiqc.relative_to(project.report_root / "data").as_posix() if multiqc else None,
    )
    out = project.report_root / "data" / "index.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    return out
