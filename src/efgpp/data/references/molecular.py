"""Genetic prediction model resources for genotype-derived molecular data (origin REFERENCE).

    predictdb-gtex-v8-expression   GTEx v8 MASHR eQTL PredictDB models   Zenodo 3518299 mashr_eqtl.tar
    predictdb-gtex-v8-splicing     GTEx v8 MASHR sQTL PredictDB models   Zenodo 3518299 mashr_sqtl.tar
    predictdb-protein              PredictDB protein models (--dataset <file.db>), Zenodo 4837327
    omicspred                      OmicsPred genetic scores (--dataset OPDxxxxxx), resolved via its REST API
    mimosa                         MIMOSA whole-blood methylation models   Zenodo 8400313 MIMOSA-Models.zip
    hibag                          a published HIBAG HLA classifier (--url/--path with its metadata)
    spliceai-resources             records SpliceAI's inputs (FASTA, built-in annotation, licence)

Software is never installed here and resources are never installed by `efgpp setup`. Nothing
is invented: unknown metadata stays null, and a dataset must be chosen explicitly.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tarfile
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, ClassVar

import polars as pl

from efgpp.data.references.base import FetchedResource, ReferenceProvider
from efgpp.resources.downloader import download
from efgpp.resources.zenodo import GB, check_space, fetch_file, record_files

OMICSPRED_API = "https://rest.omicspred.org/api"
OMICS_TYPES = {"proteomics": "protein", "metabolomics": "metabolite", "transcriptomics": "gene expression",
               "protein": "protein", "metabolite": "metabolite", "expression": "gene expression"}
MODALITY_OF = {"protein": "proteomics", "metabolite": "metabolomics", "gene expression": "expression"}


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def molecular_root(project: Any) -> Path:
    return project.resource_root / "molecular_models"


def _extract(archive: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    if archive.name.endswith(".zip"):
        with zipfile.ZipFile(archive) as z:
            z.extractall(dest)
    else:
        with tarfile.open(archive) as tar:
            tar.extractall(dest, filter="data")


class _ZenodoArchive(ReferenceProvider):
    """A single Zenodo archive with a published checksum, extracted into raw/."""

    record: ClassVar[str]
    filename: ClassVar[str]
    md5: ClassVar[str]
    subdir: ClassVar[str]
    version_name: ClassVar[str]
    modality: ClassVar[str]
    license = "CC BY 4.0"
    citation: ClassVar[str] = ""
    context: ClassVar[dict[str, Any]] = {}
    options: dict[str, Any]

    def cache(self, version: str | None = None) -> Path:
        d = molecular_root(self.project) / self.subdir
        if version:
            d = d / version
        d.mkdir(parents=True, exist_ok=True)
        return d

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "zenodo_record": self.record, "file": self.filename, "md5": self.md5,
                "stored_under": str(molecular_root(self.project) / self.subdir)}

    def _download(self, extract_factor: float = 1.2) -> tuple[Path, str, Any]:
        say = getattr(self, "options", {}).get("progress")
        cache = molecular_root(self.project) / "_downloads"
        return fetch_file(self.record, self.filename, cache, client=self.client, expected_md5=self.md5,
                          extract_factor=extract_factor, say=say)

    def _manifest(self, zf: Any, sha: str, **extra: Any) -> dict[str, Any]:
        return {"name": self.name, "provider": self.name.split("-")[0], **zf.manifest(), "version": self.version_name,
                "downloaded_at": _now(), "sha256": sha, "download_sha256": sha, "modality": self.modality,
                "license": self.license, "citation": self.citation, **self.context, **extra}


class _PredictDBGTEx(_ZenodoArchive):
    record = "3518299"
    subdir = "predictdb/gtex_v8"
    citation = "Barbeira et al. 2021, Genome Biology 22:49 (GTEx v8 MASHR models); GTEx Consortium 2020"
    context: ClassVar[dict[str, Any]] = {"genome_build": "GRCh38", "training_cohort": "GTEx v8",
                                         "training_ancestry": "mostly European (GTEx v8 donors)",
                                         "platform": "RNA-seq (GTEx v8)", "tissue": "one model database per tissue"}

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        archive, sha, zf = self._download()
        raw = staging / "raw"
        _extract(archive, raw)
        dbs = sorted(raw.rglob("*.db"))
        if not dbs:
            raise RuntimeError(f"{self.filename} contains no PredictDB .db models")
        tissues = sorted(re.sub(r"^mashr_", "", d.stem) for d in dbs)
        manifest = self._manifest(zf, sha, tissues=tissues, n_tissues=len(tissues),
                                  models_dir=str(dbs[0].parent.relative_to(staging)),
                                  archive_removed_after_verification=True)
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
        archive.unlink()
        return FetchedResource(self.name, self.version_name, [raw, staging / "manifest.json"], dbs[0].parent,
                               zf.url, self.license, genome_build="GRCh38", metadata=manifest)


class PredictDBExpressionProvider(_PredictDBGTEx):
    name = "predictdb-gtex-v8-expression"
    filename = "mashr_eqtl.tar"
    md5 = "87f3470bf2676043c748b684fb35fa7d"
    version_name = "expression_mashr"
    modality = "expression"

    def cache(self, version: str | None = None) -> Path:
        d = molecular_root(self.project) / self.subdir / "expression_mashr"
        d.mkdir(parents=True, exist_ok=True)
        return d


class PredictDBSplicingProvider(_PredictDBGTEx):
    name = "predictdb-gtex-v8-splicing"
    filename = "mashr_sqtl.tar"
    md5 = "fa9167cfd2a9699f9f58ad04781e9576"
    version_name = "splicing_mashr"
    modality = "splicing"

    def cache(self, version: str | None = None) -> Path:
        d = molecular_root(self.project) / self.subdir / "splicing_mashr"
        d.mkdir(parents=True, exist_ok=True)
        return d


class PredictDBProteinProvider(ReferenceProvider):
    """PredictDB protein models (multi-population TOPMed/MESA; Zenodo 4837327) - one .db per --dataset."""

    name = "predictdb-protein"
    license = "CC BY 4.0"
    record = "4837327"
    options: dict[str, Any]

    def cache(self, version: str | None = None) -> Path:
        d = molecular_root(self.project) / "predictdb_protein"
        if version:
            d = d / version
        d.mkdir(parents=True, exist_ok=True)
        return d

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "zenodo_record": self.record, "note": "choose a .db with --dataset"}

    def datasets(self) -> list[str]:
        return [f.key for f in record_files(self.record, self.client) if f.key.endswith(".db")]

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        opts = getattr(self, "options", {})
        dataset = opts.get("dataset")
        if not dataset:
            raise RuntimeError("choose a model database with --dataset; available: " + ", ".join(self.datasets()))
        path, sha, zf = fetch_file(self.record, dataset, staging, client=self.client, say=opts.get("progress"))
        import sqlite3

        con = sqlite3.connect(path)
        try:
            ids = [r[0] for r in con.execute("SELECT varID FROM weights LIMIT 2000")]
        finally:
            con.close()
        suffix = {m.group(1) for v in ids if (m := re.search(r"_b(3[78])$", v or ""))}
        model_build = opts.get("model_build") or ({"38": "GRCh38", "37": "GRCh37"}[suffix.pop()] if len(suffix) == 1
                                                   else None)
        if model_build is None:
            raise RuntimeError(f"{dataset}: the genome build of its variant ids cannot be read; "
                               "pass --model-build GRCh37|GRCh38 as documented by the model's authors")
        ancestry = dataset.split("_", 1)[0]
        manifest = {"name": self.name, "provider": "predictdb", **zf.manifest(), "version": Path(dataset).stem,
                    "downloaded_at": _now(), "sha256": sha, "download_sha256": sha, "genome_build": model_build,
                    "modality": "proteomics", "platform": None, "training_cohort": "TOPMed MESA",
                    "training_ancestry": ancestry, "tissue": "plasma", "license": self.license,
                    "citation": "Schubert et al. 2022 (protein prediction for trait mapping in diverse populations)"}
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
        return FetchedResource(self.name, Path(dataset).stem, [path, staging / "manifest.json"], path, zf.url,
                               self.license, genome_build=model_build, metadata=manifest)


# ------------------------------------------------------------------ OmicsPred
def omicspred_catalog_path(project: Any) -> Path:
    return molecular_root(project) / "omicspred" / "catalog" / "catalog.json"


def omicspred_refresh(project: Any, client: Any) -> dict[str, Any]:
    """Download the OmicsPred dataset catalogue (REST API) and cache it with retrieval metadata."""
    import hashlib

    url: str | None = f"{OMICSPRED_API}/dataset/all"
    results: list[dict[str, Any]] = []
    digest = hashlib.sha256()
    count = None
    version = None
    while url:  # follow the API's own pagination links
        r = client.get(url)
        r.raise_for_status()
        digest.update(r.content)
        data = r.json()
        if not isinstance(data, dict) or "results" not in data:
            raise RuntimeError("OmicsPred API answered in an unexpected format; the cached catalogue was kept")
        results += data["results"]
        count = data.get("count", count)
        version = version or r.headers.get("x-api-version")
        url = data.get("next")
    out = {"retrieved_at": _now(), "api": f"{OMICSPRED_API}/dataset/all", "sha256": digest.hexdigest(),
           "api_version": version, "count": count, "results": results}
    path = omicspred_catalog_path(project)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


def omicspred_catalog(project: Any, client: Any | None = None) -> dict[str, Any]:
    path = omicspred_catalog_path(project)
    if path.exists():
        return dict(json.loads(path.read_text(encoding="utf-8")))
    if client is None:
        raise RuntimeError("no OmicsPred catalogue yet: run `efgpp resources omicspred refresh`")
    return omicspred_refresh(project, client)


def _training(ds: dict[str, Any]) -> tuple[str | None, str | None, int | None]:
    samples = ds.get("samples_training") or []
    cohorts = sorted({c.get("name_short") for s in samples for c in (s.get("cohorts") or []) if c.get("name_short")})
    ancestry = sorted({s.get("ancestry_broad") for s in samples if s.get("ancestry_broad")})
    n = sum(int(s.get("sample_number") or 0) for s in samples) or None
    return ", ".join(cohorts) or None, ", ".join(ancestry) or None, n


def omicspred_rows(catalog: dict[str, Any], *, modality: str | None = None, cohort: str | None = None,
                   platform: str | None = None, ancestry: str | None = None) -> list[dict[str, Any]]:
    want = OMICS_TYPES.get(modality.lower(), modality.lower()) if modality else None
    rows = []
    for ds in catalog.get("results", []):
        coh, anc, n = _training(ds)
        plat = (ds.get("platform") or {}).get("name")
        if want and (ds.get("omics_type") or "").lower() != want:
            continue
        if cohort and cohort.lower() not in (coh or "").lower() and cohort.lower() not in (ds.get("name") or "").lower():
            continue
        if platform and platform.lower() not in (plat or "").lower():
            continue
        if ancestry and ancestry.lower() not in (anc or "").lower():
            continue
        rows.append({"id": ds.get("id"), "name": ds.get("name"), "omics_type": ds.get("omics_type"),
                     "platform": plat, "tissue": (ds.get("tissue") or {}).get("label"),
                     "scores": ds.get("scores_count"), "training_cohort": coh, "training_ancestry": anc,
                     "training_n": n, "method": ds.get("method_name")})
    return rows


def read_validation(path: Path) -> pl.DataFrame | None:
    """OmicsPred validation table (TSV: OMICSPRED ID, Internal_R2, <cohort>_R2, ...) when readable."""
    try:
        df = pl.read_csv(path, separator="\t", infer_schema=False)
    except Exception:  # noqa: BLE001 - some datasets publish spreadsheets; kept raw, R2 left null
        return None
    id_col = next((c for c in df.columns if c.replace(" ", "").upper() == "OMICSPREDID"), None)
    if id_col is None:
        return None
    r2_cols = [c for c in df.columns if c.upper().endswith("_R2")]
    main = next((c for c in r2_cols if c.lower().startswith("internal")), r2_cols[0] if r2_cols else None)
    if main is None:
        return None
    return df.select(pl.col(id_col).alias("feature_id"), pl.col(main).cast(pl.Float64, strict=False)
                     .alias("validation_r2"), pl.lit(main.removesuffix("_R2")).alias("validation_cohort"))


class OmicsPredProvider(ReferenceProvider):
    """OmicsPred genetic scores for one dataset (proteomics, metabolomics, transcriptomics)."""

    name = "omicspred"
    license = "CC BY 4.0 (OmicsPred genetic scores)"
    homepage = "https://www.omicspred.org/"
    options: dict[str, Any]

    def cache(self, version: str | None = None) -> Path:
        d = molecular_root(self.project) / "omicspred" / "datasets"
        if version:
            d = d / version
        d.mkdir(parents=True, exist_ok=True)
        return d

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "api": OMICSPRED_API, "note": "choose --dataset from `efgpp resources omicspred list`"}

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        from efgpp.data.predicted.models import FEATURE_SCHEMA, read_scoring_file, standardize

        opts = getattr(self, "options", {})
        dataset = opts.get("dataset")
        if not dataset:
            raise RuntimeError("choose a dataset with --dataset (see `efgpp resources omicspred list`)")
        catalog = omicspred_refresh(self.project, self.client)
        ds = next((d for d in catalog["results"] if d.get("id") == dataset), None)
        if ds is None:
            raise RuntimeError(f"{dataset} is not in the OmicsPred catalogue ({catalog['count']} datasets)")
        target = build or self.project.config.defaults.target_build
        urls = ds.get("scoring_files_urls") or {}
        key = "scoring_files_hm_38" if target == "GRCh38" and urls.get("scoring_files_hm_38") else "scoring_files"
        if not urls.get(key):
            raise RuntimeError(f"{dataset}: the OmicsPred API lists no scoring files")
        say = opts.get("progress")
        if say:
            say(f"{dataset} {ds.get('name')}: {ds.get('scores_count')} {ds.get('omics_type')} scores "
                f"({(ds.get('platform') or {}).get('name')}, {(ds.get('tissue') or {}).get('label')})")
        raw = staging / "raw"
        raw.mkdir(parents=True)
        (raw / "dataset.json").write_text(json.dumps(ds, indent=1), encoding="utf-8")
        archive = raw / f"{dataset}_{key}.zip"
        sha = download(urls[key], archive, client=self.client, timeout=600)
        validation = None
        if urls.get("validation_results"):
            vpath = raw / f"{dataset}_validation_results"
            try:
                download(urls["validation_results"], vpath, client=self.client, timeout=600)
                validation = read_validation(vpath)
            except Exception:  # noqa: BLE001 - validation metadata is optional
                validation = None
        coh, anc, n = _training(ds)
        plat = (ds.get("platform") or {}).get("name")
        tissue = (ds.get("tissue") or {}).get("label")
        pub = ds.get("publication") or {}
        modality = MODALITY_OF.get((ds.get("omics_type") or "").lower(), ds.get("omics_type"))
        frames, feats, builds = [], [], set()
        with zipfile.ZipFile(archive) as z:
            members = [m for m in z.namelist() if not m.endswith("/") and m.endswith((".txt", ".txt.gz", ".tsv"))]
            tmp = staging / "_scores"
            for m in members:
                out = tmp / Path(m).name
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_bytes(z.read(m))
                header, w = read_scoring_file(out)
                builds.add(header["_genome_build"])
                fid = header.get("omicspred_id") or out.name.split("_")[0].split(".")[0]
                frames.append(w.with_columns(pl.lit(fid).alias("feature_id"),
                                             pl.lit(header.get("trait_reported") or header.get("pgs_name"))
                                             .alias("feature_name")))
                feats.append({"feature_id": fid, "feature_name": header.get("trait_reported") or header.get("pgs_name"),
                              "gene_name": header.get("pgs_name"), "n_variants": w.height})
            shutil.rmtree(tmp, ignore_errors=True)
        if len(builds) != 1 or "unknown" in builds:
            raise RuntimeError(f"{dataset}: scoring files report genome builds {sorted(builds)}; refusing to guess")
        model_build = builds.pop()
        model_id = f"omicspred:{dataset}"
        weights = pl.concat(frames, how="vertical_relaxed")
        if validation is not None:
            weights = weights.join(validation.select("feature_id", "validation_r2", "validation_cohort"),
                                   on="feature_id", how="left")
        weights = standardize(weights, provider="omicspred", provider_dataset_id=dataset, model_id=model_id,
                              modality=modality, genome_build=model_build, tissue=tissue, sample_type=tissue,
                              platform=plat, training_cohort=coh, training_ancestry=anc, training_n=n,
                              source_publication=pub.get("doi"), model_version=f"{dataset}@{catalog['retrieved_at']}")
        features = pl.DataFrame(feats)
        if validation is not None:
            features = features.join(validation, on="feature_id", how="left")
        features = standardize(features.with_columns(pl.lit(model_id).alias("model_id")), FEATURE_SCHEMA)
        std = staging / "standardized"
        std.mkdir()
        weights.write_parquet(std / "weights.parquet")
        features.write_parquet(std / "features.parquet")
        manifest = {"name": "omicspred", "provider": "omicspred", "provider_record": dataset,
                    "version": dataset, "source_filename": archive.name, "source_url": urls[key],
                    "source_checksum": None, "sha256": sha, "download_sha256": sha,
                    "downloaded_at": _now(), "catalog_retrieved_at": catalog["retrieved_at"],
                    "genome_build": model_build, "modality": modality, "tissue": tissue, "platform": plat,
                    "training_cohort": coh, "training_ancestry": anc, "training_n": n,
                    "license": ds.get("license") or self.license, "citation": pub.get("doi"),
                    "n_scores": len(feats), "validation_r2_available": validation is not None,
                    "scoring_file_kind": key}
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
        return FetchedResource(self.name, dataset, [raw, std, staging / "manifest.json"], std, urls[key],
                               manifest["license"], genome_build=model_build, metadata=manifest)


# ------------------------------------------------------------------ MIMOSA
def rscript_for(project: Any, envs: tuple[str, ...] = ("methylation", "r")) -> Path | None:
    for env in envs:
        for base in (project.envs_dir, *(d / "envs" for d in project.legacy_software_dirs)):
            p = base / env / "bin" / "Rscript"
            if p.is_file():
                return p
    found = shutil.which("Rscript")
    return Path(found) if found else None


class MIMOSAProvider(_ZenodoArchive):
    name = "mimosa"
    record = "8400313"
    filename = "MIMOSA-Models.zip"
    md5 = "2a4be3b8af1ca1d489e8e2df8fba515e"
    subdir = "mimosa/v2"
    version_name = "whole_blood_v2"
    modality = "methylation"
    citation = ("Melton HJ, Zhang Z, Deng H, Wu L, Wu C (2023) MIMOSA: a resource consisting of improved "
                "methylome imputation models increases power to identify CpG site-phenotype associations. "
                "doi:10.5281/zenodo.8400313")
    context: ClassVar[dict[str, Any]] = {
        "tissue": "whole blood", "platform": "Illumina 450K CpG sites", "training_cohort": "GoDMC mQTL summary data",
        "training_ancestry": "European (GoDMC / Framingham Heart Study)",
        "validation_cohort": "Framingham Heart Study (test data)",
        "note": "DNA methylation prediction in whole blood; not brain, liver or tumour methylation; "
                "values are genetically predicted methylation scores, not beta values"}

    def cache(self, version: str | None = None) -> Path:
        d = molecular_root(self.project) / "mimosa" / "v2"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        from efgpp.constants import GenomeBuild
        from efgpp.data.genotype.liftover import Lifter, chain_file
        from efgpp.data.predicted.mimosa import EXPORT_SCRIPT, read_export, standardized

        rscript = rscript_for(self.project)
        if rscript is None:
            raise RuntimeError("MIMOSA models are .rds files: install R first (`efgpp setup toolkit methylation`)")
        target = build or self.project.config.defaults.target_build
        chain = None
        if target == "GRCh38":
            chain = chain_file(self.project.resource_root, GenomeBuild.GRCH37, GenomeBuild.GRCH38)
        archive, sha, zf = self._download(extract_factor=2.5)
        raw = staging / "raw"
        check_space(staging, int(3.5 * GB), "MIMOSA extraction")
        _extract(archive, raw)
        export = staging / "mimosa_export.tsv"
        structure = staging / "rds_structure.txt"
        proc = subprocess.run([str(rscript), str(EXPORT_SCRIPT), str(raw), str(export), str(structure)],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"MIMOSA conversion failed:\n{proc.stderr[-2000:]}")
        table = read_export(export)
        minimum_r2 = self.project.data.predicted.methylation.minimum_model_r2 or 0.005
        model_build = "GRCh37"  # GoDMC mQTL positions (hg19)
        lifted_from = None
        if target == "GRCh38" and chain is not None and chain.exists():
            lifter = Lifter(chain)
            pos = table.select("chromosome", "position").unique().drop_nulls()
            lifted = {(c, p): lifter.lift(str(c), int(p)) for c, p in pos.iter_rows()}
            def new_pos(s: dict[str, Any]) -> int | None:
                hit = lifted.get((s["chromosome"], s["position"]))
                # only same-chromosome, forward-strand lifts; anything else is left unplaced (excluded, reported)
                if hit is None or hit[0] != str(s["chromosome"]).removeprefix("chr") or hit[2] != "+":
                    return None
                return hit[1]

            table = table.with_columns(pl.struct("chromosome", "position").map_elements(
                new_pos, return_dtype=pl.Int64).alias("position"))
            model_build, lifted_from = "GRCh38", "GRCh37"
        weights, features = standardized(table, minimum_r2, genome_build=model_build, version=self.version_name)
        std = staging / "standardized"
        (std / "weights").mkdir(parents=True)
        for (chrom,), part in weights.partition_by("chromosome", as_dict=True).items():
            part.write_parquet(std / "weights" / f"chromosome={chrom}.parquet")
        features.write_parquet(std / "models.parquet")
        shutil.move(str(export), raw / "mimosa_export.tsv")
        shutil.move(str(structure), raw / "rds_structure.txt") if structure.exists() else None
        manifest = self._manifest(zf, sha, genome_build=model_build, lifted_from=lifted_from,
                                  minimum_model_r2=minimum_r2, n_cpgs=features.height,
                                  n_cpgs_above_threshold=int((~features.get_column("below_threshold")).sum()),
                                  selection="highest test R2 among valid models; ties ElNet, MNet, SCAD, MCP, LASSO",
                                  archive_removed_after_verification=True)
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
        archive.unlink()
        return FetchedResource(self.name, self.version_name, [raw, std, staging / "manifest.json"], std, zf.url,
                               self.license, genome_build=model_build, metadata=manifest)


# ------------------------------------------------------------------ HIBAG and SpliceAI
class HIBAGProvider(ReferenceProvider):
    """A published HIBAG classifier (.RData). Its ancestry, platform, build and loci must be stated."""

    name = "hibag"
    license = "see the classifier's publication"
    options: dict[str, Any]

    def cache(self, version: str | None = None) -> Path:
        d = molecular_root(self.project) / "hibag"
        if version:
            d = d / version
        d.mkdir(parents=True, exist_ok=True)
        return d

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "needs": "--url or --path, --ancestry, --platform, --model-build, --loci"}

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        opts = getattr(self, "options", {})
        missing = [k for k in ("ancestry", "platform", "model_build", "loci") if not opts.get(k)]
        if not (opts.get("url") or opts.get("path")) or missing:
            raise RuntimeError("a HIBAG classifier needs --url or --path and its metadata: "
                               f"missing {', '.join(missing) or 'url/path'} (never one model for everyone)")
        src = opts.get("path")
        if src:
            dest = staging / Path(src).name
            shutil.copy2(src, dest)
            from efgpp.resources.downloader import sha256_file

            sha, source = sha256_file(dest), str(src)
        else:
            dest = staging / str(opts["url"]).rsplit("/", 1)[-1]
            sha, source = download(str(opts["url"]), dest, client=self.client, timeout=600), str(opts["url"])
        version = re.sub(r"[^\w.\-]", "_", dest.stem)
        manifest = {"name": self.name, "provider": "hibag", "source_filename": dest.name, "source_url": source,
                    "sha256": sha, "download_sha256": sha, "downloaded_at": _now(), "version": version,
                    "ancestry": opts["ancestry"], "platform": opts["platform"], "genome_build": opts["model_build"],
                    "loci": opts["loci"], "modality": "hla", "license": opts.get("license"),
                    "citation": opts.get("citation")}
        (staging / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
        return FetchedResource(self.name, version, [dest, staging / "manifest.json"], dest, source,
                               opts.get("license") or self.license, genome_build=opts["model_build"], metadata=manifest)


class SpliceAIResourcesProvider(ReferenceProvider):
    """Records what SpliceAI runs with: the project FASTA, its built-in gene annotation, the licence.

    Nothing is downloaded: SpliceAI's trained models and GRCh37/GRCh38 annotation ship with the
    package (`efgpp setup toolkit spliceai`); Illumina's precomputed scores need a BaseSpace
    account and are not fetched."""

    name = "spliceai-resources"
    license = "SpliceAI code GPL-3.0; trained models CC BY-NC 4.0 (non-commercial)"

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "downloads": "none (models ship with the spliceai package)"}

    def fetch(self, *, build: str | None = None, staging: Path) -> FetchedResource:
        from efgpp.data.references.genome import installed_fasta
        from efgpp.setup.tools import ToolNotFoundError, resolve

        fasta = installed_fasta(self.project)
        if fasta is None:
            raise RuntimeError("SpliceAI needs the reference FASTA: efgpp resources install genome")
        try:
            tool = resolve(self.project, "spliceai")
        except ToolNotFoundError as exc:
            raise RuntimeError("SpliceAI is not installed: efgpp setup toolkit spliceai") from exc
        target = build or self.project.config.defaults.target_build
        info = {"name": self.name, "fasta": str(fasta), "genome_build": target,
                "annotation": f"spliceai built-in {'grch38' if target == 'GRCh38' else 'grch37'}",
                "spliceai": str(tool.path), "license": self.license, "recorded_at": _now()}
        out = staging / "spliceai_resources.json"
        out.write_text(json.dumps(info, indent=1), encoding="utf-8")
        return FetchedResource(self.name, target, [out], out, "local", self.license, genome_build=target,
                               metadata=info)
