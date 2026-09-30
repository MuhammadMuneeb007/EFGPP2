# EFGPP Data Module — Technical Document

This document maps **every operation** in the EFGPP Data layer to the **package or tool that
performs it**, the EFGPP module that calls it, and what it produces.

Paths such as `data/genotype/qc.py` are relative to `src/efgpp/`.

---

## 1. Architecture at a glance

```text
                    project.yaml   data.yaml   resources.yaml          (Pydantic)
                           │           │             │
                           ▼           ▼             ▼
  efgpp CLI (Typer + Rich) ──► DataManager / plan ──► executor ──► Snakemake | built-in
                                      │                                  │
             ┌────────────────────────┼──────────────────────────────────┤
             ▼                        ▼                                  ▼
     Adapters (Polars,          Genotype steps (PLINK 2)      Annotation / prediction
     Pandera, AnnData/Zarr)     QC · PCA · KING · ROH          (VEP, OpenCRAVAT, DuckDB
             │                        │                        lookups, MetaXcan)
             └──────────────┬─────────┴──────────────┬──────────────────┘
                            ▼                        ▼
               Registry (DuckDB + Parquet)    Provenance logs (.efgpp/logs)
                            │
                            ▼
      Availability ─► Report (Plotly + Jinja2 [+ MultiQC]) ─► DataSnapshot (YAML + Parquet)
```

EFGPP decides **what** must happen (the plan). Snakemake, or the built-in executor when Snakemake
is not installed, decides **when** and **where** (order, parallelism, restarts, SLURM).

---

## 2. Package stack and responsibilities

EFGPP itself is installed into one environment: mamba/conda (`mamba env create -f environment.yml`,
then `conda activate efgpp`) or uv (`uv sync`). The external scientific tools in §2.2 live in
separate per-purpose environments inside each project. EFGPP activates the right one per step,
so users never switch environments by hand. See README → *Install*.

### 2.1 Python packages (installed by `environment.yml` or `uv sync`)

| Package | Used for | Main EFGPP modules |
|---|---|---|
| **pydantic** | Validating `project.yaml`, `data.yaml` and `resources.yaml` (unknown keys rejected); the `Artifact` model | `config/*.py`, `data/artifacts.py` |
| **PyYAML** | Reading/writing configuration, snapshots, manifests, run records, the lock file | `config/models.py`, `data/snapshots.py`, `setup/lock.py` |
| **typer** + **rich** | Command-line interface; tables and ✓/✗/! output | `cli/*.py` |
| **polars** | Reading user tables (all columns as text so IDs keep leading zeros), typing, reshaping, Parquet I/O | `data/io.py`, all adapters, QC parsers |
| **pyarrow** | Arrow tables handed from Polars to DuckDB; Parquet engine | `data/registry.py` |
| **duckdb** | Embedded registry database; bulk joins of cohort variants against ClinVar/gnomAD/AlphaMissense files; Parquet exports of registry tables | `data/registry.py`, `data/annotation/__init__.py` |
| **pandera** (`pandera.polars`) | Contract schemas for every *standardized* table (participants, aliases, phenotype observations, covariates, events, biospecimens, variants) | `data/schemas/*.py` |
| **numpy** | Omics matrices, simulation, ancestry centroid distances, bit-packing PLINK BED files | `data/adapters/expression.py`, `data/simulate.py`, `data/genotype/ancestry.py` |
| **scipy** | Available to adapters for statistics (declared dependency) | — |
| **anndata** / **zarr** / **mudata** *(extra `omics`)* | Standardized omics storage (AnnData in a Zarr store); reading `.h5ad` and Zarr inputs | `data/adapters/expression.py` |
| **plotly** | Interactive charts in the HTML report | `reporting/charts.py` |
| **jinja2** | HTML report template | `reporting/data_report.py`, `reporting/templates/` |
| **httpx** | HTTP for downloads and reference APIs (GWAS Catalog, PGS Catalog, Open Targets, GTEx) | `resources/downloader.py`, `data/references/*.py` |
| **tenacity** | Retries with exponential back-off for downloads | `resources/downloader.py` |
| **filelock** | Cross-process lock around the DuckDB registry (DuckDB allows one writer) | `data/registry.py` |
| **psutil** | Peak memory / CPU time of external tools; total RAM detection | `data/provenance.py`, `setup/platform.py` |
| **platformdirs**, **fsspec** | Declared for cache locations and remote file systems | — |
| **hashlib** (stdlib) | SHA-256 artifact identity | `resources/checksums.py` |
| **blake3** *(extra `fast`)* | Optional fast local integrity hashing | `resources/checksums.py` |
| **sqlite3** (stdlib) | Reading PredictDB model metadata (`.db` files) | `data/predicted/metaxcan.py` |
| **phenopackets** *(extra `standards`)* | Validating exported GA4GH Phenopackets | `data/phenopackets.py` |
| **ga4gh.vrs** *(extra `standards`)* | Reserved for VRS identifiers (`vrs_id` column exists; not yet populated) | `data/schemas/artifact.py` |
| **alphagenome** *(extra `alphagenome`)* | AlphaGenome API mode | `data/annotation/alphagenome.py` |

Development: **pytest**, **pytest-cov**, **hypothesis** (property tests for checksums),
**ruff**, **mypy**, **pre-commit**.

### 2.2 External scientific tools

Each tool lives in its own small environment under `.efgpp/envs/<env>/` (Pixi preferred; if
Pixi is unavailable or fails, the fastest conda-family tool found: mamba → micromamba → conda;
set with `execution.environment_manager` / `fallback_environment_manager`; definitions in
`workflow/envs/*.yaml`). `setup/tools.py` resolves each tool
in this order: `project.yaml: execution.tools` override → `.efgpp/bin` → its environment → `PATH`.

| Tool | Environment | Used for | Called from |
|---|---|---|---|
| **PLINK 2** | `genetics` | Missingness, allele frequency, HWE, heterozygosity, LD pruning, KING kinship, sex check, QC-passed PGEN, PCA, format conversion, liftover (`--update-map`), ancestry projection (`--score`) | `data/genotype/qc.py`, `pca.py`, `conversion.py`, `harmonization.py`, `ancestry.py`, `roh.py` |
| **PLINK 1.9** | `genetics` | Runs of homozygosity (`--homozyg`, not available in PLINK 2) | `data/genotype/roh.py` |
| **FlashPCA2** | `genetics` | Optional high-scale PCA engine | `data/genotype/pca.py` |
| **bcftools / htslib (tabix, bgzip)** | `genetics` | Installed and version-checked for VCF/BCF work | `setup/manager.py` |
| **Ensembl VEP** | `annotation` | Molecular consequences (offline cache) | `data/annotation/vep.py` |
| **OpenCRAVAT** (`oc`) | `annotation` | Modular annotation aggregation | `data/annotation/opencravat.py` |
| **MetaXcan / PrediXcan** (`Predict.py`) | `metaxcan` | Genetically predicted expression from PredictDB models | `data/predicted/metaxcan.py` |
| **Snakemake** (+ `snakemake-executor-plugin-slurm`) | `core` | Scheduling, parallelism, restarts, SLURM execution | `workflow/snakemake.py` |
| **MultiQC** | `reporting` | Aggregates QC tables (`*_mqc.tsv`) next to the EFGPP report | `reporting/data_report.py` |
| **Pixi** / **mamba** / micromamba / conda | system | Creating the per-purpose environments from conda-forge + Bioconda | `setup/pixi.py`, `setup/micromamba.py` |
| **Apptainer / Docker** | system | Optional container runtime for Linux-only tools | `setup/containers.py` |

---

## 3. Lifecycle: operation → package

```text
REGISTER → INSPECT → VALIDATE → ID RESOLUTION → TIMELINE/BIOSPECIMENS → STANDARDIZE
→ GENOME BUILD → QC → DERIVE → ANNOTATE → PREDICT → AVAILABILITY → REPORT → FREEZE
```

| # | Operation | Performed by | EFGPP module | Output |
|---|---|---|---|---|
| 1 | Create project layout and configs | pathlib, Pydantic, PyYAML | `project.py` (`efgpp init`) | `project.yaml`, `data.yaml`, `resources.yaml`, folders, `workflow/envs`, profiles |
| 2 | Register a source (storage policy) | shutil (copy), os.symlink (link), nothing (reference); hashlib SHA-256 | `data/storage.py`, `data/adapters/base.py` | `artifacts` row (`REGISTERED`), copied/linked files under `data/observed/…` |
| 3 | Inspect sources | Polars (tables), sidecar parsers (PLINK/VCF/BGEN), AnnData | `adapter.inspect()` | Console table |
| 4 | Validate (collect all errors) | Polars checks per adapter; Pydantic for config | `adapter.validate()`, `data/validation.py`, `data/manager.py` | `qc/<modality>/<id>_validation.json`; status `VALIDATED` |
| 5 | Resolve participant IDs | Polars + alias files | `data/aliases.py`, `data/participants.py` | `participants`, `sample_aliases` tables |
| 6 | Timeline and biospecimens | Polars + DuckDB | `data/timeline.py`, `data/biospecimens.py` | `events`, `biospecimens`, `assays` tables |
| 7 | Standardize tables | Polars → Parquet (zstd), Pandera contracts | phenotype/covariate/clinical adapters | `phenotypes/<id>/observations.parquet`, `data/observed/covariates/…` |
| 8 | Standardize omics | AnnData + Zarr (fallback Parquet) | `data/adapters/expression.py` (+ methylation/proteomics/metabolomics) | `data/<origin>/<modality>/<id>/<id>.zarr` |
| 9 | Genome-build inference | Python: contig lengths, coordinate bounds, FASTA REF concordance (own `.fai` reader), optional rsID panel | `data/genotype/build.py` | `build_inference` metadata; validation fails when uncertain |
| 10 | Phenotype / covariate / omics QC | Polars, NumPy | `adapter.qc()` | `qc/<modality>/<id>_qc_summary.json`, `qc_runs`, `qc_metrics` |
| 11 | Genotype QC | **PLINK 2** | `data/genotype/qc.py` | `qc/genotype/<id>/sample_qc.parquet`, `variant_qc.parquet`, QC-passed PGEN in `data/derived/genotype_qc/` |
| 12 | Relatedness | **PLINK 2** `--make-king-table`; Python greedy removal | `data/genotype/qc.py`, `relatedness.py` | `data/derived/kinship/<id>_kinship.parquet` |
| 13 | QC PCA | **PLINK 2** `--pca` (or **FlashPCA2**) | `data/genotype/pca.py` | `data/derived/qc_pca/<id>_pca_<engine>.parquet` |
| 14 | Runs of homozygosity | **PLINK 1.9** `--homozyg` | `data/genotype/roh.py` | `data/derived/roh/…` |
| 15 | Ancestry | **PLINK 2** `--pca allele-wts` + `--score`; NumPy nearest centroid | `data/genotype/ancestry.py` | `data/derived/ancestry/<id>_ancestry.parquet` |
| 16 | Variant table | Polars (from `.pvar`/`.bim`/VCF) | `data/adapters/genotype.py` | `data/derived/variant_annotations/<id>/<id>_variants.parquet` |
| 17 | Consequence annotation | **Ensembl VEP** | `data/annotation/vep.py` | `<id>_vep.parquet` |
| 18 | Aggregated annotation | **OpenCRAVAT** | `data/annotation/opencravat.py` | `<id>_opencravat.parquet` |
| 19 | ClinVar / gnomAD / dbSNP / AlphaMissense | **DuckDB** `read_csv` joins on chrom/pos/ref/alt | `data/annotation/clinvar.py`, `gnomad.py`, `alphamissense.py` | `<id>_<source>.parquet` |
| 20 | AlphaGenome | Atlas table (Polars join) or **alphagenome** client | `data/annotation/alphagenome.py` | `<id>_alphagenome.parquet` |
| 21 | Predicted expression | **PLINK 2** (export VCF) → **MetaXcan** `Predict.py` | `data/predicted/metaxcan.py` | `data/predicted/expression/<tissue>/…` (origin `predicted`) |
| 22 | Availability | DuckDB + Polars | `data/availability.py` | `registry/availability.parquet`, `availability_columns.parquet` |
| 23 | Report | Plotly + Jinja2 (+ MultiQC) | `reporting/data_report.py`, `reporting/charts.py` | `reports/data/index.html` |
| 24 | Freeze snapshot | DuckDB `COPY … TO PARQUET`, PyYAML, SHA-256 re-verification | `data/snapshots.py` | `snapshots/<name>.yaml`, `snapshots/<name>/*.parquet` (read-only) |
| 25 | Lock file | importlib.metadata, registry | `setup/lock.py` | `efgpp.lock.yaml` |

---

## 4. Genotype QC in detail (PLINK 2)

All passes read the original fileset. The original is never modified. HWE is computed on all
retained samples, so no phenotype is consulted.

| Pass | PLINK 2 arguments | Parsed by | Decision (defaults in `project.yaml: genotype.qc`) |
|---|---|---|---|
| 1 Sample missingness | `--missing sample-only` | `read_plink_table(.smiss)` | fail if `F_MISS > sample_missingness` (0.02) |
| 2 Variant metrics | `--keep keep1 --missing variant-only --freq --hardy` | `.vmiss`, `.afreq`, `.hardy`, `.hardy.x` | fail if `F_MISS > 0.02`, `MAF < 0.01`, `HWE p < 1e-6` |
| 3 LD pruning | `--extract pass --autosome --indep-pairwise 200 50 0.2` | `.prune.in` | pruned set for passes 4–5 and PCA |
| 4 Heterozygosity | `--extract pruned --het` | `.het` | fail if \|F − mean\| > 3 SD |
| 5 Kinship | `--make-king-table --king-table-filter 0.0442` | `.kin0` → `classify_pairs` | duplicates (> 0.354) removed; related (> 0.0884) flagged, removed only if `remove_related: true` |
| 6 Sex check | `--check-sex` | `.sexcheck` | `PROBLEM` fails; skipped with a warning when chrX data are absent |
| 7 QC-passed genotype | `--keep final --extract pass --make-pgen` | — | new DERIVED artifact `<id>_qc` |

Every PLINK run adds `--threads N --memory M`. Its command, version, exit code, wall/CPU time
and peak memory are recorded (§6).

---

## 5. Where each kind of data is stored

| Origin | Location | Format |
|---|---|---|
| observed | `data/observed/<modality>/<id>/` (copy/link) or original path (reference) | as supplied; standardized tables in Parquet, omics in AnnData/Zarr |
| derived | `data/derived/{genotype_qc,qc_pca,kinship,roh,ancestry,variant_annotations}/` | PGEN (genotypes), Parquet (tables) |
| predicted | `data/predicted/<modality>/<tissue>/` | Parquet |
| simulated | `data/simulated/<modality>/`, ground truth in `data/simulated/truth/` | BED, CSV, Parquet, YAML |
| reference | `resources/<name>/<version>/` + `manifest.yaml` | as published |
| phenotypes | `phenotypes/<id>/observations.parquet`, `definition.yaml` | Parquet, YAML |
| QC | `qc/<modality>/` | JSON, Parquet, MultiQC `*_mqc.tsv` |
| registry | `registry/efgpp.duckdb` + Parquet mirrors | DuckDB, Parquet |

**Storage policy** (`data/storage.py`): `copy` (shutil), `link` (symlink; falls back to reference
without copying if symlinks are not permitted), `reference` (register path, size, SHA-256, mtime)
and `auto`. Under `auto`, genotype formats are referenced, files under `copy_threshold_gb` are
copied, and larger files are referenced.

**Checksums** (`resources/checksums.py`): SHA-256 for single files; a manifest hash for filesets
and directories; a labelled `sha256-sampled` for files above `full_checksum_limit_gb`.

---

## 6. Registry and provenance

DuckDB tables (`data/registry.py`):
`participants`, `sample_aliases`, `events`, `biospecimens`, `assays`, `phenotype_definitions`,
`phenotype_observations`, `artifacts`, `artifact_parents`, `resources`, `software`, `tool_runs`,
`qc_runs`, `qc_metrics`, `modality_availability`, `snapshots`.

Artifact states: `REGISTERED → VALIDATED → STANDARDIZED → QC_PASS | QC_WARN | QC_FAIL → READY`,
with `SUPERSEDED` for replaced versions. Superseded artifacts are never deleted.

Each external command run through `data/provenance.py` (subprocess + psutil) writes:

- `.efgpp/logs/RUNxxxxxx.jsonl`: structured start/end events
- `.efgpp/logs/RUNxxxxxx.log`: combined stdout/stderr
- `.efgpp/logs/RUNxxxxxx.yaml`: the run record (tool, version, command, inputs, outputs, exit code, resources)
- a `tool_runs` row in the registry

---

## 7. Reference resources

| Resource | Provider class | Access | Default source |
|---|---|---|---|
| Genome FASTA + liftOver chains | `GenomeProvider` | httpx download, gunzip, own `.fai` builder | UCSC hg19 / hg38 |
| ClinVar | `ClinVarProvider` | download VCF + index | NCBI FTP `vcf_<build>/clinvar.vcf.gz` |
| gnomAD | `GnomADProvider` | user path/URL (terabyte scale) | — |
| AlphaMissense | `AlphaMissenseProvider` | download TSV | `storage.googleapis.com/dm_alphamissense` |
| AlphaGenome | `AlphaGenomeProvider` | atlas table or API | — |
| GWAS Catalog | `GWASCatalogProvider` | REST API v2 (metadata); FTP for summary statistics | EBI |
| PGS Catalog | `PGSCatalogProvider` | REST API + scoring-file download (no scoring here) | pgscatalog.org |
| Open Targets | `OpenTargetsProvider` | GraphQL for targeted queries; bulk files via FTP | platform.opentargets.org |
| GTEx | `GTExProvider` | API v2 / listed files | gtexportal.org |
| PredictDB | `PredictDBProvider` | tar download of models | Zenodo (GTEx v8 MASHR) |
| VEP cache | `install_vep_cache` | `vep_install` | Ensembl |

`resources/manager.py` installs each version into its own directory. It refuses to overwrite a
version referenced by a frozen snapshot, and `resources update --check` only reports newer
versions.

---

## 8. Workflow and execution

| Concern | Implementation |
|---|---|
| Plan (what can run, what cannot and why) | `data/plan.py`; written to `workflow/plan.json` |
| Local execution without Snakemake | `workflow/executor.py` (topological order; skips up-to-date steps using input fingerprints in `work/steps/<id>.done`) |
| Snakemake | `workflow/snakemake.py` generates `workflow/Snakefile` + `rules/*.smk`; each rule calls `efgpp data step <id>` |
| SLURM through Snakemake | `efgpp data prepare --executor slurm` with `workflow/profiles/slurm/config.yaml` |
| Explicit SLURM scripts | `workflow/slurm.py` → `hpc/*.sbatch` + `submit_all.sh` (with `afterok` dependencies) |
| Setup | `setup/manager.py`: platform detection → Pixi (fallback mamba → micromamba → conda) → environments → PLINK 2 binary where Bioconda is unavailable → tool tests → lock file |
| Doctor | `setup/doctor.py` |

---

## 9. Phenotype-independence rules (leakage prevention)

Performed in the Data layer: genotype QC, relatedness, sex check, ancestry, QC PCA, ROH, variant
annotation, liftover, PrediXcan prediction, profiling.

Never performed in the Data layer: association-based feature selection, P-value thresholds,
target-specific PRS, hyperparameter tuning, scaling fitted on all participants, missing-value
imputation. Covariates and omics are typed and profiled only. QC PCs are labelled as QC and must
be refit inside training folds.

---

## 10. Validation status of this release

| Area | How it was verified |
|---|---|
| Registry, validation, standardization, timeline, availability, report, snapshots, migration, CLI | Automated tests (59) on Windows |
| Genotype QC, kinship, PCA | End-to-end tests with the real PLINK 2 binary |
| VEP, OpenCRAVAT parsers; ClinVar/AlphaMissense lookups | Tests on small synthetic files; the external tools were not run |
| PrediXcan | Output handling tested with a mocked run |
| Snakefile and SLURM generation | Tested; not executed |
| ROH, ancestry, AlphaGenome API, Pixi environments on Linux/macOS | Implemented, not yet exercised |
