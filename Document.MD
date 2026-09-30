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
  efgpp CLI (Typer + Rich) ──► DataManager / plan ──► built-in executor (SLURM: export)
                                      │                                  │
             ┌────────────────────────┼──────────────────────────────────┤
             ▼                        ▼                                  ▼
     Adapters (Polars,          Genotype steps (PLINK 2)      Annotation / prediction
     Pandera, AnnData/Zarr)     QC · PCA · KING · ROH          (VEP, OpenCRAVAT, DuckDB
             │                        │                        lookups, MetaXcan)
             └──────────────┬─────────┴──────────────┬──────────────────┘
                            ▼                        ▼
               Registry (DuckDB + Parquet)    Provenance logs (logs/)      
                            │
                            ▼
      Availability ─► Report (Plotly + Jinja2 [+ MultiQC]) ─► DataSnapshot (YAML + Parquet)
```

EFGPP decides **what** must happen (the plan). Local runs use EFGPP's built-in executor (errors
printed directly in the terminal). For clusters, `efgpp export slurm` writes explicit sbatch
scripts (one per step group, with `afterok` dependencies). Snakemake is not used.

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
| **pyliftover** | Coordinate liftover GRCh37↔GRCh38 with local UCSC chain files; build check confirmation | `data/genotype/liftover.py`, `harmonization.py` |
| **hashlib** (stdlib) | SHA-256 artifact identity | `resources/checksums.py` |
| **blake3** *(extra `fast`)* | Optional fast local integrity hashing | `resources/checksums.py` |
| **sqlite3** (stdlib) | Reading PredictDB model metadata (`.db` files) | `data/predicted/metaxcan.py` |
| **phenopackets** *(extra `standards`)* | Validating exported GA4GH Phenopackets | `data/phenopackets.py` |
| **ga4gh.vrs** *(extra `standards`)* | Reserved for VRS identifiers (`vrs_id` column exists; not yet populated) | `data/schemas/artifact.py` |
| **alphagenome** *(extra `alphagenome`)* | AlphaGenome API mode | `data/annotation/alphagenome.py` |

Development: **pytest**, **pytest-cov**, **hypothesis** (property tests for checksums),
**ruff**, **mypy**, **pre-commit**.

### 2.2 External scientific tools

Each tool lives in its own conda environment under `software/envs/<tool>/` in the project (Pixi preferred; if
Pixi is unavailable or fails, the fastest conda-family tool found: mamba → micromamba → conda;
set with `execution.environment_manager` / `fallback_environment_manager`; definitions in
`workflow/envs/*.yaml`). `setup/tools.py` resolves each tool
in this order: `project.yaml: execution.tools` override → `software/bin` → its environment →
legacy `.efgpp/bin` → `$EFGPP_TOOLS_HOME/bin` (only when set) → `PATH`. Nothing is installed in
the home directory or in the `efgpp` environment.

When a tool cannot be obtained from conda, `setup/installers.py` (`efgpp setup tools`, also
called automatically by `efgpp setup data`) downloads it from its official source:

| Tool | Method | Source |
|---|---|---|
| PLINK 2, PLINK 1.9 | official binary (zip) → `bin/` | cog-genomics.org download pages (build chosen by OS/CPU/AVX2) |
| FlashPCA2 | official static binary → `bin/flashpca` | github.com/gabraham/flashpca release v2.0 |
| bcftools, tabix, bgzip | source build (`./configure && make && make install`, zlib only) → `opt/samtools`, linked in `bin/` | github.com/samtools htslib + bcftools latest releases |
| MultiQC, OpenCRAVAT, GWASLab | own conda environment; fallback uv environment → `envs/<name>`, linked in `bin/` | Bioconda / PyPI |
| MetaXcan / PrediXcan | GitHub source archive → `opt/MetaXcan` + Python 3.11 environment (`envs/metaxcan`); `bin/predixcan` wrapper | github.com/hakyimlab/MetaXcan |
| Ensembl VEP | container image → `containers/ensembl-vep.sif`; `bin/vep`, `bin/vep_install` wrappers | `ensemblorg/ensembl-vep` via Apptainer/Singularity or Docker |


| Tool | Environment | Used for | Called from |
|---|---|---|---|
| **PLINK 2** | `genetics` | Missingness, allele frequency, HWE, heterozygosity, LD pruning, KING kinship, sex check, QC-passed PGEN, PCA, format conversion, liftover (`--update-map`), ancestry projection (`--score`) | `data/genotype/qc.py`, `pca.py`, `conversion.py`, `harmonization.py`, `ancestry.py`, `roh.py` |
| **PLINK 1.9** | `genetics` | Runs of homozygosity (`--homozyg`, not available in PLINK 2) | `data/genotype/roh.py` |
| **FlashPCA2** | `genetics` | Optional high-scale PCA engine | `data/genotype/pca.py` |
| **bcftools / htslib (tabix, bgzip)** | `genetics` | Installed and version-checked for VCF/BCF work | `setup/manager.py` |
| **Ensembl VEP** | `annotation` | Molecular consequences (offline cache) | `data/annotation/vep.py` |
| **OpenCRAVAT** (`oc`) | `annotation` | Modular annotation aggregation | `data/annotation/opencravat.py` |
| **MetaXcan / PrediXcan** (`Predict.py`) | `metaxcan` | Genetically predicted expression from PredictDB models | `data/predicted/metaxcan.py` |
| **Data modules** (`src/efgpp/modules/*.py`) | EFGPP env | Column inference per kind of data (IDs, phenotypes and their type, covariate roles and categorical columns, GWAS columns); editable per project | `efgpp modules export` |
| **MultiQC** | `reporting` | Aggregates QC tables (`*_mqc.tsv`) next to the EFGPP report | `reporting/data_report.py` |
| **GWASLab** (`gwaslab-python`) | `gwaslab` | GWAS summary statistics: load, `basic_check`, `infer_build`, `liftover` to GRCh38 | `data/gwas.py` → `data/gwas_gwaslab_job.py` |
| **Pixi** / **mamba** / micromamba / conda | system | Creating the per-purpose environments from conda-forge + Bioconda | `setup/pixi.py`, `setup/micromamba.py` |
| **uv** | EFGPP env | Fast isolated Python environments for tools installed by `efgpp setup tools` | `setup/installers.py` |
| **Apptainer / Docker** | system | Optional container runtime for Linux-only tools | `setup/containers.py` |

### 2.3 Toolkits (`efgpp setup toolkit`, `setup/toolkits.py`, `setup/toolkit_installer.py`)

Declarative specs installed into `<root>/envs/<env>` (conda), `<root>/opt/<name>` (repositories,
downloads) and exposed as commands in `<root>/bin`; `<root>` is `<project>/software` (or
`$EFGPP_TOOLS_HOME` with `--shared`, only when set). `efgpp setup check`
verifies commands, R packages (`requireNamespace`), pip packages and repositories.

| Toolkit | Environment / source | Contents |
|---|---|---|
| `perl` | conda-forge | perl, perl-app-cpanminus |
| `r` | conda-forge (R ≥ 4.3 + compilers), then CRAN / GitHub / Bioconductor via Rscript | bigsnpr (LDpred-2, SCT, lassosum2), bigstatsr, data.table, magrittr, remotes, devtools, optparse, pROC, glmnet, ranger, caret, SuperLearner, dplyr, Matrix, Rcpp, RcppArmadillo, R.utils, BiocManager, genio, ggplot2, xgboost, susieR, GenomicRanges; permutations, RapidoPGS, sim1000G, PANPRSnext (CRAN); lassosum, penRegSum, CTSLEB, R2BGLiMS, EBPRS (GitHub) |
| `prs` | Bioconda | plink, plink2, gcta, gemma, vcftools, gsl, armadillo, compilers |
| `prs` | official downloads | PRSice-2 2.3.5, GCTB 2.5.5, LDAK 6.3, BOLT-LMM 2.5 |
| `prs` | GitHub archives + wrappers | PRScs, PRScsx, PRSbils, LDpred-funct, PolyFun (→ prs-python); LDSC, AnnoPred, PleioPred (→ prs-py27); DBSLMM.R (→ r); SDPR (bundled GSL/MKL on `LD_LIBRARY_PATH`), NPS (`make`), CTPR (unpacked), XP-BLUP, smtpred, MTG2 (source) |
| `prs-python` | conda-forge + PyPI | Python 3.10, numpy<2, scipy, pandas, h5py, scikit-learn, statsmodels, openjdk 11; ldpred, viprs, magenpy, pandas-plink, plinkio, pyliftover, pgenlib, hail |
| `prs-py27` | conda-forge + PyPI | Python 2.7, numpy, scipy, pandas, bitarray, h5py, scikit-learn, statsmodels; plinkio |
| `simulation` | conda-forge | Python 3.11, simuPOP, msprime, tskit, stdpopsim |

Conda installs use mamba > micromamba > conda (micromamba is downloaded when none exists) with
`--override-channels -c conda-forge -c bioconda`. All six environments were solved for
`linux-64` (dry run). Manual step: DBSLMM's `dbslmm` executable (Google Drive only).
Running these methods belongs to the Representation layer; the Data layer only installs them.

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
| 9 | Genome-build inference | Ensembl reference bases (GRCh37 + GRCh38) confirmed with **pyliftover**; contig lengths; coordinate bounds; FASTA REF concordance; optional rsID panel | `data/genotype/build.py`, `liftover.py` | `build_inference` metadata; validation fails when uncertain |
| 9b | Harmonization to GRCh38 | **pyliftover** (positions) + **PLINK 2** (`--update-map --sort-vars --make-pgen`) | `data/genotype/harmonization.py` (step `harmonize.<id>`) | `data/derived/genotype_qc/<id>/liftover/<id>_GRCh38.*`, liftover report; QC uses this copy |
| 9c | GWAS summary statistics | **GWASLab** (own environment) with the project's chain files | `data/gwas.py` (step `gwas.<id>`) | `resources/gwas/<id>/<id>.GRCh38.parquet` + GWASLab report |
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

- `logs/RUNxxxxxx.jsonl`: structured start/end events
- `logs/RUNxxxxxx.log`: combined stdout/stderr
- `logs/RUNxxxxxx.yaml`: the run record (tool, version, command, inputs, outputs, exit code, resources)
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
| Local execution (default) | `workflow/executor.py` built-in executor (topological order; skips up-to-date steps using input fingerprints in `work/steps/<id>.done`) |
| SLURM scripts | `workflow/slurm.py` → `hpc/*.sbatch` + `submit_all.sh`, each step runs `efgpp data step <id>` |
| Engine | built-in only (`execution.engine: auto`/`builtin`; Snakemake removed) |
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
| Registry, validation, standardization, timeline, availability, report, snapshots, migration, CLI | Automated tests (87) on Windows |
| Genotype QC, kinship, PCA | End-to-end tests with the real PLINK 2 binary |
| Build check + liftover | Real Ensembl/UCSC check: GRCh37-derived alleles → GRCh37 (100% vs 51%, pyliftover-lifted 100%), GRCh38 → GRCh38; harmonize step tested end to end with PLINK 2 and pyliftover |
| GWASLab step | Job wiring tested with a mocked GWASLab run; GWASLab environment solved for linux-64; not yet run on Linux |
| VEP, OpenCRAVAT parsers; ClinVar/AlphaMissense lookups | Tests on small synthetic files; the external tools were not run |
| PrediXcan | Output handling tested with a mocked run |
| SLURM script generation | Tested; not executed |
| Toolkits (`efgpp setup toolkit`) | Catalogue, R install/check script, repository wrappers unit-tested; all six conda environments solved for linux-64; not yet installed on Linux |
| Direct installers (`efgpp setup tools`) | Build selection unit-tested; binary and Python-environment installs run on Windows; Linux source build, Apptainer VEP and MetaXcan environment not yet run |
| ROH, ancestry, AlphaGenome API, Pixi environments on Linux/macOS | Implemented, not yet exercised |

---

## 11. Worked example: adding an existing PLINK cohort (migraine)

The legacy EFGPP folder `/data/ascher02/uqmmune1/ANNOVAR/migraine` contains, per trait:

| File | Content | EFGPP source |
|---|---|---|
| `migraine.bed` / `.bim` / `.fam` | genotypes (PLINK 1 binary fileset) | `genotype` (GENO001) |
| `migraine.height` | phenotype: header `IID FID Height`, `1` = control, `2` = case | `phenotype` (PH001) |
| `migraine.cov` | covariates: header `FID IID <covariates…>` | `covariates` (COV001) |
| `migraine_QC.*` | the old pipeline's QC output | not needed; EFGPP runs its own QC |
| `migraine*.txt`, `.gz`, `.PRSCS`, … | GWAS summary statistics and PRS results of the old pipeline | not participant data; see "Migration" below |

The phenotype values are also column 6 of the `.fam`, but `.fam` has no header line, so the
phenotype is read from `migraine.height` (same values, with a header).

### Commands

```bash
conda activate efgpp
cd /data/ascher02/uqmmune1/EFGPP/efgpp_projects/my_project     # the EFGPP project (efgpp init . done)

M=/data/ascher02/uqmmune1/ANNOVAR/migraine

# 0. Look at the headers once: the covariate names are needed in step 3
head -3 $M/migraine.cov $M/migraine.height

# 1. Genotype: give the PLINK prefix (no extension). It is referenced in place, never copied.
efgpp data add genotype --path $M/migraine --format bed --build auto --mode reference

# 2. Phenotype (binary, PLINK 1/2 coding is detected automatically: 2 = case, 1 = control)
efgpp phenotype add --name migraine --path $M/migraine.height \
    --id-column IID --value-column Height --type binary

# 3. Covariates: every column except FID/IID is taken automatically
#    (optional: --columns <names>, --categorical <names> for codes such as sex)
efgpp data add covariates --path $M/migraine.cov --id-column IID

# 4. Check, prepare, look at the result, freeze
efgpp data validate
efgpp data prepare
efgpp data availability
efgpp data report                      # reports/data/index.html
efgpp data freeze --name migraine_v1
```

A source added wrongly is removed with `efgpp data remove COV001` and added again.

### What EFGPP does with these files

| Step | Result |
|---|---|
| `data add genotype` | `data.yaml` gets GENO001 pointing at `$M/migraine`; the three files are checksummed in place (storage mode `reference`) |
| `phenotype add` | PH001 `migraine` (binary); the small table is copied to `data/observed/phenotype/PH001/` |
| `data add covariates` | COV001; copied to `data/observed/covariates/COV001/`; no imputation or scaling |
| `data validate` | duplicate IDs, invalid phenotype codes, missing covariate columns, genome build, and how many phenotype participants have genotype data |
| `data prepare` | participant registry, standardized phenotype/covariates, PLINK 2 QC (`qc/genotype/GENO001/`, QC-passed PGEN in `data/derived/genotype_qc/GENO001/`), QC PCA, kinship, availability, report |

### Genome build

`--build auto` infers the build from the variant positions (and the reference FASTA once
`efgpp resources install genome` has run). For the migraine data, positions beyond the GRCh38
ends of chromosomes 1, 2 and 4 identify GRCh37. If the evidence is not conclusive,
`efgpp data validate` stops with "genome build could not be determined confidently". Then set it explicitly — the
UK Biobank genotypes of the original EFGPP study are GRCh37 (hg19) — either by re-adding with
`--build GRCh37` or by editing `data.yaml`:

```yaml
observed:
  genotype:
    - id: GENO001
      genome_build: GRCh37
```

### IDs

The `.fam`, `.cov` and `.height` files are matched on `IID` (never by row order). If IIDs are
only unique together with the FID, use `--sample-id-mode fid_iid` for the genotype and give the
phenotype/covariate tables an ID column built as `FID_IID`.

### Several traits

Add each trait's phenotype as another `efgpp phenotype add` (PH002, PH003, …). Re-add a genotype
only if it is a different fileset. `efgpp data availability` shows which participants have which
phenotype.

### Migration (the whole legacy folder at once)

To register everything in the legacy folder in one go (genotype filesets, `.cov`, `.height`,
GWAS summary statistics as reference resources, old PCA/PRS outputs for provenance only):

```bash
efgpp data migrate --source /data/ascher02/uqmmune1/ANNOVAR/migraine --profile legacy-efgpp --dry-run
# review migration_plan.yaml (it also lists migraine_QC.*, migraine_5/, … — remove what you do not want)
efgpp data migrate --source /data/ascher02/uqmmune1/ANNOVAR/migraine --profile legacy-efgpp --apply
```

The source folder is never modified.

---

## 12. Genome build policy: everything in GRCh38

`project.yaml: defaults.target_build: GRCh38`. Every dataset is brought to this build; the
originals are never modified and every lifted copy is a new artifact with its own lineage.

| Data | How the build is determined | How it becomes GRCh38 |
|---|---|---|
| Genotypes (PLINK/VCF/BGEN) | 1. **reference bases + pyliftover**: ~80 SNVs are compared with the GRCh37 and GRCh38 reference (Ensembl REST, batched); the build whose reference matches the alleles (~100% vs ~50%) wins, and the positions lifted with pyliftover must match the other build's reference too. 2. contig lengths (VCF), coordinate bounds (positions past chromosome ends), local FASTA, declared value. Contradictions stop validation. | step `harmonize.<id>`: pyliftover converts positions (UCSC chain), PLINK 2 writes a new PGEN; unmapped, cross-chromosome and reverse-strand variants are excluded and listed in the liftover report. QC, PCA, annotation and prediction then use the GRCh38 copy. |
| Reference panels (`ld_reference`) | declared `genome_build`, else the same evidence | lifted once into `resources/ld_reference/<name>_GRCh38` |
| GWAS summary statistics | GWASLab `infer_build` (HapMap3 positions shipped with GWASLab) or the declared build | GWASLab `liftover(to_build="38", chain_path=resources/liftover/…)` |
| Downloaded resources (ClinVar, AlphaMissense, VEP cache, genome, PGS scoring files) | — | installed for GRCh38 only; requesting GRCh37 is refused |
| Predicted expression (PredictDB GTEx v8) | models are GRCh38 | uses the GRCh38 genotype |

Chain files (`hg19ToHg38.over.chain.gz`, `hg38ToHg19.over.chain.gz`) are downloaded once from
UCSC into `resources/liftover/` and passed explicitly to pyliftover and GWASLab, so no tool
downloads anything into the home directory. `EFGPP_OFFLINE=1` (or
`defaults.online_build_check: false`) skips the Ensembl check and relies on local evidence.

### GWAS summary statistics with GWASLab

```bash
efgpp setup tools gwaslab                              # software/envs/gwaslab (never the efgpp env)
efgpp data add gwas --path <file> --trait <text> [--build auto] [--fmt auto] \
    [--col chrom=CHR --col pos=BP --col snpid=SNP --col ea=A1 --col nea=A2 --col p=P ...] [--n N]
efgpp data validate                                    # file readable, mapped columns present
efgpp data prepare                                     # step gwas.<id>
```

`--col` uses GWASLab keyword names (`snpid rsid chrom pos ea nea eaf beta se p OR z n info`).
The job (`data/gwas_gwaslab_job.py`, run by the environment's Python) does
`Sumstats(...)` → `basic_check(remove=True, remove_dup=True)` → `infer_build()` →
`liftover(...)` → Parquet, and writes a JSON report (rows loaded / after checks / after
liftover, detected build, GWASLab version) to `resources/gwas/<id>/` and `qc/gwas/`. The result
is a REFERENCE artifact (`gwas_sumstats`) in GRCh38; choosing variants or thresholds for a
target phenotype remains a Representation-layer decision. `efgpp data migrate` adds legacy GWAS
files as GWAS sources automatically.

---

## 13. Data modules and automatic column inference

Each kind of data has one readable file in `src/efgpp/modules/`; `efgpp modules export <name>`
copies it to `<project>/modules/<name>.py`, and that copy is used for the project.

| Command | Inferred (module) | Stored in data.yaml |
|---|---|---|
| `efgpp data add genotype --path <prefix>` | format (pgen/bed/bgen/vcf/bcf), sample IDs IID vs FID_IID (`genotype.py`) | `format`, `sample_id.mode`; build auto → detected and lifted to GRCh38 |
| `efgpp phenotype add --path <file>` | ID column; every phenotype column (ID, FID and covariate-looking columns skipped); type binary / continuous / multiclass; binary coding (PLINK 1/2 → 2 = case, 0/1, case/control words); name = file name for one-column files (`phenotype.py`) | one `PHxxx` per phenotype column with `value_column` = original name, `type`, `case_values`/`control_values` or `levels` |
| `efgpp data add covariates --path <file>` | ID column; variables = all other columns except FID; role (sex, age, pc, batch, other) from name and values; categorical = text, sex, batch, or ≤10 integer codes (`covariates.py`) | `variables`, `categorical`, `roles` |
| `efgpp data add gwas --path <file>` | GWASLab keyword per column (CHR/chrom, BP/pos, SNP/snpid, A1/ea, A2/nea, BETA, OR, SE, P, N, INFO, EAF, …) and whether the essentials are present (`gwas.py`) | `columns` |

Every explicit option (`--id-column`, `--value-column`, `--type`, `--columns`, `--categorical`,
`--col keyword=column`) overrides the inference.

**Original column names are kept** for feature engineering:

| Output | Columns |
|---|---|
| `phenotypes/<id>/phenotype.parquet` | `participant_id` + the phenotype column with its original name and values |
| `phenotypes/<id>/observations.parquet` | EFGPP's coded form (value_numeric 0/1, value_text) used internally |
| `data/observed/covariates/<id>/covariates.parquet` | `participant_id` + every covariate with its original name |
| `resources/gwas/<id>/<id>.GRCh38.parquet` | the original GWAS column names (values checked and lifted by GWASLab) |
| `resources/gwas/<id>/<id>.GRCh38.gwaslab.parquet` | GWASLab standard names (CHR, POS, EA, NEA, …, STATUS) |

### GWAS with study details

```bash
efgpp data add gwas --path <file> --trait <GWAS trait> [--phenotype <project phenotype> ...] \
    [--ancestry European] [--n-cases N] [--n-controls N] [--n N] [--study <name/accession>] \
    [--build auto] [--col keyword=column] [--no-run]
efgpp data gwas list [--phenotype <name>]        # which GWAS exist for a phenotype
efgpp data gwas run <GWAS id> [--force]          # (re)process with GWASLab
```

`--phenotype` links a GWAS to the project phenotype(s) it is meant for; a GWAS of one trait can
serve another phenotype (list it with several `--phenotype`). By default the GWAS runs through
GWASLab immediately (needs `efgpp setup tools gwaslab`); unchanged files are not reprocessed.
Case/control counts given as numbers become `N_CASE`/`N_CONTROL` columns; ancestry, counts,
study and target phenotypes are recorded with the artifact and in the data report.
