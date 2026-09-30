# EFGPP — Data layer

**Exploratory Framework for Genotype–Phenotype Prediction**, layer 1 of 4:

```text
DATA  →  REPRESENTATION  →  PREDICTION  →  INTERPRETATION
```

The Data module turns an arbitrary cohort (genotypes, any number of phenotypes, covariates,
omics, timelines, biospecimens) into a validated, provenance-tracked, immutable
**DataSnapshot**. The Representation layer consumes that snapshot by name; it never searches
folders.

Everything is kept in **GRCh38**: genotypes, reference panels and GWAS summary statistics in
another build are lifted automatically (pyliftover; GWASLab for GWAS files) into new artifacts —
see `Document.MD` §12.

It is **phenotype-agnostic**: a phenotype is a configuration object (`binary`, `continuous`,
`multiclass`, `ordinal`; survival/longitudinal types are declared for later). No disease or
trait is hard-coded, and a test enforces that. The one trait-like token is the legacy
`.height` file extension, which the migration profile reads as a generic phenotype file.

## Install

EFGPP uses **two kinds of environment**. You only ever activate the first one.

| Environment | Contains | Who manages it |
|---|---|---|
| **`efgpp`** (conda or a uv venv) | the EFGPP Python application: Polars, DuckDB, Pandera, Typer… | you: create it once, `conda activate efgpp` per session |
| **`<project>/software/envs/<tool>`** — one per tool | scientific tools: PLINK 2, bcftools, VEP (+ Perl), OpenCRAVAT, MetaXcan, MultiQC, GWASLab, R, PRS tools | EFGPP: created by `efgpp setup data` and switched automatically per step |

Scientific tools conflict with each other (e.g. VEP needs Perl, MetaXcan its own Python), so
each group gets its own environment. You never `conda activate` these: when a step runs,
EFGPP executes the tool inside its environment for that command only.

`efgpp setup data` builds these with **Pixi** by default (parallel, lock-file based). If Pixi
cannot be used or fails, it falls back to **mamba**, then micromamba, then conda, whichever is
installed. To use mamba directly, set it in `project.yaml`:

```yaml
execution:
  environment_manager: mamba             # pixi | mamba | micromamba | conda | system
```

### Option A — mamba (fastest; recommended on HPC clusters)

```bash
git clone https://github.com/MuhammadMuneeb007/EFGPP2.git
cd EFGPP2
mamba env create -f environment.yml      # creates the `efgpp` env; installs EFGPP in editable mode
conda activate efgpp                     # or: mamba activate efgpp
efgpp --version                          # -> efgpp 0.1.0
```

`environment.yml` takes every dependency pre-built from conda-forge in a single solve with
parallel downloads. pip only links the EFGPP source (it downloads nothing else). Measured on one
Windows workstation:

| Command | Empty package cache | Warm cache |
|---|---|---|
| `mamba env create -f environment.yml` | 111 s | 56 s |
| `conda env create -f environment.yml` | — | 121 s |
| previous version (dependencies via pip) | 225 s | — |

No mamba? `conda install -n base -c conda-forge mamba`, or install
[Miniforge](https://github.com/conda-forge/miniforge), which ships it. `micromamba create -f
environment.yml` works too, and `conda env create -f environment.yml` also works (slower).

Every new shell then only needs `conda activate efgpp`. Editable mode means `git pull` updates
EFGPP without reinstalling; re-run `mamba env update -f environment.yml` only when
dependencies change. To install into an existing conda environment instead:

```bash
conda activate <your-env>                # needs Python >= 3.11
pip install -e "/path/to/EFGPP2[omics,test]"
```

Remove the environment with `conda env remove -n efgpp`.

### Option B — uv

```bash
cd EFGPP2
uv sync                                  # creates EFGPP2/.venv with efgpp + dev tools
source .venv/bin/activate                # Windows: .venv\Scripts\activate
efgpp --version
```

Without activating, prefix commands with `uv run` (this only works inside the EFGPP2
folder), or install a user-wide command with `uv tool install --editable /path/to/EFGPP2`.

### Extras

`omics` (AnnData/Zarr/MuData; without it omics are stored as Parquet), `test` (pytest,
hypothesis), `reporting` (Kaleido/UpSetPlot), `standards` (Phenopackets validation, GA4GH VRS),
`alphagenome`. Example: `pip install -e ".[omics,test,reporting]"`.

### Toolkits: Perl, R, PRS tools, simulation

Beyond the data-layer tools, EFGPP installs complete toolkits for the software used around it
(the requirements are taken from PRSTools):

```bash
efgpp setup toolkit --list          # contents of every toolkit
efgpp setup toolkit all             # perl, r, prs-python, prs-py27, prs, simulation -> ./software
efgpp setup check                   # ✓/✗ for every tool, R package and repository
```

| Toolkit | Contents | How |
|---|---|---|
| `perl` | Perl, cpanminus | conda-forge |
| `r` | R ≥ 4.3; bigsnpr (**LDpred-2**, SCT, lassosum2), bigstatsr, data.table, glmnet, caret, SuperLearner, susieR, GenomicRanges, genio, …; lassosum, PANPRSnext, CTSLEB, RapidoPGS, EBPRS, R2BGLiMS (JAMPred), penRegSum (tlpSum), sim1000G, permutations | conda-forge, then CRAN / GitHub / Bioconductor inside the environment (compilers included) |
| `prs` | PLINK 1.9/2, GCTA, GEMMA, vcftools (conda); **GCTB 2.5.5, PRSice-2, LDAK 6.3, BOLT-LMM 2.5** (official downloads); PRScs, PRScsx, PRSbils, LDpred-funct, PolyFun, SDPR, DBSLMM, CTPR, NPS, XP-BLUP, smtpred, LDSC, AnnoPred, PleioPred, MTG2 (GitHub, with wrapper commands in `bin/`) | Bioconda + official downloads + GitHub archives |
| `prs-python` | Python 3.10: LDpred, VIPRS, magenpy, pandas-plink, pgenlib, **Hail** (+ Java 11) | conda-forge + PyPI |
| `prs-py27` | Python 2.7 for LDSC, AnnoPred, PleioPred | conda-forge |
| `simulation` | **simuPOP**, msprime, tskit, stdpopsim (R: sim1000G in `r`) | conda-forge |

Every environment was solved for `linux-64` against conda-forge + Bioconda only (EFGPP passes
`--override-channels`, so Anaconda's commercial `defaults` channel is never used). If no
mamba/micromamba/conda is installed, EFGPP downloads micromamba itself. Bioconda's `bolt-lmm`
is uninstallable (it needs an `nlopt` release that does not exist), so BOLT-LMM comes from the
official tarball. DBSLMM's `dbslmm` executable is only on Google Drive and must be downloaded by
hand into `<install root>/opt/DBSLMM/software/dbslmm`.

### Check the installation

```bash
efgpp doctor                             # what is installed / missing
pytest                                   # from the EFGPP2 folder; 87 tests
```

### If a tool cannot be installed with mamba/conda

`efgpp setup data` first builds the tool environments (Pixi, else mamba/micromamba/conda). Any
tool that is **still missing** is then downloaded from its **official source** and placed in the
project, so a missing Bioconda package never blocks you. You can also do this directly:

```bash
cd ~/efgpp_projects/my_project
efgpp setup tools                        # every missing tool
efgpp setup tools plink2 gwaslab multiqc     # only these
efgpp setup tools --force plink2         # reinstall
```

| Tool | How EFGPP gets it (Linux) | Needs |
|---|---|---|
| `plink2` | official binary from cog-genomics.org (AVX2 / AMD / x86_64 / ARM build chosen automatically) | internet |
| `plink` (1.9) | official binary from cog-genomics.org | internet |
| `flashpca2` | official static binary (GitHub release v2.0, x86_64) | internet |
| `bcftools`, `tabix`, `bgzip` | built from the official htslib + bcftools release tarballs | `gcc`/`cc`, `make`, zlib headers (e.g. `module load gcc`) |
| `multiqc`, `oc` (OpenCRAVAT), `gwaslab` | own conda environment each (fallback: uv environment from PyPI) | internet |
| `predixcan` (MetaXcan) | MetaXcan source from GitHub + its own Python 3.11 environment (uv downloads Python if needed) | internet |
| `vep` | official `ensemblorg/ensembl-vep` image with `vep` / `vep_install` wrapper scripts | Apptainer/Singularity (or Docker) |

Everything lands in the project directory, in plain sight: `software/envs/<tool>` (one
conda environment per tool — never the `efgpp` environment), commands in `software/bin/`,
sources and downloads in `software/opt/`, logs in `software/logs/`; reference data goes to
`resources/`. Nothing is written to your home directory. (`--shared` installs into
`$EFGPP_TOOLS_HOME` instead, only if you set it.) Each tool is installed with conda first
(conda-forge + Bioconda); the official download is only a fallback. EFGPP finds them automatically. To call them yourself from the shell:

```bash
eval "$(efgpp setup path)"               # puts ./software/bin on PATH
plink2 --version && multiqc --version
```

Tools you already have (e.g. `plink`, `bcftools` in your `efgpp` mamba env) are detected and
skipped.

## Quick start

Create each project **outside the EFGPP2 repository**. A project accumulates data, registry
files, logs and reports that must never be committed to the code repository.

```bash
conda activate efgpp
mkdir -p ~/efgpp_projects/my_project && cd ~/efgpp_projects/my_project
efgpp init .
efgpp setup data                         # PLINK 2, bcftools, MultiQC, GWASLab, … in software/envs/<tool>

efgpp data add genotype --path /data/genotype/cohort --format pgen --build auto --mode reference
efgpp phenotype add --name trait_a --path phenotypes.csv --id-column IID --value-column TRAIT_A --type binary
efgpp phenotype add --name trait_b --path phenotypes.csv --id-column IID --value-column TRAIT_B --type continuous
efgpp data add covariates --path covariates.csv --id-column IID --columns age sex bmi
efgpp data add expression --path expression.parquet --id-column participant_id \
    --event-column visit_id --time-column collection_date --tissue whole_blood

efgpp data inspect        # what exists
efgpp data validate       # every problem, with counts (✓ / ✗ / !)
efgpp data plan           # what can run now, what cannot, and why
efgpp data prepare        # register → validate → standardize → QC → derive → annotate → predict → availability → report
efgpp data availability   # participants per modality, grouped by origin, plus intersections
efgpp data report         # reports/data/index.html (23 sections)
efgpp data freeze --name cohort_data_v1
```

The Representation layer then does:

```python
from efgpp.project import Project
from efgpp.data.snapshots import DataSnapshot

snap = DataSnapshot.load(Project.load("my_project"), "cohort_data_v1")
sel = snap.select("trait_a", require=["genotype_qc"], include=["expression", "covariates", "qc_pca"])
sel.participants            # participant IDs with trait_a AND QC-passed genotype
sel.artifacts["genotype_qc"]  # exact registered artifacts (path, checksum, lineage, versions)
```

## Data modules (one editable file per kind of data)

How EFGPP reads each kind of data lives in its own small file under `src/efgpp/modules/`:

| File | Decides |
|---|---|
| `columns.py` | participant ID column names (IID, participant_id, eid, …), missing-value tokens |
| `genotype.py` | genotype format and sample IDs (IID or FID_IID) |
| `phenotype.py` | which columns are phenotypes; binary / continuous / multiclass; case/control coding |
| `covariates.py` | covariate roles (sex, age, PCs, batch) from names and values; categorical vs numeric |
| `gwas.py` | which GWAS column is CHR, BP, A1, A2, BETA/OR, SE, P, N, … (for GWASLab) |

`efgpp modules export covariates` copies a module into `<project>/modules/`; edit that copy and
EFGPP uses it for this project (`efgpp modules list` shows which file is active). Column names
are never changed: standardized tables keep the original names for feature engineering.

## Concepts

| Concept | Where |
|---|---|
| Five origins: `observed`, `derived`, `predicted`, `simulated`, `reference` | `efgpp.constants.Origin` |
| Artifacts (id, type, modality, origin, status, path, checksum, parents, tool, versions, …) | DuckDB `registry/efgpp.duckdb` |
| States `REGISTERED → VALIDATED → STANDARDIZED → QC_* → READY`, older versions `SUPERSEDED` (never deleted) | `data/artifacts.py` |
| Canonical `participant_id`; native IDs resolved through alias files, never by row position | `data/aliases.py` |
| participant → event → biospecimen → assay | `data/timeline.py`, `data/biospecimens.py` |
| Storage policy `copy` / `link` / `reference` / `auto` (genotype cohorts are referenced in place) | `data/storage.py` |
| Every external command: run record + `logs/RUN*.jsonl` + `.log` + `.yaml` | `data/provenance.py` |
| `efgpp.lock.yaml`: software, environments, resources, config checksums (generated) | `setup/lock.py` |

**Data vs Representation boundary.** The Data layer only does phenotype-independent work:
genotype QC (missingness, MAF, HWE on all samples, heterozygosity, KING kinship, sex check),
QC PCA, ROH, ancestry projection, annotation, PrediXcan prediction. It never selects features by
association, fits scalers or imputes values across participants. QC PCs are marked as such;
modelling PCs belong inside training folds.

## What runs where

| Component | Implementation | Status in this release |
|---|---|---|
| Registry, validation, standardization, availability, timeline, report, snapshots, migration | pure Python (Polars, DuckDB, Pandera, Pydantic) | tested on Windows |
| Genotype QC, QC PCA, kinship, conversion, liftover | PLINK 2 (`efgpp setup data` installs the official binary / Bioconda build) | QC + PCA tested end to end with the real PLINK 2 |
| ROH | PLINK 1.9 (`--homozyg` is not in PLINK 2) | implemented; not exercised in tests |
| Ancestry | PLINK 2 projection onto `resources.yaml: ld_reference` | implemented; not exercised in tests |
| VEP, OpenCRAVAT | wrapped (never re-implemented); Linux/macOS or WSL2/Docker | command building and parsers tested; the tools themselves not run here |
| ClinVar, gnomAD, dbSNP, AlphaMissense | DuckDB lookups against the official files | tested on small synthetic files |
| AlphaGenome | atlas table lookup; optional API mode | API mode is written against the published client and untested |
| PrediXcan | MetaXcan `Predict.py` + PredictDB models | output handling tested with a mocked run |
| Workflow | EFGPP's built-in executor; `efgpp export slurm` writes sbatch scripts for clusters | built-in executor tested; SLURM script generation tested, not executed |

Native Windows runs the core and PLINK 2. Bioconda has no Windows builds, so VEP, bcftools and
MetaXcan need WSL2 or Docker (`efgpp doctor` says so).

## Running on an HPC cluster

- **Run `efgpp setup data` / `efgpp setup tools` on a login node.** They download packages,
  binaries and images, and compute nodes often have no internet access.
- **Keep projects on a shared filesystem** (home or project storage, not node-local `/tmp`), so
  compute jobs see the same `software/envs` and registry.
- **Submit work** with explicit SLURM scripts: `efgpp export slurm` writes `hpc/*.sbatch` and
  `hpc/submit_all.sh` (set `execution.hpc.partition` / `account` in `project.yaml`). They call
  EFGPP through the absolute path of the `efgpp` environment's Python, so jobs need no
  `conda activate`; regenerate them if you recreate or move that environment.
- **Reuse tools you already have in conda environments** instead of letting EFGPP build new
  ones. Point to them in `project.yaml`; EFGPP runs those binaries directly:

  ```yaml
  execution:
    environment_manager: system          # do not create environments
    tools:
      plink2: /home/<user>/miniconda3/envs/genetics/bin/plink2
      bcftools: /home/<user>/miniconda3/envs/genetics/bin/bcftools
      vep: /home/<user>/miniconda3/envs/vep/bin/vep
  ```

  Then check with `efgpp doctor`.

> Automatic environment creation with Pixi/mamba on Linux is implemented but has not yet
> been exercised on a cluster; if `efgpp setup data` fails, the `tools:` override above does not
> depend on it.

## Development

```bash
uv run pytest            # unit + acceptance tests A–J (downloads the PLINK 2 binary once; set EFGPP_TEST_OFFLINE=1 to skip)
uv run ruff check src tests
uv run mypy
```
