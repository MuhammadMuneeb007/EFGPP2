# EFGPP — Data layer

**Exploratory Framework for Genotype–Phenotype Prediction**, layer 1 of 4:

```text
DATA  →  REPRESENTATION  →  PREDICTION  →  INTERPRETATION
```

The Data module turns an arbitrary cohort (genotypes, any number of phenotypes, covariates,
omics, timelines, biospecimens) into a validated, provenance-tracked, immutable
**DataSnapshot**. The Representation layer consumes that snapshot by name; it never searches
folders.

It is **phenotype-agnostic**: a phenotype is a configuration object (`binary`, `continuous`,
`multiclass`, `ordinal`; survival/longitudinal types are declared for later). No disease or
trait is hard-coded, and a test enforces that. The one trait-like token is the legacy
`.height` file extension, which the migration profile reads as a generic phenotype file.

## Install

```bash
uv sync                      # Python ≥ 3.11; creates .venv with efgpp + dev tools
uv run efgpp --version
```

Optional extras: `efgpp[omics]` (AnnData/Zarr/MuData; otherwise omics are stored as Parquet),
`efgpp[reporting]` (Kaleido/UpSetPlot), `efgpp[standards]` (Phenopackets validation, GA4GH VRS),
`efgpp[alphagenome]`.

## Quick start

```bash
mkdir my_project && cd my_project
efgpp init .
efgpp setup data                         # PLINK 2, bcftools, Snakemake, … in .efgpp/envs/<purpose>

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

## Concepts

| Concept | Where |
|---|---|
| Five origins: `observed`, `derived`, `predicted`, `simulated`, `reference` | `efgpp.constants.Origin` |
| Artifacts (id, type, modality, origin, status, path, checksum, parents, tool, versions, …) | DuckDB `registry/efgpp.duckdb` |
| States `REGISTERED → VALIDATED → STANDARDIZED → QC_* → READY`, older versions `SUPERSEDED` (never deleted) | `data/artifacts.py` |
| Canonical `participant_id`; native IDs resolved through alias files, never by row position | `data/aliases.py` |
| participant → event → biospecimen → assay | `data/timeline.py`, `data/biospecimens.py` |
| Storage policy `copy` / `link` / `reference` / `auto` (genotype cohorts are referenced in place) | `data/storage.py` |
| Every external command: run record + `.efgpp/logs/RUN*.jsonl` + `.log` + `.yaml` | `data/provenance.py` |
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
| Workflow | Snakemake (generated `workflow/Snakefile`); built-in executor when Snakemake is absent | built-in executor tested; Snakefile/SLURM generation tested, not executed |

Native Windows runs the core and PLINK 2. Bioconda has no Windows builds, so VEP, bcftools and
MetaXcan need WSL2 or Docker (`efgpp doctor` says so).

## Development

```bash
uv run pytest            # unit + acceptance tests A–J (downloads the PLINK 2 binary once; set EFGPP_TEST_OFFLINE=1 to skip)
uv run ruff check src tests
uv run mypy
```
