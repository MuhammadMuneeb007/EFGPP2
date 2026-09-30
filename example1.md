# Example 1 — add the migraine PLINK cohort

Files in `/data/ascher02/uqmmune1/ANNOVAR/migraine`:

| File | Content | Becomes |
|---|---|---|
| `migraine.bed` / `.bim` / `.fam` | genotypes (733 samples, 619,653 variants, GRCh37) | genotype `GENO001` |
| `migraine.height` | phenotype: header `IID FID Height`, `1` = control, `2` = case | phenotype `PH001` (name `migraine`) |
| `migraine.cov` | covariates: header `FID IID <covariates…>` | covariates `COV001` |
| `migraine.gz` | GWAS summary statistics: `CHR BP SNP A1 A2 N SE P OR INFO MAF` | GWAS `GWAS001` (GWASLab) |

EFGPP infers everything itself (the rules are in `src/efgpp/modules/*.py`, editable per
project with `efgpp modules export`):

- **genotype**: format (bed), sample IDs (IID), genome build (GRCh37, checked against the
  reference with pyliftover) → lifted to GRCh38 before QC.
- **phenotype**: ID column `IID`, phenotype column `Height`, type binary, PLINK coding
  (2 = case, 1 = control); named `migraine` after the file.
- **covariates**: ID column `IID`, every other column except `FID`; sex / age / PCs / batch
  recognised by name and values; categorical vs numeric decided automatically.
- **GWAS**: CHR→chrom, BP→pos, SNP→snpid, A1→ea, A2→nea, N, SE, P, OR, INFO recognised;
  GWASLab checks it, infers the build and lifts it to GRCh38.

Original column names are kept in every output file (for feature engineering later).
`migraine_QC.*` is not needed (EFGPP runs its own QC); the other `migraine*.txt`, `.PRSCS`, …
are old PRS results.

```bash
# ---------------------------------------------------------------------------
# Setup: go to your EFGPP project
# ---------------------------------------------------------------------------
conda activate efgpp
cd /data/ascher02/uqmmune1/EFGPP/EFGPP2/my_project
[ -f project.yaml ] || efgpp init .
efgpp setup tools gwaslab                  # GWASLab in its own conda env (once per project)

# ---------------------------------------------------------------------------
# 1. Genotype (referenced in place, never copied)
# ---------------------------------------------------------------------------
efgpp data add genotype --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine

# ---------------------------------------------------------------------------
# 2. Phenotype (ID column, phenotype column, type and coding inferred)
# ---------------------------------------------------------------------------
efgpp phenotype add --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.height

# ---------------------------------------------------------------------------
# 3. Covariates (ID column, covariates, sex/age/PC/batch roles, categorical inferred)
# ---------------------------------------------------------------------------
efgpp data add covariates --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.cov

# ---------------------------------------------------------------------------
# 4. GWAS for migraine: columns inferred; GWASLab runs now and stores the result
#    (add what you know about the study on the same command, for example
#     --ancestry European --n-cases <number> --n-controls <number> --study <name>)
# ---------------------------------------------------------------------------
efgpp data add gwas --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.gz \
    --trait migraine --phenotype migraine

# A GWAS of another trait used for migraine (e.g. the legacy depression GWAS):
# efgpp data add gwas --path <depression GWAS file> --trait depression --phenotype migraine

efgpp data gwas list --phenotype migraine  # the GWAS available for migraine

# ---------------------------------------------------------------------------
# 5. Check, run, inspect, freeze (run one at a time; stop if validate shows ✗)
#    validate shows: "genome build GRCh37 (... reference_bases_pyliftover)"
#                    "will be lifted GRCh37 -> GRCh38 with pyliftover"
#    prepare runs:   harmonize.GENO001 (liftover) -> genotype_qc -> pca -> availability -> report
# ---------------------------------------------------------------------------
efgpp data validate                    # every problem, with counts
efgpp data prepare                     # QC, PCA, kinship, availability, report
efgpp data availability                # who has genotype + phenotype + covariates
efgpp data report                      # reports/data/index.html
efgpp data freeze --name migraine_v1   # immutable snapshot for the next layer

# ---------------------------------------------------------------------------
# Fixing a source that was added wrongly
# ---------------------------------------------------------------------------
# efgpp data remove COV001             # also works for GWAS ids, e.g. GWAS001
# efgpp modules export covariates      # edit modules/covariates.py for this project, then re-add

# ---------------------------------------------------------------------------
# Alternative: register the whole legacy folder at once (source never modified)
# ---------------------------------------------------------------------------
# efgpp data migrate --source /data/ascher02/uqmmune1/ANNOVAR/migraine --profile legacy-efgpp --dry-run
# efgpp data migrate --source /data/ascher02/uqmmune1/ANNOVAR/migraine --profile legacy-efgpp --apply
```

Where the results are

| What | File |
|---|---|
| phenotype, original column name | `phenotypes/PH001/phenotype.parquet` (`participant_id`, `Height`) |
| covariates, original column names | `data/observed/covariates/COV001/covariates.parquet` |
| genotype lifted to GRCh38 | `data/derived/genotype_qc/GENO001/liftover/GENO001_GRCh38.*` + `…_liftover_report.parquet` |
| QC-passed genotype (GRCh38) | `data/derived/genotype_qc/GENO001/GENO001_qc.*` |
| GWAS, original column names (GRCh38) | `resources/gwas/GWAS001/GWAS001.GRCh38.parquet` |
| GWAS, GWASLab names (GRCh38) | `resources/gwas/GWAS001/GWAS001.GRCh38.gwaslab.parquet` + `GWAS001.gwaslab_report.json` |

Notes

- Samples are matched on `IID` across the `.fam`, `.cov` and `.height` files, never by row order.
- The 2 variants at position 0 are unplaced array probes: kept in the original, left out of the
  GRCh38 copy (they cannot be lifted); see the liftover report.
- `efgpp data prepare` runs everything locally and prints errors directly. For a cluster,
  `efgpp export slurm` writes sbatch scripts.
- More detail: `Document.MD`, sections 11–13.
