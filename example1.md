# Example 1 — the migraine PLINK cohort, from installation to snapshot

Files in `/data/ascher02/uqmmune1/ANNOVAR/migraine`:

| File | Content | Becomes |
|---|---|---|
| `migraine.bed` / `.bim` / `.fam` | genotypes (733 samples, 619,653 variants, GRCh37) | genotype `GENO001` |
| `migraine.height` | phenotype: header `IID FID Height`, `1` = control, `2` = case | phenotype `PH001` (name `migraine`) |
| `migraine.cov` | covariates: header `FID IID <covariates…>` | covariates `COV001` |
| `migraine.gz` | GWAS summary statistics: `CHR BP SNP A1 A2 N SE P OR INFO MAF` | GWAS `GWAS001` (GWASLab) |

Run everything on the **login node** (downloads need internet); for long steps use an
interactive job or `efgpp export slurm`. Everything is installed inside the project folder:
tools in `software/`, reference data in `resources/`, logs in `logs/` — nothing in your home
directory, and every tool has its own conda environment (never the `efgpp` one).

```bash
# ===========================================================================
# 0. EFGPP itself (once; again after every `git pull`)
# ===========================================================================
cd /data/ascher02/uqmmune1/EFGPP/EFGPP2
git pull
mamba env create -f environment.yml || mamba env update -n efgpp -f environment.yml
conda activate efgpp
efgpp --version

# ===========================================================================
# 1. Project folder (everything below is downloaded and written here)
# ===========================================================================
mkdir -p /data/ascher02/uqmmune1/EFGPP/EFGPP2/my_project
cd /data/ascher02/uqmmune1/EFGPP/EFGPP2/my_project
[ -f project.yaml ] || efgpp init .

# ===========================================================================
# 2. Tools, each in its own conda env under ./software/envs/<tool>
#    plink2, plink, bcftools/tabix/bgzip, flashpca, multiqc, oc (OpenCRAVAT),
#    vep (+perl), predixcan, gwaslab
# ===========================================================================
efgpp setup tools
# efgpp setup tools vep --force                 # reinstall one tool if it failed

# ===========================================================================
# 3. Toolkits (several GB; R packages compile for a while)
#    perl, r (bigsnpr/LDpred-2, lassosum, ...), prs (PRSice-2, PRScs, LDSC, AnnoPred, ...),
#    prs-python, prs-py27, simulation (simuPOP, msprime, ...)
# ===========================================================================
efgpp setup toolkit --list
efgpp setup toolkit all

# ===========================================================================
# 4. Reference data into ./resources (all GRCh38)
# ===========================================================================
efgpp resources install genome --build GRCh38   # FASTA + liftover chains (GRCh37 -> GRCh38)
efgpp resources install vep                     # VEP cache (large, one time)
efgpp resources install clinvar                 # clinical significance per variant
efgpp resources install alphamissense           # missense pathogenicity scores
efgpp resources install predictdb               # GTEx v8 MASHR models for PrediXcan (49 tissues)
# efgpp resources install pgs_catalog           # set pgs_catalog.score_ids in resources.yaml first
efgpp resources list

# ===========================================================================
# 5. Check that everything is installed (✓ / ✗ with the reason)
# ===========================================================================
efgpp doctor
efgpp setup check
eval "$(efgpp setup path)"                      # ./software/bin on PATH for this shell (optional)

# ===========================================================================
# 6. Add the data (columns, types, roles and genome build are inferred;
#    original column names are kept; genotype is referenced, never copied)
# ===========================================================================
efgpp data add genotype   --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine
efgpp phenotype add       --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.height
efgpp data add covariates --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.cov
efgpp data add gwas       --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.gz \
    --trait migraine --phenotype migraine --ancestry European
#   add what you know: --n-cases <number> --n-controls <number> --study <name>
#   a GWAS of another trait for migraine:
#   efgpp data add gwas --path <depression GWAS> --trait depression --phenotype migraine

efgpp data inspect                              # what was inferred
efgpp phenotype list
efgpp data gwas list --phenotype migraine

# ===========================================================================
# 7. Per-variant annotation of the genotype variants + predicted expression
#    (run inside `efgpp data prepare`)
# ===========================================================================
efgpp resources enable vep clinvar alphamissense
efgpp data predict --tissue Whole_Blood         # add more: --tissue Brain_Cortex ...

# ===========================================================================
# 8. Validate, run, inspect, freeze (stop if validate shows ✗)
#    prepare: liftover to GRCh38 -> genotype QC -> PCA/kinship -> VEP/ClinVar/AlphaMissense
#             -> PrediXcan -> availability -> report
# ===========================================================================
efgpp data validate
efgpp data plan                                 # every step, and why a step is skipped
efgpp data prepare --cores 8
efgpp data availability
efgpp data report                               # reports/data/index.html
efgpp data freeze --name migraine_v1            # immutable snapshot for the next layer
efgpp data snapshots
efgpp data verify migraine_v1

# ===========================================================================
# Fixing things
# ===========================================================================
# efgpp data remove COV001                      # a source added wrongly (also GWAS001, PH001 via phenotype remove)
# efgpp modules export covariates               # edit modules/covariates.py for this project, then re-add
# efgpp data gwas run GWAS001 --force           # re-run GWASLab
# efgpp resources enable clinvar --off          # switch an annotation off
# efgpp export slurm                            # sbatch scripts in hpc/ instead of running locally
# logs: ./software/logs/<tool>.install.log, ./logs/
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
| variant annotations (VEP, ClinVar, AlphaMissense) | `data/derived/variant_annotations/GENO001/GENO001_{vep,clinvar,alphamissense}.parquet` |
| predicted expression | `data/predicted/expression/Whole_Blood/predicted_expression_Whole_Blood.parquet` |
| report, snapshot | `reports/data/index.html`, `snapshots/migraine_v1.yaml` + `snapshots/migraine_v1/` |

Notes

- Samples are matched on `IID` across `.fam`, `.cov` and `.height`, never by row order.
- The 2 variants at position 0 are unplaced array probes: kept in the original, left out of the
  GRCh38 copy (see the liftover report).
- `migraine_QC.*` is not needed (EFGPP runs its own QC); the other `migraine*.txt`, `.PRSCS`, …
  are old PRS results.
- DBSLMM's `dbslmm` binary is on Google Drive only: put it in `./software/opt/DBSLMM/software/dbslmm`.
- More detail: `Document.MD`, sections 11–13; all install commands: `Commands.md`.
