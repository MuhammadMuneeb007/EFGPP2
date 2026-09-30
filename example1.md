# Example 1 — the migraine PLINK cohort, from a fresh Linux account to a frozen snapshot

This walks through everything: installing conda/mamba, cloning EFGPP, creating its
environment, installing the tools and reference data, adding the migraine data, running the
Data layer and freezing a snapshot.

**Input files** in `/data/ascher02/uqmmune1/ANNOVAR/migraine`:

| File | Content | Becomes |
|---|---|---|
| `migraine.bed` / `.bim` / `.fam` | genotypes (733 samples, 619,653 variants, GRCh37) | genotype `GENO001` |
| `migraine.height` | phenotype: header `IID FID Height`, `1` = control, `2` = case | phenotype `PH001` (name `migraine`) |
| `migraine.cov` | covariates: header `FID IID <covariates…>` | covariates `COV001` |
| `migraine.gz` | GWAS summary statistics: `CHR BP SNP A1 A2 N SE P OR INFO MAF` | GWAS `GWAS001` (GWASLab) |

**Folders used**

| Folder | What goes there |
|---|---|
| `/data/ascher02/uqmmune1/miniforge3` | conda + mamba (only if you do not have them yet) |
| `/data/ascher02/uqmmune1/EFGPP/EFGPP2` | the EFGPP code (git clone) |
| `/data/ascher02/uqmmune1/EFGPP/EFGPP2/my_project` | the project: `software/` (tools, one conda env each), `resources/` (reference data), `data/`, `reports/`, `snapshots/`, `logs/` |

Nothing is installed in your home directory. **Requirements:** Linux x86_64, `git`, internet
access (run the install steps on the login node), and tens of GB of free space in `/data`
(the VEP cache and the toolkits are large).

```bash
# ###########################################################################
# PART A — INSTALLATION (once)
# ###########################################################################

# ---------------------------------------------------------------------------
# A0. Long downloads: work inside tmux (or screen) so a dropped SSH connection
#     does not kill them. Re-attach later with: tmux attach -t efgpp
# ---------------------------------------------------------------------------
# tmux new -s efgpp                              # start it first, then run the rest inside it

# ---------------------------------------------------------------------------
# A1. conda + mamba
#     Already have them?  `mamba --version` prints a version -> skip to A2.
#     Otherwise install Miniforge (conda-forge's installer, includes mamba) into /data:
# ---------------------------------------------------------------------------
mamba --version || {
  cd /data/ascher02/uqmmune1
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash Miniforge3-Linux-x86_64.sh -b -p /data/ascher02/uqmmune1/miniforge3
  rm Miniforge3-Linux-x86_64.sh
}
source /data/ascher02/uqmmune1/miniforge3/etc/profile.d/conda.sh   # conda in this shell
# (optional, makes `conda` available in every new shell; edits ~/.bashrc only)
# /data/ascher02/uqmmune1/miniforge3/bin/conda init bash

# ---------------------------------------------------------------------------
# A2. Get the EFGPP code (private repository)
#     HTTPS asks for your GitHub username and a personal access token as password
#     (GitHub -> Settings -> Developer settings -> Personal access tokens, scope "repo").
#     With an SSH key on GitHub use: git clone git@github.com:MuhammadMuneeb007/EFGPP2.git
# ---------------------------------------------------------------------------
mkdir -p /data/ascher02/uqmmune1/EFGPP
cd /data/ascher02/uqmmune1/EFGPP
[ -d EFGPP2 ] || git clone https://github.com/MuhammadMuneeb007/EFGPP2.git
cd EFGPP2
git pull                                          # latest version

# ---------------------------------------------------------------------------
# A3. The `efgpp` environment (EFGPP itself only; tools get their own envs later)
#     create the first time, update after every `git pull`
# ---------------------------------------------------------------------------
mamba env create -f environment.yml || mamba env update -n efgpp -f environment.yml
conda activate efgpp
efgpp --version
efgpp --help
# optional self-test (~2 minutes):  pytest -q

# ---------------------------------------------------------------------------
# A4. Create the project (everything below is downloaded and written inside it)
# ---------------------------------------------------------------------------
mkdir -p /data/ascher02/uqmmune1/EFGPP/EFGPP2/my_project
cd /data/ascher02/uqmmune1/EFGPP/EFGPP2/my_project
[ -f project.yaml ] || efgpp init .
ls                                                # project.yaml data.yaml resources.yaml software/ resources/ ...

# ---------------------------------------------------------------------------
# A5. Data-layer tools, each in its own conda env under ./software/envs/<tool>
#     plink2, plink, bcftools/tabix/bgzip, flashpca, multiqc, oc (OpenCRAVAT),
#     vep (+perl), predixcan (MetaXcan), gwaslab
# ---------------------------------------------------------------------------
efgpp setup tools
# efgpp setup tools vep --force                   # reinstall a single tool if it failed

# ---------------------------------------------------------------------------
# A6. Toolkits (several GB; the R packages compile for a while)
#     perl | r (bigsnpr/LDpred-2, lassosum, sim1000G, ...) | prs (PRSice-2, PRScs, LDSC,
#     AnnoPred, GCTA, ...) | prs-python | prs-py27 | simulation (simuPOP, msprime, ...)
# ---------------------------------------------------------------------------
efgpp setup toolkit --list
efgpp setup toolkit all
# DBSLMM's `dbslmm` binary is on Google Drive only (see github.com/biostat0903/DBSLMM):
# put it in ./software/opt/DBSLMM/software/dbslmm

# ---------------------------------------------------------------------------
# A7. Reference data into ./resources (everything GRCh38)
# ---------------------------------------------------------------------------
efgpp resources install genome --build GRCh38     # FASTA + liftover chains (GRCh37 -> GRCh38)
efgpp resources install vep                       # VEP cache (large, one time)
efgpp resources install clinvar                   # clinical significance per variant
efgpp resources install alphamissense             # missense pathogenicity scores
efgpp resources install predictdb                 # GTEx v8 models for PrediXcan (49 tissues)
# efgpp resources install pgs_catalog             # first set pgs_catalog.score_ids in resources.yaml
efgpp resources list

# ---------------------------------------------------------------------------
# A8. Check the installation: every item ✓, or ✗ with the reason
# ---------------------------------------------------------------------------
efgpp doctor                                      # data-layer tools with versions
efgpp setup check                                 # tools, toolkits, R packages, repositories
ls software/envs                                  # one conda environment per tool

# ###########################################################################
# PART B — EVERY NEW SESSION
# ###########################################################################
source /data/ascher02/uqmmune1/miniforge3/etc/profile.d/conda.sh   # not needed after `conda init`
conda activate efgpp
cd /data/ascher02/uqmmune1/EFGPP/EFGPP2/my_project

# ###########################################################################
# PART C — THE MIGRAINE DATA
# ###########################################################################

# ---------------------------------------------------------------------------
# C1. Add the data. Columns, types, roles and genome build are inferred; original
#     column names are kept; the genotype is referenced in place, never copied.
# ---------------------------------------------------------------------------
efgpp data add genotype   --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine
efgpp phenotype add       --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.height
efgpp data add covariates --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.cov
efgpp data add gwas       --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.gz \
    --trait migraine --phenotype migraine --ancestry European
#   add what you know about the study: --n-cases <number> --n-controls <number> --study <name>
#   a GWAS of another trait used for migraine:
#   efgpp data add gwas --path <depression GWAS file> --trait depression --phenotype migraine

# ---------------------------------------------------------------------------
# C2. Look at what was inferred
# ---------------------------------------------------------------------------
efgpp data inspect
efgpp phenotype list
efgpp data gwas list --phenotype migraine
cat data.yaml

# ---------------------------------------------------------------------------
# C3. Switch on variant annotation and predicted expression (run by `prepare`)
# ---------------------------------------------------------------------------
efgpp resources enable vep clinvar alphamissense
efgpp data predict --tissue Whole_Blood           # more tissues: --tissue Brain_Cortex ...

# ---------------------------------------------------------------------------
# C4. Validate and run (stop and fix if validate shows ✗)
#     prepare: liftover to GRCh38 -> genotype QC -> PCA / kinship
#              -> VEP / ClinVar / AlphaMissense -> PrediXcan -> availability -> report
# ---------------------------------------------------------------------------
efgpp data validate
efgpp data plan                                   # every step, and why a step is skipped
efgpp data prepare --cores 8
efgpp data availability                           # who has genotype + phenotype + covariates
efgpp data report                                 # open reports/data/index.html in a browser

# ---------------------------------------------------------------------------
# C5. Freeze an immutable snapshot for the next layer
# ---------------------------------------------------------------------------
efgpp data freeze --name migraine_v1
efgpp data snapshots
efgpp data verify migraine_v1

# ###########################################################################
# PART D — IF SOMETHING GOES WRONG
# ###########################################################################
# efgpp data remove COV001                        # remove a wrongly added source (also GWAS001)
# efgpp phenotype remove PH001
# efgpp modules export covariates                 # edit modules/covariates.py for this project, re-add
# efgpp data gwas run GWAS001 --force             # re-run GWASLab
# efgpp resources enable clinvar --off            # switch an annotation off
# efgpp export slurm                              # sbatch scripts in hpc/ instead of running locally
# logs: ./software/logs/<tool>.install.log  and  ./logs/
# never `mamba install` a tool into the efgpp env; use `efgpp setup tools <tool> --force`
```

**Where the results are**

| What | File |
|---|---|
| phenotype, original column name | `phenotypes/PH001/phenotype.parquet` (`participant_id`, `Height`) |
| covariates, original column names | `data/observed/covariates/COV001/covariates.parquet` |
| genotype lifted to GRCh38 | `data/derived/genotype_qc/GENO001/liftover/GENO001_GRCh38.*` + `…_liftover_report.parquet` |
| QC-passed genotype (GRCh38) | `data/derived/genotype_qc/GENO001/GENO001_qc.*` |
| GWAS, original column names (GRCh38) | `resources/gwas/GWAS001/GWAS001.GRCh38.parquet` |
| GWAS, GWASLab names (GRCh38) | `resources/gwas/GWAS001/GWAS001.GRCh38.gwaslab.parquet` + `GWAS001.gwaslab_report.json` |
| variant annotations | `data/derived/variant_annotations/GENO001/GENO001_{vep,clinvar,alphamissense}.parquet` |
| predicted expression | `data/predicted/expression/Whole_Blood/predicted_expression_Whole_Blood.parquet` |
| report, snapshot | `reports/data/index.html`, `snapshots/migraine_v1.yaml` + `snapshots/migraine_v1/` |

**Notes**

- Samples are matched on `IID` across `.fam`, `.cov` and `.height`, never by row order.
- The 2 variants at position 0 are unplaced array probes: kept in the original, left out of the
  GRCh38 copy (see the liftover report).
- `migraine_QC.*` is not needed (EFGPP runs its own QC); the other `migraine*.txt`, `.PRSCS`, …
  are old PRS results.
- More detail: `Document.MD` (sections 11–13) and `Commands.md`.
