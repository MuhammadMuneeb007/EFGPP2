# EFGPP — install and check all tools

Run on the **login node** (the downloads need internet access).

```bash
# ---------------------------------------------------------------------------
# 1. Update EFGPP and its environment (once, and after every `git pull`)
# ---------------------------------------------------------------------------
cd /data/ascher02/uqmmune1/EFGPP/EFGPP2
git pull
mamba env update -n efgpp -f environment.yml     # first time: mamba env create -f environment.yml
conda activate efgpp
efgpp --version

# Shared tools live on /data, not in the home directory (quota). Add this line to ~/.bashrc
# so every session and every project finds them.
export EFGPP_TOOLS_HOME=/data/ascher02/uqmmune1/EFGPP/efgpp_tools

# ---------------------------------------------------------------------------
# 2. Project folder (outside the EFGPP2 code repository)
# ---------------------------------------------------------------------------
mkdir -p /data/ascher02/uqmmune1/EFGPP/efgpp_projects/my_project
cd /data/ascher02/uqmmune1/EFGPP/efgpp_projects/my_project
[ -f project.yaml ] || efgpp init .              # only the first time

# ---------------------------------------------------------------------------
# 3. Data-layer tools (PLINK 2/1.9, bcftools, FlashPCA, Snakemake, MultiQC,
#    OpenCRAVAT, MetaXcan, VEP + Perl); already-installed tools are skipped
# ---------------------------------------------------------------------------
module load gcc 2>/dev/null || true              # compiler for source builds, if needed
module load apptainer 2>/dev/null || true        # VEP container fallback, if available
efgpp setup tools                                # all tools, into ./.efgpp
# efgpp setup tools --shared                     # alternative: once for all projects ($EFGPP_TOOLS_HOME)
# efgpp setup tools bcftools --force             # retry / reinstall one tool

# ---------------------------------------------------------------------------
# 4. Toolkits: Perl, R + PRS R packages (LDpred-2 ...), PRS tools from PRSTools,
#    Python 3.10 / 2.7 method environments, simulation (simuPOP, msprime ...)
#    --shared installs once into $EFGPP_TOOLS_HOME for every project
#    (recommended: several GB, and R packages compile for a while)
# ---------------------------------------------------------------------------
efgpp setup toolkit --list                       # what each toolkit contains
efgpp setup toolkit all --shared                 # everything below in one go
# efgpp setup toolkit perl --shared              # Perl
# efgpp setup toolkit r --shared                 # R, bigsnpr (LDpred-2, SCT, lassosum2), lassosum, PANPRSnext,
#                                                #   CTSLEB, RapidoPGS, EBPRS, R2BGLiMS, penRegSum, sim1000G, ...
# efgpp setup toolkit prs --shared               # PLINK, GCTA, GCTB, GEMMA, BOLT-LMM, PRSice-2, LDAK, vcftools,
#                                                #   PRScs, PRScsx, SDPR, DBSLMM, CTPR, NPS, XP-BLUP, smtpred,
#                                                #   PRSbils, LDpred-funct, PolyFun, LDSC, AnnoPred, PleioPred
#                                                #   (also installs prs-python, prs-py27 and r)
# efgpp setup toolkit simulation --shared        # simuPOP, msprime, tskit, stdpopsim

# ---------------------------------------------------------------------------
# 5. Check what is installed
# ---------------------------------------------------------------------------
efgpp doctor                                     # data-layer tools with versions
efgpp setup check                                # every tool, toolkit, R package, repository: ✓ / ✗
# efgpp setup check --toolkits r,prs             # only some toolkits

eval "$(efgpp setup path)"                       # put EFGPP's tool folders on PATH for this shell
for t in plink2 plink bcftools tabix bgzip flashpca snakemake multiqc oc predixcan vep perl \
         R Rscript gcta64 gctb gemma bolt PRSice ldak vcftools PRScs.py SDPR ldsc.py; do
  printf '%-12s ' "$t"
  command -v "$t" >/dev/null && echo "OK       $(command -v "$t")" || echo "MISSING"
done
Rscript -e 'for (p in c("bigsnpr","lassosum","PANPRSnext","CTSLEB","RapidoPGS","EBPRS","R2BGLiMS","sim1000G"))
              cat(sprintf("%-12s %s\n", p, if (requireNamespace(p, quietly=TRUE)) "OK" else "MISSING"))'

# ---------------------------------------------------------------------------
# 6. VEP annotation data (only needed if you use VEP; large one-time download)
# ---------------------------------------------------------------------------
efgpp resources install genome --build GRCh38
efgpp resources install vep --build GRCh38
efgpp resources list

# ---------------------------------------------------------------------------
# Manual step (cannot be automated)
# ---------------------------------------------------------------------------
# DBSLMM's `dbslmm` executable is distributed on Google Drive only: download it from
# https://github.com/biostat0903/DBSLMM (README) into <install root>/opt/DBSLMM/software/dbslmm

# ---------------------------------------------------------------------------
# If something fails
# ---------------------------------------------------------------------------
# - the reason is printed next to each failed item
# - toolkit logs:        <install root>/toolkits/<name>.install.log
# - bcftools build log:  .efgpp/opt/build/build.log
# - Python tool logs:    .efgpp/envs/<name>.install.log
#   (<install root> is ./.efgpp, or $EFGPP_TOOLS_HOME with --shared)
```
