# EFGPP — install and check all tools

Run on the **login node** (the downloads need internet access). Everything is installed inside
the project directory you are in — tools in `software/`, reference data in `resources/`, logs in
`logs/` — never in your home directory and never into the `efgpp` environment itself.

```bash
# ---------------------------------------------------------------------------
# 1. Update EFGPP and its environment (once, and after every `git pull`)
#    The efgpp env only holds EFGPP itself; every tool gets its own conda env.
# ---------------------------------------------------------------------------
cd /data/ascher02/uqmmune1/EFGPP/EFGPP2
git pull
mamba env update -n efgpp -f environment.yml     # first time: mamba env create -f environment.yml
conda activate efgpp
efgpp --version

# ---------------------------------------------------------------------------
# 2. Go to the project directory (everything below is installed here)
# ---------------------------------------------------------------------------
mkdir -p /data/ascher02/uqmmune1/EFGPP/efgpp_projects/my_project
cd /data/ascher02/uqmmune1/EFGPP/efgpp_projects/my_project
[ -f project.yaml ] || efgpp init .              # only the first time

# ---------------------------------------------------------------------------
# 3. Data-layer tools: each in its OWN conda env under ./software/envs/<tool>
#    (plink2, plink, bcftools+tabix+bgzip, snakemake, multiqc, oc, vep+perl,
#    predixcan); official downloads only if conda fails; flashpca2 = official binary
# ---------------------------------------------------------------------------
efgpp setup tools                                # all tools -> ./software
# efgpp setup tools vep                          # one tool
# efgpp setup tools vep --force                  # reinstall

# ---------------------------------------------------------------------------
# 4. Toolkits (also into ./software): Perl, R + PRS R packages (LDpred-2 ...),
#    PRS tools from PRSTools, Python 3.10 / 2.7 method envs, simulation
#    (several GB; R packages compile for a while)
# ---------------------------------------------------------------------------
efgpp setup toolkit --list                       # what each toolkit contains
efgpp setup toolkit all                          # everything below in one go
# efgpp setup toolkit perl                       # Perl
# efgpp setup toolkit r                          # R, bigsnpr (LDpred-2, SCT, lassosum2), lassosum, PANPRSnext,
#                                                #   CTSLEB, RapidoPGS, EBPRS, R2BGLiMS, penRegSum, sim1000G, ...
# efgpp setup toolkit prs                        # PLINK, GCTA, GCTB, GEMMA, BOLT-LMM, PRSice-2, LDAK, vcftools,
#                                                #   PRScs, PRScsx, SDPR, DBSLMM, CTPR, NPS, XP-BLUP, smtpred,
#                                                #   PRSbils, LDpred-funct, PolyFun, LDSC, AnnoPred, PleioPred
#                                                #   (also installs prs-python, prs-py27 and r)
# efgpp setup toolkit simulation                 # simuPOP, msprime, tskit, stdpopsim

# ---------------------------------------------------------------------------
# 5. Reference data (into ./resources)
# ---------------------------------------------------------------------------
efgpp resources install genome --build GRCh38
efgpp resources install vep --build GRCh38       # VEP cache (large, one time)
efgpp resources list

# ---------------------------------------------------------------------------
# 6. Check what is installed
# ---------------------------------------------------------------------------
ls software/envs software/bin                    # one environment per tool, all commands
efgpp doctor                                     # data-layer tools with versions
efgpp setup check                                # every tool, toolkit, R package, repository: ✓ / ✗
# efgpp setup check --toolkits r,prs             # only some toolkits

eval "$(efgpp setup path)"                       # put ./software/bin on PATH for this shell
for t in plink2 plink bcftools tabix bgzip flashpca snakemake multiqc oc predixcan vep \
         perl R Rscript gcta64 gctb gemma bolt PRSice ldak vcftools PRScs.py SDPR ldsc.py; do
  printf '%-12s ' "$t"
  command -v "$t" >/dev/null && echo "OK       $(command -v "$t")" || echo "MISSING"
done
Rscript -e 'for (p in c("bigsnpr","lassosum","PANPRSnext","CTSLEB","RapidoPGS","EBPRS","R2BGLiMS","sim1000G"))
              cat(sprintf("%-12s %s\n", p, if (requireNamespace(p, quietly=TRUE)) "OK" else "MISSING"))'

# ---------------------------------------------------------------------------
# Manual step (cannot be automated)
# ---------------------------------------------------------------------------
# DBSLMM's `dbslmm` executable is distributed on Google Drive only: download it from
# https://github.com/biostat0903/DBSLMM (README) into ./software/opt/DBSLMM/software/dbslmm

# ---------------------------------------------------------------------------
# If something fails
# ---------------------------------------------------------------------------
# - the reason is printed next to each failed item
# - conda environment logs:  ./software/logs/<env>.install.log   (e.g. vep.install.log)
# - toolkit logs:            ./software/logs/<toolkit>.install.log
# - bcftools source build:   ./software/opt/build/build.log
# - never `mamba install` a tool into the efgpp env; `efgpp setup tools <tool> --force` instead
# - tools installed earlier into ./.efgpp are still found; delete .efgpp/bin and .efgpp/envs
#   once `efgpp setup check` is all ✓
```
