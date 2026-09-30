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

# ---------------------------------------------------------------------------
# 2. Project folder (outside the EFGPP2 code repository)
# ---------------------------------------------------------------------------
mkdir -p /data/ascher02/uqmmune1/EFGPP/efgpp_projects/my_project
cd /data/ascher02/uqmmune1/EFGPP/efgpp_projects/my_project
[ -f project.yaml ] || efgpp init .              # only the first time

# ---------------------------------------------------------------------------
# 3. Install every missing tool from its official source
#    (tools already in the efgpp mamba env are detected and skipped)
# ---------------------------------------------------------------------------
module load gcc 2>/dev/null || true              # compiler for the bcftools/htslib build, if needed
module load apptainer 2>/dev/null || true        # container runtime for VEP, if available
efgpp setup tools                                # all tools, into ./.efgpp
# efgpp setup tools --shared                     # alternative: once for all projects (~/.local/share/efgpp)
# efgpp setup tools bcftools --force             # retry / reinstall one tool

# ---------------------------------------------------------------------------
# 4. Check what is installed
# ---------------------------------------------------------------------------
efgpp doctor                                     # ✓ installed   ✗ required missing   ! optional missing

eval "$(efgpp setup path)"                       # put EFGPP's tool folders on PATH for this shell
for t in plink2 plink bcftools tabix bgzip flashpca snakemake multiqc oc predixcan vep; do
  printf '%-10s ' "$t"
  command -v "$t" >/dev/null && echo "OK       $(command -v "$t")" || echo "MISSING"
done

# ---------------------------------------------------------------------------
# 5. VEP annotation data (only needed if you use VEP; large one-time download)
# ---------------------------------------------------------------------------
efgpp resources install genome --build GRCh38
efgpp resources install vep --build GRCh38
efgpp resources list

# ---------------------------------------------------------------------------
# If something fails
# ---------------------------------------------------------------------------
# - the reason is printed next to each failed tool by `efgpp setup tools`
# - bcftools build log:  .efgpp/opt/build/build.log
# - Python tool logs:    .efgpp/envs/<name>.install.log
```
