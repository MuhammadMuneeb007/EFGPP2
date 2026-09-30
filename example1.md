# Example 1 — add the migraine PLINK cohort

Files in `/data/ascher02/uqmmune1/ANNOVAR/migraine`:

| File | Content | Becomes |
|---|---|---|
| `migraine.bed` / `.bim` / `.fam` | genotypes (733 samples, 619,653 variants, GRCh37) | genotype `GENO001` |
| `migraine.height` | phenotype: header `IID FID Height`, `1` = control, `2` = case | phenotype `PH001` |
| `migraine.cov` | covariates: header `FID IID <covariates…>` | covariates `COV001` |

`migraine_QC.*` is not needed (EFGPP runs its own QC). `migraine*.txt`, `.gz`, `.PRSCS`, … are
GWAS summary statistics and old PRS results, not participant data.

```bash
# ---------------------------------------------------------------------------
# Setup: go to your EFGPP project
# ---------------------------------------------------------------------------
conda activate efgpp
cd /data/ascher02/uqmmune1/EFGPP/EFGPP2/my_project
[ -f project.yaml ] || efgpp init .

# ---------------------------------------------------------------------------
# 1. Genotype: path WITHOUT extension; used where it is, never copied.
#    The build is detected automatically (positions past the GRCh38 chromosome ends -> GRCh37).
# ---------------------------------------------------------------------------
efgpp data add genotype --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine --format bed --build auto --mode reference

# ---------------------------------------------------------------------------
# 2. Phenotype: binary; PLINK 1/2 coding is detected automatically (2 = case)
# ---------------------------------------------------------------------------
efgpp phenotype add --name migraine --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.height \
    --id-column IID --value-column Height --type binary

# ---------------------------------------------------------------------------
# 3. Covariates: every column except FID/IID is taken automatically
#    (optional: --columns <names> to choose, --categorical <names> for codes such as sex)
# ---------------------------------------------------------------------------
efgpp data add covariates --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.cov --id-column IID

# ---------------------------------------------------------------------------
# 4. Check, run, inspect, freeze (run one at a time; stop if validate shows ✗)
# ---------------------------------------------------------------------------
efgpp data validate                    # every problem, with counts
efgpp data prepare                     # QC, PCA, kinship, availability, report
efgpp data availability                # who has genotype + phenotype + covariates
efgpp data report                      # reports/data/index.html
efgpp data freeze --name migraine_v1   # immutable snapshot for the next layer

# ---------------------------------------------------------------------------
# Fixing a source that was added wrongly (e.g. wrong covariate columns)
# ---------------------------------------------------------------------------
# efgpp data remove COV001
# efgpp data add covariates --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.cov --id-column IID

# ---------------------------------------------------------------------------
# Alternative: register the whole legacy folder at once (source never modified)
# ---------------------------------------------------------------------------
# efgpp data migrate --source /data/ascher02/uqmmune1/ANNOVAR/migraine --profile legacy-efgpp --dry-run
# efgpp data migrate --source /data/ascher02/uqmmune1/ANNOVAR/migraine --profile legacy-efgpp --apply
```

Notes

- Samples are matched on `IID` across the `.fam`, `.cov` and `.height` files, never by row order.
- The 2 variants at position 0 are unplaced array probes: reported as a warning and kept.
- `efgpp data prepare` runs locally with EFGPP's built-in executor and prints errors directly.
  Snakemake is only used for the cluster: `efgpp data prepare --executor slurm`.
- More detail: `Document.MD`, section 11.
