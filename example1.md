# Example 1 — add the migraine PLINK cohort

Files in `/data/ascher02/uqmmune1/ANNOVAR/migraine`:

| File | Content | Becomes |
|---|---|---|
| `migraine.bed` / `.bim` / `.fam` | genotypes | genotype `GENO001` |
| `migraine.height` | phenotype: header `IID FID Height`, `1` = control, `2` = case | phenotype `PH001` |
| `migraine.cov` | covariates: header `FID IID <covariates…>` | covariates `COV001` |

`migraine_QC.*` is not needed (EFGPP runs its own QC). `migraine*.txt`, `.gz`, `.PRSCS`, … are
GWAS summary statistics and old PRS results, not participant data.

```bash
# ---------------------------------------------------------------------------
# Setup: go to your EFGPP project
# ---------------------------------------------------------------------------
conda activate efgpp
cd /home/s4750311/efgpp_projects/my_project
[ -f project.yaml ] || efgpp init .

# ---------------------------------------------------------------------------
# 0. Look at the column names (needed for step 3)
# ---------------------------------------------------------------------------
head -3 /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.cov /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.height

# ---------------------------------------------------------------------------
# 1. Genotype: path WITHOUT extension; used where it is, never copied
#    (if validation cannot determine the build, use --build GRCh37: the UK Biobank
#     genotypes of the original EFGPP study are hg19)
# ---------------------------------------------------------------------------
efgpp data add genotype --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine --format bed --build auto --mode reference

# ---------------------------------------------------------------------------
# 2. Phenotype: binary; PLINK 1/2 coding is detected automatically (2 = case)
# ---------------------------------------------------------------------------
efgpp phenotype add --name migraine --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.height \
    --id-column IID --value-column Height --type binary

# ---------------------------------------------------------------------------
# 3. Covariates: every column of migraine.cov except FID and IID (from step 0);
#    code columns such as sex go in --categorical
#    e.g. --columns Sex Age PC1 PC2 PC3 PC4 PC5 --categorical Sex
# ---------------------------------------------------------------------------
efgpp data add covariates --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.cov --id-column IID \
    --columns COV1 COV2 COV3 --categorical COV1

# ---------------------------------------------------------------------------
# 4. Check, run, inspect, freeze
# ---------------------------------------------------------------------------
efgpp data validate                    # every problem, with counts
efgpp data prepare                     # QC, PCA, kinship, availability, report
efgpp data availability                # who has genotype + phenotype + covariates
efgpp data report                      # reports/data/index.html
efgpp data freeze --name migraine_v1   # immutable snapshot for the next layer

# ---------------------------------------------------------------------------
# Alternative: register the whole legacy folder at once (source never modified)
# ---------------------------------------------------------------------------
# efgpp data migrate --source /data/ascher02/uqmmune1/ANNOVAR/migraine --profile legacy-efgpp --dry-run
# efgpp data migrate --source /data/ascher02/uqmmune1/ANNOVAR/migraine --profile legacy-efgpp --apply
```

Samples are matched on `IID` across the `.fam`, `.cov` and `.height` files, never by row order.
More detail: `Document.MD`, section 11.
