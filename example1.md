# Example 1 — add the migraine PLINK cohort

Files in `/data/ascher02/uqmmune1/ANNOVAR/migraine`:

| File | Content | Becomes |
|---|---|---|
| `migraine.bed` / `.bim` / `.fam` | genotypes (733 samples, 619,653 variants, GRCh37) | genotype `GENO001` |
| `migraine.height` | phenotype: header `IID FID Height`, `1` = control, `2` = case | phenotype `PH001` |
| `migraine.cov` | covariates: header `FID IID <covariates…>` | covariates `COV001` |
| `migraine.gz` | GWAS summary statistics: `CHR BP SNP A1 A2 N SE P OR INFO MAF` | GWAS `GWAS001` (GWASLab) |

`migraine_QC.*` is not needed (EFGPP runs its own QC). The other `migraine*.txt`, `.PRSCS`, … are
old PRS results.

**Everything ends up in GRCh38.** The genotypes are GRCh37: EFGPP detects this (reference bases
from Ensembl for GRCh37 and GRCh38, confirmed with pyliftover), then lifts them to GRCh38 with
pyliftover before QC. The GWAS file goes through GWASLab (`basic_check`, `infer_build`,
`liftover` to GRCh38). The UCSC chain files are downloaded once into `resources/liftover/`.
The original files are never changed.

```bash
# ---------------------------------------------------------------------------
# Setup: go to your EFGPP project
# ---------------------------------------------------------------------------
conda activate efgpp
cd /data/ascher02/uqmmune1/EFGPP/EFGPP2/my_project
[ -f project.yaml ] || efgpp init .

# ---------------------------------------------------------------------------
# 0. GWASLab in its own conda environment (only needed once per project)
# ---------------------------------------------------------------------------
efgpp setup tools gwaslab

# ---------------------------------------------------------------------------
# 1. Genotype: path WITHOUT extension; used where it is, never copied.
#    Build detected automatically (GRCh37) and lifted to GRCh38 by `efgpp data prepare`.
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
# 4. GWAS summary statistics: GWASLab loads, checks, infers the build, lifts to GRCh38
# ---------------------------------------------------------------------------
efgpp data add gwas --path /data/ascher02/uqmmune1/ANNOVAR/migraine/migraine.gz --trait migraine \
    --col chrom=CHR --col pos=BP --col snpid=SNP --col ea=A1 --col nea=A2 \
    --col n=N --col se=SE --col p=P --col OR=OR --col info=INFO

# ---------------------------------------------------------------------------
# 5. Check, run, inspect, freeze (run one at a time; stop if validate shows ✗)
#    validate shows: "genome build GRCh37 (... reference_bases_pyliftover)"
#                    "will be lifted GRCh37 -> GRCh38 with pyliftover"
#    prepare runs:   harmonize.GENO001 (liftover) -> genotype_qc -> pca ... and gwas.GWAS001
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
- The 2 variants at position 0 are unplaced array probes: reported as a warning and kept
  (they cannot be lifted, so the GRCh38 copy leaves them out; see the liftover report).
- Results: lifted genotype `data/derived/genotype_qc/GENO001/liftover/GENO001_GRCh38.*` with
  `…_liftover_report.parquet`; GWAS `resources/gwas/GWAS001/GWAS001.GRCh38.parquet` with
  `GWAS001.gwaslab_report.json`.
- `efgpp data prepare` runs locally with EFGPP's built-in executor and prints errors directly.
  Snakemake is only used for the cluster: `efgpp data prepare --executor slurm`.
- More detail: `Document.MD`, section 11.
