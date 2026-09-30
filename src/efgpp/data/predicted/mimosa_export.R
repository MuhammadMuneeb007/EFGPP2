# Export MIMOSA DNA-methylation prediction models (Zenodo 8400313, MIMOSA-Models.zip) to a long TSV.
#
# Each cg*.rds is a list of five method lists in this order (upstream MIMOSA-MWAS.R):
#   ElNet, MNet, SCAD, MCP, LASSO
# and each method list holds, by position:
#   [[1]] TRUE/FALSE  a satisfactory model exists
#   [[2]] data.frame  mQTL SNPs (SNP = rsID, SNPChr, position, a1, a2, CpG, ...)
#   [[3]] weights     one per SNP row in [[2]], for allele a1
#   [[4]] test R2     (Framingham Heart Study test data)
#   [[5]] lambda
# Every method is exported (valid or not); EFGPP selects the model afterwards in Python.
# Only base R is used.
#
# Usage: Rscript mimosa_export.R <models_dir> <out.tsv> [<structure.txt>]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript mimosa_export.R <models_dir> <out.tsv> [<structure.txt>]")
models_dir <- args[[1]]
out_file <- args[[2]]
structure_file <- if (length(args) >= 3) args[[3]] else NULL
methods <- c("ElNet", "MNet", "SCAD", "MCP", "LASSO")

files <- list.files(models_dir, pattern = "^cg.*\\.rds$", full.names = TRUE, recursive = TRUE)
if (length(files) == 0) stop(paste("no cg*.rds files under", models_dir))

pick_col <- function(df, candidates) {
  hit <- names(df)[tolower(names(df)) %in% tolower(candidates)]
  if (length(hit) == 0) return(rep(NA, nrow(df)))
  as.character(df[[hit[1]]])
}

con <- file(out_file, open = "w")
writeLines(paste("cpg_id", "method", "method_index", "valid", "test_r2", "lambda", "snp", "chromosome",
                 "position", "a1", "a2", "weight", sep = "\t"), con)
wrote_structure <- FALSE
n_ok <- 0
for (f in files) {
  cpg <- sub("\\.rds$", "", basename(f))
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj) || !is.list(obj) || length(obj) < 5) {
    message("skipping ", basename(f), ": not a list of five method lists")
    next
  }
  if (!wrote_structure && !is.null(structure_file)) {
    capture.output(str(obj, max.level = 2), file = structure_file)
    wrote_structure <- TRUE
  }
  for (i in seq_len(5)) {
    m <- obj[[i]]
    valid <- isTRUE(as.logical(m[[1]]))
    r2 <- suppressWarnings(as.numeric(m[[4]]))
    r2 <- if (length(r2) == 0) NA else r2[1]
    lambda <- suppressWarnings(as.numeric(m[[5]]))
    lambda <- if (length(lambda) == 0) NA else lambda[1]
    snps <- m[[2]]
    w <- suppressWarnings(as.numeric(m[[3]]))
    if (!is.data.frame(snps) || nrow(snps) == 0 || length(w) != nrow(snps)) {
      writeLines(paste(cpg, methods[i], i, valid, r2, lambda, NA, NA, NA, NA, NA, NA, sep = "\t"), con)
      next
    }
    rows <- paste(cpg, methods[i], i, valid, r2, lambda,
                  pick_col(snps, c("SNP", "rsid", "snp")),
                  pick_col(snps, c("SNPChr", "chr", "chromosome", "CHR")),
                  pick_col(snps, c("SNPPos", "SNPpos", "pos", "position", "BP", "bp")),
                  pick_col(snps, c("a1", "A1")), pick_col(snps, c("a2", "A2")), w, sep = "\t")
    writeLines(rows, con)
  }
  n_ok <- n_ok + 1
}
close(con)
cat(sprintf("EFGPP_MIMOSA\t%d\t%d\n", length(files), n_ok))
