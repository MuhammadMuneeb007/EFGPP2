"""Optional HLA imputation with HIBAG (Bioconductor), origin DERIVED.

HIBAG predicts classical HLA alleles from SNP genotypes with a classifier trained for a given
ancestry, genotyping platform and genome build. EFGPP never picks one population model for
everyone: the classifier is installed explicitly with its metadata
(`efgpp resources install hibag --url ... --ancestry ... --platform ... --model-build ... --loci ...`)
and the genotype must be in the classifier's build (lift explicitly with `efgpp data liftover`).

Output: data/derived/hla/<GENOTYPE_ID>/hla_types.parquet
    participant_id, locus, allele1, allele2, posterior_probability, matching, classifier, training_population
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import polars as pl

from efgpp.constants import ArtifactStatus, GenomeBuild, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.io import write_parquet
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths

HIBAG_SCRIPT = r"""
suppressMessages(library(HIBAG))
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[[1]]; classifier <- args[[2]]; assembly <- args[[3]]; out <- args[[4]]
loci <- strsplit(args[[5]], ",")[[1]]
geno <- hlaBED2Geno(bed.fn = paste0(prefix, ".bed"), fam.fn = paste0(prefix, ".fam"),
                    bim.fn = paste0(prefix, ".bim"), assembly = assembly)
models <- get(load(classifier))
rows <- list()
for (locus in loci) {
  if (is.null(models[[locus]])) { message("classifier has no model for HLA-", locus); next }
  m <- hlaModelFromObj(models[[locus]])
  pred <- hlaPredict(m, geno, type = "response+prob")
  v <- pred$value
  rows[[locus]] <- data.frame(sample_id = v$sample.id, locus = locus, allele1 = v$allele1,
                              allele2 = v$allele2, posterior_probability = v$prob,
                              matching = if ("matching" %in% names(v)) v$matching else NA)
  hlaClose(m)
}
res <- do.call(rbind, rows)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
"""


def _classifier(project: Project) -> dict[str, Any] | None:
    from efgpp.data.predicted.engine import installed_resource

    return installed_resource(project, "hibag", project.data.hla.classifier)


def hla_ready(project: Project) -> str | None:
    """Why HLA imputation cannot run (None = ready)."""
    from efgpp.data.references.molecular import rscript_for

    if not project.data.hla.classifier:
        return "no HIBAG classifier chosen (hla.classifier; efgpp resources install hibag ...)"
    if _classifier(project) is None:
        return f"HIBAG classifier {project.data.hla.classifier!r} not installed (efgpp resources install hibag ...)"
    if rscript_for(project, ("hla",)) is None:
        return "HIBAG not installed (efgpp setup toolkit hla)"
    return None


def hla_genotype(project: Project, source_id: str, build: str) -> Artifact:
    """A genotype of `source_id` in the classifier's build: the source itself or an explicit liftover."""
    want = GenomeBuild.normalize(build)
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        cands = [store.latest(source_id=source_id, artifact_type="qc_genotype"),
                 *store.find(source_id=source_id, artifact_type="liftover_genotype"),
                 store.latest(source_id=source_id, artifact_type="source")]
    for a in cands:
        if a is not None and GenomeBuild.normalize(a.genome_build) == want:
            return a
    raise RuntimeError(f"the HIBAG classifier is {want.value}, but {source_id} has no genotype in that build; "
                       f"create one explicitly: efgpp data liftover <artifact id> --to {want.value}")


def run_hla(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> dict[str, Any]:
    from efgpp.data.genotype.carriers import genotype_samples, write_assays
    from efgpp.data.genotype.formats import resolve_fileset
    from efgpp.data.provenance import run_tool
    from efgpp.data.references.molecular import rscript_for

    problem = hla_ready(project)
    if problem:
        raise RuntimeError(problem)
    res = _classifier(project)
    assert res is not None
    meta = res.get("metadata") or {}
    build = str(meta.get("genome_build") or res.get("genome_build"))
    genotype = hla_genotype(project, source_id, build)
    fs = resolve_fileset(Path(genotype.path), genotype.format)
    work = project.work_root / "hla" / source_id
    work.mkdir(parents=True, exist_ok=True)
    prefix = work / "chr6"
    run_tool(project, "plink2", [*fs.plink_input_args(), "--chr", "6", "--make-bed", "--out", str(prefix),
                                 "--threads", str(threads)], step_id=step_id,
             inputs=[genotype.artifact_id])  # type: ignore[list-item]
    script = work / "hibag_predict.R"
    script.write_text(HIBAG_SCRIPT, encoding="utf-8")
    out = work / "hla_raw.tsv"
    assembly = "hg38" if GenomeBuild.normalize(build) == GenomeBuild.GRCH38 else "hg19"
    rscript = rscript_for(project, ("hla",))
    loci = ",".join(project.data.hla.loci)
    proc = subprocess.run([str(rscript), str(script), str(prefix), str(project.resolve(res["local_path"])),
                           assembly, str(out), loci], capture_output=True, text=True)
    (work / "hibag.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"HIBAG failed (log: {work / 'hibag.log'}):\n{proc.stderr[-1500:]}")
    samples = genotype_samples(project, source_id, resolve_fileset(prefix, "bed"))
    raw = pl.read_csv(out, separator="\t", infer_schema=False)
    table = raw.join(samples.select(pl.col("IID").alias("sample_id"), "participant_id"), on="sample_id", how="left") \
        .select("participant_id", pl.col("sample_id").alias("native_sample_id"), "locus", "allele1", "allele2",
                pl.col("posterior_probability").cast(pl.Float64, strict=False), "matching",
                pl.lit(res.get("version")).alias("classifier"), pl.lit(meta.get("ancestry")).alias("training_population"))
    out_dir = project.data_root / Origin.DERIVED.value / "hla" / source_id
    path = write_parquet(table, out_dir / "hla_types.parquet")
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=f"{source_id}_hla", artifact_type="hla", modality=Modality.HLA, origin=Origin.DERIVED,
            status=ArtifactStatus.READY, path=str(path), format="parquet", size=path.stat().st_size,
            checksum=str(checksum_paths([path])), participant_count=table.get_column("participant_id").n_unique(),
            feature_count=len(project.data.hla.loci), genome_build=build, temporal_type=TemporalType.STATIC,
            source_id=source_id, parent_artifact_ids=[genotype.artifact_id],  # type: ignore[list-item]
            tool="HIBAG", resource_versions={"hibag_classifier": str(res.get("version")),
                                             "classifier_sha256": str(res.get("checksum"))},
            metadata={"classifier": meta, "loci": project.data.hla.loci,
                      "note": "imputed HLA alleles with posterior probabilities; phenotype not used"},
        ))
        write_assays(reg, art, f"{source_id}_hla", Modality.HLA,
                     table.get_column("participant_id").drop_nulls().to_list())
    return {"artifact": art.artifact_id, "rows": table.height}
