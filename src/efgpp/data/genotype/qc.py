"""Phenotype-independent genotype QC with PLINK 2.

Passes (all on the same fileset, never modifying it):
  1. sample + variant missingness                          --missing
  2. variant metrics on retained samples                   --missing variant-only --freq --hardy
  3. LD pruning on retained samples / passing variants     --indep-pairwise
  4. heterozygosity and KING kinship on the pruned set     --het, --make-king-table
  5. sex check where X chromosome data allow               --check-sex
  6. QC-passed genotype                                    --keep --extract --make-pgen

Every decision is written to Parquet (sample_qc, variant_qc, kinship) with its reason.
HWE is computed on all retained samples: no phenotype is consulted.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import polars as pl

from efgpp.config.project import GenotypeQCConfig
from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.formats import (
    GenotypeFileset,
    count_variants,
    read_variants,
    resolve_fileset,
)
from efgpp.data.genotype.loader import read_plink_table, write_id_list, write_keep_file
from efgpp.data.genotype.relatedness import classify_pairs, greedy_unrelated_removal
from efgpp.data.io import write_parquet
from efgpp.data.provenance import RunRecord, run_tool
from efgpp.data.registry import Registry, to_json, utcnow
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths


@dataclass
class GenotypeInput:
    source_id: str
    fileset: GenotypeFileset
    samples: pl.DataFrame  # FID, IID, SEX, native_id, participant_id
    genome_build: str | None
    artifact_id: str


def load_genotype_input(project: Project, source_id: str) -> GenotypeInput:
    from efgpp.data.aliases import plink_native_ids
    from efgpp.data.genotype.formats import read_samples

    _, src = project.data.get_source(source_id)
    target = project.config.defaults.target_build
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        source_art = store.latest(source_id=source_id, artifact_type="source")
        if source_art is None:
            raise RuntimeError(f"{source_id} is not registered; run `efgpp data prepare`")
        lifted = store.latest(source_id=source_id, artifact_type="liftover_genotype")
        art = lifted if lifted is not None and lifted.genome_build == target else source_art
        if art.genome_build not in (None, target):
            raise RuntimeError(f"{source_id} is {art.genome_build}, not the target build {target}; "
                               f"run the harmonize.{source_id} step first")
        aliases = reg.frame(
            "SELECT native_id, participant_id FROM sample_aliases WHERE source_id = ?", [source_id]
        )
    fs = resolve_fileset(Path(art.path), art.format)
    samples = read_samples(fs)
    sid = src.sample_id  # type: ignore[attr-defined]
    samples = samples.with_columns(plink_native_ids(samples, sid.mode, sid.column, sid.separator))
    samples = samples.join(aliases, on="native_id", how="left", maintain_order="left")
    return GenotypeInput(source_id, fs, samples, art.genome_build, art.artifact_id)  # type: ignore[arg-type]


def plink_base(project: Project, fs: GenotypeFileset, threads: int) -> list[str]:
    args = [*fs.plink_input_args(), "--threads", str(threads)]
    mem = project.config.execution.hpc.default_mem_mb
    return [*args, "--memory", str(max(1000, int(mem * 0.8)))]


@dataclass
class QCRunOutputs:
    status: ArtifactStatus
    metrics: dict[str, Any]
    messages: list[str]
    sample_qc: Path
    variant_qc: Path
    kinship: Path | None
    prune_in: Path | None
    qc_prefix: Path
    runs: list[RunRecord] = field(default_factory=list)


def _maf_expr() -> pl.Expr:
    first = pl.col("ALT_FREQS").cast(pl.Utf8).str.split(",").list.first().cast(pl.Float64, strict=False)
    return pl.min_horizontal(first, 1 - first).alias("MAF")


def run_genotype_qc(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> QCRunOutputs:
    cfg: GenotypeQCConfig = project.config.genotype.qc
    gi = load_genotype_input(project, source_id)
    work = project.work_root / "genotype_qc" / source_id
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    qc_dir = project.qc_dir("genotype") / source_id
    qc_dir.mkdir(parents=True, exist_ok=True)
    derived = project.artifact_dir(Origin.DERIVED, Modality.GENOTYPE_QC, source_id)
    base = plink_base(project, gi.fileset, threads)
    runs: list[RunRecord] = []
    messages: list[str] = []

    def plink(args: list[str], check: bool = True) -> RunRecord:
        rec = run_tool(project, "plink2", [*base, *args], step_id=step_id,
                       inputs=[gi.artifact_id], check=check)
        runs.append(rec)
        return rec

    # 1. sample missingness ----------------------------------------------------
    plink(["--missing", "sample-only", "--out", str(work / "pass1")])
    smiss = read_plink_table(work / "pass1.smiss")
    samples = gi.samples.join(
        smiss.select("IID", pl.col("F_MISS").alias("sample_missingness")), on="IID", how="left"
    ).with_columns((pl.col("sample_missingness") > cfg.sample_missingness).alias("fail_missingness"))
    keep1 = samples.filter(~pl.col("fail_missingness"))
    if keep1.height == 0:
        raise RuntimeError(f"{source_id}: every sample fails the missingness threshold {cfg.sample_missingness}")
    keep1_path = write_keep_file(work / "keep1.txt", keep1)

    # 2. variant metrics on retained samples -----------------------------------
    plink(["--keep", str(keep1_path), "--missing", "variant-only", "--freq", "--hardy", "--out", str(work / "pass2")])
    vmiss = read_plink_table(work / "pass2.vmiss").select(
        pl.col("CHROM").alias("chromosome"), pl.col("ID").alias("variant_id"), pl.col("F_MISS").alias("variant_missingness")
    )
    afreq = read_plink_table(work / "pass2.afreq").select("ID", "REF", "ALT", "ALT_FREQS", _maf_expr())
    hwe_frames = [read_plink_table(p).select("ID", pl.col("P").alias("hwe_p"))
                  for p in (work / "pass2.hardy", work / "pass2.hardy.x") if p.exists()]
    hwe = pl.concat(hwe_frames) if hwe_frames else pl.DataFrame(schema={"ID": pl.Utf8, "hwe_p": pl.Float64})
    variants = (
        vmiss.join(afreq.rename({"ID": "variant_id"}), on="variant_id", how="left")
        .join(hwe.rename({"ID": "variant_id"}), on="variant_id", how="left")
        .with_columns(
            (pl.col("variant_missingness") > cfg.variant_missingness).alias("fail_missingness"),
            (pl.col("MAF").fill_null(0) < cfg.maf).alias("fail_maf"),
            (pl.col("hwe_p").is_not_null() & (pl.col("hwe_p") < cfg.hwe_p)).alias("fail_hwe"),
        )
        .with_columns((~(pl.col("fail_missingness") | pl.col("fail_maf") | pl.col("fail_hwe"))).alias("pass"))
    )
    passing_variants = variants.filter(pl.col("pass")).get_column("variant_id").to_list()
    if not passing_variants:
        messages.append("no variant passes QC")
    extract_path = write_id_list(work / "variants_pass.txt", passing_variants)

    # 3. LD pruning ---------------------------------------------------------------
    prune_in: Path | None = None
    if passing_variants:
        rec = plink(["--keep", str(keep1_path), "--extract", str(extract_path), "--autosome",
                     "--indep-pairwise", str(cfg.ld_window), str(cfg.ld_step), str(cfg.ld_r2),
                     "--out", str(work / "prune")], check=False)
        if rec.exit_code == 0 and (work / "prune.prune.in").exists():
            prune_in = derived / f"{source_id}_ld_pruned.prune.in"
            shutil.copy2(work / "prune.prune.in", prune_in)
        else:
            messages.append("LD pruning failed (see run log); heterozygosity/kinship use all passing variants")

    pruned_extract = prune_in or extract_path

    # 4. heterozygosity and kinship --------------------------------------------------
    samples = samples.with_columns(pl.lit(None, dtype=pl.Float64).alias("het_f"),
                                   pl.lit(False).alias("fail_heterozygosity"))
    rec = plink(["--keep", str(keep1_path), "--extract", str(pruned_extract), "--het", "--out", str(work / "het")], check=False)
    if rec.exit_code == 0 and (work / "het.het").exists():
        het = read_plink_table(work / "het.het").select("IID", pl.col("F").alias("_f"))
        mean, sd = het.get_column("_f").mean(), het.get_column("_f").std()
        samples = samples.join(het, on="IID", how="left").with_columns(pl.col("_f").alias("het_f")).drop("_f")
        if sd:
            samples = samples.with_columns(
                ((pl.col("het_f") - mean).abs() > cfg.heterozygosity_sd * sd).fill_null(False).alias("fail_heterozygosity")
            )
    else:
        messages.append("heterozygosity could not be computed")

    kin_path: Path | None = None
    pairs = pl.DataFrame()
    if project.config.genotype.relatedness.enabled and keep1.height >= 2:
        rec = plink(["--keep", str(keep1_path), "--extract", str(pruned_extract), "--make-king-table",
                     "--king-table-filter", str(min(cfg.kinship_cutoff, 0.0442)), "--out", str(work / "king")], check=False)
        if rec.exit_code == 0 and (work / "king.kin0").exists():
            pairs = classify_pairs(read_plink_table(work / "king.kin0"))
        else:
            messages.append("KING kinship could not be computed")
    remove_dup, remove_rel = set(), set()
    if pairs.height:
        miss = dict(samples.select("IID", "sample_missingness").iter_rows())
        dup_pairs = pairs.filter(pl.col("KINSHIP") > cfg.duplicate_kinship)
        remove_dup = greedy_unrelated_removal(dup_pairs, miss)
        rel_pairs = pairs.filter(pl.col("KINSHIP") > cfg.kinship_cutoff)
        remove_rel = greedy_unrelated_removal(rel_pairs, miss)
        kin_dir = project.artifact_dir(Origin.DERIVED, Modality.KINSHIP, source_id)
        kin_path = write_parquet(pairs, kin_dir / f"{source_id}_kinship.parquet")
    samples = samples.with_columns(
        pl.col("IID").is_in(list(remove_dup)).alias("duplicate_removed"),
        pl.col("IID").is_in(list(remove_rel)).alias("related_exclusion_suggested"),
    )

    # 5. sex check -----------------------------------------------------------------------
    samples = samples.with_columns(pl.lit(None, dtype=pl.Utf8).alias("sex_check"), pl.lit(False).alias("fail_sex"))
    if cfg.sex_check:
        rec = plink(["--keep", str(keep1_path), "--check-sex", "--out", str(work / "sex")], check=False)
        sex_file = work / "sex.sexcheck"
        if rec.exit_code == 0 and sex_file.exists():
            sx = read_plink_table(sex_file).select("IID", pl.col("STATUS").alias("_st"))
            samples = samples.join(sx, on="IID", how="left").with_columns(
                pl.col("_st").alias("sex_check"), (pl.col("_st") == "PROBLEM").fill_null(False).alias("fail_sex")
            ).drop("_st")
        else:
            messages.append("sex check not performed (needs chrX data with sex information and a PLINK 2 build supporting --check-sex)")

    # 6. final decisions and QC-passed genotype ---------------------------------------------
    fail_cols = ["fail_missingness", "fail_heterozygosity", "fail_sex", "duplicate_removed"]
    if cfg.remove_related:
        fail_cols.append("related_exclusion_suggested")
    samples = samples.with_columns(
        (~pl.any_horizontal([pl.col(c).fill_null(False) for c in fail_cols])).alias("pass"),
        pl.concat_str(
            [pl.when(pl.col(c).fill_null(False)).then(pl.lit(c.removeprefix("fail_"))) for c in fail_cols],
            separator=";", ignore_nulls=True,
        ).alias("fail_reasons"),
    )
    final_keep = samples.filter(pl.col("pass"))
    qc_prefix = derived / f"{source_id}_qc"
    if final_keep.height and passing_variants:
        keep_path = write_keep_file(work / "keep_final.txt", final_keep)
        plink(["--keep", str(keep_path), "--extract", str(extract_path), "--make-pgen", "--out", str(qc_prefix)])

    sample_qc = write_parquet(samples, qc_dir / "sample_qc.parquet")
    variant_qc = write_parquet(variants, qc_dir / "variant_qc.parquet")
    write_multiqc_table(samples, qc_dir / f"{source_id}_sample_qc_mqc.tsv", source_id)

    n = samples.height
    n_fail = n - final_keep.height
    frac = n_fail / n if n else 1.0
    metrics = {
        "samples_total": n, "samples_pass": final_keep.height, "samples_fail": n_fail,
        "sample_fail_fraction": round(frac, 4),
        "samples_fail_missingness": int(samples.get_column("fail_missingness").sum()),
        "samples_fail_heterozygosity": int(samples.get_column("fail_heterozygosity").sum()),
        "samples_fail_sex": int(samples.get_column("fail_sex").sum()),
        "duplicates_removed": len(remove_dup),
        "related_exclusion_suggested": len(remove_rel),
        "related_pairs": int(pairs.filter(pl.col("KINSHIP") > cfg.kinship_cutoff).height) if pairs.height else 0,
        "variants_total": variants.height, "variants_pass": len(passing_variants),
        "variants_fail_missingness": int(variants.get_column("fail_missingness").sum()),
        "variants_fail_maf": int(variants.get_column("fail_maf").sum()),
        "variants_fail_hwe": int(variants.get_column("fail_hwe").sum()),
        "ld_pruned_variants": _count_lines(prune_in) if prune_in else None,
    }
    if not final_keep.height or not passing_variants or frac > cfg.fail_sample_fail_fraction:
        status = ArtifactStatus.QC_FAIL
    elif frac > cfg.warn_sample_fail_fraction or messages:
        status = ArtifactStatus.QC_WARN
    else:
        status = ArtifactStatus.QC_PASS
    return QCRunOutputs(status, metrics, messages, sample_qc, variant_qc, kin_path, prune_in, qc_prefix, runs)


def _count_lines(path: Path) -> int:
    with open(path, "rb") as fh:
        return sum(1 for ln in fh if ln.strip())


def write_multiqc_table(samples: pl.DataFrame, path: Path, source_id: str) -> Path:
    """MultiQC custom-content table (picked up automatically by `multiqc qc/`)."""
    header = (
        f"# id: 'efgpp_{source_id}_sample_qc'\n# section_name: 'EFGPP genotype sample QC ({source_id})'\n"
        "# plot_type: 'table'\n"
    )
    cols = ["IID", "sample_missingness", "het_f", "sex_check", "pass", "fail_reasons"]
    body = samples.select([c for c in cols if c in samples.columns]).write_csv(separator="\t")
    path.write_text(header + body, encoding="utf-8")
    return path


def register_qc_outputs(project: Project, source_id: str, out: QCRunOutputs, step_id: str | None) -> dict[str, str]:
    """Register every QC output as a DERIVED artifact with lineage and provenance."""
    ids: dict[str, str] = {}
    last_run = out.runs[-1] if out.runs else None
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        src = store.latest(source_id=source_id, artifact_type="source")
        lifted = store.latest(source_id=source_id, artifact_type="liftover_genotype")
        if lifted is not None and lifted.genome_build == project.config.defaults.target_build:
            src = lifted  # QC ran on the genotype lifted to the target build
        parent = [src.artifact_id] if src and src.artifact_id else []
        build = src.genome_build if src else None
        common: dict[str, Any] = dict(origin=Origin.DERIVED, source_id=source_id, genome_build=build, tool="plink2",
                      tool_version=last_run.tool_version if last_run else None,
                      configuration_hash=project.config.genotype.qc.checksum(),
                      environment={"runs": [r.run_id for r in out.runs]})

        def replace(art: Artifact) -> Artifact:
            olds = [a for a in store.find(source_id=source_id, artifact_type=art.artifact_type)
                    if a.artifact_name == art.artifact_name]
            store.register(art)
            for o in olds:
                store.supersede(o.artifact_id, art.artifact_id)  # type: ignore[arg-type]
            return art

        for name, path, atype, modality in (
            ("sample_qc", out.sample_qc, "qc_table", Modality.GENOTYPE_QC),
            ("variant_qc", out.variant_qc, "qc_table", Modality.GENOTYPE_QC),
        ):
            a = replace(Artifact(artifact_name=f"{source_id}_{name}", artifact_type=atype, modality=modality,
                                 status=ArtifactStatus.READY, path=str(path), format="parquet",
                                 size=path.stat().st_size, checksum=str(checksum_paths([path])),
                                 parent_artifact_ids=parent, **common))
            ids[name] = a.artifact_id  # type: ignore[assignment]
        if out.kinship is not None:
            a = replace(Artifact(artifact_name=f"{source_id}_kinship", artifact_type="kinship",
                                 modality=Modality.KINSHIP, status=ArtifactStatus.READY, path=str(out.kinship),
                                 format="parquet", size=out.kinship.stat().st_size,
                                 checksum=str(checksum_paths([out.kinship])), parent_artifact_ids=parent, **common))
            ids["kinship"] = a.artifact_id  # type: ignore[assignment]
        if out.prune_in is not None:
            a = replace(Artifact(artifact_name=f"{source_id}_ld_pruned", artifact_type="ld_prune",
                                 modality=Modality.GENOTYPE_QC, status=ArtifactStatus.READY, path=str(out.prune_in),
                                 format="text", size=out.prune_in.stat().st_size,
                                 checksum=str(checksum_paths([out.prune_in])), feature_count=out.metrics["ld_pruned_variants"],
                                 parent_artifact_ids=parent, **common))
            ids["ld_prune"] = a.artifact_id  # type: ignore[assignment]

        pgen = out.qc_prefix.with_name(out.qc_prefix.name + ".pgen")
        if pgen.exists():
            fs = resolve_fileset(out.qc_prefix, "pgen")
            status = ArtifactStatus.READY if out.status != ArtifactStatus.QC_FAIL else ArtifactStatus.QC_FAIL
            a = replace(Artifact(
                artifact_name=f"{source_id}_qc", artifact_type="qc_genotype", modality=Modality.GENOTYPE_QC,
                status=status, path=str(out.qc_prefix), format="pgen", size=sum(m.stat().st_size for m in fs.members),
                checksum=str(checksum_paths(fs.members)), participant_count=out.metrics["samples_pass"],
                feature_count=count_variants(fs), temporal_type=TemporalType.STATIC,
                parent_artifact_ids=parent + [ids["sample_qc"], ids["variant_qc"]],
                command=" ".join(out.runs[-1].command) if out.runs else None,
                metadata={"qc_status": out.status.value, "qc_metrics": out.metrics, "messages": out.messages},
                **common,
            ))
            ids["qc_genotype"] = a.artifact_id  # type: ignore[assignment]
            # Availability of QC-passed genotype per participant.
            passed = pl.read_parquet(out.sample_qc).filter(pl.col("pass") & pl.col("participant_id").is_not_null())
            reg.execute("DELETE FROM assays WHERE source_id = ?", [f"{source_id}_QC"])
            reg.insert_frame("assays", pl.DataFrame({
                "assay_id": [f"{a.artifact_id}:{i}" for i in range(passed.height)],
                "artifact_id": [a.artifact_id] * passed.height,
                "source_id": [f"{source_id}_QC"] * passed.height,
                "modality": [Modality.GENOTYPE_QC.value] * passed.height,
                "origin": [Origin.DERIVED.value] * passed.height,
                "participant_id": passed.get_column("participant_id"),
            }))

        qc_run_id = reg.next_id("QC")
        reg.upsert("qc_runs", {
            "qc_run_id": qc_run_id, "artifact_id": ids.get("qc_genotype") or (src.artifact_id if src else None),
            "modality": "genotype", "run_id": out.runs[-1].run_id if out.runs else None,
            "status": out.status.value, "created_at": utcnow(),
            "config": project.config.genotype.qc.model_dump(), "summary": out.metrics,
        })
        for k, v in out.metrics.items():
            reg.execute("INSERT INTO qc_metrics VALUES (?, ?, 'dataset', ?, ?, ?, NULL, ?)", [
                qc_run_id, ids.get("qc_genotype"), k, float(v) if isinstance(v, int | float) else None,
                None if isinstance(v, int | float) else to_json(v), out.status.value,
            ])
        ids["qc_run"] = qc_run_id
    return ids


def placed_variants(variants: pl.DataFrame) -> tuple[pl.DataFrame, int]:
    """Drop unplaced variants (position 0, e.g. array control probes) from a variant table.
    They stay in the genotype files; they just cannot be annotated, lifted or keyed."""
    keep = variants.filter(pl.col("position").fill_null(0) >= 1)
    return keep, variants.height - keep.height


def variant_table(fs: GenotypeFileset, build: str | None, limit: int | None = None) -> pl.DataFrame:
    """Standardized variant table (section 27 keys; VRS ids are filled by annotation)."""
    v = read_variants(fs, limit=limit)
    return v.with_columns(
        pl.lit(build, dtype=pl.Utf8).alias("genome_build"),
        pl.when(pl.col("variant_id").str.starts_with("rs")).then(pl.col("variant_id")).otherwise(None).alias("rsid"),
        pl.lit(None, dtype=pl.Utf8).alias("vrs_id"),
    ).select("genome_build", "chromosome", "position", "reference", "alternate", "variant_id", "rsid", "vrs_id")
