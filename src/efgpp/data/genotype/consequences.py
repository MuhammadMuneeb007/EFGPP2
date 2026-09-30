"""Participant-specific variant consequences and burdens (all DERIVED from the genotype).

Variants are annotated once (unique cohort variants -> VEP, ClinVar, gnomAD, AlphaMissense,
SpliceAI); the annotation table is then joined to the participant carrier table:

  variant_annotation_table.parquet        one row per variant (VEP PICK / canonical primary
                                          consequence + every enabled annotation)
  participant_annotated_variants/         carrier rows + annotation (partitioned by chromosome)
  participant_consequence_counts.parquet  one row per participant: n_<consequence> (variants
                                          carried) and alt_burden_<consequence> (ALT alleles)
  participant_gene_burden.parquet         long: participant x gene counts and maxima

Counting rules
  * one consequence per variant: VEP's PICK row, else the canonical transcript, else the most
    severe consequence (Ensembl order), so a variant overlapping five transcripts counts once;
    the full transcript-level annotation stays in the VEP annotation table
  * a variant is carried when its hard-call ALT count is > 0 (dosage >= 0.5 without a hard call)
  * n_<term> counts carried variants, alt_burden_<term> sums ALT dosage (a homozygous site = 2)
  * a variant with several SO terms (missense_variant&splice_region_variant) counts for each term
  * burdens are EFGPP-derived counts, not a published burden model; no phenotype is read
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import duckdb
import polars as pl

from efgpp.constants import ArtifactStatus, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.io import write_parquet
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths

# Ensembl variant consequences, most severe first.
SEVERITY = [
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
    "stop_lost", "start_lost", "transcript_amplification", "feature_elongation", "feature_truncation",
    "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant",
    "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant", "coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification",
    "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification",
    "regulatory_region_variant", "intergenic_variant", "sequence_variant",
]
RANK = {t: i for i, t in enumerate(SEVERITY)}
# Columns always present in the per-participant counts, even when nobody carries such a variant.
REPORTED_TERMS = [
    "missense_variant", "synonymous_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
    "inframe_insertion", "inframe_deletion", "splice_acceptor_variant", "splice_donor_variant",
    "splice_region_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
    "intergenic_variant", "regulatory_region_variant",
]
SPLICE_TERMS = ["splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
                "splice_donor_5th_base_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant"]
PATHOGENIC = ("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")


def split_terms(consequence: str | None) -> list[str]:
    if not consequence:
        return []
    return [t for t in consequence.replace(",", "&").split("&") if t]


def most_severe(consequence: str | None) -> str | None:
    terms = split_terms(consequence)
    return min(terms, key=lambda t: RANK.get(t, len(RANK))) if terms else None


def pick_vep(vep: pl.DataFrame, key_map: pl.DataFrame | None = None) -> pl.DataFrame:
    """One row per variant: VEP PICK row, else canonical transcript, else most severe consequence."""
    df = vep
    if "variant_key" not in df.columns:
        uploaded = df.get_column("Uploaded_variation")
        if uploaded.drop_nulls().str.count_matches(":").min() == 3:  # sites VCF ID = variant key
            df = df.with_columns(pl.col("Uploaded_variation").alias("variant_key"))
        elif key_map is not None:
            df = df.join(key_map.select(pl.col("variant_id").alias("Uploaded_variation"), "variant_key"),
                         on="Uploaded_variation", how="left", maintain_order="left")
        else:
            raise ValueError("VEP table has no variant keys; re-run annotation")

    def col(name: str) -> pl.Expr:
        return pl.col(name) if name in df.columns else pl.lit(None, dtype=pl.Utf8)

    df = df.with_columns(
        col("Consequence").str.replace_all(",", "&").alias("_cons"),
        (col("PICK") == "1").fill_null(False).alias("_pick"),
        (col("CANONICAL") == "YES").fill_null(False).alias("_canonical"),
    ).with_columns(
        pl.col("_cons").map_elements(most_severe, return_dtype=pl.Utf8).alias("most_severe_consequence"),
    ).with_columns(
        pl.col("most_severe_consequence").replace_strict(RANK, default=len(RANK), return_dtype=pl.Int64)
        .alias("_rank"),
    )
    picked = (df.sort(["variant_key", "_pick", "_canonical", "_rank", "Feature"],
                      descending=[False, True, True, False, False], nulls_last=True)
              .unique(subset="variant_key", keep="first", maintain_order=True))
    return picked.select(
        "variant_key",
        col("Gene").alias("gene_id"), col("SYMBOL").alias("gene_symbol"), col("Feature").alias("transcript_id"),
        col("Feature_type").alias("feature_type"), col("BIOTYPE").alias("biotype"),
        pl.col("_cons").alias("consequence"), "most_severe_consequence", col("IMPACT").alias("impact"),
        col("HGVSc").alias("HGVSc"), col("HGVSp").alias("HGVSp"), col("VARIANT_CLASS").alias("variant_class"),
        pl.col("_canonical").alias("canonical"), col("Existing_variation").alias("existing_variation"),
        pl.when(pl.col("_pick")).then(pl.lit("PICK")).when(pl.col("_canonical")).then(pl.lit("canonical"))
        .otherwise(pl.lit("most_severe")).alias("picked_by"),
    )


def _latest_annotation(project: Project, source_id: str, name: str) -> Artifact | None:
    with Registry.open(project) as reg:
        return ArtifactStore(reg).latest(source_id=source_id, artifact_type=f"annotation_{name}")


def variant_annotation_table(project: Project, source_id: str) -> tuple[pl.DataFrame, list[Artifact]]:
    """One row per cohort variant with every available annotation (keyed by variant_key)."""
    from efgpp.data.annotation import cohort_variants, with_variant_key

    variants, _, _ = cohort_variants(project, source_id)
    base = with_variant_key(variants).select("variant_key", "variant_id", "rsid", "chromosome", "position",
                                             "reference", "alternate").unique("variant_key", maintain_order=True)
    used: list[Artifact] = []
    vep_art = _latest_annotation(project, source_id, "vep")
    if vep_art is not None:
        used.append(vep_art)
        base = base.join(pick_vep(pl.read_parquet(vep_art.path), base), on="variant_key", how="left", maintain_order="left")
    else:
        base = base.with_columns([pl.lit(None, dtype=pl.Utf8).alias(c) for c in (
            "gene_id", "gene_symbol", "transcript_id", "consequence", "most_severe_consequence", "impact",
            "HGVSc", "HGVSp")])
    cv = _latest_annotation(project, source_id, "clinvar")
    if cv is not None:
        used.append(cv)
        t = with_variant_key(pl.read_parquet(cv.path))
        base = base.join(t.select(
            "variant_key", pl.col("clinvar_variation_id"), pl.col("clnsig").alias("clinvar_significance"),
            pl.col("clnrevstat").alias("clinvar_review_status"), pl.col("clndn").alias("clinvar_condition"),
        ).unique("variant_key"), on="variant_key", how="left", maintain_order="left")
    gn = _latest_annotation(project, source_id, "gnomad")
    if gn is not None:
        used.append(gn)
        t = with_variant_key(pl.read_parquet(gn.path))
        af_cols = [c for c in t.columns if c == "af" or c.startswith("af_")]
        base = base.join(t.select("variant_key", *[pl.col(c).cast(pl.Float64, strict=False).alias(f"gnomad_{c}")
                                                   for c in af_cols]).unique("variant_key"),
                         on="variant_key", how="left", maintain_order="left")
    am = _latest_annotation(project, source_id, "alphamissense")
    if am is not None:
        used.append(am)
        t = with_variant_key(pl.read_parquet(am.path))
        if "transcript_id" in base.columns and "transcript_id" in t.columns:
            picked_tx = base.select("variant_key", pl.col("transcript_id").str.replace(r"\.\d+$", "").alias("_tx"))
            t = t.join(picked_tx, on="variant_key", how="left", maintain_order="left").with_columns(
                (pl.col("transcript_id").str.replace(r"\.\d+$", "") == pl.col("_tx")).fill_null(False).alias("_same"))
        else:
            t = t.with_columns(pl.lit(False).alias("_same"))
        t = t.sort(["variant_key", "_same", "alphamissense_score"], descending=[False, True, True], nulls_last=True)
        base = base.join(t.unique("variant_key", keep="first").select(
            "variant_key", "alphamissense_score", "alphamissense_class"), on="variant_key", how="left", maintain_order="left")
    sa = _latest_annotation(project, source_id, "spliceai")
    if sa is not None:
        used.append(sa)
        t = pl.read_parquet(sa.path).sort(["variant_key", "spliceai_max_score"], descending=[False, True],
                                          nulls_last=True)
        base = base.join(t.unique("variant_key", keep="first").select(
            "variant_key", "spliceai_gene", "spliceai_max_score", "spliceai_effect",
            "DS_AG", "DS_AL", "DS_DG", "DS_DL"), on="variant_key", how="left", maintain_order="left")
    for c, dtype in (("gnomad_af", pl.Float64), ("alphamissense_score", pl.Float64), ("spliceai_max_score", pl.Float64),
                     ("clinvar_significance", pl.Utf8), ("alphamissense_class", pl.Utf8)):
        if c not in base.columns:
            base = base.with_columns(pl.lit(None, dtype=dtype).alias(c))
    return base, used


def _participants(project: Project, source_id: str) -> list[str]:
    with Registry.open(project) as reg:
        return reg.frame("SELECT DISTINCT participant_id FROM assays WHERE source_id = ? ORDER BY 1",
                         [f"{source_id}_{Modality.PARTICIPANT_VARIANTS.value}"]).get_column("participant_id").to_list()


def aggregate(carrier_root: Path, annotations: pl.DataFrame, out_dir: Path, participants: list[str],
              cfg: Any) -> dict[str, Path]:
    """Join carriers to annotations and write the three participant-level tables."""
    out_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(carrier_root.rglob("*.parquet"))
    con = duckdb.connect()
    try:
        con.register("ann_arrow", annotations.to_arrow())
        con.execute("CREATE TABLE ann AS SELECT * FROM ann_arrow")
        con.execute("CREATE TABLE pw (term VARCHAR, weight DOUBLE)")
        if cfg.burden_weights:
            con.executemany("INSERT INTO pw VALUES (?, ?)", [(k, float(v)) for k, v in cfg.burden_weights.items()])
        if files:
            listing = "[" + ", ".join(f"'{f.as_posix()}'" for f in files) + "]"
            con.execute(f"CREATE VIEW carriers AS SELECT * FROM read_parquet({listing}, hive_partitioning = false)")
        else:
            con.execute("CREATE VIEW carriers AS SELECT NULL::VARCHAR AS participant_id, NULL::VARCHAR AS variant_key, "
                        "NULL::TINYINT AS hardcall_alt_count, NULL::FLOAT AS dosage, NULL::BOOLEAN AS heterozygous, "
                        "NULL::BOOLEAN AS homozygous_alt, NULL::VARCHAR AS chromosome, NULL::BIGINT AS position, "
                        "NULL::VARCHAR AS ref, NULL::VARCHAR AS alt WHERE false")
        lof = ", ".join(f"'{t}'" for t in cfg.lof_consequences)
        splice = ", ".join(f"'{t}'" for t in SPLICE_TERMS)
        patho = " OR ".join(f"a.clinvar_significance LIKE '%{p}%'" for p in PATHOGENIC)
        con.execute(f"""
            CREATE TABLE carried AS
            SELECT c.participant_id, c.variant_key, c.chromosome, c.position, c.ref, c.alt,
                   c.hardcall_alt_count AS alt_count, coalesce(c.dosage, c.hardcall_alt_count) AS dosage,
                   c.heterozygous, c.homozygous_alt,
                   a.gene_id, a.gene_symbol, a.transcript_id, a.consequence, a.most_severe_consequence, a.impact,
                   a.HGVSc, a.HGVSp, a.clinvar_significance, a.gnomad_af, a.alphamissense_score,
                   a.alphamissense_class, a.spliceai_max_score,
                   list_has_any(string_split(coalesce(a.consequence, ''), '&'), [{lof}]) AS is_lof,
                   list_has_any(string_split(coalesce(a.consequence, ''), '&'), [{splice}]) AS is_splice,
                   (list_contains(string_split(coalesce(a.consequence, ''), '&'), 'missense_variant')
                    AND a.alphamissense_score >= {float(cfg.damaging_alphamissense)}) AS is_damaging_missense,
                   coalesce(({patho}) AND a.clinvar_significance NOT LIKE '%Conflicting%', false) AS is_pathogenic,
                   a.gnomad_af < {float(cfg.rare_af)} AS is_rare,
                   coalesce(pw.weight, 0) AS burden_weight
            FROM carriers c LEFT JOIN ann a USING (variant_key)
            LEFT JOIN pw ON pw.term = a.most_severe_consequence
            WHERE coalesce(c.hardcall_alt_count, 0) > 0
               OR (c.hardcall_alt_count IS NULL AND c.dosage >= 0.5)
        """)
        annotated = out_dir / "participant_annotated_variants"
        if annotated.exists():
            import shutil

            shutil.rmtree(annotated)
        annotated.mkdir(parents=True)
        con.execute(f"COPY (SELECT * FROM carried ORDER BY chromosome, position, participant_id) TO "
                    f"'{annotated.as_posix()}' (FORMAT PARQUET, PARTITION_BY (chromosome), COMPRESSION ZSTD, "
                    f"OVERWRITE_OR_IGNORE)")
        term_long = con.execute("""
            SELECT participant_id, term, count(*) AS n, sum(dosage) AS alt_burden
            FROM (SELECT participant_id, dosage, unnest(string_split(consequence, '&')) AS term
                  FROM carried WHERE consequence IS NOT NULL)
            GROUP BY 1, 2
        """).pl()
        thresholds = sorted(cfg.spliceai_thresholds)
        splice_cols = ", ".join(
            f"count(*) FILTER (WHERE spliceai_max_score >= {t}) AS n_spliceai_ge_{str(t).replace('.', '_')}"
            for t in thresholds)
        totals = con.execute(f"""
            SELECT participant_id,
                   count(*) AS n_variants_total, sum(dosage) AS alt_alleles_total,
                   count(*) FILTER (WHERE heterozygous) AS n_heterozygous,
                   count(*) FILTER (WHERE homozygous_alt) AS n_homozygous_alt,
                   count(*) FILTER (WHERE consequence IS NULL) AS n_unannotated,
                   count(*) FILTER (WHERE is_lof) AS n_lof, coalesce(sum(dosage) FILTER (WHERE is_lof), 0) AS alt_burden_lof,
                   count(*) FILTER (WHERE is_splice) AS n_splice,
                   count(*) FILTER (WHERE is_damaging_missense) AS n_damaging_missense,
                   coalesce(sum(dosage) FILTER (WHERE is_damaging_missense), 0) AS alt_burden_damaging_missense,
                   count(*) FILTER (WHERE is_pathogenic) AS n_clinvar_pathogenic,
                   count(*) FILTER (WHERE is_rare) AS n_rare,
                   count(*) FILTER (WHERE is_rare AND (is_lof OR is_damaging_missense)) AS n_rare_damaging
                   {", " + splice_cols if splice_cols else ""}
            FROM carried GROUP BY 1
        """).pl()
        gene = con.execute("""
            SELECT participant_id, gene_id, any_value(gene_symbol) AS gene_symbol,
                   count(*) AS n_variants, sum(dosage) AS n_alt_alleles,
                   count(*) FILTER (WHERE list_contains(string_split(consequence, '&'), 'missense_variant')) AS n_missense,
                   count(*) FILTER (WHERE is_lof) AS n_lof, count(*) FILTER (WHERE is_splice) AS n_splice,
                   count(*) FILTER (WHERE is_damaging_missense) AS n_damaging_missense,
                   max(alphamissense_score) AS max_alphamissense, max(spliceai_max_score) AS max_spliceai,
                   min(gnomad_af) AS min_gnomad_af,
                   sum(dosage * burden_weight) AS efgpp_weighted_burden
            FROM carried WHERE gene_id IS NOT NULL
            GROUP BY 1, 2 ORDER BY 1, 2
        """).pl()
    finally:
        con.close()

    counts = pl.DataFrame({"participant_id": participants}, schema={"participant_id": pl.Utf8})
    counts = counts.join(totals, on="participant_id", how="left", maintain_order="left")
    terms = list(dict.fromkeys([*REPORTED_TERMS, *sorted(term_long.get_column("term").unique().to_list(),
                                                          key=lambda t: RANK.get(t, len(RANK)))]))
    if term_long.height:
        n_wide = term_long.pivot(on="term", index="participant_id", values="n")
        b_wide = term_long.pivot(on="term", index="participant_id", values="alt_burden")
        n_wide = n_wide.rename({t: f"n_{t}" for t in n_wide.columns if t != "participant_id"})
        b_wide = b_wide.rename({t: f"alt_burden_{t}" for t in b_wide.columns if t != "participant_id"})
        counts = counts.join(n_wide, on="participant_id", how="left", maintain_order="left").join(b_wide, on="participant_id", how="left", maintain_order="left")
    for t in terms:
        for prefix in ("n_", "alt_burden_"):
            if f"{prefix}{t}" not in counts.columns:
                counts = counts.with_columns(pl.lit(0).alias(f"{prefix}{t}"))
    fixed = [c for c in totals.columns if c != "participant_id"]
    ordered = ["participant_id", *fixed, *[f"n_{t}" for t in terms], *[f"alt_burden_{t}" for t in terms]]
    counts = counts.select(ordered).with_columns(
        [pl.col(c).fill_null(0) for c in ordered if c != "participant_id"])
    if not cfg.burden_weights:
        gene = gene.drop("efgpp_weighted_burden")
    return {
        "annotated": annotated,
        "counts": write_parquet(counts, out_dir / "participant_consequence_counts.parquet"),
        "gene_burden": write_parquet(gene, out_dir / "participant_gene_burden.parquet"),
    }


def run_consequences(project: Project, source_id: str, *, step_id: str | None = None,
                     threads: int = 1) -> dict[str, Any]:
    """Participant consequence counts and gene burden for one genotype source."""
    from efgpp.data.genotype.carriers import carrier_dataset, write_assays

    cd = carrier_dataset(project, source_id)
    if cd is None:
        raise RuntimeError(f"{source_id}: no participant carrier table (run the participant_variants step)")
    carriers_art, carrier_root, _sites = cd
    cfg = project.data.participant_variants
    annotations, used = variant_annotation_table(project, source_id)
    out_dir = project.artifact_dir(Origin.DERIVED, Modality.CONSEQUENCE_BURDEN, source_id)
    ann_path = write_parquet(annotations, out_dir / "variant_annotation_table.parquet")
    participants = _participants(project, source_id)
    outputs = aggregate(carrier_root, annotations, out_dir, participants, cfg)
    versions = {a.artifact_name: a.checksum or "" for a in used}
    parents = [carriers_art.artifact_id, *(a.artifact_id for a in used)]
    note = "EFGPP-derived counts from the genotype; not a published burden model; no phenotype used"
    ids: dict[str, Any] = {}
    counts = pl.read_parquet(outputs["counts"])
    gene_n = pl.read_parquet(outputs["gene_burden"]).get_column("gene_id").n_unique()
    specs = [
        ("participant_annotated_variants", Modality.PARTICIPANT_VARIANTS, [outputs["annotated"], ann_path], None),
        ("consequence_counts", Modality.CONSEQUENCE_BURDEN, [outputs["counts"]], counts.width - 1),
        ("gene_burden", Modality.GENE_BURDEN, [outputs["gene_burden"]], gene_n),
    ]
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        for name, modality, files, features in specs:
            a = store.register_replacing(Artifact(
                artifact_name=f"{source_id}_{name}", artifact_type=name, modality=modality, origin=Origin.DERIVED,
                status=ArtifactStatus.READY, path=str(files[0]),
                format="parquet_dataset" if files[0].is_dir() else "parquet",
                size=sum(f.stat().st_size for d in files for f in ([d] if d.is_file() else d.rglob("*")) if f.is_file()),
                checksum=str(checksum_paths(files, full_limit_bytes=_limit(project))), participant_count=len(participants), feature_count=features,
                genome_build=carriers_art.genome_build, temporal_type=TemporalType.STATIC, source_id=source_id,
                parent_artifact_ids=parents,  # type: ignore[arg-type]
                tool="efgpp", resource_versions=versions,
                metadata={"members": [str(f) for f in files], "note": note, "annotations": [a.artifact_name for a in used],
                          "counting": "VEP PICK / canonical / most severe consequence, one per variant"},
            ))
            ids[name] = a.artifact_id
            if modality in (Modality.CONSEQUENCE_BURDEN, Modality.GENE_BURDEN):
                write_assays(reg, a, f"{source_id}_{modality.value}", modality, participants)
    ids["participants"] = len(participants)
    return ids


def _limit(project: Project) -> int:
    """Same full-checksum size limit as snapshot verification (larger files: labelled sampled hash)."""
    return int(project.config.storage.full_checksum_limit_gb * 1024**3)
