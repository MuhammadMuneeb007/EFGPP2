"""Harmonize standardized model weights with a genotype's variants (every generic model passes here).

For each model variant the genotype site is found by chromosome + position (or by rsID when the
model has no position) and the alleles are compared:

  exact              effect = genotype ALT, other = REF        score counts the ALT allele
  allele_swap        effect = genotype REF, other = ALT        score counts the REF allele
  strand_flip        complement(effect) = ALT                  score counts the ALT allele
  strand_flip_swap   complement(effect) = REF                  score counts the REF allele

A/T and C/G variants are never strand-flipped. They are kept only when the model's alleles are
known to be on the forward strand of the reference (PredictDB GTEx v8 ids chr_pos_REF_ALT);
otherwise they are excluded and reported. Several genotype sites matching one model variant
(multiallelic duplicates), duplicate model rows, missing variants and allele mismatches are
counted per feature. PLINK scores the named allele, so no weight is ever sign-flipped.
"""

from __future__ import annotations

from dataclasses import dataclass

import polars as pl

from efgpp.constants import GenomeBuild

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
QC_COLUMNS = ["n_model_variants", "n_matched_variants", "n_missing_variants", "n_exact_matches", "n_allele_swaps",
              "n_strand_flips", "n_ambiguous_excluded", "n_multiallelic_excluded", "n_allele_mismatch",
              "n_duplicate_excluded"]


class BuildMismatchError(RuntimeError):
    pass


@dataclass
class Harmonized:
    matched: pl.DataFrame  # model_id, feature_id, geno_variant_id, counted_allele, counts_alt, weight, match_type
    qc: pl.DataFrame  # one row per model_id/feature_id
    excluded: pl.DataFrame  # model variants not used, with the reason


def check_build(model_build: str | None, genotype_build: str | None) -> None:
    m, g = GenomeBuild.normalize(model_build), GenomeBuild.normalize(genotype_build)
    if m in (GenomeBuild.UNKNOWN, GenomeBuild.AUTO):
        raise BuildMismatchError("model genome build is unknown; EFGPP never assumes a build "
                                 "(record it in the model resource manifest)")
    if g in (GenomeBuild.UNKNOWN, GenomeBuild.AUTO):
        raise BuildMismatchError("genotype genome build is unknown; run `efgpp data prepare` (build detection)")
    if m != g:
        raise BuildMismatchError(
            f"model build {m.value} != genotype build {g.value}; create a matching genotype explicitly "
            f"(`efgpp data liftover <artifact> --to {m.value}`) or install the {g.value} version of the model. "
            "EFGPP never converts silently.")


def _comp(allele: pl.Expr) -> pl.Expr:
    return allele.str.replace_all("A", "t").str.replace_all("T", "a").str.replace_all("C", "g") \
        .str.replace_all("G", "c").str.to_uppercase()


def harmonize(weights: pl.DataFrame, sites: pl.DataFrame, *, genotype_build: str | None,
              model_build: str | None, forward_strand: bool = False, match_rsid: bool = True) -> Harmonized:
    """weights: standard schema; sites: genotype variants (variant_id, chromosome, position, ref, alt)."""
    check_build(model_build, genotype_build)
    w = weights.with_row_index("_row").with_columns(
        pl.col("effect_allele").str.to_uppercase(), pl.col("other_allele").str.to_uppercase(),
        pl.col("chromosome").cast(pl.Utf8).str.replace(r"^(?i)chr", ""))
    first = w.unique(["model_id", "feature_id", "chromosome", "position", "variant_id", "effect_allele",
                      "other_allele"], keep="first", maintain_order=True)
    duplicates = w.join(first.select("_row"), on="_row", how="anti")
    w = first
    rsid = pl.col("rsid") if "rsid" in sites.columns else pl.col("variant_id")
    s = sites.select(pl.col("variant_id").alias("geno_variant_id"), rsid.alias("_rsid"),
                     pl.col("chromosome").cast(pl.Utf8).str.replace(r"^(?i)chr", "").alias("chromosome"),
                     pl.col("position").cast(pl.Int64), pl.col("ref").str.to_uppercase().alias("g_ref"),
                     pl.col("alt").str.to_uppercase().alias("g_alt"))
    by_pos = w.filter(pl.col("position").is_not_null()).join(s.drop("_rsid"), on=["chromosome", "position"],
                                                             how="inner")
    cand = by_pos
    if match_rsid:  # models without positions (rsID only)
        no_pos = w.filter(pl.col("position").is_null() & pl.col("variant_id").str.starts_with("rs"))
        by_id = no_pos.drop("chromosome", "position").join(
            s.filter(pl.col("_rsid").str.starts_with("rs")), left_on="variant_id", right_on="_rsid", how="inner")
        cand = pl.concat([by_pos, by_id.select(by_pos.columns)], how="vertical_relaxed")
    ea, oa = pl.col("effect_allele"), pl.col("other_allele")
    ambiguous = (pl.concat_str([ea, oa]).is_in(["AT", "TA", "CG", "GC"]))
    has_other = oa.is_not_null() & (oa != "")
    exact = (ea == pl.col("g_alt")) & (~has_other | (oa == pl.col("g_ref")))
    swap = (ea == pl.col("g_ref")) & (~has_other | (oa == pl.col("g_alt")))
    flip = (_comp(ea) == pl.col("g_alt")) & (~has_other | (_comp(oa) == pl.col("g_ref")))
    flip_swap = (_comp(ea) == pl.col("g_ref")) & (~has_other | (_comp(oa) == pl.col("g_alt")))
    cand = cand.with_columns(
        pl.when(ambiguous & ~pl.lit(forward_strand)).then(pl.lit("ambiguous"))
        .when(exact).then(pl.lit("exact"))
        .when(swap).then(pl.lit("allele_swap"))
        .when(~ambiguous & flip).then(pl.lit("strand_flip"))
        .when(~ambiguous & flip_swap).then(pl.lit("strand_flip_swap"))
        .otherwise(pl.lit("mismatch")).alias("match_type"))
    good = cand.filter(pl.col("match_type").is_in(["exact", "allele_swap", "strand_flip", "strand_flip_swap"]))
    n_good = good.group_by("_row").len("n_sites")
    good = good.join(n_good, on="_row")
    matched = good.filter(pl.col("n_sites") == 1).with_columns(
        pl.col("match_type").is_in(["exact", "strand_flip"]).alias("counts_alt"),
    ).with_columns(
        pl.when(pl.col("counts_alt")).then(pl.col("g_alt")).otherwise(pl.col("g_ref")).alias("counted_allele"))
    multi_rows = good.filter(pl.col("n_sites") > 1).select("_row").unique()
    used = set(matched.get_column("_row").to_list()) | set(multi_rows.get_column("_row").to_list())
    cand_rows = cand.group_by("_row").agg(
        (pl.col("match_type") == "ambiguous").any().alias("amb"), pl.len().alias("n"))
    reasons = w.select("_row", "model_id", "feature_id", "chromosome", "position", "variant_id", "effect_allele",
                       "other_allele").join(cand_rows, on="_row", how="left", maintain_order="left").with_columns(
        pl.when(pl.col("_row").is_in(multi_rows.get_column("_row").to_list())).then(pl.lit("multiallelic"))
        .when(pl.col("n").is_null()).then(pl.lit("missing"))
        .when(pl.col("amb")).then(pl.lit("ambiguous"))
        .otherwise(pl.lit("allele_mismatch")).alias("reason"))
    excluded = pl.concat([
        reasons.filter(~pl.col("_row").is_in(list(used)) | pl.col("_row").is_in(multi_rows.get_column("_row").to_list()))
        .drop("amb", "n"),
        duplicates.select("_row", "model_id", "feature_id", "chromosome", "position", "variant_id", "effect_allele",
                          "other_allele").with_columns(pl.lit("duplicate").alias("reason")),
    ], how="vertical_relaxed")
    matched_out = matched.select("model_id", "feature_id", "geno_variant_id", "counted_allele", "counts_alt",
                                 "weight", "match_type", "chromosome", "position")
    return Harmonized(matched_out, model_qc(weights, matched_out, excluded), excluded.drop("_row"))


def model_qc(weights: pl.DataFrame, matched: pl.DataFrame, excluded: pl.DataFrame) -> pl.DataFrame:
    keys = ["model_id", "feature_id"]
    total = weights.group_by(keys).agg(pl.len().alias("n_model_variants"),
                                       pl.col("validation_r2").first().alias("validation_r2"))
    m = matched.group_by(keys).agg(
        pl.len().alias("n_matched_variants"),
        (pl.col("match_type") == "exact").sum().alias("n_exact_matches"),
        (pl.col("match_type") == "allele_swap").sum().alias("n_allele_swaps"),
        pl.col("match_type").is_in(["strand_flip", "strand_flip_swap"]).sum().alias("n_strand_flips"))
    e = excluded.group_by(keys).agg(
        (pl.col("reason") == "missing").sum().alias("n_missing_variants"),
        (pl.col("reason") == "ambiguous").sum().alias("n_ambiguous_excluded"),
        (pl.col("reason") == "multiallelic").sum().alias("n_multiallelic_excluded"),
        (pl.col("reason") == "allele_mismatch").sum().alias("n_allele_mismatch"),
        (pl.col("reason") == "duplicate").sum().alias("n_duplicate_excluded"))
    qc = total.join(m, on=keys, how="left", maintain_order="left").join(e, on=keys, how="left", maintain_order="left").with_columns(
        [pl.col(c).fill_null(0).cast(pl.Int64) for c in QC_COLUMNS if c != "n_model_variants"])
    return qc.with_columns(
        (pl.col("n_matched_variants") / (pl.col("n_model_variants") - pl.col("n_duplicate_excluded"))
         .clip(lower_bound=1)).alias("coverage_fraction")
    ).select(*keys, *QC_COLUMNS, "coverage_fraction", "validation_r2").sort(keys)


def apply_coverage(qc: pl.DataFrame, minimum: float, allow_low: bool,
                   minimum_r2: float | None = None) -> pl.DataFrame:
    """Add status (OK | LOW_COVERAGE | NO_MATCH), below_threshold and `predict` (whether scored)."""
    r2_low = (pl.col("validation_r2") <= minimum_r2).fill_null(False) if minimum_r2 is not None else pl.lit(False)
    return qc.with_columns(
        pl.when(pl.col("n_matched_variants") == 0).then(pl.lit("NO_MATCH"))
        .when(pl.col("coverage_fraction") < minimum).then(pl.lit("LOW_COVERAGE"))
        .otherwise(pl.lit("OK")).alias("status"),
        r2_low.alias("below_threshold"),
    ).with_columns(
        ((pl.col("status") == "OK") | ((pl.col("status") == "LOW_COVERAGE") & pl.lit(allow_low))).alias("predict"))
