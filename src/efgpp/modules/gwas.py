"""GWAS module: map the columns of a GWAS summary-statistics file to GWASLab keywords.

Used by `efgpp data add gwas --path <file>` when --col is not given. The original column
names are kept in the output (`<id>.GRCh38.parquet`); the GWASLab-standard copy
(`<id>.GRCh38.gwaslab.parquet`) uses GWASLab's names (CHR, POS, EA, NEA, ...).

Each GWASLab keyword lists the header names it accepts (case-insensitive), first match wins.
Add your consortium's spellings here.
"""

from __future__ import annotations

from efgpp.modules.columns import normalize

# GWASLab keyword -> accepted header names (lower case)
GWAS_COLUMNS: dict[str, list[str]] = {
    "snpid": ["snp", "snpid", "markername", "marker", "variant_id", "variantid", "id", "snp_id", "cptid"],
    "rsid": ["rsid", "rs_id", "rs", "rsnumber", "dbsnp_id"],
    "chrom": ["chr", "chrom", "chromosome", "#chrom", "#chr", "hg19chrc", "chr_name"],
    "pos": ["bp", "pos", "position", "base_pair_location", "bp_hg19", "bp_hg38", "genpos", "chr_position"],
    "ea": ["a1", "ea", "effect_allele", "allele1", "alt", "tested_allele", "effectallele", "inc_allele"],
    "nea": ["a2", "nea", "other_allele", "non_effect_allele", "allele2", "ref", "otherallele", "dec_allele"],
    "eaf": ["eaf", "effect_allele_frequency", "freq", "freq1", "a1freq", "af", "frq", "a1_freq", "freq_a1"],
    "beta": ["beta", "b", "effect", "estimate", "logor", "log_or"],
    "OR": ["or", "odds_ratio", "oddsratio"],
    "se": ["se", "stderr", "standard_error", "sebeta", "se_beta"],
    "z": ["z", "zscore", "z_score", "zstat"],
    "p": ["p", "pval", "p_value", "pvalue", "p-value", "p.value", "p_bolt_lmm", "p_bolt_lmm_inf"],
    "mlog10p": ["mlog10p", "neg_log_10_p_value", "log10p", "-log10p", "lp"],
    "n": ["n", "neff", "n_total", "total_n", "samplesize", "sample_size", "nobs", "obs_ct"],
    "ncase": ["n_case", "ncase", "n_cases", "ncases"],
    "ncontrol": ["n_control", "ncontrol", "n_controls", "ncontrols"],
    "info": ["info", "imputation_info", "info_score", "r2", "rsq"],
}


def infer_gwas_columns(header: list[str]) -> dict[str, str]:
    """GWASLab keyword -> original column name, for every column that can be recognized."""
    lookup: dict[str, str] = {}
    for col in header:
        lookup.setdefault(normalize(col), col)
    mapping: dict[str, str] = {}
    used: set[str] = set()
    for keyword, names in GWAS_COLUMNS.items():
        for name in names:
            original = lookup.get(name)
            if original is not None and original not in used:
                mapping[keyword] = original
                used.add(original)
                break
    return mapping


def missing_essentials(mapping: dict[str, str]) -> list[str]:
    """What GWASLab needs at least: position (CHR+POS) or a variant ID, both alleles, one
    statistic and a p-value (or something to derive it from)."""
    missing = []
    if not ({"chrom", "pos"} <= mapping.keys() or {"snpid"} & mapping.keys() or {"rsid"} & mapping.keys()):
        missing.append("chrom+pos or a variant ID")
    if not {"ea", "nea"} <= mapping.keys():
        missing.append("effect and other allele")
    if not ({"beta", "OR", "z"} & mapping.keys()):
        missing.append("beta, OR or z")
    if not ({"p", "mlog10p", "z"} & mapping.keys() or {"beta", "se"} <= mapping.keys()):
        missing.append("p-value")
    return missing
