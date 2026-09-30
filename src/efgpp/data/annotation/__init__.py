"""Variant annotation (VEP, OpenCRAVAT) and reference-score lookups (ClinVar, gnomAD,
AlphaMissense, AlphaGenome). Annotations are phenotype-independent reference knowledge
attached to the cohort's variants; they are DERIVED artifacts with resource versions.
"""

from __future__ import annotations

import gzip
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

KEY = ["chromosome", "position", "reference", "alternate"]


def cohort_variants(project: Project, source_id: str, qc_only: bool = True) -> tuple[pl.DataFrame, Artifact, str | None]:
    """Variants of a genotype source to annotate. Returns (variants, parent artifact, genome build).

    When the participant carrier table exists its sites are used: REF/ALT verified against the
    reference FASTA and keyed chromosome:position:REF:ALT. Otherwise the standardized variant
    table (QC-passing variants when QC has run; PLINK's provisional REF = A2)."""
    from efgpp.data.genotype.carriers import carrier_dataset

    carriers = carrier_dataset(project, source_id)
    if carriers is not None and carriers[2].exists():
        art, _root, sites = carriers
        v = pl.read_parquet(sites).filter(pl.col("included"))
        return v.select("chromosome", "position", "reference", "alternate", "variant_id", "variant_key",
                        "rsid"), art, art.genome_build
    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        vt = store.latest(source_id=source_id, artifact_type="variant_table")
        vqc = next((a for a in store.find(source_id=source_id, artifact_type="qc_table")
                    if a.artifact_name.endswith("variant_qc")), None) if qc_only else None
    if vt is None:
        raise RuntimeError(f"{source_id}: no variant table; run `efgpp data prepare` first")
    variants = pl.read_parquet(vt.path)
    if vqc is not None:
        passing = pl.read_parquet(vqc.path).filter(pl.col("pass")).select("variant_id")
        variants = variants.join(passing, on="variant_id", how="semi")
    return with_variant_key(variants), vt, vt.genome_build


def key_expr(chrom: str = "chromosome", pos: str = "position", ref: str = "reference",
             alt: str = "alternate") -> pl.Expr:
    """Canonical variant key chromosome:position:REF:ALT."""
    return pl.concat_str([pl.col(chrom).cast(pl.Utf8).str.replace(r"^(?i)chr", ""), pl.col(pos).cast(pl.Utf8),
                          pl.col(ref), pl.col(alt)], separator=":")


def with_variant_key(df: pl.DataFrame) -> pl.DataFrame:
    if "variant_key" in df.columns:
        return df
    return df.with_columns(key_expr().alias("variant_key"))


def write_sites_vcf(variants: pl.DataFrame, path: Path, build: str | None) -> Path:
    """Sites-only VCF (no samples) used as annotation input."""
    path.parent.mkdir(parents=True, exist_ok=True)
    v = variants.filter(pl.col("reference").is_not_null() & pl.col("alternate").is_not_null())
    v = v.sort(["chromosome", "position"])
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        if build:
            fh.write(f"##reference={build}\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        # ID = canonical variant key, so annotations join back to participants unambiguously
        for chrom, pos, key, ref, alt in with_variant_key(v).select(
                "chromosome", "position", "variant_key", "reference", "alternate").iter_rows():
            fh.write(f"{chrom}\t{pos}\t{key}\t{ref}\t{alt}\t.\t.\t.\n")
    return path


def header_line_index(path: Path, marker: str = "#CHROM") -> int:
    """0-based index of the column-header line in a (gzipped) VCF/TSV with comment lines."""
    opener = gzip.open if path.name.endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        for i, line in enumerate(fh):
            if line.startswith(marker):
                return i
            if not line.startswith("#"):
                break
    raise ValueError(f"{path}: header line starting with {marker!r} not found")


def lookup_table(
    path: Path,
    variants: pl.DataFrame,
    *,
    select_sql: str,
    chrom_col: str = "CHROM",
    pos_col: str = "POS",
    ref_col: str = "REF",
    alt_col: str = "ALT",
    header_marker: str = "#CHROM",
) -> pl.DataFrame:
    """Join cohort variants to a large delimited reference file with DuckDB.

    `select_sql` selects extra columns from alias `r` (e.g. "r.am_pathogenicity").
    Multi-allelic ALT fields are split so each allele matches independently.
    """
    skip = header_line_index(path, header_marker)
    con = duckdb.connect()
    try:
        con.register("cohort", variants.select(KEY).unique().to_arrow())
        compression = "gzip" if path.name.endswith((".gz", ".bgz")) else "none"
        con.execute(
            f"""
            CREATE VIEW raw AS
            SELECT * FROM read_csv('{path.as_posix()}', delim='	', skip={skip}, header=true,
                                   all_varchar=true, compression='{compression}', quote='',
                                   ignore_errors=true)
            """
        )
        cols = [c[0] for c in con.execute("DESCRIBE raw").fetchall()]

        def actual(name: str) -> str:
            for c in cols:
                if c.lstrip("#").upper() == name.lstrip("#").upper():
                    return c
            raise KeyError(f"{path.name}: column {name!r} not found in {cols[:12]}")

        chrom, pos, ref, alt = (actual(c) for c in (chrom_col, pos_col, ref_col, alt_col))
        con.execute(
            f"""
            CREATE VIEW r AS
            SELECT *, upper(_alt_raw) AS _alt FROM (
                SELECT *, regexp_replace("{chrom}", '^chr', '') AS _chrom,
                       TRY_CAST("{pos}" AS BIGINT) AS _pos, upper("{ref}") AS _ref,
                       unnest(string_split("{alt}", ',')) AS _alt_raw
                FROM raw)
            """
        )
        sql = f"""
            SELECT c.chromosome, c.position, c.reference, c.alternate, {select_sql}
            FROM cohort c
            JOIN r ON r._chrom = c.chromosome AND r._pos = c.position
                  AND r._ref = c.reference AND r._alt = c.alternate
        """
        return con.execute(sql).pl()
    finally:
        con.close()


def info_field(name: str) -> str:
    """SQL extracting one INFO key from a VCF INFO column."""
    return f"regexp_extract(r.INFO, '(?:^|;){name}=([^;]*)', 1) AS {name.lower()}"


def register_annotation(
    project: Project,
    *,
    source_id: str,
    name: str,
    table: pl.DataFrame,
    parent: Artifact,
    tool: str,
    tool_version: str | None,
    resource_versions: dict[str, str],
    genome_build: str | None,
    extra_files: list[Path] | None = None,
    metadata: dict[str, Any] | None = None,
) -> str:
    out_dir = project.artifact_dir(Origin.DERIVED, Modality.VARIANT_ANNOTATIONS, source_id)
    if "variant_key" not in table.columns and set(KEY) <= set(table.columns):
        table = table.with_columns(key_expr().alias("variant_key"))
    out = write_parquet(table, out_dir / f"{source_id}_{name}.parquet", project.config.storage.compression)
    files = [out, *(extra_files or [])]
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=f"{source_id}_{name}", artifact_type=f"annotation_{name}",
            modality=Modality.VARIANT_ANNOTATIONS, origin=Origin.DERIVED, status=ArtifactStatus.READY,
            path=str(out), format="parquet", size=sum(f.stat().st_size for f in files),
            checksum=str(checksum_paths(files)), feature_count=table.height, genome_build=genome_build,
            temporal_type=TemporalType.STATIC, source_id=source_id,
            parent_artifact_ids=[parent.artifact_id],  # type: ignore[list-item]
            tool=tool, tool_version=tool_version, resource_versions=resource_versions,
            metadata={**(metadata or {}), "members": [str(f) for f in files]},
        ))
    return art.artifact_id  # type: ignore[return-value]


def resource_version(project: Project, name: str) -> tuple[str | None, Path | None]:
    """Latest installed (version, local path) of a reference resource."""
    with Registry.open(project) as reg:
        row = reg.one("SELECT version, local_path FROM resources WHERE name = ? ORDER BY download_date DESC LIMIT 1", [name])
    if row is None:
        return None, None
    return row["version"], Path(row["local_path"]) if row["local_path"] else None
