"""Generic genetic scores: prediction = sum over variants of counted-allele dosage x weight.

Backends
  plink2   `plink2 --score` (primary): dosage-aware, low memory, reads PGEN/BED/BGEN directly.
           The genotype is first given unique chromosome:position:REF:ALT variant ids once
           (cached keyed PGEN) so model variants are matched unambiguously. Features are scored
           in batches; variants whose ALT is counted and those whose REF is counted go into two
           score files, so every weight is applied to exactly the allele the model names.
  python   numpy on PLINK .bed filesets (same arithmetic; used without PLINK 2 and in tests).

Missing genotypes contribute 0 (no mean imputation). An intercept is added only when a model
supplies one. The genotype matrix is never loaded whole: only the matched variants are read.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
from scipy import sparse

from efgpp.data.genotype.carriers import decode_bed
from efgpp.data.genotype.formats import (
    GenotypeFileset,
    read_samples,
    read_variants,
    resolve_fileset,
)
from efgpp.project import Project


def feature_keys(matched: pl.DataFrame) -> list[tuple[str, str]]:
    return list(matched.select("model_id", "feature_id").unique(maintain_order=True).iter_rows())


# ------------------------------------------------------------------ python backend
def bed_sites(fs: GenotypeFileset) -> pl.DataFrame:
    """Sites of a .bed fileset with row-unique ids (#<row>); .bim A2 = ref, A1 = alt as PLINK reads it."""
    v = read_variants(fs).with_row_index("_row")
    return v.select(pl.concat_str([pl.lit("#"), pl.col("_row").cast(pl.Utf8)]).alias("variant_id"),
                    pl.col("variant_id").alias("rsid"), "chromosome", "position",
                    pl.col("reference").alias("ref"), pl.col("alternate").alias("alt"))


def score_python(fs: GenotypeFileset, matched: pl.DataFrame, intercepts: dict[tuple[str, str], float] | None = None
                 ) -> pl.DataFrame:
    """Scores from a .bed fileset for `matched` rows whose geno_variant_id is `#<bim row>`."""
    if fs.format != "bed":
        raise ValueError("the Python scorer reads PLINK .bed filesets only; install PLINK 2 for other formats")
    samples = read_samples(fs)
    n = samples.height
    keys = feature_keys(matched)
    col = {k: i for i, k in enumerate(keys)}
    rows = matched.with_columns(pl.col("geno_variant_id").str.slice(1).cast(pl.Int64).alias("_row"))
    needed = sorted(set(rows.get_column("_row").to_list()))
    pos = {r: i for i, r in enumerate(needed)}
    bpv = (n + 3) // 4
    a1 = np.zeros((len(needed), n), dtype=np.int8)
    with open(fs.prefix.with_name(fs.prefix.name + ".bed"), "rb") as fh:
        if fh.read(3) != b"\x6c\x1b\x01":
            raise ValueError("not a SNP-major .bed file")
        for r in needed:
            fh.seek(3 + r * bpv)
            a1[pos[r]] = decode_bed(np.frombuffer(fh.read(bpv), dtype=np.uint8).reshape(1, bpv), n)[0]
    missing = a1 < 0
    scores = np.zeros((n, len(keys)))
    for counts_alt in (True, False):
        sub = rows.filter(pl.col("counts_alt") == counts_alt)
        if sub.height == 0:
            continue
        # .bim: alt = A1, so counting ALT = A1 count; counting REF = 2 - A1 count
        dos = a1.astype(np.float64) if counts_alt else (2 - a1).astype(np.float64)
        dos[missing] = 0.0
        w = sparse.csr_matrix((sub.get_column("weight").to_numpy(),
                               ([pos[r] for r in sub.get_column("_row").to_list()],
                                [col[(m, f)] for m, f in sub.select("model_id", "feature_id").iter_rows()])),
                              shape=(len(needed), len(keys)))
        scores += np.asarray(w.T.dot(dos).T)
    for (m, f), b in (intercepts or {}).items():
        if (m, f) in col:
            scores[:, col[(m, f)]] += b
    out = samples.select("FID", "IID") if "FID" in samples.columns else samples.select("IID")
    return out.with_columns([pl.Series(f, scores[:, i]) for i, (_m, f) in enumerate(keys)])


# ------------------------------------------------------------------ plink2 backend
def keyed_genotype(project: Project, fs: GenotypeFileset, cache: Path, threads: int = 1,
                   step_id: str | None = None) -> GenotypeFileset:
    """PGEN copy with unique chromosome:position:REF:ALT ids (made once, then reused)."""
    from efgpp.data.provenance import run_tool

    prefix = cache / "keyed"
    if not prefix.with_suffix(".pgen").exists():
        cache.mkdir(parents=True, exist_ok=True)
        run_tool(project, "plink2", [*fs.plink_input_args(), "--set-all-var-ids", "@:#:$r:$a",
                                     "--new-id-max-allele-len", "1000", "missing", "--rm-dup", "force-first",
                                     "--make-pgen", "--out", str(prefix), "--threads", str(threads)], step_id=step_id)
    return resolve_fileset(prefix, "pgen")


def keyed_sites(fs: GenotypeFileset) -> pl.DataFrame:
    v = read_variants(fs)
    return v.select("variant_id", "chromosome", "position", pl.col("reference").alias("ref"),
                    pl.col("alternate").alias("alt"))


def score_command(fs: GenotypeFileset, score_file: Path, extract: Path, n_features: int, out: Path,
                  threads: int = 1) -> list[str]:
    return [*fs.plink_input_args(), "--extract", str(extract),
            "--score", str(score_file), "1", "2", "header-read", "no-mean-imputation", "cols=+scoresums",
            "--score-col-nums", f"3-{2 + n_features}", "--out", str(out), "--threads", str(threads)]


def score_plink2(project: Project, fs: GenotypeFileset, matched: pl.DataFrame, work: Path, *, threads: int = 1,
                 batch_size: int = 500, step_id: str | None = None,
                 intercepts: dict[tuple[str, str], float] | None = None) -> pl.DataFrame:
    """Scores with `plink2 --score` on a keyed genotype (geno_variant_id = keyed variant id)."""
    from efgpp.data.provenance import run_tool

    keys = feature_keys(matched)
    names = {k: f"F{i}" for i, k in enumerate(keys)}
    work.mkdir(parents=True, exist_ok=True)
    total: pl.DataFrame | None = None
    for b in range(0, len(keys), batch_size):
        batch = keys[b:b + batch_size]
        sub = matched.join(pl.DataFrame(batch, schema=["model_id", "feature_id"], orient="row"),
                           on=["model_id", "feature_id"], how="semi").with_columns(
            pl.struct("model_id", "feature_id").map_elements(lambda s: names[(s["model_id"], s["feature_id"])],
                                                             return_dtype=pl.Utf8).alias("_col"))
        cols = [names[k] for k in batch]
        for counts_alt in (True, False):
            part = sub.filter(pl.col("counts_alt") == counts_alt)
            if part.height == 0:
                continue
            wide = part.pivot(on="_col", index=["geno_variant_id", "counted_allele"], values="weight",
                              aggregate_function="sum")
            wide = wide.with_columns([pl.lit(0.0).alias(c) for c in cols if c not in wide.columns]) \
                .select(pl.col("geno_variant_id").alias("ID"), pl.col("counted_allele").alias("A1"),
                        *[pl.col(c).fill_null(0.0) for c in cols])
            tag = f"b{b // batch_size}_{'alt' if counts_alt else 'ref'}"
            score_file, extract = work / f"{tag}.score.tsv", work / f"{tag}.extract.txt"
            wide.write_csv(score_file, separator="\t")
            extract.write_text("\n".join(wide.get_column("ID").to_list()) + "\n", encoding="utf-8")
            run_tool(project, "plink2", score_command(fs, score_file, extract, len(cols), work / tag, threads),
                     step_id=step_id)
            sscore = pl.read_csv(work / f"{tag}.sscore", separator="\t", infer_schema_length=0)
            id_cols = [c for c in ("#FID", "FID", "#IID", "IID") if c in sscore.columns]
            got = sscore.select(
                *[pl.col(c).alias(c.lstrip("#")) for c in id_cols],
                *[pl.col(f"{c}_SUM").cast(pl.Float64).alias(c) for c in cols if f"{c}_SUM" in sscore.columns])
            got = got.with_columns([pl.lit(0.0).alias(c) for c in cols if c not in got.columns])
            ids = [c.lstrip("#") for c in id_cols]
            if total is None:
                total = got.select(*ids, *[pl.lit(0.0).alias(n) for n in names.values()])
            total = total.join(got, on=ids, how="left", suffix="_new").with_columns(
                [(pl.col(c) + pl.col(f"{c}_new").fill_null(0.0)).alias(c) for c in cols]).drop(
                [f"{c}_new" for c in cols])
    if total is None:
        samples = read_samples(fs)
        total = samples.select("IID")
    for (m, f), bval in (intercepts or {}).items():
        if (m, f) in names:
            total = total.with_columns((pl.col(names[(m, f)]) + bval).alias(names[(m, f)]))
    return total.rename({names[k]: k[1] for k in keys if names[k] in total.columns})
