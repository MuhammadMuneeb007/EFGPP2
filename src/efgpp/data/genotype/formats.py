"""Genotype file formats: detection, fileset membership, sample and variant listing.

Genotype matrices are never loaded here; only the small sidecar files (.psam/.fam,
.pvar/.bim, VCF headers, BGEN sample blocks) are read. Heavy operations are PLINK2's job.
"""

from __future__ import annotations

import gzip
import io
import struct
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

GENOTYPE_EXTENSIONS = {
    "pgen": (".pgen", ".pvar", ".psam"),
    "bed": (".bed", ".bim", ".fam"),
}


@dataclass
class GenotypeFileset:
    format: str
    prefix: Path  # PLINK prefix, or the main file for BGEN/VCF/BCF
    members: list[Path] = field(default_factory=list)
    missing: list[Path] = field(default_factory=list)

    @property
    def main(self) -> Path:
        return self.members[0] if self.members else self.prefix

    @property
    def complete(self) -> bool:
        return not self.missing

    def plink_input_args(self) -> list[str]:
        """Arguments that make PLINK2 read this fileset."""
        if self.format == "pgen":
            pvar = next((m for m in self.members if ".pvar" in m.name), None)
            if pvar is not None and pvar.name.endswith(".zst"):
                return [
                    "--pgen", str(self.prefix.with_suffix(".pgen")),
                    "--pvar", str(pvar),
                    "--psam", str(self.prefix.with_suffix(".psam")),
                ]
            return ["--pfile", str(self.prefix)]
        if self.format == "bed":
            return ["--bfile", str(self.prefix)]
        if self.format == "bgen":
            args = ["--bgen", str(self.main), "ref-first"]
            sample = _bgen_sample_file(self.main)
            if sample:
                args += ["--sample", str(sample)]
            return args
        if self.format == "vcf":
            return ["--vcf", str(self.main)]
        if self.format == "bcf":
            return ["--bcf", str(self.main)]
        raise ValueError(f"unsupported genotype format {self.format!r}")


def _strip_known(path: Path) -> Path:
    name = path.name
    for suffix in (".pvar.zst", ".pgen", ".pvar", ".psam", ".bed", ".bim", ".fam"):
        if name.endswith(suffix):
            return path.with_name(name[: -len(suffix)])
    return path


def detect_format(path: Path) -> str | None:
    """Detect a genotype format from a path or PLINK prefix."""
    name = path.name.lower()
    if name.endswith((".vcf", ".vcf.gz", ".vcf.bgz")):
        return "vcf"
    if name.endswith(".bcf"):
        return "bcf"
    if name.endswith(".bgen"):
        return "bgen"
    prefix = _strip_known(path)
    if prefix.with_name(prefix.name + ".pgen").exists():
        return "pgen"
    if prefix.with_name(prefix.name + ".bed").exists():
        return "bed"
    return None


def _bgen_sample_file(bgen: Path) -> Path | None:
    for candidate in (bgen.with_suffix(".sample"), Path(str(bgen) + ".sample")):
        if candidate.exists():
            return candidate
    return None


def resolve_fileset(path: Path, fmt: str = "auto") -> GenotypeFileset:
    detected = detect_format(path) if fmt in ("auto", None) else fmt
    if detected is None:
        raise FileNotFoundError(f"cannot find a genotype fileset at {path}")
    fmt = detected
    if fmt in GENOTYPE_EXTENSIONS:
        prefix = _strip_known(path)
        members: list[Path] = []
        missing: list[Path] = []
        for ext in GENOTYPE_EXTENSIONS[fmt]:
            f = prefix.with_name(prefix.name + ext)
            if ext == ".pvar" and not f.exists():
                zst = prefix.with_name(prefix.name + ".pvar.zst")
                f = zst if zst.exists() else f
            (members if f.exists() else missing).append(f)
        return GenotypeFileset(fmt, prefix, members, missing)
    main = path
    members = [main] if main.exists() else []
    missing = [] if main.exists() else [main]
    if fmt == "bgen":
        for extra in (_bgen_sample_file(main), Path(str(main) + ".bgi")):
            if extra is not None and extra.exists():
                members.append(extra)
    else:
        for idx in (".tbi", ".csi"):
            f = Path(str(main) + idx)
            if f.exists():
                members.append(f)
    return GenotypeFileset(fmt, main, members, missing)


# --------------------------------------------------------------------------- samples
def _read_ws_table(path: Path, has_header: bool) -> list[list[str]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        lines = [ln.split() for ln in fh if ln.strip()]
    return lines[1:] if has_header else lines


def read_samples(fs: GenotypeFileset) -> pl.DataFrame:
    """Return a frame with columns FID (nullable), IID and SEX (nullable) in file order."""
    if fs.format == "pgen":
        psam = fs.prefix.with_name(fs.prefix.name + ".psam")
        with open(psam, encoding="utf-8") as fh:
            first = fh.readline()
        if first.startswith("#"):
            header = first.lstrip("#").split()
            rows = _read_ws_table(psam, has_header=True)
        else:  # headerless .psam behaves like .fam
            header = ["FID", "IID", "PAT", "MAT", "SEX", "PHENO1"]
            rows = _read_ws_table(psam, has_header=False)
        df = pl.DataFrame([dict(zip(header, r, strict=False)) for r in rows]) if rows else None
        if df is None:
            return pl.DataFrame(schema={"FID": pl.Utf8, "IID": pl.Utf8, "SEX": pl.Utf8})
        if "FID" not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("FID"))
        if "SEX" not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("SEX"))
        return df.select("FID", "IID", "SEX")
    if fs.format == "bed":
        rows = _read_ws_table(fs.prefix.with_name(fs.prefix.name + ".fam"), has_header=False)
        return pl.DataFrame(
            {"FID": [r[0] for r in rows], "IID": [r[1] for r in rows], "SEX": [r[4] for r in rows]},
            schema={"FID": pl.Utf8, "IID": pl.Utf8, "SEX": pl.Utf8},
        )
    if fs.format == "vcf":
        ids = vcf_samples(fs.main)
    elif fs.format == "bgen":
        ids = bgen_samples(fs.main)
    else:  # BCF is binary; sample listing needs bcftools
        raise NotImplementedError("listing BCF samples requires bcftools (`bcftools query -l`)")
    return pl.DataFrame(
        {"FID": [None] * len(ids), "IID": ids, "SEX": [None] * len(ids)},
        schema={"FID": pl.Utf8, "IID": pl.Utf8, "SEX": pl.Utf8},
    )


def vcf_header(path: Path) -> tuple[list[str], list[str]]:
    """Return (meta lines, column header fields) of a plain or gzip/BGZF VCF."""
    opener = gzip.open if path.name.endswith((".gz", ".bgz")) else open
    meta: list[str] = []
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        for line in fh:
            if line.startswith("##"):
                meta.append(line.rstrip("\n"))
            elif line.startswith("#CHROM"):
                return meta, line.rstrip("\n").lstrip("#").split("\t")
            else:
                break
    raise ValueError(f"{path} has no #CHROM header line")


def vcf_samples(path: Path) -> list[str]:
    _, cols = vcf_header(path)
    return cols[9:]


def bgen_samples(path: Path) -> list[str]:
    """Sample IDs from a .sample sidecar or the BGEN v1.2 embedded sample block."""
    sample = _bgen_sample_file(path)
    if sample is not None:
        rows = _read_ws_table(sample, has_header=False)
        return [r[1] for r in rows[2:]]  # two header lines: names + types
    with open(path, "rb") as fh:
        offset, header_len = struct.unpack("<II", fh.read(8))
        fh.seek(4 + header_len - 4)
        (flags,) = struct.unpack("<I", fh.read(4))
        if not flags >> 31 & 1:
            raise ValueError(f"{path} has no embedded sample IDs and no .sample file")
        fh.seek(4 + header_len)
        _block_len, n = struct.unpack("<II", fh.read(8))
        ids = []
        for _ in range(n):
            (length,) = struct.unpack("<H", fh.read(2))
            ids.append(fh.read(length).decode())
        return ids


def bgen_variant_count(path: Path) -> int:
    with open(path, "rb") as fh:
        fh.read(4)
        _header_len, m = struct.unpack("<II", fh.read(8))
    return int(m)


# -------------------------------------------------------------------------- variants
VARIANT_SCHEMA = {
    "chromosome": pl.Utf8,
    "position": pl.Int64,
    "variant_id": pl.Utf8,
    "reference": pl.Utf8,
    "alternate": pl.Utf8,
}


def _pvar_text(path: Path) -> io.StringIO:
    if path.name.endswith(".zst"):
        raise NotImplementedError(".pvar.zst needs PLINK2 to decompress (plink2 --make-just-pvar)")
    return io.StringIO(path.read_text(encoding="utf-8"))


def read_variants(fs: GenotypeFileset, limit: int | None = None) -> pl.DataFrame:
    """Variant table with harmonized columns. For PLINK .bim, A2 is REF and A1 is ALT."""
    if fs.format == "bed":
        bim = fs.prefix.with_name(fs.prefix.name + ".bim")
        names = ["chromosome", "variant_id", "cm", "position", "alternate", "reference"]
        df = pl.read_csv(bim, separator="\t", has_header=False, n_rows=limit, infer_schema=False)
        if df.width == 6:
            df.columns = names
        else:  # space-delimited .bim
            df = _ws_frame(bim, names, limit)
        return _harmonize_variants(df)
    if fs.format == "pgen":
        pvar = next(m for m in fs.members if ".pvar" in m.name)
        text = _pvar_text(pvar)
        header = None
        body = []
        for line in text:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split()
                continue
            if line.strip():
                body.append(line.split())
                if limit and len(body) >= limit:
                    break
        header = header or ["CHROM", "ID", "CM", "POS", "ALT", "REF"]  # headerless = .bim layout
        df = pl.DataFrame(body, schema=header[: len(body[0])] if body else header, orient="row")
        rename = {"CHROM": "chromosome", "POS": "position", "ID": "variant_id", "REF": "reference", "ALT": "alternate"}
        df = df.rename({k: v for k, v in rename.items() if k in df.columns})
        return _harmonize_variants(df)
    if fs.format == "vcf":
        return _harmonize_variants(pl.DataFrame(list(iter_vcf_sites(fs.main, limit)), schema=list(VARIANT_SCHEMA), orient="row"))
    if fs.format == "bgen":
        raise NotImplementedError("BGEN variant listing requires the .bgi index or PLINK2 export")
    raise NotImplementedError(f"variant listing for {fs.format} requires bcftools")


def _ws_frame(path: Path, names: list[str], limit: int | None) -> pl.DataFrame:
    rows = _read_ws_table(path, has_header=False)
    if limit:
        rows = rows[:limit]
    return pl.DataFrame(rows, schema=names, orient="row")


def _harmonize_variants(df: pl.DataFrame) -> pl.DataFrame:
    for col in VARIANT_SCHEMA:
        if col not in df.columns:
            df = df.with_columns(pl.lit(None).alias(col))
    return df.select(
        pl.col("chromosome").cast(pl.Utf8).str.replace(r"^(?i)chr", ""),
        pl.col("position").cast(pl.Int64),
        pl.col("variant_id").cast(pl.Utf8),
        pl.col("reference").cast(pl.Utf8).str.to_uppercase(),
        pl.col("alternate").cast(pl.Utf8).str.to_uppercase(),
    )


def iter_vcf_sites(path: Path, limit: int | None = None) -> Iterator[tuple[str, int, str, str, str]]:
    opener = gzip.open if path.name.endswith((".gz", ".bgz")) else open
    n = 0
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            yield f[0], int(f[1]), f[2], f[3], f[4]
            n += 1
            if limit and n >= limit:
                return


def count_variants(fs: GenotypeFileset) -> int | None:
    """Count variants cheaply where possible (text sidecars); None when too costly."""
    if fs.format == "bgen":
        return bgen_variant_count(fs.main)
    if fs.format == "bed":
        with open(fs.prefix.with_name(fs.prefix.name + ".bim"), "rb") as fh:
            return sum(1 for ln in fh if ln.strip())
    if fs.format == "pgen":
        pvar = next((m for m in fs.members if ".pvar" in m.name), None)
        if pvar is None or pvar.name.endswith(".zst"):
            return None
        with open(pvar, "rb") as fh:
            return sum(1 for ln in fh if ln.strip() and not ln.startswith(b"#"))
    return None
