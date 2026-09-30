"""GENOTYPE module: what a genotype path is and how its samples are identified.

Used by `efgpp data add genotype --path <prefix or file>` when --format / --sample-id-mode
are not given.

Rules:
  * format: .pgen/.pvar/.psam -> pgen; .bed/.bim/.fam -> bed; .bgen; .vcf/.vcf.gz; .bcf
  * sample IDs: IID when IIDs are unique; FID_IID when IIDs repeat across families
  * genome build: `auto` - detected from the data (reference bases + pyliftover) and lifted to
    the project's target build (GRCh38) by the harmonize step; the files are never modified.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from efgpp.data.genotype.formats import detect_format, read_samples, resolve_fileset


@dataclass
class GenotypeGuess:
    format: str
    prefix: Path
    sample_id_mode: str  # iid | fid_iid
    samples: int
    note: str = ""


def infer_genotype(path: Path) -> GenotypeGuess:
    fmt = detect_format(path)
    if fmt is None:
        raise FileNotFoundError(f"no PLINK/BGEN/VCF genotype found at {path}")
    fs = resolve_fileset(path, fmt)
    try:
        samples = read_samples(fs)
    except NotImplementedError:
        return GenotypeGuess(fmt, fs.prefix, "iid", 0, "sample list needs bcftools (BCF)")
    iids = samples.get_column("IID")
    if iids.is_duplicated().any() and samples.get_column("FID").null_count() == 0:
        return GenotypeGuess(fmt, fs.prefix, "fid_iid", samples.height, "IIDs repeat across families: using FID_IID")
    return GenotypeGuess(fmt, fs.prefix, "iid", samples.height)
