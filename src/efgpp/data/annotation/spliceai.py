"""SpliceAI (Illumina): variant-level splice-disruption scores.

SpliceAI predicts, for each variant, the probability of acceptor/donor gain/loss (delta scores
DS_AG, DS_AL, DS_DG, DS_DL) and the positions of those changes (DP_*). It is NOT the same as
GTEx sQTL-predicted participant splicing (predicted.splicing). Raw scores are always kept; the
configurable thresholds (0.2 / 0.5 / 0.8) are only used for reporting counts, never to label a
variant pathogenic or benign.

Runs in its own environment (software/envs/spliceai, TensorFlow inside) on the sites VCF of
the cohort, with the project's GRCh38 FASTA.

Licence: SpliceAI source code GPL-3.0; the trained models (shipped with the package) are
CC BY-NC 4.0 - academic and non-commercial use.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl

from efgpp.data.annotation import (
    cohort_variants,
    register_annotation,
    with_variant_key,
    write_sites_vcf,
)
from efgpp.data.provenance import run_tool
from efgpp.project import Project

LICENSE = ("SpliceAI (Illumina) source code GPL-3.0; trained models CC BY-NC 4.0 "
           "(academic and non-commercial use only)")
SCORE_FIELDS = ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]
POSITION_FIELDS = ["DP_AG", "DP_AL", "DP_DG", "DP_DL"]
EFFECTS = {"DS_AG": "acceptor_gain", "DS_AL": "acceptor_loss", "DS_DG": "donor_gain", "DS_DL": "donor_loss"}


def parse_spliceai_info(value: str) -> list[dict[str, object]]:
    """SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL[,...] -> one dict per gene."""
    out: list[dict[str, object]] = []
    for entry in value.split(","):
        f = entry.split("|")
        if len(f) < 10:
            continue
        scores = {k: _num(v) for k, v in zip(SCORE_FIELDS, f[2:6], strict=True)}
        positions = {k: _int(v) for k, v in zip(POSITION_FIELDS, f[6:10], strict=True)}
        valid = {k: v for k, v in scores.items() if v is not None}
        best = max(valid, key=lambda k: valid[k]) if valid else None
        out.append({"spliceai_allele": f[0], "spliceai_gene": f[1], **scores, **positions,
                    "spliceai_max_score": valid[best] if best else None,
                    "spliceai_effect": EFFECTS[best] if best else None})
    return out


def _num(v: str) -> float | None:
    try:
        return float(v)
    except ValueError:
        return None


def _int(v: str) -> int | None:
    try:
        return int(v)
    except ValueError:
        return None


def parse_spliceai_vcf(path: Path) -> pl.DataFrame:
    """One row per variant x gene from a SpliceAI-annotated VCF (ID column = variant key)."""
    rows = []
    opener = gzip.open if path.name.endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 8)
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if "SpliceAI" not in info:
                continue
            base = {"variant_key": f[2], "chromosome": f[0].removeprefix("chr"), "position": int(f[1]),
                    "reference": f[3], "alternate": f[4]}
            rows += [{**base, **r} for r in parse_spliceai_info(info["SpliceAI"])]
    schema = {"variant_key": pl.Utf8, "chromosome": pl.Utf8, "position": pl.Int64, "reference": pl.Utf8,
              "alternate": pl.Utf8, "spliceai_allele": pl.Utf8, "spliceai_gene": pl.Utf8,
              **dict.fromkeys(SCORE_FIELDS, pl.Float64), **dict.fromkeys(POSITION_FIELDS, pl.Int64),
              "spliceai_max_score": pl.Float64, "spliceai_effect": pl.Utf8}
    return pl.DataFrame(rows, schema=schema)


def _with_contigs(sites: Path, fasta: Path | None) -> None:
    """Add ##contig lines (lengths from the FASTA index) so pysam/SpliceAI accept the VCF."""
    if fasta is None or not Path(str(fasta) + ".fai").exists():
        return
    lengths = {}
    for line in Path(str(fasta) + ".fai").read_text(encoding="utf-8").splitlines():
        name, length = line.split("\t")[:2]
        lengths[name.removeprefix("chr")] = length
    lines = sites.read_text(encoding="utf-8").splitlines(keepends=True)
    used = dict.fromkeys(ln.split("\t", 1)[0] for ln in lines if not ln.startswith("#"))
    contigs = [f"##contig=<ID={c},length={lengths[c.removeprefix('chr')]}>\n" for c in used
               if c.removeprefix("chr") in lengths]
    sites.write_text(lines[0] + "".join(contigs) + "".join(lines[1:]), encoding="utf-8")


def run_spliceai(project: Project, source_id: str, *, step_id: str | None = None, threads: int = 1) -> str:
    from efgpp.data.references.genome import installed_fasta

    variants, parent, build = cohort_variants(project, source_id)
    fasta = installed_fasta(project)
    if fasta is None:
        raise RuntimeError("SpliceAI needs the reference FASTA: efgpp resources install genome")
    if build not in ("GRCh37", "GRCh38"):
        raise RuntimeError(f"{source_id}: genome build {build!r}; SpliceAI supports GRCh37/GRCh38 annotations")
    work = project.work_root / "spliceai" / source_id
    sites = write_sites_vcf(with_variant_key(variants), work / "sites.vcf", build)
    _with_contigs(sites, fasta)
    out = work / "spliceai.vcf"
    distance = str(project.resources.spliceai.distance)
    rec = run_tool(project, "spliceai", ["-I", str(sites), "-O", str(out), "-R", str(fasta),
                                         "-A", "grch38" if build == "GRCh38" else "grch37", "-D", distance],
                   step_id=step_id, inputs=[parent.artifact_id])  # type: ignore[list-item]
    table = parse_spliceai_vcf(out)
    thresholds = project.data.participant_variants.spliceai_thresholds
    per_variant = table.group_by("variant_key").agg(pl.col("spliceai_max_score").max())
    return register_annotation(
        project, source_id=source_id, name="spliceai", table=table, parent=parent, tool="spliceai",
        tool_version=rec.tool_version, resource_versions={"spliceai": rec.tool_version or "unknown",
                                                          "annotation": f"spliceai built-in {build}"},
        genome_build=build,
        metadata={"license": LICENSE, "distance": int(distance), "command": " ".join(rec.command),
                  "variants_scored": per_variant.height,
                  "at_or_above": {str(t): int((per_variant.get_column("spliceai_max_score") >= t).sum())
                                  for t in thresholds},
                  "note": "variant-level splice-disruption prediction; not participant splicing"},
    )
