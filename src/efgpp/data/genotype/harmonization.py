"""Genome harmonization: chain-file liftover (UCSC chain format).

Liftover always creates a *new* artifact. The original coordinates are never replaced.
Variants that fail to map, map to a different chromosome, or map to the reverse strand
(which would need allele complementing) are excluded and reported, never guessed.
"""

from __future__ import annotations

import bisect
import gzip
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from efgpp.constants import ArtifactStatus, GenomeBuild, Modality, Origin, TemporalType
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.genotype.build import normalize_chrom
from efgpp.data.genotype.formats import count_variants, read_variants, resolve_fileset
from efgpp.data.genotype.loader import write_id_list
from efgpp.data.io import write_parquet
from efgpp.data.provenance import run_tool
from efgpp.data.registry import Registry
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths


@dataclass(frozen=True)
class Block:
    t_start: int  # 0-based, inclusive
    t_end: int  # exclusive
    q_chrom: str
    q_start: int  # 0-based position on the query strand
    q_strand: str
    q_size: int
    score: int


class ChainMap:
    """Point liftover over UCSC chain blocks."""

    def __init__(self, blocks: dict[str, list[Block]]) -> None:
        self.blocks = {c: sorted(b, key=lambda x: x.t_start) for c, b in blocks.items()}
        self.starts = {c: [b.t_start for b in bl] for c, bl in self.blocks.items()}
        self.max_len = {c: max((b.t_end - b.t_start for b in bl), default=0) for c, bl in self.blocks.items()}

    @classmethod
    def from_file(cls, path: Path) -> ChainMap:
        opener = gzip.open if path.name.endswith(".gz") else open
        blocks: dict[str, list[Block]] = {}
        with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[operator]
            header = None
            t_pos = q_pos = 0
            for line in fh:
                parts = line.split()
                if not parts:
                    continue
                if parts[0] == "chain":
                    header = parts
                    t_pos, q_pos = int(parts[5]), int(parts[10])
                    continue
                assert header is not None, "malformed chain file"
                size = int(parts[0])
                t_chrom = normalize_chrom(header[2])
                blocks.setdefault(t_chrom, []).append(Block(
                    t_pos, t_pos + size, normalize_chrom(header[7]), q_pos, header[9], int(header[8]), int(header[1]),
                ))
                if len(parts) == 3:
                    t_pos += size + int(parts[1])
                    q_pos += size + int(parts[2])
        return cls(blocks)

    def map(self, chrom: str, pos: int) -> tuple[str, int, str] | None:
        """Map a 1-based position; returns (chrom, 1-based position, strand) or None."""
        c = normalize_chrom(chrom)
        starts = self.starts.get(c)
        if not starts:
            return None
        p0 = pos - 1
        i = bisect.bisect_right(starts, p0) - 1
        best: Block | None = None
        while i >= 0 and starts[i] > p0 - self.max_len[c] - 1:
            b = self.blocks[c][i]
            if b.t_start <= p0 < b.t_end and (best is None or b.score > best.score):
                best = b
            i -= 1
        if best is None:
            return None
        q = best.q_start + (p0 - best.t_start)
        if best.q_strand == "-":
            q = best.q_size - q - 1
        return best.q_chrom, q + 1, best.q_strand


def liftover_variants(variants: pl.DataFrame, chain: ChainMap) -> pl.DataFrame:
    rows: list[tuple[str, str | None, int | None, str]] = []
    for vid, chrom, pos in variants.select("variant_id", "chromosome", "position").iter_rows():
        m = chain.map(chrom, pos)
        if m is None:
            rows.append((vid, None, None, "unmapped"))
        elif m[0] != normalize_chrom(chrom):
            rows.append((vid, m[0], m[1], "chromosome_changed"))
        elif m[2] == "-":
            rows.append((vid, m[0], m[1], "reverse_strand"))
        else:
            rows.append((vid, m[0], m[1], "mapped"))
    return pl.DataFrame(rows, schema={"variant_id": pl.Utf8, "new_chromosome": pl.Utf8,
                                      "new_position": pl.Int64, "status": pl.Utf8}, orient="row")


def run_liftover(project: Project, artifact_id: str, target_build: str, *, step_id: str | None = None,
                 threads: int = 1) -> str:
    with Registry.open(project) as reg:
        src = ArtifactStore(reg).get(artifact_id)
    source_build = GenomeBuild.normalize(src.genome_build)
    target = GenomeBuild.normalize(target_build)
    if source_build not in (GenomeBuild.GRCH37, GenomeBuild.GRCH38):
        raise RuntimeError(f"{artifact_id}: source genome build is {src.genome_build!r}; set it before lifting over")
    if source_build == target:
        raise RuntimeError(f"{artifact_id} is already {target.value}")
    key = f"{source_build.value}->{target.value}"
    chain_path = project.resources.liftover.chain_files.get(key)
    if not chain_path:
        raise RuntimeError(f"no chain file configured for {key} (resources.yaml: liftover.chain_files)")
    chain = ChainMap.from_file(project.resolve(chain_path))
    fs = resolve_fileset(Path(src.path), src.format)
    variants = read_variants(fs)
    mapping = liftover_variants(variants, chain)

    out_dir = project.artifact_dir(Origin.DERIVED, Modality.GENOTYPE_QC, src.source_id or "lifted", "liftover")
    report = write_parquet(mapping, out_dir / f"{src.artifact_name}_{target.value}_liftover_report.parquet")
    mapped = mapping.filter(pl.col("status") == "mapped")
    work = project.work_root / "liftover" / src.artifact_name
    work.mkdir(parents=True, exist_ok=True)
    exclude = write_id_list(work / "exclude.txt", mapping.filter(pl.col("status") != "mapped").get_column("variant_id").to_list())
    update = work / "update_map.txt"
    update.write_text("".join(f"{v}\t{p}\n" for v, p in mapped.select("variant_id", "new_position").iter_rows()), encoding="utf-8")
    prefix = out_dir / f"{src.artifact_name}_{target.value}"
    rec = run_tool(project, "plink2", [*fs.plink_input_args(), "--threads", str(threads), "--exclude", str(exclude),
                                       "--update-map", str(update), "--sort-vars", "--make-pgen", "--out", str(prefix)],
                   step_id=step_id, inputs=[artifact_id])
    lifted = resolve_fileset(prefix, "pgen")
    counts = {s: int(n) for s, n in mapping.group_by("status").len().iter_rows()}
    with Registry.open(project) as reg:
        art = ArtifactStore(reg).register_replacing(Artifact(
            artifact_name=f"{src.artifact_name}_{target.value}", artifact_type="liftover_genotype",
            modality=src.modality, origin=Origin.DERIVED, status=ArtifactStatus.READY, path=str(prefix),
            format="pgen", size=sum(m.stat().st_size for m in lifted.members),
            checksum=str(checksum_paths(lifted.members)), participant_count=src.participant_count,
            feature_count=count_variants(lifted), genome_build=target.value, temporal_type=TemporalType.STATIC,
            source_id=src.source_id, parent_artifact_ids=[artifact_id], tool="plink2", tool_version=rec.tool_version,
            resource_versions={"chain": Path(chain_path).name}, command=" ".join(rec.command),
            metadata={"from_build": source_build.value, "liftover_counts": counts, "report": str(report)},
        ))
    return art.artifact_id  # type: ignore[return-value]
