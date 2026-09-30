"""Stand-alone GWASLab job, executed with the Python of the `gwaslab` conda environment
(software/envs/gwaslab) - never inside the efgpp environment. It imports only the standard
library, GWASLab and pandas/pyarrow, so it runs without EFGPP installed.

    gwaslab-python gwas_gwaslab_job.py job.json

job.json: {input, fmt, build ("19"/"38"/"99"), columns {gwaslab keyword: column}, other [columns],
           n, ncase, ncontrol, target ("38"), chains {"19->38": path, "38->19": path},
           output (original column names), output_gwaslab (GWASLab names), report}

Two outputs:
  output          the original column names of the input file (values checked and lifted by
                  GWASLab; rows removed by basic_check are gone). Columns GWASLab adds are kept
                  with GWASLab's names (e.g. N_CASE / N_CONTROL when given as numbers).
  output_gwaslab  GWASLab's standard names (SNPID, CHR, POS, EA, NEA, BETA, SE, P, ..., STATUS)
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

# GWASLab keyword -> the column name GWASLab gives it
GWASLAB_NAMES = {
    "snpid": "SNPID", "rsid": "rsID", "chrom": "CHR", "pos": "POS", "ea": "EA", "nea": "NEA", "eaf": "EAF",
    "beta": "BETA", "OR": "OR", "se": "SE", "z": "Z", "p": "P", "mlog10p": "MLOG10P", "n": "N",
    "ncase": "N_CASE", "ncontrol": "N_CONTROL", "info": "INFO",
}


def _pyliftover(ss, chain: str) -> None:  # type: ignore[no-untyped-def]
    """Lift CHR/POS with pyliftover (same chain file). Variants that do not map to the same
    chromosome on the forward strand are removed, as GWASLab does."""
    from pyliftover import LiftOver

    lo = LiftOver(chain)
    data = ss.data
    ucsc = {"23": "X", "24": "Y", "25": "X", "26": "M", "MT": "M"}
    new_pos, keep = [], []
    for c, p in zip(data["CHR"].astype(str), data["POS"], strict=True):
        cc = ucsc.get(c, c)
        hits = lo.convert_coordinate(f"chr{cc}", int(p) - 1)
        ok = bool(hits) and hits[0][0] == f"chr{cc}" and hits[0][2] == "+"
        keep.append(ok)
        new_pos.append(int(hits[0][1]) + 1 if ok else 0)
    data = data.assign(POS=new_pos)[keep].copy()
    data["POS"] = data["POS"].astype("Int64")
    ss.data = data


def main(job_path: str) -> int:
    job = json.loads(Path(job_path).read_text(encoding="utf-8"))
    import importlib.metadata as md

    import gwaslab as gl

    report: dict = {"gwaslab_version": md.version("gwaslab"), "input": job["input"], "started": time.time()}
    columns = dict(job.get("columns") or {})
    kwargs: dict = dict(columns)
    for key in ("n", "ncase", "ncontrol"):
        if job.get(key) is not None and key not in columns:
            kwargs[key] = int(job[key])  # constant sample sizes become N / N_CASE / N_CONTROL
    if job.get("other"):
        kwargs["other"] = list(job["other"])  # keep every unmapped column
    fmt = job.get("fmt") or "auto"
    build = job.get("build") or "99"
    try:
        ss = gl.Sumstats(job["input"], fmt=None if fmt in ("none", "auto") and columns else fmt, build=build, **kwargs)
    except Exception as exc:  # noqa: BLE001 - retry with explicit columns only
        if not columns:
            raise
        report["fmt_failed"] = f"{type(exc).__name__}: {exc}"
        ss = gl.Sumstats(job["input"], build=build, **kwargs)
    report["rows_loaded"] = int(len(ss.data))

    # Standardize and QC: IDs, chromosomes, positions, alleles, statistics; remove bad/duplicated rows.
    ss.basic_check(remove=True, remove_dup=True)
    report["rows_after_basic_check"] = int(len(ss.data))

    if build in ("99", "auto", None):
        ss.infer_build()
    report["build_detected"] = str(ss.build)
    if str(ss.build) not in ("19", "38"):
        raise SystemExit(f"GWASLab could not determine the genome build of {job['input']}; "
                         "set `build` for this GWAS in data.yaml")

    target = job["target"]
    report["build_final"] = str(ss.build)
    if str(ss.build) != target:
        chain = job["chains"][f"{ss.build}->{target}"]
        try:
            ss.liftover(from_build=str(ss.build), to_build=target, chain_path=chain)
            report["liftover_engine"] = "gwaslab"
            report["build_final"] = str(ss.build)
        except ImportError as exc:  # GWASLab without its optional `sumstats-liftover` package
            report["liftover_engine"] = f"pyliftover (GWASLab liftover unavailable: {exc})"
            _pyliftover(ss, chain)
            report["build_final"] = target
        report["lifted"] = f"{report['build_detected']}->{target}"
        report["rows_after_liftover"] = int(len(ss.data))
    else:
        report["lifted"] = None

    data = ss.data
    report["rows_final"] = int(len(data))
    report["columns_gwaslab"] = [str(c) for c in data.columns]
    Path(job["output_gwaslab"]).parent.mkdir(parents=True, exist_ok=True)
    data.to_parquet(job["output_gwaslab"], index=False)

    # Same rows and values, original column names.
    back = {GWASLAB_NAMES[k]: v for k, v in columns.items() if k in GWASLAB_NAMES}
    original = data.drop(columns=[c for c in ("STATUS",) if c in data.columns]).rename(columns=back)
    report["columns_original"] = [str(c) for c in original.columns]
    original.to_parquet(job["output"], index=False)

    report["finished"] = time.time()
    Path(job["report"]).write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
