"""Stand-alone GWASLab job, executed with the Python of the `gwaslab` conda environment
(software/envs/gwaslab) - never inside the efgpp environment. It imports only the standard
library, GWASLab and pandas/pyarrow, so it runs without EFGPP installed.

    gwaslab-python gwas_gwaslab_job.py job.json

job.json: {input, fmt, build ("19"/"38"/"99"), columns {gwaslab keyword: column}, n,
           target ("38"), chains {"19->38": path, "38->19": path}, output, report}
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path


def main(job_path: str) -> int:
    job = json.loads(Path(job_path).read_text(encoding="utf-8"))
    import importlib.metadata as md

    import gwaslab as gl

    report: dict = {"gwaslab_version": md.version("gwaslab"), "input": job["input"], "started": time.time()}
    kwargs = dict(job.get("columns") or {})
    if job.get("n") is not None:
        kwargs["n"] = int(job["n"])
    fmt = job.get("fmt") or "auto"
    build = job.get("build") or "99"
    try:
        ss = gl.Sumstats(job["input"], fmt=None if fmt == "none" else fmt, build=build, **kwargs)
    except Exception as exc:  # noqa: BLE001 - retry with explicit columns only
        if not kwargs:
            raise
        report["fmt_auto_failed"] = f"{type(exc).__name__}: {exc}"
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
    if str(ss.build) != target:
        chain = job["chains"][f"{ss.build}->{target}"]
        ss.liftover(from_build=str(ss.build), to_build=target, chain_path=chain)
        report["lifted"] = f"{report['build_detected']}->{target}"
        report["rows_after_liftover"] = int(len(ss.data))
    else:
        report["lifted"] = None
    report["build_final"] = str(ss.build)

    data = ss.data
    report["rows_final"] = int(len(data))
    report["columns"] = [str(c) for c in data.columns]
    out = Path(job["output"])
    out.parent.mkdir(parents=True, exist_ok=True)
    data.to_parquet(out, index=False)
    report["finished"] = time.time()
    Path(job["report"]).write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
