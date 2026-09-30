"""Toolkit catalogue and installer mechanics (no network, no conda)."""

from __future__ import annotations

import io
import zipfile
from pathlib import Path

import pytest

from efgpp.setup import toolkit_installer as ti
from efgpp.setup.toolkits import ALIASES, TOOLKITS, Repo, resolve_names


def test_aliases_and_dependencies_resolve_in_install_order() -> None:
    order = resolve_names(["prs"])
    assert order.index("prs-python") < order.index("prs")
    assert order.index("prs-py27") < order.index("prs")
    assert order.index("r") < order.index("prs")
    assert resolve_names(["all"]) == resolve_names(ALIASES["all"])
    assert resolve_names(["ldpred2"]) == ["r"]
    with pytest.raises(KeyError):
        resolve_names(["nope"])


def test_catalogue_is_consistent() -> None:
    envs = {t.env for t in TOOLKITS.values() if t.env}
    for tk in TOOLKITS.values():
        for repo in tk.repos:
            for _cmd, (interp, _rel) in repo.commands.items():
                assert interp in {"", *envs}, (repo.name, interp)
        for d in tk.downloads:
            assert d.kind in ("zip", "tgz", "gz", "raw") and d.members
    r = TOOLKITS["r"]
    assert "r-bigsnpr" in r.conda  # LDpred-2, SCT, lassosum2
    assert {"tshmak/lassosum", "andrewhaoyu/CTSLEB", "pjnewcombe/R2BGLiMS"} <= set(r.r_github)
    assert "perl" in TOOLKITS["perl"].conda
    assert {"simupop", "msprime", "stdpopsim"} <= {c.split("=")[0] for c in TOOLKITS["simulation"].conda}


def test_r_script_install_and_check_modes() -> None:
    install = ti.R_SCRIPT.format(cran=ti.CRAN, install="TRUE", cran_pkgs=ti._r_vector(["sim1000G"]),
                                 bioc_pkgs="", gh_pkgs=ti._r_vector(["tshmak/lassosum"]), archive_pkgs="")
    check = ti.R_SCRIPT.format(cran=ti.CRAN, install="FALSE", cran_pkgs=ti._r_vector(["sim1000G"]),
                               bioc_pkgs="", gh_pkgs="", archive_pkgs="")
    assert "install_mode <- TRUE" in install and '"tshmak/lassosum"' in install
    assert "install_mode <- FALSE" in check
    assert 'sprintf("EFGPP_R\\t%s\\t%s\\t%s\\n"' in check  # escapes reach R intact


def _fake_repo_zip(files: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, data in files.items():
            z.writestr(f"SDPR-main/{name}", data)
    return buf.getvalue()


def test_repo_install_creates_wrappers(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    payload = _fake_repo_zip({"SDPR": b"#!/bin/sh\necho sdpr\n", "tool.py": b"print(1)\n", "gsl/lib/x.so": b""})
    monkeypatch.setattr(ti, "download", lambda url, dest, **_: dest.parent.mkdir(parents=True, exist_ok=True)
                        or dest.write_bytes(payload))
    repo = Repo("SDPR", "eldronzhou/SDPR", commands={"SDPR": ("", "SDPR"), "tool.py": ("prs-python", "tool.py")},
                ld_library_path=["gsl/lib"])
    result = ti.ToolkitResult("prs")
    ti.install_repo(tmp_path, repo, tmp_path / "log.txt", lambda _m: None, result, [])
    assert (tmp_path / "opt" / "SDPR" / "SDPR").exists()
    sdpr = (tmp_path / "bin" / "SDPR").read_text()
    assert "LD_LIBRARY_PATH=" in sdpr and str(tmp_path / "opt" / "SDPR" / "gsl" / "lib") in sdpr
    py = (tmp_path / "bin" / "tool.py").read_text()
    assert str(tmp_path / "envs" / "prs-python" / "bin" / "python") in py
    assert [i.status for i in result.items] == ["installed"]


def test_check_reports_missing_items(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("EFGPP_TOOLS_HOME", raising=False)
    monkeypatch.setattr(ti, "install_root", lambda project, shared: tmp_path / ("shared" if shared else "proj"))
    res = ti.check_toolkit(None, "simulation")
    assert not res.ok and res.items[0].kind == "conda" and res.items[0].status == "missing"
    (tmp_path / "proj" / "envs" / "simulation" / "conda-meta").mkdir(parents=True)
    assert ti.check_toolkit(None, "simulation").ok


def test_genotype_derived_toolkits() -> None:
    from efgpp.setup.installers import CONDA_TOOLS
    from efgpp.setup.toolkits import TOOLKITS, resolve_names

    assert resolve_names(["predicted-omics"]) == ["metaxcan", "methylation"]
    assert TOOLKITS["metaxcan"].tools == ["predixcan", "plink2"] and TOOLKITS["spliceai"].tools == ["spliceai"]
    assert "bioconductor-hibag" in TOOLKITS["hla"].conda and "HIBAG" in TOOLKITS["hla"].r_bioc
    assert "r-bedmatrix" in TOOLKITS["methylation"].conda
    env, conda, pip, _ = CONDA_TOOLS["spliceai"]
    assert env == "spliceai" and "tensorflow>=2.10,<2.16" in conda  # TensorFlow only in its own environment
    env, conda, pip, _ = CONDA_TOOLS["predixcan"]
    assert "python=3.11" in conda and {"sqlalchemy", "patsy"} <= set(conda) and "bgen-reader" in pip


def test_r_archive_packages_and_cxx_override() -> None:
    from efgpp.setup.toolkit_installer import CRAN, R_SCRIPT, _r_vector
    from efgpp.setup.toolkits import TOOLKITS

    r = TOOLKITS["r"]
    assert r.r_archive == ["hapsim@0.31", "sim1000G@1.40"] and "sim1000G" not in r.r_cran
    assert {"r-gmp", "r-partitions"} <= set(r.conda)  # libgmp for permutations -> partitions
    script = R_SCRIPT.format(cran=CRAN, install="TRUE", cran_pkgs=_r_vector(r.r_cran), bioc_pkgs="",
                             gh_pkgs=_r_vector(r.r_github), archive_pkgs=_r_vector(r.r_archive))
    assert f"{CRAN}/src/contrib/Archive/%s/%s_%s.tar.gz" in script and '"hapsim@0.31", "sim1000G@1.40"' in script


def test_tool_in_software_bin_runs_in_its_env(tmp_path) -> None:  # type: ignore[no-untyped-def]
    import os

    import pytest

    from efgpp.setup.tools import _owning_env

    env = tmp_path / "envs" / "vep"
    (env / "conda-meta").mkdir(parents=True)
    (env / "bin").mkdir()
    real = env / "bin" / "vep"
    real.write_text("#!/bin/sh\n")
    link = tmp_path / "bin" / "vep"
    link.parent.mkdir()
    try:
        os.symlink(real, link)
    except OSError:
        pytest.skip("symlinks not permitted")
    assert _owning_env(link) == env.resolve()  # vep_install then sees the env's tabix/bgzip/perl
    assert _owning_env(real) is None
