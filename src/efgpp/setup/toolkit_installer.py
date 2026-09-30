"""Installing and checking toolkits (see efgpp.setup.toolkits for what each contains).

Layout under the install root (<project>/software; $EFGPP_TOOLS_HOME only with --shared):

    envs/<env>/          conda environments (mamba > micromamba > conda; micromamba is
                         downloaded automatically when none is installed)
    opt/<repo>/          GitHub repositories and unpacked downloads
    bin/                 commands: links to environment executables and wrapper scripts
    toolkits/<name>.yaml what was installed, with per-item status
    logs/<name>.install.log  full installer output
"""

from __future__ import annotations

import os
import shutil
import stat
import subprocess
import tarfile
import zipfile
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from efgpp.project import Project
from efgpp.resources.downloader import download
from efgpp.setup.installers import InstallError, install_root
from efgpp.setup.micromamba import find_conda_tool
from efgpp.setup.platform import detect
from efgpp.setup.toolkits import CRAN, TOOLKITS, Download, Repo, Toolkit, resolve_names

MICROMAMBA_URL = "https://micro.mamba.pm/api/micromamba/{platform}/latest"
# Only community channels: never Anaconda's commercial `defaults` channel from a user's .condarc.
CHANNELS = ["--override-channels", "-c", "conda-forge", "-c", "bioconda"]

Progress = Callable[[str], None]


@dataclass
class ItemResult:
    kind: str  # conda | pip | R | repo | download | command
    name: str
    status: str  # ok | installed | failed | missing | manual
    detail: str = ""


@dataclass
class ToolkitResult:
    toolkit: str
    items: list[ItemResult] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not any(i.status in ("failed", "missing") for i in self.items)

    def add(self, kind: str, name: str, status: str, detail: str = "") -> None:
        self.items.append(ItemResult(kind, name, status, detail))


# ---------------------------------------------------------------------- helpers
def _executable(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _log_run(cmd: list[str], log: Path, cwd: Path | None = None, env: dict[str, str] | None = None) -> int:
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a", encoding="utf-8") as fh:
        fh.write(f"\n$ {' '.join(cmd)}\n")
        fh.flush()
        return subprocess.run(cmd, cwd=cwd, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode


def _tail(log: Path, n: int = 12) -> str:
    return "".join(log.read_text(encoding="utf-8", errors="replace").splitlines(True)[-n:])


def env_prefix(root: Path, env: str) -> Path:
    return root / "envs" / env


def env_bin(root: Path, env: str) -> Path:
    return env_prefix(root, env) / "bin"


def _link(target: Path, bin_dir: Path, name: str) -> Path:
    bin_dir.mkdir(parents=True, exist_ok=True)
    link = bin_dir / name
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(target)
    return link


def _wrapper(bin_dir: Path, name: str, body: str) -> Path:
    bin_dir.mkdir(parents=True, exist_ok=True)
    path = bin_dir / name
    if path.exists() or path.is_symlink():
        path.unlink()
    path.write_text("#!/bin/sh\n" + body + "\n", encoding="utf-8")
    _executable(path)
    return path


# ------------------------------------------------------------------- conda envs
def conda_tool(project: Project | None, root: Path, say: Progress) -> tuple[str, Path]:
    """mamba > micromamba > conda; downloads micromamba into <root>/bin when none exists."""
    pref = project.config.execution.fallback_environment_manager if project else "mamba"
    found = find_conda_tool(project, pref if pref != "system" else "mamba") if project else None
    if found is None:
        for kind in ("mamba", "micromamba", "conda"):
            exe = shutil.which(kind)
            if exe:
                return kind, Path(exe)
        local = root / "bin" / "micromamba"
        if local.exists():
            return "micromamba", local
        info = detect()
        platform = {("linux", "x86_64"): "linux-64", ("linux", "aarch64"): "linux-aarch64",
                    ("darwin", "x86_64"): "osx-64", ("darwin", "aarch64"): "osx-arm64"}.get((info.os, info.arch))
        if platform is None:
            raise InstallError(f"no micromamba build for {info.os}/{info.arch}")
        say("downloading micromamba (no conda/mamba found)")
        archive = root / "downloads" / "micromamba.tar.bz2"
        download(MICROMAMBA_URL.format(platform=platform), archive)
        with tarfile.open(archive, "r:bz2") as tar:
            member = tar.getmember("bin/micromamba")
            data = tar.extractfile(member).read()  # type: ignore[union-attr]
        local.parent.mkdir(parents=True, exist_ok=True)
        local.write_bytes(data)
        _executable(local)
        return "micromamba", local
    return found  # type: ignore[return-value]


def ensure_conda_env(tool: tuple[str, Path], prefix: Path, packages: list[str], log: Path, say: Progress) -> None:
    kind, exe = tool
    exists = (prefix / "conda-meta").exists()
    say(f"{'updating' if exists else 'creating'} environment {prefix.name} with {kind} ({len(packages)} packages)")
    verb = "install" if exists else "create"
    env = dict(os.environ)
    if kind == "micromamba":
        env.setdefault("MAMBA_ROOT_PREFIX", str(prefix.parent.parent / "micromamba-root"))
    code = _log_run([str(exe), verb, "-y", "-p", str(prefix), *CHANNELS, *packages], log, env=env)
    if code != 0:
        raise InstallError(f"{kind} {verb} failed for {prefix.name} (log: {log})\n{_tail(log)}")


# ----------------------------------------------------------------------- R / pip
R_SCRIPT = r"""
options(repos = c(CRAN = "{cran}"), timeout = 600, warn = 1)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
install_mode <- {install}
report <- function(pkg, ok, msg = "") cat(sprintf("EFGPP_R\t%s\t%s\t%s\n", pkg, if (ok) "ok" else "failed", gsub("[\r\n\t]", " ", msg)))
have <- function(pkg) requireNamespace(pkg, quietly = TRUE)
try_install <- function(expr) tryCatch({{ expr; "" }}, error = function(e) conditionMessage(e))
for (p in c({cran_pkgs})) {{
  msg <- if (install_mode && !have(p)) try_install(install.packages(p)) else ""
  report(p, have(p), msg)
}}
for (p in c({bioc_pkgs})) {{
  msg <- ""
  if (install_mode && !have(p)) {{
    if (!have("BiocManager")) install.packages("BiocManager")
    msg <- try_install(BiocManager::install(p, update = FALSE, ask = FALSE))
  }}
  report(p, have(p), msg)
}}
for (r in c({gh_pkgs})) {{
  p <- sub("@.*$", "", basename(r))
  msg <- ""
  if (install_mode && !have(p)) {{
    if (!have("remotes")) install.packages("remotes")
    msg <- try_install(remotes::install_github(r, upgrade = "never", dependencies = TRUE))
  }}
  report(p, have(p), msg)
}}
"""


def _r_vector(items: list[str]) -> str:
    return ", ".join(f'"{i}"' for i in items)


def install_r_packages(root: Path, tk: Toolkit, log: Path, say: Progress, result: ToolkitResult,
                       check_only: bool = False) -> None:
    rscript = env_bin(root, tk.env) / "Rscript"  # type: ignore[arg-type]
    if not rscript.exists():
        for p in [*tk.r_cran, *tk.r_bioc, *tk.r_github]:
            result.add("R", p, "missing", "R environment not installed")
        return
    cran, bioc, gh = tk.r_cran, tk.r_bioc, tk.r_github
    if check_only:  # an empty install list turns the script into a pure check
        say("checking R packages")
    else:
        say(f"installing R packages ({len(cran) + len(bioc) + len(gh)}; compiling may take a while)")
    script = R_SCRIPT.format(cran=CRAN, install="FALSE" if check_only else "TRUE", cran_pkgs=_r_vector(cran),
                             bioc_pkgs=_r_vector(bioc), gh_pkgs=_r_vector(gh))
    script_path = log.parent / f"{tk.name}_packages.R"
    script_path.write_text(script, encoding="utf-8")
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([str(env_bin(root, tk.env)), env.get("PATH", "")])  # type: ignore[arg-type]
    proc = subprocess.run([str(rscript), str(script_path)], capture_output=True, text=True, env=env)
    with open(log, "a", encoding="utf-8") as fh:
        fh.write(proc.stdout + proc.stderr)
    seen = set()
    for line in proc.stdout.splitlines():
        if line.startswith("EFGPP_R\t"):
            _, pkg, status, msg = (line.split("\t") + [""])[:4]
            seen.add(pkg)
            result.add("R", pkg, "ok" if status == "ok" else ("missing" if check_only else "failed"), msg[:200])
    for p in [*cran, *bioc, *(sub.split("/")[-1].split("@")[0] for sub in gh)]:
        if p not in seen:
            result.add("R", p, "failed", f"R did not report (log: {log})")
    # Core packages installed through conda are checked by the same mechanism.
    core = [c.split("=")[0].removeprefix("r-") for c in tk.conda if c.startswith("r-") and not c.startswith("r-base")]
    if core and check_only:
        names = {"bigsnpr": "bigsnpr", "data.table": "data.table", "proc": "pROC", "superlearner": "SuperLearner",
                 "r.utils": "R.utils", "rcpparmadillo": "RcppArmadillo", "biocmanager": "BiocManager",
                 "matrix": "Matrix", "rcpp": "Rcpp", "susier": "susieR"}
        check = "; ".join(f'cat("{names.get(c, c)}", requireNamespace("{names.get(c, c)}", quietly=TRUE), "\\n")' for c in core)
        out = subprocess.run([str(rscript), "-e", check], capture_output=True, text=True, env=env).stdout
        for line in out.splitlines():
            parts = line.split()
            if len(parts) == 2:
                result.add("R", parts[0], "ok" if parts[1] == "TRUE" else "missing", "conda")


def install_pip(root: Path, tk: Toolkit, log: Path, say: Progress, result: ToolkitResult) -> None:
    py = env_bin(root, tk.env) / "python"  # type: ignore[arg-type]
    if not tk.pip:
        return
    say(f"installing pip packages into {tk.env}: {', '.join(tk.pip)}")
    for pkg in tk.pip:
        code = _log_run([str(py), "-m", "pip", "install", "--no-input", pkg], log)
        result.add("pip", pkg, "installed" if code == 0 else "failed", "" if code == 0 else f"see {log}")


# ------------------------------------------------------------ repos / downloads
def _archive_url(repo: Repo, commit: str | None = None) -> str:
    return f"https://github.com/{repo.github}/archive/{commit or repo.ref or 'HEAD'}.zip"


def resolve_commit(github: str, ref: str | None) -> str | None:
    """The commit a branch/tag points to now (GitHub API), so installs are pinned, never floating HEAD."""
    import httpx

    try:
        r = httpx.get(f"https://api.github.com/repos/{github}/commits/{ref or 'HEAD'}", timeout=30,
                      headers={"Accept": "application/vnd.github.sha"}, follow_redirects=True)
        return r.text.strip() if r.status_code == 200 and len(r.text.strip()) == 40 else None
    except httpx.HTTPError:
        return None


def write_source_pin(dest: Path, github: str, ref: str | None, commit: str | None, sha256: str) -> Path:
    """software/opt/<name>/EFGPP_SOURCE.json: repository, commit, archive SHA-256 (read by the lock file)."""
    import json
    from datetime import UTC, datetime

    pin = dest / "EFGPP_SOURCE.json"
    pin.write_text(json.dumps({"github": github, "ref": ref, "commit": commit, "archive_sha256": sha256,
                               "installed_at": datetime.now(UTC).isoformat(timespec="seconds")}, indent=1),
                   encoding="utf-8")
    return pin


def install_repo(root: Path, repo: Repo, log: Path, say: Progress, result: ToolkitResult, build_path: list[str]) -> None:
    dest = root / "opt" / repo.name
    commit = resolve_commit(repo.github, repo.ref)
    say(f"downloading {repo.github}@{commit[:10] if commit else repo.ref or 'HEAD'}")
    archive = root / "downloads" / f"{repo.name}-{(commit or repo.ref or 'HEAD')[:12]}.zip"
    archive_sha = download(_archive_url(repo, commit), archive)
    shutil.rmtree(dest, ignore_errors=True)
    staging = root / "opt" / f".{repo.name}.tmp"
    shutil.rmtree(staging, ignore_errors=True)
    with zipfile.ZipFile(archive) as z:
        z.extractall(staging)
        modes = {i.filename: (i.external_attr >> 16) for i in z.infolist()}
    top = next(staging.iterdir())
    for rel, mode in modes.items():  # keep executable bits (e.g. prebuilt SDPR)
        f = staging / rel
        if mode and f.is_file():
            f.chmod(mode & 0o777 or 0o644)
    top.rename(dest)
    shutil.rmtree(staging, ignore_errors=True)
    write_source_pin(dest, repo.github, repo.ref, commit, archive_sha)
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([*build_path, env.get("PATH", "")])
    for cmd in repo.build:
        say(f"{repo.name}: {cmd}")
        if _log_run(["sh", "-c", cmd], log, cwd=dest, env=env) != 0:
            result.add("repo", repo.name, "failed", f"`{cmd}` failed (log: {log})")
            return
    for command, (interp, rel) in repo.commands.items():
        target = dest / rel
        if not target.exists():
            result.add("command", command, "failed", f"{rel} not found in {repo.github}")
            continue
        lib = ""
        if repo.ld_library_path:
            dirs = ":".join(str(dest / d) for d in repo.ld_library_path)
            lib = f'LD_LIBRARY_PATH="{dirs}${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}" '
        if interp == "":
            _executable(target)
            _wrapper(root / "bin", command, f'{lib}exec "{target}" "$@"')
        elif interp == "r":
            _wrapper(root / "bin", command, f'exec "{env_bin(root, "r") / "Rscript"}" "{target}" "$@"')
        else:
            _wrapper(root / "bin", command, f'{lib}exec "{env_bin(root, interp) / "python"}" "{target}" "$@"')
    result.add("repo", repo.name, "manual" if repo.note and not repo.commands else "installed",
               repo.note or f"{repo.github} -> {dest}")


def install_download(root: Path, d: Download, say: Progress, result: ToolkitResult) -> None:
    dest = root / "opt" / d.name
    dest.mkdir(parents=True, exist_ok=True)
    say(f"downloading {d.url}")
    archive = root / "downloads" / d.url.rsplit("/", 1)[-1]
    download(d.url, archive)
    if d.kind == "zip":
        with zipfile.ZipFile(archive) as z:
            for member in z.namelist():
                base = member.rsplit("/", 1)[-1]
                if base in d.members:
                    out = dest / base
                    out.write_bytes(z.read(member))
                    _executable(out)
    elif d.kind == "tgz":
        with tarfile.open(archive) as tar:
            for info in tar.getmembers():
                base = info.name.rsplit("/", 1)[-1]
                if info.isfile() and base in d.members:
                    data = tar.extractfile(info).read()  # type: ignore[union-attr]
                    (dest / base).write_bytes(data)
    elif d.kind == "raw":
        shutil.copy2(archive, dest / archive.name)
    for member, command in d.members.items():
        f = dest / member
        if not f.exists():
            result.add("download", d.name, "failed", f"{member} not found in {d.url}")
            return
        _executable(f)
        if f.suffix == ".R":
            continue  # helper scripts stay in opt/<name>
        _link(f, root / "bin", command)
    result.add("download", d.name, "installed", d.note or d.url)


# ------------------------------------------------------------------ orchestration
def install_toolkit(project: Project | None, name: str, *, shared: bool = False,
                    progress: Progress | None = None) -> ToolkitResult:
    say = progress or (lambda _m: None)
    tk = TOOLKITS[name]
    info = detect()
    if info.is_windows:
        raise InstallError("toolkits target Linux (and macOS); on Windows use WSL2")
    root = install_root(project, shared)
    log = root / "logs" / f"{name}.install.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("", encoding="utf-8")
    result = ToolkitResult(name)
    if tk.env:
        try:
            tool = conda_tool(project, root, say)
            ensure_conda_env(tool, env_prefix(root, tk.env), tk.conda, log, say)
            result.add("conda", tk.env, "installed", f"{len(tk.conda)} packages via {tool[0]}")
        except InstallError as exc:
            result.add("conda", tk.env, "failed", str(exc))
            write_manifest(root, result)
            return result
        for exe in tk.expose:
            target = env_bin(root, tk.env) / exe
            if target.exists():
                _link(target, root / "bin", exe)
        install_pip(root, tk, log, say, result)
        if tk.r_cran or tk.r_bioc or tk.r_github:
            install_r_packages(root, tk, log, say, result)
    if tk.tools:
        from efgpp.setup.installers import install_tools

        for tool_name, status, detail in install_tools(project, tk.tools, shared=shared, progress=say):
            result.add("tool", tool_name, "installed" if status in ("installed", "ok") else status, detail)
    build_path = [str(env_bin(root, tk.env))] if tk.env else []
    for repo in tk.repos:
        try:
            install_repo(root, repo, log, say, result, build_path)
        except Exception as exc:  # noqa: BLE001
            result.add("repo", repo.name, "failed", str(exc))
    for d in tk.downloads:
        try:
            install_download(root, d, say, result)
        except Exception as exc:  # noqa: BLE001
            result.add("download", d.name, "failed", str(exc))
    write_manifest(root, result)
    return result


def install_toolkits(project: Project | None, names: list[str], *, shared: bool = False,
                     progress: Progress | None = None) -> list[ToolkitResult]:
    return [install_toolkit(project, n, shared=shared, progress=progress) for n in resolve_names(names)]


def write_manifest(root: Path, result: ToolkitResult) -> Path:
    path = root / "toolkits" / f"{result.toolkit}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    data: dict[str, Any] = {"toolkit": result.toolkit, "ok": result.ok,
                            "items": [i.__dict__ for i in result.items]}
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return path


def check_toolkit(project: Project | None, name: str, *, shared: bool = False) -> ToolkitResult:
    """Verify every item of a toolkit without installing anything."""
    tk = TOOLKITS[name]
    result = ToolkitResult(name)
    from efgpp.setup.tools import shared_root

    roots = [install_root(project, False)]
    if project is not None:
        roots += list(project.legacy_software_dirs)
    shared_home = shared_root()
    if shared_home is not None:
        roots = [shared_home] if shared else [*roots, shared_home]
    root = next((r for r in roots if tk.env and (env_prefix(r, tk.env) / "conda-meta").exists()), roots[0])
    if tk.tools:
        from efgpp.setup.tools import available

        for tool in tk.tools:
            result.add("tool", tool, "ok" if available(project, tool) else "missing", "")
    if tk.env:
        ok = (env_prefix(root, tk.env) / "conda-meta").exists()
        result.add("conda", tk.env, "ok" if ok else "missing", str(env_prefix(root, tk.env)))
    for cmd in [*tk.checks, *tk.expose, *(c for r in tk.repos for c in r.commands),
                *(c for d in tk.downloads for c in d.members.values() if not c.endswith(".R"))]:
        hit = next((r / "bin" / cmd for r in roots if (r / "bin" / cmd).exists()), None)
        in_env = tk.env and (env_bin(root, tk.env) / cmd).exists()
        result.add("command", cmd, "ok" if hit or in_env else "missing", str(hit or ""))
    for repo in tk.repos:
        present = any((r / "opt" / repo.name).exists() for r in roots)
        result.add("repo", repo.name, ("manual" if repo.note and present else "ok") if present else "missing",
                   repo.note if present else "")
    if tk.env and (tk.r_cran or tk.r_bioc or tk.r_github or any(c.startswith("r-") for c in tk.conda)):
        log = root / "logs" / f"{name}.check.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text("", encoding="utf-8")
        install_r_packages(root, tk, log, lambda _m: None, result, check_only=True)
    if tk.env and tk.pip:
        py = env_bin(root, tk.env) / "python"
        for pkg in tk.pip:
            if not py.exists():
                result.add("pip", pkg, "missing", "environment not installed")
                continue
            code = subprocess.run([str(py), "-m", "pip", "show", "-q", pkg], capture_output=True).returncode
            result.add("pip", pkg, "ok" if code == 0 else "missing")
    # de-duplicate command rows (a command can be both exposed and checked)
    seen, items = set(), []
    for i in result.items:
        key = (i.kind, i.name)
        if key not in seen:
            seen.add(key)
            items.append(i)
    result.items = items
    return result
