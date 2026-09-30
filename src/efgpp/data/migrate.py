"""Migration of existing (legacy) project folders into an EFGPP project.

    efgpp data migrate --source /old/project --profile legacy-efgpp --dry-run   -> migration_plan.yaml
    efgpp data migrate --source /old/project --profile legacy-efgpp --apply

The source directory is only ever read. Participant-level files become sources in
data.yaml (stored by reference by default); GWAS summary statistics and annotation tables
are registered as reference resources; legacy PCA/PRS outputs are registered for
provenance only (target-specific PRS belongs to the Representation layer).
"""

from __future__ import annotations

import gzip
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from efgpp.config.data import CovariateSource, GenotypeSource, OmicsSource, PhenotypeSource
from efgpp.constants import ArtifactStatus, Modality, Origin, PhenotypeType, StorageMode
from efgpp.data.artifacts import Artifact, ArtifactStore
from efgpp.data.registry import Registry, utcnow
from efgpp.project import Project
from efgpp.resources.checksums import checksum_paths

PROFILES = ("legacy-efgpp", "generic")
GWAS_COLUMNS = {"CHR", "BP", "SNP", "P", "A1", "A2", "BETA", "OR", "POS", "PVAL", "P_VALUE", "RSID"}
PRS_SUFFIXES = (".profile", ".sscore", ".all_score", ".prs", ".all.score")
PHENO_SUFFIXES = (".height", ".pheno", ".phen")
SKIP_DIRS = {".git", "__pycache__", ".ipynb_checkpoints", "node_modules"}


@dataclass
class PlanItem:
    kind: str  # genotype | phenotype | covariates | gwas | annotations | pca | prs | expression | unknown
    path: str
    format: str | None = None
    details: dict[str, Any] = field(default_factory=dict)
    action: str = ""


def _header(path: Path, n: int = 1) -> list[str]:
    opener = gzip.open if path.name.endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[operator]
            return [fh.readline() for _ in range(n)]
    except (OSError, EOFError, UnicodeDecodeError):
        return []


def _columns(path: Path) -> list[str]:
    head = _header(path)
    return re.split(r"[\s,]+", head[0].strip().lstrip("#")) if head and head[0].strip() else []


def _values(path: Path, column: str, limit: int = 2000) -> set[str]:
    lines = _header(path, limit + 1)
    cols = re.split(r"[\s,]+", lines[0].strip()) if lines else []
    if column not in cols:
        return set()
    i = cols.index(column)
    out = set()
    for ln in lines[1:]:
        parts = re.split(r"[\s,]+", ln.strip())
        if len(parts) > i:
            out.add(parts[i])
    return out


def _guess_type(values: set[str]) -> PhenotypeType:
    vals = {v for v in values if v not in ("", "NA", "-9", "nan")}
    if vals and vals <= {"0", "1", "2"}:
        return PhenotypeType.BINARY
    return PhenotypeType.CONTINUOUS


def scan(source: Path, profile: str = "legacy-efgpp") -> list[PlanItem]:
    if profile not in PROFILES:
        raise ValueError(f"unknown profile {profile!r}; choose from {PROFILES}")
    source = source.resolve()
    items: list[PlanItem] = []
    seen_prefixes: set[Path] = set()
    files = sorted(p for p in source.rglob("*") if p.is_file() and not (set(p.relative_to(source).parts) & SKIP_DIRS))
    for f in files:
        name = f.name.lower()
        rel = f.relative_to(source).as_posix()
        if name.endswith((".bed", ".pgen")):
            prefix = f.with_suffix("")
            companions = (".bim", ".fam") if name.endswith(".bed") else (".pvar", ".psam")
            if prefix not in seen_prefixes and all(prefix.with_suffix(c).exists() for c in companions):
                seen_prefixes.add(prefix)
                items.append(PlanItem("genotype", str(prefix), "bed" if name.endswith(".bed") else "pgen",
                                      {"relative": rel}, "add genotype source (reference)"))
            continue
        if name.endswith((".bim", ".fam", ".pvar", ".psam", ".log", ".nosex", ".eigenval", ".tbi", ".csi")):
            continue
        if name.endswith((".vcf", ".vcf.gz", ".bcf", ".bgen")):
            fmt = "bgen" if name.endswith(".bgen") else ("bcf" if name.endswith(".bcf") else "vcf")
            items.append(PlanItem("genotype", str(f), fmt, {"relative": rel}, "add genotype source (reference)"))
        elif name.endswith(".cov"):
            cols = _columns(f)
            variables = [c for c in cols if c not in ("FID", "IID")]
            items.append(PlanItem("covariates", str(f), "txt", {"id_column": "IID" if "IID" in cols else cols[0],
                                                                 "variables": variables}, "add covariate source"))
        elif name.endswith(PHENO_SUFFIXES):
            cols = _columns(f)
            value_cols = [c for c in cols if c not in ("FID", "IID")]
            value = value_cols[0] if value_cols else None
            items.append(PlanItem("phenotype", str(f), "txt", {
                "name": f.parent.name if f.parent != source else f.stem,
                "id_column": "IID" if "IID" in cols else (cols[0] if cols else None),
                "value_column": value, "type": _guess_type(_values(f, value)).value if value else None,
            }, "add phenotype"))
        elif name.endswith(".eigenvec"):
            items.append(PlanItem("pca", str(f), "eigenvec", {"relative": rel}, "register legacy PCA (provenance only)"))
        elif name.endswith(PRS_SUFFIXES):
            items.append(PlanItem("prs", str(f), "text", {"relative": rel},
                                  "register legacy PRS (provenance only; PRS belongs to Representation)"))
        elif name.endswith((".h5ad",)) or re.search(r"(expr|rna|transcript)", name):
            items.append(PlanItem("expression", str(f), None, {"relative": rel}, "add expression source (review metadata)"))
        elif "annot" in name:
            items.append(PlanItem("annotations", str(f), "tsv", {"relative": rel}, "register reference annotation table"))
        elif name.endswith((".gz", ".txt", ".tsv", ".csv", ".assoc", ".glm.logistic", ".glm.linear")) and \
                len(GWAS_COLUMNS & {c.upper() for c in _columns(f)}) >= 3:
            items.append(PlanItem("gwas", str(f), "sumstats", {"columns": _columns(f)[:12]},
                                  "register GWAS summary statistics as reference resource"))
        else:
            items.append(PlanItem("unknown", str(f), None, {"relative": rel}, "ignored (review)"))
    _dedupe_names(items)
    return items


def _dedupe_names(items: list[PlanItem]) -> None:
    counts: dict[str, int] = {}
    for it in items:
        if it.kind != "phenotype":
            continue
        base = re.sub(r"\W+", "_", str(it.details["name"])).strip("_") or "phenotype"
        counts[base] = counts.get(base, 0) + 1
        it.details["name"] = base if counts[base] == 1 else f"{base}_{counts[base]}"


def write_plan(project: Project, source: Path, profile: str, mode: StorageMode, items: list[PlanItem]) -> Path:
    summary: dict[str, int] = {}
    for it in items:
        summary[it.kind] = summary.get(it.kind, 0) + 1
    plan = {"source": str(source.resolve()), "profile": profile, "mode": mode.value,
            "created_at": utcnow().isoformat(timespec="seconds"), "summary": summary,
            "items": [asdict(i) for i in items]}
    path = project.path("migration_plan.yaml")
    path.write_text(yaml.safe_dump(plan, sort_keys=False), encoding="utf-8")
    return path


def apply(project: Project, source: Path, profile: str, mode: StorageMode, plan_path: Path | None = None) -> dict[str, list[str]]:
    plan_path = plan_path or project.path("migration_plan.yaml")
    if plan_path.exists():
        plan = yaml.safe_load(plan_path.read_text(encoding="utf-8"))
        if Path(plan["source"]) != source.resolve():
            raise RuntimeError(f"{plan_path} was made for {plan['source']}; re-run with --dry-run")
        items = [PlanItem(**i) for i in plan["items"]]
    else:
        items = scan(source, profile)
    obs = project.data.observed
    known_paths = {s.path for _, s in project.data.iter_sources()}
    done: dict[str, list[str]] = {}

    def note(kind: str, what: str) -> None:
        done.setdefault(kind, []).append(what)

    for it in items:
        if it.path in known_paths:
            continue
        if it.kind == "genotype":
            sid = project.data.next_id("GENO")
            obs.genotype.append(GenotypeSource(id=sid, path=it.path, format=it.format or "auto", mode=mode))
            note("genotype", sid)
        elif it.kind == "phenotype" and it.details.get("value_column") and it.details.get("type"):
            sid = project.data.next_id("PH")
            obs.phenotypes.append(PhenotypeSource(
                id=sid, name=it.details["name"], path=it.path, format="txt", mode=mode,
                participant_id_column=it.details["id_column"], value_column=it.details["value_column"],
                type=PhenotypeType(it.details["type"])))
            note("phenotype", sid)
        elif it.kind == "covariates" and it.details.get("variables"):
            sid = project.data.next_id("COV")
            obs.covariates.append(CovariateSource(id=sid, path=it.path, format="txt", mode=mode,
                                                  participant_id_column=it.details["id_column"],
                                                  variables=it.details["variables"]))
            note("covariates", sid)
        elif it.kind == "expression":
            sid = project.data.next_id("RNA")
            obs.expression.append(OmicsSource(id=sid, path=it.path, mode=mode))
            note("expression", sid)
        project.data = type(project.data).model_validate(project.data.model_dump())
        obs = project.data.observed
    project.save_data_config()

    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        for it in items:
            p = Path(it.path)
            if it.kind in ("gwas", "annotations"):
                checksum = str(checksum_paths([p]))
                rid = f"RES_legacy_{it.kind}_{p.name}"
                reg.upsert("resources", {
                    "resource_id": rid, "name": f"legacy_{it.kind}", "version": p.name, "release_date": None,
                    "genome_build": None, "source": f"migrated from {source}", "license": "unknown (legacy)",
                    "checksum": checksum, "download_date": utcnow().isoformat(timespec="seconds"),
                    "local_path": str(p), "metadata": it.details,
                })
                note(it.kind, rid)
            elif it.kind in ("pca", "prs"):
                art = store.register_replacing(Artifact(
                    artifact_name=f"legacy_{it.kind}_{p.stem}", artifact_type=f"legacy_{it.kind}",
                    modality=Modality.QC_PCA if it.kind == "pca" else "legacy_prs", origin=Origin.DERIVED,
                    status=ArtifactStatus.REGISTERED, path=str(p), format=it.format or "text", size=p.stat().st_size,
                    checksum=str(checksum_paths([p])), source_id=f"LEGACY_{it.kind.upper()}", storage_mode="reference",
                    tool="legacy", metadata={"migrated_from": str(source), "note": it.action},
                ))
                note(it.kind, art.artifact_id or "")
    log = project.path("migration_applied.yaml")
    log.write_text(yaml.safe_dump({"source": str(source), "profile": profile, "applied_at": utcnow().isoformat(),
                                   "registered": done}, sort_keys=False), encoding="utf-8")
    return done
