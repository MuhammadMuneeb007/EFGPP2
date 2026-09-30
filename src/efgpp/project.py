"""The EFGPP project: a cohort rooted at a directory containing project.yaml."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from efgpp.config import DataConfig, ProjectConfig, ResourcesConfig, load_model, write_yaml
from efgpp.config.models import file_checksum
from efgpp.constants import Modality, Origin

PROJECT_FILE = "project.yaml"
DATA_FILE = "data.yaml"
RESOURCES_FILE = "resources.yaml"
LOCK_FILE = "efgpp.lock.yaml"

DATA_YAML_HEADER = """EFGPP data.yaml - every participant-level source of this cohort.
Edit by hand or with `efgpp data add ...` / `efgpp phenotype add ...`.
Phenotypes are configuration objects: add as many as you need."""

# Directory layout created by `efgpp init` (section 10 of the specification).
PROJECT_DIRS = (
    # Installed software, in plain sight: one conda environment per tool, commands, sources.
    "software/bin",
    "software/envs",
    "software/opt",
    "software/downloads",
    "logs",
    # Internal bookkeeping only (locks, caches, temporary downloads).
    ".efgpp/cache",
    ".efgpp/downloads",
    ".efgpp/locks",
    "registry",
    *(f"data/observed/{m}" for m in (
        "genotype", "phenotype", "covariates", "expression",
        "methylation", "proteomics", "metabolomics", "clinical",
    )),
    *(f"data/derived/{m}" for m in (
        "genotype_qc", "ancestry", "qc_pca", "kinship", "roh", "hla",
        "genotype_imputation", "variant_annotations",
    )),
    *(f"data/predicted/{m}" for m in ("expression", "proteomics", "metabolomics")),
    *(f"data/simulated/{m}" for m in ("genotype", "phenotype", "expression", "covariates", "truth")),
    *(f"resources/{r}" for r in (
        "genomes", "vep", "opencravat", "clinvar", "gnomad", "alphamissense", "alphagenome",
        "gtex", "predictdb", "omicspred", "gwas_catalog", "pgs_catalog", "opentargets",
        "ld_reference", "ontologies", "pathways",
    )),
    *(f"qc/{m}" for m in (
        "genotype", "phenotype", "covariates", "expression", "methylation",
        "proteomics", "metabolomics", "clinical",
    )),
    "phenotypes",
    "snapshots",
    "work",
    "reports/data",
    "workflow/rules",
    "workflow/envs",
    "workflow/profiles/local",
    "workflow/profiles/slurm",
)

# Where each (origin, modality) pair lives under data/.
_MODALITY_DIR = {
    Modality.GENOTYPE_QC: "genotype_qc",
    Modality.QC_PCA: "qc_pca",
    Modality.VARIANTS: "variant_annotations",
    Modality.VARIANT_ANNOTATIONS: "variant_annotations",
    Modality.PHENOTYPE: "phenotype",
}


class ProjectNotFoundError(RuntimeError):
    pass


@dataclass
class Project:
    root: Path
    config: ProjectConfig = field(default_factory=ProjectConfig)
    data: DataConfig = field(default_factory=DataConfig)
    resources: ResourcesConfig = field(default_factory=ResourcesConfig)

    # ------------------------------------------------------------------ loading
    @classmethod
    def find_root(cls, start: Path | None = None) -> Path:
        here = (start or Path.cwd()).resolve()
        for candidate in (here, *here.parents):
            if (candidate / PROJECT_FILE).exists():
                return candidate
        raise ProjectNotFoundError(
            f"no {PROJECT_FILE} found in {here} or its parents; run `efgpp init .` first"
        )

    @classmethod
    def load(cls, root: Path | str | None = None) -> Project:
        root_path = cls.find_root(Path(root) if root else None)
        return cls(
            root=root_path,
            config=load_model(ProjectConfig, root_path / PROJECT_FILE),
            data=load_model(DataConfig, root_path / DATA_FILE),
            resources=load_model(ResourcesConfig, root_path / RESOURCES_FILE),
        )

    def save_data_config(self) -> None:
        # Re-validate before writing so a bad edit can never reach disk.
        self.data = DataConfig.model_validate(self.data.model_dump())
        write_yaml(self.root / DATA_FILE, self.data.to_yaml_dict(), header=DATA_YAML_HEADER)

    def save_project_config(self) -> None:
        write_yaml(self.root / PROJECT_FILE, self.config.to_yaml_dict())

    def save_resources_config(self) -> None:
        write_yaml(self.root / RESOURCES_FILE, self.resources.to_yaml_dict())

    # -------------------------------------------------------------------- paths
    def path(self, *parts: str) -> Path:
        return self.root.joinpath(*parts)

    def resolve(self, path: str | Path) -> Path:
        """Resolve a user path: absolute paths as-is, relative ones against the project root."""
        p = Path(path).expanduser()
        return p if p.is_absolute() else (self.root / p)

    def relative(self, path: Path) -> str:
        """Portable path string: relative to the root when inside it, absolute otherwise."""
        try:
            return path.resolve().relative_to(self.root.resolve()).as_posix()
        except ValueError:
            return path.resolve().as_posix()

    @property
    def registry_path(self) -> Path:
        return self.resolve(self.config.registry.database)

    @property
    def data_root(self) -> Path:
        return self.resolve(self.config.storage.generated_data_root)

    @property
    def resource_root(self) -> Path:
        return self.resolve(self.config.storage.resource_root)

    @property
    def work_root(self) -> Path:
        return self.resolve(self.config.storage.work_root)

    @property
    def report_root(self) -> Path:
        return self.resolve(self.config.storage.report_root)

    @property
    def logs_dir(self) -> Path:
        return self.path("logs")

    @property
    def software_dir(self) -> Path:
        """Everything EFGPP installs for this project: software/{bin,envs,opt,downloads}."""
        return self.path("software")

    @property
    def legacy_software_dirs(self) -> list[Path]:
        """Where earlier EFGPP versions installed tools (still searched)."""
        return [self.path(".efgpp")]

    @property
    def bin_dir(self) -> Path:
        return self.software_dir / "bin"

    @property
    def envs_dir(self) -> Path:
        return self.software_dir / "envs"

    def artifact_dir(self, origin: Origin, modality: Modality | str, *parts: str) -> Path:
        mod = Modality(modality) if not isinstance(modality, Modality) else modality
        folder = _MODALITY_DIR.get(mod, mod.value)
        d = self.data_root / origin.value / folder
        for p in parts:
            d = d / p
        d.mkdir(parents=True, exist_ok=True)
        return d

    def qc_dir(self, modality: str) -> Path:
        d = self.path("qc", modality)
        d.mkdir(parents=True, exist_ok=True)
        return d

    def phenotype_dir(self, phenotype_id: str) -> Path:
        d = self.path("phenotypes", phenotype_id)
        d.mkdir(parents=True, exist_ok=True)
        return d

    def config_checksums(self) -> dict[str, str | None]:
        return {
            name: file_checksum(self.root / name)
            for name in (PROJECT_FILE, DATA_FILE, RESOURCES_FILE)
        }


def init_project(root: Path, name: str | None = None, force: bool = False) -> Project:
    """Create the EFGPP directory layout and default configuration files."""
    from efgpp.workflow import write_workflow_templates

    root = root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if (root / PROJECT_FILE).exists() and not force:
        raise FileExistsError(f"{root / PROJECT_FILE} already exists (use --force to re-initialise)")
    for d in PROJECT_DIRS:
        (root / d).mkdir(parents=True, exist_ok=True)

    config = ProjectConfig()
    config.project.name = name or root.name
    project = Project(root=root, config=config)
    project.save_project_config()
    if not (root / DATA_FILE).exists() or force:
        project.save_data_config()
    if not (root / RESOURCES_FILE).exists() or force:
        project.save_resources_config()
    write_workflow_templates(project)

    gitignore = root / ".gitignore"
    if not gitignore.exists():
        gitignore.write_text(
            "# EFGPP generated content\n.efgpp/\nsoftware/\nlogs/\nwork/\nregistry/*.duckdb*\n*.tmp\n",
            encoding="utf-8",
        )

    from efgpp.data.registry import Registry

    with Registry.open(project) as reg:
        reg.set_meta("project_name", config.project.name)
    return project
