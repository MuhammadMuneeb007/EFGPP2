"""Toolkits: complete, reproducible installs of the software used around EFGPP.

A toolkit is a declarative spec (conda packages, R packages, pip packages, GitHub
repositories, direct downloads, exposed commands). `efgpp setup toolkit <name>` installs it;
`efgpp setup check` verifies every item.

    perl        Perl (and Ensembl VEP through the `annotation` tools)
    r           R 4.3 + LDpred-2/SCT/lassosum2 (bigsnpr), lassosum, PANPRSnext, CTSLEB,
                RapidoPGS, EBPRS, R2BGLiMS (JAMPred), penRegSum (tlpSum), sim1000G, ...
    prs         PRS/GWAS binaries: PLINK 1.9/2, GCTA, GCTB 2.5, GEMMA, BOLT-LMM, vcftools, PRSice-2,
                LDAK + source repositories: PRScs, PRScsx, SDPR, DBSLMM, CTPR, NPS, XP-BLUP,
                smtpred, PRSbils, LDpred-funct, PolyFun, MTG2 (source)
    prs-python  Python 3.10 environment for PRScs/PRScsx/PRSbils/LDpred/VIPRS/PolyFun/Hail
    prs-py27    Python 2.7 environment for LDSC, AnnoPred, PleioPred
    simulation  simuPOP, msprime, tskit, stdpopsim
    metaxcan    MetaXcan / PrediXcan (individual-level, Python 3) pinned to a resolved commit
    methylation R >= 4.3 with data.table, dplyr, optparse, BEDMatrix (MIMOSA model conversion)
    spliceai    SpliceAI (Illumina) with TensorFlow, in its own environment
    hla         HIBAG + SNPRelate, gdsfmt, SeqArray (Bioconductor, via conda / BiocManager)
    predicted-omics = metaxcan + methylation (+ PLINK 2 for generic scores); models are resources

The package lists mirror what PRSTools needs (see Document.MD). Toolkits are phenotype-
independent infrastructure; running PRS methods belongs to the Representation layer.
"""

from __future__ import annotations

from dataclasses import dataclass, field

CRAN = "https://cloud.r-project.org"


@dataclass(frozen=True)
class Repo:
    """A GitHub repository fetched as an archive (no git needed) into <root>/opt/<name>."""

    name: str
    github: str  # owner/repo
    ref: str | None = None  # tag/branch; None = default branch
    # command name -> (interpreter env or "" for executables, path inside the repo)
    commands: dict[str, tuple[str, str]] = field(default_factory=dict)
    build: list[str] = field(default_factory=list)  # shell commands run inside the repo (env on PATH)
    ld_library_path: list[str] = field(default_factory=list)  # repo-relative library dirs for wrappers
    lib_env: str | None = None  # conda environment whose lib/ the wrapper puts on LD_LIBRARY_PATH
    note: str = ""


@dataclass(frozen=True)
class Download:
    """A single official file: zip/gz archive or raw executable."""

    name: str
    url: str
    kind: str  # zip | tgz | gz | raw
    members: dict[str, str] = field(default_factory=dict)  # file in archive -> command name
    note: str = ""


@dataclass(frozen=True)
class Toolkit:
    name: str
    description: str
    env: str | None = None  # conda environment (under <root>/envs/<env>)
    conda: list[str] = field(default_factory=list)
    pip: list[str] = field(default_factory=list)
    r_cran: list[str] = field(default_factory=list)
    r_bioc: list[str] = field(default_factory=list)
    r_github: list[str] = field(default_factory=list)  # owner/repo[@ref]
    r_archive: list[str] = field(default_factory=list)  # pkg@version from the CRAN archive (removed from CRAN)
    repos: list[Repo] = field(default_factory=list)
    downloads: list[Download] = field(default_factory=list)
    expose: list[str] = field(default_factory=list)  # env executables linked into <root>/bin
    checks: list[str] = field(default_factory=list)  # commands that must resolve afterwards
    requires: list[str] = field(default_factory=list)  # other toolkits installed first
    tools: list[str] = field(default_factory=list)  # data-layer tools (own conda env each) installed with it


R_BUILD_TOOLS = ["c-compiler", "cxx-compiler", "fortran-compiler", "make", "pkg-config", "gsl", "zlib"]

TOOLKITS: dict[str, Toolkit] = {t.name: t for t in (
    Toolkit(
        "perl", "Perl interpreter (VEP itself: `efgpp setup tools vep`)",
        env="perl", conda=["perl", "perl-app-cpanminus"], expose=["perl", "cpanm"], checks=["perl"],
    ),
    Toolkit(
        "r", "R (>= 4.3) with every R package used by PRSTools (LDpred-2, SCT, lassosum, PANPRS, CTSLEB, ...)",
        env="r",
        conda=["r-base>=4.3", *R_BUILD_TOOLS, "r-bigsnpr", "r-bigstatsr", "r-data.table", "r-magrittr",
               "r-remotes", "r-devtools", "r-optparse", "r-proc", "r-glmnet", "r-ranger", "r-caret",
               "r-superlearner", "r-dplyr", "r-matrix", "r-rcpp", "r-rcpparmadillo", "r-r.utils",
               "r-biocmanager", "r-genio", "r-ggplot2", "r-xgboost", "r-susier",
               "bioconductor-genomicranges",
               # GMP-based packages precompiled (permutations -> partitions needs libgmp);
               # readr/stringr for sim1000G
               "r-gmp", "r-partitions", "r-readr", "r-stringr", "r-mass"],
        r_cran=["permutations", "RapidoPGS", "PANPRSnext"],
        # sim1000G was removed from CRAN (2025-06-12, its dependency hapsim was archived)
        r_archive=["hapsim@0.31", "sim1000G@1.40"],
        r_github=["tshmak/lassosum", "jpattee/penRegSum", "andrewhaoyu/CTSLEB", "pjnewcombe/R2BGLiMS",
                  "cran/EBPRS"],
        expose=["R", "Rscript"],
        checks=["Rscript"],
    ),
    Toolkit(
        "prs", "PRS / GWAS / heritability binaries and method repositories",
        env="prs",
        conda=["plink", "plink2", "gcta", "gemma", "vcftools", "gsl", "armadillo", "make", "c-compiler",
               "cxx-compiler", "fortran-compiler", "zlib"],
        downloads=[
            Download("prsice", "https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip",
                     "zip", {"PRSice_linux": "PRSice", "PRSice.R": "PRSice.R"}),
            Download("gctb", "https://cnsgenomics.com/software/gctb/download/gctb_2.5.5_Linux.zip",
                     "zip", {"gctb": "gctb"}, note="GCTB 2.5.5 (SBayesR / SBayesRC)"),
            Download("ldak", "https://github.com/dougspeed/LDAK/raw/main/ldak6.3.linux", "raw", {"ldak6.3.linux": "ldak"}),
            # Bioconda's bolt-lmm cannot be installed (it needs an nlopt release that does not exist).
            Download("bolt-lmm", "https://alkesgroup.broadinstitute.org/BOLT-LMM/downloads/BOLT-LMM_v2.5.tar.gz",
                     "tgz", {"bolt": "bolt"}, note="BOLT-LMM v2.5 official static build"),
        ],
        repos=[
            Repo("PRScs", "getian107/PRScs", commands={"PRScs.py": ("prs-python", "PRScs.py")}),
            Repo("PRScsx", "getian107/PRScsx", commands={"PRScsx.py": ("prs-python", "PRScsx.py")}),
            Repo("PRSbils", "styvon/PRSbils", commands={"PRSbils.py": ("prs-python", "PRSbils.py")}),
            Repo("LDpred-funct", "carlaml/LDpred-funct",
                 commands={"ldpredfunct.py": ("prs-python", "ldpredfunct.py")}),
            Repo("polyfun", "omerwe/polyfun", commands={"polyfun.py": ("prs-python", "polyfun.py"),
                                                        "polypred.py": ("prs-python", "polypred.py")}),
            Repo("SDPR", "eldronzhou/SDPR", commands={"SDPR": ("", "SDPR")}, ld_library_path=["gsl/lib", "MKL/lib"],
                 note="prebuilt Linux binary with bundled GSL/MKL libraries"),
            Repo("DBSLMM", "biostat0903/DBSLMM", commands={"DBSLMM.R": ("r", "software/DBSLMM.R")},
                 note="the `dbslmm` executable is distributed on Google Drive only: download it from the "
                      "DBSLMM README into opt/DBSLMM/software/dbslmm"),
            Repo("CTPR", "wonilchung/CTPR", build=["tar -xzf ctpr_v1.1.tar.gz"],
                 commands={"ctpr": ("", "ctpr")}, lib_env="ctpr",
                 note="prebuilt Linux binary; libarmadillo.so.9 + OpenBLAS from software/envs/ctpr"),
            Repo("NPS", "sgchun/nps", ref="1.1.1", build=["make"],
                 commands={"nps-run_all_chroms.sh": ("", "run_all_chroms.sh")}),
            Repo("XP-BLUP", "tanglab/XP-BLUP", commands={"xpblup.sh": ("", "xpblup.sh")}),
            Repo("smtpred", "uqrmaie1/smtpred", commands={"smtpred.py": ("prs-python", "smtpred.py")}),
            Repo("ldsc", "bulik/ldsc", commands={"ldsc.py": ("prs-py27", "ldsc.py"),
                                                 "munge_sumstats.py": ("prs-py27", "munge_sumstats.py")}),
            Repo("AnnoPred", "yiminghu/AnnoPred", commands={"AnnoPred.py": ("prs-py27", "AnnoPred.py")}),
            Repo("PleioPred", "yiminghu/PleioPred", commands={"PleioPred.py": ("prs-py27", "PleioPred.py")}),
            Repo("mtg2", "honglee0707/mtg2",
                 note="source only; its Makefile needs Intel ifort + static MKL (make.sh). Build where Intel "
                      "oneAPI is available, or put the author's mtg2 binary in software/bin"),
        ],
        expose=["plink", "plink2", "gcta64", "gemma", "vcftools"],
        checks=["plink", "plink2", "gcta64", "gemma", "bolt", "PRSice", "gctb", "ldak"],
        requires=["prs-python", "prs-py27", "r", "ctpr-libs"],
    ),
    Toolkit(
        "ctpr-libs", "Shared libraries of the prebuilt CTPR binary (Armadillo 9, OpenBLAS)",
        env="ctpr", conda=["armadillo=9.900", "libopenblas"],
    ),
    Toolkit(
        "prs-python", "Python 3.10 for PRScs, PRScsx, PRSbils, LDpred, VIPRS, PolyFun, Hail (Java 11)",
        env="prs-python",
        conda=["python=3.10", "numpy<2", "scipy", "pandas", "h5py", "scikit-learn", "statsmodels",
               "bitarray", "tqdm", "pyarrow", "networkx", "openjdk=11", "c-compiler", "cxx-compiler"],
        pip=["ldpred", "viprs", "magenpy", "pandas-plink", "plinkio", "pyliftover", "pgenlib", "hail"],
        expose=[],
        checks=[],
    ),
    Toolkit(
        "prs-py27", "Python 2.7 for LDSC, AnnoPred and PleioPred (they do not run on Python 3)",
        env="prs-py27",
        conda=["python=2.7", "pip", "numpy", "scipy", "pandas", "bitarray", "h5py", "scikit-learn", "statsmodels",
               "c-compiler"],
        pip=["plinkio"],
    ),
    Toolkit(
        "simulation", "Genotype/phenotype simulation: simuPOP, msprime, tskit, stdpopsim (+ R sim1000G via `r`)",
        env="simulation",
        conda=["python=3.11", "simupop", "msprime", "tskit", "stdpopsim", "numpy", "pandas", "scipy"],
        checks=[],
    ),
)}

TOOLKITS.update({t.name: t for t in (
    Toolkit(
        "metaxcan", "MetaXcan / PrediXcan (individual-level prediction, Python 3; never S-PrediXcan here)",
        tools=["predixcan", "plink2"],
    ),
    Toolkit(
        "methylation", "R for MIMOSA model conversion (weights are converted once; no MWAS stack needed)",
        env="methylation",
        conda=["r-base>=4.3", "r-data.table", "r-dplyr", "r-optparse", "r-bedmatrix"],
        checks=[],
    ),
    Toolkit(
        "spliceai", "SpliceAI (Illumina) + TensorFlow in software/envs/spliceai (models CC BY-NC 4.0)",
        tools=["spliceai"],
    ),
    Toolkit(
        "hla", "HIBAG HLA imputation (Bioconductor) with SNPRelate, gdsfmt, SeqArray",
        env="hla",
        conda=["r-base>=4.3", "r-biocmanager", "bioconductor-hibag", "bioconductor-snprelate",
               "bioconductor-gdsfmt", "bioconductor-seqarray"],
        r_bioc=["HIBAG", "SNPRelate", "gdsfmt", "SeqArray"],
    ),
)})

ALIASES = {"all": ["perl", "r", "prs-python", "prs-py27", "prs", "simulation"], "ldpred2": ["r"],
           "predicted-omics": ["metaxcan", "methylation"]}


def resolve_names(names: list[str]) -> list[str]:
    """Expand aliases and dependencies; order so requirements install first."""
    wanted: list[str] = []

    def add(n: str) -> None:
        for m in ALIASES.get(n, [n]):
            if m not in TOOLKITS:
                raise KeyError(f"unknown toolkit {m!r}; choose from {', '.join([*TOOLKITS, *ALIASES])}")
            for dep in TOOLKITS[m].requires:
                add(dep)
            if m not in wanted:
                wanted.append(m)

    for n in names:
        add(n)
    return wanted
