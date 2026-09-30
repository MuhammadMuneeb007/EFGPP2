"""Built-in simulation of participant-level data (origin: simulated).

A deliberately simple generator (independent loci in Hardy-Weinberg equilibrium, additive
liability model) so the Data module is simulation-ready and testable. The ground truth is
always written separately to data/simulated/truth/. Population-genetic realism
(`engine: stdpopsim`) is left to a future release.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
import yaml

from efgpp.config.data import CovariateSource, GenotypeSource, OmicsSource, PhenotypeSource
from efgpp.constants import Modality, Origin, PhenotypeType, StorageMode
from efgpp.project import Project

# BED genotype codes by count of the A1 (first .bim allele) allele: 2 -> 00, 1 -> 10, 0 -> 11.
_BED_CODE = {2: 0b00, 1: 0b10, 0: 0b11}
_MISSING = 0b01


def write_bed(prefix: Path, dosage: np.ndarray, iids: list[str], chroms: list[str], positions: list[int],
              a1: list[str], a2: list[str], sexes: list[int] | None = None,
              variant_ids: list[str] | None = None) -> Path:
    """Write a PLINK 1 BED/BIM/FAM fileset. `dosage` is samples x variants counts of A1
    (0/1/2, or -1 for missing)."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    n, m = dosage.shape
    code = np.full(dosage.shape, _MISSING, dtype=np.uint8)
    for count, bits in _BED_CODE.items():
        code[dosage == count] = bits
    nbytes = (n + 3) // 4
    padded = np.zeros((nbytes * 4, m), dtype=np.uint8)
    padded[:n] = code
    packed = (padded[0::4] | (padded[1::4] << 2) | (padded[2::4] << 4) | (padded[3::4] << 6)).astype(np.uint8)
    with open(str(prefix) + ".bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        fh.write(packed.T.tobytes())  # variant-major
    ids = variant_ids or [f"{c}:{p}:{r}:{a}" for c, p, r, a in zip(chroms, positions, a2, a1, strict=True)]
    with open(str(prefix) + ".bim", "w", encoding="utf-8", newline="\n") as fh:
        for c, vid, p, x, y in zip(chroms, ids, positions, a1, a2, strict=True):
            fh.write(f"{c}\t{vid}\t0\t{p}\t{x}\t{y}\n")
    sexes = sexes or [0] * n
    with open(str(prefix) + ".fam", "w", encoding="utf-8", newline="\n") as fh:
        for iid, sex in zip(iids, sexes, strict=True):
            fh.write(f"{iid}\t{iid}\t0\t0\t{sex}\t-9\n")
    return prefix


def read_bed_dosage(prefix: Path, n: int, m: int) -> np.ndarray:
    """Decode a BED file into A1 counts (-1 = missing). Used for tests and tiny cohorts."""
    raw = np.fromfile(str(prefix) + ".bed", dtype=np.uint8)[3:]
    nbytes = (n + 3) // 4
    raw = raw.reshape(m, nbytes)
    codes = np.stack([(raw >> s) & 0b11 for s in (0, 2, 4, 6)], axis=2).reshape(m, nbytes * 4)[:, :n].T
    out = np.full(codes.shape, -1, dtype=np.int8)
    out[codes == 0b00] = 2
    out[codes == 0b10] = 1
    out[codes == 0b11] = 0
    return out


@dataclass
class SimulatedCohort:
    iids: list[str]
    dosage: np.ndarray
    variant_ids: list[str]
    freqs: np.ndarray


BASES = np.array(list("ACGT"))


def simulate_genotypes(n: int, m: int, n_chrom: int, seed: int, prefix: Path) -> SimulatedCohort:
    rng = np.random.default_rng(seed)
    freqs = rng.uniform(0.05, 0.5, m)
    dosage = rng.binomial(2, freqs, size=(n, m)).astype(np.int8)
    chroms = [str(1 + (j * n_chrom) // m) for j in range(m)]
    positions: list[int] = []
    last: dict[str, int] = {}
    for c in chroms:
        last[c] = last.get(c, 10_000) + int(rng.integers(1_000, 20_000))
        positions.append(last[c])
    ref_idx = rng.integers(0, 4, m)
    alt_idx = (ref_idx + rng.integers(1, 4, m)) % 4
    a2, a1 = BASES[ref_idx].tolist(), BASES[alt_idx].tolist()
    iids = [f"SIM{i + 1:05d}" for i in range(n)]
    vids = [f"{c}:{p}:{r}:{a}" for c, p, r, a in zip(chroms, positions, a2, a1, strict=True)]
    sexes = rng.integers(1, 3, n).tolist()
    write_bed(prefix, dosage, iids, chroms, positions, a1, a2, sexes, vids)
    return SimulatedCohort(iids, dosage, vids, freqs)


def simulate_phenotype(cohort: SimulatedCohort, ptype: PhenotypeType, h2: float, n_causal: int,
                       prevalence: float, seed: int) -> tuple[np.ndarray, dict[str, object]]:
    rng = np.random.default_rng(seed)
    m = cohort.dosage.shape[1]
    causal = np.sort(rng.choice(m, size=min(n_causal, m), replace=False))
    x = cohort.dosage[:, causal].astype(float)
    x = (x - x.mean(0)) / np.where(x.std(0) > 0, x.std(0), 1)
    beta = rng.normal(0, 1, len(causal))
    g = x @ beta
    g = (g - g.mean()) / (g.std() or 1)
    liability = np.sqrt(h2) * g + np.sqrt(1 - h2) * rng.normal(0, 1, len(g))
    if ptype == PhenotypeType.BINARY:
        values = (liability > np.quantile(liability, 1 - prevalence)).astype(float)
    elif ptype == PhenotypeType.CONTINUOUS:
        values = liability
    else:
        values = np.digitize(liability, np.quantile(liability, [1 / 3, 2 / 3])).astype(float)
    truth = {
        "type": ptype.value, "heritability": h2, "prevalence": prevalence if ptype == PhenotypeType.BINARY else None,
        "causal_variants": [cohort.variant_ids[i] for i in causal], "effects": beta.tolist(),
    }
    return values, truth


def run_simulation(project: Project) -> list[str]:
    """Generate enabled simulated modalities and add them to data.yaml as sources."""
    cfg = project.data.simulation
    added: list[str] = []
    if not cfg.genotype.enabled:
        return added
    g = cfg.genotype
    gdir = project.artifact_dir(Origin.SIMULATED, Modality.GENOTYPE)
    truth_dir = project.artifact_dir(Origin.SIMULATED, Modality.TRUTH)
    cohort = simulate_genotypes(g.n_participants, g.n_variants, g.n_chromosomes, g.seed, gdir / "simulated")
    obs = project.data.observed

    def upsert(lst: list, src: object) -> None:  # type: ignore[type-arg]
        for i, s in enumerate(lst):
            if s.id == src.id:  # type: ignore[attr-defined]
                lst[i] = src
                return
        lst.append(src)
        added.append(src.id)  # type: ignore[attr-defined]

    upsert(obs.genotype, GenotypeSource(id="SIMGENO", path=project.relative(gdir / "simulated"), format="bed",
                                        mode=StorageMode.REFERENCE, origin=Origin.SIMULATED, genome_build="GRCh38"))
    truth: dict[str, object] = {"engine": g.engine, "seed": g.seed, "phenotypes": {}}
    if cfg.phenotype.enabled:
        p = cfg.phenotype
        cols = {"participant_id": cohort.iids}
        for k in range(p.count):
            values, t = simulate_phenotype(cohort, p.type, p.heritability, p.n_causal, p.prevalence, p.seed + k)
            name = f"sim_trait_{k + 1}"
            cols[name] = [str(int(v)) if p.type != PhenotypeType.CONTINUOUS else f"{v:.6f}" for v in values]
            truth["phenotypes"][name] = t  # type: ignore[index]
        pdir = project.artifact_dir(Origin.SIMULATED, Modality.PHENOTYPE)
        ppath = pdir / "simulated_phenotypes.csv"
        pl.DataFrame(cols).write_csv(ppath)
        for k in range(p.count):
            name = f"sim_trait_{k + 1}"
            upsert(obs.phenotypes, PhenotypeSource(
                id=f"SIMPH{k + 1:03d}", name=name, path=project.relative(ppath), value_column=name, type=p.type,
                origin=Origin.SIMULATED, mode=StorageMode.REFERENCE,
                levels=["0", "1", "2"] if p.type == PhenotypeType.ORDINAL else None,
            ))
    if cfg.covariates.enabled:
        rng = np.random.default_rng(cfg.covariates.seed)
        cdir = project.artifact_dir(Origin.SIMULATED, Modality.COVARIATES)
        cpath = cdir / "simulated_covariates.csv"
        pl.DataFrame({"participant_id": cohort.iids,
                      "age": rng.normal(55, 10, len(cohort.iids)).round(1),
                      "sex": rng.integers(0, 2, len(cohort.iids))}).write_csv(cpath)
        upsert(obs.covariates, CovariateSource(id="SIMCOV", path=project.relative(cpath), variables=["age", "sex"],
                                               origin=Origin.SIMULATED, mode=StorageMode.REFERENCE))
    if cfg.expression.enabled:
        rng = np.random.default_rng(cfg.expression.seed)
        edir = project.artifact_dir(Origin.SIMULATED, Modality.EXPRESSION)
        epath = edir / "simulated_expression.parquet"
        n_f = cfg.expression.n_features
        cis = rng.choice(cohort.dosage.shape[1], size=n_f)
        X = cohort.dosage[:, cis] * rng.normal(0.5, 0.2, n_f) + rng.normal(0, 1, (len(cohort.iids), n_f))
        frame = pl.DataFrame({"participant_id": cohort.iids}).hstack(
            pl.DataFrame(X, schema=[f"GENE{j + 1:04d}" for j in range(n_f)], orient="row"))
        frame.write_parquet(epath)
        truth["expression_cis_variants"] = {f"GENE{j + 1:04d}": cohort.variant_ids[c] for j, c in enumerate(cis)}
        upsert(obs.expression, OmicsSource(id="SIMRNA", path=project.relative(epath), origin=Origin.SIMULATED,
                                           tissue="simulated", measurement_type="simulated",
                                           feature_id_system="simulated", normalization="none",
                                           mode=StorageMode.REFERENCE))
    (truth_dir / "truth.yaml").write_text(yaml.safe_dump(truth, sort_keys=False), encoding="utf-8")
    project.save_data_config()
    return added
