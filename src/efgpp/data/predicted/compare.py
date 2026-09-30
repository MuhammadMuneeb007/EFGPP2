"""QC only: observed vs genetically predicted data of the same modality (neither replaces the other).

For every observed standardized matrix and every predicted (wide) matrix of one modality:
shared participants, shared features (Ensembl version suffixes ignored) and the per-feature
Pearson correlation across participants who have both values.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

from efgpp.constants import Modality, Origin
from efgpp.data.artifacts import ArtifactStore
from efgpp.data.registry import Registry
from efgpp.project import Project


def _base(feature: str) -> str:
    return feature.split(".")[0] if feature.startswith("ENS") else feature


def compare_observed_predicted(project: Project, modality: Modality, min_participants: int = 10
                               ) -> list[dict[str, Any]]:
    from efgpp.data.adapters.expression import load_omics_artifact

    with Registry.open(project) as reg:
        store = ArtifactStore(reg)
        observed = [a for a in store.find(modality=modality.value, artifact_type="standardized")
                    if a.origin == Origin.OBSERVED]
        predicted = [a for a in store.find(modality=modality.value, artifact_type="predicted")
                     if a.origin == Origin.PREDICTED and a.format == "parquet"]
    out: list[dict[str, Any]] = []
    for o in observed:
        m = load_omics_artifact(Path(o.path))
        obs = pl.DataFrame(m.X, schema=[_base(f) for f in m.features], orient="row").with_columns(
            m.obs.get_column("participant_id")).filter(pl.col("participant_id").is_not_null())
        obs = obs.group_by("participant_id", maintain_order=True).first()  # first measurement per participant
        for p in predicted:
            pred = pl.read_parquet(p.path).filter(pl.col("participant_id").is_not_null())
            pred = pred.rename({c: _base(c) for c in pred.columns if c not in ("participant_id", "native_id")})
            features = sorted(set(obs.columns) & set(pred.columns) - {"participant_id", "native_id"})
            joined = obs.select("participant_id", *features).join(
                pred.select("participant_id", *features), on="participant_id", how="inner", suffix="__pred")
            rs = []
            for f in features:
                a = joined.get_column(f).cast(pl.Float64).to_numpy()
                b = joined.get_column(f"{f}__pred").cast(pl.Float64).to_numpy()
                ok = ~(np.isnan(a) | np.isnan(b))
                if ok.sum() >= min_participants and a[ok].std() > 0 and b[ok].std() > 0:
                    rs.append(float(np.corrcoef(a[ok], b[ok])[0, 1]))
            out.append({"observed": o.source_id, "predicted": p.artifact_name, "tissue": p.tissue,
                        "shared_participants": joined.height, "shared_features": len(features),
                        "features_correlated": len(rs), "median_r": float(np.median(rs)) if rs else None,
                        "features_r_above_0.1": sum(r > 0.1 for r in rs)})
    return out
