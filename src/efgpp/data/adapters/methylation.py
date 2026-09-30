"""Processed methylation matrices (participant x CpG)."""

from __future__ import annotations

from typing import ClassVar

import numpy as np

from efgpp.constants import Modality
from efgpp.data.adapters.expression import OmicsAdapter, OmicsMatrix
from efgpp.data.validation import ValidationReport


class MethylationAdapter(OmicsAdapter):
    modality = Modality.METHYLATION
    recommended_metadata: ClassVar[tuple[str, ...]] = (
        "tissue", "platform", "feature_id_system", "genome_build", "measurement_type", "normalization",
    )

    def modality_checks(self, m: OmicsMatrix) -> ValidationReport:
        rep = super().modality_checks(m)
        mt = (self.source.measurement_type or "").lower()
        if mt in ("beta", "beta_value") and m.X.size:
            finite = m.X[~np.isnan(m.X)]
            n_out = int(((finite < 0) | (finite > 1)).sum())
            if n_out:
                rep.fail(f"{n_out:,} beta values outside [0, 1]", n_out)
            else:
                rep.ok("beta values within [0, 1]")
        return rep
