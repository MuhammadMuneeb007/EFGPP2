"""Processed proteomics matrices (participant x protein)."""

from __future__ import annotations

from typing import ClassVar

from efgpp.constants import Modality
from efgpp.data.adapters.expression import OmicsAdapter


class ProteomicsAdapter(OmicsAdapter):
    modality = Modality.PROTEOMICS
    recommended_metadata: ClassVar[tuple[str, ...]] = (
        "tissue", "platform", "feature_id_system", "units", "normalization",
    )
