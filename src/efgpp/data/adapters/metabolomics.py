"""Processed metabolomics matrices (participant x metabolite)."""

from __future__ import annotations

from typing import ClassVar

from efgpp.constants import Modality
from efgpp.data.adapters.expression import OmicsAdapter


class MetabolomicsAdapter(OmicsAdapter):
    modality = Modality.METABOLOMICS
    recommended_metadata: ClassVar[tuple[str, ...]] = (
        "tissue", "platform", "feature_id_system", "units", "normalization",
    )
