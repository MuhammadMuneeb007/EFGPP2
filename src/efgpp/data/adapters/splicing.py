"""Processed splicing matrices (participant x intron/junction/event, e.g. PSI or LeafCutter ratios)."""

from __future__ import annotations

from typing import ClassVar

from efgpp.constants import Modality
from efgpp.data.adapters.expression import OmicsAdapter


class SplicingAdapter(OmicsAdapter):
    modality = Modality.SPLICING
    recommended_metadata: ClassVar[tuple[str, ...]] = (
        "tissue", "platform", "feature_id_system", "genome_build", "measurement_type", "normalization",
    )
