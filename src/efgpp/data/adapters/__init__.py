"""Modality adapters implementing the DataAdapter interface."""

from __future__ import annotations

from efgpp.config.data import SourceBase
from efgpp.constants import Modality
from efgpp.data.adapters.base import AdapterContext, DataAdapter, QCResult
from efgpp.data.adapters.clinical import ClinicalAdapter
from efgpp.data.adapters.covariates import CovariateAdapter
from efgpp.data.adapters.expression import ExpressionAdapter, OmicsAdapter
from efgpp.data.adapters.genotype import GenotypeAdapter
from efgpp.data.adapters.metabolomics import MetabolomicsAdapter
from efgpp.data.adapters.methylation import MethylationAdapter
from efgpp.data.adapters.phenotype import PhenotypeAdapter
from efgpp.data.adapters.proteomics import ProteomicsAdapter

ADAPTERS: dict[Modality, type[DataAdapter]] = {
    Modality.GENOTYPE: GenotypeAdapter,
    Modality.PHENOTYPE: PhenotypeAdapter,
    Modality.COVARIATES: CovariateAdapter,
    Modality.EXPRESSION: ExpressionAdapter,
    Modality.METHYLATION: MethylationAdapter,
    Modality.PROTEOMICS: ProteomicsAdapter,
    Modality.METABOLOMICS: MetabolomicsAdapter,
    Modality.CLINICAL: ClinicalAdapter,
}


def adapter_for(ctx: AdapterContext, modality: Modality, source: SourceBase) -> DataAdapter:
    return ADAPTERS[modality](ctx, source)


__all__ = [
    "ADAPTERS",
    "AdapterContext",
    "ClinicalAdapter",
    "CovariateAdapter",
    "DataAdapter",
    "ExpressionAdapter",
    "GenotypeAdapter",
    "MetabolomicsAdapter",
    "MethylationAdapter",
    "OmicsAdapter",
    "PhenotypeAdapter",
    "ProteomicsAdapter",
    "QCResult",
    "adapter_for",
]
