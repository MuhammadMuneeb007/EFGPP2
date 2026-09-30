"""Reference knowledge providers (origin: reference)."""

from __future__ import annotations

from efgpp.data.references.base import FetchedResource, ReferenceProvider
from efgpp.data.references.genome import GenomeProvider
from efgpp.data.references.gtex import GTExProvider
from efgpp.data.references.gwas_catalog import GWASCatalogProvider
from efgpp.data.references.molecular import (
    HIBAGProvider,
    MIMOSAProvider,
    OmicsPredProvider,
    PredictDBExpressionProvider,
    PredictDBProteinProvider,
    PredictDBSplicingProvider,
    SpliceAIResourcesProvider,
)
from efgpp.data.references.opentargets import OpenTargetsProvider
from efgpp.data.references.pgs_catalog import PGSCatalogProvider
from efgpp.data.references.predictdb import PredictDBProvider
from efgpp.data.references.variants import (
    AlphaGenomeProvider,
    AlphaMissenseProvider,
    ClinVarProvider,
    GnomADProvider,
)

PROVIDERS: dict[str, type[ReferenceProvider]] = {
    p.name: p
    for p in (
        GenomeProvider, ClinVarProvider, GnomADProvider, AlphaMissenseProvider, AlphaGenomeProvider,
        GWASCatalogProvider, PGSCatalogProvider, OpenTargetsProvider, GTExProvider, PredictDBProvider,
        PredictDBExpressionProvider, PredictDBSplicingProvider, PredictDBProteinProvider, OmicsPredProvider,
        MIMOSAProvider, HIBAGProvider, SpliceAIResourcesProvider,
    )
}

__all__ = ["PROVIDERS", "FetchedResource", "ReferenceProvider"]
