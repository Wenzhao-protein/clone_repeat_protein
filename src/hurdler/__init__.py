"""HURDLER repeat-protein cloning toolkit."""

from .constants import PLASMIDS
from .design import DesignRequest, DesignResult, design_construct
from .index import PatternIndex
from .plasmid_reference import PlasmidProfile, PlasmidReference, VectorCutScheme
from .protein_index import ProteinPatternIndex
from .vector_design import (
    CompatibilityQuery,
    DesignRequestV2,
    DesignResultV2,
    DesignSelection,
    design_construct_v2,
    design_query,
)
from .matching import match_module, query_all_plasmids
from .rules import LEGACY_OPTIMIZED_V1, RuleProfile

__all__ = [
    "LEGACY_OPTIMIZED_V1",
    "PLASMIDS",
    "DesignRequest",
    "DesignResult",
    "PatternIndex",
    "ProteinPatternIndex",
    "PlasmidReference",
    "PlasmidProfile",
    "VectorCutScheme",
    "CompatibilityQuery",
    "DesignSelection",
    "DesignRequestV2",
    "DesignResultV2",
    "RuleProfile",
    "match_module",
    "query_all_plasmids",
    "design_construct",
    "design_query",
    "design_construct_v2",
]

__version__ = "0.6.0"
