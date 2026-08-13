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
from .progress import DesignProgressEvent
from .exact_dna_design import (
    EXACT_DNA_SCHEMA_VERSION,
    ExactDNAQuery,
    ExactDNAResult,
    ExactDNASelection,
    confirm_exact_dna_route,
    load_exact_dna_enzyme_catalog,
    query_exact_dna,
    write_exact_dna_outputs,
)

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
    "DesignProgressEvent",
    "EXACT_DNA_SCHEMA_VERSION",
    "ExactDNAQuery",
    "ExactDNASelection",
    "ExactDNAResult",
    "RuleProfile",
    "match_module",
    "query_all_plasmids",
    "design_construct",
    "design_query",
    "design_construct_v2",
    "query_exact_dna",
    "confirm_exact_dna_route",
    "load_exact_dna_enzyme_catalog",
    "write_exact_dna_outputs",
]

__version__ = "0.6.0"
