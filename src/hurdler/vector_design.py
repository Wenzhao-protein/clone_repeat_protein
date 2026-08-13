"""Shared controller for the annotation-aware HURDLER designer v2."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import time
from dataclasses import asdict, dataclass, field, fields, replace
from pathlib import Path
from typing import Any, Callable, Mapping, Protocol, Sequence

from .constants import validate_protein_sequence
from .design import (
    GBLOCK_MAX_BP,
    GENERATION_SCHEDULE,
    ConfirmedBoundary,
    _locked_construct,
    _score_fragments,
    analyze_repeat_sequence,
    confirm_repeat_boundary,
)
from .ga_optimization import (
    GAPopulationState,
    GA_SCORE_PROFILE,
    adaptive_copy_search,
    adjust_ga_score_profile_from_idt,
    genetic_refine_dna,
    load_restriction_sites,
)
from .dna_assembly import IDT_FRIENDLY_CLAMP, _type_iis_flank, load_enzyme_catalog
from .io import sha256_file, utc_now, write_json_atomic
from .optimization import (
    diversify_codons,
    load_codon_weights,
    recognition_site_count,
    reverse_complement,
    translate_dna,
)
from .paths import ProjectPaths
from .plasmid_reference import (
    PLASMID_REFERENCE_VERSION,
    PlasmidReferenceDatabase,
    SilencingDecision,
    VectorCutScheme,
    decide_cutter_silencing,
    load_plasmid_reference,
    retained_backbone_contains_site,
)
from .protein_index import (
    PROTEIN_INDEX_VERSION,
    ProteinPatternIndex,
    enumerate_protein_solutions,
)
from .progress import DesignRunControl, ProgressCallback, emit_progress


DESIGN_SCHEMA_VERSION_V2 = "vector-aware-hurdler-designer-v2"
INPUT_MODES = ("split", "full")
VALIDATION_MODES = ("none", "api", "batch")
ASSEMBLY_STRATEGIES = (
    "exact_reused_secondary_rdl",
    "legacy_adaptive_max",
    "single_exact",
)


class ComplexityScorer(Protocol):
    def score(self, name: str, sequence: str) -> dict[str, Any]: ...


CheckpointCallback = Callable[[Mapping[str, Any]], None]


class IDTScoringError(RuntimeError):
    """Fatal interactive scoring failure that must not become a GA rejection."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        response_sha256: str = "",
        invalid_rules: Sequence[str] = (),
    ) -> None:
        super().__init__(message)
        self.code = str(code)
        self.response_sha256 = str(response_sha256)
        self.invalid_rules = tuple(str(value) for value in invalid_rules)


@dataclass(frozen=True)
class CompatibilityQuery:
    schema_version: str
    input_mode: str
    sequence_id: str
    n_cap: str = ""
    repeat_module: str = ""
    c_cap: str = ""
    repeat_copies: int = 2
    full_protein_sequence: str = ""
    repeat_region_start: int | None = None
    repeat_region_end: int | None = None
    repeat_period: int | None = None
    site_i_allowlist: tuple[str, ...] = ()
    site_ii_allowlist: tuple[str, ...] = ()
    site_iii_allowlist: tuple[str, ...] = ()
    plasmid_allowlist: tuple[str, ...] = ()
    allow_left_cutter_in_hurdler_pair: bool = False
    allow_right_cutter_in_hurdler_pair: bool = False
    max_restoration_length_bp: int | None = None

    def __post_init__(self) -> None:
        if self.schema_version != DESIGN_SCHEMA_VERSION_V2:
            raise ValueError(f"schema_version must be {DESIGN_SCHEMA_VERSION_V2!r}")
        if self.input_mode not in INPUT_MODES:
            raise ValueError(f"input_mode must be one of {INPUT_MODES}")
        for name in ("n_cap", "repeat_module", "c_cap", "full_protein_sequence"):
            value = str(getattr(self, name) or "").replace(" ", "").replace("\n", "").upper()
            if value:
                value = validate_protein_sequence(value)
            object.__setattr__(self, name, value)
        for name in ("site_i_allowlist", "site_ii_allowlist", "site_iii_allowlist", "plasmid_allowlist"):
            value = getattr(self, name)
            normalized = tuple(part.strip() for part in value.split(",") if part.strip()) if isinstance(value, str) else tuple(str(part).strip() for part in value if str(part).strip())
            object.__setattr__(self, name, normalized)
        if self.max_restoration_length_bp is not None:
            if isinstance(self.max_restoration_length_bp, bool) or not isinstance(
                self.max_restoration_length_bp, int
            ):
                raise ValueError("max_restoration_length_bp must be a non-negative integer or None")
            if self.max_restoration_length_bp < 0:
                raise ValueError("max_restoration_length_bp must be a non-negative integer or None")
        if self.input_mode == "split":
            if not self.repeat_module:
                raise ValueError("split input requires repeat_module")
            if int(self.repeat_copies) < 2:
                raise ValueError("repeat_copies must be at least two")
        else:
            if not self.full_protein_sequence:
                raise ValueError("full input requires full_protein_sequence")
            coordinates = (self.repeat_region_start, self.repeat_region_end, self.repeat_period)
            if any(value is not None for value in coordinates) and not all(value is not None for value in coordinates):
                raise ValueError("full input repeat coordinates must be provided together")

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "CompatibilityQuery":
        allowed = {item.name for item in fields(cls)}
        unknown = sorted(set(payload) - allowed)
        if unknown:
            raise ValueError("Unknown v2 CompatibilityQuery fields: " + ", ".join(unknown))
        return cls(**dict(payload))


@dataclass(frozen=True)
class DesignSelection:
    candidate_id: str
    profile_id: str
    scheme_id: str
    site_iii_enzyme: str
    advanced_pair_per_round: bool = False


@dataclass(frozen=True)
class DesignRequestV2:
    schema_version: str
    query: CompatibilityQuery
    selection: DesignSelection
    validation_mode: str = "none"
    assembly_strategy: str = "exact_reused_secondary_rdl"
    max_repeat_copies: int | None = None
    population_size: int = 16
    mutation_rate: float = 0.08
    crossover_rate: float = 0.75
    elite_fraction: float = 0.15
    seed: int = 42
    generation_schedule: tuple[int, ...] = GENERATION_SCHEDULE
    score_weights: Mapping[str, float] = field(default_factory=lambda: dict(GA_SCORE_PROFILE))
    auto_adjust_weights_from_idt: bool = True
    minimum_secondary_copies: int = 1
    maximum_secondary_copies: int | None = None
    max_idt_feedback_rounds: int = 100
    generations_per_feedback_round: int = 10
    elite_seed_count: int = 10
    auto_adjust_ga_parameters_from_idt: bool = True
    max_population_size: int = 256
    max_mutation_rate: float = 0.35
    max_crossover_rate: float = 0.95
    max_weight_multiplier: float = 1024.0
    max_purchase_bp: int = GBLOCK_MAX_BP
    ga_workers: int = 1

    def __post_init__(self) -> None:
        if self.schema_version != DESIGN_SCHEMA_VERSION_V2:
            raise ValueError(f"schema_version must be {DESIGN_SCHEMA_VERSION_V2!r}")
        if self.validation_mode not in VALIDATION_MODES:
            raise ValueError(f"validation_mode must be one of {VALIDATION_MODES}")
        if self.assembly_strategy not in ASSEMBLY_STRATEGIES:
            raise ValueError(f"assembly_strategy must be one of {ASSEMBLY_STRATEGIES}")
        if self.max_repeat_copies is not None and self.assembly_strategy != "legacy_adaptive_max":
            raise ValueError(
                "max_repeat_copies is only valid with assembly_strategy='legacy_adaptive_max'; "
                "exact RDL uses query.repeat_copies as the immutable final target"
            )
        if int(self.population_size) < 2:
            raise ValueError("population_size must be at least two")
        if isinstance(self.ga_workers, bool) or not isinstance(self.ga_workers, int):
            raise ValueError("ga_workers must be a positive integer")
        if self.ga_workers < 1:
            raise ValueError("ga_workers must be a positive integer")
        if isinstance(self.minimum_secondary_copies, bool) or not isinstance(
            self.minimum_secondary_copies, int
        ):
            raise ValueError("minimum_secondary_copies must be a positive integer")
        if int(self.minimum_secondary_copies) < 1:
            raise ValueError("minimum_secondary_copies must be at least one")
        if self.maximum_secondary_copies is not None:
            if isinstance(self.maximum_secondary_copies, bool) or not isinstance(
                self.maximum_secondary_copies, int
            ):
                raise ValueError("maximum_secondary_copies must be a positive integer or None")
            if self.maximum_secondary_copies < self.minimum_secondary_copies:
                raise ValueError(
                    "maximum_secondary_copies cannot be smaller than minimum_secondary_copies"
                )
        if not 1 <= int(self.max_idt_feedback_rounds) <= 1000:
            raise ValueError("max_idt_feedback_rounds must be between 1 and 1000")
        if not 1 <= int(self.generations_per_feedback_round) <= 1000:
            raise ValueError("generations_per_feedback_round must be between 1 and 1000")
        if not 1 <= int(self.elite_seed_count) <= int(self.max_population_size):
            raise ValueError("elite_seed_count must be between one and max_population_size")
        if int(self.max_population_size) < int(self.population_size):
            raise ValueError("max_population_size cannot be smaller than population_size")
        if not float(self.mutation_rate) <= float(self.max_mutation_rate) <= 1.0:
            raise ValueError("max_mutation_rate must be between mutation_rate and one")
        if not float(self.crossover_rate) <= float(self.max_crossover_rate) <= 1.0:
            raise ValueError("max_crossover_rate must be between crossover_rate and one")
        if float(self.max_weight_multiplier) < 1.0:
            raise ValueError("max_weight_multiplier must be at least one")
        schedule = tuple(sorted({int(value) for value in self.generation_schedule}))
        if not schedule or schedule[-1] != 100:
            raise ValueError("generation_schedule must terminate at 100 generations")
        object.__setattr__(self, "generation_schedule", schedule)
        if int(self.max_purchase_bp) > 3000:
            raise ValueError("The limit applies per purchase fragment and cannot exceed 3000 bp")

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "DesignRequestV2":
        allowed = {item.name for item in fields(cls)}
        unknown = sorted(set(payload) - allowed)
        if unknown:
            raise ValueError("Unknown v2 DesignRequest fields: " + ", ".join(unknown))
        values = dict(payload)
        values["query"] = CompatibilityQuery.from_dict(values["query"])
        values["selection"] = DesignSelection(**values["selection"])
        return cls(**values)


@dataclass
class DesignResultV2:
    schema_version: str
    status: str
    message: str
    request: dict[str, Any]
    boundary_analysis: dict[str, Any] | None = None
    confirmed_boundary: dict[str, Any] | None = None
    protein_candidates: list[dict[str, Any]] = field(default_factory=list)
    vector_routes: list[dict[str, Any]] = field(default_factory=list)
    selected_route: dict[str, Any] | None = None
    primary_fragments: list[dict[str, Any]] = field(default_factory=list)
    secondary_fragments: list[dict[str, Any]] = field(default_factory=list)
    restoration_segments: list[dict[str, Any]] = field(default_factory=list)
    stop_rescue_records: list[dict[str, Any]] = field(default_factory=list)
    cloning_steps: list[dict[str, Any]] = field(default_factory=list)
    optimization_attempts: list[dict[str, Any]] = field(default_factory=list)
    ga_elite_candidates: list[dict[str, Any]] = field(default_factory=list)
    ga_parameter_history: list[dict[str, Any]] = field(default_factory=list)
    idt_feedback_history: list[dict[str, Any]] = field(default_factory=list)
    maximum_secondary_evidence: dict[str, Any] = field(default_factory=dict)
    rdl_plan: dict[str, Any] = field(default_factory=dict)
    intermediate_validations: list[dict[str, Any]] = field(default_factory=list)
    assembly_steps: list[dict[str, Any]] = field(default_factory=list)
    idt_audit: list[dict[str, Any]] = field(default_factory=list)
    final_plasmid: dict[str, Any] | None = None
    final_protein_sequence: str = ""
    final_dna_sequence: str = ""
    termination_reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class DesignRouteUniverse:
    """Unranked routes for one protein/boundary and cutter-reuse policy.

    RE role, Site-III, plasmid, and restoration-length selections are applied
    later by :func:`filter_route_universe`.  Keeping these cheap filters out of
    enumeration lets interactive clients update support counts without
    rebuilding the protein index or plasmid schemes.
    """

    boundary_analysis: dict[str, Any] | None = None
    confirmed_boundary: dict[str, Any] | None = None
    protein_candidates: list[dict[str, Any]] = field(default_factory=list)
    vector_routes: list[dict[str, Any]] = field(default_factory=list)
    final_protein_sequence: str = ""


def bundled_protein_index_dir() -> Path:
    return ProjectPaths.discover().root / "data" / "artifacts" / PROTEIN_INDEX_VERSION


def _query_boundary(query: CompatibilityQuery) -> tuple[dict[str, Any] | None, ConfirmedBoundary | None, list[str]]:
    if query.input_mode == "split":
        module = query.repeat_module
        units = [module for _ in range(max(2, int(query.repeat_copies)))]
        boundary = ConfirmedBoundary(
            repeat_region_start=len(query.n_cap) + 1,
            repeat_region_end=len(query.n_cap) + len(module) * len(units),
            period=len(module),
            repeat_count=len(units),
            unit_sequences=tuple(units),
            middle_unit_index=(len(units) - 1) // 2 + 1,
            middle_unit_start=len(query.n_cap) + ((len(units) - 1) // 2) * len(module) + 1,
            middle_unit_end=len(query.n_cap) + (((len(units) - 1) // 2) + 1) * len(module),
            middle_module=module,
            consensus_module=module,
            position_conservation=tuple(1.0 for _ in module),
            fixed_positions_1based=tuple(range(1, len(module) + 1)),
            variable_ranges_1based=(),
            n_terminal_flank=query.n_cap,
            c_terminal_flank=query.c_cap,
            confirmation_token="split-input-explicit",
        )
        return None, boundary, [module]
    analysis = analyze_repeat_sequence(query.full_protein_sequence, sequence_id=query.sequence_id)
    if query.repeat_region_start is None:
        return analysis.to_dict(), None, []
    boundary = confirm_repeat_boundary(
        query.full_protein_sequence,
        start=int(query.repeat_region_start),
        end=int(query.repeat_region_end),
        period=int(query.repeat_period),
    )
    # Full-sequence mode deliberately searches the exact heterogeneous repeat
    # region as a whole.  This exposes RE opportunities created by variant
    # units and never replaces them with a consensus or middle unit.
    region = query.full_protein_sequence[int(query.repeat_region_start) - 1:int(query.repeat_region_end)]
    return analysis.to_dict(), boundary, [region]


def _allow(value: str, allowlist: Sequence[str]) -> bool:
    return not allowlist or value in set(allowlist)


def _candidate_id(row: Mapping[str, Any]) -> str:
    token = "|".join(str(row.get(key, "")) for key in ("pair_id", "pattern_key", "direction", "site_i_position", "site_ii_position", "site_i_enzyme", "site_ii_enzyme"))
    return hashlib.sha256(token.encode()).hexdigest()[:20]


def _protein_candidates(
    units: Sequence[str], query: CompatibilityQuery, index: ProteinPatternIndex
) -> list[dict[str, Any]]:
    per_unit = [enumerate_protein_solutions(unit, index) for unit in units]
    if not per_unit or any(not values for values in per_unit):
        return []
    common_pairs = set(row["pair_id"] for row in per_unit[0])
    for values in per_unit[1:]:
        common_pairs &= {row["pair_id"] for row in values}
    rows: list[dict[str, Any]] = []
    # A full-sequence route retains one representative match per exact unit;
    # every residue remains unchanged in later construction.
    for row in per_unit[0]:
        if row["pair_id"] not in common_pairs:
            continue
        if not _allow(str(row["site_i_enzyme"]), query.site_i_allowlist) or not _allow(str(row["site_ii_enzyme"]), query.site_ii_allowlist):
            continue
        site_iii = [value for value in str(row["site_iii_enzymes"]).split(",") if value]
        allowed_iii = [value for value in site_iii if _allow(value, query.site_iii_allowlist)]
        if not allowed_iii:
            continue
        expanded = dict(row)
        expanded["candidate_id"] = _candidate_id(row)
        expanded["site_iii_options"] = allowed_iii
        expanded["matched_exact_unit_count"] = len(units)
        expanded["all_exact_units_supported"] = True
        expanded["unit_candidate_ids_json"] = json.dumps(
            {
                unit: _candidate_id(next(item for item in values if item["pair_id"] == row["pair_id"]))
                for unit, values in zip(units, per_unit, strict=True)
            },
            sort_keys=True,
        )
        rows.append(expanded)
    return rows


def _scheme_route(
    database: PlasmidReferenceDatabase,
    scheme: VectorCutScheme,
    candidate: Mapping[str, Any],
    query: CompatibilityQuery,
) -> dict[str, Any] | None:
    if not scheme.valid:
        return None
    profile = database.profile(scheme.profile_id)
    if query.plasmid_allowlist and profile.profile_id not in query.plasmid_allowlist:
        return None
    for key in ("site_i_recognition_site", "site_ii_recognition_site"):
        if retained_backbone_contains_site(scheme, str(candidate[key])):
            return None
    left_reused = str(candidate["site_i_enzyme"]) in set(scheme.left_cutter.enzyme_aliases) or str(candidate["site_ii_enzyme"]) in set(scheme.left_cutter.enzyme_aliases)
    right_reused = str(candidate["site_i_enzyme"]) in set(scheme.right_cutter.enzyme_aliases) or str(candidate["site_ii_enzyme"]) in set(scheme.right_cutter.enzyme_aliases)
    if left_reused and not query.allow_left_cutter_in_hurdler_pair:
        return None
    if right_reused and not query.allow_right_cutter_in_hurdler_pair:
        return None
    decisions: list[SilencingDecision] = []
    if left_reused:
        decisions.append(decide_cutter_silencing(database, profile.profile_id, scheme.left_cutter))
    if right_reused:
        decisions.append(decide_cutter_silencing(database, profile.profile_id, scheme.right_cutter))
    if any(not item.allowed for item in decisions):
        return None
    return {
        "candidate_id": candidate["candidate_id"],
        "profile_id": profile.profile_id,
        "reference_id": profile.reference_id,
        "scheme_id": scheme.scheme_id,
        "cut_scheme": f"{scheme.left_location}/{scheme.right_location}",
        "left_cutter": scheme.left_cutter.canonical_enzyme,
        "right_cutter": scheme.right_cutter.canonical_enzyme,
        "left_cutter_site": scheme.left_cutter.recognition_site,
        "right_cutter_site": scheme.right_cutter.recognition_site,
        "left_cutter_aliases": list(scheme.left_cutter.enzyme_aliases),
        "right_cutter_aliases": list(scheme.right_cutter.enzyme_aliases),
        "site_i_enzyme": candidate["site_i_enzyme"],
        "site_ii_enzyme": candidate["site_ii_enzyme"],
        "site_iii_options": candidate["site_iii_options"],
        "cutter_reuse": left_reused or right_reused,
        "left_cutter_reused": left_reused,
        "right_cutter_reused": right_reused,
        "silencing_decisions": [asdict(item) for item in decisions],
        "left_restoration_sequence": scheme.left_restoration_sequence,
        "right_restoration_sequence": scheme.right_restoration_sequence,
        "restoration_length_bp": len(scheme.left_restoration_sequence) + len(scheme.right_restoration_sequence),
        "retained_backbone_sha256": scheme.retained_backbone_sha256,
        "protein_rule_profile": "legacy-optimized-v1",
        "protein_index_version": PROTEIN_INDEX_VERSION,
        "plasmid_reference_version": PLASMID_REFERENCE_VERSION,
    }


def _query_request_payload(query: CompatibilityQuery) -> dict[str, Any]:
    payload = asdict(query)
    for key in ("site_i_allowlist", "site_ii_allowlist", "site_iii_allowlist", "plasmid_allowlist"):
        payload[key] = list(payload[key])
    return payload


def build_route_universe(
    query: CompatibilityQuery,
    *,
    protein_index: ProteinPatternIndex | None = None,
    protein_index_dir: str | Path | None = None,
    plasmid_database: PlasmidReferenceDatabase | None = None,
    plasmid_reference_path: str | Path | None = None,
) -> DesignRouteUniverse:
    """Enumerate all selectable routes while retaining cutter-policy constraints."""

    unrestricted = replace(
        query,
        site_i_allowlist=(),
        site_ii_allowlist=(),
        site_iii_allowlist=(),
        plasmid_allowlist=(),
        max_restoration_length_bp=None,
    )
    analysis, boundary, units = _query_boundary(unrestricted)
    if boundary is None:
        return DesignRouteUniverse(
            boundary_analysis=analysis,
        )
    index = protein_index or ProteinPatternIndex.load(protein_index_dir or bundled_protein_index_dir())
    candidates = _protein_candidates(units, unrestricted, index)
    database = plasmid_database or load_plasmid_reference(plasmid_reference_path)
    routes: list[dict[str, Any]] = []
    for candidate in candidates:
        for scheme in database.schemes:
            route = _scheme_route(database, scheme, candidate, unrestricted)
            if route is not None:
                routes.append(route)
    target = (
        query.n_cap + query.repeat_module * int(query.repeat_copies) + query.c_cap
        if query.input_mode == "split"
        else query.full_protein_sequence
    )
    return DesignRouteUniverse(
        boundary_analysis=analysis,
        confirmed_boundary=boundary.to_dict(),
        protein_candidates=candidates,
        vector_routes=routes,
        final_protein_sequence=target,
    )


def filter_route_universe(
    universe: DesignRouteUniverse,
    query: CompatibilityQuery,
) -> DesignResultV2:
    """Apply interactive selections and rank the remaining annotation-safe routes."""

    request_payload = _query_request_payload(query)
    if universe.confirmed_boundary is None:
        return DesignResultV2(
            DESIGN_SCHEMA_VERSION_V2,
            "needs_boundary_confirmation",
            "Confirm or edit the 1-based inclusive repeat boundary and period; the full input sequence will not be homogenized.",
            request_payload,
            boundary_analysis=universe.boundary_analysis,
        )

    candidates: list[dict[str, Any]] = []
    for candidate in universe.protein_candidates:
        if not _allow(str(candidate["site_i_enzyme"]), query.site_i_allowlist):
            continue
        if not _allow(str(candidate["site_ii_enzyme"]), query.site_ii_allowlist):
            continue
        site_iii_options = [
            value
            for value in candidate["site_iii_options"]
            if _allow(str(value), query.site_iii_allowlist)
        ]
        if not site_iii_options:
            continue
        selected = dict(candidate)
        selected["site_iii_options"] = site_iii_options
        candidates.append(selected)

    if not candidates:
        return DesignResultV2(
            DESIGN_SCHEMA_VERSION_V2,
            "no_hurdler_pair_match",
            "No allowed protein-level Site-I/Site-II RE pair matches every required repeat unit.",
            request_payload,
            boundary_analysis=universe.boundary_analysis,
            confirmed_boundary=universe.confirmed_boundary,
            final_protein_sequence=universe.final_protein_sequence,
        )

    candidates_by_id = {row["candidate_id"]: row for row in candidates}
    routes_before_restoration: list[dict[str, Any]] = []
    for route in universe.vector_routes:
        candidate = candidates_by_id.get(route["candidate_id"])
        if candidate is None:
            continue
        if query.plasmid_allowlist and route["profile_id"] not in query.plasmid_allowlist:
            continue
        selected = dict(route)
        selected["site_iii_options"] = list(candidate["site_iii_options"])
        routes_before_restoration.append(selected)

    routes = routes_before_restoration
    if query.max_restoration_length_bp is not None:
        routes = [
            row
            for row in routes
            if int(row["restoration_length_bp"]) <= query.max_restoration_length_bp
        ]

    # Restoration filtering precedes this final fallback.  A short cutter-reuse
    # route therefore remains eligible when every unmodified route is too long.
    if any(not row["cutter_reuse"] for row in routes):
        routes = [row for row in routes if not row["cutter_reuse"]]
    routes.sort(key=lambda row: (row["cutter_reuse"], row["restoration_length_bp"], row["profile_id"], row["scheme_id"], row["candidate_id"]))
    for rank, row in enumerate(routes, 1):
        row["rank"] = rank
    status = "compatible_unoptimized" if routes else "no_vector_route"
    if routes:
        message = "At least one protein RE pair and annotation-safe plasmid cut route is available; select a route before optimization."
    elif query.max_restoration_length_bp is not None and routes_before_restoration:
        message = (
            "Routes exist, but all require restoration longer than "
            f"{query.max_restoration_length_bp} bp."
        )
    else:
        message = "Protein-level pairs exist, but none is clean on a retained annotated plasmid backbone under the selected cutter policy."
    return DesignResultV2(
        DESIGN_SCHEMA_VERSION_V2,
        status,
        message,
        request_payload,
        boundary_analysis=universe.boundary_analysis,
        confirmed_boundary=universe.confirmed_boundary,
        protein_candidates=candidates,
        vector_routes=routes,
        final_protein_sequence=universe.final_protein_sequence,
    )


def design_query(
    query: CompatibilityQuery,
    *,
    protein_index: ProteinPatternIndex | None = None,
    protein_index_dir: str | Path | None = None,
    plasmid_database: PlasmidReferenceDatabase | None = None,
    plasmid_reference_path: str | Path | None = None,
) -> DesignResultV2:
    universe = build_route_universe(
        query,
        protein_index=protein_index,
        protein_index_dir=protein_index_dir,
        plasmid_database=plasmid_database,
        plasmid_reference_path=plasmid_reference_path,
    )
    return filter_route_universe(universe, query)


def _actual_fragments(
    sequence_id: str,
    dna: str,
    route: Mapping[str, Any],
    *,
    max_purchase_bp: int,
) -> list[dict[str, Any]]:
    prefix = str(route.get("left_restoration_sequence", "")) + str(route.get("left_cutter_site", ""))
    suffix = str(route.get("right_restoration_sequence", "")) + str(route.get("right_cutter_site", ""))
    available = max_purchase_bp - len(prefix) - len(suffix)
    if available <= 0:
        raise ValueError("Restoration/adaptor sequences alone exceed max_purchase_bp")
    chunks = [dna[start:start + available] for start in range(0, len(dna), available)] or [""]
    rows = []
    for ordinal, chunk in enumerate(chunks, 1):
        purchase = (prefix if ordinal == 1 else "") + chunk + (suffix if ordinal == len(chunks) else "")
        rows.append({
            "fragment_id": f"{sequence_id}_{ordinal:02d}",
            "fragment_index": ordinal,
            "purchase_sequence": purchase,
            "purchase_length_bp": len(purchase),
            "purchase_sha256": hashlib.sha256(purchase.encode()).hexdigest(),
            "product_type": "gBlock" if 125 <= len(purchase) <= 3000 else "duplexed_ultramer" if 20 <= len(purchase) < 125 else "unsupported",
            "product_length_valid": 20 <= len(purchase) <= 3000,
            "left_restoration_bp": len(prefix) if ordinal == 1 else 0,
            "right_restoration_bp": len(suffix) if ordinal == len(chunks) else 0,
            "core_start_bp": len(prefix) if ordinal == 1 else 0,
            "core_end_bp": (len(prefix) if ordinal == 1 else 0) + len(chunk),
            "core_length_bp": len(chunk),
        })
    return rows


def _route_and_candidate(
    query_result: DesignResultV2, selection: DesignSelection
) -> tuple[dict[str, Any], dict[str, Any]]:
    routes = [row for row in query_result.vector_routes if row["candidate_id"] == selection.candidate_id and row["profile_id"] == selection.profile_id and row["scheme_id"] == selection.scheme_id and selection.site_iii_enzyme in row["site_iii_options"]]
    if not routes:
        raise ValueError("DesignSelection is not present in the current vector-route table")
    candidates = [row for row in query_result.protein_candidates if row["candidate_id"] == selection.candidate_id]
    if not candidates:
        raise ValueError("Selected protein candidate is missing")
    route, candidate = dict(routes[0]), dict(candidates[0])
    route["site_iii_enzyme"] = selection.site_iii_enzyme
    return route, candidate


def _batch_fragments(path: Path, fragments: Sequence[Mapping[str, Any]]) -> dict[str, str]:
    path.mkdir(parents=True, exist_ok=True)
    csv_path, tsv_path, fasta_path = path / "idt_bulk_input.csv", path / "idt_bulk_input.tsv", path / "idt_bulk_input.fasta"
    for destination, delimiter in ((csv_path, ","), (tsv_path, "\t")):
        with destination.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter=delimiter)
            writer.writerow(["Name", "Sequence"])
            writer.writerows((row["fragment_id"], row["purchase_sequence"]) for row in fragments)
    fasta_path.write_text("".join(f">{row['fragment_id']}\n{row['purchase_sequence']}\n" for row in fragments))
    return {"csv": str(csv_path), "tsv": str(tsv_path), "fasta": str(fasta_path)}


def _circular_site_count(sequence: str, motif: str) -> int:
    motif = str(motif).upper()
    if not motif:
        return 0
    patterns = {motif, reverse_complement(motif)}
    extended = sequence.upper() + sequence.upper()[: len(motif) - 1]
    return len(
        {
            position
            for pattern in patterns
            for position in range(len(sequence))
            if extended.startswith(pattern, position)
        }
    )


def _simulate_circular_vector(
    database: PlasmidReferenceDatabase,
    route: Mapping[str, Any],
    cds: str,
    protein: str,
    candidate: Mapping[str, Any],
) -> dict[str, Any]:
    scheme = next(item for item in database.schemes if item.scheme_id == route["scheme_id"])
    profile = database.profile(str(route["profile_id"]))
    reference = database.reference(profile.reference_id)
    donor_arc = scheme.left_restoration_sequence + cds + scheme.right_restoration_sequence
    final = scheme.retained_backbone_sequence + donor_arc
    cds_start = len(scheme.retained_backbone_sequence) + len(scheme.left_restoration_sequence)
    if translate_dna(final[cds_start:cds_start + len(cds)]) != protein:
        raise AssertionError("Circular-vector simulation changed the target translation")
    mcs = (profile.mcs_start, profile.mcs_end)
    missing_protected: list[str] = []
    for feature in reference.features:
        if not feature.protected:
            continue
        if any(start < mcs[1] and end > mcs[0] for start, end in feature.intervals):
            # MCS is the explicit modification exception; broad features that
            # span it are validated at their retained flanks by restoration.
            continue
        original = "".join(reference.sequence[start:end] for start, end in feature.intervals)
        if original and original not in final and reverse_complement(original) not in final:
            missing_protected.append(f"{feature.feature_class}:{feature.label}")
    if missing_protected:
        raise AssertionError("Protected feature sequence lost during assembly: " + ", ".join(missing_protected))
    return {
        "profile_id": profile.profile_id,
        "reference_id": reference.reference_id,
        "source_sequence_sha256": reference.sequence_sha256,
        "scheme_id": scheme.scheme_id,
        "final_plasmid_sequence": final,
        "final_plasmid_length_bp": len(final),
        "final_plasmid_sha256": hashlib.sha256(final.encode()).hexdigest(),
        "circular": True,
        "coordinate_orientation": "expression_direction_rotated_at_right_cut",
        "cds_start_0based": cds_start,
        "cds_end_0based_exclusive": cds_start + len(cds),
        "cds_translation_exact": True,
        "retained_backbone_sha256": hashlib.sha256(scheme.retained_backbone_sequence.encode()).hexdigest(),
        "restoration_exact": True,
        "protected_feature_sequences_preserved": True,
        "missing_protected_features": [],
        "selected_site_i_count": _circular_site_count(final, str(candidate["site_i_recognition_site"])),
        "selected_site_ii_count": _circular_site_count(final, str(candidate["site_ii_recognition_site"])),
        "left_cutter_count": _circular_site_count(final, scheme.left_cutter.recognition_site),
        "right_cutter_count": _circular_site_count(final, scheme.right_cutter.recognition_site),
    }


def _evaluate_split_copy_count(
    request: DesignRequestV2,
    candidate: Mapping[str, Any],
    route: Mapping[str, Any],
    copy_count: int,
    generations: int,
    *,
    codon_weights: dict[str, float],
    recognition_sites: tuple[str, ...],
    score_profile: dict[str, float],
    idt_scorer: ComplexityScorer | None,
    aggregate_audit: list[dict[str, Any]],
    progress_callback: ProgressCallback | None = None,
    run_control: DesignRunControl | None = None,
) -> dict[str, Any]:
    query = replace(request.query, repeat_copies=int(copy_count))
    _analysis, boundary, _units = _query_boundary(query)
    assert boundary is not None
    protein = query.n_cap + query.repeat_module * int(copy_count) + query.c_cap
    try:
        dna, locked_positions, site_limits, _ = _locked_construct(
            protein, boundary, candidate, codon_weights
        )
    except Exception as exc:
        return {
            "passed": False,
            "terminal": False,
            "error": f"{type(exc).__name__}: {exc}",
            "ga_local_constraints_passed": False,
            "selected_pair_re_site_excess": None,
        }
    refined, metrics = genetic_refine_dna(
        dna,
        locked_positions=locked_positions,
        selected_site_limits=site_limits,
        recognition_sites=recognition_sites,
        codon_weights=codon_weights,
        seed=request.seed + int(copy_count) * 1000 + int(generations),
        population_size=request.population_size,
        generations=int(generations),
        score_profile=score_profile,
        mutation_rate=request.mutation_rate,
        crossover_rate=request.crossover_rate,
        elite_fraction=request.elite_fraction,
        progress_callback=progress_callback,
        progress_context={"fragment_kind": "legacy_whole_construct", "copies": copy_count},
        ga_workers=request.ga_workers,
        run_control=run_control,
    )
    if translate_dna(refined) != protein:
        raise AssertionError("Adaptive GA changed the exact target protein")
    fragments = _actual_fragments(
        f"{query.sequence_id}_{copy_count}copies",
        refined,
        route,
        max_purchase_bp=request.max_purchase_bp,
    )
    local_pass = bool(metrics["ga_local_constraints_passed"]) and all(
        row["product_length_valid"] for row in fragments
    )
    idt_pass = False
    scored_fragments = fragments
    audits: list[dict[str, Any]] = []
    if local_pass and request.validation_mode == "api":
        if idt_scorer is None:
            raise RuntimeError("API validation mode requires a configured scorer")
        idt_pass, scored_fragments, audits = _score_fragments(
            fragments,
            idt_scorer,
            safe_point=run_control.safe_point if run_control is not None else None,
        )
        for row in audits:
            row.update({"repeat_copies": copy_count, "ga_generations": generations})
        aggregate_audit.extend(audits)
        if not idt_pass and request.auto_adjust_weights_from_idt:
            for row in audits:
                updated, _changes = adjust_ga_score_profile_from_idt(score_profile, row)
                score_profile.clear()
                score_profile.update(updated)
    passed = local_pass and (
        request.validation_mode == "batch" or (request.validation_mode == "api" and idt_pass)
    )
    scores = [row.get("idt_complexity_score") for row in scored_fragments]
    return {
        "passed": passed,
        "terminal": False,
        "copies": copy_count,
        "generations": generations,
        "dna_sequence": refined,
        "protein_sequence": protein,
        "boundary": boundary.to_dict(),
        "fragments": scored_fragments,
        "ga_score": metrics["ga_score"],
        "ga_local_constraints_passed": local_pass,
        "selected_pair_re_site_excess": metrics["selected_pair_re_site_excess"],
        "repeated_re_site_excess": metrics["repeated_re_site_excess"],
        "idt_request_attempted": bool(audits),
        "idt_api_called": bool(audits),
        "idt_status": "passed" if idt_pass else "batch_not_called" if request.validation_mode == "batch" else "rejected" if audits else "not_scored_local_failure",
        "idt_explicit_pass": idt_pass if request.validation_mode == "api" else None,
        "idt_complexity_score": max((float(value) for value in scores if isinstance(value, (int, float))), default=None),
        "ga_score_profile_after_idt_json": json.dumps(score_profile, sort_keys=True),
    }


def _minimum_locked_copy_count(candidate: Mapping[str, Any], module_length: int) -> int:
    """Return the smallest repeat count containing both selected 3-AA windows."""
    farthest = max(int(candidate["site_i_position"]), int(candidate["site_ii_position"])) + 3
    return max(1, math.ceil(farthest / int(module_length)))


def _secondary_adapters(
    route: Mapping[str, Any],
    *,
    project_root: Path,
) -> tuple[str, str, dict[str, Any]]:
    """Create disposable Site-III adapters for the selected Site-I/II ends."""
    geometries, _plasmids = load_enzyme_catalog(
        project_root / "data" / "reference_output",
        artifact_dir=project_root / "output",
    )
    names = (
        str(route["site_i_enzyme"]),
        str(route["site_ii_enzyme"]),
        str(route["site_iii_enzyme"]),
    )
    missing = [name for name in names if name not in geometries]
    if missing:
        raise ValueError("Missing restriction-enzyme cut geometry: " + ", ".join(missing))
    site_i, site_ii, site_iii = (geometries[name] for name in names)
    left = _type_iis_flank(
        site_iii,
        left=True,
        overhang_sequence=site_i.ovhgseq,
        clamp_sequence=IDT_FRIENDLY_CLAMP,
    )
    right = _type_iis_flank(
        site_iii,
        left=False,
        overhang_sequence=site_ii.ovhgseq,
        clamp_sequence=IDT_FRIENDLY_CLAMP,
    )
    return left, right, {
        "site_iii_enzyme": site_iii.enzyme,
        "left_adapter_sequence": left,
        "right_adapter_sequence": right,
        "left_exposed_overhang": site_i.ovhgseq,
        "right_exposed_overhang": site_ii.ovhgseq,
        "adapters_removed_after_digest": True,
    }


def _exact_split_boundary(query: CompatibilityQuery, copies: int) -> ConfirmedBoundary:
    """Describe an exact split-input repeat block, including a one-copy primary.

    ``CompatibilityQuery`` deliberately requires at least two target repeats so
    HURDLER can scan the doubled module.  The exact RDL construction search has
    a different lower bound: its primary seed may contain one physical copy
    after the protein-level route has already been selected.  Building this
    boundary directly prevents that construction detail from weakening the
    public query validation rule.
    """
    copies = int(copies)
    if copies < 1:
        raise ValueError("An exact split boundary requires at least one repeat")
    module = query.repeat_module
    period = len(module)
    start = len(query.n_cap) + 1
    end = len(query.n_cap) + period * copies
    middle_index = (copies - 1) // 2
    middle_start = start + middle_index * period
    protein = query.n_cap + module * copies + query.c_cap
    return ConfirmedBoundary(
        repeat_region_start=start,
        repeat_region_end=end,
        period=period,
        repeat_count=copies,
        unit_sequences=tuple(module for _ in range(copies)),
        middle_unit_index=middle_index + 1,
        middle_unit_start=middle_start,
        middle_unit_end=middle_start + period - 1,
        middle_module=module,
        consensus_module=module,
        position_conservation=tuple(1.0 for _ in module),
        fixed_positions_1based=tuple(range(1, period + 1)),
        variable_ranges_1based=(),
        n_terminal_flank=query.n_cap,
        c_terminal_flank=query.c_cap,
        confirmation_token=hashlib.sha256(
            f"{protein}|{start}|{end}|{period}".encode()
        ).hexdigest(),
    )


def _strict_idt_audit(audits: Sequence[Mapping[str, Any]]) -> None:
    """Raise when a live response cannot be classified by the score-sum policy."""
    for row in audits:
        status = str(row.get("idt_status") or "unknown")
        response_sha = str(row.get("response_sha256") or "")
        if status == "api_failure":
            error_type = str(row.get("idt_error_type") or "IDT API failure")
            raise IDTScoringError(
                "idt_api_error",
                f"IDT complexity scoring failed after retry handling: {error_type}",
                response_sha256=response_sha,
            )
        score = row.get("idt_complexity_score")
        numeric = (
            isinstance(score, (int, float))
            and not isinstance(score, bool)
            and math.isfinite(float(score))
        )
        if not bool(row.get("idt_score_complete")) or not numeric:
            try:
                invalid = json.loads(str(row.get("idt_invalid_score_names_json") or "[]"))
            except (TypeError, ValueError, json.JSONDecodeError):
                invalid = []
            invalid = invalid if isinstance(invalid, list) else []
            detail = ", ".join(str(value) for value in invalid) or "missing aggregate score"
            raise IDTScoringError(
                "idt_score_error",
                f"IDT returned an unclassifiable complexity score ({detail})",
                response_sha256=response_sha,
                invalid_rules=tuple(str(value) for value in invalid),
            )


def _adapt_ga_parameters_from_idt(
    request: DesignRequestV2,
    score: float,
    *,
    population_size: int,
    mutation_rate: float,
    crossover_rate: float,
) -> tuple[int, float, float, dict[str, Any]]:
    """Apply the frozen score-tier parameter policy after an IDT rejection."""
    if score < 10:
        return population_size, mutation_rate, crossover_rate, {"tier": "passed"}
    if score < 20:
        population_factor, mutation_factor, crossover_delta, tier = 1.25, 1.10, 0.02, "10_to_20"
    elif score < 50:
        population_factor, mutation_factor, crossover_delta, tier = 1.50, 1.25, 0.05, "20_to_50"
    else:
        population_factor, mutation_factor, crossover_delta, tier = 2.0, 1.50, 0.10, "50_or_more"
    proposed_population = int(math.ceil(population_size * population_factor / 4.0) * 4)
    updated_population = min(int(request.max_population_size), proposed_population)
    updated_mutation = min(float(request.max_mutation_rate), mutation_rate * mutation_factor)
    updated_crossover = min(float(request.max_crossover_rate), crossover_rate + crossover_delta)
    return updated_population, updated_mutation, updated_crossover, {
        "tier": tier,
        "population_factor": population_factor,
        "mutation_factor": mutation_factor,
        "crossover_delta": crossover_delta,
    }


def _finite_number(value: Any) -> float | None:
    if (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    ):
        return float(value)
    return None


def _idt_feedback_guidance(
    audits: Sequence[Mapping[str, Any]],
    fragments: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Translate IDT rule geometry into core-DNA GA targets.

    API coordinates refer to the submitted purchase fragment, which can carry
    disposable Type-IIS adapters or vector-restoration sequence.  Guidance is
    clipped and translated to the optimizable protein-coding core so a reported
    location cannot accidentally mutate an adapter or protected vector base.
    """
    fragment_by_id = {str(row.get("fragment_id")): row for row in fragments}
    target_ranges: set[tuple[int, int]] = set()
    repeat_windows: list[dict[str, Any]] = []
    terminal_gc_windows: list[dict[str, Any]] = []
    avoid_segments: list[str] = []
    rule_targets: list[dict[str, Any]] = []
    repeat_threshold = 100.0

    for audit in audits:
        fragment = fragment_by_id.get(str(audit.get("fragment_id")))
        if fragment is None:
            continue
        purchase = str(fragment.get("purchase_sequence") or "").upper()
        core_start = int(fragment.get("core_start_bp", 0))
        core_end = int(fragment.get("core_end_bp", len(purchase)))
        core_length = max(0, core_end - core_start)
        try:
            details = json.loads(str(audit.get("idt_rule_details_json") or "[]"))
        except (TypeError, ValueError, json.JSONDecodeError):
            details = []
        for detail in details if isinstance(details, list) else []:
            if not isinstance(detail, dict):
                continue
            score = _finite_number(detail.get("score"))
            if score is None or score <= 0:
                continue
            name = str(detail.get("name") or "unnamed_rule")
            lowered = name.lower()
            threshold = _finite_number(detail.get("threshold_value"))
            repeated_segment = str(detail.get("repeated_segment") or "").upper()
            locations = [
                int(value)
                for key in ("forward_locations", "reverse_locations")
                for value in detail.get(key, [])
                if isinstance(value, (int, float)) and not isinstance(value, bool)
            ]
            detail_ranges: list[tuple[int, int]] = []

            if repeated_segment and set(repeated_segment) <= set("ACGT"):
                core = purchase[core_start:core_end]
                if repeated_segment in core or reverse_complement(repeated_segment) in core:
                    avoid_segments.append(repeated_segment)
                for location in locations:
                    clipped_start = max(core_start, int(location)) - core_start
                    clipped_end = min(core_end, int(location) + len(repeated_segment)) - core_start
                    if clipped_end > clipped_start:
                        detail_ranges.append((clipped_start, clipped_end))

            if "windowed repeat" in lowered:
                start_index = int(_finite_number(detail.get("start_index")) or 0)
                window_length = int(
                    _finite_number(detail.get("threshold_window_length")) or 90
                )
                clipped_start = max(core_start, start_index) - core_start
                clipped_end = min(core_end, start_index + window_length) - core_start
                if clipped_end > clipped_start:
                    window_threshold = threshold if threshold is not None else 88.0
                    repeat_windows.append(
                        {
                            "start": clipped_start,
                            "end": clipped_end,
                            "threshold": window_threshold,
                            "rule": name,
                        }
                    )
                    detail_ranges.append((clipped_start, clipped_end))

            if "overall repeat" in lowered and threshold is not None:
                repeat_threshold = min(repeat_threshold, threshold)
                # This global rule has no locations.  The complete coding core
                # is its actionable target; adapters remain immutable here.
                if core_length:
                    detail_ranges.append((0, core_length))

            terminal_end_value = int(
                _finite_number(detail.get("terminal_end")) or 0
            )
            if (
                "terminal" in lowered or terminal_end_value in {3, 5}
            ) and "gc" in lowered:
                terminal_end = terminal_end_value
                purchase_start, purchase_end = (
                    (0, min(60, len(purchase)))
                    if terminal_end == 5
                    else (max(0, len(purchase) - 60), len(purchase))
                )
                clipped_start = max(core_start, purchase_start) - core_start
                clipped_end = min(core_end, purchase_end) - core_start
                if clipped_end > clipped_start:
                    terminal_gc_windows.append(
                        {
                            "start": clipped_start,
                            "end": clipped_end,
                            "threshold": threshold if threshold is not None else 68.0,
                            "terminal_end": terminal_end,
                            "rule": name,
                        }
                    )
                    detail_ranges.append((clipped_start, clipped_end))

            target_ranges.update(detail_ranges)
            rule_targets.append(
                {
                    "name": name,
                    "score": score,
                    "actual_value": detail.get("actual_value"),
                    "threshold_value": threshold,
                    "repeated_segment": repeated_segment,
                    "core_ranges": [list(value) for value in detail_ranges],
                }
            )

    # Latest concrete segments are most useful; cap retained guidance so a
    # long 100-round run cannot accumulate an unbounded obsolete blacklist.
    unique_segments = list(dict.fromkeys(avoid_segments))[-16:]
    return {
        "schema_version": "idt-structured-ga-feedback-v1",
        "repeat_coverage_threshold": repeat_threshold,
        "repeat_windows": repeat_windows[-8:],
        "terminal_gc_windows": terminal_gc_windows[-4:],
        "avoid_segments": unique_segments,
        "target_ranges": [list(value) for value in sorted(target_ranges)][-32:],
        "hotspot_mutation_rate": 0.65,
        "repeat_aware_steps": 40_000 if repeat_threshold < 100.0 else 0,
        "rule_targets": rule_targets,
    }


def _merge_idt_feedback_guidance(
    previous: Mapping[str, Any],
    current: Mapping[str, Any],
) -> dict[str, Any]:
    """Retain recently discovered hotspots while adding the newest rules.

    IDT reports one worst repeat window or repeated segment at a time.  If the
    next round targets only that newest item, an older repaired hotspot may
    reappear.  Keep a bounded, de-duplicated working set across feedback
    rounds so the search makes cumulative progress.
    """
    def unique_rows(key: str, limit: int) -> list[Any]:
        rows = [*previous.get(key, []), *current.get(key, [])]
        deduplicated: dict[str, Any] = {}
        for row in rows:
            deduplicated[json.dumps(row, sort_keys=True)] = row
        return list(deduplicated.values())[-limit:]

    previous_threshold = _finite_number(previous.get("repeat_coverage_threshold"))
    current_threshold = _finite_number(current.get("repeat_coverage_threshold"))
    thresholds = [value for value in (previous_threshold, current_threshold) if value is not None]
    return {
        "schema_version": "idt-structured-ga-feedback-v1",
        "repeat_coverage_threshold": min(thresholds) if thresholds else 100.0,
        "repeat_windows": unique_rows("repeat_windows", 16),
        "terminal_gc_windows": unique_rows("terminal_gc_windows", 8),
        "avoid_segments": unique_rows("avoid_segments", 32),
        "target_ranges": unique_rows("target_ranges", 64),
        "hotspot_mutation_rate": max(
            float(previous.get("hotspot_mutation_rate", 0.0)),
            float(current.get("hotspot_mutation_rate", 0.0)),
        ),
        "repeat_aware_steps": max(
            int(previous.get("repeat_aware_steps", 0)),
            int(current.get("repeat_aware_steps", 0)),
        ),
        "rule_targets": unique_rows("rule_targets", 64),
    }


def _rdl_purchase_fragments(
    request: DesignRequestV2,
    route: Mapping[str, Any],
    *,
    fragment_kind: str,
    copies: int,
    dna_sequence: str,
    secondary_adapters: tuple[str, str],
) -> tuple[list[dict[str, Any]], bool]:
    if fragment_kind == "primary":
        fragments = _actual_fragments(
            f"{request.query.sequence_id}_primary_{copies}copies",
            dna_sequence,
            route,
            max_purchase_bp=request.max_purchase_bp,
        )
        return fragments, len(fragments) == 1
    left_adapter, right_adapter = secondary_adapters
    purchase = left_adapter + dna_sequence + right_adapter
    length = len(purchase)
    return [{
        "fragment_id": f"{request.query.sequence_id}_secondary_{copies}copies",
        "fragment_index": 1,
        "purchase_sequence": purchase,
        "purchase_length_bp": length,
        "purchase_sha256": hashlib.sha256(purchase.encode()).hexdigest(),
        "product_type": "gBlock" if 125 <= length <= 3000 else "duplexed_ultramer" if 20 <= length < 125 else "unsupported",
        "product_length_valid": 20 <= length <= request.max_purchase_bp,
        "left_adapter": left_adapter,
        "right_adapter": right_adapter,
        "adapter_bases_removed_after_digest": True,
        "core_start_bp": len(left_adapter),
        "core_end_bp": len(left_adapter) + len(dna_sequence),
        "core_length_bp": len(dna_sequence),
    }], True


def _rdl_fragment_attempt(
    request: DesignRequestV2,
    candidate: Mapping[str, Any],
    route: Mapping[str, Any],
    *,
    fragment_kind: str,
    copies: int,
    generations: int,
    codon_weights: dict[str, float],
    recognition_sites: tuple[str, ...],
    score_profile: dict[str, float],
    idt_scorer: ComplexityScorer | None,
    aggregate_audit: list[dict[str, Any]],
    secondary_adapters: tuple[str, str] = ("", ""),
    population_state: GAPopulationState | None = None,
    population_size: int | None = None,
    mutation_rate: float | None = None,
    crossover_rate: float | None = None,
    feedback_round: int = 1,
    excluded_idt_sha256: set[str] | None = None,
    idt_feedback_guidance: dict[str, Any] | None = None,
    progress_callback: ProgressCallback | None = None,
    run_control: DesignRunControl | None = None,
) -> dict[str, Any]:
    """Optimize one exact primary or reusable secondary purchase sequence."""
    started = time.monotonic()
    module = request.query.repeat_module
    if fragment_kind == "primary":
        boundary = _exact_split_boundary(request.query, copies)
        protein = request.query.n_cap + module * int(copies) + request.query.c_cap
        try:
            dna, locked_positions, site_limits, _absolute = _locked_construct(
                protein, boundary, candidate, codon_weights
            )
        except Exception as exc:
            return {
                "passed": False,
                "error": f"{type(exc).__name__}: {exc}",
                "fragment_kind": fragment_kind,
                "copies": copies,
                "generations": generations,
            }
    elif fragment_kind == "secondary":
        protein = module * int(copies)
        site_limits = {
            str(candidate.get("site_i_recognition_site") or ""): 0,
            str(candidate.get("site_ii_recognition_site") or ""): 0,
        }
        site_limits = {site: limit for site, limit in site_limits.items() if site}
        dna = diversify_codons(
            protein,
            codon_weights=codon_weights,
            site_limits=site_limits,
        )
        locked_positions = set()
    else:
        raise ValueError(f"Unknown RDL fragment kind: {fragment_kind}")

    active_population = int(population_size or request.population_size)
    active_mutation = float(request.mutation_rate if mutation_rate is None else mutation_rate)
    active_crossover = float(request.crossover_rate if crossover_rate is None else crossover_rate)
    refined, metrics = genetic_refine_dna(
        dna,
        locked_positions=locked_positions,
        selected_site_limits=site_limits,
        recognition_sites=recognition_sites,
        codon_weights=codon_weights,
        seed=request.seed + (1 if fragment_kind == "primary" else 2) * 1_000_000 + copies * 1000 + feedback_round,
        population_size=active_population,
        generations=int(generations),
        score_profile=score_profile,
        idt_feedback_guidance=idt_feedback_guidance,
        mutation_rate=active_mutation,
        crossover_rate=active_crossover,
        elite_fraction=request.elite_fraction,
        population_state=population_state,
        elite_seed_count=request.elite_seed_count,
        capture_population_state=True,
        progress_callback=progress_callback,
        progress_context={
            "fragment_kind": fragment_kind,
            "copies": int(copies),
            "feedback_round": int(feedback_round),
            "max_feedback_rounds": int(request.max_idt_feedback_rounds),
            "population_size": active_population,
            "mutation_rate": active_mutation,
            "crossover_rate": active_crossover,
        },
        ga_workers=request.ga_workers,
        run_control=run_control,
    )
    if translate_dna(refined) != protein:
        raise AssertionError("RDL fragment GA changed the exact target protein")

    next_population_state = metrics.pop("ga_population_state")
    elite_candidates = list(metrics.pop("ga_elite_candidates", []))
    excluded = excluded_idt_sha256 or set()
    selected_row: dict[str, Any] | None = None
    fragments: list[dict[str, Any]] = []
    single_purchase = False
    for elite in elite_candidates:
        elite_dna = str(elite["dna_sequence"])
        candidate_fragments, candidate_single = _rdl_purchase_fragments(
            request,
            route,
            fragment_kind=fragment_kind,
            copies=copies,
            dna_sequence=elite_dna,
            secondary_adapters=secondary_adapters,
        )
        locally_valid = (
            bool(elite.get("ga_local_constraints_passed"))
            and candidate_single
            and all(bool(row["product_length_valid"]) for row in candidate_fragments)
        )
        already_scored = any(
            str(row.get("purchase_sha256") or "") in excluded
            for row in candidate_fragments
        )
        if locally_valid and not (request.validation_mode == "api" and already_scored):
            selected_row = elite
            refined = elite_dna
            fragments = candidate_fragments
            single_purchase = candidate_single
            break
    if not fragments:
        fragments, single_purchase = _rdl_purchase_fragments(
            request,
            route,
            fragment_kind=fragment_kind,
            copies=copies,
            dna_sequence=refined,
            secondary_adapters=secondary_adapters,
        )

    local_pass = (
        bool(selected_row and selected_row["ga_local_constraints_passed"])
        and single_purchase
        and all(bool(row["product_length_valid"]) for row in fragments)
    )
    idt_pass = False
    audits: list[dict[str, Any]] = []
    scored_fragments = fragments
    if local_pass and request.validation_mode == "api":
        if idt_scorer is None:
            raise RuntimeError("API validation mode requires a configured IDT scorer")
        emit_progress(
            progress_callback,
            stage="idt",
            status="request_started",
            fragment_kind=fragment_kind,
            copies=int(copies),
            generations=int(generations),
            feedback_round=int(feedback_round),
            max_feedback_rounds=int(request.max_idt_feedback_rounds),
            population_size=active_population,
            mutation_rate=active_mutation,
            crossover_rate=active_crossover,
            elapsed_seconds=time.monotonic() - started,
        )
        idt_pass, scored_fragments, audits = _score_fragments(
            fragments,
            idt_scorer,
            safe_point=run_control.safe_point if run_control is not None else None,
        )
        for row in audits:
            row.update({
                "fragment_kind": fragment_kind,
                "repeat_copies": int(copies),
                "ga_generations": int(generations),
                "feedback_round": int(feedback_round),
            })
        aggregate_audit.extend(audits)
        _strict_idt_audit(audits)
        positive_rules: list[str] = []
        for row in audits:
            try:
                names = json.loads(str(row.get("idt_positive_score_names_json") or "[]"))
            except (TypeError, ValueError, json.JSONDecodeError):
                names = []
            if isinstance(names, list):
                positive_rules.extend(str(value) for value in names)
        emit_progress(
            progress_callback,
            stage="idt",
            status="request_completed",
            fragment_kind=fragment_kind,
            copies=int(copies),
            generations=int(generations),
            feedback_round=int(feedback_round),
            max_feedback_rounds=int(request.max_idt_feedback_rounds),
            population_size=active_population,
            mutation_rate=active_mutation,
            crossover_rate=active_crossover,
            idt_score=max(
                (float(row["idt_complexity_score"]) for row in audits),
                default=None,
            ),
            idt_positive_rules=tuple(sorted(set(positive_rules))),
            elapsed_seconds=time.monotonic() - started,
            details={
                "passed": bool(idt_pass),
                "scores": [row.get("idt_complexity_score") for row in scored_fragments],
            },
        )
    elif local_pass and request.validation_mode == "batch":
        idt_pass = False

    passed = local_pass and (
        request.validation_mode == "batch"
        or (request.validation_mode == "api" and idt_pass)
    )
    scores = [row.get("idt_complexity_score") for row in scored_fragments]
    return {
        "passed": bool(passed),
        "fragment_kind": fragment_kind,
        "copies": int(copies),
        "generations": int(generations),
        "dna_sequence": refined,
        "protein_sequence": protein,
        "fragments": scored_fragments,
        "ga_score": float(selected_row["ga_score"] if selected_row else metrics["ga_score"]),
        "ga_local_constraints_passed": bool(local_pass),
        "selected_pair_re_site_excess": int(selected_row["selected_pair_re_site_excess"] if selected_row else metrics["selected_pair_re_site_excess"]),
        "repeated_re_site_excess": int(selected_row["repeated_re_site_excess"] if selected_row else metrics["repeated_re_site_excess"]),
        "feedback_round": int(feedback_round),
        "ga_total_generations": int(next_population_state.total_generations),
        "ga_population_size": active_population,
        "ga_mutation_rate": active_mutation,
        "ga_crossover_rate": active_crossover,
        "population_state": next_population_state,
        "elite_candidates": elite_candidates,
        "novel_candidate_available": bool(selected_row),
        "idt_request_attempted": bool(audits),
        "idt_status": "passed" if idt_pass else "batch_not_called" if request.validation_mode == "batch" else "rejected" if audits else "not_scored_local_failure",
        "idt_explicit_pass": idt_pass if request.validation_mode == "api" else None,
        "idt_complexity_score": max(
            (float(value) for value in scores if isinstance(value, (int, float))),
            default=None,
        ),
        "ga_score_profile_after_idt_json": json.dumps(score_profile, sort_keys=True),
        "elapsed_seconds": time.monotonic() - started,
    }


def _run_fragment_schedule(
    request: DesignRequestV2,
    candidate: Mapping[str, Any],
    route: Mapping[str, Any],
    *,
    fragment_kind: str,
    copies: int,
    codon_weights: dict[str, float],
    recognition_sites: tuple[str, ...],
    score_profile: dict[str, float],
    idt_scorer: ComplexityScorer | None,
    aggregate_audit: list[dict[str, Any]],
    parameter_history: list[dict[str, Any]],
    feedback_history: list[dict[str, Any]],
    secondary_adapters: tuple[str, str],
    progress_callback: ProgressCallback | None,
    run_control: DesignRunControl | None = None,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    attempts: list[dict[str, Any]] = []
    population_state: GAPopulationState | None = None
    population_size = int(request.population_size)
    mutation_rate = float(request.mutation_rate)
    crossover_rate = float(request.crossover_rate)
    scored_hashes: set[str] = set()
    guidance: dict[str, Any] = {}
    rounds = 1 if request.validation_mode == "batch" else int(request.max_idt_feedback_rounds)
    base_profile = {
        **GA_SCORE_PROFILE,
        **{str(key): float(value) for key, value in request.score_weights.items()},
    }
    for feedback_round in range(1, rounds + 1):
        if run_control is not None:
            run_control.safe_point()
        audit_start = len(aggregate_audit)
        result = _rdl_fragment_attempt(
            request,
            candidate,
            route,
            fragment_kind=fragment_kind,
            copies=copies,
            generations=int(request.generations_per_feedback_round),
            codon_weights=codon_weights,
            recognition_sites=recognition_sites,
            score_profile=score_profile,
            idt_scorer=idt_scorer,
            aggregate_audit=aggregate_audit,
            secondary_adapters=secondary_adapters,
            population_state=population_state,
            population_size=population_size,
            mutation_rate=mutation_rate,
            crossover_rate=crossover_rate,
            feedback_round=feedback_round,
            excluded_idt_sha256=scored_hashes,
            idt_feedback_guidance=guidance,
            progress_callback=progress_callback,
            run_control=run_control,
        )
        population_state = result["population_state"]
        compact = {
            key: value
            for key, value in result.items()
            if key not in {
                "dna_sequence", "protein_sequence", "fragments",
                "population_state", "elite_candidates",
            }
        }
        attempts.append(compact)
        compact["elite_candidates"] = [
            {key: value for key, value in elite.items() if key != "dna_sequence"}
            for elite in result["elite_candidates"]
        ]
        for fragment in result.get("fragments", []):
            if result.get("idt_request_attempted"):
                scored_hashes.add(str(fragment.get("purchase_sha256") or ""))
        round_audits = aggregate_audit[audit_start:]
        feedback_row = {
            "fragment_kind": fragment_kind,
            "repeat_copies": int(copies),
            "feedback_round": int(feedback_round),
            "max_feedback_rounds": int(request.max_idt_feedback_rounds),
            "total_generations": int(result["ga_total_generations"]),
            "novel_candidate_available": bool(result["novel_candidate_available"]),
            "idt_status": result["idt_status"],
            "idt_complexity_score": result["idt_complexity_score"],
            "candidate_sha256": hashlib.sha256(str(result["dna_sequence"]).encode()).hexdigest(),
            "guidance_applied_json": json.dumps(guidance, sort_keys=True),
        }
        feedback_history.append(feedback_row)
        if result.get("passed"):
            return result, attempts
        old_parameters = {
            "population_size": population_size,
            "mutation_rate": mutation_rate,
            "crossover_rate": crossover_rate,
        }
        adjustments: list[dict[str, Any]] = []
        score = result.get("idt_complexity_score")
        if request.validation_mode == "api" and isinstance(score, (int, float)):
            next_guidance = _idt_feedback_guidance(
                round_audits,
                result.get("fragments", []),
            )
            next_guidance = _merge_idt_feedback_guidance(guidance, next_guidance)
            if request.auto_adjust_weights_from_idt:
                for row in round_audits:
                    updated, changes = adjust_ga_score_profile_from_idt(score_profile, row)
                    for key, value in updated.items():
                        ceiling = float(base_profile.get(key, value)) * float(request.max_weight_multiplier)
                        updated[key] = min(float(value), ceiling)
                    score_profile.clear()
                    score_profile.update(updated)
                    adjustments.extend(changes)
            if request.auto_adjust_ga_parameters_from_idt:
                population_size, mutation_rate, crossover_rate, tier = _adapt_ga_parameters_from_idt(
                    request,
                    float(score),
                    population_size=population_size,
                    mutation_rate=mutation_rate,
                    crossover_rate=crossover_rate,
                )
            else:
                tier = {"tier": "disabled"}
            guidance = next_guidance
        else:
            tier = {"tier": "no_numeric_idt_rejection"}
        parameter_row = {
            **feedback_row,
            **{f"old_{key}": value for key, value in old_parameters.items()},
            "new_population_size": population_size,
            "new_mutation_rate": mutation_rate,
            "new_crossover_rate": crossover_rate,
            "parameter_tier": tier["tier"],
            "parameter_policy_json": json.dumps(tier, sort_keys=True),
            "weight_adjustments_json": json.dumps(adjustments, sort_keys=True),
            "score_weights_json": json.dumps(score_profile, sort_keys=True),
            "next_guidance_json": json.dumps(guidance, sort_keys=True),
        }
        parameter_history.append(parameter_row)
        emit_progress(
            progress_callback,
            stage="feedback",
            status="parameters_adjusted" if score is not None else "no_novel_candidate",
            fragment_kind=fragment_kind,
            copies=int(copies),
            feedback_round=int(feedback_round),
            max_feedback_rounds=int(request.max_idt_feedback_rounds),
            generations=int(request.generations_per_feedback_round),
            generation=int(request.generations_per_feedback_round),
            ga_score=float(result["ga_score"]),
            population_size=population_size,
            mutation_rate=mutation_rate,
            crossover_rate=crossover_rate,
            idt_score=float(score) if isinstance(score, (int, float)) else None,
            details={"parameter_tier": tier["tier"]},
        )
    return None, attempts


def _assemble_reused_secondary(
    request: DesignRequestV2,
    candidate: Mapping[str, Any],
    *,
    primary_copies: int,
    secondary_copies: int,
    rdl_rounds: int,
    primary_dna: str,
    secondary_dna: str,
) -> tuple[str, list[dict[str, Any]]]:
    module_length_bp = len(request.query.repeat_module) * 3
    insert_at = (len(request.query.n_cap) + int(primary_copies) * len(request.query.repeat_module)) * 3
    current = primary_dna
    validations: list[dict[str, Any]] = []
    limits = {
        str(candidate.get("site_i_recognition_site") or ""): 1,
        str(candidate.get("site_ii_recognition_site") or ""): 0,
    }
    limits = {site: limit for site, limit in limits.items() if site}
    for round_index in range(1, int(rdl_rounds) + 1):
        current = current[:insert_at] + secondary_dna + current[insert_at:]
        insert_at += len(secondary_dna)
        copies = int(primary_copies) + round_index * int(secondary_copies)
        expected = request.query.n_cap + request.query.repeat_module * copies + request.query.c_cap
        observed = translate_dna(current)
        site_counts = {site: recognition_site_count(current, site) for site in limits}
        passed = observed == expected and all(site_counts[site] <= limits[site] for site in limits)
        validations.append({
            "round": round_index,
            "primary_copies": int(primary_copies),
            "secondary_copies": int(secondary_copies),
            "result_copy_count": copies,
            "expected_length_bp": len(expected) * 3,
            "observed_length_bp": len(current),
            "translation_exact": observed == expected,
            "selected_site_counts": site_counts,
            "selected_site_limits": limits,
            "selected_pair_clean": all(site_counts[site] <= limits[site] for site in limits),
            "same_secondary_sha256": hashlib.sha256(secondary_dna.encode()).hexdigest(),
            "passed": passed,
        })
        if not passed:
            return current, validations
    if len(secondary_dna) != int(secondary_copies) * module_length_bp:
        raise AssertionError("Secondary DNA length does not equal its exact module count")
    return current, validations


def _design_exact_reused_secondary_rdl(
    request: DesignRequestV2,
    result: DesignResultV2,
    route: dict[str, Any],
    candidate: dict[str, Any],
    *,
    idt_scorer: ComplexityScorer | None,
    progress_callback: ProgressCallback | None,
    checkpoint_callback: CheckpointCallback | None,
    run_control: DesignRunControl | None = None,
) -> DesignResultV2:
    """Build an exact target from one primary and one reusable secondary."""
    started = time.monotonic()
    root = ProjectPaths.discover().root
    codon_weights = load_codon_weights(root / "data" / "reference_output" / "codon_usage.csv")
    recognition_sites = load_restriction_sites(root / "data" / "reference_output" / "restriction_enzyme.csv")
    target_copies = int(request.query.repeat_copies)
    module = request.query.repeat_module
    module_bp = len(module) * 3
    minimum_primary = _minimum_locked_copy_count(candidate, len(module))
    if target_copies < minimum_primary:
        result.status = "optimization_failed"
        result.message = (
            f"The selected route needs at least {minimum_primary} primary repeats, "
            f"but the exact target contains {target_copies}."
        )
        result.termination_reason = "target_below_selected_route_geometry_minimum"
        return result

    profile = {
        **GA_SCORE_PROFILE,
        **{str(key): float(value) for key, value in request.score_weights.items()},
    }
    emit_progress(
        progress_callback,
        stage="rdl_plan",
        status="started",
        message=f"Planning an exact {target_copies}-copy construct",
        copies=target_copies,
        elapsed_seconds=0.0,
    )
    if run_control is not None:
        run_control.safe_point()

    prefix_length = len(str(route.get("left_restoration_sequence", ""))) + len(str(route.get("left_cutter_site", "")))
    suffix_length = len(str(route.get("right_restoration_sequence", ""))) + len(str(route.get("right_cutter_site", "")))
    full_primary_bp = (
        len(request.query.n_cap + module * target_copies + request.query.c_cap) * 3
        + prefix_length
        + suffix_length
    )
    if full_primary_bp <= request.max_purchase_bp:
        direct, direct_attempts = _run_fragment_schedule(
            request,
            candidate,
            route,
            fragment_kind="primary",
            copies=target_copies,
            codon_weights=codon_weights,
            recognition_sites=recognition_sites,
            score_profile=profile,
            idt_scorer=idt_scorer,
            aggregate_audit=result.idt_audit,
            parameter_history=result.ga_parameter_history,
            feedback_history=result.idt_feedback_history,
            secondary_adapters=("", ""),
            progress_callback=progress_callback,
            run_control=run_control,
        )
        result.optimization_attempts.extend(
            [{**row, "component": "direct_primary"} for row in direct_attempts]
        )
        if direct is not None:
            result.final_protein_sequence = str(direct["protein_sequence"])
            result.final_dna_sequence = str(direct["dna_sequence"])
            result.primary_fragments = list(direct["fragments"])
            result.ga_elite_candidates = [
                {
                    **elite,
                    "fragment_kind": "primary",
                    "repeat_copies": target_copies,
                }
                for elite in direct["elite_candidates"]
            ]
            result.rdl_plan = {
                "strategy": "exact_reused_secondary_rdl",
                "target_repeat_copies": target_copies,
                "primary_repeat_copies": target_copies,
                "secondary_repeat_copies": 0,
                "secondary_reuse_count": 0,
                "equation": f"{target_copies} + 0 x 0 = {target_copies}",
                "initial_vector_insertions": 1,
                "rdl_rounds": 0,
                "final_copy_count_exact": True,
                "minimum_secondary_copies": int(request.minimum_secondary_copies),
                "minimum_secondary_bypassed_by_single_purchase": True,
            }
            result.status = "idt_accepted" if request.validation_mode == "api" else "optimized_unvalidated_batch"
            result.termination_reason = "whole_exact_target_single_purchase"
            result.message = (
                f"The exact {target_copies}-copy target fits one accepted primary purchase; no RDL round is required."
                if request.validation_mode == "api"
                else f"The exact {target_copies}-copy target fits one locally validated primary purchase; IDT was not called."
            )
            result.final_plasmid = _simulate_circular_vector(
                load_plasmid_reference(), route, result.final_dna_sequence,
                result.final_protein_sequence, candidate,
            )
            result.cloning_steps = [
                {"step": 1, "stage": "vector_digest", "action": f"Open {route['profile_id']} with {route['left_cutter']}/{route['right_cutter']}."},
                {"step": 2, "stage": "primary_installation", "action": f"Install the exact {target_copies}-copy primary once."},
                {"step": 3, "stage": "sequence_qc", "action": "Verify exact translation, selected-pair cleanliness, restoration, and protected annotations."},
            ]
            emit_progress(
                progress_callback,
                stage="rdl_plan",
                status="completed",
                message=result.message,
                copies=target_copies,
                elapsed_seconds=time.monotonic() - started,
            )
            return result

    left_adapter, right_adapter, adapter_evidence = _secondary_adapters(route, project_root=root)
    capacity = (request.max_purchase_bp - len(left_adapter) - len(right_adapter)) // module_bp
    minimum_secondary = int(request.minimum_secondary_copies)
    if capacity < minimum_secondary:
        result.status = "optimization_failed"
        result.message = (
            f"Site-III adapters permit at most {capacity} complete modules, below the "
            f"required minimum secondary size N={minimum_secondary}."
        )
        result.termination_reason = "minimum_secondary_exceeds_capacity"
        result.maximum_secondary_evidence = {
            "mathematical_capacity_copies": int(capacity),
            "required_minimum_copies": minimum_secondary,
            "maximum_verified_copies": 0,
            "proof": "minimum_secondary_exceeds_capacity",
            "adapter_evidence": adapter_evidence,
        }
        return result

    secondary_cache: dict[int, dict[str, Any]] = {}
    evaluated_secondary: dict[int, dict[str, Any] | None] = {}
    physical_exact_route_capacity = min(int(capacity), target_copies - minimum_primary)
    exact_route_capacity = physical_exact_route_capacity
    if request.maximum_secondary_copies is not None:
        exact_route_capacity = min(exact_route_capacity, int(request.maximum_secondary_copies))
    if exact_route_capacity < minimum_secondary:
        result.status = "optimization_failed"
        result.message = (
            f"The exact {target_copies}-copy target and selected primary geometry leave room for "
            f"at most {exact_route_capacity} modules in a reusable secondary, below N={minimum_secondary}."
        )
        result.termination_reason = "minimum_secondary_exceeds_exact_target_route"
        result.maximum_secondary_evidence = {
            "mathematical_capacity_copies": int(capacity),
            "exact_route_capacity_copies": int(exact_route_capacity),
            "required_minimum_copies": minimum_secondary,
            "maximum_verified_copies": 0,
            "proof": "minimum_secondary_exceeds_exact_target_route",
            "adapter_evidence": adapter_evidence,
        }
        return result

    def evaluate_secondary(copies: int) -> dict[str, Any] | None:
        copies = int(copies)
        if copies in evaluated_secondary:
            return evaluated_secondary[copies]
        accepted, attempts = _run_fragment_schedule(
            request,
            candidate,
            route,
            fragment_kind="secondary",
            copies=copies,
            codon_weights=codon_weights,
            recognition_sites=recognition_sites,
            score_profile=profile,
            idt_scorer=idt_scorer,
            aggregate_audit=result.idt_audit,
            parameter_history=result.ga_parameter_history,
            feedback_history=result.idt_feedback_history,
            secondary_adapters=(left_adapter, right_adapter),
            progress_callback=progress_callback,
            run_control=run_control,
        )
        result.optimization_attempts.extend(
            [
                {**row, "component": "maximum_secondary_search"}
                for row in attempts
            ]
        )
        evaluated_secondary[copies] = accepted
        if accepted is not None:
            secondary_cache[copies] = accepted
            if checkpoint_callback is not None:
                fragment = dict(accepted["fragments"][0])
                checkpoint_callback(
                    {
                        "event": "accepted_secondary",
                        "sequence_id": request.query.sequence_id,
                        "repeat_copies": copies,
                        "core_sequence": str(accepted["dna_sequence"]),
                        "core_length_bp": len(str(accepted["dna_sequence"])),
                        "purchase_sequence": str(fragment["purchase_sequence"]),
                        "purchase_length_bp": int(fragment["purchase_length_bp"]),
                        "purchase_sha256": str(fragment["purchase_sha256"]),
                        "idt_complexity_score": accepted.get("idt_complexity_score"),
                        "idt_status": accepted.get("idt_status"),
                        "idt_response_sha256": fragment.get("idt_response_sha256", ""),
                        "validation_mode": request.validation_mode,
                        "ga_total_generations": int(accepted["ga_total_generations"]),
                        "ga_score": float(accepted["ga_score"]),
                        "ga_weights": dict(profile),
                        "tested_secondary_copies": sorted(evaluated_secondary),
                        "tested_lengths": sorted(evaluated_secondary),
                        "failed_secondary_copies": sorted(
                            key for key, value in evaluated_secondary.items() if value is None
                        ),
                        "failure_reason": (
                            "Some tested copy counts exhausted the configured GA/IDT feedback rounds"
                            if any(value is None for value in evaluated_secondary.values())
                            else ""
                        ),
                        "query_fingerprint": hashlib.sha256(
                            json.dumps(asdict(request), sort_keys=True).encode()
                        ).hexdigest(),
                        "route_fingerprint": hashlib.sha256(
                            json.dumps(route, sort_keys=True).encode()
                        ).hexdigest(),
                    }
                )
        return accepted

    # The user-specified floor is a hard gate.  Only after N itself passes do
    # we use a binary probe followed by exact one-copy boundary advancement.
    minimum_result = evaluate_secondary(minimum_secondary)
    if minimum_result is None:
        result.status = "no_accepted_repeat_construct"
        result.message = (
            f"The required {minimum_secondary}-module secondary did not pass "
            f"within {request.max_idt_feedback_rounds} GA/IDT feedback rounds."
            if request.validation_mode == "api"
            else f"The required {minimum_secondary}-module secondary did not pass local Batch-mode constraints."
        )
        result.termination_reason = "minimum_secondary_not_idt_accepted"
        result.maximum_secondary_evidence = {
            "mathematical_capacity_copies": int(capacity),
            "exact_route_capacity_copies": int(exact_route_capacity),
            "required_minimum_copies": minimum_secondary,
            "maximum_verified_copies": 0,
            "proof": result.termination_reason,
            "feedback_round_limit": int(request.max_idt_feedback_rounds),
            "adapter_evidence": adapter_evidence,
        }
        return result

    low = minimum_secondary
    high = exact_route_capacity
    while low < high:
        if run_control is not None:
            run_control.safe_point()
        midpoint = (low + high + 1) // 2
        if evaluate_secondary(midpoint) is not None:
            low = midpoint
        else:
            high = midpoint - 1
    maximum_secondary = low
    best_secondary = secondary_cache[maximum_secondary]
    while maximum_secondary < exact_route_capacity:
        if run_control is not None:
            run_control.safe_point()
        next_copy = maximum_secondary + 1
        next_result = evaluate_secondary(next_copy)
        if next_result is None:
            break
        maximum_secondary = next_copy
        best_secondary = next_result
    if maximum_secondary == exact_route_capacity:
        search_reason = (
            "user_bounded_secondary_limit_reached"
            if request.maximum_secondary_copies is not None
            and exact_route_capacity < physical_exact_route_capacity
            else "exact_route_capacity_reached"
        )
    else:
        suffix = (
            f"failed_after_{request.max_idt_feedback_rounds}_feedback_rounds"
            if request.validation_mode == "api"
            else "failed_local_batch_validation"
        )
        search_reason = f"copy_{maximum_secondary + 1}_{suffix}"
    result.maximum_secondary_evidence = {
        "mathematical_capacity_copies": int(capacity),
        "physical_exact_route_capacity_copies": int(physical_exact_route_capacity),
        "exact_route_capacity_copies": int(exact_route_capacity),
        "requested_maximum_secondary_copies": request.maximum_secondary_copies,
        "required_minimum_copies": minimum_secondary,
        "maximum_verified_copies": int(maximum_secondary),
        "validation_mode": request.validation_mode,
        "proof": search_reason,
        "next_copy_failed_after_max_feedback_rounds": search_reason
        == f"copy_{maximum_secondary + 1}_failed_after_{request.max_idt_feedback_rounds}_feedback_rounds",
        "feedback_round_limit_per_copy": int(request.max_idt_feedback_rounds),
        "generations_per_feedback_round": int(request.generations_per_feedback_round),
        "adapter_evidence": adapter_evidence,
        "idt_verified": request.validation_mode == "api" and maximum_secondary > 0,
    }
    secondary_cache[int(maximum_secondary)] = best_secondary

    largest_route_secondary = min(int(maximum_secondary), target_copies - minimum_primary)
    primary_available_bp = request.max_purchase_bp - prefix_length - suffix_length
    if primary_available_bp <= 0:
        result.status = "optimization_failed"
        result.message = "Primary restoration/adaptor sequences leave no purchase capacity."
        result.termination_reason = "primary_adapter_capacity_below_one_module"
        return result

    # Search equations in increasing RDL-round order.  Within one round count,
    # prefer the longest already orderable secondary.  This implements the
    # promised "maximum secondary, minimum rounds" behavior while retaining a
    # deterministic fallback when that equation's independently optimized
    # primary is rejected.
    equations: list[tuple[int, int, int]] = []
    for rounds in range(1, target_copies - minimum_primary + 1):
        for secondary_copies in range(largest_route_secondary, minimum_secondary - 1, -1):
            primary_copies = target_copies - rounds * secondary_copies
            if primary_copies < minimum_primary:
                continue
            primary_core_bp = len(
                request.query.n_cap
                + module * primary_copies
                + request.query.c_cap
            ) * 3
            if primary_core_bp > primary_available_bp:
                continue
            equations.append((rounds, secondary_copies, primary_copies))

    for rounds, secondary_copies, primary_copies in equations:
        if run_control is not None:
            run_control.safe_point()
        equation_context = {
            "rdl_candidate_rounds": int(rounds),
            "rdl_candidate_secondary_copies": int(secondary_copies),
            "rdl_candidate_primary_copies": int(primary_copies),
        }
        secondary = secondary_cache.get(secondary_copies)
        if secondary is None:
            secondary = evaluate_secondary(secondary_copies)
            if secondary is None:
                continue
            secondary_cache[secondary_copies] = secondary

        primary, attempts = _run_fragment_schedule(
            request,
            candidate,
            route,
            fragment_kind="primary",
            copies=primary_copies,
            codon_weights=codon_weights,
            recognition_sites=recognition_sites,
            score_profile=profile,
            idt_scorer=idt_scorer,
            aggregate_audit=result.idt_audit,
            parameter_history=result.ga_parameter_history,
            feedback_history=result.idt_feedback_history,
            secondary_adapters=(left_adapter, right_adapter),
            progress_callback=progress_callback,
            run_control=run_control,
        )
        result.optimization_attempts.extend(
            [
                {**row, "component": "rdl_primary_candidate", **equation_context}
                for row in attempts
            ]
        )
        if primary is None:
            continue

        final_dna, validations = _assemble_reused_secondary(
            request,
            candidate,
            primary_copies=primary_copies,
            secondary_copies=secondary_copies,
            rdl_rounds=rounds,
            primary_dna=str(primary["dna_sequence"]),
            secondary_dna=str(secondary["dna_sequence"]),
        )
        result.intermediate_validations.extend(validations)
        final_protein = request.query.n_cap + module * target_copies + request.query.c_cap
        if not validations or not all(row["passed"] for row in validations):
            continue
        if translate_dna(final_dna) != final_protein:
            raise AssertionError("The selected RDL route did not produce the exact final protein")

        result.final_protein_sequence = final_protein
        result.final_dna_sequence = final_dna
        result.primary_fragments = [
            {**row, "fragment_role": "primary", "module_copies": int(primary_copies)}
            for row in primary["fragments"]
        ]
        result.secondary_fragments = [
            {
                **row,
                "fragment_role": "reusable_secondary",
                "module_copies": int(secondary_copies),
                "reuse_count": int(rounds),
            }
            for row in secondary["fragments"]
        ]
        result.ga_elite_candidates = [
            *[
                {
                    **elite,
                    "fragment_kind": "primary",
                    "repeat_copies": int(primary_copies),
                }
                for elite in primary["elite_candidates"]
            ],
            *[
                {
                    **elite,
                    "fragment_kind": "secondary",
                    "repeat_copies": int(secondary_copies),
                }
                for elite in secondary["elite_candidates"]
            ],
        ]
        result.rdl_plan = {
            "strategy": "exact_reused_secondary_rdl",
            "target_repeat_copies": target_copies,
            "minimum_secondary_copies": minimum_secondary,
            "maximum_secondary_copies": request.maximum_secondary_copies,
            "minimum_secondary_satisfied": int(secondary_copies) >= minimum_secondary,
            "minimum_secondary_bypassed_by_single_purchase": False,
            "minimum_primary_copies_for_selected_geometry": minimum_primary,
            "primary_repeat_copies": int(primary_copies),
            "secondary_repeat_copies": int(secondary_copies),
            "secondary_reuse_count": int(rounds),
            "equation": f"{primary_copies} + {rounds} x {secondary_copies} = {target_copies}",
            "initial_vector_insertions": 1,
            "rdl_rounds": int(rounds),
            "final_copy_count_exact": True,
            "secondary_purchase_sha256": result.secondary_fragments[0]["purchase_sha256"],
            "same_secondary_reused_every_round": True,
            "selected_pair_fixed_every_round": True,
        }
        result.status = "idt_accepted" if request.validation_mode == "api" else "optimized_unvalidated_batch"
        result.termination_reason = "exact_target_rdl_route_verified"
        result.message = (
            f"Exact {target_copies}-copy route verified: {primary_copies} primary + "
            f"{rounds} x {secondary_copies}-copy reusable secondary."
        )
        result.final_plasmid = _simulate_circular_vector(
            load_plasmid_reference(), route, final_dna, final_protein, candidate
        )
        result.cloning_steps = [
            {"step": 1, "stage": "vector_digest", "action": f"Open {route['profile_id']} with {route['left_cutter']}/{route['right_cutter']}; retain the annotated long backbone."},
            {"step": 2, "stage": "primary_installation", "action": f"Install the independently optimized {primary_copies}-copy primary containing both caps."},
            *[
                {
                    "step": 2 + round_index,
                    "stage": "rdl_secondary_round",
                    "action": (
                        f"RDL round {round_index}: reuse the same {secondary_copies}-copy secondary "
                        f"with {route['site_i_enzyme']}/{route['site_ii_enzyme']}/{route['site_iii_enzyme']} "
                        f"to reach {primary_copies + round_index * secondary_copies} copies."
                    ),
                }
                for round_index in range(1, rounds + 1)
            ],
            {"step": 3 + rounds, "stage": "sequence_qc", "action": "Verify exact 25-copy translation, selected-pair cleanliness, restoration, and protected annotations."},
        ]
        emit_progress(
            progress_callback,
            stage="rdl_plan",
            status="completed",
            message=result.message,
            copies=target_copies,
            elapsed_seconds=time.monotonic() - started,
            details=result.rdl_plan,
        )
        return result

    result.status = "no_accepted_repeat_construct"
    result.message = (
        f"No primary + reusable-secondary equation produced the exact {target_copies}-copy target "
        "after the required 100-generation attempts."
    )
    result.termination_reason = "all_exact_rdl_equations_exhausted"
    emit_progress(
        progress_callback,
        stage="rdl_plan",
        status="failed",
        message=result.message,
        copies=target_copies,
        elapsed_seconds=time.monotonic() - started,
    )
    return result


def design_construct_v2(
    request: DesignRequestV2,
    *,
    protein_index_dir: str | Path | None = None,
    plasmid_reference_path: str | Path | None = None,
    idt_scorer: ComplexityScorer | None = None,
    progress_callback: ProgressCallback | None = None,
    checkpoint_callback: CheckpointCallback | None = None,
    run_control: DesignRunControl | None = None,
) -> DesignResultV2:
    started = time.monotonic()
    emit_progress(
        progress_callback,
        stage="design",
        status="started",
        message="Validating the confirmed HURDLER route",
        elapsed_seconds=0.0,
    )
    if run_control is not None:
        run_control.safe_point()
    result = design_query(request.query, protein_index_dir=protein_index_dir, plasmid_reference_path=plasmid_reference_path)
    result.request = asdict(request)
    if result.status != "compatible_unoptimized":
        return result
    route, candidate = _route_and_candidate(result, request.selection)
    result.selected_route = route
    result.restoration_segments = [
        {"side": side, "sequence": route[f"{side}_restoration_sequence"], "length_bp": len(route[f"{side}_restoration_sequence"])}
        for side in ("left", "right") if route[f"{side}_restoration_sequence"]
    ]
    result.stop_rescue_records = [item for item in route["silencing_decisions"] if item["status"] == "stop_rescue_then_silence"]
    if request.validation_mode == "none":
        result.status = "compatible_unoptimized"
        result.message = "The route is selected, but no codon optimization was requested; no purchase DNA was emitted."
        return result
    if request.query.input_mode == "split" and request.assembly_strategy == "exact_reused_secondary_rdl":
        try:
            return _design_exact_reused_secondary_rdl(
                request,
                result,
                route,
                candidate,
                idt_scorer=idt_scorer,
                progress_callback=progress_callback,
                checkpoint_callback=checkpoint_callback,
                run_control=run_control,
            )
        except IDTScoringError as exc:
            result.status = exc.code
            result.message = str(exc)
            result.termination_reason = exc.code
            result.final_dna_sequence = ""
            result.final_protein_sequence = ""
            result.primary_fragments = []
            result.secondary_fragments = []
            result.idt_feedback_history.append(
                {
                    "status": exc.code,
                    "message": str(exc),
                    "response_sha256": exc.response_sha256,
                    "invalid_score_names": list(exc.invalid_rules),
                }
            )
            emit_progress(
                progress_callback,
                stage="idt",
                status="failed",
                message=str(exc),
                elapsed_seconds=time.monotonic() - started,
                details={
                    "error_code": exc.code,
                    "response_sha256": exc.response_sha256,
                    "invalid_score_names": list(exc.invalid_rules),
                },
            )
            return result
    root = ProjectPaths.discover().root
    weights = load_codon_weights(root / "data" / "reference_output" / "codon_usage.csv")
    recognition_sites = load_restriction_sites(root / "data" / "reference_output" / "restriction_enzyme.csv")
    if request.query.input_mode == "split" and request.assembly_strategy == "legacy_adaptive_max":
        upper = int(request.max_repeat_copies)
        if upper < 2:
            raise ValueError("max_repeat_copies must be at least two")
        adaptive_profile = {
            **GA_SCORE_PROFILE,
            **{str(key): float(value) for key, value in request.score_weights.items()},
        }
        best_copies, best, trace, reason = adaptive_copy_search(
            2,
            upper,
            short_generations=10,
            generation_schedule=request.generation_schedule,
            evaluate=lambda copies, generations: _evaluate_split_copy_count(
                request,
                candidate,
                route,
                copies,
                generations,
                codon_weights=weights,
                recognition_sites=recognition_sites,
                score_profile=adaptive_profile,
                idt_scorer=idt_scorer,
                aggregate_audit=result.idt_audit,
                progress_callback=progress_callback,
                run_control=run_control,
            ),
            progress_callback=progress_callback,
            progress_context={"fragment_kind": "legacy_whole_construct"},
            run_control=run_control,
        )
        result.optimization_attempts = trace
        result.termination_reason = reason
        if best_copies < 2 or best is None:
            result.status = "no_accepted_repeat_construct"
            result.message = "No repeat construct with at least two copies passed through the required 100-generation attempt."
            return result
        result.confirmed_boundary = dict(best["boundary"])
        result.final_protein_sequence = str(best["protein_sequence"])
        result.final_dna_sequence = str(best["dna_sequence"])
        result.primary_fragments = list(best["fragments"])
        result.selected_route.update(
            {
                "maximum_verified_repeat_copies": int(best_copies),
                "maximum_search_upper_bound": upper,
                "maximum_proof": reason,
                "next_copy_failed_at_100": reason == f"copy_{best_copies + 1}_failed_at_100",
                "final_ga_weights": dict(adaptive_profile),
            }
        )
        result.status = (
            "optimized_unvalidated_batch"
            if request.validation_mode == "batch"
            else "idt_accepted"
        )
        result.message = (
            f"Maximum verified repeat count is {best_copies}; Bulk Input files are unvalidated and no API/order call was made."
            if request.validation_mode == "batch"
            else f"Maximum verified repeat count is {best_copies}; every accepted purchase fragment passed live IDT score sum <10."
        )
        database = load_plasmid_reference(plasmid_reference_path)
        result.final_plasmid = _simulate_circular_vector(
            database,
            route,
            result.final_dna_sequence,
            result.final_protein_sequence,
            candidate,
        )
        result.cloning_steps = [
            {"step": 1, "stage": "vector_digest", "action": f"Digest {route['profile_id']} with {route['left_cutter']} and {route['right_cutter']}; retain the annotated long backbone."},
            {"step": 2, "stage": "primary_seed", "action": "Insert an independently optimized seed containing at least two repeat copies and every restoration/stop-rescue segment."},
            {"step": 3, "stage": "secondary_rounds", "action": f"Grow one repeat copy at a time with {route['site_i_enzyme']}/{route['site_ii_enzyme']}/{route['site_iii_enzyme']}; fixed-pair mode is active."},
            {"step": 4, "stage": "maximum_proof", "action": reason},
            {"step": 5, "stage": "sequence_qc", "action": "Verify exact translation, selected-pair site counts, restored vector sequence and protected annotations."},
        ]
        return result
    boundary_payload = result.confirmed_boundary or {}
    boundary = ConfirmedBoundary(
        **{
            **boundary_payload,
            "unit_sequences": tuple(boundary_payload["unit_sequences"]),
            "position_conservation": tuple(boundary_payload["position_conservation"]),
            "fixed_positions_1based": tuple(boundary_payload["fixed_positions_1based"]),
            "variable_ranges_1based": tuple(tuple(value) for value in boundary_payload["variable_ranges_1based"]),
        }
    )
    protein = result.final_protein_sequence
    try:
        dna, locked_positions, site_limits, _ = _locked_construct(protein, boundary, candidate, weights)
    except Exception as exc:
        result.status = "optimization_failed"
        result.message = f"Could not place the selected locked windows without changing translation: {type(exc).__name__}: {exc}"
        result.termination_reason = "locked_window_geometry_failure"
        return result
    profile = {**GA_SCORE_PROFILE, **{str(key): float(value) for key, value in request.score_weights.items()}}
    fragments: list[dict[str, Any]] = []
    for generations in request.generation_schedule:
        if run_control is not None:
            run_control.safe_point()
        refined, metrics = genetic_refine_dna(
            dna,
            locked_positions=locked_positions,
            selected_site_limits=site_limits,
            recognition_sites=recognition_sites,
            codon_weights=weights,
            seed=request.seed + generations,
            population_size=request.population_size,
            generations=generations,
            score_profile=profile,
            mutation_rate=request.mutation_rate,
            crossover_rate=request.crossover_rate,
            elite_fraction=request.elite_fraction,
            ga_workers=request.ga_workers,
            run_control=run_control,
        )
        if translate_dna(refined) != protein:
            raise AssertionError("GA changed the exact target protein")
        fragments = _actual_fragments(request.query.sequence_id, refined, route, max_purchase_bp=request.max_purchase_bp)
        local_pass = bool(metrics["ga_local_constraints_passed"]) and all(row["product_length_valid"] for row in fragments)
        audit: list[dict[str, Any]] = []
        idt_pass = False
        if local_pass and request.validation_mode == "api":
            if idt_scorer is None:
                raise RuntimeError("API validation mode requires a configured live or explicit mock IDT scorer")
            idt_pass, fragments, audit = _score_fragments(
                fragments,
                idt_scorer,
                safe_point=run_control.safe_point if run_control is not None else None,
            )
            result.idt_audit.extend(audit)
            if not idt_pass and request.auto_adjust_weights_from_idt:
                for row in audit:
                    profile, _ = adjust_ga_score_profile_from_idt(profile, row)
        elif local_pass and request.validation_mode == "batch":
            idt_pass = False
        result.optimization_attempts.append({
            "generations": generations,
            "dna_sha256": hashlib.sha256(refined.encode()).hexdigest(),
            "translation_exact": True,
            "local_constraints_passed": local_pass,
            "selected_pair_re_site_excess": metrics["selected_pair_re_site_excess"],
            "ga_weights_json": json.dumps(profile, sort_keys=True),
            "idt_passed": idt_pass if request.validation_mode == "api" else None,
        })
        dna = refined
        if local_pass and (request.validation_mode == "batch" or idt_pass):
            result.final_dna_sequence = refined
            result.primary_fragments = fragments
            result.status = "optimized_unvalidated_batch" if request.validation_mode == "batch" else "idt_accepted"
            result.message = (
                "Optimized exact sequences were exported for IDT Bulk Input; no HTTP validation or ordering was performed."
                if request.validation_mode == "batch"
                else "Every exact purchase fragment passed the live IDT score-sum <10 policy."
            )
            result.termination_reason = "batch_export_ready" if request.validation_mode == "batch" else "all_fragments_idt_accepted"
            break
    if not result.final_dna_sequence:
        result.status = "idt_rejected" if request.validation_mode == "api" else "optimization_failed"
        result.message = "No acceptable exact construct was obtained after the required 100-generation attempt."
        result.termination_reason = "failed_after_100_generations"
    else:
        database = load_plasmid_reference(plasmid_reference_path)
        result.final_plasmid = _simulate_circular_vector(
            database, route, result.final_dna_sequence, result.final_protein_sequence, candidate
        )
    result.cloning_steps = [
        {"step": 1, "stage": "vector_digest", "action": f"Digest {route['profile_id']} with {route['left_cutter']} and {route['right_cutter']}; retain the annotated long backbone."},
        {"step": 2, "stage": "primary_insert", "action": "Ligate the primary fragment, including every recorded restoration/stop-rescue segment."},
        {"step": 3, "stage": "hurdler_expansion", "action": f"Use Site I {route['site_i_enzyme']}, Site II {route['site_ii_enzyme']}, and Site III {route['site_iii_enzyme']}; preserve the selected pair by default across rounds."},
        {"step": 4, "stage": "sequence_qc", "action": "Verify the assembled CDS and all protected vector annotations by sequencing/digest simulation."},
    ]
    return result


def write_design_outputs_v2(result: DesignResultV2, output_dir: str | Path) -> dict[str, str]:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    if result.final_plasmid:
        from .design_artifacts import write_assembly_artifacts

        paths.update(write_assembly_artifacts(result, destination))
    summary = destination / "design_summary.json"
    write_json_atomic(result.to_dict(), summary)
    paths["design_summary_json"] = str(summary)
    for name, rows in (
        ("protein_candidates.csv", result.protein_candidates),
        ("vector_routes.csv", result.vector_routes),
        ("primary_fragments.csv", result.primary_fragments),
        ("secondary_fragments.csv", result.secondary_fragments),
        ("restoration_segments.csv", result.restoration_segments),
        ("stop_rescue_records.csv", result.stop_rescue_records),
        ("cloning_steps.csv", result.cloning_steps),
        ("ga_audit.csv", result.optimization_attempts),
        ("ga_elite_candidates.csv", result.ga_elite_candidates),
        ("ga_parameter_history.csv", result.ga_parameter_history),
        ("idt_feedback_history.csv", result.idt_feedback_history),
        ("rdl_intermediate_validations.csv", result.intermediate_validations),
        ("assembly_steps.csv", result.assembly_steps),
    ):
        path = destination / name
        columns = sorted({key for row in rows for key in row})
        with path.open("w", newline="") as handle:
            if columns:
                writer = csv.DictWriter(handle, fieldnames=columns)
                writer.writeheader()
                writer.writerows({key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list, tuple)) else value for key, value in row.items()} for row in rows)
        paths[name.replace(".", "_")] = str(path)
    if result.final_dna_sequence:
        fasta = destination / "optimized_construct.fasta"
        fasta.write_text(f">{result.request['query']['sequence_id']}\n{result.final_dna_sequence}\n")
        paths["optimized_construct_fasta"] = str(fasta)
        protein_fasta = destination / "optimized_construct_protein.fasta"
        protein_fasta.write_text(
            f">{result.request['query']['sequence_id']}_translated\n{result.final_protein_sequence}\n"
        )
        paths["optimized_construct_protein_fasta"] = str(protein_fasta)
    if result.ga_elite_candidates:
        elite_fasta = destination / "ga_elite_candidates.fasta"
        seen_elites: set[str] = set()
        records: list[str] = []
        for row in result.ga_elite_candidates:
            sequence = str(row.get("dna_sequence") or "")
            sequence_sha = str(row.get("dna_sha256") or hashlib.sha256(sequence.encode()).hexdigest())
            if not sequence or sequence_sha in seen_elites:
                continue
            seen_elites.add(sequence_sha)
            records.append(
                f">{row.get('fragment_kind', 'fragment')}_copies{row.get('repeat_copies', 'na')}_"
                f"rank{row.get('rank', 'na')}_{sequence_sha[:12]}\n{sequence}\n"
            )
        elite_fasta.write_text("".join(records))
        paths["ga_elite_candidates_fasta"] = str(elite_fasta)
    if result.rdl_plan:
        rdl_json = destination / "rdl_plan.json"
        write_json_atomic(result.rdl_plan, rdl_json)
        paths["rdl_plan_json"] = str(rdl_json)
        rdl_csv = destination / "rdl_plan.csv"
        with rdl_csv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=sorted(result.rdl_plan))
            writer.writeheader()
            writer.writerow({
                key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list, tuple)) else value
                for key, value in result.rdl_plan.items()
            })
        paths["rdl_plan_csv"] = str(rdl_csv)
    if result.final_plasmid:
        plasmid_fasta = destination / "final_plasmid.fasta"
        plasmid_fasta.write_text(
            f">{result.request['query']['sequence_id']}_final_circular_plasmid\n"
            f"{result.final_plasmid['final_plasmid_sequence']}\n"
        )
        paths["final_plasmid_fasta"] = str(plasmid_fasta)
        plasmid_genbank = destination / "final_plasmid.gb"
        final_step_paths = sorted(
            Path(value)
            for key, value in paths.items()
            if key.startswith("step") and key.endswith("_plasmid_gb")
        )
        if not final_step_paths:
            raise AssertionError("Accepted design did not produce a stepwise plasmid GenBank")
        plasmid_genbank.write_bytes(final_step_paths[-1].read_bytes())
        paths["final_plasmid_genbank"] = str(plasmid_genbank)
    all_fragments = [*result.primary_fragments, *result.secondary_fragments]
    if all_fragments and result.status in {"optimized_unvalidated_batch", "idt_accepted"}:
        paths.update({f"idt_bulk_{key}": value for key, value in _batch_fragments(destination, all_fragments).items()})
    audit = destination / "idt_audit.jsonl"
    audit.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in result.idt_audit))
    paths["idt_audit_jsonl"] = str(audit)
    manifest = {
        "schema_version": DESIGN_SCHEMA_VERSION_V2,
        "created_at": utc_now(),
        "status": result.status,
        "max_restoration_length_bp": result.request.get("query", {}).get(
            "max_restoration_length_bp"
        ),
        "minimum_secondary_copies": result.request.get("minimum_secondary_copies"),
        "maximum_secondary_copies": result.request.get("maximum_secondary_copies"),
        "ga_workers": result.request.get("ga_workers", 1),
        "query_fingerprint": hashlib.sha256(
            json.dumps(result.request.get("query", {}), sort_keys=True).encode()
        ).hexdigest(),
        "route_fingerprint": hashlib.sha256(
            json.dumps(result.selected_route or {}, sort_keys=True).encode()
        ).hexdigest(),
        "protein_index_version": PROTEIN_INDEX_VERSION,
        "plasmid_reference_version": PLASMID_REFERENCE_VERSION,
        "credentials_persisted": False,
        "ordering_performed": False,
        "files": {key: {"path": Path(value).name, "sha256": sha256_file(value)} for key, value in paths.items()},
    }
    manifest_path = destination / "run_manifest.json"
    write_json_atomic(manifest, manifest_path)
    paths["run_manifest_json"] = str(manifest_path)
    return paths
