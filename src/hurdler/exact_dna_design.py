"""Interactive, exact-sequence HURDLER design for arbitrary DNA and arrays.

Unlike the protein designer, this workflow never changes the requested target
sequence.  One-base latent restriction sites may be activated transiently by
donor-derived bases, but every complete route ends with the exact target DNA.
"""

from __future__ import annotations

import csv
import hashlib
import heapq
import io
import json
import math
import tempfile
import time
import zipfile
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Mapping, Protocol, Sequence

import pandas as pd
from Bio import SeqIO

from .constants import PLASMIDS
from .dna_assembly import (
    DEFAULT_GBLOCK_MAX_BP,
    DEFAULT_OLIGO_MAX_BP,
    DEFAULT_OLIGO_MIN_BP,
    PRIMER_PAIR_CORE_THRESHOLD_BP,
    EnzymeGeometry,
    TargetRecord,
    _product_type,
    _select_type_iis,
    _type_iis_flank,
    enumerate_active_latent_pairs,
    load_enzyme_catalog,
    plan_target,
    reverse_complement,
    scan_re_sites,
)
from .exact_dna_verification import verify_exact_dna_assembly
from .idt import IDT_SCORE_POLICY
from .io import sha256_file, utc_now, write_json_atomic
from .paths import ProjectPaths
from .plasmid_reference import (
    PLASMID_REFERENCE_VERSION,
    PlasmidReferenceDatabase,
    VectorCutScheme,
    decide_cutter_silencing,
    load_plasmid_reference,
    retained_backbone_contains_site,
)
from .progress import ProgressCallback, emit_progress


EXACT_DNA_SCHEMA_VERSION = "exact-dna-hurdler-designer-v1"
EXACT_DNA_INPUT_MODES = ("array", "exact")
EXACT_DNA_VALIDATION_MODES = ("none", "api", "batch")
EXACT_DNA_PAIR_POLICY = "fixed_then_variable"
DEFAULT_MAX_STATES = 10_000
DEFAULT_TIMEOUT_SECONDS = 600
DEFAULT_PATHS_PER_STATE = 3
DEFAULT_MAX_COMPLETE_ROUTES = 25


class ComplexityScorer(Protocol):
    def score(self, name: str, sequence: str) -> dict[str, Any]: ...


def _sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def parse_exact_dna_input(value: str, *, sequence_id: str = "") -> tuple[str, str]:
    """Parse one raw DNA sequence or one FASTA record without base changes."""
    text = str(value or "").strip()
    if not text:
        raise ValueError("DNA sequence is empty")
    resolved_id = str(sequence_id or "").strip()
    if text.startswith(">"):
        records = list(SeqIO.parse(io.StringIO(text), "fasta"))
        if len(records) != 1:
            raise ValueError("Exact-DNA input must contain exactly one FASTA record")
        if not resolved_id:
            resolved_id = str(records[0].id)
        text = str(records[0].seq)
    normalized = "".join(text.split()).upper()
    invalid = sorted(set(normalized) - set("ACGT"))
    if invalid:
        raise ValueError(
            "Exact-DNA input accepts A/C/G/T only and does not silently convert bases: "
            + "".join(invalid)
        )
    return resolved_id or "exact_dna_target", normalized


def _normalize_allowlist(value: Sequence[str] | str) -> tuple[str, ...]:
    if isinstance(value, str):
        values = value.split(",")
    else:
        values = value
    return tuple(dict.fromkeys(part.strip() for part in values if part.strip()))


def load_exact_dna_enzyme_catalog(
    *,
    reference_dir: str | Path | None = None,
    artifact_dir: str | Path | None = None,
) -> dict[str, EnzymeGeometry]:
    """Load broad Site-I/II roles plus the maintained Site-III adapter pool.

    Historical exact-DNA screening treated every maintained conventional RE
    as eligible for either HURDLER role. The optimized protein tables narrow
    those roles, but contain the maintained Type-IIS adapter list. Combining
    these sources preserves exact-DNA geometry without changing legacy callers
    of :func:`load_enzyme_catalog`.
    """
    paths = ProjectPaths.discover()
    reference = Path(reference_dir) if reference_dir else paths.reference_output
    role_artifacts = Path(artifact_dir) if artifact_dir else paths.root / "data/artifacts"
    geometries, _legacy = load_enzyme_catalog(reference, artifact_dir=role_artifacts)
    if not any(item.is_type_iis for item in geometries.values()):
        adapter_geometries, _unused = load_enzyme_catalog(
            reference, artifact_dir=paths.output
        )
        geometries.update(
            {
                name: geometry
                for name, geometry in adapter_geometries.items()
                if geometry.is_type_iis and geometry.site_iii_eligible
            }
        )
    return geometries


@dataclass(frozen=True)
class ExactDNAQuery:
    schema_version: str
    input_mode: str
    sequence_id: str = ""
    exact_dna: str = ""
    repeat_unit: str = ""
    spacer: str = ""
    repeat_copies: int = 4
    site_i_allowlist: tuple[str, ...] = ()
    site_ii_allowlist: tuple[str, ...] = ()
    site_iii_allowlist: tuple[str, ...] = ()
    plasmid_allowlist: tuple[str, ...] = ()
    allow_left_cutter_in_hurdler_pair: bool = False
    allow_right_cutter_in_hurdler_pair: bool = False
    pair_policy: str = EXACT_DNA_PAIR_POLICY
    max_purchase_bp: int = DEFAULT_GBLOCK_MAX_BP
    max_states: int = DEFAULT_MAX_STATES
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS
    paths_per_state: int = DEFAULT_PATHS_PER_STATE
    max_complete_routes: int = DEFAULT_MAX_COMPLETE_ROUTES

    def __post_init__(self) -> None:
        if self.schema_version != EXACT_DNA_SCHEMA_VERSION:
            raise ValueError(f"schema_version must be {EXACT_DNA_SCHEMA_VERSION!r}")
        if self.input_mode not in EXACT_DNA_INPUT_MODES:
            raise ValueError(f"input_mode must be one of {EXACT_DNA_INPUT_MODES}")
        if self.pair_policy != EXACT_DNA_PAIR_POLICY:
            raise ValueError(f"pair_policy must be {EXACT_DNA_PAIR_POLICY!r}")
        for name in (
            "site_i_allowlist",
            "site_ii_allowlist",
            "site_iii_allowlist",
            "plasmid_allowlist",
        ):
            object.__setattr__(self, name, _normalize_allowlist(getattr(self, name)))
        if not DEFAULT_OLIGO_MIN_BP <= int(self.max_purchase_bp) <= DEFAULT_GBLOCK_MAX_BP:
            raise ValueError("max_purchase_bp must be between 20 and 3000")
        if int(self.max_states) < 1 or int(self.timeout_seconds) < 1:
            raise ValueError("Search state and timeout limits must be positive")
        if int(self.paths_per_state) < 1 or int(self.max_complete_routes) < 1:
            raise ValueError("Route retention limits must be positive")
        if self.input_mode == "array":
            _id, unit = parse_exact_dna_input(self.repeat_unit, sequence_id=self.sequence_id)
            spacer = ""
            if self.spacer:
                _unused, spacer = parse_exact_dna_input(self.spacer, sequence_id="spacer")
            if int(self.repeat_copies) < 2:
                raise ValueError("repeat_copies must be at least two")
            object.__setattr__(self, "repeat_unit", unit)
            object.__setattr__(self, "spacer", spacer)
            object.__setattr__(self, "sequence_id", _id)
        else:
            resolved_id, sequence = parse_exact_dna_input(
                self.exact_dna, sequence_id=self.sequence_id
            )
            object.__setattr__(self, "exact_dna", sequence)
            object.__setattr__(self, "sequence_id", resolved_id)

    @property
    def array_unit(self) -> str:
        return self.repeat_unit + self.spacer

    def array_sequence(self, copies: int) -> str:
        if copies < 1:
            raise ValueError("copies must be positive")
        return self.spacer.join([self.repeat_unit] * copies)

    @property
    def target_sequence(self) -> str:
        if self.input_mode == "array":
            return self.array_sequence(int(self.repeat_copies))
        return self.exact_dna

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ExactDNAQuery":
        allowed = {item.name for item in fields(cls)}
        unknown = sorted(set(payload) - allowed)
        if unknown:
            raise ValueError("Unknown ExactDNAQuery fields: " + ", ".join(unknown))
        return cls(**dict(payload))


@dataclass(frozen=True)
class ExactDNASelection:
    route_id: str
    validation_mode: str = "none"
    plasmid_profile: str = ""
    cut_scheme_id: str = ""
    site_i_enzyme: str = ""
    site_ii_enzyme: str = ""

    def __post_init__(self) -> None:
        if self.validation_mode not in EXACT_DNA_VALIDATION_MODES:
            raise ValueError(
                f"validation_mode must be one of {EXACT_DNA_VALIDATION_MODES}"
            )


@dataclass
class ExactDNAResult:
    schema_version: str
    status: str
    message: str
    query: dict[str, Any]
    target_sequence: str
    target_sequence_sha256: str
    target_length_bp: int
    direct_purchase_product: str = ""
    restriction_hits: list[dict[str, Any]] = field(default_factory=list)
    pair_candidates: list[dict[str, Any]] = field(default_factory=list)
    route_candidates: list[dict[str, Any]] = field(default_factory=list)
    selected_route: dict[str, Any] | None = None
    seed: dict[str, Any] | None = None
    transitions: list[dict[str, Any]] = field(default_factory=list)
    purchase_fragments: list[dict[str, Any]] = field(default_factory=list)
    latent_transitions: list[dict[str, Any]] = field(default_factory=list)
    restoration_segments: list[dict[str, Any]] = field(default_factory=list)
    cloning_steps: list[dict[str, Any]] = field(default_factory=list)
    idt_audit: list[dict[str, Any]] = field(default_factory=list)
    route_attempts: list[dict[str, Any]] = field(default_factory=list)
    whole_target_idt: dict[str, Any] = field(default_factory=dict)
    final_insert_sequence: str = ""
    final_plasmid: dict[str, Any] | None = None
    independent_verification: dict[str, Any] = field(default_factory=dict)
    search_summary: dict[str, Any] = field(default_factory=dict)
    termination_reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _query_payload(query: ExactDNAQuery) -> dict[str, Any]:
    payload = asdict(query)
    for key in (
        "site_i_allowlist",
        "site_ii_allowlist",
        "site_iii_allowlist",
        "plasmid_allowlist",
    ):
        payload[key] = list(payload[key])
    return payload


def _seed_evidence(sequence: str) -> dict[str, Any] | None:
    length = len(sequence)
    if not DEFAULT_OLIGO_MIN_BP <= length <= DEFAULT_GBLOCK_MAX_BP:
        return None
    if length < PRIMER_PAIR_CORE_THRESHOLD_BP:
        return {
            "seed_sequence": sequence,
            "seed_length_bp": length,
            "product_type": "duplexed_seed_oligo_pair",
            "purchase_sequence": sequence,
            "secondary_purchase_sequence": reverse_complement(sequence),
            "purchase_sequence_count": 2,
            "purchase_length_bp": length * 2,
            "purchase_sha256": _sha(
                "duplexed_seed_oligo_pair|" + sequence + "|" + reverse_complement(sequence)
            ),
        }
    product = _product_type(length)
    if product is None:
        return None
    return {
        "seed_sequence": sequence,
        "seed_length_bp": length,
        "product_type": product,
        "purchase_sequence": sequence,
        "secondary_purchase_sequence": "",
        "purchase_sequence_count": 1,
        "purchase_length_bp": length,
        "purchase_sha256": _sha(product + "|" + sequence),
    }


def _primer_fragment_lengths_valid(row: Mapping[str, Any]) -> bool:
    if row.get("product_type") != "annealed_sticky_end_primer_pair":
        length = int(row.get("purchase_length_bp", 0))
        return DEFAULT_OLIGO_MIN_BP <= length <= DEFAULT_GBLOCK_MAX_BP
    forward = str(row.get("primer_forward_5to3", ""))
    reverse = str(row.get("primer_reverse_5to3", ""))
    return all(
        DEFAULT_OLIGO_MIN_BP <= len(sequence) <= DEFAULT_OLIGO_MAX_BP
        for sequence in (forward, reverse)
    )


def _pair_key(edge: Mapping[str, Any]) -> str:
    return f"{edge['site_i_enzyme']}|{edge['site_ii_enzyme']}"


def _edge_rank(edge: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        int(edge["step_count"]),
        int(edge["total_purchase_bp"]),
        -int(edge["replacement_length_bp"]),
        _pair_key(edge),
        str(edge["edge_id"]),
    )


def _enumerate_molecular_edges(
    sequence: str,
    query: ExactDNAQuery,
    geometries: dict[str, EnzymeGeometry],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], bool]:
    """Return locally validated reverse edges for one exact result state."""
    state_id = "state_" + _sha(sequence)[:16]
    allowed_i = set(query.site_i_allowlist)
    allowed_ii = set(query.site_ii_allowlist)
    allowed_iii = set(query.site_iii_allowlist)
    scan_geometries = {
        name: geometry
        for name, geometry in geometries.items()
        if (
            geometry.site_i_eligible and (not allowed_i or name in allowed_i)
        )
        or (
            geometry.site_ii_eligible and (not allowed_ii or name in allowed_ii)
        )
        or (
            geometry.site_iii_eligible and (not allowed_iii or name in allowed_iii)
        )
    }
    hits = scan_re_sites(state_id, sequence, scan_geometries)
    schemes = enumerate_active_latent_pairs(hits, None)
    schemes = [
        scheme
        for scheme in schemes
        if (
            not allowed_i or str(scheme["site_i"].enzyme) in allowed_i
        )
        and (
            not allowed_ii or str(scheme["site_ii"].enzyme) in allowed_ii
        )
    ]
    if not schemes:
        return [], [hit.to_dict() for hit in hits], False
    prepared = []
    for scheme in schemes:
        item = dict(scheme)
        # A single placeholder column lets the battle-tested molecular
        # simulator run without applying the historical plasmid mask.
        item["compatible_plasmids"] = ["pUC18"]
        prepared.append(item)
    route_limit = max(100, len(prepared) * 2)
    planned = [
        plan_target(
            TargetRecord(state_id, sequence),
            geometries,
            pd.DataFrame(),
            require_idt=False,
            top_routes=route_limit,
            max_purchase_bp=ceiling,
            _precomputed_hits=hits,
            _precomputed_schemes=prepared,
        )
        for ceiling in dict.fromkeys((int(query.max_purchase_bp), 60))
    ]
    routes = pd.concat(
        [item["routes"] for item in planned if not item["routes"].empty],
        ignore_index=True,
        sort=False,
    ) if any(not item["routes"].empty for item in planned) else pd.DataFrame()
    if routes.empty:
        return [], [hit.to_dict() for hit in hits], False
    all_steps = pd.concat(
        [item["steps"] for item in planned if not item["steps"].empty],
        ignore_index=True,
        sort=False,
    ) if any(not item["steps"].empty for item in planned) else pd.DataFrame()
    all_fragments = pd.concat(
        [item["fragments"] for item in planned if not item["fragments"].empty],
        ignore_index=True,
        sort=False,
    ) if any(not item["fragments"].empty for item in planned) else pd.DataFrame()
    steps_by_route = {
        route_id: rows.to_dict("records")
        for route_id, rows in all_steps.groupby("route_id", sort=False)
    } if not all_steps.empty else {}
    fragments_by_route = {
        route_id: rows.to_dict("records")
        for route_id, rows in all_fragments.groupby("route_id", sort=False)
    } if not all_fragments.empty else {}
    edges: list[dict[str, Any]] = []
    for row in routes.to_dict("records"):
        if not bool(row.get("local_constraints_passed")):
            continue
        fragments = fragments_by_route.get(str(row["route_id"]), [])
        if not fragments or not all(_primer_fragment_lengths_valid(item) for item in fragments):
            continue
        start, end = int(row["replacement_start"]), int(row["replacement_end"])
        precursor = sequence[:start] + sequence[end:]
        if _seed_evidence(precursor) is None and len(precursor) < DEFAULT_OLIGO_MIN_BP:
            continue
        edge_payload = {
            "result_sha256": _sha(sequence),
            "precursor_sha256": _sha(precursor),
            "route_id": row["route_id"],
            "pair": [row["site_i_enzyme"], row["site_ii_enzyme"]],
            "replacement": [start, end],
        }
        edge_id = "edge_" + _sha(
            json.dumps(edge_payload, sort_keys=True, separators=(",", ":"))
        )[:20]
        edge = {
            **row,
            "edge_id": edge_id,
            "result_sequence": sequence,
            "result_sequence_sha256": _sha(sequence),
            "precursor_sequence": precursor,
            "precursor_sequence_sha256": _sha(precursor),
            "donor_core_sequence": sequence[start:end],
            "donor_core_sha256": _sha(sequence[start:end]),
            "steps": steps_by_route.get(str(row["route_id"]), []),
            "fragments": fragments,
        }
        if edge["precursor_sequence"] + edge["donor_core_sequence"] == sequence:
            # This equality only holds for terminal insertions.  General
            # intervals are checked by coordinates below.
            pass
        rebuilt = precursor[:start] + sequence[start:end] + precursor[start:]
        if rebuilt != sequence:
            continue
        edges.append(edge)
    deduplicated: dict[tuple[Any, ...], dict[str, Any]] = {}
    for edge in sorted(edges, key=_edge_rank):
        key = (
            edge["precursor_sequence_sha256"],
            _pair_key(edge),
            int(edge["replacement_start"]),
            int(edge["replacement_end"]),
            str(edge["direction"]),
        )
        deduplicated.setdefault(key, edge)
    truncated = any(len(item["routes"]) >= route_limit for item in planned)
    return list(deduplicated.values()), [hit.to_dict() for hit in hits], truncated


def _path_rank(path: Mapping[str, Any]) -> tuple[Any, ...]:
    edges = path["edges"]
    return (
        sum(int(edge["step_count"]) for edge in edges),
        len({fragment["purchase_sha256"] for edge in edges for fragment in edge["fragments"]}),
        sum(int(edge["total_purchase_bp"]) for edge in edges),
        len({_pair_key(edge) for edge in edges}),
        str(path["seed"]["purchase_sha256"]),
    )


def _path_pair_signature(path: Mapping[str, Any]) -> tuple[str, ...]:
    return tuple(sorted({_pair_key(edge) for edge in path["edges"]}))


def _retain_paths_per_pair(
    paths: Sequence[Mapping[str, Any]],
    *,
    paths_per_pair: int,
) -> list[dict[str, Any]]:
    """Retain representative paths without letting one RE pair crowd out others."""
    grouped: dict[tuple[str, ...], list[dict[str, Any]]] = {}
    for path in sorted(paths, key=_path_rank):
        key = _path_pair_signature(path)
        values = grouped.setdefault(key, [])
        if len(values) < paths_per_pair:
            values.append(dict(path))
    return [
        path
        for key in sorted(grouped)
        for path in grouped[key]
    ]


def _array_paths(
    query: ExactDNAQuery,
    geometries: dict[str, EnzymeGeometry],
    *,
    fixed_pair: bool,
    deadline: float,
    progress_callback: ProgressCallback | None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    target_copies = int(query.repeat_copies)
    sequence_by_copies = {
        copies: query.array_sequence(copies)
        for copies in range(1, target_copies + 1)
    }
    edges_by_recipient: dict[int, list[tuple[int, dict[str, Any]]]] = {}
    states_explored = 0
    truncated = False
    target_hits: list[dict[str, Any]] = []
    for result_copies in range(2, target_copies + 1):
        if time.monotonic() >= deadline or states_explored >= int(query.max_states):
            return [], {
                "complete": False,
                "reason": "search_budget_exhausted",
                "states_explored": states_explored,
                "target_hits": target_hits,
            }
        sequence = sequence_by_copies[result_copies]
        edges, hits, was_truncated = _enumerate_molecular_edges(sequence, query, geometries)
        states_explored += 1
        truncated = truncated or was_truncated
        if result_copies == target_copies:
            target_hits = hits
        for edge in edges:
            precursor = str(edge["precursor_sequence"])
            recipient_copies = next(
                (
                    copies
                    for copies, candidate in sequence_by_copies.items()
                    if copies < result_copies and precursor == candidate
                ),
                0,
            )
            if recipient_copies < 1:
                continue
            edge["recipient_copy_count"] = recipient_copies
            edge["result_copy_count"] = result_copies
            edge["donor_copy_count"] = result_copies - recipient_copies
            edges_by_recipient.setdefault(recipient_copies, []).append((result_copies, edge))
        emit_progress(
            progress_callback,
            stage="molecular_search",
            status="state_completed",
            copies=result_copies,
            message=f"Scanned exact {result_copies}-copy array state",
            details={"edge_count": len(edges)},
        )
    state_paths: dict[tuple[int, str], list[dict[str, Any]]] = {}
    # Match complete-route-v2: establish the shortest exact purchasable seed,
    # then prove every subsequent nucleotide through HURDLER growth.
    for copies in range(1, target_copies):
        seed = _seed_evidence(sequence_by_copies[copies])
        if seed is not None:
            state_paths[(copies, "")] = [{"seed": seed, "edges": []}]
            break
    for recipient in range(1, target_copies):
        source_items = [item for key, values in state_paths.items() if key[0] == recipient for item in values]
        for path in source_items:
            for result_copies, edge in edges_by_recipient.get(recipient, []):
                pair = _pair_key(edge)
                existing_pairs = {_pair_key(item) for item in path["edges"]}
                if fixed_pair and existing_pairs and pair not in existing_pairs:
                    continue
                new_path = {"seed": path["seed"], "edges": [*path["edges"], edge]}
                key = (result_copies, pair if fixed_pair else "variable")
                values = [*state_paths.get(key, []), new_path]
                state_paths[key] = sorted(values, key=_path_rank)[: int(query.paths_per_state)]
    complete = [
        path for (copies, _pair), values in state_paths.items()
        if copies == target_copies for path in values if path["edges"]
    ]
    retained_complete = _retain_paths_per_pair(
        complete,
        paths_per_pair=int(query.paths_per_state),
    )
    return retained_complete, {
        "complete": not truncated,
        "reason": "route_found" if complete else "exhausted",
        "states_explored": states_explored,
        "target_hits": target_hits,
        "edge_count": sum(len(values) for values in edges_by_recipient.values()),
    }


def _arbitrary_paths(
    query: ExactDNAQuery,
    geometries: dict[str, EnzymeGeometry],
    *,
    fixed_pair: bool,
    deadline: float,
    progress_callback: ProgressCallback | None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    target = query.target_sequence
    serial = 0
    heap: list[tuple[int, int, int, str, list[dict[str, Any]], str]] = [
        (len(target), 0, serial, target, [], "")
    ]
    retained: dict[tuple[str, str], int] = {}
    edge_cache: dict[str, tuple[list[dict[str, Any]], list[dict[str, Any]], bool]] = {}
    complete: list[dict[str, Any]] = []
    states_explored = 0
    target_hits: list[dict[str, Any]] = []
    truncated = False
    complete_per_pair: dict[tuple[str, ...], int] = {}
    while heap:
        if time.monotonic() >= deadline or states_explored >= int(query.max_states):
            return sorted(complete, key=_path_rank), {
                "complete": False,
                "reason": "search_budget_exhausted",
                "states_explored": states_explored,
                "target_hits": target_hits,
            }
        _length, _steps, _serial, sequence, reverse_edges, required_pair = heapq.heappop(heap)
        key = (_sha(sequence), required_pair if fixed_pair else "variable")
        count = retained.get(key, 0)
        if count >= int(query.paths_per_state):
            continue
        retained[key] = count + 1
        states_explored += 1
        if reverse_edges:
            seed = _seed_evidence(sequence)
            if seed is not None:
                candidate = {"seed": seed, "edges": list(reversed(reverse_edges))}
                signature = _path_pair_signature(candidate)
                if complete_per_pair.get(signature, 0) < int(query.paths_per_state):
                    complete.append(candidate)
                    complete_per_pair[signature] = complete_per_pair.get(signature, 0) + 1
                continue
        state_sha = _sha(sequence)
        if state_sha not in edge_cache:
            edge_cache[state_sha] = _enumerate_molecular_edges(sequence, query, geometries)
        edges, hits, was_truncated = edge_cache[state_sha]
        truncated = truncated or was_truncated
        if sequence == target:
            target_hits = hits
        for edge in edges:
            pair = _pair_key(edge)
            if fixed_pair and required_pair and pair != required_pair:
                continue
            precursor = str(edge["precursor_sequence"])
            serial += 1
            next_pair = pair if fixed_pair else "variable"
            heapq.heappush(
                heap,
                (
                    len(precursor),
                    sum(int(item["step_count"]) for item in reverse_edges) + int(edge["step_count"]),
                    serial,
                    precursor,
                    [*reverse_edges, edge],
                    next_pair,
                ),
            )
        emit_progress(
            progress_callback,
            stage="molecular_search",
            status="state_completed",
            message=f"Explored exact precursor state {states_explored}",
            details={"state_length_bp": len(sequence), "edge_count": len(edges)},
        )
    return _retain_paths_per_pair(
        complete,
        paths_per_pair=int(query.paths_per_state),
    ), {
        "complete": not truncated and not heap,
        "reason": "route_found" if complete else "exhausted",
        "states_explored": states_explored,
        "target_hits": target_hits,
    }


def _feature_preservation(
    database: PlasmidReferenceDatabase,
    scheme: VectorCutScheme,
    final_sequence: str,
) -> tuple[bool, list[str]]:
    profile = database.profile(scheme.profile_id)
    reference = database.reference(profile.reference_id)
    missing: list[str] = []
    for feature in reference.features:
        if not feature.protected:
            continue
        if any(
            start < profile.mcs_end and end > profile.mcs_start
            for start, end in feature.intervals
        ):
            continue
        segments = [reference.sequence[start:end] for start, end in feature.intervals]
        if any(
            segment
            and segment not in (final_sequence + final_sequence[: max(0, len(segment) - 1)])
            and reverse_complement(segment)
            not in (final_sequence + final_sequence[: max(0, len(segment) - 1)])
            for segment in segments
        ):
            missing.append(f"{feature.feature_class}:{feature.label}")
    return not missing, missing


def _motif_positions(sequence: str, motif: str) -> set[int]:
    positions: set[int] = set()
    for oriented in {motif.upper(), reverse_complement(motif.upper())}:
        start = sequence.upper().find(oriented)
        while start >= 0:
            positions.add(start)
            start = sequence.upper().find(oriented, start + 1)
    return positions


def _selected_pair_excess_counts(
    final_sequence: str,
    target_sequence: str,
    pair_rows: Mapping[str, Mapping[str, str]],
) -> dict[str, int]:
    excess: dict[str, int] = {}
    for pair in pair_rows.values():
        for enzyme_key, site_key in (
            ("site_i_enzyme", "site_i_recognition_site"),
            ("site_ii_enzyme", "site_ii_recognition_site"),
        ):
            enzyme = pair[enzyme_key]
            motif = pair[site_key]
            excess[enzyme] = max(
                excess.get(enzyme, 0),
                len(_motif_positions(final_sequence, motif))
                - len(_motif_positions(target_sequence, motif)),
            )
    return excess


def _annotation_routes_for_path(
    path: Mapping[str, Any],
    query: ExactDNAQuery,
    database: PlasmidReferenceDatabase,
) -> list[dict[str, Any]]:
    pair_rows: dict[str, dict[str, str]] = {}
    for edge in path["edges"]:
        pair_rows.setdefault(
            _pair_key(edge),
            {
                "site_i_enzyme": str(edge["site_i_enzyme"]),
                "site_ii_enzyme": str(edge["site_ii_enzyme"]),
                "site_i_recognition_site": _GEOMETRY_CONTEXT[
                    str(edge["site_i_enzyme"])
                ].recognition_site,
                "site_ii_recognition_site": _GEOMETRY_CONTEXT[
                    str(edge["site_ii_enzyme"])
                ].recognition_site,
            },
        )
    routes: list[dict[str, Any]] = []
    for scheme in database.schemes:
        if not scheme.valid or scheme.left_cutter is None or scheme.right_cutter is None:
            continue
        profile = database.profile(scheme.profile_id)
        if query.plasmid_allowlist and profile.profile_id not in set(query.plasmid_allowlist):
            continue
        if any(
            retained_backbone_contains_site(scheme, pair[site_key])
            for pair in pair_rows.values()
            for site_key in ("site_i_recognition_site", "site_ii_recognition_site")
        ):
            continue
        used_enzymes = {
            pair[key] for pair in pair_rows.values()
            for key in ("site_i_enzyme", "site_ii_enzyme")
        }
        left_reused = bool(used_enzymes & set(scheme.left_cutter.enzyme_aliases))
        right_reused = bool(used_enzymes & set(scheme.right_cutter.enzyme_aliases))
        if left_reused and not query.allow_left_cutter_in_hurdler_pair:
            continue
        if right_reused and not query.allow_right_cutter_in_hurdler_pair:
            continue
        final = (
            scheme.retained_backbone_sequence
            + scheme.left_restoration_sequence
            + query.target_sequence
            + scheme.right_restoration_sequence
        )
        selected_pair_excess = _selected_pair_excess_counts(
            final, query.target_sequence, pair_rows
        )
        # Exact-DNA design never silently edits vector restoration sequence.
        # A reused cutter is accepted only if ligation naturally destroys its
        # additional site; otherwise this is not yet a complete exact route.
        if any(value > 0 for value in selected_pair_excess.values()):
            continue
        decisions = []
        if left_reused:
            decisions.append(
                decide_cutter_silencing(
                    database, profile.profile_id, scheme.left_cutter,
                    junction_destroyed=True,
                )
            )
        if right_reused:
            decisions.append(
                decide_cutter_silencing(
                    database, profile.profile_id, scheme.right_cutter,
                    junction_destroyed=True,
                )
            )
        if any(not decision.allowed for decision in decisions):
            continue
        preserved, missing = _feature_preservation(database, scheme, final)
        if not preserved:
            continue
        route_payload = {
            "target_sha256": _sha(query.target_sequence),
            "edge_ids": [edge["edge_id"] for edge in path["edges"]],
            "profile_id": profile.profile_id,
            "scheme_id": scheme.scheme_id,
        }
        routes.append(
            {
                "route_id": "exact_route_" + _sha(
                    json.dumps(route_payload, sort_keys=True, separators=(",", ":"))
                )[:20],
                "profile_id": profile.profile_id,
                "reference_id": profile.reference_id,
                "scheme_id": scheme.scheme_id,
                "cut_scheme": f"{scheme.left_location}/{scheme.right_location}",
                "left_cutter": scheme.left_cutter.canonical_enzyme,
                "right_cutter": scheme.right_cutter.canonical_enzyme,
                "left_cutter_reused": left_reused,
                "right_cutter_reused": right_reused,
                "cutter_reuse": left_reused or right_reused,
                "silencing_decisions": [asdict(item) for item in decisions],
                "left_restoration_sequence": scheme.left_restoration_sequence,
                "right_restoration_sequence": scheme.right_restoration_sequence,
                "restoration_length_bp": len(scheme.left_restoration_sequence)
                + len(scheme.right_restoration_sequence),
                "retained_backbone_sha256": scheme.retained_backbone_sha256,
                "pair_mode": "fixed" if len(pair_rows) == 1 else "pair_changes",
                "pairs": list(pair_rows.values()),
                "pair_change_count": max(0, len(pair_rows) - 1),
                "hurdler_step_count": sum(int(edge["step_count"]) for edge in path["edges"]),
                "transition_count": len(path["edges"]),
                "unique_purchase_count": len(
                    {fragment["purchase_sha256"] for edge in path["edges"] for fragment in edge["fragments"]}
                    | {str(path["seed"]["purchase_sha256"])}
                ),
                "total_purchase_bp": int(path["seed"]["purchase_length_bp"])
                + sum(int(edge["total_purchase_bp"]) for edge in path["edges"]),
                "seed": path["seed"],
                "edges": path["edges"],
                "protected_features_preserved": preserved,
                "missing_protected_features": missing,
                "selected_pair_excess_sites": selected_pair_excess,
                "plasmid_reference_version": PLASMID_REFERENCE_VERSION,
            }
        )
    if any(not route["cutter_reuse"] for route in routes):
        routes = [route for route in routes if not route["cutter_reuse"]]
    return routes


# The geometry context is scoped to a synchronous query and avoids serializing
# Biopython enzyme objects into every path record.
_GEOMETRY_CONTEXT: dict[str, EnzymeGeometry] = {}


def _public_route(route: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in route.items() if key != "edges"}


def _retain_annotated_route_groups(
    routes: Sequence[Mapping[str, Any]],
    *,
    routes_per_group: int,
) -> list[dict[str, Any]]:
    """Keep each discovered RE-pair/plasmid/cut-scheme group visible."""
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for route in routes:
        pair_signature = tuple(
            sorted(
                (str(pair["site_i_enzyme"]), str(pair["site_ii_enzyme"]))
                for pair in route["pairs"]
            )
        )
        key = (
            pair_signature,
            str(route["profile_id"]),
            str(route["scheme_id"]),
            str(route["pair_mode"]),
        )
        values = grouped.setdefault(key, [])
        if len(values) < routes_per_group:
            values.append(dict(route))
    return [route for key in sorted(grouped, key=str) for route in grouped[key]]


def query_exact_dna(
    query: ExactDNAQuery,
    *,
    reference_dir: str | Path | None = None,
    artifact_dir: str | Path | None = None,
    plasmid_database: PlasmidReferenceDatabase | None = None,
    plasmid_reference_path: str | Path | None = None,
    progress_callback: ProgressCallback | None = None,
) -> ExactDNAResult:
    """Find complete molecular routes without making any IDT HTTP request."""
    paths = ProjectPaths.discover()
    geometries = load_exact_dna_enzyme_catalog(
        reference_dir=reference_dir or paths.reference_output,
        # Keep the broad conventional-RE pool; the exact loader adds the
        # maintained Site-III adapters from ``output/`` when needed.
        artifact_dir=artifact_dir or paths.root / "data" / "artifacts",
    )
    global _GEOMETRY_CONTEXT
    _GEOMETRY_CONTEXT = geometries
    unknown_i = sorted(set(query.site_i_allowlist) - set(geometries))
    unknown_ii = sorted(set(query.site_ii_allowlist) - set(geometries))
    unknown_iii = sorted(set(query.site_iii_allowlist) - set(geometries))
    invalid_iii = sorted(
        name for name in query.site_iii_allowlist
        if name in geometries and not geometries[name].site_iii_eligible
    )
    if unknown_i or unknown_ii or unknown_iii:
        raise ValueError(
            "Unknown restriction enzymes: "
            + ", ".join(unknown_i + unknown_ii + unknown_iii)
        )
    if invalid_iii:
        raise ValueError("Enzymes are not eligible for Site III: " + ", ".join(invalid_iii))
    unknown_plasmids = sorted(set(query.plasmid_allowlist) - set(PLASMIDS))
    if unknown_plasmids:
        raise ValueError("Unknown plasmid profiles: " + ", ".join(unknown_plasmids))
    database = plasmid_database or load_plasmid_reference(plasmid_reference_path)
    target = query.target_sequence
    started = time.monotonic()
    deadline = started + int(query.timeout_seconds)
    emit_progress(
        progress_callback,
        stage="molecular_search",
        status="started",
        message="Searching exact seed-to-target HURDLER routes",
        details={"target_length_bp": len(target), "input_mode": query.input_mode},
    )
    search = _array_paths if query.input_mode == "array" else _arbitrary_paths
    fixed_paths, fixed_summary = search(
        query,
        geometries,
        fixed_pair=True,
        deadline=deadline,
        progress_callback=progress_callback,
    )
    search_summaries = {"fixed_pair": fixed_summary}
    annotated: list[dict[str, Any]] = []
    for path in fixed_paths:
        annotated.extend(_annotation_routes_for_path(path, query, database))
    pair_mode = "fixed"
    if not annotated and time.monotonic() < deadline:
        variable_paths, variable_summary = search(
            query,
            geometries,
            fixed_pair=False,
            deadline=deadline,
            progress_callback=progress_callback,
        )
        pair_mode = "pair_changes"
        search_summaries["pair_changes"] = variable_summary
        for path in variable_paths:
            annotated.extend(_annotation_routes_for_path(path, query, database))
    annotated.sort(
        key=lambda route: (
            route["pair_mode"] != "fixed",
            bool(route["cutter_reuse"]),
            int(route["hurdler_step_count"]),
            int(route["unique_purchase_count"]),
            int(route["total_purchase_bp"]),
            int(route["restoration_length_bp"]),
            str(route["profile_id"]),
            str(route["scheme_id"]),
            str(route["route_id"]),
        )
    )
    annotated = _retain_annotated_route_groups(
        annotated,
        routes_per_group=int(query.max_complete_routes),
    )
    target_hits = fixed_summary.get("target_hits", [])
    pair_candidates: dict[tuple[str, str], dict[str, Any]] = {}
    for route in annotated:
        for pair in route["pairs"]:
            key = (pair["site_i_enzyme"], pair["site_ii_enzyme"])
            row = pair_candidates.setdefault(
                key,
                {
                    **pair,
                    "supported_plasmids": set(),
                    "route_count": 0,
                    "pair_mode": pair_mode,
                },
            )
            row["supported_plasmids"].add(route["profile_id"])
            row["route_count"] += 1
    public_pairs = []
    for row in pair_candidates.values():
        row = dict(row)
        row["supported_plasmids"] = sorted(row["supported_plasmids"])
        public_pairs.append(row)
    exhaustive = all(bool(item.get("complete")) for item in search_summaries.values())
    if annotated:
        status = "hurdler_compatible_molecular"
        message = (
            "At least one exact seed-to-target route is molecularly verified; "
            "IDT has not been called."
        )
        termination = "complete_route_found"
    elif not exhaustive:
        status = "search_incomplete"
        message = "The route search reached its state or time budget; compatibility is unclassified."
        termination = "search_budget_exhausted"
    elif _seed_evidence(target) is not None:
        status = "direct_purchase_only"
        message = "The exact target fits a direct synthesis product, but no complete HURDLER growth route was found."
        termination = "no_hurdler_transition"
    else:
        status = "no_complete_hurdler_route"
        message = "The complete exact route search finished without an annotation-safe seed-to-target route."
        termination = "search_exhausted"
    result = ExactDNAResult(
        schema_version=EXACT_DNA_SCHEMA_VERSION,
        status=status,
        message=message,
        query=_query_payload(query),
        target_sequence=target,
        target_sequence_sha256=_sha(target),
        target_length_bp=len(target),
        direct_purchase_product=(
            str(_seed_evidence(target)["product_type"])
            if _seed_evidence(target) is not None else ""
        ),
        restriction_hits=target_hits,
        pair_candidates=sorted(
            public_pairs,
            key=lambda row: (row["site_i_enzyme"], row["site_ii_enzyme"]),
        ),
        route_candidates=[_public_route(route) for route in annotated],
        search_summary={
            **search_summaries,
            "pair_search_used": pair_mode,
            "elapsed_seconds": time.monotonic() - started,
            "route_count": len(annotated),
            "complete": exhaustive,
        },
        termination_reason=termination,
    )
    # Private route bodies are retained only in memory for confirmation.  They
    # are intentionally absent from public tables to keep Colab output small.
    setattr(result, "_route_bodies", {route["route_id"]: route for route in annotated})
    emit_progress(
        progress_callback,
        stage="molecular_search",
        status="completed" if status != "search_incomplete" else "incomplete",
        message=message,
        elapsed_seconds=time.monotonic() - started,
        details={"route_count": len(annotated), "status": status},
    )
    return result


def _circular_slice(sequence: str, start: int, end: int) -> str:
    length = len(sequence)
    start %= length
    end %= length
    if start == end:
        return sequence
    return sequence[start:end] if start < end else sequence[start:] + sequence[:end]


def _primary_purchase_row(
    route: Mapping[str, Any],
    query: Mapping[str, Any],
    database: PlasmidReferenceDatabase,
    geometries: Mapping[str, EnzymeGeometry],
) -> dict[str, Any]:
    """Design the actual vector-ready primary purchase molecule."""
    scheme = next(item for item in database.schemes if item.scheme_id == route["scheme_id"])
    if scheme.left_cutter is None or scheme.right_cutter is None:
        raise ValueError("Selected cut scheme has no vector cutters")
    profile = database.profile(str(route["profile_id"]))
    reference = database.reference(profile.reference_id)
    oriented = (
        reference.sequence
        if profile.expression_strand == 1
        else reverse_complement(reference.sequence)
    )
    left_overhang = _circular_slice(
        oriented,
        min(scheme.left_cutter.top_cut_oriented, scheme.left_cutter.bottom_cut_oriented),
        max(scheme.left_cutter.top_cut_oriented, scheme.left_cutter.bottom_cut_oriented),
    )
    right_overhang = _circular_slice(
        oriented,
        min(scheme.right_cutter.top_cut_oriented, scheme.right_cutter.bottom_cut_oriented),
        max(scheme.right_cutter.top_cut_oriented, scheme.right_cutter.bottom_cut_oriented),
    )
    core = (
        scheme.left_restoration_sequence
        + str(route["seed"]["seed_sequence"])
        + scheme.right_restoration_sequence
    )
    if len(core) < PRIMER_PAIR_CORE_THRESHOLD_BP:
        forward = left_overhang + core
        reverse = right_overhang + reverse_complement(core)
        product_type = "annealed_sticky_end_primer_pair"
        purchase = forward
        secondary = reverse
        left_adapter = right_adapter = ""
        left_adapter_enzyme = right_adapter_enzyme = "not_required"
        fragment_hash_input = f"{product_type}|{forward}|{reverse}"
    else:
        selected_geometries = {
            name: geometry
            for name, geometry in geometries.items()
            if not query.get("site_iii_allowlist")
            or name in set(query.get("site_iii_allowlist", []))
            or not geometry.site_iii_eligible
        }
        left_geometry = _select_type_iis(
            int(scheme.left_cutter.overhang_length), core, selected_geometries
        )
        right_geometry = _select_type_iis(
            int(scheme.right_cutter.overhang_length), core, selected_geometries
        )
        if left_geometry is None or right_geometry is None:
            raise ValueError("No selected Site-III enzyme can prepare the primary insert ends")
        left_adapter = _type_iis_flank(
            left_geometry, left=True, overhang_sequence=left_overhang
        )
        right_adapter = _type_iis_flank(
            right_geometry, left=False, overhang_sequence=right_overhang
        )
        purchase = left_adapter + core + right_adapter
        if len(purchase) > int(query["max_purchase_bp"]):
            raise ValueError("The vector-ready primary fragment exceeds max_purchase_bp")
        product_type = _product_type(len(purchase))
        if product_type is None:
            raise ValueError("The vector-ready primary fragment has no supported purchase format")
        secondary = ""
        forward = reverse = ""
        left_adapter_enzyme = left_geometry.enzyme
        right_adapter_enzyme = right_geometry.enzyme
        fragment_hash_input = f"{product_type}|{purchase}"
    digest = left_overhang + core + right_overhang
    purchase_sha = _sha(fragment_hash_input)
    return {
        "fragment_id": "primary_" + purchase_sha[:16],
        "stage": "primary_seed",
        "transition_index": 0,
        "target_start": 0,
        "target_end": len(core),
        "core_sequence": core,
        "core_length_bp": len(core),
        "purchase_sequence": purchase,
        "secondary_purchase_sequence": secondary,
        "primer_forward_5to3": forward,
        "primer_reverse_5to3": reverse,
        "purchase_sequence_count": 2 if secondary else 1,
        "purchase_length_bp": len(purchase) + len(secondary),
        "purchase_sha256": purchase_sha,
        "product_type": product_type,
        "left_adapter_enzyme": left_adapter_enzyme,
        "right_adapter_enzyme": right_adapter_enzyme,
        "left_adapter_sequence": left_adapter,
        "right_adapter_sequence": right_adapter,
        "digest_fragment_sequence": digest,
        "left_digest_overhang": left_overhang,
        "right_digest_overhang": right_overhang,
        "vector_left_cutter": scheme.left_cutter.canonical_enzyme,
        "vector_right_cutter": scheme.right_cutter.canonical_enzyme,
        "restoration_length_bp": len(scheme.left_restoration_sequence)
        + len(scheme.right_restoration_sequence),
        "idt_policy": IDT_SCORE_POLICY,
        "idt_status": "not_applicable_primer_pair_under_90bp" if secondary else "not_run",
        "idt_score": None,
        "idt_response_sha256": "",
    }


def _purchase_rows(
    route: Mapping[str, Any],
    query: Mapping[str, Any],
    database: PlasmidReferenceDatabase,
    geometries: Mapping[str, EnzymeGeometry],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    rows.append(_primary_purchase_row(route, query, database, geometries))
    for transition_index, edge in enumerate(route["edges"], start=1):
        for fragment in edge["fragments"]:
            row = dict(fragment)
            row["stage"] = "hurdler_donor"
            row["transition_index"] = transition_index
            rows.append(row)
    return rows


def _score_one_purchase(
    scorer: ComplexityScorer,
    name: str,
    sequence: str,
) -> dict[str, Any]:
    try:
        result = scorer.score(name, sequence)
    except Exception as exc:
        raise RuntimeError("idt_api_error") from exc
    score = result.get("idt_complexity_score")
    score_complete = result.get("idt_score_complete")
    if score_complete is not True or not isinstance(score, (int, float)) or isinstance(score, bool) or not math.isfinite(float(score)):
        invalid = str(result.get("idt_invalid_score_names_json", "[]"))
        raise ValueError(f"idt_score_error: invalid or incomplete rule scores: {invalid}")
    expected_pass = float(score) < 10.0
    explicit_pass = result.get("idt_explicit_pass")
    if not isinstance(explicit_pass, bool) or explicit_pass != expected_pass:
        raise ValueError("idt_score_error: IDT pass flag disagrees with strict score-sum policy")
    return result


def _simulate_final_plasmid(
    route: Mapping[str, Any],
    target: str,
    database: PlasmidReferenceDatabase,
) -> dict[str, Any]:
    scheme = next(item for item in database.schemes if item.scheme_id == route["scheme_id"])
    final = (
        scheme.retained_backbone_sequence
        + scheme.left_restoration_sequence
        + target
        + scheme.right_restoration_sequence
    )
    insert_start = len(scheme.retained_backbone_sequence) + len(scheme.left_restoration_sequence)
    insert = final[insert_start : insert_start + len(target)]
    if insert != target:
        raise AssertionError("Final vector simulation changed the exact target DNA")
    pair_rows = {
        f"{item['site_i_enzyme']}|{item['site_ii_enzyme']}": item
        for item in route["pairs"]
    }
    selected_pair_excess = _selected_pair_excess_counts(final, target, pair_rows)
    if any(value > 0 for value in selected_pair_excess.values()):
        raise AssertionError("Final vector contains an excess selected-pair RE site")
    preserved, missing = _feature_preservation(database, scheme, final)
    if not preserved:
        raise AssertionError("Protected plasmid features were lost: " + ", ".join(missing))
    return {
        "profile_id": route["profile_id"],
        "reference_id": route["reference_id"],
        "scheme_id": route["scheme_id"],
        "final_plasmid_sequence": final,
        "final_plasmid_length_bp": len(final),
        "final_plasmid_sha256": _sha(final),
        "circular": True,
        "insert_start_0based": insert_start,
        "insert_end_0based_exclusive": insert_start + len(target),
        "insert_sequence_sha256": _sha(insert),
        "insert_exact": True,
        "retained_backbone_sha256": scheme.retained_backbone_sha256,
        "restoration_exact": True,
        "protected_feature_sequences_preserved": True,
        "selected_pair_excess_sites": selected_pair_excess,
    }


def _latent_transition_rows(
    edges: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Normalize the transient one-base activations used by a complete route."""
    rows: list[dict[str, Any]] = []
    for transition_index, edge in enumerate(edges, start=1):
        donor_positions = {
            int(value)
            for value in json.loads(
                str(edge.get("donor_derived_mutation_positions_json", "[]"))
            )
        }
        orientation = {
            str(edge["site_i_enzyme"]): str(edge["site_i_orientation"]),
            str(edge["site_ii_enzyme"]): str(edge["site_ii_orientation"]),
        }
        role = {
            str(edge["site_i_enzyme"]): "Site I",
            str(edge["site_ii_enzyme"]): "Site II",
        }
        for step in edge["steps"]:
            states = json.loads(str(step.get("latent_state_transition_json", "{}")))
            for enzyme, state in states.items():
                mismatch = state.get("mismatch_position")
                donor_derived = mismatch is not None and int(mismatch) in donor_positions
                final_state = str(state.get("final", ""))
                rows.append(
                    {
                        "transition_index": transition_index,
                        "edge_id": edge["edge_id"],
                        "step_number": step.get("step_number"),
                        "enzyme_role": role.get(str(enzyme), ""),
                        "enzyme": enzyme,
                        "orientation": orientation.get(str(enzyme), ""),
                        "mechanism": edge["mechanism"],
                        "mismatch_position_0based": mismatch,
                        "mismatch_position_1based": (
                            int(mismatch) + 1 if mismatch is not None else None
                        ),
                        "target_base": state.get("target_base"),
                        "temporary_active_base": state.get("temporary_active_base"),
                        "intermediate_state": state.get("intermediate"),
                        "final_state": final_state,
                        "activation_base_source": (
                            "donor" if donor_derived else "retained_recipient"
                        ),
                        "donor_derived": donor_derived,
                        "final_target_restored": (
                            final_state == "latent" and bool(edge["final_sequence_exact"])
                        ),
                        "final_target_exact": bool(edge["final_sequence_exact"]),
                    }
                )
    return rows


def confirm_exact_dna_route(
    query_result: ExactDNAResult,
    selection: ExactDNASelection,
    *,
    idt_scorer: ComplexityScorer | None = None,
    plasmid_database: PlasmidReferenceDatabase | None = None,
    plasmid_reference_path: str | Path | None = None,
    progress_callback: ProgressCallback | None = None,
) -> ExactDNAResult:
    """Confirm one route, optionally score its exact purchase sequences."""
    if query_result.status != "hurdler_compatible_molecular":
        raise ValueError("Only a molecularly compatible query can confirm a route")
    bodies = getattr(query_result, "_route_bodies", {})
    if selection.route_id not in bodies:
        raise ValueError("Selected route is absent or belongs to a stale query")
    route = bodies[selection.route_id]
    expected_values = {
        "plasmid_profile": (selection.plasmid_profile, route["profile_id"]),
        "cut_scheme_id": (selection.cut_scheme_id, route["scheme_id"]),
    }
    selected_pairs = {
        (item["site_i_enzyme"], item["site_ii_enzyme"])
        for item in route["pairs"]
    }
    if selection.site_i_enzyme or selection.site_ii_enzyme:
        requested_pair = (selection.site_i_enzyme, selection.site_ii_enzyme)
        if not all(requested_pair) or requested_pair not in selected_pairs:
            raise ValueError("Selected Site-I/Site-II enzymes do not match the route")
    for label, (requested, observed) in expected_values.items():
        if requested and requested != observed:
            raise ValueError(f"Selected {label} does not match the route")
    result = ExactDNAResult(**{
        key: value for key, value in asdict(query_result).items()
        if key in {item.name for item in fields(ExactDNAResult)}
    })
    result.selected_route = _public_route(route)
    result.seed = dict(route["seed"])
    result.transitions = [
        {key: value for key, value in edge.items() if key not in {"steps", "fragments", "result_sequence", "precursor_sequence"}}
        for edge in route["edges"]
    ]
    result.latent_transitions = _latent_transition_rows(route["edges"])
    result.restoration_segments = [
        {
            "side": "left",
            "sequence": route["left_restoration_sequence"],
            "length_bp": len(route["left_restoration_sequence"]),
        },
        {
            "side": "right",
            "sequence": route["right_restoration_sequence"],
            "length_bp": len(route["right_restoration_sequence"]),
        },
    ]
    database = plasmid_database or load_plasmid_reference(plasmid_reference_path)
    geometries = load_exact_dna_enzyme_catalog()
    try:
        result.purchase_fragments = _purchase_rows(
            route, query_result.query, database, geometries
        )
        verification = verify_exact_dna_assembly(
            query=query_result.query,
            target_sequence=query_result.target_sequence,
            route=route,
            primary_fragment=result.purchase_fragments[0],
            database=database,
            geometries=geometries,
        )
    except Exception as exc:
        result.status = "independent_assembly_validation_failed"
        result.message = f"Independent assembly validation failed: {exc}"
        result.termination_reason = "independent_verifier_exception"
        result.independent_verification = {
            "passed": False,
            "verifier": "independent-circular-digest-ligation-v1",
            "errors": [str(exc)],
        }
        return result
    result.independent_verification = {
        "passed": bool(verification["passed"]),
        "verifier": verification["verifier"],
        "errors": list(verification["errors"]),
        "cloning_step_count": len(verification["steps"]),
        "protected_features_preserved": bool(
            verification["protected_features_preserved"]
        ),
        "selected_pair_excess_sites": dict(
            verification["selected_pair_excess_sites"]
        ),
    }
    if not verification["passed"]:
        result.status = "independent_assembly_validation_failed"
        result.message = (
            "The route planner found a candidate, but the independent digest/ligation "
            "verifier rejected it: " + "; ".join(verification["errors"])
        )
        result.termination_reason = "independent_verifier_disagreement"
        return result
    setattr(result, "_assembly_bundle", verification)
    result.cloning_steps = [
        {
            "step": int(step["step"]),
            "stage": str(step["stage"]),
            "input_plasmid": str(step["input_plasmid"]),
            "insert": f"step{int(step['step']):02d}_insert.gb",
            "purchase_fragment_ids": ";".join(step["purchase_fragment_ids"]),
            "restriction_enzymes": "/".join(step["enzymes"]),
            "site_iii_enzymes": "/".join(step.get("site_iii_enzymes", [])),
            "insert_core_length_bp": len(str(step["insert_core_sequence"])),
            "result_insert_length_bp": len(str(step["result_insert_sequence"])),
            "output_plasmid": str(step["output_plasmid"]),
            "insert_reused": False,
        }
        for step in verification["steps"]
    ]
    fragment_use_counts: dict[str, int] = {}
    for step in result.cloning_steps:
        fragment_id = str(step["purchase_fragment_ids"])
        fragment_use_counts[fragment_id] = fragment_use_counts.get(fragment_id, 0) + 1
    seen_fragments: set[str] = set()
    for step in result.cloning_steps:
        fragment_id = str(step["purchase_fragment_ids"])
        step["insert_reused"] = (
            fragment_use_counts[fragment_id] > 1 and fragment_id in seen_fragments
        )
        seen_fragments.add(fragment_id)
    result.final_insert_sequence = str(verification["final_insert_sequence"])
    insert_start = len(str(verification["retained_backbone_sequence"])) + len(
        route["left_restoration_sequence"]
    )
    result.final_plasmid = {
        "profile_id": route["profile_id"],
        "reference_id": route["reference_id"],
        "scheme_id": route["scheme_id"],
        "final_plasmid_sequence": verification["final_plasmid_sequence"],
        "final_plasmid_length_bp": len(verification["final_plasmid_sequence"]),
        "final_plasmid_sha256": _sha(verification["final_plasmid_sequence"]),
        "circular": True,
        "insert_start_0based": insert_start,
        "insert_end_0based_exclusive": insert_start + len(query_result.target_sequence),
        "insert_sequence_sha256": _sha(query_result.target_sequence),
        "insert_exact": True,
        "restoration_exact": True,
        "protected_feature_sequences_preserved": True,
        "selected_pair_excess_sites": verification["selected_pair_excess_sites"],
        "independently_verified": True,
    }
    if selection.validation_mode == "none":
        result.status = "hurdler_compatible_molecular"
        result.message = "The selected exact route is confirmed molecularly; IDT was not called."
        result.termination_reason = "selected_route_molecularly_verified"
        return result
    if selection.validation_mode == "batch":
        result.status = "bulk_export_unvalidated"
        result.message = "Exact purchase sequences are ready for Bulk Input; no IDT HTTP request was made."
        result.termination_reason = "bulk_export_ready"
        return result
    if idt_scorer is None:
        result.status = "idt_api_error"
        result.message = "Live IDT validation was requested without a configured scorer."
        result.termination_reason = "missing_idt_scorer"
        return result
    emit_progress(
        progress_callback,
        stage="idt",
        status="started",
        message="Scoring the exact whole target and selected purchase sequences",
    )
    try:
        if PRIMER_PAIR_CORE_THRESHOLD_BP <= len(query_result.target_sequence) <= DEFAULT_GBLOCK_MAX_BP:
            try:
                whole = _score_one_purchase(
                    idt_scorer,
                    query_result.query["sequence_id"] + "|whole_target_diagnostic",
                    query_result.target_sequence,
                )
                result.whole_target_idt = {
                    "diagnostic_only": True,
                    "score": whole["idt_complexity_score"],
                    "status": whole["idt_status"],
                    "response_sha256": whole.get("idt_response_sha256", ""),
                    "positive_rules": whole.get("idt_positive_score_names_json", "[]"),
                }
            except RuntimeError:
                result.whole_target_idt = {
                    "diagnostic_only": True,
                    "status": "diagnostic_api_error",
                }
            except ValueError:
                result.whole_target_idt = {
                    "diagnostic_only": True,
                    "status": "diagnostic_score_error",
                }
        elif len(query_result.target_sequence) > DEFAULT_GBLOCK_MAX_BP:
            result.whole_target_idt = {
                "diagnostic_only": True,
                "status": "not_scored_target_over_3000bp",
            }
        scored_by_sha: dict[str, dict[str, Any]] = {}
        all_passed = True
        scored_purchase_count = 0
        unscored_purchase_count = 0
        for ordinal, fragment in enumerate(result.purchase_fragments, start=1):
            product = str(fragment.get("product_type", ""))
            if product in {"annealed_sticky_end_primer_pair", "duplexed_seed_oligo_pair"}:
                fragment.update(
                    {
                        "idt_status": "unscored_primer_pair",
                        "idt_score": None,
                        "idt_response_sha256": "",
                        "idt_accepted": None,
                    }
                )
                unscored_purchase_count += 1
                continue
            sequence = str(fragment.get("purchase_sequence", ""))
            sequence_sha = _sha(sequence)
            if sequence_sha not in scored_by_sha:
                scored_by_sha[sequence_sha] = _score_one_purchase(
                    idt_scorer,
                    str(fragment.get("fragment_id", f"purchase_{ordinal}")),
                    sequence,
                )
            scored = scored_by_sha[sequence_sha]
            scored_purchase_count += 1
            accepted = bool(scored["idt_explicit_pass"])
            all_passed = all_passed and accepted
            fragment.update(
                {
                    "idt_policy": IDT_SCORE_POLICY,
                    "idt_status": scored["idt_status"],
                    "idt_score": scored["idt_complexity_score"],
                    "idt_response_sha256": scored.get("idt_response_sha256", ""),
                    "idt_scored_sequence_sha256": scored.get("idt_scored_sequence_sha256", sequence_sha),
                    "idt_positive_score_names_json": scored.get("idt_positive_score_names_json", "[]"),
                    "idt_rule_details_json": scored.get("idt_rule_details_json", "[]"),
                    "idt_accepted": accepted,
                }
            )
            result.idt_audit.append(
                {
                    "fragment_id": fragment.get("fragment_id", f"purchase_{ordinal}"),
                    "length_bp": len(sequence),
                    "dna_sha256": sequence_sha,
                    "score": scored["idt_complexity_score"],
                    "accepted": accepted,
                    "positive_rules": scored.get("idt_positive_score_names_json", "[]"),
                    "response_sha256": scored.get("idt_response_sha256", ""),
                    "policy": IDT_SCORE_POLICY,
                }
            )
        if not all_passed:
            result.status = "idt_rejected_route"
            result.message = "The molecular route remains valid, but at least one exact purchase sequence failed IDT scoring."
            result.termination_reason = "purchase_fragment_idt_rejected"
        elif unscored_purchase_count:
            result.status = "idt_unscored_primer_route"
            result.message = (
                "At least one short primer-pair purchase is outside the IDT gBlocks scorer; "
                "the molecular route remains valid but is not labeled IDT accepted."
            )
            result.termination_reason = "route_contains_unscored_primer_pairs"
        elif scored_purchase_count:
            result.status = "idt_accepted_route"
            result.message = "Every scored exact purchase sequence passed the strict IDT score-sum threshold."
            result.termination_reason = "all_scored_purchases_idt_accepted"
        else:
            result.status = "idt_score_error"
            result.message = "The selected route contains no purchase sequence eligible for a valid IDT decision."
            result.termination_reason = "no_idt_eligible_purchase_sequence"
    except RuntimeError as exc:
        result.status = "idt_api_error"
        result.message = "IDT API validation failed without changing the molecular compatibility result."
        result.termination_reason = str(exc)
    except ValueError as exc:
        result.status = "idt_score_error"
        result.message = str(exc)
        result.termination_reason = "invalid_idt_score_structure"
    emit_progress(
        progress_callback,
        stage="idt",
        status="completed" if result.status == "idt_accepted_route" else "failed",
        message=result.message,
        details={"status": result.status},
    )
    return result


def confirm_best_exact_dna_route(
    query_result: ExactDNAResult,
    selection: ExactDNASelection,
    *,
    idt_scorer: ComplexityScorer | None = None,
    plasmid_database: PlasmidReferenceDatabase | None = None,
    plasmid_reference_path: str | Path | None = None,
    progress_callback: ProgressCallback | None = None,
) -> ExactDNAResult:
    """Try alternate exact routes after IDT rejection without changing the target.

    Alternatives remain within the user's confirmed RE-pair, plasmid profile,
    and vector cut scheme. Authentication and malformed-score failures are
    systemic, so they stop the retry loop immediately.
    """
    if selection.validation_mode != "api":
        return confirm_exact_dna_route(
            query_result,
            selection,
            idt_scorer=idt_scorer,
            plasmid_database=plasmid_database,
            plasmid_reference_path=plasmid_reference_path,
            progress_callback=progress_callback,
        )
    pair = (selection.site_i_enzyme, selection.site_ii_enzyme)
    candidates = [
        row for row in query_result.route_candidates
        if (not selection.plasmid_profile or row["profile_id"] == selection.plasmid_profile)
        and (not selection.cut_scheme_id or row["scheme_id"] == selection.cut_scheme_id)
        and (
            not all(pair)
            or any(
                (item["site_i_enzyme"], item["site_ii_enzyme"]) == pair
                for item in row["pairs"]
            )
        )
    ]
    candidates.sort(key=lambda row: row["route_id"] != selection.route_id)
    if not candidates:
        raise ValueError("No current route matches the confirmed RE/plasmid selection")
    attempts: list[dict[str, Any]] = []
    last: ExactDNAResult | None = None
    for ordinal, route in enumerate(candidates, start=1):
        emit_progress(
            progress_callback,
            stage="idt_route_selection",
            status="attempt_started",
            message=f"Trying exact molecular route {ordinal}/{len(candidates)}",
            details={"route_id": route["route_id"]},
        )
        attempt = confirm_exact_dna_route(
            query_result,
            ExactDNASelection(
                route_id=str(route["route_id"]),
                validation_mode="api",
                plasmid_profile=str(route["profile_id"]),
                cut_scheme_id=str(route["scheme_id"]),
                site_i_enzyme=pair[0],
                site_ii_enzyme=pair[1],
            ),
            idt_scorer=idt_scorer,
            plasmid_database=plasmid_database,
            plasmid_reference_path=plasmid_reference_path,
            progress_callback=progress_callback,
        )
        failed_fragments = [
            str(row.get("fragment_id", ""))
            for row in attempt.purchase_fragments
            if row.get("idt_accepted") is not True
        ]
        attempts.append(
            {
                "attempt": ordinal,
                "route_id": route["route_id"],
                "status": attempt.status,
                "failed_fragment_ids": failed_fragments,
                "termination_reason": attempt.termination_reason,
            }
        )
        attempt.route_attempts = list(attempts)
        last = attempt
        if attempt.status == "idt_accepted_route":
            return attempt
        if attempt.status in {"idt_api_error", "idt_score_error"}:
            break
    assert last is not None
    last.route_attempts = attempts
    last.message += " Bulk Input files will be exported for manual validation."
    return last


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    columns = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        if not columns:
            return
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, sort_keys=True)
                    if isinstance(value, (dict, list, tuple, set)) else value
                    for key, value in row.items()
                }
            )


def _bulk_records(fragments: Sequence[Mapping[str, Any]]) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    seen_sequences: set[str] = set()
    for index, row in enumerate(fragments, start=1):
        name = str(row.get("fragment_id", f"purchase_{index}"))
        product = str(row.get("product_type", ""))
        candidates = (
            [
                (name + "_forward", str(row.get("purchase_sequence", ""))),
                (name + "_reverse", str(row.get("secondary_purchase_sequence", ""))),
            ]
            if product == "duplexed_seed_oligo_pair" else
            [
                (name + "_forward", str(row.get("primer_forward_5to3", ""))),
                (name + "_reverse", str(row.get("primer_reverse_5to3", ""))),
            ]
            if product == "annealed_sticky_end_primer_pair" else
            [(name, str(row.get("purchase_sequence", "")))]
        )
        for candidate in candidates:
            if candidate[1] and candidate[1] not in seen_sequences:
                seen_sequences.add(candidate[1])
                records.append(candidate)
    return records


def _user_purchase_rows(
    result: ExactDNAResult,
) -> list[dict[str, Any]]:
    step_usage: dict[str, list[int]] = {}
    for step in result.cloning_steps:
        for fragment_id in str(step.get("purchase_fragment_ids", "")).split(";"):
            if fragment_id:
                step_usage.setdefault(fragment_id, []).append(int(step["step"]))
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, fragment in enumerate(result.purchase_fragments, start=1):
        fragment_id = str(fragment.get("fragment_id", f"purchase_{index}"))
        product = str(fragment.get("product_type", ""))
        candidates = (
            [
                (fragment_id + "_forward", str(fragment.get("purchase_sequence", ""))),
                (fragment_id + "_reverse", str(fragment.get("secondary_purchase_sequence", ""))),
            ]
            if product in {"duplexed_seed_oligo_pair", "annealed_sticky_end_primer_pair"}
            else [(fragment_id, str(fragment.get("purchase_sequence", "")))]
        )
        for name, sequence in candidates:
            if not sequence or sequence in seen:
                continue
            seen.add(sequence)
            rows.append(
                {
                    "Name": name,
                    "Sequence": sequence,
                    "Length_bp": len(sequence),
                    "Product_type": product,
                    "Used_in_cloning_steps": ";".join(
                        str(value) for value in sorted(set(step_usage.get(fragment_id, [])))
                    ),
                    "IDT_status": fragment.get("idt_status", "not_run"),
                    "IDT_score": fragment.get("idt_score", ""),
                }
            )
    return rows


def write_exact_dna_outputs(result: ExactDNAResult, output_dir: str | Path) -> dict[str, str]:
    """Write a concise cloning package plus a separate technical audit ZIP."""
    from .exact_dna_artifacts import write_exact_dna_genbanks

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    route = result.selected_route or {}
    pair_text = ";".join(
        f"{row['site_i_enzyme']}/{row['site_ii_enzyme']}"
        for row in route.get("pairs", [])
    )
    summary_row = {
        "sequence_id": result.query["sequence_id"],
        "status": result.status,
        "message": result.message,
        "target_length_bp": result.target_length_bp,
        "cloning_step_count": len(result.cloning_steps),
        "unique_purchase_sequence_count": len(_bulk_records(result.purchase_fragments)),
        "plasmid_profile": route.get("profile_id", ""),
        "cut_scheme": route.get("cut_scheme", ""),
        "restriction_enzyme_pairs": pair_text,
        "independent_assembly_verification": (
            "passed" if result.independent_verification.get("passed") else "not_passed"
        ),
        "idt_decision": (
            "all_purchase_fragments_accepted"
            if result.status == "idt_accepted_route"
            else "manual_bulk_validation_required"
        ),
        "order_submitted": False,
    }
    summary_csv = destination / "cloning_summary.csv"
    _write_rows(summary_csv, [summary_row])
    paths["cloning_summary.csv"] = str(summary_csv)
    summary_md = destination / "cloning_summary.md"
    summary_md.write_text(
        "# Exact-DNA HURDLER cloning summary\n\n"
        f"- Status: `{result.status}`\n"
        f"- Cloning steps: **{len(result.cloning_steps)}**\n"
        f"- Plasmid: **{route.get('profile_id', 'not selected')}**\n"
        f"- RE pair(s): **{pair_text or 'not selected'}**\n"
        f"- IDT: **{summary_row['idt_decision']}**\n\n"
        "No order was submitted. Follow `cloning_steps.csv` in numeric order.\n"
    )
    paths["cloning_summary.md"] = str(summary_md)
    steps_path = destination / "cloning_steps.csv"
    _write_rows(steps_path, result.cloning_steps)
    paths["cloning_steps.csv"] = str(steps_path)

    purchase_rows = _user_purchase_rows(result)
    purchase_csv = destination / "purchase_fragments.csv"
    _write_rows(purchase_csv, purchase_rows)
    paths["purchase_fragments.csv"] = str(purchase_csv)
    purchase_fasta = destination / "purchase_fragments.fasta"
    purchase_fasta.write_text(
        "".join(f">{row['Name']}\n{row['Sequence']}\n" for row in purchase_rows)
    )
    paths["purchase_fragments.fasta"] = str(purchase_fasta)

    paths.update(write_exact_dna_genbanks(result, destination, include_manifest=False))
    if result.final_insert_sequence:
        final_insert = destination / "final_exact_insert.fasta"
        final_insert.write_text(
            f">{result.query['sequence_id']}_final_exact_insert\n{result.final_insert_sequence}\n"
        )
        paths["final_insert_fasta"] = str(final_insert)
    if result.final_plasmid:
        plasmid = destination / "final_plasmid.fasta"
        plasmid.write_text(
            f">{result.query['sequence_id']}_final_circular_plasmid\n"
            f"{result.final_plasmid['final_plasmid_sequence']}\n"
        )
        paths["final_plasmid_fasta"] = str(plasmid)
    records = [(str(row["Name"]), str(row["Sequence"])) for row in purchase_rows]
    if records and result.status == "idt_accepted_route":
        order_csv = destination / "order_ready_fragments.csv"
        _write_rows(order_csv, purchase_rows)
        paths["order_ready_fragments.csv"] = str(order_csv)
        order_fasta = destination / "order_ready_fragments.fasta"
        order_fasta.write_text(
            "".join(f">{name}\n{sequence}\n" for name, sequence in records)
        )
        paths["order_ready_fragments.fasta"] = str(order_fasta)
    elif records:
        for filename, delimiter in (
            ("idt_bulk_input.csv", ","),
            ("idt_bulk_input.tsv", "\t"),
        ):
            path = destination / filename
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter=delimiter)
                writer.writerow(["Name", "Sequence"])
                writer.writerows(records)
            paths[filename] = str(path)
        fasta = destination / "idt_bulk_input.fasta"
        fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))
        paths["idt_bulk_input.fasta"] = str(fasta)

    audit_zip = destination.parent / f"{destination.name}_technical_audit.zip"
    with tempfile.TemporaryDirectory(prefix="hurdler_exact_audit_") as temporary:
        audit_root = Path(temporary)
        write_json_atomic(result.to_dict(), audit_root / "exact_dna_design_summary.json")
        for filename, rows in (
            ("restriction_hits.csv", result.restriction_hits),
            ("pair_candidates.csv", result.pair_candidates),
            ("route_candidates.csv", result.route_candidates),
            ("transitions.csv", result.transitions),
            ("latent_transitions.csv", result.latent_transitions),
            ("restoration_segments.csv", result.restoration_segments),
            ("route_attempts.csv", result.route_attempts),
        ):
            _write_rows(audit_root / filename, rows)
        (audit_root / "idt_audit.jsonl").write_text(
            "".join(json.dumps(row, sort_keys=True) + "\n" for row in result.idt_audit)
        )
        bundle = getattr(result, "_assembly_bundle", None)
        if bundle:
            write_json_atomic(bundle, audit_root / "independent_assembly_verification.json")
        write_json_atomic(
            {
                "schema_version": EXACT_DNA_SCHEMA_VERSION,
                "created_at": utc_now(),
                "status": result.status,
                "target_sequence_sha256": result.target_sequence_sha256,
                "plasmid_reference_version": PLASMID_REFERENCE_VERSION,
                "idt_policy": IDT_SCORE_POLICY,
                "credentials_persisted": False,
                "ordering_performed": False,
                "main_files": {
                    key: {"path": Path(value).name, "sha256": sha256_file(value)}
                    for key, value in paths.items()
                },
            },
            audit_root / "run_manifest.json",
        )
        with zipfile.ZipFile(audit_zip, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(audit_root.iterdir()):
                archive.write(path, path.name)
    paths["technical_audit_zip"] = str(audit_zip)
    return paths
