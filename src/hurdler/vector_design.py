"""Shared controller for the annotation-aware HURDLER designer v2."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass, field, fields, replace
from pathlib import Path
from typing import Any, Mapping, Protocol, Sequence

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
    GA_SCORE_PROFILE,
    adaptive_copy_search,
    adjust_ga_score_profile_from_idt,
    genetic_refine_dna,
    load_restriction_sites,
)
from .io import sha256_file, utc_now, write_json_atomic
from .optimization import diversify_codons, load_codon_weights, reverse_complement, translate_dna
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


DESIGN_SCHEMA_VERSION_V2 = "vector-aware-hurdler-designer-v2"
INPUT_MODES = ("split", "full")
VALIDATION_MODES = ("none", "api", "batch")


class ComplexityScorer(Protocol):
    def score(self, name: str, sequence: str) -> dict[str, Any]: ...


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
    max_repeat_copies: int | None = None
    population_size: int = 16
    mutation_rate: float = 0.08
    crossover_rate: float = 0.75
    elite_fraction: float = 0.15
    seed: int = 42
    generation_schedule: tuple[int, ...] = GENERATION_SCHEDULE
    score_weights: Mapping[str, float] = field(default_factory=lambda: dict(GA_SCORE_PROFILE))
    auto_adjust_weights_from_idt: bool = True
    max_purchase_bp: int = GBLOCK_MAX_BP

    def __post_init__(self) -> None:
        if self.schema_version != DESIGN_SCHEMA_VERSION_V2:
            raise ValueError(f"schema_version must be {DESIGN_SCHEMA_VERSION_V2!r}")
        if self.validation_mode not in VALIDATION_MODES:
            raise ValueError(f"validation_mode must be one of {VALIDATION_MODES}")
        if int(self.population_size) < 2:
            raise ValueError("population_size must be at least two")
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
    idt_audit: list[dict[str, Any]] = field(default_factory=list)
    final_plasmid: dict[str, Any] | None = None
    final_protein_sequence: str = ""
    final_dna_sequence: str = ""
    termination_reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


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


def design_query(
    query: CompatibilityQuery,
    *,
    protein_index: ProteinPatternIndex | None = None,
    protein_index_dir: str | Path | None = None,
    plasmid_database: PlasmidReferenceDatabase | None = None,
    plasmid_reference_path: str | Path | None = None,
) -> DesignResultV2:
    analysis, boundary, units = _query_boundary(query)
    request_payload = asdict(query)
    for key in ("site_i_allowlist", "site_ii_allowlist", "site_iii_allowlist", "plasmid_allowlist"):
        request_payload[key] = list(request_payload[key])
    if boundary is None:
        return DesignResultV2(
            DESIGN_SCHEMA_VERSION_V2,
            "needs_boundary_confirmation",
            "Confirm or edit the 1-based inclusive repeat boundary and period; the full input sequence will not be homogenized.",
            request_payload,
            boundary_analysis=analysis,
        )
    index = protein_index or ProteinPatternIndex.load(protein_index_dir or bundled_protein_index_dir())
    candidates = _protein_candidates(units, query, index)
    if not candidates:
        return DesignResultV2(
            DESIGN_SCHEMA_VERSION_V2,
            "no_hurdler_pair_match",
            "No allowed protein-level Site-I/Site-II RE pair matches every required repeat unit.",
            request_payload,
            boundary_analysis=analysis,
            confirmed_boundary=boundary.to_dict(),
        )
    database = plasmid_database or load_plasmid_reference(plasmid_reference_path)
    routes: list[dict[str, Any]] = []
    for candidate in candidates:
        for scheme in database.schemes:
            route = _scheme_route(database, scheme, candidate, query)
            if route is not None:
                routes.append(route)
    # Cutter reuse is a final fallback, never co-ranked with unmodified routes.
    if any(not row["cutter_reuse"] for row in routes):
        routes = [row for row in routes if not row["cutter_reuse"]]
    routes.sort(key=lambda row: (row["cutter_reuse"], row["restoration_length_bp"], row["profile_id"], row["scheme_id"], row["candidate_id"]))
    for rank, row in enumerate(routes, 1):
        row["rank"] = rank
    status = "compatible_unoptimized" if routes else "no_vector_route"
    message = (
        "At least one protein RE pair and annotation-safe plasmid cut route is available; select a route before optimization."
        if routes
        else "Protein-level pairs exist, but none is clean on a retained annotated plasmid backbone under the selected cutter policy."
    )
    target = query.n_cap + query.repeat_module * int(query.repeat_copies) + query.c_cap if query.input_mode == "split" else query.full_protein_sequence
    return DesignResultV2(
        DESIGN_SCHEMA_VERSION_V2,
        status,
        message,
        request_payload,
        boundary_analysis=analysis,
        confirmed_boundary=boundary.to_dict(),
        protein_candidates=candidates,
        vector_routes=routes,
        final_protein_sequence=target,
    )


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
        idt_pass, scored_fragments, audits = _score_fragments(fragments, idt_scorer)
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


def design_construct_v2(
    request: DesignRequestV2,
    *,
    protein_index_dir: str | Path | None = None,
    plasmid_reference_path: str | Path | None = None,
    idt_scorer: ComplexityScorer | None = None,
) -> DesignResultV2:
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
    root = ProjectPaths.discover().root
    weights = load_codon_weights(root / "data" / "reference_output" / "codon_usage.csv")
    recognition_sites = load_restriction_sites(root / "data" / "reference_output" / "restriction_enzyme.csv")
    if request.query.input_mode == "split" and request.max_repeat_copies is not None:
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
            ),
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
            idt_pass, fragments, audit = _score_fragments(fragments, idt_scorer)
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
    if result.final_plasmid:
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO

        plasmid_fasta = destination / "final_plasmid.fasta"
        plasmid_fasta.write_text(
            f">{result.request['query']['sequence_id']}_final_circular_plasmid\n"
            f"{result.final_plasmid['final_plasmid_sequence']}\n"
        )
        paths["final_plasmid_fasta"] = str(plasmid_fasta)
        record = SeqRecord(
            Seq(str(result.final_plasmid["final_plasmid_sequence"])),
            id=f"{result.request['query']['sequence_id']}_v2",
            name="HURDLER_v2",
            description="annotation-aware HURDLER design; design files only",
        )
        record.annotations.update({"molecule_type": "DNA", "topology": "circular"})
        record.features.append(
            SeqFeature(
                FeatureLocation(
                    int(result.final_plasmid["cds_start_0based"]),
                    int(result.final_plasmid["cds_end_0based_exclusive"]),
                    strand=1,
                ),
                type="CDS",
                qualifiers={
                    "label": ["optimized_repeat_protein_CDS"],
                    "translation": [result.final_protein_sequence],
                },
            )
        )
        plasmid_genbank = destination / "final_plasmid.gb"
        SeqIO.write(record, plasmid_genbank, "genbank")
        paths["final_plasmid_genbank"] = str(plasmid_genbank)
    all_fragments = [*result.primary_fragments, *result.secondary_fragments]
    if result.status == "optimized_unvalidated_batch":
        paths.update({f"idt_bulk_{key}": value for key, value in _batch_fragments(destination, all_fragments).items()})
    audit = destination / "idt_audit.jsonl"
    audit.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in result.idt_audit))
    paths["idt_audit_jsonl"] = str(audit)
    manifest = {
        "schema_version": DESIGN_SCHEMA_VERSION_V2,
        "created_at": utc_now(),
        "status": result.status,
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
