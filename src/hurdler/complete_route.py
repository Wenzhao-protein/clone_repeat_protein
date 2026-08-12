"""End-to-end repeat-array assembly from a purchasable exact seed.

The legacy arbitrary-DNA planner proves one replacement against an assumed
recipient.  This module treats exact copy counts as graph states and only
reports a target after a path from a purchasable seed has been verified.  A
single plasmid is fixed along a path; the Site-I/Site-II pair may change.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from .constants import PLASMIDS
from .dna_assembly import (
    DEFAULT_GBLOCK_MAX_BP,
    DNA_COMPLETE_ROUTE_VERSION,
    PRIMER_PAIR_CORE_THRESHOLD_BP,
    EnzymeGeometry,
    TargetRecord,
    _product_type,
    _select_type_iis,
    enumerate_active_latent_pairs,
    exact_repeat_gain_for_interval,
    load_enzyme_catalog,
    plan_target,
    reverse_complement,
    scan_re_sites,
    validate_dna,
)
from .idt import IDT_SCORE_POLICY, IDTComplexityScorer
from .io import sha256_file, utc_now, write_json_atomic


TARGET_COPY_COUNTS = (2, 4, 8, 16, 32)
MAX_COPY_COUNT = max(TARGET_COPY_COUNTS)


@dataclass(frozen=True)
class SeedEvidence:
    copy_count: int
    core_sequence: str
    product_type: str
    purchase_sequence: str
    secondary_purchase_sequence: str
    purchase_sequence_count: int
    purchase_length_bp: int
    purchase_sha256: str
    idt_status: str
    idt_score: float | None
    idt_response_sha256: str
    accepted: bool

    def to_dict(self) -> dict[str, Any]:
        return dict(self.__dict__)


@dataclass(frozen=True)
class TransitionEvidence:
    recipient_copy_count: int
    result_copy_count: int
    donor_copy_count: int
    plasmid: str
    site_i_enzyme: str
    site_ii_enzyme: str
    route_id: str
    hurdle_steps: int
    purchase_sha256s: tuple[str, ...]
    total_purchase_bp: int
    maximum_idt_score: float | None
    whole_target_idt_status: str
    whole_target_idt_score: float | None
    whole_target_idt_response_sha256: str
    final_target_exact: bool
    all_purchase_fragments_accepted: bool
    route_row: dict[str, Any]
    step_rows: tuple[dict[str, Any], ...]
    fragment_rows: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class _Path:
    plasmid: str
    copy_count: int
    transitions: tuple[TransitionEvidence, ...]
    purchase_sha256s: frozenset[str]
    total_purchase_bp: int
    maximum_idt_score: float
    pair_change_count: int
    last_pair: tuple[str, str] | None
    last_fragment_ceiling_bp: int | None

    @property
    def hurdle_steps(self) -> int:
        return sum(edge.hurdle_steps for edge in self.transitions)

    def ranking_key(self) -> tuple[Any, ...]:
        return (
            1 + self.hurdle_steps,  # seed installation plus HURDLER cycles
            len(self.purchase_sha256s),
            self.total_purchase_bp,
            self.pair_change_count,
            self.maximum_idt_score,
            self.plasmid,
            tuple(
                (edge.site_i_enzyme, edge.site_ii_enzyme, edge.route_id)
                for edge in self.transitions
            ),
        )


TransitionProvider = Callable[[int, int], list[TransitionEvidence]]


def _safe_idt_score(
    scorer: IDTComplexityScorer | None,
    name: str,
    sequence: str,
    *,
    require_idt: bool,
) -> dict[str, Any]:
    if not require_idt:
        return {
            "idt_status": "not_required",
            "idt_explicit_pass": True,
            "idt_complexity_score": None,
            "idt_response_sha256": "",
        }
    if scorer is None:
        return {
            "idt_status": "api_unclassified",
            "idt_explicit_pass": False,
            "idt_complexity_score": None,
            "idt_response_sha256": "",
        }
    try:
        return scorer.score(name, sequence)
    except Exception:
        return {
            "idt_status": "api_failure",
            "idt_explicit_pass": False,
            "idt_complexity_score": None,
            "idt_response_sha256": "",
        }


def _whole_target_evidence(
    scorer: IDTComplexityScorer | None,
    name: str,
    sequence: str,
    *,
    require_idt: bool,
) -> dict[str, Any]:
    if len(sequence) < PRIMER_PAIR_CORE_THRESHOLD_BP:
        return {
            "idt_status": "not_applicable_primer_pair_under_90bp",
            "idt_explicit_pass": True,
            "idt_complexity_score": None,
            "idt_response_sha256": "",
        }
    return _safe_idt_score(
        scorer, name, sequence, require_idt=require_idt
    )


def find_shortest_purchasable_seed(
    unit_sequence: str,
    *,
    max_copy_count: int = MAX_COPY_COUNT,
    idt_scorer: IDTComplexityScorer | None = None,
    require_idt: bool = True,
) -> SeedEvidence | None:
    """Return the smallest exact subarray that satisfies the purchase policy."""
    unit = validate_dna(unit_sequence)
    for copies in range(1, max_copy_count + 1):
        core = unit * copies
        if len(core) < PRIMER_PAIR_CORE_THRESHOLD_BP:
            forward = core
            reverse = reverse_complement(core)
            purchase_key = f"duplexed_seed_oligo_pair|{forward}|{reverse}"
            return SeedEvidence(
                copy_count=copies,
                core_sequence=core,
                product_type="duplexed_seed_oligo_pair",
                purchase_sequence=forward,
                secondary_purchase_sequence=reverse,
                purchase_sequence_count=2,
                purchase_length_bp=len(forward) + len(reverse),
                purchase_sha256=hashlib.sha256(purchase_key.encode()).hexdigest(),
                idt_status="not_applicable_primer_pair_under_90bp",
                idt_score=None,
                idt_response_sha256="",
                accepted=True,
            )
        product_type = _product_type(len(core))
        if product_type is None or len(core) > DEFAULT_GBLOCK_MAX_BP:
            continue
        result = _safe_idt_score(
            idt_scorer,
            f"seed|copies={copies}",
            core,
            require_idt=require_idt,
        )
        if result.get("idt_explicit_pass") is not True:
            continue
        return SeedEvidence(
            copy_count=copies,
            core_sequence=core,
            product_type=product_type,
            purchase_sequence=core,
            secondary_purchase_sequence="",
            purchase_sequence_count=1,
            purchase_length_bp=len(core),
            purchase_sha256=hashlib.sha256(
                f"{product_type}|{core}".encode()
            ).hexdigest(),
            idt_status=str(result.get("idt_status", "passed")),
            idt_score=(
                float(result["idt_complexity_score"])
                if result.get("idt_complexity_score") is not None else None
            ),
            idt_response_sha256=str(result.get("idt_response_sha256", "")),
            accepted=True,
        )
    return None


def _default_transition_provider(
    base: TargetRecord,
    geometries: dict[str, EnzymeGeometry],
    plasmids: pd.DataFrame,
    *,
    idt_scorer: IDTComplexityScorer | None,
    require_idt: bool,
    top_routes: int,
    purchase_ceilings: tuple[int, ...] = (DEFAULT_GBLOCK_MAX_BP, 200, 60, 55),
) -> TransitionProvider:
    unit = validate_dna(base.unit_sequence)
    cache: dict[tuple[int, int], list[TransitionEvidence]] = {}
    diagnostics: dict[tuple[int, int], dict[str, Any]] = {}
    prepared_results: set[int] = set()

    def prepare_result(result_copies: int) -> None:
        """Scan one final state once, then split routes by exact unit gain."""
        if result_copies in prepared_results:
            return
        prepared_results.add(result_copies)
        edge_target = TargetRecord(
            target_id=f"{base.element_id or base.target_id}|copies={result_copies}",
            sequence=unit * result_copies,
            cohort=base.cohort,
            architecture="exact_tandem",
            source_url=base.source_url,
            source_accession=base.source_accession,
            unit_sequence=unit,
            copy_count=result_copies,
            notes=base.notes,
            source_database=base.source_database,
            element_id=base.element_id,
        )
        edge_sequence = unit * result_copies
        precomputed_hits = scan_re_sites(
            edge_target.target_id, edge_sequence, geometries
        )
        precomputed_schemes = enumerate_active_latent_pairs(
            precomputed_hits, plasmids
        )
        exact_schemes = [
            scheme
            for scheme in precomputed_schemes
            if exact_repeat_gain_for_interval(
                edge_target,
                int(scheme["replacement_start"]),
                int(scheme["replacement_end"]),
            ) is not None
        ]
        long_fragment_possible = any(
            _select_type_iis(
                scheme["left_hit"].ovhg, edge_sequence, geometries
            ) is not None
            and _select_type_iis(
                scheme["right_hit"].ovhg, edge_sequence, geometries
            ) is not None
            for scheme in exact_schemes
        )
        initial_ceilings = (
            tuple(
                ceiling
                for ceiling in dict.fromkeys(purchase_ceilings)
                if ceiling >= PRIMER_PAIR_CORE_THRESHOLD_BP or ceiling == 60
            )
            if long_fragment_possible else
            (60,)
        )

        def plan_ceiling(ceiling: int) -> dict[str, pd.DataFrame]:
            return plan_target(
                edge_target,
                geometries,
                plasmids,
                idt_scorer=None,
                require_idt=False,
                top_routes=top_routes,
                max_purchase_bp=ceiling,
                enumerate_repeat_unit_gains=True,
                _precomputed_hits=precomputed_hits,
                _precomputed_schemes=precomputed_schemes,
            )

        planned_by_ceiling = [
            plan_ceiling(ceiling) for ceiling in initial_ceilings
        ]
        expected_gains = {
            int(gain)
            for scheme in exact_schemes
            if (
                gain := exact_repeat_gain_for_interval(
                    edge_target,
                    int(scheme["replacement_start"]),
                    int(scheme["replacement_end"]),
                )
            ) is not None
        }
        passing_gains: set[int] = set()
        for result in planned_by_ceiling:
            result_routes = result["routes"]
            if not result_routes.empty:
                passing_gains.update(
                    result_routes.loc[
                        result_routes.local_constraints_passed,
                        "repeat_unit_gain",
                    ].dropna().astype(int).tolist()
                )
        # A 55-bp layout can remove a selected site spanning a 60-bp
        # fragment boundary.  It cannot improve the primary ranking when all
        # exact gains already have a passing 60-bp or long-fragment route.
        if 55 in purchase_ceilings and expected_gains - passing_gains:
            planned_by_ceiling.append(plan_ceiling(55))
        routes = pd.concat(
            [result["routes"] for result in planned_by_ceiling if not result["routes"].empty],
            ignore_index=True,
            sort=False,
        ) if any(not result["routes"].empty for result in planned_by_ceiling) else pd.DataFrame()
        steps = pd.concat(
            [result["steps"] for result in planned_by_ceiling if not result["steps"].empty],
            ignore_index=True,
            sort=False,
        ) if any(not result["steps"].empty for result in planned_by_ceiling) else pd.DataFrame()
        fragments = pd.concat(
            [result["fragments"] for result in planned_by_ceiling if not result["fragments"].empty],
            ignore_index=True,
            sort=False,
        ) if any(not result["fragments"].empty for result in planned_by_ceiling) else pd.DataFrame()
        summary = planned_by_ceiling[0]["summary"].iloc[0]
        for gain in range(1, result_copies):
            recipient_copies = result_copies - gain
            key = (recipient_copies, result_copies)
            gain_routes = (
                routes.loc[routes.repeat_unit_gain.eq(gain)].copy()
                if not routes.empty and "repeat_unit_gain" in routes.columns
                else pd.DataFrame()
            )
            diagnostics[key] = {
                "failure_reason": (
                    "" if bool(gain_routes.get("passed", pd.Series(dtype=bool)).any())
                    else str(summary.failure_reason)
                ),
                "candidate_pair_count": int(len(gain_routes)),
                "whole_target_idt_status": "pending_independent_score",
                "failed_idt_routes": int(
                    gain_routes.status.eq("failed_idt").sum()
                ) if not gain_routes.empty else 0,
                "local_constraint_routes": int(
                    gain_routes.local_constraints_passed.sum()
                ) if not gain_routes.empty else 0,
            }
            results: list[TransitionEvidence] = []
            if not gain_routes.empty:
                ordered_routes = gain_routes.loc[gain_routes.passed].sort_values(
                    [
                        "step_count",
                        "unique_fragment_count",
                        "total_purchase_bp",
                        "site_i_enzyme",
                        "site_ii_enzyme",
                        "plasmid",
                        "direction",
                        "site_i_position",
                        "site_ii_position",
                    ],
                    kind="mergesort",
                )
                for route in ordered_routes.itertuples(index=False):
                    route_steps = steps.loc[steps.route_id.eq(route.route_id)]
                    route_fragments = fragments.loc[
                        fragments.route_id.eq(route.route_id)
                    ]
                    results.append(
                        TransitionEvidence(
                            recipient_copy_count=recipient_copies,
                            result_copy_count=result_copies,
                            donor_copy_count=gain,
                            plasmid=str(route.plasmid),
                            site_i_enzyme=str(route.site_i_enzyme),
                            site_ii_enzyme=str(route.site_ii_enzyme),
                            route_id=str(route.route_id),
                            hurdle_steps=int(route.step_count),
                            purchase_sha256s=tuple(
                                route_fragments.purchase_sha256.astype(str).tolist()
                            ),
                            total_purchase_bp=int(route.total_purchase_bp),
                            maximum_idt_score=(
                                float(route.maximum_idt_score)
                                if pd.notna(route.maximum_idt_score) else None
                            ),
                            whole_target_idt_status="pending_independent_score",
                            whole_target_idt_score=None,
                            whole_target_idt_response_sha256="",
                            final_target_exact=bool(route.final_sequence_exact),
                            all_purchase_fragments_accepted=bool(
                                route.local_constraints_passed
                            ),
                            route_row=route._asdict(),
                            step_rows=tuple(route_steps.to_dict("records")),
                            fragment_rows=tuple(
                                route_fragments.to_dict("records")
                            ),
                        )
                    )
            cache[key] = results

    def provide(recipient_copies: int, result_copies: int) -> list[TransitionEvidence]:
        key = (recipient_copies, result_copies)
        if key in cache:
            return cache[key]
        if not 1 <= recipient_copies < result_copies <= MAX_COPY_COUNT:
            cache[key] = []
            return []
        prepare_result(result_copies)
        return cache.get(key, [])

    provide.diagnostics = diagnostics  # type: ignore[attr-defined]
    return provide


def _extend_path(path: _Path, edge: TransitionEvidence) -> _Path:
    pair = (edge.site_i_enzyme, edge.site_ii_enzyme)
    scores = [path.maximum_idt_score]
    if edge.maximum_idt_score is not None:
        scores.append(edge.maximum_idt_score)
    return _Path(
        plasmid=path.plasmid,
        copy_count=edge.result_copy_count,
        transitions=(*path.transitions, edge),
        purchase_sha256s=frozenset(
            set(path.purchase_sha256s).union(edge.purchase_sha256s)
        ),
        total_purchase_bp=path.total_purchase_bp + edge.total_purchase_bp,
        maximum_idt_score=max(scores),
        pair_change_count=(
            path.pair_change_count
            + int(path.last_pair is not None and path.last_pair != pair)
        ),
        last_pair=pair,
        last_fragment_ceiling_bp=int(
            edge.route_row.get("fragment_purchase_ceiling_bp", 0) or 0
        ),
    )


def _validate_path_purchases(
    path: _Path,
    *,
    idt_scorer: IDTComplexityScorer | None,
    require_idt: bool,
) -> tuple[_Path, bool, list[str]]:
    """Score only purchases belonging to a candidate complete path."""
    validated_edges: list[TransitionEvidence] = []
    scores = [path.maximum_idt_score]
    failures: list[str] = []
    for edge in path.transitions:
        fragment_rows: list[dict[str, Any]] = []
        edge_scores: list[float] = []
        edge_accepted = True
        for original in edge.fragment_rows:
            fragment = dict(original)
            if fragment.get("product_type") == "annealed_sticky_end_primer_pair":
                result = {
                    "idt_status": "not_applicable_primer_pair_under_90bp",
                    "idt_explicit_pass": True,
                    "idt_complexity_score": None,
                    "idt_response_sha256": "",
                    "idt_positive_score_names_json": "[]",
                    "idt_rule_details_json": "[]",
                }
            else:
                result = _safe_idt_score(
                    idt_scorer,
                    str(fragment.get("fragment_id", edge.route_id)),
                    str(fragment.get("purchase_sequence", "")),
                    require_idt=require_idt,
                )
            accepted = result.get("idt_explicit_pass") is True
            edge_accepted = edge_accepted and accepted
            if result.get("idt_complexity_score") is not None:
                score = float(result["idt_complexity_score"])
                edge_scores.append(score)
                scores.append(score)
            fragment.update(
                {
                    "idt_policy": IDT_SCORE_POLICY,
                    "idt_status": result.get("idt_status"),
                    "idt_score": result.get("idt_complexity_score"),
                    "idt_response_sha256": result.get(
                        "idt_response_sha256", ""
                    ),
                    "idt_scored_sequence_sha256": result.get(
                        "idt_scored_sequence_sha256", ""
                    ),
                    "idt_positive_score_names_json": result.get(
                        "idt_positive_score_names_json", "[]"
                    ),
                    "idt_rule_details_json": result.get(
                        "idt_rule_details_json", "[]"
                    ),
                    "purchase_accepted": accepted,
                }
            )
            if not accepted:
                positive_rules = str(
                    result.get("idt_positive_score_names_json", "[]")
                )
                failures.append(
                    f"{fragment.get('fragment_id', 'fragment')}:"
                    f"{result.get('idt_status', 'api_unclassified')}:"
                    f"positive_rules={positive_rules}"
                )
            fragment_rows.append(fragment)
        route_row = dict(edge.route_row)
        route_row.update(
            {
                "maximum_idt_score": max(edge_scores) if edge_scores else None,
                "all_fragments_idt_passed": edge_accepted,
                "all_purchase_fragments_accepted": edge_accepted,
                "passed": bool(
                    route_row.get("local_constraints_passed", True)
                    and edge_accepted
                ),
                "status": "passed" if edge_accepted else "failed_idt",
            }
        )
        validated_edges.append(
            replace(
                edge,
                maximum_idt_score=max(edge_scores) if edge_scores else None,
                all_purchase_fragments_accepted=edge_accepted,
                route_row=route_row,
                fragment_rows=tuple(fragment_rows),
            )
        )
    validated = replace(
        path,
        transitions=tuple(validated_edges),
        maximum_idt_score=max(scores) if scores else 0.0,
    )
    return validated, not failures, failures


def search_complete_repeat_routes(
    base: TargetRecord,
    geometries: dict[str, EnzymeGeometry],
    plasmids: pd.DataFrame,
    *,
    idt_scorer: IDTComplexityScorer | None = None,
    require_idt: bool = True,
    target_copy_counts: tuple[int, ...] = TARGET_COPY_COUNTS,
    max_copy_count: int = MAX_COPY_COUNT,
    top_routes_per_transition: int = 10,
    paths_per_state: int = 2,
    max_paths_per_plasmid_state: int = 80,
    candidate_paths_per_fragmentation: int = 10,
    transition_provider: TransitionProvider | None = None,
) -> dict[str, pd.DataFrame]:
    """Find lexicographically optimal complete paths to all requested copies."""
    unit = validate_dna(base.unit_sequence)
    if not unit:
        raise ValueError("Complete repeat routing requires unit_sequence")
    seed = find_shortest_purchasable_seed(
        unit,
        max_copy_count=max_copy_count,
        idt_scorer=idt_scorer,
        require_idt=require_idt,
    )
    source_database = base.source_database or (
        base.target_id.split("|", 1)[0] if "|" in base.target_id else ""
    )
    element_id = base.element_id or base.target_id.split("|copies=", 1)[0]
    def failed_rows(reason: str) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for copies in target_copy_counts:
            target_sequence = unit * copies
            whole = _whole_target_evidence(
                idt_scorer,
                f"{source_database}|{element_id}|copies={copies}|whole_target",
                target_sequence,
                require_idt=require_idt,
            )
            rows.append(
                {
                    "version": DNA_COMPLETE_ROUTE_VERSION,
                    "source_database": source_database,
                    "element_id": element_id,
                    "unit_sequence": unit,
                    "unit_length_bp": len(unit),
                    "target_copy_count": copies,
                    "target_id": f"{source_database}|{element_id}|copies={copies}",
                    "target_length_bp": len(target_sequence),
                    "target_sequence_sha256": hashlib.sha256(
                        target_sequence.encode()
                    ).hexdigest(),
                    "seed_copy_count": seed.copy_count if seed else None,
                    "seed_product_type": seed.product_type if seed else "",
                    "complete_route_id": "",
                    "complete_route_verified": False,
                    "final_target_exact": False,
                    "hurdler_compatible": False,
                    "plasmid": "",
                    "experimental_step_count": None,
                    "hurdler_step_count": None,
                    "pair_change_count": None,
                    "unique_purchase_count": None,
                    "total_purchase_bp": None,
                    "maximum_idt_score": None,
                    "whole_target_idt_status": whole.get("idt_status"),
                    "whole_target_idt_score": whole.get(
                        "idt_complexity_score"
                    ),
                    "whole_target_idt_response_sha256": whole.get(
                        "idt_response_sha256", ""
                    ),
                    "whole_target_idt_positive_score_names_json": whole.get(
                        "idt_positive_score_names_json", "[]"
                    ),
                    "whole_target_idt_rule_details_json": whole.get(
                        "idt_rule_details_json", "[]"
                    ),
                    "idt_evidence_tier": "not_applicable",
                    "fragment_rescued_by_hurdler": False,
                    "failure_reason": reason,
                }
            )
        return rows

    if seed is None:
        return {
            "targets": pd.DataFrame(failed_rows("seed_unavailable")),
            "selected_routes": pd.DataFrame(),
            "transitions": pd.DataFrame(),
            "steps": pd.DataFrame(),
            "fragments": pd.DataFrame(),
            "seeds": pd.DataFrame(),
            "candidate_routes": pd.DataFrame(),
        }
    # Safe global rejection: every active/latent occurrence present in a
    # shorter exact array is also present in the 32-copy array.  If that
    # superset has no vector-compatible pair, none of the requested states can
    # support a HURDLER growth edge.  This avoids 31 redundant scans for the
    # dominant negative class without changing scientific outcomes.
    possible_repeat_gains: set[int] | None = None
    if transition_provider is None:
        superset_hits = scan_re_sites(
            f"{element_id}|preflight", unit * max_copy_count, geometries
        )
        superset_schemes = enumerate_active_latent_pairs(superset_hits, plasmids)
        superset_target = TargetRecord(
            target_id=f"{element_id}|preflight",
            sequence=unit * max_copy_count,
            unit_sequence=unit,
            copy_count=max_copy_count,
        )
        exact_growth_schemes = [
            scheme
            for scheme in superset_schemes
            if exact_repeat_gain_for_interval(
                superset_target,
                int(scheme["replacement_start"]),
                int(scheme["replacement_end"]),
            ) is not None
        ]
        possible_repeat_gains = {
            int(gain)
            for scheme in exact_growth_schemes
            if (
                gain := exact_repeat_gain_for_interval(
                    superset_target,
                    int(scheme["replacement_start"]),
                    int(scheme["replacement_end"]),
                )
            ) is not None
        }
        if not exact_growth_schemes:
            reason = (
                "no_active_latent_pair"
                if not superset_schemes else
                "no_exact_repeat_gain_pair"
            )
            return {
                "targets": pd.DataFrame(failed_rows(reason)),
                "selected_routes": pd.DataFrame(),
                "transitions": pd.DataFrame(),
                "steps": pd.DataFrame(),
                "fragments": pd.DataFrame(),
                "seeds": pd.DataFrame(
                    [
                        {
                            "version": DNA_COMPLETE_ROUTE_VERSION,
                            "source_database": source_database,
                            "element_id": element_id,
                            **seed.to_dict(),
                        }
                    ]
                ),
                "candidate_routes": pd.DataFrame(),
            }
        # Do not preflight Type IIS adapter availability here.  Every exact
        # growth interval is also evaluated with 60/55-bp fragmentation,
        # where each core is supplied as a complementary sticky-end primer
        # pair and therefore needs no Type IIS adapter.  Long-fragment routes
        # still enforce adapter availability inside ``plan_target()``.
    provider = transition_provider or _default_transition_provider(
        base,
        geometries,
        plasmids,
        idt_scorer=idt_scorer,
        require_idt=require_idt,
        top_routes=top_routes_per_transition,
    )
    states: dict[int, list[_Path]] = {seed.copy_count: []}
    for plasmid in PLASMIDS:
        if plasmid in plasmids.columns:
            states[seed.copy_count].append(
                _Path(
                    plasmid=plasmid,
                    copy_count=seed.copy_count,
                    transitions=(),
                    purchase_sha256s=frozenset({seed.purchase_sha256}),
                    total_purchase_bp=seed.purchase_length_bp,
                    maximum_idt_score=seed.idt_score or 0.0,
                    pair_change_count=0,
                    last_pair=None,
                    last_fragment_ceiling_bp=None,
                )
            )
    for result_copies in range(seed.copy_count + 1, max_copy_count + 1):
        candidates: list[_Path] = []
        for recipient_copies in sorted(states):
            if recipient_copies >= result_copies:
                continue
            if (
                possible_repeat_gains is not None
                and result_copies - recipient_copies
                not in possible_repeat_gains
            ):
                continue
            edges = provider(recipient_copies, result_copies)
            if not edges:
                continue
            for path in states[recipient_copies]:
                for edge in edges:
                    if edge.plasmid != path.plasmid:
                        continue
                    if not edge.final_target_exact or not edge.all_purchase_fragments_accepted:
                        continue
                    candidates.append(_extend_path(path, edge))
        grouped: dict[
            tuple[str, tuple[str, str] | None, int | None], list[_Path]
        ] = {}
        for candidate in candidates:
            grouped.setdefault(
                (
                    candidate.plasmid,
                    candidate.last_pair,
                    candidate.last_fragment_ceiling_bp,
                ),
                [],
            ).append(candidate)
        retained: list[_Path] = []
        for group in grouped.values():
            retained.extend(sorted(group, key=lambda item: item.ranking_key())[:paths_per_state])
        if retained:
            by_plasmid: dict[str, list[_Path]] = {}
            for path in retained:
                by_plasmid.setdefault(path.plasmid, []).append(path)
            states[result_copies] = [
                path
                for plasmid in sorted(by_plasmid)
                for path in sorted(
                    by_plasmid[plasmid], key=lambda item: item.ranking_key()
                )[:max_paths_per_plasmid_state]
            ]

    target_rows: list[dict[str, Any]] = []
    route_rows: list[dict[str, Any]] = []
    transition_rows: list[dict[str, Any]] = []
    step_rows: list[dict[str, Any]] = []
    fragment_rows: list[dict[str, Any]] = []
    candidate_route_rows: list[dict[str, Any]] = []
    seed_row = {
        "version": DNA_COMPLETE_ROUTE_VERSION,
        "source_database": source_database,
        "element_id": element_id,
        **seed.to_dict(),
    }
    provider_diagnostics = getattr(provider, "diagnostics", {})
    for target_copies in target_copy_counts:
        eligible = sorted(
            [path for path in states.get(target_copies, []) if path.transitions],
            key=lambda item: item.ranking_key(),
        )
        by_fragmentation: dict[tuple[int, ...], list[_Path]] = {}
        for path in eligible:
            signature = tuple(
                int(
                    edge.route_row.get("fragment_purchase_ceiling_bp", 0) or 0
                )
                for edge in path.transitions
            )
            by_fragmentation.setdefault(signature, []).append(path)
        candidates_to_validate = sorted(
            [
                path
                for paths in by_fragmentation.values()
                for path in paths[:candidate_paths_per_fragmentation]
            ],
            key=lambda item: item.ranking_key(),
        )
        validated_candidates: list[_Path] = []
        for candidate_rank, candidate in enumerate(
            candidates_to_validate, start=1
        ):
            validated, accepted, failures = _validate_path_purchases(
                candidate,
                idt_scorer=idt_scorer,
                require_idt=require_idt,
            )
            candidate_route_rows.append(
                {
                    "version": DNA_COMPLETE_ROUTE_VERSION,
                    "source_database": source_database,
                    "element_id": element_id,
                    "target_copy_count": target_copies,
                    "candidate_rank_before_idt": candidate_rank,
                    "plasmid": candidate.plasmid,
                    "hurdler_step_count": candidate.hurdle_steps,
                    "fragmentation_signature_json": json.dumps(
                        [
                            edge.route_row.get(
                                "fragment_purchase_ceiling_bp"
                            )
                            for edge in candidate.transitions
                        ]
                    ),
                    "transition_route_ids_json": json.dumps(
                        [edge.route_id for edge in candidate.transitions]
                    ),
                    "all_purchase_fragments_accepted": accepted,
                    "maximum_idt_score": validated.maximum_idt_score,
                    "failure_reasons_json": json.dumps(failures),
                }
            )
            if accepted:
                validated_candidates.append(validated)
        best = (
            min(validated_candidates, key=lambda item: item.ranking_key())
            if validated_candidates else None
        )
        target_sequence = unit * target_copies
        target_sha = hashlib.sha256(target_sequence.encode()).hexdigest()
        whole_target = _whole_target_evidence(
            idt_scorer,
            f"{source_database}|{element_id}|copies={target_copies}|whole_target",
            target_sequence,
            require_idt=require_idt,
        )
        complete = best is not None
        last_edge = best.transitions[-1] if best else None
        full_route_id = ""
        evidence_tier = "not_applicable"
        if best:
            full_route_payload = {
                "source_database": source_database,
                "element_id": element_id,
                "target_copy_count": target_copies,
                "plasmid": best.plasmid,
                "edges": [edge.route_id for edge in best.transitions],
            }
            full_route_id = "complete_" + hashlib.sha256(
                json.dumps(full_route_payload, sort_keys=True).encode()
            ).hexdigest()[:20]
            products = {seed.product_type}
            scored = seed.idt_response_sha256 != ""
            unscored = seed.product_type == "duplexed_seed_oligo_pair"
            for edge in best.transitions:
                for fragment in edge.fragment_rows:
                    products.add(str(fragment.get("product_type", "")))
                    scored = scored or bool(fragment.get("idt_response_sha256"))
                    unscored = unscored or str(fragment.get("product_type")) == (
                        "annealed_sticky_end_primer_pair"
                    )
            evidence_tier = (
                "mixed_scored_and_primer"
                if scored and unscored else
                "all_long_fragments_live_idt_passed"
                if scored else
                "primer_only_unscored"
            )
            route_rows.append(
                {
                    "version": DNA_COMPLETE_ROUTE_VERSION,
                    "complete_route_id": full_route_id,
                    "source_database": source_database,
                    "element_id": element_id,
                    "target_copy_count": target_copies,
                    "target_id": f"{source_database}|{element_id}|copies={target_copies}",
                    "seed_copy_count": seed.copy_count,
                    "plasmid": best.plasmid,
                    "experimental_step_count": 1 + best.hurdle_steps,
                    "hurdler_step_count": best.hurdle_steps,
                    "unique_purchase_count": len(best.purchase_sha256s),
                    "total_purchase_bp": best.total_purchase_bp,
                    "pair_change_count": best.pair_change_count,
                    "maximum_idt_score": best.maximum_idt_score,
                    "transition_count": len(best.transitions),
                    "transition_route_ids_json": json.dumps(
                        [edge.route_id for edge in best.transitions]
                    ),
                    "final_target_exact": True,
                    "target_sequence_sha256": target_sha,
                    "final_sequence_sha256": target_sha,
                    "idt_evidence_tier": evidence_tier,
                }
            )
            for edge_index, edge in enumerate(best.transitions, start=1):
                transition_rows.append(
                    {
                        "version": DNA_COMPLETE_ROUTE_VERSION,
                        "complete_route_id": full_route_id,
                        "transition_index": edge_index,
                        "source_database": source_database,
                        "element_id": element_id,
                        **{
                            key: value for key, value in edge.__dict__.items()
                            if key not in {"route_row", "step_rows", "fragment_rows"}
                        },
                    }
                )
                for row in edge.step_rows:
                    step_rows.append(
                        {
                            **row,
                            "version": DNA_COMPLETE_ROUTE_VERSION,
                            "complete_route_id": full_route_id,
                            "transition_index": edge_index,
                        }
                    )
                for row in edge.fragment_rows:
                    fragment_rows.append(
                        {
                            **row,
                            "version": DNA_COMPLETE_ROUTE_VERSION,
                            "complete_route_id": full_route_id,
                            "transition_index": edge_index,
                        }
                    )
        whole_failed = whole_target.get("idt_status") == "failed"
        failure_reason = ""
        if not complete:
            target_candidate_audit = [
                row for row in candidate_route_rows
                if int(row["target_copy_count"]) == target_copies
            ]
            diagnostics = [
                value
                for (recipient, result), value in provider_diagnostics.items()
                if result == target_copies and recipient < result
            ]
            if target_candidate_audit and not any(
                bool(row["all_purchase_fragments_accepted"])
                for row in target_candidate_audit
            ):
                failure_reason = "purchase_or_idt_failure"
            elif diagnostics and all(
                int(value.get("candidate_pair_count", 0)) == 0
                for value in diagnostics
            ):
                failure_reason = "no_active_latent_pair"
            elif any(
                value.get("whole_target_idt_status")
                in {"api_failure", "api_unclassified"}
                or int(value.get("failed_idt_routes", 0)) > 0
                for value in diagnostics
            ):
                failure_reason = "purchase_or_idt_failure"
            elif diagnostics:
                failure_reason = "vector_or_digest_failure"
            else:
                failure_reason = "no_reachable_precursor_state"
        target_rows.append(
            {
                "version": DNA_COMPLETE_ROUTE_VERSION,
                "source_database": source_database,
                "element_id": element_id,
                "unit_sequence": unit,
                "unit_length_bp": len(unit),
                "target_copy_count": target_copies,
                "target_id": f"{source_database}|{element_id}|copies={target_copies}",
                "target_length_bp": len(target_sequence),
                "target_sequence_sha256": target_sha,
                "seed_copy_count": seed.copy_count,
                "seed_product_type": seed.product_type,
                "complete_route_id": full_route_id,
                "complete_route_verified": complete,
                "final_target_exact": complete,
                "hurdler_compatible": complete,
                "plasmid": best.plasmid if best else "",
                "experimental_step_count": 1 + best.hurdle_steps if best else None,
                "hurdler_step_count": best.hurdle_steps if best else None,
                "pair_change_count": best.pair_change_count if best else None,
                "unique_purchase_count": len(best.purchase_sha256s) if best else None,
                "total_purchase_bp": best.total_purchase_bp if best else None,
                "maximum_idt_score": best.maximum_idt_score if best else None,
                "whole_target_idt_status": whole_target.get("idt_status"),
                "whole_target_idt_score": whole_target.get(
                    "idt_complexity_score"
                ),
                "whole_target_idt_response_sha256": whole_target.get(
                    "idt_response_sha256", ""
                ),
                "whole_target_idt_positive_score_names_json": whole_target.get(
                    "idt_positive_score_names_json", "[]"
                ),
                "whole_target_idt_rule_details_json": whole_target.get(
                    "idt_rule_details_json", "[]"
                ),
                "idt_evidence_tier": evidence_tier,
                "fragment_rescued_by_hurdler": bool(complete and whole_failed),
                "failure_reason": failure_reason,
            }
        )
    return {
        "targets": pd.DataFrame(target_rows),
        "selected_routes": pd.DataFrame(route_rows),
        "transitions": pd.DataFrame(transition_rows),
        "steps": pd.DataFrame(step_rows),
        "fragments": pd.DataFrame(fragment_rows),
        "seeds": pd.DataFrame([seed_row]),
        "candidate_routes": pd.DataFrame(candidate_route_rows),
    }


def plan_complete_route_catalog(
    catalog: str | Path,
    reference_dir: str | Path,
    output_dir: str | Path,
    *,
    artifact_dir: str | Path | None = None,
    idt_scorer: IDTComplexityScorer | None = None,
    require_idt: bool = True,
    shard_index: int = 0,
    shard_count: int = 1,
    limit_elements: int | None = None,
) -> dict[str, pd.DataFrame]:
    frame = pd.read_parquet(catalog) if Path(catalog).suffix != ".csv" else pd.read_csv(catalog)
    real = frame.loc[frame.cohort.eq("real_element_derived")].copy()
    if real.empty:
        raise ValueError("Complete repeat routing requires real_element_derived rows")
    element_columns = ["source_database", "element_id", "unit_sequence"]
    missing = [column for column in element_columns if column not in real.columns]
    if missing:
        raise ValueError(f"Catalog lacks complete-route columns: {missing}")
    elements = (
        real.sort_values(["source_database", "element_id", "copy_count"], kind="mergesort")
        .drop_duplicates(["source_database", "element_id"])
        .reset_index(drop=True)
    )
    elements = elements.loc[
        [index % shard_count == shard_index for index in range(len(elements))]
    ]
    if limit_elements is not None:
        elements = elements.head(limit_elements)
    geometries, plasmids = load_enzyme_catalog(reference_dir, artifact_dir=artifact_dir)
    collected = {
        name: [] for name in (
            "targets", "selected_routes", "transitions", "steps", "fragments",
            "seeds", "candidate_routes"
        )
    }
    for row in elements.to_dict("records"):
        base = TargetRecord(
            target_id=str(row["element_id"]),
            sequence=validate_dna(str(row["unit_sequence"])),
            cohort="real_element_derived",
            architecture="exact_tandem",
            source_url=str(row.get("source_url", "")),
            source_accession=str(row.get("source_accession", "")),
            unit_sequence=validate_dna(str(row["unit_sequence"])),
            copy_count=1,
            notes=str(row.get("notes", "")),
            source_database=str(row["source_database"]),
            element_id=str(row["element_id"]),
        )
        result = search_complete_repeat_routes(
            base,
            geometries,
            plasmids,
            idt_scorer=idt_scorer,
            require_idt=require_idt,
        )
        for name, table in result.items():
            if not table.empty:
                collected[name].append(table)
    merged = {
        name: pd.concat(tables, ignore_index=True, sort=False) if tables else pd.DataFrame()
        for name, tables in collected.items()
    }
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    for name, table in merged.items():
        table.to_parquet(root / f"complete_route_{name}.parquet", index=False)
    write_json_atomic(
        {
            "version": DNA_COMPLETE_ROUTE_VERSION,
            "created_at": utc_now(),
            "catalog": str(Path(catalog).absolute()),
            "catalog_sha256": sha256_file(catalog),
            "shard_index": shard_index,
            "shard_count": shard_count,
            "element_rows": len(elements),
            "output_rows": {name: len(table) for name, table in merged.items()},
            "idt_required": require_idt,
            "idt_policy": IDT_SCORE_POLICY,
        },
        root / "complete_route_manifest.json",
    )
    return merged


def _read_complete_table(root: Path, name: str) -> pd.DataFrame:
    path = root / f"complete_route_{name}.parquet"
    if not path.is_file():
        return pd.DataFrame()
    return pd.read_parquet(path)


def _write_selected_purchase_fasta(
    fragments: pd.DataFrame,
    seeds: pd.DataFrame,
    destination: Path,
) -> None:
    """Write every unique selected seed/donor purchase in 5'-to-3' form."""
    lines: list[str] = []
    if not seeds.empty:
        for row in seeds.drop_duplicates("purchase_sha256").itertuples(index=False):
            lines.extend(
                [
                    f">seed_{row.purchase_sha256}_forward product={row.product_type}",
                    str(row.purchase_sequence),
                ]
            )
            if str(row.secondary_purchase_sequence):
                lines.extend(
                    [
                        f">seed_{row.purchase_sha256}_reverse product={row.product_type}",
                        str(row.secondary_purchase_sequence),
                    ]
                )
    if not fragments.empty:
        for row in fragments.drop_duplicates("purchase_sha256").itertuples(index=False):
            if row.product_type == "annealed_sticky_end_primer_pair":
                forward = getattr(
                    row, "primer_forward_5to3", row.purchase_sequence
                )
                reverse = getattr(
                    row,
                    "primer_reverse_5to3",
                    getattr(row, "secondary_purchase_sequence", ""),
                )
                lines.extend(
                    [
                        f">{row.fragment_id}_forward product={row.product_type}",
                        str(forward),
                        f">{row.fragment_id}_reverse product={row.product_type}",
                        str(reverse),
                    ]
                )
            else:
                lines.extend(
                    [
                        f">{row.fragment_id} product={row.product_type}",
                        str(row.purchase_sequence),
                    ]
                )
    destination.write_text("\n".join(lines) + ("\n" if lines else ""))


def build_element_matrix(targets: pd.DataFrame) -> pd.DataFrame:
    """Build one strict five-target row per public element."""
    required = {
        "source_database",
        "element_id",
        "unit_sequence",
        "unit_length_bp",
        "target_copy_count",
        "complete_route_verified",
        "final_target_exact",
    }
    missing = sorted(required - set(targets.columns))
    if missing:
        raise ValueError(f"Complete target table lacks columns: {missing}")
    real = targets.loc[targets.source_database.ne("Synthetic")].copy()
    if real.duplicated(["source_database", "element_id", "target_copy_count"]).any():
        raise ValueError("Duplicate element/copy-count rows in complete target table")
    rows: list[dict[str, Any]] = []
    for (source, element_id), group in real.groupby(
        ["source_database", "element_id"], sort=True
    ):
        observed = set(group.target_copy_count.astype(int))
        missing_copies = sorted(set(TARGET_COPY_COUNTS) - observed)
        if missing_copies:
            raise ValueError(
                f"{source}/{element_id} lacks copy counts: {missing_copies}"
            )
        first = group.sort_values("target_copy_count").iloc[0]
        row: dict[str, Any] = {
            "version": DNA_COMPLETE_ROUTE_VERSION,
            "source_database": source,
            "element_id": element_id,
            "unit_sequence": first.unit_sequence,
            "unit_sequence_sha256": hashlib.sha256(
                str(first.unit_sequence).encode()
            ).hexdigest(),
            "unit_length_bp": int(first.unit_length_bp),
            "seed_copy_count": (
                int(first.seed_copy_count) if pd.notna(first.seed_copy_count) else None
            ),
        }
        successes: list[int] = []
        eligibility: list[bool] = []
        for copies in TARGET_COPY_COUNTS:
            record = group.loc[group.target_copy_count.eq(copies)].iloc[0]
            passed = bool(
                record.complete_route_verified and record.final_target_exact
            )
            row[f"copy_{copies}_complete"] = passed
            row[f"copy_{copies}_failure_reason"] = str(record.failure_reason)
            row[f"copy_{copies}_route_id"] = str(record.complete_route_id)
            reviewer_eligible = bool(
                getattr(record, "reviewer_eligible", True)
            )
            row[f"copy_{copies}_reviewer_eligible"] = reviewer_eligible
            eligibility.append(reviewer_eligible)
            if passed:
                successes.append(copies)
        row["successful_target_count"] = len(successes)
        row["all_five_complete"] = len(successes) == len(TARGET_COPY_COUNTS)
        row["maximum_verified_copy_count"] = max(successes, default=0)
        row["reviewer_eligible_all_targets"] = all(eligibility)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["source_database", "element_id"], kind="mergesort"
    ).reset_index(drop=True)


def finalize_complete_route_shards(
    shard_dirs: Iterable[str | Path],
    output_dir: str | Path,
    *,
    expected_public_elements: int | None = None,
    expected_real_targets: int | None = None,
) -> dict[str, pd.DataFrame]:
    """Merge complete-route shards and write plotting/reviewer artifacts."""
    roots = [Path(path) for path in shard_dirs]
    if not roots:
        raise ValueError("At least one complete-route shard is required")
    names = (
        "targets", "selected_routes", "transitions", "steps", "fragments",
        "seeds", "candidate_routes"
    )
    merged: dict[str, pd.DataFrame] = {}
    metric_rows: list[dict[str, Any]] = []
    for root in roots:
        manifest_path = root / "complete_route_manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(manifest_path)
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("version") != DNA_COMPLETE_ROUTE_VERSION:
            raise ValueError(f"Unexpected complete-route version in {manifest_path}")
        if manifest.get("idt_required") is not True:
            raise ValueError(f"Production shard lacks mandatory live IDT mode: {root}")
        audit_path = root / "idt_audit.jsonl"
        if not audit_path.is_file():
            raise FileNotFoundError(
                f"Production shard lacks its IDT audit artifact: {audit_path}"
            )
        audit_records = 0
        if audit_path.is_file():
            with audit_path.open() as handle:
                audit_records = sum(1 for line in handle if line.strip())
        metric_rows.append(
            {
                "version": manifest.get("version"),
                "shard_index": int(manifest["shard_index"]),
                "shard_count": int(manifest["shard_count"]),
                "element_rows": int(manifest["element_rows"]),
                "target_rows": int(manifest.get("output_rows", {}).get("targets", 0)),
                "selected_route_rows": int(
                    manifest.get("output_rows", {}).get("selected_routes", 0)
                ),
                "idt_audit_path": str(audit_path.absolute()),
                "idt_audit_records": audit_records,
                "idt_audit_sha256": (
                    sha256_file(audit_path) if audit_path.is_file() else ""
                ),
                "status": "passed",
            }
        )
    metrics = pd.DataFrame(metric_rows).sort_values("shard_index")
    if metrics.shard_index.duplicated().any():
        raise ValueError("Duplicate complete-route shard indices")
    declared_shards = set(metrics.shard_count.astype(int))
    if len(declared_shards) != 1:
        raise ValueError("Shard manifests disagree on shard_count")
    expected_shards = declared_shards.pop()
    if set(metrics.shard_index.astype(int)) != set(range(expected_shards)):
        raise ValueError("Complete-route shard set is incomplete")
    for name in names:
        tables = [_read_complete_table(root, name) for root in roots]
        nonempty = [table for table in tables if not table.empty]
        merged[name] = (
            pd.concat(nonempty, ignore_index=True, sort=False)
            if nonempty else pd.DataFrame()
        )
    targets = merged["targets"]
    routes = merged["selected_routes"]
    transitions = merged["transitions"]
    steps = merged["steps"]
    fragments = merged["fragments"]
    seeds = merged["seeds"]
    if targets.target_id.duplicated().any():
        raise ValueError("Complete production has duplicate target IDs")
    expected_copy_rows = (
        targets.groupby(["source_database", "element_id"])
        .target_copy_count.agg(lambda values: set(map(int, values)))
    )
    if not expected_copy_rows.map(
        lambda observed: observed == set(TARGET_COPY_COUNTS)
    ).all():
        raise ValueError("Every public element must retain all five target copy counts")
    exact_lengths = (
        targets.unit_length_bp.astype(int)
        * targets.target_copy_count.astype(int)
    )
    if not targets.target_length_bp.astype(int).eq(exact_lengths).all():
        raise ValueError("A target length is not unit_length * target_copy_count")
    expected_target_hashes = targets.apply(
        lambda row: hashlib.sha256(
            (str(row.unit_sequence) * int(row.target_copy_count)).encode()
        ).hexdigest(),
        axis=1,
    )
    if not targets.target_sequence_sha256.eq(expected_target_hashes).all():
        raise ValueError("A target SHA does not match its exact tandem sequence")
    scored_targets = targets.loc[
        targets.whole_target_idt_status.isin(["passed", "failed"])
    ]
    if not scored_targets.empty:
        if scored_targets.whole_target_idt_response_sha256.fillna("").eq("").any():
            raise ValueError("A classified intact target lacks a live-IDT response hash")
        if scored_targets.whole_target_idt_score.isna().any():
            raise ValueError("A classified intact target lacks a numeric IDT score")
    passed = targets.loc[targets.complete_route_verified]
    if not passed.final_target_exact.all():
        raise ValueError("A complete route does not match its exact target")
    if not passed.hurdler_step_count.astype(int).ge(1).all():
        raise ValueError("A complete route must contain at least one HURDLER growth step")
    if passed.complete_route_id.eq("").any():
        raise ValueError("A passing target lacks a complete route ID")
    if not routes.empty and routes.complete_route_id.duplicated().any():
        raise ValueError("Selected complete route IDs are not unique")
    route_ids = set(routes.complete_route_id) if not routes.empty else set()
    if set(passed.complete_route_id) != route_ids:
        raise ValueError("Passing targets and selected routes are not one-to-one")
    if not routes.empty:
        if not routes.final_target_exact.all():
            raise ValueError("A selected route lacks exact final-sequence proof")
        if not routes.target_sequence_sha256.eq(routes.final_sequence_sha256).all():
            raise ValueError("A selected route final SHA differs from its target SHA")
    if not transitions.empty:
        plasmid_counts = transitions.groupby("complete_route_id").plasmid.nunique()
        if not plasmid_counts.eq(1).all():
            raise ValueError("A complete route changes plasmid during assembly")
    if not steps.empty:
        if not steps.unintended_cut_count.fillna(0).astype(int).eq(0).all():
            raise ValueError("A selected route contains an unintended selected-enzyme cut")
        if not steps.double_strand_source_verified.fillna(False).all():
            raise ValueError("A selected route lacks double-strand source proof")
    if not fragments.empty:
        long_fragments = fragments.loc[
            fragments.product_type.ne("annealed_sticky_end_primer_pair")
        ]
        if not long_fragments.empty:
            if not long_fragments.purchase_accepted.fillna(False).all():
                raise ValueError("A selected long purchase fragment was not accepted")
            if long_fragments.idt_response_sha256.fillna("").eq("").any():
                raise ValueError("A selected long fragment lacks live-IDT evidence")
            if not long_fragments.idt_score.astype(float).lt(10).all():
                raise ValueError("A selected long fragment has IDT score >=10")
        primer_fragments = fragments.loc[
            fragments.product_type.eq("annealed_sticky_end_primer_pair")
        ]
        if not primer_fragments.empty and not primer_fragments.core_length_bp.astype(int).lt(
            PRIMER_PAIR_CORE_THRESHOLD_BP
        ).all():
            raise ValueError("A >=90-bp core was incorrectly emitted as primers")
    if not seeds.empty:
        long_seeds = seeds.loc[
            seeds.product_type.ne("duplexed_seed_oligo_pair")
        ]
        if not long_seeds.empty:
            if long_seeds.idt_response_sha256.fillna("").eq("").any():
                raise ValueError("A selected long seed lacks live-IDT evidence")
            if not long_seeds.idt_score.astype(float).lt(10).all():
                raise ValueError("A selected long seed has IDT score >=10")
    invalid_rescue = targets.fragment_rescued_by_hurdler & ~targets.complete_route_verified
    if invalid_rescue.any():
        raise ValueError("An incompatible target was marked as HURDLER-rescued")
    targets["reviewer_eligible"] = ~targets.whole_target_idt_status.isin(
        ["api_failure", "api_unclassified", "scored_unclassified"]
    )
    element_matrix = build_element_matrix(targets)
    if expected_public_elements is not None and len(element_matrix) != expected_public_elements:
        raise ValueError(
            f"Expected {expected_public_elements} elements, observed {len(element_matrix)}"
        )
    if expected_real_targets is not None and len(targets) != expected_real_targets:
        raise ValueError(
            f"Expected {expected_real_targets} targets, observed {len(targets)}"
        )
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    targets.to_parquet(output / "production_target_analysis.parquet", index=False)
    targets.to_csv(output / "production_target_analysis.csv", index=False)
    element_matrix.to_parquet(output / "production_element_matrix.parquet", index=False)
    element_matrix.to_csv(output / "production_element_matrix.csv", index=False)
    routes.to_parquet(output / "production_selected_routes.parquet", index=False)
    routes.to_csv(output / "production_selected_routes.csv", index=False)
    metrics.to_parquet(output / "production_run_metrics.parquet", index=False)
    metrics[[
        "shard_index", "shard_count", "status", "element_rows",
        "target_rows", "selected_route_rows", "idt_audit_records",
    ]].to_csv(output / "production_status_matrix.csv", index=False)
    metrics[[
        "shard_index", "idt_audit_path", "idt_audit_records",
        "idt_audit_sha256",
    ]].to_csv(output / "production_idt_audit_index.csv", index=False)
    for name in ("transitions", "steps", "fragments", "seeds", "candidate_routes"):
        merged[name].to_parquet(output / f"production_{name}.parquet", index=False)
    for name in ("transitions", "steps", "fragments", "seeds"):
        merged[name].to_csv(output / f"production_{name}.csv", index=False)
    _write_selected_purchase_fasta(
        fragments, seeds, output / "production_selected_purchases.fasta"
    )
    public_targets = targets.loc[targets.source_database.ne("Synthetic")].copy()
    reviewer_targets = public_targets.loc[public_targets.reviewer_eligible]
    public_elements = element_matrix
    eligibility_by_element = (
        public_targets.groupby(["source_database", "element_id"])
        .reviewer_eligible.all()
    )
    eligible_element_keys = set(eligibility_by_element.loc[eligibility_by_element].index)
    eligible_elements = public_elements.loc[
        [
            (row.source_database, row.element_id) in eligible_element_keys
            for row in public_elements.itertuples(index=False)
        ]
    ]
    all_five = int(eligible_elements.all_five_complete.sum())
    target_passes = int(reviewer_targets.complete_route_verified.sum())
    primer_only = int(
        reviewer_targets.idt_evidence_tier.eq("primer_only_unscored").sum()
    )
    live_scored = int(
        reviewer_targets.idt_evidence_tier.isin(
            ["all_long_fragments_live_idt_passed", "mixed_scored_and_primer"]
        ).sum()
    )
    rescued = int(reviewer_targets.fragment_rescued_by_hurdler.sum())
    headline = {
        "version": DNA_COMPLETE_ROUTE_VERSION,
        "created_at": utc_now(),
        "unique_public_elements": int(len(public_elements)),
        "elements_all_five_complete": all_five,
        "real_exact_targets": int(len(public_targets)),
        "reviewer_eligible_exact_targets": int(len(reviewer_targets)),
        "reviewer_eligible_elements": int(len(eligible_elements)),
        "real_exact_targets_complete": target_passes,
        "primer_only_complete_targets": primer_only,
        "live_idt_scored_complete_targets": live_scored,
        "valid_hurdler_rescues": rescued,
        "api_unclassified_targets": int(
            public_targets.whole_target_idt_status.isin(
                ["api_failure", "api_unclassified", "scored_unclassified"]
            ).sum()
        ),
        "legacy_final_step_success_fraction": 0.5367,
        "legacy_value_is_reviewer_eligible": False,
    }
    headline["reviewer_response_text"] = (
        "We tested "
        f"{headline['unique_public_elements']} unique public DNA elements as exact "
        "2-, 4-, 8-, 16- and 32-copy arrays. "
        f"{all_five} elements admitted complete in-silico HURDLER routes from a "
        "purchasable seed at all five lengths. Across "
        f"{headline['real_exact_targets']} exact targets, "
        f"{headline['reviewer_eligible_exact_targets']} had classifiable API evidence; "
        f"among those, {target_passes} completed the full route with a "
        "sequence-identical final product."
    )
    write_json_atomic(headline, output / "production_headline_summary.json")
    (output / "reviewer_response.md").write_text(
        str(headline["reviewer_response_text"]) + "\n"
    )
    write_json_atomic(
        {
            "version": DNA_COMPLETE_ROUTE_VERSION,
            "created_at": utc_now(),
            "shard_count": expected_shards,
            "rows": {
                "targets": len(targets),
                "elements": len(element_matrix),
                "selected_routes": len(routes),
                "transitions": len(merged["transitions"]),
                "steps": len(merged["steps"]),
                "fragments": len(merged["fragments"]),
            },
        },
        output / "production_finalize_manifest.json",
    )
    return {
        **merged,
        "element_matrix": element_matrix,
        "run_metrics": metrics,
    }
