"""Deterministic genetic refinement of synonymous HURDLER constructs."""

from __future__ import annotations

import json
import hashlib
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd

from .idt import (
    IDTComplexityScorer,
    credentials_available,
    screen_gblock_sequences,
    summarize_complexity_response,
    write_cached_response,
)
from .optimization import (
    GENETIC_CODE,
    _construct_metrics,
    _repair_site_limits,
    codon_adaptation_index,
    gc_window_extrema,
    hairpin_proxy,
    load_codon_weights,
    recognition_site_count,
    repeated_nmer_count,
    reverse_complement,
    translate_dna,
)
from .progress import ProgressCallback, emit_progress


GA_SCORE_PROFILE = {
    "selected_re_site_excess": 1_000_000_000,
    "gc_window_violation": 1_000_000_000,
    "repeated_re_site_excess": 10_000,
    "repeated_14mer": 250,
    "repeated_13mer": 100,
    "repeated_8mer": 5,
    "hairpin_10mer_proxy": 25,
    "homopolymer_excess": 250,
    "terminal_repeat_proxy": 100,
    "gc_window_soft_violation": 100,
    "negative_log_cai": 50,
}

IDT_FEEDBACK_MULTIPLIER = 2.0
GA_RE_SITE_POLICY = "nonselected-re-sites-soft-score-selected-sites-hard-v2"


@dataclass
class GAPopulationState:
    """Serializable warm-start state for a continuous synonymous GA run."""

    protein_sequence: str
    elite_sequences: tuple[str, ...] = ()
    total_generations: int = 0
    rng_state: dict[str, Any] = field(default_factory=dict)

    def to_dict(self, *, include_sequences: bool = True) -> dict[str, Any]:
        payload = {
            "protein_sequence": self.protein_sequence,
            "total_generations": int(self.total_generations),
            "rng_state": self.rng_state,
            "elite_count": len(self.elite_sequences),
            "elite_sha256": [
                hashlib.sha256(sequence.encode()).hexdigest()
                for sequence in self.elite_sequences
            ],
        }
        if include_sequences:
            payload["elite_sequences"] = list(self.elite_sequences)
        return payload


def load_restriction_sites(path: str | Path) -> tuple[str, ...]:
    frame = pd.read_csv(path)
    if "site" not in frame:
        raise ValueError("Restriction-enzyme table must contain a site column")
    sites = {
        str(site).upper()
        for site in frame.site.dropna()
        if set(str(site).upper()) <= set("ACGT")
    }
    # A recognition site and its reverse complement describe the same physical
    # site because recognition_site_count already searches both orientations.
    sites = {min(site, reverse_complement(site)) for site in sites}
    sites = sorted(sites)
    return tuple(sites)


def repeated_re_site_excess(dna: str, recognition_sites: tuple[str, ...]) -> int:
    """Count every RE-site occurrence after the first occurrence per site."""
    return sum(max(0, recognition_site_count(dna, site) - 1) for site in recognition_sites)


def homopolymer_excess(dna: str, allowed_run: int = 6) -> int:
    excess = 0
    run = 0
    previous = ""
    for base in dna:
        run = run + 1 if base == previous else 1
        previous = base
        excess += int(run > allowed_run)
    return excess


def terminal_repeat_proxy(dna: str, k: int = 10, terminal_window: int = 30) -> int:
    """Count terminal k-mers that recur elsewhere in either orientation."""
    if len(dna) < k:
        return 0
    counts: dict[str, int] = {}
    for start in range(len(dna) - k + 1):
        kmer = dna[start : start + k]
        canonical = min(kmer, reverse_complement(kmer))
        counts[canonical] = counts.get(canonical, 0) + 1
    starts = list(range(0, min(terminal_window, len(dna) - k + 1)))
    starts.extend(
        range(max(0, len(dna) - terminal_window), len(dna) - k + 1)
    )
    terminal_kmers = {
        min(dna[start : start + k], reverse_complement(dna[start : start + k]))
        for start in starts
    }
    return sum(max(0, counts[kmer] - 1) for kmer in terminal_kmers)


def adjust_ga_score_profile_from_idt(
    score_profile: dict[str, float],
    idt_summary: dict[str, Any],
) -> tuple[dict[str, float], list[dict[str, Any]]]:
    """Raise weights implicated by positive scores in an IDT-rejected sum.

    ``IsViolated`` is diagnostic under ``idt-rule-score-sum-lt10-v1``.  Every
    positive rule score contributes to the rejection total and therefore
    participates in deterministic feedback.
    """
    adjusted = dict(score_profile)
    try:
        details = json.loads(str(idt_summary.get("idt_rule_details_json", "[]")))
    except (TypeError, ValueError, json.JSONDecodeError):
        details = []
    requested: dict[str, list[dict[str, Any]]] = {}
    for detail in details if isinstance(details, list) else []:
        if not isinstance(detail, dict):
            continue
        score = detail.get("score")
        if (
            not isinstance(score, (int, float))
            or isinstance(score, bool)
            or not math.isfinite(float(score))
            or float(score) <= 0
        ):
            continue
        name = str(detail.get("name", "unnamed_rule"))
        lowered = name.lower()
        keys: set[str] = set()
        if "gc" in lowered:
            keys.add("gc_window_soft_violation")
        if "hairpin" in lowered or "palindrome" in lowered:
            keys.add("hairpin_10mer_proxy")
        if "homopolymer" in lowered:
            keys.add("homopolymer_excess")
        if "terminal" in lowered:
            keys.update({"terminal_repeat_proxy", "repeated_13mer"})
        if "repeat" in lowered or "ssa" in lowered:
            keys.update({"repeated_8mer", "repeated_13mer", "repeated_14mer"})
        if "restriction site" in lowered:
            keys.add("repeated_re_site_excess")
        if not keys:
            keys.update({"repeated_13mer", "hairpin_10mer_proxy"})
        for key in keys:
            requested.setdefault(key, []).append(detail)
    adjustments: list[dict[str, Any]] = []
    for key, reasons in sorted(requested.items()):
        old = float(adjusted[key])
        new = old * IDT_FEEDBACK_MULTIPLIER
        adjusted[key] = new
        adjustments.append(
            {
                "score_component": key,
                "old_weight": old,
                "new_weight": new,
                "multiplier": IDT_FEEDBACK_MULTIPLIER,
                "idt_reasons": [str(item.get("name", "unnamed_rule")) for item in reasons],
                "idt_scores": [item.get("score") for item in reasons],
                "idt_actual_values": [item.get("actual_value") for item in reasons],
            }
        )
    return adjusted, adjustments


def ga_sequence_metrics(
    dna: str,
    codon_weights: dict[str, float],
    recognition_sites: tuple[str, ...],
    selected_site_limits: dict[str, int],
    score_profile: dict[str, float] | None = None,
) -> dict[str, Any]:
    profile = dict(GA_SCORE_PROFILE if score_profile is None else score_profile)
    gc_min, gc_max = gc_window_extrema(dna, 50)
    selected_excess = sum(
        max(0, recognition_site_count(dna, site) - limit)
        for site, limit in selected_site_limits.items()
        if site
    )
    re_repeat_excess = repeated_re_site_excess(dna, recognition_sites)
    repeat8 = repeated_nmer_count(dna, 8)
    repeat13 = repeated_nmer_count(dna, 13)
    repeat14 = repeated_nmer_count(dna, 14)
    hairpins = hairpin_proxy(dna, 10)
    homopolymers = homopolymer_excess(dna)
    terminal_repeats = terminal_repeat_proxy(dna)
    cai = codon_adaptation_index(dna, codon_weights)
    gc_violation = max(0.0, 0.25 - gc_min) + max(0.0, gc_max - 0.75)
    gc_soft_violation = max(0.0, 0.35 - gc_min) + max(0.0, gc_max - 0.65)
    score = (
        profile["selected_re_site_excess"] * selected_excess
        + profile["gc_window_violation"] * gc_violation
        + profile["gc_window_soft_violation"] * gc_soft_violation
        + profile["repeated_re_site_excess"] * re_repeat_excess
        + profile["repeated_14mer"] * repeat14
        + profile["repeated_13mer"] * repeat13
        + profile["repeated_8mer"] * repeat8
        + profile["hairpin_10mer_proxy"] * hairpins
        + profile["homopolymer_excess"] * homopolymers
        + profile["terminal_repeat_proxy"] * terminal_repeats
        - profile["negative_log_cai"] * math.log(max(cai, 1e-12))
    )
    return {
        "ga_score": float(score),
        "selected_re_site_excess": int(selected_excess),
        "selected_pair_re_site_excess": int(selected_excess),
        "repeated_re_site_excess": int(re_repeat_excess),
        "repeat_8mer_count": int(repeat8),
        "repeat_13mer_count": int(repeat13),
        "repeat_14mer_count": int(repeat14),
        "hairpin_10mer_proxy": int(hairpins),
        "homopolymer_excess": int(homopolymers),
        "terminal_repeat_proxy": int(terminal_repeats),
        "gc_50bp_min": float(gc_min),
        "gc_50bp_max": float(gc_max),
        "ga_gc_bounds_passed": bool(gc_violation == 0),
        "gc_50bp_soft_violation": float(gc_soft_violation),
        "cai_kazusa_83333": float(cai),
        # Live IDT, not a duplicated local GC cutoff, is the orderability
        # gate. Locally, only the frozen selected Site-I/Site-II counts are
        # hard constraints.
        "ga_local_constraints_passed": bool(selected_excess == 0),
    }


def genetic_refine_dna(
    dna: str,
    *,
    locked_positions: set[int],
    selected_site_limits: dict[str, int],
    recognition_sites: tuple[str, ...],
    codon_weights: dict[str, float],
    seed: int,
    population_size: int = 16,
    generations: int = 20,
    mutation_rate: float = 0.015,
    crossover_rate: float = 1.0,
    elite_fraction: float = 0.125,
    score_profile: dict[str, float] | None = None,
    population_state: GAPopulationState | None = None,
    elite_seed_count: int = 6,
    capture_population_state: bool = False,
    progress_callback: ProgressCallback | None = None,
    progress_context: dict[str, Any] | None = None,
) -> tuple[str, dict[str, Any]]:
    """Refine synonymous codons with the RE-repeat term in every fitness call."""
    if translate_dna(dna) == "":
        raise ValueError("DNA cannot be empty")
    protein = translate_dna(dna)
    codons = [dna[index : index + 3] for index in range(0, len(dna), 3)]
    repaired = _repair_site_limits(
        dna,
        protein,
        locked_positions,
        selected_site_limits,
        codon_weights,
    )
    rng = np.random.default_rng(seed)
    if population_state is not None:
        if population_state.protein_sequence != protein:
            raise ValueError("Warm-start GA state encodes a different protein")
        if population_state.rng_state:
            rng.bit_generator.state = population_state.rng_state
    cache: dict[str, dict[str, Any]] = {}
    profile = dict(GA_SCORE_PROFILE if score_profile is None else score_profile)
    context = dict(progress_context or {})
    started = time.monotonic()
    emit_progress(
        progress_callback,
        stage="ga",
        status="started",
        generations=int(generations),
        generation=0,
        elapsed_seconds=0.0,
        **context,
    )

    def metrics(sequence: str) -> dict[str, Any]:
        if sequence not in cache:
            cache[sequence] = ga_sequence_metrics(
                sequence,
                codon_weights,
                recognition_sites,
                selected_site_limits,
                score_profile=profile,
            )
        return cache[sequence]

    def mutate(sequence: str, rate: float) -> str:
        child = [sequence[index : index + 3] for index in range(0, len(sequence), 3)]
        for position, aa in enumerate(protein):
            if position in locked_positions or rng.random() >= rate:
                continue
            choices = GENETIC_CODE[aa]
            weights = np.array([codon_weights.get(codon, 1.0) for codon in choices], dtype=float)
            weights /= weights.sum()
            child[position] = str(rng.choice(choices, p=weights))
        return "".join(child)

    elite_seed_count = int(elite_seed_count)
    if elite_seed_count < 1:
        raise ValueError("elite_seed_count must be at least one")
    warm_sequences: list[str] = []
    if population_state is not None:
        for sequence in population_state.elite_sequences:
            if len(sequence) != len(repaired) or translate_dna(sequence) != protein:
                continue
            if any(
                sequence[position * 3 : position * 3 + 3]
                != dna[position * 3 : position * 3 + 3]
                for position in locked_positions
            ):
                continue
            warm_sequences.append(sequence)
    population = list(dict.fromkeys([*warm_sequences, repaired]))
    seed_pool = tuple(population)
    while len(population) < population_size:
        parent = seed_pool[len(population) % len(seed_pool)]
        population.append(mutate(parent, max(mutation_rate, 0.03)))
    initial_metrics = metrics(dna)
    best = repaired
    for _generation in range(generations):
        population = sorted(
            set(population),
            key=lambda sequence: (metrics(sequence)["ga_score"], sequence),
        )
        if metrics(population[0])["ga_score"] < metrics(best)["ga_score"]:
            best = population[0]
        parents = population[: max(2, min(elite_seed_count, len(population)))]
        elite_count = max(1, min(len(parents), round(population_size * float(elite_fraction))))
        next_population = parents[:elite_count]
        while len(next_population) < population_size:
            first, second = rng.choice(parents, size=2, replace=True)
            first_codons = [first[index : index + 3] for index in range(0, len(first), 3)]
            second_codons = [second[index : index + 3] for index in range(0, len(second), 3)]
            if rng.random() < float(crossover_rate):
                mask = rng.random(len(codons)) < 0.5
                child = "".join(
                    first_codons[position] if mask[position] else second_codons[position]
                    for position in range(len(codons))
                )
            else:
                child = str(first)
            next_population.append(mutate(child, mutation_rate))
        population = next_population
        current = metrics(best)
        emit_progress(
            progress_callback,
            stage="ga",
            status="running",
            generations=int(generations),
            generation=int(_generation + 1),
            ga_score=float(current["ga_score"]),
            selected_pair_re_site_excess=int(
                current.get("selected_pair_re_site_excess", current["selected_re_site_excess"])
            ),
            elapsed_seconds=time.monotonic() - started,
            **context,
        )
    ranked_population = sorted(
        set([best, *population]),
        key=lambda sequence: (metrics(sequence)["ga_score"], sequence),
    )
    best = ranked_population[0]
    best = _repair_site_limits(
        best,
        protein,
        locked_positions,
        selected_site_limits,
        codon_weights,
    )
    final_metrics = metrics(best)
    if translate_dna(best) != protein:
        raise AssertionError("Genetic refinement changed translation")
    for position in locked_positions:
        if best[position * 3 : position * 3 + 3] != dna[position * 3 : position * 3 + 3]:
            raise AssertionError("Genetic refinement changed a locked HURDLER codon")
    final_metrics.update(
        {
            "ga_initial_score": float(initial_metrics["ga_score"]),
            "ga_initial_repeated_re_site_excess": int(initial_metrics["repeated_re_site_excess"]),
            "ga_repeated_re_site_excess_removed": int(
                initial_metrics["repeated_re_site_excess"] - final_metrics["repeated_re_site_excess"]
            ),
            "ga_improved": bool(final_metrics["ga_score"] < initial_metrics["ga_score"]),
            "ga_population_size": population_size,
            "ga_generations": generations,
            "ga_seed": seed,
            "ga_mutation_rate": float(mutation_rate),
            "ga_crossover_rate": float(crossover_rate),
            "ga_elite_fraction": float(elite_fraction),
            "ga_score_profile_json": json.dumps(profile, sort_keys=True),
        }
    )
    if capture_population_state:
        state_candidates = set([best, *ranked_population])
        fill_attempts = 0
        fill_limit = max(50, elite_seed_count * 50)
        fill_parents = tuple(state_candidates) or (repaired,)
        while len(state_candidates) < elite_seed_count and fill_attempts < fill_limit:
            parent = fill_parents[fill_attempts % len(fill_parents)]
            state_candidates.add(mutate(parent, max(float(mutation_rate), 0.10)))
            fill_attempts += 1
        final_ranked = sorted(
            state_candidates,
            key=lambda sequence: (metrics(sequence)["ga_score"], sequence),
        )[:elite_seed_count]
        elite_rows = []
        for rank, sequence in enumerate(final_ranked, start=1):
            row_metrics = metrics(sequence)
            elite_rows.append(
                {
                    "rank": rank,
                    "dna_sequence": sequence,
                    "dna_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                    "ga_score": float(row_metrics["ga_score"]),
                    "ga_local_constraints_passed": bool(
                        row_metrics["ga_local_constraints_passed"]
                    ),
                    "selected_pair_re_site_excess": int(
                        row_metrics["selected_pair_re_site_excess"]
                    ),
                    "repeated_re_site_excess": int(
                        row_metrics["repeated_re_site_excess"]
                    ),
                }
            )
        next_state = GAPopulationState(
            protein_sequence=protein,
            elite_sequences=tuple(row["dna_sequence"] for row in elite_rows),
            total_generations=(
                int(population_state.total_generations) if population_state else 0
            )
            + int(generations),
            rng_state=json.loads(json.dumps(rng.bit_generator.state)),
        )
        final_metrics["ga_elite_candidates"] = elite_rows
        final_metrics["ga_population_state"] = next_state
    emit_progress(
        progress_callback,
        stage="ga",
        status="completed",
        generations=int(generations),
        generation=int(generations),
        ga_score=float(final_metrics["ga_score"]),
        selected_pair_re_site_excess=int(final_metrics["selected_pair_re_site_excess"]),
        elapsed_seconds=time.monotonic() - started,
        **context,
    )
    return best, final_metrics


def adaptive_copy_search(
    minimum_copies: int,
    maximum_copies: int,
    *,
    short_generations: int,
    generation_schedule: tuple[int, ...],
    evaluate: Callable[[int, int], dict[str, Any]],
    progress_callback: ProgressCallback | None = None,
    progress_context: dict[str, Any] | None = None,
) -> tuple[int, dict[str, Any] | None, list[dict[str, Any]], str]:
    """Find the longest passing repeat count using the requested two-stage route.

    Stage one uses only ``short_generations`` in a binary search.  Stage two
    starts one copy above that result, advances exactly one module at a time,
    and raises the generation budget through a deterministic schedule ending
    at 100.  The first copy count that still fails at 100 terminates the search.
    """
    if minimum_copies < 1 or maximum_copies < minimum_copies:
        return 0, None, [], "no_feasible_copy_range"
    if short_generations < 1 or short_generations > 100:
        raise ValueError("short_generations must be between 1 and 100")
    schedule = tuple(
        sorted({short_generations, *(int(value) for value in generation_schedule), 100})
    )
    if schedule[0] < 1 or schedule[-1] > 100:
        raise ValueError("generation_schedule values must be between 1 and 100")

    cache: dict[tuple[int, int], dict[str, Any]] = {}
    trace: list[dict[str, Any]] = []
    context = dict(progress_context or {})
    started = time.monotonic()

    def run(copies: int, generations: int, phase: str) -> dict[str, Any]:
        key = (copies, generations)
        cached = key in cache
        emit_progress(
            progress_callback,
            stage="copy_search",
            status="attempt_started",
            copies=int(copies),
            phase=phase,
            generations=int(generations),
            generation=0,
            elapsed_seconds=time.monotonic() - started,
            details={"cached": cached},
            **context,
        )
        if not cached:
            cache[key] = evaluate(copies, generations)
        result = cache[key]
        trace.append(
            {
                "phase": phase,
                "copies": copies,
                "generations": generations,
                "passed": bool(result.get("passed", False)),
                "terminal": bool(result.get("terminal", False)),
                "cached": cached,
                "ga_score": result.get("ga_score"),
                "repeated_re_site_excess": result.get("repeated_re_site_excess"),
                "selected_re_site_excess": result.get("selected_re_site_excess"),
                "selected_pair_re_site_excess": result.get(
                    "selected_pair_re_site_excess",
                    result.get("selected_re_site_excess"),
                ),
                "ga_local_constraints_passed": result.get(
                    "ga_local_constraints_passed"
                ),
                "idt_request_attempted": result.get(
                    "idt_request_attempted", False
                ),
                "idt_api_called": result.get("idt_api_called", False),
                "idt_cache_hit": result.get("idt_cache_hit", False),
                "idt_status": result.get("idt_status", "not_scored_local_failure"),
                "idt_explicit_pass": result.get("idt_explicit_pass"),
                "idt_violation_count": result.get("idt_violation_count"),
                "idt_complexity_score": result.get("idt_complexity_score"),
                "idt_score_complete": result.get("idt_score_complete"),
                "idt_score_policy": result.get("idt_score_policy"),
                "idt_violation_names_json": result.get(
                    "idt_violation_names_json", "[]"
                ),
                "idt_rule_scores_json": result.get("idt_rule_scores_json", "{}"),
                "idt_rule_details_json": result.get("idt_rule_details_json", "[]"),
                "idt_scored_sequence_sha256": result.get(
                    "idt_scored_sequence_sha256", ""
                ),
                "idt_response_sha256": result.get("idt_response_sha256", ""),
                "idt_positive_score_names_json": result.get(
                    "idt_positive_score_names_json", "[]"
                ),
                "idt_invalid_score_names_json": result.get(
                    "idt_invalid_score_names_json", "[]"
                ),
                "candidate_dna_length_bp": len(
                    str(result.get("dna_sequence", ""))
                ),
                "candidate_dna_sha256": hashlib.sha256(
                    str(result.get("dna_sequence", "")).encode()
                ).hexdigest()
                if result.get("dna_sequence")
                else "",
                "ga_score_profile_before_idt_json": result.get(
                    "ga_score_profile_before_idt_json", "{}"
                ),
                "ga_score_profile_after_idt_json": result.get(
                    "ga_score_profile_after_idt_json", "{}"
                ),
                "idt_feedback_adjustments_json": result.get(
                    "idt_feedback_adjustments_json", "[]"
                ),
                "error": result.get("error", ""),
                "idt_error": result.get("idt_error", ""),
            }
        )
        emit_progress(
            progress_callback,
            stage="copy_search",
            status="attempt_completed",
            copies=int(copies),
            phase=phase,
            generations=int(generations),
            generation=int(generations),
            ga_score=(
                float(result["ga_score"])
                if isinstance(result.get("ga_score"), (int, float))
                else None
            ),
            selected_pair_re_site_excess=result.get(
                "selected_pair_re_site_excess", result.get("selected_re_site_excess")
            ),
            elapsed_seconds=time.monotonic() - started,
            details={
                "passed": bool(result.get("passed", False)),
                "idt_status": result.get("idt_status", ""),
                "cached": cached,
            },
            **context,
        )
        return result

    best_copies = 0
    best_result: dict[str, Any] | None = None
    low = minimum_copies
    high = maximum_copies
    while low <= high:
        middle = (low + high) // 2
        result = run(middle, short_generations, "binary_short")
        if result.get("passed"):
            best_copies = middle
            best_result = result
            low = middle + 1
        else:
            high = middle - 1

    next_copies = best_copies + 1 if best_copies else minimum_copies
    for copies in range(next_copies, maximum_copies + 1):
        passing_result = None
        for generations in schedule:
            result = run(copies, generations, "linear_escalation")
            if result.get("passed"):
                passing_result = result
                break
        if passing_result is None:
            # Even deterministic construction failures are materialized at
            # every requested budget through 100 generations. This gives all
            # maxima the same machine-checkable next-copy proof.
            reason = f"copy_{copies}_failed_at_100"
            return best_copies, best_result, trace, reason
        best_copies = copies
        best_result = passing_result
    return best_copies, best_result, trace, "reached_local_upper_bound"


def _locked_positions_and_limits(payload: dict[str, Any]) -> tuple[set[int], dict[str, int]]:
    locked_positions = {
        int(payload["site_i_position"]) + offset for offset in range(3)
    } | {
        int(payload["site_ii_position"]) + offset for offset in range(3)
    }
    selected_limits = {
        str(payload.get("site_i_recognition_site", "")): 1,
        str(payload.get("site_ii_recognition_site", "")): 0,
    }
    selected_limits = {
        site: limit
        for site, limit in selected_limits.items()
        if site and site != "nan"
    }
    return locked_positions, selected_limits


def _adaptive_seed(base_seed: int, payload: dict[str, Any], copies: int) -> int:
    token = (
        f"{base_seed}|{payload.get('module_id')}|"
        f"{payload.get('fragment_limit_bp')}|{copies}"
    ).encode()
    return int.from_bytes(hashlib.sha256(token).digest()[:4], "big")


def _score_idt_candidate(
    scorer: IDTComplexityScorer,
    request_name: str,
    dna: str,
) -> dict[str, Any]:
    """Score one exact GA candidate without losing the surrounding shard.

    Authentication, rate-limit, transport, and response-schema failures are
    scientific non-passes. They remain retryable at the next generation
    budget and are distinguished from a genuine IDT score rejection.
    """
    try:
        summary = scorer.score(request_name, dna)
        summary["idt_request_attempted"] = True
        return summary
    except Exception as exc:  # network/API failure is auditable, never a pass
        return {
            "idt_request_attempted": True,
            "idt_api_called": False,
            "idt_cache_hit": False,
            "idt_status": "api_failure",
            "idt_explicit_pass": None,
            "idt_violation_count": None,
            "idt_violation_names_json": "[]",
            "idt_positive_score_names_json": "[]",
            "idt_invalid_score_names_json": "[]",
            "idt_rule_scores_json": "{}",
            "idt_rule_details_json": "[]",
            "idt_complexity_score": None,
            "idt_score_complete": False,
            "idt_score_policy": None,
            "idt_scored_sequence_length_bp": len(dna),
            "idt_scored_sequence_sha256": hashlib.sha256(dna.encode()).hexdigest(),
            "idt_response_sha256": "",
            "idt_error": f"{type(exc).__name__}: {str(exc)[:500]}",
        }


def _adaptive_refine_payload(
    payload: dict[str, Any],
    *,
    codon_weights: dict[str, float],
    recognition_sites: tuple[str, ...],
    seed: int,
    population_size: int,
    short_generations: int,
    generation_schedule: tuple[int, ...],
    idt_scorer: IDTComplexityScorer | None = None,
    require_idt_orderable: bool = False,
) -> dict[str, Any]:
    original_dna = payload.get("dna_sequence")
    pre_adaptive_verified = int(payload.get("verified_max_copies", 0) or 0)
    # The adaptive search must be allowed to recover copies that the legacy
    # deterministic optimizer dropped. The fragment-length ceiling, rather
    # than that earlier result, is therefore the scientific upper bound.
    local_upper = int(payload.get("mathematical_max_copies", 0) or 0)
    unit = str(payload["unit_sequence"])
    payload["pre_adaptive_verified_max_copies"] = pre_adaptive_verified
    payload["adaptive_search_upper_bound_copies"] = local_upper
    payload["dna_sequence_pre_adaptive"] = original_dna
    payload["ga_re_site_policy"] = GA_RE_SITE_POLICY
    geometric_minimum = max(
        1,
        math.ceil(
            (max(int(payload["site_i_position"]), int(payload["site_ii_position"])) + 3)
            / len(unit)
        ),
    )
    known_orderable_copies = int(payload.get("known_orderable_copies", 0) or 0)
    known_orderable_dna = payload.get("known_orderable_dna_sequence")
    known_orderable_pre_ga = payload.get("known_orderable_dna_pre_ga")
    if known_orderable_copies:
        if not isinstance(known_orderable_dna, str) or not known_orderable_dna:
            raise ValueError("known_orderable_copies requires known_orderable_dna_sequence")
        if not geometric_minimum <= known_orderable_copies <= local_upper:
            raise ValueError("known orderable copies must be within the adaptive search range")
        if translate_dna(known_orderable_dna) != unit * known_orderable_copies:
            raise ValueError("known orderable DNA does not encode the requested repeat count")
    # A HURDLER repeat construct is only scientifically verified from two
    # complete repeats onward.  Single-copy successes may be useful debug
    # probes but must never be reported as a maximum repeat count.
    minimum = max(2, geometric_minimum, known_orderable_copies)
    payload["adaptive_search_minimum_copies"] = minimum
    payload["adaptive_known_orderable_copies"] = known_orderable_copies
    payload["adaptive_orderability_gate"] = bool(require_idt_orderable)
    locked_positions, selected_limits = _locked_positions_and_limits(payload)
    build_cache: dict[int, str | Exception] = {}
    score_profiles: dict[int, dict[str, float]] = {}

    def evaluate(copies: int, generations: int) -> dict[str, Any]:
        if copies == known_orderable_copies and isinstance(known_orderable_dna, str):
            profile = score_profiles.setdefault(copies, dict(GA_SCORE_PROFILE))
            profile_before = dict(profile)
            metrics = ga_sequence_metrics(
                known_orderable_dna,
                codon_weights,
                recognition_sites,
                selected_limits,
                score_profile=profile_before,
            )
            metrics.update(
                ga_initial_score=float(metrics["ga_score"]),
                ga_initial_repeated_re_site_excess=int(
                    metrics["repeated_re_site_excess"]
                ),
                ga_repeated_re_site_excess_removed=0,
                ga_improved=False,
                ga_population_size=population_size,
                ga_generations=generations,
                ga_seed=_adaptive_seed(seed, payload, copies),
                ga_score_profile_json=json.dumps(profile_before, sort_keys=True),
            )
            # Repeated non-selected RE sites and GC are soft GA objectives.
            # Only selected HURDLER-site counts are local hard constraints;
            # IDT decides whether the optimized sequence is orderable.
            local_passed = bool(metrics["ga_local_constraints_passed"])
            idt_summary: dict[str, Any] = {
                "idt_request_attempted": False,
                "idt_api_called": False,
                "idt_cache_hit": False,
                "idt_status": "not_scored_local_failure",
                "idt_explicit_pass": None,
                "idt_violation_count": None,
                "idt_violation_names_json": "[]",
                "idt_rule_scores_json": "{}",
                "idt_rule_details_json": "[]",
                "idt_scored_sequence_sha256": "",
            }
            adjustments: list[dict[str, Any]] = []
            if local_passed and idt_scorer is not None:
                request_name = (
                    f"{payload['module_id']}|cap={payload['fragment_limit_bp']}|"
                    f"copies={copies}|known-lower-bound"
                )
                idt_summary = _score_idt_candidate(
                    idt_scorer, request_name, known_orderable_dna
                )
                passed = idt_summary.get("idt_explicit_pass") is True
                if not passed:
                    profile, adjustments = adjust_ga_score_profile_from_idt(
                        profile_before, idt_summary
                    )
                    score_profiles[copies] = profile
            elif local_passed and require_idt_orderable:
                raise RuntimeError(
                    "Adaptive orderability search requires an IDT complexity scorer"
                )
            else:
                passed = local_passed and not require_idt_orderable
            return {
                "passed": bool(passed),
                "terminal": False,
                "dna_sequence": known_orderable_dna,
                "dna_sequence_pre_ga": (
                    known_orderable_pre_ga
                    if isinstance(known_orderable_pre_ga, str)
                    and known_orderable_pre_ga
                    else known_orderable_dna
                ),
                "metrics": metrics,
                "ga_score": metrics["ga_score"],
                "ga_local_constraints_passed": local_passed,
                "repeated_re_site_excess": metrics["repeated_re_site_excess"],
                "selected_re_site_excess": metrics["selected_re_site_excess"],
                "selected_pair_re_site_excess": metrics.get(
                    "selected_pair_re_site_excess",
                    metrics["selected_re_site_excess"],
                ),
                "ga_score_profile_before_idt_json": json.dumps(
                    profile_before, sort_keys=True
                ),
                "ga_score_profile_after_idt_json": json.dumps(
                    score_profiles[copies], sort_keys=True
                ),
                "idt_feedback_adjustments_json": json.dumps(
                    adjustments, sort_keys=True
                ),
                "idt_summary": idt_summary,
                **idt_summary,
            }
        if copies not in build_cache:
            try:
                build_cache[copies] = str(
                    _construct_metrics(
                        unit,
                        copies,
                        payload,
                        codon_weights,
                        validate_hard_constraints=False,
                    )["dna_sequence"]
                )
            except Exception as exc:  # auditable hard construction failure
                build_cache[copies] = exc
        built = build_cache[copies]
        if isinstance(built, Exception):
            return {
                "passed": False,
                "terminal": True,
                "error": f"{type(built).__name__}: {built}",
            }
        profile = score_profiles.setdefault(copies, dict(GA_SCORE_PROFILE))
        profile_before = dict(profile)
        refined, metrics = genetic_refine_dna(
            built,
            locked_positions=locked_positions,
            selected_site_limits=selected_limits,
            recognition_sites=recognition_sites,
            codon_weights=codon_weights,
            seed=_adaptive_seed(seed, payload, copies),
            population_size=population_size,
            generations=generations,
            score_profile=profile_before,
        )
        # Keep GC and repeated non-selected RE sites in the fitness function,
        # not in the local hard gate. The selected site limits remain hard
        # through ``ga_local_constraints_passed``; live IDT decides whether
        # the completed exact DNA is orderable.
        local_passed = bool(metrics["ga_local_constraints_passed"])
        idt_summary: dict[str, Any] = {
            "idt_request_attempted": False,
            "idt_api_called": False,
            "idt_cache_hit": False,
            "idt_status": "not_scored_local_failure",
            "idt_explicit_pass": None,
            "idt_violation_count": None,
            "idt_violation_names_json": "[]",
            "idt_rule_scores_json": "{}",
            "idt_rule_details_json": "[]",
            "idt_scored_sequence_sha256": "",
        }
        adjustments: list[dict[str, Any]] = []
        if local_passed and idt_scorer is not None:
            request_name = (
                f"{payload['module_id']}|cap={payload['fragment_limit_bp']}|"
                f"copies={copies}|generations={generations}"
            )
            idt_summary = _score_idt_candidate(idt_scorer, request_name, refined)
            passed = idt_summary.get("idt_explicit_pass") is True
            if not passed:
                profile, adjustments = adjust_ga_score_profile_from_idt(
                    profile_before, idt_summary
                )
                score_profiles[copies] = profile
        elif local_passed and require_idt_orderable:
            raise RuntimeError(
                "Adaptive orderability search requires an IDT complexity scorer"
            )
        else:
            passed = local_passed and not require_idt_orderable
        profile_after = score_profiles[copies]
        return {
            "passed": bool(passed),
            "terminal": False,
            "dna_sequence": refined,
            "dna_sequence_pre_ga": built,
            "metrics": metrics,
            "ga_score": metrics["ga_score"],
            "ga_local_constraints_passed": local_passed,
            "repeated_re_site_excess": metrics["repeated_re_site_excess"],
            "selected_re_site_excess": metrics["selected_re_site_excess"],
            "selected_pair_re_site_excess": metrics.get(
                "selected_pair_re_site_excess",
                metrics["selected_re_site_excess"],
            ),
            "ga_score_profile_before_idt_json": json.dumps(
                profile_before, sort_keys=True
            ),
            "ga_score_profile_after_idt_json": json.dumps(
                profile_after, sort_keys=True
            ),
            "idt_feedback_adjustments_json": json.dumps(
                adjustments, sort_keys=True
            ),
            "idt_summary": idt_summary,
            **idt_summary,
        }

    best_copies, best, trace, stop_reason = adaptive_copy_search(
        minimum,
        local_upper,
        short_generations=short_generations,
        generation_schedule=generation_schedule,
        evaluate=evaluate,
    )
    payload["adaptive_search_trace_json"] = json.dumps(trace, sort_keys=True)
    payload["adaptive_search_evaluations"] = len(trace)
    payload["adaptive_short_generations"] = short_generations
    payload["adaptive_generation_schedule"] = ",".join(
        str(value) for value in generation_schedule
    )
    payload["adaptive_stop_reason"] = stop_reason
    reported_maximum: int | None = best_copies if best_copies >= 2 else None
    payload["adaptive_verified_max_copies"] = reported_maximum
    payload["verified_max_copies"] = reported_maximum
    payload["adaptive_maximum_proof_status"] = (
        "capacity_limit_reached"
        if stop_reason == "reached_local_upper_bound" and reported_maximum is not None
        else "next_copy_failed_at_100"
        if stop_reason.endswith("_failed_at_100") and reported_maximum is not None
        else "no_accepted_repeat_construct"
    )
    idt_trace = [
        item
        for item in trace
        if item.get("idt_request_attempted") or item.get("idt_api_called")
    ]
    payload["adaptive_idt_scored_evaluations"] = len(idt_trace)
    payload["adaptive_idt_cache_hits"] = sum(
        bool(item.get("idt_cache_hit")) for item in idt_trace
    )
    payload["adaptive_idt_status_counts_json"] = json.dumps(
        {
            str(status): sum(item.get("idt_status") == status for item in idt_trace)
            for status in sorted({str(item.get("idt_status")) for item in idt_trace})
        },
        sort_keys=True,
    )
    if best is None:
        payload["dna_sequence_pre_ga"] = original_dna
        payload["dna_sequence"] = None
        payload["dna_length"] = 0
        payload["ga_status"] = "no_accepted_repeat_construct"
        payload["ga_adaptive_constraints_passed"] = False
        payload["adaptive_orderable_passed"] = False
        payload["optimization_status"] = "no_accepted_repeat_construct"
        payload["failure_reason"] = stop_reason
        if idt_trace:
            last_idt = next(
                result
                for result in reversed(trace)
                if result.get("idt_request_attempted")
                or result.get("idt_api_called")
            )
            for key in (
                "idt_explicit_pass",
                "idt_violation_count",
                "idt_violation_names_json",
                "idt_rule_scores_json",
                "idt_rule_details_json",
                "idt_positive_score_names_json",
                "idt_complexity_score",
                "idt_score_complete",
                "idt_score_policy",
                "idt_scored_sequence_sha256",
                "idt_response_sha256",
                "idt_error",
            ):
                payload[key] = last_idt.get(key)
            payload["idt_status"] = "failed_no_orderable_candidate"
            payload["idt_api_called"] = any(
                bool(item.get("idt_api_called")) for item in idt_trace
            )
        else:
            payload["idt_status"] = "not_scored_no_local_candidate"
            payload["idt_api_called"] = False
        weight_rows = [
            item
            for item in trace
            if item.get("ga_score_profile_after_idt_json") not in (None, "{}")
        ]
        payload["final_ga_weights_json"] = (
            weight_rows[-1]["ga_score_profile_after_idt_json"]
            if weight_rows
            else json.dumps(GA_SCORE_PROFILE, sort_keys=True)
        )
        return payload
    metrics = dict(best["metrics"])
    metrics["ga_adaptive_constraints_passed"] = True
    payload.update(metrics)
    payload["dna_sequence_pre_ga"] = best["dna_sequence_pre_ga"]
    payload["dna_sequence"] = best["dna_sequence"]
    payload["dna_length"] = len(str(best["dna_sequence"]))
    payload["ga_status"] = "passed"
    payload.update(best["idt_summary"])
    payload["final_ga_weights_json"] = best.get(
        "ga_score_profile_after_idt_json", metrics.get("ga_score_profile_json", "{}")
    )
    payload["adaptive_orderable_passed"] = bool(
        best.get("idt_explicit_pass") is True
        if require_idt_orderable
        else True
    )
    payload["optimization_status"] = (
        "passed_adaptive_maximum"
        if best_copies == local_upper
        else "passed_adaptive_reduced"
    )
    payload["failure_reason"] = ""
    return payload


def refine_construct_table(
    constructs_path: str | Path,
    codon_usage_path: str | Path,
    restriction_sites_path: str | Path,
    output_dir: str | Path,
    *,
    shard_index: int = 0,
    shard_count: int = 1,
    population_size: int = 16,
    generations: int = 20,
    seed: int = 42,
    use_idt: bool = False,
    adaptive_copy_search_enabled: bool = False,
    short_generations: int = 10,
    generation_schedule: tuple[int, ...] = (10, 20, 40, 60, 80, 100),
) -> pd.DataFrame:
    if adaptive_copy_search_enabled and not use_idt:
        raise ValueError(
            "Adaptive maximum-copy search requires --use-idt so orderability "
            "participates in every length decision"
        )
    if adaptive_copy_search_enabled and not credentials_available():
        raise RuntimeError(
            "Adaptive maximum-copy search requires available IDT credentials"
        )
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    idt_scorer = (
        IDTComplexityScorer(destination / "idt_optimization_responses.jsonl")
        if adaptive_copy_search_enabled and use_idt
        else None
    )
    frame = pd.read_parquet(constructs_path)
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("shard_index must satisfy 0 <= shard_index < shard_count")
    frame = frame.iloc[shard_index::shard_count].copy()
    if adaptive_copy_search_enabled and {
        "module_id",
        "fragment_limit_bp",
    }.issubset(frame.columns):
        frame = frame.sort_values(
            ["module_id", "fragment_limit_bp"], kind="mergesort"
        ).reset_index(drop=True)
    codon_weights = load_codon_weights(codon_usage_path)
    recognition_sites = load_restriction_sites(restriction_sites_path)
    rows: list[dict[str, Any]] = []
    idt_inputs: list[dict[str, str]] = []
    idt_row_positions: list[int] = []
    known_orderable_by_module: dict[str, dict[str, Any]] = {}
    for row_number, row in enumerate(frame.to_dict(orient="records")):
        payload = dict(row)
        payload["ga_re_site_policy"] = GA_RE_SITE_POLICY
        module_key = str(payload.get("module_id", ""))
        known = known_orderable_by_module.get(module_key)
        if adaptive_copy_search_enabled and known is not None:
            payload.update(known)
        dna = payload.get("dna_sequence")
        if not isinstance(dna, str) or not dna:
            payload.update(
                ga_status="not_applicable_no_construct",
                idt_status="not_applicable_no_construct",
                final_passed=False,
            )
            if adaptive_copy_search_enabled:
                payload.update(
                    pre_adaptive_verified_max_copies=0,
                    adaptive_search_upper_bound_copies=int(
                        payload.get("mathematical_max_copies", 0) or 0
                    ),
                    adaptive_search_minimum_copies=2,
                    adaptive_verified_max_copies=None,
                    verified_max_copies=None,
                    adaptive_search_trace_json="[]",
                    adaptive_search_evaluations=0,
                    adaptive_stop_reason="no_accepted_repeat_construct",
                    adaptive_maximum_proof_status="no_accepted_repeat_construct",
                    ga_adaptive_constraints_passed=False,
                    adaptive_orderable_passed=False,
                    adaptive_orderability_gate=True,
                    adaptive_idt_scored_evaluations=0,
                    adaptive_idt_cache_hits=0,
                    adaptive_idt_status_counts_json="{}",
                )
            rows.append(payload)
            continue
        if adaptive_copy_search_enabled:
            payload = _adaptive_refine_payload(
                payload,
                codon_weights=codon_weights,
                recognition_sites=recognition_sites,
                seed=seed,
                population_size=population_size,
                short_generations=short_generations,
                generation_schedule=generation_schedule,
                idt_scorer=idt_scorer,
                require_idt_orderable=True,
            )
            if (
                payload.get("ga_status") == "passed"
                and payload.get("idt_status") == "passed"
                and int(payload.get("verified_max_copies", 0) or 0) >= 2
                and isinstance(payload.get("dna_sequence"), str)
            ):
                known_orderable_by_module[module_key] = {
                    "known_orderable_copies": int(payload["verified_max_copies"]),
                    "known_orderable_dna_sequence": payload["dna_sequence"],
                    "known_orderable_dna_pre_ga": payload.get(
                        "dna_sequence_pre_ga", payload["dna_sequence"]
                    ),
                    "known_orderable_source_fragment_limit_bp": int(
                        payload.get("fragment_limit_bp", 0)
                    ),
                }
        else:
            locked_positions, selected_limits = _locked_positions_and_limits(payload)
            refined, metrics = genetic_refine_dna(
                dna,
                locked_positions=locked_positions,
                selected_site_limits=selected_limits,
                recognition_sites=recognition_sites,
                codon_weights=codon_weights,
                seed=seed + shard_index * 1_000_003 + row_number,
                population_size=population_size,
                generations=generations,
            )
            payload.update(metrics)
            payload["dna_sequence_pre_ga"] = dna
            payload["dna_sequence"] = refined
            payload["ga_status"] = (
                "passed" if metrics["ga_local_constraints_passed"] else "failed"
            )
        if not adaptive_copy_search_enabled:
            payload["idt_status"] = "not_requested"
            payload["idt_api_called"] = False
        payload["final_passed"] = False
        rows.append(payload)
        final_dna = payload.get("dna_sequence")
        if (
            use_idt
            and not adaptive_copy_search_enabled
            and payload["ga_status"] == "passed"
            and isinstance(final_dna, str)
        ):
            idt_inputs.append(
                {
                    "Name": f"{payload['module_id']}|cap={payload['fragment_limit_bp']}",
                    "Sequence": final_dna,
                }
            )
            idt_row_positions.append(len(rows) - 1)

    if use_idt and not adaptive_copy_search_enabled:
        if not credentials_available():
            for position in idt_row_positions:
                rows[position]["idt_status"] = "not_run_missing_credentials"
        elif idt_inputs:
            response = screen_gblock_sequences(idt_inputs)
            if (
                isinstance(response, list)
                and response
                and all(isinstance(item, list) for item in response)
                and len(response) != len(idt_inputs)
            ):
                raise RuntimeError(
                    "IDT response cardinality did not match the submitted sequence batch"
                )
            write_cached_response(
                response,
                destination / "idt_complexity_response.json",
                sequences=idt_inputs,
            )
            # The public Swagger response is an untyped object. Until explicit
            # per-sequence fields are present, conservatively inspect a result
            # matched by Name (falling back to the batch response) and never
            # infer success from HTTP 200 alone.
            for batch_index, (position, item) in enumerate(
                zip(idt_row_positions, idt_inputs, strict=True)
            ):
                summary = summarize_complexity_response(
                    response,
                    name=item["Name"],
                    sequence_index=batch_index,
                )
                rows[position].update(summary)
                rows[position]["idt_api_called"] = True
                rows[position]["idt_scored_sequence_length_bp"] = len(item["Sequence"])
                rows[position]["idt_scored_sequence_sha256"] = hashlib.sha256(
                    item["Sequence"].encode()
                ).hexdigest()
                rows[position]["idt_scored_sequence_unchanged"] = bool(
                    item["Sequence"] == rows[position]["dna_sequence"]
                )
                rows[position]["idt_raw_response_file"] = str(destination / "idt_complexity_response.json")
    for payload in rows:
        payload["final_passed"] = bool(
            payload.get("ga_status") == "passed" and payload.get("idt_status") == "passed"
        )
        if payload.get("ga_status") != "passed":
            payload["final_status"] = "failed_local_optimization"
        elif payload.get("idt_status") == "passed":
            payload["final_status"] = "passed_local_and_idt"
        else:
            payload["final_status"] = f"passed_local_idt_{payload.get('idt_status')}"
    result = pd.DataFrame(rows)
    result.to_parquet(destination / "optimized_constructs_ga.parquet", index=False)
    result.to_csv(destination / "optimized_constructs_ga.csv", index=False)
    with (destination / "optimized_constructs_ga.fasta").open("w") as handle:
        for row in result.itertuples(index=False):
            dna = getattr(row, "dna_sequence", None)
            if not isinstance(dna, str) or not dna:
                continue
            handle.write(f">{row.module_id}|cap={row.fragment_limit_bp}|copies={row.verified_max_copies}|idt={row.idt_status}\n")
            for start in range(0, len(dna), 80):
                handle.write(dna[start : start + 80] + "\n")
    idt_called_series = result.get(
        "idt_api_called", pd.Series(False, index=result.index)
    ).fillna(False).astype(bool)
    idt_unchanged_series = result.get(
        "idt_scored_sequence_unchanged", pd.Series(False, index=result.index)
    ).fillna(False).astype(bool)
    validation = {
        "rows": len(result),
        "ga_passed": int(result.ga_status.eq("passed").sum()),
        "adaptive_copy_search": adaptive_copy_search_enabled,
        "adaptive_copy_search_passed": int(
            result.get("ga_adaptive_constraints_passed", pd.Series(dtype=bool))
            .fillna(False)
            .sum()
        ),
        "adaptive_orderability_gate": bool(adaptive_copy_search_enabled),
        "adaptive_orderable_passed": int(
            result.get("adaptive_orderable_passed", pd.Series(dtype=bool))
            .fillna(False)
            .sum()
        ),
        "adaptive_idt_http_calls": int(idt_scorer.api_calls if idt_scorer else 0),
        "adaptive_idt_request_attempts": int(
            idt_scorer.api_attempts if idt_scorer else 0
        ),
        "adaptive_idt_cache_hits": int(idt_scorer.cache_hits if idt_scorer else 0),
        "final_passed": int(result.final_passed.sum()),
        "idt_requested": use_idt,
        "idt_credentials_available": credentials_available(),
        "idt_api_called_rows": int(result.get("idt_api_called", pd.Series(dtype=bool)).fillna(False).sum()),
        "idt_status_counts": {
            str(status): int(count)
            for status, count in result.get("idt_status", pd.Series(dtype=str))
            .value_counts(dropna=False)
            .items()
        },
        "idt_rule_violation_total": int(
            pd.to_numeric(
                result.get("idt_violation_count", pd.Series(0, index=result.index)),
                errors="coerce",
            )
            .fillna(0)
            .sum()
        ),
        "idt_scored_sequence_changed_rows": int(
            (idt_called_series & ~idt_unchanged_series).sum()
        ),
        "ga_score_profile": GA_SCORE_PROFILE,
    }
    (destination / "ga_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    return result
