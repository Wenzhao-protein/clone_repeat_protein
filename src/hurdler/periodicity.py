"""Infer primitive repeat units from complete protein sequences.

The implementation deliberately combines three independent signals:

* a Fourier concentration score on the 20-amino-acid one-hot encoding;
* residue identity and BLOSUM62-positive self-similarity at a candidate lag;
* phase-wise conservation across the inferred repeat copies.

Source annotations are accepted as priors, but they are never silently treated
as the primitive sequence period.  Coordinates exposed by this module are
1-based and inclusive; internal calculations are 0-based and half-open.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np
from Bio.Align import substitution_matrices

from .constants import AMINO_ACIDS, validate_protein_sequence

BOUNDARY_METHOD_VERSION = "spectral-secondary-structure-v4-middle-unit"
MODULE_SELECTION_POLICY = "repeat-region-middle-unit-tie-earlier-v1"

_AA_TO_INDEX = {amino_acid: index for index, amino_acid in enumerate(AMINO_ACIDS)}
_BLOSUM62 = substitution_matrices.load("BLOSUM62")


@dataclass(frozen=True)
class RepeatCandidate:
    """One candidate period and its best supported contiguous repeat block."""

    period: int
    local_start: int
    local_end: int
    repeat_count: int
    adjacent_identity: float
    adjacent_positive_fraction: float
    identity_enrichment: float
    phase_conservation: float
    fixed_fraction: float
    spectral_concentration: float
    coverage_fraction: float
    prior_harmonic: bool
    prior_overlap: bool
    prior_boundary_aligned: bool
    score: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RepeatBoundaryResult:
    """Selected primitive repeat block and auditable supporting measurements."""

    full_sequence: str
    full_sequence_origin: int
    full_sequence_sha256: str
    period: int
    repeat_region_start: int
    repeat_region_end: int
    first_module_start: int
    first_module_end: int
    first_module_sequence: str
    selected_module_index: int
    selected_module_start: int
    selected_module_end: int
    selected_module_sequence: str
    selected_module_policy: str
    repeat_count: int
    unit_sequences: tuple[str, ...]
    consensus_sequence: str
    position_conservation: tuple[float, ...]
    fixed_mask: str
    fixed_positions: tuple[int, ...]
    variable_positions: tuple[int, ...]
    variable_ranges: tuple[tuple[int, int], ...]
    left_flank_sequence: str
    right_flank_sequence: str
    score: float
    confidence: str
    selection_reason: str
    prior_period: int | None
    prior_unit_start: int | None
    prior_unit_end: int | None
    harmonic_ratio: float | None
    boundary_method_version: str = BOUNDARY_METHOD_VERSION

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["unit_sequences_json"] = json.dumps(self.unit_sequences)
        payload["position_conservation_json"] = json.dumps(self.position_conservation)
        payload["fixed_positions_json"] = json.dumps(self.fixed_positions)
        payload["variable_positions_json"] = json.dumps(self.variable_positions)
        payload["variable_ranges_json"] = json.dumps(self.variable_ranges)
        for field in (
            "unit_sequences",
            "position_conservation",
            "fixed_positions",
            "variable_positions",
            "variable_ranges",
        ):
            payload.pop(field)
        return payload


def _composition_identity_baseline(sequence: str) -> float:
    counts = np.bincount([_AA_TO_INDEX[amino_acid] for amino_acid in sequence], minlength=20)
    frequencies = counts / counts.sum()
    return float(np.square(frequencies).sum())


def _pair_scores(first: str, second: str) -> tuple[float, float]:
    identities = 0
    positive = 0
    for left, right in zip(first, second, strict=True):
        identities += left == right
        positive += float(_BLOSUM62[left, right]) > 0
    length = len(first)
    return identities / length, positive / length


def _position_statistics(
    units: Iterable[str], *, fixed_threshold: float
) -> tuple[str, tuple[float, ...], str, tuple[int, ...], tuple[int, ...], tuple[tuple[int, int], ...]]:
    unit_list = list(units)
    if not unit_list or len({len(unit) for unit in unit_list}) != 1:
        raise ValueError("Repeat units must be a non-empty equal-length collection")
    consensus: list[str] = []
    conservation: list[float] = []
    for column in zip(*unit_list, strict=True):
        counts: dict[str, int] = {}
        for amino_acid in column:
            counts[amino_acid] = counts.get(amino_acid, 0) + 1
        best = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
        consensus.append(best[0])
        conservation.append(best[1] / len(unit_list))
    fixed = tuple(index + 1 for index, value in enumerate(conservation) if value >= fixed_threshold)
    variable = tuple(index + 1 for index, value in enumerate(conservation) if value < fixed_threshold)
    ranges: list[tuple[int, int]] = []
    for position in variable:
        if not ranges or position != ranges[-1][1] + 1:
            ranges.append((position, position))
        else:
            ranges[-1] = (ranges[-1][0], position)
    mask = "".join("F" if value >= fixed_threshold else "V" for value in conservation)
    return "".join(consensus), tuple(conservation), mask, fixed, variable, tuple(ranges)


def _spectral_concentration(units: list[str]) -> float:
    """Return the fraction of non-DC Fourier power at repeat harmonics.

    For ``q`` equal-length copies, an exact tandem repeat has non-DC energy
    only at FFT bins divisible by ``q``.  Variable residues spread energy away
    from those bins, producing a naturally bounded 0--1 score.
    """
    sequence = "".join(units)
    encoded = np.zeros((len(sequence), 20), dtype=np.float64)
    encoded[np.arange(len(sequence)), [_AA_TO_INDEX[value] for value in sequence]] = 1.0
    encoded -= encoded.mean(axis=0, keepdims=True)
    power = np.square(np.abs(np.fft.fft(encoded, axis=0))).sum(axis=1)
    power[0] = 0.0
    total = float(power.sum())
    if total <= 0:
        return 0.0
    copy_count = len(units)
    harmonic_bins = np.arange(copy_count, len(sequence), copy_count)
    return float(power[harmonic_bins].sum() / total)


def _candidate_score(
    *,
    identity: float,
    positive: float,
    identity_enrichment: float,
    phase_conservation: float,
    spectral: float,
    coverage: float,
    repeat_count: int,
    prior_harmonic: bool,
    prior_overlap: bool,
    prior_boundary_aligned: bool,
) -> float:
    # Two-copy phase conservation is inflated by construction, so identity
    # enrichment and independent third/fourth copies carry explicit weight.
    copy_support = min(1.0, max(0.0, (repeat_count - 2) / 3))
    score = (
        0.24 * identity_enrichment
        + 0.18 * positive
        + 0.18 * phase_conservation
        + 0.18 * spectral
        + 0.08 * min(1.0, coverage * 2)
        + 0.08 * copy_support
        + 0.02 * float(prior_harmonic)
        + 0.02 * float(prior_overlap)
        + 0.02 * float(prior_boundary_aligned)
    )
    # A pair alone is useful evidence, but not enough to drive an aggressive
    # harmonic split unless it is nearly exact.
    if repeat_count == 2 and identity < 0.8:
        score -= 0.08
    return float(max(0.0, min(1.0, score)))


def _best_block_for_period(
    sequence: str,
    period: int,
    *,
    fixed_threshold: float,
    prior_period: int | None,
    prior_local_start: int | None,
    prior_local_end: int | None,
) -> RepeatCandidate | None:
    baseline = _composition_identity_baseline(sequence)
    best: RepeatCandidate | None = None
    # Every absolute boundary belongs to exactly one of these phases.  A
    # contiguous sub-block of the phase tiling can therefore start anywhere.
    for phase in range(min(period, len(sequence) - 2 * period + 1)):
        tiled_units = [
            sequence[start : start + period]
            for start in range(phase, len(sequence) - period + 1, period)
        ]
        if len(tiled_units) < 2:
            continue
        for first_index in range(len(tiled_units) - 1):
            # Cap exhaustive spans to keep pathological low-period sequences
            # inexpensive while retaining the longest supported region.
            for end_index in range(first_index + 2, len(tiled_units) + 1):
                units = tiled_units[first_index:end_index]
                local_start = phase + first_index * period
                local_end = local_start + len(units) * period
                prior_overlap = bool(
                    prior_local_start is not None
                    and prior_local_end is not None
                    and local_start < prior_local_end
                    and local_end > prior_local_start
                )
                if prior_local_start is not None and prior_local_end is not None and not prior_overlap:
                    continue
                pair_values = [_pair_scores(left, right) for left, right in zip(units, units[1:])]
                identity = float(np.mean([value[0] for value in pair_values]))
                positive = float(np.mean([value[1] for value in pair_values]))
                enrichment = max(0.0, (identity - baseline) / max(1e-9, 1.0 - baseline))
                consensus, conservation, _, fixed, _, _ = _position_statistics(
                    units, fixed_threshold=fixed_threshold
                )
                del consensus
                phase_conservation = float(np.mean(conservation))
                fixed_fraction = len(fixed) / period
                # Fast rejection prevents thousands of FFTs on unsupported
                # spans, especially for p=1..5 in long natural proteins.
                if identity < max(0.12, baseline + 0.035) and positive < 0.52:
                    continue
                prior_harmonic = bool(
                    prior_period
                    and prior_period >= period
                    and abs(prior_period / period - round(prior_period / period)) < 0.08
                )
                prior_boundary_aligned = bool(
                    prior_local_start is not None
                    and (prior_local_start - local_start) % period == 0
                )
                spectral = _spectral_concentration(units)
                coverage = (local_end - local_start) / len(sequence)
                score = _candidate_score(
                    identity=identity,
                    positive=positive,
                    identity_enrichment=enrichment,
                    phase_conservation=phase_conservation,
                    spectral=spectral,
                    coverage=coverage,
                    repeat_count=len(units),
                    prior_harmonic=prior_harmonic,
                    prior_overlap=prior_overlap,
                    prior_boundary_aligned=prior_boundary_aligned,
                )
                candidate = RepeatCandidate(
                    period=period,
                    local_start=local_start,
                    local_end=local_end,
                    repeat_count=len(units),
                    adjacent_identity=identity,
                    adjacent_positive_fraction=positive,
                    identity_enrichment=enrichment,
                    phase_conservation=phase_conservation,
                    fixed_fraction=fixed_fraction,
                    spectral_concentration=spectral,
                    coverage_fraction=coverage,
                    prior_harmonic=prior_harmonic,
                    prior_overlap=prior_overlap,
                    prior_boundary_aligned=prior_boundary_aligned,
                    score=score,
                )
                if best is None or (
                    candidate.score,
                    candidate.repeat_count,
                    candidate.coverage_fraction,
                    -candidate.local_start,
                ) > (
                    best.score,
                    best.repeat_count,
                    best.coverage_fraction,
                    -best.local_start,
                ):
                    best = candidate
    return best


def _quick_lag_score(sequence: str, period: int) -> float:
    """Cheap pre-screen used before exhaustive boundary optimization."""
    if len(sequence) < 2 * period:
        return 0.0
    identity, positive = _pair_scores(sequence[:-period], sequence[period:])
    baseline = _composition_identity_baseline(sequence)
    enrichment = max(0.0, (identity - baseline) / max(1e-9, 1.0 - baseline))
    return 0.65 * enrichment + 0.35 * positive


def scan_repeat_periods(
    full_sequence: str,
    *,
    minimum_period: int = 3,
    maximum_period: int = 120,
    fixed_threshold: float = 0.8,
    prior_period: int | None = None,
    prior_unit_start: int | None = None,
    full_sequence_origin: int = 1,
) -> list[RepeatCandidate]:
    """Score candidate primitive periods on a complete protein sequence."""
    sequence = validate_protein_sequence(full_sequence)
    # A single-residue run has no sequence information from which to infer a
    # biological module boundary, even though every numerical lag is perfect.
    if len(set(sequence)) == 1:
        return []
    if minimum_period < 1:
        raise ValueError("minimum_period must be positive")
    if not 0.5 <= fixed_threshold <= 1.0:
        raise ValueError("fixed_threshold must be between 0.5 and 1.0")
    upper = min(maximum_period, len(sequence) // 2, prior_period or maximum_period)
    prior_local_start = None
    prior_local_end = None
    if prior_unit_start is not None and prior_period is not None:
        prior_local_start = prior_unit_start - full_sequence_origin
        prior_local_end = prior_local_start + prior_period
    all_periods = list(range(minimum_period, upper + 1))
    if prior_period is not None and prior_local_start is not None:
        margin = max(prior_period * 3, 30)
        window_start = max(0, prior_local_start - margin)
        window_end = min(len(sequence), (prior_local_end or prior_local_start) + margin)
        local_window = sequence[window_start:window_end]
        quick = sorted(
            ((_quick_lag_score(local_window, period), period) for period in all_periods),
            reverse=True,
        )
        selected_periods = {period for _, period in quick[:24]}
        selected_periods.update(
            period
            for period in all_periods
            if abs(prior_period / period - round(prior_period / period)) < 0.08
        )
        selected_periods.add(prior_period)
        all_periods = sorted(selected_periods)
    candidates = [
        candidate
        for period in all_periods
        if (
            candidate := _best_block_for_period(
                sequence,
                period,
                fixed_threshold=fixed_threshold,
                prior_period=prior_period,
                prior_local_start=prior_local_start,
                prior_local_end=prior_local_end,
            )
        )
        is not None
    ]
    return sorted(candidates, key=lambda value: (-value.score, value.period))


def select_primitive_candidate(
    candidates: list[RepeatCandidate], *, prior_period: int | None
) -> tuple[RepeatCandidate | None, str]:
    if not candidates:
        raise ValueError("No supported repeat period found")
    strongest = candidates[0]
    # A primitive harmonic must be independently supported by >=3 copies and
    # remain close to the best overall evidence.  This is what prevents a
    # longest-exact-repeat method from returning two or more true modules.
    primitive_pool = [
        candidate
        for candidate in candidates
        if candidate.repeat_count >= 4
        and candidate.score >= 0.64
        and candidate.adjacent_identity >= 0.35
        and candidate.adjacent_positive_fraction >= 0.60
        and candidate.spectral_concentration >= 0.60
        and (
            prior_period is None
            or candidate.period * candidate.repeat_count >= 2 * prior_period
        )
    ]
    if prior_period:
        harmonic_pool = [
            candidate
            for candidate in primitive_pool
            if candidate.period <= prior_period
            and abs(prior_period / candidate.period - round(prior_period / candidate.period)) < 0.08
        ]
        if harmonic_pool:
            selected = min(harmonic_pool, key=lambda value: (value.period, -value.score))
            if selected.period < prior_period:
                return selected, "smallest independently supported harmonic of source-annotated unit"
    if primitive_pool:
        selected = min(primitive_pool, key=lambda value: (value.period, -value.score))
        if selected.period < strongest.period:
            return selected, "smallest independently supported period near strongest spectral peak"
    if prior_period:
        source_period_pool = [
            candidate
            for candidate in candidates
            if abs(candidate.period / prior_period - 1.0) <= 0.08
            and candidate.score >= 0.80
            and candidate.adjacent_identity >= 0.80
        ]
        if source_period_pool:
            return max(source_period_pool, key=lambda value: value.score), (
                "source period confirmed by strong full-sequence self-similarity"
            )
        return None, (
            "no independently supported sequence subperiod; source-annotated unit retained"
        )
    if strongest.repeat_count >= 3 and strongest.score >= 0.64:
        return strongest, "strongest combined spectral and self-similarity evidence"
    raise ValueError("No independently supported repeat period found")


def materialize_repeat_boundary(
    full_sequence: str,
    selected: RepeatCandidate,
    *,
    selection_reason: str,
    full_sequence_origin: int = 1,
    prior_period: int | None = None,
    prior_unit_start: int | None = None,
    fixed_threshold: float = 0.8,
) -> RepeatBoundaryResult:
    """Build coordinates, flanks, and fixed/variable calls from a candidate."""
    sequence = validate_protein_sequence(full_sequence)
    local_start = selected.local_start
    local_end = selected.local_end
    units = tuple(
        sequence[position : position + selected.period]
        for position in range(local_start, local_end, selected.period)
    )
    consensus, conservation, mask, fixed, variable, variable_ranges = _position_statistics(
        units, fixed_threshold=fixed_threshold
    )
    absolute_start = full_sequence_origin + local_start
    absolute_end = full_sequence_origin + local_end - 1
    selected_index_zero_based = (len(units) - 1) // 2
    selected_start = absolute_start + selected_index_zero_based * selected.period
    harmonic_ratio = prior_period / selected.period if prior_period else None
    if selected.score >= 0.72 and selected.repeat_count >= 3:
        confidence = "high"
    elif selected.score >= 0.58:
        confidence = "medium"
    else:
        confidence = "low"
    return RepeatBoundaryResult(
        full_sequence=sequence,
        full_sequence_origin=full_sequence_origin,
        full_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
        period=selected.period,
        repeat_region_start=absolute_start,
        repeat_region_end=absolute_end,
        first_module_start=absolute_start,
        first_module_end=absolute_start + selected.period - 1,
        first_module_sequence=units[0],
        selected_module_index=selected_index_zero_based + 1,
        selected_module_start=selected_start,
        selected_module_end=selected_start + selected.period - 1,
        selected_module_sequence=units[selected_index_zero_based],
        selected_module_policy=MODULE_SELECTION_POLICY,
        repeat_count=len(units),
        unit_sequences=units,
        consensus_sequence=consensus,
        position_conservation=conservation,
        fixed_mask=mask,
        fixed_positions=fixed,
        variable_positions=variable,
        variable_ranges=variable_ranges,
        left_flank_sequence=sequence[:local_start],
        right_flank_sequence=sequence[local_end:],
        score=selected.score,
        confidence=confidence,
        selection_reason=selection_reason,
        prior_period=prior_period,
        prior_unit_start=prior_unit_start,
        prior_unit_end=(prior_unit_start + prior_period - 1)
        if prior_unit_start is not None and prior_period is not None
        else None,
        harmonic_ratio=harmonic_ratio,
    )


def infer_repeat_boundaries(
    full_sequence: str,
    *,
    full_sequence_origin: int = 1,
    prior_unit_sequence: str | None = None,
    prior_unit_start: int | None = None,
    minimum_period: int = 3,
    maximum_period: int = 120,
    fixed_threshold: float = 0.8,
) -> tuple[RepeatBoundaryResult, list[RepeatCandidate]]:
    """Infer a primitive repeat block and return the complete candidate audit."""
    sequence = validate_protein_sequence(full_sequence)
    prior_sequence = validate_protein_sequence(prior_unit_sequence) if prior_unit_sequence else None
    prior_period = len(prior_sequence) if prior_sequence else None
    candidates = scan_repeat_periods(
        sequence,
        minimum_period=minimum_period,
        maximum_period=maximum_period,
        fixed_threshold=fixed_threshold,
        prior_period=prior_period,
        prior_unit_start=prior_unit_start,
        full_sequence_origin=full_sequence_origin,
    )
    selected, reason = select_primitive_candidate(candidates, prior_period=prior_period)
    if selected is None:
        if prior_period is None or prior_unit_start is None:
            raise ValueError("A source boundary is required when sequence-period evidence is insufficient")
        local_start = prior_unit_start - full_sequence_origin
        local_end = local_start + prior_period
        if local_start < 0 or local_end > len(sequence):
            raise ValueError("Source unit falls outside the complete protein sequence")
        selected = RepeatCandidate(
            period=prior_period,
            local_start=local_start,
            local_end=local_end,
            repeat_count=1,
            adjacent_identity=0.0,
            adjacent_positive_fraction=0.0,
            identity_enrichment=0.0,
            phase_conservation=0.0,
            fixed_fraction=0.0,
            spectral_concentration=0.0,
            coverage_fraction=prior_period / len(sequence),
            prior_harmonic=True,
            prior_overlap=True,
            prior_boundary_aligned=True,
            score=0.0,
        )
    result = materialize_repeat_boundary(
        sequence,
        selected,
        selection_reason=reason,
        full_sequence_origin=full_sequence_origin,
        prior_period=prior_period,
        prior_unit_start=prior_unit_start,
        fixed_threshold=fixed_threshold,
    )
    return result, candidates
