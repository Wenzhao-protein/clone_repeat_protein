"""Interactive, sequence-only HURDLER construct design workflow.

The versioned protein pattern index is the sole authority for HURDLER
compatibility.  DNA-level calculations only optimize and audit a selected
protein-derived solution; they can never create a new compatibility result.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Protocol

import pandas as pd
from Bio import Restriction, SeqIO

from .constants import PLASMIDS, RULE_PROFILE_NAME, validate_protein_sequence
from .ga_optimization import (
    GA_SCORE_PROFILE,
    adjust_ga_score_profile_from_idt,
    genetic_refine_dna,
    load_restriction_sites,
)
from .index import PatternIndex
from .io import sha256_file, utc_now, write_json_atomic
from .matching import enumerate_module_solutions
from .optimization import (
    _candidate_rank,
    diversify_codons,
    load_codon_weights,
    recognition_site_count,
    reverse_complement,
    translate_dna,
)
from .paths import ProjectPaths
from .periodicity import RepeatCandidate, scan_repeat_periods


DESIGN_SCHEMA_VERSION = "interactive-hurdler-designer-v1"
BOUNDARY_METHOD_VERSION = "sequence-self-periodicity-confirmed-v1"
GENERATION_SCHEDULE = (10, 20, 40, 60, 80, 100)
GBLOCK_MIN_BP = 125
GBLOCK_MAX_BP = 3000
ASSEMBLY_OVERLAP_BP = 30

PLASMID_MCS: dict[str, tuple[str, str]] = {
    "pGEX-4T-1": ("BamHI", "EagI"),
    "pMAL-c5X": ("NdeI", "HindIII"),
    "pET-21a(+)": ("NdeI", "XhoI"),
    "pET-28a(+)": ("BamHI", "XhoI"),
    "pET-28a(+)_start_codon": ("NcoI", "XhoI"),
    "pCold_I": ("NdeI", "BspMI"),
    "pUC18": ("EcoRI", "HindIII"),
    "pQE-3": ("BamHI", "HindIII"),
}


class ComplexityScorer(Protocol):
    def score(self, name: str, sequence: str) -> dict[str, Any]: ...


def bundled_index_dir() -> Path:
    """Return the committed full lookup used by clone-and-run workflows."""
    return ProjectPaths.discover().root / "data" / "artifacts" / RULE_PROFILE_NAME


def parse_protein_input(value: str) -> tuple[str, str]:
    """Parse one raw amino-acid sequence or one FASTA record."""
    text = str(value).strip()
    if not text:
        raise ValueError("Protein input is empty")
    if text.startswith(">"):
        lines = text.splitlines()
        headers = [line for line in lines if line.startswith(">")]
        if len(headers) != 1:
            raise ValueError("Exactly one FASTA record is required")
        identifier = headers[0][1:].strip().split()[0] or "interactive_input"
        sequence = "".join(line.strip() for line in lines[1:] if not line.startswith(">"))
    else:
        identifier = "interactive_input"
        sequence = text
    return identifier, validate_protein_sequence(sequence)


@dataclass(frozen=True)
class BoundaryCandidate:
    period: int
    repeat_region_start: int
    repeat_region_end: int
    repeat_count: int
    score: float
    adjacent_identity: float
    phase_conservation: float
    coverage_fraction: float
    evidence: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class BoundaryAnalysis:
    sequence_id: str
    full_protein_sequence: str
    candidates: tuple[BoundaryCandidate, ...]
    proposed_start: int | None
    proposed_end: int | None
    proposed_period: int | None
    proposed_confidence: str

    def to_dict(self) -> dict[str, object]:
        return {
            **{key: value for key, value in asdict(self).items() if key != "candidates"},
            "candidates": [candidate.to_dict() for candidate in self.candidates],
        }


@dataclass(frozen=True)
class ConfirmedBoundary:
    repeat_region_start: int
    repeat_region_end: int
    period: int
    repeat_count: int
    unit_sequences: tuple[str, ...]
    middle_unit_index: int
    middle_unit_start: int
    middle_unit_end: int
    middle_module: str
    consensus_module: str
    position_conservation: tuple[float, ...]
    fixed_positions_1based: tuple[int, ...]
    variable_ranges_1based: tuple[tuple[int, int], ...]
    n_terminal_flank: str
    c_terminal_flank: str
    confirmation_token: str

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["unit_sequences"] = list(self.unit_sequences)
        payload["position_conservation"] = list(self.position_conservation)
        payload["fixed_positions_1based"] = list(self.fixed_positions_1based)
        payload["variable_ranges_1based"] = [list(value) for value in self.variable_ranges_1based]
        return payload


@dataclass(frozen=True)
class DesignRequest:
    """Strict public request schema shared by the CLI and notebook."""

    full_protein_sequence: str
    target_repeat_copies: int
    plasmid: str
    sequence_id: str = "interactive_input"
    confirmed_repeat_start: int | None = None
    confirmed_repeat_end: int | None = None
    confirmed_period: int | None = None
    confirmation_token: str | None = None
    site_i_allowlist: tuple[str, ...] = ()
    site_ii_allowlist: tuple[str, ...] = ()
    site_iii_allowlist: tuple[str, ...] = ()
    selected_candidate_id: str | None = None
    optimize: bool = False
    seed: int = 42
    population_size: int = 16
    max_purchase_bp: int = GBLOCK_MAX_BP

    def __post_init__(self) -> None:
        raw_protein_input = str(self.full_protein_sequence)
        parsed_id, sequence = parse_protein_input(raw_protein_input)
        object.__setattr__(self, "full_protein_sequence", sequence)
        requested_id = str(self.sequence_id or "interactive_input")
        if requested_id == "interactive_input" and raw_protein_input.lstrip().startswith(">"):
            requested_id = parsed_id
        object.__setattr__(self, "sequence_id", requested_id)
        for name in ("site_i_allowlist", "site_ii_allowlist", "site_iii_allowlist"):
            value = getattr(self, name)
            if isinstance(value, str):
                normalized = tuple(part.strip() for part in value.split(",") if part.strip())
            else:
                normalized = tuple(str(part).strip() for part in value if str(part).strip())
            object.__setattr__(self, name, normalized)
        if self.plasmid not in PLASMIDS:
            raise ValueError(f"Unsupported plasmid {self.plasmid!r}; choose one of {', '.join(PLASMIDS)}")
        if int(self.target_repeat_copies) < 2:
            raise ValueError("target_repeat_copies must be at least 2")
        if int(self.population_size) < 2:
            raise ValueError("population_size must be at least 2")
        if not GBLOCK_MIN_BP <= int(self.max_purchase_bp) <= GBLOCK_MAX_BP:
            raise ValueError("max_purchase_bp must be between 125 and 3000")
        coordinates = (
            self.confirmed_repeat_start,
            self.confirmed_repeat_end,
            self.confirmed_period,
        )
        if any(value is not None for value in coordinates) and not all(
            value is not None for value in coordinates
        ):
            raise ValueError("Confirmed start, end, and period must be provided together")

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "DesignRequest":
        allowed = {item.name for item in fields(cls)}
        unknown = sorted(set(payload) - allowed)
        if unknown:
            raise ValueError("Unknown DesignRequest fields: " + ", ".join(unknown))
        return cls(**dict(payload))

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        for name in ("site_i_allowlist", "site_ii_allowlist", "site_iii_allowlist"):
            payload[name] = list(payload[name])
        return payload


@dataclass
class DesignResult:
    schema_version: str
    status: str
    message: str
    request: dict[str, object]
    boundary_analysis: dict[str, object]
    confirmed_boundary: dict[str, object] | None = None
    target_protein_sequence: str | None = None
    candidates: list[dict[str, object]] = field(default_factory=list)
    selected_candidate: dict[str, object] | None = None
    optimization: dict[str, object] | None = None
    purchase_fragments: list[dict[str, object]] = field(default_factory=list)
    cloning_steps: list[dict[str, object]] = field(default_factory=list)
    final_plasmid: dict[str, object] | None = None
    idt_audit: list[dict[str, object]] = field(default_factory=list)
    output_files: dict[str, str] = field(default_factory=dict)

    @property
    def orderable(self) -> bool:
        return self.status == "orderable_design_complete"

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _exact_short_period_candidates(sequence: str) -> list[BoundaryCandidate]:
    """Find maximal exact period-1/2 blocks, including homopolymers."""
    candidates: list[BoundaryCandidate] = []
    for period in (1, 2):
        if len(sequence) < 2 * period:
            continue
        comparisons = [sequence[index] == sequence[index - period] for index in range(period, len(sequence))]
        run_start: int | None = None
        for offset, agrees in enumerate([*comparisons, False]):
            if agrees and run_start is None:
                run_start = offset
            if agrees or run_start is None:
                continue
            run_end = offset
            if run_end - run_start >= period:
                local_start = run_start
                local_end = period + run_end
                length = local_end - local_start
                length -= length % period
                local_end = local_start + length
                repeat_count = length // period
                if repeat_count >= 2:
                    candidates.append(
                        BoundaryCandidate(
                            period=period,
                            repeat_region_start=local_start + 1,
                            repeat_region_end=local_end,
                            repeat_count=repeat_count,
                            score=min(1.0, 0.82 + 0.03 * min(repeat_count, 6)),
                            adjacent_identity=1.0,
                            phase_conservation=1.0,
                            coverage_fraction=length / len(sequence),
                            evidence=f"exact_period_{period}_self_alignment",
                        )
                    )
            run_start = None
    return candidates


def _convert_period_candidate(candidate: RepeatCandidate, sequence_length: int) -> BoundaryCandidate:
    return BoundaryCandidate(
        period=int(candidate.period),
        repeat_region_start=int(candidate.local_start) + 1,
        repeat_region_end=int(candidate.local_end),
        repeat_count=int(candidate.repeat_count),
        score=float(candidate.score),
        adjacent_identity=float(candidate.adjacent_identity),
        phase_conservation=float(candidate.phase_conservation),
        coverage_fraction=(int(candidate.local_end) - int(candidate.local_start)) / sequence_length,
        evidence="self_alignment_periodicity_conservation",
    )


def analyze_repeat_sequence(
    protein_input: str,
    *,
    sequence_id: str | None = None,
    maximum_period: int = 120,
) -> BoundaryAnalysis:
    """Return fast sequence-only period candidates for mandatory confirmation."""
    parsed_id, sequence = parse_protein_input(protein_input)
    resolved_id = sequence_id or parsed_id
    raw: list[BoundaryCandidate] = _exact_short_period_candidates(sequence)
    if len(sequence) >= 6:
        try:
            inferred = scan_repeat_periods(
                sequence,
                minimum_period=3,
                maximum_period=min(int(maximum_period), max(3, len(sequence) // 2)),
            )
        except ValueError:
            inferred = []
        raw.extend(_convert_period_candidate(item, len(sequence)) for item in inferred)
    deduplicated: dict[tuple[int, int, int], BoundaryCandidate] = {}
    for candidate in raw:
        key = (candidate.period, candidate.repeat_region_start, candidate.repeat_region_end)
        old = deduplicated.get(key)
        if old is None or candidate.score > old.score:
            deduplicated[key] = candidate
    candidates = sorted(
        deduplicated.values(),
        key=lambda item: (-item.score, item.period, -item.repeat_count, item.repeat_region_start),
    )
    selected = candidates[0] if candidates else None
    confidence = (
        "high" if selected and selected.score >= 0.72
        else "medium" if selected and selected.score >= 0.58
        else "low" if selected
        else "not_inferred"
    )
    return BoundaryAnalysis(
        sequence_id=str(resolved_id),
        full_protein_sequence=sequence,
        candidates=tuple(candidates),
        proposed_start=selected.repeat_region_start if selected else None,
        proposed_end=selected.repeat_region_end if selected else None,
        proposed_period=selected.period if selected else None,
        proposed_confidence=confidence,
    )


def boundary_confirmation_token(sequence: str, start: int, end: int, period: int) -> str:
    normalized = validate_protein_sequence(sequence)
    return hashlib.sha256(f"{normalized}|{start}|{end}|{period}".encode()).hexdigest()


def confirm_repeat_boundary(
    sequence: str,
    *,
    start: int,
    end: int,
    period: int,
    expected_token: str | None = None,
) -> ConfirmedBoundary:
    """Validate 1-based inclusive coordinates and select the earlier middle unit."""
    protein = validate_protein_sequence(sequence)
    start = int(start)
    end = int(end)
    period = int(period)
    if start < 1 or end > len(protein) or start > end:
        raise ValueError("Confirmed repeat coordinates fall outside the protein sequence")
    if period < 1:
        raise ValueError("Confirmed period must be positive")
    region = protein[start - 1 : end]
    if len(region) % period:
        raise ValueError("Confirmed repeat-region length must be divisible by the period")
    repeat_count = len(region) // period
    if repeat_count < 2:
        raise ValueError("Confirmed repeat region must contain at least two complete units")
    units = tuple(region[offset : offset + period] for offset in range(0, len(region), period))
    middle_index = (repeat_count - 1) // 2
    middle_start = start + middle_index * period
    consensus_chars: list[str] = []
    conservation: list[float] = []
    for position in range(period):
        counts: dict[str, int] = {}
        for unit in units:
            counts[unit[position]] = counts.get(unit[position], 0) + 1
        residue, count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
        consensus_chars.append(residue)
        conservation.append(count / repeat_count)
    fixed_positions = tuple(
        position + 1 for position, value in enumerate(conservation) if value >= 0.8
    )
    variable_positions = [
        position + 1 for position, value in enumerate(conservation) if value < 0.8
    ]
    variable_ranges: list[tuple[int, int]] = []
    for position in variable_positions:
        if not variable_ranges or position != variable_ranges[-1][1] + 1:
            variable_ranges.append((position, position))
        else:
            variable_ranges[-1] = (variable_ranges[-1][0], position)
    token = boundary_confirmation_token(protein, start, end, period)
    if expected_token is not None and expected_token != token:
        raise ValueError("Boundary confirmation is stale because the sequence or coordinates changed")
    return ConfirmedBoundary(
        repeat_region_start=start,
        repeat_region_end=end,
        period=period,
        repeat_count=repeat_count,
        unit_sequences=units,
        middle_unit_index=middle_index + 1,
        middle_unit_start=middle_start,
        middle_unit_end=middle_start + period - 1,
        middle_module=units[middle_index],
        consensus_module="".join(consensus_chars),
        position_conservation=tuple(conservation),
        fixed_positions_1based=fixed_positions,
        variable_ranges_1based=tuple(variable_ranges),
        n_terminal_flank=protein[: start - 1],
        c_terminal_flank=protein[end:],
        confirmation_token=token,
    )


def role_enzyme_options(index: PatternIndex) -> dict[str, tuple[str, ...]]:
    pair_table = index.pair_table
    site_iii = {
        enzyme
        for value in pair_table["site_iii_enzymes"].fillna("").astype(str)
        for enzyme in value.split(",")
        if enzyme
    }
    return {
        "site_i": tuple(sorted(pair_table["site_i_enzyme"].astype(str).unique())),
        "site_ii": tuple(sorted(pair_table["site_ii_enzyme"].astype(str).unique())),
        "site_iii": tuple(sorted(site_iii)),
    }


def _validate_allowlist(
    requested: Iterable[str],
    *,
    role: str,
    options: Mapping[str, tuple[str, ...]],
) -> set[str] | None:
    values = {str(value) for value in requested if str(value)}
    if not values:
        return None
    allowed = set(options[role])
    invalid = sorted(values - allowed)
    if invalid:
        other_roles = {
            enzyme: [name for name, pool in options.items() if enzyme in pool]
            for enzyme in invalid
        }
        details = "; ".join(
            f"{enzyme} (available roles: {', '.join(other_roles[enzyme]) or 'none'})"
            for enzyme in invalid
        )
        raise ValueError(f"Invalid {role.replace('_', ' ').title()} whitelist: {details}")
    return values


def _python_scalar(value: Any) -> Any:
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def enumerate_design_candidates(
    module: str,
    plasmid: str,
    index: PatternIndex,
    *,
    site_i_allowlist: Iterable[str] = (),
    site_ii_allowlist: Iterable[str] = (),
    site_iii_allowlist: Iterable[str] = (),
) -> list[dict[str, object]]:
    """Enumerate and deterministically rank every allowed frozen-index route."""
    if plasmid not in PLASMIDS:
        raise ValueError(f"Unsupported plasmid: {plasmid}")
    options = role_enzyme_options(index)
    allowed_i = _validate_allowlist(site_i_allowlist, role="site_i", options=options)
    allowed_ii = _validate_allowlist(site_ii_allowlist, role="site_ii", options=options)
    allowed_iii = _validate_allowlist(site_iii_allowlist, role="site_iii", options=options)
    rows: list[dict[str, object]] = []
    for raw in enumerate_module_solutions(module, index):
        row = {key: _python_scalar(value) for key, value in raw.items()}
        if row.get("plasmid") != plasmid:
            continue
        if allowed_i is not None and row.get("site_i_enzyme") not in allowed_i:
            continue
        if allowed_ii is not None and row.get("site_ii_enzyme") not in allowed_ii:
            continue
        enzymes = str(row.get("site_iii_enzymes") or "").split(",")
        sites = str(row.get("site_iii_sites") or "").split(",")
        for position, enzyme in enumerate(enzymes):
            if not enzyme or (allowed_iii is not None and enzyme not in allowed_iii):
                continue
            expanded = dict(row)
            expanded["site_iii_enzyme"] = enzyme
            expanded["site_iii_recognition_site"] = sites[position] if position < len(sites) else ""
            expanded["site_iii_ovhg"] = expanded.get("site_ii_ovhg")
            token = "|".join(
                str(expanded.get(key, ""))
                for key in (
                    "plasmid", "pattern_key", "site_i_position", "site_ii_position",
                    "candidate_pair_id", "site_i_enzyme", "site_ii_enzyme", "site_iii_enzyme",
                )
            )
            expanded["candidate_id"] = hashlib.sha256(token.encode()).hexdigest()[:20]
            expanded["rule_profile"] = RULE_PROFILE_NAME
            expanded["scan_copy_count"] = 2
            rows.append(expanded)
    rows.sort(
        key=lambda row: (
            _candidate_rank(row),
            str(row.get("site_iii_enzyme", "")),
            str(row.get("candidate_id", "")),
        )
    )
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def _locked_construct(
    target_protein: str,
    boundary: ConfirmedBoundary,
    candidate: Mapping[str, object],
    codon_weights: dict[str, float],
) -> tuple[str, set[int], dict[str, int], dict[str, int]]:
    repeat_offset = len(boundary.n_terminal_flank)
    locks: dict[int, str] = {}
    locked_codons: dict[int, str] = {}
    absolute_positions: dict[str, int] = {}
    for role, position_key, sequence_key in (
        ("site_i", "site_i_position", "site_i_9mer_bp"),
        ("site_ii", "site_ii_position", "site_ii_9mer_bp_mutated"),
    ):
        start = repeat_offset + int(candidate[position_key])
        dna9 = str(candidate[sequence_key])
        if start < 0 or start + 3 > len(target_protein):
            raise ValueError(
                f"Target repeat count is too short for the selected {role.replace('_', ' ')} window"
            )
        if translate_dna(dna9) != target_protein[start : start + 3]:
            raise ValueError(f"Selected {role.replace('_', ' ')} window does not match the target protein")
        for offset in range(3):
            old = locked_codons.setdefault(start + offset, dna9[offset * 3 : offset * 3 + 3])
            if old != dna9[offset * 3 : offset * 3 + 3]:
                raise ValueError("Selected Site-I and Site-II locked codons conflict")
        locks[start] = dna9
        absolute_positions[role] = start
    site_limits = {
        str(candidate.get("site_i_recognition_site") or ""): 1,
        str(candidate.get("site_ii_recognition_site") or ""): 0,
    }
    site_limits = {site: limit for site, limit in site_limits.items() if site}
    dna = diversify_codons(target_protein, locks, codon_weights, site_limits)
    if translate_dna(dna) != target_protein:
        raise AssertionError("Initial construct changed the requested protein")
    return dna, set(locked_codons), site_limits, absolute_positions


def _balanced_core_intervals(
    length: int,
    *,
    prefix_length: int,
    suffix_length: int,
    max_purchase_bp: int,
    overlap: int = ASSEMBLY_OVERLAP_BP,
) -> list[tuple[int, int]]:
    if length <= 0:
        raise ValueError("Construct DNA is empty")
    if prefix_length + length + suffix_length <= max_purchase_bp:
        return [(0, length)]
    for count in range(2, 10000):
        total_core = length + overlap * (count - 1)
        span = math.ceil(total_core / count)
        intervals: list[tuple[int, int]] = []
        start = 0
        for index in range(count):
            end = length if index == count - 1 else min(length, start + span)
            intervals.append((start, end))
            if end >= length:
                break
            start = end - overlap
        if len(intervals) != count or intervals[-1][1] != length:
            continue
        purchase_lengths = [
            end - start + (prefix_length if index == 0 else 0) + (suffix_length if index == count - 1 else 0)
            for index, (start, end) in enumerate(intervals)
        ]
        if max(purchase_lengths) <= max_purchase_bp and min(end - start for start, end in intervals) >= overlap:
            return intervals
    raise ValueError("Could not split the construct into purchaseable fragments")


def plan_purchase_fragments(
    target_id: str,
    dna: str,
    plasmid: str,
    *,
    max_purchase_bp: int = GBLOCK_MAX_BP,
) -> list[dict[str, object]]:
    """Create exact MCS-adapted fragments and verify overlap reconstruction."""
    left_enzyme, right_enzyme = PLASMID_MCS[plasmid]
    left_site = str(getattr(Restriction, left_enzyme).site)
    right_site = str(getattr(Restriction, right_enzyme).site)
    intervals = _balanced_core_intervals(
        len(dna),
        prefix_length=len(left_site),
        suffix_length=len(right_site),
        max_purchase_bp=max_purchase_bp,
    )
    fragments: list[dict[str, object]] = []
    reconstructed = ""
    previous_end = 0
    for index, (start, end) in enumerate(intervals):
        core = dna[start:end]
        prefix = left_site if index == 0 else ""
        suffix = right_site if index == len(intervals) - 1 else ""
        purchase = prefix + core + suffix
        length = len(purchase)
        product_type = (
            "gBlock" if GBLOCK_MIN_BP <= length <= GBLOCK_MAX_BP
            else "duplexed_ultramer" if 20 <= length < GBLOCK_MIN_BP
            else "unsupported"
        )
        overlap_with_previous = max(0, previous_end - start) if index else 0
        reconstructed += core[overlap_with_previous:]
        previous_end = end
        fragments.append(
            {
                "fragment_id": f"{target_id}_fragment_{index + 1:02d}",
                "fragment_index": index + 1,
                "core_start_0based": start,
                "core_end_exclusive": end,
                "core_length_bp": len(core),
                "overlap_with_previous_bp": overlap_with_previous,
                "left_adapter": prefix,
                "right_adapter": suffix,
                "purchase_sequence": purchase,
                "purchase_length_bp": length,
                "product_type": product_type,
                "product_length_valid": product_type != "unsupported",
                "purchase_sha256": hashlib.sha256(purchase.encode()).hexdigest(),
                "idt_score_policy": "idt-rule-score-sum-lt10-v1",
            }
        )
    if reconstructed != dna:
        raise AssertionError("Purchase-fragment overlap reconstruction changed the optimized CDS")
    return fragments


def _score_fragments(
    fragments: list[dict[str, object]],
    scorer: ComplexityScorer,
    *,
    safe_point: Callable[[], None] | None = None,
) -> tuple[bool, list[dict[str, object]], list[dict[str, object]]]:
    scored: list[dict[str, object]] = []
    audit: list[dict[str, object]] = []
    all_passed = True
    for fragment in fragments:
        if safe_point is not None:
            safe_point()
        sequence = str(fragment["purchase_sequence"])
        expected_sha = hashlib.sha256(sequence.encode()).hexdigest()
        try:
            summary = dict(scorer.score(str(fragment["fragment_id"]), sequence))
        except Exception as exc:
            summary = {
                "idt_status": "api_failure",
                "idt_explicit_pass": None,
                "idt_complexity_score": None,
                "idt_score_complete": False,
                "idt_rule_details_json": "[]",
                "idt_positive_score_names_json": "[]",
                "idt_violation_names_json": "[]",
                "idt_scored_sequence_sha256": expected_sha,
                "idt_response_sha256": "",
                "idt_error_type": type(exc).__name__,
            }
        if safe_point is not None:
            safe_point()
        observed_sha = str(summary.get("idt_scored_sequence_sha256") or expected_sha)
        if observed_sha != expected_sha:
            raise ValueError("IDT response DNA hash does not match the exact purchase fragment")
        score = summary.get("idt_complexity_score")
        numeric_score = (
            isinstance(score, (int, float))
            and not isinstance(score, bool)
            and math.isfinite(float(score))
        )
        passed = (
            summary.get("idt_explicit_pass") is True
            and numeric_score
            and float(score) < 10.0
        )
        all_passed = all_passed and passed and bool(fragment["product_length_valid"])
        row = dict(fragment)
        row.update(
            {
                "idt_status": summary.get("idt_status", "unknown"),
                "idt_explicit_pass": passed,
                "idt_complexity_score": summary.get("idt_complexity_score"),
                "idt_positive_score_names_json": summary.get("idt_positive_score_names_json", "[]"),
                "idt_violation_names_json": summary.get("idt_violation_names_json", "[]"),
                "idt_response_sha256": summary.get("idt_response_sha256", ""),
                "idt_scored_sequence_sha256": observed_sha,
                "idt_score_complete": summary.get("idt_score_complete", False),
                "idt_invalid_score_names_json": summary.get(
                    "idt_invalid_score_names_json", "[]"
                ),
            }
        )
        scored.append(row)
        audit.append(
            {
                "fragment_id": fragment["fragment_id"],
                "request_length_bp": len(sequence),
                "dna_sha256": expected_sha,
                "idt_status": summary.get("idt_status", "unknown"),
                "idt_explicit_pass": passed,
                "idt_complexity_score": summary.get("idt_complexity_score"),
                "idt_rule_details_json": summary.get("idt_rule_details_json", "[]"),
                "idt_positive_score_names_json": summary.get("idt_positive_score_names_json", "[]"),
                "idt_violation_names_json": summary.get("idt_violation_names_json", "[]"),
                "response_sha256": summary.get("idt_response_sha256", ""),
                "idt_error_type": summary.get("idt_error_type", ""),
                "idt_score_complete": summary.get("idt_score_complete", False),
                "idt_invalid_score_names_json": summary.get(
                    "idt_invalid_score_names_json", "[]"
                ),
            }
        )
    return all_passed, scored, audit


def _plasmid_fasta_path(plasmid: str) -> Path:
    name = "pET-28a(+)" if plasmid == "pET-28a(+)_start_codon" else plasmid
    return ProjectPaths.discover().root / "data" / "reference_input" / "plasmids" / f"{name}.fa"


def simulate_vector_assembly(plasmid: str, cds: str, target_protein: str) -> dict[str, object]:
    """Replace the shortest circular MCS arc and verify the exact CDS feature."""
    fasta_path = _plasmid_fasta_path(plasmid)
    record = next(SeqIO.parse(fasta_path, "fasta"))
    vector = str(record.seq).upper()
    left_enzyme, right_enzyme = PLASMID_MCS[plasmid]
    left_site = str(getattr(Restriction, left_enzyme).site)
    right_site = str(getattr(Restriction, right_enzyme).site)
    left_positions = [index for index in range(len(vector)) if vector.startswith(left_site, index)]
    right_positions = [index for index in range(len(vector)) if vector.startswith(right_site, index)]
    if len(left_positions) != 1 or len(right_positions) != 1:
        raise ValueError("Maintained MCS enzymes must each have one exact plasmid recognition site")
    left_start = left_positions[0]
    right_start = right_positions[0]
    left_end = left_start + len(left_site)
    right_end = right_start + len(right_site)
    forward_gap = (right_start - left_end) % len(vector)
    reverse_gap = (left_start - right_end) % len(vector)
    if forward_gap <= reverse_gap:
        final_sequence = left_site + cds + vector[right_start:] + vector[:left_start]
        cds_start = len(left_site)
        cds_sequence_in_fasta = cds
        strand = "+"
        removed_bp = forward_gap
    else:
        coding_reverse = reverse_complement(cds)
        final_sequence = right_site + coding_reverse + vector[left_start:] + vector[:right_start]
        cds_start = len(right_site)
        cds_sequence_in_fasta = coding_reverse
        strand = "-"
        removed_bp = reverse_gap
    cds_end = cds_start + len(cds)
    recovered = final_sequence[cds_start:cds_end]
    recovered_coding = recovered if strand == "+" else reverse_complement(recovered)
    if recovered != cds_sequence_in_fasta or recovered_coding != cds:
        raise AssertionError("Final-plasmid assembly did not preserve the optimized CDS")
    if translate_dna(recovered_coding) != target_protein:
        raise AssertionError("Final-plasmid CDS does not translate to the requested protein")
    return {
        "plasmid": plasmid,
        "source_fasta": str(fasta_path.relative_to(ProjectPaths.discover().root)),
        "source_vector_sha256": hashlib.sha256(vector.encode()).hexdigest(),
        "mcs_left_enzyme": left_enzyme,
        "mcs_right_enzyme": right_enzyme,
        "mcs_left_site": left_site,
        "mcs_right_site": right_site,
        "removed_mcs_arc_bp": removed_bp,
        "insert_orientation_in_fasta": strand,
        "cds_start_1based": cds_start + 1,
        "cds_end_1based": cds_end,
        "cds_strand": strand,
        "final_plasmid_sequence": final_sequence,
        "final_plasmid_length_bp": len(final_sequence),
        "final_plasmid_sha256": hashlib.sha256(final_sequence.encode()).hexdigest(),
        "cds_sequence_sha256": hashlib.sha256(cds.encode()).hexdigest(),
        "assembly_sequence_exact": True,
        "translation_exact": True,
    }


def circular_diagnostic_digest(
    sequence: str,
    enzyme_sites: Mapping[str, str],
) -> dict[str, list[int]]:
    """Return deterministic circular fragment sizes for sequence-level QC."""
    result: dict[str, list[int]] = {}
    for enzyme, raw_site in enzyme_sites.items():
        site = str(raw_site).upper()
        if not site or set(site) - set("ACGT"):
            result[str(enzyme)] = []
            continue
        reverse = reverse_complement(site)
        positions: set[int] = set()
        for pattern in {site, reverse}:
            start = sequence.find(pattern)
            while start >= 0:
                positions.add(start)
                start = sequence.find(pattern, start + 1)
        ordered = sorted(positions)
        if not ordered:
            result[str(enzyme)] = []
        elif len(ordered) == 1:
            result[str(enzyme)] = [len(sequence)]
        else:
            result[str(enzyme)] = sorted(
                [
                    ordered[index + 1] - ordered[index]
                    for index in range(len(ordered) - 1)
                ]
                + [len(sequence) - ordered[-1] + ordered[0]]
            )
    return result


def _cloning_steps(
    request: DesignRequest,
    candidate: Mapping[str, object],
    fragments: list[dict[str, object]],
    *,
    optimized: bool,
) -> list[dict[str, object]]:
    left_mcs, right_mcs = PLASMID_MCS[request.plasmid]
    if not optimized:
        return [
            {
                "step": 1,
                "stage": "topology_draft",
                "action": "Confirm the repeat boundary and selected HURDLER route before sequence ordering.",
                "status": "not_orderable_not_for_purchase",
            },
            {
                "step": 2,
                "stage": "planned_vector_cloning",
                "action": f"Planned directional insertion into {request.plasmid} through {left_mcs}/{right_mcs}.",
                "status": "not_simulated_without_optimized_dna",
            },
            {
                "step": 3,
                "stage": "planned_hurdler_cycle",
                "action": (
                    f"Candidate cycle uses Site I {candidate['site_i_enzyme']}, "
                    f"Site II {candidate['site_ii_enzyme']}, and Site III {candidate['site_iii_enzyme']}."
                ),
                "status": "design_level_only",
            },
        ]
    steps: list[dict[str, object]] = []
    if len(fragments) > 1:
        steps.append(
            {
                "step": len(steps) + 1,
                "stage": "insert_assembly",
                "action": f"Assemble {len(fragments)} sequence-verified overlapping purchase fragments.",
                "expected_product": "optimized full-length CDS",
                "orientation": "5prime_to_3prime_coding",
            }
        )
    steps.extend(
        [
            {
                "step": len(steps) + 1,
                "stage": "vector_digest",
                "action": f"Open {request.plasmid} at its maintained {left_mcs}/{right_mcs} MCS pair.",
                "expected_product": "directional linearized backbone",
                "orientation": "defined_by_nonidentical_MCS_ends",
            },
            {
                "step": len(steps) + 2,
                "stage": "insert_ligation",
                "action": "Ligate the exact reconstructed CDS between the maintained MCS ends.",
                "expected_product": "final expression plasmid",
                "orientation": "verified_by_sequence_simulation",
            },
            {
                "step": len(steps) + 3,
                "stage": "hurdler_expansion_route",
                "action": (
                    f"For subsequent repeat expansion, use {candidate['site_i_enzyme']} (Site I), "
                    f"{candidate['site_ii_enzyme']} (Site II), and {candidate['site_iii_enzyme']} "
                    "(Site III); the Site-II recognition sequence is silent-mutated after joining."
                ),
                "expected_product": "directional scar-compatible expanded repeat cassette",
                "orientation": str(candidate.get("direction", "")),
            },
            {
                "step": len(steps) + 4,
                "stage": "diagnostic_digest",
                "action": (
                    f"Verify the MCS pair and the selected Site-I/Site-II internal-site pattern; "
                    f"Site III {candidate['site_iii_enzyme']} is disposable donor-adapter sequence."
                ),
                "expected_product": "sequence-consistent diagnostic fragment pattern",
                "orientation": "not_applicable",
            },
        ]
    )
    return steps


def _optimize_candidates(
    request: DesignRequest,
    boundary: ConfirmedBoundary,
    target_protein: str,
    candidates: list[dict[str, object]],
    scorer: ComplexityScorer,
) -> tuple[
    dict[str, object] | None,
    dict[str, object],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    root = ProjectPaths.discover().root
    codon_weights = load_codon_weights(root / "data" / "reference_output" / "codon_usage.csv")
    recognition_sites = load_restriction_sites(
        root / "data" / "reference_output" / "restriction_enzyme.csv"
    )
    if request.selected_candidate_id:
        selected = [row for row in candidates if row["candidate_id"] == request.selected_candidate_id]
        if not selected:
            raise ValueError("selected_candidate_id is not present in the allowed candidate table")
        candidates = selected
    attempts: list[dict[str, object]] = []
    aggregate_audit: list[dict[str, object]] = []
    last_weights = dict(GA_SCORE_PROFILE)
    failure_reasons: list[str] = []
    for candidate in candidates:
        try:
            dna, locked_positions, site_limits, absolute_positions = _locked_construct(
                target_protein, boundary, candidate, codon_weights
            )
        except Exception as exc:
            failure_reasons.append(f"{candidate['candidate_id']}: {type(exc).__name__}: {exc}")
            continue
        profile = dict(GA_SCORE_PROFILE)
        for generations in GENERATION_SCHEDULE:
            refined, metrics = genetic_refine_dna(
                dna,
                locked_positions=locked_positions,
                selected_site_limits=site_limits,
                recognition_sites=recognition_sites,
                codon_weights=codon_weights,
                seed=int(request.seed) + int(candidate["rank"]) * 1000 + generations,
                population_size=int(request.population_size),
                generations=generations,
                score_profile=profile,
            )
            if translate_dna(refined) != target_protein:
                raise AssertionError("GA changed the full target protein")
            local_passed = bool(metrics["ga_local_constraints_passed"])
            scored_fragments: list[dict[str, object]] = []
            attempt_audit: list[dict[str, object]] = []
            idt_passed = False
            adjustments: list[dict[str, Any]] = []
            if local_passed:
                fragments = plan_purchase_fragments(
                    request.sequence_id,
                    refined,
                    request.plasmid,
                    max_purchase_bp=int(request.max_purchase_bp),
                )
                idt_passed, scored_fragments, attempt_audit = _score_fragments(fragments, scorer)
                for record in attempt_audit:
                    aggregate_audit.append(
                        {
                            **record,
                            "candidate_id": candidate["candidate_id"],
                            "ga_generations": generations,
                        }
                    )
                if not idt_passed:
                    for fragment, audit_row in zip(scored_fragments, attempt_audit):
                        summary = {
                            "idt_rule_details_json": audit_row.get("idt_rule_details_json", "[]"),
                            "idt_positive_score_names_json": audit_row.get(
                                "idt_positive_score_names_json", "[]"
                            ),
                        }
                        profile, changes = adjust_ga_score_profile_from_idt(profile, summary)
                        for change in changes:
                            change["fragment_id"] = fragment["fragment_id"]
                        adjustments.extend(changes)
            attempt = {
                "candidate_id": candidate["candidate_id"],
                "candidate_rank": candidate["rank"],
                "ga_generations": generations,
                "ga_local_constraints_passed": local_passed,
                "ga_score": metrics["ga_score"],
                "selected_pair_re_site_excess": metrics["selected_pair_re_site_excess"],
                "repeated_re_site_excess": metrics["repeated_re_site_excess"],
                "dna_sha256": hashlib.sha256(refined.encode()).hexdigest(),
                "fragment_count": len(scored_fragments),
                "all_fragments_idt_passed": idt_passed,
                "ga_weights_before_json": metrics["ga_score_profile_json"],
                "ga_weights_after_json": json.dumps(profile, sort_keys=True),
                "idt_feedback_adjustments_json": json.dumps(adjustments, sort_keys=True),
            }
            attempts.append(attempt)
            dna = refined
            last_weights = dict(profile)
            if local_passed and idt_passed:
                selected_candidate = dict(candidate)
                selected_candidate.update(
                    {
                        "absolute_site_i_aa_position_0based": absolute_positions["site_i"],
                        "absolute_site_ii_aa_position_0based": absolute_positions["site_ii"],
                    }
                )
                result = {
                    "status": "passed",
                    "optimized_dna": refined,
                    "optimized_dna_sha256": hashlib.sha256(refined.encode()).hexdigest(),
                    "optimized_length_bp": len(refined),
                    "translation_exact": True,
                    "selected_pair_re_site_excess": metrics["selected_pair_re_site_excess"],
                    "scheme_wide_repeated_site_excess": metrics["repeated_re_site_excess"],
                    "ga_metrics": metrics,
                    "final_ga_weights": profile,
                    "attempts": attempts,
                    "termination_reason": "all_purchase_fragments_idt_score_sum_below_10",
                }
                return selected_candidate, result, scored_fragments, aggregate_audit
        failure_reasons.append(f"{candidate['candidate_id']}: no accepted construct after 100 generations")
    return None, {
        "status": "failed",
        "optimized_dna": None,
        "final_ga_weights": last_weights,
        "attempts": attempts,
        "failure_reasons": failure_reasons,
        "termination_reason": "all_allowed_candidates_exhausted",
    }, [], aggregate_audit


def _safe_request_payload(request: DesignRequest) -> dict[str, object]:
    payload = request.to_dict()
    payload["full_protein_sequence_sha256"] = hashlib.sha256(
        request.full_protein_sequence.encode()
    ).hexdigest()
    return payload


def design_construct(
    request: DesignRequest,
    *,
    index: PatternIndex | None = None,
    index_dir: str | Path | None = None,
    idt_scorer: ComplexityScorer | None = None,
) -> DesignResult:
    """Run inference, mandatory confirmation, HURDLER query, and optional GA."""
    analysis = analyze_repeat_sequence(
        request.full_protein_sequence,
        sequence_id=request.sequence_id,
    )
    base = DesignResult(
        schema_version=DESIGN_SCHEMA_VERSION,
        status="needs_boundary_confirmation",
        message="Confirm or edit the 1-based inclusive repeat boundary and period before HURDLER analysis.",
        request=_safe_request_payload(request),
        boundary_analysis=analysis.to_dict(),
    )
    if request.confirmed_repeat_start is None:
        return base
    boundary = confirm_repeat_boundary(
        request.full_protein_sequence,
        start=int(request.confirmed_repeat_start),
        end=int(request.confirmed_repeat_end),
        period=int(request.confirmed_period),
        expected_token=request.confirmation_token,
    )
    target_protein = (
        boundary.n_terminal_flank
        + boundary.middle_module * int(request.target_repeat_copies)
        + boundary.c_terminal_flank
    )
    resolved_index = index or PatternIndex.load(index_dir or bundled_index_dir())
    candidates = enumerate_design_candidates(
        boundary.middle_module,
        request.plasmid,
        resolved_index,
        site_i_allowlist=request.site_i_allowlist,
        site_ii_allowlist=request.site_ii_allowlist,
        site_iii_allowlist=request.site_iii_allowlist,
    )
    for candidate in candidates:
        minimum_copies = math.ceil(
            (max(int(candidate["site_i_position"]), int(candidate["site_ii_position"])) + 3)
            / len(boundary.middle_module)
        )
        candidate["minimum_target_repeat_copies_for_locked_windows"] = minimum_copies
        candidate["requested_target_geometry_supported"] = (
            int(request.target_repeat_copies) >= minimum_copies
        )
    base.confirmed_boundary = boundary.to_dict()
    base.target_protein_sequence = target_protein
    base.candidates = candidates
    if not candidates:
        base.status = "hurdler_incompatible"
        base.message = "No frozen-index HURDLER solution passed the selected plasmid and enzyme whitelists."
        return base
    default_candidate = candidates[0]
    if request.selected_candidate_id:
        matches = [row for row in candidates if row["candidate_id"] == request.selected_candidate_id]
        if not matches:
            raise ValueError("selected_candidate_id is not present in the allowed candidate table")
        default_candidate = matches[0]
    base.selected_candidate = default_candidate
    if not request.optimize:
        base.status = "not_orderable_not_for_purchase"
        base.message = (
            "HURDLER compatibility is confirmed, but codon optimization and live IDT fragment scoring "
            "were not run. This topology draft is not a purchase design."
        )
        base.cloning_steps = _cloning_steps(
            request, default_candidate, [], optimized=False
        )
        return base
    if idt_scorer is None:
        raise RuntimeError("Codon optimization requires a configured live or explicit mock IDT scorer")
    selected, optimization, fragments, audit = _optimize_candidates(
        request, boundary, target_protein, candidates, idt_scorer
    )
    base.optimization = optimization
    base.idt_audit = audit
    if selected is None:
        base.status = "optimization_failed_not_orderable"
        base.message = "No allowed route produced selected-pair-clean, IDT-orderable purchase fragments."
        return base
    dna = str(optimization["optimized_dna"])
    final_plasmid = simulate_vector_assembly(request.plasmid, dna, target_protein)
    left_mcs, right_mcs = PLASMID_MCS[request.plasmid]
    final_plasmid["diagnostic_digest_fragments_bp"] = circular_diagnostic_digest(
        str(final_plasmid["final_plasmid_sequence"]),
        {
            left_mcs: str(getattr(Restriction, left_mcs).site),
            right_mcs: str(getattr(Restriction, right_mcs).site),
            str(selected["site_i_enzyme"]): str(selected.get("site_i_recognition_site") or ""),
            str(selected["site_ii_enzyme"]): str(selected.get("site_ii_recognition_site") or ""),
        },
    )
    base.status = "orderable_design_complete"
    base.message = "The exact purchase fragments passed live IDT scoring and the final plasmid simulation."
    base.selected_candidate = selected
    base.purchase_fragments = fragments
    base.final_plasmid = final_plasmid
    base.cloning_steps = _cloning_steps(request, selected, fragments, optimized=True)
    base.cloning_steps[-1]["expected_fragment_sizes_bp_json"] = json.dumps(
        final_plasmid["diagnostic_digest_fragments_bp"], sort_keys=True
    )
    return base


def _write_csv(rows: list[dict[str, object]], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        destination.write_text("")
        return
    columns = sorted({key for row in rows for key in row})
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, sort_keys=True)
                    if isinstance(value, (dict, list, tuple))
                    else value
                    for key, value in row.items()
                }
            )


def _fasta_record(identifier: str, sequence: str, description: str = "") -> str:
    header = f">{identifier}" + (f" {description}" if description else "")
    lines = [sequence[index : index + 80] for index in range(0, len(sequence), 80)]
    return "\n".join([header, *lines, ""])


def _summary_row(result: DesignResult) -> dict[str, object]:
    boundary = result.confirmed_boundary or {}
    selected = result.selected_candidate or {}
    optimization = result.optimization or {}
    return {
        "schema_version": result.schema_version,
        "status": result.status,
        "sequence_id": result.request.get("sequence_id"),
        "plasmid": result.request.get("plasmid"),
        "input_protein_sha256": result.request.get("full_protein_sequence_sha256"),
        "middle_module": boundary.get("middle_module"),
        "module_length_aa": boundary.get("period"),
        "target_repeat_copies": result.request.get("target_repeat_copies"),
        "target_protein_length_aa": len(result.target_protein_sequence or ""),
        "hurdler_candidate_count": len(result.candidates),
        "selected_candidate_id": selected.get("candidate_id"),
        "site_i_enzyme": selected.get("site_i_enzyme"),
        "site_ii_enzyme": selected.get("site_ii_enzyme"),
        "site_iii_enzyme": selected.get("site_iii_enzyme"),
        "optimization_status": optimization.get("status"),
        "optimized_dna_sha256": optimization.get("optimized_dna_sha256"),
        "purchase_fragment_count": len(result.purchase_fragments),
        "orderable": result.orderable,
        "message": result.message,
    }


def write_design_outputs(result: DesignResult, output_dir: str | Path) -> dict[str, str]:
    """Write the authoritative design bundle without credential material."""
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    paths["design_summary_json"] = write_json_atomic(result.to_dict(), destination / "design_summary.json")
    paths["design_summary_csv"] = destination / "design_summary.csv"
    _write_csv([_summary_row(result)], paths["design_summary_csv"])
    paths["hurdler_candidates_csv"] = destination / "hurdler_candidates.csv"
    _write_csv(result.candidates, paths["hurdler_candidates_csv"])
    paths["cloning_steps_csv"] = destination / "cloning_steps.csv"
    _write_csv(result.cloning_steps, paths["cloning_steps_csv"])
    plan_lines = [
        "# HURDLER cloning plan",
        "",
        f"Status: `{result.status}`",
        "",
        result.message,
        "",
        "This is a sequence-level design. Reaction volumes, enzyme amounts, temperatures, and times are intentionally not prescribed.",
        "",
    ]
    for step in result.cloning_steps:
        plan_lines.append(f"{step.get('step')}. **{step.get('stage')}** — {step.get('action')}")
    paths["cloning_plan_md"] = destination / "cloning_plan.md"
    paths["cloning_plan_md"].write_text("\n".join(plan_lines) + "\n")
    if result.orderable:
        dna = str((result.optimization or {})["optimized_dna"])
        paths["optimized_construct_fasta"] = destination / "optimized_construct.fasta"
        paths["optimized_construct_fasta"].write_text(
            _fasta_record(str(result.request["sequence_id"]), dna, "optimized exact CDS")
        )
        paths["purchase_fragments_csv"] = destination / "purchase_fragments.csv"
        _write_csv(result.purchase_fragments, paths["purchase_fragments_csv"])
        paths["purchase_fragments_fasta"] = destination / "purchase_fragments.fasta"
        paths["purchase_fragments_fasta"].write_text(
            "".join(
                _fasta_record(
                    str(row["fragment_id"]),
                    str(row["purchase_sequence"]),
                    f"{row['product_type']} IDT_score={row.get('idt_complexity_score')}",
                )
                for row in result.purchase_fragments
            )
        )
        plasmid = result.final_plasmid or {}
        paths["final_plasmid_fasta"] = destination / "final_plasmid.fasta"
        paths["final_plasmid_fasta"].write_text(
            _fasta_record(
                f"{result.request['sequence_id']}_final_plasmid",
                str(plasmid["final_plasmid_sequence"]),
                f"CDS={plasmid['cds_start_1based']}..{plasmid['cds_end_1based']} strand={plasmid['cds_strand']}",
            )
        )
    else:
        paths["fragment_topology_csv"] = destination / "fragment_topology.csv"
        _write_csv(
            [
                {
                    "status": "not_orderable_not_for_purchase",
                    "plasmid": result.request.get("plasmid"),
                    "mcs_left_enzyme": PLASMID_MCS[str(result.request["plasmid"])][0],
                    "mcs_right_enzyme": PLASMID_MCS[str(result.request["plasmid"])][1],
                    "target_cds_length_bp": 3 * len(result.target_protein_sequence or ""),
                }
            ],
            paths["fragment_topology_csv"],
        )
    paths["idt_audit_jsonl"] = destination / "idt_audit.jsonl"
    if not paths["idt_audit_jsonl"].exists():
        paths["idt_audit_jsonl"].write_text(
            "".join(json.dumps(row, sort_keys=True) + "\n" for row in result.idt_audit)
        )
    file_hashes = {
        name: {"path": path.name, "sha256": sha256_file(path), "bytes": path.stat().st_size}
        for name, path in paths.items()
        if path.is_file()
    }
    manifest = {
        "schema_version": DESIGN_SCHEMA_VERSION,
        "created_at": utc_now(),
        "status": result.status,
        "rule_profile": RULE_PROFILE_NAME,
        "index_dir": "data/artifacts/legacy-optimized-v1",
        "input_protein_sha256": result.request.get("full_protein_sequence_sha256"),
        "credential_material_persisted": False,
        "files": file_hashes,
    }
    paths["run_manifest_json"] = write_json_atomic(manifest, destination / "run_manifest.json")
    result.output_files = {name: str(path) for name, path in paths.items()}
    return result.output_files
