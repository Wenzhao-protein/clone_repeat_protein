"""Deterministic codon diversification coupled to HURDLER matches."""

from __future__ import annotations

import math
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import pandas as pd

from .constants import DEFAULT_FRAGMENT_LIMITS_BP, PLASMIDS, validate_protein_sequence
from .index import PatternIndex
from .matching import enumerate_module_solutions, materialize_best_solution, query_all_plasmids

GENETIC_CODE = {
    "A": ("GCG", "GCC", "GCA", "GCT"), "C": ("TGC", "TGT"),
    "D": ("GAT", "GAC"), "E": ("GAA", "GAG"), "F": ("TTT", "TTC"),
    "G": ("GGC", "GGT", "GGA", "GGG"), "H": ("CAT", "CAC"),
    "I": ("ATT", "ATC", "ATA"), "K": ("AAA", "AAG"),
    "L": ("CTG", "TTA", "TTG", "CTT", "CTC", "CTA"), "M": ("ATG",),
    "N": ("AAC", "AAT"), "P": ("CCG", "CCT", "CCA", "CCC"),
    "Q": ("CAG", "CAA"), "R": ("CGC", "CGT", "CGG", "CGA", "AGA", "AGG"),
    "S": ("TCG", "TCA", "AGC", "TCT", "TCC", "AGT"),
    "T": ("ACC", "ACG", "ACT", "ACA"), "V": ("GTG", "GTT", "GTC", "GTA"),
    "W": ("TGG",), "Y": ("TAT", "TAC"),
}
CODON_TO_AA = {codon: aa for aa, codons in GENETIC_CODE.items() for codon in codons}


def translate_dna(dna: str) -> str:
    if len(dna) % 3:
        raise ValueError("DNA length must be divisible by three")
    try:
        return "".join(CODON_TO_AA[dna[index : index + 3]] for index in range(0, len(dna), 3))
    except KeyError as exc:
        raise ValueError(f"Unsupported codon: {exc.args[0]}") from exc


def repeated_nmer_count(sequence: str, n: int) -> int:
    if len(sequence) < n:
        return 0
    counts = Counter(sequence[index : index + n] for index in range(len(sequence) - n + 1))
    return sum(count - 1 for count in counts.values())


def hairpin_proxy(sequence: str, n: int = 10) -> int:
    complement = str.maketrans("ACGT", "TGCA")
    windows = {sequence[index : index + n] for index in range(max(0, len(sequence) - n + 1))}
    return sum(1 for window in windows if window.translate(complement)[::-1] in windows)


def gc_window_extrema(sequence: str, window: int = 50) -> tuple[float, float]:
    """Return the minimum and maximum GC fraction over fixed-size windows."""
    if not sequence:
        raise ValueError("DNA sequence cannot be empty")
    if len(sequence) <= window:
        value = (sequence.count("G") + sequence.count("C")) / len(sequence)
        return value, value
    fractions = [
        (chunk.count("G") + chunk.count("C")) / window
        for start in range(len(sequence) - window + 1)
        for chunk in (sequence[start : start + window],)
    ]
    return min(fractions), max(fractions)


def load_codon_weights(path: str | Path) -> dict[str, float]:
    """Load Kazusa 83333 frequencies as within-amino-acid relative weights."""
    frame = pd.read_csv(path)
    if not {"codon", "frequency"}.issubset(frame.columns):
        raise ValueError("Codon table must contain codon and frequency columns")
    frequency = {str(row.codon).upper(): float(row.frequency) for row in frame.itertuples(index=False)}
    weights: dict[str, float] = {}
    for aa, codons in GENETIC_CODE.items():
        maximum = max(frequency[codon] for codon in codons)
        for codon in codons:
            weights[codon] = frequency[codon] / maximum
    return weights


def codon_adaptation_index(dna: str, weights: dict[str, float]) -> float:
    codons = [dna[start : start + 3] for start in range(0, len(dna), 3)]
    if not codons:
        return 0.0
    return math.exp(sum(math.log(max(weights[codon], 1e-12)) for codon in codons) / len(codons))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def recognition_site_count(dna: str, site: str) -> int:
    site = site.upper()
    if not site:
        return 0
    return len(_site_occurrences(dna, site))


def _site_occurrences(dna: str, site: str) -> list[tuple[int, int]]:
    """Return overlapping forward/reverse occurrences as half-open intervals."""
    site = site.upper()
    if not site:
        return []
    reverse = reverse_complement(site)
    patterns = (site,) if reverse == site else (site, reverse)
    occurrences: list[tuple[int, int]] = []
    for pattern in patterns:
        start = dna.find(pattern)
        while start >= 0:
            occurrences.append((start, start + len(pattern)))
            # Advance one base so self-overlapping recognition sites are kept.
            start = dna.find(pattern, start + 1)
    return occurrences


def _site_limit_excess(dna: str, site_limits: dict[str, int]) -> int:
    return sum(
        max(0, recognition_site_count(dna, site) - limit)
        for site, limit in site_limits.items()
    )


def _repair_site_limits(
    dna: str,
    protein: str,
    locked_positions: set[int],
    site_limits: dict[str, int],
    codon_weights: dict[str, float] | None,
) -> str:
    """Remove look-ahead restriction sites with translation-preserving edits."""
    maximum_steps = max(20, _site_limit_excess(dna, site_limits) * 10)
    for _step in range(maximum_steps):
        current_excess = _site_limit_excess(dna, site_limits)
        if current_excess == 0:
            return dna
        affected_positions: set[int] = set()
        for site, limit in site_limits.items():
            occurrences = _site_occurrences(dna, site)
            if len(occurrences) <= limit:
                continue
            for start, end in occurrences:
                affected_positions.update(range(start // 3, (end - 1) // 3 + 1))
        best: tuple[tuple[float, ...], str] | None = None
        for position in sorted(affected_positions):
            if position in locked_positions:
                continue
            old_codon = dna[position * 3 : position * 3 + 3]
            for preference, codon in enumerate(GENETIC_CODE[protein[position]]):
                if codon == old_codon:
                    continue
                candidate = dna[: position * 3] + codon + dna[position * 3 + 3 :]
                excess = _site_limit_excess(candidate, site_limits)
                if excess >= current_excess:
                    continue
                usage_penalty = (
                    -math.log(max(codon_weights.get(codon, 1.0), 1e-12))
                    if codon_weights is not None
                    else preference * 0.05
                )
                score = (
                    float(excess),
                    float(repeated_nmer_count(candidate, 14)),
                    float(repeated_nmer_count(candidate, 13)),
                    float(repeated_nmer_count(candidate, 8)),
                    usage_penalty,
                    float(preference),
                )
                if best is None or score < best[0]:
                    best = (score, candidate)
        if best is None:
            return dna
        dna = best[1]
    return dna


def unexpected_selected_sites(dna: str, solution: dict[str, object]) -> int:
    """Count excess sites for the selected Site-I/Site-II enzyme pair only."""
    site_i = str(solution.get("site_i_recognition_site", ""))
    site_ii = str(solution.get("site_ii_recognition_site", ""))
    unexpected = max(0, recognition_site_count(dna, site_i) - 1)
    unexpected += recognition_site_count(dna, site_ii)
    return unexpected


def _new_duplicate_windows(sequence: str, n: int, appended: int = 3) -> int:
    """Count newly completed n-mers that already occurred upstream."""
    if len(sequence) < n:
        return 0
    first_start = max(0, len(sequence) - n - appended + 1)
    duplicate_count = 0
    for start in range(first_start, len(sequence) - n + 1):
        window = sequence[start : start + n]
        if window in sequence[:start]:
            duplicate_count += 1
    return duplicate_count


def _new_duplicate_windows_from_seen(
    sequence: str,
    n: int,
    earliest_end_by_window: dict[str, int],
    appended: int = 3,
) -> int:
    """Incremental equivalent of `_new_duplicate_windows` for appended DNA."""
    if len(sequence) < n:
        return 0
    first_start = max(0, len(sequence) - n - appended + 1)
    duplicate_count = 0
    for start in range(first_start, len(sequence) - n + 1):
        window = sequence[start : start + n]
        if earliest_end_by_window.get(window, len(sequence) + 1) <= start:
            duplicate_count += 1
    return duplicate_count


def _update_seen_windows(
    sequence: str,
    n: int,
    earliest_end_by_window: dict[str, int],
    appended: int = 3,
) -> None:
    if len(sequence) < n:
        return
    first_start = max(0, len(sequence) - n - appended + 1)
    for start in range(first_start, len(sequence) - n + 1):
        window = sequence[start : start + n]
        earliest_end_by_window.setdefault(window, start + n)


def _new_site_occurrence_count(sequence: str, site: str, appended: int = 3) -> int:
    """Count forward/reverse site occurrences completed by appended bases."""
    site = site.upper()
    reverse = reverse_complement(site)
    patterns = (site,) if reverse == site else (site, reverse)
    count = 0
    for pattern in patterns:
        first_start = max(0, len(sequence) - len(pattern) - appended + 1)
        for start in range(first_start, len(sequence) - len(pattern) + 1):
            count += int(sequence.startswith(pattern, start))
    return count


def diversify_codons(
    protein: str,
    locked: dict[int, str] | None = None,
    codon_weights: dict[str, float] | None = None,
    site_limits: dict[str, int] | None = None,
) -> str:
    """Greedy repeat-aware reverse translation with locked and banned sites."""
    protein = validate_protein_sequence(protein)
    locked = locked or {}
    site_limits = {site.upper(): int(limit) for site, limit in (site_limits or {}).items() if site}
    locked_codons: dict[int, str] = {}
    for start, dna9 in locked.items():
        if len(dna9) != 9:
            raise ValueError("Locked HURDLER windows must be 9 bp")
        for offset in range(3):
            locked_codons[start + offset] = dna9[offset * 3 : offset * 3 + 3]
    dna = ""
    site_counts = {site: 0 for site in site_limits}
    seen_windows: dict[int, dict[str, int]] = {n: {} for n in (8, 13, 14)}

    def update_incremental_state(sequence: str) -> None:
        for site in site_limits:
            site_counts[site] += _new_site_occurrence_count(sequence, site)
        for n, seen in seen_windows.items():
            _update_seen_windows(sequence, n, seen)

    for position, aa in enumerate(protein):
        if position in locked_codons:
            codon = locked_codons[position]
            if CODON_TO_AA.get(codon) != aa:
                raise ValueError(f"Locked codon {codon} does not encode {aa} at position {position}")
            dna += codon
            update_incremental_state(dna)
            continue
        best = None
        for preference, codon in enumerate(GENETIC_CODE[aa]):
            candidate = dna + codon
            site_penalty = sum(
                1_000_000
                * max(
                    0,
                    site_counts[site] + _new_site_occurrence_count(candidate, site) - limit,
                )
                for site, limit in site_limits.items()
            )
            repeat_penalty = (
                250 * _new_duplicate_windows_from_seen(candidate, 14, seen_windows[14])
                + 150 * _new_duplicate_windows_from_seen(candidate, 13, seen_windows[13])
                + 25 * _new_duplicate_windows_from_seen(candidate, 8, seen_windows[8])
            )
            gc_window = candidate[-60:]
            gc = (gc_window.count("G") + gc_window.count("C")) / len(gc_window)
            gc_penalty = max(0.0, abs(gc - 0.5) - 0.18) * 100
            usage_penalty = (
                -math.log(max(codon_weights.get(codon, 1.0), 1e-12))
                if codon_weights is not None
                else preference * 0.05
            )
            score = site_penalty + repeat_penalty + gc_penalty + usage_penalty
            item = (score, preference, codon)
            if best is None or item < best:
                best = item
        dna += best[2]
        update_incremental_state(dna)
    if site_limits:
        dna = _repair_site_limits(
            dna,
            protein,
            set(locked_codons),
            site_limits,
            codon_weights,
        )
    if translate_dna(dna) != protein:
        raise AssertionError("Codon diversification changed the protein sequence")
    return dna


def _candidate_rank(row: dict[str, object]) -> tuple[object, ...]:
    """Stable legacy-compatible preference order after optimizability."""
    return (
        -(float(row.get("site_i_codon_usage_freq", 0)) + float(row.get("site_ii_codon_usage_freq", 0))),
        -float(row.get("orthogonality", 0)),
        str(row.get("site_i_enzyme", "")),
        str(row.get("site_ii_enzyme", "")),
        PLASMIDS.index(str(row["plasmid"])),
        int(row.get("site_i_position", 0)),
        int(row.get("site_ii_position", 0)),
        int(row.get("candidate_pair_id", row.get("best_pair_id", 0))),
    )


_OPTIMIZATION_SIGNATURE_FIELDS = (
    "site_i_position",
    "site_ii_position",
    "site_i_9mer_bp",
    "site_ii_9mer_bp_mutated",
    "site_i_recognition_site",
    "site_ii_recognition_site",
    "site_iii_sites",
)


def _optimization_signature(row: dict[str, object]) -> tuple[object, ...]:
    """Fields that completely determine construct DNA and hard-site checks."""
    return tuple(row.get(field) for field in _OPTIMIZATION_SIGNATURE_FIELDS)


def _unique_optimization_candidates(
    ordered_candidates: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Keep the stable first representative of each equivalent DNA problem.

    Plasmid and enzyme-pair bookkeeping can create many candidates with the
    same locked codons and the same recognition-site constraints.  Their
    construct optimization is identical.  Deduplicating only after applying
    ``_candidate_rank`` preserves the exact reporting winner while avoiding
    repeated work and duplicate work split across process groups.
    """
    seen: set[tuple[object, ...]] = set()
    unique: list[dict[str, object]] = []
    for candidate in ordered_candidates:
        signature = _optimization_signature(candidate)
        if signature in seen:
            continue
        seen.add(signature)
        unique.append(candidate)
    return unique


def _construct_metrics(
    unit: str,
    copies: int,
    solution: dict[str, object],
    codon_weights: dict[str, float],
    *,
    validate_hard_constraints: bool = True,
) -> dict[str, Any]:
    """Build and validate one synonymously diversified construct candidate."""
    if copies < 1:
        raise ValueError("A construct must contain at least one module copy")
    protein = unit * copies
    locks = {
        int(solution["site_i_position"]): str(solution["site_i_9mer_bp"]),
        int(solution["site_ii_position"]): str(solution["site_ii_9mer_bp_mutated"]),
    }
    if max(locks) + 3 > len(protein):
        raise ValueError("Construct is too short to contain the locked HURDLER windows")
    site_limits = {
        str(solution.get("site_i_recognition_site", "")): 1,
        str(solution.get("site_ii_recognition_site", "")): 0,
    }
    dna = diversify_codons(protein, locks, codon_weights, site_limits)
    gc_min, gc_max = gc_window_extrema(dna, 50)
    extra_sites = unexpected_selected_sites(dna, solution)
    if validate_hard_constraints and extra_sites:
        raise ValueError(f"{extra_sites} unexpected selected-enzyme recognition sites")
    if validate_hard_constraints and (gc_min < 0.25 or gc_max > 0.75):
        raise ValueError(f"50-bp GC range {gc_min:.3f}--{gc_max:.3f} outside 0.25--0.75")
    return {
        "dna_sequence": dna,
        "dna_length": len(dna),
        "gc_fraction": (dna.count("G") + dna.count("C")) / len(dna),
        "gc_50bp_min": gc_min,
        "gc_50bp_max": gc_max,
        "cai_kazusa_83333": codon_adaptation_index(dna, codon_weights),
        "unexpected_selected_sites": extra_sites,
        "selected_pair_re_site_excess": extra_sites,
        "site_iii_re_site_count": sum(
            recognition_site_count(dna, site)
            for site in str(solution.get("site_iii_sites", "")).split(",")
            if site and site != "nan"
        ),
        "repeat_8mer_count": repeated_nmer_count(dna, 8),
        "repeat_13mer_count": repeated_nmer_count(dna, 13),
        "repeat_14mer_count": repeated_nmer_count(dna, 14),
        "hairpin_10mer_proxy": hairpin_proxy(dna, 10),
    }


def _locked_site_excess_lower_bound(solution: dict[str, object]) -> int:
    """Return unavoidable selected-site excess wholly inside locked codons."""
    required = {
        "site_i_position",
        "site_ii_position",
        "site_i_9mer_bp",
        "site_ii_9mer_bp_mutated",
    }
    if not required.issubset(solution):
        return 0
    locked_codons: dict[int, str] = {}
    for start_key, sequence_key in (
        ("site_i_position", "site_i_9mer_bp"),
        ("site_ii_position", "site_ii_9mer_bp_mutated"),
    ):
        start = int(solution[start_key])
        dna9 = str(solution[sequence_key])
        for offset in range(3):
            codon = dna9[offset * 3 : offset * 3 + 3]
            old = locked_codons.setdefault(start + offset, codon)
            if old != codon:
                return 1
    segments: list[str] = []
    current: list[str] = []
    previous: int | None = None
    for position, codon in sorted(locked_codons.items()):
        if previous is not None and position != previous + 1:
            segments.append("".join(current))
            current = []
        current.append(codon)
        previous = position
    if current:
        segments.append("".join(current))
    site_limits = {
        str(solution.get("site_i_recognition_site", "")): 1,
        str(solution.get("site_ii_recognition_site", "")): 0,
    }
    return sum(
        max(0, sum(recognition_site_count(segment, site) for segment in segments) - limit)
        for site, limit in site_limits.items()
        if site
    )


def _maximum_verified_construct(
    unit: str,
    mathematical_max: int,
    candidates: list[dict[str, object]],
    codon_weights: dict[str, float],
    workers: int = 1,
) -> tuple[int, dict[str, object] | None, dict[str, Any] | None, list[str]]:
    """Find the largest validated construct, then apply stable solution ranking.

    For a fixed candidate, generated shorter constructs are prefixes of longer
    constructs.  The hard checks (extra restriction sites and GC extrema) are
    therefore monotone, so binary search avoids a potentially thousand-step
    linear descent for 1-AA units.
    """
    ordered_candidates = _unique_optimization_candidates(
        sorted(candidates, key=_candidate_rank)
    )
    if workers > 1 and len(ordered_candidates) > 1:
        active_workers = min(workers, len(ordered_candidates))
        group_size = math.ceil(len(ordered_candidates) / active_workers)
        groups = [
            ordered_candidates[start : start + group_size]
            for start in range(0, len(ordered_candidates), group_size)
        ]
        arguments = [
            (unit, mathematical_max, group, codon_weights)
            for group in groups
        ]
        with ProcessPoolExecutor(max_workers=active_workers) as executor:
            group_results = list(executor.map(_maximum_verified_construct_group, arguments))
        best_copies = max(result[0] for result in group_results)
        finalists = [
            result for result in group_results if result[0] == best_copies and result[1] is not None
        ]
        errors = [error for result in group_results for error in result[3]][:50]
        if not finalists:
            return 0, None, None, errors
        finalists.sort(key=lambda result: _candidate_rank(result[1]))
        return finalists[0][0], finalists[0][1], finalists[0][2], errors

    best_copies = 0
    passing: list[tuple[dict[str, object], dict[str, Any]]] = []
    errors: list[str] = []
    evaluation_cache: dict[
        tuple[tuple[object, ...], int],
        tuple[dict[str, Any] | None, str | None],
    ] = {}
    for candidate in ordered_candidates:
        if _locked_site_excess_lower_bound(candidate):
            if len(errors) < 50:
                errors.append(
                    f"{candidate.get('plasmid')}:{candidate.get('site_i_enzyme')}/"
                    f"{candidate.get('site_ii_enzyme')}: selected RE site is unavoidable inside locked windows"
                )
            continue
        minimum = max(
            1,
            math.ceil((max(int(candidate["site_i_position"]), int(candidate["site_ii_position"])) + 3) / len(unit)),
        )
        if minimum > mathematical_max or mathematical_max < best_copies:
            continue
        optimization_signature = _optimization_signature(candidate)

        def evaluate(copies: int) -> tuple[dict[str, Any] | None, str | None]:
            cache_key = (optimization_signature, copies)
            if cache_key not in evaluation_cache:
                try:
                    evaluation_cache[cache_key] = (
                        _construct_metrics(unit, copies, candidate, codon_weights),
                        None,
                    )
                except Exception as exc:  # retained as an auditable per-candidate reason
                    evaluation_cache[cache_key] = (None, f"{type(exc).__name__}: {exc}")
            return evaluation_cache[cache_key]

        maximum_metrics, maximum_error = evaluate(mathematical_max)
        if maximum_metrics is not None:
            candidate_copies = mathematical_max
            metrics = maximum_metrics
        else:
            low = minimum
            high = mathematical_max - 1
            candidate_copies = 0
            metrics = None
            while low <= high:
                middle = (low + high) // 2
                middle_metrics, _ = evaluate(middle)
                if middle_metrics is not None:
                    candidate_copies = middle
                    metrics = middle_metrics
                    low = middle + 1
                else:
                    high = middle - 1
            if maximum_error and len(errors) < 50:
                errors.append(
                    f"{candidate.get('plasmid')}:{candidate.get('site_i_enzyme')}/"
                    f"{candidate.get('site_ii_enzyme')}: {maximum_error}"
                )
        if candidate_copies > best_copies:
            best_copies = candidate_copies
            passing = [(candidate, metrics)] if metrics is not None else []
        elif candidate_copies == best_copies and candidate_copies > 0 and metrics is not None:
            passing.append((candidate, metrics))
        if best_copies == mathematical_max:
            # Candidates are pre-sorted: the first full-length pass is the
            # stable winner, and no later candidate can improve copy count.
            break
    if not passing:
        return 0, None, None, errors
    passing.sort(key=lambda item: _candidate_rank(item[0]))
    solution, metrics = passing[0]
    return best_copies, solution, metrics, errors


def _maximum_verified_construct_group(
    arguments: tuple[str, int, list[dict[str, object]], dict[str, float]],
) -> tuple[int, dict[str, object] | None, dict[str, Any] | None, list[str]]:
    unit, mathematical_max, candidates, codon_weights = arguments
    return _maximum_verified_construct(
        unit,
        mathematical_max,
        candidates,
        codon_weights,
        workers=1,
    )


def optimize_module_catalog(
    catalog_path: str | Path,
    index_dir: str | Path,
    output_dir: str | Path,
    *,
    fragment_limits: tuple[int, ...] = DEFAULT_FRAGMENT_LIMITS_BP,
    external_deduction_bp: int = 0,
    codon_usage_path: str | Path | None = None,
    shard_index: int = 0,
    shard_count: int = 1,
    limit: int | None = None,
    workers: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    catalog = pd.read_parquet(catalog_path) if Path(catalog_path).suffix == ".parquet" else pd.read_csv(catalog_path)
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("shard_index must satisfy 0 <= shard_index < shard_count")
    if workers < 1:
        raise ValueError("workers must be positive")
    catalog = catalog.iloc[shard_index::shard_count].copy()
    if limit is not None:
        if limit < 1:
            raise ValueError("limit must be positive")
        catalog = catalog.head(limit).copy()
    index = PatternIndex.load(index_dir)
    if external_deduction_bp < 0:
        raise ValueError("external_deduction_bp must be non-negative")
    codon_weights = load_codon_weights(codon_usage_path) if codon_usage_path is not None else {
        codon: 1.0 for codon in CODON_TO_AA
    }
    result_rows: list[dict[str, object]] = []
    construct_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str]] = []

    for module in catalog.itertuples(index=False):
        unit = validate_protein_sequence(module.unit_sequence)
        matches = query_all_plasmids(unit, index, expand_short=True)
        solutions = [materialize_best_solution(match, index) for match in matches]
        for solution in solutions:
            solution["module_id"] = module.module_id
            solution["collection"] = module.collection
            solution["family"] = module.family
            solution["evidence_tier"] = module.evidence_tier
            solution["in_designed_primary100"] = bool(
                getattr(module, "in_designed_primary100", False)
            )
            result_rows.append(solution)
        all_candidates = enumerate_module_solutions(unit, index, expand_short=True)
        for candidate in all_candidates:
            candidate["module_id"] = module.module_id
            candidate["collection"] = module.collection
            candidate["family"] = module.family
            candidate["evidence_tier"] = module.evidence_tier
            candidate["in_designed_primary100"] = bool(
                getattr(module, "in_designed_primary100", False)
            )
            candidate_rows.append(candidate)
        successful = all_candidates
        ranked = sorted(successful, key=_candidate_rank)

        for cap in fragment_limits:
            available_bp = max(0, cap - external_deduction_bp)
            mathematical_max = available_bp // (3 * len(unit))
            row: dict[str, object] = {
                "module_id": module.module_id,
                "collection": module.collection,
                "family": module.family,
                "evidence_tier": module.evidence_tier,
                "in_designed_primary100": bool(getattr(module, "in_designed_primary100", False)),
                "unit_sequence": unit,
                "unit_length": len(unit),
                "fragment_limit_bp": cap,
                "external_deduction_bp": external_deduction_bp,
                "available_coding_bp": available_bp,
                "mathematical_max_copies": mathematical_max,
                "verified_max_copies": 0,
                "optimization_status": "no_hurdler_solution" if not ranked else "pending",
                "failure_reason": (
                    "No compatible legacy-optimized-v1 HURDLER solution across the eight plasmids"
                    if not ranked
                    else ""
                ),
            }
            if ranked and mathematical_max:
                verified, best, metrics, errors = _maximum_verified_construct(
                    unit, mathematical_max, ranked, codon_weights, workers=workers
                )
                if best is not None and metrics is not None:
                    row.update(
                        {
                            "verified_max_copies": verified,
                            "optimization_status": "passed" if verified == mathematical_max else "passed_reduced",
                            "plasmid": best["plasmid"],
                            "direction": best["direction"],
                            "site_i_position": best["site_i_position"],
                            "site_ii_position": best["site_ii_position"],
                            "site_i_enzyme": best["site_i_enzyme"],
                            "site_ii_enzyme": best["site_ii_enzyme"],
                            "site_iii_enzymes": best["site_iii_enzymes"],
                            **metrics,
                        }
                    )
                    row["failure_reason"] = ""
                    fasta_records.append((f"{module.module_id}|cap={cap}|copies={verified}", str(metrics["dna_sequence"])))
                else:
                    row["optimization_status"] = "failed"
                    row["failure_reason"] = " | ".join(errors) if errors else "No candidate passed construct validation"
            construct_rows.append(row)

    results = pd.DataFrame(result_rows)
    constructs = pd.DataFrame(construct_rows)
    candidates = pd.DataFrame(candidate_rows)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    results.to_parquet(destination / "module_hurdler_results.parquet", index=False)
    results.to_csv(destination / "module_hurdler_results.csv", index=False)
    candidates.to_parquet(destination / "module_hurdler_candidates.parquet", index=False)
    constructs.to_parquet(destination / "optimized_constructs.parquet", index=False)
    constructs.to_csv(destination / "optimized_constructs.csv", index=False)
    with (destination / "optimized_constructs.fasta").open("w") as handle:
        for name, dna in fasta_records:
            handle.write(f">{name}\n")
            for start in range(0, len(dna), 80):
                handle.write(dna[start : start + 80] + "\n")
    return results, constructs
