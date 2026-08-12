"""Plan exact-DNA HURDLER assembly through active and latent RE sites.

The protein workflow may redesign synonymous codons.  This module never
changes the requested target DNA.  It searches exact recognition sites and
one-base latent sites, constructs disposable Type-IIS donor adapters, and
audits a linear digest/ligation route whose final insert is byte-for-byte the
requested target sequence.
"""

from __future__ import annotations

import bisect
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator

import numpy as np
import pandas as pd
from Bio import Restriction
from Bio.Seq import Seq
from Bio import SeqIO

from .constants import PLASMIDS
from .idt import IDT_SCORE_POLICY, IDTComplexityScorer
from .io import sha256_file, utc_now, write_json_atomic


DNA_ASSEMBLY_VERSION = "arbitrary-dna-active-latent-v1"
DNA_COMPLETE_ROUTE_VERSION = "arbitrary-dna-complete-route-v2"
DEFAULT_CLAMP = "GCGCGC"
DEFAULT_OLIGO_MIN_BP = 20
DEFAULT_OLIGO_MAX_BP = 200
DEFAULT_GBLOCK_MIN_BP = 125
DEFAULT_GBLOCK_MAX_BP = 3000
PRIMER_PAIR_CORE_THRESHOLD_BP = 90
DEFAULT_TOP_ROUTES = 10
DNA_ALPHABET = frozenset("ACGT")


@dataclass(frozen=True)
class EnzymeGeometry:
    enzyme: str
    canonical_enzyme: str
    recognition_site: str
    ovhg: int
    ovhgseq: str
    top_cut_offset: int
    bottom_cut_offset: int
    elucidate: str
    is_type_iis: bool
    site_i_eligible: bool
    site_ii_eligible: bool
    site_iii_eligible: bool
    methylation_compatible: bool
    ligation_ok: bool
    no_star_activity: bool

    @property
    def overhang_length(self) -> int:
        return abs(self.ovhg)


@dataclass(frozen=True)
class RestrictionHit:
    target_id: str
    enzyme: str
    canonical_enzyme: str
    site: str
    orientation: str
    start: int
    end: int
    state: str
    observed: str
    mismatch_index: int | None
    mismatch_position: int | None
    active_base: str | None
    target_base: str | None
    top_cut: int
    bottom_cut: int
    ovhg: int
    ovhgseq: str
    site_i_eligible: bool
    site_ii_eligible: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class TargetRecord:
    target_id: str
    sequence: str
    cohort: str = "user_input"
    architecture: str = "unspecified"
    source_url: str = ""
    source_accession: str = ""
    unit_sequence: str = ""
    copy_count: int | None = None
    notes: str = ""
    source_database: str = ""
    element_id: str = ""
    synthetic_unit_length_bp: int | None = None
    gc_target: float | None = None
    synthetic_replicate: int | None = None

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["length_bp"] = len(self.sequence)
        payload["sequence_sha256"] = hashlib.sha256(self.sequence.encode()).hexdigest()
        return payload


def validate_dna(sequence: str, *, allow_empty: bool = False) -> str:
    normalized = "".join(str(sequence).split()).upper().replace("U", "T")
    if not normalized and not allow_empty:
        raise ValueError("DNA sequence is empty")
    invalid = sorted(set(normalized) - DNA_ALPHABET)
    if invalid:
        raise ValueError(f"Unsupported DNA symbols: {''.join(invalid)}")
    return normalized


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


def _bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _canonical_name(row: pd.Series) -> str:
    aliases = str(row.get("isoschizomers", "") or "")
    names = sorted(part for part in aliases.split("__") if part)
    return names[0] if names else str(row["enzyme"])


def _cut_offsets(enzyme_name: str, orientation: str = "+") -> tuple[int, int, str]:
    enzyme = getattr(Restriction, enzyme_name)
    elucidate = str(enzyme.elucidate())
    clean_index = 0
    top_cut = None
    bottom_cut = None
    for character in elucidate:
        if character == "^":
            top_cut = clean_index
        elif character == "_":
            bottom_cut = clean_index
        else:
            clean_index += 1
    if top_cut is None or bottom_cut is None:
        raise ValueError(f"Could not parse cut markers for {enzyme_name}: {elucidate}")
    if orientation == "-":
        top_cut, bottom_cut = clean_index - bottom_cut, clean_index - top_cut
    return int(top_cut), int(bottom_cut), elucidate


def load_enzyme_catalog(
    reference_dir: str | Path,
    *,
    artifact_dir: str | Path | None = None,
) -> tuple[dict[str, EnzymeGeometry], pd.DataFrame]:
    """Load maintained Site-I/II/III pools and plasmid compatibility."""
    reference_root = Path(reference_dir)
    artifact_root = Path(artifact_dir) if artifact_dir else reference_root.parents[1] / "output"
    base = pd.read_csv(reference_root / "restriction_enzyme.csv").rename(columns={"name": "enzyme"})
    site_i_path = artifact_root / "selected_site_i_enzymes.csv"
    site_ii_path = artifact_root / "selected_site_ii_enzymes.csv"
    site_iii_path = artifact_root / "selected_site_iii_enzymes.csv"
    if site_i_path.is_file() and site_ii_path.is_file():
        site_i_names = set(pd.read_csv(site_i_path)["enzyme"].astype(str))
        site_ii_names = set(pd.read_csv(site_ii_path)["enzyme"].astype(str))
    else:
        available = pd.read_csv(reference_root / "available_restriction_enzyme.csv")
        site_i_names = site_ii_names = set(available["name"].astype(str))
    type_rows: list[pd.DataFrame] = []
    if site_iii_path.is_file():
        third = pd.read_csv(site_iii_path).copy()
        third["ovhgseq"] = ""
        third["elucidate"] = ""
        third["isoschizomers"] = third["enzyme"]
        type_rows.append(third)
    combined = pd.concat([base, *type_rows], ignore_index=True, sort=False)
    combined = combined.drop_duplicates(subset=["enzyme"], keep="last")
    geometries: dict[str, EnzymeGeometry] = {}
    for _, row in combined.iterrows():
        name = str(row["enzyme"])
        try:
            enzyme = getattr(Restriction, name)
            top_cut, bottom_cut, elucidate = _cut_offsets(name)
            site = str(enzyme.site)
            ovhg = int(enzyme.ovhg)
            ovhgseq = str(enzyme.ovhgseq)
            is_type_iis = bool(
                int(enzyme.fst5) < 0
                or int(enzyme.fst5) > len(site)
                or int(enzyme.fst3) < -len(site)
                or int(enzyme.fst3) > len(site)
            )
        except (AttributeError, TypeError, ValueError):
            continue
        if not site or set(site) - DNA_ALPHABET or ovhg == 0:
            continue
        geometries[name] = EnzymeGeometry(
            enzyme=name,
            canonical_enzyme=_canonical_name(row),
            recognition_site=site,
            ovhg=ovhg,
            ovhgseq=ovhgseq,
            top_cut_offset=top_cut,
            bottom_cut_offset=bottom_cut,
            elucidate=elucidate,
            is_type_iis=is_type_iis,
            site_i_eligible=name in site_i_names and not is_type_iis,
            site_ii_eligible=name in site_ii_names and not is_type_iis,
            site_iii_eligible=is_type_iis,
            methylation_compatible=_bool(row.get("methylation_compatible", True)),
            ligation_ok=_bool(row.get("ligation_ok", True)),
            no_star_activity=_bool(row.get("no_star_activity", True)),
        )
    plasmids = pd.read_csv(reference_root / "plasmid_digest_check.csv", index_col=0)
    plasmids.index = plasmids.index.astype(str)
    for column in plasmids.columns:
        plasmids[column] = plasmids[column].map(_bool)
    return geometries, plasmids


def _oriented_geometry(geometry: EnzymeGeometry, orientation: str) -> tuple[str, int, int, str]:
    if orientation == "+":
        return (
            geometry.recognition_site,
            geometry.top_cut_offset,
            geometry.bottom_cut_offset,
            geometry.ovhgseq,
        )
    site = reverse_complement(geometry.recognition_site)
    top, bottom, _ = _cut_offsets(geometry.enzyme, "-")
    return site, top, bottom, reverse_complement(geometry.ovhgseq)


def scan_re_sites(
    target_id: str,
    sequence: str,
    geometries: dict[str, EnzymeGeometry],
    *,
    include_type_iis: bool = False,
) -> list[RestrictionHit]:
    """Return all active and exact one-base latent hits on both orientations."""
    dna = validate_dna(sequence)
    hits: dict[tuple[Any, ...], RestrictionHit] = {}
    for geometry in geometries.values():
        if geometry.is_type_iis != include_type_iis:
            continue
        if not include_type_iis and not (geometry.site_i_eligible or geometry.site_ii_eligible):
            continue
        for orientation in ("+", "-"):
            site, top_offset, bottom_offset, ovhgseq = _oriented_geometry(geometry, orientation)
            if orientation == "-" and site == geometry.recognition_site:
                continue
            width = len(site)
            for start in range(0, len(dna) - width + 1):
                observed = dna[start : start + width]
                differences = [index for index, (left, right) in enumerate(zip(observed, site)) if left != right]
                if len(differences) > 1:
                    continue
                state = "active" if not differences else "latent"
                mismatch_index = differences[0] if differences else None
                mismatch_position = start + mismatch_index if mismatch_index is not None else None
                hit = RestrictionHit(
                    target_id=target_id,
                    enzyme=geometry.enzyme,
                    canonical_enzyme=geometry.canonical_enzyme,
                    site=site,
                    orientation=orientation,
                    start=start,
                    end=start + width,
                    state=state,
                    observed=observed,
                    mismatch_index=mismatch_index,
                    mismatch_position=mismatch_position,
                    active_base=site[mismatch_index] if mismatch_index is not None else None,
                    target_base=observed[mismatch_index] if mismatch_index is not None else None,
                    top_cut=start + top_offset,
                    bottom_cut=start + bottom_offset,
                    ovhg=geometry.ovhg,
                    ovhgseq=ovhgseq,
                    site_i_eligible=geometry.site_i_eligible,
                    site_ii_eligible=geometry.site_ii_eligible,
                )
                key = (
                    geometry.canonical_enzyme,
                    orientation,
                    start,
                    state,
                    mismatch_position,
                    geometry.ovhg,
                )
                hits.setdefault(key, hit)
    return sorted(
        hits.values(),
        key=lambda item: (item.start, item.end, item.canonical_enzyme, item.orientation, item.state),
    )


def _nearest_hits(sorted_hits: list[RestrictionHit], position: int, neighbours: int = 3) -> Iterator[RestrictionHit]:
    starts = [item.start for item in sorted_hits]
    center = bisect.bisect_left(starts, position)
    lower = max(0, center - neighbours)
    upper = min(len(sorted_hits), center + neighbours + 1)
    yield from sorted_hits[lower:upper]


def _plasmid_options(
    enzyme_i: str,
    enzyme_ii: str,
    plasmid_compatibility: pd.DataFrame,
) -> list[str]:
    if enzyme_i not in plasmid_compatibility.index or enzyme_ii not in plasmid_compatibility.index:
        return []
    return [
        plasmid
        for plasmid in PLASMIDS
        if plasmid in plasmid_compatibility.columns
        and bool(plasmid_compatibility.loc[enzyme_i, plasmid])
        and bool(plasmid_compatibility.loc[enzyme_ii, plasmid])
    ]


def enumerate_active_latent_pairs(
    hits: list[RestrictionHit],
    plasmid_compatibility: pd.DataFrame,
) -> list[dict[str, object]]:
    """Enumerate spatially useful Site-I/Site-II pairs with at least one latent site.

    To avoid a quadratic expansion in long low-complexity DNA, every latent
    hit is paired with the nearest Site-I/II occurrences for each distinct
    enzyme.  These are the only pairs that can improve the linear step-count
    objective without adding a longer replaceable interval.
    """
    site_i_by_enzyme: dict[str, list[RestrictionHit]] = {}
    site_ii_by_enzyme: dict[str, list[RestrictionHit]] = {}
    for hit in hits:
        if hit.site_i_eligible:
            site_i_by_enzyme.setdefault(hit.canonical_enzyme, []).append(hit)
        if hit.site_ii_eligible:
            site_ii_by_enzyme.setdefault(hit.canonical_enzyme, []).append(hit)
    candidates: dict[tuple[Any, ...], dict[str, object]] = {}
    latent_hits = [hit for hit in hits if hit.state == "latent"]
    for latent in latent_hits:
        role_pairs: list[tuple[RestrictionHit, RestrictionHit]] = []
        if latent.site_i_eligible:
            for enzyme, other_hits in site_ii_by_enzyme.items():
                if enzyme != latent.canonical_enzyme:
                    role_pairs.extend((latent, other) for other in _nearest_hits(other_hits, latent.start))
        if latent.site_ii_eligible:
            for enzyme, other_hits in site_i_by_enzyme.items():
                if enzyme != latent.canonical_enzyme:
                    role_pairs.extend((other, latent) for other in _nearest_hits(other_hits, latent.start))
        for site_i, site_ii in role_pairs:
            if site_i.canonical_enzyme == site_ii.canonical_enzyme:
                continue
            left, right = sorted((site_i, site_ii), key=lambda item: (item.start, item.end))
            replacement_start = max(left.top_cut, left.bottom_cut)
            replacement_end = min(right.top_cut, right.bottom_cut)
            if replacement_start >= replacement_end:
                continue
            latent_positions = [
                hit.mismatch_position
                for hit in (site_i, site_ii)
                if hit.state == "latent" and hit.mismatch_position is not None
            ]
            donor_derived = [
                position for position in latent_positions if replacement_start <= int(position) < replacement_end
            ]
            if not donor_derived:
                continue
            plasmids = _plasmid_options(site_i.enzyme, site_ii.enzyme, plasmid_compatibility)
            if not plasmids:
                continue
            mechanism = "latent+latent" if site_i.state == site_ii.state == "latent" else "active+latent"
            key = (
                site_i.canonical_enzyme,
                site_ii.canonical_enzyme,
                site_i.start,
                site_ii.start,
                site_i.orientation,
                site_ii.orientation,
            )
            candidates[key] = {
                "site_i": site_i,
                "site_ii": site_ii,
                "left_hit": left,
                "right_hit": right,
                "replacement_start": replacement_start,
                "replacement_end": replacement_end,
                "replacement_length_bp": replacement_end - replacement_start,
                "donor_derived_mutation_positions": donor_derived,
                "mutation_between_overhangs": True,
                "mechanism": mechanism,
                "compatible_plasmids": plasmids,
            }
    return sorted(
        candidates.values(),
        key=lambda item: (
            int(item["replacement_length_bp"]),
            str(item["site_i"].canonical_enzyme),
            str(item["site_ii"].canonical_enzyme),
            int(item["replacement_start"]),
        ),
    )


def _select_type_iis(
    overhang: int,
    sequence: str,
    geometries: dict[str, EnzymeGeometry],
) -> EnzymeGeometry | None:
    candidates = [
        geometry
        for geometry in geometries.values()
        if geometry.site_iii_eligible
        and geometry.ovhg == overhang
        and geometry.recognition_site not in sequence
        and reverse_complement(geometry.recognition_site) not in sequence
        and geometry.ligation_ok
        and geometry.no_star_activity
    ]
    return sorted(
        candidates,
        key=lambda item: (
            not item.methylation_compatible,
            len(item.recognition_site),
            item.canonical_enzyme,
            item.enzyme,
        ),
    )[0] if candidates else None


def _type_iis_flank(
    geometry: EnzymeGeometry,
    *,
    left: bool,
    overhang_sequence: str,
) -> str:
    """Build a disposable adapter that exposes the requested cohesive end.

    For an outward-cutting enzyme the bases between its top and bottom cuts
    determine the donor overhang.  They are adapter bases, not target bases,
    and are removed or become the transient cohesive end after digestion.
    """
    overhang = validate_dna(overhang_sequence)
    if len(overhang) != geometry.overhang_length:
        raise ValueError(
            f"{geometry.enzyme} requires a {geometry.overhang_length}-nt overhang, got {len(overhang)}"
        )
    spacer_length = max(0, geometry.top_cut_offset - len(geometry.recognition_site))
    spacer = ("ACGT" * (math.ceil(spacer_length / 4)))[:spacer_length]
    if left:
        return DEFAULT_CLAMP + geometry.recognition_site + spacer + overhang
    return reverse_complement(overhang + spacer + geometry.recognition_site) + DEFAULT_CLAMP


def _product_type(length: int) -> str | None:
    if DEFAULT_OLIGO_MIN_BP <= length <= DEFAULT_OLIGO_MAX_BP:
        return "duplexed_ultramer"
    if DEFAULT_GBLOCK_MIN_BP <= length <= DEFAULT_GBLOCK_MAX_BP:
        return "gblock"
    return None


def _split_interval(
    sequence: str,
    *,
    left_adapter: str,
    right_adapter: str,
    max_purchase_bp: int = DEFAULT_GBLOCK_MAX_BP,
) -> list[tuple[int, int, str]]:
    max_core = max_purchase_bp - len(left_adapter) - len(right_adapter)
    if max_core < DEFAULT_OLIGO_MIN_BP:
        raise ValueError("Adapters leave no synthesis capacity for a target fragment")
    if not sequence:
        return []
    step_count = math.ceil(len(sequence) / max_core)
    base, remainder = divmod(len(sequence), step_count)
    lengths = [base + (1 if index < remainder else 0) for index in range(step_count)]
    fragments: list[tuple[int, int, str]] = []
    cursor = 0
    for length in lengths:
        end = cursor + length
        fragments.append((cursor, end, sequence[cursor:end]))
        cursor = end
    return fragments


def _active_site_positions(sequence: str, geometries: Iterable[EnzymeGeometry]) -> dict[str, list[int]]:
    positions: dict[str, list[int]] = {}
    for geometry in geometries:
        found: set[int] = set()
        for site in {geometry.recognition_site, reverse_complement(geometry.recognition_site)}:
            start = sequence.find(site)
            while start >= 0:
                found.add(start)
                start = sequence.find(site, start + 1)
        positions[geometry.canonical_enzyme] = sorted(found)
    return positions


def _temporary_active_sequence(
    sequence: str,
    hits: Iterable[RestrictionHit],
    *,
    replacement_start: int | None = None,
    replacement_end: int | None = None,
    assembled_length: int | None = None,
    direction: str = "left_to_right",
) -> str:
    mutable = list(sequence)
    for hit in hits:
        if hit.state == "latent" and hit.mismatch_position is not None and hit.active_base is not None:
            position = int(hit.mismatch_position)
            if replacement_start is not None and replacement_end is not None and assembled_length is not None:
                removed_length = replacement_end - replacement_start
                if position >= replacement_end:
                    position += assembled_length - removed_length
                elif replacement_start <= position < replacement_end:
                    if direction == "left_to_right":
                        if position >= replacement_start + assembled_length:
                            continue
                    else:
                        assembled_source_start = replacement_end - assembled_length
                        if position < assembled_source_start:
                            continue
                        position = replacement_start + (position - assembled_source_start)
            if 0 <= position < len(mutable):
                mutable[position] = str(hit.active_base)
    return "".join(mutable)


def _route_sort_key(route: dict[str, object]) -> tuple[Any, ...]:
    return (
        not bool(route["passed"]),
        int(route["step_count"]),
        int(route["unique_fragment_count"]),
        int(route["total_purchase_bp"]),
        float(route.get("maximum_idt_score") if route.get("maximum_idt_score") is not None else math.inf),
        str(route["site_i_enzyme"]),
        str(route["site_ii_enzyme"]),
        str(route["plasmid"]),
        str(route["direction"]),
    )


def is_exact_repeat_gain_route(
    target: TargetRecord,
    replacement_start: int,
    replacement_end: int,
    *,
    repeat_unit_gain: int = 1,
) -> bool:
    """Return whether a replacement expands an exact tandem array by full units.

    Removing the proposed donor interval from the final target must recover the
    exact shorter precursor. This prevents a formally valid short repair inside
    an already complete array from being reported as repeat-array expansion.
    """
    if repeat_unit_gain < 1:
        raise ValueError("repeat_unit_gain must be at least one")
    unit = validate_dna(target.unit_sequence, allow_empty=True)
    if not unit or target.copy_count is None:
        raise ValueError(
            "Exact repeat-gain validation requires unit_sequence and copy_count"
        )
    copies = int(target.copy_count)
    if copies <= repeat_unit_gain:
        raise ValueError("copy_count must exceed repeat_unit_gain")
    final_sequence = validate_dna(target.sequence)
    if final_sequence != unit * copies:
        raise ValueError(
            "Exact repeat-gain validation requires sequence == unit_sequence * copy_count"
        )
    precursor = final_sequence[:replacement_start] + final_sequence[replacement_end:]
    return precursor == unit * (copies - repeat_unit_gain)


def exact_repeat_gain_for_interval(
    target: TargetRecord,
    replacement_start: int,
    replacement_end: int,
) -> int | None:
    """Return the number of exact tandem units removed by an interval.

    This is the batched counterpart to :func:`is_exact_repeat_gain_route`.
    It allows the complete-route planner to scan one final copy-number state
    once and then reuse its molecular routes for every shorter precursor.
    """
    unit = validate_dna(target.unit_sequence, allow_empty=True)
    if not unit or target.copy_count is None:
        return None
    copies = int(target.copy_count)
    removed_length = replacement_end - replacement_start
    if removed_length <= 0 or removed_length % len(unit):
        return None
    gain = removed_length // len(unit)
    if gain < 1 or gain >= copies:
        return None
    final_sequence = validate_dna(target.sequence)
    if final_sequence != unit * copies:
        return None
    precursor = final_sequence[:replacement_start] + final_sequence[replacement_end:]
    return gain if precursor == unit * (copies - gain) else None


def _score_idt_safely(
    scorer: IDTComplexityScorer | None,
    name: str,
    sequence: str,
) -> dict[str, Any]:
    if scorer is None:
        return {
            "idt_status": "not_run",
            "idt_explicit_pass": None,
            "idt_complexity_score": None,
            "idt_response_sha256": "",
            "idt_scored_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        }
    try:
        return scorer.score(name, sequence)
    except Exception:
        # API/library exception text is intentionally not propagated because
        # upstream responses can echo account context.  The sequence hash is
        # sufficient for a missing-only rerun and contains no credential.
        return {
            "idt_status": "api_failure",
            "idt_explicit_pass": None,
            "idt_complexity_score": None,
            "idt_response_sha256": "",
            "idt_scored_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        }


def plan_target(
    target: TargetRecord,
    geometries: dict[str, EnzymeGeometry],
    plasmid_compatibility: pd.DataFrame,
    *,
    idt_scorer: IDTComplexityScorer | None = None,
    require_idt: bool = True,
    top_routes: int = DEFAULT_TOP_ROUTES,
    max_purchase_bp: int = DEFAULT_GBLOCK_MAX_BP,
    min_replacement_bp: int = 0,
    required_repeat_unit_gain: int = 0,
    enumerate_repeat_unit_gains: bool = False,
    _breakpoint_attempts: tuple[int, ...] = (),
    _precomputed_hits: list[RestrictionHit] | None = None,
    _precomputed_schemes: list[dict[str, object]] | None = None,
) -> dict[str, pd.DataFrame]:
    """Plan and audit the best linear routes for one immutable target DNA."""
    sequence = validate_dna(target.sequence)
    if min_replacement_bp < 0:
        raise ValueError("min_replacement_bp must be non-negative")
    if required_repeat_unit_gain < 0:
        raise ValueError("required_repeat_unit_gain must be non-negative")
    if required_repeat_unit_gain:
        # Validate exact-tandem metadata before doing the expensive hit scan.
        is_exact_repeat_gain_route(
            target,
            0,
            len(sequence),
            repeat_unit_gain=required_repeat_unit_gain,
        )
    if len(sequence) < PRIMER_PAIR_CORE_THRESHOLD_BP:
        whole_target_idt = {
            "idt_status": "not_applicable_primer_pair_under_90bp",
            "idt_explicit_pass": True,
            "idt_complexity_score": None,
            "idt_response_sha256": "",
            "idt_scored_sequence_sha256": "",
        }
    else:
        whole_target_idt = _score_idt_safely(
            idt_scorer,
            f"{target.target_id}|whole_target",
            sequence,
        )
    # When the intact target is explicitly rejected by the IDT score policy,
    # move breakpoints into the preferred short-fragment overlap.  Only the
    # disposable adapters and target partition change; the concatenated donor
    # cores remain the immutable target interval.
    effective_max_purchase_bp = max_purchase_bp
    breakpoint_policy = "fewest_steps_up_to_3000bp"
    if require_idt and whole_target_idt.get("idt_explicit_pass") is not True:
        effective_max_purchase_bp = min(max_purchase_bp, DEFAULT_OLIGO_MAX_BP)
        breakpoint_policy = f"idt_rejected_whole_target_retry_at_{effective_max_purchase_bp}bp"
    breakpoint_attempts = (*_breakpoint_attempts, effective_max_purchase_bp)
    hits = (
        _precomputed_hits
        if _precomputed_hits is not None else
        scan_re_sites(target.target_id, sequence, geometries)
    )
    schemes = (
        _precomputed_schemes
        if _precomputed_schemes is not None else
        enumerate_active_latent_pairs(hits, plasmid_compatibility)
    )
    hit_rows = [hit.to_dict() for hit in hits]
    route_rows: list[dict[str, object]] = []
    step_rows: list[dict[str, object]] = []
    fragment_rows: list[dict[str, object]] = []
    for scheme_index, scheme in enumerate(schemes):
        site_i: RestrictionHit = scheme["site_i"]  # type: ignore[assignment]
        site_ii: RestrictionHit = scheme["site_ii"]  # type: ignore[assignment]
        left_hit: RestrictionHit = scheme["left_hit"]  # type: ignore[assignment]
        right_hit: RestrictionHit = scheme["right_hit"]  # type: ignore[assignment]
        selected_geometries = [geometries[site_i.enzyme], geometries[site_ii.enzyme]]
        replacement_start = int(scheme["replacement_start"])
        replacement_end = int(scheme["replacement_end"])
        if replacement_end - replacement_start < min_replacement_bp:
            continue
        repeat_unit_gain = exact_repeat_gain_for_interval(
            target, replacement_start, replacement_end
        )
        if required_repeat_unit_gain and repeat_unit_gain != required_repeat_unit_gain:
            continue
        if enumerate_repeat_unit_gains and repeat_unit_gain is None:
            continue
        donor_core = sequence[replacement_start:replacement_end]
        # Cores below 90 bp are purchased as two complementary oligos that
        # already expose the selected sticky ends.  Likewise, a deliberately
        # narrow fragmentation ceiling (60/55 bp in complete-route recovery)
        # guarantees primer-pair fragments.  Neither case needs disposable
        # Type IIS adapters, so requiring an adapter here would create a
        # systematic false negative before the primer branch is reached.
        primer_only_fragmentation = (
            len(donor_core) < PRIMER_PAIR_CORE_THRESHOLD_BP
            or effective_max_purchase_bp < PRIMER_PAIR_CORE_THRESHOLD_BP
        )
        adapter_left: EnzymeGeometry | None = None
        adapter_right: EnzymeGeometry | None = None
        if primer_only_fragmentation:
            left_flank = ""
            right_flank = ""
        else:
            adapter_left = _select_type_iis(
                left_hit.ovhg, sequence, geometries
            )
            adapter_right = _select_type_iis(
                right_hit.ovhg, sequence, geometries
            )
            if adapter_left is None or adapter_right is None:
                continue
            left_flank = _type_iis_flank(
                adapter_left,
                left=True,
                overhang_sequence=left_hit.ovhgseq,
            )
            right_flank = _type_iis_flank(
                adapter_right,
                left=False,
                overhang_sequence=right_hit.ovhgseq,
            )
        try:
            raw_fragments = _split_interval(
                donor_core,
                left_adapter=left_flank,
                right_adapter=right_flank,
                max_purchase_bp=effective_max_purchase_bp,
            )
        except ValueError:
            continue
        if not raw_fragments:
            continue
        active_target = _active_site_positions(sequence, selected_geometries)
        # The molecular identity includes the role, strand and cut geometry of
        # both sites.  Earlier IDs omitted positions/orientations and collided
        # when two latent occurrences exchanged Site-I/Site-II roles.
        pair_payload = {
            "target_id": target.target_id,
            "site_i": {
                "enzyme": site_i.enzyme,
                "position": site_i.start,
                "orientation": site_i.orientation,
                "top_cut": site_i.top_cut,
                "bottom_cut": site_i.bottom_cut,
            },
            "site_ii": {
                "enzyme": site_ii.enzyme,
                "position": site_ii.start,
                "orientation": site_ii.orientation,
                "top_cut": site_ii.top_cut,
                "bottom_cut": site_ii.bottom_cut,
            },
            "replacement": [replacement_start, replacement_end],
            "site_iii": [
                adapter_left.enzyme if adapter_left is not None else "not_required_primer_pair",
                adapter_right.enzyme if adapter_right is not None else "not_required_primer_pair",
            ],
        }
        pair_id = hashlib.sha256(
            json.dumps(pair_payload, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()[:16]
        for plasmid in scheme["compatible_plasmids"]:  # type: ignore[union-attr]
            for direction in ("left_to_right", "right_to_left"):
                route_identity = {
                    "pair_id": pair_id,
                    "plasmid": plasmid,
                    "direction": direction,
                    "fragment_purchase_ceiling_bp": effective_max_purchase_bp,
                }
                route_hash = hashlib.sha256(
                    json.dumps(
                        route_identity, sort_keys=True, separators=(",", ":")
                    ).encode()
                ).hexdigest()[:16]
                route_id = f"{target.target_id}|route={route_hash}"
                ordered = raw_fragments if direction == "left_to_right" else list(reversed(raw_fragments))
                local_fragments: list[dict[str, object]] = []
                route_steps: list[dict[str, object]] = []
                assembled_parts: list[str] = []
                route_valid = True
                idt_scores: list[float] = []
                for step_number, (relative_start, relative_end, core) in enumerate(ordered, start=1):
                    is_primer_pair = len(core) < PRIMER_PAIR_CORE_THRESHOLD_BP
                    if is_primer_pair:
                        primer_forward = left_hit.ovhgseq + core
                        primer_reverse = right_hit.ovhgseq + reverse_complement(core)
                        purchase = primer_forward
                        secondary_purchase = primer_reverse
                        purchase_length = len(primer_forward) + len(primer_reverse)
                        product_type = "annealed_sticky_end_primer_pair"
                        fragment_hash_input = (
                            f"{product_type}|{primer_forward}|{primer_reverse}"
                        )
                        local_left_adapter = ""
                        local_right_adapter = ""
                        local_left_adapter_enzyme = "not_required"
                        local_right_adapter_enzyme = "not_required"
                        local_site_iii_left = "not_required_primer_pair"
                        local_site_iii_right = "not_required_primer_pair"
                        donor_digest_sequence = left_hit.ovhgseq + core + right_hit.ovhgseq
                        donor_type_iis_sites = {}
                        expected_type_iis_sites = 0
                        actual_type_iis_sites = 0
                        idt_summary: dict[str, Any] = {
                            "idt_status": "not_applicable_primer_pair_under_90bp",
                            "idt_explicit_pass": True,
                            "idt_complexity_score": None,
                            "idt_response_sha256": "",
                            "idt_scored_sequence_sha256": "",
                        }
                    else:
                        purchase = left_flank + core + right_flank
                        secondary_purchase = ""
                        purchase_length = len(purchase)
                        product_type = _product_type(purchase_length)
                        if product_type is None:
                            route_valid = False
                            break
                        primer_forward = ""
                        primer_reverse = ""
                        fragment_hash_input = f"{product_type}|{purchase}"
                        local_left_adapter = left_flank
                        local_right_adapter = right_flank
                        local_left_adapter_enzyme = adapter_left.enzyme
                        local_right_adapter_enzyme = adapter_right.enzyme
                        local_site_iii_left = adapter_left.enzyme
                        local_site_iii_right = adapter_right.enzyme
                        donor_digest_sequence = core
                        donor_type_iis_sites = _active_site_positions(
                            purchase,
                            {adapter_left, adapter_right},
                        )
                        expected_type_iis_sites = 2
                        actual_type_iis_sites = sum(
                            len(positions) for positions in donor_type_iis_sites.values()
                        )
                        idt_summary = {
                            "idt_status": "not_run",
                            "idt_explicit_pass": None,
                            "idt_complexity_score": None,
                            "idt_response_sha256": "",
                        }
                    fragment_sha = hashlib.sha256(fragment_hash_input.encode()).hexdigest()
                    fragment_id = f"frag_{fragment_sha[:16]}"
                    donor_selected_sites = _active_site_positions(
                        donor_digest_sequence, selected_geometries
                    )
                    donor_unintended_selected_cuts = sum(
                        len(positions) for positions in donor_selected_sites.values()
                    )
                    if donor_unintended_selected_cuts or actual_type_iis_sites != expected_type_iis_sites:
                        route_valid = False
                    local_fragments.append(
                        {
                            "target_id": target.target_id,
                            "route_id": route_id,
                            "fragment_id": fragment_id,
                            "step_number": step_number,
                            "target_start": replacement_start + relative_start,
                            "target_end": replacement_start + relative_end,
                            "core_sequence": core,
                            "core_length_bp": len(core),
                            "purchase_sequence": purchase,
                            "secondary_purchase_sequence": secondary_purchase,
                            "primer_forward_5to3": primer_forward,
                            "primer_reverse_5to3": primer_reverse,
                            "primer_pair_exposes_sticky_ends": is_primer_pair,
                            "purchase_sequence_count": 2 if is_primer_pair else 1,
                            "purchase_length_bp": purchase_length,
                            "purchase_sha256": fragment_sha,
                            "product_type": product_type,
                            "left_adapter_enzyme": local_left_adapter_enzyme,
                            "right_adapter_enzyme": local_right_adapter_enzyme,
                            "left_adapter_sequence": local_left_adapter,
                            "right_adapter_sequence": local_right_adapter,
                            "digest_fragment_sequence": donor_digest_sequence,
                            "digest_fragment_sha256": hashlib.sha256(donor_digest_sequence.encode()).hexdigest(),
                            "left_digest_overhang": left_hit.ovhgseq,
                            "right_digest_overhang": right_hit.ovhgseq,
                            "donor_selected_re_sites_json": json.dumps(donor_selected_sites, sort_keys=True),
                            "donor_unintended_selected_cut_count": donor_unintended_selected_cuts,
                            "donor_type_iis_sites_json": json.dumps(donor_type_iis_sites, sort_keys=True),
                            "donor_type_iis_site_count": actual_type_iis_sites,
                            "idt_policy": IDT_SCORE_POLICY,
                            "idt_status": idt_summary.get("idt_status"),
                            "idt_score": idt_summary.get("idt_complexity_score"),
                            "idt_response_sha256": idt_summary.get("idt_response_sha256", ""),
                            "idt_scored_sequence_sha256": idt_summary.get("idt_scored_sequence_sha256", ""),
                        }
                    )
                    if direction == "left_to_right":
                        assembled_parts.append(core)
                    else:
                        assembled_parts.insert(0, core)
                    assembled = "".join(assembled_parts)
                    is_final = step_number == len(ordered)
                    if direction == "left_to_right":
                        intermediate = sequence[:replacement_start] + assembled + sequence[replacement_end:]
                    else:
                        intermediate = sequence[:replacement_start] + assembled + sequence[replacement_end:]
                    if not is_final:
                        intermediate = _temporary_active_sequence(
                            intermediate,
                            (left_hit, right_hit),
                            replacement_start=replacement_start,
                            replacement_end=replacement_end,
                            assembled_length=len(assembled),
                            direction=direction,
                        )
                    active_positions = _active_site_positions(intermediate, selected_geometries)
                    unintended = sum(
                        max(0, len(values) - 1) for values in active_positions.values()
                    )
                    if unintended:
                        route_valid = False
                    route_steps.append(
                        {
                            "target_id": target.target_id,
                            "route_id": route_id,
                            "step_number": step_number,
                            "direction": direction,
                            "plasmid": plasmid,
                            "fragment_id": fragment_id,
                            "recipient_insert_length_bp": len(intermediate) - len(core),
                            "donor_core_length_bp": len(core),
                            "result_insert_length_bp": len(intermediate),
                            "site_i_enzyme": site_i.enzyme,
                            "site_ii_enzyme": site_ii.enzyme,
                            "site_iii_left": local_site_iii_left,
                            "site_iii_right": local_site_iii_right,
                            "site_i_top_cut": site_i.top_cut,
                            "site_i_bottom_cut": site_i.bottom_cut,
                            "site_ii_top_cut": site_ii.top_cut,
                            "site_ii_bottom_cut": site_ii.bottom_cut,
                            "site_i_overhang": site_i.ovhgseq,
                            "site_ii_overhang": site_ii.ovhgseq,
                            "site_i_target_state": site_i.state,
                            "site_ii_target_state": site_ii.state,
                            "latent_state_transition_json": json.dumps(
                                {
                                    hit.canonical_enzyme: {
                                        "target_before": "latent",
                                        "intermediate": "active" if not is_final else "not_applicable",
                                        "final": "latent",
                                        "mismatch_position": hit.mismatch_position,
                                        "target_base": hit.target_base,
                                        "temporary_active_base": hit.active_base,
                                    }
                                    for hit in (site_i, site_ii)
                                    if hit.state == "latent"
                                },
                                sort_keys=True,
                            ),
                            "donor_top_strand_interval_json": json.dumps(
                                [replacement_start + relative_start, replacement_start + relative_end]
                            ),
                            "donor_bottom_strand_interval_json": json.dumps(
                                [
                                    len(sequence) - (replacement_start + relative_end),
                                    len(sequence) - (replacement_start + relative_start),
                                ]
                            ),
                            "double_strand_source_verified": True,
                            "active_sites_json": json.dumps(active_positions, sort_keys=True),
                            "unintended_cut_count": unintended,
                            "is_final_step": is_final,
                            "intermediate_sequence_sha256": hashlib.sha256(intermediate.encode()).hexdigest(),
                        }
                    )
                concatenated = "".join(item[2] for item in raw_fragments)
                final_sequence = sequence[:replacement_start] + concatenated + sequence[replacement_end:]
                final_exact = final_sequence == sequence
                route_valid = route_valid and final_exact
                unique_fragments = {str(row["purchase_sha256"]) for row in local_fragments}
                route = {
                    "version": DNA_ASSEMBLY_VERSION,
                    "target_id": target.target_id,
                    "route_id": route_id,
                    "scheme_index": scheme_index,
                    "passed": bool(route_valid and not require_idt),
                    "local_constraints_passed": bool(route_valid),
                    "status": (
                        "pending_idt" if route_valid and require_idt
                        else "passed" if route_valid
                        else "failed_route_constraints"
                    ),
                    "mechanism": scheme["mechanism"],
                    "growth_mode": "growth_front" if direction == "left_to_right" else "terminal_only",
                    "direction": direction,
                    "plasmid": plasmid,
                    "site_i_enzyme": site_i.enzyme,
                    "site_ii_enzyme": site_ii.enzyme,
                    "site_i_state": site_i.state,
                    "site_ii_state": site_ii.state,
                    "site_i_position": site_i.start,
                    "site_ii_position": site_ii.start,
                    "site_i_orientation": site_i.orientation,
                    "site_ii_orientation": site_ii.orientation,
                    "site_i_overhang": site_i.ovhgseq,
                    "site_ii_overhang": site_ii.ovhgseq,
                    "site_iii_left": (
                        adapter_left.enzyme
                        if adapter_left is not None else
                        "not_required_primer_pair"
                    ),
                    "site_iii_right": (
                        adapter_right.enzyme
                        if adapter_right is not None else
                        "not_required_primer_pair"
                    ),
                    "replacement_start": replacement_start,
                    "replacement_end": replacement_end,
                    "replacement_length_bp": replacement_end - replacement_start,
                    "repeat_unit_gain": repeat_unit_gain,
                    "donor_derived_mutation_positions_json": json.dumps(
                        scheme["donor_derived_mutation_positions"]
                    ),
                    "mutation_between_overhangs": True,
                    "target_active_sites_json": json.dumps(active_target, sort_keys=True),
                    "step_count": len(local_fragments),
                    "unique_fragment_count": len(unique_fragments),
                    "purchase_sequence_count": sum(
                        int(row["purchase_sequence_count"]) for row in local_fragments
                    ),
                    "total_purchase_bp": sum(int(row["purchase_length_bp"]) for row in local_fragments),
                    "fragment_purchase_ceiling_bp": effective_max_purchase_bp,
                    "breakpoint_policy": breakpoint_policy,
                    "maximum_idt_score": None,
                    "all_fragments_idt_passed": False,
                    "all_purchase_fragments_accepted": False,
                    "final_sequence_exact": final_exact,
                    "target_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                    "final_sequence_sha256": hashlib.sha256(final_sequence.encode()).hexdigest(),
                }
                route_rows.append(route)
                step_rows.extend(route_steps)
                fragment_rows.extend(local_fragments)
    # Rank on local constraints first.  Only retained routes are sent to IDT;
    # scoring every discarded plasmid/direction duplicate would multiply live
    # API traffic without changing the top-10 route set.
    route_rows.sort(
        key=lambda route: (
            not bool(route["local_constraints_passed"]),
            int(route["step_count"]),
            int(route["unique_fragment_count"]),
            int(route["total_purchase_bp"]),
            str(route["site_i_enzyme"]),
            str(route["site_ii_enzyme"]),
            str(route["plasmid"]),
            str(route["direction"]),
        )
    )
    if enumerate_repeat_unit_gains:
        retained_rows: list[dict[str, object]] = []
        per_gain_counts: dict[int, int] = {}
        for row in route_rows:
            gain = int(row["repeat_unit_gain"])
            if per_gain_counts.get(gain, 0) >= top_routes:
                continue
            retained_rows.append(row)
            per_gain_counts[gain] = per_gain_counts.get(gain, 0) + 1
        route_rows = retained_rows
    else:
        route_rows = route_rows[:top_routes]
    selected_route_ids = {str(row["route_id"]) for row in route_rows}
    step_rows = [row for row in step_rows if str(row["route_id"]) in selected_route_ids]
    fragment_rows = [row for row in fragment_rows if str(row["route_id"]) in selected_route_ids]
    if require_idt:
        summaries_by_sha: dict[str, dict[str, Any]] = {}
        for fragment in fragment_rows:
            if fragment.get("product_type") == "annealed_sticky_end_primer_pair":
                continue
            purchase_sha = str(fragment["purchase_sha256"])
            if purchase_sha not in summaries_by_sha:
                summaries_by_sha[purchase_sha] = _score_idt_safely(
                    idt_scorer,
                    str(fragment["fragment_id"]),
                    str(fragment["purchase_sequence"]),
                )
            idt_summary = summaries_by_sha[purchase_sha]
            fragment["idt_status"] = idt_summary.get("idt_status")
            fragment["idt_score"] = idt_summary.get("idt_complexity_score")
            fragment["idt_response_sha256"] = idt_summary.get("idt_response_sha256", "")
            fragment["idt_scored_sequence_sha256"] = idt_summary.get(
                "idt_scored_sequence_sha256", ""
            )
        fragments_by_route: dict[str, list[dict[str, object]]] = {}
        for fragment in fragment_rows:
            fragments_by_route.setdefault(str(fragment["route_id"]), []).append(fragment)
        for route in route_rows:
            route_fragments = fragments_by_route.get(str(route["route_id"]), [])
            explicit_pass = bool(route_fragments) and all(
                fragment.get("idt_status")
                in {"passed", "not_applicable_primer_pair_under_90bp"}
                for fragment in route_fragments
            )
            scores = [
                float(fragment["idt_score"])
                for fragment in route_fragments
                if fragment.get("idt_score") is not None
            ]
            route["all_fragments_idt_passed"] = explicit_pass
            route["all_purchase_fragments_accepted"] = explicit_pass
            route["maximum_idt_score"] = max(scores) if scores else None
            route["passed"] = bool(route["local_constraints_passed"] and explicit_pass)
            route["status"] = (
                "passed" if route["passed"]
                else "failed_idt" if route["local_constraints_passed"]
                else "failed_route_constraints"
            )
        route_rows.sort(key=_route_sort_key)
    passed_routes = [row for row in route_rows if bool(row["passed"])]
    next_purchase_ceiling = (
        60 if effective_max_purchase_bp > 60
        else 55 if effective_max_purchase_bp > 55
        else None
    )
    if (
        require_idt
        and idt_scorer is not None
        and not passed_routes
        and next_purchase_ceiling is not None
        and any(bool(row.get("local_constraints_passed")) for row in route_rows)
    ):
        return plan_target(
            target,
            geometries,
            plasmid_compatibility,
            idt_scorer=idt_scorer,
            require_idt=True,
            top_routes=top_routes,
            max_purchase_bp=next_purchase_ceiling,
            min_replacement_bp=min_replacement_bp,
            required_repeat_unit_gain=required_repeat_unit_gain,
            enumerate_repeat_unit_gains=enumerate_repeat_unit_gains,
            _breakpoint_attempts=breakpoint_attempts,
        )
    best = passed_routes[0] if passed_routes else (route_rows[0] if route_rows else None)
    best_fragments_passed = bool(best and best.get("all_fragments_idt_passed"))
    whole_target_failed = whole_target_idt.get("idt_explicit_pass") is False
    summary = {
        **target.to_dict(),
        "version": DNA_ASSEMBLY_VERSION,
        "active_hit_count": sum(hit.state == "active" for hit in hits),
        "latent_hit_count": sum(hit.state == "latent" for hit in hits),
        "candidate_pair_count": len(schemes),
        "route_count_retained": len(route_rows),
        "hurdler_compatible": bool(passed_routes),
        "status": "passed" if passed_routes else "no_complete_route",
        "failure_reason": "" if passed_routes else (
            "no_active_latent_pair" if not schemes else "no_route_passed_vector_idt_and_digest_constraints"
        ),
        "best_route_id": str(best["route_id"]) if best else "",
        "best_plasmid": str(best["plasmid"]) if best else "",
        "best_site_i": str(best["site_i_enzyme"]) if best else "",
        "best_site_ii": str(best["site_ii_enzyme"]) if best else "",
        "best_step_count": int(best["step_count"]) if best else None,
        "best_unique_fragment_count": int(best["unique_fragment_count"]) if best else None,
        "whole_target_idt_status": whole_target_idt.get("idt_status"),
        "whole_target_idt_score": whole_target_idt.get("idt_complexity_score"),
        "whole_target_idt_response_sha256": whole_target_idt.get("idt_response_sha256", ""),
        "whole_target_idt_scored_sequence_sha256": whole_target_idt.get(
            "idt_scored_sequence_sha256", ""
        ),
        "fragment_rescued_by_hurdler": bool(
            whole_target_failed and passed_routes and best_fragments_passed
        ),
        "breakpoint_policy": breakpoint_policy,
        "fragment_purchase_ceiling_bp": effective_max_purchase_bp,
        "breakpoint_attempts_json": json.dumps(breakpoint_attempts),
        "candidate_repeat_unit_gains_json": json.dumps(
            sorted(
                {
                    int(row["repeat_unit_gain"])
                    for row in route_rows
                    if row.get("repeat_unit_gain") is not None
                }
            )
        ),
        "idt_required": require_idt,
    }
    return {
        "summary": pd.DataFrame([summary]),
        "hits": pd.DataFrame(hit_rows),
        "routes": pd.DataFrame(route_rows),
        "steps": pd.DataFrame(step_rows),
        "fragments": pd.DataFrame(fragment_rows),
    }


def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _optional_text(value: object, default: str = "") -> str:
    return default if value is None or pd.isna(value) else str(value)


def read_target_catalog(
    *,
    fasta: str | Path | None = None,
    table: str | Path | None = None,
) -> list[TargetRecord]:
    records: list[TargetRecord] = []
    if fasta is not None:
        for record in SeqIO.parse(str(fasta), "fasta"):
            records.append(TargetRecord(target_id=str(record.id), sequence=validate_dna(str(record.seq))))
    if table is not None:
        frame = _read_table(Path(table))
        if "target_id" not in frame.columns:
            raise ValueError("Target table requires target_id")
        for row in frame.to_dict("records"):
            unit = validate_dna(_optional_text(row.get("unit_sequence", "")), allow_empty=True)
            sequence_value = _optional_text(row.get("sequence", ""))
            copy_value = row.get("copy_count")
            if sequence_value:
                sequence = validate_dna(sequence_value)
            elif unit and pd.notna(copy_value):
                sequence = unit * int(copy_value)
            else:
                raise ValueError(f"{row['target_id']} requires sequence or unit_sequence+copy_count")
            records.append(
                TargetRecord(
                    target_id=str(row["target_id"]),
                    sequence=sequence,
                    cohort=_optional_text(row.get("cohort"), "user_input"),
                    architecture=_optional_text(row.get("architecture"), "unspecified"),
                    source_url=_optional_text(row.get("source_url")),
                    source_accession=_optional_text(row.get("source_accession")),
                    unit_sequence=unit,
                    copy_count=int(copy_value) if pd.notna(copy_value) else None,
                    notes=_optional_text(row.get("notes")),
                    source_database=_optional_text(row.get("source_database")),
                    element_id=_optional_text(row.get("element_id")),
                    synthetic_unit_length_bp=(
                        int(row["synthetic_unit_length_bp"])
                        if pd.notna(row.get("synthetic_unit_length_bp")) else None
                    ),
                    gc_target=(
                        float(row["gc_target"])
                        if pd.notna(row.get("gc_target")) else None
                    ),
                    synthetic_replicate=(
                        int(row["synthetic_replicate"])
                        if pd.notna(row.get("synthetic_replicate")) else None
                    ),
                )
            )
    unique: dict[str, TargetRecord] = {}
    for record in records:
        if record.target_id in unique and unique[record.target_id].sequence != record.sequence:
            raise ValueError(f"Duplicate target_id with different sequence: {record.target_id}")
        unique[record.target_id] = record
    return list(unique.values())


def target_records_from_frame(frame: pd.DataFrame) -> list[TargetRecord]:
    """Validate an in-memory catalog without materializing a secret-prone temp file."""
    if "target_id" not in frame.columns:
        raise ValueError("Target table requires target_id")
    records: list[TargetRecord] = []
    for row in frame.to_dict("records"):
        unit = validate_dna(_optional_text(row.get("unit_sequence", "")), allow_empty=True)
        sequence_value = _optional_text(row.get("sequence", ""))
        copy_value = row.get("copy_count")
        if sequence_value:
            sequence = validate_dna(sequence_value)
        elif unit and pd.notna(copy_value):
            sequence = unit * int(copy_value)
        else:
            raise ValueError(f"{row['target_id']} requires sequence or unit_sequence+copy_count")
        records.append(
            TargetRecord(
                target_id=str(row["target_id"]),
                sequence=sequence,
                cohort=_optional_text(row.get("cohort"), "user_input"),
                architecture=_optional_text(row.get("architecture"), "unspecified"),
                source_url=_optional_text(row.get("source_url")),
                source_accession=_optional_text(row.get("source_accession")),
                unit_sequence=unit,
                copy_count=int(copy_value) if pd.notna(copy_value) else None,
                notes=_optional_text(row.get("notes")),
                source_database=_optional_text(row.get("source_database")),
                element_id=_optional_text(row.get("element_id")),
                synthetic_unit_length_bp=(
                    int(row["synthetic_unit_length_bp"])
                    if pd.notna(row.get("synthetic_unit_length_bp")) else None
                ),
                gc_target=(
                    float(row["gc_target"])
                    if pd.notna(row.get("gc_target")) else None
                ),
                synthetic_replicate=(
                    int(row["synthetic_replicate"])
                    if pd.notna(row.get("synthetic_replicate")) else None
                ),
            )
        )
    unique: dict[str, TargetRecord] = {}
    for record in records:
        if record.target_id in unique and unique[record.target_id].sequence != record.sequence:
            raise ValueError(f"Duplicate target_id with different sequence: {record.target_id}")
        unique[record.target_id] = record
    return list(unique.values())


def build_synthetic_factorial(seed: int = 42) -> pd.DataFrame:
    """Build the frozen 900-case architecture/GC/length benchmark."""
    generator = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for unit_length in (6, 12, 24, 48, 96):
        for gc_target in (0.30, 0.50, 0.70):
            probabilities = np.array([(1 - gc_target) / 2, gc_target / 2, gc_target / 2, (1 - gc_target) / 2])
            for copy_count in (2, 4, 8, 16, 32):
                for replicate in range(3):
                    motif_a = "".join(generator.choice(list("ACGT"), unit_length, p=probabilities))
                    motif_b = "".join(generator.choice(list("ACGT"), unit_length, p=probabilities))
                    spacer_length = max(4, unit_length // 4)
                    spacer = "".join(generator.choice(list("ACGT"), spacer_length, p=probabilities))
                    target_length = unit_length * copy_count
                    control = "".join(generator.choice(list("ACGT"), target_length, p=probabilities))
                    architectures = {
                        "exact_tandem": motif_a * copy_count,
                        "fixed_spacer": (motif_a + spacer) * copy_count,
                        "alternating_ab": "".join(motif_a if index % 2 == 0 else motif_b for index in range(copy_count)),
                        "nonrepetitive_control": control,
                    }
                    for architecture, sequence in architectures.items():
                        target_id = f"syn_u{unit_length}_gc{int(gc_target*100)}_c{copy_count}_r{replicate}_{architecture}"
                        rows.append(
                            {
                                "target_id": target_id,
                                "sequence": sequence,
                                "cohort": "synthetic_factorial",
                                "architecture": architecture,
                                "unit_sequence": motif_a if architecture != "nonrepetitive_control" else "",
                                "copy_count": copy_count,
                                "gc_target": gc_target,
                                "synthetic_unit_length_bp": unit_length,
                                "synthetic_replicate": replicate,
                                "source_database": "Synthetic",
                                "element_id": f"syn_u{unit_length}_gc{int(gc_target*100)}_r{replicate}",
                                "seed": seed,
                                "source_url": "",
                                "source_accession": "",
                                "notes": "deterministic synthetic benchmark",
                            }
                        )
    frame = pd.DataFrame(rows)
    if len(frame) != 900 or frame.target_id.nunique() != 900:
        raise AssertionError("Synthetic factorial must contain exactly 900 stable IDs")
    return frame


def expand_element_inventory(
    elements: pd.DataFrame,
    *,
    copy_counts: tuple[int, ...] = (2, 4, 8, 16, 32),
) -> pd.DataFrame:
    """Expand exact public elements into the frozen real-derived array series."""
    required = {"element_id", "element_sequence", "source_url"}
    missing = sorted(required - set(elements.columns))
    if missing:
        raise ValueError(f"Element inventory is missing columns: {', '.join(missing)}")
    rows: list[dict[str, object]] = []
    for element in elements.to_dict("records"):
        sequence = validate_dna(str(element["element_sequence"]))
        spacer = validate_dna(str(element.get("spacer_sequence", "") or ""), allow_empty=True)
        source_name = str(element.get("source_database", "public_element"))
        for copies in copy_counts:
            target_sequence = (sequence + spacer) * copies
            rows.append(
                {
                    "target_id": f"{source_name}|{element['element_id']}|copies={copies}",
                    "sequence": target_sequence,
                    "cohort": "real_element_derived",
                    "architecture": "fixed_spacer" if spacer else "exact_tandem",
                    "source_url": str(element["source_url"]),
                    "source_accession": str(element.get("source_accession", element["element_id"])),
                    "unit_sequence": sequence + spacer,
                    "copy_count": copies,
                    "notes": str(element.get("notes", "") or ""),
                    "element_sequence": sequence,
                    "element_length_bp": len(sequence),
                    "spacer_sequence": spacer,
                    "source_database": source_name,
                    "element_id": str(element["element_id"]),
                }
            )
    result = pd.DataFrame(rows)
    if result.target_id.duplicated().any():
        raise ValueError("Expanded element inventory produced duplicate target IDs")
    return result


def build_target_corpus(
    output: str | Path,
    *,
    source_tables: Iterable[str | Path] = (),
    include_synthetic: bool = True,
    seed: int = 42,
) -> pd.DataFrame:
    frames = [_read_table(Path(path)) for path in source_tables]
    frames = [
        expand_element_inventory(frame) if "element_sequence" in frame.columns and "sequence" not in frame.columns else frame
        for frame in frames
    ]
    if include_synthetic:
        frames.append(build_synthetic_factorial(seed))
    if not frames:
        raise ValueError("No corpus source was selected")
    combined = pd.concat(frames, ignore_index=True, sort=False)
    records = target_records_from_frame(combined)
    frame = pd.DataFrame(record.to_dict() for record in records)
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.suffix.lower() == ".csv":
        frame.to_csv(destination, index=False)
    else:
        frame.to_parquet(destination, index=False)
    return frame


def _write_frame(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)


def _write_fragment_fasta(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if frame.empty:
        path.write_text("")
        return
    sort_columns = [column for column in ("fragment_id", "route_id") if column in frame.columns]
    unique = frame.sort_values(sort_columns).drop_duplicates("fragment_id")
    lines: list[str] = []
    for row in unique.itertuples(index=False):
        if row.product_type == "annealed_sticky_end_primer_pair":
            lines.extend(
                [
                    f">{row.fragment_id}_forward product={row.product_type} orientation=5to3",
                    str(row.primer_forward_5to3),
                    f">{row.fragment_id}_reverse product={row.product_type} orientation=5to3",
                    str(row.primer_reverse_5to3),
                ]
            )
        else:
            lines.extend(
                [
                    f">{row.fragment_id} product={row.product_type} length={row.purchase_length_bp}",
                    str(row.purchase_sequence),
                ]
            )
    path.write_text("\n".join(lines) + "\n")


def unique_fragment_inventory(frame: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate by complete purchase DNA plus product type and count reuse."""
    columns = [
        "fragment_id",
        "purchase_sha256",
        "product_type",
        "purchase_sequence",
        "secondary_purchase_sequence",
        "primer_forward_5to3",
        "primer_reverse_5to3",
        "primer_pair_exposes_sticky_ends",
        "purchase_sequence_count",
        "purchase_length_bp",
        "occurrence_count",
        "target_count",
        "route_count",
    ]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    grouped = frame.groupby(
        [
            "fragment_id",
            "purchase_sha256",
            "product_type",
            "purchase_sequence",
            "secondary_purchase_sequence",
            "primer_forward_5to3",
            "primer_reverse_5to3",
            "primer_pair_exposes_sticky_ends",
            "purchase_sequence_count",
            "purchase_length_bp",
        ],
        dropna=False,
        as_index=False,
    ).agg(
        occurrence_count=("fragment_id", "size"),
        target_count=("target_id", "nunique"),
        route_count=("route_id", "nunique"),
    )
    return grouped[columns].sort_values(["purchase_length_bp", "fragment_id"]).reset_index(drop=True)


def plan_target_catalog(
    catalog: str | Path,
    reference_dir: str | Path,
    output_dir: str | Path,
    *,
    artifact_dir: str | Path | None = None,
    idt_scorer: IDTComplexityScorer | None = None,
    require_idt: bool = True,
    shard_index: int = 0,
    shard_count: int = 1,
) -> dict[str, pd.DataFrame]:
    targets = read_target_catalog(table=catalog)
    selected = [record for index, record in enumerate(targets) if index % shard_count == shard_index]
    geometries, plasmids = load_enzyme_catalog(reference_dir, artifact_dir=artifact_dir)
    collected = {name: [] for name in ("summary", "hits", "routes", "steps", "fragments")}
    for target in selected:
        result = plan_target(
            target,
            geometries,
            plasmids,
            idt_scorer=idt_scorer,
            require_idt=require_idt,
        )
        for name, frame in result.items():
            if not frame.empty:
                collected[name].append(frame)
    merged = {
        name: pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
        for name, frames in collected.items()
    }
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    for name, frame in merged.items():
        _write_frame(frame, root / f"dna_assembly_{name}.parquet")
    inventory = unique_fragment_inventory(merged["fragments"])
    _write_frame(inventory, root / "dna_unique_fragment_inventory.parquet")
    inventory.to_csv(root / "dna_unique_fragment_inventory.csv", index=False)
    _write_fragment_fasta(inventory, root / "dna_unique_fragments.fasta")
    manifest = {
        "version": DNA_ASSEMBLY_VERSION,
        "created_at": utc_now(),
        "catalog": str(Path(catalog).absolute()),
        "catalog_sha256": sha256_file(catalog),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "target_rows": len(selected),
        "output_rows": {name: len(frame) for name, frame in merged.items()},
        "idt_required": require_idt,
        "idt_policy": IDT_SCORE_POLICY,
    }
    write_json_atomic(manifest, root / "run_manifest.json")
    return merged


def finalize_target_plans(
    shard_dirs: Iterable[str | Path],
    output_dir: str | Path,
) -> dict[str, pd.DataFrame]:
    roots = [Path(path) for path in shard_dirs]
    merged: dict[str, pd.DataFrame] = {}
    for name in ("summary", "hits", "routes", "steps", "fragments"):
        frames = [pd.read_parquet(root / f"dna_assembly_{name}.parquet") for root in roots]
        frame = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
        if name == "summary" and not frame.empty:
            if frame.target_id.duplicated().any():
                raise ValueError("Duplicate target IDs across finalized shards")
        merged[name] = frame
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    for name, frame in merged.items():
        _write_frame(frame, root / f"dna_assembly_{name}.parquet")
        frame.to_csv(root / f"dna_assembly_{name}.csv", index=False)
    fragments = merged["fragments"]
    inventory = unique_fragment_inventory(fragments)
    _write_frame(inventory, root / "dna_unique_fragment_inventory.parquet")
    inventory.to_csv(root / "dna_unique_fragment_inventory.csv", index=False)
    _write_fragment_fasta(inventory, root / "dna_unique_fragments.fasta")
    summary = merged["summary"]
    real = summary[summary.cohort.isin(["real_construct", "real_element_derived"])] if not summary.empty else summary
    reviewer = {
        "version": DNA_ASSEMBLY_VERSION,
        "real_examples_tested": int(len(real)),
        "real_examples_hurdler_compatible": int(real.hurdler_compatible.fillna(False).sum()) if not real.empty else 0,
        "real_success_fraction": float(real.hurdler_compatible.fillna(False).mean()) if not real.empty else None,
        "synthetic_examples_tested": int(summary.cohort.eq("synthetic_factorial").sum()) if not summary.empty else 0,
        "synthetic_examples_hurdler_compatible": int(
            summary.loc[summary.cohort.eq("synthetic_factorial"), "hurdler_compatible"].fillna(False).sum()
        ) if not summary.empty else 0,
    }
    reviewer["response_text"] = (
        "We added a separate notebook that tested "
        f"{reviewer['real_examples_tested']} exact DNA targets. "
        f"{reviewer['real_examples_hurdler_compatible']} targets admitted at least one complete "
        "in-silico HURDLER assembly route on one of the eight maintained plasmids."
    )
    write_json_atomic(reviewer, root / "reviewer_response_numbers.json")
    (root / "reviewer_response.md").write_text(str(reviewer["response_text"]) + "\n")
    write_json_atomic(
        {
            "version": DNA_ASSEMBLY_VERSION,
            "created_at": utc_now(),
            "shard_dirs": [str(path.absolute()) for path in roots],
            "rows": {name: len(frame) for name, frame in merged.items()},
        },
        root / "finalize_manifest.json",
    )
    return merged


def plot_dna_assembly_report(summary: pd.DataFrame, output_dir: str | Path) -> list[Path]:
    import matplotlib.pyplot as plt
    import seaborn as sns

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    frame = summary.copy()
    if frame.empty:
        raise ValueError("Cannot plot an empty DNA assembly summary")
    frame["compatibility"] = np.where(frame.hurdler_compatible, "compatible", "incompatible")
    order = sorted(frame.cohort.dropna().unique())
    counts = frame.groupby(["cohort", "compatibility"]).size().rename("count").reset_index()
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.8), facecolor="white")
    sns.barplot(data=counts, x="cohort", y="count", hue="compatibility", order=order, ax=axes[0], palette={"compatible": "#4B2E83", "incompatible": "#B7A57A"})
    scatter = frame[frame.hurdler_compatible & frame.best_step_count.notna()]
    sns.scatterplot(data=scatter, x="length_bp", y="best_step_count", hue="cohort", style="architecture", ax=axes[1], palette="colorblind")
    axes[0].set_title("Exact-DNA HURDLER compatibility")
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", rotation=25)
    axes[1].set_title("Best linear assembly route")
    axes[1].set_xlabel("Target length (bp)")
    axes[1].set_ylabel("Digest–ligation steps")
    sns.despine()
    figure.tight_layout()
    outputs = []
    for suffix in ("png", "pdf"):
        path = destination / f"dna_assembly_summary.{suffix}"
        figure.savefig(path, dpi=300, facecolor="white")
        outputs.append(path)
    plt.close(figure)
    return outputs


def plot_graphical_abstract_panel(output_dir: str | Path) -> list[Path]:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(8, 3.6), facecolor="white")
    axis.set_xlim(0, 10)
    axis.set_ylim(0, 4)
    axis.axis("off")
    colors = ["#4B2E83", "#B7A57A", "#2D7DD2", "#F45D01"]
    for index, x_value in enumerate((0.5, 2.0, 3.5)):
        axis.add_patch(FancyBboxPatch((x_value, 2.35), 1.1, 0.55, boxstyle="round,pad=0.08", facecolor=colors[index], edgecolor="none"))
    axis.text(2.35, 3.25, "Protein repeat modules", ha="center", fontsize=11)
    for index, x_value in enumerate((0.5, 1.8, 3.1, 4.4)):
        axis.add_patch(FancyBboxPatch((x_value, 0.75), 1.0, 0.5, boxstyle="round,pad=0.05", facecolor=colors[index % len(colors)], edgecolor="none"))
    axis.text(2.7, 0.25, "Regulatory and repetitive DNA arrays", ha="center", fontsize=11)
    axis.add_patch(FancyArrowPatch((5.6, 2.0), (7.0, 2.0), arrowstyle="-|>", mutation_scale=18, linewidth=2, color="#333333"))
    axis.text(6.3, 2.35, "active/latent\nRE geometry", ha="center", fontsize=9)
    axis.add_patch(FancyBboxPatch((7.2, 1.25), 2.2, 1.5, boxstyle="round,pad=0.12", facecolor="#ECE9F2", edgecolor="#4B2E83", linewidth=1.5))
    axis.text(8.3, 2.25, "HURDLER", ha="center", va="center", fontsize=15, color="#4B2E83", weight="bold")
    axis.text(8.3, 1.65, "scarless, exact-DNA\nstepwise assembly", ha="center", va="center", fontsize=9)
    figure.tight_layout()
    outputs = []
    for suffix in ("png", "pdf", "svg"):
        path = destination / f"graphical_abstract_dna_application.{suffix}"
        figure.savefig(path, dpi=300, facecolor="white", bbox_inches="tight")
        outputs.append(path)
    plt.close(figure)
    return outputs


def plot_active_latent_transition(output_dir: str | Path) -> list[Path]:
    """Render the exact-DNA state transition used in the methods notebook."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(2, 1, figsize=(8.2, 4.8), facecolor="white")
    examples = [
        ("active + latent", "GAATTC", "GAGCTT", "AAGCTT"),
        ("latent + latent", "GAATTA", "GAGCTT", "GAATTC  …  AAGCTT"),
    ]
    for axis, (label, left, target_right, temporary) in zip(axes, examples):
        axis.set_xlim(0, 12)
        axis.set_ylim(0, 2.2)
        axis.axis("off")
        axis.text(0.1, 1.85, label, weight="bold", color="#4B2E83")
        boxes = [
            (0.3, "final target", f"5′ {left} … {target_right} 3′", "#ECE9F2"),
            (4.4, "temporary active", f"5′ {temporary} 3′", "#FFF3CD"),
            (8.3, "donor-restored final", f"5′ {left} … {target_right} 3′", "#D9EAD3"),
        ]
        for x_value, title, sequence, color in boxes:
            axis.add_patch(
                FancyBboxPatch(
                    (x_value, 0.45),
                    3.1,
                    0.9,
                    boxstyle="round,pad=0.08",
                    facecolor=color,
                    edgecolor="#333333",
                    linewidth=0.8,
                )
            )
            axis.text(x_value + 1.55, 1.08, title, ha="center", fontsize=8)
            axis.text(x_value + 1.55, 0.7, sequence, ha="center", family="monospace", fontsize=8)
        for start, end in ((3.45, 4.3), (7.55, 8.15)):
            axis.add_patch(
                FancyArrowPatch(
                    (start, 0.9),
                    (end, 0.9),
                    arrowstyle="-|>",
                    mutation_scale=12,
                    color="#4B2E83",
                )
            )
    figure.tight_layout()
    outputs: list[Path] = []
    for suffix in ("png", "pdf", "svg"):
        path = destination / f"active_latent_state_transition.{suffix}"
        figure.savefig(path, dpi=300, facecolor="white", bbox_inches="tight")
        outputs.append(path)
    plt.close(figure)
    return outputs


def write_figure_manifest(
    figures: Iterable[str | Path],
    destination: str | Path,
    *,
    source_notebook: str = "notebooks/tasks/08_long_repetitive_dna_assembly.ipynb",
) -> pd.DataFrame:
    rows = []
    for value in figures:
        path = Path(value)
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing or empty figure: {path}")
        rows.append(
            {
                "version": DNA_ASSEMBLY_VERSION,
                "figure": str(path.absolute()),
                "filename": path.name,
                "format": path.suffix.lstrip(".").lower(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "source_notebook": source_notebook,
                "status": "passed",
            }
        )
    frame = pd.DataFrame(rows)
    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() == ".json":
        output.write_text(frame.to_json(orient="records", indent=2) + "\n")
    else:
        frame.to_csv(output, index=False)
    return frame
