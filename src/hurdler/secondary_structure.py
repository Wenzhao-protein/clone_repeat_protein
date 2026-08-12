"""Secondary-structure evidence for primitive repeat-module boundaries.

Sequence periodicity alone cannot distinguish a biological repeat from a
harmonic of that repeat.  This module maps residue-level DSSP annotations (or
an author-published residue template) onto the complete construct and scores
whether the same H/E/C topology recurs at each proposed boundary.

Unknown residues are represented by ``?`` and never count as agreement.
Coordinates returned to callers are one-based and inclusive.
"""

from __future__ import annotations

import json
import math
import re
import shutil
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from Bio import Align
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.DSSP import DSSP

from .periodicity import RepeatCandidate

SECONDARY_STRUCTURE_METHOD_VERSION = "dssp-joint-periodicity-v1"


@dataclass(frozen=True)
class SecondaryStructureSupport:
    """Structural support for one sequence-derived period candidate."""

    period: int
    local_start: int
    local_end: int
    known_fraction: float
    state_agreement: float
    state_enrichment: float
    transition_agreement: float
    transition_enrichment: float
    phase_conservation: float
    boundary_transition_fraction: float
    states_observed: int
    transitions_observed: int
    informative: bool
    score: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def collapse_dssp_state(state: str) -> str:
    """Collapse DSSP's eight states to helix, strand, or coil."""
    value = str(state or "-").strip() or "-"
    if value in {"H", "G", "I"}:
        return "H"
    if value in {"E", "B"}:
        return "E"
    return "C"


def _optional_dssp_float(value: object) -> float:
    """Convert DSSP numeric fields while preserving explicit missing values."""
    if value is None or str(value).strip().upper() in {"", "NA", "N/A", "NULL", "NONE"}:
        return float("nan")
    return float(value)


def _agreement_enrichment(values: list[str], shifted: list[str]) -> tuple[float, float]:
    pairs = [(left, right) for left, right in zip(values, shifted, strict=True) if "?" not in {left, right}]
    if not pairs:
        return 0.0, 0.0
    agreement = sum(left == right for left, right in pairs) / len(pairs)
    observed = [value for pair in pairs for value in pair]
    frequencies = {state: observed.count(state) / len(observed) for state in set(observed)}
    baseline = sum(value * value for value in frequencies.values())
    enrichment = max(0.0, (agreement - baseline) / max(1e-9, 1.0 - baseline))
    return float(agreement), float(enrichment)


def score_secondary_structure_period(
    secondary_structure: str,
    *,
    period: int,
    local_start: int,
    local_end: int,
) -> SecondaryStructureSupport:
    """Score repetition of H/E/C states and their transition pattern."""
    if period < 1 or local_start < 0 or local_end > len(secondary_structure):
        raise ValueError("Candidate coordinates fall outside secondary-structure sequence")
    if local_end - local_start < 2 * period:
        raise ValueError("At least two complete periods are required")
    block = secondary_structure[local_start:local_end].upper()
    if set(block) - {"H", "E", "C", "?"}:
        raise ValueError("Secondary structure must use only H, E, C, and ?")
    known_fraction = sum(value != "?" for value in block) / len(block)
    state_agreement, state_enrichment = _agreement_enrichment(
        list(block[:-period]), list(block[period:])
    )

    transitions: list[str] = []
    transition_count = 0
    for index in range(1, len(block)):
        left, right = block[index - 1], block[index]
        if "?" in {left, right}:
            transitions.append("?")
        else:
            changed = left != right
            transitions.append("1" if changed else "0")
            transition_count += int(changed)
    if len(transitions) > period:
        transition_agreement, transition_enrichment = _agreement_enrichment(
            transitions[:-period], transitions[period:]
        )
    else:
        transition_agreement, transition_enrichment = 0.0, 0.0

    copy_count = len(block) // period
    phase_values: list[float] = []
    for phase in range(period):
        column = [block[phase + copy_index * period] for copy_index in range(copy_count)]
        column = [value for value in column if value != "?"]
        if len(column) >= 2:
            phase_values.append(max(column.count(state) for state in set(column)) / len(column))
    phase_conservation = float(np.mean(phase_values)) if phase_values else 0.0

    boundary_values: list[bool] = []
    for offset in range(period, copy_count * period, period):
        left, right = block[offset - 1], block[offset]
        if "?" not in {left, right}:
            boundary_values.append(left != right)
    boundary_transition_fraction = (
        sum(boundary_values) / len(boundary_values) if boundary_values else 0.0
    )
    states_observed = len(set(block) - {"?"})
    transitions_per_copy = transition_count / max(1, copy_count)
    informative = bool(
        known_fraction >= 0.70
        and states_observed >= 2
        and transitions_per_copy >= 1.0
        and len(phase_values) >= max(3, period // 2)
    )
    score = (
        0.27 * state_enrichment
        + 0.27 * transition_enrichment
        + 0.20 * phase_conservation
        + 0.16 * boundary_transition_fraction
        + 0.10 * known_fraction
    )
    if not informative:
        score *= 0.5
    return SecondaryStructureSupport(
        period=period,
        local_start=local_start,
        local_end=local_end,
        known_fraction=float(known_fraction),
        state_agreement=state_agreement,
        state_enrichment=state_enrichment,
        transition_agreement=transition_agreement,
        transition_enrichment=transition_enrichment,
        phase_conservation=phase_conservation,
        boundary_transition_fraction=float(boundary_transition_fraction),
        states_observed=states_observed,
        transitions_observed=transition_count,
        informative=informative,
        score=float(max(0.0, min(1.0, score))),
    )


def score_candidates_with_secondary_structure(
    candidates: Iterable[RepeatCandidate], secondary_structure: str
) -> list[SecondaryStructureSupport]:
    return [
        score_secondary_structure_period(
            secondary_structure,
            period=candidate.period,
            local_start=candidate.local_start,
            local_end=candidate.local_end,
        )
        for candidate in candidates
    ]


def select_joint_sequence_structure_candidate(
    candidates: list[RepeatCandidate],
    supports: list[SecondaryStructureSupport],
    *,
    prior_period: int | None,
) -> tuple[RepeatCandidate | None, SecondaryStructureSupport | None, str]:
    """Choose the smallest harmonic supported by both sequence and structure.

    The H/E/C alphabet is intentionally not allowed to create a repeat on its
    own.  A structural call must coincide with a sequence candidate and span
    at least four copies.  Strong structure permits sequence variability, but
    still requires BLOSUM-positive and spectral evidence.
    """
    if len(candidates) != len(supports):
        raise ValueError("Sequence candidates and structural supports must be paired")
    paired = list(zip(candidates, supports, strict=True))
    eligible: list[tuple[RepeatCandidate, SecondaryStructureSupport, float]] = []
    for candidate, support in paired:
        harmonic = bool(
            prior_period is None
            or (
                candidate.period <= prior_period
                and abs(prior_period / candidate.period - round(prior_period / candidate.period))
                < 0.08
            )
        )
        sequence_supported = bool(
            candidate.score >= 0.55
            and candidate.adjacent_positive_fraction >= 0.55
            and candidate.spectral_concentration >= 0.52
        )
        structure_supported = bool(
            support.informative
            and support.score >= 0.62
            and support.state_agreement >= 0.75
            and support.transition_agreement >= 0.65
            and support.phase_conservation >= 0.75
        )
        enough_span = bool(
            candidate.repeat_count >= 4
            and (
                prior_period is None
                or candidate.period * candidate.repeat_count >= 2 * prior_period
            )
        )
        if harmonic and sequence_supported and structure_supported and enough_span:
            joint = 0.55 * candidate.score + 0.45 * support.score
            eligible.append((candidate, support, joint))
    if eligible:
        candidate, support, joint = min(
            eligible, key=lambda value: (value[0].period, -value[2])
        )
        return candidate, support, (
            "smallest harmonic jointly supported by amino-acid and residue-level "
            f"secondary-structure periodicity (joint_score={joint:.4f})"
        )
    return None, None, (
        "no harmonic jointly supported by sequence and residue-level secondary structure; "
        "source-annotated unit retained"
    )


def _align_annotations_to_full_sequence(
    full_sequence: str, structure_sequence: str, annotations: list[dict[str, object]]
) -> tuple[str, list[dict[str, object]], float]:
    if len(structure_sequence) != len(annotations):
        raise ValueError("Structure sequence and residue annotations have different lengths")
    if structure_sequence == full_sequence:
        mapped = []
        for index, annotation in enumerate(annotations):
            mapped.append({**annotation, "full_sequence_index": index, "full_sequence_position": index + 1})
        return "".join(str(value["ss3"]) for value in mapped), mapped, 1.0
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(full_sequence, structure_sequence)[0]
    ss = ["?"] * len(full_sequence)
    mapped: list[dict[str, object]] = []
    matches = 0
    aligned_count = 0
    for (full_start, full_end), (structure_start, structure_end) in zip(
        alignment.aligned[0], alignment.aligned[1], strict=True
    ):
        block_length = min(full_end - full_start, structure_end - structure_start)
        for offset in range(block_length):
            full_index = int(full_start + offset)
            structure_index = int(structure_start + offset)
            annotation = annotations[structure_index]
            ss[full_index] = str(annotation["ss3"])
            mapped.append(
                {
                    **annotation,
                    "full_sequence_index": full_index,
                    "full_sequence_position": full_index + 1,
                }
            )
            aligned_count += 1
            matches += full_sequence[full_index] == structure_sequence[structure_index]
    identity = matches / aligned_count if aligned_count else 0.0
    return "".join(ss), mapped, float(identity)


def dssp_annotations_for_structure(
    structure_path: str | Path,
    full_sequence: str,
    *,
    chain_id: str | None = "A",
    dssp_executable: str | Path | None = None,
    minimum_sequence_identity: float = 0.75,
) -> tuple[str, pd.DataFrame, dict[str, object]]:
    """Run DSSP and align its residue annotations to a complete sequence."""
    path = Path(structure_path)
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
        file_type = "MMCIF"
    else:
        parser = PDBParser(QUIET=True)
        file_type = "PDB"
    structure = parser.get_structure(path.stem, path)
    executable = str(
        dssp_executable
        or shutil.which("mkdssp")
        or shutil.which("dssp")
        or "mkdssp"
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="parse error at line 1:.*")
        dssp = DSSP(structure[0], str(path), dssp=executable, file_type=file_type)
    annotations_by_chain: dict[str, list[dict[str, object]]] = {}
    for key in dssp.keys():
        chain, residue_id = key
        value = dssp[key]
        raw_state = str(value[2] or "-")
        annotations_by_chain.setdefault(str(chain), []).append(
            {
                "structure_chain": str(chain),
                "auth_residue_number": int(residue_id[1]),
                "insertion_code": str(residue_id[2]).strip(),
                "amino_acid": str(value[1]).upper(),
                "dssp8": raw_state,
                "ss3": collapse_dssp_state(raw_state),
                "relative_accessibility": _optional_dssp_float(value[3]),
                "phi": _optional_dssp_float(value[4]),
                "psi": _optional_dssp_float(value[5]),
            }
        )
    if chain_id:
        requested = str(chain_id)
        if requested not in annotations_by_chain:
            raise ValueError(f"DSSP returned no residues for chain {requested} in {path}")
        candidate_chains = [requested]
    else:
        candidate_chains = sorted(annotations_by_chain)
    if not candidate_chains:
        raise ValueError(f"DSSP returned no protein chains in {path}")
    aligned_options = []
    for candidate_chain in candidate_chains:
        annotations = annotations_by_chain[candidate_chain]
        structure_sequence = "".join(str(value["amino_acid"]) for value in annotations)
        ss_value, mapped_value, identity_value = _align_annotations_to_full_sequence(
            full_sequence, structure_sequence, annotations
        )
        aligned_options.append(
            (identity_value, len(mapped_value) / len(full_sequence), candidate_chain, ss_value, mapped_value)
        )
    identity, _, selected_chain, ss, mapped = max(aligned_options)
    if identity < minimum_sequence_identity:
        raise ValueError(
            f"DSSP/full-sequence identity {identity:.3f} is below "
            f"{minimum_sequence_identity:.3f}; probable chain mismatch"
        )
    frame = pd.DataFrame(mapped)
    metadata = {
        "secondary_structure_method": "DSSP 4 via Biopython",
        "secondary_structure_method_version": SECONDARY_STRUCTURE_METHOD_VERSION,
        "secondary_structure_source_path": str(path),
        "secondary_structure_chain": selected_chain,
        "secondary_structure_known_fraction": (len(frame) / len(full_sequence)),
        "structure_sequence_identity": identity,
        "dssp_executable": executable,
    }
    return ss, frame, metadata


def parse_dhr_secondary_structure_table(text: str) -> pd.DataFrame:
    """Parse DHR repeat-length and H/L residue counts from Supplementary Table 2."""
    start_marker = "Supplementary Table 2 | Global geometric parameters."
    end_marker = "For each design, repeat size, secondary structure length"
    if start_marker not in text or end_marker not in text:
        raise ValueError("DHR Supplementary Table 2 markers were not found")
    table = text.split(start_marker, 1)[1].split(end_marker, 1)[0]
    rows: list[dict[str, object]] = []
    pattern = re.compile(r"^\s*(\d{1,2})\s+(\d{2})\s+(\d{1,2})\s+(\d)\s+(\d{1,2})\s+(\d)\b")
    for line in table.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        number, length, helix1, loop1, helix2, loop2 = map(int, match.groups())
        if helix1 + loop1 + helix2 + loop2 != length:
            raise ValueError(f"DHR{number} secondary-structure lengths do not sum to repeat length")
        ss3 = "H" * helix1 + "C" * loop1 + "H" * helix2 + "C" * loop2
        rows.append(
            {
                "source_accession": f"DHR{number}",
                "author_repeat_length": length,
                "author_helix1_length": helix1,
                "author_loop1_length": loop1,
                "author_helix2_length": helix2,
                "author_loop2_length": loop2,
                "author_repeat_ss3": ss3,
            }
        )
    frame = pd.DataFrame(rows).sort_values("source_accession").reset_index(drop=True)
    if len(frame) != 83 or frame.source_accession.nunique() != 83:
        raise ValueError(f"Expected 83 DHR topology rows, found {len(frame)}")
    return frame


def map_author_template_to_full_sequence(
    full_sequence: str,
    *,
    repeat_start: int,
    full_sequence_origin: int,
    ss3_template: str,
) -> tuple[str, pd.DataFrame]:
    """Tile a published residue-level repeat template around a known boundary."""
    period = len(ss3_template)
    anchor = repeat_start - full_sequence_origin
    ss = ["?"] * len(full_sequence)
    rows: list[dict[str, object]] = []
    first_copy = math.floor(-anchor / period)
    last_copy = math.ceil((len(full_sequence) - anchor) / period)
    for copy_index in range(first_copy, last_copy + 1):
        start = anchor + copy_index * period
        for phase, state in enumerate(ss3_template):
            index = start + phase
            if 0 <= index < len(full_sequence):
                ss[index] = state
                rows.append(
                    {
                        "full_sequence_index": index,
                        "full_sequence_position": full_sequence_origin + index,
                        "amino_acid": full_sequence[index],
                        "ss3": state,
                        "dssp8": None,
                        "template_copy_index": copy_index,
                        "template_phase": phase + 1,
                    }
                )
    return "".join(ss), pd.DataFrame(rows).sort_values("full_sequence_index")


def secondary_structure_json(secondary_structure: str) -> str:
    """Stable JSON representation used in Parquet/CSV catalog fields."""
    return json.dumps({"alphabet": "HEC?", "states": secondary_structure}, sort_keys=True)
