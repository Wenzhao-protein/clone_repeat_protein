"""Strict DSSP/Foldseek boundary inference for designed repeat proteins."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any, Iterable

import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB.Polypeptide import protein_letters_3to1_extended

from .constants import validate_protein_sequence
from .periodicity import MODULE_SELECTION_POLICY

DESIGNED_BOUNDARY_METHOD = "biotite-dssp-foldseek-strict-dual-evidence-v1"
DESIGNED_CORPUS_VERSION = "expanded-middle-repeatsdb-foldseek-v1"
DEFAULT_FOLDSEEK = Path("/net/software/utils/foldseek")

_BLOSUM62 = substitution_matrices.load("BLOSUM62")


@dataclass(frozen=True)
class DualEvidenceThresholds:
    dssp_state_agreement: float = 0.70
    dssp_transition_agreement: float = 0.70
    minimum_transitions_per_unit: int = 2
    foldseek_3di_identity: float = 0.35
    foldseek_tm_score: float = 0.50
    foldseek_lddt: float = 0.50
    foldseek_coverage: float = 0.70
    minimum_complete_copies: int = 2


def _pairwise_agreement(sequence: str, lag: int, unknown: str = "?X-") -> float:
    pairs = [
        (left, right)
        for left, right in zip(sequence[:-lag], sequence[lag:], strict=True)
        if left not in unknown and right not in unknown
    ]
    return (
        sum(left == right for left, right in pairs) / len(pairs)
        if pairs
        else 0.0
    )


def _transition_string(sequence: str) -> str:
    if not sequence:
        return ""
    return "".join(
        "1" if sequence[index] != sequence[index - 1] else "0"
        for index in range(1, len(sequence))
    )


def _median_transitions_per_unit(dssp8: str, period: int) -> float:
    counts = [
        sum(unit[index] != unit[index - 1] for index in range(1, len(unit)))
        for start in range(0, len(dssp8) - period + 1, period)
        if len(unit := dssp8[start : start + period]) == period
    ]
    return float(median(counts)) if counts else 0.0


def _blosum_positive_fraction(sequence: str, lag: int) -> float:
    values: list[bool] = []
    for left, right in zip(sequence[:-lag], sequence[lag:], strict=True):
        try:
            values.append(float(_BLOSUM62[left, right]) > 0)
        except (IndexError, KeyError):
            continue
    return sum(values) / len(values) if values else 0.0


def lag_evidence(
    amino_acids: str,
    dssp8: str,
    foldseek_3di: str,
    period: int,
) -> dict[str, float | int]:
    if not (len(amino_acids) == len(dssp8) == len(foldseek_3di)):
        raise ValueError("AA, DSSP and Foldseek 3Di strings must have equal length")
    if not 1 <= period <= len(amino_acids) // 2:
        raise ValueError("Candidate period must allow at least two copies")
    dssp_transitions = _transition_string(dssp8)
    return {
        "period": period,
        "dssp_state_agreement": _pairwise_agreement(dssp8, period),
        "dssp_transition_agreement": _pairwise_agreement(
            dssp_transitions, period
        ),
        "dssp_median_transitions_per_unit": _median_transitions_per_unit(
            dssp8, period
        ),
        "foldseek_3di_identity": _pairwise_agreement(foldseek_3di, period),
        "aa_identity": _pairwise_agreement(amino_acids, period, unknown=""),
        "aa_blosum_positive_fraction": _blosum_positive_fraction(
            amino_acids, period
        ),
    }


def _local_peak_periods(frame: pd.DataFrame, column: str) -> list[int]:
    values = dict(zip(frame.period.astype(int), frame[column].astype(float), strict=True))
    return [
        period
        for period in sorted(values)
        if values[period] >= values.get(period - 1, -math.inf)
        and values[period] >= values.get(period + 1, -math.inf)
    ]


def scan_dual_evidence_periods(
    amino_acids: str,
    dssp8: str,
    foldseek_3di: str,
    *,
    minimum_period: int = 6,
    thresholds: DualEvidenceThresholds = DualEvidenceThresholds(),
) -> pd.DataFrame:
    """Return candidate lags whose independent DSSP and 3Di peaks agree."""
    maximum = len(amino_acids) // 2
    rows = [
        lag_evidence(amino_acids, dssp8, foldseek_3di, period)
        for period in range(minimum_period, maximum + 1)
    ]
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    dssp_peaks = _local_peak_periods(frame, "dssp_transition_agreement")
    foldseek_peaks = _local_peak_periods(frame, "foldseek_3di_identity")
    matched: set[int] = set()
    for dssp_period in dssp_peaks:
        for foldseek_period in foldseek_peaks:
            tolerance = max(2, round(0.05 * min(dssp_period, foldseek_period)))
            if abs(dssp_period - foldseek_period) > tolerance:
                continue
            lo, hi = sorted((dssp_period, foldseek_period))
            window = frame.loc[frame.period.between(lo, hi)].copy()
            window["_joint"] = (
                window.dssp_state_agreement
                * window.dssp_transition_agreement
                * window.foldseek_3di_identity
            )
            matched.add(
                int(window.sort_values(["_joint", "period"], ascending=[False, True]).iloc[0].period)
            )
    frame["dssp_peak"] = frame.period.isin(dssp_peaks)
    frame["foldseek_3di_peak"] = frame.period.isin(foldseek_peaks)
    frame["dual_peak_matched"] = frame.period.isin(matched)
    frame["lag_thresholds_passed"] = (
        frame.dual_peak_matched
        & frame.dssp_state_agreement.ge(thresholds.dssp_state_agreement)
        & frame.dssp_transition_agreement.ge(thresholds.dssp_transition_agreement)
        & frame.dssp_median_transitions_per_unit.ge(
            thresholds.minimum_transitions_per_unit
        )
        & frame.foldseek_3di_identity.ge(thresholds.foldseek_3di_identity)
    )
    frame["joint_lag_score"] = (
        frame.dssp_state_agreement
        * frame.dssp_transition_agreement
        * frame.foldseek_3di_identity
    ) ** (1 / 3)
    return frame.sort_values(
        ["lag_thresholds_passed", "period", "joint_lag_score"],
        ascending=[False, True, False],
        kind="mergesort",
    ).reset_index(drop=True)


def _unit_pair_string_metrics(
    dssp8: str, foldseek_3di: str, left: int, period: int
) -> dict[str, float]:
    right = left + period
    dssp_left, dssp_right = dssp8[left:right], dssp8[right : right + period]
    di_left, di_right = foldseek_3di[left:right], foldseek_3di[right : right + period]
    return {
        "dssp": _pairwise_agreement(dssp_left + dssp_right, period),
        "transition": _pairwise_agreement(
            _transition_string(dssp_left) + _transition_string(dssp_right),
            max(1, period - 1),
        ),
        "three_di": _pairwise_agreement(di_left + di_right, period),
    }


def select_maximal_repeat_block(
    dssp8: str,
    foldseek_3di: str,
    period: int,
    *,
    thresholds: DualEvidenceThresholds = DualEvidenceThresholds(),
) -> dict[str, Any] | None:
    """Select the longest phase-consistent chain of passing adjacent copies."""
    candidates: list[dict[str, Any]] = []
    for phase in range(period):
        starts = list(range(phase, len(dssp8) - 2 * period + 1, period))
        run: list[int] = []
        for start in starts:
            metrics = _unit_pair_string_metrics(dssp8, foldseek_3di, start, period)
            passed = (
                metrics["dssp"] >= thresholds.dssp_state_agreement
                and metrics["transition"] >= thresholds.dssp_transition_agreement
                and metrics["three_di"] >= thresholds.foldseek_3di_identity
            )
            if passed:
                run.append(start)
            else:
                if run:
                    block_start = run[0]
                    candidates.append(
                        {
                            "start": block_start,
                            "end": run[-1] + 2 * period,
                            "copy_count": len(run) + 1,
                            "phase": phase,
                            "pair_metrics": [
                                _unit_pair_string_metrics(
                                    dssp8, foldseek_3di, value, period
                                )
                                for value in run
                            ],
                        }
                    )
                    run = []
        if run:
            candidates.append(
                {
                    "start": run[0],
                    "end": run[-1] + 2 * period,
                    "copy_count": len(run) + 1,
                    "phase": phase,
                    "pair_metrics": [
                        _unit_pair_string_metrics(dssp8, foldseek_3di, value, period)
                        for value in run
                    ],
                }
            )
    if not candidates:
        return None
    for candidate in candidates:
        candidate["string_joint_score"] = float(
            np.mean(
                [
                    (item["dssp"] * item["transition"] * item["three_di"])
                    ** (1 / 3)
                    for item in candidate["pair_metrics"]
                ]
            )
        )
    return sorted(
        candidates,
        key=lambda value: (
            -int(value["copy_count"]),
            -float(value["string_joint_score"]),
            int(value["start"]),
        ),
    )[0]


def _biotite_chain_and_dssp(
    structure_path: Path,
    *,
    chain_id: str | None,
    dssp_executable: Path,
) -> tuple[Any, str, str, str]:
    try:
        import biotite.structure as struc
        from biotite.application.dssp import DsspApp
        from biotite.structure.io import load_structure
    except ImportError as exc:
        raise RuntimeError("Biotite is required for designed boundary inference") from exc

    atoms = load_structure(str(structure_path), model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    if atoms.array_length() == 0:
        raise ValueError(f"No amino-acid atoms in {structure_path}")
    chain_ids = sorted(set(str(value) for value in atoms.chain_id))
    if chain_id and chain_id in chain_ids:
        selected_chain = chain_id
    else:
        selected_chain = max(
            chain_ids,
            key=lambda value: struc.get_residue_count(atoms[atoms.chain_id == value]),
        )
    chain = atoms[atoms.chain_id == selected_chain]
    residue_starts = struc.get_residue_starts(chain, add_exclusive_stop=True)
    complete_backbone = []
    for index in range(len(residue_starts) - 1):
        atom_names = set(
            str(value)
            for value in chain.atom_name[
                residue_starts[index] : residue_starts[index + 1]
            ]
        )
        complete_backbone.append({"N", "CA", "C", "O"}.issubset(atom_names))
    if not any(complete_backbone):
        raise ValueError("No residue has a complete N/CA/C/O backbone for DSSP")
    first_complete = complete_backbone.index(True)
    last_complete = len(complete_backbone) - 1 - complete_backbone[::-1].index(True)
    if not all(complete_backbone[first_complete : last_complete + 1]):
        raise ValueError("Internal residues lack a complete backbone for DSSP")
    # mkdssp legitimately omits incomplete terminal residues from its output.
    # Trim only such termini so the eight-state string remains residue-mapped.
    chain = chain[
        int(residue_starts[first_complete]) : int(residue_starts[last_complete + 1])
    ]
    _, residue_names = struc.get_residues(chain)
    amino_acids = "".join(
        protein_letters_3to1_extended.get(str(name).upper(), "X")
        for name in residue_names
    )
    dssp_chain = chain.copy()
    # Long mmCIF asym IDs can make mkdssp choose an incompatible legacy DSSP
    # representation.  Relabel the already-selected single chain solely for
    # the Biotite DsspApp call; residue order and coordinates are unchanged.
    dssp_chain.chain_id[:] = "A"
    dssp8 = "".join(
        str(value) if str(value) not in {"", "-"} else "C"
        for value in DsspApp.annotate_sse(
            dssp_chain, bin_path=str(dssp_executable)
        )
    )
    if len(amino_acids) != len(dssp8):
        raise ValueError("Biotite residue and DSSP lengths differ")
    return chain, validate_protein_sequence(amino_acids), dssp8, selected_chain


def foldseek_descriptor(
    structure_path: str | Path,
    *,
    foldseek_binary: str | Path = DEFAULT_FOLDSEEK,
    chain_id: str | None = None,
) -> tuple[str, str, str]:
    """Return the exact Foldseek AA and 3Di strings for one chain."""
    executable = Path(foldseek_binary)
    if not executable.is_file():
        raise FileNotFoundError(f"Foldseek binary not found: {executable}")
    with tempfile.TemporaryDirectory(prefix="hurdler_foldseek_descriptor_") as temporary:
        output = Path(temporary) / "descriptor"
        subprocess.run(
            [
                str(executable),
                "structureto3didescriptor",
                str(structure_path),
                str(output),
                "--threads",
                "1",
                "--chain-name-mode",
                "1",
                "-v",
                "0",
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=180,
        )
        rows = []
        for line in output.read_text().splitlines():
            values = line.rstrip("\0").split("\t")
            if len(values) >= 3:
                rows.append((values[0], values[1], values[2]))
    if not rows:
        raise ValueError(f"Foldseek returned no descriptor for {structure_path}")
    if chain_id:
        matched = [row for row in rows if row[0].rsplit("_", 1)[-1] == chain_id]
        if matched:
            rows = matched
    name, amino_acids, three_di = max(rows, key=lambda value: len(value[1]))
    return validate_protein_sequence(amino_acids), three_di, name


def project_foldseek_descriptor(
    biotite_amino_acids: str,
    foldseek_amino_acids: str,
    foldseek_3di: str,
    *,
    minimum_coverage: float = 0.98,
    minimum_identity: float = 0.97,
) -> tuple[str, dict[str, float | int]]:
    """Project a Foldseek descriptor onto the Biotite residue coordinates.

    Foldseek can retain a few terminal/linker residues that Biotite excludes
    from the amino-acid atom array.  A global residue alignment makes that
    coordinate difference explicit.  Descriptor gaps are represented by
    ``?`` and are ignored by agreement calculations; low-coverage or
    sequence-discordant structures remain strict failures.
    """
    if len(foldseek_amino_acids) != len(foldseek_3di):
        raise ValueError("Foldseek AA and 3Di descriptor lengths differ")
    if biotite_amino_acids == foldseek_amino_acids:
        return foldseek_3di, {
            "foldseek_biotite_alignment_coverage": 1.0,
            "foldseek_biotite_alignment_identity": 1.0,
            "foldseek_biotite_unmapped_residues": 0,
        }
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -1.0
    alignment = aligner.align(biotite_amino_acids, foldseek_amino_acids)[0]
    projected = ["?"] * len(biotite_amino_acids)
    mapped = 0
    identical = 0
    for (biotite_start, biotite_end), (foldseek_start, foldseek_end) in zip(
        alignment.aligned[0], alignment.aligned[1], strict=True
    ):
        biotite_length = int(biotite_end - biotite_start)
        foldseek_length = int(foldseek_end - foldseek_start)
        if biotite_length != foldseek_length:
            raise ValueError("Unequal aligned residue blocks in descriptor projection")
        for offset in range(biotite_length):
            biotite_index = int(biotite_start) + offset
            foldseek_index = int(foldseek_start) + offset
            projected[biotite_index] = foldseek_3di[foldseek_index]
            mapped += 1
            identical += (
                biotite_amino_acids[biotite_index]
                == foldseek_amino_acids[foldseek_index]
            )
    coverage = mapped / len(biotite_amino_acids) if biotite_amino_acids else 0.0
    identity = identical / mapped if mapped else 0.0
    if coverage < minimum_coverage or identity < minimum_identity:
        raise ValueError(
            "Biotite/Foldseek residue alignment is inadequate "
            f"(coverage={coverage:.4f}, identity={identity:.4f})"
        )
    return "".join(projected), {
        "foldseek_biotite_alignment_coverage": coverage,
        "foldseek_biotite_alignment_identity": identity,
        "foldseek_biotite_unmapped_residues": len(biotite_amino_acids) - mapped,
    }


def map_structure_positions_to_full_sequence(
    structure_amino_acids: str,
    full_amino_acids: str,
    *,
    minimum_coverage: float = 0.98,
    minimum_identity: float = 0.97,
) -> tuple[list[int | None], dict[str, float | int]]:
    """Map observed structure residues to zero-based full-sequence positions."""
    if structure_amino_acids == full_amino_acids:
        return list(range(len(structure_amino_acids))), {
            "structure_full_alignment_coverage": 1.0,
            "structure_full_alignment_identity": 1.0,
            "structure_full_unmapped_residues": 0,
        }
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -1.0
    alignment = aligner.align(structure_amino_acids, full_amino_acids)[0]
    mapping: list[int | None] = [None] * len(structure_amino_acids)
    mapped = 0
    identical = 0
    for (structure_start, structure_end), (full_start, full_end) in zip(
        alignment.aligned[0], alignment.aligned[1], strict=True
    ):
        structure_length = int(structure_end - structure_start)
        full_length = int(full_end - full_start)
        if structure_length != full_length:
            raise ValueError("Unequal aligned blocks in full-sequence mapping")
        for offset in range(structure_length):
            structure_index = int(structure_start) + offset
            full_index = int(full_start) + offset
            mapping[structure_index] = full_index
            mapped += 1
            identical += (
                structure_amino_acids[structure_index]
                == full_amino_acids[full_index]
            )
    coverage = mapped / len(structure_amino_acids) if structure_amino_acids else 0.0
    identity = identical / mapped if mapped else 0.0
    if coverage < minimum_coverage or identity < minimum_identity:
        raise ValueError(
            "Structure cannot be reliably mapped into full designed sequence "
            f"(coverage={coverage:.4f}, identity={identity:.4f})"
        )
    return mapping, {
        "structure_full_alignment_coverage": coverage,
        "structure_full_alignment_identity": identity,
        "structure_full_unmapped_residues": len(structure_amino_acids) - mapped,
    }


def _save_unit_fragments(chain: Any, block: dict[str, Any], period: int, directory: Path) -> list[Path]:
    import biotite.structure as struc
    from biotite.structure.io import save_structure

    directory.mkdir(parents=True, exist_ok=True)
    residue_starts = struc.get_residue_starts(chain, add_exclusive_stop=True)
    paths: list[Path] = []
    for index in range(int(block["copy_count"])):
        residue_start = int(block["start"]) + index * period
        residue_end = residue_start + period
        atom_start = int(residue_starts[residue_start])
        atom_end = int(residue_starts[residue_end])
        fragment = chain[atom_start:atom_end].copy()
        fragment.chain_id[:] = "A"
        path = directory / f"unit_{index:03d}.pdb"
        save_structure(str(path), fragment)
        paths.append(path)
    return paths


def validate_foldseek_fragments(
    chain: Any,
    block: dict[str, Any],
    period: int,
    *,
    foldseek_binary: str | Path = DEFAULT_FOLDSEEK,
) -> dict[str, Any]:
    """Validate finalized adjacent unit structures with global TM-align."""
    with tempfile.TemporaryDirectory(prefix="hurdler_foldseek_units_") as temporary:
        root = Path(temporary)
        fragments = _save_unit_fragments(chain, block, period, root / "units")
        output = root / "alignments.tsv"
        subprocess.run(
            [
                str(foldseek_binary),
                "easy-search",
                str(root / "units"),
                str(root / "units"),
                str(output),
                str(root / "work"),
                "--exhaustive-search",
                "1",
                "--alignment-type",
                "1",
                "--max-seqs",
                str(max(100, len(fragments) ** 2)),
                "-e",
                "1000",
                "--format-output",
                "query,target,qtmscore,ttmscore,alntmscore,lddt,alnlen,qlen,tlen,qstart,qend,tstart,tend",
                "--threads",
                "1",
                "-v",
                "1",
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=300,
        )
        columns = [
            "query",
            "target",
            "qtmscore",
            "ttmscore",
            "alntmscore",
            "lddt",
            "alnlen",
            "qlen",
            "tlen",
            "qstart",
            "qend",
            "tstart",
            "tend",
        ]
        alignments = pd.read_csv(output, sep="\t", names=columns)
    pair_rows = []
    for index in range(len(fragments) - 1):
        left, right = f"unit_{index:03d}", f"unit_{index + 1:03d}"
        matched = alignments.loc[
            alignments["query"].astype(str).str.contains(left, regex=False)
            & alignments["target"].astype(str).str.contains(right, regex=False)
        ]
        if matched.empty:
            matched = alignments.loc[
                alignments["query"].astype(str).str.contains(right, regex=False)
                & alignments["target"].astype(str).str.contains(left, regex=False)
            ]
        if matched.empty:
            continue
        value = matched.sort_values("alntmscore", ascending=False).iloc[0]
        pair_rows.append(
            {
                "left_index": index,
                "right_index": index + 1,
                "qtmscore": float(value.qtmscore),
                "ttmscore": float(value.ttmscore),
                "alntmscore": float(value.alntmscore),
                "lddt": float(value.lddt),
                "coverage": float(value.alnlen) / min(float(value.qlen), float(value.tlen)),
            }
        )
    if len(pair_rows) != len(fragments) - 1:
        raise ValueError("Foldseek did not return every adjacent unit alignment")
    return {
        "foldseek_pair_alignments": pair_rows,
        "foldseek_median_min_tm": float(
            median(min(row["qtmscore"], row["ttmscore"]) for row in pair_rows)
        ),
        "foldseek_median_lddt": float(median(row["lddt"] for row in pair_rows)),
        "foldseek_median_coverage": float(
            median(row["coverage"] for row in pair_rows)
        ),
    }


def _resolve_structure(row: dict[str, Any]) -> tuple[Path, str]:
    for column, source in (
        ("author_structure_path", "author_or_pdb"),
        ("pdb_structure_path", "author_or_pdb"),
        ("alphafold_structure_path", "alphafolddb"),
        ("af3_structure_path", "alphafold3"),
        ("structure_path", "provided"),
    ):
        value = row.get(column)
        if isinstance(value, str) and value and Path(value).is_file():
            return Path(value), source
    raise FileNotFoundError("No author/PDB, AlphaFoldDB or AlphaFold3 structure is available")


def _mafft_alignment(units: list[str], mafft_binary: str | Path) -> list[str]:
    if len(units) < 2:
        raise ValueError("At least two units are required for MAFFT")
    with tempfile.TemporaryDirectory(prefix="hurdler_mafft_") as temporary:
        fasta = Path(temporary) / "units.fasta"
        fasta.write_text(
            "".join(f">unit_{index:03d}\n{sequence}\n" for index, sequence in enumerate(units))
        )
        result = subprocess.run(
            [str(mafft_binary), "--auto", "--quiet", str(fasta)],
            check=True,
            capture_output=True,
            text=True,
            timeout=300,
        )
        aligned = Path(temporary) / "aligned.fasta"
        aligned.write_text(result.stdout)
        records = list(AlignIO.read(aligned, "fasta"))
    return [str(record.seq) for record in records]


def infer_designed_boundary(
    row: dict[str, Any],
    *,
    dssp_executable: str | Path,
    foldseek_binary: str | Path = DEFAULT_FOLDSEEK,
    mafft_binary: str | Path = "mafft",
    fixed_threshold: float = 0.8,
    thresholds: DualEvidenceThresholds = DualEvidenceThresholds(),
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Infer one designed module or return a strict, auditable exclusion."""
    payload = dict(row)
    structure_path, structure_source = _resolve_structure(payload)
    chain_id = str(payload.get("source_chain") or "") or None
    chain, dssp_aa, dssp8, selected_chain = _biotite_chain_and_dssp(
        structure_path,
        chain_id=chain_id,
        dssp_executable=Path(dssp_executable),
    )
    foldseek_aa, three_di, foldseek_name = foldseek_descriptor(
        structure_path, foldseek_binary=foldseek_binary, chain_id=selected_chain
    )
    three_di, descriptor_projection = project_foldseek_descriptor(
        dssp_aa, foldseek_aa, three_di
    )
    candidates = scan_dual_evidence_periods(
        dssp_aa, dssp8, three_di, thresholds=thresholds
    )
    # A short internal helix/strand texture can pass the pointwise thresholds
    # without being a complete repeat module.  First find the dominant
    # recurrent scale using contiguous block span and DSSP/3Di recurrence,
    # then choose the smallest fully validated period among that scale and its
    # harmonics.  This rejects the trivial local alternatives returned by
    # whole-chain self-search while preserving primitive-period selection.
    blocks: dict[int, dict[str, Any]] = {}
    for candidate in candidates.loc[candidates.lag_thresholds_passed].itertuples(
        index=False
    ):
        period = int(candidate.period)
        block = select_maximal_repeat_block(
            dssp8, three_di, period, thresholds=thresholds
        )
        if block is None or int(block["copy_count"]) < thresholds.minimum_complete_copies:
            continue
        block["recurrence_composite"] = (
            (int(block["end"]) - int(block["start"]))
            / len(dssp8)
            * float(block["string_joint_score"])
        )
        blocks[period] = block
    best_recurrence = max(
        (float(block["recurrence_composite"]) for block in blocks.values()),
        default=0.0,
    )
    shortlisted_periods = {
        period
        for period, block in blocks.items()
        if float(block["recurrence_composite"]) >= 0.85 * best_recurrence
    }
    structural_by_period: dict[int, dict[str, Any]] = {}
    accepted: tuple[pd.Series, dict[str, Any], dict[str, Any]] | None = None
    for candidate in candidates.loc[
        candidates.period.isin(shortlisted_periods)
    ].sort_values(["period", "joint_lag_score"], ascending=[True, False]).itertuples(
        index=False
    ):
        period = int(candidate.period)
        block = blocks[period]
        structural = validate_foldseek_fragments(
            chain, block, period, foldseek_binary=foldseek_binary
        )
        structural["foldseek_global_thresholds_passed"] = bool(
            structural["foldseek_median_min_tm"] >= thresholds.foldseek_tm_score
            and structural["foldseek_median_lddt"] >= thresholds.foldseek_lddt
            and structural["foldseek_median_coverage"] >= thresholds.foldseek_coverage
        )
        structural_by_period[period] = structural
        if structural["foldseek_global_thresholds_passed"]:
            accepted = (pd.Series(candidate._asdict()), block, structural)
            break
    candidate_output = candidates.copy()
    candidate_output["repeat_block_start"] = candidate_output.period.map(
        lambda value: blocks.get(int(value), {}).get("start")
    )
    candidate_output["repeat_block_end"] = candidate_output.period.map(
        lambda value: blocks.get(int(value), {}).get("end")
    )
    candidate_output["repeat_block_copy_count"] = candidate_output.period.map(
        lambda value: blocks.get(int(value), {}).get("copy_count")
    )
    candidate_output["repeat_block_string_joint_score"] = candidate_output.period.map(
        lambda value: blocks.get(int(value), {}).get("string_joint_score")
    )
    candidate_output["repeat_block_recurrence_composite"] = candidate_output.period.map(
        lambda value: blocks.get(int(value), {}).get("recurrence_composite")
    )
    candidate_output["dominant_recurrence_shortlist"] = candidate_output.period.isin(
        shortlisted_periods
    )
    for column in (
        "foldseek_median_min_tm",
        "foldseek_median_lddt",
        "foldseek_median_coverage",
        "foldseek_global_thresholds_passed",
    ):
        candidate_output[column] = candidate_output.period.map(
            lambda value, key=column: structural_by_period.get(int(value), {}).get(key)
        )
    candidate_output.insert(0, "module_id", str(payload["module_id"]))
    if accepted is None:
        payload.update(
            {
                "corpus_version": DESIGNED_CORPUS_VERSION,
                "boundary_method_version": DESIGNED_BOUNDARY_METHOD,
                "boundary_refinement_status": "boundary_ambiguous_strict_dual_evidence",
                "strict_dual_evidence_passed": False,
                "structure_path": str(structure_path),
                "structure_source_type": structure_source,
                "structure_sequence_sha256": hashlib.sha256(dssp_aa.encode()).hexdigest(),
                "foldseek_descriptor_name": foldseek_name,
            }
        )
        return payload, candidate_output, pd.DataFrame(), pd.DataFrame()

    candidate, block, structural = accepted
    period = int(candidate.period)
    block_start = int(block["start"])
    copy_count = int(block["copy_count"])
    units = [
        dssp_aa[block_start + index * period : block_start + (index + 1) * period]
        for index in range(copy_count)
    ]
    aligned_units = _mafft_alignment(units, mafft_binary)
    middle_index = (copy_count - 1) // 2
    full_sequence = validate_protein_sequence(
        str(payload.get("full_sequence") or dssp_aa)
    )
    structure_to_full, full_mapping_audit = map_structure_positions_to_full_sequence(
        dssp_aa, full_sequence
    )
    unit_coordinates: list[tuple[int, int]] = []
    for index, unit in enumerate(units):
        structure_start = block_start + index * period
        positions = structure_to_full[structure_start : structure_start + period]
        if any(position is None for position in positions):
            raise ValueError("Inferred repeat crosses an unmapped structure residue")
        mapped_positions = [int(position) for position in positions if position is not None]
        if any(
            right != left + 1
            for left, right in zip(
                mapped_positions[:-1], mapped_positions[1:], strict=True
            )
        ):
            raise ValueError("Inferred repeat crosses an insertion in the full sequence")
        full_unit = full_sequence[mapped_positions[0] : mapped_positions[-1] + 1]
        if full_unit != unit:
            raise ValueError("Inferred structure repeat differs from the full sequence")
        unit_coordinates.append((mapped_positions[0] + 1, mapped_positions[-1] + 1))
    selected_start, selected_end = unit_coordinates[middle_index]
    selected_sequence = units[middle_index]

    unit_rows = []
    for index, (unit, aligned) in enumerate(zip(units, aligned_units, strict=True)):
        start, end = unit_coordinates[index]
        unit_rows.append(
            {
                "module_id": payload["module_id"],
                "repeat_index": index + 1,
                "unit_start": start,
                "unit_end": end,
                "unit_sequence": unit,
                "aligned_unit_sequence": aligned,
                "is_selected_module": index == middle_index,
                "boundary_source": DESIGNED_BOUNDARY_METHOD,
            }
        )
    position_rows = []
    for position, column in enumerate(zip(*aligned_units, strict=True), start=1):
        counts = {value: column.count(value) for value in sorted(set(column))}
        consensus, count = sorted(
            counts.items(), key=lambda item: (-item[1], item[0] == "-", item[0])
        )[0]
        conservation = count / len(aligned_units)
        position_rows.append(
            {
                "module_id": payload["module_id"],
                "module_position": position,
                "consensus_amino_acid": consensus,
                "conservation": conservation,
                "fixed": conservation >= fixed_threshold,
                "variants_json": json.dumps(sorted(counts)),
            }
        )
    variable_positions = [
        row["module_position"] for row in position_rows if not row["fixed"]
    ]
    variable_ranges: list[list[int]] = []
    for position in variable_positions:
        if not variable_ranges or position != variable_ranges[-1][1] + 1:
            variable_ranges.append([position, position])
        else:
            variable_ranges[-1][1] = position

    payload.update(
        {
            "collection": "designed_all",
            "module_type": "Designed",
            "corpus_version": DESIGNED_CORPUS_VERSION,
            "boundary_method": DESIGNED_BOUNDARY_METHOD,
            "boundary_method_version": DESIGNED_BOUNDARY_METHOD,
            "boundary_refinement_status": "strict_dual_evidence_passed",
            "strict_dual_evidence_passed": True,
            "unit_sequence": selected_sequence,
            "unit_length": period,
            "unit_start": selected_start,
            "unit_end": selected_end,
            "selected_module_sequence": selected_sequence,
            "selected_module_start": selected_start,
            "selected_module_end": selected_end,
            "selected_module_index": middle_index + 1,
            "selected_module_count": copy_count,
            "selected_module_policy": MODULE_SELECTION_POLICY,
            "module_selection_policy": MODULE_SELECTION_POLICY,
            "repeat_region_start": unit_coordinates[0][0],
            "repeat_region_end": unit_coordinates[-1][1],
            "repeat_count": copy_count,
            "period": period,
            "primitive_period": period,
            "full_sequence": full_sequence,
            "full_sequence_origin": 1,
            "full_sequence_sha256": hashlib.sha256(full_sequence.encode()).hexdigest(),
            "unit_sequences_json": json.dumps(units),
            "aligned_unit_sequences_json": json.dumps(aligned_units),
            "source_unit_coordinates_json": json.dumps(
                [[row["unit_start"], row["unit_end"]] for row in unit_rows]
            ),
            "fixed_positions_json": json.dumps(
                [row["module_position"] for row in position_rows if row["fixed"]]
            ),
            "variable_positions_json": json.dumps(variable_positions),
            "variable_ranges_json": json.dumps(variable_ranges),
            "dssp_state_agreement": float(candidate.dssp_state_agreement),
            "dssp_transition_agreement": float(candidate.dssp_transition_agreement),
            "dssp_median_transitions_per_unit": float(
                candidate.dssp_median_transitions_per_unit
            ),
            "foldseek_3di_identity": float(candidate.foldseek_3di_identity),
            **structural,
            "foldseek_pair_alignments_json": json.dumps(
                structural["foldseek_pair_alignments"], sort_keys=True
            ),
            "structure_path": str(structure_path),
            "structure_source_type": structure_source,
            "structure_chain": selected_chain,
            "structure_sequence_sha256": hashlib.sha256(dssp_aa.encode()).hexdigest(),
            "foldseek_descriptor_name": foldseek_name,
            **descriptor_projection,
            **full_mapping_audit,
            "selection_reason": (
                "smallest complete period passing independent DSSP, Foldseek 3Di, "
                "adjacent-copy TM-score, LDDT and coverage thresholds"
            ),
        }
    )
    payload.pop("foldseek_pair_alignments", None)
    return (
        payload,
        candidate_output,
        pd.DataFrame(unit_rows),
        pd.DataFrame(position_rows),
    )


def infer_designed_catalog(
    input_path: str | Path,
    output_path: str | Path,
    *,
    candidates_path: str | Path,
    units_path: str | Path,
    positions_path: str | Path,
    exclusions_path: str | Path,
    dssp_executable: str | Path,
    foldseek_binary: str | Path = DEFAULT_FOLDSEEK,
    mafft_binary: str | Path = "mafft",
    shard_index: int = 0,
    shard_count: int = 1,
) -> pd.DataFrame:
    source_path = Path(input_path)
    source = (
        pd.read_parquet(source_path)
        if source_path.suffix == ".parquet"
        else pd.read_csv(source_path)
    )
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("shard_index must satisfy 0 <= shard_index < shard_count")
    if "module_id" not in source:
        raise ValueError("Designed inventory must contain module_id")
    source = source.sort_values("module_id", kind="mergesort").reset_index(drop=True)
    source = source.iloc[shard_index::shard_count].copy()
    accepted: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    candidate_frames: list[pd.DataFrame] = []
    unit_frames: list[pd.DataFrame] = []
    position_frames: list[pd.DataFrame] = []
    for source_row in source.to_dict(orient="records"):
        try:
            row, candidates, units, positions = infer_designed_boundary(
                source_row,
                dssp_executable=dssp_executable,
                foldseek_binary=foldseek_binary,
                mafft_binary=mafft_binary,
            )
            candidate_frames.append(candidates)
            if bool(row.get("strict_dual_evidence_passed")):
                accepted.append(row)
                unit_frames.append(units)
                position_frames.append(positions)
            else:
                exclusions.append(row)
        except Exception as exc:
            exclusions.append(
                {
                    **source_row,
                    "corpus_version": DESIGNED_CORPUS_VERSION,
                    "boundary_method_version": DESIGNED_BOUNDARY_METHOD,
                    "boundary_refinement_status": "designed_boundary_error",
                    "strict_dual_evidence_passed": False,
                    "boundary_error": f"{type(exc).__name__}: {exc}",
                }
            )
    accepted_frame = pd.DataFrame(accepted)
    accepted_mappings = accepted_frame.copy()
    if not accepted_frame.empty:
        accepted_frame = accepted_frame.sort_values(
            ["unit_sequence", "evidence_tier", "module_id"], kind="mergesort"
        ).drop_duplicates("unit_sequence", keep="first")
    outputs = (
        (Path(output_path), accepted_frame),
        (Path(candidates_path), pd.concat(candidate_frames, ignore_index=True) if candidate_frames else pd.DataFrame()),
        (Path(units_path), pd.concat(unit_frames, ignore_index=True) if unit_frames else pd.DataFrame()),
        (Path(positions_path), pd.concat(position_frames, ignore_index=True) if position_frames else pd.DataFrame()),
        (Path(exclusions_path), pd.DataFrame(exclusions)),
    )
    for path, frame in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".parquet":
            frame.to_parquet(path, index=False)
            frame.to_csv(path.with_suffix(".csv"), index=False)
        else:
            frame.to_csv(path, index=False)
    mapping_path = Path(output_path).with_name(
        Path(output_path).stem + "_source_mappings.parquet"
    )
    mapping_path.parent.mkdir(parents=True, exist_ok=True)
    accepted_mappings.to_parquet(mapping_path, index=False)
    accepted_mappings.to_csv(mapping_path.with_suffix(".csv"), index=False)
    manifest = {
        "corpus_version": DESIGNED_CORPUS_VERSION,
        "boundary_method": DESIGNED_BOUNDARY_METHOD,
        "shard_index": shard_index,
        "shard_count": shard_count,
        "input_rows": len(source),
        "accepted_source_rows": len(accepted_mappings),
        "accepted_unique_middle_units": len(accepted_frame),
        "exclusion_rows": len(exclusions),
        "foldseek_binary": str(Path(foldseek_binary).resolve()),
        "mafft_binary": str(mafft_binary),
        "dssp_executable": str(Path(dssp_executable).resolve()),
    }
    Path(output_path).with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return accepted_frame.reset_index(drop=True)


def finalize_designed_catalog(
    mapping_paths: Iterable[str | Path],
    exclusion_paths: Iterable[str | Path],
    output_path: str | Path,
    *,
    candidate_paths: Iterable[str | Path] = (),
    unit_paths: Iterable[str | Path] = (),
    position_paths: Iterable[str | Path] = (),
) -> pd.DataFrame:
    """Deduplicate strict-pass designed modules across resumable shards."""
    mappings_files = [Path(path) for path in mapping_paths]
    if not mappings_files:
        raise ValueError("At least one designed source-mapping shard is required")
    mappings = pd.concat(
        [pd.read_parquet(path) for path in mappings_files], ignore_index=True
    )
    catalog = mappings.sort_values(
        ["unit_sequence", "evidence_tier", "module_id"], kind="mergesort"
    ).drop_duplicates("unit_sequence", keep="first")
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_parquet(destination, index=False)
    catalog.to_csv(destination.with_suffix(".csv"), index=False)
    mapping_destination = destination.with_name(
        destination.stem + "_source_mappings.parquet"
    )
    mappings.to_parquet(mapping_destination, index=False)
    mappings.to_csv(mapping_destination.with_suffix(".csv"), index=False)
    exclusion_files = [Path(path) for path in exclusion_paths]
    exclusions = (
        pd.concat(
            [
                pd.read_parquet(path)
                if path.suffix == ".parquet"
                else pd.read_csv(path)
                for path in exclusion_files
            ],
            ignore_index=True,
        )
        if exclusion_files
        else pd.DataFrame()
    )
    exclusion_destination = destination.with_name(
        destination.stem + "_exclusions.csv"
    )
    exclusions.to_csv(exclusion_destination, index=False)
    audit_tables = {
        "boundary_candidates": [Path(path) for path in candidate_paths],
        "unit_alignment": [Path(path) for path in unit_paths],
        "position_variability": [Path(path) for path in position_paths],
    }
    audit_row_counts: dict[str, int] = {}
    for suffix, files in audit_tables.items():
        frames = [pd.read_parquet(path) for path in files]
        frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        audit_row_counts[suffix] = len(frame)
        audit_destination = destination.with_name(
            destination.stem + f"_{suffix}.parquet"
        )
        frame.to_parquet(audit_destination, index=False)
        frame.to_csv(audit_destination.with_suffix(".csv"), index=False)
    destination.with_suffix(".manifest.json").write_text(
        json.dumps(
            {
                "corpus_version": DESIGNED_CORPUS_VERSION,
                "boundary_method": DESIGNED_BOUNDARY_METHOD,
                "mapping_shards": [str(path.resolve()) for path in mappings_files],
                "exclusion_shards": [str(path.resolve()) for path in exclusion_files],
                "candidate_shards": [
                    str(path.resolve()) for path in audit_tables["boundary_candidates"]
                ],
                "unit_shards": [
                    str(path.resolve()) for path in audit_tables["unit_alignment"]
                ],
                "position_shards": [
                    str(path.resolve()) for path in audit_tables["position_variability"]
                ],
                "accepted_source_rows": len(mappings),
                "accepted_unique_middle_units": len(catalog),
                "exclusion_rows": len(exclusions),
                "audit_row_counts": audit_row_counts,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return catalog.reset_index(drop=True)
