"""Natural and designed repeat-module curation."""

from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import date
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests
from Bio import Align
from Bio.PDB.Polypeptide import protein_letters_3to1_extended

from .constants import AMINO_ACIDS, validate_protein_sequence
from .io import sha256_file, utc_now, write_json_atomic
from .periodicity import (
    BOUNDARY_METHOD_VERSION,
    MODULE_SELECTION_POLICY,
    RepeatCandidate,
    infer_repeat_boundaries,
    materialize_repeat_boundary,
    select_primitive_candidate,
)
from .secondary_structure import (
    SECONDARY_STRUCTURE_METHOD_VERSION,
    SecondaryStructureSupport,
    select_joint_sequence_structure_candidate,
)
from .schemas import ModuleRecord

REPEATSDB_API = "https://repeatsdb.org/api/production"


def _parse_fasta(text: str) -> tuple[str, int]:
    """Return an annotation FASTA sequence and its absolute coordinate origin."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    header = next((line for line in lines if line.startswith(">")), "")
    match = re.search(r"\brange\s+(-?\d+):(-?\d+)", header)
    origin = int(match.group(1)) if match else 1
    sequence = "".join(line for line in lines if not line.startswith(">"))
    return validate_protein_sequence(sequence), origin


def _middle_unit(loci: list[dict[str, object]], region_index: int) -> dict[str, object] | None:
    units = [
        locus
        for locus in loci
        if locus.get("type") == "unit"
        and locus.get("parent") is not None
        and int(locus["parent"]) == region_index
    ]
    if not units:
        return None
    units.sort(key=lambda unit: int(unit["start"]))
    region = loci[region_index]
    region_midpoint = (int(region["start"]) + int(region["end"])) / 2
    return min(
        units,
        key=lambda unit: (
            abs((int(unit["start"]) + int(unit["end"])) / 2 - region_midpoint),
            int(unit["start"]),
        ),
    )


def _recorded_unit_copies_from_row(
    row: dict[str, object],
    *,
    default_sequence: str,
    default_start: int,
    default_end: int,
) -> tuple[list[str], list[tuple[int, int]]]:
    try:
        sequences = [
            validate_protein_sequence(str(value))
            for value in json.loads(str(row.get("unit_sequences_json") or "[]"))
        ]
    except (TypeError, ValueError, json.JSONDecodeError):
        sequences = []
    if not sequences:
        return [default_sequence], [(default_start, default_end)]
    try:
        coordinates = [
            (int(value[0]), int(value[1]))
            for value in json.loads(str(row.get("source_unit_coordinates_json") or "[]"))
        ]
    except (TypeError, ValueError, json.JSONDecodeError, IndexError):
        coordinates = []
    if len(coordinates) != len(sequences):
        region_start = int(row.get("repeat_region_start") or default_start)
        period = int(row.get("period") or row.get("primitive_period") or len(sequences[0]))
        coordinates = [
            (region_start + index * period, region_start + (index + 1) * period - 1)
            for index in range(len(sequences))
        ]
    return sequences, coordinates


def _selected_middle_unit_from_row(
    row: dict[str, object],
    *,
    default_sequence: str,
    default_start: int,
    default_end: int,
) -> tuple[str, int, int, int, int]:
    """Return the repeat-region middle unit, preferring recorded real copies.

    The returned indices are one-based. Natural units may have unequal lengths,
    so their recorded source coordinates take precedence over an equal-period
    reconstruction.
    """
    sequences, coordinates = _recorded_unit_copies_from_row(
        row,
        default_sequence=default_sequence,
        default_start=default_start,
        default_end=default_end,
    )
    region_start = int(row.get("repeat_region_start") or coordinates[0][0])
    region_end = int(row.get("repeat_region_end") or coordinates[-1][1])
    region_midpoint = (region_start + region_end) / 2
    selected_index = min(
        range(len(sequences)),
        key=lambda index: (
            abs((coordinates[index][0] + coordinates[index][1]) / 2 - region_midpoint),
            coordinates[index][0],
        ),
    )
    start, end = coordinates[selected_index]
    return sequences[selected_index], start, end, selected_index + 1, len(sequences)


def fetch_natural_modules(
    output_path: str | Path,
    *,
    per_class: int = 20,
    timeout: int = 30,
    page_size: int = 100,
    boundary_workers: int = 1,
) -> pd.DataFrame:
    """Select 100 structurally annotated units, balanced over five classes."""
    if not 1 <= page_size <= 100:
        raise ValueError("RepeatsDB v4 page_size must be between 1 and 100")
    session = requests.Session()
    candidates: dict[int, list[dict[str, object]]] = defaultdict(list)
    # Query each class independently.  This is both much faster than walking
    # the full AlphaFold-heavy collection and makes empty classes explicit.
    for requested_class in range(1, 6):
        skip = 0
        total = None
        candidate_target = max(per_class * 10, per_class + 25)
        while (total is None or skip < total) and len(candidates[requested_class]) < candidate_target:
            response = session.get(
                f"{REPEATSDB_API}/annotations",
                params={
                    "limit": page_size,
                    "skip": skip,
                    "chain.source": "PDB",
                    "region.classes": str(requested_class),
                },
                timeout=timeout,
            )
            response.raise_for_status()
            payload = response.json()
            total = int(payload["count"])
            items = payload.get("items", {})
            iterable = items.values() if isinstance(items, dict) else items
            for item in iterable:
                content = item.get("content", {})
                chain = content.get("chain", {})
                if "PDB" not in str(chain.get("source", "")).upper():
                    continue
                loci = content.get("loci", [])
                for region_index, region in enumerate(loci):
                    if region.get("type") != "region" or not region.get("class"):
                        continue
                    class_id = int(str(region["class"]).split(".")[0])
                    if class_id != requested_class:
                        continue
                    unit = _middle_unit(loci, region_index)
                    if unit is None:
                        continue
                    classification = str(region["class"])
                    candidates[class_id].append(
                        {
                            "uuid": item["uuid"],
                            "pdb": str(chain.get("structure", "")).lower(),
                            "chain": str(chain.get("id", "")),
                            "class_id": class_id,
                            "classification": classification,
                            "topology": ".".join(classification.split(".")[:2]),
                            "unit_start": int(unit["start"]),
                            "unit_end": int(unit["end"]),
                            "reviewed": bool(content.get("reviewed", False)),
                        }
                    )
            skip += page_size

    records: list[ModuleRecord] = []
    seen_sequences: set[str] = set()
    seen_regions: set[tuple[str, str, str]] = set()
    fetched: dict[str, tuple[str, int, str]] = {}

    def materialize(candidate: dict[str, object]) -> ModuleRecord | None:
        uuid = str(candidate["uuid"])
        if uuid not in fetched:
            response = session.get(
                f"{REPEATSDB_API}/annotations/{uuid}", params={"format": "fasta"}, timeout=timeout
            )
            response.raise_for_status()
            fasta_sha256 = hashlib.sha256(response.content).hexdigest()
            try:
                chain_sequence, origin = _parse_fasta(response.text)
            except ValueError:
                # Ambiguous residues are explicitly ineligible for the
                # standard-20-AA catalog, but should not abort the run.
                fetched[uuid] = ("", 1, fasta_sha256)
                return None
            fetched[uuid] = (chain_sequence, origin, fasta_sha256)
        chain_sequence, origin, fasta_sha256 = fetched[uuid]
        if not chain_sequence:
            return None
        start = int(candidate["unit_start"])
        end = int(candidate["unit_end"])
        local_start = start - origin
        local_end = end - origin + 1
        if local_start < 0 or local_end > len(chain_sequence) or end < start:
            return None
        unit_sequence = chain_sequence[local_start:local_end]
        if not unit_sequence or set(unit_sequence) - set(AMINO_ACIDS):
            return None
        identity = (str(candidate["pdb"]), str(candidate["chain"]), str(candidate["classification"]))
        if unit_sequence in seen_sequences or identity in seen_regions:
            return None
        seen_sequences.add(unit_sequence)
        seen_regions.add(identity)
        return ModuleRecord(
            module_id=f"natural_{candidate['pdb']}_{candidate['chain']}_{start}_{end}",
            collection="natural100",
            family=f"RepeatsDB class {candidate['class_id']}",
            unit_sequence=unit_sequence,
            evidence_tier="A",
            source_name="RepeatsDB v4 / RCSB PDB",
            source_url=f"https://repeatsdb.org/structure/{candidate['pdb']}",
            source_accession=str(candidate["pdb"]),
            source_chain=str(candidate["chain"]),
            unit_start=start,
            unit_end=end,
            reviewed=bool(candidate["reviewed"]),
            boundary_method="RepeatsDB unit coordinates; internal unit preferred",
            retrieved_date=date.today().isoformat(),
            source_sha256=fasta_sha256,
            license_name="RepeatsDB terms; underlying PDB structure",
            citation="RepeatsDB v4 (Nucleic Acids Research, 2025)",
            full_sequence=chain_sequence,
            full_sequence_origin=origin,
            full_sequence_sha256=hashlib.sha256(chain_sequence.encode()).hexdigest(),
            source_annotation_id=uuid,
            notes=(
                f"classification={candidate['classification']}; topology={candidate['topology']}; "
                f"fasta_origin={origin}; retrieved={date.today().isoformat()}"
            ),
        )

    ranked_by_class: dict[int, list[dict[str, object]]] = {
        class_id: sorted(
            candidates[class_id],
            key=lambda row: (not row["reviewed"], row["topology"], row["pdb"], row["chain"], row["unit_start"]),
        )
        for class_id in range(1, 6)
    }
    for class_id in range(1, 6):
        for candidate in ranked_by_class[class_id]:
            if sum(record.family == f"RepeatsDB class {class_id}" for record in records) >= per_class:
                break
            record = materialize(candidate)
            if record is not None:
                records.append(record)

    # Some RepeatsDB releases have no entries in every top-level class (v4
    # currently returns zero PDB annotations for class 1).  Fill the deficit
    # from the least represented topology, preserving only PDB-backed units.
    selected_ids = {record.module_id for record in records}
    topology_counts: dict[str, int] = defaultdict(int)
    for record in records:
        topology = re.search(r"topology=([^;]+)", record.notes)
        topology_counts[topology.group(1) if topology else record.family] += 1
    remaining = [
        candidate
        for rows in ranked_by_class.values()
        for candidate in rows
        if (
            f"natural_{candidate['pdb']}_{candidate['chain']}_"
            f"{candidate['unit_start']}_{candidate['unit_end']}"
        )
        not in selected_ids
    ]
    while len(records) < per_class * 5:
        remaining.sort(
            key=lambda row: (
                topology_counts[str(row["topology"])],
                not row["reviewed"],
                row["class_id"],
                row["topology"],
                row["pdb"],
                row["chain"],
                row["unit_start"],
            )
        )
        added = False
        for index, candidate in enumerate(remaining):
            module_id = (
                f"natural_{candidate['pdb']}_{candidate['chain']}_"
                f"{candidate['unit_start']}_{candidate['unit_end']}"
            )
            if module_id in selected_ids:
                continue
            record = materialize(candidate)
            remaining.pop(index)
            if record is None:
                break
            records.append(record)
            selected_ids.add(record.module_id)
            topology_counts[str(candidate["topology"])] += 1
            added = True
            break
        if not added and not remaining:
            break

    if len(records) != per_class * 5:
        counts = {class_id: sum(record.family == f"RepeatsDB class {class_id}" for record in records) for class_id in range(1, 6)}
        raise RuntimeError(f"Could not collect {per_class * 5} unique PDB-backed natural units: {counts}")
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    frame = refine_module_boundaries(
        pd.DataFrame(record.to_dict() for record in records),
        audit_path=destination.with_name(destination.stem + "_boundary_audit.parquet"),
        candidates_path=destination.with_name(destination.stem + "_period_candidates.parquet"),
        unit_alignment_path=destination.with_name(destination.stem + "_unit_alignment.parquet"),
        position_variability_path=destination.with_name(
            destination.stem + "_position_variability.parquet"
        ),
        workers=boundary_workers,
    )
    frame.to_parquet(destination, index=False)
    frame.to_csv(destination.with_suffix(".csv"), index=False)
    return frame


def longest_adjacent_exact_repeat(sequence: str, minimum: int = 6, maximum: int = 120) -> tuple[str, int] | None:
    """Find the longest adjacent exact repeat, used only for auditable designed constructs."""
    normalized = validate_protein_sequence(sequence)
    upper = min(maximum, len(normalized) // 2)
    for length in range(upper, minimum - 1, -1):
        for start in range(0, len(normalized) - 2 * length + 1):
            unit = normalized[start : start + length]
            if unit == normalized[start + length : start + 2 * length]:
                return unit, start
    return None


def parse_fasta_modules(
    fasta_paths: Iterable[str | Path],
    *,
    family: str,
    evidence_tier: str,
    source_url: str,
    exclusions_path: str | Path | None = None,
) -> pd.DataFrame:
    """Extract exact adjacent repeat units from public designed FASTA files."""
    rows: list[dict[str, object]] = []
    exclusions: list[dict[str, str]] = []
    seen: set[str] = set()
    for path_like in fasta_paths:
        path = Path(path_like)
        name = None
        sequence_parts: list[str] = []
        for line in path.read_text().splitlines() + [">"]:
            if line.startswith(">"):
                if name is not None:
                    raw_sequence = "".join(sequence_parts).rstrip("*")
                    try:
                        sequence = validate_protein_sequence(raw_sequence)
                    except ValueError as exc:
                        exclusions.append({"source_record": name, "reason": str(exc)})
                        name = line[1:].strip() or None
                        sequence_parts = []
                        continue
                    repeat = longest_adjacent_exact_repeat(sequence)
                    if repeat is not None and repeat[0] not in seen:
                        unit, start = repeat
                        seen.add(unit)
                        identifier = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")
                        rows.append(
                            ModuleRecord(
                                module_id=f"designed_{identifier}",
                                collection="designed_all",
                                family=family,
                                unit_sequence=unit,
                                evidence_tier=evidence_tier,
                                source_name=path.name,
                                source_url=source_url,
                                source_accession=name,
                                unit_start=start + 1,
                                unit_end=start + len(unit),
                                boundary_method="longest pair of adjacent exact internal repeats",
                                retrieved_date=date.today().isoformat(),
                                source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                license_name="See source repository or supplement",
                                citation="See source_url and source_accession",
                                full_sequence=sequence,
                                full_sequence_origin=1,
                                full_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                                notes=f"exact_adjacent_repeat; file_sha256={hashlib.sha256(path.read_bytes()).hexdigest()}",
                            ).to_dict()
                        )
                    elif repeat is None:
                        exclusions.append({"source_record": name, "reason": "no adjacent exact repeat of at least 6AA"})
                name = line[1:].strip() or None
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())
    if exclusions_path is not None:
        destination = Path(exclusions_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(exclusions).to_csv(destination, index=False)
    return pd.DataFrame(rows)


def parse_dhr_supplement(path: str | Path) -> pd.DataFrame:
    """Extract DHR1--83 units using author-provided repeat lengths.

    Supplementary Table 2 defines the repeat length and Table 16 defines the
    complete experimentally tested construct.  A record is accepted only if
    two adjacent exact units of that author-provided length are present.
    """
    source = Path(path)
    text = source.read_text(errors="replace")
    table2_starts = [match.start() for match in re.finditer(r"Supplementary Table 2", text)]
    table3_starts = [match.start() for match in re.finditer(r"Supplementary Table 3", text)]
    table16_starts = [match.start() for match in re.finditer(r"Supplementary Table 16", text)]
    reference_starts = [match.start() for match in re.finditer(r"Supplementary References", text)]
    if not all((table2_starts, table3_starts, table16_starts, reference_starts)):
        raise ValueError("Could not locate DHR Supplementary Tables 2 and 16")
    length_section = text[table2_starts[-1] : table3_starts[-1]]
    repeat_lengths = {
        int(number): int(length)
        for number, length in re.findall(r"^\s*(\d{1,2})\s+(\d{2,3})\s+", length_section, re.MULTILINE)
        if 1 <= int(number) <= 83
    }
    sequence_section = text[table16_starts[-1] : reference_starts[-1]]
    sequences: dict[int, str] = {}
    current: int | None = None
    for raw_line in sequence_section.splitlines():
        line = raw_line.replace("\f", "")
        start = re.match(r"^\s*DHR(\d+)\s+([A-Z]+)\s*$", line)
        if start:
            current = int(start.group(1))
            sequences[current] = start.group(2)
            continue
        continuation = re.match(r"^\s+([A-Z]+)\s*$", line)
        if continuation and current is not None:
            sequences[current] += continuation.group(1)
    if set(repeat_lengths) != set(range(1, 84)) or set(sequences) != set(range(1, 84)):
        raise ValueError(
            f"Expected DHR1--83; lengths={len(repeat_lengths)}, sequences={len(sequences)}"
        )

    structure_validated = {4, 5, 7, 8, 10, 14, 18, 49, 53, 54, 64, 71, 76, 79, 81}
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    rows: list[dict[str, object]] = []
    for number in range(1, 84):
        sequence = validate_protein_sequence(sequences[number])
        repeat_length = repeat_lengths[number]
        matches = [
            position
            for position in range(len(sequence) - 2 * repeat_length + 1)
            if sequence[position : position + repeat_length]
            == sequence[position + repeat_length : position + 2 * repeat_length]
        ]
        if not matches:
            raise ValueError(f"DHR{number} has no exact adjacent units of author length {repeat_length}")
        position = matches[0]
        unit = sequence[position : position + repeat_length]
        rows.append(
            ModuleRecord(
                module_id=f"designed_DHR{number}",
                collection="designed_all",
                family="DHR",
                unit_sequence=unit,
                evidence_tier="A" if number in structure_validated else "B",
                source_name=source.name,
                source_url="https://www.nature.com/articles/nature16162",
                source_accession=f"DHR{number}",
                unit_start=position + 1,
                unit_end=position + repeat_length,
                reviewed=True,
                boundary_method="author repeat length (Supplementary Table 2) plus adjacent exact-unit verification",
                retrieved_date=date.today().isoformat(),
                source_sha256=source_hash,
                license_name="Nature author manuscript supplementary material",
                citation="Brunette et al., Nature 2015, doi:10.1038/nature16162",
                full_sequence=sequence,
                full_sequence_origin=1,
                full_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                notes=(
                    f"full_construct_length={len(sequence)}; repeat_length={repeat_length}; "
                    "all 83 constructs experimentally tested; A marks deposited crystal structures"
                ),
            ).to_dict()
        )
    return pd.DataFrame(rows)


def _pdb_chain_sequences(path: Path) -> dict[str, str]:
    chains: dict[str, list[str]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    model_seen = False
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                if model_seen:
                    break
                model_seen = True
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM") or line[12:16].strip() != "CA" or line[16] not in {" ", "A"}:
                continue
            chain = line[21].strip() or "_"
            residue_key = (chain, line[22:27])
            if residue_key in seen:
                continue
            seen.add(residue_key)
            amino_acid = protein_letters_3to1_extended.get(line[17:20].strip().upper())
            if amino_acid:
                chains[chain].append(amino_acid)
    return {chain: "".join(sequence) for chain, sequence in chains.items() if sequence}


def parse_pdb_exact_modules(
    pdb_paths: Iterable[str | Path],
    *,
    family: str,
    source_url: str,
    experimental_accessions: Iterable[str] = (),
    structure_accessions: Iterable[str] = (),
) -> pd.DataFrame:
    """Extract exact adjacent repeat units from author-supplied PDB models."""
    experimental = {str(value).upper() for value in experimental_accessions}
    structures = {str(value).upper() for value in structure_accessions}
    rows: list[dict[str, object]] = []
    for path_like in pdb_paths:
        path = Path(path_like)
        accession = re.sub(r"_(design|model|xtal).*$", "", path.stem, flags=re.IGNORECASE)
        chains = _pdb_chain_sequences(path)
        if not chains:
            continue
        chain, sequence = max(chains.items(), key=lambda item: (len(item[1]), item[0]))
        repeat = longest_adjacent_exact_repeat(sequence)
        if repeat is None:
            continue
        unit, position = repeat
        upper_accession = accession.upper()
        evidence = "A" if upper_accession in structures else "B" if upper_accession in experimental else "C"
        rows.append(
            ModuleRecord(
                module_id=f"designed_{accession}",
                collection="designed_all",
                family=family,
                unit_sequence=unit,
                evidence_tier=evidence,
                source_name=path.name,
                source_url=source_url,
                source_accession=accession,
                source_chain=chain,
                unit_start=position + 1,
                unit_end=position + len(unit),
                reviewed=True,
                boundary_method="two adjacent exact internal repeats in author-supplied PDB model",
                retrieved_date=date.today().isoformat(),
                source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                license_name="CC BY 4.0 article; supplementary-data terms apply",
                citation="Huddy et al., Nature 2024, doi:10.1038/s41586-024-07188-4",
                full_sequence=sequence,
                full_sequence_origin=1,
                full_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                notes=f"model_chain_length={len(sequence)}; exact_repeat_start={position + 1}",
            ).to_dict()
        )
    return pd.DataFrame(rows)


def parse_module_manifest(path: str | Path) -> pd.DataFrame:
    """Read author-supplied or manually reviewed exact repeat-unit boundaries."""
    source = Path(path)
    frame = pd.read_csv(source)
    required = {"module_id", "unit_sequence", "family", "evidence_tier", "source_url"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Designed manifest is missing columns: {sorted(missing)}")
    frame["unit_sequence"] = frame["unit_sequence"].map(validate_protein_sequence)
    if not frame["evidence_tier"].isin(["A", "B", "C"]).all():
        raise ValueError("Designed evidence_tier must be A, B, or C")
    frame["collection"] = "designed_all"
    frame["source_name"] = frame.get("source_name", source.name)
    frame["source_accession"] = frame.get("source_accession", frame["module_id"])
    frame["boundary_method"] = frame.get("boundary_method", "author-supplied or manually reviewed boundary")
    frame["retrieved_date"] = frame.get("retrieved_date", date.today().isoformat())
    frame["download_date"] = frame.get("download_date", frame["retrieved_date"])
    frame["source_sha256"] = frame.get("source_sha256", hashlib.sha256(source.read_bytes()).hexdigest())
    frame["license_name"] = frame.get("license_name", "See primary source")
    frame["citation"] = frame.get("citation", "See source_url")
    return frame


def _infer_boundary_worker(
    payload: tuple[dict[str, object], float, int, int]
) -> tuple[dict[str, object], object | None, list[object], str | None]:
    """Pickle-safe worker for one full-protein boundary inference."""
    source_row, fixed_threshold, minimum_period, maximum_period = payload
    row = dict(source_row)
    if (
        row.get("boundary_refinement_status") == "refined"
        and row.get("boundary_method_version") == BOUNDARY_METHOD_VERSION
    ):
        return row, None, [], "already_refined"
    prior_sequence = validate_protein_sequence(str(row["unit_sequence"]))
    row["prior_unit_sequence"] = prior_sequence
    row["prior_unit_length"] = len(prior_sequence)
    row["prior_unit_start"] = row.get("unit_start")
    row["prior_unit_end"] = row.get("unit_end")
    full_value = row.get("full_sequence")
    full_sequence = "" if full_value is None or pd.isna(full_value) else str(full_value).strip()
    if not full_sequence:
        row["boundary_refinement_status"] = "missing_full_sequence"
        row["boundary_method_version"] = "source-boundary-unrefined"
        return row, None, [], "missing_full_sequence"
    origin_value = row.get("full_sequence_origin")
    origin = 1 if origin_value is None or pd.isna(origin_value) else int(origin_value)
    start_value = row.get("unit_start")
    prior_start = None if start_value is None or pd.isna(start_value) else int(start_value)
    inference_origin = origin
    inference_prior_start = prior_start
    mapping_value = row.get("full_sequence_auth_mapping_json")
    if isinstance(mapping_value, str) and mapping_value.startswith("[") and prior_start is not None:
        mapping = [_auth_coordinate(value) for value in json.loads(mapping_value)]
        mapped_indices = [
            index
            for index, coordinate in enumerate(mapping)
            if coordinate is not None
            and prior_start <= coordinate <= int(row.get("unit_end") or prior_start)
        ]
        if mapped_indices:
            inference_origin = 1
            inference_prior_start = min(mapped_indices) + 1
            row["periodicity_scan_coordinate_system"] = "rcsb_entity_poly_seq_index"
    try:
        result, candidates = infer_repeat_boundaries(
            full_sequence,
            full_sequence_origin=inference_origin,
            prior_unit_sequence=prior_sequence,
            prior_unit_start=inference_prior_start,
            minimum_period=minimum_period,
            maximum_period=maximum_period,
            fixed_threshold=fixed_threshold,
        )
    except (ValueError, IndexError) as exc:
        row["boundary_refinement_status"] = "no_supported_period"
        row["boundary_refinement_error"] = str(exc)
        row["boundary_method_version"] = BOUNDARY_METHOD_VERSION
        return row, None, [], "no_supported_period"
    return row, result, candidates, None


def refine_module_boundaries(
    frame: pd.DataFrame,
    *,
    audit_path: str | Path | None = None,
    candidates_path: str | Path | None = None,
    unit_alignment_path: str | Path | None = None,
    position_variability_path: str | Path | None = None,
    fixed_threshold: float = 0.8,
    minimum_period: int = 3,
    maximum_period: int = 120,
    workers: int = 1,
) -> pd.DataFrame:
    """Replace source candidates with the middle inferred primitive unit.

    The source boundary remains in ``prior_*`` columns.  Records without a
    complete protein sequence are retained but explicitly marked unrefined;
    they are never presented as frequency-derived boundaries.
    """
    refined_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    unit_rows: list[dict[str, object]] = []
    position_rows: list[dict[str, object]] = []
    if workers < 1:
        raise ValueError("workers must be positive")
    payloads = [
        (row, fixed_threshold, minimum_period, maximum_period)
        for row in frame.to_dict(orient="records")
    ]
    if workers == 1:
        outcomes = map(_infer_boundary_worker, payloads)
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        outcomes = executor.map(_infer_boundary_worker, payloads)
    try:
        for row, result, candidates, status in outcomes:
            if status is not None:
                refined_rows.append(row)
                continue
            assert result is not None
            result_payload = result.to_dict()
            row.update(result_payload)
            row["unit_sequence"] = result.selected_module_sequence
            row["unit_start"] = result.selected_module_start
            row["unit_end"] = result.selected_module_end
            row["unit_length"] = result.period
            row["primitive_period"] = result.period
            row["periodicity_score"] = result.score
            row["periodicity_confidence"] = result.confidence
            row["boundary_method"] = (
                f"{BOUNDARY_METHOD_VERSION}; source boundary retained in prior_* columns"
            )
            source_fallback = result.selection_reason.startswith(
                "no independently supported sequence subperiod"
            )
            row["boundary_refinement_status"] = (
                "source_prior_fallback" if source_fallback else "refined"
            )
            row["fixed_variable_assessment_status"] = (
                "insufficient_inferred_copies" if source_fallback else "inferred_from_repeat_copies"
            )
            refined_rows.append(row)
            for repeat_index, unit_sequence in enumerate(result.unit_sequences, start=1):
                unit_start = result.repeat_region_start + (repeat_index - 1) * result.period
                unit_rows.append(
                    {
                        "module_id": row["module_id"],
                        "repeat_index": repeat_index,
                        "unit_start": unit_start,
                        "unit_end": unit_start + result.period - 1,
                        "unit_sequence": unit_sequence,
                        "is_first_module": repeat_index == 1,
                        "is_selected_module": repeat_index == result.selected_module_index,
                        "module_selection_policy": result.selected_module_policy,
                    }
                )
            for position, (consensus_amino_acid, conservation) in enumerate(
                zip(result.consensus_sequence, result.position_conservation, strict=True), start=1
            ):
                variants = sorted({unit[position - 1] for unit in result.unit_sequences})
                position_rows.append(
                    {
                        "module_id": row["module_id"],
                        "module_position": position,
                        "consensus_amino_acid": consensus_amino_acid,
                        "conservation": conservation,
                        "fixed": conservation >= fixed_threshold,
                        "variants_json": json.dumps(variants),
                    }
                )
            for rank, candidate in enumerate(candidates, start=1):
                candidate_rows.append(
                    {
                        "module_id": row["module_id"],
                        "candidate_rank": rank,
                        **candidate.to_dict(),
                    }
                )
    finally:
        if workers != 1:
            executor.shutdown()
    refined = pd.DataFrame(refined_rows)
    if audit_path is not None:
        destination = Path(audit_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        refined.to_parquet(destination, index=False)
        refined.to_csv(destination.with_suffix(".csv"), index=False)
    if candidates_path is not None:
        destination = Path(candidates_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        candidates_frame = pd.DataFrame(candidate_rows)
        candidates_frame.to_parquet(destination, index=False)
        candidates_frame.to_csv(destination.with_suffix(".csv"), index=False)
    if unit_alignment_path is not None:
        destination = Path(unit_alignment_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        units_frame = pd.DataFrame(unit_rows)
        units_frame.to_parquet(destination, index=False)
        units_frame.to_csv(destination.with_suffix(".csv"), index=False)
    if position_variability_path is not None:
        destination = Path(position_variability_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        positions_frame = pd.DataFrame(position_rows)
        positions_frame.to_parquet(destination, index=False)
        positions_frame.to_csv(destination.with_suffix(".csv"), index=False)
    return refined


def reselect_module_boundaries(
    frame: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    secondary_structure_support: pd.DataFrame | None = None,
    output_path: str | Path | None = None,
    unit_alignment_path: str | Path | None = None,
    position_variability_path: str | Path | None = None,
    fixed_threshold: float = 0.8,
) -> pd.DataFrame:
    """Apply the current selector to a previously computed candidate audit.

    This makes selector revisions cheap and transparent: Fourier and
    self-similarity measurements remain frozen, while the versioned selection
    rule can be compared without rescanning complete sequences.
    """
    candidate_fields = tuple(RepeatCandidate.__dataclass_fields__)
    support_fields = tuple(SecondaryStructureSupport.__dataclass_fields__)
    grouped = {key: value for key, value in candidates.groupby("module_id", sort=False)}
    support_grouped = (
        {
            key: value
            for key, value in secondary_structure_support.groupby("module_id", sort=False)
        }
        if secondary_structure_support is not None
        else {}
    )
    refined_rows: list[dict[str, object]] = []
    unit_rows: list[dict[str, object]] = []
    position_rows: list[dict[str, object]] = []
    for source_row in frame.to_dict(orient="records"):
        row = dict(source_row)
        module_id = str(row["module_id"])
        prior_sequence = validate_protein_sequence(str(row["prior_unit_sequence"]))
        prior_period = len(prior_sequence)
        source_candidates = (
            [
                RepeatCandidate(**{field: record[field] for field in candidate_fields})
                for record in grouped[module_id]
                .sort_values("candidate_rank")
                .to_dict(orient="records")
            ]
            if module_id in grouped
            else []
        )
        selected_support: SecondaryStructureSupport | None = None
        if not source_candidates:
            selected = None
            reason = (
                "no amino-acid periodicity candidate passed the scan floor; "
                "source-annotated unit retained"
            )
        elif secondary_structure_support is not None:
            source_support_frame = support_grouped.get(module_id)
            support_by_key: dict[tuple[int, int, int], SecondaryStructureSupport] = {}
            if source_support_frame is not None:
                for record in source_support_frame.to_dict(orient="records"):
                    support = SecondaryStructureSupport(
                        **{field: record[field] for field in support_fields}
                    )
                    support_by_key[(support.period, support.local_start, support.local_end)] = support
            paired_candidates: list[RepeatCandidate] = []
            paired_supports: list[SecondaryStructureSupport] = []
            for candidate in source_candidates:
                key = (candidate.period, candidate.local_start, candidate.local_end)
                if key in support_by_key:
                    paired_candidates.append(candidate)
                    paired_supports.append(support_by_key[key])
            if paired_candidates:
                selected, selected_support, reason = select_joint_sequence_structure_candidate(
                    paired_candidates,
                    paired_supports,
                    prior_period=prior_period,
                )
            else:
                selected = None
                reason = (
                    "no residue-level secondary-structure scores matched the sequence candidates; "
                    "source-annotated unit retained"
                )
        else:
            selected, reason = select_primitive_candidate(
                source_candidates, prior_period=prior_period
            )
        if selected is None:
            # Do not convert an absolute PDB author coordinate to a sequence
            # offset by subtraction.  Missing residues and insertion codes
            # make that transformation invalid for natural structures.  The
            # already validated source unit is retained verbatim instead.
            prior_start = int(row["prior_unit_start"])
            prior_end = int(row["prior_unit_end"])
            (
                selected_sequence,
                selected_start,
                selected_end,
                selected_index,
                selected_count,
            ) = _selected_middle_unit_from_row(
                row,
                default_sequence=prior_sequence,
                default_start=prior_start,
                default_end=prior_end,
            )
            recorded_sequences, recorded_coordinates = _recorded_unit_copies_from_row(
                row,
                default_sequence=prior_sequence,
                default_start=prior_start,
                default_end=prior_end,
            )
            row.update(
                {
                    "unit_sequence": selected_sequence,
                    "unit_start": selected_start,
                    "unit_end": selected_end,
                    "unit_length": len(selected_sequence),
                    "period": len(selected_sequence),
                    "primitive_period": len(selected_sequence),
                    "first_module_sequence": row.get("first_module_sequence", prior_sequence),
                    "first_module_start": row.get("first_module_start", prior_start),
                    "first_module_end": row.get("first_module_end", prior_end),
                    "selected_module_sequence": selected_sequence,
                    "selected_module_start": selected_start,
                    "selected_module_end": selected_end,
                    "selected_module_index": selected_index,
                    "selected_module_count": selected_count,
                    "selected_module_policy": MODULE_SELECTION_POLICY,
                    "module_selection_policy": MODULE_SELECTION_POLICY,
                    "selection_reason": reason,
                    "harmonic_ratio": 1.0,
                    "boundary_method_version": BOUNDARY_METHOD_VERSION,
                    "boundary_method": (
                        f"{BOUNDARY_METHOD_VERSION}; exact source boundary retained because "
                        "joint sequence/secondary-structure support was insufficient"
                    ),
                    "boundary_refinement_status": "source_prior_fallback",
                    "fixed_variable_assessment_status": row.get(
                        "fixed_variable_assessment_status", "source_annotation"
                    ),
                    "boundary_reselected_from_candidates": True,
                    "secondary_structure_evidence_used": (
                        secondary_structure_support is not None
                    ),
                    "secondary_structure_method_version": (
                        SECONDARY_STRUCTURE_METHOD_VERSION
                        if secondary_structure_support is not None
                        else None
                    ),
                    "secondary_structure_selected_support": False,
                }
            )
            refined_rows.append(row)
            for repeat_index, (sequence, coordinates) in enumerate(
                zip(recorded_sequences, recorded_coordinates, strict=True), start=1
            ):
                unit_rows.append(
                    {
                        "module_id": module_id,
                        "repeat_index": repeat_index,
                        "unit_start": coordinates[0],
                        "unit_end": coordinates[1],
                        "unit_sequence": sequence,
                        "is_first_module": repeat_index == 1,
                        "is_selected_module": repeat_index == selected_index,
                        "module_selection_policy": MODULE_SELECTION_POLICY,
                        "boundary_source": "source_prior_fallback",
                    }
                )
            continue
        result = materialize_repeat_boundary(
            str(row["full_sequence"]),
            selected,
            selection_reason=reason,
            full_sequence_origin=int(row.get("full_sequence_origin") or 1),
            prior_period=prior_period,
            prior_unit_start=int(row["prior_unit_start"]),
            fixed_threshold=fixed_threshold,
        )
        row.update(result.to_dict())
        row["unit_sequence"] = result.selected_module_sequence
        row["unit_start"] = result.selected_module_start
        row["unit_end"] = result.selected_module_end
        row["unit_length"] = result.period
        row["primitive_period"] = result.period
        row["periodicity_score"] = result.score
        row["periodicity_confidence"] = result.confidence
        row["boundary_method"] = (
            f"{BOUNDARY_METHOD_VERSION}; source boundary retained in prior_* columns"
        )
        source_fallback = "source-annotated unit retained" in reason or reason.startswith(
            "no independently supported sequence subperiod"
        )
        row["boundary_refinement_status"] = (
            "source_prior_fallback" if source_fallback else "refined"
        )
        row["fixed_variable_assessment_status"] = (
            "insufficient_inferred_copies" if source_fallback else "inferred_from_repeat_copies"
        )
        row["boundary_reselected_from_candidates"] = True
        row["secondary_structure_evidence_used"] = secondary_structure_support is not None
        row["secondary_structure_method_version"] = (
            SECONDARY_STRUCTURE_METHOD_VERSION if secondary_structure_support is not None else None
        )
        row["secondary_structure_selected_support"] = selected_support is not None
        if selected_support is not None:
            for field, value in selected_support.to_dict().items():
                row[f"secondary_structure_{field}"] = value
        refined_rows.append(row)
        for repeat_index, unit_sequence in enumerate(result.unit_sequences, start=1):
            unit_start = result.repeat_region_start + (repeat_index - 1) * result.period
            unit_rows.append(
                {
                    "module_id": module_id,
                    "repeat_index": repeat_index,
                    "unit_start": unit_start,
                    "unit_end": unit_start + result.period - 1,
                    "unit_sequence": unit_sequence,
                    "is_first_module": repeat_index == 1,
                    "is_selected_module": repeat_index == result.selected_module_index,
                    "module_selection_policy": result.selected_module_policy,
                }
            )
        for position, (consensus_amino_acid, conservation) in enumerate(
            zip(result.consensus_sequence, result.position_conservation, strict=True), start=1
        ):
            variants = sorted({unit[position - 1] for unit in result.unit_sequences})
            position_rows.append(
                {
                    "module_id": module_id,
                    "module_position": position,
                    "consensus_amino_acid": consensus_amino_acid,
                    "conservation": conservation,
                    "fixed": conservation >= fixed_threshold,
                    "variants_json": json.dumps(variants),
                }
            )
    refined = pd.DataFrame(refined_rows)
    for path_like, output_frame in (
        (output_path, refined),
        (unit_alignment_path, pd.DataFrame(unit_rows)),
        (position_variability_path, pd.DataFrame(position_rows)),
    ):
        if path_like is None:
            continue
        destination = Path(path_like)
        destination.parent.mkdir(parents=True, exist_ok=True)
        output_frame.to_parquet(destination, index=False)
        output_frame.to_csv(destination.with_suffix(".csv"), index=False)
    return refined


def _auth_coordinate(value: object) -> int | None:
    match = re.match(r"^-?\d+", str(value))
    return int(match.group()) if match else None


def _fetch_natural_annotation_source(
    payload: tuple[str, str, str, int]
) -> tuple[str, dict[str, object]]:
    annotation_id, pdb_id, chain_id, timeout = payload
    response = requests.get(f"{REPEATSDB_API}/annotations/{annotation_id}", timeout=timeout)
    response.raise_for_status()
    # RepeatsDB records author chain IDs.  The RCSB instance REST endpoint,
    # however, requires label_asym_id.  Case-changing fallbacks are unsafe:
    # large complexes can contain distinct author chains ``w`` and ``W``.
    graph_query = """
    query($id: String!) {
      entry(entry_id: $id) {
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            entity_id auth_asym_ids asym_ids
          }
          entity_poly { pdbx_seq_one_letter_code_can }
        }
      }
    }
    """
    graph_response = requests.post(
        "https://data.rcsb.org/graphql",
        json={"query": graph_query, "variables": {"id": pdb_id.upper()}},
        timeout=timeout,
    )
    graph_response.raise_for_status()
    entry = graph_response.json().get("data", {}).get("entry")
    if not entry:
        raise RuntimeError(f"RCSB entry not found for {pdb_id}")
    selected_entity = None
    label_chain = None
    for entity in entry.get("polymer_entities", []):
        identifiers = entity["rcsb_polymer_entity_container_identifiers"]
        author_chains = list(identifiers.get("auth_asym_ids") or [])
        label_chains = list(identifiers.get("asym_ids") or [])
        if chain_id in author_chains:
            index = author_chains.index(chain_id)
            label_chain = label_chains[index] if index < len(label_chains) else label_chains[0]
            selected_entity = entity
            break
        if chain_id in label_chains:
            label_chain = chain_id
            selected_entity = entity
            break
    if selected_entity is None or label_chain is None:
        raise RuntimeError(
            f"RCSB auth/label chain mapping not found for {pdb_id} chain {chain_id}"
        )
    instance_response = requests.get(
        f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/"
        f"{pdb_id.upper()}/{label_chain}",
        timeout=timeout,
    )
    instance_response.raise_for_status()
    identifiers = instance_response.json()[
        "rcsb_polymer_entity_instance_container_identifiers"
    ]
    entity_id = str(identifiers["entity_id"])
    canonical_sequence = re.sub(
        r"\s+",
        "",
        str(selected_entity["entity_poly"]["pdbx_seq_one_letter_code_can"]),
    ).upper()
    auth_mapping = list(identifiers["auth_to_entity_poly_seq_mapping"])
    if len(canonical_sequence) != len(auth_mapping):
        raise ValueError(
            f"RCSB sequence/mapping length mismatch for {pdb_id}/{chain_id}: "
            f"{len(canonical_sequence)} != {len(auth_mapping)}"
        )
    return annotation_id, {
        "annotation": response.json(),
        "canonical_sequence": canonical_sequence,
        "auth_mapping": auth_mapping,
        "entity_id": entity_id,
        "auth_chain_id": chain_id,
        "label_chain_id": label_chain,
        "rcsb_instance_url": instance_response.url,
        "rcsb_entity_url": (
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/"
            f"{pdb_id.upper()}/{entity_id}"
        ),
    }


def _align_units_to_reference(units: list[str], reference_index: int) -> list[str]:
    if not 0 <= reference_index < len(units):
        raise ValueError("reference_index is outside the repeat-unit list")
    reference = units[reference_index]
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -3.0
    aligner.extend_gap_score = -0.5
    aligned_units: list[str] = []
    for index, unit in enumerate(units):
        if index == reference_index:
            aligned_units.append(reference)
            continue
        alignment = aligner.align(reference, unit)[0]
        mapped = ["-"] * len(reference)
        for (reference_start, reference_end), (unit_start, unit_end) in zip(
            alignment.aligned[0], alignment.aligned[1], strict=True
        ):
            for offset in range(reference_end - reference_start):
                mapped[int(reference_start) + offset] = unit[int(unit_start) + offset]
        aligned_units.append("".join(mapped))
    return aligned_units


def apply_natural_middle_unit_annotations(
    frame: pd.DataFrame,
    *,
    output_path: str | Path | None = None,
    unit_alignment_path: str | Path | None = None,
    position_variability_path: str | Path | None = None,
    workers: int = 12,
    timeout: int = 30,
    fixed_threshold: float = 0.8,
) -> pd.DataFrame:
    """Use the middle RepeatsDB unit and align all units in its region.

    This is the explicit fallback for natural proteins whose divergent repeat
    sequences do not support a shorter equal-period Fourier call.  The full
    sequence scan remains recorded; the database annotation is not presented
    as a spectral discovery. The selected unit has the midpoint closest to the
    repeat-region midpoint; an exact tie is resolved toward the earlier unit.
    """
    source_payloads = [
        (
            str(row.source_annotation_id),
            str(row.source_accession),
            str(row.source_chain),
            timeout,
        )
        for row in frame.drop_duplicates("source_annotation_id").itertuples(index=False)
    ]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        annotations = dict(
            executor.map(
                _fetch_natural_annotation_source,
                source_payloads,
            )
        )
    rows: list[dict[str, object]] = []
    unit_rows: list[dict[str, object]] = []
    position_rows: list[dict[str, object]] = []
    for source_row in frame.to_dict(orient="records"):
        row = dict(source_row)
        source = annotations[str(row["source_annotation_id"])]
        annotation = source["annotation"]
        loci = annotation["content"]["loci"]
        prior_start = int(row["prior_unit_start"])
        prior_end = int(row["prior_unit_end"])
        matching = [
            (index, locus)
            for index, locus in enumerate(loci)
            if locus.get("type") == "unit"
            and int(locus["start"]) == prior_start
            and int(locus["end"]) == prior_end
        ]
        if len(matching) != 1 or matching[0][1].get("parent") is None:
            row["boundary_refinement_status"] = "source_annotation_unit_not_found"
            rows.append(row)
            continue
        parent_index = int(matching[0][1]["parent"])
        region = loci[parent_index]
        region_units = sorted(
            [
                locus
                for locus in loci
                if locus.get("type") == "unit"
                and locus.get("parent") is not None
                and int(locus["parent"]) == parent_index
            ],
            key=lambda value: int(value["start"]),
        )
        full_sequence = str(source["canonical_sequence"])
        auth_mapping_raw = list(source["auth_mapping"])
        auth_mapping = [_auth_coordinate(value) for value in auth_mapping_raw]
        if any(value is None for value in auth_mapping):
            row["boundary_refinement_status"] = "rcsb_auth_mapping_unparseable"
            rows.append(row)
            continue
        auth_mapping_int = [int(value) for value in auth_mapping if value is not None]
        origin = auth_mapping_int[0]
        sequences: list[str] = []
        coordinate_rows: list[tuple[int, int]] = []
        for locus in region_units:
            start = int(locus["start"])
            end = int(locus["end"])
            indices = [
                index
                for index, coordinate in enumerate(auth_mapping_int)
                if start <= coordinate <= end
            ]
            if not indices:
                sequences = []
                break
            sequence = "".join(full_sequence[index] for index in indices)
            if not sequence or set(sequence) - set(AMINO_ACIDS):
                sequences = []
                break
            sequences.append(sequence)
            coordinate_rows.append((start, end))
        if len(sequences) < 2:
            row["boundary_refinement_status"] = "source_annotation_units_outside_sequence"
            rows.append(row)
            continue
        region_start = int(region["start"])
        region_end = int(region["end"])
        region_midpoint = (region_start + region_end) / 2
        selected_index_zero_based = min(
            range(len(coordinate_rows)),
            key=lambda index: (
                abs(
                    (coordinate_rows[index][0] + coordinate_rows[index][1]) / 2
                    - region_midpoint
                ),
                coordinate_rows[index][0],
            ),
        )
        aligned_units = _align_units_to_reference(sequences, selected_index_zero_based)
        first_sequence = sequences[0]
        selected_sequence = sequences[selected_index_zero_based]
        conservation: list[float] = []
        consensus: list[str] = []
        variants_by_position: list[list[str]] = []
        for column in zip(*aligned_units, strict=True):
            counts: dict[str, int] = defaultdict(int)
            for amino_acid in column:
                counts[amino_acid] += 1
            best_amino_acid, best_count = sorted(
                counts.items(), key=lambda item: (-item[1], item[0] == "-", item[0])
            )[0]
            consensus.append(best_amino_acid)
            conservation.append(best_count / len(aligned_units))
            variants_by_position.append(sorted(counts))
        fixed_positions = tuple(
            index + 1 for index, value in enumerate(conservation) if value >= fixed_threshold
        )
        variable_positions = tuple(
            index + 1 for index, value in enumerate(conservation) if value < fixed_threshold
        )
        variable_ranges: list[tuple[int, int]] = []
        for position in variable_positions:
            if not variable_ranges or position != variable_ranges[-1][1] + 1:
                variable_ranges.append((position, position))
            else:
                variable_ranges[-1] = (variable_ranges[-1][0], position)
        first_start, first_end = coordinate_rows[0]
        selected_start, selected_end = coordinate_rows[selected_index_zero_based]
        region_indices = [
            index
            for index, coordinate in enumerate(auth_mapping_int)
            if region_start <= coordinate <= region_end
        ]
        if not region_indices:
            row["boundary_refinement_status"] = "source_region_outside_rcsb_mapping"
            rows.append(row)
            continue
        local_region_start = min(region_indices)
        local_region_end = max(region_indices) + 1
        initial_source_unit_sequence = row.get("prior_unit_sequence", row.get("unit_sequence"))
        initial_source_unit_start = row.get("prior_unit_start", row.get("unit_start"))
        initial_source_unit_end = row.get("prior_unit_end", row.get("unit_end"))
        original_module_id = str(row["module_id"])
        row["module_id"] = (
            f"natural_{row['source_accession']}_{row['source_chain']}_"
            f"{selected_start}_{selected_end}"
        )
        row.update(
            {
                "full_sequence": full_sequence,
                "full_sequence_origin": origin,
                "full_sequence_sha256": hashlib.sha256(full_sequence.encode()).hexdigest(),
                "full_sequence_auth_mapping_json": json.dumps(auth_mapping_raw),
                "full_sequence_source": "RCSB polymer entity canonical sequence",
                "rcsb_entity_id": source["entity_id"],
                "rcsb_auth_chain_id": source["auth_chain_id"],
                "rcsb_label_chain_id": source["label_chain_id"],
                "rcsb_instance_url": source["rcsb_instance_url"],
                "rcsb_entity_url": source["rcsb_entity_url"],
                "source_candidate_module_id": original_module_id,
                "unit_sequence": selected_sequence,
                "unit_length": len(selected_sequence),
                "unit_start": selected_start,
                "unit_end": selected_end,
                "prior_unit_sequence": selected_sequence,
                "prior_unit_length": len(selected_sequence),
                "prior_unit_start": selected_start,
                "prior_unit_end": selected_end,
                "initial_source_candidate_unit_sequence": initial_source_unit_sequence,
                "initial_source_candidate_unit_start": initial_source_unit_start,
                "initial_source_candidate_unit_end": initial_source_unit_end,
                "primitive_period": len(selected_sequence),
                "period": len(selected_sequence),
                "first_module_sequence": first_sequence,
                "first_module_start": first_start,
                "first_module_end": first_end,
                "selected_module_sequence": selected_sequence,
                "selected_module_start": selected_start,
                "selected_module_end": selected_end,
                "selected_module_index": selected_index_zero_based + 1,
                "selected_module_count": len(sequences),
                "selected_module_policy": MODULE_SELECTION_POLICY,
                "module_selection_policy": MODULE_SELECTION_POLICY,
                "repeat_region_start": region_start,
                "repeat_region_end": region_end,
                "repeat_count": len(sequences),
                "unit_sequences_json": json.dumps(sequences),
                "aligned_unit_sequences_json": json.dumps(aligned_units),
                "source_unit_coordinates_json": json.dumps(coordinate_rows),
                "consensus_sequence": "".join(consensus),
                "position_conservation_json": json.dumps(conservation),
                "fixed_mask": "".join(
                    "F" if value >= fixed_threshold else "V" for value in conservation
                ),
                "fixed_positions_json": json.dumps(fixed_positions),
                "variable_positions_json": json.dumps(variable_positions),
                "variable_ranges_json": json.dumps(variable_ranges),
                "left_flank_sequence": full_sequence[:local_region_start],
                "right_flank_sequence": full_sequence[local_region_end:],
                "boundary_refinement_status": "source_annotation_middle_unit",
                "fixed_variable_assessment_status": "aligned_source_annotated_units",
                "periodicity_confidence": "source_annotation",
                "selection_reason": (
                    "middle RepeatsDB unit selected from the annotated repeat region"
                ),
                "boundary_method": (
                    f"{BOUNDARY_METHOD_VERSION}; middle unit and repeat region from RepeatsDB; "
                    "all annotated units aligned to the selected middle unit for fixed/variable positions"
                ),
                "boundary_method_version": BOUNDARY_METHOD_VERSION,
                "source_region_start": region_start,
                "source_region_end": region_end,
            }
        )
        rows.append(row)
        for repeat_index, (sequence, aligned, coordinates) in enumerate(
            zip(sequences, aligned_units, coordinate_rows, strict=True), start=1
        ):
            unit_rows.append(
                {
                    "module_id": row["module_id"],
                    "repeat_index": repeat_index,
                    "unit_start": coordinates[0],
                    "unit_end": coordinates[1],
                    "unit_sequence": sequence,
                    "aligned_to_selected_sequence": aligned,
                    "is_first_module": repeat_index == 1,
                    "is_selected_module": repeat_index == selected_index_zero_based + 1,
                    "module_selection_policy": MODULE_SELECTION_POLICY,
                    "boundary_source": "RepeatsDB",
                }
            )
        for position, (amino_acid, value, variants) in enumerate(
            zip(consensus, conservation, variants_by_position, strict=True), start=1
        ):
            position_rows.append(
                {
                    "module_id": row["module_id"],
                    "module_position": position,
                    "consensus_amino_acid": amino_acid,
                    "conservation": value,
                    "fixed": value >= fixed_threshold,
                    "variants_json": json.dumps(variants),
                    "boundary_source": "RepeatsDB alignment",
                    "reference_repeat_index": selected_index_zero_based + 1,
                }
            )
    refined = pd.DataFrame(rows)
    for path_like, output_frame in (
        (output_path, refined),
        (unit_alignment_path, pd.DataFrame(unit_rows)),
        (position_variability_path, pd.DataFrame(position_rows)),
    ):
        if path_like is None:
            continue
        destination = Path(path_like)
        destination.parent.mkdir(parents=True, exist_ok=True)
        output_frame.to_parquet(destination, index=False)
        output_frame.to_csv(destination.with_suffix(".csv"), index=False)
    return refined


def merge_module_catalogs(frames: Iterable[pd.DataFrame], output_path: str | Path) -> pd.DataFrame:
    frame = pd.concat(list(frames), ignore_index=True)
    finalized_statuses = {
        "source_annotation_middle_unit",
        "strict_dual_evidence_passed",
    }
    already_final = (
        "boundary_refinement_status" in frame
        and frame["boundary_refinement_status"].isin(finalized_statuses).all()
    )
    if not already_final and (
        "boundary_method_version" not in frame
        or not frame["boundary_method_version"].eq(BOUNDARY_METHOD_VERSION).all()
    ):
        frame = refine_module_boundaries(frame)
    frame["retrieved_date"] = frame.get("retrieved_date", date.today().isoformat())
    frame["download_date"] = frame.get("download_date", frame["retrieved_date"])
    frame["download_date"] = frame["download_date"].fillna(frame["retrieved_date"])
    frame["unit_sequence"] = frame["unit_sequence"].map(validate_protein_sequence)
    frame["unit_length"] = frame["unit_sequence"].str.len()
    frame["sequence_sha256"] = frame["unit_sequence"].map(lambda value: hashlib.sha256(value.encode()).hexdigest())
    mappings = frame.sort_values(["unit_sequence", "collection", "evidence_tier", "family", "module_id"]).copy()
    canonical_by_collection_sequence = {
        (row.collection, row.unit_sequence): row.module_id
        for row in mappings.drop_duplicates(["collection", "unit_sequence"], keep="first").itertuples(index=False)
    }
    mappings["canonical_module_id"] = [
        canonical_by_collection_sequence[(collection, sequence)]
        for collection, sequence in zip(mappings["collection"], mappings["unit_sequence"])
    ]
    frame = mappings.drop_duplicates(["collection", "unit_sequence"], keep="first").copy()
    designed = frame[frame["collection"].eq("designed_all")].copy()
    designed["_dhr_number"] = designed["source_accession"].astype(str).str.extract(r"DHR\D*(\d+)", expand=False).fillna("999999").astype(int)
    designed["_evidence_rank"] = designed["evidence_tier"].map({"A": 0, "B": 1, "C": 2}).fillna(3)
    family_upper = designed["family"].str.upper()
    designed["_primary_rank"] = 3
    designed.loc[family_upper.str.contains("DHR") & designed["_dhr_number"].between(1, 83), "_primary_rank"] = 0
    designed.loc[family_upper.str.contains("THR") & designed["evidence_tier"].isin(["A", "B"]), "_primary_rank"] = 1
    designed.loc[family_upper.str.contains("THR") & designed["_primary_rank"].eq(3), "_primary_rank"] = 2
    primary_ids = set(
        designed.sort_values(["_primary_rank", "_dhr_number", "_evidence_rank", "module_id"]).head(100)["module_id"]
    )
    frame["in_designed_primary100"] = frame["module_id"].isin(primary_ids)
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(destination, index=False)
    frame.to_csv(destination.with_suffix(".csv"), index=False)
    mappings_path = destination.with_name(destination.stem + "_source_mappings.parquet")
    mappings.to_parquet(mappings_path, index=False)
    mappings.to_csv(destination.with_name(destination.stem + "_source_mappings.csv"), index=False)
    manifest = {
        "created_at": utc_now(),
        "rows": int(len(frame)),
        "source_mapping_rows": int(len(mappings)),
        "collection_counts": {
            str(key): int(value)
            for key, value in frame["collection"].value_counts().sort_index().items()
        },
        "deduplication_key": ["collection", "unit_sequence"],
        "duplicate_source_rows_collapsed": int(len(mappings) - len(frame)),
        "catalog_sha256": sha256_file(destination),
        "source_mappings_sha256": sha256_file(mappings_path),
    }
    write_json_atomic(manifest, destination.with_suffix(".manifest.json"))
    return frame
