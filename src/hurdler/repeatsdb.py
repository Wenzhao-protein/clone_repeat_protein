"""Direct RepeatsDB module acquisition for the expanded natural corpus.

Natural repeat boundaries in this module are database observations, never
periodicity inferences.  RepeatsDB uses two public annotation layouts:

* explicit ``region``/``unit`` records in ``content.loci``;
* mapped ``RepeatsDB-*`` feature loci, most often on AlphaFoldDB entries.

Both are normalized to the same row contract before selecting one longest
repeat region per biological protein and its earlier middle unit.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import requests

from .constants import AMINO_ACIDS, validate_protein_sequence

REPEATSDB_API = "https://repeatsdb.org/api/production"
NATURAL_CORPUS_VERSION = "expanded-middle-repeatsdb-foldseek-v1"
NATURAL_BOUNDARY_METHOD = "repeatsdb-direct-middle-unit-v1"
NATURAL_SELECTION_POLICY = "one-protein-longest-region-middle-unit-tie-earlier-v1"


def _mkdir_stable(path: Path) -> None:
    """Create a shared directory despite transient parallel-filesystem races."""
    for attempt in range(5):
        try:
            path.mkdir(parents=True, exist_ok=True)
            return
        except FileExistsError:
            if path.is_dir():
                return
            if attempt == 4:
                raise
            time.sleep(0.1 * (attempt + 1))


def _integer(value: object) -> int:
    match = re.match(r"^-?\d+", str(value))
    if match is None:
        raise ValueError(f"Coordinate is not an integer: {value!r}")
    return int(match.group())


def _uniprot_accessions(content: dict[str, Any]) -> tuple[str, ...]:
    values = {str(value) for value in content.get("features_uniprot", []) if value}
    for key, value in (content.get("features") or {}).items():
        if str(key).lower().startswith("uniprot-"):
            accession = (
                value.get("uniprot_id")
                if isinstance(value, dict)
                else str(key).split("-", 1)[-1]
            )
            if accession:
                values.add(str(accession))
    return tuple(sorted(values))


def _normalized_locus(value: dict[str, Any]) -> dict[str, Any]:
    return {
        **value,
        "type": str(value.get("type", "")),
        "start": _integer(value["start"]),
        "end": _integer(value["end"]),
    }


def annotation_repeat_regions(annotation: dict[str, Any]) -> list[dict[str, Any]]:
    """Normalize every directly annotated repeat region in one API record."""
    content = annotation.get("content") or {}
    chain = content.get("chain") or {}
    annotation_uuid = str(annotation.get("uuid") or annotation.get("_id") or "")
    uniprot_accessions = set(_uniprot_accessions(content))
    if "ALPHAFOLD" in str(chain.get("source") or "").upper() and chain.get(
        "structure"
    ):
        # AlphaFoldDB structure accessions are canonical UniProt accessions,
        # even when the optional features_uniprot array is empty.
        uniprot_accessions.add(str(chain["structure"]))
    base = {
        "annotation_uuid": annotation_uuid,
        "structure_accession": str(chain.get("structure") or ""),
        "structure_chain": str(chain.get("id") or ""),
        "structure_source": str(chain.get("source") or ""),
        "reviewed": bool(content.get("reviewed", False)),
        "uniprot_accessions_json": json.dumps(sorted(uniprot_accessions)),
        "annotation_updated": str((annotation.get("updated") or {}).get("at") or ""),
    }
    loci = [_normalized_locus(value) for value in (content.get("loci") or [])]
    regions: list[dict[str, Any]] = []
    for region_index, region in enumerate(loci):
        if region["type"] != "region":
            continue
        units = sorted(
            [
                value
                for value in loci
                if value["type"] == "unit"
                and value.get("parent") is not None
                and int(value["parent"]) == region_index
            ],
            key=lambda value: (value["start"], value["end"]),
        )
        if len(units) < 2:
            continue
        insertions = sorted(
            [
                value
                for value in loci
                if value["type"] == "insertion"
                and value.get("parent") is not None
                and int(value["parent"]) == region_index
            ],
            key=lambda value: (value["start"], value["end"]),
        )
        regions.append(
            {
                **base,
                "region_locator": f"content.loci[{region_index}]",
                "annotation_schema": "content_loci",
                "region_start": region["start"],
                "region_end": region["end"],
                "classification": str(region.get("class") or ""),
                "unit_coordinates_json": json.dumps(
                    [[unit["start"], unit["end"]] for unit in units]
                ),
                "insertion_coordinates_json": json.dumps(
                    [[item["start"], item["end"]] for item in insertions]
                ),
                "unit_count": len(units),
            }
        )
    if regions:
        return regions

    # Mapped AlphaFoldDB annotations can have no top-level region.  Each
    # RepeatsDB feature represents one mapped source annotation and therefore
    # one candidate region.  Insertions are retained but never treated as units.
    for feature_name, feature in sorted((content.get("features") or {}).items()):
        if not str(feature_name).startswith("RepeatsDB-") or not isinstance(feature, dict):
            continue
        feature_loci = [
            _normalized_locus(value) for value in (feature.get("loci") or [])
        ]
        units = sorted(
            [value for value in feature_loci if value["type"] == "unit"],
            key=lambda value: (value["start"], value["end"]),
        )
        if len(units) < 2:
            continue
        insertions = sorted(
            [value for value in feature_loci if value["type"] == "insertion"],
            key=lambda value: (value["start"], value["end"]),
        )
        region_loci = [*units, *insertions]
        feature_chain = feature.get("chain") or {}
        regions.append(
            {
                **base,
                "region_locator": f"content.features[{feature_name!r}]",
                "annotation_schema": "mapped_repeatsdb_feature",
                "region_start": min(value["start"] for value in region_loci),
                "region_end": max(value["end"] for value in region_loci),
                "classification": str(feature.get("class") or ""),
                "mapped_feature_name": str(feature_name),
                "mapped_from_structure": str(feature_chain.get("structure") or ""),
                "mapped_from_chain": str(feature_chain.get("id") or ""),
                "unit_coordinates_json": json.dumps(
                    [[unit["start"], unit["end"]] for unit in units]
                ),
                "insertion_coordinates_json": json.dumps(
                    [[item["start"], item["end"]] for item in insertions]
                ),
                "unit_count": len(units),
            }
        )
    return regions


def select_longest_region_per_protein(frame: pd.DataFrame) -> pd.DataFrame:
    """Select exactly one longest annotated region for each protein key."""
    if frame.empty:
        return frame.copy()
    ranked = frame.copy()
    ranked["region_span"] = ranked.region_end.astype(int) - ranked.region_start.astype(int) + 1
    ranked["_source_rank"] = ranked.structure_source.astype(str).str.upper().map(
        lambda value: 0 if "PDB" in value else 1 if "ALPHAFOLD" in value else 2
    )
    ranked = ranked.sort_values(
        [
            "protein_key",
            "region_span",
            "unit_count",
            "reviewed",
            "_source_rank",
            "region_start",
            "annotation_uuid",
            "region_locator",
        ],
        ascending=[True, False, False, False, True, True, True, True],
        kind="mergesort",
    )
    selected = ranked.drop_duplicates("protein_key", keep="first").drop(
        columns=["_source_rank"]
    )
    return selected.reset_index(drop=True)


def _parse_fasta(text: str) -> tuple[str, int]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    header = next((line for line in lines if line.startswith(">")), "")
    match = re.search(r"\brange\s+(-?\d+):(-?\d+)", header)
    origin = int(match.group(1)) if match else 1
    sequence = validate_protein_sequence(
        "".join(line for line in lines if not line.startswith(">"))
    )
    return sequence, origin


def _slice_units(
    sequence: str, origin: int, coordinates: Iterable[Iterable[int]]
) -> list[str]:
    units: list[str] = []
    for start_value, end_value in coordinates:
        start, end = int(start_value), int(end_value)
        local_start, local_end = start - origin, end - origin + 1
        if local_start < 0 or local_end > len(sequence) or local_end <= local_start:
            raise ValueError(
                f"Annotated unit {start}:{end} is outside FASTA range "
                f"{origin}:{origin + len(sequence) - 1}"
            )
        units.append(validate_protein_sequence(sequence[local_start:local_end]))
    return units


def _slice_pdb_units_with_rcsb_mapping(
    row: dict[str, Any], coordinates: list[list[int]], timeout: int
) -> tuple[list[str], str, int, str]:
    """Recover PDB units when author coordinates are not contiguous FASTA indices."""
    # Local import avoids a module-import cycle while reusing the repository's
    # tested author-chain to entity-poly mapping implementation.
    from .modules import _auth_coordinate, _fetch_natural_annotation_source

    _, source = _fetch_natural_annotation_source(
        (
            str(row["annotation_uuid"]),
            str(row["structure_accession"]),
            str(row["structure_chain"]),
            timeout,
        )
    )
    sequence = validate_protein_sequence(str(source["canonical_sequence"]))
    mapping = [_auth_coordinate(value) for value in source["auth_mapping"]]
    units: list[str] = []
    for start_value, end_value in coordinates:
        start, end = int(start_value), int(end_value)
        indices = [
            index
            for index, coordinate in enumerate(mapping)
            if coordinate is not None and start <= coordinate <= end
        ]
        if not indices:
            raise ValueError(f"RCSB mapping contains no residues for {start}:{end}")
        units.append(validate_protein_sequence("".join(sequence[index] for index in indices)))
    origin = min(value for value in mapping if value is not None)
    return units, sequence, origin, json.dumps(source["auth_mapping"])


def _request_json(url: str, *, params: dict[str, Any] | None, timeout: int) -> Any:
    error: Exception | None = None
    for attempt in range(4):
        try:
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, ValueError) as exc:
            error = exc
            if attempt < 3:
                time.sleep(0.5 * (2**attempt))
    assert error is not None
    raise error


def enumerate_repeatsdb_annotations(
    *,
    sources: tuple[str, ...] = ("PDB", "AlphaFoldDB"),
    page_size: int = 100,
    timeout: int = 30,
    max_annotations: int | None = None,
) -> list[dict[str, Any]]:
    """Walk all requested public API pages in stable source/skip order."""
    if not 1 <= page_size <= 100:
        raise ValueError("RepeatsDB page_size must be between 1 and 100")
    annotations: list[dict[str, Any]] = []
    for source in sources:
        skip = 0
        total: int | None = None
        while total is None or skip < total:
            payload = _request_json(
                f"{REPEATSDB_API}/annotations",
                params={"limit": page_size, "skip": skip, "chain.source": source},
                timeout=timeout,
            )
            total = int(payload["count"])
            items = payload.get("items") or []
            page = list(items.values()) if isinstance(items, dict) else list(items)
            annotations.extend(page)
            if max_annotations is not None and len(annotations) >= max_annotations:
                return annotations[:max_annotations]
            skip += page_size
    return annotations


def write_annotation_inventory(
    output_path: str | Path,
    *,
    timeout: int = 30,
    max_annotations: int | None = None,
) -> pd.DataFrame:
    """Materialize the API metadata once so CPU FASTA shards do not rewalk it."""
    annotations = enumerate_repeatsdb_annotations(
        timeout=timeout, max_annotations=max_annotations
    )
    rows = []
    for annotation in annotations:
        content = annotation.get("content") or {}
        chain = content.get("chain") or {}
        rows.append(
            {
                "annotation_uuid": str(
                    annotation.get("uuid") or annotation.get("_id") or ""
                ),
                "structure_source": str(chain.get("source") or ""),
                "structure_accession": str(chain.get("structure") or ""),
                "structure_chain": str(chain.get("id") or ""),
                "annotation_json": json.dumps(
                    annotation, sort_keys=True, separators=(",", ":")
                ),
            }
        )
    frame = pd.DataFrame(rows).sort_values(
        ["structure_source", "annotation_uuid"], kind="mergesort"
    )
    if frame.annotation_uuid.duplicated().any():
        raise ValueError("RepeatsDB inventory contains duplicate annotation UUIDs")
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(destination, index=False)
    frame.to_csv(destination.with_suffix(".csv.gz"), index=False, compression="gzip")
    destination.with_suffix(".manifest.json").write_text(
        json.dumps(
            {
                "corpus_version": NATURAL_CORPUS_VERSION,
                "annotation_count": len(frame),
                "source_counts": frame.structure_source.value_counts().to_dict(),
                "api": REPEATSDB_API,
                "retrieved_date": date.today().isoformat(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return frame.reset_index(drop=True)


def _materialize_annotation(
    payload: tuple[dict[str, Any], Path | None, int]
) -> tuple[str, str, int, str, str, str, str]:
    annotation, cache_dir, timeout = payload
    annotation_uuid = str(annotation.get("uuid") or annotation.get("_id"))
    cache_path = cache_dir / f"{annotation_uuid}.fasta" if cache_dir else None
    cache_metadata_path = (
        cache_path.with_suffix(".metadata.json") if cache_path is not None else None
    )
    sequence_source = "repeatsdb_fasta"
    try:
        if cache_path is not None and cache_path.exists() and cache_path.stat().st_size:
            text = cache_path.read_text()
            sequence_source = "cached_sequence"
            if cache_metadata_path is not None and cache_metadata_path.is_file():
                try:
                    sequence_source = str(
                        json.loads(cache_metadata_path.read_text()).get(
                            "full_sequence_source", sequence_source
                        )
                    )
                except (ValueError, json.JSONDecodeError):
                    pass
            if sequence_source == "cached_sequence":
                header = text.splitlines()[0] if text.splitlines() else ""
                if re.search(r"\brange\s+-?\d+:-?\d+", header):
                    sequence_source = "repeatsdb_fasta"
                elif header.startswith((">sp|", ">tr|", ">UniRef")):
                    sequence_source = "uniprot_canonical_fallback"
                elif header.startswith(">RCSB fallback"):
                    sequence_source = "rcsb_entity_poly_fallback"
        else:
            response = requests.get(
                f"{REPEATSDB_API}/annotations/{annotation_uuid}",
                params={"format": "fasta"},
                timeout=timeout,
            )
            response.raise_for_status()
            text = response.text
        sequence, origin = _parse_fasta(text)
    except (requests.RequestException, ValueError) as primary_error:
        content = annotation.get("content") or {}
        chain = content.get("chain") or {}
        structure_source = str(chain.get("source") or "")
        accession = str(chain.get("structure") or "")
        try:
            if "ALPHAFOLD" in structure_source.upper() and accession:
                response = requests.get(
                    f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
                    timeout=timeout,
                )
                response.raise_for_status()
                text = response.text
                sequence, origin = _parse_fasta(text)
                sequence_source = "uniprot_canonical_fallback"
            elif "PDB" in structure_source.upper() and accession:
                from .modules import _fetch_natural_annotation_source

                _, source = _fetch_natural_annotation_source(
                    (
                        annotation_uuid,
                        accession,
                        str(chain.get("id") or ""),
                        timeout,
                    )
                )
                sequence = validate_protein_sequence(
                    str(source["canonical_sequence"])
                )
                origin = 1
                text = f">RCSB fallback {accession}\n{sequence}\n"
                sequence_source = "rcsb_entity_poly_fallback"
            else:
                raise ValueError(
                    f"No sequence fallback for source {structure_source!r}"
                )
        except Exception as fallback_error:
            return (
                annotation_uuid,
                "",
                0,
                "",
                "",
                "sequence_unavailable",
                f"primary={type(primary_error).__name__}: {primary_error}; "
                f"fallback={type(fallback_error).__name__}: {fallback_error}",
            )
    if cache_path is not None and not cache_path.exists():
        _mkdir_stable(cache_path.parent)
        cache_path.write_text(text)
        assert cache_metadata_path is not None
        cache_metadata_path.write_text(
            json.dumps(
                {
                    "annotation_uuid": annotation_uuid,
                    "full_sequence_source": sequence_source,
                    "content_sha256": hashlib.sha256(text.encode()).hexdigest(),
                },
                sort_keys=True,
            )
            + "\n"
        )
    return (
        annotation_uuid,
        sequence,
        origin,
        hashlib.sha256(text.encode()).hexdigest(),
        hashlib.sha256(sequence.encode()).hexdigest(),
        sequence_source,
        "",
    )


def build_natural_corpus(
    output_path: str | Path,
    *,
    mappings_path: str | Path | None = None,
    exclusions_path: str | Path | None = None,
    cache_dir: str | Path | None = None,
    workers: int = 12,
    timeout: int = 30,
    max_annotations: int | None = None,
    shard_index: int = 0,
    shard_count: int = 1,
    annotation_inventory_path: str | Path | None = None,
) -> pd.DataFrame:
    """Build the no-cap, one-protein RepeatsDB-direct natural corpus."""
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("shard_index must satisfy 0 <= shard_index < shard_count")
    if annotation_inventory_path is not None:
        inventory = pd.read_parquet(annotation_inventory_path)
        if "annotation_json" not in inventory:
            raise ValueError("RepeatsDB inventory lacks annotation_json")
        annotations = [json.loads(value) for value in inventory.annotation_json]
        if max_annotations is not None:
            annotations = annotations[:max_annotations]
    else:
        annotations = enumerate_repeatsdb_annotations(
            timeout=timeout, max_annotations=max_annotations
        )
    candidates: list[dict[str, Any]] = []
    annotation_by_uuid: dict[str, dict[str, Any]] = {}
    for annotation in annotations:
        annotation_uuid = str(annotation.get("uuid") or annotation.get("_id") or "")
        annotation_by_uuid[annotation_uuid] = annotation
        candidates.extend(annotation_repeat_regions(annotation))
    candidate_frame = pd.DataFrame(candidates)
    if candidate_frame.empty:
        raise RuntimeError("RepeatsDB enumeration returned no regions with at least two units")

    annotation_ids = sorted(set(candidate_frame.annotation_uuid.astype(str)))
    retained_ids = set(annotation_ids[shard_index::shard_count])
    candidate_frame = candidate_frame.loc[
        candidate_frame.annotation_uuid.astype(str).isin(retained_ids)
    ].copy()
    cache = Path(cache_dir) if cache_dir is not None else None
    if cache is not None:
        _mkdir_stable(cache)
    unique_annotations = [
        annotation_by_uuid[value]
        for value in sorted(retained_ids)
    ]
    payloads = [(annotation, cache, timeout) for annotation in unique_annotations]
    if workers == 1:
        materialized = list(map(_materialize_annotation, payloads))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            materialized = list(executor.map(_materialize_annotation, payloads))
    sequence_by_uuid = {
        uuid: {
            "full_sequence": sequence,
            "full_sequence_origin": origin,
            "source_fasta_sha256": fasta_sha,
            "full_sequence_sha256": sequence_sha,
            "full_sequence_source": sequence_source,
            "materialization_error": materialization_error,
        }
        for (
            uuid,
            sequence,
            origin,
            fasta_sha,
            sequence_sha,
            sequence_source,
            materialization_error,
        ) in materialized
    }
    for column in (
        "full_sequence",
        "full_sequence_origin",
        "source_fasta_sha256",
        "full_sequence_sha256",
        "full_sequence_source",
        "materialization_error",
    ):
        candidate_frame[column] = candidate_frame.annotation_uuid.map(
            lambda uuid, key=column: sequence_by_uuid[str(uuid)][key]
        )
    materialization_failures = candidate_frame.loc[
        candidate_frame.materialization_error.astype(str).ne("")
    ].copy()
    exclusions: list[dict[str, Any]] = [
        {
            "annotation_uuid": row.annotation_uuid,
            "protein_key": "",
            "status": "source_sequence_unavailable",
            "reason": row.materialization_error,
        }
        for row in materialization_failures.drop_duplicates(
            "annotation_uuid"
        ).itertuples(index=False)
    ]
    candidate_frame = candidate_frame.loc[
        candidate_frame.materialization_error.astype(str).eq("")
    ].copy()
    if candidate_frame.empty:
        raise RuntimeError("No RepeatsDB candidate had a retrievable full sequence")
    candidate_frame["uniprot_accessions"] = candidate_frame.uniprot_accessions_json.map(
        json.loads
    )
    candidate_frame["protein_key"] = candidate_frame.apply(
        lambda row: (
            f"uniprot:{row.uniprot_accessions[0]}"
            if row.uniprot_accessions
            else f"sha256:{row.full_sequence_sha256}"
        ),
        axis=1,
    )
    selected = select_longest_region_per_protein(candidate_frame)

    rows: list[dict[str, Any]] = []
    for row in selected.to_dict(orient="records"):
        coordinates = json.loads(str(row["unit_coordinates_json"]))
        try:
            units = _slice_units(
                str(row["full_sequence"]),
                int(row["full_sequence_origin"]),
                coordinates,
            )
        except ValueError as direct_error:
            if "PDB" in str(row["structure_source"]).upper():
                try:
                    units, mapped_sequence, mapped_origin, auth_mapping_json = (
                        _slice_pdb_units_with_rcsb_mapping(row, coordinates, timeout)
                    )
                    row["full_sequence"] = mapped_sequence
                    row["full_sequence_origin"] = mapped_origin
                    row["full_sequence_sha256"] = hashlib.sha256(
                        mapped_sequence.encode()
                    ).hexdigest()
                    row["full_sequence_auth_mapping_json"] = auth_mapping_json
                    row["sequence_coordinate_method"] = "rcsb_author_to_entity_poly_mapping"
                except (ValueError, RuntimeError, requests.RequestException) as mapping_error:
                    exclusions.append(
                        {
                            "annotation_uuid": row["annotation_uuid"],
                            "protein_key": row["protein_key"],
                            "status": "source_coordinate_extraction_failed",
                            "reason": f"direct={direct_error}; rcsb={mapping_error}",
                        }
                    )
                    continue
            else:
                exclusions.append(
                    {
                        "annotation_uuid": row["annotation_uuid"],
                        "protein_key": row["protein_key"],
                        "status": "source_coordinate_extraction_failed",
                        "reason": str(direct_error),
                    }
                )
                continue
        middle_index = (len(units) - 1) // 2
        start, end = map(int, coordinates[middle_index])
        middle = units[middle_index]
        accession = str(row["structure_accession"])
        chain = str(row["structure_chain"])
        rows.append(
            {
                **row,
                "module_id": (
                    f"natural_{accession}_{chain}_{start}_{end}_"
                    f"{str(row['annotation_uuid'])[:8]}"
                ),
                "collection": "natural_all",
                "module_type": "Natural",
                "family": f"RepeatsDB {row.get('classification') or 'unclassified'}",
                "unit_sequence": middle,
                "unit_length": len(middle),
                "unit_start": start,
                "unit_end": end,
                "selected_module_sequence": middle,
                "selected_module_start": start,
                "selected_module_end": end,
                "selected_module_index": middle_index + 1,
                "selected_module_count": len(units),
                "unit_sequences_json": json.dumps(units),
                "source_unit_coordinates_json": json.dumps(coordinates),
                "repeat_region_start": int(row["region_start"]),
                "repeat_region_end": int(row["region_end"]),
                "repeat_count": len(units),
                "evidence_tier": "A" if "PDB" in str(row["structure_source"]).upper() else "B",
                "source_name": "RepeatsDB v4",
                "source_url": f"https://repeatsdb.org/structure/{accession}",
                "source_accession": accession,
                "source_chain": chain,
                "source_annotation_id": str(row["annotation_uuid"]),
                "source_sha256": str(row["source_fasta_sha256"]),
                "reviewed": bool(row["reviewed"]),
                "retrieved_date": date.today().isoformat(),
                "download_date": date.today().isoformat(),
                "citation": "RepeatsDB v4",
                "license_name": "RepeatsDB terms; underlying structure terms apply",
                "boundary_method": NATURAL_BOUNDARY_METHOD,
                "boundary_method_version": NATURAL_BOUNDARY_METHOD,
                "boundary_refinement_status": "source_annotation_middle_unit",
                "module_selection_policy": NATURAL_SELECTION_POLICY,
                "selected_module_policy": NATURAL_SELECTION_POLICY,
                "corpus_version": NATURAL_CORPUS_VERSION,
                "periodicity_confidence": "source_annotation",
                "selection_reason": (
                    "longest RepeatsDB region for protein; earlier middle annotated unit"
                ),
            }
        )
    selected_rows = pd.DataFrame(rows)
    if selected_rows.empty:
        raise RuntimeError("No selected RepeatsDB region could be materialized")

    selected_rows = selected_rows.sort_values(
        ["unit_sequence", "reviewed", "evidence_tier", "protein_key", "module_id"],
        ascending=[True, False, True, True, True],
        kind="mergesort",
    )
    selected_rows["canonical_module_id"] = selected_rows.groupby(
        "unit_sequence", sort=False
    ).module_id.transform("first")
    catalog = selected_rows.drop_duplicates("unit_sequence", keep="first").copy()

    selected_region_keys = {
        (str(value["annotation_uuid"]), str(value["region_locator"]))
        for value in selected.to_dict(orient="records")
    }
    candidate_frame["selected_for_protein"] = [
        (str(uuid), str(locator)) in selected_region_keys
        for uuid, locator in zip(
            candidate_frame.annotation_uuid, candidate_frame.region_locator, strict=True
        )
    ]
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_parquet(destination, index=False)
    catalog.to_csv(destination.with_suffix(".csv"), index=False)
    mappings_destination = Path(mappings_path) if mappings_path else destination.with_name(
        destination.stem + "_source_mappings.parquet"
    )
    mappings_destination.parent.mkdir(parents=True, exist_ok=True)
    selected_rows.to_parquet(mappings_destination, index=False)
    selected_rows.to_csv(mappings_destination.with_suffix(".csv"), index=False)
    inventory_destination = destination.with_name(
        destination.stem + "_region_inventory.parquet"
    )
    candidate_frame.to_parquet(inventory_destination, index=False)
    candidate_frame.to_csv(inventory_destination.with_suffix(".csv"), index=False)
    exclusion_destination = Path(exclusions_path) if exclusions_path else destination.with_name(
        destination.stem + "_exclusions.csv"
    )
    pd.DataFrame(
        exclusions,
        columns=["annotation_uuid", "protein_key", "status", "reason"],
    ).to_csv(exclusion_destination, index=False)
    manifest = {
        "corpus_version": NATURAL_CORPUS_VERSION,
        "boundary_method": NATURAL_BOUNDARY_METHOD,
        "selection_policy": NATURAL_SELECTION_POLICY,
        "annotation_shard_index": shard_index,
        "annotation_shard_count": shard_count,
        "enumerated_annotation_count": len(annotations),
        "retained_annotation_count": len(retained_ids),
        "materialized_annotation_count": len(unique_annotations)
        - materialization_failures.annotation_uuid.nunique(),
        "sequence_unavailable_annotation_count": materialization_failures.annotation_uuid.nunique(),
        "full_sequence_source_counts": selected_rows.full_sequence_source.value_counts().to_dict(),
        "selected_protein_rows": len(selected_rows),
        "unique_middle_unit_rows": len(catalog),
        "exclusion_rows": len(exclusions),
    }
    destination.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return catalog.reset_index(drop=True)


def finalize_natural_corpus(
    mapping_paths: Iterable[str | Path],
    output_path: str | Path,
    *,
    region_inventory_paths: Iterable[str | Path] = (),
    exclusion_paths: Iterable[str | Path] = (),
    annotation_inventory_path: str | Path | None = None,
) -> pd.DataFrame:
    """Globally resolve protein and exact-unit duplicates across CPU shards."""
    paths = [Path(path) for path in mapping_paths]
    if not paths:
        raise ValueError("At least one natural source-mapping shard is required")
    mappings = pd.concat(
        [pd.read_parquet(path) for path in paths], ignore_index=True
    )
    protein_rows = select_longest_region_per_protein(mappings)
    if protein_rows.protein_key.duplicated().any():
        raise AssertionError("Natural finalization did not yield one row per protein")
    protein_rows = protein_rows.sort_values(
        ["unit_sequence", "reviewed", "evidence_tier", "protein_key", "module_id"],
        ascending=[True, False, True, True, True],
        kind="mergesort",
    )
    protein_rows["canonical_module_id"] = protein_rows.groupby(
        "unit_sequence", sort=False
    ).module_id.transform("first")
    catalog = protein_rows.drop_duplicates("unit_sequence", keep="first").copy()
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_parquet(destination, index=False)
    catalog.to_csv(destination.with_suffix(".csv"), index=False)
    mapping_destination = destination.with_name(
        destination.stem + "_source_mappings.parquet"
    )
    protein_rows.to_parquet(mapping_destination, index=False)
    protein_rows.to_csv(mapping_destination.with_suffix(".csv"), index=False)
    region_files = [Path(path) for path in region_inventory_paths]
    region_rows = (
        pd.concat(
            [pd.read_parquet(path) for path in region_files], ignore_index=True
        )
        if region_files
        else pd.DataFrame()
    )
    if not region_rows.empty:
        selected_region_keys = {
            (str(row.annotation_uuid), str(row.region_locator))
            for row in protein_rows.itertuples(index=False)
        }
        region_rows["selected_for_protein_global"] = [
            (str(uuid), str(locator)) in selected_region_keys
            for uuid, locator in zip(
                region_rows.annotation_uuid,
                region_rows.region_locator,
                strict=True,
            )
        ]
        region_destination = destination.with_name(
            destination.stem + "_all_region_source_mappings.parquet"
        )
        region_rows.to_parquet(region_destination, index=False)
        region_rows.to_csv(region_destination.with_suffix(".csv.gz"), index=False, compression="gzip")
    exclusion_files = [Path(path) for path in exclusion_paths]
    exclusions = (
        pd.concat([pd.read_csv(path) for path in exclusion_files], ignore_index=True)
        if exclusion_files
        else pd.DataFrame(columns=["annotation_uuid", "protein_key", "status", "reason"])
    )
    exclusion_destination = destination.with_name(
        destination.stem + "_exclusions.csv"
    )
    exclusions.to_csv(exclusion_destination, index=False)
    annotation_rows = None
    annotation_source_counts: dict[str, int] = {}
    if annotation_inventory_path is not None:
        annotation_inventory = pd.read_parquet(
            annotation_inventory_path,
            columns=["annotation_uuid", "structure_source"],
        )
        annotation_rows = len(annotation_inventory)
        annotation_source_counts = {
            str(key): int(value)
            for key, value in annotation_inventory.structure_source.value_counts().items()
        }
    manifest = {
        "corpus_version": NATURAL_CORPUS_VERSION,
        "mapping_shards": [str(path.resolve()) for path in paths],
        "region_inventory_shards": [str(path.resolve()) for path in region_files],
        "exclusion_shards": [str(path.resolve()) for path in exclusion_files],
        "annotation_inventory": (
            str(Path(annotation_inventory_path).resolve())
            if annotation_inventory_path is not None
            else None
        ),
        "annotation_inventory_rows": annotation_rows,
        "annotation_source_counts": annotation_source_counts,
        "protein_rows_before_global_selection": len(mappings),
        "one_per_protein_rows": len(protein_rows),
        "unique_middle_unit_rows": len(catalog),
        "all_annotated_region_rows": len(region_rows),
        "coordinate_exclusion_rows": len(exclusions),
    }
    destination.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return catalog.reset_index(drop=True)
