#!/usr/bin/env python3
"""Expand the designed inventory with exact, structure-verified PDB constructs."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from hurdler.constants import validate_protein_sequence
from hurdler.structural_repeats import DESIGNED_CORPUS_VERSION


PDB_FAMILIES = {
    "designed_ankyrin": [
        "1N0Q", "1N0R", "4GPM", "4HQD", "4GMR", "4HB5",
    ],
    "designed_armadillo": [
        "4HXT", "4DB6", "4DB8", "4DB9", "4DBA", "4V3O", "4V3Q",
        "4V3R", "4PLQ", "4PLR", "4PLS", "4D49", "4D4E", "5AEI",
        "5MFB", "5MFC", "5MFD", "5MFE", "5MFF", "5MFG", "5MFH",
        "5MFI", "5MFJ", "5MFK", "5MFL", "5MFM", "5MFN", "5MFO",
        "7QNP", "7R0R",
    ],
    "designed_leucine_rich_repeat": [
        "4PQ8", "4PSJ", "4R58", "4R5C", "4R5D", "4R6F", "4R6G", "4R6J",
    ],
    "designed_tetratricopeptide_repeat": ["1NA0", "1NA3", "2AVP", "3KD7"],
    "designed_ice_binding_iTHR": ["9MG8", "9D01"],
}


def _session() -> requests.Session:
    client = requests.Session()
    retry = Retry(
        total=4,
        connect=4,
        read=4,
        status=4,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
        respect_retry_after_header=True,
    )
    client.mount("https://", HTTPAdapter(max_retries=retry))
    return client


def _get_json(client: requests.Session, url: str) -> dict[str, Any]:
    response = client.get(url, timeout=60)
    response.raise_for_status()
    return response.json()


def _choose_polymer_entity(
    client: requests.Session, accession: str
) -> tuple[str, dict[str, Any], dict[str, Any]]:
    entry = _get_json(
        client, f"https://data.rcsb.org/rest/v1/core/entry/{accession}"
    )
    identifiers = entry.get("rcsb_entry_container_identifiers") or {}
    entity_ids = identifiers.get("polymer_entity_ids") or []
    if not entity_ids:
        raise ValueError("RCSB entry has no polymer entity")
    candidates = []
    for entity_id in entity_ids:
        entity = _get_json(
            client,
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{accession}/{entity_id}",
        )
        polymer = entity.get("entity_poly") or {}
        if polymer.get("rcsb_entity_polymer_type") != "Protein":
            continue
        description = str(
            (entity.get("rcsb_polymer_entity") or {}).get("pdbx_description")
            or ""
        )
        sequence = re.sub(
            r"\s+", "", str(polymer.get("pdbx_seq_one_letter_code_can") or "")
        )
        keyword_score = sum(
            token in description.lower()
            for token in (
                "design", "engineer", "repeat", "armadillo", "ankyrin",
                "leucine", "tetratricopeptide", "ithr", "ice-binding",
            )
        )
        exclusion_score = sum(
            token in description.lower()
            for token in ("green fluorescent", "lysozyme", "peptide ligand")
        )
        candidates.append(
            (keyword_score - 2 * exclusion_score, len(sequence), str(entity_id), entity)
        )
    if not candidates:
        raise ValueError("RCSB entry has no protein entity")
    _, _, entity_id, entity = max(candidates, key=lambda item: (item[0], item[1]))
    return entity_id, entity, entry


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-inventory", type=Path, required=True)
    parser.add_argument("--structure-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = pd.read_parquet(args.base_inventory)
    structure_dir = args.structure_dir.resolve()
    structure_dir.mkdir(parents=True, exist_ok=True)
    client = _session()
    rows: list[dict[str, Any]] = []
    exclusions: list[dict[str, Any]] = []
    for family, accessions in PDB_FAMILIES.items():
        for accession in accessions:
            try:
                entity_id, entity, entry = _choose_polymer_entity(client, accession)
                polymer = entity["entity_poly"]
                sequence = validate_protein_sequence(
                    re.sub(
                        r"\s+",
                        "",
                        str(polymer["pdbx_seq_one_letter_code_can"]),
                    )
                )
                identifiers = entity.get(
                    "rcsb_polymer_entity_container_identifiers"
                ) or {}
                chains = identifiers.get("auth_asym_ids") or identifiers.get(
                    "asym_ids"
                ) or []
                structure_path = structure_dir / f"{accession}.cif"
                response = client.get(
                    f"https://files.rcsb.org/download/{accession}.cif", timeout=120
                )
                response.raise_for_status()
                structure_path.write_bytes(response.content)
                description = str(
                    (entity.get("rcsb_polymer_entity") or {}).get(
                        "pdbx_description"
                    )
                    or accession
                )
                citation = entry.get("rcsb_primary_citation") or {}
                doi = str(citation.get("pdbx_database_id_DOI") or "")
                rows.append(
                    {
                        "module_id": f"designed_pdb_{accession}",
                        "canonical_module_id": f"designed_pdb_{accession}",
                        "collection": "designed_all",
                        "module_type": "Designed",
                        "family": family,
                        "unit_sequence": sequence,
                        "unit_length": len(sequence),
                        "prior_unit_sequence": "",
                        "prior_unit_length": None,
                        "full_sequence": sequence,
                        "full_sequence_length": len(sequence),
                        "full_sequence_origin": 1,
                        "full_sequence_sha256": hashlib.sha256(
                            sequence.encode()
                        ).hexdigest(),
                        "full_sequence_source": "RCSB polymer entity canonical sequence",
                        "evidence_tier": "A",
                        "source_name": "RCSB Protein Data Bank",
                        "source_url": f"https://www.rcsb.org/structure/{accession}",
                        "source_accession": accession,
                        "source_chain": str(chains[0]) if chains else "",
                        "source_annotation_id": f"PDB:{accession}:entity:{entity_id}",
                        "source_sha256": hashlib.sha256(response.content).hexdigest(),
                        "citation": doi or description,
                        "license_name": "wwPDB data terms",
                        "reviewed": True,
                        "retrieved_date": date.today().isoformat(),
                        "download_date": date.today().isoformat(),
                        "notes": description,
                        "author_structure_path": str(structure_path),
                        "structure_candidate_paths_json": json.dumps(
                            [str(structure_path)]
                        ),
                        "structure_inventory_status": "author_or_pdb_structure_available",
                        "boundary_method": "pending strict DSSP/Foldseek inference",
                        "boundary_refinement_status": "pending_strict_dual_evidence",
                        "corpus_version": DESIGNED_CORPUS_VERSION,
                    }
                )
            except Exception as exc:
                exclusions.append(
                    {
                        "source_accession": accession,
                        "family": family,
                        "status": "designed_pdb_acquisition_failed",
                        "reason": f"{type(exc).__name__}: {exc}",
                    }
                )
    additions = pd.DataFrame(rows)
    combined = pd.concat([base, additions], ignore_index=True, sort=False)
    combined = combined.sort_values("module_id", kind="mergesort").drop_duplicates(
        "module_id", keep="first"
    )
    destination = args.output.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(destination, index=False)
    combined.to_csv(destination.with_suffix(".csv"), index=False)
    additions.to_parquet(
        destination.with_name(destination.stem + "_pdb_additions.parquet"), index=False
    )
    additions.to_csv(
        destination.with_name(destination.stem + "_pdb_additions.csv"), index=False
    )
    pd.DataFrame(
        exclusions,
        columns=["source_accession", "family", "status", "reason"],
    ).to_csv(
        destination.with_name(destination.stem + "_pdb_exclusions.csv"), index=False
    )
    destination.with_suffix(".manifest.json").write_text(
        json.dumps(
            {
                "corpus_version": DESIGNED_CORPUS_VERSION,
                "base_rows": len(base),
                "requested_pdb_rows": sum(map(len, PDB_FAMILIES.values())),
                "added_pdb_rows": len(additions),
                "exclusion_rows": len(exclusions),
                "combined_rows": len(combined),
                "families": {key: len(value) for key, value in PDB_FAMILIES.items()},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
