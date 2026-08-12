#!/usr/bin/env python3
"""Build residue-level secondary-structure evidence and jointly select modules."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import requests

from hurdler.modules import merge_module_catalogs, reselect_module_boundaries
from hurdler.periodicity import BOUNDARY_METHOD_VERSION, MODULE_SELECTION_POLICY, RepeatCandidate
from hurdler.secondary_structure import (
    SECONDARY_STRUCTURE_METHOD_VERSION,
    dssp_annotations_for_structure,
    map_author_template_to_full_sequence,
    parse_dhr_secondary_structure_table,
    score_candidates_with_secondary_structure,
    secondary_structure_json,
)


def _write_frame(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    frame.to_csv(path.with_suffix(".csv"), index=False)


def _fetch_rcsb_structure(
    accession: str, cache_dir: Path, timeout: int, *, cache_key: str = ""
) -> Path:
    safe_key = re.sub(r"[^A-Za-z0-9_.-]+", "_", cache_key).strip("_")
    filename = f"{accession.lower()}_{safe_key}.cif" if safe_key else f"{accession.lower()}.cif"
    destination = cache_dir / filename
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    response = None
    for attempt in range(3):
        try:
            response = requests.get(
                f"https://files.rcsb.org/download/{accession.upper()}.cif",
                timeout=timeout,
            )
            response.raise_for_status()
            break
        except requests.RequestException:
            if attempt == 2:
                raise
            time.sleep(0.5 * (2**attempt))
    assert response is not None
    destination.write_bytes(response.content)
    return destination


def _resolve_local_structure(row: dict[str, object], structure_dir: Path) -> Path:
    accession = str(row["source_accession"])
    component_match = re.search(r"_([AB])comp$", accession)
    component_stripped = re.sub(r"_(?:A|B)comp$", "", accession)
    component_specific = (
        f"{component_stripped}_design_{component_match.group(1)}comp.pdb"
        if component_match
        else ""
    )
    names = [
        str(row.get("source_name") or ""),
        f"{accession}_design.pdb",
        f"{accession}_xtal.pdb",
        component_specific,
        f"{component_stripped}_design.pdb",
        f"{component_stripped}_design_Acomp.pdb",
        f"{component_stripped}_design_Bcomp.pdb",
    ]
    for name in names:
        candidate = structure_dir / name
        if name and candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"No local structure for {row['module_id']}; tried {', '.join(names)}"
    )


def _split_candidate_blocks(frame: pd.DataFrame) -> list[pd.DataFrame]:
    """Split historical candidate rows into independent scan results.

    Some early curation artifacts contain more than one whole-protein scan for
    the same module ID (for example, before and after a construct-sequence
    correction).  Each scan starts again at ``candidate_rank == 1``.  Mixing
    those blocks can attach coordinates from an obsolete sequence to the
    current catalog row.
    """
    if frame.empty:
        return []
    block_ids = frame["candidate_rank"].astype(int).eq(1).cumsum()
    return [block.copy() for _, block in frame.groupby(block_ids, sort=False)]


def _candidate_records_for_row(
    row: dict[str, object], candidate_frame: pd.DataFrame
) -> list[dict[str, object]]:
    """Select the scan block whose top call matches the current catalog row."""
    module_frame = candidate_frame.loc[
        candidate_frame["module_id"].astype(str).eq(str(row["module_id"]))
    ]
    blocks = _split_candidate_blocks(module_frame)
    if not blocks:
        return []

    origin = int(row.get("full_sequence_origin") or 1)
    full_length = len(str(row["full_sequence"]))
    expected_period = int(row.get("period") or row.get("unit_length") or 0)
    expected_start = int(row.get("repeat_region_start") or row.get("unit_start")) - origin
    expected_end = int(row.get("repeat_region_end") or row.get("unit_end")) - origin + 1

    ranked: list[tuple[tuple[int, int, int, int], pd.DataFrame]] = []
    for block in blocks:
        first = block.iloc[0]
        all_in_bounds = bool(
            block["local_start"].astype(int).ge(0).all()
            and block["local_end"].astype(int).le(full_length).all()
        )
        exact_top_call = bool(
            int(first["period"]) == expected_period
            and int(first["local_start"]) == expected_start
            and int(first["local_end"]) == expected_end
        )
        top_period_match = int(int(first["period"]) == expected_period)
        coordinate_distance = abs(int(first["local_start"]) - expected_start) + abs(
            int(first["local_end"]) - expected_end
        )
        # Maximize the first three components and then minimize distance.
        rank = (
            int(exact_top_call),
            int(all_in_bounds),
            top_period_match,
            -coordinate_distance,
        )
        ranked.append((rank, block))
    rank, selected = max(ranked, key=lambda value: value[0])
    if not rank[1]:
        raise ValueError(
            f"No candidate scan for {row['module_id']} fits its {full_length}-residue sequence"
        )
    return selected.sort_values("candidate_rank").to_dict(orient="records")


def _annotate_one(
    payload: tuple[
        dict[str, object],
        list[dict[str, object]],
        dict[str, dict[str, object]],
        str,
        str,
        str,
        int,
    ]
) -> dict[str, object]:
    (
        row,
        candidate_records,
        dhr_templates,
        structure_dir_string,
        cache_dir_string,
        dssp_executable,
        timeout,
    ) = payload
    module_id = str(row["module_id"])
    try:
        family = str(row.get("family") or "")
        accession = str(row.get("source_accession") or "")
        full_sequence = str(row["full_sequence"])
        if family.upper() == "DHR" and accession in dhr_templates:
            template = dhr_templates[accession]
            ss3, residues = map_author_template_to_full_sequence(
                full_sequence,
                repeat_start=int(row["prior_unit_start"]),
                full_sequence_origin=int(row.get("full_sequence_origin") or 1),
                ss3_template=str(template["author_repeat_ss3"]),
            )
            metadata = {
                "secondary_structure_method": (
                    "author residue template (helix1-loop1-helix2-loop2)"
                ),
                "secondary_structure_method_version": SECONDARY_STRUCTURE_METHOD_VERSION,
                "secondary_structure_source_type": "author_supplementary_table",
                "secondary_structure_source": (
                    "Brunette et al. Nature 2015 Supplementary Table 2"
                ),
                "secondary_structure_known_fraction": len(residues) / len(full_sequence),
                "structure_sequence_identity": 1.0,
                **template,
            }
        else:
            if str(row.get("collection")) == "natural100":
                structure_path = _fetch_rcsb_structure(
                    accession,
                    Path(cache_dir_string),
                    timeout,
                    cache_key=module_id,
                )
                source_type = "RCSB_experimental_structure_DSSP"
            else:
                structure_path = _resolve_local_structure(
                    row, Path(structure_dir_string)
                )
                source_type = "author_structure_model_DSSP"
            chain_value = row.get("source_chain")
            chain_id = (
                None
                if chain_value is None
                or pd.isna(chain_value)
                or not str(chain_value).strip()
                else str(chain_value).strip()
            )
            ss3, residues, metadata = dssp_annotations_for_structure(
                structure_path,
                full_sequence,
                chain_id=chain_id,
                dssp_executable=dssp_executable,
            )
            metadata["secondary_structure_source_type"] = source_type
            metadata["secondary_structure_source_sha256"] = hashlib.sha256(
                structure_path.read_bytes()
            ).hexdigest()

        candidates = [
            RepeatCandidate(
                **{
                    field: record[field]
                    for field in RepeatCandidate.__dataclass_fields__
                }
            )
            for record in candidate_records
        ]
        supports = score_candidates_with_secondary_structure(candidates, ss3)
        selected_candidate_rows = [
            {
                **candidate_record,
                "candidate_scan_full_sequence_length": len(full_sequence),
                "candidate_scan_full_sequence_sha256": hashlib.sha256(
                    full_sequence.encode()
                ).hexdigest(),
            }
            for candidate_record in candidate_records
        ]
        support_rows = [
            {
                "module_id": module_id,
                "candidate_rank": candidate_record.get("candidate_rank"),
                "sequence_score": candidate_record.get("score"),
                "sequence_adjacent_identity": candidate_record.get("adjacent_identity"),
                "sequence_positive_fraction": candidate_record.get(
                    "adjacent_positive_fraction"
                ),
                "sequence_spectral_concentration": candidate_record.get(
                    "spectral_concentration"
                ),
                **support.to_dict(),
            }
            for candidate_record, support in zip(candidate_records, supports, strict=True)
        ]
        residue_rows = residues.to_dict(orient="records")
        for residue in residue_rows:
            residue["module_id"] = module_id
            residue["secondary_structure_source_type"] = metadata[
                "secondary_structure_source_type"
            ]
        return {
            "module_id": module_id,
            "status": "passed",
            "secondary_structure": ss3,
            "secondary_structure_json": secondary_structure_json(ss3),
            "metadata": metadata,
            "selected_candidate_rows": selected_candidate_rows,
            "support_rows": support_rows,
            "residue_rows": residue_rows,
            "error": None,
        }
    except Exception as error:  # preserved in the audit instead of aborting a shard
        return {
            "module_id": module_id,
            "status": "secondary_structure_unavailable",
            "secondary_structure": None,
            "secondary_structure_json": None,
            "metadata": {},
            "selected_candidate_rows": [],
            "support_rows": [],
            "residue_rows": [],
            "error": f"{type(error).__name__}: {error}",
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--structure-dir", type=Path, required=True)
    parser.add_argument("--structure-cache", type=Path, required=True)
    parser.add_argument("--dhr-supplement", type=Path, required=True)
    parser.add_argument("--dssp-executable", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--allow-missing", action="store_true")
    args = parser.parse_args()

    source = pd.read_parquet(args.catalog)
    candidates = pd.read_parquet(args.candidates)
    dhr_frame = parse_dhr_secondary_structure_table(
        args.dhr_supplement.read_text(errors="replace")
    )
    dhr_frame["secondary_structure_source_sha256"] = hashlib.sha256(
        args.dhr_supplement.read_bytes()
    ).hexdigest()
    dhr_frame["secondary_structure_source_path"] = str(args.dhr_supplement)
    dhr_templates = {
        str(row["source_accession"]): row
        for row in dhr_frame.to_dict(orient="records")
    }
    payloads = [
        (
            row,
            _candidate_records_for_row(row, candidates),
            dhr_templates,
            str(args.structure_dir),
            str(args.structure_cache),
            str(args.dssp_executable),
            args.timeout,
        )
        for row in source.to_dict(orient="records")
    ]
    if args.workers == 1:
        results = [_annotate_one(payload) for payload in payloads]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            results = list(executor.map(_annotate_one, payloads))

    audit_rows: list[dict[str, object]] = []
    support_rows: list[dict[str, object]] = []
    selected_candidate_rows: list[dict[str, object]] = []
    residue_rows: list[dict[str, object]] = []
    metadata_by_module: dict[str, dict[str, object]] = {}
    for result in results:
        module_id = str(result["module_id"])
        metadata = dict(result["metadata"])
        metadata_by_module[module_id] = {
            **metadata,
            "full_secondary_structure": result["secondary_structure"],
            "full_secondary_structure_json": result["secondary_structure_json"],
            "secondary_structure_status": result["status"],
            "secondary_structure_error": result["error"],
        }
        audit_rows.append(
            {
                "module_id": module_id,
                "status": result["status"],
                "error": result["error"],
                "full_sequence_length": len(
                    str(source.loc[source.module_id.eq(module_id), "full_sequence"].iloc[0])
                ),
                "annotated_residue_count": len(result["residue_rows"]),
                "candidate_count": len(result["support_rows"]),
                **metadata,
            }
        )
        support_rows.extend(result["support_rows"])
        selected_candidate_rows.extend(result["selected_candidate_rows"])
        residue_rows.extend(result["residue_rows"])

    source_with_structure = source.copy()
    for index, row in source_with_structure.iterrows():
        for key, value in metadata_by_module[str(row.module_id)].items():
            source_with_structure.at[index, key] = value
    audit = pd.DataFrame(audit_rows)
    supports = pd.DataFrame(support_rows)
    selected_candidates = pd.DataFrame(selected_candidate_rows)
    residues = pd.DataFrame(residue_rows)
    prefix = args.output.with_suffix("")
    _write_frame(audit, prefix.with_name(prefix.name + "_secondary_structure_audit.parquet"))
    _write_frame(
        supports,
        prefix.with_name(prefix.name + "_secondary_structure_candidates.parquet"),
    )
    _write_frame(
        selected_candidates,
        prefix.with_name(prefix.name + "_sequence_candidates_selected.parquet"),
    )
    _write_frame(
        residues,
        prefix.with_name(prefix.name + "_secondary_structure_residues.parquet"),
    )
    _write_frame(
        dhr_frame,
        prefix.with_name(prefix.name + "_dhr_author_secondary_structure.parquet"),
    )
    natural_source_annotation_policy = bool(
        source.collection.astype(str).eq("natural100").all()
    )
    # RepeatsDB unit boundaries are themselves structure-derived and use PDB
    # author coordinates, while sequence candidates use entity-poly indices.
    # Keep those authoritative natural boundaries in the final catalog and
    # expose shorter joint candidates for manual review until a 3D
    # superposition confirms that the apparent harmonic is a complete unit.
    selection_supports = supports.iloc[0:0].copy() if natural_source_annotation_policy else supports
    refined = reselect_module_boundaries(
        source_with_structure,
        candidates,
        secondary_structure_support=selection_supports,
        output_path=prefix.with_name(prefix.name + "_boundary_audit.parquet"),
        unit_alignment_path=prefix.with_name(prefix.name + "_unit_alignment.parquet"),
        position_variability_path=prefix.with_name(
            prefix.name + "_position_variability.parquet"
        ),
    )
    catalog = merge_module_catalogs([refined], args.output)
    passed_count = int(audit.status.eq("passed").sum())
    validation = {
        "boundary_method_version": BOUNDARY_METHOD_VERSION,
        "secondary_structure_method_version": SECONDARY_STRUCTURE_METHOD_VERSION,
        "source_rows": len(source),
        "catalog_rows": len(catalog),
        "secondary_structure_passed": passed_count,
        "secondary_structure_unavailable": int(len(audit) - passed_count),
        "jointly_selected": int(
            refined.secondary_structure_selected_support.fillna(False).sum()
        ),
        "natural_source_annotation_policy": natural_source_annotation_policy,
        "source_annotation_fallback": int(
            refined.boundary_refinement_status.eq("source_prior_fallback").sum()
        ),
        "all_boundaries_resolved": bool(
            refined.boundary_refinement_status.isin(
                ["refined", "source_prior_fallback"]
            ).all()
        ),
        "module_selection_policy": MODULE_SELECTION_POLICY,
        "all_selected_modules_are_middle_policy": bool(
            refined.selected_module_policy.eq(MODULE_SELECTION_POLICY).all()
        ),
        "all_unit_sequences_are_selected_module": bool(
            refined.unit_sequence.eq(refined.selected_module_sequence).all()
            and refined.unit_start.eq(refined.selected_module_start).all()
            and refined.unit_end.eq(refined.selected_module_end).all()
        ),
        "output": str(args.output),
    }
    is_natural = bool(source.collection.astype(str).eq("natural100").all())
    catalog_count_valid = len(catalog) == 100 if is_natural else len(catalog) >= 100
    validation["passed"] = bool(
        catalog_count_valid
        and validation["catalog_rows"] <= validation["source_rows"]
        and validation["all_boundaries_resolved"]
        and validation["all_selected_modules_are_middle_policy"]
        and validation["all_unit_sequences_are_selected_module"]
        and (args.allow_missing or validation["secondary_structure_unavailable"] == 0)
    )
    validation_path = prefix.with_name(prefix.name + "_validation.json")
    validation_path.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0 if validation["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
