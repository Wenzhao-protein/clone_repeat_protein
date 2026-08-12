"""Build the public, spreadsheet-friendly repeat-module result catalog.

The exporter intentionally projects the normalized scientific artifacts into
one row per unique middle module.  It exposes only validated maximum-copy DNA
constructs and never copies raw IDT responses, adaptive traces, credentials, or
machine-specific paths into the public CSV.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pandas.api.types import is_object_dtype, is_string_dtype

from .constants import RULE_PROFILE_NAME, validate_protein_sequence
from .idt import IDT_SCORE_POLICY
from .module_experiments import CORPUS_VERSION
from .optimization import translate_dna

PUBLIC_RESULT_FILENAME = "natural_designed_repeat_protein_hurdler_idt.csv"
# GitHub rejects regular Git blobs above 100,000,000 bytes (decimal MB).
GITHUB_BLOB_LIMIT_BYTES = 100_000_000

_CATALOG_COLUMNS = (
    "module_id",
    "collection",
    "module_type",
    "family",
    "unit_sequence",
    "unit_length",
    "full_sequence",
    "full_sequence_length",
    "source_name",
    "source_accession",
    "source_url",
    "source_chain",
    "source_annotation_id",
    "annotation_uuid",
    "structure_accession",
    "structure_chain",
    "structure_source",
    "uniprot_accessions",
    "protein_key",
    "evidence_tier",
    "reviewed",
    "citation",
    "license_name",
    "retrieved_date",
    "download_date",
    "source_sha256",
    "full_sequence_sha256",
    "sequence_sha256",
    "unit_start",
    "unit_end",
    "selected_module_index",
    "selected_module_count",
    "repeat_region_start",
    "repeat_region_end",
    "repeat_count",
    "boundary_method",
    "boundary_method_version",
    "boundary_refinement_status",
    "module_selection_policy",
    "sequence_coordinate_method",
    "strict_dual_evidence_passed",
    "dssp_state_agreement",
    "dssp_transition_agreement",
    "dssp_median_transitions_per_unit",
    "foldseek_3di_identity",
    "foldseek_median_min_tm",
    "foldseek_median_lddt",
    "foldseek_median_coverage",
    "structure_source_type",
    "corpus_version",
)

_COMPATIBILITY_COLUMNS = (
    "module_id",
    "collection",
    "unit_sequence",
    "hurdler_compatible",
    "compatible_plasmid_count",
    "compatible_plasmids_json",
    "candidate_solution_count",
    "selected_solution_json",
    "selected_plasmid",
    "selected_site_i_enzyme",
    "selected_site_ii_enzyme",
    "selected_site_i_recognition_site",
    "selected_site_ii_recognition_site",
    "selected_direction",
    "selected_site_i_position",
    "selected_site_ii_position",
    "selected_candidate_pair_id",
    "compatibility_method_version",
    "rules_version",
)

_MAXIMUM_COLUMNS = (
    "module_id",
    "collection",
    "unit_sequence",
    "fragment_limit_bp",
    "mathematical_max_copies",
    "verified_max_copies",
    "adaptive_maximum_proof_status",
    "adaptive_stop_reason",
    "optimization_status",
    "failure_reason",
    "final_ga_weights_json",
    "idt_status",
    "idt_complexity_score",
    "idt_score_policy",
    "idt_positive_score_names_json",
    "idt_rule_details_json",
    "idt_response_sha256",
    "idt_scored_sequence_sha256",
    "selected_pair_re_site_excess",
    "dna_sequence",
    "dna_length",
    "validation_passed",
    "validation_reasons_json",
)

_SOURCE_MAPPING_COLUMNS = (
    "module_id",
    "source_name",
    "source_accession",
    "source_url",
    "structure_accession",
    "uniprot_accessions",
    "uniprot_accessions_json",
    "annotation_uuid",
    "source_annotation_id",
    "citation",
)


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_projected(
    path: str | Path,
    requested: Sequence[str],
    *,
    required: Iterable[str],
) -> pd.DataFrame:
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(source)
    if source.suffix == ".parquet":
        available = set(pq.read_schema(source).names)
        columns = [column for column in requested if column in available]
        frame = pd.read_parquet(source, columns=columns)
    elif source.suffix == ".csv":
        available = set(pd.read_csv(source, nrows=0).columns)
        columns = [column for column in requested if column in available]
        frame = pd.read_csv(source, usecols=columns, low_memory=False)
    else:
        raise ValueError(f"Expected Parquet or CSV input: {source}")
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{source} is missing required columns: {missing}")
    return frame


def _collection_label(value: Any, module_type: Any = "") -> str:
    text = str(value or "").strip().lower()
    kind = str(module_type or "").strip().lower()
    if kind == "natural" or text.startswith("natural"):
        return "Natural"
    if kind == "designed" or text.startswith("designed"):
        return "Designed"
    raise ValueError(f"Unknown module collection: {value!r}")


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def _truthy(value: Any) -> bool:
    if _is_missing(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "passed"}
    return bool(value)


def _compact_text(value: Any) -> str:
    if _is_missing(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def _spreadsheet_safe(value: Any) -> str:
    text = _compact_text(value)
    if text and text[0] in "=+-@":
        return "'" + text
    return text


def _decoded_values(value: Any) -> list[str]:
    if _is_missing(value):
        return []
    if isinstance(value, (list, tuple, set, np.ndarray)):
        output: list[str] = []
        for item in value:
            output.extend(_decoded_values(item))
        return output
    if isinstance(value, dict):
        output = []
        for item in value.values():
            output.extend(_decoded_values(item))
        return output
    text = _compact_text(value)
    if not text:
        return []
    if text[:1] in "[{(" and text[-1:] in "]})":
        for decoder in (json.loads, ast.literal_eval):
            try:
                decoded = decoder(text)
            except (ValueError, SyntaxError, json.JSONDecodeError):
                continue
            return _decoded_values(decoded)
    return [text]


def _joined_unique(values: Iterable[Any]) -> str:
    unique = {
        _spreadsheet_safe(item)
        for value in values
        for item in _decoded_values(value)
        if _compact_text(item)
    }
    return "|".join(sorted(unique, key=lambda item: item.casefold()))


def _aggregate_source_mappings(paths: Sequence[str | Path]) -> pd.DataFrame:
    frames = [
        _read_projected(path, _SOURCE_MAPPING_COLUMNS, required={"module_id"})
        for path in paths
    ]
    if not frames:
        return pd.DataFrame(columns=["module_id", "source_mapping_count"])
    mappings = pd.concat(frames, ignore_index=True, sort=False)
    output: list[dict[str, Any]] = []
    aggregate_fields = {
        "source_name": "all_source_names",
        "source_accession": "all_source_accessions",
        "source_url": "all_source_urls",
        "structure_accession": "all_structure_accessions",
        "uniprot_accessions": "all_uniprot_accessions",
        "uniprot_accessions_json": "all_uniprot_accessions_json",
        "annotation_uuid": "all_annotation_ids",
        "source_annotation_id": "all_source_annotation_ids",
        "citation": "all_citations",
    }
    for module_id, group in mappings.groupby("module_id", sort=False):
        row: dict[str, Any] = {
            "module_id": str(module_id),
            "source_mapping_count": int(len(group)),
        }
        for source, destination in aggregate_fields.items():
            row[destination] = (
                _joined_unique(group[source]) if source in group else ""
            )
        output.append(row)
    aggregated = pd.DataFrame(output)
    if "all_uniprot_accessions_json" in aggregated:
        aggregated["all_uniprot_accessions"] = [
            _joined_unique([left, right])
            for left, right in zip(
                aggregated["all_uniprot_accessions"],
                aggregated.pop("all_uniprot_accessions_json"),
                strict=True,
            )
        ]
    return aggregated


def _parse_json_object(value: Any) -> dict[str, Any]:
    if _is_missing(value) or not str(value).strip():
        return {}
    try:
        decoded = json.loads(str(value))
    except json.JSONDecodeError:
        return {}
    return decoded if isinstance(decoded, dict) else {}


def _rule_reason_summary(value: Any) -> str:
    if _is_missing(value) or not str(value).strip():
        return ""
    try:
        rules = json.loads(str(value))
    except json.JSONDecodeError:
        return "invalid_rule_details"
    if not isinstance(rules, list):
        return "invalid_rule_details"
    summaries: list[str] = []
    for rule in rules:
        if not isinstance(rule, dict):
            continue
        score = rule.get("score")
        try:
            numeric_score = float(score)
        except (TypeError, ValueError):
            numeric_score = math.nan
        if not (rule.get("is_violated") or (math.isfinite(numeric_score) and numeric_score > 0)):
            continue
        name = _compact_text(rule.get("name")) or "unnamed_rule"
        fields = [name]
        if math.isfinite(numeric_score):
            fields.append(f"score={numeric_score:g}")
        if not _is_missing(rule.get("actual_value")):
            fields.append(f"actual={_compact_text(rule.get('actual_value'))}")
        if not _is_missing(rule.get("threshold_value")):
            fields.append(f"threshold={_compact_text(rule.get('threshold_value'))}")
        summaries.append(":".join(fields))
    return "|".join(summaries)


def _compact_weight_mapping(value: Any) -> str:
    mapping = _parse_json_object(value)
    if not mapping:
        return ""
    return ";".join(
        f"{_compact_text(key)}={_compact_text(mapping[key])}"
        for key in sorted(mapping)
    )


def _repo_relative(path: str | Path, repository_root: str | Path) -> str:
    resolved = Path(path).resolve()
    root = Path(repository_root).resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError as exc:
        raise ValueError(f"Public artifact is outside the repository: {resolved}") from exc


def _validate_accepted_construct(row: pd.Series, capacity: int) -> tuple[str, int, str]:
    dna = re.sub(r"\s+", "", str(row.get("dna_sequence") or "")).upper()
    copies = int(row["verified_max_copies"])
    unit = validate_protein_sequence(str(row["unit_sequence"]))
    if copies < 2:
        raise ValueError("Validated repeat maximum must contain at least two copies")
    if not dna or len(dna) > capacity or len(dna) != int(row["dna_length"]):
        raise ValueError("Validated construct has an invalid DNA length")
    if translate_dna(dna) != unit * copies:
        raise ValueError("Validated construct does not translate to unit_sequence × copies")
    if int(row["selected_pair_re_site_excess"]) != 0:
        raise ValueError("Validated construct has selected-pair excess RE sites")
    score = float(row["idt_complexity_score"])
    if not math.isfinite(score) or score >= 10:
        raise ValueError("Validated construct does not have an IDT score below 10")
    if str(row["idt_score_policy"]) != IDT_SCORE_POLICY:
        raise ValueError("Validated construct uses the wrong IDT score policy")
    dna_sha = hashlib.sha256(dna.encode()).hexdigest()
    if str(row["idt_scored_sequence_sha256"]) != dna_sha:
        raise ValueError("Validated construct DNA does not match its IDT request hash")
    response_sha = str(row["idt_response_sha256"])
    if not re.fullmatch(r"[0-9a-f]{64}", response_sha):
        raise ValueError("Validated construct is missing a valid IDT response hash")
    if row.get("adaptive_maximum_proof_status") not in {
        "capacity_limit_reached",
        "next_copy_failed_at_100",
    }:
        raise ValueError("Validated construct is missing a maximum-copy proof")
    return dna, copies, dna_sha


def _capacity_frame(maximum: pd.DataFrame, capacity: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    subset = maximum.loc[maximum.fragment_limit_bp.astype(int).eq(capacity)]
    for source in subset.to_dict(orient="records"):
        row = pd.Series(source)
        accepted = _truthy(row.get("validation_passed"))
        dna = ""
        copies: int | str = ""
        dna_sha = ""
        if accepted:
            dna, copies, dna_sha = _validate_accepted_construct(row, capacity)
        status = "idt_accepted" if accepted else _compact_text(
            row.get("optimization_status")
        ) or "no_accepted_repeat_construct"
        positive_names = _joined_unique(
            [row.get("idt_positive_score_names_json")]
        )
        rows.append(
            {
                "collection": row["collection"],
                "module_id": row["module_id"],
                f"cap{capacity}_result_status": status,
                f"cap{capacity}_mathematical_max_copies": row.get(
                    "mathematical_max_copies"
                ),
                f"cap{capacity}_maximum_verified_copies": copies,
                f"cap{capacity}_maximum_proof": _compact_text(
                    row.get("adaptive_maximum_proof_status")
                ),
                f"cap{capacity}_stop_reason": _spreadsheet_safe(
                    row.get("adaptive_stop_reason")
                ),
                f"cap{capacity}_optimization_status": _compact_text(
                    row.get("optimization_status")
                ),
                f"cap{capacity}_failure_reason": _spreadsheet_safe(
                    row.get("failure_reason")
                ),
                f"cap{capacity}_final_ga_weights": _compact_weight_mapping(
                    row.get("final_ga_weights_json")
                ),
                f"cap{capacity}_idt_passed": accepted,
                f"cap{capacity}_idt_status": _compact_text(row.get("idt_status")),
                f"cap{capacity}_idt_score_sum": (
                    row.get("idt_complexity_score") if accepted else ""
                ),
                f"cap{capacity}_idt_positive_rules": positive_names,
                f"cap{capacity}_idt_rule_reasons": _rule_reason_summary(
                    row.get("idt_rule_details_json")
                ),
                f"cap{capacity}_idt_response_sha256": (
                    _compact_text(row.get("idt_response_sha256")) if accepted else ""
                ),
                f"cap{capacity}_selected_pair_excess_sites": row.get(
                    "selected_pair_re_site_excess"
                ),
                f"cap{capacity}_validation_reasons": _compact_text(
                    row.get("validation_reasons_json")
                ),
                f"cap{capacity}_idt_accepted_dna": dna,
                f"cap{capacity}_idt_accepted_dna_length_bp": len(dna) if dna else "",
                f"cap{capacity}_idt_accepted_dna_sha256": dna_sha,
            }
        )
    return pd.DataFrame(rows)


def export_module_results(
    catalog_path: str | Path,
    source_mapping_paths: Sequence[str | Path],
    compatibility_path: str | Path,
    maximum_results_path: str | Path,
    output_path: str | Path,
    *,
    repository_root: str | Path,
    generated_at_utc: str | None = None,
) -> pd.DataFrame:
    """Export one manually searchable row per active middle-repeat module."""
    catalog = _read_projected(
        catalog_path,
        _CATALOG_COLUMNS,
        required={"module_id", "collection", "unit_sequence", "unit_length", "full_sequence"},
    )
    catalog["collection"] = [
        _collection_label(collection, module_type)
        for collection, module_type in zip(
            catalog.collection,
            catalog.get("module_type", pd.Series([""] * len(catalog))),
            strict=True,
        )
    ]
    catalog["unit_sequence"] = catalog.unit_sequence.map(validate_protein_sequence)
    catalog["full_sequence"] = catalog.full_sequence.map(validate_protein_sequence)
    if not catalog.unit_length.astype(int).eq(catalog.unit_sequence.str.len()).all():
        raise ValueError("Catalog unit_length does not match unit_sequence")
    if catalog.duplicated(["collection", "module_id"]).any():
        raise ValueError("Catalog has duplicate collection/module_id rows")
    if catalog.duplicated(["collection", "unit_sequence"]).any():
        raise ValueError("Catalog is not deduplicated by exact middle-module sequence")
    if "corpus_version" in catalog and not catalog.corpus_version.fillna(
        CORPUS_VERSION
    ).eq(CORPUS_VERSION).all():
        raise ValueError("Catalog contains a non-active corpus version")

    compatibility = _read_projected(
        compatibility_path,
        _COMPATIBILITY_COLUMNS,
        required={"module_id", "collection", "unit_sequence", "hurdler_compatible"},
    )
    compatibility["collection"] = compatibility.collection.map(_collection_label)
    if compatibility.duplicated(["collection", "module_id"]).any():
        raise ValueError("Compatibility input has duplicate module rows")
    if len(compatibility) != len(catalog):
        raise ValueError(
            f"Catalog/compatibility row mismatch: {len(catalog)} != {len(compatibility)}"
        )

    maximum = _read_projected(
        maximum_results_path,
        _MAXIMUM_COLUMNS,
        required={
            "module_id",
            "collection",
            "unit_sequence",
            "fragment_limit_bp",
            "validation_passed",
        },
    )
    maximum["collection"] = maximum.collection.map(_collection_label)
    if maximum.duplicated(["collection", "module_id", "fragment_limit_bp"]).any():
        raise ValueError("Maximum-copy input has duplicate module/capacity rows")
    compatible_keys = set(
        compatibility.loc[
            compatibility.hurdler_compatible.astype(bool),
            ["collection", "module_id"],
        ].itertuples(index=False, name=None)
    )
    expected = {
        (collection, module_id, capacity)
        for collection, module_id in compatible_keys
        for capacity in (1800, 3000)
    }
    observed = {
        (str(row.collection), str(row.module_id), int(row.fragment_limit_bp))
        for row in maximum.itertuples(index=False)
    }
    if expected != observed:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        raise ValueError(
            "Maximum-copy table is incomplete: "
            f"missing={len(missing)} {missing[:5]}, extra={len(extra)} {extra[:5]}"
        )

    selected = compatibility.get(
        "selected_solution_json", pd.Series([""] * len(compatibility))
    ).map(_parse_json_object)
    compatibility["selected_site_iii_enzymes"] = selected.map(
        lambda item: _compact_text(item.get("site_iii_enzymes"))
    )
    compatibility["selected_site_iii_recognition_sites"] = selected.map(
        lambda item: _compact_text(item.get("site_iii_sites"))
    )
    compatibility["selected_site_i_overhang"] = selected.map(
        lambda item: item.get("site_i_ovhg", "")
    )
    compatibility["selected_site_ii_overhang"] = selected.map(
        lambda item: item.get("site_ii_ovhg", "")
    )
    compatibility["compatible_plasmids"] = compatibility.get(
        "compatible_plasmids_json", pd.Series([""] * len(compatibility))
    ).map(lambda value: _joined_unique([value]))
    compatibility["hurdler_status"] = compatibility.hurdler_compatible.map(
        lambda value: "compatible" if _truthy(value) else "incompatible"
    )

    provenance = _aggregate_source_mappings(source_mapping_paths)
    output = catalog.merge(
        provenance, on="module_id", how="left", validate="one_to_one"
    )
    compatibility_drop = ["unit_sequence", "selected_solution_json", "compatible_plasmids_json"]
    output = output.merge(
        compatibility.drop(columns=compatibility_drop, errors="ignore"),
        on=["collection", "module_id"],
        how="left",
        validate="one_to_one",
    )
    for capacity in (1800, 3000):
        output = output.merge(
            _capacity_frame(maximum, capacity),
            on=["collection", "module_id"],
            how="left",
            validate="one_to_one",
        )
        status_column = f"cap{capacity}_result_status"
        pass_column = f"cap{capacity}_idt_passed"
        output[status_column] = output[status_column].fillna(
            "not_applicable_hurdler_incompatible"
        )
        output[pass_column] = output[pass_column].fillna(False).astype(bool)
        for column in output.columns:
            if column.startswith(f"cap{capacity}_") and (
                is_object_dtype(output[column].dtype)
                or is_string_dtype(output[column].dtype)
            ):
                output[column] = output[column].fillna("")

    output["display_name"] = [
        _spreadsheet_safe(accession) or _spreadsheet_safe(structure) or str(module_id)
        for accession, structure, module_id in zip(
            output.get("source_accession", pd.Series([""] * len(output))),
            output.get("structure_accession", pd.Series([""] * len(output))),
            output.module_id,
            strict=True,
        )
    ]
    output["full_protein_length_aa"] = output.full_sequence.str.len()
    output["middle_module_number_one_based"] = pd.to_numeric(
        output.get("selected_module_index"), errors="coerce"
    ).add(1).astype("Int64")
    output["coordinate_system"] = "1-based inclusive"
    output["record_status"] = [
        (
            "hurdler_incompatible"
            if not _truthy(compatible)
            else "idt_accepted_both_capacities"
            if _truthy(cap1800) and _truthy(cap3000)
            else "idt_accepted_partial"
            if _truthy(cap1800) or _truthy(cap3000)
            else "no_accepted_repeat_construct"
        )
        for compatible, cap1800, cap3000 in zip(
            output.hurdler_compatible,
            output.cap1800_idt_passed,
            output.cap3000_idt_passed,
            strict=True,
        )
    ]
    output["search_terms"] = [
        _joined_unique(values)
        for values in zip(
            output.display_name,
            output.module_id,
            output.collection,
            output.family,
            output.get("source_name", pd.Series([""] * len(output))),
            output.get("source_accession", pd.Series([""] * len(output))),
            output.get("all_source_accessions", pd.Series([""] * len(output))),
            output.get("all_uniprot_accessions", pd.Series([""] * len(output))),
            output.get("all_structure_accessions", pd.Series([""] * len(output))),
            output.get("all_annotation_ids", pd.Series([""] * len(output))),
            output.get("selected_plasmid", pd.Series([""] * len(output))),
            output.get("selected_site_i_enzyme", pd.Series([""] * len(output))),
            output.get("selected_site_ii_enzyme", pd.Series([""] * len(output))),
            strict=True,
        )
    ]

    generated = generated_at_utc or datetime.now(timezone.utc).replace(
        microsecond=0
    ).isoformat()
    input_paths = [
        Path(catalog_path),
        *map(Path, source_mapping_paths),
        Path(compatibility_path),
        Path(maximum_results_path),
    ]
    bundle_hash = hashlib.sha256(
        "\n".join(f"{_repo_relative(path, repository_root)}:{_sha256(path)}" for path in input_paths).encode()
    ).hexdigest()
    output["corpus_version"] = CORPUS_VERSION
    output["hurdler_rules_version"] = RULE_PROFILE_NAME
    output["idt_policy"] = IDT_SCORE_POLICY
    output["generated_at_utc"] = generated
    output["input_bundle_sha256"] = bundle_hash
    rename = {
        "unit_sequence": "middle_module_sequence_aa",
        "unit_length": "middle_module_length_aa",
        "full_sequence": "full_protein_sequence_aa",
        "selected_module_index": "middle_module_index_zero_based",
        "unit_start": "middle_module_start",
        "unit_end": "middle_module_end",
        "boundary_refinement_status": "boundary_status",
    }
    output = output.rename(columns=rename)
    output = output.drop(
        columns=[
            "module_type",
            "full_sequence_length",
            "source_annotation_id",
            "annotation_uuid",
            "uniprot_accessions",
            "download_date",
            "source_sha256",
            "full_sequence_sha256",
            "sequence_sha256",
            "middle_module_index_zero_based",
            "selected_module_count",
            "boundary_method_version",
            "module_selection_policy",
            "all_source_names",
            "all_source_accessions",
            "all_source_urls",
            "all_structure_accessions",
            "all_annotation_ids",
            "all_source_annotation_ids",
            "all_citations",
            "selected_candidate_pair_id",
            "compatibility_method_version",
            "rules_version",
            "hurdler_status",
            "cap1800_idt_positive_rules",
            "cap3000_idt_positive_rules",
            "cap1800_validation_reasons",
            "cap3000_validation_reasons",
            "cap1800_optimization_status",
            "cap3000_optimization_status",
            "coordinate_system",
        ],
        errors="ignore",
    )
    leading = [
        "display_name",
        "search_terms",
        "collection",
        "module_id",
        "family",
        "record_status",
        "middle_module_sequence_aa",
        "middle_module_length_aa",
        "full_protein_sequence_aa",
        "full_protein_length_aa",
        "source_accession",
        "source_url",
        "all_uniprot_accessions",
        "all_structure_accessions",
        "all_annotation_ids",
        "source_mapping_count",
        "hurdler_compatible",
        "selected_plasmid",
        "selected_site_i_enzyme",
        "selected_site_ii_enzyme",
        "cap1800_maximum_verified_copies",
        "cap1800_idt_passed",
        "cap1800_idt_accepted_dna",
        "cap3000_maximum_verified_copies",
        "cap3000_idt_passed",
        "cap3000_idt_accepted_dna",
    ]
    ordered = [column for column in leading if column in output]
    ordered.extend(column for column in output.columns if column not in ordered)
    output = output[ordered].sort_values(
        ["collection", "family", "middle_module_length_aa", "display_name", "module_id"],
        kind="mergesort",
    ).reset_index(drop=True)

    text_fields = [
        "display_name",
        "source_name",
        "source_accession",
        "source_url",
        "citation",
        "license_name",
        "all_source_names",
        "all_source_accessions",
        "all_source_urls",
        "all_citations",
    ]
    for column in text_fields:
        if column in output:
            output[column] = output[column].map(_spreadsheet_safe)

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(
        destination,
        index=False,
        encoding="utf-8",
        quoting=csv.QUOTE_MINIMAL,
        lineterminator="\n",
    )
    if destination.stat().st_size >= GITHUB_BLOB_LIMIT_BYTES:
        destination.unlink()
        raise ValueError(
            "Public result CSV exceeds GitHub's 100-MiB blob limit; "
            "remove nonessential projected columns before publishing"
        )
    return output
