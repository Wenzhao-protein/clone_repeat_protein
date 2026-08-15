"""Audit whether every component in a selected complete route can be ordered.

This is deliberately downstream of the molecular route search.  It does not
invent a new route or change an exact target: it normalizes the selected seed
and donor purchases into explicit complementary oligo pairs or gBlocks and
then applies the product-specific acceptance rules.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .dna_assembly import (
    DEFAULT_GBLOCK_MAX_BP,
    DEFAULT_GBLOCK_MIN_BP,
    DEFAULT_OLIGO_MAX_BP,
    DEFAULT_OLIGO_MIN_BP,
    reverse_complement,
    validate_dna,
)
from .idt import IDT_ORDERABILITY_THRESHOLD, IDT_SCORE_POLICY, IDTComplexityScorer
from .io import sha256_file, utc_now, write_json_atomic


PURCHASE_ORDERABILITY_VERSION = "complete-route-purchase-orderability-v1"
PURCHASE_POLICY = "idt-custom-dna-oligo-pair-or-gblock-v1"
STANDARD_PRIMER_MAX_BP = 89


def _purchase_id(product: str, sequences: Iterable[str]) -> str:
    payload = product + "|" + "|".join(sequences)
    return "purchase_" + hashlib.sha256(payload.encode()).hexdigest()[:20]


def _numeric_score(result: dict[str, Any]) -> float | None:
    value = result.get("idt_complexity_score")
    if isinstance(value, bool):
        return None
    try:
        score = float(value)
    except (TypeError, ValueError):
        return None
    return score if math.isfinite(score) else None


def classify_double_stranded_purchase(
    name: str,
    sequence: str,
    *,
    idt_scorer: IDTComplexityScorer | None,
) -> dict[str, Any]:
    """Normalize one exact dsDNA purchase and apply its product policy."""
    core = validate_dna(sequence)
    length = len(core)
    base: dict[str, Any] = {
        "component_name": str(name),
        "core_sequence": core,
        "core_length_bp": length,
        "idt_policy": "not_applicable_to_oligo_pair",
        "idt_status": "not_required",
        "idt_score": None,
        "idt_response_sha256": "",
        "failure_reason": "",
    }
    # Prefer a live-screened gBlock throughout its supported interval.  The
    # 125--200 overlap could also be purchased as Ultramers, but treating it as
    # a gBlock gives the production claim stronger, product-specific evidence.
    if DEFAULT_OLIGO_MIN_BP <= length < DEFAULT_GBLOCK_MIN_BP:
        reverse = reverse_complement(core)
        product = (
            "complementary_standard_primer_pair"
            if length <= STANDARD_PRIMER_MAX_BP
            else "complementary_ultramer_pair"
        )
        base.update(
            {
                "product_class": product,
                "purchase_sequence_1": core,
                "purchase_sequence_2": reverse,
                "purchase_sequence_count": 2,
                "maximum_single_oligo_length_nt": length,
                "purchase_id": _purchase_id(product, (core, reverse)),
                "orderable": True,
            }
        )
        return base
    if DEFAULT_GBLOCK_MIN_BP <= length <= DEFAULT_GBLOCK_MAX_BP:
        product = "gblock"
        base.update(
            {
                "product_class": product,
                "purchase_sequence_1": core,
                "purchase_sequence_2": "",
                "purchase_sequence_count": 1,
                "maximum_single_oligo_length_nt": None,
                "purchase_id": _purchase_id(product, (core,)),
                "idt_policy": IDT_SCORE_POLICY,
            }
        )
        if idt_scorer is None:
            base.update(
                {
                    "idt_status": "api_unclassified",
                    "orderable": False,
                    "failure_reason": "live_idt_score_required",
                }
            )
            return base
        try:
            result = idt_scorer.score(str(name), core)
        except Exception as exc:  # production output records class, never secrets
            base.update(
                {
                    "idt_status": "api_failure",
                    "orderable": False,
                    "failure_reason": f"idt_api_error:{type(exc).__name__}",
                }
            )
            return base
        score = _numeric_score(result)
        response_sha = str(result.get("idt_response_sha256") or "")
        passed = (
            score is not None
            and score < IDT_ORDERABILITY_THRESHOLD
            and result.get("idt_explicit_pass") is True
            and bool(response_sha)
        )
        base.update(
            {
                "idt_status": str(result.get("idt_status") or "scored_unclassified"),
                "idt_score": score,
                "idt_response_sha256": response_sha,
                "orderable": passed,
                "failure_reason": "" if passed else "gblock_not_idt_accepted",
            }
        )
        return base
    base.update(
        {
            "product_class": "unsupported_length",
            "purchase_sequence_1": core,
            "purchase_sequence_2": "",
            "purchase_sequence_count": 0,
            "maximum_single_oligo_length_nt": None,
            "purchase_id": _purchase_id("unsupported_length", (core,)),
            "orderable": False,
            "failure_reason": "outside_supported_oligo_and_gblock_ranges",
        }
    )
    return base


def classify_existing_primer_pair(
    name: str,
    forward: str,
    reverse: str,
    *,
    core_sequence: str,
) -> dict[str, Any]:
    """Validate the two actual 5'-to-3' oligos emitted by the route planner."""
    first = validate_dna(forward)
    second = validate_dna(reverse)
    core = validate_dna(core_sequence)
    lengths = (len(first), len(second))
    accepted = all(
        DEFAULT_OLIGO_MIN_BP <= length <= DEFAULT_OLIGO_MAX_BP
        for length in lengths
    )
    maximum = max(lengths)
    product = (
        "sticky_end_standard_primer_pair"
        if maximum <= STANDARD_PRIMER_MAX_BP
        else "sticky_end_ultramer_pair"
    )
    return {
        "component_name": str(name),
        "core_sequence": core,
        "core_length_bp": len(core),
        "product_class": product,
        "purchase_sequence_1": first,
        "purchase_sequence_2": second,
        "purchase_sequence_count": 2,
        "maximum_single_oligo_length_nt": maximum,
        "purchase_id": _purchase_id(product, (first, second)),
        "idt_policy": "not_applicable_to_oligo_pair",
        "idt_status": "not_required",
        "idt_score": None,
        "idt_response_sha256": "",
        "orderable": accepted,
        "failure_reason": "" if accepted else "oligo_length_out_of_range",
    }


def _read_shards(raw_root: Path, name: str) -> pd.DataFrame:
    paths = sorted(raw_root.glob(f"shard_*/complete_route_{name}.parquet"))
    frames = [pd.read_parquet(path) for path in paths]
    nonempty = [frame for frame in frames if not frame.empty]
    return pd.concat(nonempty, ignore_index=True, sort=False) if nonempty else pd.DataFrame()


def _validate_shards(raw_root: Path, expected_shards: int | None) -> list[Path]:
    manifests = sorted(raw_root.glob("shard_*/complete_route_manifest.json"))
    if not manifests:
        raise FileNotFoundError(f"No complete-route shards found under {raw_root}")
    indices: list[int] = []
    declared_counts: set[int] = set()
    for path in manifests:
        payload = json.loads(path.read_text())
        indices.append(int(payload["shard_index"]))
        declared_counts.add(int(payload["shard_count"]))
        if payload.get("idt_required") is not True:
            raise ValueError(f"Input shard was not a live-IDT production shard: {path}")
    if len(declared_counts) != 1:
        raise ValueError("Input shard manifests disagree on shard_count")
    declared = declared_counts.pop()
    expected = declared if expected_shards is None else int(expected_shards)
    if declared != expected or set(indices) != set(range(expected)):
        raise ValueError(
            f"Incomplete shard set: declared={declared}, expected={expected}, "
            f"observed={len(indices)}"
        )
    return manifests


def audit_complete_route_purchase_orderability(
    raw_root: str | Path,
    output_dir: str | Path,
    *,
    idt_scorer: IDTComplexityScorer | None,
    expected_shards: int | None = None,
    expected_routes: int | None = None,
    expected_elements: int | None = None,
) -> dict[str, pd.DataFrame]:
    """Audit every selected complete route and write compact production tables."""
    raw = Path(raw_root)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    manifests = _validate_shards(raw, expected_shards)
    targets = _read_shards(raw, "targets")
    seeds = _read_shards(raw, "seeds")
    fragments = _read_shards(raw, "fragments")
    passed = targets.loc[targets.complete_route_verified.fillna(False)].copy()
    if passed.complete_route_id.duplicated().any():
        raise ValueError("Selected complete route IDs are not unique")
    if expected_routes is not None and len(passed) != int(expected_routes):
        raise ValueError(f"Expected {expected_routes} selected routes, observed {len(passed)}")
    element_keys = passed[["source_database", "element_id"]].drop_duplicates()
    if expected_elements is not None and len(element_keys) != int(expected_elements):
        raise ValueError(
            f"Expected {expected_elements} elements with routes, observed {len(element_keys)}"
        )

    seed_by_element = seeds.drop_duplicates(["source_database", "element_id"]).set_index(
        ["source_database", "element_id"]
    )
    component_rows: list[dict[str, Any]] = []
    seed_classification: dict[tuple[str, str], dict[str, Any]] = {}
    for key in map(tuple, element_keys.itertuples(index=False, name=None)):
        if key not in seed_by_element.index:
            raise ValueError(f"Selected element lacks seed evidence: {key}")
        seed = seed_by_element.loc[key]
        classified = classify_double_stranded_purchase(
            f"seed|{key[0]}|{key[1]}",
            str(seed.core_sequence),
            idt_scorer=idt_scorer,
        )
        seed_classification[key] = classified

    passed_ids = set(passed.complete_route_id.astype(str))
    selected_fragments = fragments.loc[
        fragments.complete_route_id.astype(str).isin(passed_ids)
    ].copy()
    fragment_classification: dict[tuple[str, str], dict[str, Any]] = {}
    for fragment in selected_fragments.itertuples(index=False):
        key = (str(fragment.complete_route_id), str(fragment.fragment_id))
        if str(fragment.product_type) == "annealed_sticky_end_primer_pair":
            classified = classify_existing_primer_pair(
                f"donor|{key[0]}|{key[1]}",
                str(fragment.primer_forward_5to3),
                str(fragment.primer_reverse_5to3),
                core_sequence=str(fragment.core_sequence),
            )
        else:
            classified = classify_double_stranded_purchase(
                f"donor|{key[0]}|{key[1]}",
                str(fragment.purchase_sequence),
                idt_scorer=idt_scorer,
            )
        fragment_classification[key] = classified

    fragments_by_route = {
        str(route_id): frame
        for route_id, frame in selected_fragments.groupby("complete_route_id", sort=False)
    }
    route_rows: list[dict[str, Any]] = []
    for target in passed.itertuples(index=False):
        route_id = str(target.complete_route_id)
        element_key = (str(target.source_database), str(target.element_id))
        seed_component = dict(seed_classification[element_key])
        occurrence = {
            "route_id": route_id,
            "target_id": str(target.target_id),
            "source_database": element_key[0],
            "element_id": element_key[1],
            "target_copy_count": int(target.target_copy_count),
            "component_role": "seed",
            "component_index": 0,
            **seed_component,
        }
        component_rows.append(occurrence)
        route_components = [occurrence]
        route_fragments = fragments_by_route.get(route_id)
        if route_fragments is None or route_fragments.empty:
            raise ValueError(f"Selected route has no donor purchase fragments: {route_id}")
        for index, fragment in enumerate(route_fragments.itertuples(index=False), start=1):
            key = (route_id, str(fragment.fragment_id))
            item = {
                "route_id": route_id,
                "target_id": str(target.target_id),
                "source_database": element_key[0],
                "element_id": element_key[1],
                "target_copy_count": int(target.target_copy_count),
                "component_role": "donor",
                "component_index": index,
                **fragment_classification[key],
            }
            component_rows.append(item)
            route_components.append(item)
        failures = sorted(
            {str(row["failure_reason"]) for row in route_components if not row["orderable"]}
        )
        product_counts = pd.Series(
            [row["product_class"] for row in route_components]
        ).value_counts()
        route_rows.append(
            {
                "version": PURCHASE_ORDERABILITY_VERSION,
                "purchase_policy": PURCHASE_POLICY,
                "source_database": element_key[0],
                "element_id": element_key[1],
                "target_id": str(target.target_id),
                "target_copy_count": int(target.target_copy_count),
                "complete_route_id": route_id,
                "component_occurrence_count": len(route_components),
                "unique_purchase_count": len(
                    {str(row["purchase_id"]) for row in route_components}
                ),
                "standard_primer_pair_count": int(
                    sum(
                        int(value)
                        for key, value in product_counts.items()
                        if "standard_primer_pair" in key
                    )
                ),
                "ultramer_pair_count": int(
                    sum(
                        int(value)
                        for key, value in product_counts.items()
                        if "ultramer_pair" in key
                    )
                ),
                "gblock_count": int(product_counts.get("gblock", 0)),
                "every_component_orderable": not failures,
                "failure_reasons": ";".join(failures),
            }
        )

    components = pd.DataFrame(component_rows)
    routes = pd.DataFrame(route_rows).sort_values(
        ["source_database", "element_id", "target_copy_count"]
    )
    if len(routes) != len(passed) or routes.complete_route_id.duplicated().any():
        raise ValueError("Route audit cardinality differs from selected-route input")
    unique_purchase_rows: list[dict[str, Any]] = []
    for purchase_id, frame in components.groupby("purchase_id", sort=True):
        representative = frame.iloc[0]
        varying = [
            column
            for column in (
                "product_class", "purchase_sequence_1", "purchase_sequence_2",
                "orderable", "idt_score", "idt_response_sha256"
            )
            if frame[column].fillna("").astype(str).nunique() != 1
        ]
        if varying:
            raise ValueError(f"Purchase identity has inconsistent fields: {purchase_id}: {varying}")
        unique_purchase_rows.append(
            {
                "version": PURCHASE_ORDERABILITY_VERSION,
                "purchase_policy": PURCHASE_POLICY,
                "purchase_id": purchase_id,
                "product_class": representative.product_class,
                "purchase_sequence_1": representative.purchase_sequence_1,
                "purchase_sequence_2": representative.purchase_sequence_2,
                "purchase_sequence_count": int(representative.purchase_sequence_count),
                "maximum_single_oligo_length_nt": representative.maximum_single_oligo_length_nt,
                "core_length_bp": int(representative.core_length_bp),
                "route_count": int(frame.route_id.nunique()),
                "occurrence_count": int(len(frame)),
                "orderable": bool(representative.orderable),
                "idt_policy": representative.idt_policy,
                "idt_status": representative.idt_status,
                "idt_score": representative.idt_score,
                "idt_response_sha256": representative.idt_response_sha256,
                "failure_reason": representative.failure_reason,
            }
        )
    purchases = pd.DataFrame(unique_purchase_rows).sort_values(
        ["product_class", "purchase_id"]
    )
    element_rows: list[dict[str, Any]] = []
    for (source, element), frame in routes.groupby(
        ["source_database", "element_id"], sort=True
    ):
        element_rows.append(
            {
                "version": PURCHASE_ORDERABILITY_VERSION,
                "source_database": source,
                "element_id": element,
                "found_route_count": int(len(frame)),
                "orderable_route_count": int(frame.every_component_orderable.sum()),
                "every_found_route_orderable": bool(frame.every_component_orderable.all()),
                "any_orderable_route": bool(frame.every_component_orderable.any()),
                "maximum_orderable_copy_count": (
                    int(frame.loc[frame.every_component_orderable, "target_copy_count"].max())
                    if frame.every_component_orderable.any() else None
                ),
            }
        )
    elements = pd.DataFrame(element_rows)
    source_summary = (
        routes.groupby("source_database", sort=True)
        .agg(
            found_routes=("complete_route_id", "size"),
            all_components_orderable=("every_component_orderable", "sum"),
            elements_with_routes=("element_id", "nunique"),
        )
        .reset_index()
    )
    source_summary["orderable_fraction_among_found_routes"] = (
        source_summary.all_components_orderable / source_summary.found_routes
    )

    for name, table in {
        "route_purchase_orderability": routes,
        "element_purchase_orderability": elements,
        "unique_purchase_components": purchases,
        "source_purchase_orderability": source_summary,
    }.items():
        table.to_parquet(output / f"{name}.parquet", index=False)
        table.to_csv(output / f"{name}.csv", index=False)

    passed_routes = int(routes.every_component_orderable.sum())
    passed_elements = int(elements.every_found_route_orderable.sum())
    product_counts = purchases.product_class.value_counts().sort_index().to_dict()
    summary = {
        "version": PURCHASE_ORDERABILITY_VERSION,
        "purchase_policy": PURCHASE_POLICY,
        "created_at": utc_now(),
        "input": {
            "complete_route_run": raw.parent.name,
            "complete_route_table_set": raw.name,
            "shard_count": len(manifests),
            "manifest_set_sha256": hashlib.sha256(
                "".join(sha256_file(path) for path in manifests).encode()
            ).hexdigest(),
        },
        "criteria": {
            "custom_dna_oligo_length_nt_inclusive": [
                DEFAULT_OLIGO_MIN_BP, DEFAULT_OLIGO_MAX_BP
            ],
            "standard_primer_reporting_cutoff_nt": STANDARD_PRIMER_MAX_BP,
            "gblock_length_bp_inclusive": [
                DEFAULT_GBLOCK_MIN_BP, DEFAULT_GBLOCK_MAX_BP
            ],
            "gblock_live_idt_policy": IDT_SCORE_POLICY,
            "gblock_score_threshold_exclusive": IDT_ORDERABILITY_THRESHOLD,
        },
        "found_routes": int(len(routes)),
        "routes_every_component_orderable": passed_routes,
        "orderable_fraction_among_found_routes": passed_routes / len(routes),
        "elements_with_found_routes": int(len(elements)),
        "elements_every_found_route_orderable": passed_elements,
        "element_fraction": passed_elements / len(elements),
        "component_occurrences": int(len(components)),
        "unique_purchase_components": int(len(purchases)),
        "unique_purchase_product_counts": {key: int(value) for key, value in product_counts.items()},
        "unique_gblocks_live_idt_accepted": int(
            purchases.product_class.eq("gblock").mul(purchases.orderable).sum()
        ),
        "failure_counts": {
            str(key): int(value)
            for key, value in components.loc[~components.orderable, "failure_reason"]
            .value_counts().items()
        },
        "output_sha256": {},
    }
    summary_path = output / "purchase_orderability_summary.json"
    write_json_atomic(summary, summary_path)
    output_paths = sorted(output.glob("*.csv")) + sorted(output.glob("*.parquet"))
    summary["output_sha256"] = {
        path.name: sha256_file(path) for path in output_paths
    }
    write_json_atomic(summary, summary_path)
    return {
        "routes": routes,
        "elements": elements,
        "components": components,
        "purchases": purchases,
        "source_summary": source_summary,
    }
