"""User-facing, annotation-rich molecular records for exact-DNA designs."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping, Sequence

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from .design_artifacts import (
    _annotate_enzyme_sites,
    _base_record,
    _source_feature,
)
from .plasmid_reference import PlasmidReferenceDatabase, load_plasmid_reference

if TYPE_CHECKING:  # pragma: no cover
    from .exact_dna_design import ExactDNAResult


def _add_exact_target_features(
    record,
    *,
    start: int,
    sequence: str,
    query: Mapping[str, Any],
    label: str,
) -> None:
    record.features.append(
        SeqFeature(
            FeatureLocation(start, start + len(sequence), strand=1),
            type="regulatory" if query.get("input_mode") == "array" else "misc_feature",
            qualifiers={
                "label": [label],
                "feature_kind": ["exact_target_DNA"],
                "exact_seq_preserved": ["true"],
            },
        )
    )
    if query.get("input_mode") != "array":
        return
    unit = str(query.get("repeat_unit", ""))
    spacer = str(query.get("spacer", ""))
    if not unit:
        return
    cursor = 0
    ordinal = 1
    while cursor + len(unit) <= len(sequence):
        if sequence[cursor : cursor + len(unit)] == unit:
            record.features.append(
                SeqFeature(
                    FeatureLocation(start + cursor, start + cursor + len(unit), strand=1),
                    type="repeat_region",
                    qualifiers={
                        "label": [f"repeat unit {ordinal}"],
                        "repeat_number": [str(ordinal)],
                    },
                )
            )
            ordinal += 1
            cursor += len(unit)
            if spacer and sequence[cursor : cursor + len(spacer)] == spacer:
                record.features.append(
                    SeqFeature(
                        FeatureLocation(start + cursor, start + cursor + len(spacer), strand=1),
                        type="misc_feature",
                        qualifiers={"label": [f"spacer after repeat {ordinal - 1}"]},
                    )
                )
                cursor += len(spacer)
        else:
            cursor += 1


def _used_enzymes(result: "ExactDNAResult") -> list[dict[str, str]]:
    route = result.selected_route or {}
    rows: list[dict[str, str]] = []
    if route:
        rows.extend(
            [
                {
                    "enzyme": str(route["left_cutter"]),
                    "role": "vector left cutter",
                    "recognition_site": "",
                },
                {
                    "enzyme": str(route["right_cutter"]),
                    "role": "vector right cutter",
                    "recognition_site": "",
                },
            ]
        )
        for pair in route.get("pairs", []):
            rows.extend(
                [
                    {
                        "enzyme": str(pair["site_i_enzyme"]),
                        "role": "Site I",
                        "recognition_site": str(pair["site_i_recognition_site"]),
                    },
                    {
                        "enzyme": str(pair["site_ii_enzyme"]),
                        "role": "Site II",
                        "recognition_site": str(pair["site_ii_recognition_site"]),
                    },
                ]
            )
    for fragment in result.purchase_fragments:
        for side in ("left", "right"):
            enzyme = str(fragment.get(f"{side}_adapter_enzyme", ""))
            if enzyme and enzyme != "not_required":
                rows.append(
                    {"enzyme": enzyme, "role": f"Site III {side} adapter", "recognition_site": ""}
                )
    unique: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for row in rows:
        key = (row["enzyme"], row["role"])
        if key not in seen:
            seen.add(key)
            unique.append(row)
    return unique


def _add_cohesive_end_features(record, fragment: Mapping[str, Any]) -> None:
    left = str(fragment.get("left_digest_overhang", ""))
    right = str(fragment.get("right_digest_overhang", ""))
    if left:
        record.features.append(
            SeqFeature(
                FeatureLocation(0, min(len(left), len(record)), strand=1),
                type="misc_feature",
                qualifiers={
                    "label": [f"left cohesive end {left}"],
                    "generated_by": [str(fragment.get("left_adapter_enzyme") or fragment.get("vector_left_cutter", ""))],
                },
            )
        )
    if right:
        record.features.append(
            SeqFeature(
                FeatureLocation(max(0, len(record) - len(right)), len(record), strand=1),
                type="misc_feature",
                qualifiers={
                    "label": [f"right cohesive end {right}"],
                    "generated_by": [str(fragment.get("right_adapter_enzyme") or fragment.get("vector_right_cutter", ""))],
                },
            )
        )


def _remove_public_hash_metadata(record) -> None:
    """Keep reproducibility hashes in the audit ZIP, not user-facing GenBank."""
    for key in list(record.annotations):
        if "hash" in key.lower() or "sha" in key.lower():
            del record.annotations[key]
    for feature in record.features:
        for key in list(feature.qualifiers):
            if "hash" in key.lower() or "sha" in key.lower():
                del feature.qualifiers[key]


def build_exact_dna_assembly_records(
    result: "ExactDNAResult",
    *,
    plasmid_database: PlasmidReferenceDatabase | None = None,
    plasmid_reference_path: str | Path | None = None,
) -> tuple[list[tuple[str, Any]], list[dict[str, Any]]]:
    bundle = getattr(result, "_assembly_bundle", None)
    if not bundle or not bundle.get("passed") or not result.selected_route:
        return [], []
    database = plasmid_database or load_plasmid_reference(plasmid_reference_path)
    route = result.selected_route
    profile = database.profile(str(route["profile_id"]))
    reference = database.reference(profile.reference_id)
    scheme = next(item for item in database.schemes if item.scheme_id == route["scheme_id"])
    enzymes = _used_enzymes(result)
    origin = int(bundle["rotation_origin_0based"])
    retained_length = len(str(bundle["retained_backbone_sequence"]))

    step00 = _base_record(
        str(bundle["initial_plasmid_sequence"]),
        "step00_plasmid",
        f"Initial {profile.profile_id} plasmid rotated to the cloning origin",
        circular=True,
    )
    step00.features[0].qualifiers.update(
        {"profile_id": [profile.profile_id], "assembly_step": ["0"], "role": ["initial_plasmid"]}
    )
    for feature in reference.features:
        mapped = _source_feature(
            feature,
            reference_length=len(reference.sequence),
            expression_strand=profile.expression_strand,
            origin=origin,
            retained_length=None,
        )
        if mapped is not None:
            step00.features.append(mapped)
    _annotate_enzyme_sites(step00, enzymes, circular=True, step_role="initial_plasmid")
    records: list[tuple[str, Any]] = [("step00_plasmid.gb", step00)]
    manifest: list[dict[str, Any]] = [
        {
            "step": 0,
            "molecule": "plasmid",
            "file": "step00_plasmid.gb",
            "role": "initial_plasmid",
            "length_bp": len(step00),
            "sequence_sha256": hashlib.sha256(str(step00.seq).encode()).hexdigest(),
            "cloning_region_start_0based": 0,
            "cloning_region_end_0based_exclusive": len(step00),
        }
    ]

    fragments = {str(row["fragment_id"]): row for row in result.purchase_fragments}
    for step in bundle["steps"]:
        number = int(step["step"])
        fragment_id = str(step["purchase_fragment_ids"][0])
        fragment = fragments[fragment_id]
        product_type = str(fragment.get("product_type", ""))
        primer_pair = product_type in {
            "annealed_sticky_end_primer_pair",
            "duplexed_seed_oligo_pair",
        }
        insert_sequence = (
            str(step["insert_sequence"])
            if primer_pair
            else str(fragment.get("purchase_sequence", step["insert_sequence"]))
        )
        insert = _base_record(
            insert_sequence,
            f"step{number:02d}_insert",
            (
                f"Assembly-ready annealed insert for cloning step {number}"
                if primer_pair
                else f"Purchased insert before Site-III preparation for cloning step {number}"
            ),
            circular=False,
        )
        insert.features[0].qualifiers.update(
            {
                "assembly_step": [str(number)],
                "purchase_fragment_id": [fragment_id],
                "molecule_role": ["assembly_ready_insert"],
                "preparation": [
                    "anneal complementary oligos"
                    if primer_pair
                    else "digest disposable Site-III adapters before ligation"
                ],
            }
        )
        left = len(
            str(
                fragment.get(
                    "left_digest_overhang" if primer_pair else "left_adapter_sequence",
                    "",
                )
            )
        )
        core = str(step["insert_core_sequence"])
        target_feature_start = left
        target_feature_sequence = core
        if number == 1:
            left_restoration = str(route.get("left_restoration_sequence", ""))
            right_restoration = str(route.get("right_restoration_sequence", ""))
            if left_restoration:
                insert.features.append(
                    SeqFeature(
                        FeatureLocation(left, left + len(left_restoration), strand=1),
                        type="misc_feature",
                        qualifiers={"label": ["left vector restoration segment"]},
                    )
                )
            if right_restoration:
                right_start = left + len(core) - len(right_restoration)
                insert.features.append(
                    SeqFeature(
                        FeatureLocation(right_start, right_start + len(right_restoration), strand=1),
                        type="misc_feature",
                        qualifiers={"label": ["right vector restoration segment"]},
                    )
                )
            target_feature_start += len(left_restoration)
            target_feature_sequence = str(result.seed.get("seed_sequence", "")) if result.seed else core
        _add_exact_target_features(
            insert,
            start=target_feature_start,
            sequence=target_feature_sequence,
            query=result.query,
            label=f"step {number} inserted exact core",
        )
        if primer_pair:
            _add_cohesive_end_features(insert, fragment)
        else:
            if left:
                insert.features.append(
                    SeqFeature(
                        FeatureLocation(0, left, strand=1),
                        type="misc_feature",
                        qualifiers={
                            "label": [f"disposable left Site-III adapter ({fragment.get('left_adapter_enzyme', '')})"],
                            "disposable": ["true"],
                        },
                    )
                )
            right_length = len(str(fragment.get("right_adapter_sequence", "")))
            if right_length:
                insert.features.append(
                    SeqFeature(
                        FeatureLocation(len(insert) - right_length, len(insert), strand=1),
                        type="misc_feature",
                        qualifiers={
                            "label": [f"disposable right Site-III adapter ({fragment.get('right_adapter_enzyme', '')})"],
                            "disposable": ["true"],
                        },
                    )
                )
        _annotate_enzyme_sites(insert, enzymes, circular=False, step_role="insert")
        for latent in step.get("latent_sites", []):
            insert.features[0].qualifiers.setdefault("latent_site_used", []).append(
                f"{latent['role']}:{latent['enzyme']}:{latent['target_base']}->{latent['temporary_active_base']}"
            )
            mismatch = int(latent["mismatch_position_in_result_0based"])
            fragment_start = int(fragment.get("target_start", -1))
            fragment_end = int(fragment.get("target_end", -1))
            if fragment_start <= mismatch < fragment_end:
                position = left + mismatch - fragment_start
                insert.features.append(
                    SeqFeature(
                        FeatureLocation(position, position + 1, strand=1),
                        type="misc_feature",
                        qualifiers={
                            "label": [f"latent {latent['enzyme']} activation base"],
                            "enzyme": [str(latent["enzyme"])],
                            "role": [str(latent["role"])],
                            "target_base": [str(latent["target_base"])],
                            "temporary_base": [str(latent["temporary_active_base"])],
                            "donor_derived": ["true"],
                        },
                    )
                )
        insert_name = f"step{number:02d}_insert.gb"
        records.append((insert_name, insert))
        manifest.append(
            {
                "step": number,
                "molecule": "insert",
                "file": insert_name,
                "role": str(step["stage"]),
                "purchase_fragment_id": fragment_id,
                "length_bp": len(insert),
                "sequence_sha256": hashlib.sha256(str(insert.seq).encode()).hexdigest(),
                "cloning_region_start_0based": target_feature_start,
                "cloning_region_end_0based_exclusive": (
                    target_feature_start + len(target_feature_sequence)
                ),
            }
        )

        plasmid_sequence = str(step["plasmid_sequence"])
        plasmid_name = f"step{number:02d}_plasmid.gb"
        plasmid = _base_record(
            plasmid_sequence,
            f"step{number:02d}_plasmid",
            f"Circular {profile.profile_id} plasmid after cloning step {number}",
            circular=True,
        )
        plasmid.features[0].qualifiers.update(
            {
                "profile_id": [profile.profile_id],
                "assembly_step": [str(number)],
                "independent_verify": ["passed"],
            }
        )
        for feature in reference.features:
            mapped = _source_feature(
                feature,
                reference_length=len(reference.sequence),
                expression_strand=profile.expression_strand,
                origin=origin,
                retained_length=retained_length,
            )
            if mapped is not None:
                plasmid.features.append(mapped)
        cursor = retained_length
        if scheme.left_restoration_sequence:
            plasmid.features.append(
                SeqFeature(
                    FeatureLocation(cursor, cursor + len(scheme.left_restoration_sequence), strand=1),
                    type="misc_feature",
                    qualifiers={"label": ["left vector restoration segment"]},
                )
            )
            cursor += len(scheme.left_restoration_sequence)
        current_insert = str(step["result_insert_sequence"])
        _add_exact_target_features(
            plasmid,
            start=cursor,
            sequence=current_insert,
            query=result.query,
            label=f"exact insert after cloning step {number}",
        )
        for latent in step.get("latent_sites", []):
            position = cursor + int(latent["mismatch_position_in_intermediate_0based"])
            plasmid.features.append(
                SeqFeature(
                    FeatureLocation(position, position + 1, strand=1),
                    type="misc_feature",
                    qualifiers={
                        "label": [f"temporarily active latent {latent['enzyme']} site"],
                        "enzyme": [str(latent["enzyme"])],
                        "role": [str(latent["role"])],
                        "target_base": [str(latent["target_base"])],
                        "temporary_base": [str(latent["temporary_active_base"])],
                    },
                )
            )
        _annotate_enzyme_sites(plasmid, enzymes, circular=True, step_role="assembled_plasmid")
        records.append((plasmid_name, plasmid))
        manifest.append(
            {
                "step": number,
                "molecule": "plasmid",
                "file": plasmid_name,
                "role": "assembled_plasmid",
                "length_bp": len(plasmid),
                "insert_length_bp": len(current_insert),
                "independent_verification": "passed",
                "sequence_sha256": hashlib.sha256(str(plasmid.seq).encode()).hexdigest(),
                "cloning_region_start_0based": cursor,
                "cloning_region_end_0based_exclusive": cursor + len(current_insert),
            }
        )
    return records, manifest


def write_exact_dna_genbanks(
    result: "ExactDNAResult",
    output_dir: str | Path,
    *,
    plasmid_database: PlasmidReferenceDatabase | None = None,
    plasmid_reference_path: str | Path | None = None,
    include_manifest: bool = False,
    export_maps: bool = False,
) -> dict[str, str]:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    records, manifest = build_exact_dna_assembly_records(
        result,
        plasmid_database=plasmid_database,
        plasmid_reference_path=plasmid_reference_path,
    )
    paths: dict[str, str] = {}
    for filename, record in records:
        _remove_public_hash_metadata(record)
        path = destination / filename
        SeqIO.write(record, path, "genbank")
        paths[filename] = str(path)
    if records:
        # The final target is clearer as a separate record than the final donor.
        target_record = _base_record(
            result.final_insert_sequence,
            "final_exact_insert",
            "Final exact user-requested DNA insert",
            circular=False,
        )
        _add_exact_target_features(
            target_record,
            start=0,
            sequence=result.final_insert_sequence,
            query=result.query,
            label="final exact target DNA",
        )
        target_path = destination / "final_exact_insert.gb"
        _remove_public_hash_metadata(target_record)
        SeqIO.write(target_record, target_path, "genbank")
        paths["final_exact_insert.gb"] = str(target_path)
        final_plasmid_source = next(
            record for filename, record in records if filename == f"step{len(result.cloning_steps):02d}_plasmid.gb"
        )
        final_path = destination / "final_plasmid.gb"
        SeqIO.write(final_plasmid_source, final_path, "genbank")
        paths["final_plasmid.gb"] = str(final_path)
    if include_manifest:
        manifest_path = destination / "assembly_step_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        paths["assembly_step_manifest.json"] = str(manifest_path)
        legacy_path = destination / "step_genbank_index.json"
        legacy_path.write_text(manifest_path.read_text())
        paths["step_genbank_index.json"] = str(legacy_path)
        csv_path = destination / "assembly_step_manifest.csv"
        columns = sorted({key for row in manifest for key in row})
        with csv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows({key: row.get(key, "") for key in columns} for row in manifest)
        paths["assembly_step_manifest.csv"] = str(csv_path)
    if export_maps and records:
        try:
            import matplotlib.pyplot as plt
            from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
        except ImportError:  # pragma: no cover - optional outside notebook installs
            return paths
        map_dir = destination / "maps"
        map_dir.mkdir(exist_ok=True)
        for filename, record in records:
            views = (
                ("circular", "linear")
                if record.annotations.get("topology") == "circular"
                else ("linear",)
            )
            for view in views:
                def feature_filter(feature):
                    if feature.type == "source":
                        return False
                    if view != "circular":
                        return True
                    qualifiers = feature.qualifiers
                    if qualifiers.get("feature_kind", [""])[0] in {
                        "restriction_site",
                        "exact_target_DNA",
                    }:
                        return True
                    if feature.type in {"regulatory", "repeat_region"}:
                        return True
                    return qualifiers.get("feature_class", [""])[0] in {
                        "antibiotic_resistance",
                        "origin",
                        "replication_origin",
                        "promoter",
                        "terminator",
                        "operator",
                    }

                translator = BiopythonTranslator(features_filters=(feature_filter,))
                if view == "circular":
                    graphic = translator.translate_record(
                        record, record_class=CircularGraphicRecord
                    )
                    figure, axis = plt.subplots(figsize=(8, 8))
                    graphic.plot(ax=axis)
                else:
                    graphic = translator.translate_record(record)
                    figure, axis = plt.subplots(figsize=(13, 3.6))
                    graphic.plot(ax=axis, figure_width=13)
                axis.set_title(f"{filename} · {view}")
                figure.tight_layout()
                try:
                    for suffix in ("png", "svg"):
                        path = map_dir / f"{Path(filename).stem}.{view}.{suffix}"
                        figure.savefig(
                            path,
                            dpi=200 if suffix == "png" else None,
                            bbox_inches="tight",
                            facecolor="white",
                        )
                        paths[f"map_{Path(filename).stem}_{view}_{suffix}"] = str(path)
                finally:
                    plt.close(figure)
    return paths
