"""Auditable GenBank timelines and static maps for HURDLER designs.

The optimization controller intentionally stores plain dictionaries.  This
module turns an accepted result into molecular records only at export time, so
the JSON result remains portable and every GenBank file can be regenerated.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, Mapping, Sequence

from Bio import Restriction, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .dna_assembly import _cut_offsets
from .optimization import reverse_complement, translate_dna
from .plasmid_reference import PlasmidFeature, load_plasmid_reference

if TYPE_CHECKING:  # pragma: no cover - imported only for static analysis
    from .vector_design import DesignResultV2


def _sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def _hash_parts(digest: str) -> list[str]:
    """Keep a SHA256 lossless across GenBank's line-wrapping parser."""
    value = str(digest)
    return [value[:32], value[32:]]


def _safe_identifier(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._-")
    return cleaned or "sequence"


def write_secondary_checkpoint(payload: Mapping[str, Any], output_zip: str | Path) -> Path:
    """Atomically save the latest accepted secondary without any credentials."""
    destination = Path(output_zip)
    destination.parent.mkdir(parents=True, exist_ok=True)
    core = str(payload.get("core_sequence") or "")
    purchase = str(payload.get("purchase_sequence") or "")
    score = payload.get("idt_complexity_score")
    accepted = bool(
        core
        and purchase
        and payload.get("event") == "accepted_secondary"
        and payload.get("validation_mode") == "api"
        and isinstance(score, (int, float))
        and not isinstance(score, bool)
        and float(score) < 10
    )
    public = {
        key: value
        for key, value in payload.items()
        if key not in {"core_sequence", "purchase_sequence"}
    }
    public["accepted_secondary_available"] = accepted
    with tempfile.TemporaryDirectory(prefix="hurdler_checkpoint_") as temporary:
        root = Path(temporary)
        (root / "checkpoint.json").write_text(json.dumps(public, indent=2, sort_keys=True) + "\n")
        if accepted:
            sequence_id = _safe_identifier(str(payload.get("sequence_id") or "sequence"))
            copies = int(payload["repeat_copies"])
            (root / "best_secondary_core.fasta").write_text(
                f">{sequence_id}_secondary_{copies}copies_core\n{core}\n"
            )
            (root / "best_secondary_purchase.fasta").write_text(
                f">{sequence_id}_secondary_{copies}copies_purchase\n{purchase}\n"
            )
        temporary_zip = destination.with_suffix(destination.suffix + ".tmp")
        with zipfile.ZipFile(temporary_zip, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(root.iterdir()):
                archive.write(path, path.name)
        temporary_zip.replace(destination)
    return destination


def timestamped_results_archive(
    output_dir: str | Path,
    archive_dir: str | Path,
    *,
    sequence_id: str,
    timestamp: datetime | None = None,
) -> Path:
    """Create the final UTC-stamped ZIP, suitable for optional Drive sync."""
    moment = timestamp or datetime.now(timezone.utc)
    stamp = moment.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    destination = Path(archive_dir)
    destination.mkdir(parents=True, exist_ok=True)
    stem = destination / f"hurdler_{_safe_identifier(sequence_id)}_{stamp}_results"
    archive = shutil.make_archive(str(stem), "zip", root_dir=Path(output_dir))
    return Path(archive)


def _parts_from_positions(positions: Iterable[int], *, strand: int = 1):
    ordered = sorted(set(int(value) for value in positions))
    if not ordered:
        return None
    runs: list[tuple[int, int]] = []
    start = previous = ordered[0]
    for value in ordered[1:]:
        if value != previous + 1:
            runs.append((start, previous + 1))
            start = value
        previous = value
    runs.append((start, previous + 1))
    locations = [FeatureLocation(start, end, strand=strand or None) for start, end in runs]
    return locations[0] if len(locations) == 1 else CompoundLocation(locations, operator="join")


def _source_feature(
    feature: PlasmidFeature,
    *,
    reference_length: int,
    expression_strand: int,
    origin: int,
    retained_length: int | None,
) -> SeqFeature | None:
    output_positions: list[int] = []
    for start, end in feature.intervals:
        for physical in range(int(start), int(end)):
            oriented = physical if expression_strand == 1 else reference_length - 1 - physical
            output = (oriented - origin) % reference_length
            if retained_length is None or output < retained_length:
                output_positions.append(output)
    strand = int(feature.strand) * int(expression_strand)
    location = _parts_from_positions(output_positions, strand=strand)
    if location is None:
        return None
    expected = sum(end - start for start, end in feature.intervals)
    qualifiers: dict[str, list[str]] = {
        "label": [feature.label],
        "feature_id": [feature.feature_id],
        "feature_class": [feature.feature_class],
        "source_provider": list(feature.sources),
        "protected": [str(bool(feature.protected)).lower()],
    }
    if len(output_positions) != expected:
        qualifiers["clipped_by_assembly"] = ["true"]
    return SeqFeature(location, type=feature.feature_type, qualifiers=qualifiers)


def _site_occurrences(sequence: str, motif: str, *, circular: bool) -> list[tuple[int, str]]:
    motif = motif.upper()
    if not motif:
        return []
    patterns = [(motif, "+")]
    reverse = reverse_complement(motif)
    if reverse != motif:
        patterns.append((reverse, "-"))
    extended = sequence + sequence[: len(motif) - 1] if circular else sequence
    hits: set[tuple[int, str]] = set()
    for pattern, strand in patterns:
        for start in range(len(sequence)):
            if extended.startswith(pattern, start):
                hits.add((start, strand))
    return sorted(hits)


def _site_location(start: int, length: int, sequence_length: int, strand: int):
    end = start + length
    if end <= sequence_length:
        return FeatureLocation(start, end, strand=strand)
    return CompoundLocation(
        [
            FeatureLocation(start, sequence_length, strand=strand),
            FeatureLocation(0, end - sequence_length, strand=strand),
        ],
        operator="join",
    )


def _enzyme_geometry(enzyme: str, fallback_site: str = "") -> dict[str, Any]:
    try:
        restriction = getattr(Restriction, enzyme)
        site = str(restriction.site)
        top, bottom, elucidate = _cut_offsets(enzyme)
        return {
            "recognition_site": site,
            "top_cut_offset": top,
            "bottom_cut_offset": bottom,
            "sticky_end": str(restriction.ovhgseq or ""),
            "overhang_length": int(restriction.ovhg),
            "elucidate": elucidate,
        }
    except (AttributeError, TypeError, ValueError):
        return {
            "recognition_site": fallback_site,
            "top_cut_offset": 0,
            "bottom_cut_offset": 0,
            "sticky_end": "",
            "overhang_length": 0,
            "elucidate": fallback_site,
        }


def _annotate_enzyme_sites(
    record: SeqRecord,
    enzymes: Sequence[Mapping[str, str]],
    *,
    circular: bool,
    step_role: str,
) -> list[dict[str, Any]]:
    audit: list[dict[str, Any]] = []
    for item in enzymes:
        enzyme = str(item["enzyme"])
        role = str(item["role"])
        geometry = _enzyme_geometry(enzyme, str(item.get("recognition_site", "")))
        motif = str(geometry["recognition_site"])
        hits = _site_occurrences(str(record.seq), motif, circular=circular)
        audit.append({
            "enzyme": enzyme,
            "role": role,
            "recognition_site": motif,
            "occurrence_count": len(hits),
            "molecule_role": step_role,
        })
        for ordinal, (start, orientation) in enumerate(hits, 1):
            strand = 1 if orientation == "+" else -1
            top_offset = int(geometry["top_cut_offset"])
            bottom_offset = int(geometry["bottom_cut_offset"])
            if strand == -1:
                top_offset, bottom_offset = len(motif) - bottom_offset, len(motif) - top_offset
            record.features.append(
                SeqFeature(
                    _site_location(start, len(motif), len(record), strand),
                    type="misc_feature",
                    qualifiers={
                        "label": [f"{role}: {enzyme}"],
                        "feature_kind": ["restriction_site"],
                        "enzyme": [enzyme],
                        "step_role": [role],
                        "recognition_sequence": [motif],
                        "orientation": [orientation],
                        "top_cut_0based": [str((start + top_offset) % len(record))],
                        "bottom_cut_0based": [str((start + bottom_offset) % len(record))],
                        "sticky_end": [str(geometry["sticky_end"])],
                        "overhang_length": [str(geometry["overhang_length"])],
                        "elucidate": [str(geometry["elucidate"])],
                        "occurrence": [str(ordinal)],
                    },
                )
            )
    return audit


def _base_record(sequence: str, record_id: str, description: str, *, circular: bool) -> SeqRecord:
    record = SeqRecord(Seq(sequence), id=record_id[:16], name=record_id[:16], description=description)
    record.annotations.update(
        {
            "molecule_type": "DNA",
            "topology": "circular" if circular else "linear",
            "data_file_division": "SYN",
        }
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(0, len(sequence), strand=1),
            type="source",
            qualifiers={
                "label": [record_id],
                "molecule_sha256": _hash_parts(_sha(sequence)),
                "molecule_topology": ["circular" if circular else "linear"],
                "design_only": ["true"],
            },
        )
    )
    return record


def _add_coding_features(
    record: SeqRecord,
    *,
    cds_start: int,
    dna: str,
    n_cap_aa: int,
    module_aa: int,
    copies: int,
    c_cap_aa: int,
    label: str,
) -> None:
    protein = translate_dna(dna)
    record.features.append(
        SeqFeature(
            FeatureLocation(cds_start, cds_start + len(dna), strand=1),
            type="CDS",
            qualifiers={
                "label": [label],
                "translation": [protein],
                "codon_start": ["1"],
                "transl_table": ["11"],
                "repeat_copy_count": [str(copies)],
                "dna_sha256": _hash_parts(_sha(dna)),
            },
        )
    )
    cursor = cds_start
    if n_cap_aa:
        end = cursor + n_cap_aa * 3
        record.features.append(SeqFeature(FeatureLocation(cursor, end, strand=1), type="misc_feature", qualifiers={"label": ["N-terminal cap"]}))
        cursor = end
    for index in range(copies):
        end = cursor + module_aa * 3
        record.features.append(
            SeqFeature(
                FeatureLocation(cursor, end, strand=1),
                type="repeat_region",
                qualifiers={"label": [f"repeat unit {index + 1}"], "repeat_number": [str(index + 1)]},
            )
        )
        cursor = end
    if c_cap_aa:
        record.features.append(SeqFeature(FeatureLocation(cursor, cursor + c_cap_aa * 3, strand=1), type="misc_feature", qualifiers={"label": ["C-terminal cap"]}))


def _fragment_record(
    fragment: Mapping[str, Any],
    *,
    record_id: str,
    fragment_kind: str,
    copies: int,
    query: Mapping[str, Any],
    enzymes: Sequence[Mapping[str, str]],
    reuse_round: int = 1,
) -> tuple[SeqRecord, dict[str, Any]]:
    purchase = str(fragment["purchase_sequence"])
    core_start = int(fragment["core_start_bp"])
    core_end = int(fragment["core_end_bp"])
    record = _base_record(purchase, record_id, f"HURDLER {fragment_kind} purchase DNA", circular=False)
    record.annotations["fragment_kind"] = fragment_kind
    record.annotations["purchase_sha256"] = str(fragment["purchase_sha256"])
    record.features[0].qualifiers.update({
        "fragment_kind": [fragment_kind],
        "purchase_sha256": _hash_parts(str(fragment["purchase_sha256"])),
        "repeat_copy_count": [str(copies)],
        "reuse_round": [str(reuse_round)],
    })
    if core_start:
        record.features.append(
            SeqFeature(
                FeatureLocation(0, core_start),
                type="misc_feature",
                qualifiers={"label": ["disposable left adapter/restoration"], "disposable": ["true"]},
            )
        )
    if core_end < len(purchase):
        record.features.append(
            SeqFeature(
                FeatureLocation(core_end, len(purchase)),
                type="misc_feature",
                qualifiers={"label": ["disposable right adapter/restoration"], "disposable": ["true"]},
            )
        )
    core = purchase[core_start:core_end]
    if fragment_kind == "primary":
        _add_coding_features(
            record,
            cds_start=core_start,
            dna=core,
            n_cap_aa=len(str(query["n_cap"])),
            module_aa=len(str(query["repeat_module"])),
            copies=copies,
            c_cap_aa=len(str(query["c_cap"])),
            label="primary repeat-protein insert",
        )
    else:
        _add_coding_features(
            record,
            cds_start=core_start,
            dna=core,
            n_cap_aa=0,
            module_aa=len(str(query["repeat_module"])),
            copies=copies,
            c_cap_aa=0,
            label="reusable secondary repeat insert",
        )
        record.annotations["reuse_round"] = reuse_round
    site_audit = _annotate_enzyme_sites(record, enzymes, circular=False, step_role="insert")
    return record, {
        "sequence_sha256": _sha(purchase),
        "length_bp": len(purchase),
        "core_sha256": _sha(core),
        "cloning_region_start_0based": core_start,
        "cloning_region_end_0based_exclusive": core_end,
        "core_length_bp": len(core),
        "translation": translate_dna(core),
        "site_audit": site_audit,
    }


def build_assembly_records(
    result: "DesignResultV2",
    *,
    plasmid_reference_path: str | Path | None = None,
) -> tuple[list[tuple[str, SeqRecord]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Build every starting, insert, and intermediate plasmid record."""
    if not result.final_plasmid or not result.final_dna_sequence or not result.selected_route:
        return [], [], []
    request = result.request
    query = request["query"]
    route = result.selected_route
    candidate = next(
        row for row in result.protein_candidates if row["candidate_id"] == route["candidate_id"]
    )
    database = load_plasmid_reference(plasmid_reference_path)
    profile = database.profile(str(route["profile_id"]))
    reference = database.reference(profile.reference_id)
    scheme = next(row for row in database.schemes if row.scheme_id == route["scheme_id"])
    if scheme.left_cutter is None or scheme.right_cutter is None:
        raise ValueError("Selected cut scheme has no vector cutters")
    site_iii = str(request["selection"]["site_iii_enzyme"])
    enzymes = [
        {"enzyme": scheme.left_cutter.canonical_enzyme, "role": "vector left cutter", "recognition_site": scheme.left_cutter.recognition_site},
        {"enzyme": scheme.right_cutter.canonical_enzyme, "role": "vector right cutter", "recognition_site": scheme.right_cutter.recognition_site},
        {"enzyme": str(candidate["site_i_enzyme"]), "role": "Site I", "recognition_site": str(candidate["site_i_recognition_site"])},
        {"enzyme": str(candidate["site_ii_enzyme"]), "role": "Site II", "recognition_site": str(candidate["site_ii_recognition_site"])},
        {"enzyme": site_iii, "role": "Site III", "recognition_site": ""},
    ]
    # Avoid duplicate features when a vector cutter is intentionally reused.
    unique_enzymes: list[dict[str, str]] = []
    for enzyme in enzymes:
        if (enzyme["enzyme"], enzyme["role"]) not in {
            (row["enzyme"], row["role"]) for row in unique_enzymes
        }:
            unique_enzymes.append(enzyme)

    length = len(reference.sequence)
    oriented = reference.sequence if profile.expression_strand == 1 else reverse_complement(reference.sequence)
    origin = int(scheme.right_cutter.top_cut_oriented)
    initial = oriented[origin:] + oriented[:origin]
    step00 = _base_record(initial, "step00_plasmid", f"Initial {profile.profile_id} plasmid, rotated to cloning origin", circular=True)
    step00.features[0].qualifiers.update({
        "profile_id": [profile.profile_id],
        "source_sha256": _hash_parts(reference.sequence_sha256),
        "rotation_origin_bp": [str(origin)],
        "correctness_qc": ["circular rotation of the complete source plasmid"],
    })
    for feature in reference.features:
        mapped = _source_feature(
            feature,
            reference_length=length,
            expression_strand=profile.expression_strand,
            origin=origin,
            retained_length=None,
        )
        if mapped is not None:
            step00.features.append(mapped)
    initial_site_audit = _annotate_enzyme_sites(step00, unique_enzymes, circular=True, step_role="initial_plasmid")

    records: list[tuple[str, SeqRecord]] = [("step00_plasmid.gb", step00)]
    manifest: list[dict[str, Any]] = [{
        "step": 0,
        "molecule": "plasmid",
        "file": "step00_plasmid.gb",
        "length_bp": len(step00),
        "sequence_sha256": _sha(str(step00.seq)),
        "copy_count": 0,
        "translation_exact": True,
        "site_audit": initial_site_audit,
        "cloning_region_start_0based": len(scheme.retained_backbone_sequence),
        "cloning_region_end_0based_exclusive": len(step00),
    }]
    translations: list[dict[str, Any]] = []

    primary_fragment = result.primary_fragments[0]
    primary_copies = int(result.rdl_plan.get("primary_repeat_copies", query["repeat_copies"]))
    primary_record, primary_meta = _fragment_record(
        primary_fragment,
        record_id="step01_insert",
        fragment_kind="primary",
        copies=primary_copies,
        query=query,
        enzymes=unique_enzymes,
    )
    records.append(("step01_insert.gb", primary_record))
    manifest.append({"step": 1, "molecule": "insert", "file": "step01_insert.gb", "copy_count": primary_copies, **primary_meta})

    primary_core = str(primary_fragment["purchase_sequence"])[
        int(primary_fragment["core_start_bp"]):int(primary_fragment["core_end_bp"])
    ]
    secondary_core = ""
    secondary_copies = 0
    rounds = int(result.rdl_plan.get("secondary_reuse_count", 0))
    if rounds:
        secondary_fragment = result.secondary_fragments[0]
        secondary_core = str(secondary_fragment["purchase_sequence"])[
            int(secondary_fragment["core_start_bp"]):int(secondary_fragment["core_end_bp"])
        ]
        secondary_copies = int(result.rdl_plan["secondary_repeat_copies"])

    current_cds = primary_core
    module_bp = len(str(query["repeat_module"])) * 3
    insertion_offset = (len(str(query["n_cap"])) + primary_copies * len(str(query["repeat_module"]))) * 3
    for step in range(1, rounds + 2):
        if step > 1:
            insert_record, insert_meta = _fragment_record(
                result.secondary_fragments[0],
                record_id=f"step{step:02d}_insert",
                fragment_kind="secondary",
                copies=secondary_copies,
                query=query,
                enzymes=unique_enzymes,
                reuse_round=step - 1,
            )
            insert_name = f"step{step:02d}_insert.gb"
            records.append((insert_name, insert_record))
            manifest.append({
                "step": step,
                "molecule": "insert",
                "file": insert_name,
                "copy_count": secondary_copies,
                "reuse_round": step - 1,
                "reused_purchase_sha256": str(result.secondary_fragments[0]["purchase_sha256"]),
                **insert_meta,
            })
            current_cds = current_cds[:insertion_offset] + secondary_core + current_cds[insertion_offset:]
            insertion_offset += len(secondary_core)

        current_copies = primary_copies + max(0, step - 1) * secondary_copies
        plasmid_sequence = (
            scheme.retained_backbone_sequence
            + scheme.left_restoration_sequence
            + current_cds
            + scheme.right_restoration_sequence
        )
        cds_start = len(scheme.retained_backbone_sequence) + len(scheme.left_restoration_sequence)
        plasmid_name = f"step{step:02d}_plasmid.gb"
        plasmid = _base_record(
            plasmid_sequence,
            f"step{step:02d}_plasmid",
            f"HURDLER plasmid after assembly step {step}",
            circular=True,
        )
        plasmid.features[0].qualifiers.update({
            "profile_id": [profile.profile_id],
            "assembly_step": [str(step)],
            "repeat_copy_count": [str(current_copies)],
            "cds_sha256": _hash_parts(_sha(current_cds)),
            "correctness_qc": ["translation exact; molecular simulation exact"],
        })
        for feature in reference.features:
            mapped = _source_feature(
                feature,
                reference_length=length,
                expression_strand=profile.expression_strand,
                origin=origin,
                retained_length=len(scheme.retained_backbone_sequence),
            )
            if mapped is not None:
                plasmid.features.append(mapped)
        _add_coding_features(
            plasmid,
            cds_start=cds_start,
            dna=current_cds,
            n_cap_aa=len(str(query["n_cap"])),
            module_aa=len(str(query["repeat_module"])),
            copies=current_copies,
            c_cap_aa=len(str(query["c_cap"])),
            label=f"repeat-protein CDS ({current_copies} copies)",
        )
        site_audit = _annotate_enzyme_sites(plasmid, unique_enzymes, circular=True, step_role="assembled_plasmid")
        translated = translate_dna(current_cds)
        expected = str(query["n_cap"]) + str(query["repeat_module"]) * current_copies + str(query["c_cap"])
        translation_exact = translated == expected
        if not translation_exact:
            raise AssertionError(f"Step {step} GenBank translation does not match its protein")
        records.append((plasmid_name, plasmid))
        manifest.append({
            "step": step,
            "molecule": "plasmid",
            "file": plasmid_name,
            "length_bp": len(plasmid),
            "sequence_sha256": _sha(plasmid_sequence),
            "copy_count": current_copies,
            "translation_exact": translation_exact,
            "protein_sha256": _sha(translated),
            "site_audit": site_audit,
            "cloning_region_start_0based": cds_start,
            "cloning_region_end_0based_exclusive": cds_start + len(current_cds),
        })
        translations.append({
            "step": step,
            "copy_count": current_copies,
            "cds_length_bp": len(current_cds),
            "protein_length_aa": len(translated),
            "translation": translated,
            "translation_exact": translation_exact,
            "dna_sha256": _sha(current_cds),
            "protein_sha256": _sha(translated),
        })

    final_step = next(row for row in reversed(manifest) if row["molecule"] == "plasmid")
    final_matches = final_step["sequence_sha256"] == str(result.final_plasmid["final_plasmid_sha256"])
    if not final_matches:
        raise AssertionError("Step timeline final plasmid differs from molecular simulation")
    final_step["matches_final_plasmid"] = True
    return records, manifest, translations


def write_assembly_artifacts(
    result: "DesignResultV2",
    output_dir: str | Path,
    *,
    plasmid_reference_path: str | Path | None = None,
    export_maps: bool = True,
) -> dict[str, str]:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    records, manifest, translations = build_assembly_records(
        result, plasmid_reference_path=plasmid_reference_path
    )
    if not records:
        return {}
    paths: dict[str, str] = {}
    for filename, record in records:
        path = destination / filename
        SeqIO.write(record, path, "genbank")
        paths[filename.replace(".", "_")] = str(path)
    result.assembly_steps = manifest
    manifest_path = destination / "assembly_step_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    paths["assembly_step_manifest_json"] = str(manifest_path)
    import csv

    csv_path = destination / "assembly_step_manifest.csv"
    with csv_path.open("w", newline="") as handle:
        fields = sorted({key for row in manifest for key in row if key != "site_audit"})
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fields} for row in manifest)
    paths["assembly_step_manifest_csv"] = str(csv_path)
    translation_path = destination / "step_translations.csv"
    with translation_path.open("w", newline="") as handle:
        fields = sorted({key for row in translations for key in row})
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(translations)
    paths["step_translations_csv"] = str(translation_path)
    if export_maps:
        try:
            from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
            import matplotlib.pyplot as plt
        except ImportError:  # optional outside notebook installs
            return paths
        map_dir = destination / "maps"
        map_dir.mkdir(exist_ok=True)
        for filename, record in records:
            translator = BiopythonTranslator(features_filters=(lambda feature: feature.type != "source",))
            views = ("circular", "linear") if record.annotations.get("topology") == "circular" else ("linear",)
            for view in views:
                if view == "circular":
                    graphic = translator.translate_record(record, record_class=CircularGraphicRecord)
                    figure, axis = plt.subplots(figsize=(8, 8))
                    graphic.plot(ax=axis)
                else:
                    graphic = translator.translate_record(record)
                    figure, axis = plt.subplots(figsize=(12, 3))
                    graphic.plot(ax=axis, figure_width=12)
                axis.set_title(f"{filename} · {view}")
                figure.tight_layout()
                for suffix in ("png", "svg"):
                    path = map_dir / f"{Path(filename).stem}.{view}.{suffix}"
                    figure.savefig(path, dpi=200, bbox_inches="tight")
                    paths[f"map_{Path(filename).stem}_{view}_{suffix}"] = str(path)
                plt.close(figure)
    return paths
