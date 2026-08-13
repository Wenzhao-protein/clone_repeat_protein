from __future__ import annotations

import hashlib
import json
import zipfile
from datetime import datetime, timezone

from Bio import SeqIO

from hurdler.design_artifacts import (
    build_step00_plasmid_record,
    timestamped_results_archive,
    write_secondary_checkpoint,
)
from hurdler.optimization import translate_dna
from hurdler.vector_design import (
    DESIGN_SCHEMA_VERSION_V2,
    CompatibilityQuery,
    DesignRequestV2,
    DesignSelection,
    design_construct_v2,
    design_query,
    write_design_outputs_v2,
)


def _accepted_direct_design():
    query = CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        input_mode="split",
        sequence_id="artifact_roundtrip",
        repeat_module="ACDEFGHIKLMNPQRSTVWY",
        repeat_copies=3,
    )
    route = design_query(query).vector_routes[0]
    request = DesignRequestV2(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=query,
        selection=DesignSelection(
            route["candidate_id"], route["profile_id"], route["scheme_id"],
            route["site_iii_options"][0],
        ),
        validation_mode="batch",
        population_size=4,
        generations_per_feedback_round=1,
        minimum_secondary_copies=100,
    )
    return design_construct_v2(request)


def test_step_genbanks_roundtrip_translation_annotations_and_final_identity(tmp_path):
    result = _accepted_direct_design()
    paths = write_design_outputs_v2(result, tmp_path)
    step00 = SeqIO.read(tmp_path / "step00_plasmid.gb", "genbank")
    insert = SeqIO.read(tmp_path / "step01_insert.gb", "genbank")
    assembled = SeqIO.read(tmp_path / "step01_plasmid.gb", "genbank")
    assert step00.annotations["topology"] == "circular"
    assert insert.annotations["topology"] == "linear"
    assert assembled.annotations["topology"] == "circular"
    assert any(feature.qualifiers.get("protected") == ["true"] for feature in step00.features)
    assert any(feature.qualifiers.get("feature_kind") == ["restriction_site"] for feature in insert.features)
    cds = next(
        feature for feature in assembled.features
        if feature.type == "CDS" and "repeat-protein CDS" in feature.qualifiers.get("label", [""])[0]
    )
    coding_dna = str(cds.extract(assembled.seq))
    assert translate_dna(coding_dna) == cds.qualifiers["translation"][0]
    assert cds.qualifiers["translation"][0] == result.final_protein_sequence
    assert str(assembled.seq) == result.final_plasmid["final_plasmid_sequence"]
    assert hashlib.sha256(str(assembled.seq).encode()).hexdigest() == result.final_plasmid["final_plasmid_sha256"]
    assert (tmp_path / "final_plasmid.gb").read_bytes() == (tmp_path / "step01_plasmid.gb").read_bytes()
    manifest = json.loads((tmp_path / "assembly_step_manifest.json").read_text())
    assert manifest[-1]["matches_final_plasmid"] is True
    assert result.assembly_steps == manifest
    assert paths["assembly_step_manifest_json"].endswith("assembly_step_manifest.json")
    assert (tmp_path / "maps" / "step00_plasmid.circular.png").stat().st_size > 0
    assert (tmp_path / "maps" / "step00_plasmid.linear.svg").stat().st_size > 0
    assert (tmp_path / "maps" / "step01_insert.linear.png").stat().st_size > 0

    candidate = next(
        row for row in result.protein_candidates
        if row["candidate_id"] == result.selected_route["candidate_id"]
    )
    preview, preview_row = build_step00_plasmid_record(
        result.selected_route,
        candidate,
        result.request["selection"]["site_iii_enzyme"],
    )
    assert str(preview.seq) == str(step00.seq)
    assert preview_row["sequence_sha256"] == hashlib.sha256(str(step00.seq).encode()).hexdigest()
    assert preview_row["site_audit"] == result.assembly_steps[0]["site_audit"]


def test_checkpoint_never_emits_fake_accepted_fasta_and_final_zip_is_timestamped(tmp_path):
    checkpoint = tmp_path / "checkpoint.zip"
    write_secondary_checkpoint(
        {"event": "heartbeat", "sequence_id": "sample", "failure_reason": "no accepted candidate"},
        checkpoint,
    )
    with zipfile.ZipFile(checkpoint) as archive:
        assert archive.namelist() == ["checkpoint.json"]

    payload = {
        "event": "accepted_secondary",
        "sequence_id": "sample id",
        "repeat_copies": 4,
        "core_sequence": "GCT" * 4,
        "purchase_sequence": "AA" + "GCT" * 4 + "TT",
        "purchase_sha256": hashlib.sha256(("AA" + "GCT" * 4 + "TT").encode()).hexdigest(),
        "idt_complexity_score": 2.5,
        "validation_mode": "api",
        "ga_weights": {"repeated_re_site": 12.0},
    }
    write_secondary_checkpoint(payload, checkpoint)
    with zipfile.ZipFile(checkpoint) as archive:
        assert set(archive.namelist()) == {
            "checkpoint.json", "best_secondary_core.fasta", "best_secondary_purchase.fasta"
        }
        public = json.loads(archive.read("checkpoint.json"))
        assert "core_sequence" not in public and "purchase_sequence" not in public
        assert public["accepted_secondary_available"] is True

    output = tmp_path / "results"
    output.mkdir()
    (output / "proof.txt").write_text("verified")
    archive = timestamped_results_archive(
        output,
        tmp_path / "archives",
        sequence_id="sample id",
        timestamp=datetime(2026, 8, 12, 12, 34, 56, tzinfo=timezone.utc),
    )
    assert archive.name == "hurdler_sample_id_20260812T123456Z_results.zip"
    with zipfile.ZipFile(archive) as handle:
        assert handle.read("proof.txt") == b"verified"
