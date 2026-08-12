import hashlib
import json

import pandas as pd
import pytest

from hurdler.module_results import export_module_results


def _write_fixture_inputs(tmp_path):
    dna = "GAATTCGAGTTC"  # EFEF
    dna_sha = hashlib.sha256(dna.encode()).hexdigest()
    response_sha = "a" * 64
    catalog = pd.DataFrame(
        [
            {
                "module_id": "natural_fixture",
                "collection": "natural_all",
                "module_type": "natural",
                "family": "RepeatsDB 4.1",
                "unit_sequence": "EF",
                "unit_length": 2,
                "full_sequence": "EFEFEF",
                "source_name": "=unsafe label",
                "source_accession": "P12345",
                "source_url": "https://repeatsdb.org/structure/P12345",
                "structure_accession": "P12345",
                "uniprot_accessions": "['P12345']",
                "unit_start": 3,
                "unit_end": 4,
                "selected_module_index": 1,
                "selected_module_count": 3,
                "repeat_region_start": 1,
                "repeat_region_end": 6,
                "repeat_count": 3,
                "boundary_method": "repeatsdb-direct-middle-unit-v1",
                "boundary_refinement_status": "source_annotation_middle_unit",
                "corpus_version": "expanded-middle-repeatsdb-foldseek-v1",
            },
            {
                "module_id": "designed_fixture",
                "collection": "designed_all",
                "module_type": "designed",
                "family": "DHR",
                "unit_sequence": "ACDEFG",
                "unit_length": 6,
                "full_sequence": "ACDEFGACDEFG",
                "source_name": "DHR supplement",
                "source_accession": "DHRX",
                "source_url": "https://example.org/dhrx",
                "unit_start": 1,
                "unit_end": 6,
                "selected_module_index": 0,
                "selected_module_count": 2,
                "repeat_region_start": 1,
                "repeat_region_end": 12,
                "repeat_count": 2,
                "boundary_method": "biotite-dssp-foldseek-strict-dual-evidence-v1",
                "boundary_refinement_status": "strict_dual_evidence_passed",
                "corpus_version": "expanded-middle-repeatsdb-foldseek-v1",
            },
        ]
    )
    compatibility = pd.DataFrame(
        [
            {
                "module_id": "natural_fixture",
                "collection": "Natural",
                "unit_sequence": "EF",
                "hurdler_compatible": True,
                "compatible_plasmid_count": 1,
                "compatible_plasmids_json": '["pUC18"]',
                "candidate_solution_count": 1,
                "selected_solution_json": json.dumps(
                    {
                        "site_iii_enzymes": "BsaI,BsmBI",
                        "site_iii_sites": "GGTCTC,CGTCTC",
                        "site_i_ovhg": -4,
                        "site_ii_ovhg": -4,
                    }
                ),
                "selected_plasmid": "pUC18",
                "selected_site_i_enzyme": "EcoRI",
                "selected_site_ii_enzyme": "BsaI",
                "selected_site_i_recognition_site": "GAATTC",
                "selected_site_ii_recognition_site": "GGTCTC",
                "selected_direction": "left",
                "selected_site_i_position": 0,
                "selected_site_ii_position": 3,
                "rules_version": "legacy-optimized-v1",
            },
            {
                "module_id": "designed_fixture",
                "collection": "Designed",
                "unit_sequence": "ACDEFG",
                "hurdler_compatible": False,
                "compatible_plasmid_count": 0,
                "compatible_plasmids_json": "[]",
                "candidate_solution_count": 0,
                "selected_solution_json": "",
            },
        ]
    )
    maximum_rows = []
    for capacity in (1800, 3000):
        maximum_rows.append(
            {
                "module_id": "natural_fixture",
                "collection": "Natural",
                "unit_sequence": "EF",
                "fragment_limit_bp": capacity,
                "mathematical_max_copies": capacity // 6,
                "verified_max_copies": 2,
                "adaptive_maximum_proof_status": "next_copy_failed_at_100",
                "adaptive_stop_reason": "copy_3_failed_at_100",
                "optimization_status": "passed_adaptive_reduced",
                "failure_reason": "",
                "final_ga_weights_json": '{"repeated_8mer":80}',
                "idt_status": "passed",
                "idt_complexity_score": 5.0,
                "idt_score_policy": "idt-rule-score-sum-lt10-v1",
                "idt_positive_score_names_json": '["Repeat Length"]',
                "idt_rule_details_json": json.dumps(
                    [
                        {
                            "name": "Repeat Length",
                            "score": 5.0,
                            "is_violated": True,
                            "actual_value": 14,
                            "threshold_value": 13,
                        }
                    ]
                ),
                "idt_response_sha256": response_sha,
                "idt_scored_sequence_sha256": dna_sha,
                "selected_pair_re_site_excess": 0,
                "dna_sequence": dna,
                "dna_length": len(dna),
                "validation_passed": True,
                "validation_reasons_json": "[]",
            }
        )

    paths = {
        "catalog": tmp_path / "catalog.parquet",
        "mapping": tmp_path / "mappings.parquet",
        "compatibility": tmp_path / "compatibility.parquet",
        "maximum": tmp_path / "maximum.parquet",
    }
    catalog.to_parquet(paths["catalog"], index=False)
    catalog[[
        "module_id",
        "source_name",
        "source_accession",
        "source_url",
        "structure_accession",
        "uniprot_accessions",
    ]].to_parquet(paths["mapping"], index=False)
    compatibility.to_parquet(paths["compatibility"], index=False)
    pd.DataFrame(maximum_rows).to_parquet(paths["maximum"], index=False)
    return paths, pd.DataFrame(maximum_rows), dna


def test_export_module_results_is_one_row_per_module_and_hides_unaccepted_dna(tmp_path):
    paths, _, dna = _write_fixture_inputs(tmp_path)
    destination = tmp_path / "data/results/natural_designed_repeat_protein_hurdler_idt.csv"
    exported = export_module_results(
        paths["catalog"],
        [paths["mapping"]],
        paths["compatibility"],
        paths["maximum"],
        destination,
        repository_root=tmp_path,
        generated_at_utc="2026-08-11T00:00:00+00:00",
    )

    assert len(exported) == 2
    assert exported.module_id.is_unique
    natural = exported.set_index("module_id").loc["natural_fixture"]
    designed = exported.set_index("module_id").loc["designed_fixture"]
    assert natural.cap1800_idt_accepted_dna == dna
    assert natural.cap3000_idt_accepted_dna == dna
    assert natural.cap1800_maximum_verified_copies == 2
    assert natural.record_status == "idt_accepted_both_capacities"
    assert designed.record_status == "hurdler_incompatible"
    assert designed.cap1800_idt_accepted_dna == ""
    assert natural.source_name == "'=unsafe label"
    assert "P12345" in natural.search_terms
    assert "EcoRI" in natural.search_terms
    reread = pd.read_csv(destination, keep_default_na=False)
    assert len(reread) == 2
    assert not reread.to_csv(index=False).count("/net/scratch")


def test_export_module_results_rejects_missing_capacity(tmp_path):
    paths, maximum, _ = _write_fixture_inputs(tmp_path)
    maximum.loc[maximum.fragment_limit_bp.eq(1800)].to_parquet(
        paths["maximum"], index=False
    )
    with pytest.raises(ValueError, match="Maximum-copy table is incomplete"):
        export_module_results(
            paths["catalog"],
            [paths["mapping"]],
            paths["compatibility"],
            paths["maximum"],
            tmp_path / "result.csv",
            repository_root=tmp_path,
        )
