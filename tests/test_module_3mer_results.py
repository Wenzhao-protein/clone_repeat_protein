import json

import pandas as pd
import pytest

from hurdler.constants import PLASMIDS
from hurdler.module_3mer_results import build_module_3mer_results


def _fixture_index(tmp_path):
    root = tmp_path / "index"
    root.mkdir()
    (root / "metadata.json").write_text(
        json.dumps({"rule_profile": {"name": "legacy-optimized-v1"}})
    )
    pd.DataFrame(
        [{"site_i_enzyme": "EnzI", "site_i_3mer_aa": "EFG"}]
    ).to_parquet(root / "site_i_variants.parquet", index=False)
    pd.DataFrame(
        [{"site_ii_enzyme": "EnzII", "site_ii_3mer_aa": "ACD"}]
    ).to_parquet(root / "site_ii_variants.parquet", index=False)
    pd.DataFrame(
        [
            {
                "site_i_enzyme": "EnzI",
                "site_ii_enzyme": "EnzII",
                "plasmid_mask": 1 << PLASMIDS.index("pUC18"),
            }
        ]
    ).to_parquet(root / "enzyme_pairs.parquet", index=False)
    return root


def _source():
    return pd.DataFrame(
        [
            {
                "module_id": "natural_1abc_A_1_6_fixture",
                "display_name": "1abc",
                "source_accession": "1abc",
                "collection": "Natural",
                "family": "fixture",
                "middle_module_sequence_aa": "ACDEFG",
                "middle_module_length_aa": 6,
                "hurdler_compatible": True,
                "selected_plasmid": "pUC18",
                "selected_site_i_enzyme": "EnzI",
                "selected_site_ii_enzyme": "EnzII",
                "selected_direction": "left",
                "selected_site_i_position": 3,
                "selected_site_ii_position": 0,
                "corpus_version": "expanded-middle-repeatsdb-foldseek-v1",
                "hurdler_rules_version": "legacy-optimized-v1",
            },
            {
                "module_id": "designed_fixture",
                "display_name": "fixture",
                "source_accession": "fixture",
                "collection": "Designed",
                "family": "DHR",
                "middle_module_sequence_aa": "AAAAAA",
                "middle_module_length_aa": 6,
                "hurdler_compatible": False,
                "selected_plasmid": pd.NA,
                "selected_site_i_enzyme": pd.NA,
                "selected_site_ii_enzyme": pd.NA,
                "selected_direction": pd.NA,
                "selected_site_i_position": pd.NA,
                "selected_site_ii_position": pd.NA,
                "corpus_version": "expanded-middle-repeatsdb-foldseek-v1",
                "hurdler_rules_version": "legacy-optimized-v1",
            },
        ]
    )


def test_build_module_3mer_results_recovers_both_selected_windows(tmp_path):
    result = build_module_3mer_results(
        _source(), _fixture_index(tmp_path), expected_counts=None
    ).set_index("sequence_id")

    compatible = result.loc["natural_1abc_A_1_6_fixture"]
    assert compatible["selected_module_sequence_aa"] == "ACDEFG"
    assert compatible["selected_re_pair"] == "EnzI / EnzII"
    assert compatible["site_i_3mer_aa"] == "EFG"
    assert compatible["site_ii_3mer_aa"] == "ACD"
    assert compatible["three_mer_aa_pair"] == "EFG / ACD"
    assert compatible["direction"] == "left"

    incompatible = result.loc["designed_fixture"]
    assert not incompatible["hurdler_compatible"]
    assert incompatible["selected_re_pair"] == ""
    assert incompatible["three_mer_aa_pair"] == ""


def test_build_module_3mer_results_rejects_lookup_mismatch(tmp_path):
    source = _source()
    source.loc[0, "selected_site_i_enzyme"] = "WrongEnzyme"
    with pytest.raises(ValueError, match="absent from lookup"):
        build_module_3mer_results(source, _fixture_index(tmp_path))


def test_build_module_3mer_results_checks_production_counts(tmp_path):
    with pytest.raises(ValueError, match="Unexpected production collection counts"):
        build_module_3mer_results(
            _source(),
            _fixture_index(tmp_path),
            expected_counts={"Natural": 25_913, "Designed": 182},
        )
