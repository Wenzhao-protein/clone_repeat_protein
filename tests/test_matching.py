import numpy as np
import pandas as pd

from hurdler.constants import PLASMIDS
from hurdler.index import PatternIndex, decode_pattern, encode_pattern
from hurdler.matching import expand_short_module, match_module, materialize_best_solution, query_all_plasmids
from hurdler.short_screen import SHORT_MOTIF_COUNTS, SHORT_MOTIF_TOTAL, iter_motifs


def make_index(key, plasmid_bit=0):
    pair_table = pd.DataFrame(
        [
            {
                "pair_id": 0,
                "site_i_enzyme": "EcoRI",
                "site_ii_enzyme": "HindIII",
                "site_i_ovhg": -4,
                "site_ii_ovhg": -4,
                "site_iii_enzymes": "BsaI",
                "orthogonality": 2.0,
                "plasmid_mask": 1 << plasmid_bit,
            }
        ]
    )
    site_i = pd.DataFrame(
        [{"site_i_enzyme": "EcoRI", "site_i_3mer_aa": "AAA", "site_i_9mer_bp": "GCTGCTGCT"}]
    )
    site_ii = pd.DataFrame(
        [
            {
                "site_ii_enzyme": "HindIII",
                "site_ii_3mer_aa": "CCC",
                "site_ii_9mer_bp_original": "TGTTGTTGT",
                "site_ii_9mer_bp_mutated": "TGCTGCTGC",
            }
        ]
    )
    return PatternIndex(
        keys=np.array([key], dtype=np.uint32),
        plasmid_masks=np.array([1 << plasmid_bit], dtype=np.uint8),
        solution_counts=np.array([[1 if bit == plasmid_bit else 0 for bit in range(8)]], dtype=np.uint16),
        best_pair_ids=np.array([[0 if bit == plasmid_bit else 65535 for bit in range(8)]], dtype=np.uint16),
        pair_table=pair_table,
        site_i_table=site_i,
        site_ii_table=site_ii,
        metadata={"schema_version": 1},
    )


def test_short_module_expands_to_minimum_six():
    assert expand_short_module("VLA") == ("VLAVLA", 2)
    assert expand_short_module("VL") == ("VLVLVL", 3)
    assert expand_short_module("A") == ("AAAAAA", 6)


def test_short_motif_space_is_exponential_and_not_deduplicated():
    assert SHORT_MOTIF_COUNTS == {1: 20, 2: 400, 3: 8_000, 4: 160_000, 5: 3_200_000}
    assert SHORT_MOTIF_TOTAL == 3_368_420
    assert len(list(iter_motifs(3))) == 20**3
    assert len(list(iter_motifs(5, "AC"))) == 20**3


def test_right_match_uses_inclusive_distance_five():
    index = make_index(encode_pattern("AAA", "CCC", "right"))
    result = match_module("AAAXXCCC".replace("X", "G"), PLASMIDS[0], index, expand_short=False)
    assert result.success
    assert result.site_i_position == 0
    assert result.site_ii_position == 5
    assert materialize_best_solution(result, index)["site_iii_enzymes"] == "BsaI"


def test_left_match_reverses_site_roles():
    index = make_index(encode_pattern("AAA", "CCC", "left"))
    result = match_module("CCCGGAAA", PLASMIDS[0], index, expand_short=False)
    assert result.success
    assert result.site_i_position == 5
    assert result.site_ii_position == 0


def test_plasmid_mask_is_enforced():
    index = make_index(encode_pattern("AAA", "CCC", "right"), plasmid_bit=1)
    assert not match_module("AAAGGCCC", PLASMIDS[0], index, expand_short=False).success
    assert match_module("AAAGGCCC", PLASMIDS[1], index, expand_short=False).success


def test_batched_query_matches_independent_queries():
    index = make_index(encode_pattern("AAA", "CCC", "right"), plasmid_bit=1)
    batched = query_all_plasmids("AAAGGCCC", index, expand_short=False)
    independent = [match_module("AAAGGCCC", plasmid, index, expand_short=False) for plasmid in PLASMIDS]
    assert [result.to_dict() for result in batched] == [result.to_dict() for result in independent]


def test_sparse_index_matches_brute_force_golden_random_cases():
    """Property-style golden check against an implementation without searchsorted."""
    rng = np.random.default_rng(42)
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    modules = [
        "".join(rng.choice(alphabet, size=int(rng.integers(6, 31))))
        for _case in range(80)
    ]
    forced_keys = {
        encode_pattern(module[:3], (module * 2)[5:8], direction)
        for module in modules
        for direction in ("left", "right")
    }
    random_keys = set(int(value) for value in rng.choice(20**6 * 2, size=300, replace=False))
    keys = np.array(sorted(forced_keys | random_keys), dtype=np.uint32)
    masks = rng.integers(1, 256, size=len(keys), dtype=np.uint8)
    index = PatternIndex(
        keys=keys,
        plasmid_masks=masks,
        solution_counts=np.ones((len(keys), 8), dtype=np.uint16),
        best_pair_ids=np.zeros((len(keys), 8), dtype=np.uint16),
        pair_table=pd.DataFrame(),
        site_i_table=pd.DataFrame(),
        site_ii_table=pd.DataFrame(),
        metadata={"schema_version": 1},
    )
    brute = {int(key): int(mask) for key, mask in zip(keys, masks)}

    for module in modules:
        sequence = module * 2
        module_length = len(module)
        for plasmid_bit, plasmid in enumerate(PLASMIDS):
            expected = None
            for left_position in range(len(sequence) - 2):
                left = sequence[left_position : left_position + 3]
                maximum_right = min(len(sequence) - 3, left_position + module_length - 1)
                for right_position in range(left_position + 5, maximum_right + 1):
                    right = sequence[right_position : right_position + 3]
                    candidates = (
                        (left, right, "right", left_position, right_position),
                        (right, left, "left", right_position, left_position),
                    )
                    for site_i, site_ii, direction, site_i_position, site_ii_position in candidates:
                        key = encode_pattern(site_i, site_ii, direction)
                        if brute.get(key, 0) & (1 << plasmid_bit):
                            expected = (key, site_i_position, site_ii_position, direction)
                            break
                    if expected is not None:
                        break
                if expected is not None:
                    break
            observed = match_module(module, plasmid, index, expand_short=False)
            assert observed.success == (expected is not None)
            if expected is not None:
                assert (
                    observed.pattern_key,
                    observed.site_i_position,
                    observed.site_ii_position,
                    observed.direction,
                ) == expected
                assert decode_pattern(observed.pattern_key)[2] == observed.direction
