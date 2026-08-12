import pytest

from hurdler.index import decode_3mer, decode_pattern, encode_3mer, encode_pattern
from hurdler.rules import LEGACY_OPTIMIZED_V1


@pytest.mark.parametrize("three_mer", ["AAA", "ACD", "NEQ", "YYY"])
def test_3mer_round_trip(three_mer):
    assert decode_3mer(encode_3mer(three_mer)) == three_mer


@pytest.mark.parametrize("direction", ["left", "right"])
def test_pattern_round_trip(direction):
    key = encode_pattern("NEQ", "IQA", direction)
    assert decode_pattern(key) == ("NEQ", "IQA", direction)


def test_legacy_distance_boundary_is_inclusive_at_five():
    assert LEGACY_OPTIMIZED_V1.distance_is_valid(5, 6)
    assert not LEGACY_OPTIMIZED_V1.distance_is_valid(4, 6)
    assert not LEGACY_OPTIMIZED_V1.distance_is_valid(6, 6)


def test_invalid_3mer_is_rejected():
    with pytest.raises(ValueError):
        encode_3mer("ABZ")
