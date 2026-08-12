import pytest

from hurdler.optimization import (
    _new_duplicate_windows,
    _new_duplicate_windows_from_seen,
    _new_site_occurrence_count,
    _site_occurrences,
    _update_seen_windows,
    _locked_site_excess_lower_bound,
    _maximum_verified_construct,
    diversify_codons,
    repeated_nmer_count,
    translate_dna,
)
from hurdler.rate import legacy_random_modules


def test_legacy_random_sequence_order_is_frozen():
    iterator = legacy_random_modules(min_length=7, max_length=7, tests_per_plasmid=1, seed=42)
    assert next(iterator) == (7, "pGEX-4T-1", 0, "EAKIIFE")
    assert next(iterator) == (7, "pMAL-c5X", 0, "VDWQCAD")


def test_codon_diversification_preserves_translation():
    protein = "NEQIQAVIDAGALPALVQLLSSP"
    dna = diversify_codons(protein)
    assert len(dna) == len(protein) * 3
    assert translate_dna(dna) == protein


def test_codon_diversification_avoids_selected_recognition_sites():
    dna = diversify_codons("AAAA", site_limits={"GCAGCA": 0})
    assert "GCAGCA" not in dna
    assert translate_dna(dna) == "AAAA"


def test_locked_window_is_translation_checked():
    with pytest.raises(ValueError):
        diversify_codons("AAA", {0: "TGGTGGTGG"})


def test_repeat_count_counts_only_duplicates():
    assert repeated_nmer_count("AAAAAAAA", 4) == 4
    assert repeated_nmer_count("ACGT", 4) == 0


def test_incremental_repeat_windows_equal_full_prefix_scan():
    sequence = ""
    seen = {n: {} for n in (8, 13, 14)}
    for codon in ("GCT", "GCC", "GCA", "GCG") * 8:
        sequence += codon
        for n in seen:
            assert _new_duplicate_windows_from_seen(sequence, n, seen[n]) == _new_duplicate_windows(
                sequence, n
            )
            _update_seen_windows(sequence, n, seen[n])


def test_incremental_site_occurrences_equal_full_scan():
    sequence = ""
    counts = {"GGTCTC": 0, "GAATTC": 0, "AGGAG": 0}
    for codon in ("GGT", "CTC", "GAA", "TTC", "CTA", "GGA", "GAG") * 3:
        sequence += codon
        for site in counts:
            counts[site] += _new_site_occurrence_count(sequence, site)
            assert counts[site] == len(_site_occurrences(sequence, site))


def test_maximum_construct_falls_back_to_shorter_copy_count(monkeypatch):
    candidate = {
        "plasmid": "pGEX-4T-1",
        "site_i_position": 0,
        "site_ii_position": 5,
        "site_i_enzyme": "A",
        "site_ii_enzyme": "B",
    }

    def fake_metrics(unit, copies, solution, weights):
        if copies > 7:
            raise ValueError("synthetic length failure")
        return {"dna_sequence": "GCT" * (len(unit) * copies)}

    monkeypatch.setattr("hurdler.optimization._construct_metrics", fake_metrics)
    copies, solution, metrics, errors = _maximum_verified_construct("ACDEFG", 20, [candidate], {})
    assert copies == 7
    assert solution == candidate
    assert metrics is not None
    assert errors


def test_maximum_construct_deduplicates_equivalent_optimization_signatures(monkeypatch):
    base = {
        "plasmid": "pGEX-4T-1",
        "site_i_position": 0,
        "site_ii_position": 5,
        "site_i_9mer_bp": "GCTGCTGCT",
        "site_ii_9mer_bp_mutated": "GCTGCTGCC",
        "site_i_recognition_site": "GGTCTC",
        "site_ii_recognition_site": "GAGACC",
        "site_iii_sites": "",
    }
    candidates = [
        {**base, "site_i_enzyme": "A", "site_ii_enzyme": "B"},
        {**base, "site_i_enzyme": "C", "site_ii_enzyme": "D"},
    ]
    calls = 0

    def always_fails(unit, copies, solution, weights):
        nonlocal calls
        calls += 1
        raise ValueError("synthetic failure")

    monkeypatch.setattr("hurdler.optimization._construct_metrics", always_fails)
    copies, solution, metrics, errors = _maximum_verified_construct("ACDEFG", 10, candidates, {})
    assert copies == 0
    assert solution is None and metrics is None and errors
    # The second bookkeeping-only duplicate is removed before evaluation.
    # Full length plus the same three binary-search points run only once.
    assert calls == 4


def test_locked_site_lower_bound_prunes_only_unavoidable_sites(monkeypatch):
    candidate = {
        "plasmid": "pGEX-4T-1",
        "site_i_position": 0,
        "site_ii_position": 5,
        "site_i_enzyme": "A",
        "site_ii_enzyme": "B",
        "site_i_9mer_bp": "GGTCTCGCT",
        "site_ii_9mer_bp_mutated": "GCTGCTGCT",
        "site_i_recognition_site": "GGTCTC",
        "site_ii_recognition_site": "GAGACC",
        "site_iii_sites": "GCTGCT",
    }
    assert _locked_site_excess_lower_bound(candidate) > 0

    def must_not_run(*args, **kwargs):
        raise AssertionError("an impossible locked-window candidate was evaluated")

    monkeypatch.setattr("hurdler.optimization._construct_metrics", must_not_run)
    copies, solution, metrics, errors = _maximum_verified_construct("ALAALA", 10, [candidate], {})
    assert copies == 0 and solution is None and metrics is None
    assert errors and "unavoidable" in errors[0]
