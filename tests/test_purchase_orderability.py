from __future__ import annotations

import hashlib

from hurdler.purchase_orderability import (
    classify_double_stranded_purchase,
    classify_existing_primer_pair,
)


class _Scorer:
    def __init__(self, score: float):
        self.score_value = score
        self.calls: list[tuple[str, str]] = []

    def score(self, name: str, sequence: str):
        self.calls.append((name, sequence))
        passed = self.score_value < 10
        return {
            "idt_status": "passed" if passed else "failed",
            "idt_explicit_pass": passed,
            "idt_complexity_score": self.score_value,
            "idt_response_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        }


def test_exact_ds_purchase_uses_oligo_pairs_below_gblock_minimum():
    scorer = _Scorer(0)
    standard = classify_double_stranded_purchase("short", "ACGT" * 10, idt_scorer=scorer)
    ultramer = classify_double_stranded_purchase("long", "ACGT" * 30, idt_scorer=scorer)
    assert standard["product_class"] == "complementary_standard_primer_pair"
    assert ultramer["product_class"] == "complementary_ultramer_pair"
    assert standard["orderable"] and ultramer["orderable"]
    assert scorer.calls == []


def test_gblock_requires_live_idt_score_strictly_below_ten():
    sequence = "ACGT" * 55
    accepted = classify_double_stranded_purchase(
        "accepted", sequence, idt_scorer=_Scorer(9.9)
    )
    rejected = classify_double_stranded_purchase(
        "rejected", sequence, idt_scorer=_Scorer(10.0)
    )
    assert accepted["product_class"] == "gblock"
    assert accepted["orderable"]
    assert not rejected["orderable"]
    assert rejected["failure_reason"] == "gblock_not_idt_accepted"


def test_sticky_end_pair_checks_the_actual_two_oligos():
    accepted = classify_existing_primer_pair(
        "donor", "A" * 25, "C" * 64, core_sequence="ACGT" * 8
    )
    rejected = classify_existing_primer_pair(
        "donor", "A" * 19, "C" * 64, core_sequence="ACGT" * 8
    )
    assert accepted["orderable"]
    assert accepted["maximum_single_oligo_length_nt"] == 64
    assert not rejected["orderable"]
    assert rejected["failure_reason"] == "oligo_length_out_of_range"
