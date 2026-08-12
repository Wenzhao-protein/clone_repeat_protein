from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_designed_citation_guide_accounts_for_all_182_modules():
    with (ROOT / "data/results/designed_repeat_protein_citations.csv").open(
        newline=""
    ) as handle:
        citations = list(csv.DictReader(handle))
    with (ROOT / "data/results/natural_designed_repeat_protein_hurdler_idt.csv").open(
        newline=""
    ) as handle:
        designed = [
            row
            for row in csv.DictReader(handle)
            if row["collection"].strip().lower() == "designed"
        ]

    assert len(designed) == 182
    assert len(citations) == 19
    assert sum(int(row["designed_module_rows"]) for row in citations) == len(designed)
    assert sum(row["publication_status"] == "peer_reviewed" for row in citations) == 17
    assert sum(row["publication_status"] == "pdb_only_to_be_published" for row in citations) == 2
    assert all(row["doi_or_record_url"].startswith("https://") for row in citations)
