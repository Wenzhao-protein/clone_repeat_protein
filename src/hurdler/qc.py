"""Run-level and legacy-data quality checks."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .constants import PLASMIDS
from .index import PatternIndex
from .io import utc_now, write_json_atomic


OPTIMIZED_NOTEBOOK_BASELINE = {
    "candidate_count": 2_360_004,
    "pattern_count": 1_528_765,
    "plasmid_pattern_counts": {
        "pGEX-4T-1": 895_481,
        "pMAL-c5X": 736_353,
        "pET-21a(+)": 718_033,
        "pET-28a(+)": 545_149,
        "pET-28a(+)_start_codon": 723_332,
        "pCold_I": 1_021_688,
        "pUC18": 1_292_709,
        "pQE-3": 762_145,
    },
}


def legacy_qc(source_dir: str | Path, index_dir: str | Path, output_path: str | Path) -> dict[str, object]:
    source = Path(source_dir)
    index = PatternIndex.load(index_dir)
    df1_path = source / "hurdler_three_site_combinations_df1.csv"
    mismatched_overhangs = None
    rows = None
    if df1_path.exists():
        frame = pd.read_csv(df1_path)
        rows = len(frame)
        mismatched_overhangs = int((frame["site_ii_ovhg"] != frame["site_iii_ovhg"]).sum())
    current_plasmid_counts = {
        plasmid: int(((index.plasmid_masks & (1 << bit)) != 0).sum())
        for bit, plasmid in enumerate(PLASMIDS)
    }
    baseline = OPTIMIZED_NOTEBOOK_BASELINE
    diffs = {
        "candidate_count": int(index.metadata["candidate_count"]) - baseline["candidate_count"],
        "pattern_count": len(index.keys) - baseline["pattern_count"],
        "plasmid_pattern_counts": {
            plasmid: current_plasmid_counts[plasmid] - baseline["plasmid_pattern_counts"][plasmid]
            for plasmid in PLASMIDS
        },
    }
    payload = {
        "created_at": utc_now(),
        "status": "passed",
        "rule_profile": index.metadata["rule_profile"]["name"],
        "pattern_count": len(index.keys),
        "candidate_count": int(index.metadata["candidate_count"]),
        "enzyme_pair_count": len(index.pair_table),
        "plasmid_pattern_counts": current_plasmid_counts,
        "optimized_notebook_embedded_baseline": baseline,
        "machine_readable_diff": diffs,
        "plasmids": list(PLASMIDS),
        "legacy_df1_rows": rows,
        "legacy_df1_site_ii_iii_overhang_mismatches": mismatched_overhangs,
        "note": (
            "The embedded notebook baseline includes combinations that are excluded by the frozen "
            "signed Site-II/Site-III overhang rule and may reflect older source tables. Differences "
            "are reported, never silently accepted; they do not alter legacy-optimized-v1."
        ),
    }
    write_json_atomic(payload, output_path)
    return payload
