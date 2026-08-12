"""Deterministic legacy-compatible Monte Carlo success rates."""

from __future__ import annotations

import random
from pathlib import Path

import pandas as pd

from .constants import AMINO_ACIDS, PLASMIDS
from .index import PatternIndex
from .matching import match_module
from .rules import RuleProfile, LEGACY_OPTIMIZED_V1


def legacy_random_modules(
    min_length: int = 7,
    max_length: int = 60,
    tests_per_plasmid: int = 1000,
    seed: int = 42,
):
    """Yield modules in the exact optimized-notebook loop order."""
    generator = random.Random(seed)
    for length in range(min_length, max_length + 1):
        for plasmid in PLASMIDS:
            for test_index in range(tests_per_plasmid):
                module = "".join(generator.choice(AMINO_ACIDS) for _ in range(length))
                yield length, plasmid, test_index, module


def run_success_rate(
    index_dir: str | Path,
    output_path: str | Path,
    *,
    min_length: int = 7,
    max_length: int = 60,
    tests_per_plasmid: int = 1000,
    seed: int = 42,
    rules: RuleProfile = LEGACY_OPTIMIZED_V1,
) -> pd.DataFrame:
    index = PatternIndex.load(index_dir)
    counts = {(length, plasmid): 0 for length in range(min_length, max_length + 1) for plasmid in PLASMIDS}
    for length, plasmid, _test_index, module in legacy_random_modules(min_length, max_length, tests_per_plasmid, seed):
        if match_module(module, plasmid, index, rules=rules, expand_short=False).success:
            counts[(length, plasmid)] += 1
    rows = [
        {
            "module_length": length,
            "plasmid": plasmid,
            "successes": counts[(length, plasmid)],
            "tests": tests_per_plasmid,
            "success_rate": counts[(length, plasmid)] / tests_per_plasmid,
            "method": "monte_carlo_legacy",
            "seed": seed,
        }
        for length in range(min_length, max_length + 1)
        for plasmid in PLASMIDS
    ]
    frame = pd.DataFrame(rows)
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.suffix == ".parquet":
        frame.to_parquet(destination, index=False)
    else:
        frame.to_csv(destination, index=False)
    return frame
