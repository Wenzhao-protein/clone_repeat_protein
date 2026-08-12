#!/usr/bin/env python3
"""Run the historical 1--60 AA HURDLER success landscape with 16 workers.

The exhaustive 1--5 AA screen is written to one Parquet file.  The seeded
6--60 AA Monte Carlo observations are written to a second Parquet file.  No
worker writes a shard: the parent process is the sole Parquet writer.

This workflow deliberately reconstructs the plasmid-specific pattern sets
used by ``hurdler_success_rate_analysis.ipynb``.  The scan-copy count is an
explicit 2-or-3 parameter; both modes preserve the notebook's exact
``re.search()`` semantics and distance criterion (3-mer start distance
``5 <= d < module_length``).  It does not use the later package-wide sparse
lookup, whose enzyme/pattern population differs from the historical
success-landscape notebook.
"""

from __future__ import annotations

import argparse
import bisect
import itertools
import json
import multiprocessing as mp
import os
import random
import sys
import time
from pathlib import Path
from typing import Iterable, Iterator, NamedTuple, Sequence

import duckdb
import matplotlib
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from hurdler.constants import AA_TO_INT, AMINO_ACIDS, PLASMIDS, THREE_MER_SPACE
from hurdler.io import sha256_file
from hurdler.matching import expand_short_module


DEFAULT_REPO = Path(__file__).resolve().parents[1]
DEFAULT_SHORT_OUTPUT = Path(
    "/net/scratch/wendai/projects/hurdler/clone_repeat_protein/studies/"
    "hurdler_validation/step02_success_landscape/runs/"
    "run06_three_copy_16core/raw/short_motifs_1_5.parquet"
)
DEFAULT_RANDOM_OUTPUT = DEFAULT_SHORT_OUTPUT.with_name("random_modules_6_60.parquet")
DEFAULT_FIGURE_DIR = Path(
    "/home/wendai/projects/hurdler/clone_repeat_protein/studies/"
    "hurdler_validation/step02_success_landscape/figures/scan_3x"
)
SHORT_COUNTS = {length: 20**length for length in range(1, 6)}
SHORT_TOTAL = sum(SHORT_COUNTS.values())
SIX_AA_SEED = 420006
LEGACY_SEED = 42
SUCCESS_RULE_PROFILE = "historical-notebook-success-v1"


class HistoricalPatternIndex(NamedTuple):
    """Compact historical ``left_3mer -> right_3mer`` plasmid masks."""

    masks: np.ndarray
    three_mer_present: np.ndarray
    source_hashes: dict[str, str]


_HISTORICAL_INDEX: HistoricalPatternIndex | None = None
_SCAN_COPIES = 2


def _short_schema() -> pa.Schema:
    fields = [
        pa.field("motif", pa.string(), nullable=False),
        pa.field("motif_length", pa.int8(), nullable=False),
        pa.field("expansion_copies", pa.int8(), nullable=False),
        pa.field("expanded_module", pa.string(), nullable=False),
        pa.field("expanded_length", pa.int8(), nullable=False),
    ]
    for plasmid in PLASMIDS:
        fields.append(pa.field(f"{plasmid}_success", pa.bool_(), nullable=False))
    fields.extend(
        [
            pa.field("success_mask", pa.uint8(), nullable=False),
            pa.field("successful_plasmids", pa.int8(), nullable=False),
            pa.field("any_success", pa.bool_(), nullable=False),
        ]
    )
    return pa.schema(fields)


RANDOM_SCHEMA = pa.schema(
    [
        pa.field("module_length", pa.int16(), nullable=False),
        pa.field("plasmid", pa.string(), nullable=False),
        pa.field("test_index", pa.int32(), nullable=False),
        pa.field("module", pa.string(), nullable=False),
        pa.field("seed", pa.int64(), nullable=False),
        pa.field("success", pa.bool_(), nullable=False),
    ]
)


def _metadata(schema: pa.Schema, **values: object) -> pa.Schema:
    encoded = {str(key).encode(): str(value).encode() for key, value in values.items()}
    return schema.with_metadata(encoded)


def _encode_three_mer(sequence: str) -> int:
    """Encode a standard amino-acid 3-mer in the historical alphabet."""
    return (
        AA_TO_INT[sequence[0]] * len(AMINO_ACIDS) ** 2
        + AA_TO_INT[sequence[1]] * len(AMINO_ACIDS)
        + AA_TO_INT[sequence[2]]
    )


def _boolean(value: object) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def build_historical_pattern_index(repo: Path = DEFAULT_REPO) -> HistoricalPatternIndex:
    """Rebuild the exact plasmid pattern population used by the old notebook.

    The 64-million-entry uint8 vector is a dense 3-mer-pair table whose bits
    identify compatible plasmids.  It is built in memory and inherited by fork
    workers, so it does not create another persistent result file.
    """
    source_paths = {
        "site_i": repo / "output/hurdler_site_i_data.csv",
        "site_ii": repo / "output/hurdler_site_ii_data.csv",
        "pairing": repo / "output/site_i_site_ii_pairing_matrix.csv",
        "site_ii_enzymes": repo / "output/selected_site_ii_enzymes.csv",
        "site_iii_enzymes": repo / "output/selected_site_iii_enzymes.csv",
        "plasmids": repo / "data/reference_output/plasmid_digest_check.csv",
    }
    missing = [str(path) for path in source_paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing historical success inputs: {missing}")

    site_i = pd.read_csv(source_paths["site_i"])
    site_ii = pd.read_csv(source_paths["site_ii"])
    pairing = pd.read_csv(source_paths["pairing"], index_col=0)
    pairing = pairing.apply(lambda column: column.map(_boolean))
    site_ii_enzymes = pd.read_csv(source_paths["site_ii_enzymes"]).set_index("enzyme")
    site_iii_enzymes = pd.read_csv(source_paths["site_iii_enzymes"])
    plasmids = pd.read_csv(source_paths["plasmids"], index_col=0)
    if tuple(plasmids.columns) != PLASMIDS:
        raise ValueError("Historical plasmid table order differs from the maintained order")

    site_i = site_i.loc[site_i["enzyme"].isin(pairing.index)].copy()
    site_ii = site_ii.loc[site_ii["enzyme"].isin(pairing.columns)].copy()
    site_i["code"] = site_i["3mer_aa"].map(_encode_three_mer)
    site_ii["code"] = site_ii["3mer_aa"].map(_encode_three_mer)
    site_i_codes = {
        enzyme: np.unique(group["code"].to_numpy(dtype=np.int32))
        for enzyme, group in site_i.groupby("enzyme", sort=False)
    }
    # The notebook treated any non-'right' value as the left-search ordering.
    site_ii["historical_direction"] = np.where(
        site_ii["search_direction"].eq("right"), "right", "left"
    )
    site_ii_codes = {
        (enzyme, direction): np.unique(group["code"].to_numpy(dtype=np.int32))
        for (enzyme, direction), group in site_ii.groupby(
            ["enzyme", "historical_direction"], sort=False
        )
    }

    masks = np.zeros(THREE_MER_SPACE * THREE_MER_SPACE, dtype=np.uint8)
    for bit, plasmid in enumerate(PLASMIDS):
        for enzyme_ii in pairing.columns:
            if enzyme_ii not in plasmids.index or not _boolean(
                plasmids.loc[enzyme_ii, plasmid]
            ):
                continue
            overhang = int(site_ii_enzymes.loc[enzyme_ii, "ovhg"])
            if not site_iii_enzymes["ovhg"].eq(overhang).any():
                continue
            compatible_i = [
                enzyme_i
                for enzyme_i in pairing.index[pairing[enzyme_ii]]
                if enzyme_i in site_i_codes
                and enzyme_i in plasmids.index
                and _boolean(plasmids.loc[enzyme_i, plasmid])
            ]
            if not compatible_i:
                continue
            i_codes = np.unique(np.concatenate([site_i_codes[name] for name in compatible_i]))
            for direction in ("right", "left"):
                ii_codes = site_ii_codes.get((enzyme_ii, direction))
                if ii_codes is None:
                    continue
                if direction == "right":
                    flat_indices = (
                        i_codes[:, None] * THREE_MER_SPACE + ii_codes[None, :]
                    ).ravel()
                else:
                    flat_indices = (
                        ii_codes[:, None] * THREE_MER_SPACE + i_codes[None, :]
                    ).ravel()
                masks[flat_indices] |= np.uint8(1 << bit)

    nonzero = np.flatnonzero(masks)
    present = np.zeros(THREE_MER_SPACE, dtype=np.bool_)
    present[nonzero // THREE_MER_SPACE] = True
    present[nonzero % THREE_MER_SPACE] = True
    return HistoricalPatternIndex(
        masks=masks,
        three_mer_present=present,
        source_hashes={name: sha256_file(path) for name, path in source_paths.items()},
    )


def _set_historical_index(index: HistoricalPatternIndex) -> None:
    global _HISTORICAL_INDEX
    _HISTORICAL_INDEX = index


def _require_historical_index() -> HistoricalPatternIndex:
    if _HISTORICAL_INDEX is None:  # pragma: no cover - defensive multiprocessing guard
        raise RuntimeError("Worker historical pattern index was not initialized")
    return _HISTORICAL_INDEX


def historical_success_mask(
    module: str,
    index: HistoricalPatternIndex,
    *,
    scan_copies: int = 2,
) -> int:
    """Return successful plasmid bits using the parameterized repeated matcher.

    The old notebook used ``re.search()`` once for each direction-specific
    pattern.  Therefore, if the earliest non-overlapping occurrence of a
    pattern has an invalid span, a later occurrence of that same pattern is
    not reconsidered.  This otherwise surprising detail is retained so the
    committed 7--60AA results can be reproduced exactly.
    """
    if scan_copies not in (2, 3):
        raise ValueError("scan_copies must be 2 or 3")
    module_length = len(module)
    scanned_sequence = module * scan_copies
    codes = [
        _encode_three_mer(scanned_sequence[position : position + 3])
        for position in range(len(scanned_sequence) - 2)
    ]
    # This is the original helper's explicit "two different 3-mers" gate.
    if np.unique(np.asarray(codes, dtype=np.int16)[index.three_mer_present[codes]]).size < 2:
        return 0

    positions: dict[int, list[int]] = {}
    for position, code in enumerate(codes):
        positions.setdefault(code, []).append(position)
    first_distance: dict[int, int | None] = {}
    successes = 0
    all_plasmids = (1 << len(PLASMIDS)) - 1
    for left_position, left_code in enumerate(codes):
        right_stop = min(len(codes), left_position + module_length)
        for right_position in range(left_position + 5, right_stop):
            right_code = codes[right_position]
            key = left_code * THREE_MER_SPACE + right_code
            candidate_mask = int(index.masks[key]) & ~successes
            if not candidate_mask:
                continue
            if key not in first_distance:
                distance: int | None = None
                right_positions = positions[right_code]
                for possible_left in positions[left_code]:
                    location = bisect.bisect_left(right_positions, possible_left + 3)
                    if location < len(right_positions):
                        distance = right_positions[location] - possible_left
                        break
                first_distance[key] = distance
            distance = first_distance[key]
            if distance is not None and 5 <= distance < module_length:
                successes |= candidate_mask
                if successes == all_plasmids:
                    return successes
    return successes


def _set_scan_copies(scan_copies: int) -> None:
    global _SCAN_COPIES
    if scan_copies not in (2, 3):
        raise ValueError("scan_copies must be 2 or 3")
    _SCAN_COPIES = scan_copies


def _short_tasks(max_length: int = 5) -> list[tuple[int, str]]:
    """Use similarly sized prefix tasks without changing motif enumeration."""
    tasks: list[tuple[int, str]] = []
    for length in range(1, max_length + 1):
        if length <= 2:
            tasks.append((length, ""))
        elif length <= 4:
            tasks.extend((length, prefix) for prefix in AMINO_ACIDS)
        else:
            tasks.extend((length, "".join(prefix)) for prefix in itertools.product(AMINO_ACIDS, repeat=2))
    return tasks


def _motifs(length: int, prefix: str) -> Iterator[str]:
    for suffix in itertools.product(AMINO_ACIDS, repeat=length - len(prefix)):
        yield prefix + "".join(suffix)


def _screen_short_task(task: tuple[int, str]) -> pa.Table:
    length, prefix = task
    rows: list[dict[str, object]] = []
    index = _require_historical_index()
    for motif in _motifs(length, prefix):
        expanded, expansion_copies = expand_short_module(motif)
        mask = historical_success_mask(expanded, index, scan_copies=_SCAN_COPIES)
        row: dict[str, object] = {
            "motif": motif,
            "motif_length": length,
            "expansion_copies": expansion_copies,
            "expanded_module": expanded,
            "expanded_length": len(expanded),
        }
        for bit, plasmid in enumerate(PLASMIDS):
            row[f"{plasmid}_success"] = bool(mask & (1 << bit))
        row["success_mask"] = mask
        row["successful_plasmids"] = mask.bit_count()
        row["any_success"] = bool(mask)
        rows.append(row)
    return pa.Table.from_pylist(rows, schema=_short_schema())


def _random_tasks(
    min_length: int,
    max_length: int,
    tests: int,
) -> Iterable[tuple[int, str, int, list[str]]]:
    """Generate inputs in the frozen notebook order before parallel matching."""
    if min_length <= 6 <= max_length:
        generator = random.Random(SIX_AA_SEED)
        for plasmid in PLASMIDS:
            modules = [
                "".join(generator.choice(AMINO_ACIDS) for _ in range(6))
                for _ in range(tests)
            ]
            yield 6, plasmid, SIX_AA_SEED, modules
    legacy_minimum = max(7, min_length)
    if legacy_minimum <= max_length:
        generator = random.Random(LEGACY_SEED)
        # Advance the legacy generator exactly as the original 7--60 loop did
        # when a restricted range is requested for a smoke test.
        for length in range(7, max_length + 1):
            for plasmid in PLASMIDS:
                modules = [
                    "".join(generator.choice(AMINO_ACIDS) for _ in range(length))
                    for _ in range(tests)
                ]
                if length >= legacy_minimum:
                    yield length, plasmid, LEGACY_SEED, modules


def _screen_random_task(task: tuple[int, str, int, list[str]]) -> pa.Table:
    length, plasmid, seed, modules = task
    index = _require_historical_index()
    plasmid_bit = PLASMIDS.index(plasmid)
    rows: list[dict[str, object]] = []
    for test_index, module in enumerate(modules):
        success = bool(
            historical_success_mask(module, index, scan_copies=_SCAN_COPIES)
            & (1 << plasmid_bit)
        )
        rows.append(
            {
                "module_length": length,
                "plasmid": plasmid,
                "test_index": test_index,
                "module": module,
                "seed": seed,
                "success": success,
            }
        )
    return pa.Table.from_pylist(rows, schema=RANDOM_SCHEMA)


def _ordered_parallel_map(function, tasks: Sequence | Iterable, workers: int):
    if workers == 1:
        yield from map(function, tasks)
        return
    context = mp.get_context("fork")
    with context.Pool(processes=workers) as pool:
        yield from pool.imap(function, tasks, chunksize=1)


def _write_tables_atomic(
    destination: Path,
    tables: Iterable[pa.Table],
    schema: pa.Schema,
) -> int:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.partial-{os.getpid()}")
    writer: pq.ParquetWriter | None = None
    rows = 0
    try:
        writer = pq.ParquetWriter(temporary, schema, compression="zstd")
        for table in tables:
            writer.write_table(table.cast(schema))
            rows += table.num_rows
        writer.close()
        writer = None
        temporary.replace(destination)
    except BaseException:
        if writer is not None:
            writer.close()
        temporary.unlink(missing_ok=True)
        raise
    return rows


def run_short(
    index: HistoricalPatternIndex,
    output: Path,
    workers: int,
    max_length: int = 5,
    scan_copies: int = 2,
) -> int:
    _set_historical_index(index)
    _set_scan_copies(scan_copies)
    tasks = _short_tasks(max_length)
    schema = _metadata(
        _short_schema(),
        artifact="exhaustive-short-motif-screen",
        rules=SUCCESS_RULE_PROFILE,
        matching_sequence="+".join(["module"] * scan_copies),
        scan_copies=scan_copies,
        distance="5<=d<effective_module_length",
        regex_semantics="first-nonoverlapping-match-per-pattern",
        alphabet=AMINO_ACIDS,
        workers=workers,
        max_length=max_length,
        expected_rows=sum(20**length for length in range(1, max_length + 1)),
    )
    return _write_tables_atomic(
        output,
        _ordered_parallel_map(_screen_short_task, tasks, workers),
        schema,
    )


def run_random(
    index: HistoricalPatternIndex,
    output: Path,
    workers: int,
    min_length: int = 6,
    max_length: int = 60,
    tests: int = 1000,
    scan_copies: int = 2,
) -> int:
    _set_historical_index(index)
    _set_scan_copies(scan_copies)
    schema = _metadata(
        RANDOM_SCHEMA,
        artifact="seeded-random-module-screen",
        rules=SUCCESS_RULE_PROFILE,
        matching_sequence="+".join(["module"] * scan_copies),
        scan_copies=scan_copies,
        distance="5<=d<module_length",
        regex_semantics="first-nonoverlapping-match-per-pattern",
        alphabet=AMINO_ACIDS,
        workers=workers,
        min_length=min_length,
        max_length=max_length,
        tests_per_length_plasmid=tests,
        seed_6aa=SIX_AA_SEED,
        seed_7_60aa=LEGACY_SEED,
    )
    return _write_tables_atomic(
        output,
        _ordered_parallel_map(
            _screen_random_task,
            _random_tasks(min_length, max_length, tests),
            workers,
        ),
        schema,
    )


def _sql_path(path: Path) -> str:
    return str(path).replace("'", "''")


def validate_outputs(
    short_path: Path,
    random_path: Path,
    *,
    golden_path: Path | None = DEFAULT_REPO / "output/hurdler_success_rate_7_60aa_per_plasmid.csv",
    short_max_length: int = 5,
    random_min_length: int = 6,
    random_max_length: int = 60,
    tests: int = 1000,
) -> dict[str, object]:
    for path in (short_path, random_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty result file: {path}")
        pq.ParquetFile(path).scan_contents()
    connection = duckdb.connect()
    short = _sql_path(short_path)
    random_path_sql = _sql_path(random_path)
    short_counts = connection.execute(
        f"SELECT motif_length, count(*) AS row_count, count(DISTINCT motif) AS unique_motifs, "
        f"min(expanded_length) min_expanded, max(expanded_length) max_expanded "
        f"FROM read_parquet('{short}') GROUP BY motif_length ORDER BY motif_length"
    ).fetchdf()
    expected_short = {length: 20**length for length in range(1, short_max_length + 1)}
    observed_lengths = set(short_counts["motif_length"].astype(int))
    if observed_lengths != set(expected_short):
        raise ValueError(f"Unexpected short lengths: {sorted(observed_lengths)}")
    for row in short_counts.itertuples(index=False):
        expected = expected_short[int(row.motif_length)]
        if int(row.row_count) != expected or int(row.unique_motifs) != expected:
            raise ValueError(f"Invalid exhaustive count: {row}")
        if int(row.min_expanded) < 6:
            raise ValueError(f"Short expansion below 6 AA: {row}")

    random_counts = connection.execute(
        f"SELECT module_length, plasmid, count(*) AS row_count, "
        f"count(DISTINCT test_index) AS test_count, "
        f"min(length(module)) min_length, max(length(module)) max_length "
        f"FROM read_parquet('{random_path_sql}') GROUP BY module_length, plasmid "
        f"ORDER BY module_length, plasmid"
    ).fetchdf()
    expected_groups = (random_max_length - random_min_length + 1) * len(PLASMIDS)
    if len(random_counts) != expected_groups:
        raise ValueError(f"Expected {expected_groups} random groups, found {len(random_counts)}")
    if set(random_counts["plasmid"]) != set(PLASMIDS):
        raise ValueError("Random result does not contain all maintained plasmids")
    if not (
        random_counts["row_count"].eq(tests)
        & random_counts["test_count"].eq(tests)
    ).all():
        raise ValueError("Random group count or test-index uniqueness is invalid")
    if not (
        random_counts["min_length"].eq(random_counts["module_length"])
        & random_counts["max_length"].eq(random_counts["module_length"])
    ).all():
        raise ValueError("Random module sequence length does not match its group")

    golden_mismatches = None
    if golden_path is not None and random_min_length <= 7 and random_max_length >= 60:
        if not golden_path.is_file():
            raise FileNotFoundError(f"Missing historical golden result: {golden_path}")
        golden = pd.read_csv(golden_path)[
            ["module_length", "plasmid", "n_success", "n_total"]
        ]
        observed = connection.execute(
            f"SELECT module_length, plasmid, "
            f"sum(CASE WHEN success THEN 1 ELSE 0 END)::BIGINT AS n_success, "
            f"count(*)::BIGINT AS n_total FROM read_parquet('{random_path_sql}') "
            f"WHERE module_length BETWEEN 7 AND 60 GROUP BY module_length, plasmid"
        ).fetchdf()
        comparison = golden.merge(
            observed,
            on=["module_length", "plasmid"],
            how="outer",
            suffixes=("_golden", "_observed"),
            indicator=True,
        )
        golden_mismatches = int(
            (
                comparison["_merge"].ne("both")
                | comparison["n_success_golden"].ne(comparison["n_success_observed"])
                | comparison["n_total_golden"].ne(comparison["n_total_observed"])
            ).sum()
        )
        if golden_mismatches:
            raise ValueError(
                f"Historical 7--60AA golden regression has {golden_mismatches} mismatched groups"
            )
    connection.close()
    return {
        "passed": True,
        "short_rows": int(short_counts["row_count"].sum()),
        "short_expected": sum(expected_short.values()),
        "short_counts": {
            str(int(row.motif_length)): int(row.row_count)
            for row in short_counts.itertuples(index=False)
        },
        "random_rows": int(random_counts["row_count"].sum()),
        "random_expected": expected_groups * tests,
        "historical_7_60_golden_groups": 432 if golden_mismatches is not None else None,
        "historical_7_60_golden_mismatches": golden_mismatches,
        "short_sha256": sha256_file(short_path),
        "random_sha256": sha256_file(random_path),
    }


def success_summary(short_path: Path, random_path: Path) -> pd.DataFrame:
    connection = duckdb.connect()
    short = _sql_path(short_path)
    random_path_sql = _sql_path(random_path)
    parts: list[str] = []
    for order, plasmid in enumerate(PLASMIDS):
        column = f'{plasmid}_success'.replace('"', '""')
        label = plasmid.replace("'", "''")
        parts.append(
            f"SELECT motif_length AS module_length, '{label}' AS plasmid, {order} AS plasmid_order, "
            f"sum(CASE WHEN \"{column}\" THEN 1 ELSE 0 END)::BIGINT successes, "
            f"count(*)::BIGINT tests, 'exhaustive_expanded' AS method "
            f"FROM read_parquet('{short}') GROUP BY motif_length"
        )
    parts.append(
        f"SELECT module_length, plasmid, CASE plasmid "
        + " ".join(
            f"WHEN '{plasmid.replace(chr(39), chr(39) * 2)}' THEN {order}"
            for order, plasmid in enumerate(PLASMIDS)
        )
        + " END AS plasmid_order, sum(CASE WHEN success THEN 1 ELSE 0 END)::BIGINT successes, "
        + "count(*)::BIGINT tests, 'monte_carlo' AS method "
        + f"FROM read_parquet('{random_path_sql}') GROUP BY module_length, plasmid"
    )
    frame = connection.execute(" UNION ALL ".join(parts)).fetchdf()
    connection.close()
    frame["success_rate"] = frame["successes"] / frame["tests"]
    return frame.sort_values(["plasmid_order", "module_length"]).reset_index(drop=True)


def compare_scan_copy_results(
    two_short: Path,
    two_random: Path,
    three_short: Path,
    three_random: Path,
    output: Path,
) -> pd.DataFrame:
    """Compare 2-copy and 3-copy scans on identical exhaustive/random inputs."""
    connection = duckdb.connect()
    identity_queries = [
        (two_short, three_short, "motif, motif_length, expanded_module, expanded_length"),
        (two_random, three_random, "module_length, plasmid, test_index, module, seed"),
    ]
    for left, right, columns in identity_queries:
        left_sql = _sql_path(left)
        right_sql = _sql_path(right)
        differences = connection.execute(
            f"SELECT count(*) FROM ("
            f"(SELECT {columns} FROM read_parquet('{left_sql}') EXCEPT ALL "
            f" SELECT {columns} FROM read_parquet('{right_sql}')) UNION ALL "
            f"(SELECT {columns} FROM read_parquet('{right_sql}') EXCEPT ALL "
            f" SELECT {columns} FROM read_parquet('{left_sql}')))"
        ).fetchone()[0]
        if differences:
            raise ValueError(
                f"2-copy and 3-copy observations differ in {differences} identity rows"
            )
    connection.close()

    two = success_summary(two_short, two_random).rename(
        columns={"successes": "successes_2x", "success_rate": "success_rate_2x"}
    )
    three = success_summary(three_short, three_random).rename(
        columns={"successes": "successes_3x", "success_rate": "success_rate_3x"}
    )
    expected_rows = len(two)
    comparison = two.merge(
        three[
            [
                "module_length",
                "plasmid",
                "tests",
                "successes_3x",
                "success_rate_3x",
            ]
        ],
        on=["module_length", "plasmid", "tests"],
        how="outer",
        validate="one_to_one",
    )
    if len(comparison) != expected_rows or comparison.isna().any().any():
        raise ValueError("Incomplete 2-copy/3-copy comparison")
    comparison["success_delta"] = (
        comparison["successes_3x"] - comparison["successes_2x"]
    )
    comparison["success_rate_delta"] = (
        comparison["success_rate_3x"] - comparison["success_rate_2x"]
    )
    comparison["improved"] = comparison["success_delta"].gt(0)
    comparison = comparison.sort_values(["plasmid_order", "module_length"])
    output.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(output, index=False)
    by_length = (
        comparison.groupby("module_length", sort=True)
        .agg(
            tests_per_plasmid=("tests", "first"),
            improved_plasmids=("improved", "sum"),
            added_successes_all_plasmids=("success_delta", "sum"),
            mean_success_rate_2x=("success_rate_2x", "mean"),
            mean_success_rate_3x=("success_rate_3x", "mean"),
            mean_success_rate_delta=("success_rate_delta", "mean"),
            maximum_plasmid_rate_delta=("success_rate_delta", "max"),
        )
        .reset_index()
    )
    improved_names = (
        comparison.loc[comparison["improved"]]
        .groupby("module_length")["plasmid"]
        .agg(lambda values: ";".join(values))
    )
    by_length["improved_plasmid_names"] = by_length["module_length"].map(
        improved_names
    ).fillna("")
    by_length["any_improvement"] = by_length["improved_plasmids"].gt(0)
    by_length.to_csv(output.with_name("scan_copy_improvement_by_length.csv"), index=False)
    return comparison


def plot_success_curve(
    summary: pd.DataFrame,
    figure_dir: Path,
    *,
    file_stem: str = "success_rate_1_60",
) -> list[Path]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    # Retain the original get_re_sites.ipynb typography/line settings,
    # plasmid/color order, exact text, and 50AA upper limit on a user-requested
    # near-square 6 x 5 inch canvas. Only the 1--5AA values and 3-copy scan are new.
    original_plasmid_order = (
        "pET-28a(+)",
        "pET-28a(+)_start_codon",
        "pGEX-4T-1",
        "pMAL-c5X",
        "pUC18",
        "pQE-3",
        "pCold_I",
        "pET-21a(+)",
    )
    displayed = summary.loc[summary["module_length"].le(50)]
    with plt.style.context("default"):
        fig, axis = plt.subplots(figsize=(6, 5))
        for plasmid in original_plasmid_order:
            group = displayed.loc[displayed["plasmid"] == plasmid].sort_values(
                "module_length"
            )
            axis.plot(
                group["module_length"],
                group["success_rate"] * 100,
                marker="o",
                label=plasmid,
            )
        axis.set_xlabel("Sequence Length")
        axis.set_ylabel("Probability (%)")
        axis.set_title("3-mer Probability vs Sequence Length")
        axis.set_xlim(1, 50)
        axis.legend(title="Plasmid")
        axis.grid(True)
    outputs = [figure_dir / f"{file_stem}.png", figure_dir / f"{file_stem}.pdf"]
    for path in outputs:
        fig.savefig(path, dpi=600)
    plt.close(fig)
    return outputs


def benchmark(
    index: HistoricalPatternIndex,
    worker_counts: Sequence[int],
    *,
    scan_copies: int = 2,
) -> dict[str, object]:
    """Benchmark an identical 160,000-motif workload and verify equality."""
    task = (5, "A")
    tasks = [task]
    # Split the same 160,000 motifs into 20 equal prefix tasks so multiple
    # workers can participate without changing the tested sequence set.
    parallel_tasks = [(5, "A" + aa) for aa in AMINO_ACIDS]
    _set_historical_index(index)
    _set_scan_copies(scan_copies)
    results: dict[str, object] = {}
    reference: pd.DataFrame | None = None
    reference_seconds: float | None = None
    for workers in worker_counts:
        started = time.perf_counter()
        selected_tasks = tasks if workers == 1 else parallel_tasks
        tables = list(_ordered_parallel_map(_screen_short_task, selected_tasks, workers))
        frame = pa.concat_tables(tables).to_pandas().sort_values("motif").reset_index(drop=True)
        elapsed = time.perf_counter() - started
        if reference is None:
            reference = frame
            reference_seconds = elapsed
            exact = True
        else:
            exact = frame.equals(reference)
        if not exact:
            raise ValueError(f"Multiprocessing result mismatch at {workers} workers")
        results[str(workers)] = {
            "seconds": round(elapsed, 6),
            "rows": len(frame),
            "exact_match": exact,
            "speedup_vs_serial": round(reference_seconds / elapsed, 3) if reference_seconds else 1.0,
        }
    return {"benchmark": results, "motifs": 160000, "passed": True}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-dir", type=Path, default=DEFAULT_REPO)
    parser.add_argument("--short-output", type=Path, default=DEFAULT_SHORT_OUTPUT)
    parser.add_argument("--random-output", type=Path, default=DEFAULT_RANDOM_OUTPUT)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--short-max-length", type=int, default=5)
    parser.add_argument("--random-min-length", type=int, default=6)
    parser.add_argument("--random-max-length", type=int, default=60)
    parser.add_argument("--tests", type=int, default=1000)
    parser.add_argument("--scan-copies", type=int, choices=(2, 3), default=3)
    parser.add_argument("--golden-path", type=Path)
    parser.add_argument("--compare-two-short", type=Path)
    parser.add_argument("--compare-two-random", type=Path)
    parser.add_argument("--compare-three-short", type=Path)
    parser.add_argument("--compare-three-random", type=Path)
    parser.add_argument("--comparison-output", type=Path)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--benchmark-workers", default="1,16")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    comparison_paths = (
        args.compare_two_short,
        args.compare_two_random,
        args.compare_three_short,
        args.compare_three_random,
        args.comparison_output,
    )
    if any(path is not None for path in comparison_paths):
        if not all(path is not None for path in comparison_paths):
            raise ValueError("All five comparison paths must be supplied together")
        comparison = compare_scan_copy_results(*comparison_paths)
        improved = comparison.loc[comparison["improved"]]
        print(
            json.dumps(
                {
                    "comparison_rows": len(comparison),
                    "improved_length_plasmid_rows": len(improved),
                    "improved_module_lengths": sorted(
                        improved["module_length"].astype(int).unique().tolist()
                    ),
                    "maximum_success_rate_delta": float(
                        comparison["success_rate_delta"].max()
                    ),
                    "output": str(args.comparison_output),
                    "by_length_output": str(
                        args.comparison_output.with_name(
                            "scan_copy_improvement_by_length.csv"
                        )
                    ),
                },
                indent=2,
            )
        )
        return 0
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    if not 1 <= args.short_max_length <= 5:
        raise ValueError("--short-max-length must be in 1..5")
    if not 6 <= args.random_min_length <= args.random_max_length <= 60:
        raise ValueError("Random length range must be within 6..60")
    index = build_historical_pattern_index(args.repo_dir)
    if args.benchmark:
        worker_counts = [int(value) for value in args.benchmark_workers.split(",")]
        print(
            json.dumps(
                benchmark(index, worker_counts, scan_copies=args.scan_copies), indent=2
            )
        )
        return 0
    started = time.time()
    if not args.validate_only:
        short_rows = run_short(
            index,
            args.short_output,
            args.workers,
            args.short_max_length,
            args.scan_copies,
        )
        random_rows = run_random(
            index,
            args.random_output,
            args.workers,
            args.random_min_length,
            args.random_max_length,
            args.tests,
            args.scan_copies,
        )
        print(json.dumps({"short_rows_written": short_rows, "random_rows_written": random_rows}))
    validation = validate_outputs(
        args.short_output,
        args.random_output,
        short_max_length=args.short_max_length,
        random_min_length=args.random_min_length,
        random_max_length=args.random_max_length,
        tests=args.tests,
        golden_path=(
            args.golden_path
            if args.golden_path is not None
            else args.repo_dir / "output/hurdler_success_rate_7_60aa_per_plasmid.csv"
            if args.scan_copies == 2
            else None
        ),
    )
    summary = success_summary(args.short_output, args.random_output)
    figures = plot_success_curve(
        summary,
        args.figure_dir,
        file_stem=(
            "success_rate_1_60" if args.scan_copies == 2 else "success_rate_1_60_scan_3x"
        ),
    )
    validation.update(
        {
            "elapsed_seconds": round(time.time() - started, 3),
            "figures": [str(path) for path in figures],
            "rules": SUCCESS_RULE_PROFILE,
            "matching_sequence": "+".join(["module"] * args.scan_copies),
            "scan_copies": args.scan_copies,
            "source_hashes": index.source_hashes,
            "plasmids": list(PLASMIDS),
        }
    )
    print(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
