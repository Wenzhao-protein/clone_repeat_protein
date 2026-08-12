"""Exhaustive 1--5 amino-acid motif screening."""

from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Iterable

import pandas as pd

from .constants import AMINO_ACIDS, PLASMIDS
from .index import PatternIndex
from .matching import materialize_best_solution, query_all_plasmids
from .rules import RuleProfile, LEGACY_OPTIMIZED_V1

SHORT_MOTIF_COUNTS = {length: 20**length for length in range(1, 6)}
SHORT_MOTIF_TOTAL = sum(SHORT_MOTIF_COUNTS.values())


def iter_motifs(length: int, prefix: str = "") -> Iterable[str]:
    if not 1 <= length <= 5:
        raise ValueError("Short motif length must be between 1 and 5")
    if len(prefix) > length or any(aa not in AMINO_ACIDS for aa in prefix):
        raise ValueError(f"Invalid prefix {prefix!r} for motif length {length}")
    for suffix in itertools.product(AMINO_ACIDS, repeat=length - len(prefix)):
        yield prefix + "".join(suffix)


def screen_short_shard(
    index_dir: str | Path,
    output_dir: str | Path,
    *,
    length: int,
    prefix: str = "",
    rules: RuleProfile = LEGACY_OPTIMIZED_V1,
) -> dict[str, object]:
    index = PatternIndex.load(index_dir)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    tag = f"k{length}_{prefix or 'all'}"
    motif_rows: list[dict[str, object]] = []
    hit_rows: list[dict[str, object]] = []

    for motif in iter_motifs(length, prefix):
        results = query_all_plasmids(motif, index, rules=rules, expand_short=True)
        row: dict[str, object] = {
            "motif": motif,
            "motif_length": length,
            "expansion_copies": results[0].expansion_copies,
            "expanded_module": results[0].effective_module,
            "expanded_length": results[0].effective_length,
        }
        success_mask = 0
        for bit, result in enumerate(results):
            row[f"{result.plasmid}_success"] = result.success
            row[f"{result.plasmid}_solution_count"] = result.solution_count
            if result.success:
                success_mask |= 1 << bit
                hit_rows.append(materialize_best_solution(result, index))
        row["success_mask"] = success_mask
        row["successful_plasmids"] = success_mask.bit_count()
        row["any_success"] = bool(success_mask)
        motif_rows.append(row)

    motif_frame = pd.DataFrame(motif_rows)
    hit_frame = pd.DataFrame(hit_rows)
    motif_path = destination / "motif_shards" / f"short_motifs_{tag}.parquet"
    hit_path = (
        destination
        / "short_motif_hits.parquet"
        / f"motif_length={length}"
        / f"prefix={prefix or 'all'}"
        / f"part-{tag}.parquet"
    )
    motif_path.parent.mkdir(parents=True, exist_ok=True)
    hit_path.parent.mkdir(parents=True, exist_ok=True)
    motif_frame.to_parquet(motif_path, index=False)
    if hit_frame.empty:
        hit_frame = pd.DataFrame(columns=["module", "plasmid", "success"])
    hit_frame.to_parquet(hit_path, index=False)
    return {
        "tag": tag,
        "motif_count": len(motif_frame),
        "hit_count": len(hit_frame),
        "motif_path": str(motif_path),
        "hit_path": str(hit_path),
        "expected_count": 20 ** (length - len(prefix)),
    }


def finalize_short_results(output_dir: str | Path) -> dict[str, object]:
    """Merge motif shards, validate exhaustive counts, and write 5x8 summary."""
    try:
        import duckdb
    except ImportError as exc:  # pragma: no cover - dependency guard
        raise RuntimeError("DuckDB is required to finalize short-motif shards") from exc

    root = Path(output_dir)
    shard_glob = str(root / "motif_shards" / "short_motifs_k*.parquet")
    shards = sorted((root / "motif_shards").glob("short_motifs_k*.parquet"))
    if len(shards) != 404:
        raise ValueError(f"Expected 404 motif shards, found {len(shards)}")
    combined = root / "short_motifs.parquet"
    summary_path = root / "short_success_summary.csv"
    validation_path = root / "short_validation.json"
    connection = duckdb.connect()
    escaped_glob = shard_glob.replace("'", "''")
    escaped_combined = str(combined).replace("'", "''")
    connection.execute(
        f"COPY (SELECT * FROM read_parquet('{escaped_glob}', union_by_name=true) "
        f"ORDER BY motif_length, motif) TO '{escaped_combined}' (FORMAT PARQUET, COMPRESSION ZSTD)"
    )

    counts = connection.execute(
        f"SELECT motif_length, count(*) AS rows, count(DISTINCT motif) AS unique_motifs, "
        f"min(expanded_length) AS min_expanded, max(expanded_length) AS max_expanded "
        f"FROM read_parquet('{escaped_combined}') GROUP BY motif_length ORDER BY motif_length"
    ).fetchdf()
    expected = SHORT_MOTIF_COUNTS
    for row in counts.itertuples(index=False):
        if int(row.rows) != expected[int(row.motif_length)] or int(row.unique_motifs) != int(row.rows):
            raise ValueError(f"Invalid exhaustive count for {row.motif_length}AA: {row}")
        if int(row.min_expanded) < 6:
            raise ValueError(f"Expanded module below 6AA for {row.motif_length}AA")

    union_parts: list[str] = []
    for plasmid in PLASMIDS:
        column = plasmid.replace('"', '""') + "_success"
        label = plasmid.replace("'", "''")
        union_parts.append(
            f"SELECT motif_length AS module_length, '{label}' AS plasmid, "
            f"sum(CASE WHEN \"{column}\" THEN 1 ELSE 0 END)::BIGINT AS successes, "
            f"count(*)::BIGINT AS tests FROM read_parquet('{escaped_combined}') GROUP BY motif_length"
        )
    summary = connection.execute(" UNION ALL ".join(union_parts)).fetchdf()
    summary["success_rate"] = summary["successes"] / summary["tests"]
    summary["method"] = "exhaustive_expanded"
    summary = summary.sort_values(["module_length", "plasmid"]).reset_index(drop=True)
    summary.to_csv(summary_path, index=False)
    validation = {
        "shard_count": len(shards),
        "total_motifs": int(counts["rows"].sum()),
        "expected_total": sum(expected.values()),
        "lengths": counts.to_dict(orient="records"),
        "passed": int(counts["rows"].sum()) == sum(expected.values()),
    }
    validation_path.write_text(json.dumps(validation, indent=2, default=int) + "\n")
    connection.close()
    return {
        **validation,
        "motifs": str(combined),
        "hits": str(root / "short_motif_hits.parquet"),
        "summary": str(summary_path),
    }


def summarize_short_results(paths: list[str | Path], output_path: str | Path) -> pd.DataFrame:
    frames = [pd.read_parquet(path) for path in paths]
    frame = pd.concat(frames, ignore_index=True)
    rows: list[dict[str, object]] = []
    for length, group in frame.groupby("motif_length", sort=True):
        for plasmid in PLASMIDS:
            successes = int(group[f"{plasmid}_success"].sum())
            rows.append(
                {
                    "module_length": int(length),
                    "plasmid": plasmid,
                    "successes": successes,
                    "tests": int(len(group)),
                    "success_rate": successes / len(group),
                    "method": "exhaustive_expanded",
                }
            )
    summary = pd.DataFrame(rows)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, index=False)
    return summary
