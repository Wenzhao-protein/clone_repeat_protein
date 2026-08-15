#!/usr/bin/env python3
"""Live-IDT smoke using the actual gBlock selected by production routes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from hurdler.idt import (
    IDTComplexityScorer,
    clear_idt_secret_environment,
    configure_idt_credentials,
)
from hurdler.purchase_orderability import classify_double_stranded_purchase


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--credential-path", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected = None
    for path in sorted(args.raw_root.glob("shard_*/complete_route_seeds.parquet")):
        frame = pd.read_parquet(path)
        candidates = frame.loc[frame.core_sequence.str.len().gt(200)]
        if not candidates.empty:
            selected = candidates.iloc[0]
            break
    if selected is None:
        raise RuntimeError("Production input does not contain a representative gBlock seed")
    configure_idt_credentials(
        mode="path",
        path=args.credential_path,
        headless=True,
        include_path_in_status=False,
    )
    try:
        scorer = IDTComplexityScorer(args.output_dir / "idt_audit.jsonl")
        result = classify_double_stranded_purchase(
            "production-purchase-smoke",
            str(selected.core_sequence),
            idt_scorer=scorer,
        )
    finally:
        clear_idt_secret_environment()
    if result["orderable"] is not True:
        raise RuntimeError(
            f"Representative production gBlock did not pass: {result['failure_reason']}"
        )
    summary = {
        "status": "passed",
        "product_class": result["product_class"],
        "length_bp": result["core_length_bp"],
        "sequence_sha256": hashlib.sha256(
            str(selected.core_sequence).encode()
        ).hexdigest(),
        "idt_status": result["idt_status"],
        "idt_score": result["idt_score"],
        "idt_response_sha256": result["idt_response_sha256"],
        "credential_contents_recorded": False,
    }
    (args.output_dir / "smoke_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
