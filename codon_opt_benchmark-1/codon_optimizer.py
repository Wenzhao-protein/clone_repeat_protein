"""Codon optimization script ported from the notebook workflow.

The script runs a genetic algorithm to minimize synthesis/cloning pain points:
- repeated n-mers in the coding region
- repeated n-mers when flanked by provided context
- high GC windows
- simple hairpin-like reverse-complement windows

Usage (example):
    python codon_optimizer.py --tasks batch_optimization_tasks.csv \
        --output-dir ./codon_opt_results --pop-size 400 --max-iter 50

Input CSV columns (required/optional):
- aa_seq (required): amino-acid sequence to optimize
- name (optional): identifier; falls back to row index
- prefix (optional): DNA to prepend before the optimized region
- suffix (optional): DNA to append after the optimized region
"""

from __future__ import annotations

import argparse
import csv
import functools
import math
import os
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from sko.GA import GA
from sko.tools import set_run_mode


CODON_TABLE: Dict[str, List[str]] = {
    "A": ["GCT", "GCA", "GCC", "GCG"],
    "C": ["TGT", "TGC"],
    "D": ["GAC", "GAT"],
    "E": ["GAG", "GAA"],
    "F": ["TTC", "TTT"],
    "G": ["GGA", "GGG", "GGT", "GGC"],
    "H": ["CAC", "CAT"],
    "I": ["ATC", "ATT"],
    "K": ["AAG", "AAA"],
    "L": ["CTC", "CTT", "TTG", "TTA", "CTG"],
    "M": ["ATG"],
    "N": ["AAT", "AAC"],
    "P": ["CCC", "CCT", "CCA", "CCG"],
    "Q": ["CAA", "CAG"],
    "R": ["CGG", "CGC", "CGT"],
    "S": ["TCG", "TCA", "TCC", "AGT", "TCT", "AGC"],
    "T": ["ACA", "ACT", "ACG", "ACC"],
    "V": ["GTA", "GTC", "GTT", "GTG"],
    "W": ["TGG"],
    "Y": ["TAC", "TAT"],
}


REGISTRATION_SITES: Tuple[str, ...] = ("AAGCTT", "GAGCTC", "TGTACA")  # HindIII, SacI, BsrGI


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run codon optimization GA")
    parser.add_argument("--tasks", required=True, type=Path, help="CSV with aa_seq column")
    parser.add_argument("--output-dir", type=Path, default=Path("codon_opt_results"), help="Folder for outputs")
    parser.add_argument("--pop-size", type=int, default=800, help="GA population size (ipynb default 800)")
    parser.add_argument("--max-iter", type=int, default=1000, help="GA generations (ipynb default 1000)")
    parser.add_argument("--mutation-prob", type=float, default=0.0015, help="GA mutation probability (ipynb default 0.0015)")
    parser.add_argument(
        "--mutation-grid",
        type=str,
        default=None,
        help="Comma-separated mutation probs to sweep; overrides --mutation-prob",
    )
    parser.add_argument("--run-mode", choices=["common", "multithreading", "multiprocessing"], default="multiprocessing")
    parser.add_argument("--verbose", action="store_true", help="Save every individual per generation")
    return parser.parse_args()


def load_tasks(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if "aa_seq" not in df.columns:
        raise ValueError("Input CSV must contain an 'aa_seq' column")
    df["name"] = df.get("name", pd.Series(df.index.astype(str)))
    df["prefix"] = df.get("prefix", "")
    df["suffix"] = df.get("suffix", "")
    return df


def reverse_translate_vector(aa_seq: str, vector: Iterable[float]) -> str:
    # sanitize GA parameters to avoid NaNs/infs during translation
    vals = list(vector)
    if len(vals) != len(aa_seq):
        raise ValueError(f"Length mismatch: got {len(vals)} params for {len(aa_seq)} aa")
    cleaned = np.nan_to_num(np.asarray(vals, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)

    dna_parts: List[str] = []
    for aa, val in zip(aa_seq, cleaned):
        codons = CODON_TABLE.get(aa)
        if not codons:
            raise ValueError(f"Unsupported amino acid '{aa}'")
        idx = int(round(val)) % len(codons)
        dna_parts.append(codons[idx])
    return "".join(dna_parts)


def nmer_repeat_count(seq: str, n: int) -> int:
    if len(seq) <= n:
        return 0
    windows = [seq[i : i + n] for i in range(len(seq) - n + 1)]
    return max(0, len(windows) - len(set(windows)))


def count_high_gc_windows(seq: str, window: int, threshold: float) -> int:
    if len(seq) < window:
        return 0
    hits = 0
    for i in range(len(seq) - window + 1):
        frag = seq[i : i + window]
        gc_ratio = (frag.count("G") + frag.count("C")) / window
        if gc_ratio >= threshold:
            hits += 1
    return hits


def hairpin_count(seq: str, n: int = 10) -> int:
    if len(seq) <= n:
        return 0
    windows = [seq[i : i + n] for i in range(len(seq) - n + 1)]
    rc_windows = [str(Seq(frag).reverse_complement()) for frag in set(windows)]
    return len(set(windows)) + len(set(rc_windows)) - len(set(windows + rc_windows))


def registration_bonus(seq: str) -> float:
    return sum(0.1 for site in REGISTRATION_SITES if site in seq)


def objective_func(params: Iterable[float], aa_seq: str, prefix: str, suffix: str) -> float:
    try:
        core_dna = reverse_translate_vector(aa_seq, params)
    except Exception:
        # Guard against bad population members
        return 1e9

    full_dna = prefix + core_dna + suffix

    score = 0.0
    for n in (8, 13):
        score += nmer_repeat_count(core_dna, n)

    score += 0.05 * nmer_repeat_count(full_dna, 14)
    score += registration_bonus(core_dna)
    score += 0.02 * count_high_gc_windows(core_dna, 20, 0.9)
    score += 0.02 * count_high_gc_windows(core_dna, 60, 0.68)
    score += 0.02 * hairpin_count(full_dna, 10)
    return score


def make_objective(aa_seq: str, prefix: str, suffix: str):
    # functools.partial keeps the objective picklable for multiprocessing.
    return functools.partial(objective_func, aa_seq=aa_seq, prefix=prefix, suffix=suffix)


def run_ga_for_sequence(
    name: str,
    aa_seq: str,
    prefix: str,
    suffix: str,
    pop_size: int,
    max_iter: int,
    mutation_prob: float,
    run_mode: str,
    verbose: bool,
    output_dir: Path,
):
    objective = make_objective(aa_seq, prefix, suffix)
    set_run_mode(objective, run_mode)

    ub = [len(CODON_TABLE[aa]) - 1 for aa in aa_seq]
    ga = GA(
        func=objective,
        n_dim=len(aa_seq),
        size_pop=pop_size,
        max_iter=1,
        prob_mut=mutation_prob,
        lb=[0] * len(aa_seq),
        ub=ub,
        precision=1,
    )
    best_score_overall = math.inf
    best_dna_overall = ""

    # Prepare output files (overwrite if exist) and write headers once.
    best_path = output_dir / f"{name}_best_per_iteration.csv"
    best_file = best_path.open("w", newline="")
    best_writer = csv.writer(best_file)
    best_writer.writerow(["iteration", "aa_seq", "dna_seq", "score"])

    verbose_file = None
    verbose_writer = None
    if verbose:
        verbose_path = output_dir / f"{name}_all_per_iteration.csv"
        verbose_file = verbose_path.open("w", newline="")
        verbose_writer = csv.writer(verbose_file)
        verbose_writer.writerow(["iteration", "idx", "score"])

    for gen in range(1, max_iter + 1):
        best_x, best_y = ga.run(1)
        best_y_val = float(np.asarray(best_y).reshape(-1)[0])
        best_dna_core = reverse_translate_vector(aa_seq, best_x)
        best_dna_full = prefix + best_dna_core + suffix

        if best_y_val < best_score_overall:
            best_score_overall = best_y_val
            best_dna_overall = best_dna_full

        best_writer.writerow((gen, aa_seq, best_dna_full, best_y_val))
        best_file.flush()

        if verbose:
            scores = np.asarray(ga.Y).reshape(-1)
            for idx, (candidate, score) in enumerate(zip(ga.X, scores)):
                #dna = prefix + reverse_translate_vector(aa_seq, candidate) + suffix
                verbose_writer.writerow((gen, idx, float(score)))
            verbose_file.flush()

        # Clear GA caches to control memory
        if hasattr(ga, "all_history_Y"):
            ga.all_history_Y.clear()
        if hasattr(ga, "all_history_FitV"):
            ga.all_history_FitV.clear()

    best_file.close()
    if verbose_file:
        verbose_file.close()

    return {
        "name": name,
        "aa_seq": aa_seq,
        "core_dna": best_dna_overall[len(prefix) : len(prefix) + 3 * len(aa_seq)] if best_dna_overall else "",
        "best_dna": best_dna_overall,
        "score": best_score_overall,
    }


def main() -> None:
    args = parse_args()
    tasks = load_tasks(args.tasks)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.mutation_grid:
        mutation_list = [float(x) for x in args.mutation_grid.split(",") if x.strip()]
    else:
        mutation_list = [args.mutation_prob]

    for mut in mutation_list:
        sweep_dir = args.output_dir if len(mutation_list) == 1 else args.output_dir / f"mut_{mut}"
        sweep_dir.mkdir(parents=True, exist_ok=True)

        results: List[dict] = []
        for _, row in tasks.iterrows():
            name = str(row["name"])
            aa_seq = str(row["aa_seq"]).strip()
            prefix = str(row.get("prefix", ""))
            suffix = str(row.get("suffix", ""))

            if not aa_seq:
                raise ValueError(f"Empty aa_seq for task '{name}'")

            run_result = run_ga_for_sequence(
                name=name,
                aa_seq=aa_seq,
                prefix=prefix,
                suffix=suffix,
                pop_size=args.pop_size,
                max_iter=args.max_iter,
                mutation_prob=mut,
                run_mode=args.run_mode,
                verbose=args.verbose,
                output_dir=sweep_dir,
            )

            results.append({
                "name": run_result["name"],
                "aa_seq": run_result["aa_seq"],
                "core_dna": run_result["core_dna"],
                "best_dna_with_context": run_result["best_dna"],
                "score": run_result["score"],
                "mutation_prob": mut,
            })

        results_df = pd.DataFrame(results)
        results_path = sweep_dir / "codon_opt_results.csv"
        results_df.to_csv(results_path, index=False)
        print(f"Saved results to {results_path}")


if __name__ == "__main__":
    main()