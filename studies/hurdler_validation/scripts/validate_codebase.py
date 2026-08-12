#!/usr/bin/env python3
"""Digs-side environment, CLI, source, kernel, and artifact smoke checks."""

from __future__ import annotations

import argparse
import compileall
import hashlib
import importlib
import json
import subprocess
from pathlib import Path

import nbformat
from nbclient import NotebookClient


IMPORTS = [
    "Bio",
    "duckdb",
    "h5py",
    "hurdler",
    "matplotlib",
    "numpy",
    "openpyxl",
    "pandas",
    "papermill",
    "PIL",
    "plotly",
    "pyarrow",
    "scipy",
    "seaborn",
    "sko",
]


def checked(command: list[str], cwd: Path) -> str:
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(command)}\n{result.stderr}")
    return result.stdout


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    env_bin = Path("/home/wendai/.conda/envs/hurdler/bin")

    imported = {}
    for name in IMPORTS:
        module = importlib.import_module(name)
        imported[name] = getattr(module, "__version__", "imported")
    pip_check = checked([str(env_bin / "python"), "-m", "pip", "check"], repo).strip()
    if not compileall.compile_dir(repo / "src", quiet=1):
        raise RuntimeError("src compilation failed")
    if not compileall.compile_dir(repo / "studies" / "hurdler_validation" / "scripts", quiet=1):
        raise RuntimeError("study-script compilation failed")
    python_sources = sorted(repo.rglob("*.py"))
    for source in python_sources:
        compile(source.read_text(), str(source), "exec")
    shell_scripts = sorted(repo.rglob("*.sh"))
    for script in shell_scripts:
        checked(["/bin/bash", "-n", str(script)], repo)
    pytest_output = checked([str(env_bin / "python"), "-m", "pytest", "-q", str(repo / "tests")], repo).strip()

    cli_help = {}
    for arguments in (
        ["--help"],
        ["reference", "--help"],
        ["lookup", "--help"],
        ["query", "--help"],
        ["screen-short", "--help"],
        ["success-rate", "--help"],
        ["curate-modules", "--help"],
        ["merge-module-catalogs", "--help"],
        ["infer-boundaries", "--help"],
        ["designed-inventory", "--help"],
        ["validate-designed-structures", "--help"],
        ["infer-designed-boundaries", "--help"],
        ["module-compatibility", "--help"],
        ["adaptive-copy-search", "--help"],
        ["optimize-modules", "--help"],
        ["refine-ga", "--help"],
        ["validate-run", "--help"],
    ):
        cli_help[" ".join(arguments)] = checked([str(env_bin / "hurdler"), *arguments], repo).splitlines()[0]

    book = nbformat.v4.new_notebook(
        cells=[nbformat.v4.new_code_cell("import hurdler, pandas as pd; assert hurdler.__version__")]
    )
    book.metadata.kernelspec = {"display_name": "HURDLER", "language": "python", "name": "hurdler"}
    NotebookClient(book, kernel_name="hurdler", timeout=120).execute(cwd=str(repo))

    payload = {
        "passed": True,
        "imports": imported,
        "pip_check": pip_check,
        "compiled_src": True,
        "compiled_study_scripts": True,
        "python_sources_syntax_checked": [str(path.relative_to(repo)) for path in python_sources],
        "shell_scripts_checked": [str(path.relative_to(repo)) for path in shell_scripts],
        "source_sha256": {
            str(path.relative_to(repo)): sha256(path)
            for path in [*python_sources, *shell_scripts]
        },
        "environment_sha256": {
            relative: sha256(repo / relative)
            for relative in (
                "envs/hurdler.yml",
                "envs/hurdler-linux-64.lock",
                "envs/hurdler-pip.lock",
            )
        },
        "pytest": pytest_output,
        "cli_help": cli_help,
        "kernel_smoke": "passed",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
