"""Compatibility entry point for run validation."""

from __future__ import annotations

from .cli import main as cli_main


def main() -> int:
    return cli_main(["validate-run"])


if __name__ == "__main__":
    raise SystemExit(main())
