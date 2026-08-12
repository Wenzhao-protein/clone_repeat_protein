"""Compatibility entry point for the maintained reference/index pipeline."""

from __future__ import annotations

from .cli import main as cli_main


def main() -> int:
    return cli_main(["lookup", "build"])


if __name__ == "__main__":
    raise SystemExit(main())
