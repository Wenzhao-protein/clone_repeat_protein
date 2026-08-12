"""Compatibility entry point for deterministic success-rate analysis."""

from .rate import legacy_random_modules, run_success_rate

__all__ = ["legacy_random_modules", "run_success_rate"]
