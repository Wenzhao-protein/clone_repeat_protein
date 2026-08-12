"""Versioned scientific rule profiles."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from .constants import DEFAULT_RANDOM_SEED, RULE_PROFILE_NAME, SCHEMA_VERSION


@dataclass(frozen=True)
class RuleProfile:
    name: str = RULE_PROFILE_NAME
    schema_version: int = SCHEMA_VERSION
    orthogonality_threshold: int = 1
    missing_fidelity_is_compatible: bool = True
    require_signed_site_ii_iii_overhang_match: bool = True
    plasmid_sites: tuple[str, ...] = ("site_i", "site_ii")
    minimum_start_distance: int = 5
    minimum_is_inclusive: bool = True
    maximum_distance_is_module_length_exclusive: bool = True
    double_module_before_matching: bool = True
    random_seed: int = DEFAULT_RANDOM_SEED
    description: str = "Frozen behavior of the latest optimized HURDLER notebook."

    def distance_is_valid(self, distance: int, module_length: int) -> bool:
        lower = distance >= self.minimum_start_distance
        upper = distance < module_length
        return lower and upper

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["plasmid_sites"] = list(self.plasmid_sites)
        return payload


LEGACY_OPTIMIZED_V1 = RuleProfile()


def load_rule_profile(path: str | Path | None = None) -> RuleProfile:
    if path is None:
        return LEGACY_OPTIMIZED_V1
    payload = json.loads(Path(path).read_text())
    if payload.get("name") != RULE_PROFILE_NAME:
        raise ValueError(f"Unsupported rule profile: {payload.get('name')!r}")
    payload["plasmid_sites"] = tuple(payload.get("plasmid_sites", ("site_i", "site_ii")))
    return RuleProfile(**payload)
