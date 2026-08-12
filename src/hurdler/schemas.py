"""Small typed records used at public HURDLER boundaries."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class MatchResult:
    module: str
    effective_module: str
    original_length: int
    effective_length: int
    expansion_copies: int
    plasmid: str
    success: bool
    site_i_3mer: str | None = None
    site_ii_3mer: str | None = None
    direction: str | None = None
    site_i_position: int | None = None
    site_ii_position: int | None = None
    pattern_key: int | None = None
    solution_count: int = 0
    best_pair_id: int | None = None

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ModuleRecord:
    module_id: str
    collection: str
    family: str
    unit_sequence: str
    evidence_tier: str
    source_name: str
    source_url: str
    source_accession: str = ""
    source_chain: str = ""
    unit_start: int | None = None
    unit_end: int | None = None
    reviewed: bool | None = None
    boundary_method: str = ""
    retrieved_date: str = ""
    download_date: str = ""
    source_sha256: str = ""
    license_name: str = ""
    citation: str = ""
    notes: str = ""
    full_sequence: str = ""
    full_sequence_origin: int = 1
    full_sequence_sha256: str = ""
    source_annotation_id: str = ""

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["download_date"] = self.download_date or self.retrieved_date
        payload["unit_length"] = len(self.unit_sequence)
        return payload
