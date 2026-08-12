"""Stable constants shared by HURDLER workflows."""

from __future__ import annotations

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INT = {aa: index for index, aa in enumerate(AMINO_ACIDS)}
INT_TO_AA = dict(enumerate(AMINO_ACIDS))
THREE_MER_SPACE = len(AMINO_ACIDS) ** 3

PLASMIDS = (
    "pGEX-4T-1",
    "pMAL-c5X",
    "pET-21a(+)",
    "pET-28a(+)",
    "pET-28a(+)_start_codon",
    "pCold_I",
    "pUC18",
    "pQE-3",
)

RULE_PROFILE_NAME = "legacy-optimized-v1"
SCHEMA_VERSION = 1
DEFAULT_RANDOM_SEED = 42
DEFAULT_FRAGMENT_LIMITS_BP = (1800, 3000)


def validate_protein_sequence(sequence: str, *, allow_empty: bool = False) -> str:
    """Normalize and validate a sequence against the 20-residue alphabet."""
    normalized = "".join(sequence.split()).upper()
    if not normalized and not allow_empty:
        raise ValueError("Protein sequence is empty")
    invalid = sorted(set(normalized) - set(AMINO_ACIDS))
    if invalid:
        raise ValueError(f"Unsupported amino acids: {''.join(invalid)}")
    return normalized
