from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from hurdler.data_deposition import (
    HIDDEN_METADATA_NOTE,
    build_directory_structure,
)


REPO = Path(__file__).resolve().parents[1]
DEPOSITION = REPO / "constructs_and_sequencing_result"


def test_deposition_directory_structure_is_current():
    stored = (DEPOSITION / "directory_structure.txt").read_text(encoding="utf-8")
    assert stored == build_directory_structure(DEPOSITION)
    assert stored.splitlines()[0] == HIDDEN_METADATA_NOTE
    assert "\uf026" not in stored and "#Uf026" not in stored
    assert ".clc" not in "\n".join(stored.splitlines()[2:])
    assert "README.md" not in "\n".join(stored.splitlines()[2:])


def test_deposition_paths_are_ascii_and_design_genbanks_parse():
    non_ascii = [
        str(path.relative_to(DEPOSITION))
        for path in DEPOSITION.rglob("*")
        if not path.name.isascii()
    ]
    assert non_ascii == []
    design_files = sorted((DEPOSITION / "extended_armrp/1_fragment").glob("*.gbk"))
    design_files += sorted((DEPOSITION / "extended_armrp/2_plasmid_construct").glob("*.gbk"))
    design_files += sorted((DEPOSITION / "single_chain_armrp_assembly/1_fragment").glob("*.gbk"))
    design_files += sorted((DEPOSITION / "single_chain_armrp_assembly/2_plasmid_construct").glob("*.gbk"))
    assert len(design_files) == 27
    for path in design_files:
        record = SeqIO.read(path, "genbank")
        assert len(record.seq) > 0, path


def test_deposition_retains_repository_only_clc_evidence():
    triangle = DEPOSITION / "single_chain_armrp_assembly/3_sequencing_result/dArmRP_triangle_step6"
    assert {".acl", ".clcinfo", ".orderlist2"}.issubset(
        {path.name for path in triangle.iterdir()}
    )
    assert len(list(triangle.glob("*.clc"))) == 9
