"""Build an auditable designed-protein structure inventory and AF3 recovery plan."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

from .constants import validate_protein_sequence
from .structural_repeats import DESIGNED_CORPUS_VERSION

AF3_CONTAINER = Path(
    "/net/software/containers/versions/af3/af3-open_20260729.sif"
)
AF3_WRAPPER = Path(
    "/home/wendai/projects/method_dev/pipelines/src/blocks/scripts/run_alphafold_AF.py"
)


def _read_table(path: str | Path) -> pd.DataFrame:
    source = Path(path)
    return pd.read_parquet(source) if source.suffix == ".parquet" else pd.read_csv(source)


def _identifier(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _structure_aliases(path: Path) -> set[str]:
    stem = path.stem
    aliases = {
        _identifier(stem),
        _identifier(
            re.sub(
                r"(?i)(?:^|_)(design|xtal|cryo|fit|model|relaxed|rank_?\d+)(?:_|$)",
                "_",
                stem,
            )
        ),
    }
    for suffix in (
        "_design",
        "_xtal",
        "_cryo_fit",
        "_model",
        "_relaxed",
        "_rank_001",
    ):
        if suffix in stem.lower():
            aliases.add(_identifier(stem[: stem.lower().index(suffix)]))
    return {value for value in aliases if value}


def build_designed_structure_inventory(
    catalog_path: str | Path,
    structure_roots: Iterable[str | Path],
    output_path: str | Path,
    *,
    af3_output_root: str | Path,
    af3_task_path: str | Path,
) -> pd.DataFrame:
    """Match supplied structures and emit missing-only, fixed-seed AF3 tasks."""
    catalog = _read_table(catalog_path)
    required = {"module_id", "full_sequence", "source_accession"}
    missing = sorted(required - set(catalog.columns))
    if missing:
        raise ValueError(f"Designed catalog is missing columns: {missing}")
    if "collection" in catalog:
        catalog = catalog.loc[
            catalog.collection.astype(str).str.startswith("designed")
        ].copy()
    if catalog.empty:
        raise ValueError("Designed catalog contains no designed collection rows")
    structure_paths: list[Path] = []
    for root_value in structure_roots:
        root = Path(root_value)
        if root.is_file() and root.suffix.lower() in {".pdb", ".cif", ".mmcif"}:
            structure_paths.append(root.resolve())
        elif root.is_dir():
            structure_paths.extend(
                path.resolve()
                for path in root.rglob("*")
                if path.is_file()
                and path.suffix.lower() in {".pdb", ".cif", ".mmcif"}
            )
    by_alias: dict[str, list[Path]] = {}
    for path in sorted(set(structure_paths)):
        for alias in _structure_aliases(path):
            by_alias.setdefault(alias, []).append(path)

    af3_root = Path(af3_output_root).resolve()
    af3_input_root = af3_root.parent / "af3_inputs"
    af3_input_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for source in catalog.to_dict(orient="records"):
        sequence = validate_protein_sequence(str(source["full_sequence"]))
        identifiers = {
            _identifier(source["module_id"]),
            _identifier(source["source_accession"]),
            _identifier(str(source["module_id"]).removeprefix("designed_")),
        }
        matches = sorted(
            {
                path
                for identifier in identifiers
                for path in by_alias.get(identifier, [])
            },
            key=lambda value: ("xtal" not in value.stem.lower(), str(value)),
        )
        module_id = str(source["module_id"])
        af3_directory = af3_root / module_id
        row = {
            **source,
            "corpus_version": DESIGNED_CORPUS_VERSION,
            "full_sequence": sequence,
            "full_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "full_sequence_length": len(sequence),
            "author_structure_path": str(matches[0]) if matches else "",
            "structure_candidate_paths_json": json.dumps(
                [str(path) for path in matches]
            ),
            "structure_inventory_status": (
                "author_or_pdb_structure_available"
                if matches
                else "missing_structure_af3_requested"
            ),
            "af3_structure_path": str(af3_directory / f"{module_id}_model.cif"),
            "af3_output_dir": str(af3_directory),
            "af3_input_json": str(af3_input_root / f"{module_id}.json"),
            "af3_seed": 42,
            "af3_diffusion_samples": 1,
        }
        rows.append(row)
    inventory = pd.DataFrame(rows).sort_values("module_id", kind="mergesort")
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    inventory.to_parquet(destination, index=False)
    inventory.to_csv(destination.with_suffix(".csv"), index=False)

    missing_rows = inventory.loc[
        inventory.structure_inventory_status.eq("missing_structure_af3_requested")
    ].copy()
    task_path = Path(af3_task_path)
    task_path.parent.mkdir(parents=True, exist_ok=True)
    commands = []
    for row in missing_rows.itertuples(index=False):
        input_json = Path(row.af3_input_json)
        input_json.write_text(
            json.dumps(
                {
                    "name": row.module_id,
                    "modelSeeds": [42],
                    "sequences": [
                        {
                            "protein": {
                                "id": "A",
                                "sequence": row.full_sequence,
                                "unpairedMsa": "",
                                "pairedMsa": "",
                                "templates": [],
                            }
                        }
                    ],
                    "dialect": "alphafold3",
                    "version": 1,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        command = " ".join(
            [
                "apptainer exec --nv",
                "--bind /home:/home",
                "--bind /net/scratch:/net/scratch",
                "--bind /net/databases:/databases",
                str(AF3_CONTAINER),
                "python",
                str(AF3_WRAPPER),
                f"--json_path={input_json}",
                "--num_diffusion_samples=1",
                "--run_data_pipeline=False",
                "--run_inference=True",
                "--skip_existing=True",
                "--flash_attention_implementation=xla",
                "--db_dir=/databases",
                "--model_dir=/databases/alphafold",
                f"--jax_compilation_cache_dir={af3_root / 'jax_cache'}",
                f"--output_dir={af3_root}",
            ]
        )
        if "\n" in command:
            raise AssertionError("AF3 task command must occupy exactly one line")
        commands.append(command)
    task_path.write_text("\n".join(commands) + ("\n" if commands else ""))
    if commands:
        longest_position = int(
            missing_rows.full_sequence_length.astype(int).to_numpy().argmax()
        )
        task_path.with_name("smoke_tasks.txt").write_text(
            commands[longest_position] + "\n"
        )
    missing_manifest = task_path.with_name("af3_missing_manifest.csv")
    missing_rows[
        [
            "module_id",
            "family",
            "source_accession",
            "full_sequence_length",
            "full_sequence_sha256",
            "af3_seed",
            "af3_diffusion_samples",
            "af3_output_dir",
            "af3_input_json",
            "af3_structure_path",
        ]
    ].to_csv(missing_manifest, index=False)
    task_index = pd.DataFrame(
        {
            "task_index": range(len(missing_rows)),
            "module_id": missing_rows.module_id.tolist(),
            "expected_output": missing_rows.af3_structure_path.tolist(),
            "status": "missing",
        }
    )
    task_index.to_csv(task_path.with_name("task_index.csv"), index=False)
    destination.with_suffix(".manifest.json").write_text(
        json.dumps(
            {
                "corpus_version": DESIGNED_CORPUS_VERSION,
                "input_catalog": str(Path(catalog_path).resolve()),
                "structure_roots": [str(Path(root).resolve()) for root in structure_roots],
                "input_rows": len(inventory),
                "matched_structure_rows": int(
                    inventory.structure_inventory_status.eq(
                        "author_or_pdb_structure_available"
                    ).sum()
                ),
                "af3_missing_rows": len(missing_rows),
                "af3_seed": 42,
                "af3_diffusion_samples": 1,
                "af3_container": str(AF3_CONTAINER),
                "af3_wrapper": str(AF3_WRAPPER),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return inventory.reset_index(drop=True)


def validate_af3_outputs(
    inventory_path: str | Path, output_path: str | Path
) -> pd.DataFrame:
    """Validate AF3 CIF sequence identity and retain confidence provenance."""
    import biotite.structure as struc
    from biotite.structure.io import load_structure
    from Bio.PDB.Polypeptide import protein_letters_3to1_extended

    inventory = _read_table(inventory_path)
    rows = []
    for source in inventory.loc[
        inventory.structure_inventory_status.eq("missing_structure_af3_requested")
    ].to_dict(orient="records"):
        expected = validate_protein_sequence(str(source["full_sequence"]))
        structure_path = Path(str(source["af3_structure_path"]))
        payload = {
            "module_id": source["module_id"],
            "expected_sequence_sha256": hashlib.sha256(expected.encode()).hexdigest(),
            "expected_length": len(expected),
            "af3_structure_path": str(structure_path),
            "af3_structure_exists": structure_path.is_file(),
            "af3_seed": int(source.get("af3_seed", 42)),
            "af3_diffusion_samples": int(source.get("af3_diffusion_samples", 1)),
        }
        if not structure_path.is_file():
            payload.update(
                af3_validation_status="missing_structure",
                observed_sequence_sha256="",
                observed_length=0,
                exact_sequence_match=False,
                ranking_score=None,
                ptm=None,
                confidence_file="",
            )
            rows.append(payload)
            continue
        try:
            atoms = load_structure(str(structure_path), model=1)
            atoms = atoms[struc.filter_amino_acids(atoms)]
            chain_ids = sorted(set(str(value) for value in atoms.chain_id))
            chain_id = max(
                chain_ids,
                key=lambda value: struc.get_residue_count(
                    atoms[atoms.chain_id == value]
                ),
            )
            _, residue_names = struc.get_residues(
                atoms[atoms.chain_id == chain_id]
            )
            observed = validate_protein_sequence(
                "".join(
                    protein_letters_3to1_extended.get(str(name).upper(), "X")
                    for name in residue_names
                )
            )
            confidence_file = structure_path.with_name(
                structure_path.name.replace("_model.cif", "_summary_confidences.json")
            )
            confidence = (
                json.loads(confidence_file.read_text())
                if confidence_file.is_file()
                else {}
            )
            exact = observed == expected
            payload.update(
                af3_validation_status=(
                    "passed_exact_sequence" if exact else "failed_sequence_mismatch"
                ),
                observed_sequence_sha256=hashlib.sha256(observed.encode()).hexdigest(),
                observed_length=len(observed),
                exact_sequence_match=exact,
                structure_chain=chain_id,
                structure_sha256=hashlib.sha256(structure_path.read_bytes()).hexdigest(),
                ranking_score=confidence.get("ranking_score"),
                ptm=confidence.get("ptm"),
                confidence_file=str(confidence_file) if confidence_file.is_file() else "",
                confidence_sha256=(
                    hashlib.sha256(confidence_file.read_bytes()).hexdigest()
                    if confidence_file.is_file()
                    else ""
                ),
            )
        except Exception as exc:
            payload.update(
                af3_validation_status="failed_parse",
                observed_sequence_sha256="",
                observed_length=0,
                exact_sequence_match=False,
                ranking_score=None,
                ptm=None,
                confidence_file="",
                validation_error=f"{type(exc).__name__}: {exc}",
            )
        rows.append(payload)
    result = pd.DataFrame(rows).sort_values("module_id", kind="mergesort")
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(destination, index=False)
    result.to_csv(destination.with_suffix(".csv"), index=False)
    destination.with_suffix(".manifest.json").write_text(
        json.dumps(
            {
                "corpus_version": DESIGNED_CORPUS_VERSION,
                "inventory": str(Path(inventory_path).resolve()),
                "rows": len(result),
                "status_counts": result.af3_validation_status.value_counts().to_dict(),
                "all_exact_when_present": bool(
                    result.loc[result.af3_structure_exists, "exact_sequence_match"].all()
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return result.reset_index(drop=True)
