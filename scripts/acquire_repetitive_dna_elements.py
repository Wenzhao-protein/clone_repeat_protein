#!/usr/bin/env python3
"""Acquire traceable public DNA/RNA elements for exact-DNA array tests.

Raw downloads belong in the study scratch mirror.  Compact inventories and
explicit exclusions may be promoted into the study tree.  This script never
uses IDT credentials and never sends sequences to IDT.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from hurdler.dna_assembly import expand_element_inventory, validate_dna  # noqa: E402


CRISPR_URL = (
    "https://crisprcas.i2bc.paris-saclay.fr/"
    "Home/DownloadFile?filename=dr_34.zip"
)
APTAMER_URL = "https://aptamer.ribocentre.org/apidata/sequences_cleaned.json"
RFAM_TEMPLATE = "https://rfam.org/family/{accession}/alignment/fastau"
RFAM_REGULATORY_FAMILIES = (
    "RF00050",  # FMN riboswitch
    "RF00059",  # TPP riboswitch
    "RF00162",  # SAM riboswitch
    "RF00167",  # purine riboswitch
    "RF00168",  # lysine riboswitch
    "RF00234",  # glmS ribozyme
    "RF00379",  # M-box riboswitch
    "RF00504",  # glycine riboswitch
    "RF00522",  # preQ1 riboswitch
    "RF01051",  # fluoride riboswitch
    "RF01734",  # c-di-GMP-I riboswitch
    "RF01750",  # ZTP riboswitch
    "RF01786",  # c-di-GMP-II riboswitch
    "RF01831",  # c-di-AMP riboswitch
    "RF03072",  # guanidine-I riboswitch
)


def session() -> requests.Session:
    client = requests.Session()
    retry = Retry(
        total=4,
        connect=4,
        read=4,
        status=4,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
        respect_retry_after_header=True,
    )
    client.mount("https://", HTTPAdapter(max_retries=retry))
    client.headers["User-Agent"] = "HURDLER academic reproducibility workflow"
    return client


def fetch(client: requests.Session, url: str, cache_path: Path) -> bytes:
    if cache_path.is_file() and cache_path.stat().st_size:
        return cache_path.read_bytes()
    response = client.get(url, timeout=120)
    response.raise_for_status()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_bytes(response.content)
    return response.content


def normalize_public_sequence(value: object) -> str:
    raw = re.sub(r"[\s.-]", "", str(value)).upper().replace("U", "T")
    return validate_dna(raw)


def parse_crispr(payload: bytes) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    exclusions: list[dict[str, object]] = []
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        fasta = archive.read("direct_repeat_id.fsa").decode()
    for record in SeqIO.parse(io.StringIO(fasta), "fasta"):
        try:
            sequence = normalize_public_sequence(record.seq)
        except ValueError as exc:
            exclusions.append(
                {"source_database": "CRISPRCasdb", "source_accession": record.id, "reason": str(exc)}
            )
            continue
        rows.append(
            {
                "element_id": f"CRISPRCasdb:{record.id}",
                "element_sequence": sequence,
                "source_database": "CRISPRCasdb",
                "source_accession": record.id,
                "source_url": CRISPR_URL,
                "notes": "CRISPRCasdb direct-repeat entry",
            }
        )
    return rows, exclusions


def parse_aptamers(payload: bytes) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    body = json.loads(payload)
    entries = body.get("Sheet1", body) if isinstance(body, dict) else body
    rows: list[dict[str, object]] = []
    exclusions: list[dict[str, object]] = []
    for entry in entries:
        accession = str(entry.get("ID", "missing_id"))
        try:
            sequence = normalize_public_sequence(entry.get("Sequence", ""))
        except ValueError as exc:
            exclusions.append(
                {"source_database": "Ribocentre_Aptamer", "source_accession": accession, "reason": str(exc)}
            )
            continue
        rows.append(
            {
                "element_id": f"Ribocentre_Aptamer:{accession}",
                "element_sequence": sequence,
                "source_database": "Ribocentre_Aptamer",
                "source_accession": accession,
                "source_url": APTAMER_URL,
                "notes": str(entry.get("Link to PubMed Entry", "") or "curated aptamer sequence"),
            }
        )
    return rows, exclusions


def parse_rfam(
    client: requests.Session,
    cache_dir: Path,
    families: tuple[str, ...],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    rows: list[dict[str, object]] = []
    exclusions: list[dict[str, object]] = []
    raw_records: list[dict[str, object]] = []
    downloads: list[dict[str, object]] = []
    for accession in families:
        url = RFAM_TEMPLATE.format(accession=accession)
        try:
            payload = fetch(client, url, cache_dir / "rfam" / f"{accession}.fasta")
        except requests.RequestException:
            exclusions.append(
                {"source_database": "Rfam", "source_accession": accession, "reason": "download_failed"}
            )
            continue
        downloads.append(
            {
                "source": "Rfam",
                "accession": accession,
                "url": url,
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
        )
        records = list(SeqIO.parse(io.StringIO(payload.decode()), "fasta"))
        if not records:
            exclusions.append(
                {"source_database": "Rfam", "source_accession": accession, "reason": "empty_seed_alignment"}
            )
            continue
        # All seed sequences remain in source mappings; one exact central seed
        # is the experiment element so very large families do not dominate it.
        selected = records[(len(records) - 1) // 2]
        for record in records:
            raw_records.append(
                {
                    "source_database": "Rfam",
                    "family_accession": accession,
                    "sequence_accession": record.id,
                    "selected_for_experiment": record.id == selected.id,
                }
            )
        try:
            sequence = normalize_public_sequence(selected.seq)
        except ValueError as exc:
            exclusions.append(
                {"source_database": "Rfam", "source_accession": accession, "reason": str(exc)}
            )
            continue
        rows.append(
            {
                "element_id": f"Rfam:{accession}:{selected.id}",
                "element_sequence": sequence,
                "source_database": "Rfam",
                "source_accession": f"{accession}:{selected.id}",
                "source_url": url,
                "notes": "earlier-middle sequence in the official Rfam seed alignment",
            }
        )
    return rows, exclusions, raw_records, downloads


def deduplicate_elements(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = frame.copy()
    frame["element_sha256"] = frame.element_sequence.map(
        lambda sequence: hashlib.sha256(str(sequence).encode()).hexdigest()
    )
    mappings = frame.copy()
    selected = (
        frame.sort_values(["source_database", "element_id"])
        .drop_duplicates(["source_database", "element_sha256"])
        .reset_index(drop=True)
    )
    selected["representative_source_element_id"] = selected["element_id"]
    selected["element_id"] = selected.apply(
        lambda row: f"{row['source_database']}:{row['element_sha256'][:16]}",
        axis=1,
    )
    selected["element_length_bp"] = selected.element_sequence.str.len()
    return selected, mappings


def write_table(frame: pd.DataFrame, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(stem.with_suffix(".parquet"), index=False)
    frame.to_csv(stem.with_suffix(".csv"), index=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--skip-crispr", action="store_true")
    parser.add_argument("--skip-aptamer", action="store_true")
    parser.add_argument("--skip-rfam", action="store_true")
    parser.add_argument("--rfam-family", action="append", default=[])
    args = parser.parse_args()

    client = session()
    elements: list[dict[str, object]] = []
    exclusions: list[dict[str, object]] = []
    rfam_mappings: list[dict[str, object]] = []
    raw_manifest: list[dict[str, object]] = []
    if not args.skip_crispr:
        payload = fetch(client, CRISPR_URL, args.cache_dir / "crisprcasdb" / "dr_34.zip")
        rows, rejected = parse_crispr(payload)
        elements.extend(rows)
        exclusions.extend(rejected)
        raw_manifest.append({"source": "CRISPRCasdb", "url": CRISPR_URL, "sha256": hashlib.sha256(payload).hexdigest()})
    if not args.skip_aptamer:
        payload = fetch(client, APTAMER_URL, args.cache_dir / "ribocentre" / "sequences_cleaned.json")
        rows, rejected = parse_aptamers(payload)
        elements.extend(rows)
        exclusions.extend(rejected)
        raw_manifest.append({"source": "Ribocentre_Aptamer", "url": APTAMER_URL, "sha256": hashlib.sha256(payload).hexdigest()})
    if not args.skip_rfam:
        families = tuple(args.rfam_family) or RFAM_REGULATORY_FAMILIES
        rows, rejected, mappings, downloads = parse_rfam(
            client, args.cache_dir, families
        )
        elements.extend(rows)
        exclusions.extend(rejected)
        rfam_mappings.extend(mappings)
        raw_manifest.extend(downloads)

    if not elements:
        raise RuntimeError("No public elements were acquired; see source exclusions")
    inventory, mappings = deduplicate_elements(pd.DataFrame(elements))
    targets = expand_element_inventory(inventory)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_table(inventory, args.output_dir / "public_element_inventory")
    write_table(mappings, args.output_dir / "public_element_source_mappings")
    exclusion_frame = pd.DataFrame(
        exclusions,
        columns=["source_database", "source_accession", "reason"],
    )
    rfam_mapping_frame = pd.DataFrame(
        rfam_mappings,
        columns=[
            "source_database",
            "family_accession",
            "sequence_accession",
            "selected_for_experiment",
        ],
    )
    write_table(exclusion_frame, args.output_dir / "public_element_exclusions")
    write_table(rfam_mapping_frame, args.output_dir / "rfam_seed_source_mappings")
    write_table(targets, args.output_dir / "real_element_derived_targets")
    manifest = {
        "version": "arbitrary-dna-active-latent-v1",
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "source_downloads": raw_manifest,
        "source_element_rows": len(mappings),
        "unique_element_rows": len(inventory),
        "derived_target_rows": len(targets),
        "exclusion_rows": len(exclusions),
        "source_counts": inventory.source_database.value_counts().to_dict(),
        "copy_counts": [2, 4, 8, 16, 32],
    }
    (args.output_dir / "public_element_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
