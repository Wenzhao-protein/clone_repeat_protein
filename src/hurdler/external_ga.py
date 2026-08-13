"""Credential-safe, reproducible local/Slurm GA run bundles."""

from __future__ import annotations

import hashlib
import json
import re
import shlex
import tempfile
import zipfile
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .idt import IDT_SCORE_POLICY
from .vector_design import DesignRequestV2


DEFAULT_REPOSITORY_URL = "https://github.com/Wenzhao-protein/clone_repeat_protein.git"
_SAFE_SLURM_VALUE = re.compile(r"^[A-Za-z0-9_.-]+$")
_SAFE_COMMIT = re.compile(r"^[0-9a-f]{40}$")
_SAFE_WALLTIME = re.compile(r"^(?:[0-9]+-)?[0-9]{1,2}:[0-9]{2}:[0-9]{2}$")


def _sha_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _safe_optional(value: str, name: str) -> str:
    value = str(value).strip()
    if value and not _SAFE_SLURM_VALUE.fullmatch(value):
        raise ValueError(f"{name} contains unsafe characters")
    return value


@dataclass(frozen=True)
class ExternalGAResources:
    """Structured compute settings; arbitrary scheduler directives are forbidden."""

    worker_cpus: int = 16
    memory_gb: int = 32
    walltime: str = "24:00:00"
    partition: str = "cpu"
    account: str = ""
    qos: str = ""
    constraint: str = ""
    conda_environment: str = "hurdler"
    result_directory: str = "results"

    def __post_init__(self) -> None:
        if isinstance(self.worker_cpus, bool) or not 1 <= int(self.worker_cpus) <= 1024:
            raise ValueError("worker_cpus must be an integer between 1 and 1024")
        if isinstance(self.memory_gb, bool) or not 1 <= int(self.memory_gb) <= 1_048_576:
            raise ValueError("memory_gb must be a positive integer")
        if not _SAFE_WALLTIME.fullmatch(str(self.walltime)):
            raise ValueError("walltime must use HH:MM:SS or D-HH:MM:SS")
        for name in ("partition", "account", "qos", "constraint", "conda_environment"):
            _safe_optional(getattr(self, name), name)
        if not self.partition:
            raise ValueError("partition cannot be empty")
        if not self.conda_environment:
            raise ValueError("conda_environment cannot be empty")
        if not str(self.result_directory).strip():
            raise ValueError("result_directory cannot be empty")


def _slurm_header(resources: ExternalGAResources) -> str:
    directives = [
        "#SBATCH --job-name=hurdler-ga",
        f"#SBATCH --cpus-per-task={int(resources.worker_cpus)}",
        f"#SBATCH --mem={int(resources.memory_gb)}G",
        f"#SBATCH --time={resources.walltime}",
        f"#SBATCH --partition={resources.partition}",
        "#SBATCH --output=slurm-%j.out",
        "#SBATCH --error=slurm-%j.err",
    ]
    if resources.account:
        directives.append(f"#SBATCH --account={resources.account}")
    if resources.qos:
        directives.append(f"#SBATCH --qos={resources.qos}")
    if resources.constraint:
        directives.append(f"#SBATCH --constraint={resources.constraint}")
    return "\n".join(directives)


def _run_script(
    *,
    request_json: str,
    request_sha256: str,
    environment_sha256: str,
    repository_url: str,
    repository_commit: str,
    resources: ExternalGAResources,
    credential_path: str,
    auth_method: str,
    validation_mode: str,
) -> str:
    idt_arguments = ""
    if validation_mode == "api":
        idt_arguments = ' --idt-credential-file "$IDT_CREDENTIAL_FILE"'
        if auth_method != "auto":
            idt_arguments += f" --auth-method {shlex.quote(auth_method)}"
    return f'''#!/usr/bin/env bash
{_slurm_header(resources)}
set -euo pipefail

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

BUNDLE_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REQUEST_FILE="$BUNDLE_DIR/request.json"
ENVIRONMENT_FILE="$BUNDLE_DIR/envs/hurdler.yml"
REQUEST_SHA256={shlex.quote(request_sha256)}
ENVIRONMENT_SHA256={shlex.quote(environment_sha256)}
REPOSITORY_URL={shlex.quote(repository_url)}
REPOSITORY_COMMIT={shlex.quote(repository_commit)}
CONDA_ENVIRONMENT={shlex.quote(resources.conda_environment)}
GA_WORKERS={int(resources.worker_cpus)}
REQUESTED_MEMORY_GB={int(resources.memory_gb)}
RESULT_DIRECTORY_SETTING={shlex.quote(resources.result_directory)}
IDT_CREDENTIAL_FILE={shlex.quote(credential_path)}
IDT_AUTH_METHOD={shlex.quote(auth_method)}

case "$RESULT_DIRECTORY_SETTING" in
  /*) RESULT_DIR="$RESULT_DIRECTORY_SETTING" ;;
  *) RESULT_DIR="$BUNDLE_DIR/$RESULT_DIRECTORY_SETTING" ;;
esac
case "$IDT_CREDENTIAL_FILE" in
  '~/'*) IDT_CREDENTIAL_FILE="$HOME/${{IDT_CREDENTIAL_FILE#\\~/}}" ;;
esac
CHECKOUT_DIR="${{HURDLER_REPO_DIR:-$BUNDLE_DIR/work/repository}}"

sha256_of() {{
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{{print $1}}';
  else python - "$1" <<'PY'
import hashlib, pathlib, sys
print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())
PY
  fi
}}

write_embedded_request() {{
  cat <<'HURDLER_REQUEST_JSON'
{request_json.rstrip()}
HURDLER_REQUEST_JSON
}}

activate_environment() {{
  if ! command -v conda >/dev/null 2>&1; then
    echo "conda is unavailable; run 'bash run_ga.sh setup' after installing Miniconda/Mambaforge" >&2
    return 2
  fi
  eval "$(conda shell.bash hook)"
  conda activate "$CONDA_ENVIRONMENT"
}}

verify_bundle() {{
  test -f "$REQUEST_FILE" || {{ echo "Missing $REQUEST_FILE" >&2; return 2; }}
  test -f "$ENVIRONMENT_FILE" || {{ echo "Missing $ENVIRONMENT_FILE" >&2; return 2; }}
  test "$(sha256_of "$REQUEST_FILE")" = "$REQUEST_SHA256" || {{ echo "request.json SHA256 mismatch" >&2; return 2; }}
  test "$(sha256_of "$ENVIRONMENT_FILE")" = "$ENVIRONMENT_SHA256" || {{ echo "environment YAML SHA256 mismatch" >&2; return 2; }}
  local embedded
  embedded="$(mktemp)"
  write_embedded_request > "$embedded"
  test "$(sha256_of "$embedded")" = "$REQUEST_SHA256" || {{ rm -f "$embedded"; echo "embedded request SHA256 mismatch" >&2; return 2; }}
  cmp -s "$embedded" "$REQUEST_FILE" || {{ rm -f "$embedded"; echo "embedded request differs from request.json" >&2; return 2; }}
  rm -f "$embedded"
}}

verify_cpu_allocation() {{
  local available
  if test -n "${{SLURM_CPUS_PER_TASK:-}}"; then available="$SLURM_CPUS_PER_TASK"
  else available="$(python - <<'PY'
import os
print(os.cpu_count() or 1)
PY
)"; fi
  test "$available" -ge "$GA_WORKERS" || {{ echo "Allocated CPU count $available is below ga_workers=$GA_WORKERS" >&2; return 2; }}
}}

setup_environment() {{
  verify_bundle
  command -v git >/dev/null 2>&1 || {{ echo "git is required" >&2; return 2; }}
  command -v conda >/dev/null 2>&1 || {{ echo "conda is required" >&2; return 2; }}
  mkdir -p "$(dirname "$CHECKOUT_DIR")"
  if test -d "$CHECKOUT_DIR/.git"; then
    git -C "$CHECKOUT_DIR" fetch origin "$REPOSITORY_COMMIT"
  else
    git clone "$REPOSITORY_URL" "$CHECKOUT_DIR"
  fi
  git -C "$CHECKOUT_DIR" checkout --detach "$REPOSITORY_COMMIT"
  eval "$(conda shell.bash hook)"
  if conda env list --json | python -c 'import json,sys; n=sys.argv[1]; raise SystemExit(0 if any(p.rsplit("/",1)[-1] == n for p in json.load(sys.stdin)["envs"]) else 1)' "$CONDA_ENVIRONMENT"; then
    conda env update -n "$CONDA_ENVIRONMENT" -f "$ENVIRONMENT_FILE" --prune
  else
    conda env create -n "$CONDA_ENVIRONMENT" -f "$ENVIRONMENT_FILE"
  fi
  conda activate "$CONDA_ENVIRONMENT"
  python -m pip install -e "$CHECKOUT_DIR"
  echo "Environment $CONDA_ENVIRONMENT is ready at commit $REPOSITORY_COMMIT"
}}

preflight() {{
  verify_bundle
  activate_environment
  verify_cpu_allocation
  test -d "$CHECKOUT_DIR/.git" || {{ echo "Checkout missing: $CHECKOUT_DIR; run setup or set HURDLER_REPO_DIR" >&2; return 2; }}
  test "$(git -C "$CHECKOUT_DIR" rev-parse HEAD)" = "$REPOSITORY_COMMIT" || {{ echo "Checkout is not at frozen commit $REPOSITORY_COMMIT" >&2; return 2; }}
  mkdir -p "$RESULT_DIR"
  test -w "$RESULT_DIR" || {{ echo "Result directory is not writable: $RESULT_DIR" >&2; return 2; }}
  if test "$REQUESTED_MEMORY_GB" -lt 1; then echo "Invalid requested memory" >&2; return 2; fi
  if test -z "${{SLURM_JOB_ID:-}}" && test -r /proc/meminfo; then
    local available_memory_gb
    available_memory_gb="$(awk '/MemAvailable:/ {{printf "%d", $2 / 1024 / 1024}}' /proc/meminfo)"
    if test "$available_memory_gb" -lt "$REQUESTED_MEMORY_GB"; then
      echo "WARNING: approximately $available_memory_gb GB RAM is currently available; this bundle requests $REQUESTED_MEMORY_GB GB." >&2
    else
      echo "Local memory preflight: approximately $available_memory_gb GB available for a $REQUESTED_MEMORY_GB GB request."
    fi
  fi
'''+(f'''  hurdler idt-preflight --idt-credential-file "$IDT_CREDENTIAL_FILE"{'' if auth_method == 'auto' else ' --auth-method ' + shlex.quote(auth_method)}
''' if validation_mode == "api" else '''  echo "Batch/compatibility mode: IDT credential and network preflight omitted."
''')+'''  echo "Preflight passed."
}

run_design() {
  preflight
  local started output checkpoint progress archives
  started="$(date -u +%Y%m%dT%H%M%SZ)"
  output="$RESULT_DIR/run_$started"
  checkpoint="$RESULT_DIR/checkpoint_latest.zip"
  progress="$RESULT_DIR/progress_$started.jsonl"
  archives="$RESULT_DIR/archives"
  mkdir -p "$output" "$archives"
  set +e
  hurdler design-construct \
    --request "$REQUEST_FILE" \
    --output-dir "$output" \
    --progress-jsonl "$progress" \
    --checkpoint-zip "$checkpoint" \
    --checkpoint-interval-seconds 180 \
    --final-archive-dir "$archives" \
    --fail-on-nonaccepted'''+idt_arguments+'''
  local code=$?
  set -e
  echo "HURDLER exit code: $code; output: $output"
  return "$code"
}

submit_job() {
  command -v sbatch >/dev/null 2>&1 || { echo "sbatch is unavailable" >&2; return 2; }
  verify_bundle
  sbatch "$0" __slurm_run
}

show_status() {
  local job_id="${1:-}"
  test -n "$job_id" || { echo "Usage: bash run_ga.sh status JOB_ID" >&2; return 2; }
  if command -v squeue >/dev/null 2>&1; then squeue -j "$job_id" || true; fi
  if command -v sacct >/dev/null 2>&1; then sacct -j "$job_id" --format=JobID,JobName,Partition,State,Elapsed,ExitCode,AllocCPUS,MaxRSS || true; fi
  echo "Expected Slurm logs: $BUNDLE_DIR/slurm-$job_id.out and $BUNDLE_DIR/slurm-$job_id.err"
}

case "${1:-}" in
  setup) setup_environment ;;
  preflight) preflight ;;
  run) run_design ;;
  submit) submit_job ;;
  status) shift; show_status "$@" ;;
  __slurm_run) run_design ;;
  *) echo "Usage: bash run_ga.sh {setup|preflight|run|submit|status JOB_ID}" >&2; exit 2 ;;
esac
'''


def _readme(
    resources: ExternalGAResources,
    *,
    repository_commit: str,
    validation_mode: str,
    credential_path: str,
) -> str:
    idt = (
        f"Live IDT scoring is enabled. Create `{credential_path}` on a filesystem visible "
        "to the compute node, keep it outside the checkout, and run `chmod 600` on it."
        if validation_mode == "api"
        else "This request does not use the live IDT API; no credential argument is passed."
    )
    return f"""# External HURDLER GA run

This bundle freezes the complete Colab request and repository commit
`{repository_commit}`. It requests {resources.worker_cpus} worker CPUs,
{resources.memory_gb} GB total memory and {resources.walltime} walltime on the
Slurm `{resources.partition}` partition. GA fitness is parallel; random choices
and all IDT calls remain ordered in the main process.

## New environment

Install Conda/Miniforge, unpack this ZIP on a compute-node-visible filesystem,
then run:

```bash
bash run_ga.sh setup
bash run_ga.sh preflight
bash run_ga.sh run
```

`setup` clones the repository, checks out the frozen commit, creates or updates
the `{resources.conda_environment}` environment from `envs/hurdler.yml`, and
performs an editable install. It never submits a job.

## Existing environment or checkout

```bash
conda activate {resources.conda_environment}
export HURDLER_REPO_DIR=/shared/path/to/clone_repeat_protein
bash run_ga.sh preflight
bash run_ga.sh run
```

The checkout must be exactly at the frozen commit. For Slurm:

```bash
bash run_ga.sh submit
bash run_ga.sh status JOB_ID
```

Repository, credential and output paths must be visible from compute nodes.
An offline cluster may receive a separately copied frozen checkout by setting
`HURDLER_REPO_DIR`, but live IDT mode still needs outbound access to the IDT
OAuth and complexity endpoints.

## IDT credential file

{idt}

Create the external file from the blank template:

```bash
mkdir -p ~/.config/hurdler
cp idt.env.example ~/.config/hurdler/idt.env
chmod 600 ~/.config/hurdler/idt.env
```

Use either `IDT_ACCESS_TOKEN` or all four password-grant fields, never both.
The bundle contains only the external path and authentication method. It never
contains a token, password, client secret, uploaded credential contents or an
environment export of a secret. IDT is used only to score completed DNA.

## Recovery and outputs

Progress is appended to JSONL, the best accepted secondary is checkpointed at
least every 180 seconds (and immediately when improved), and the final result
is written to a UTC-stamped ZIP. A non-accepted design exits nonzero. The
criterion for API acceptance is `{IDT_SCORE_POLICY}`.
"""


def create_external_ga_bundle(
    request: DesignRequestV2,
    destination_dir: str | Path,
    *,
    repository_commit: str,
    environment_file: str | Path,
    resources: ExternalGAResources | None = None,
    repository_url: str = DEFAULT_REPOSITORY_URL,
    idt_credential_path: str = "~/.config/hurdler/idt.env",
    auth_method: str = "auto",
    timestamp: datetime | None = None,
) -> Path:
    """Write one portable ZIP without reading or copying any credential file."""
    resources = resources or ExternalGAResources()
    commit = str(repository_commit).strip().lower()
    if not _SAFE_COMMIT.fullmatch(commit):
        raise ValueError("repository_commit must be a full 40-character Git SHA")
    if auth_method not in {"auto", "password", "access_token"}:
        raise ValueError("auth_method must be auto, password, or access_token")
    if request.validation_mode == "api" and not str(idt_credential_path).strip():
        raise ValueError("Live IDT mode requires an external credential path")
    environment_path = Path(environment_file)
    if not environment_path.is_file():
        raise FileNotFoundError(environment_path)
    request = replace(request, ga_workers=int(resources.worker_cpus))
    request_json = _canonical_json(asdict(request))
    request_bytes = request_json.encode()
    environment_bytes = environment_path.read_bytes()
    request_sha = _sha_bytes(request_bytes)
    environment_sha = _sha_bytes(environment_bytes)
    script = _run_script(
        request_json=request_json,
        request_sha256=request_sha,
        environment_sha256=environment_sha,
        repository_url=repository_url,
        repository_commit=commit,
        resources=resources,
        credential_path=(str(idt_credential_path) if request.validation_mode == "api" else ""),
        auth_method=auth_method,
        validation_mode=request.validation_mode,
    )
    example = (
        "# Use either the token line OR all four password-grant lines.\n"
        "IDT_ACCESS_TOKEN=\n"
        "IDT_CLIENT_ID=\n"
        "IDT_CLIENT_SECRET=\n"
        "IDT_USERNAME=\n"
        "IDT_PASSWORD=\n"
    )
    readme = _readme(
        resources,
        repository_commit=commit,
        validation_mode=request.validation_mode,
        credential_path=str(idt_credential_path),
    )
    files: dict[str, bytes] = {
        "run_ga.sh": script.encode(),
        "request.json": request_bytes,
        "envs/hurdler.yml": environment_bytes,
        "idt.env.example": example.encode(),
        "README_external_ga.md": readme.encode(),
    }
    manifest = {
        "schema_version": "hurdler-external-ga-bundle-v1",
        "created_at_utc": (timestamp or datetime.now(timezone.utc)).astimezone(timezone.utc).isoformat(),
        "repository_url": repository_url,
        "repository_commit": commit,
        "request_sha256": request_sha,
        "environment_sha256": environment_sha,
        "validation_mode": request.validation_mode,
        "idt_credential_path": str(idt_credential_path) if request.validation_mode == "api" else None,
        "idt_auth_method": auth_method if request.validation_mode == "api" else None,
        "credentials_embedded": False,
        "resources": asdict(resources),
        "files": {name: {"sha256": _sha_bytes(data), "size_bytes": len(data)} for name, data in files.items()},
    }
    files["bundle_manifest.json"] = _canonical_json(manifest).encode()
    moment = timestamp or datetime.now(timezone.utc)
    stamp = moment.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    sequence_id = re.sub(r"[^A-Za-z0-9._-]+", "_", request.query.sequence_id).strip("._-") or "sequence"
    destination = Path(destination_dir)
    destination.mkdir(parents=True, exist_ok=True)
    bundle = destination / f"hurdler_{sequence_id}_{stamp}_external_ga.zip"
    with tempfile.NamedTemporaryFile(dir=destination, suffix=".zip", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name, data in files.items():
                info = zipfile.ZipInfo(name)
                info.date_time = (1980, 1, 1, 0, 0, 0)
                info.external_attr = (0o755 if name == "run_ga.sh" else 0o644) << 16
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, data)
        temporary.replace(bundle)
    finally:
        temporary.unlink(missing_ok=True)
    return bundle
