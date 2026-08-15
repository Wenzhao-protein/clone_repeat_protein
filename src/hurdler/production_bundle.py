"""Portable Digs/taskrunner production bundle generation.

No command in this module submits work.  Generated scripts require an
explicit ``submit`` action, preserve taskrunner state, and never clean up
files automatically.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import stat
import subprocess
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping

from .io import sha256_file, utc_now, write_json_atomic


PRODUCTION_BUNDLE_SCHEMA = "hurdler-production-bundle-v2"
_SAFE_OPTION = re.compile(r"^[A-Za-z0-9_.:/+@-]*$")


@dataclass(frozen=True)
class WorkflowSpec:
    workflow_id: str
    title: str
    default_shards: int
    cpu_per_task: int
    memory: str
    walltime: str
    array_throttle: int
    requires_idt: bool = False
    gpu_per_task: int = 0


WORKFLOWS: dict[str, WorkflowSpec] = {
    item.workflow_id: item
    for item in (
        WorkflowSpec("success-landscape", "Success landscape 1-60 AA", 1, 16, "32G", "02:00:00", 1),
        WorkflowSpec("repeatsdb-natural", "RepeatsDB natural extraction", 128, 8, "16G", "02:00:00", 16),
        WorkflowSpec("designed-structure", "Designed DSSP/Foldseek", 1, 4, "8G", "02:00:00", 16),
        WorkflowSpec("missing-af3", "Missing-only AlphaFold3", 1, 4, "32G", "08:00:00", 4, gpu_per_task=1),
        WorkflowSpec("module-stage1", "Stage-1 module compatibility", 128, 1, "8G", "02:00:00", 16),
        WorkflowSpec("module-stage2", "Stage-2 adaptive GA/IDT", 1, 1, "8G", "02:00:00", 16, requires_idt=True),
        WorkflowSpec("exact-dna-routes", "Exact-DNA complete routes", 512, 1, "8G", "02:00:00", 16, requires_idt=True),
        WorkflowSpec("exact-dna-purchase", "Exact-DNA purchase audit", 1, 1, "8G", "00:30:00", 1, requires_idt=True),
        WorkflowSpec("reports", "Papermill/report generation", 1, 1, "8G", "02:00:00", 1),
    )
}


@dataclass(frozen=True)
class ClusterProfile:
    repo_root: str
    scratch_root: str
    conda_prefix: str
    taskrunner: str = "/net/software/taskrunner/taskrunner"
    container: str = ""
    foldseek: str = "/net/software/utils/foldseek"
    mkdssp: str = "mkdssp"
    mafft: str = "mafft"
    af3_runner: str = ""
    partition: str = "cpu"
    account: str = ""
    qos: str = ""
    constraint: str = ""
    cpu_per_task: int = 1
    memory: str = "8G"
    walltime: str = "02:00:00"
    array_throttle: int = 16
    gpu: str = ""
    idt_env_path: str = "~/.config/hurdler/idt.env"

    def validate(self) -> None:
        for name in ("repo_root", "scratch_root", "conda_prefix", "taskrunner"):
            value = getattr(self, name)
            if not Path(value).expanduser().is_absolute():
                raise ValueError(f"cluster_profile.{name} must be an absolute path")
        for name in ("container", "af3_runner"):
            value = getattr(self, name)
            if value and not Path(value).expanduser().is_absolute():
                raise ValueError(f"cluster_profile.{name} must be an absolute path when set")
        if self.cpu_per_task < 1 or self.array_throttle < 1:
            raise ValueError("CPU and array throttle must be positive")
        for name in ("partition", "account", "qos", "constraint", "memory", "walltime", "gpu"):
            if not _SAFE_OPTION.fullmatch(getattr(self, name)):
                raise ValueError(f"Unsafe cluster option: {name}")


@dataclass(frozen=True)
class ProductionBundleRequest:
    workflow_id: str
    parameter_version: str
    repo_commit: str
    cluster_profile: ClusterProfile
    inputs: tuple[dict[str, Any], ...] = ()
    shard_count: int | None = None
    scientific_parameters: dict[str, Any] = field(default_factory=dict)
    random_seed: int = 42
    idt_mode: str = "batch"
    output_dir: str = ""
    scratch_dir: str = ""
    finalization_strategy: str = "validate_then_finalize"
    resume_strategy: str = "missing_only"

    def validate(self) -> WorkflowSpec:
        if self.workflow_id not in WORKFLOWS:
            raise ValueError(f"Unknown workflow: {self.workflow_id}")
        if self.idt_mode not in {"external_path", "batch"}:
            raise ValueError("idt_mode must be external_path or batch")
        spec = WORKFLOWS[self.workflow_id]
        if spec.requires_idt and self.idt_mode == "external_path" and not self.cluster_profile.idt_env_path:
            raise ValueError("Live IDT workflows require an external credential path")
        if self.shard_count is not None and self.shard_count < 1:
            raise ValueError("shard_count must be positive")
        self.cluster_profile.validate()
        for field_name in ("output_dir", "scratch_dir"):
            if not Path(getattr(self, field_name)).expanduser().is_absolute():
                raise ValueError(f"{field_name} must be an absolute path")
        for item in self.inputs:
            if "path" not in item or "sha256" not in item:
                raise ValueError("Every input requires path and sha256")
            if not Path(str(item["path"])).expanduser().is_absolute():
                raise ValueError("Production input paths must be absolute")
            if item.get("kind", "file") not in {"file", "directory"}:
                raise ValueError("Production input kind must be file or directory")
            if len(str(item["sha256"])) != 64:
                raise ValueError("Production input SHA256 must contain 64 hexadecimal characters")
        return spec

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ProductionBundleRequest":
        data = dict(payload)
        data["cluster_profile"] = ClusterProfile(**data["cluster_profile"])
        data["inputs"] = tuple(data.get("inputs", ()))
        return cls(**data)

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        if self.idt_mode == "batch":
            payload["cluster_profile"]["idt_env_path"] = ""
        return payload


def _q(value: object) -> str:
    return shlex.quote(str(value))


def directory_manifest_sha256(path: str | Path) -> str:
    """Hash a directory as sorted ``sha256sum`` lines without archiving it."""
    root = Path(path).expanduser().absolute()
    if not root.is_dir():
        raise NotADirectoryError(root)
    digest = hashlib.sha256()
    for candidate in sorted((item for item in root.rglob("*") if item.is_file()), key=lambda item: item.relative_to(root).as_posix()):
        relative = candidate.relative_to(root).as_posix()
        digest.update(f"{sha256_file(candidate)}  ./{relative}\n".encode())
    return digest.hexdigest()


def _hurdler(request: ProductionBundleRequest) -> str:
    return str(Path(request.cluster_profile.conda_prefix) / "bin" / "hurdler")


def _python(request: ProductionBundleRequest) -> str:
    return str(Path(request.cluster_profile.conda_prefix) / "bin" / "python")


def _input(request: ProductionBundleRequest, index: int = 0) -> str:
    try:
        return str(request.inputs[index]["path"])
    except IndexError as exc:
        raise ValueError(f"Workflow {request.workflow_id} requires input {index + 1}") from exc


def _live_idt_args(request: ProductionBundleRequest) -> list[str]:
    if request.idt_mode != "external_path":
        return []
    return ["--use-idt", "--credential-mode", "path", "--credential-path", request.cluster_profile.idt_env_path]


def _table_row_count(path: str | Path, *, compatible_only: bool = False) -> int:
    """Read only metadata/the compatibility flag when dynamic sharding is needed."""
    source = Path(path)
    if source.suffix.lower() == ".parquet":
        if compatible_only:
            import pandas as pd

            frame = pd.read_parquet(source, columns=["hurdler_compatible"])
            return int(frame.hurdler_compatible.astype(bool).sum())
        import pyarrow.parquet as pq

        return int(pq.ParquetFile(source).metadata.num_rows)
    if source.suffix.lower() in {".csv", ".tsv"}:
        import pandas as pd

        if compatible_only:
            frame = pd.read_csv(
                source, sep="\t" if source.suffix.lower() == ".tsv" else ",",
                usecols=["hurdler_compatible"],
            )
            return int(frame.hurdler_compatible.astype(bool).sum())
        with source.open("rb") as handle:
            return max(0, sum(1 for _ in handle) - 1)
    raise ValueError(f"Dynamic production sharding requires CSV/TSV/Parquet: {source}")


def _resolved_shard_count(request: ProductionBundleRequest, spec: WorkflowSpec) -> int:
    if request.shard_count is not None:
        return request.shard_count
    if request.workflow_id in {"designed-structure", "missing-af3"}:
        return max(1, _table_row_count(_input(request)))
    if request.workflow_id == "module-stage2":
        modules_per_task = int(request.scientific_parameters.get("modules_per_task", 4))
        if modules_per_task < 1:
            raise ValueError("modules_per_task must be positive")
        compatible = _table_row_count(_input(request), compatible_only=True)
        return max(1, math.ceil(compatible / modules_per_task))
    return spec.default_shards


def _command_rows(request: ProductionBundleRequest) -> tuple[list[str], list[dict[str, Any]], str]:
    spec = request.validate()
    shards = _resolved_shard_count(request, spec)
    hurdler = _hurdler(request)
    python = _python(request)
    repo = request.cluster_profile.repo_root
    scratch = request.scratch_dir
    output = request.output_dir
    params = request.scientific_parameters
    commands: list[list[str]] = []
    outputs: list[str] = []
    finalize: list[str]

    if request.workflow_id == "success-landscape":
        commands.append([
            python, str(Path(repo) / "scripts" / "run_success_landscape_single_files.py"),
            "--repo-dir", repo,
            "--short-output", str(Path(scratch) / "short_motifs_1_5.parquet"),
            "--random-output", str(Path(scratch) / "random_modules_6_60.parquet"),
            "--figure-dir", str(Path(output) / "figures"),
            "--workers", str(request.cluster_profile.cpu_per_task),
            "--scan-copies", str(params.get("scan_copies", 2)),
            "--short-max-length", "5", "--random-min-length", "6",
            "--random-max-length", "60", "--tests", str(params.get("tests", 1000)),
        ])
        outputs.append(str(Path(scratch) / "random_modules_6_60.parquet"))
        finalize = [python, str(Path(repo) / "scripts" / "run_success_landscape_single_files.py"),
                    "--repo-dir", repo, "--short-output", str(Path(scratch) / "short_motifs_1_5.parquet"),
                    "--random-output", str(Path(scratch) / "random_modules_6_60.parquet"),
                    "--figure-dir", str(Path(output) / "figures"), "--workers", str(request.cluster_profile.cpu_per_task),
                    "--scan-copies", str(params.get("scan_copies", 2)), "--validate-only"]
    elif request.workflow_id == "repeatsdb-natural":
        inventory = _input(request)
        for shard in range(shards):
            shard_root = Path(scratch) / f"shard_{shard:05d}"
            commands.append([
                hurdler, "curate-modules", "--all-repeatsdb", "--one-per-protein",
                "--annotation-inventory", inventory, "--natural-output", str(shard_root / "natural.parquet"),
                "--natural-mappings-output", str(shard_root / "natural_source_mappings.parquet"),
                "--natural-exclusions-output", str(shard_root / "natural_exclusions.csv"),
                "--natural-cache-dir", str(Path(scratch) / "fasta_cache"), "--natural-workers", str(request.cluster_profile.cpu_per_task),
                "--natural-shard-index", str(shard), "--natural-shard-count", str(shards),
            ])
            outputs.append(str(shard_root / "natural.parquet"))
        finalize = [python, str(Path(repo) / "scripts" / "finalize_v2_workflow.py"),
                    "--workflow", request.workflow_id, "--scratch-dir", scratch, "--output-dir", output,
                    "--input", inventory]
    elif request.workflow_id == "designed-structure":
        inventory = _input(request)
        for shard in range(shards):
            shard_root = Path(scratch) / f"shard_{shard:05d}"
            commands.append([
                hurdler, "infer-designed-boundaries", "--input", inventory,
                "--output", str(shard_root / "designed.parquet"),
                "--candidates-output", str(shard_root / "candidates.parquet"),
                "--units-output", str(shard_root / "units.parquet"),
                "--positions-output", str(shard_root / "positions.parquet"),
                "--exclusions-output", str(shard_root / "exclusions.parquet"),
                "--dssp-engine", "biotite",
                "--mkdssp", request.cluster_profile.mkdssp, "--foldseek", request.cluster_profile.foldseek,
                "--mafft", request.cluster_profile.mafft, "--shard-index", str(shard), "--shard-count", str(shards),
            ])
            outputs.append(str(shard_root / "designed.parquet"))
        finalize = [python, str(Path(repo) / "scripts" / "finalize_v2_workflow.py"),
                    "--workflow", request.workflow_id, "--scratch-dir", scratch, "--output-dir", output,
                    "--input", inventory]
    elif request.workflow_id == "missing-af3":
        if not request.cluster_profile.af3_runner:
            raise ValueError("missing-af3 requires cluster_profile.af3_runner")
        if not request.cluster_profile.gpu:
            raise ValueError("missing-af3 requires cluster_profile.gpu (the Digs GPU type)")
        inventory = _input(request)
        for shard in range(shards):
            destination = str(Path(scratch) / f"shard_{shard:05d}")
            commands.append([request.cluster_profile.af3_runner, "--inventory", inventory,
                             "--shard-index", str(shard), "--shard-count", str(shards),
                             "--seed", str(request.random_seed), "--diffusion-samples", "1", "--output-dir", destination])
            outputs.append(destination)
        finalize = [python, str(Path(repo) / "scripts" / "finalize_v2_workflow.py"),
                    "--workflow", request.workflow_id, "--scratch-dir", scratch, "--output-dir", output,
                    "--input", inventory]
    elif request.workflow_id == "module-stage1":
        catalog = _input(request)
        index_dir = str(params.get("index_dir", Path(repo) / "data/artifacts/legacy-optimized-v1"))
        for shard in range(shards):
            shard_root = Path(scratch) / f"shard_{shard:05d}"
            suffix = f"shard-{shard:05d}-of-{shards:05d}"
            commands.append([hurdler, "module-compatibility", "--catalog", catalog,
                             "--index-dir", index_dir, "--output-dir", str(shard_root),
                             "--shard-index", str(shard), "--shard-count", str(shards)])
            outputs.append(str(shard_root / f"module_compatibility_{suffix}.parquet"))
        finalize = [python, str(Path(repo) / "scripts" / "finalize_v2_workflow.py"),
                    "--workflow", request.workflow_id, "--scratch-dir", scratch, "--output-dir", output,
                    "--input", catalog]
    elif request.workflow_id == "module-stage2":
        module_inputs = _input(request)
        for shard in range(shards):
            shard_root = Path(scratch) / f"shard_{shard:05d}"
            idt_arguments = (
                ["--credential-path", request.cluster_profile.idt_env_path]
                if request.idt_mode == "external_path"
                else ["--idt-batch"]
            )
            command = [hurdler, "adaptive-copy-search", "--compatibility", module_inputs,
                       "--output-dir", str(shard_root), "--shard-index", str(shard),
                       "--shard-count", str(shards), "--seed", str(request.random_seed),
                       "--short-generations", "10", "--generation-schedule", "10", "20", "40", "60", "80", "100",
                       "--idt-policy", "idt-rule-score-sum-lt10-v1", *idt_arguments]
            commands.append(command)
            outputs.append(str(shard_root / "optimized_constructs_ga.parquet"))
        finalize = [python, str(Path(repo) / "scripts" / "finalize_v2_workflow.py"),
                    "--workflow", request.workflow_id, "--scratch-dir", scratch, "--output-dir", output,
                    "--input", module_inputs]
    elif request.workflow_id == "exact-dna-routes":
        catalog = _input(request)
        for shard in range(shards):
            shard_root = Path(scratch) / f"shard_{shard:05d}"
            commands.append([hurdler, "dna-assembly", "plan-complete", "--catalog", catalog,
                             "--reference-dir", str(Path(repo) / "data/reference_output"),
                             "--artifact-dir", str(Path(repo) / "data/artifacts"), "--output-dir", str(shard_root),
                             "--shard-index", str(shard), "--shard-count", str(shards), *_live_idt_args(request)])
            outputs.append(str(shard_root / "complete_route_manifest.json"))
        finalize = [python, str(Path(repo) / "scripts" / "finalize_v2_workflow.py"),
                    "--workflow", request.workflow_id, "--scratch-dir", scratch, "--output-dir", output,
                    "--input", catalog,
                    "--expected-elements", str(params.get("expected_elements", 29042)),
                    "--expected-targets", str(params.get("expected_targets", 145210))]
    elif request.workflow_id == "exact-dna-purchase":
        raw = _input(request)
        purchase_output = str(Path(scratch) / "production")
        commands.append([hurdler, "dna-assembly", "audit-purchases", "--raw-root", raw,
                         "--output-dir", purchase_output,
                         "--expected-shards", str(params.get("expected_shards", 512)),
                         "--expected-routes", str(params.get("expected_routes", 15535)),
                         "--expected-elements", str(params.get("expected_elements", 3129)),
                         *_live_idt_args(request)])
        outputs.append(str(Path(purchase_output) / "purchase_orderability_summary.parquet"))
        finalize = [python, str(Path(repo) / "scripts" / "finalize_v2_workflow.py"),
                    "--workflow", request.workflow_id, "--scratch-dir", scratch, "--output-dir", output,
                    "--expected-elements", str(params.get("expected_elements", 3129))]
    else:  # reports
        notebook = _input(request)
        executed = str(Path(output) / (Path(notebook).stem + "_executed.ipynb"))
        commands.append([str(Path(request.cluster_profile.conda_prefix) / "bin" / "papermill"), notebook, executed])
        outputs.append(executed)
        finalize = [str(Path(request.cluster_profile.conda_prefix) / "bin" / "jupyter"), "nbconvert", "--to", "html", executed,
                    "--output-dir", output]

    environment_prefix = [
        "/usr/bin/env", "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1", "NUMEXPR_NUM_THREADS=1",
    ]
    rendered = [
        " ".join(_q(token) for token in [*environment_prefix, *command])
        for command in commands
    ]
    rows = [
        {"task_id": index + 1, "shard_index": index, "shard_count": len(rendered), "expected_output": outputs[index], "command": command}
        for index, command in enumerate(rendered)
    ]
    return rendered, rows, " ".join(_q(token) for token in finalize)


def _script_header() -> str:
    return "#!/usr/bin/env bash\nset -euo pipefail\nIFS=$'\\n\\t'\n"


def _submission_script(request: ProductionBundleRequest, spec: WorkflowSpec) -> str:
    profile = request.cluster_profile
    conda_hurdler = Path(profile.conda_prefix) / "bin" / "hurdler"
    options = ["--partition", profile.partition, "--cpu", str(profile.cpu_per_task),
               "--mem", profile.memory, "--time", profile.walltime,
               "--array-throttle", str(profile.array_throttle)]
    if profile.gpu:
        options += ["--gpu", profile.gpu]
    if profile.container:
        options += ["--apptainer", profile.container]
    option_text = " ".join(_q(item) for item in options)
    slurm_environment = ""
    if profile.account:
        slurm_environment += f"export SBATCH_ACCOUNT={_q(profile.account)}\n"
    if profile.qos:
        slurm_environment += f"export SBATCH_QOS={_q(profile.qos)}\n"
    if profile.constraint:
        slurm_environment += f"export SBATCH_CONSTRAINT={_q(profile.constraint)}\n"
    credential_check = ""
    if request.idt_mode == "external_path":
        credential_check = f'''\ncredential={_q(profile.idt_env_path)}
credential="${{credential/#\\~/$HOME}}"
[[ -f "$credential" ]] || {{ echo "Missing external IDT credential file" >&2; exit 2; }}
[[ ! "$credential" -ef "$repo_root/config/idt.env.example" ]] || {{ echo "Credential must be outside the repository" >&2; exit 2; }}
mode=$(stat -c '%a' "$credential")
(( (8#$mode & 8#077) == 0 )) || {{ echo "IDT credential file must be mode 600 or stricter" >&2; exit 2; }}
[[ $(stat -c '%u' "$credential") -eq $(id -u) ]] || {{ echo "IDT credential file must belong to the current user" >&2; exit 2; }}
'''
    input_checks = ""
    for item in request.inputs:
        input_checks += f"\ninput_path={_q(item['path'])}\n"
        if item.get("kind", "file") == "directory":
            input_checks += (
                '[[ -d "$input_path" ]] || { echo "Missing production input directory: $input_path" >&2; exit 2; }\n'
                'observed_input_sha=$(cd "$input_path" && find . -type f -print0 | LC_ALL=C sort -z | '
                "xargs -0 -r sha256sum | sha256sum | cut -d' ' -f1)\n"
                f"[[ \"$observed_input_sha\" == {_q(item['sha256'])} ]] || "
                '{ echo "Production input directory checksum mismatch: $input_path" >&2; exit 2; }\n'
            )
        else:
            input_checks += (
                '[[ -f "$input_path" ]] || { echo "Missing production input: $input_path" >&2; exit 2; }\n'
                f"[[ $(sha256sum \"$input_path\" | cut -d' ' -f1) == {_q(item['sha256'])} ]] || "
                '{ echo "Production input checksum mismatch: $input_path" >&2; exit 2; }\n'
            )
    return _script_header() + f'''bundle_dir=$(cd -- "$(dirname -- "$0")" && pwd)
repo_root={_q(profile.repo_root)}
scratch_dir={_q(request.scratch_dir)}
output_dir={_q(request.output_dir)}
state="$bundle_dir/production.state"
taskrunner={_q(profile.taskrunner)}
action="${{1:-help}}"
{slurm_environment}

verify() {{
  cd "$bundle_dir"
  sha256sum --check SHA256SUMS
  [[ $(git -C "$repo_root" rev-parse HEAD) == {_q(request.repo_commit)} ]] || {{ echo "Repository commit mismatch" >&2; exit 2; }}
  [[ -x "$taskrunner" ]] || {{ echo "Taskrunner is not executable: $taskrunner" >&2; exit 2; }}
  [[ -x {_q(conda_hurdler)} ]] || {{ echo "HURDLER is not installed in the requested conda prefix" >&2; exit 2; }}
  mkdir -p "$scratch_dir" "$output_dir"
  [[ -d "$scratch_dir" && -w "$scratch_dir" ]] || {{ echo "Scratch directory is not writable" >&2; exit 2; }}
  [[ -d "$output_dir" && -w "$output_dir" ]] || {{ echo "Output directory is not writable" >&2; exit 2; }}
  while IFS= read -r command; do
    [[ "$command" == /* ]] || {{ echo "Task does not start with an absolute executable: $command" >&2; exit 2; }}
    executable=${{command%% *}}
    [[ -x "$executable" ]] || {{ echo "Task executable is not available: $executable" >&2; exit 2; }}
  done < "$bundle_dir/tasks.txt"
{input_checks}{credential_check}}}

case "$action" in
  preflight)
    verify
    ;;
  add)
    verify
    [[ ! -e "$state" ]] || {{ echo "Refusing to overwrite taskrunner state: $state" >&2; exit 2; }}
    "$taskrunner" --state-file="$state" add "$bundle_dir/tasks.txt"
    "$taskrunner" --state-file="$state" submit {option_text} --dry-run > "$bundle_dir/production_dry_run.sh"
    ;;
  submit)
    verify
    [[ -f "$state" ]] || {{ echo "Run '$0 add' first" >&2; exit 2; }}
    "$taskrunner" --state-file="$state" submit {option_text} | tee "$bundle_dir/submission_output.txt"
    grep -Eo '[0-9]+' "$bundle_dir/submission_output.txt" | tail -1 > "$bundle_dir/job_id.txt" || true
    "$taskrunner" --state-file="$state" status | tee "$bundle_dir/submission_status.txt"
    ;;
  status)
    [[ -f "$state" ]] || {{ echo "No taskrunner state exists" >&2; exit 2; }}
    "$taskrunner" --state-file="$state" status
    ;;
  *)
    echo "Usage: $0 preflight|add|submit|status"
    ;;
esac
'''


def _simple_script(command: str) -> str:
    return _script_header() + 'bundle_dir=$(cd -- "$(dirname -- "$0")" && pwd)\n' + command + "\n"


def _write_checksums(bundle: Path) -> None:
    members = sorted(
        path for path in bundle.rglob("*")
        if path.is_file() and path.name not in {"SHA256SUMS", "production.state"}
    )
    text = "".join(
        f"{sha256_file(path)}  {path.relative_to(bundle)}\n" for path in members
    )
    (bundle / "SHA256SUMS").write_text(text)


def build_production_bundle(request: ProductionBundleRequest, output_dir: str | Path) -> Path:
    spec = request.validate()
    bundle = Path(output_dir).expanduser().absolute()
    if bundle.exists() and any(bundle.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty bundle directory: {bundle}")
    bundle.mkdir(parents=True, exist_ok=True)
    tasks, rows, finalize = _command_rows(request)
    (bundle / "tasks.txt").write_text("\n".join(tasks) + "\n")
    (bundle / "smoke_tasks.txt").write_text(tasks[0] + "\n")
    with (bundle / "task_index.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    request_payload = request.to_dict()
    write_json_atomic(request_payload, bundle / "request.json")
    write_json_atomic(request_payload["cluster_profile"], bundle / "cluster_profile.json")
    shutil_source = Path(request.cluster_profile.repo_root) / "envs" / "hurdler.yml"
    env_dir = bundle / "envs"
    env_dir.mkdir()
    if shutil_source.exists():
        (env_dir / "hurdler.yml").write_bytes(shutil_source.read_bytes())
    else:
        (env_dir / "hurdler.yml").write_text("name: hurdler\ndependencies:\n  - python=3.11\n")
    output_paths = [row["expected_output"] for row in rows]
    write_json_atomic(
        {
            "schema_version": PRODUCTION_BUNDLE_SCHEMA,
            "created_at": utc_now(),
            "workflow": request.workflow_id,
            "repo_commit": request.repo_commit,
            "task_count": len(tasks),
            "status": "bundle_created_not_submitted",
            "expected_outputs": output_paths,
            "cleanup_performed": False,
        },
        bundle / "run_manifest.json",
    )
    (bundle / "submit_digs.sh").write_text(_submission_script(request, spec))
    (bundle / "preflight.sh").write_text(_simple_script('exec "$bundle_dir/submit_digs.sh" preflight'))
    (bundle / "status.sh").write_text(_simple_script('exec "$bundle_dir/submit_digs.sh" status'))
    (bundle / "finalize.sh").write_text(_simple_script(finalize))
    (bundle / "resume_missing.sh").write_text(_simple_script(
        f'{_q(_python(request))} {_q(Path(request.cluster_profile.repo_root) / "scripts" / "create_missing_v2_tasks.py")} '
        ' --task-index "$bundle_dir/task_index.csv" --output "$bundle_dir/missing_tasks.txt"'
    ))
    (bundle / "README.md").write_text(
        f"# HURDLER production bundle: {request.workflow_id}\n\n"
        "This bundle never submits automatically. Run `bash preflight.sh`, then "
        "`bash submit_digs.sh add`, inspect `production_dry_run.sh`, and explicitly run "
        "`bash submit_digs.sh submit`. Use `status.sh`, `resume_missing.sh`, and "
        "`finalize.sh` afterwards. No cleanup command is executed by this bundle.\n"
    )
    for name in ("submit_digs.sh", "preflight.sh", "status.sh", "resume_missing.sh", "finalize.sh"):
        path = bundle / name
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
    _write_checksums(bundle)
    return bundle


def validate_production_bundle(path: str | Path) -> dict[str, Any]:
    bundle = Path(path).expanduser().absolute()
    required = {
        "request.json", "cluster_profile.json", "tasks.txt", "task_index.csv",
        "preflight.sh", "submit_digs.sh", "status.sh", "resume_missing.sh",
        "finalize.sh", "run_manifest.json", "SHA256SUMS", "README.md",
        "envs/hurdler.yml",
    }
    missing = sorted(name for name in required if not (bundle / name).is_file())
    if missing:
        raise FileNotFoundError(f"Bundle is missing: {', '.join(missing)}")
    subprocess.run(["sha256sum", "--check", "SHA256SUMS"], cwd=bundle, check=True, capture_output=True, text=True)
    for script in ("preflight.sh", "submit_digs.sh", "status.sh", "resume_missing.sh", "finalize.sh"):
        subprocess.run(["bash", "-n", str(bundle / script)], check=True)
    request = ProductionBundleRequest.from_dict(json.loads((bundle / "request.json").read_text()))
    request.validate()
    tasks = [line for line in (bundle / "tasks.txt").read_text().splitlines() if line.strip()]
    if not tasks or any(not line.startswith("/") for line in tasks):
        raise ValueError("Every task must start with an absolute executable")
    with (bundle / "task_index.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != len(tasks):
        raise ValueError("Task index and task file counts differ")
    if any(row["command"] != tasks[index] for index, row in enumerate(rows)):
        raise ValueError("Task index commands differ from tasks.txt")
    return {"status": "passed", "workflow": request.workflow_id, "task_count": len(tasks)}
