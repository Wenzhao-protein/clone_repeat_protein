from __future__ import annotations

import hashlib
import json
import subprocess
import zipfile
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import pytest

import hurdler.cli as cli
from hurdler.external_ga import ExternalGAResources, create_external_ga_bundle
from hurdler.progress import DesignProgressEvent
from hurdler.vector_design import (
    DESIGN_SCHEMA_VERSION_V2,
    CompatibilityQuery,
    DesignRequestV2,
    DesignResultV2,
    DesignSelection,
)


COMMIT = "2e1b40cdf2647df91ca13f7774769b0c92c7d49c"


def _request(*, validation_mode: str = "api") -> DesignRequestV2:
    query = CompatibilityQuery(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        input_mode="split",
        sequence_id="external_demo",
        n_cap="M",
        repeat_module="ACDEFGHIK",
        c_cap="G",
        repeat_copies=8,
        max_restoration_length_bp=100,
    )
    return DesignRequestV2(
        schema_version=DESIGN_SCHEMA_VERSION_V2,
        query=query,
        selection=DesignSelection("candidate", "pET-28a(+)", "scheme", "BsaI"),
        validation_mode=validation_mode,
        generation_schedule=(10, 20, 40, 60, 80, 100),
        score_weights={"repeated_re_site_excess": 12345.0},
        minimum_secondary_copies=3,
        maximum_secondary_copies=9,
    )


def test_external_bundle_is_reproducible_complete_and_secret_free(tmp_path: Path):
    environment = tmp_path / "hurdler.yml"
    environment.write_text("name: hurdler\ndependencies:\n  - python=3.11\n")
    resources = ExternalGAResources(
        worker_cpus=16,
        memory_gb=32,
        walltime="24:00:00",
        partition="cpu",
        account="lab_account",
        qos="normal",
        constraint="avx2",
    )
    bundle = create_external_ga_bundle(
        _request(),
        tmp_path,
        repository_commit=COMMIT,
        environment_file=environment,
        resources=resources,
        idt_credential_path="~/.config/hurdler/idt.env",
        auth_method="auto",
        timestamp=datetime(2026, 8, 12, 12, 0, tzinfo=timezone.utc),
    )
    assert bundle.name == "hurdler_external_demo_20260812T120000Z_external_ga.zip"
    extracted = tmp_path / "extracted"
    with zipfile.ZipFile(bundle) as archive:
        assert set(archive.namelist()) == {
            "run_ga.sh",
            "request.json",
            "envs/hurdler.yml",
            "idt.env.example",
            "README_external_ga.md",
            "bundle_manifest.json",
        }
        archive.extractall(extracted)
    request_bytes = (extracted / "request.json").read_bytes()
    request = json.loads(request_bytes)
    manifest = json.loads((extracted / "bundle_manifest.json").read_text())
    script = (extracted / "run_ga.sh").read_text()
    assert request["ga_workers"] == 16
    assert request["query"]["max_restoration_length_bp"] == 100
    assert request["maximum_secondary_copies"] == 9
    assert request["score_weights"]["repeated_re_site_excess"] == 12345.0
    assert hashlib.sha256(request_bytes).hexdigest() == manifest["request_sha256"]
    assert manifest["credentials_embedded"] is False
    assert "#SBATCH --cpus-per-task=16" in script
    assert "#SBATCH --mem=32G" in script
    assert "#SBATCH --account=lab_account" in script
    assert "OMP_NUM_THREADS=1" in script
    assert "hurdler idt-preflight" in script
    assert "--idt-credential-file" in script
    assert "IDT_CLIENT_SECRET=" not in script
    assert "IDT_PASSWORD=" not in script
    assert "IDT_ACCESS_TOKEN=" not in script
    subprocess.run(["bash", "-n", extracted / "run_ga.sh"], check=True)


def test_batch_bundle_omits_all_idt_cli_arguments(tmp_path: Path):
    environment = tmp_path / "hurdler.yml"
    environment.write_text("name: hurdler\n")
    bundle = create_external_ga_bundle(
        _request(validation_mode="batch"),
        tmp_path,
        repository_commit=COMMIT,
        environment_file=environment,
        timestamp=datetime(2026, 8, 12, tzinfo=timezone.utc),
    )
    with zipfile.ZipFile(bundle) as archive:
        script = archive.read("run_ga.sh").decode()
        manifest = json.loads(archive.read("bundle_manifest.json"))
    assert "hurdler idt-preflight" not in script
    assert "--idt-credential-file" not in script
    assert manifest["idt_credential_path"] is None


@pytest.mark.parametrize(
    "values",
    [
        {"worker_cpus": 0},
        {"memory_gb": 0},
        {"walltime": "tomorrow"},
        {"partition": "cpu\n#SBATCH --mail-user=x"},
        {"account": "lab account"},
    ],
)
def test_external_resource_values_reject_unsafe_scheduler_input(values):
    with pytest.raises(ValueError):
        ExternalGAResources(**values)


def test_live_bundle_requires_only_an_external_path_not_credential_contents(tmp_path: Path):
    environment = tmp_path / "hurdler.yml"
    environment.write_text("name: hurdler\n")
    with pytest.raises(ValueError, match="external credential path"):
        create_external_ga_bundle(
            _request(),
            tmp_path,
            repository_commit=COMMIT,
            environment_file=environment,
            idt_credential_path="",
        )


def test_ga_workers_schema_is_strict_and_defaults_to_serial():
    request = _request(validation_mode="batch")
    assert request.ga_workers == 1
    assert DesignRequestV2.from_dict(asdict(request)).ga_workers == 1
    with pytest.raises(ValueError, match="ga_workers"):
        DesignRequestV2(**{**asdict(request), "ga_workers": True})


def test_generated_runner_setup_run_submit_and_status_with_fake_tools(tmp_path: Path):
    environment = tmp_path / "hurdler.yml"
    environment.write_text("name: hurdler\n")
    bundle = create_external_ga_bundle(
        _request(validation_mode="batch"),
        tmp_path,
        repository_commit=COMMIT,
        environment_file=environment,
        resources=ExternalGAResources(worker_cpus=1, memory_gb=1),
        timestamp=datetime(2026, 8, 12, tzinfo=timezone.utc),
    )
    root = tmp_path / "runner"
    with zipfile.ZipFile(bundle) as archive:
        archive.extractall(root)
    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir()
    log = tmp_path / "calls.log"

    def executable(name: str, text: str) -> None:
        path = fake_bin / name
        path.write_text(text)
        path.chmod(0o755)

    executable(
        "git",
        f"""#!/usr/bin/env bash
echo "git $*" >> {log}
if test "${{1:-}}" = clone; then mkdir -p "${{@: -1}}/.git"; fi
case " $* " in *" rev-parse HEAD "*) echo {COMMIT};; esac
""",
    )
    executable(
        "conda",
        f"""#!/usr/bin/env bash
echo "conda $*" >> {log}
if test "${{1:-}} ${{2:-}}" = "shell.bash hook"; then echo ':'; fi
if test "${{1:-}} ${{2:-}}" = "env list"; then echo '{{"envs":[]}}'; fi
""",
    )
    executable(
        "python",
        f"""#!/usr/bin/env bash
echo "python $*" >> {log}
if test "${{1:-}} ${{2:-}}" = "-m pip"; then exit 0; fi
exec /usr/bin/python3 "$@"
""",
    )
    executable(
        "hurdler",
        f"""#!/usr/bin/env bash
echo "hurdler $*" >> {log}
while test "$#" -gt 0; do
  if test "$1" = --output-dir; then shift; mkdir -p "$1"; fi
  shift || true
done
""",
    )
    executable("sbatch", f"#!/usr/bin/env bash\necho \"sbatch $*\" >> {log}\necho 'Submitted batch job 123'\n")
    executable("squeue", f"#!/usr/bin/env bash\necho \"squeue $*\" >> {log}\necho RUNNING\n")
    executable("sacct", f"#!/usr/bin/env bash\necho \"sacct $*\" >> {log}\necho COMPLETED\n")
    runner_env = {"PATH": f"{fake_bin}:{Path('/usr/bin')}"}
    script = root / "run_ga.sh"
    for arguments in (("setup",), ("preflight",), ("run",), ("submit",), ("status", "123")):
        completed = subprocess.run(
            ["bash", script, *arguments],
            cwd=root,
            env=runner_env,
            text=True,
            capture_output=True,
        )
        assert completed.returncode == 0, completed.stderr
    calls = log.read_text()
    assert "git clone" in calls
    assert "conda env create" in calls
    assert "python -m pip install -e" in calls
    assert "hurdler design-construct" in calls
    assert "--progress-jsonl" in calls and "--checkpoint-zip" in calls
    assert "sbatch" in calls and "squeue -j 123" in calls and "sacct -j 123" in calls


def test_design_construct_cli_writes_progress_checkpoint_archive_and_failure_code(
    monkeypatch, tmp_path: Path
):
    request = _request(validation_mode="batch")
    request_path = tmp_path / "request.json"
    request_path.write_text(json.dumps(asdict(request)))

    status = {"value": "optimized_unvalidated_batch"}

    def design(request_value, **kwargs):
        kwargs["progress_callback"](
            DesignProgressEvent(stage="ga", status="running", copies=3, generation=1)
        )
        kwargs["checkpoint_callback"](
            {
                "event": "accepted_secondary",
                "validation_mode": "batch",
                "sequence_id": request_value.query.sequence_id,
                "repeat_copies": 3,
                "core_sequence": "ATG" * 3,
                "purchase_sequence": "ATG" * 3,
                "idt_complexity_score": None,
            }
        )
        return DesignResultV2(
            schema_version=DESIGN_SCHEMA_VERSION_V2,
            status=status["value"],
            message="mock",
            request=asdict(request_value),
        )

    def write_outputs(result, destination):
        destination.mkdir(parents=True, exist_ok=True)
        path = destination / "design_summary.json"
        path.write_text(json.dumps(result.to_dict()))
        return {"design_summary_json": str(path)}

    monkeypatch.setattr(cli, "design_construct_v2", design)
    monkeypatch.setattr(cli, "write_design_outputs_v2", write_outputs)
    output = tmp_path / "output"
    progress = tmp_path / "progress.jsonl"
    checkpoint = tmp_path / "checkpoint.zip"
    archives = tmp_path / "archives"
    arguments = [
        "design-construct", "--request", str(request_path), "--output-dir", str(output),
        "--progress-jsonl", str(progress), "--checkpoint-zip", str(checkpoint),
        "--checkpoint-interval-seconds", "180", "--final-archive-dir", str(archives),
        "--fail-on-nonaccepted",
    ]
    assert cli.main(arguments) == 0
    assert json.loads(progress.read_text().splitlines()[0])["stage"] == "ga"
    with zipfile.ZipFile(checkpoint) as archive:
        assert archive.namelist() == ["checkpoint.json"]
    assert list(archives.glob("hurdler_external_demo_*_results.zip"))

    status["value"] = "no_accepted_repeat_construct"
    assert cli.main(arguments) == 1
