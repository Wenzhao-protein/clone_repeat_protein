#!/usr/bin/env python3
"""Wait for Stage 2, recover missing shards, and submit finalization.

Completion is defined from readable Parquet row counts, never from Slurm or
taskrunner status alone.  Every recovery cycle waits for NFS directory caches
to settle and then regenerates the missing list from the original task index.
"""

from __future__ import annotations

import json
import re
import subprocess
import time
from pathlib import Path


REPO = Path("/home/wendai/projects/hurdler/clone_repeat_protein")
RUN_DIR = (
    REPO
    / "studies/hurdler_validation/step04_module_optimization/runs/"
    "run103_expanded_stage2_adaptive"
)
TASK_DIR = RUN_DIR / "taskfiles"
TASKRUNNER = Path("/net/software/taskrunner/taskrunner")
PYTHON = Path("/home/wendai/.conda/envs/hurdler/bin/python")
PRIMARY_JOB_ID = 17506750
POLL_SECONDS = 30
NFS_SETTLE_POLLS = 4
MAX_RECOVERY_CYCLES = 8


def _run(command: list[str], *, capture: bool = False) -> subprocess.CompletedProcess:
    result = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=capture,
        cwd=REPO,
    )
    if capture and result.stdout:
        print(result.stdout, end="", flush=True)
    return result


def _job_is_present(job_id: int) -> bool:
    result = subprocess.run(
        ["squeue", "-h", "-j", str(job_id), "-o", "%A"],
        check=True,
        text=True,
        capture_output=True,
    )
    return bool(result.stdout.strip())


def _wait_for_job(job_id: int) -> None:
    print(f"waiting_for_job={job_id}", flush=True)
    while _job_is_present(job_id):
        time.sleep(POLL_SECONDS)
    print(f"job_finished={job_id}", flush=True)


def _settle_nfs() -> None:
    for poll in range(1, NFS_SETTLE_POLLS + 1):
        print(f"nfs_cache_settle={poll}/{NFS_SETTLE_POLLS}", flush=True)
        time.sleep(POLL_SECONDS)


def _taskrunner_status(state: Path) -> str:
    result = _run(
        [str(TASKRUNNER), f"--state-file={state}", "status"],
        capture=True,
    )
    return result.stdout


def _submit_taskfile(
    tasks: Path,
    state: Path,
    dry_run: Path,
    *,
    cpu: int,
    memory: str,
    time_limit: str,
    throttle: int,
) -> int:
    if state.exists():
        raise FileExistsError(f"Refusing to overwrite taskrunner state: {state}")
    _run([str(TASKRUNNER), f"--state-file={state}", "add", str(tasks)])
    options = [
        "submit",
        "--partition",
        "cpu",
        "--cpu",
        str(cpu),
        "--mem",
        memory,
        "--time",
        time_limit,
        "--array-throttle",
        str(throttle),
    ]
    preview = _run(
        [str(TASKRUNNER), f"--state-file={state}", *options, "--dry-run"],
        capture=True,
    )
    dry_run.write_text(preview.stdout)
    _run([str(TASKRUNNER), f"--state-file={state}", *options])
    status = _taskrunner_status(state)
    match = re.search(r"Slurm JobID:\s+(\d+)", status)
    if match is None:
        raise RuntimeError(f"Taskrunner did not report a Slurm job ID for {state}")
    return int(match.group(1))


def _assert_taskrunner_passed(state: Path) -> None:
    status = _taskrunner_status(state)
    failed = re.search(r"\n\s*Failed:\s+(\d+)", status)
    if failed and int(failed.group(1)):
        raise RuntimeError(f"Taskrunner state has failed tasks: {state}")


def _create_missing(cycle: int) -> tuple[Path, Path, int]:
    tasks = TASK_DIR / f"recovery_authoritative_cycle{cycle}_tasks.txt"
    index = TASK_DIR / f"recovery_authoritative_cycle{cycle}_index.csv"
    result = _run(
        [
            str(PYTHON),
            str(REPO / "scripts/create_missing_stage2_tasks.py"),
            "--task-index",
            str(TASK_DIR / "task_index.csv"),
            "--tasks",
            str(TASK_DIR / "tasks.txt"),
            "--output-taskfile",
            str(tasks),
            "--output-index",
            str(index),
        ],
        capture=True,
    )
    payload = json.loads(result.stdout)
    return tasks, index, int(payload["missing_tasks"])


def main() -> int:
    _wait_for_job(PRIMARY_JOB_ID)
    _settle_nfs()

    for cycle in range(1, MAX_RECOVERY_CYCLES + 1):
        tasks, _index, missing = _create_missing(cycle)
        print(f"recovery_cycle={cycle} missing_tasks={missing}", flush=True)
        if missing == 0:
            break
        state = TASK_DIR / f"recovery_authoritative_cycle{cycle}.state"
        job_id = _submit_taskfile(
            tasks,
            state,
            TASK_DIR / f"recovery_authoritative_cycle{cycle}_dry_run.sh",
            cpu=1,
            memory="4G",
            time_limit="02:00:00",
            throttle=64,
        )
        _wait_for_job(job_id)
        _settle_nfs()
        # A failed/cancelled element is evidence for the next missing-only
        # cycle, not a reason to abandon the recovery loop.
        _taskrunner_status(state)
    else:
        raise RuntimeError("Stage-2 recovery did not converge within eight cycles")

    final_state = TASK_DIR / "finalize_authoritative.state"
    final_job = _submit_taskfile(
        TASK_DIR / "finalize_task.txt",
        final_state,
        TASK_DIR / "finalize_authoritative_dry_run.sh",
        cpu=2,
        memory="32G",
        time_limit="02:00:00",
        throttle=1,
    )
    _wait_for_job(final_job)
    _assert_taskrunner_passed(final_state)
    print(json.dumps({"status": "passed", "finalize_job_id": final_job}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
