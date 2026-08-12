"""Credential-safe client for IDT SciTools complexity screening."""

from __future__ import annotations

import json
import hashlib
import getpass
import math
import os
import shlex
import stat
import time
from copy import deepcopy
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

TOKEN_URL = "https://www.idtdna.com/Identityserver/connect/token"
GBLOCK_COMPLEXITY_URL = "https://www.idtdna.com/restapi/v1/Complexities/ScreenGblockSequences"
REQUIRED_CREDENTIAL_ENV = (
    "IDT_CLIENT_ID",
    "IDT_CLIENT_SECRET",
    "IDT_USERNAME",
    "IDT_PASSWORD",
)
ACCESS_TOKEN_ENV = "IDT_ACCESS_TOKEN"
ALLOWED_CREDENTIAL_ENV = REQUIRED_CREDENTIAL_ENV + (ACCESS_TOKEN_ENV,)

IDT_SCORE_POLICY = "idt-rule-score-sum-lt10-v1"
IDT_ORDERABILITY_THRESHOLD = 10.0
IDT_CREDENTIAL_PATH = Path.home() / ".config" / "hurdler" / "idt.env"
_EPHEMERAL_CREDENTIAL_FLAG = "HURDLER_IDT_EPHEMERAL_CREDENTIALS"


def _repository_root() -> Path | None:
    """Return the repository root when called from an installed checkout."""
    current = Path(__file__).absolute()
    for candidate in current.parents:
        if (candidate / "pyproject.toml").is_file() and (candidate / "src" / "hurdler").is_dir():
            return candidate
    return None


def validate_idt_credential_path(
    path: str | Path,
    *,
    repository_root: str | Path | None = None,
) -> Path:
    """Validate an external, owner-only IDT credential file.

    The path itself is safe to report, but its contents and hash are not.  On
    POSIX systems the file must be owned by the current user and have no
    group/other permission bits.
    """
    credential_path = Path(path).expanduser().absolute()
    if not credential_path.is_file():
        raise FileNotFoundError(f"IDT credential file not found: {credential_path}")
    root_value = Path(repository_root).expanduser().absolute() if repository_root else _repository_root()
    if root_value is not None:
        try:
            credential_path.relative_to(root_value)
        except ValueError:
            pass
        else:
            raise ValueError("IDT credential files must be stored outside the repository")
    file_stat = credential_path.stat()
    if os.name == "posix":
        if stat.S_IMODE(file_stat.st_mode) & 0o077:
            raise PermissionError("IDT credential file must have mode 600 or stricter")
        if file_stat.st_uid != os.getuid():
            raise PermissionError("IDT credential file must be owned by the current user")
    return credential_path


def _parse_credential_lines(lines: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[7:].lstrip()
        if "=" not in line:
            continue
        name, raw_value = line.split("=", 1)
        name = name.strip()
        if name not in ALLOWED_CREDENTIAL_ENV:
            continue
        tokens = shlex.split(raw_value, comments=True, posix=True)
        value = tokens[0] if tokens else ""
        if value:
            values[name] = value
    return values


def _validate_credential_values(values: dict[str, str], auth_method: str | None = None) -> str:
    password_complete = all(values.get(name) for name in REQUIRED_CREDENTIAL_ENV)
    password_partial = any(values.get(name) for name in REQUIRED_CREDENTIAL_ENV)
    token_present = bool(values.get(ACCESS_TOKEN_ENV))
    if password_complete and token_present:
        raise RuntimeError("Provide either IDT password-grant fields or IDT_ACCESS_TOKEN, not both")
    if auth_method == "password" and not password_complete:
        raise RuntimeError("Password authentication requires all four IDT OAuth fields")
    if auth_method == "access_token" and not token_present:
        raise RuntimeError("Access-token authentication requires IDT_ACCESS_TOKEN")
    if token_present:
        if password_partial:
            raise RuntimeError("Do not mix IDT_ACCESS_TOKEN with password-grant fields")
        return "access_token"
    if password_complete:
        return "password"
    missing = [name for name in REQUIRED_CREDENTIAL_ENV if not values.get(name)]
    raise RuntimeError(
        "IDT credentials require IDT_ACCESS_TOKEN or all password-grant fields; missing: "
        + ", ".join(missing)
    )


def load_idt_credentials(
    path: str | Path = IDT_CREDENTIAL_PATH,
    *,
    auth_method: str | None = None,
    repository_root: str | Path | None = None,
) -> Path:
    """Load a validated external env file without returning or logging values."""
    credential_path = validate_idt_credential_path(path, repository_root=repository_root)
    values = _parse_credential_lines(credential_path.read_text().splitlines())
    _validate_credential_values(values, auth_method)
    clear_idt_secret_environment()
    os.environ.update(values)
    os.environ.pop(_EPHEMERAL_CREDENTIAL_FLAG, None)
    return credential_path


def prompt_idt_credentials(
    *,
    auth_method: str = "password",
    headless: bool = False,
    prompt: Any = getpass.getpass,
) -> dict[str, object]:
    """Read credentials invisibly for an interactive notebook kernel."""
    if headless:
        raise RuntimeError("Manual IDT credentials are disabled in headless/Papermill execution; use path mode")
    if auth_method not in {"password", "access_token"}:
        raise ValueError("auth_method must be 'password' or 'access_token'")
    names = REQUIRED_CREDENTIAL_ENV if auth_method == "password" else (ACCESS_TOKEN_ENV,)
    values = {name: str(prompt(f"{name}: ")).strip() for name in names}
    resolved = _validate_credential_values(values, auth_method)
    clear_idt_secret_environment()
    os.environ.update(values)
    os.environ[_EPHEMERAL_CREDENTIAL_FLAG] = "1"
    return {
        "credential_mode": "manual",
        "auth_method": resolved,
        "required_fields_complete": True,
    }


def configure_idt_credentials(
    *,
    mode: str = "path",
    path: str | Path = IDT_CREDENTIAL_PATH,
    auth_method: str | None = None,
    headless: bool = False,
    repository_root: str | Path | None = None,
    prompt: Any = getpass.getpass,
    include_path_in_status: bool = True,
) -> dict[str, object]:
    """Configure path or manual credentials and return only non-secret status."""
    if mode == "manual":
        return prompt_idt_credentials(
            auth_method=auth_method or "password", headless=headless, prompt=prompt
        )
    if mode != "path":
        raise ValueError("credential mode must be 'path' or 'manual'")
    credential_path = load_idt_credentials(
        path, auth_method=auth_method, repository_root=repository_root
    )
    resolved = "access_token" if os.environ.get(ACCESS_TOKEN_ENV) else "password"
    status: dict[str, object] = {
        "credential_mode": "path",
        "auth_method": resolved,
        "required_fields_complete": True,
        "permission_check": "passed",
    }
    if include_path_in_status:
        status["credential_path"] = str(credential_path)
    return status


def configure_idt_credentials_from_values(
    values: dict[str, str],
    *,
    auth_method: str | None = None,
) -> dict[str, object]:
    """Configure manually entered notebook credentials without persisting them.

    The returned status contains no credential values. Callers must clear the
    password widgets immediately after this function returns.
    """
    filtered = {
        name: str(value).strip()
        for name, value in values.items()
        if name in ALLOWED_CREDENTIAL_ENV and str(value).strip()
    }
    resolved = _validate_credential_values(filtered, auth_method)
    clear_idt_secret_environment()
    os.environ.update(filtered)
    os.environ[_EPHEMERAL_CREDENTIAL_FLAG] = "1"
    return {
        "credential_mode": "manual",
        "auth_method": resolved,
        "required_fields_complete": True,
    }


def configure_idt_credentials_from_bytes(
    payload: bytes,
    *,
    auth_method: str | None = None,
    maximum_bytes: int = 64 * 1024,
) -> dict[str, object]:
    """Parse an uploaded env file entirely in memory and retain no file data."""
    if len(payload) > maximum_bytes:
        raise ValueError("Uploaded credential env file is unexpectedly large")
    try:
        text = bytes(payload).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("Uploaded credential env file must be UTF-8 text") from exc
    values = _parse_credential_lines(text.splitlines())
    resolved = _validate_credential_values(values, auth_method)
    clear_idt_secret_environment()
    os.environ.update(values)
    os.environ[_EPHEMERAL_CREDENTIAL_FLAG] = "1"
    # Drop all direct references before returning a non-secret status.
    values.clear()
    text = ""
    return {
        "credential_mode": "upload",
        "auth_method": resolved,
        "required_fields_complete": True,
        "upload_retained": False,
    }


def clear_idt_secret_environment(*, keep_access_token: bool = False) -> None:
    """Remove IDT login material from the current process environment."""
    for name in REQUIRED_CREDENTIAL_ENV:
        os.environ.pop(name, None)
    if not keep_access_token:
        os.environ.pop(ACCESS_TOKEN_ENV, None)
    os.environ.pop(_EPHEMERAL_CREDENTIAL_FLAG, None)


def audit_artifacts_for_idt_secrets(
    paths: list[str | Path],
    credential_path: str | Path,
    *,
    repository_root: str | Path | None = None,
) -> dict[str, object]:
    """Scan text artifacts for exact credential values without reporting them."""
    source = validate_idt_credential_path(
        credential_path, repository_root=repository_root
    )
    values = _parse_credential_lines(source.read_text().splitlines())
    findings: list[dict[str, str]] = []
    scanned = 0
    for raw_path in paths:
        path = Path(raw_path)
        candidates = [path] if path.is_file() else (
            [item for item in path.rglob("*") if item.is_file()] if path.is_dir() else []
        )
        for candidate in candidates:
            try:
                content = candidate.read_text(errors="ignore")
            except OSError:
                continue
            scanned += 1
            for field, secret in values.items():
                if len(secret) >= 4 and secret in content:
                    findings.append({"artifact": str(candidate), "credential_field": field})
    return {
        "scanned_file_count": scanned,
        "finding_count": len(findings),
        "findings": findings,
        "passed": not findings,
    }


def _retrying_session() -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=4,
        connect=4,
        read=4,
        status=4,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"POST"}),
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def credentials_available() -> bool:
    return bool(os.environ.get(ACCESS_TOKEN_ENV)) or all(os.environ.get(name) for name in REQUIRED_CREDENTIAL_ENV)


def get_access_token(session: requests.Session | None = None, timeout: int = 30) -> str:
    supplied = os.environ.get(ACCESS_TOKEN_ENV)
    if supplied:
        if os.environ.get(_EPHEMERAL_CREDENTIAL_FLAG) == "1":
            clear_idt_secret_environment()
        return supplied
    missing = [name for name in REQUIRED_CREDENTIAL_ENV if not os.environ.get(name)]
    if missing:
        raise RuntimeError(f"Missing IDT OAuth environment variables: {', '.join(missing)}")
    client = session or _retrying_session()
    try:
        response = client.post(
            TOKEN_URL,
            auth=(os.environ["IDT_CLIENT_ID"], os.environ["IDT_CLIENT_SECRET"]),
            data={
                "grant_type": "password",
                "scope": "test",
                "username": os.environ["IDT_USERNAME"],
                "password": os.environ["IDT_PASSWORD"],
            },
            timeout=timeout,
        )
        response.raise_for_status()
        token = response.json().get("access_token")
    except (requests.RequestException, ValueError) as exc:
        raise RuntimeError("IDT OAuth authentication failed; credentials and response were redacted") from exc
    finally:
        # The scorer keeps only the bearer token in memory.  Login material is
        # removed after every OAuth attempt for both path and manual modes.
        clear_idt_secret_environment()
    if not token:
        raise RuntimeError("IDT OAuth response did not contain an access token")
    return str(token)


def screen_gblock_sequences(
    sequences: list[dict[str, str]],
    *,
    session: requests.Session | None = None,
    access_token: str | None = None,
    timeout: int = 120,
) -> Any:
    """Call the documented gBlocks complexity endpoint in one batch."""
    if not sequences:
        return []
    payload = [{"Name": str(item["Name"]), "Sequence": str(item["Sequence"]).upper()} for item in sequences]
    client = session or _retrying_session()
    token = access_token or get_access_token(client, timeout=min(timeout, 30))
    try:
        response = client.post(
            GBLOCK_COMPLEXITY_URL,
            headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json"},
            json=payload,
            timeout=timeout,
        )
        response.raise_for_status()
        return response.json()
    except (requests.RequestException, ValueError) as exc:
        raise RuntimeError("IDT complexity scoring failed; request authorization and response were redacted") from exc


def _flatten(value: Any, prefix: str = "") -> list[tuple[str, Any]]:
    rows: list[tuple[str, Any]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            child_prefix = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(_flatten(child, child_prefix))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            rows.extend(_flatten(child, f"{prefix}[{index}]"))
    else:
        rows.append((prefix.lower(), value))
    return rows


def _named_result(value: Any, name: str) -> Any | None:
    if isinstance(value, dict):
        lowered = {str(key).lower(): child for key, child in value.items()}
        if str(lowered.get("name", "")) == name:
            return value
        for child in value.values():
            found = _named_result(child, name)
            if found is not None:
                return found
    elif isinstance(value, list):
        for child in value:
            found = _named_result(child, name)
            if found is not None:
                return found
    return None


def _idt_rules(value: Any) -> list[dict[str, Any]]:
    """Return rule records carrying IDT's explicit ``IsViolated`` field."""
    rules: list[dict[str, Any]] = []
    if isinstance(value, dict):
        if isinstance(value.get("IsViolated"), bool):
            rules.append(value)
        else:
            for child in value.values():
                rules.extend(_idt_rules(child))
    elif isinstance(value, list):
        for child in value:
            rules.extend(_idt_rules(child))
    return rules


def summarize_complexity_response(
    response: Any,
    *,
    name: str | None = None,
    sequence_index: int | None = None,
) -> dict[str, Any]:
    """Summarize IDT gBlocks scores using the frozen score-only policy.

    The live endpoint returns one ordered rule list per submitted sequence.
    Orderability is based on the sum of every finite numeric rule ``Score``;
    ``IsViolated`` is retained only as diagnostic evidence.  An empty,
    position-matched rule list has score zero.  A returned rule with a missing
    or non-finite score makes the result unclassified rather than silently
    contributing zero.
    """
    selected = _named_result(response, name) if name is not None else None
    matched_by_name = selected is not None
    selected_by_index = False
    if selected is None and sequence_index is not None:
        if (
            isinstance(response, list)
            and response
            and all(isinstance(item, list) for item in response)
            and 0 <= sequence_index < len(response)
        ):
            selected = response[sequence_index]
            selected_by_index = True
    flattened = _flatten(selected if selected is not None else response)
    numeric_scores = {
        key: float(value)
        for key, value in flattened
        if "score" in key and isinstance(value, (int, float)) and not isinstance(value, bool)
    }
    pass_values: list[bool] = []
    issue_values: list[bool] = []
    for key, value in flattened:
        if not isinstance(value, bool):
            continue
        terminal = key.rsplit(".", 1)[-1]
        if any(word in terminal for word in ("pass", "accepted", "orderable", "synthesizable")):
            pass_values.append(value)
        if any(word in terminal for word in ("complex", "issue", "problem", "warning", "error")):
            issue_values.append(value)
    selected_value = selected if selected is not None else response
    rules = _idt_rules(selected_value)
    violations = [str(rule.get("Name", "unnamed_rule")) for rule in rules if rule["IsViolated"]]
    rule_scores: dict[str, float] = {}
    invalid_score_names: list[str] = []
    for index, rule in enumerate(rules):
        rule_name = str(rule.get("Name", f"rule_{index}"))
        score = rule.get("Score")
        if (
            isinstance(score, (int, float))
            and not isinstance(score, bool)
            and math.isfinite(float(score))
        ):
            rule_scores[rule_name] = rule_scores.get(rule_name, 0.0) + float(score)
        else:
            invalid_score_names.append(rule_name)
    rule_details = []
    for index, rule in enumerate(rules):
        threshold = rule.get("ThresholdOutput")
        threshold = threshold if isinstance(threshold, dict) else {}
        rule_details.append(
            {
                "name": str(rule.get("Name", f"rule_{index}")),
                "is_violated": bool(rule["IsViolated"]),
                "actual_value": rule.get("ActualValue"),
                "score": rule.get("Score"),
                "display_text": str(rule.get("DisplayText") or ""),
                "repeated_segment": str(rule.get("RepeatedSegment") or "").upper(),
                "forward_locations": [
                    int(value)
                    for value in (rule.get("ForwardLocations") or [])
                    if isinstance(value, (int, float)) and not isinstance(value, bool)
                ],
                "reverse_locations": [
                    int(value)
                    for value in (rule.get("ReverseLocations") or [])
                    if isinstance(value, (int, float)) and not isinstance(value, bool)
                ],
                "start_index": rule.get("StartIndex"),
                "terminal_end": rule.get("TerminalEnd"),
                "length": rule.get("Length"),
                "minimum_repeat_length": rule.get("MinimumRepeatLength"),
                "gc_percentage": rule.get("GCPercentage"),
                "repeat_percentage": rule.get("RepeatPercentage"),
                "rank": rule.get("Rank"),
                "threshold_type": threshold.get("ThresholdType"),
                "threshold_value": threshold.get("Value"),
                "threshold_window_length": threshold.get("WindowLength"),
                "threshold_min_length": threshold.get("MinLength"),
                "threshold_max_length": threshold.get("MaxLength"),
                "threshold_min_percentage": threshold.get("MinPercentage"),
                "threshold_max_percentage": threshold.get("MaxPercentage"),
                "threshold_quantity": threshold.get("Quantity"),
            }
        )
    score_total: float | None = None
    explicit_pass: bool | None = None
    score_complete = False
    if selected_by_index and isinstance(selected, list):
        score_complete = not invalid_score_names
        if score_complete:
            score_total = float(sum(rule_scores.values()))
            explicit_pass = score_total < IDT_ORDERABILITY_THRESHOLD
    elif rules:
        score_complete = not invalid_score_names
        if score_complete:
            score_total = float(sum(rule_scores.values()))
            explicit_pass = score_total < IDT_ORDERABILITY_THRESHOLD
    elif pass_values:
        # Retain compatibility with older named mock responses.  Production
        # gBlocks calls always take the score-list route above.
        explicit_pass = all(pass_values)
    elif issue_values:
        explicit_pass = not any(issue_values)
    positive_score_names = sorted(
        name for name, score in rule_scores.items() if score > 0
    )
    return {
        "idt_status": "passed" if explicit_pass is True else "failed" if explicit_pass is False else "scored_unclassified",
        "idt_explicit_pass": explicit_pass,
        "idt_numeric_scores_json": json.dumps(numeric_scores, sort_keys=True),
        "idt_rule_scores_json": json.dumps(rule_scores, sort_keys=True),
        "idt_rule_details_json": json.dumps(rule_details, sort_keys=True),
        "idt_rule_count": len(rules),
        "idt_violation_count": len(violations),
        "idt_violation_names_json": json.dumps(violations),
        "idt_positive_score_names_json": json.dumps(positive_score_names),
        "idt_invalid_score_names_json": json.dumps(invalid_score_names),
        "idt_score_complete": score_complete,
        "idt_complexity_score": score_total,
        "idt_score_aggregation": "sum",
        "idt_orderability_threshold": IDT_ORDERABILITY_THRESHOLD,
        "idt_score_policy": IDT_SCORE_POLICY,
        "idt_result_matched_by_name": matched_by_name,
        "idt_result_selected_by_index": selected_by_index,
    }


class IDTComplexityScorer:
    """Stateful, credential-safe scorer used inside adaptive optimization.

    OAuth is acquired once per process, exact-DNA results are cached by SHA256,
    and every real HTTP response is appended to a JSONL audit without storing
    the submitted sequence itself.
    """

    def __init__(self, audit_path: str | Path, *, timeout: int = 120) -> None:
        self.audit_path = Path(audit_path)
        self.audit_path.parent.mkdir(parents=True, exist_ok=True)
        # Every Stage-2 shard has an audit artifact, including shards that
        # never reach a locally valid DNA candidate.  This keeps finalization
        # cardinality exact without inventing an API response.
        self.audit_path.touch(exist_ok=True)
        self.timeout = timeout
        self.session = _retrying_session()
        self.access_token: str | None = None
        self.cache: dict[str, dict[str, Any]] = {}
        self.api_attempts = 0
        self.api_calls = 0
        self.cache_hits = 0
        if self.audit_path.is_file():
            for line in self.audit_path.read_text().splitlines():
                try:
                    record = json.loads(line)
                    sequence_sha = str(record["request"]["sha256"])
                    summary = record["summary"]
                except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                    continue
                if isinstance(summary, dict) and sequence_sha:
                    self.cache[sequence_sha] = deepcopy(summary)

    def _screen(self, item: dict[str, str]) -> Any:
        if self.access_token is None:
            self.access_token = get_access_token(
                self.session, timeout=min(self.timeout, 30)
            )
        try:
            return screen_gblock_sequences(
                [item],
                session=self.session,
                access_token=self.access_token,
                timeout=self.timeout,
            )
        except requests.HTTPError as exc:
            if exc.response is None or exc.response.status_code != 401:
                raise
            self.access_token = get_access_token(
                self.session, timeout=min(self.timeout, 30)
            )
            return screen_gblock_sequences(
                [item],
                session=self.session,
                access_token=self.access_token,
                timeout=self.timeout,
            )

    def score(self, name: str, sequence: str) -> dict[str, Any]:
        normalized = str(sequence).upper()
        sequence_sha = hashlib.sha256(normalized.encode()).hexdigest()
        if sequence_sha in self.cache:
            self.cache_hits += 1
            cached = deepcopy(self.cache[sequence_sha])
            cached["idt_cache_hit"] = True
            cached["idt_request_name"] = str(name)
            return cached

        item = {"Name": str(name), "Sequence": normalized}
        self.api_attempts += 1
        response = self._screen(item)
        if not (
            isinstance(response, list)
            and len(response) == 1
            and isinstance(response[0], list)
        ):
            raise RuntimeError(
                "IDT single-sequence response was not one ordered rule list"
            )
        summary = summarize_complexity_response(
            response, name=str(name), sequence_index=0
        )
        response_sha = hashlib.sha256(
            json.dumps(
                response, sort_keys=True, separators=(",", ":"), default=str
            ).encode()
        ).hexdigest()
        summary.update(
            {
                "idt_api_called": True,
                "idt_cache_hit": False,
                "idt_request_name": str(name),
                "idt_scored_sequence_length_bp": len(normalized),
                "idt_scored_sequence_sha256": sequence_sha,
                "idt_scored_sequence_unchanged": True,
                "idt_raw_response_file": str(self.audit_path),
                "idt_response_sha256": response_sha,
            }
        )
        record = {
            "endpoint": GBLOCK_COMPLEXITY_URL,
            "retrieved_unix_time": time.time(),
            "request": {
                "name": str(name),
                "length_bp": len(normalized),
                "sha256": sequence_sha,
            },
            "summary": summary,
            "response_sha256": response_sha,
            "response": response,
        }
        with self.audit_path.open("a") as handle:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
        self.api_calls += 1
        self.cache[sequence_sha] = deepcopy(summary)
        return summary


def write_cached_response(
    response: Any,
    destination: str | Path,
    *,
    sequences: list[dict[str, str]] | None = None,
) -> Path:
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "endpoint": GBLOCK_COMPLEXITY_URL,
        "retrieved_unix_time": time.time(),
        "request_sequences": [
            {
                "name": str(item["Name"]),
                "length_bp": len(str(item["Sequence"])),
                "sha256": hashlib.sha256(str(item["Sequence"]).encode()).hexdigest(),
            }
            for item in (sequences or [])
        ],
        "response": response,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path
