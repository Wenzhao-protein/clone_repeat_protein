from __future__ import annotations

import os
from pathlib import Path

import pytest

from hurdler.idt import (
    ACCESS_TOKEN_ENV,
    REQUIRED_CREDENTIAL_ENV,
    audit_artifacts_for_idt_secrets,
    clear_idt_secret_environment,
    configure_idt_credentials,
    configure_idt_credentials_from_bytes,
    configure_idt_credentials_from_values,
    get_access_token,
    load_idt_credentials,
    validate_idt_credential_path,
)


@pytest.fixture(autouse=True)
def clear_credentials_after_test():
    clear_idt_secret_environment()
    yield
    clear_idt_secret_environment()


def _secure(path: Path, text: str) -> Path:
    path.write_text(text)
    path.chmod(0o600)
    return path


def test_path_mode_parses_export_comments_quotes_and_access_token(tmp_path):
    path = _secure(
        tmp_path / "external.env",
        '# comment\nexport IDT_ACCESS_TOKEN="temporary token" # comment\n',
    )
    status = configure_idt_credentials(
        mode="path",
        path=path,
        auth_method="access_token",
        repository_root=tmp_path / "different-repository",
    )
    assert status["credential_mode"] == "path"
    assert status["auth_method"] == "access_token"
    assert os.environ[ACCESS_TOKEN_ENV] == "temporary token"
    assert "temporary token" not in str(status)


def test_password_and_token_formats_are_mutually_exclusive(tmp_path):
    path = _secure(
        tmp_path / "mixed.env",
        "\n".join(
            [
                "IDT_CLIENT_ID=id",
                "IDT_CLIENT_SECRET=secret",
                "IDT_USERNAME=user",
                "IDT_PASSWORD=password",
                "IDT_ACCESS_TOKEN=token",
            ]
        ),
    )
    with pytest.raises(RuntimeError, match="either"):
        load_idt_credentials(path, repository_root=tmp_path / "repo")


def test_repository_internal_and_permissive_files_are_rejected(tmp_path):
    repository = tmp_path / "repo"
    repository.mkdir()
    internal = _secure(repository / "idt.env", "IDT_ACCESS_TOKEN=token\n")
    with pytest.raises(ValueError, match="outside"):
        validate_idt_credential_path(internal, repository_root=repository)

    external = tmp_path / "world-readable.env"
    external.write_text("IDT_ACCESS_TOKEN=token\n")
    external.chmod(0o644)
    if os.name == "posix":
        with pytest.raises(PermissionError, match="600"):
            validate_idt_credential_path(external, repository_root=repository)


def test_manual_mode_is_in_memory_and_fails_headless():
    with pytest.raises(RuntimeError, match="headless"):
        configure_idt_credentials(mode="manual", auth_method="access_token", headless=True)

    supplied = iter(["token-value"])
    status = configure_idt_credentials(
        mode="manual",
        auth_method="access_token",
        headless=False,
        prompt=lambda _label: next(supplied),
    )
    assert status == {
        "credential_mode": "manual",
        "auth_method": "access_token",
        "required_fields_complete": True,
    }
    assert "token-value" not in str(status)
    assert get_access_token() == "token-value"
    assert ACCESS_TOKEN_ENV not in os.environ


def test_password_login_material_is_removed_after_oauth_attempt(monkeypatch):
    for name in REQUIRED_CREDENTIAL_ENV:
        monkeypatch.setenv(name, f"value-for-{name}")

    class Response:
        def raise_for_status(self):
            return None

        def json(self):
            return {"access_token": "bearer"}

    class Session:
        def post(self, *args, **kwargs):
            return Response()

    assert get_access_token(Session()) == "bearer"
    assert all(name not in os.environ for name in REQUIRED_CREDENTIAL_ENV)


def test_artifact_audit_reports_field_name_without_echoing_secret(tmp_path):
    credentials = _secure(tmp_path / "credentials.env", "IDT_ACCESS_TOKEN=do-not-copy-this\n")
    clean = tmp_path / "clean.log"
    clean.write_text("authentication passed\n")
    leaked = tmp_path / "bad.log"
    leaked.write_text("do-not-copy-this\n")
    report = audit_artifacts_for_idt_secrets(
        [clean, leaked], credentials, repository_root=tmp_path / "repo"
    )
    assert not report["passed"] and report["finding_count"] == 1
    assert report["findings"][0]["credential_field"] == "IDT_ACCESS_TOKEN"
    assert "do-not-copy-this" not in str(report)


def test_manual_values_and_uploaded_env_are_memory_only():
    status = configure_idt_credentials_from_values(
        {"IDT_ACCESS_TOKEN": "manual-value"}, auth_method="access_token"
    )
    assert status["credential_mode"] == "manual"
    assert "manual-value" not in str(status)
    clear_idt_secret_environment()

    payload = bytearray(b"IDT_ACCESS_TOKEN=uploaded-value\n")
    status = configure_idt_credentials_from_bytes(bytes(payload), auth_method="access_token")
    for index in range(len(payload)):
        payload[index] = 0
    assert status["credential_mode"] == "upload"
    assert status["upload_retained"] is False
    assert "uploaded-value" not in str(status)


def test_path_status_can_omit_the_actual_path(tmp_path):
    path = _secure(tmp_path / "private.env", "IDT_ACCESS_TOKEN=token\n")
    status = configure_idt_credentials(
        mode="path",
        path=path,
        auth_method="access_token",
        repository_root=tmp_path / "repo",
        include_path_in_status=False,
    )
    assert "credential_path" not in status
