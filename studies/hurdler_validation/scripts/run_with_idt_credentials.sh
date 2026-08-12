#!/usr/bin/env bash
set -euo pipefail

idt_credential_file="${IDT_CREDENTIAL_FILE:-${XDG_CONFIG_HOME:-${HOME}/.config}/hurdler/idt.env}"
if [[ -f "$idt_credential_file" ]]; then
    idt_mode="$(stat -c '%a' "$idt_credential_file")"
    if (( 10#$idt_mode % 100 != 0 )); then
        echo "Refusing IDT credential file with group/other permissions: $idt_credential_file" >&2
        exit 77
    fi
    set -a
    # shellcheck disable=SC1090
    source "$idt_credential_file"
    set +a
fi

exec "$@"
