#!/usr/bin/env bash
set -euo pipefail

host="127.0.0.1"
port="2718"
open_browser=true

usage() {
    cat <<'EOF'
Usage: ./scripts/start_hurdler_web.sh [--host HOST] [--port PORT] [--no-browser]

Start the local annotation-aware HURDLER Marimo designer.  The default page is
http://127.0.0.1:2718 and opens in the local browser automatically.
EOF
}

while (($#)); do
    case "$1" in
        --host)
            [[ $# -ge 2 ]] || { echo "--host requires a value" >&2; exit 2; }
            host="$2"
            shift 2
            ;;
        --port)
            [[ $# -ge 2 ]] || { echo "--port requires a value" >&2; exit 2; }
            port="$2"
            shift 2
            ;;
        --no-browser)
            open_browser=false
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

python_bin="${HURDLER_PYTHON:-python}"
command -v "$python_bin" >/dev/null 2>&1 || {
    echo "Python was not found. Activate the HURDLER environment first." >&2
    exit 1
}
"$python_bin" -c 'import hurdler, marimo' >/dev/null 2>&1 || {
    echo 'HURDLER/Marimo is not installed in this Python environment.' >&2
    echo 'Run: python -m pip install -e ".[notebooks,optimization]"' >&2
    exit 1
}

command=("$python_bin" -m hurdler web --host "$host" --port "$port")
if [[ "$open_browser" == false ]]; then
    command+=(--no-browser)
fi
echo "Starting HURDLER designer at http://${host}:${port}"
exec "${command[@]}"
