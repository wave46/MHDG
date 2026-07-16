#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
Usage:
  regression_tests/regression.sh help
  regression_tests/regression.sh --help
  regression_tests/regression.sh --settings FILE check-data

Available commands:
  help        Show this help text.
  check-data  Validate an external bundle without modifying it.

Options:
  --settings FILE  Local settings file containing MHDG_REGRESSION_DATA_ROOT.
  -h, --help       Show this help text.

Solver execution and result comparison are not implemented yet.
EOF
}

settings_file=""
command=""

while (($# > 0)); do
  case "$1" in
    help|-h|--help)
      usage
      exit 0
      ;;
    --settings)
      if (($# < 2)); then
        echo "error: --settings requires a file" >&2
        exit 2
      fi
      settings_file=$2
      shift 2
      ;;
    --settings=*)
      settings_file=${1#*=}
      shift
      ;;
    check-data)
      if [[ -n "$command" ]]; then
        echo "error: only one command may be specified" >&2
        exit 2
      fi
      command=$1
      shift
      ;;
    *)
      echo "error: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$command" ]]; then
  usage
  exit 0
fi

if [[ -z "$settings_file" ]]; then
  echo "error: check-data requires --settings FILE" >&2
  exit 2
fi

python_command=${PYTHON:-python3}
exec "$python_command" "$SCRIPT_DIR/tools/check_bundle.py" \
  --settings "$settings_file" \
  --cases "$SCRIPT_DIR/cases"
