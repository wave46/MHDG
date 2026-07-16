#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
Usage:
  regression_tests/regression.sh help
  regression_tests/regression.sh --help
  regression_tests/regression.sh bundle create --case CASE --source DIR --output DIR
  regression_tests/regression.sh --settings FILE check-data
  regression_tests/regression.sh --settings FILE prepare CASE WORKFLOW --layout LAYOUT

Available commands:
  help           Show this help text.
  bundle create  Create and validate a bundle from prepared case files.
  check-data     Validate an external bundle without modifying it.
  prepare        Create an isolated run directory without executing the solver.

Options:
  --settings FILE  Local bundle, run, executable, and launcher settings.
  -h, --help       Show this help text.

Solver execution and result comparison are not implemented yet.
EOF
}

settings_file=""
command=""
bundle_arguments=()
prepare_arguments=()

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
    bundle)
      if [[ -n "$command" ]]; then
        echo "error: only one command may be specified" >&2
        exit 2
      fi
      if (($# < 2)) || [[ "$2" != "create" ]]; then
        echo "error: expected 'bundle create'" >&2
        exit 2
      fi
      command="bundle-create"
      shift 2
      bundle_arguments=("$@")
      break
      ;;
    prepare)
      if [[ -n "$command" ]]; then
        echo "error: only one command may be specified" >&2
        exit 2
      fi
      command="prepare"
      shift
      prepare_arguments=("$@")
      break
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

python_command=${PYTHON:-python3}

case "$command" in
  check-data)
    if [[ -z "$settings_file" ]]; then
      echo "error: check-data requires --settings FILE" >&2
      exit 2
    fi
    exec "$python_command" "$SCRIPT_DIR/tools/check_bundle.py" \
      --settings "$settings_file" \
      --cases "$SCRIPT_DIR/cases"
    ;;
  bundle-create)
    if [[ -n "$settings_file" ]]; then
      echo "error: bundle create does not use --settings" >&2
      exit 2
    fi
    exec "$python_command" "$SCRIPT_DIR/tools/create_bundle.py" \
      --cases "$SCRIPT_DIR/cases" \
      "${bundle_arguments[@]}"
    ;;
  prepare)
    if [[ -z "$settings_file" ]]; then
      echo "error: prepare requires --settings FILE" >&2
      exit 2
    fi
    exec "$python_command" "$SCRIPT_DIR/tools/prepare_run.py" \
      --settings "$settings_file" \
      --cases "$SCRIPT_DIR/cases" \
      --layouts "$SCRIPT_DIR/layouts.json" \
      "${prepare_arguments[@]}"
    ;;
esac
