#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<'EOF'
Usage:
  regression_tests/regression.sh help
  regression_tests/regression.sh --help
  regression_tests/regression.sh bundle create --case CASE --source DIR --output DIR
  regression_tests/regression.sh --settings FILE bundle promote SUITE_SUMMARY --output DIR --bundle-version VERSION
  regression_tests/regression.sh --settings FILE build [--jobs N]
  regression_tests/regression.sh --settings FILE check-data
  regression_tests/regression.sh --settings FILE prepare CASE WORKFLOW --layout LAYOUT
  regression_tests/regression.sh --settings FILE run CASE WORKFLOW [WORKFLOW ...] --layout LAYOUT
  regression_tests/regression.sh compare RUN_DIRECTORY
  regression_tests/regression.sh --settings FILE suite SUITE [--run-only] [--resume]
  regression_tests/regression.sh golden-check [SUITE] [--build] [--build-jobs N]

Available commands:
  help           Show this help text.
  bundle create  Create and validate a bundle from prepared case files.
  bundle promote Create a complete golden bundle from a passing canonical run.
  build          Build clean serial and parallel regression executables.
  check-data     Validate an external bundle without modifying it.
  prepare        Create an isolated run directory without executing the solver.
  run            Prepare and execute one or more isolated workflows.
  compare        Compare a completed fixed-mesh run with its reference.
  suite          Run every workflow/layout in a tracked suite.
  golden-check   Check executables against the configured golden bundle.

Options:
  --settings FILE  Local bundle, run, executable, and launcher settings.
  -h, --help       Show this help text.

Suites: warm, warm_parallelism, cold_matrix.
EOF
}

settings_file=""
command=""
bundle_arguments=()
promotion_arguments=()
build_arguments=()
prepare_arguments=()
run_arguments=()
compare_arguments=()
suite_arguments=()
golden_arguments=()

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
    build)
      if [[ -n "$command" ]]; then
        echo "error: only one command may be specified" >&2
        exit 2
      fi
      command="build"
      shift
      build_arguments=("$@")
      break
      ;;
    bundle)
      if [[ -n "$command" ]]; then
        echo "error: only one command may be specified" >&2
        exit 2
      fi
      if (($# < 2)) || [[ "$2" != "create" && "$2" != "promote" ]]; then
        echo "error: expected 'bundle create' or 'bundle promote'" >&2
        exit 2
      fi
      command="bundle-$2"
      shift 2
      if [[ "$command" == "bundle-create" ]]; then
        bundle_arguments=("$@")
      else
        promotion_arguments=("$@")
      fi
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
    run)
      if [[ -n "$command" ]]; then
        echo "error: only one command may be specified" >&2
        exit 2
      fi
      command="run"
      shift
      run_arguments=("$@")
      break
      ;;
    compare)
      if [[ -n "$command" ]]; then
        echo "error: only one command may be specified" >&2
        exit 2
      fi
      command="compare"
      shift
      compare_arguments=("$@")
      break
      ;;
    suite)
      if [[ -n "$command" ]]; then
        echo "error: only one command may be specified" >&2
        exit 2
      fi
      command="suite"
      shift
      suite_arguments=("$@")
      break
      ;;
    golden-check)
      if [[ -n "$command" ]]; then
        echo "error: only one command may be specified" >&2
        exit 2
      fi
      command="golden-check"
      shift
      golden_arguments=("$@")
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
  bundle-promote)
    if [[ -z "$settings_file" ]]; then
      echo "error: bundle promote requires --settings FILE" >&2
      exit 2
    fi
    exec "$python_command" "$SCRIPT_DIR/tools/promote_bundle.py" \
      --settings "$settings_file" \
      --cases "$SCRIPT_DIR/cases" \
      "${promotion_arguments[@]}"
    ;;
  build)
    if [[ -z "$settings_file" ]]; then
      echo "error: build requires --settings FILE" >&2
      exit 2
    fi
    exec "$python_command" "$SCRIPT_DIR/tools/build_solver.py" \
      --settings "$settings_file" \
      "${build_arguments[@]}"
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
  run)
    if [[ -z "$settings_file" ]]; then
      echo "error: run requires --settings FILE" >&2
      exit 2
    fi
    exec "$python_command" "$SCRIPT_DIR/tools/run_case.py" \
      --settings "$settings_file" \
      --cases "$SCRIPT_DIR/cases" \
      --layouts "$SCRIPT_DIR/layouts.json" \
      "${run_arguments[@]}"
    ;;
  compare)
    if [[ -n "$settings_file" ]]; then
      echo "error: compare does not use --settings" >&2
      exit 2
    fi
    exec "$python_command" "$SCRIPT_DIR/tools/compare_run.py" \
      --cases "$SCRIPT_DIR/cases" \
      --tolerances "$SCRIPT_DIR/tolerances.json" \
      "${compare_arguments[@]}"
    ;;
  suite)
    if [[ -z "$settings_file" ]]; then
      echo "error: suite requires --settings FILE" >&2
      exit 2
    fi
    exec "$python_command" "$SCRIPT_DIR/tools/run_suite.py" \
      --settings "$settings_file" \
      --cases "$SCRIPT_DIR/cases" \
      --layouts "$SCRIPT_DIR/layouts.json" \
      --suites "$SCRIPT_DIR/suites.json" \
      --tolerances "$SCRIPT_DIR/tolerances.json" \
      "${suite_arguments[@]}"
    ;;
  golden-check)
    if [[ -z "$settings_file" ]]; then
      settings_file=${MHDG_REGRESSION_GOLDEN_SETTINGS:-$SCRIPT_DIR/golden.local.env}
    fi
    if [[ ! -f "$settings_file" ]]; then
      echo "error: golden settings file not found: $settings_file" >&2
      echo "copy settings.example.env to golden.local.env or set MHDG_REGRESSION_GOLDEN_SETTINGS" >&2
      exit 2
    fi
    golden_suite=warm
    if ((${#golden_arguments[@]} > 0)) && [[ "${golden_arguments[0]}" != -* ]]; then
      golden_suite=${golden_arguments[0]}
      golden_arguments=("${golden_arguments[@]:1}")
    fi
    exec "$python_command" "$SCRIPT_DIR/tools/run_suite.py" \
      --settings "$settings_file" \
      --cases "$SCRIPT_DIR/cases" \
      --layouts "$SCRIPT_DIR/layouts.json" \
      --suites "$SCRIPT_DIR/suites.json" \
      --tolerances "$SCRIPT_DIR/tolerances.json" \
      --require-bundle-class golden \
      "$golden_suite" \
      "${golden_arguments[@]}"
    ;;
esac
