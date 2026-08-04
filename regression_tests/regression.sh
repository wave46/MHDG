#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PYTHON_COMMAND=${PYTHON:-python3}

usage() {
  cat <<'EOF'
Usage:
  regression_tests/regression.sh help
  regression_tests/regression.sh bundle create --case CASE --source DIR --output DIR
  regression_tests/regression.sh bundle validate --settings FILE
  regression_tests/regression.sh bundle promote SUITE_SUMMARY [SUITE_SUMMARY ...] --settings FILE --output DIR --bundle-version VERSION
  regression_tests/regression.sh build --settings FILE [--jobs N]
  regression_tests/regression.sh prepare CASE WORKFLOW --layout LAYOUT --settings FILE
  regression_tests/regression.sh run CASE WORKFLOW [WORKFLOW ...] --layout LAYOUT --settings FILE
  regression_tests/regression.sh compare RUN_DIRECTORY
  regression_tests/regression.sh suite run SUITE --settings FILE [--run-only] [--resume]
  regression_tests/regression.sh suite compare SUITE_SUMMARY
  regression_tests/regression.sh suite check [SUITE] [--settings FILE] [--build]
  regression_tests/regression.sh golden update CASE --settings FILE --run-id ID --output DIR --bundle-version VERSION [--only COMPONENT] [--bootstrap-candidate] [--retry-failed]
  regression_tests/regression.sh golden status WORKSPACE

Commands:
  bundle create    Create a validated candidate bundle from prepared files.
  bundle validate  Validate a configured external bundle without changing it.
  bundle promote   Create a golden bundle from accepted suite results.
  build            Build clean serial and parallel regression executables.
  prepare          Prepare one isolated workflow without running the solver.
  run              Prepare and execute one or more workflows.
  compare          Compare one completed run using its recorded policy.
  suite run        Execute every workflow-layout cell in a tracked suite.
  suite compare    Recompare saved suite results without running the solver.
  suite check      Run a suite against a required golden bundle (default: warm).
  golden update    Start or continue an ordered golden-reference update.
  golden status    Show persisted golden-update state.

Suites: warm, impurity_scalar_baseline, impurity_mixture, impurity_references,
        initialization_smoke, stored_field_compatibility, race, cold,
        warm_parallelism, race_matrix, cold_matrix, diverted_warm,
        diverted_warm_parallelism, diverted_cold_adaptive,
        diverted_race_matrix.
EOF
}

fail() {
  echo "error: $*" >&2
  usage >&2
  exit 2
}

run_python() {
  local tool=$1
  shift
  exec "$PYTHON_COMMAND" "$SCRIPT_DIR/tools/$tool" "$@"
}

run_bundle_command() {
  local action=${1:-}
  if (($# > 0)); then
    shift
  fi
  case "$action" in
    create)
      run_python create_bundle.py --cases "$SCRIPT_DIR/cases" "$@"
      ;;
    validate)
      run_python check_bundle.py --cases "$SCRIPT_DIR/cases" "$@"
      ;;
    promote)
      run_python promote_bundle.py --cases "$SCRIPT_DIR/cases" "$@"
      ;;
    help|-h|--help)
      usage
      ;;
    *)
      fail "expected 'bundle create', 'bundle validate', or 'bundle promote'"
      ;;
  esac
}

run_suite_command() {
  local action=${1:-}
  if (($# > 0)); then
    shift
  fi
  case "$action" in
    run)
      run_python run_suite.py \
        --cases "$SCRIPT_DIR/cases" \
        --layouts "$SCRIPT_DIR/layouts.json" \
        --suites "$SCRIPT_DIR/suites.json" \
        --tolerances "$SCRIPT_DIR/tolerances.json" \
        "$@"
      ;;
    compare)
      run_python verify_suite.py \
        --cases "$SCRIPT_DIR/cases" \
        --tolerances "$SCRIPT_DIR/tolerances.json" \
        "$@"
      ;;
    check)
      run_golden_suite "$@"
      ;;
    help|-h|--help)
      usage
      ;;
    *)
      fail "expected 'suite run', 'suite compare', or 'suite check'"
      ;;
  esac
}

run_golden_suite() {
  local suite=warm
  if (($# > 0)) && [[ "$1" != -* ]]; then
    suite=$1
    shift
  fi
  local settings=${MHDG_REGRESSION_GOLDEN_SETTINGS:-$SCRIPT_DIR/golden.local.env}
  run_python run_suite.py \
    --settings "$settings" \
    --cases "$SCRIPT_DIR/cases" \
    --layouts "$SCRIPT_DIR/layouts.json" \
    --suites "$SCRIPT_DIR/suites.json" \
    --tolerances "$SCRIPT_DIR/tolerances.json" \
    --require-bundle-class golden \
    "$suite" \
    "$@"
}

run_golden_command() {
  local action=${1:-}
  case "$action" in
    update|status)
      run_python golden_update.py "$@"
      ;;
    help|-h|--help)
      run_python golden_update.py --help
      ;;
    *)
      fail "expected 'golden update' or 'golden status'"
      ;;
  esac
}

if (($# == 0)); then
  usage
  exit 0
fi

command=$1
shift
case "$command" in
  help|-h|--help)
    usage
    ;;
  bundle)
    run_bundle_command "$@"
    ;;
  build)
    run_python build_solver.py "$@"
    ;;
  prepare)
    run_python prepare_run.py \
      --cases "$SCRIPT_DIR/cases" \
      --layouts "$SCRIPT_DIR/layouts.json" \
      "$@"
    ;;
  run)
    run_python run_case.py \
      --cases "$SCRIPT_DIR/cases" \
      --layouts "$SCRIPT_DIR/layouts.json" \
      "$@"
    ;;
  compare)
    run_python compare.py \
      --cases "$SCRIPT_DIR/cases" \
      --tolerances "$SCRIPT_DIR/tolerances.json" \
      "$@"
    ;;
  suite)
    run_suite_command "$@"
    ;;
  golden)
    run_golden_command "$@"
    ;;
  *)
    fail "unknown command: $command"
    ;;
esac
