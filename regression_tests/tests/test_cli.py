"""Smoke checks for the public package, without running the scientific solver."""

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from tests.fixtures.harness import create_harness, run_command  # noqa: E402


def test_help_and_discovery_need_no_machine_setup():
    overview = run_command()
    assert overview.returncode == 0, overview.stderr
    assert "usage:" in overview.stdout
    # -S excludes installed scientific packages: even a fresh Python can show help.
    help_result = subprocess.run(
        [sys.executable, "-B", "-S", "-m", "regression_tests", "check", "--help"],
        cwd=ROOT.parent, capture_output=True, text=True,
    )
    assert help_result.returncode == 0, help_result.stderr
    assert "--cases" not in help_result.stdout
    environment = {"MHDG_REGRESSION_SETTINGS": "/nonexistent/settings.env"}
    for listing in ("cases", "suites", "layouts", "workflows"):
        result = run_command("list", listing, environment=environment)
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip()
    case = run_command("list", "cases").stdout.split(":", 1)[0]
    result = run_command("list", "workflows", case, environment=environment)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip()


def test_usage_and_runtime_failures_have_distinct_exit_codes(tmp_path):
    usage = run_command("check", "--cases", str(tmp_path))
    assert usage.returncode == 2
    assert "unrecognized arguments" in usage.stderr
    missing = run_command("check", "--settings", str(tmp_path / "missing.env"))
    assert missing.returncode == 1
    assert missing.stderr.startswith("error: cannot read settings file")
    assert "Traceback" not in missing.stderr


def test_debug_prepare_and_bundle_validation_use_discovered_catalogs(tmp_path):
    harness = create_harness(tmp_path)
    result = run_command("bundle", "validate", str(harness.bundle))
    assert result.returncode == 0, result.stderr
    assert "bundle valid:" in result.stdout
    result = run_command(
        "prepare", "legacy_case", "cold_step_adaptive", "--layout", "serial_omp1",
        "--run-id", "cli-smoke", "--settings", str(harness.settings),
    )
    assert result.returncode == 0, result.stderr
    run = harness.run_directory("cold_step_adaptive", "serial_omp1", "cli-smoke")
    plan = json.loads((run / "run_plan.json").read_text())
    assert plan["stages"]
    assert "run prepared:" in result.stdout
    assert not (run / "run_metadata.json").exists()
