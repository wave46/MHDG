"""Build selection, object isolation and provenance with a fake make process."""

import json
from pathlib import Path
import os
import subprocess
import sys
from types import SimpleNamespace

import pytest

from regression_tests import build, cli, config
from regression_tests.support import BundleError
from regression_tests.files import file_identity
from regression_tests.tests.fixtures.harness import create_harness


@pytest.fixture
def setup(tmp_path, monkeypatch):
    repository = tmp_path / "repository"
    library = repository / "lib"
    library.mkdir(parents=True)
    (repository / "test").mkdir()
    (repository / "test" / build.RUNTIME_FILE).write_text("fixture nodes\n")
    (repository / ".gitignore").write_text("*.o\nMHDG-*\n")
    fake_bin = repository / "fake-bin"
    fake_bin.mkdir()
    make = fake_bin / "make"
    make.write_text(FAKE_MAKE)
    make.chmod(0o755)
    log = tmp_path / "make.log"
    environment = {
        "PATH": f"{fake_bin}:{os.environ['PATH']}",
        "MHDG_FAKE_LOG": str(log), "MHDG_FAKE_ENV": "loaded",
        "MHDG_FIXTURE_PYTHON": sys.executable,
        "MHDG_FIXTURE_CODE": (f"import sys; sys.path.insert(0, {str(Path(__file__).resolve().parents[2])!r}); "
                              "from regression_tests.tests.fixtures.solutions import finalize_solver_outputs; "
                              "finalize_solver_outputs(*sys.argv[1:4], dirty=sys.argv[4]=='true')"),
    }
    for name, value in environment.items():
        monkeypatch.setenv(name, value)
    script = library / "environment setup.sh"
    script.write_text('test -z "$MHDG_OPTIONAL_UNSET"\nexport MHDG_FAKE_ENV=sourced\n')
    for command in (
        ["git", "init", "-q"], ["git", "add", "."],
        ["git", "-c", "user.name=Regression Test", "-c",
         "user.email=regression@example.invalid", "commit", "-qm", "fixture"],
    ):
        subprocess.run(command, cwd=repository, check=True, capture_output=True)
    return SimpleNamespace(
        repository=repository, log=log, script=script,
        settings={"MHDG_REGRESSION_BUILD_ROOT": str(tmp_path / "builds")},
    )


@pytest.mark.parametrize("variants,existing_objects,source_script", [
    (("parallel", "serial", "parallel"), True, True),
    (("parallel",), False, False),
    (("serial",), True, False),
    (("serial", "parallel", "gamma", "gamma"), True, False),
])
def test_required_variants_cleaning_and_manifest_selection(setup, variants, existing_objects, source_script):
    if existing_objects:
        (setup.repository / "lib/stale.o").touch()
    if source_script:
        setup.settings["MHDG_ENVIRONMENT_SCRIPT"] = str(setup.script)
    requirements = [("NGammaTiTeNeutralGamma" if v == "gamma" else config.DEFAULT_MODEL,
                     "mpi" if v == "parallel" else "serial") for v in variants]
    result = build.build_solver(setup.settings, setup.repository, jobs=3, requirements=requirements)
    metadata = json.loads(result.metadata_path.read_text())
    commands = metadata["commands"]
    compiles = [command for command in commands if command[1] != "clean"]
    assert len(compiles) == len(set(variants))
    assert sum(command == ["make", "clean"] for command in commands) == (
        int(existing_objects) + len(compiles) - 1
    )
    if existing_objects:
        assert commands[0] == ["make", "clean"]
        assert not (setup.repository / "lib/stale.o").exists()
    if len(compiles) == 2:
        assert commands[commands.index(compiles[1]) - 1] == ["make", "clean"]
    for command in compiles:
        assert {"-j3", "DIM=2D", "COMPTYPE=opt"} <= set(command)
    assert {(next(a[4:] for a in c if a.startswith("MDL=")),
             "mpi" if "MODE=parall" in c else "serial") for c in compiles} == set(requirements)
    assert metadata["status"] == "completed"
    assert not metadata["repository"]["dirty"]
    assert bool(metadata["environment_script"]) == source_script
    expected_env = "sourced" if source_script else "loaded"
    assert all(line.startswith(expected_env + "|") for line in setup.log.read_text().splitlines())
    for variant, executable in result.executables.items():
        assert metadata["artifacts"][variant] == {"path": str(executable), **file_identity(executable)}
    assert (result.path / "bin" / build.RUNTIME_FILE).read_bytes() == (
        setup.repository / "test" / build.RUNTIME_FILE
    ).read_bytes()

    # An explicit partial build replaces the default build as a whole.
    old_settings = setup.repository / "machine.json"
    old_settings.write_text(json.dumps({"defaults": {"build": "unavailable-old-build.json"}}))
    selected = config.settings(old_settings, build_manifest=result.metadata_path)
    assert selected["MHDG_EXECUTABLES"] == {key: str(path) for key, path in result.executables.items()}
    missing = ("NGammaTiTeNeutralGamma", "mpi")
    with pytest.raises(BundleError, match="no build selected"):
        config.runtime_settings(selected, {missing})
    removed = next(iter(metadata["artifacts"]))
    del metadata["artifacts"][removed]
    if metadata["artifacts"]:
        result.metadata_path.write_text(json.dumps(metadata))
        with pytest.raises(BundleError, match="no build selected"):
            config.runtime_settings(config.build_settings(result.metadata_path), {tuple(removed.split("/"))})


def test_failed_build_keeps_log_without_a_completed_manifest(setup, monkeypatch):
    monkeypatch.setenv("MHDG_FAKE_FAIL", "parall")
    with pytest.raises(BundleError, match="command failed; see"):
        build.build_solver(setup.settings, setup.repository, requirements={(config.DEFAULT_MODEL, "serial"), (config.DEFAULT_MODEL, "mpi")})
    root = setup.repository.parent / "builds"
    assert not list(root.rglob("build_metadata.json"))
    assert "failed compiler" in next(root.rglob("*-mpi.log")).read_text()
    assert len(list(root.rglob("bin/MHDG-*"))) == 1  # Earlier serial copy is only partial evidence.


@pytest.mark.parametrize("profile", [False, True, "gamma"])
def test_check_builds_only_its_layouts_and_runs_from_manifest(setup, tmp_path, monkeypatch, profile):
    harness_root = tmp_path / "harness"
    harness_root.mkdir()
    harness = create_harness(harness_root)
    from regression_tests.bundles import create_bundle
    diverted = tmp_path / "diverted"
    create_bundle("diverted_case", harness.source, diverted, Path(__file__).resolve().parents[1] / "cases")
    settings = tmp_path / "machine.json"
    settings.write_text(json.dumps({
        "run_root": str(harness.run_root), "build_root": setup.settings["MHDG_REGRESSION_BUILD_ROOT"],
        "mpi_launcher": str(harness.mpi_launcher),
        "defaults": {"build": "nonexistent-old-build.json", "bundles": {"legacy_case": str(harness.bundle), "diverted_case": str(diverted)}},
    }))
    real_build = build.build_solver
    calls = []
    def build_once(values, repository, jobs, **kwargs):
        calls.append(kwargs["requirements"])
        return real_build(values, setup.repository, jobs, **kwargs)
    monkeypatch.setattr(build, "build_solver", build_once)
    assert cli.main([
        "check", *(["neutralgamma"] if profile == "gamma" else ["routine-extended"] if profile else ["warm", "--case", "legacy_case"]),
        "--build", "--build-jobs", "2", "--allow-candidate",
        "--run-only", "--settings", str(settings), "--run-id", "new-build",
    ]) == 0
    manifest = next((tmp_path / "builds").rglob("build_metadata.json"))
    requirements = ({("NGammaTiTeNeutralGamma", "serial")} if profile == "gamma" else
                    {(config.DEFAULT_MODEL, "serial"), (config.DEFAULT_MODEL, "mpi")} if profile else
                    {(config.DEFAULT_MODEL, "mpi")})
    assert calls == [requirements]
    assert set(json.loads(manifest.read_text())["artifacts"]) == {config.build_key(*item) for item in requirements}
    case = "diverted_case" if profile is True else "legacy_case"
    suite = "neutralgamma" if profile == "gamma" else "warm"
    summary = json.loads((harness.run_root / f"suites/{suite}/{case}/new-build/suite_summary.json").read_text())
    assert summary["execution_inputs"]["build_manifest"]["path"] == str(manifest)


FAKE_MAKE = """#!/usr/bin/env bash
set -euo pipefail
printf '%s|%s\n' "$MHDG_FAKE_ENV" "$*" >> "$MHDG_FAKE_LOG"
if [[ "$1" == "--version" ]]; then
    printf 'fixture make 1.0\n'
    exit 0
fi
if [[ "$1" == "clean" ]]; then
    rm -f MHDG-* *.o
    exit 0
fi
for arg in "$@"; do
    if [[ "$arg" == "MODE=${MHDG_FAKE_FAIL:-unset}" ]]; then
        echo 'failed compiler'
        exit 7
    fi
    case "$arg" in
      MDL=*) model=${arg#*=};;
      MHDG_BUILD_ID=*) build_id=${arg#*=};;
      MHDG_GIT_COMMIT=*) revision=${arg#*=};;
      MHDG_GIT_DIRTY=*) dirty=${arg#*=};;
    esac
done
target=${!#}
printf '#!/bin/bash\\nprintf fixture > outputs/result.h5\\n%q -B -c %q %q %q %q %q\\n' "$MHDG_FIXTURE_PYTHON" "$MHDG_FIXTURE_CODE" "$model" "$build_id" "$revision" "$dirty" > "$target"
chmod +x "$target"
touch built.o
"""


def test_profile_rejects_missing_model_before_running(tmp_path, monkeypatch, capsys):
    from regression_tests import suites
    from regression_tests.bundles import create_bundle
    from regression_tests.tests.fixtures.harness import REGRESSION_ROOT
    harness = create_harness(tmp_path)
    diverted = tmp_path / "diverted"
    create_bundle("diverted_case", harness.source, diverted, REGRESSION_ROOT / "cases")
    machine = json.loads(harness.settings.read_text())
    machine["defaults"]["bundles"]["diverted_case"] = str(diverted)
    harness.settings.write_text(json.dumps(machine))
    monkeypatch.setattr(suites, "run_profile", lambda *a, **kw: pytest.fail("must preflight all builds"))
    assert cli.main(["check", "full", "--settings", str(harness.settings), "--allow-candidate"]) == 1
    assert "NGammaTiTeNeutralGamma/serial" in capsys.readouterr().err
    assert not harness.run_root.exists()
