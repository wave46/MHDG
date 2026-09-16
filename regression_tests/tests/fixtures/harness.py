"""Build one reusable external bundle, executable set, and settings file."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from regression_tests.bundles import create_bundle
from regression_tests.config import settings as resolve_settings
from regression_tests.files import file_identity
from regression_tests.tests.fixtures.case_data import write_case_source, write_catalog
from regression_tests.tests.fixtures.solutions import write_solver_output


REGRESSION_ROOT = Path(__file__).resolve().parents[2]
EXIT_SOLVER = "#!/usr/bin/env bash\nexit 0\n"
MPI_LAUNCHER = """#!/usr/bin/env bash
set -euo pipefail
if [[ "$1" == '--version' ]]; then
  printf 'mpirun (Open MPI) test launcher\n'
  exit 0
fi
test "$1" = '--bind-to'
test "$2" = 'core'
test "$3" = '--map-by'
[[ "$4" == slot:PE=* ]]
shift 4
test "$1" = '-n'
if [[ -d outputs ]]; then
  printf '%s\n' "$2" > outputs/mpi_ranks.txt
fi
shift 2
exec "$@"
"""


@dataclass(frozen=True)
class HarnessFixture:
    root: Path
    source: Path
    catalog: Path
    bundle: Path
    run_root: Path
    settings: Path
    serial_executable: Path
    parallel_executable: Path
    mpi_launcher: Path
    runtime_file: Path
    environment_script: Path | None = None
    build_manifest: Path | None = None

    def install_solver(self, contents: str, *variants: str) -> None:
        selected = variants or ("serial", "parallel")
        for variant in selected:
            path = (
                self.serial_executable
                if variant == "serial"
                else self.parallel_executable
            )
            write_executable(path, contents)
        self.write_build_record()

    @property
    def values(self):
        return resolve_settings(self.settings, case="legacy_case")

    def write_build_record(self):
        def record(path):
            return {"path": path.name, **file_identity(path)}
        self.build_manifest.write_text(json.dumps({
            "schema_version": 3, "status": "completed", "build_id": "fixture-build",
            "repository": {"revision": "fixture-revision"},
            "profile": {"dimension": "2D"},
            "artifacts": {"NGammaTiTeNeutral/serial": record(self.serial_executable),
                          "NGammaTiTeNeutral/mpi": record(self.parallel_executable)},
            "runtime_files": {self.runtime_file.name: record(self.runtime_file)},
        }))

    def run_directory(
        self,
        workflow: str,
        layout: str,
        run_id: str,
    ) -> Path:
        return self.run_root / "legacy_case" / workflow / layout / run_id

    def set_bundle_class(self, bundle_class: str) -> None:
        path = self.bundle / "manifest.json"
        manifest = json.loads(path.read_text(encoding="utf-8"))
        manifest["bundle_class"] = bundle_class
        path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def create_harness(
    root: Path,
    *,
    solver: str = EXIT_SOLVER,
    include_provenance: bool = False,
) -> HarnessFixture:
    """Create the external inputs needed by prepare, run, and suite tests."""
    catalog = write_catalog(root / "catalog")
    source = write_case_source(root / "source")
    bundle = root / "bundle"
    create_bundle("legacy_case", source, bundle, catalog / "cases")

    binary_directory = root / "bin"
    binary_directory.mkdir()
    serial = write_executable(binary_directory / "serial", solver)
    parallel = write_executable(binary_directory / "parallel", solver)
    write_solver_output(path=binary_directory / "seed.h5")
    launcher = write_executable(binary_directory / "mpirun", MPI_LAUNCHER)
    runtime = binary_directory / "positionFeketeNodesTri2D.h5"
    runtime.write_text("synthetic Fekete nodes\n", encoding="utf-8")

    environment_script = None
    if include_provenance:
        environment_script = binary_directory / "environment setup.sh"
        environment_script.write_text("export MHDG_TEST_ENV=loaded\n", encoding="utf-8")
    build_manifest = binary_directory / "build_metadata.json"
    run_root = root / "runs"
    settings = root / "settings.json"
    settings.write_text(json.dumps({
        "run_root": "runs", "mpi_launcher": "bin/mpirun",
        **({"environment_script": "bin/environment setup.sh"} if environment_script else {}),
        "defaults": {"build": "bin/build_metadata.json", "bundles": {"legacy_case": "bundle"}},
    }))
    harness = HarnessFixture(root, source, catalog, bundle, run_root, settings, serial, parallel,
                             launcher, runtime, environment_script, build_manifest)
    harness.write_build_record()
    return harness


def run_command(
    *arguments: str,
    environment: dict[str, str] | None = None,
    catalog: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    command = ["-m", "regression_tests"] if catalog is None else [
        "-c", f"from pathlib import Path; from regression_tests import cli; cli.ROOT=Path({str(catalog)!r}); raise SystemExit(cli.main())"]
    return subprocess.run(
        [sys.executable, "-B", *command, *arguments],
        cwd=REGRESSION_ROOT.parent,
        check=False,
        capture_output=True,
        env={
            **os.environ,
            **(environment or {}),
        },
        text=True,
    )


def write_executable(path: Path, contents: str) -> Path:
    path.write_text(contents, encoding="utf-8")
    path.chmod(0o755)
    return path
