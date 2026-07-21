"""Build one reusable external bundle, executable set, and settings file."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from bundle.creation import create_bundle
from tests.fixtures.case_data import write_case_source


REGRESSION_ROOT = Path(__file__).resolve().parents[2]
EXIT_SOLVER = "#!/usr/bin/env bash\nexit 0\n"
MPI_LAUNCHER = """#!/usr/bin/env bash
set -euo pipefail
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
    reference_writer: Callable[[Path], None] | None = None,
    include_provenance: bool = False,
) -> HarnessFixture:
    """Create the external inputs needed by prepare, run, and suite tests."""
    source = write_case_source(root / "source", reference_writer)
    bundle = root / "bundle"
    create_bundle("legacy_case", source, bundle, REGRESSION_ROOT / "cases")

    binary_directory = root / "bin"
    binary_directory.mkdir()
    serial = write_executable(binary_directory / "serial", solver)
    parallel = write_executable(binary_directory / "parallel", solver)
    launcher = write_executable(binary_directory / "mpirun", MPI_LAUNCHER)
    runtime = binary_directory / "positionFeketeNodesTri2D.h5"
    runtime.write_text("synthetic Fekete nodes\n", encoding="utf-8")

    environment_script = None
    build_manifest = None
    extra_settings = ""
    if include_provenance:
        environment_script = binary_directory / "environment setup.sh"
        environment_script.write_text(
            "export MHDG_TEST_ENV=loaded\n",
            encoding="utf-8",
        )
        build_manifest = binary_directory / "build_metadata.json"
        build_manifest.write_text("{}\n", encoding="utf-8")
        extra_settings = (
            f"MHDG_ENVIRONMENT_SCRIPT={environment_script}\n"
            f"MHDG_SOLVER_REVISION={'a' * 40}\n"
            "MHDG_BUILD_DESCRIPTION=synthetic-test-build\n"
            f"MHDG_BUILD_MANIFEST={build_manifest}\n"
        )

    run_root = root / "runs"
    settings = root / "settings.env"
    settings.write_text(
        "MHDG_REGRESSION_SETTINGS_VERSION=1\n"
        f"MHDG_REGRESSION_DATA_ROOT={bundle}\n"
        f"MHDG_REGRESSION_RUN_ROOT={run_root}\n"
        f"MHDG_SERIAL_EXECUTABLE={serial}\n"
        f"MHDG_PARALLEL_EXECUTABLE={parallel}\n"
        f"MHDG_MPI_LAUNCHER={launcher}\n"
        f"{extra_settings}",
        encoding="utf-8",
    )
    return HarnessFixture(
        root,
        source,
        bundle,
        run_root,
        settings,
        serial,
        parallel,
        launcher,
        runtime,
        environment_script,
        build_manifest,
    )


def run_command(
    *arguments: str,
    environment: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(REGRESSION_ROOT / "regression.sh"), *arguments],
        check=False,
        capture_output=True,
        env={
            **os.environ,
            "PYTHON": sys.executable,
            **(environment or {}),
        },
        text=True,
    )


def write_executable(path: Path, contents: str) -> Path:
    path.write_text(contents, encoding="utf-8")
    path.chmod(0o755)
    return path
