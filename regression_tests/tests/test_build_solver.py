from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from build_solver import BuildError, build_solver  # noqa: E402
from check_bundle import read_settings  # noqa: E402


class BuildSolverTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.repository = self.root / "repository"
        self.build_root = self.root / "builds"
        self.make_log = self.root / "make.log"
        self._create_repository()

        self.settings = self.root / "settings.env"
        self.settings.write_text(
            "MHDG_REGRESSION_SETTINGS_VERSION=1\n"
            f"MHDG_REGRESSION_BUILD_ROOT={self.build_root}\n"
            f"MHDG_ENVIRONMENT_SCRIPT={self.environment_script}\n",
            encoding="utf-8",
        )

    def test_builds_clean_variants_and_generates_provenance(self) -> None:
        result = build_solver(self.settings, self.repository, jobs=3)

        self.assertTrue(result.executables["serial"].is_file())
        self.assertTrue(result.executables["parallel"].is_file())
        self.assertTrue((result.path / "bin/positionFeketeNodesTri2D.h5").is_file())
        self.assertTrue((result.path / "logs/serial.log").is_file())
        self.assertTrue((result.path / "logs/parallel.log").is_file())

        invocations = self.make_log.read_text(encoding="utf-8").splitlines()
        self.assertEqual(
            invocations,
            [
                "loaded|clean",
                "loaded|-j3 MODE=serial COMPTYPE=opt "
                "MHDG-NGammaTiTeNeutral-serial-2D",
                "loaded|clean",
                "loaded|-j3 MODE=parall COMPTYPE=opt "
                "MHDG-NGammaTiTeNeutral-parall-2D",
                "loaded|--version",
            ],
        )

        metadata = json.loads(result.metadata_path.read_text(encoding="utf-8"))
        self.assertEqual(metadata["status"], "completed")
        self.assertFalse(metadata["repository"]["dirty"])
        self.assertEqual(metadata["profile"]["jobs"], 3)
        self.assertEqual(len(metadata["artifacts"]["serial"]["sha256"]), 64)

        settings = read_settings(result.settings_path)
        self.assertEqual(
            settings["MHDG_SERIAL_EXECUTABLE"], str(result.executables["serial"])
        )
        self.assertEqual(
            settings["MHDG_PARALLEL_EXECUTABLE"],
            str(result.executables["parallel"]),
        )
        self.assertEqual(
            settings["MHDG_BUILD_MANIFEST"], str(result.metadata_path)
        )

    def test_rejects_nonpositive_configured_job_count(self) -> None:
        self.settings.write_text(
            self.settings.read_text(encoding="utf-8")
            + "MHDG_REGRESSION_BUILD_JOBS=0\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(BuildError, "positive integer"):
            build_solver(self.settings, self.repository)

    def _create_repository(self) -> None:
        lib = self.repository / "lib"
        test = self.repository / "test"
        fake_bin = self.repository / "fake-bin"
        (lib / "Make.inc").mkdir(parents=True)
        test.mkdir()
        fake_bin.mkdir()

        (test / "positionFeketeNodesTri2D.h5").write_text(
            "synthetic nodes\n", encoding="utf-8"
        )
        fake_make = fake_bin / "make"
        fake_make.write_text(FAKE_MAKE, encoding="utf-8")
        fake_make.chmod(0o755)

        self.environment_script = lib / "Make.inc" / "environment setup.sh"
        self.environment_script.write_text(
            'test -z "$MHDG_OPTIONAL_UNSET"\n'
            f'export PATH="{fake_bin}:$PATH"\n'
            f'export MHDG_FAKE_LOG="{self.make_log}"\n'
            "export MHDG_FAKE_ENV=loaded\n",
            encoding="utf-8",
        )

        _run(["git", "init", "-q"], self.repository)
        _run(["git", "add", "."], self.repository)
        _run(
            [
                "git",
                "-c",
                "user.name=Regression Test",
                "-c",
                "user.email=regression@example.invalid",
                "commit",
                "-qm",
                "synthetic repository",
            ],
            self.repository,
        )


FAKE_MAKE = """#!/usr/bin/env bash
set -euo pipefail
printf '%s|%s\n' "$MHDG_FAKE_ENV" "$*" >> "$MHDG_FAKE_LOG"
if [[ "$1" == "--version" ]]; then
    printf 'synthetic make 1.0\n'
    exit 0
fi
if [[ "$1" == "clean" ]]; then
    rm -f MHDG-*
    exit 0
fi
target=${!#}
printf '#!/usr/bin/env bash\nexit 0\n' > "$target"
chmod +x "$target"
"""


def _run(command: list[str], cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True, capture_output=True, text=True)


if __name__ == "__main__":
    unittest.main()
