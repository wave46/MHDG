from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from check_bundle import BundleError, validate_bundle_root  # noqa: E402
from create_bundle import create_bundle  # noqa: E402


REQUIRED_FILES = (
    "mesh.msh",
    "geometry.geo",
    "equilibrium.h5",
    "current_density.h5",
    "param.txt",
    "transport_model.nml",
    "restart.h5",
    "reference_mpi4_omp4.h5",
)


class BundleCreationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.source = self.root / "prepared"
        self.source.mkdir()
        for filename in REQUIRED_FILES:
            (self.source / filename).write_text(
                f"synthetic {filename}\n", encoding="utf-8"
            )
        self.output = self.root / "bundle"

    def test_creates_and_validates_bundle(self) -> None:
        summary = create_bundle(
            "legacy_fixed", self.source, self.output, REGRESSION_ROOT / "cases"
        )

        self.assertEqual(summary.artifact_count, 8)
        self.assertEqual(summary.checked_cases, ["legacy_fixed"])
        self.assertTrue(
            (self.output / "case_data" / "legacy_fixed" / "geometry.geo").is_file()
        )
        validate_bundle_root(self.output, REGRESSION_ROOT / "cases")

    def test_public_command_copies_symlink_target(self) -> None:
        equilibrium = self.root / "shared_equilibrium.h5"
        equilibrium.write_text("shared equilibrium\n", encoding="utf-8")
        (self.source / "equilibrium.h5").unlink()
        (self.source / "equilibrium.h5").symlink_to(equilibrium)

        completed = subprocess.run(
            [
                str(REGRESSION_ROOT / "regression.sh"),
                "bundle",
                "create",
                "--case",
                "legacy_fixed",
                "--source",
                str(self.source),
                "--output",
                str(self.output),
            ],
            check=False,
            capture_output=True,
            env={**os.environ, "PYTHON": sys.executable},
            text=True,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("recorded 8 artifacts", completed.stdout)
        copied = self.output / "case_data" / "legacy_fixed" / "equilibrium.h5"
        self.assertFalse(copied.is_symlink())
        self.assertEqual(copied.read_text(encoding="utf-8"), "shared equilibrium\n")

    def test_available_optional_cold_file_is_recorded(self) -> None:
        filename = "param_cold_fixed_time_init.txt"
        (self.source / filename).write_text(
            "synthetic cold parameters\n", encoding="utf-8"
        )

        summary = create_bundle(
            "legacy_fixed", self.source, self.output, REGRESSION_ROOT / "cases"
        )

        manifest = json.loads(
            (self.output / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(summary.artifact_count, 9)
        self.assertIn(
            "cold_fixed_time_init_parameters",
            manifest["case_data"]["legacy_fixed"]["roles"],
        )
        self.assertTrue(
            (self.output / "case_data" / "legacy_fixed" / filename).is_file()
        )

    def test_missing_geometry_leaves_no_output(self) -> None:
        (self.source / "geometry.geo").unlink()

        with self.assertRaisesRegex(
            BundleError, "required file missing.*geometry"
        ):
            create_bundle(
                "legacy_fixed", self.source, self.output, REGRESSION_ROOT / "cases"
            )

        self.assertFalse(self.output.exists())

    def test_existing_output_is_not_replaced(self) -> None:
        self.output.mkdir()
        marker = self.output / "keep.txt"
        marker.write_text("keep\n", encoding="utf-8")

        with self.assertRaisesRegex(BundleError, "output already exists"):
            create_bundle(
                "legacy_fixed", self.source, self.output, REGRESSION_ROOT / "cases"
            )

        self.assertEqual(marker.read_text(encoding="utf-8"), "keep\n")

    def test_external_data_id_may_differ_from_case_id(self) -> None:
        contract_root = self.root / "contract"
        case_dir = contract_root / "cases"
        case_dir.mkdir(parents=True)
        shutil.copytree(REGRESSION_ROOT / "schemas", contract_root / "schemas")

        case = json.loads(
            (REGRESSION_ROOT / "cases" / "legacy_fixed.json").read_text(
                encoding="utf-8"
            )
        )
        case["external_data_id"] = "fixture_data"
        (case_dir / "legacy_fixed.json").write_text(
            json.dumps(case, indent=2) + "\n", encoding="utf-8"
        )

        summary = create_bundle(
            "legacy_fixed", self.source, self.output, case_dir
        )

        manifest = json.loads(
            (self.output / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(summary.checked_cases, ["legacy_fixed"])
        self.assertEqual(
            manifest["case_data"]["fixture_data"]["case_id"], "legacy_fixed"
        )


if __name__ == "__main__":
    unittest.main()
