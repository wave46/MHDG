from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from check_bundle import (  # noqa: E402
    BundleError,
    load_case_definition,
    required_case_roles,
    validate_bundle,
)


class BundleValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.bundle = self.root / "bundle"
        self.bundle.mkdir()
        self.settings = self.root / "settings.env"
        self.settings.write_text(
            "MHDG_REGRESSION_SETTINGS_VERSION=1\n"
            f"MHDG_REGRESSION_DATA_ROOT={self.bundle}\n",
            encoding="utf-8",
        )

        roles = {
            "mesh": "legacy_mesh",
            "geometry": "legacy_geometry",
            "equilibrium_magnetic_field": "legacy_magnetic_field",
            "equilibrium_current_density": "legacy_current_density",
            "warm_parameters": "legacy_warm_parameters",
            "transport_configuration": "legacy_transport_configuration",
            "warm_restart": "legacy_warm_restart",
            "warm_reference": "legacy_warm_reference_mpi4_omp4",
        }
        media_types = {
            "legacy_mesh": "application/x-gmsh",
            "legacy_geometry": "text/plain",
            "legacy_magnetic_field": "application/x-hdf5",
            "legacy_current_density": "application/x-hdf5",
            "legacy_warm_parameters": "text/plain",
            "legacy_transport_configuration": "text/plain",
            "legacy_warm_restart": "application/x-hdf5",
            "legacy_warm_reference_mpi4_omp4": "application/x-hdf5",
        }

        artifacts = {}
        for artifact_id, media_type in media_types.items():
            relative_path = f"files/{artifact_id}.dat"
            artifact_path = self.bundle / relative_path
            artifact_path.parent.mkdir(parents=True, exist_ok=True)
            contents = f"synthetic {artifact_id}\n".encode()
            artifact_path.write_bytes(contents)
            artifacts[artifact_id] = {
                "path": relative_path,
                "sha256": hashlib.sha256(contents).hexdigest(),
                "size_bytes": len(contents),
                "media_type": media_type,
            }

        self.manifest = {
            "schema_version": 1,
            "bundle_id": "synthetic_bundle",
            "bundle_version": "1.0.0",
            "artifacts": artifacts,
            "case_data": {
                "legacy_fixed": {
                    "case_id": "legacy_fixed",
                    "roles": roles,
                }
            },
        }
        self._write_manifest()

    def _write_manifest(self) -> None:
        (self.bundle / "manifest.json").write_text(
            json.dumps(self.manifest, indent=2) + "\n", encoding="utf-8"
        )

    def test_valid_bundle(self) -> None:
        summary = validate_bundle(self.settings, REGRESSION_ROOT / "cases")
        self.assertEqual(summary.artifact_count, 8)
        self.assertEqual(summary.verified_artifact_count, 8)
        self.assertEqual(summary.checked_cases, ["legacy_fixed"])

    def test_cold_roles_are_required_only_by_cold_workflow(self) -> None:
        case = load_case_definition("legacy_fixed", REGRESSION_ROOT / "cases")
        role = "cold_fixed_time_init_parameters"

        self.assertNotIn(role, required_case_roles(case))
        self.assertIn(role, required_case_roles(case, "cold_fixed"))
        self.assertNotIn(
            "adaptive_initial_mesh", required_case_roles(case, "cold_fixed")
        )
        self.assertIn(
            "adaptive_initial_mesh", required_case_roles(case, "cold_adaptive")
        )

    def test_public_check_data_command(self) -> None:
        completed = subprocess.run(
            [
                str(REGRESSION_ROOT / "regression.sh"),
                "--settings",
                str(self.settings),
                "check-data",
            ],
            check=False,
            capture_output=True,
            env={**os.environ, "PYTHON": sys.executable},
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("bundle valid: synthetic_bundle version 1.0.0", completed.stdout)

    def test_checksum_mismatch_fails(self) -> None:
        self.manifest["artifacts"]["legacy_mesh"]["sha256"] = "0" * 64
        self._write_manifest()
        with self.assertRaisesRegex(BundleError, "sha256 does not match"):
            validate_bundle(self.settings, REGRESSION_ROOT / "cases")

    def test_size_mismatch_fails(self) -> None:
        self.manifest["artifacts"]["legacy_mesh"]["size_bytes"] += 1
        self._write_manifest()
        with self.assertRaisesRegex(BundleError, "size_bytes"):
            validate_bundle(self.settings, REGRESSION_ROOT / "cases")

    def test_manifest_schema_error_fails(self) -> None:
        del self.manifest["bundle_version"]
        self._write_manifest()
        with self.assertRaisesRegex(BundleError, "bundle_version.*required"):
            validate_bundle(self.settings, REGRESSION_ROOT / "cases")

    def test_parent_path_fails(self) -> None:
        self.manifest["artifacts"]["legacy_mesh"]["path"] = "../mesh.msh"
        self._write_manifest()
        with self.assertRaisesRegex(BundleError, "legacy_mesh.path"):
            validate_bundle(self.settings, REGRESSION_ROOT / "cases")

    def test_symlink_escape_fails(self) -> None:
        outside = self.root / "outside"
        outside.mkdir()
        (outside / "mesh.msh").write_text("outside\n", encoding="utf-8")
        (self.bundle / "escape").symlink_to(outside, target_is_directory=True)
        artifact = self.manifest["artifacts"]["legacy_mesh"]
        artifact["path"] = "escape/mesh.msh"
        self._write_manifest()
        with self.assertRaisesRegex(BundleError, "resolves outside"):
            validate_bundle(self.settings, REGRESSION_ROOT / "cases")

    def test_missing_optional_artifact_warns(self) -> None:
        self.manifest["artifacts"]["supplemental_notes"] = {
            "path": "files/notes.txt",
            "sha256": "0" * 64,
            "size_bytes": 0,
            "media_type": "text/plain",
            "optional": True,
        }
        self.manifest["case_data"]["legacy_fixed"]["roles"][
            "supplemental_notes"
        ] = "supplemental_notes"
        self._write_manifest()

        summary = validate_bundle(self.settings, REGRESSION_ROOT / "cases")

        self.assertEqual(summary.artifact_count, 9)
        self.assertEqual(summary.verified_artifact_count, 8)
        self.assertIn(
            "optional artifact unavailable: supplemental_notes", summary.warnings
        )

    def test_missing_required_case_file_role_fails(self) -> None:
        del self.manifest["case_data"]["legacy_fixed"]["roles"]["geometry"]
        self._write_manifest()
        with self.assertRaisesRegex(BundleError, "missing required artifact roles"):
            validate_bundle(self.settings, REGRESSION_ROOT / "cases")

    def test_case_id_mismatch_fails(self) -> None:
        case_data = self.manifest["case_data"]["legacy_fixed"]
        case_data["case_id"] = "historical_feature"
        self._write_manifest()
        with self.assertRaisesRegex(BundleError, "expected legacy_fixed"):
            validate_bundle(self.settings, REGRESSION_ROOT / "cases")


if __name__ == "__main__":
    unittest.main()
