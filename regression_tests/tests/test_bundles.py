from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from bundle.creation import create_bundle  # noqa: E402
from bundle.validation import validate_bundle_root  # noqa: E402
from support.errors import BundleError  # noqa: E402
from tests.fixtures.case_data import write_case_source  # noqa: E402
from tests.fixtures.harness import run_command  # noqa: E402


class BundleWorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.source = write_case_source(self.root / "source")
        self.bundle = self.root / "bundle"

    def test_public_creation_copies_links_and_produces_a_valid_bundle(self) -> None:
        equilibrium = self.root / "shared_equilibrium.h5"
        equilibrium.write_text("shared equilibrium\n", encoding="utf-8")
        (self.source / "equilibrium.h5").unlink()
        (self.source / "equilibrium.h5").symlink_to(equilibrium)

        completed = run_command(
            "bundle",
            "create",
            "--case",
            "legacy_case",
            "--source",
            str(self.source),
            "--output",
            str(self.bundle),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        copied = self.bundle / "case_data/legacy_case/equilibrium.h5"
        self.assertFalse(copied.is_symlink())
        self.assertEqual(copied.read_text(encoding="utf-8"), "shared equilibrium\n")
        summary = validate_bundle_root(self.bundle, REGRESSION_ROOT / "cases")
        self.assertEqual(summary.checked_cases, ["legacy_case"])

    def test_checksum_mismatch_is_rejected(self) -> None:
        self._create_bundle()
        manifest = self._read_manifest()
        manifest["artifacts"]["legacy_mesh"]["sha256"] = "0" * 64
        self._write_manifest(manifest)

        with self.assertRaisesRegex(BundleError, "sha256 does not match"):
            validate_bundle_root(self.bundle, REGRESSION_ROOT / "cases")

    def test_artifact_cannot_escape_the_bundle(self) -> None:
        self._create_bundle()
        outside = self.root / "outside"
        outside.mkdir()
        (outside / "mesh.msh").write_text("outside\n", encoding="utf-8")
        (self.bundle / "escape").symlink_to(outside, target_is_directory=True)
        manifest = self._read_manifest()
        manifest["artifacts"]["legacy_mesh"]["path"] = "escape/mesh.msh"
        self._write_manifest(manifest)

        with self.assertRaisesRegex(BundleError, "resolves outside"):
            validate_bundle_root(self.bundle, REGRESSION_ROOT / "cases")

    def test_missing_optional_artifact_is_reported_as_a_warning(self) -> None:
        self._create_bundle()
        manifest = self._read_manifest()
        artifact_id = "legacy_cold_fixed_time_init_parameters"
        artifact = manifest["artifacts"][artifact_id]
        artifact["optional"] = True
        self._write_manifest(manifest)
        (self.bundle / artifact["path"]).unlink()

        summary = validate_bundle_root(self.bundle, REGRESSION_ROOT / "cases")

        self.assertIn(
            f"optional artifact unavailable: {artifact_id}",
            summary.warnings,
        )

    def test_creation_failures_do_not_replace_existing_data(self) -> None:
        (self.source / "geometry.geo").unlink()
        with self.assertRaisesRegex(BundleError, "required file missing.*geometry"):
            self._create_bundle()
        self.assertFalse(self.bundle.exists())

        self.bundle.mkdir()
        marker = self.bundle / "keep.txt"
        marker.write_text("keep\n", encoding="utf-8")
        with self.assertRaisesRegex(BundleError, "output already exists"):
            self._create_bundle()
        self.assertEqual(marker.read_text(encoding="utf-8"), "keep\n")

    def _create_bundle(self) -> None:
        create_bundle(
            "legacy_case",
            self.source,
            self.bundle,
            REGRESSION_ROOT / "cases",
        )

    def _read_manifest(self) -> dict:
        return json.loads((self.bundle / "manifest.json").read_text(encoding="utf-8"))

    def _write_manifest(self, manifest: dict) -> None:
        (self.bundle / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )
