from __future__ import annotations

import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

import golden_update  # noqa: E402
from support.errors import BundleError  # noqa: E402
from tests.fixtures.harness import create_harness, run_command  # noqa: E402


class GoldenUpdateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        harness_root = self.root / "harness"
        harness_root.mkdir()
        self.harness = create_harness(harness_root)
        self.harness.set_bundle_class("golden")
        self.campaigns = self.root / "golden_campaigns.json"
        self.campaigns.write_text(
            (REGRESSION_ROOT / "golden_campaigns.json").read_text(
                encoding="utf-8"
            ),
            encoding="utf-8",
        )
        schemas = self.root / "schemas"
        schemas.mkdir()
        (schemas / "golden-campaigns.schema.json").write_text(
            (
                REGRESSION_ROOT / "schemas/golden-campaigns.schema.json"
            ).read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        self.workspace = self.root / "campaign"
        self.output = self.root / "new-golden"
        self.build = SimpleNamespace(
            path=self.root / "build",
            settings_path=self.harness.settings,
            metadata_path=self.root / "build/metadata.json",
        )
        self.build.metadata_path.parent.mkdir()
        self.build.metadata_path.write_text("{}\n", encoding="utf-8")
        self.build_solver = self._patch("build_solver", return_value=self.build)
        self.run_suite = self._patch("run_suite", side_effect=self._passing_suite)
        self.compare_pairs = self._patch(
            "compare_layout_pairs",
            return_value=[{"status": "passed"}],
        )
        self.check_balance = self._patch(
            "check_balance_diagnostics",
            return_value={"status": "passed", "failures": [], "runs": []},
        )
        old_report = self.root / "old-golden.json"
        old_report.write_text("{}\n", encoding="utf-8")
        self.verify_suite = self._patch(
            "verify_suite",
            return_value=(
                old_report,
                {
                    "status": "failed",
                    "results": [{"convergence_status": "passed"}],
                },
            ),
        )
        self.promote_bundle = self._patch(
            "promote_bundle",
            side_effect=lambda *args, **kwargs: self._create_candidate(args[2]),
        )
        self.promote_mapped_bundle = self._patch(
            "promote_mapped_bundle",
            side_effect=lambda *args, **kwargs: self._create_candidate(args[3]),
        )
        self.publish_campaign_bundle = self._patch(
            "publish_campaign_bundle",
            side_effect=lambda *args, **kwargs: self._create_candidate(args[1]),
        )
        self.validate_published = self._patch("_validate_published")

    def test_update_runs_through_then_publishes_after_campaign_acceptance(self) -> None:
        suites = []

        def run_suite(*args, **kwargs):
            suite_id = args[1]
            suites.append((suite_id, kwargs["compare"], kwargs["resume"]))
            return self._passing_suite(*args, **kwargs)

        self.run_suite.side_effect = run_suite
        self.assertEqual(self._run(), 0)
        state = self._state()
        self.assertEqual(state["status"], "awaiting_acceptance")
        self.assertEqual(
            [suite for suite, _, _ in suites],
            [
                "cold_matrix",
                "warm",
                "impurity_restart_producers",
                "impurity_references",
                "neutral_sources_in_elements_restart_producer",
                "neutral_sources_in_elements_warm",
                "neutral_feature_references",
                "initialization_smoke",
                "stored_field_compatibility",
                "warm_parallelism",
                "race_matrix",
                "neutral_feature_race_matrix",
                "neutral_sources_in_elements_race_matrix",
                "balance_diagnostics_cold",
                "neutral_sources_in_elements_cold_adaptive",
                "warm",
                "impurity_mixture",
                "neutral_features_warm",
            ],
        )
        self.assertEqual(state["acceptance"]["status"], "pending")
        self.assertEqual(
            state["acceptance"]["required_stages"],
            [
                "cold_matrix",
                "warm_reference",
                "impurity_references",
                "neutral_sources_in_elements_reference",
                "neutral_feature_references",
            ],
        )
        self.assertFalse(self.output.exists())
        self.assertEqual(self._run("--accept", "campaign"), 0)
        self.assertEqual(self.promote_bundle.call_count, 1)
        self.assertEqual(self.promote_mapped_bundle.call_count, 6)
        calls = self.run_suite.call_args_list
        self.assertEqual(
            [call.args[7] for call in calls],
            ["golden", *("candidate" for _ in range(17))],
        )
        self.assertEqual(
            [Path(call.args[0]).name for call in calls[1:]],
            [
                "cold_matrix.env",
                "warm_reference.env",
                "impurity_restarts.env",
                "impurity_references.env",
                "neutral_sources_in_elements_restart.env",
                "neutral_sources_in_elements_reference.env",
                *("neutral_feature_references.env" for _ in range(11)),
            ],
        )
        self.assertEqual(self.check_balance.call_count, 3)
        self.assertEqual(
            self.promote_bundle.call_args.kwargs["matrix_warm_roles"],
            ("warm_restart",),
        )
        promotions = {
            Path(call.args[3]).name: call
            for call in self.promote_mapped_bundle.call_args_list
        }
        source_mapping = promotions[
            "neutral_sources_in_elements_reference"
        ].args[2]
        self.assertEqual(
            source_mapping,
            [
                {
                    "workflow": "warm_neutral_sources_in_elements",
                    "roles": ["warm_neutral_sources_in_elements_reference"],
                }
            ],
        )
        restart_mapping = promotions[
            "neutral_sources_in_elements_restart"
        ].args[2]
        self.assertEqual(
            restart_mapping,
            [
                {
                    "workflow": "bootstrap_neutral_sources_in_elements",
                    "roles": ["warm_neutral_sources_in_elements_restart"],
                }
            ],
        )
        impurity_restart_mappings = promotions["impurity_restarts"].args[2]
        self.assertEqual(
            [mapping["workflow"] for mapping in impurity_restart_mappings],
            [
                "bootstrap_impurity_off",
                "bootstrap_impurity_n",
                "bootstrap_impurity_nw",
            ],
        )
        neutral_mappings = promotions["neutral_feature_references"].args[2]
        self.assertEqual(
            [mapping["roles"] for mapping in neutral_mappings],
            [
                ["warm_neutral_pressure_reference"],
                ["warm_neutral_perpendicular_reference"],
                ["warm_neutral_limiter_fixed_reference"],
                ["warm_neutral_limiter_ti_reference"],
            ],
        )
        state = self._state()
        self.assertEqual(state["status"], "published")
        self.assertEqual(
            state["verification_candidate"]["root"],
            state["active_bundle"],
        )
        self.assertEqual(state["acceptance"]["status"], "accepted")
        self.assertIsNotNone(state["stages"][0].get("accepted_utc"))
        self.assertTrue(self.output.exists())
        self.assertEqual(self.publish_campaign_bundle.call_count, 1)
        published_files = dict(self.publish_campaign_bundle.call_args.args[4])
        self.assertIn("campaign.json", published_files)
        self.assertIn("build/build_metadata.json", published_files)
        self.assertIn(
            "stages/verify_impurity_mixture/suite_summary.json",
            published_files,
        )
        self.assertIn(
            "stages/verify_neutral_features/suite_summary.json",
            published_files,
        )
        self.assertIn(
            "stages/balance_diagnostics_evidence/"
            "balance_diagnostics_report.json",
            published_files,
        )
        self.assertEqual(self._run(), 0)
        self.assertEqual(self.publish_campaign_bundle.call_count, 1)
        self.validate_published.assert_called_once()
        completed = run_command("golden", "status", str(self.workspace))
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("cold_matrix: completed", completed.stdout)

    def test_partial_warm_update_warns_and_preserves_other_components(self) -> None:
        self.assertEqual(self._run("--only", "warm"), 0)
        state = self._state()
        self.assertEqual(state["stages"][0]["status"], "skipped")
        self.assertEqual(state["stages"][1]["status"], "completed")
        self.assertEqual(state["status"], "awaiting_acceptance")
        self.assertEqual(len(state["warnings"]), 1)

        self.assertEqual(
            self._run("--only", "warm", "--accept", "campaign"),
            0,
        )
        self.assertEqual(
            [call.args[1] for call in self.run_suite.call_args_list],
            [
                "warm",
                "warm",
                "impurity_mixture",
                "neutral_features_warm",
            ],
        )
        state = self._state()
        self.assertEqual(state["status"], "published")
        self.assertEqual(state["stages"][2]["status"], "skipped")
        self.assertEqual(state["stages"][3]["status"], "skipped")
        self.assertEqual(self.promote_bundle.call_count, 0)
        self.assertEqual(self.promote_mapped_bundle.call_count, 1)

    def test_interrupted_suite_uses_existing_resume(self) -> None:
        calls = []

        def run_suite(*args, **kwargs):
            calls.append(kwargs["resume"])
            if len(calls) == 1:
                raise BundleError("synthetic interruption")
            return self._passing_suite(*args, **kwargs)

        self.run_suite.side_effect = run_suite
        with redirect_stderr(StringIO()):
            self.assertEqual(self._run(), 1)
            self.assertEqual(self._run(), 0)

        self.assertEqual(calls[:2], [False, True])
        self.assertEqual(len(calls), 19)
        self.assertEqual(self._state()["stages"][0]["status"], "completed")
        self.assertEqual(self._state()["status"], "awaiting_acceptance")

    def test_failed_stage_retry_preserves_evidence_and_completed_stages(self) -> None:
        attempted_run_ids = []
        failed_once = False

        def run_suite(*args, **kwargs):
            nonlocal failed_once
            run_id = args[6]
            attempted_run_ids.append(run_id)
            if args[1] == "stored_field_compatibility" and not failed_once:
                failed_once = True
                return self._suite_summary(run_id, "failed")
            return self._passing_suite(*args, **kwargs)

        self.run_suite.side_effect = run_suite
        with redirect_stderr(StringIO()):
            self.assertEqual(self._run(), 1)

        failed_state = self._state()
        failed_stage = next(
            stage
            for stage in failed_state["stages"]
            if stage["id"] == "stored_field_compatibility"
        )
        failed_summary = failed_stage["summary"]
        self.assertEqual(failed_state["status"], "failed")
        self.assertEqual(failed_stage["status"], "failed")
        self.assertEqual(len(attempted_run_ids), 9)

        self.assertEqual(self._run("--retry-failed"), 0)
        state = self._state()
        retried_stage = next(
            stage
            for stage in state["stages"]
            if stage["id"] == "stored_field_compatibility"
        )
        self.assertEqual(state["status"], "awaiting_acceptance")
        self.assertEqual(retried_stage["status"], "completed")
        self.assertEqual(retried_stage["retry_count"], 1)
        self.assertEqual(len(retried_stage["failed_attempts"]), 1)
        self.assertEqual(
            retried_stage["failed_attempts"][0]["summary"], failed_summary
        )
        self.assertEqual(
            attempted_run_ids[9],
            "golden-test-stored_field_compatibility-retry-1",
        )
        self.assertEqual(
            attempted_run_ids.count("golden-test-cold_matrix"), 1
        )
        retry_call = next(
            call
            for call in self.run_suite.call_args_list
            if call.args[6] == "golden-test-stored_field_compatibility-retry-1"
        )
        self.assertFalse(retry_call.kwargs["resume"])

        self.assertEqual(self._run("--accept", "campaign"), 0)
        published_files = dict(self.publish_campaign_bundle.call_args.args[4])
        self.assertIn(
            "stages/stored_field_compatibility/failed_attempts/001/"
            "suite_summary.json",
            published_files,
        )

    def test_retry_requires_one_failed_stage(self) -> None:
        with redirect_stderr(StringIO()) as errors:
            self.assertEqual(self._run("--retry-failed"), 1)
        self.assertIn("no single failed stage", errors.getvalue())

    def test_reference_stage_rejects_nonconverged_output_before_promotion(
        self,
    ) -> None:
        old_report = self.root / "nonconverged-reference.json"
        old_report.write_text("{}\n", encoding="utf-8")
        self.verify_suite.return_value = (
            old_report,
            {
                "status": "failed",
                "results": [{"convergence_status": "failed"}],
            },
        )

        with redirect_stderr(StringIO()):
            self.assertEqual(self._run(), 1)

        state = self._state()
        failed = next(
            stage for stage in state["stages"] if stage["status"] == "failed"
        )
        self.assertEqual(failed["id"], "warm_reference")
        self.assertEqual(self.promote_bundle.call_count, 1)
        self.assertEqual(self.promote_mapped_bundle.call_count, 0)

    def test_balance_checker_failure_stops_campaign(self) -> None:
        self.check_balance.return_value = {
            "status": "failed",
            "failures": ["synthetic balance failure"],
            "runs": [],
        }

        with redirect_stderr(StringIO()):
            self.assertEqual(self._run(), 1)

        state = self._state()
        failed = next(
            stage for stage in state["stages"] if stage["status"] == "failed"
        )
        self.assertEqual(failed["id"], "neutral_sources_race_evidence")
        self.assertIn("balance_diagnostics_report", failed)

    def test_retry_from_rewinds_to_preceding_candidate(self) -> None:
        failed_once = False

        def run_suite(*args, **kwargs):
            nonlocal failed_once
            if args[1] == "impurity_mixture" and not failed_once:
                failed_once = True
                return self._suite_summary(args[6], "failed")
            return self._passing_suite(*args, **kwargs)

        self.run_suite.side_effect = run_suite
        with redirect_stderr(StringIO()):
            self.assertEqual(self._run(), 1)

        document = json.loads(self.campaigns.read_text(encoding="utf-8"))
        stages = document["campaigns"]["legacy_case"]["stages"]
        corrected = next(
            stage for stage in stages if stage["id"] == "impurity_references"
        )
        corrected["role_mappings"][0]["roles"].append(
            "warm_impurity_off_restart"
        )
        self.campaigns.write_text(
            json.dumps(document, indent=2) + "\n", encoding="utf-8"
        )

        calls_before_retry = self.run_suite.call_count
        self.assertEqual(self._run("--retry-from", "impurity_references"), 0)
        self.assertEqual(self.run_suite.call_count - calls_before_retry, 15)

        state = self._state()
        impurity_restarts = state["stages"][2]
        impurity_references = state["stages"][3]
        verification = next(
            stage
            for stage in state["stages"]
            if stage["id"] == "verify_impurity_mixture"
        )
        self.assertEqual(impurity_restarts.get("retry_count", 0), 0)
        self.assertEqual(impurity_references["retry_count"], 1)
        self.assertEqual(verification["retry_count"], 1)
        self.assertEqual(len(impurity_references["archived_attempts"]), 1)
        self.assertEqual(
            impurity_references["archived_attempts"][0]["declaration"]
            ["role_mappings"][0]["roles"],
            ["warm_impurity_off_reference"],
        )
        self.assertTrue(
            impurity_references["candidate"].endswith(
                "candidates/impurity_references-retry-1"
            )
        )
        self.assertEqual(
            self.run_suite.call_args_list[calls_before_retry].args[6],
            "golden-test-impurity_references-retry-1",
        )
        self.assertEqual(len(state["campaign_amendments"]), 1)
        self.assertEqual(state["status"], "awaiting_acceptance")

        self.assertEqual(self._run("--accept", "campaign"), 0)
        published_files = dict(self.publish_campaign_bundle.call_args.args[4])
        self.assertIn(
            "stages/impurity_references/archived_attempts/001/"
            "suite_summary.json",
            published_files,
        )

    def test_changed_source_settings_reject_resume(self) -> None:
        self.assertEqual(self._run(), 0)
        with self.harness.settings.open("a", encoding="utf-8") as stream:
            stream.write("# changed\n")
        with redirect_stderr(StringIO()) as errors:
            self.assertEqual(self._run(), 1)

        self.assertIn("campaign inputs changed", errors.getvalue())

    def test_unrelated_shared_catalog_change_allows_acceptance(self) -> None:
        self.assertEqual(self._run(), 0)
        run_count = self.run_suite.call_count
        document = json.loads(self.campaigns.read_text(encoding="utf-8"))
        diverted = document["campaigns"]["diverted_case"]
        diverted["stages"][0]["acceptance_required"] = False
        self.campaigns.write_text(
            json.dumps(document, indent=2) + "\n", encoding="utf-8"
        )

        self.assertEqual(self._run("--accept", "campaign"), 0)

        state = self._state()
        refresh = state["campaign_catalog_refreshes"][0]
        self.assertEqual(state["status"], "published")
        self.assertEqual(self.run_suite.call_count, run_count)
        self.assertEqual(refresh["case_id"], "legacy_case")
        self.assertNotEqual(
            refresh["previous_campaign_catalog"]["sha256"],
            refresh["campaign_catalog"]["sha256"],
        )
        self.assertEqual(
            state["inputs"]["campaign_catalog"],
            refresh["campaign_catalog"],
        )

    def test_selected_case_catalog_change_still_rejects_acceptance(self) -> None:
        self.assertEqual(self._run(), 0)
        document = json.loads(self.campaigns.read_text(encoding="utf-8"))
        legacy = document["campaigns"]["legacy_case"]
        legacy["stages"][0]["acceptance_required"] = False
        self.campaigns.write_text(
            json.dumps(document, indent=2) + "\n", encoding="utf-8"
        )

        with redirect_stderr(StringIO()) as errors:
            self.assertEqual(self._run("--accept", "campaign"), 1)

        self.assertIn("campaign inputs changed", errors.getvalue())

    def test_candidate_bootstrap_must_be_explicit(self) -> None:
        self.harness.set_bundle_class("candidate")
        self.assertEqual(self._run("--bootstrap-candidate"), 0)

        self.assertEqual(self.run_suite.call_args.args[7], "candidate")
        self.assertEqual(
            self._state()["inputs"]["source_bundle_class"],
            "candidate",
        )

    def test_candidate_bootstrap_is_rejected_by_default(self) -> None:
        self.harness.set_bundle_class("candidate")

        with redirect_stderr(StringIO()) as errors:
            self.assertEqual(self._run(), 1)

        self.assertIn("requires bundle_class=golden", errors.getvalue())

    def _arguments(self, *extra: str) -> list[str]:
        return [
            "update",
            "legacy_case",
            "--settings",
            str(self.harness.settings),
            "--run-id",
            "golden-test",
            "--workspace",
            str(self.workspace),
            "--output",
            str(self.output),
            "--bundle-version",
            "golden-test-1",
            "--campaigns",
            str(self.campaigns),
            *extra,
        ]

    def _run(self, *extra: str) -> int:
        with redirect_stdout(StringIO()):
            return golden_update.main(self._arguments(*extra))

    def _passing_suite(self, *args, **kwargs):
        return self._suite_summary(args[6], "passed")

    def _suite_summary(self, run_id: str, status: str):
        path = self.root / f"{run_id}-summary.json"
        summary = {
            "status": status,
            "execution_inputs": {
                "build_manifest": golden_update._file_record(
                    self.build.metadata_path,
                )
            },
            "results": [],
        }
        path.write_text(json.dumps(summary) + "\n", encoding="utf-8")
        return path, summary

    def _create_candidate(self, path: Path) -> None:
        path.mkdir(parents=True)
        (path / "manifest.json").write_text("{}\n", encoding="utf-8")

    def _patch(self, name: str, **kwargs):
        patcher = patch.object(golden_update, name, **kwargs)
        self.addCleanup(patcher.stop)
        return patcher.start()

    def _state(self) -> dict:
        return json.loads(
            (self.workspace / golden_update.STATE_FILE).read_text(encoding="utf-8")
        )


if __name__ == "__main__":
    unittest.main()
