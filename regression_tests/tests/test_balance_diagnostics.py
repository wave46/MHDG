from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from check_balance_diagnostics import check_suite  # noqa: E402


class BalanceDiagnosticsCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)

    def test_checks_schema_identities_and_terminal_history(self) -> None:
        summary, solution = self._write_suite()

        report = check_suite(summary)

        self.assertEqual(report["status"], "passed")
        history = report["runs"][0]["stages"][0]["terminal_history"]
        self.assertEqual(len(history), 2)
        self.assertEqual(history[-1]["time_iteration"], 2)
        self.assertEqual(history[-1]["newton_iteration"], 3)

        with h5py.File(solution, "r+") as handle:
            handle["diagnostics/particles/conservation/total"][()] = 8.5
        failed = check_suite(summary)
        self.assertEqual(failed["status"], "failed")
        self.assertTrue(
            any(
                "particles/conservation/total" in failure
                for failure in failed["failures"]
            )
        )

    def _write_suite(self) -> tuple[Path, Path]:
        run = self.root / "run"
        stage = run / "stages" / "time_init"
        stage.mkdir(parents=True)
        solution = stage / "solution.h5"
        self._write_solution(solution)
        (stage / "stdout.log").write_text(
            self._terminal_block(1, 2) + self._terminal_block(2, 3),
            encoding="utf-8",
        )
        (run / "run_metadata.json").write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "stages": [
                        {
                            "stage_id": "time_init",
                            "status": "completed",
                            "run_directory": str(stage),
                            "selected_hdf5": str(solution),
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        summary = self.root / "suite_summary.json"
        summary.write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "results": [
                        {
                            "workflow_id": "cold_fixed_balance_diagnostics",
                            "layout_id": "mpi4_omp4",
                            "run_status": "completed",
                            "run_directory": str(run),
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        return summary, solution

    def _write_solution(self, path: Path) -> None:
        particle_groups = {
            "content": (10.0, 2.0, 12.0),
            "temporal": (5.0, 1.0, 6.0),
            "volume": (2.0, -1.5, 0.5),
            "boundary_physical_inward": (-3.0, -0.5, -3.5),
            "physical_imbalance": (6.0, 3.0, 9.0),
            "hdg_tau_inward": (1.0, 0.5, 1.5),
            "conservation": (5.0, 2.5, 7.5),
        }
        components = {
            "plasma/volume": {
                "ionization": 3.0,
                "recombination": -1.0,
                "other_source": 0.0,
            },
            "plasma/boundary_inward": {
                "parallel": -4.0,
                "diffusion": 1.0,
                "pinch": 0.0,
            },
            "neutral/volume": {
                "ionization": -3.0,
                "recombination": 1.0,
                "other_source": 0.0,
                "puff_source": 1.0,
                "pump_source": -0.5,
            },
            "neutral/boundary_inward": {
                "limited_diffusion": -1.0,
                "limited_pressure": 0.5,
            },
        }
        wall_groups = {
            "plasma_particles": {
                "diffusion_inward": -1.0,
                "stabilization_inward": 1.0,
                "residual": 0.0,
            },
            "neutral": {
                "source_inward": 4.2,
                "physical_flux_inward": 4.1,
                "stabilization_inward": 0.1,
                "residual": 0.0,
                "recycled_plasma_inward": 3.8,
                "puff_source": 0.5,
                "pump_sink": 0.1,
            },
            "neutral/recycled_plasma_components": {
                "parallel_source": 4.0,
                "diffusion_source": -0.2,
                "pinch_source": 0.0,
            },
            "neutral/physical_flux_components": {
                "limited_diffusion_inward": 4.0,
                "limited_pressure_inward": 0.1,
            },
        }

        with h5py.File(path, "w") as handle:
            handle.create_dataset(
                "simulation_parameters/switches/balance_diagnostics_mode",
                data="detailed",
            )
            for group, values in particle_groups.items():
                prefix = f"diagnostics/particles/{group}"
                units = "particles" if group == "content" else "particles/s"
                handle.create_dataset(f"{prefix}/units", data=units)
                for species, value in zip(("plasma", "neutral", "total"), values):
                    handle.create_dataset(f"{prefix}/{species}", data=value)
            for group, values in components.items():
                prefix = f"diagnostics/particles/components/{group}"
                handle.create_dataset(f"{prefix}/units", data="particles/s")
                for name, value in values.items():
                    handle.create_dataset(f"{prefix}/{name}", data=value)
            handle.create_dataset(
                "diagnostics/particles/exchange/units", data="particles/s"
            )
            handle.create_dataset(
                "diagnostics/particles/exchange/charge_exchange", data=7.0
            )
            for group, values in wall_groups.items():
                prefix = f"diagnostics/wall_closure/{group}"
                handle.create_dataset(f"{prefix}/units", data="particles/s")
                for name, value in values.items():
                    handle.create_dataset(f"{prefix}/{name}", data=value)

    @staticmethod
    def _terminal_block(time_iteration: int, newton_iteration: int) -> str:
        return f"""Time iteration = {time_iteration}
NR iteration: {newton_iteration}
Particle diagnostics (detailed)
  Content [particles]
      plasma n 1.000E+01
      neutral nn 2.000E+00
      total n+nn 1.200E+01
  plasma conservation [particles/s]
      temporal 5.000E+00
      volume 2.000E+00
      boundary physical, inward -3.000E+00
      physical imbalance 6.000E+00
      HDG tau inward 1.000E+00
      discrete residual 5.000E+00
  neutral conservation [particles/s]
      temporal 1.000E+00
      volume -1.500E+00
      boundary physical, inward -5.000E-01
      physical imbalance 3.000E+00
      HDG tau inward 5.000E-01
      discrete residual 2.500E+00
  total conservation [particles/s]
      temporal 6.000E+00
      volume 5.000E-01
      boundary physical, inward -3.500E+00
      physical imbalance 9.000E+00
      HDG tau inward 1.500E+00
      discrete residual 7.500E+00
  Wall / HDG boundary-condition residuals [particles/s]
    Plasma density BC: diffusion + stabilization = 0
      residual 0.000E+00
    Neutral density BC:
      residual 0.000E+00
"""


if __name__ == "__main__":
    unittest.main()
