"""Check completed MPI restart startup and optionally an active one-rank ME run.

Supply an external prepared ME parameter file and a restart with an embedded
mesh. All fixture changes and decomposition output stay in temporary directories.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile
import sys

import h5py


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--parameters", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--active-steps", type=int, default=0,
                        help="also advance the supplied bootstrap with its prepared fields/controls")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    if args.active_steps < 0:
        parser.error("--active-steps must be nonnegative")
    text = args.parameters.read_text()
    launcher = shlex.split(os.environ.get("MPIEXEC", "mpirun.openmpi"))
    environment = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    for ranks in (1, 2):
        with tempfile.TemporaryDirectory(prefix="mhdg-completed-restart-") as name:
            scratch = Path(name)
            # Reject the obsolete projection mode before reading any inputs.
            # This empty directory also makes it impossible to start a solve.
            run = subprocess.run([*launcher, "-n", str(ranks), str(args.solver.resolve()),
                                  "new_mesh", "restart", "old_mesh"],
                                 cwd=scratch, env=environment, text=True,
                                 capture_output=True, timeout=30)
            output = run.stdout + run.stderr
            assert run.returncode != 0, output
            assert "mesh-to-mesh restart projection is unsupported" in output, output
            assert "STARTING" not in output and "exiting improperly" not in output, output
            print(f"Unsupported projection rejected before input loading: {ranks} MPI rank(s) PASS")
            shutil.copy2(args.checkpoint, scratch / "restart.h5")
            with h5py.File(scratch / "restart.h5", "r+") as checkpoint:
                parameters = checkpoint["simulation_parameters"]
                parameters["switches/ME"][...] = 1
                parameters["time/Current_time_step_number"][...] = 1
                # A static bootstrap has no controller history. Supply it when
                # turning this temporary fixture into a completed ME restart.
                for key in ("feedback_integral_error", "feedback_previous_error"):
                    if key not in parameters["physics"]:
                        parameters["physics"].create_dataset(key, data=[0.0])
            changes = {"ME": ".true.", "nts": "1", "readMeshFromSol": ".true.",
                       "field_path": repr(str(scratch / "absent_equilibrium")),
                       "jtor_path": repr(str(scratch / "absent_current"))}
            fixture = text
            for key, value in changes.items():
                fixture, count = re.subn(rf"(?im)^(\s*{key}\s*=).*?$",
                                        lambda match: match[1] + " " + value, fixture)
                assert count == 1, key
            (scratch / "param.txt").write_text(fixture)
            (scratch / "positionFeketeNodesTri2D.h5").symlink_to(
                repo / "test/positionFeketeNodesTri2D.h5")
            run = subprocess.run([*launcher, "-n", str(ranks), str(args.solver.resolve()),
                                  str(scratch / "mesh"), str(scratch / "restart")],
                                 cwd=scratch, env=environment, text=True,
                                 capture_output=True, timeout=60)
            output = run.stdout + run.stderr
            assert run.returncode == 0 and "ME final step already reached:" in output, output
            assert "Magnetic field loaded from file:" not in output, output
            assert "Time iteration" not in output, output
            assert "HDF5-DIAG" not in output, output
        print(f"Completed ME restart startup: {ranks} MPI rank(s) PASS; zero timesteps")

    if args.active_steps:
        sys.path.insert(0, str(repo))
        from regression_tests.compare import validate_me_completed_steps
        from regression_tests.parameters import read_selected_input_values
        import numpy as np

        with tempfile.TemporaryDirectory(prefix="mhdg-active-mpi1-") as name:
            scratch = Path(name)
            (scratch / "outputs").mkdir()
            (scratch / "positionFeketeNodesTri2D.h5").symlink_to(repo / "test/positionFeketeNodesTri2D.h5")
            fixture = text
            for key, value in {"nts": str(args.active_steps),
                               "save_folder": repr(str(scratch / "outputs") + "/")}.items():
                fixture, count = re.subn(rf"(?im)^(\s*{key}\s*=).*?$",
                                        lambda match: match[1] + " " + value, fixture)
                assert count == 1, key
            (scratch / "param.txt").write_text(fixture)
            run = subprocess.run([*launcher, "-n", "1", str(args.solver.resolve()),
                                  str(scratch / "mesh"), str(args.checkpoint.resolve().with_suffix(""))],
                                 cwd=scratch, env=environment, text=True, capture_output=True, timeout=300)
            assert run.returncode == 0, run.stdout + run.stderr
            assert "HDF5-DIAG" not in run.stdout + run.stderr, run.stdout + run.stderr
            (scratch / "stdout.log").write_text(run.stdout)
            paths = sorted((scratch / "outputs").glob("*.h5"))
            field = read_selected_input_values(args.parameters, {"field_path"})["field_path"]
            with h5py.File(field + "_0001.h5") as handle:
                dt = float(np.asarray(handle["dt"]).item())
            report = validate_me_completed_steps(scratch,
                {"moving_equilibrium_inputs": {"time_step_s": dt}, "parameter_overrides": {"nts": args.active_steps}},
                {"hdf5_outputs": [str(path) for path in paths]}, 2e-4)
            assert report["status"] == "passed" and report["complete"], report
            with h5py.File(report["steps"][0]["output"]) as initial, h5py.File(args.checkpoint) as source:
                for key in ("Nintfaces", "Nextfaces"):
                    assert np.array_equal(initial[f"mesh/{key}"], source[f"mesh/{key}"]), key
            print(f"Active ME restart: 1 MPI rank PASS; {args.active_steps} converged advances and preserved boundary counts")


if __name__ == "__main__":
    main()
