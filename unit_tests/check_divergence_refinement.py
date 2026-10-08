"""Exercise production recovery and history projection without a plasma solve.

Run `make check-divergence-refinement` from lib/. Parallel builds exercise one
and two MPI ranks with distinct local checkpoint values and shared retry limits.
The generic affine mesh is built in memory; no case data or Gmsh is required.
"""
from pathlib import Path
import os
import re
import shlex
import subprocess
import sys
import tempfile


def main():
    executable = Path(sys.argv[1]).resolve()
    repository = Path(__file__).resolve().parents[1]
    environment = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    with tempfile.TemporaryDirectory(prefix="mhdg-divergence-") as directory:
        scratch = Path(directory)
        (scratch / "positionFeketeNodesTri2D.h5").symlink_to(
            repository / "test/positionFeketeNodesTri2D.h5"
        )
        commands = [[str(executable), "1"]]
        if sys.argv[2] == "parall":
            launcher = shlex.split(os.environ.get("MPIEXEC", "mpirun.openmpi"))
            commands.append([*launcher, "-n", "2", str(executable), "2"])
        for command in commands:
            subprocess.run(command, cwd=scratch, env=environment, check=True, timeout=60)
        template = (repository / "test/param_initial.txt").read_text()
        # Keep the omitted-parameter check when the example declares the limit.
        template = re.sub(
            r"(?im)^\s*max_divergence_refinements\s*=.*\n", "", template
        )
        for limit in [None, 0, 1, 4, -1]:
            parameters = template
            if limit is not None:
                parameters, count = re.subn(
                    r"(?im)^(\s*div_adapt\s*=.*)$",
                    lambda match: f"{match[1]}\n    max_divergence_refinements = {limit}",
                    parameters,
                )
                assert count == 1
            (scratch / "param.txt").write_text(parameters)
            run = subprocess.run([str(executable), "1", "input"], cwd=scratch,
                                 env=environment, text=True, capture_output=True, timeout=20)
            output = run.stdout + run.stderr
            if limit == -1:
                passed = run.returncode != 0 and "max_divergence_refinements must be nonnegative" in output
            else:
                expected = 2 if limit is None else limit
                passed = run.returncode == 0 and f"DIV_REFINEMENT_LIMIT {expected}\n" in output
            if not passed:
                raise SystemExit(f"FAIL: divergence refinement input limit={limit}\n{output}")
        print("Divergence refinement input checks: default, 0/1/4 overrides, negative rejection PASS")


if __name__ == "__main__":
    main()
