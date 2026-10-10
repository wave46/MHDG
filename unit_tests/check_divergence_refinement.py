"""Exercise production refinement budgets, recovery and history projection without a plasma solve.

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
        # Exercise omitted defaults even when the example declares the settings.
        template = re.sub(
            r"(?im)^\s*max_(?:divergence|oscillation)_refinements\s*=.*\n", "", template
        )
        template = re.sub(r"(?im)^\s*check_neutral_positivity\s*=.*\n", "", template)
        cases = [({}, None)]
        for key, limits in [("divergence", [0, 1, 4, -1]), ("oscillation", [0, 1, 10, -1])]:
            for limit in limits:
                parameter = f"max_{key}_refinements"
                cases.append(({parameter: limit}, parameter if limit < 0 else None))
        cases.append(({"max_divergence_refinements": 2, "max_oscillation_refinements": 6}, None))
        cases.extend([({"check_neutral_positivity": ".true."}, None),
                      ({"check_neutral_positivity": ".false."}, None)])
        for overrides, invalid_parameter in cases:
            parameters = template
            if overrides:
                assignments = "\n".join(f"    {key} = {value}" for key, value in overrides.items())
                parameters, count = re.subn(
                    r"(?im)^(\s*div_adapt\s*=.*)$",
                    lambda match: f"{match[1]}\n{assignments}",
                    parameters,
                )
                assert count == 1
            (scratch / "param.txt").write_text(parameters)
            run = subprocess.run([str(executable), "1", "input"], cwd=scratch,
                                 env=environment, text=True, capture_output=True, timeout=20)
            output = run.stdout + run.stderr
            if invalid_parameter:
                passed = run.returncode != 0 and f"{invalid_parameter} must be nonnegative" in output
            else:
                expected_div = overrides.get("max_divergence_refinements", 4)
                expected_osc = overrides.get("max_oscillation_refinements", 10)
                expected_neutral = "F" if overrides.get("check_neutral_positivity") == ".false." else "T"
                passed = (run.returncode == 0
                          and f"DIV_REFINEMENT_LIMIT {expected_div}\n" in output
                          and f"OSC_REFINEMENT_LIMIT {expected_osc}\n" in output
                          and f"CHECK_NEUTRAL_POSITIVITY {expected_neutral}\n" in output)
            if not passed:
                raise SystemExit(f"FAIL: refinement input overrides={overrides}\n{output}")
        print("Refinement input checks: defaults, overrides, retry limits and neutral positivity policy PASS")


if __name__ == "__main__":
    main()
