"""Exercise the real input reader and mesh transfer, without a plasma solve."""
from pathlib import Path
import math
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
    parameter = "projection_recovery_distance_ratio"
    templates = ["param_reference.txt", "param_initial.txt", "param_diffred.txt", "param_steady.txt"]
    template = (repository / "test/param_initial.txt").read_text()
    omitted = re.sub(rf"(?mi)^\s*{parameter}\s*=[^\n]*\n", "", template)
    cases = [("omitted", omitted, 0.01)]
    cases.extend(("reference_" + name, (repository / "test" / name).read_text(), 0.03)
                 for name in templates)
    for name, value, expected in [("explicit", "0.03", 0.03), ("zero", "0", 0.0),
                                  ("negative", "-0.01", None), ("nan", "NaN", None),
                                  ("infinity", "Infinity", None)]:
        source = omitted.replace("&ADAPT_LST", f"&ADAPT_LST\n  {parameter}={value}", 1)
        cases.append((name, source, expected))
    with tempfile.TemporaryDirectory(prefix="mhdg-projection-recovery-") as directory:
        scratch = Path(directory)
        for name, source, expected in cases:
            (scratch / "param.txt").write_text(source)
            run = subprocess.run([str(executable), "input"], cwd=scratch, env=environment,
                                 text=True, capture_output=True, timeout=30)
            output = run.stdout + run.stderr
            match = re.search(r"^PROJECTION_RECOVERY_RATIO\s+(\S+)", output, re.M)
            if expected is None:
                passed = run.returncode != 0 and f"{parameter} must be finite and nonnegative" in output
            else:
                passed = (run.returncode == 0 and match is not None
                          and math.isclose(float(match[1]), expected, abs_tol=1.e-14))
            if not passed:
                raise SystemExit(f"FAIL: input {name}\n{output}")

        (scratch / "positionFeketeNodesTri2D.h5").symlink_to(repository / "test/positionFeketeNodesTri2D.h5")
        launchers = [[]]
        if sys.argv[2] == "parall":
            launchers.append([*shlex.split(os.environ.get("MPIEXEC", "mpirun.openmpi")), "-n", "2"])
        projection_cases = [("default", "0.005", True), ("default", "0.02", False),
                            ("0.03", "0.02", True), ("0.03", "0.04", False),
                            ("0", "0", True), ("0", "0.002", False)]
        for launcher in launchers:
            for degree in (1, 8):
                for ratio, gap, accepted in projection_cases:
                    run = subprocess.run([*launcher, str(executable), ratio, gap, str(degree)],
                                         cwd=scratch, env=environment, text=True,
                                         capture_output=True, timeout=45)
                    output = run.stdout + run.stderr
                    if accepted:
                        passed = run.returncode == 0 and "PROJECTION_RECOVERY_PASS" in output
                    else:
                        passed = (run.returncode != 0 and "Adaptive projection has unmatched points" in output
                                  and "Projection recovery distance ratio (shortest old edge):" in output)
                    if not passed:
                        raise SystemExit(f"FAIL: projection {launcher=} P{degree} {ratio=} {gap=}\n{output}")
        print(f"Projection recovery: {len(cases)} input checks and "
              f"{len(launchers)*2*len(projection_cases)} P1/P8 transfer checks PASS")


if __name__ == "__main__":
    main()
