"""Exercise the production ME restart reader and input validation on small HDF5 files."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def parameters(template, *, nts=100, me=True):
    values = {"ME": ".true." if me else ".false.", "nts": str(nts),
              "steady": ".false.", "target_variable": "1", "printint": "0", "puff": "456"}
    for key, value in values.items():
        template, count = re.subn(rf"(?im)^(\s*{key}\s*=).*?$", rf"\g<1> {value}", template)
        assert count == 1, key
    return template


def checkpoint(path, step, saved_me):
    with h5py.File(path, "w") as f:
        f.create_group("magnetic")
        p = f.create_group("simulation_parameters")
        p.create_dataset("model", data=np.array([b"N-Gamma-Ti-Te-Neutral"], dtype="S100"))
        p.create_dataset("switches/ME", data=[int(saved_me)])
        p.create_dataset("time/Current_time_step_number", data=[step])
        p.create_dataset("time/Current_time", data=[step * 2.0])
        p.create_dataset("physics/puff", data=[123.0])
        p.create_dataset("physics/feedback_integral_error", data=[12.5])
        p.create_dataset("physics/feedback_previous_error", data=[25.0])
        for name, size in [("u", 5), ("u_tilde", 5), ("q", 10)]:
            f.create_dataset(f"solution/{name}", data=np.ones(size))


def main():
    executable = Path(sys.argv[1]).resolve()
    repo = Path(__file__).resolve().parents[1]
    template = (repo / "test/param_initial.txt").read_text()
    cases = [
        ("static_bootstrap", 1, False, 100, True, 0),
        ("me_initial", 0, True, 100, True, 0),
        ("me_first_step", 1, True, 100, True, 1),
        ("me_middle", 50, True, 100, True, 50),
        ("me_penultimate", 99, True, 100, True, 99),
        ("me_final", 100, True, 100, True, 100),
        ("non_me_unchanged", 100, True, 2, False, 0),
        ("beyond_final", 100, True, 99, True, "outside [0, nts"),
        ("negative_step", -1, True, 100, True, "outside [0, nts"),
        ("zero_final", 1, False, 0, True, "ME final step nts must be positive"),
    ]
    with tempfile.TemporaryDirectory(prefix="mhdg-me-restart-") as scratch:
        for name, step, saved_me, nts, me, expected in cases:
            directory = Path(scratch) / name
            directory.mkdir()
            (directory / "param.txt").write_text(parameters(template, nts=nts, me=me))
            checkpoint(directory / "checkpoint.h5", step, saved_me)
            run = subprocess.run([str(executable), "checkpoint"], cwd=directory,
                                 text=True, capture_output=True, timeout=20)
            output = run.stdout + run.stderr
            if isinstance(expected, str):
                passed = run.returncode != 0 and expected in output and "ME_RESTART_STATE" not in output
            else:
                match = re.search(r"^ME_RESTART_STATE\s+(.*)$", output, re.M)
                passed = run.returncode == 0 and match is not None and "HDF5-DIAG" not in output
                if passed:
                    values = tuple(float(x) for x in match.group(1).split())
                    restored = expected > 0
                    passed = values == (expected, expected, expected*2.0, expected*2.0,
                                        123.0 if me else 456.0,
                                        12.5 if restored else -11.0, 25.0 if restored else -12.0)
            if not passed:
                raise SystemExit(f"FAIL: {name}\n{output}")
    print(f"ME restart reader checks: {len(cases)} PASS")


if __name__ == "__main__":
    main()
