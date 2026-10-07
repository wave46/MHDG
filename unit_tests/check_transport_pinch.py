"""Check production pinch input, runtime rows, lifecycle, and HDF5 metadata."""

from pathlib import Path
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def main():
    executable = Path(sys.argv[1]).resolve()
    valid = [
        ("omitted", "", [1]),
        ("particles", "pinch_equations=1", [1]),
        ("all", "pinch_equations=1,2,3,4", [1, 2, 3, 4]),
        ("ion_only", "pinch_equations=3", [3]),
        ("momentum_only", "pinch_equations=2", [2]),
        ("mixed", "pinch_equations=1,3,4", [1, 3, 4]),
        ("disabled", "pinch_equations=0", []),
        ("unordered", "pinch_equations=4,1,3", [1, 3, 4]),
        ("duplicates", "pinch_equations=3,3", [3]),
        ("zero_slots", "pinch_equations=0,4,0,2", [2, 4]),
    ]
    invalid = [
        ("negative", "pinch_equations=-1", "entries must be in [0,4]"),
        ("neutral", "pinch_equations=5", "entries must be in [0,4]"),
        ("overlong", "pinch_equations=1,2,3,4,0", "Invalid TRANSPORT_MODEL_1D_LST"),
    ]
    direct = [("negative", "entries must be in [0,4]"),
              ("neutral", "entries must be in [0,4]"),
              ("overlong", "accepts at most four entries")]
    with tempfile.TemporaryDirectory(prefix="mhdg-transport-pinch-") as scratch:
        for name, entries, expected in valid + invalid:
            directory = Path(scratch) / name
            directory.mkdir()
            (directory / "transport_model.nml").write_text(
                "&TRANSPORT_MODEL_1D_LST\n"
                " pinch_model=3, vpinch_const_phys=-8\n"
                f" {entries}\n/\n"
            )
            run = subprocess.run([str(executable)], cwd=directory,
                                 capture_output=True, text=True, timeout=15)
            output = run.stdout + run.stderr
            if isinstance(expected, str):
                if run.returncode == 0 or expected not in output:
                    raise SystemExit(f"FAIL: {name}\n{output}")
                continue
            if run.returncode != 0 or "transport pinch checks: PASS" not in output:
                raise SystemExit(f"FAIL: {name}\n{output}")
            with h5py.File(directory / "transport_pinch.h5") as handle:
                pinch = np.zeros((2, 6))  # HDF5 Fortran array dimensions are reversed.
                for equation in expected:
                    pinch[0, equation - 1] = -4.0  # -8 m/s * reference time/length.
                np.testing.assert_array_equal(handle["pinch"][()], pinch,
                                              err_msg=name)
                selected = expected + [0] * (4 - len(expected))
                np.testing.assert_array_equal(
                    handle["transport_1d/params/pinch_equations"][()], selected,
                    err_msg=name,
                )
        for argument, expected in direct:
            run = subprocess.run([str(executable), argument], cwd=scratch,
                                 capture_output=True, text=True, timeout=15)
            if run.returncode == 0 or expected not in run.stdout + run.stderr:
                raise SystemExit(f"FAIL: direct {argument}\n{run.stdout}{run.stderr}")
    print(f"transport pinch selection checks: {len(valid) + len(invalid) + len(direct)} PASS")


if __name__ == "__main__":
    main()
