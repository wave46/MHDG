"""Check band input, normalization, application/lifecycle, and HDF5 metadata."""

from pathlib import Path
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def main():
    executable = Path(sys.argv[1]).resolve()
    valid = [
        ("omitted", "", []),
        ("all", "band_equations=1,2,3,4", [1, 2, 3, 4]),
        ("mixed", "band_equations=4,0,1,4", [1, 4]),
        ("ion_only", "band_equations=3", [3]),
        ("explicit_width", "band_equations=2, band_transition_width=.0005", [2]),
    ]
    invalid = [
        ("negative", "band_equations=-1", "entries must be in [0,4]"),
        ("neutral", "band_equations=5", "entries must be in [0,4]"),
        ("overlong", "band_equations=1,2,3,4,0", "Invalid TRANSPORT_MODEL_1D_LST"),
        ("missing_width", "band_equations=1, band_width=0", "band_width must be positive"),
        ("negative_target", "band_equations=1, band_diff_phys=-.01,0,0,0", "finite and nonnegative"),
        ("nan_target", "band_equations=1, band_diff_phys=NaN,0,0,0", "finite and nonnegative"),
        ("nan_center", "band_equations=1, band_center=NaN", "coordinates must be finite"),
        ("zero_shoulder", "band_equations=1, band_transition_width=0", "positive or -1"),
        ("negative_shoulder", "band_equations=1, band_transition_width=-2", "positive or -1"),
        ("across_axis", "band_equations=1, band_center=.001", "landmarks must be nonnegative"),
        ("unresolved", "band_equations=1, band_width=1.e-30", "strictly increasing"),
    ]
    direct = [
        ("negative", "entries must be in [0,4]"),
        ("neutral", "entries must be in [0,4]"),
        ("overlong", "accepts at most four entries"),
        ("nan", "coordinates must be finite"),
        ("infinity", "finite and nonnegative"),
        ("no_geometry", "requires the magnetic geometry cache"),
    ]
    with tempfile.TemporaryDirectory(prefix="mhdg-transport-band-") as scratch:
        for name, entries, expected in valid + invalid:
            directory = Path(scratch) / name
            directory.mkdir()
            (directory / "transport_model.nml").write_text(
                "&TRANSPORT_MODEL_1D_LST\n"
                " band_center=1, band_width=.004\n"
                " band_diff_phys=.016,.032,.048,.064\n"
                f" {entries}\n/\n"
            )
            run = subprocess.run([str(executable)], cwd=directory,
                                 capture_output=True, text=True, timeout=15)
            output = run.stdout + run.stderr
            if isinstance(expected, str):
                if run.returncode == 0 or expected not in output:
                    raise SystemExit(f"FAIL: {name}\n{output}")
                continue
            if run.returncode != 0 or "transport band checks: PASS" not in output:
                raise SystemExit(f"FAIL: {name}\n{output}")
            with h5py.File(directory / "transport_band.h5") as handle:
                params = handle["transport_1d/params"]
                if not expected:
                    assert not any(key.startswith("band_") for key in params)
                    continue
                np.testing.assert_array_equal(params["band_equations"][()],
                                              expected + [0] * (4 - len(expected)))
                np.testing.assert_allclose(params["band_diff"][()], [.002, .004, .006, .008])
                transition = .0005 if name == "explicit_width" else .001
                np.testing.assert_allclose(params["band_transition_width"][()], transition)
                np.testing.assert_allclose(params["band_rho"][()],
                                          [.998 - transition, .998, 1., 1.002, 1.002 + transition])
                for equation in expected:
                    np.testing.assert_allclose(handle["used"][equation - 1, 3], .002 * equation)
                    np.testing.assert_allclose(handle["used"][equation - 1, 7], .002 * equation)
        for argument, expected in direct:
            run = subprocess.run([str(executable), argument], cwd=scratch,
                                 capture_output=True, text=True, timeout=15)
            if run.returncode == 0 or expected not in run.stdout + run.stderr:
                raise SystemExit(f"FAIL: direct {argument}\n{run.stdout}{run.stderr}")
    print(f"transport band input/runtime/output checks: {len(valid) + len(invalid) + len(direct)} PASS")


if __name__ == "__main__":
    main()
