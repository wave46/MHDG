"""Check diffusion-shape inputs, runtime/lifecycle, and HDF5 metadata."""

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
        ("coordinate", "diffusion_coordinate='unknown'", "diffusion_coordinate must be"),
    ]
    invalid += [
        ("far_negative", "far_sol_min_phys=-1,0,0,0", "finite and nonnegative"),
        ("far_nan", "far_sol_min_phys=NaN,0,0,0", "finite and nonnegative"),
        ("far_zero_width", "far_sol_min_phys=1, far_sol_width=0", "far_sol_width must be positive"),
        ("far_nan_start", "far_sol_min_phys=1, far_sol_start=NaN", "coordinates and PFR scale must be finite"),
        ("far_zero_scale", "far_sol_min_phys=1, far_sol_pfr_scale=0", "far_sol_pfr_scale must be positive"),
        ("far_unresolved", "far_sol_min_phys=1, far_sol_width=1.e-30", "strictly increasing"),
        ("far_inside", "far_sol_min_phys=1, far_sol_start=.99", "outside the separatrix"),
        ("far_pfr_axis", "far_sol_min_phys=1, far_sol_pfr_scale=10", "Mirrored PFR endpoints"),
        ("far_omp_inside", "far_sol_min_phys=1, diffusion_coordinate='omp_distance', far_sol_start=-.01", "outside the separatrix"),
    ]
    far_cases = [
        ("far_rho", "", None),
        ("far_narrow", "far_sol_start=1.002, far_sol_width=.004, far_sol_pfr_scale=.5", None),
        ("far_overlap", "band_equations=1,2,3,4, band_diff_phys=.1,.1,.1,.1, band_center=1.17, band_width=.08", None),
        ("far_omp", "diffusion_coordinate='omp_distance', far_sol_start=.04, far_sol_width=.02", None),
        ("far_omp_reversed", "diffusion_coordinate='omp_distance', far_sol_start=.04, far_sol_width=.02", None),
        ("far_no_geometry", "", "require the magnetic geometry cache"),
        ("far_omp_outside", "diffusion_coordinate='omp_distance', far_sol_start=.04, far_sol_width=.02", "outside the outboard equilibrium grid"),
        ("far_omp_unprepared", "diffusion_coordinate='omp_distance', far_sol_start=.04, far_sol_width=.02", "must be prepared for the current equilibrium"),
        ("far_direct_negative", "", "finite and nonnegative"),
        ("far_direct_nan", "", "finite and nonnegative"),
    ]
    direct = [
        ("negative", "entries must be in [0,4]"),
        ("neutral", "entries must be in [0,4]"),
        ("overlong", "accepts at most four entries"),
        ("nan", "coordinates must be finite"),
        ("infinity", "finite and nonnegative"),
        ("no_geometry", "require the magnetic geometry cache"),
    ]
    omp_cases = [
        ("omp", None),
        ("omp_reversed", None),
        ("omp_outside", "outside the outboard equilibrium grid"),
        ("omp_past_axis", "outside the outboard equilibrium grid"),
        ("omp_no_geometry", "require the magnetic equilibrium and geometry cache"),
        ("omp_bad_scale", "finite positive equilibrium length scale"),
        ("omp_no_crossing", "Could not locate the OMP separatrix"),
        ("omp_unresolved", "Mapped OMP band boundaries must be strictly increasing"),
        ("omp_unprepared", "must be prepared for the current equilibrium"),
        ("omp_stale", "must be prepared for the current equilibrium"),
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
                    assert not any(key.startswith(("band_", "far_sol_")) for key in params)
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
        for argument, expected in omp_cases:
            directory = Path(scratch) / argument
            directory.mkdir()
            (directory / "transport_model.nml").write_text(
                "&TRANSPORT_MODEL_1D_LST\n diffusion_coordinate='omp_distance'\n"
                " band_equations=1,2,3,4, band_center=0, band_width=.02\n"
                " band_diff_phys=.016,.032,.048,.064\n/\n"
            )
            run = subprocess.run([str(executable), argument], cwd=directory,
                                 capture_output=True, text=True, timeout=15)
            output = run.stdout + run.stderr
            if expected:
                if run.returncode == 0 or expected not in output:
                    raise SystemExit(f"FAIL: {argument}\n{output}")
                continue
            if run.returncode != 0 or "transport OMP band checks: PASS" not in output:
                raise SystemExit(f"FAIL: {argument}\n{output}")
            assert output.count("Transport band: OMP distance [m] -> rho_pol_norm") == 4, output
            with h5py.File(directory / "transport_band_omp.h5") as handle:
                params = handle["transport_1d/params"]
                s = 2.0 + np.array([-.015, -.01, 0., .01, .015]) / .5
                expected_rho = np.sqrt((s * s + .1 * s**3) / 4.8)
                np.testing.assert_allclose(params["band_rho"][()], expected_rho, atol=1e-10)
                np.testing.assert_allclose(params["band_center"][()], 0.)
                np.testing.assert_allclose(params["band_width"][()], .02)
                np.testing.assert_allclose(params["band_transition_width"][()], .005)
                assert b"omp_distance" in params["band_coordinate"][()].tobytes()
        for argument, entries, error in far_cases:
            directory = Path(scratch) / argument
            directory.mkdir()
            (directory / "transport_model.nml").write_text(
                "&TRANSPORT_MODEL_1D_LST\n far_sol_min_phys=16,0,32,64\n"
                f" {entries}\n/\n"
            )
            run = subprocess.run([str(executable), argument], cwd=directory,
                                 capture_output=True, text=True, timeout=15)
            output = run.stdout + run.stderr
            if error:
                if run.returncode == 0 or error not in output:
                    raise SystemExit(f"FAIL: {argument}\n{output}")
                continue
            if run.returncode != 0 or "transport far-SOL checks: PASS" not in output:
                raise SystemExit(f"FAIL: {argument}\n{output}")
            with h5py.File(directory / "transport_far.h5") as handle:
                params = handle["transport_1d/params"]
                np.testing.assert_allclose(params["far_sol_min"][()], [2, 0, 4, 8])
                onset, width, scale = 1.15, .05, 1/3
                if argument == "far_narrow":
                    onset, width, scale = 1.002, .004, .5
                rho_edges = np.array([onset, onset + width])
                coordinate = b"rho_pol_norm"
                if argument.startswith("far_omp"):
                    onset, width = .04, .02
                    radius = 2.0 + np.array([onset, onset + width]) / .5
                    rho_edges = np.sqrt((radius**2 + .1*radius**3) / 4.8)
                    coordinate = b"omp_distance"
                np.testing.assert_allclose(params["far_sol_start"][()], onset)
                np.testing.assert_allclose(params["far_sol_width"][()], width)
                np.testing.assert_allclose(params["far_sol_pfr_scale"][()], scale)
                np.testing.assert_allclose(params["far_sol_rho"][()], rho_edges, atol=1e-10)
                np.testing.assert_allclose(params["far_sol_pfr_rho"][()], 1-scale*(rho_edges-1), atol=1e-10)
                assert coordinate in params["far_sol_coordinate"][()].tobytes()
    count = len(valid) + len(invalid) + len(direct) + len(omp_cases) + len(far_cases)
    print(f"transport diffusion-shape input/runtime/output checks: {count} PASS")


if __name__ == "__main__":
    main()
