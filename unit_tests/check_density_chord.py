"""Check production chord reporting, normalization, snapshots and input parsing."""
from pathlib import Path
import os
import re
import shlex
import subprocess
import sys
import tempfile

import h5py
import numpy as np

executable = Path(sys.argv[1]).resolve()
template = (Path(__file__).resolve().parents[1] / "test/param_initial.txt").read_text()
template = re.sub(r"(?im)^\s*(?:density_chord|line_integrated_density_diagnostics)\s*=.*\n", "", template)
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")


def check_report(output, integral, lcfs_average, reinitialized=False):
    lines = [line for line in output.splitlines() if line.startswith("Line-integrated density")]
    expected = [[integral, 5.994, lcfs_average]] * 2
    if reinitialized:
        expected.append([5.9952, 11.9904, 23.9808])
    assert len(lines) == len(expected), output
    pattern = (r"Line-integrated density: raw =\s*(\S+) \[m\^-2\]; "
               r"average \(in-mesh chord\) =\s*(\S+) \[m\^-3\]; "
               r"average \(inside-LCFS chord\) =\s*(\S+) \[m\^-3\]")
    for line, values in zip(lines, expected):
        match = re.fullmatch(pattern, line)
        assert match, output
        np.testing.assert_allclose([float(v) for v in match.groups()], values, rtol=5.1e-4)


with tempfile.TemporaryDirectory(prefix="mhdg-density-") as scratch:
    root = Path(scratch)
    for mode, expected in [("measure", None), ("reset", None), ("steady", None),
                           ("different", None), ("same", None), ("reversed", None),
                           ("missing", None), ("partial", None),
                           ("miss", "does not cross"), ("coverage", "mesh does not cover")]:
        result = subprocess.run([str(executable), mode], cwd=root, env=env,
                                capture_output=True, text=True, timeout=30)
        output = result.stdout + result.stderr
        if expected:
            assert result.returncode != 0 and expected in output, output
        else:
            assert result.returncode == 0 and "DENSITY_CHORD_PASS" in output, output
            assert output.count("Warning: diagnostic density_chord differs") == int(mode == "different"), output
            if mode == "different":
                assert "; each uses its own chord" in output, output
            with h5py.File(root / "density.h5") as stream:
                d = stream["density_chord"]
                expected_values = {"measurement_time_seconds": .25,
                                   "full_integral_m-2": 2.997, "mesh_length_m": .5, "lcfs_length_m": .25}
                assert set(d) == {*expected_values, "endpoints_m"}, set(d)
                for key, value in expected_values.items():
                    np.testing.assert_allclose(d[key][...], value, rtol=1e-12)
                check_report(output, expected_values["full_integral_m-2"], 11.988,
                             reinitialized=mode == "reset")
                np.testing.assert_allclose(d["endpoints_m"][...].reshape(-1),
                                           [.5, -.25, .5, .25])
            if mode == "reset":
                with h5py.File(root / "density_reset.h5") as stream:
                    assert "density_chord" not in stream, list(stream)
                with h5py.File(root / "density_reinitialized.h5") as stream:
                    d = stream["density_chord"]
                    for key, value in {"measurement_time_seconds": .75, "full_integral_m-2": 5.9952,
                                       "mesh_length_m": .5, "lcfs_length_m": .25}.items():
                        np.testing.assert_allclose(d[key][...], value, rtol=1e-12)
                    np.testing.assert_allclose(d["endpoints_m"][...].reshape(-1), [.5, -.4, .5, .4])

    cases = [({}, "DENSITY_INPUT F"),
             ({"INPUT_LST": "density_chord=.5,-.25,.5,.25",
               "UTILS_LST": "line_integrated_density_diagnostics=.true."}, "DENSITY_INPUT T"),
             ({"UTILS_LST": "line_integrated_density_diagnostics=.true."}, "requires density_chord"),
             ({"INPUT_LST": "density_chord=.5,-.25"}, "requires four finite"),
             ({"INPUT_LST": "density_chord=.5,0,.5,0"}, "endpoints must be distinct"),
             ({"INPUT_LST": "density_chord=.5,NaN,.5,1"}, "requires four finite")]
    for additions, expected in cases:
        text = re.sub(r"(?im)^(\s*printint\s*=).*?$", r"\g<1> 0", template)
        for block, values in additions.items():
            text = text.replace("&" + block, "&" + block + "\n " + values + "\n", 1)
        (root / "param.txt").write_text(text)
        result = subprocess.run([str(executable), "input"], cwd=root, env=env,
                                capture_output=True, text=True, timeout=15)
        output = result.stdout + result.stderr
        assert expected in output, output
        assert (result.returncode == 0) == expected.startswith("DENSITY_INPUT"), output
    for reference in ("param_initial.txt", "param_reference.txt", "param_steady.txt", "param_diffred.txt"):
        (root / "param.txt").write_text((Path(__file__).resolve().parents[1] / "test" / reference).read_text())
        result = subprocess.run([str(executable), "input"], cwd=root, env=env,
                                capture_output=True, text=True, timeout=15)
        assert result.returncode == 0 and "DENSITY_INPUT F" in result.stdout, result.stdout + result.stderr
    if len(sys.argv) > 2 and sys.argv[2] == "parall":
        launcher = shlex.split(os.environ.get("MPIEXEC", "mpirun.openmpi"))
        for threads in (1, 2):
            directory = root / f"mpi_threads_{threads}"
            directory.mkdir()
            result = subprocess.run(launcher + ["-np", "2", str(executable), "different"],
                                    cwd=directory, env=dict(env, OMP_NUM_THREADS=str(threads)),
                                    capture_output=True, text=True, timeout=30)
            output = result.stdout + result.stderr
            assert result.returncode == 0 and output.count("DENSITY_CHORD_PASS") == 1, output
            check_report(output, 2.997, 11.988)
            assert output.count("Warning: diagnostic density_chord differs") == 1, output
            assert "; each uses its own chord" in output, output
            with h5py.File(directory / "density.h5") as stream:
                d = stream["density_chord"]
                assert set(d) == {"measurement_time_seconds", "full_integral_m-2", "mesh_length_m",
                                  "lcfs_length_m", "endpoints_m"}, set(d)
                np.testing.assert_allclose(d["full_integral_m-2"][...], 2.997, rtol=1e-12)
                np.testing.assert_allclose(d["mesh_length_m"][...], .5, rtol=1e-12)
                np.testing.assert_allclose(d["lcfs_length_m"][...], .25, rtol=1e-12)
                np.testing.assert_allclose(d["endpoints_m"][...].reshape(-1), [.5, -.25, .5, .25])
print("Density chord: independent reporting/feedback, mismatch warning, normalizations, ownership, snapshots and inputs PASS")
