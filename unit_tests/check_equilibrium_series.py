"""Check physical sampling and legacy/trajectory orientation with tiny HDF5 inputs."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def inputs(directory):
    r = np.array([2., 2.1, 2.3, 2.6])
    z = np.array([-1., -.1, .5])
    times = np.array([0., .08, .2, .56, .8])
    t, y, x = np.meshgrid(times, z, r, indexing="ij")
    fields = {"psi": 2*x*x + 3*y + 7*t, "Bphi": 5+x-y+11*t,
              "Br": .2+x+2*y+3*t, "Bz": -.4+2*x-y-4*t,
              "Jtor": 1+3*x+4*y+13*t}
    with h5py.File(directory / "trajectory.h5", "w") as f:
        f["time"], f["grid/r"], f["grid/z"] = times, r, z
        f["equilibrium/psiSep"] = 40 + 17*times
        for name, values in fields.items():
            f[f"equilibrium/{name}"] = values
    with h5py.File(directory / "stationary.h5", "w") as f:
        x, y = np.meshgrid(r, z, indexing="ij")
        f["r2D"], f["z2D"] = x, y
        f["psiSep"] = 40.  # Legacy rectangular-grid files also use scalar datasets.
        names = {"psi": "flux2D", "Bphi": "Bphi2D", "Br": "Br2D", "Bz": "Bz2D", "Jtor": "Jtor"}
        for name, values in fields.items():
            f[names[name]] = values[0].T


def main():
    executable = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="mhdg-equilibrium-series-") as scratch:
        directory = Path(scratch)
        inputs(directory)
        cases = [("trajectory.h5", "trajectory", None), ("stationary", "stationary", None)]
        with h5py.File(directory / "trajectory.h5") as source, h5py.File(directory / "single.h5", "w") as f:
            f["grid/r"], f["grid/z"] = source["grid/r"][:], source["grid/z"][:]
            f["time"] = [3.2]  # A fixed input timestamp does not restrict the solver clock.
            for name in ("psi", "Bphi", "Br", "Bz", "Jtor", "psiSep"):
                f[f"equilibrium/{name}"] = source[f"equilibrium/{name}"][:1]
        cases.append(("single.h5", "stationary", None))
        shutil.copyfile(directory / "trajectory.h5", directory / "optional.h5")
        with h5py.File(directory / "optional.h5", "r+") as f:
            for name in ("Br", "Bz", "Jtor"):
                del f[f"equilibrium/{name}"]
        cases.append(("optional.h5", "optional", None))
        invalid = {
            "times": ("time", [0., .08, .08, .56, .8], "times must be finite and strictly increasing"),
            "nonfinite_times": ("time", [0., .08, np.nan, .56, .8], "times must be finite and strictly increasing"),
            "shape": ("equilibrium/Bphi", np.ones((5, 4, 3)), "shape differs from grid/time"),
            "missing_pair": ("equilibrium/Br", None, "Br and Bz must either both be present"),
            "nonfinite_frame": ("equilibrium/psi", np.full((5, 3, 4), np.nan), "Cannot read finite equilibrium frame"),
        }
        for name, (dataset, values, expected) in invalid.items():
            filename = f"{name}.h5"
            shutil.copyfile(directory / "trajectory.h5", directory / filename)
            with h5py.File(directory / filename, "r+") as f:
                del f[dataset]
                if values is not None:
                    f[dataset] = values
            cases.append((filename, "invalid", expected))
        for filename, mode, expected in cases:
            run = subprocess.run([str(executable), filename, mode], cwd=directory,
                                 capture_output=True, text=True, timeout=20)
            output = run.stdout + run.stderr
            if expected is None:
                passed = run.returncode == 0 and "EQUILIBRIUM_SERIES_PASS" in output and "HDF5-DIAG" not in output
            else:
                passed = run.returncode != 0 and expected in output and "HDF5-DIAG" not in output
            if not passed:
                raise SystemExit(f"FAIL: {filename}\n{output}")
    print(f"Equilibrium series checks: {len(cases)} PASS")


if __name__ == "__main__":
    main()
