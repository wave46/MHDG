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


def numbered_inputs(directory):
    directory.mkdir()
    r = np.array([2., 2.1, 2.3, 2.6])
    z = np.array([-1., -.1, .5])
    jr, jz = np.array([2., 2.3, 2.6]), np.array([-1., .5])
    names = {"psi": "flux2D", "Bphi": "Bphi2D", "Br": "Br2D", "Bz": "Bz2D", "Jtor": "Jtor"}
    with h5py.File(directory.parent / "trajectory.h5") as source:
        for frame, time in enumerate(source["time"][:]):
            with h5py.File(directory / f"field_{frame:04d}.h5", "w") as f:
                f["r2D"], f["z2D"] = np.meshgrid(r, z, indexing="ij")
                f["time"], f["dt"] = time, 0. if frame == 0 else time-source["time"][frame-1]
                f["psiSep"] = source["equilibrium/psiSep"][frame]
                for name, legacy in names.items():
                    f[legacy] = source[f"equilibrium/{name}"][frame].T
            with h5py.File(directory / f"current_{frame:04d}.h5", "w") as f:
                x, y = np.meshgrid(jr, jz, indexing="ij")
                f["r2D"], f["z2D"], f["time"] = x, y, time
                f["Jtor"] = 1+3*x+4*y+13*time


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
        cases = [(name, mode, expected, None) for name, mode, expected in cases]
        numbered_inputs(directory / "numbered")
        cases.extend([
            ("numbered/field", "numbered", None, None),
            ("numbered/field_0000.h5", "numbered", None, None),
            ("numbered/field", "numbered", None, "numbered/current_0000.h5"),
            ("stationary.h5", "stationary", None, "numbered/current_0000.h5"),
        ])
        for variant, dataset, expected in [
            ("field_grid", "r2D", "Numbered magnetic R,Z grid differs"),
            ("current_grid", "r2D", "Numbered Jtor R,Z grid differs"),
            ("current_time", "time", "Paired equilibrium/Jtor timestamps differ"),
            ("missing_frame", None, "Cannot read numbered equilibrium metadata"),
        ]:
            shutil.copytree(directory / "numbered", directory / variant)
            if variant == "missing_frame":
                (directory / variant / "field_0004.h5").unlink()
            else:
                prefix = "field" if variant == "field_grid" else "current"
                with h5py.File(directory / variant / f"{prefix}_0001.h5", "r+") as f:
                    f[dataset][...] = f[dataset][...] + .01
            cases.append((f"{variant}/field", "numbered_invalid", expected, f"{variant}/current"))
        for filename, mode, expected, current in cases:
            command = [str(executable), filename, mode]
            if current is not None:
                command.append(current)
            run = subprocess.run(command, cwd=directory,
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
