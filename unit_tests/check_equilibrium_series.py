"""Check physical sampling and legacy/trajectory orientation with tiny HDF5 inputs."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.convert_equilibrium import convert_equilibria


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


def conversion_checks(directory):
    """Use the existing numerical fixtures for Python/Fortran conversion equivalence."""
    fields = [directory / "numbered" / f"field_{index:04d}.h5" for index in range(5)]
    for index, path in enumerate(fields):
        with h5py.File(path, "r+") as f:
            f["r2D"].attrs["units"] = "m"
            f["flux2D"].attrs["units"] = "Wb"
            f["flux2D"].attrs["normalization"] = "primitive flux"
            f.attrs["source_frame"] = index
    originals = [path.read_bytes() for path in fields]
    convert_equilibria(fields, directory / "converted.h5")
    convert_equilibria([directory / "stationary.h5"], directory / "converted_fixed.h5")
    with h5py.File(directory / "converted.h5") as converted, h5py.File(directory / "trajectory.h5") as expected:
        for name in ("time", "grid/r", "grid/z", "equilibrium/psi", "equilibrium/Bphi", "equilibrium/Br",
                     "equilibrium/Bz", "equilibrium/Jtor", "equilibrium/psiSep"):
            np.testing.assert_array_equal(converted[name][:], expected[name][:])
        assert converted["grid/r"].attrs["units"] == "m"
        assert converted["equilibrium/psi"].attrs["normalization"] == "primitive flux"
        assert converted["provenance/source_attributes/0004/magnetic"].attrs["source_frame"] == 4
    assert originals == [path.read_bytes() for path in fields]
    published = (directory / "converted.h5").read_bytes()
    try:
        convert_equilibria(fields, directory / "converted.h5")
    except FileExistsError:
        pass
    else:
        raise AssertionError("Converter overwrote an existing file")
    assert published == (directory / "converted.h5").read_bytes()

    # A separate same-grid Jtor series is sufficient when it is absent in the field files.
    paired = directory / "paired"
    shutil.copytree(directory / "numbered", paired)
    for index in range(5):
        with h5py.File(paired / f"field_{index:04d}.h5", "r+") as field, h5py.File(paired / f"current_{index:04d}.h5", "w") as current:
            current["r2D"], current["z2D"] = field["r2D"][:], field["z2D"][:]
            current["time"], current["Jtor"] = field["time"][...], field["Jtor"][:]
            current["Jtor"].attrs["units"] = "A/m²"
            del field["Jtor"]
    command = [sys.executable, str(Path(__file__).resolve().parents[1] / "tools/convert_equilibrium.py"),
               "--field", "paired/field_0000.h5", "--jtor", "paired/current", "--last-frame", "4",
               "--output", "converted_paired.h5"]
    run = subprocess.run(command, cwd=directory, text=True, capture_output=True, timeout=20)
    assert run.returncode == 0, run.stdout + run.stderr
    with h5py.File(directory / "converted_paired.h5") as converted:
        assert converted["equilibrium/Jtor"].attrs["units"] == "A/m²"

    for problem, expected in [("grid", "R,Z grid differs"), ("pair_grid", "Magnetic/Jtor R,Z grids differ"),
                              ("time", "timestamps differ")]:
        bad = directory / f"conversion_bad_{problem}"
        shutil.copytree(paired, bad)
        with h5py.File(bad / ("field_0004.h5" if problem in {"grid", "pair_grid"} else "current_0004.h5"), "r+") as f:
            name = "r2D" if problem in {"grid", "pair_grid"} else "time"
            f[name][...] = f[name][...] + .01
        if problem == "grid":
            with h5py.File(bad / "current_0004.h5", "r+") as f:
                f["r2D"][...] = f["r2D"][...] + .01
        output = directory / f"rejected_{problem}.h5"
        try:
            convert_equilibria([bad / path.name for path in fields], output,
                               [bad / f"current_{index:04d}.h5" for index in range(5)])
        except ValueError as error:
            assert expected in str(error), error
        else:
            raise AssertionError(f"Converter accepted inconsistent {problem}")
        assert not output.exists() and not list(directory.glob(f".{output.name}.*.tmp"))
    print("Equilibrium conversion: values, attributes, sources, CLI, rejection and publication PASS")
    return [("converted.h5", "trajectory", None, None), ("converted_fixed.h5", "stationary", None, None),
            ("converted_paired.h5", "trajectory", None, None)]


def scale_checks(directory):
    """Preserve one declared flux convention through numbered input and conversion."""
    alpha = -1 / (2*np.pi)
    scaled = directory / 'scaled'
    shutil.copytree(directory / 'numbered', scaled)
    fields = [scaled / f'field_{index:04d}.h5' for index in range(5)]
    for path in fields:
        with h5py.File(path, 'r+') as f:
            f['poloidal_scale'] = [alpha]
            f['poloidal_scale'].attrs['definition'] = 'Br=-alpha/R*dpsi/dZ; Bz=alpha/R*dpsi/dR'
    convert_equilibria(fields, directory / 'scaled.h5')
    with h5py.File(directory / 'scaled.h5') as f:
        assert f['equilibrium/poloidal_scale'].shape == ()
        assert f['equilibrium/poloidal_scale'][()] == alpha
        assert 'definition' in f['equilibrium/poloidal_scale'].attrs
    cases = [('scaled/field', 'numbered_scale', None, None), ('scaled.h5', 'scale', None, None)]
    for name, value in [('nonfinite', np.nan), ('zero', 0.), ('vector', [alpha, alpha])]:
        path = directory / f'scale_{name}.h5'
        shutil.copyfile(directory / 'scaled.h5', path)
        with h5py.File(path, 'r+') as f:
            del f['equilibrium/poloidal_scale']
            f['equilibrium/poloidal_scale'] = value
        cases.append((path.name, 'invalid', 'poloidal_scale must be a finite nonzero scalar', None))
    with h5py.File(fields[1], 'r+') as f:
        f['poloidal_scale'][...] = -alpha
    rejected = directory / 'rejected_scale.h5'
    try:
        convert_equilibria(fields, rejected)
    except ValueError as error:
        assert 'poloidal_scale differs between frames' in str(error), error
    else:
        raise AssertionError('Converter accepted changing flux conventions')
    assert not rejected.exists() and not list(directory.glob('.rejected_scale.h5.*.tmp'))
    shutil.copytree(scaled, directory / 'changed_scale')
    cases.append(('changed_scale/field', 'numbered_invalid', 'Numbered poloidal_scale differs', None))
    with h5py.File(fields[1], 'r+') as f:
        f['poloidal_scale'][...] = alpha
    return cases


def heating_inputs(directory):
    r, z = np.array([1.6, 1.9, 2.4]), np.array([-.8, .2])
    times = np.array([0., .05, .17, .45, .8])
    t, y, x = np.meshgrid(times, z, r, indexing="ij")
    with h5py.File(directory / "heating.h5", "w") as f:
        f["time"], f["grid/r"], f["grid/z"] = times, r, z
        f["heating/ions"] = 20+2*x+3*y+7*t
        f["heating/electrons"] = 30+x*x+y-4*t
    with h5py.File(directory / "heating.h5") as source, h5py.File(directory / "heating_single.h5", "w") as f:
        f["time"], f["grid/r"], f["grid/z"] = [3.2], r, z
        f["heating/ions"] = source["heating/ions"][:1]
        f["heating/electrons"] = source["heating/electrons"][:1]
    with h5py.File(directory / "heating.h5") as source, h5py.File(directory / "heating_legacy.h5", "w") as f:
        f["r2D"], f["z2D"] = np.meshgrid(r, z, indexing="ij")
        f["heating2D_ions"] = source["heating/ions"][0].T
        f["heating2D_electrons"] = source["heating/electrons"][0].T
        # Stale/omitted legacy ip,jp datasets do not determine the measured grid.
        f["ip"], f["jp"] = [999], [999]
    cases = [("heating.h5", "heating", None, None), ("heating_single.h5", "heating_fixed", None, None),
             ("heating_legacy", "heating_fixed", None, None)]
    for problem, name, values, expected in [
        ("times", "time", [0., .05, .05, .45, .8], "Heating times must be finite and strictly increasing"),
        ("shape", "heating/electrons", np.ones((5, 3, 2)), "Heating dataset shape differs from grid/time"),
        ("nonfinite", "heating/ions", np.full((5, 2, 3), np.nan), "Cannot read finite heating frame"),
    ]:
        path = directory / f"heating_bad_{problem}.h5"
        shutil.copyfile(directory / "heating.h5", path)
        with h5py.File(path, "r+") as f:
            del f[name]
            f[name] = values
        cases.append((path.name, "heating_invalid", expected, None))
    return cases


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
        cases.extend(conversion_checks(directory))
        cases.extend(scale_checks(directory))
        cases.extend(heating_inputs(directory))
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
