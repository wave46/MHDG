"""Exercise the production ME restart reader and input validation on small HDF5 files."""
from pathlib import Path
import os
import re
import shlex
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def parameters(template, *, nts=100, me=True, full_paths=False):
    values = {"ME": ".true." if me else ".false.", "nts": str(nts),
              "steady": ".false.", "target_variable": "1", "printint": "0", "puff": "456",
              "field_path": "'field_0000.h5'" if full_paths else "'field'",
              "jtor_path": "'current_0000.h5'" if full_paths else "'current'",
              "zeff_path": "'controls.h5'"}
    for key, value in values.items():
        template, count = re.subn(rf"(?im)^(\s*{key}\s*=).*?$", rf"\g<1> {value}", template)
        assert count == 1, key
    if full_paths:
        template = re.sub(r"(?im)^\s*(?:(?:field|jtor)_dimensions[^\n]*|"
                          r"(?:puff|target_density|target_density_xpr|zeff)_dimension\s*=[^\n]*)\n", "", template)
    else:
        template = template.replace("&INPUT_LST", "&INPUT_LST\n"
            "    field_dimensions = 1000, 1000\n    jtor_dimensions = 1000, 1000\n"
            "    puff_dimension = 1000\n    target_density_dimension = 1000\n"
            "    target_density_xpr_dimension = 1000\n    zeff_dimension = 1000", 1)
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
        for name, size in [("u", 15), ("u_tilde", 30), ("q", 30)]:
            f.create_dataset(f"solution/{name}", data=np.ones(size))
        mesh = f.create_group("mesh")
        for name, value in {"elemType": 0, "Ndim": 2, "Nnodes": 3, "Nelems": 1,
                            "Nnodesperelem": 3, "Nnodesperface": 2, "Nextfaces": 3, "Nfaces": 3}.items():
            mesh[name] = [value]
        mesh["X"] = np.array([[2., 0.], [3., 0.], [2., 1.]]).T
        mesh["T"] = np.array([[1, 2, 3]], dtype=np.int32).T
        mesh["Tb"] = np.array([[1, 2], [2, 3], [3, 1]], dtype=np.int32).T
        mesh["boundaryFlag"] = np.array([1, 2, 3], dtype=np.int32)


def equilibrium_frames(directory, me):
    # HDF5/Python (nr,nz) is read as Fortran (nz,nr), like legacy WEST files.
    r, z = np.meshgrid([0., .1, .3, .5], [0., .25, .5], indexing="ij")
    for frame in (0, 1):
        for prefix in ("field", "current"):
            with h5py.File(directory / f"{prefix}_{frame:04d}.h5", "w") as f:
                f["r2D"], f["z2D"] = r, z
                f["time"], f["dt"] = [.08 * frame], [.08 * frame]
                if prefix == "field":
                    f["flux2D"] = r**2 + z**2
                    f["psiSep"] = [1.]
                    f["Br2D"], f["Bz2D"], f["Bphi2D"] = r, z, np.full_like(r, 1+frame)
                else:
                    f["Jtor"] = 1 + 2*r + 3*z + 10*frame
            if not me and frame == 0:
                (directory / f"{prefix}.h5").symlink_to(f"{prefix}_0000.h5")


def check_physical_restart(executable, directory, physical, *, moving):
    # Reuse the checkpoint and trajectory fixtures already exercised above.
    cases = [("old_me", 23, .04, None, None), ("prescribed", 23, .04, True, None)] if moving else [
        ("steady", 1, 1000., False, None), ("pseudo", 99, 1000., None, (False, True)),
        ("fixed", 1, .04, True, None), ("complete", 23, .08, True, None),
        ("nonfinite", 1, np.nan, True, None)]
    for name, step, seconds, marker, flags in cases:
        text = physical.replace("dt0=.02, t_final=.08", "dt0=.02, t_start=.02, t_final=.08")
        if name == "prescribed":
            text = text.replace("time_mode='fixed'", "time_mode='prescribed'")
        (directory / "param.txt").write_text(text)
        checkpoint(directory / "checkpoint.h5", step, saved_me=moving)
        with h5py.File(directory / "checkpoint.h5", "r+") as f:
            p = f["simulation_parameters"]
            p["time/Current_time"][...] = seconds/.01
            p["adimensionalization/time_scale"] = [.01]
            if marker is not None:
                p["time/Physical_time_evolution"] = [int(marker)]
                p["time/Current_time_seconds"] = [seconds]
            if flags is not None:
                p["switches/steady"], p["switches/psdtime"] = [int(flags[0])], [int(flags[1])]
        run = subprocess.run([str(executable), "physical-restart"], cwd=directory,
                             text=True, capture_output=True, timeout=20)
        output = run.stdout + run.stderr
        if name == "nonfinite":
            if run.returncode == 0 or "Checkpoint physical time must be finite" not in output:
                raise SystemExit(f"FAIL: corrupt physical checkpoint\n{output}")
            continue
        match = re.search(r"^PHYSICAL_RESTART_STATE\s+(.*)$", output, re.M)
        if run.returncode or match is None or "HDF5-DIAG" in output:
            raise SystemExit(f"FAIL: physical restart {name}\n{output}")
        values = match.group(1).split()
        initial = name in {"steady", "pseudo"}
        accepted_step = 0 if initial else step
        accepted = .02 if initial else seconds
        target = .08 if name in {"prescribed", "complete"} else accepted+.02
        expected = [accepted, accepted, target, accepted/.08 if moving else 0.]
        numeric = [float(x) for x in values[2:6]]
        if ([int(x) for x in values[:2]] != [accepted_step]*2 or
                not np.allclose(numeric, expected, rtol=0., atol=1.e-14) or
                (values[6] == "T") != (name == "complete") or
                [float(x) for x in values[7:10]] != ([456., 0., 0.] if initial else [123., 12.5, 25.]) or
                float(values[10]) < accepted_step):
            raise SystemExit(f"FAIL: physical restart {name}: {values}")
    print(f"Physical checkpoint {'moving' if moving else 'fixed'}: {len(cases)} PASS")


def main():
    executable = Path(sys.argv[1]).resolve()
    repo = Path(__file__).resolve().parents[1]
    template = (repo / "test/param_initial.txt").read_text()
    clock = subprocess.run([str(executable), "time-clock"], text=True, capture_output=True, timeout=20)
    if clock.returncode or "PHYSICAL_CLOCK_PASS" not in clock.stdout:
        raise SystemExit(f"FAIL: physical clock/history\n{clock.stdout}{clock.stderr}")
    print("Physical endpoint selection and history growth: PASS")
    check_time_input(executable, template)
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
        ("me_full_filenames", 50, True, 100, True, 50),
        ("non_me_full_filenames", 100, True, 2, False, 0),
        ("field_shape_mismatch", 0, True, 100, True, "Equilibrium dataset shape differs from grid"),
        ("current_shape_mismatch", 0, True, 100, True, "Jtor dataset shape differs from its R,Z grid"),
        ("missing_full_filename", 0, True, 100, True, "Equilibrium input file not found: missing.h5"),
        ("trajectory_me", 50, True, 100, True, 50),
        ("trajectory_fixed", 100, True, 2, False, 0),
        ("control_shape_mismatch", 0, True, 100, True, "Control dataset length differs from time"),
        ("control_rank_mismatch", 0, True, 100, True, "Control values and time must be rank one"),
        ("control_too_short", 0, True, 100, True, "Invalid control series length"),
        ("control_nonfinite", 0, True, 100, True, "Invalid control values or timestamps"),
    ]
    with tempfile.TemporaryDirectory(prefix="mhdg-me-restart-") as scratch:
        for name, step, saved_me, nts, me, expected in cases:
            directory = Path(scratch) / name
            directory.mkdir()
            text = parameters(template, nts=nts, me=me, full_paths="full_filename" in name)
            if name == "missing_full_filename":
                text = text.replace("'field_0000.h5'", "'missing.h5'")
            if name.startswith("trajectory_"):
                text = re.sub(r"(?im)^(\s*field_path\s*=).*?$", r"\g<1> 'trajectory.h5'", text)
                # Embedded current needs no separate path; exercise its namelist default.
                text = re.sub(r"(?im)^\s*jtor_path\s*=.*?\n", "", text)
            (directory / "param.txt").write_text(text)
            checkpoint(directory / "checkpoint.h5", step, saved_me)
            with h5py.File(directory / "controls.h5", "w") as f:
                times = np.array([0., .5, .56])
                f["time"] = times
                f["puff"] = 10 + 20*times
                f["zeff"] = 2 + 3*times
                f["impurity_concentration"] = .1 + .01*times
                f["target_density"] = np.ones_like(times)
                for key, value in {"x_lower": 2.1, "x_upper": 2.5, "y_lower": .1, "y_upper": .1}.items():
                    f[key] = [value]
                if name == "control_shape_mismatch":
                    del f["puff"]
                    f["puff"] = np.ones(2)
                elif name == "control_nonfinite":
                    f["puff"][1] = np.nan
                elif name == "control_rank_mismatch":
                    del f["time"]
                    f["time"] = times.reshape(1, 3)
                elif name == "control_too_short":
                    del f["time"], f["target_density"]
                    f["time"], f["target_density"] = [0.], [1.]
            mesh_file = directory / "inputs/unit_mesh_P1.msh"
            mesh_file.parent.mkdir()
            mesh_file.write_text("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n3\n"
                "1 2 0 0\n2 3 0 0\n3 2 1 0\n$EndNodes\n$Elements\n4\n"
                "1 1 2 1 1 1 2\n2 1 2 2 2 2 3\n3 1 2 3 3 3 1\n"
                "4 2 2 4 4 1 2 3\n$EndElements\n")
            equilibrium_frames(directory, me)
            if name.startswith("trajectory_"):
                frames = (0, 1) if me else (0,)
                with h5py.File(directory / "trajectory.h5", "w") as output:
                    output["time"] = [.08*frame for frame in frames]
                    with h5py.File(directory / "field_0000.h5") as field:
                        output["grid/r"], output["grid/z"] = field["r2D"][:, 0], field["z2D"][0, :]
                    for new, old in {"psi": "flux2D", "Bphi": "Bphi2D", "Br": "Br2D", "Bz": "Bz2D",
                                     "psiSep": "psiSep", "Jtor": "Jtor"}.items():
                        values = []
                        for frame in frames:
                            prefix = "current" if new == "Jtor" else "field"
                            with h5py.File(directory / f"{prefix}_{frame:04d}.h5") as source:
                                plane = source[old][...]
                                values.append(float(plane[0]) if new == "psiSep" else plane.T)
                        output[f"equilibrium/{new}"] = np.array(values)
            if name == "missing_full_filename":
                (directory / "missing_0000.h5").symlink_to("field_0000.h5")
            if name in {"field_shape_mismatch", "current_shape_mismatch"}:
                prefix, dataset = ("field", "Bphi2D") if name.startswith("field") else ("current", "Jtor")
                with h5py.File(directory / f"{prefix}_0000.h5", "r+") as f:
                    del f[dataset]
                    f[dataset] = np.ones((4, 2))
            for frame, power in [(0, 0.), (1, 7.)]:
                with h5py.File(directory / f"heat_{frame:04d}.h5", "w") as heat:
                    heat['ip'], heat['jp'] = [999], [999]  # Shapes, not legacy hints, are authoritative.
                    heat['r2D'] = [[0., 0.], [.5, .5]]
                    heat['z2D'] = [[0., .5], [0., .5]]
                    heat['heating2D_ions'] = np.full((2, 2), power)
                    heat['heating2D_electrons'] = np.full((2, 2), 2 * power)
            with h5py.File(directory / "heat_trajectory.h5", "w") as heat:
                ht, hz, hr = np.meshgrid([0., .04, .16], [0., .5], [0., .2, .5], indexing="ij")
                heat["time"], heat["grid/r"], heat["grid/z"] = [0., .04, .16], [0., .2, .5], [0., .5]
                heat["heating/ions"] = hr+hz+87.5*ht
                heat["heating/electrons"] = 2*(hr+hz+87.5*ht)
            full_names = "full_filename" in name
            arguments = [str(executable), "inputs/unit_mesh_P1.msh" if full_names else "inputs/unit_mesh_P1",
                         "checkpoint.h5" if full_names else "checkpoint"]
            run = subprocess.run(arguments, cwd=directory,
                                 text=True, capture_output=True, timeout=20)
            output = run.stdout + run.stderr
            if isinstance(expected, str):
                passed = run.returncode != 0 and expected in " ".join(output.split()) and "ME_RESTART_STATE" not in output
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
            if name in {"me_initial", "me_full_filenames", "trajectory_me", "trajectory_fixed"}:
                physical = re.sub(r"(?im)^\s*ME\s*=.*\n", "", text)
                physical = re.sub(r"(?is)&TIME_LST\b.*?\n/", "&TIME_LST\n time_mode='fixed', "
                                  "time_units='seconds', dt0=.02, t_final=.08\n/", physical)
                (directory / "param.txt").write_text(physical)
                kind = "fixed" if name == "trajectory_fixed" else "moving"
                sampling = subprocess.run([str(executable), "time-sampling", kind], cwd=directory,
                                          text=True, capture_output=True, timeout=20)
                if sampling.returncode or "PHYSICAL_SAMPLING_PASS" not in sampling.stdout:
                    raise SystemExit(f"FAIL: physical sampling {name}\n{sampling.stdout}{sampling.stderr}")
                if name == "trajectory_fixed":
                    controls = subprocess.run([str(executable), "time-controls"], cwd=directory,
                                              text=True, capture_output=True, timeout=20)
                    if (controls.returncode or "PHYSICAL_CONTROLS_PASS" not in controls.stdout or
                            "FEEDBACK_TIMING_PASS" not in controls.stdout):
                        raise SystemExit(f"FAIL: physical controls\n{controls.stdout}{controls.stderr}")
                if name in {"trajectory_fixed", "trajectory_me"}:
                    check_physical_restart(executable, directory, physical, moving=me)


    print(f"ME restart reader checks: {len(cases)} PASS")
    print("Physical sampling: numbered stem/anchor, trajectory, fixed field and changed-node heating/current PASS")
    print("Fixed-field controls: constants, signals, derivative timing, retries and checkpoint continuation PASS")
    check_nodal_input(executable, sys.argv[2] if len(sys.argv) > 2 else "serial")


def check_time_input(executable, template):
    fixed = "time_mode='fixed', time_units='seconds', dt0=.02, t_start=.12, t_final=.56"
    cases = [
        (fixed, None),
        ("time_mode='fixed', dt0=145560.40756914119, t_start=873362.4454148472, "
         "t_final=4075691.411935954", None),
        ("time_mode='prescribed', time_units='seconds', t_start=.12, tfi=.56, "
         "dt_reduce_on_failure=.true.", None),
        ("time_mode='adaptive', time_units='seconds', dt0=.02, t_start=.12, t_final=.56, "
         "dt_reduce_on_failure=.false., dt_grow_on_easy=.true.", None),
        (fixed + ", tfi=.50", "tfi and t_final specify different final times"),
        (fixed.replace(".56", "NaN"), "t_final must be finite"),
        (fixed.replace(".56", ".10"), "t_final must be greater than t_start"),
        (fixed + ", dt_reduce_on_failure=.true.", "Fixed time_mode does not support timestep adaptation"),
        (fixed + ", tis=2", "Physical time modes currently require BDF1"),
        (fixed.replace("seconds", "milliseconds"), "time_units must be normalized or seconds"),
    ]
    with tempfile.TemporaryDirectory(prefix="mhdg-time-input-") as scratch:
        directory = Path(scratch)
        for values, error in cases:
            text = parameters(template, me=False)
            text = re.sub(r"(?is)&TIME_LST\b.*?\n/", "&TIME_LST\n " + values + "\n/", text)
            (directory / "param.txt").write_text(text)
            run = subprocess.run([str(executable), "time-input"], cwd=directory,
                                 text=True, capture_output=True, timeout=20)
            output = run.stdout + run.stderr
            if error:
                if run.returncode == 0 or error not in " ".join(output.split()):
                    raise SystemExit(f"FAIL: time input {values}\n{output}")
                continue
            match = re.search(r"^TIME_INPUT_STATE\s+(.*)$", output, re.M)
            if run.returncode or not match:
                raise SystemExit(f"FAIL: time input {values}\n{output}")
            mode, units, *state = match.group(1).split()
            dt, start, final, scale, reduce, grow, steps = map(float, state)
            np.testing.assert_allclose([start*scale, final*scale], [.12, .56], rtol=1e-14)
            np.testing.assert_allclose(dt*scale, 0. if mode == "prescribed" else .02, rtol=1e-14)
            assert (reduce, grow, steps) == (float(mode == "prescribed"), float(mode == "adaptive"), 0.)
        # Explicit initialization uses physical start time and capacity, without nts.
        text = re.sub(r"(?is)&TIME_LST\b.*?\n/", "&TIME_LST\n " + fixed + "\n/", text)
        (directory / "param.txt").write_text(text)
        run = subprocess.run([str(executable), "time-input-init"], cwd=directory,
                             text=True, capture_output=True, timeout=20)
        match = re.search(r"^TIME_INIT_STATE\s+(.*)$", run.stdout, re.M)
        if run.returncode or not match:
            raise SystemExit(f"FAIL: physical clock initialization\n{run.stdout}{run.stderr}")
        np.testing.assert_allclose(list(map(float, match[1].split())), [.12, .02, 0., 16.], rtol=1e-14)
    print(f"Time input checks: {len(cases)} cases and physical initialization PASS")


def check_nodal_input(executable, build_mode):
    environment = {**os.environ, "OMP_NUM_THREADS": "4", "OPENBLAS_NUM_THREADS": "1"}
    with tempfile.TemporaryDirectory(prefix="mhdg-nodal-input-") as scratch:
        directory = Path(scratch)
        nodes = np.arange(1., 7.)
        for frame in (0, 1):
            with h5py.File(directory / f"nodal_{frame:04d}.h5", "w") as f:
                f["Br"] = nodes + .1 + 10*frame
                f["Bz"] = 2*nodes + .2 + 10*frame
                f["Bt"] = 3*nodes + .3 + 10*frame
                f["flux"] = 10 - nodes + 100*frame
                f["psiSep"] = [20. + 100*frame]
                f["time"], f["dt"] = [.08*frame], [.08*frame]
        (directory / "stationary.h5").symlink_to("nodal_0000.h5")
        commands = [(1, "stationary", "fixed"), (1, "nodal", "moving")]
        if build_mode == "parall":
            commands += [(4, "stationary.h5", "fixed"), (4, "nodal_0000.h5", "moving")]
        for ranks, filename, mode in commands:
            command = [str(executable), "nodal", filename, mode]
            if ranks > 1:
                launcher = shlex.split(os.environ.get("MPIEXEC", "mpirun.openmpi"))
                command = [*launcher, "-n", str(ranks), *command]
            run = subprocess.run(command, cwd=directory, env=environment,
                                 capture_output=True, text=True, timeout=30)
            output = run.stdout + run.stderr
            if run.returncode != 0 or f"GLOBAL_NODAL_INPUT_PASS {ranks}" not in output or "HDF5-DIAG" in output:
                raise SystemExit(f"FAIL: nodal {ranks} ranks {mode}\n{output}")
        with h5py.File(directory / "nodal_0000.h5", "r+") as f:
            del f["Bz"]
            f["Bz"] = np.ones(5)
        run = subprocess.run([str(executable), "nodal", "stationary.h5", "fixed"],
                             cwd=directory, env=environment, capture_output=True, text=True, timeout=30)
        if run.returncode == 0 or "Nodal equilibrium length differs from global mesh" not in run.stdout + run.stderr:
            raise SystemExit(f"FAIL: nodal mesh-size mismatch\n{run.stdout}{run.stderr}")
    print(f"Global nodal equilibrium checks: {len(commands)} layouts/modes and mesh-size rejection PASS")


if __name__ == "__main__":
    main()
