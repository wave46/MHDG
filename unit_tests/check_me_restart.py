"""Exercise the production ME restart reader and input validation on small HDF5 files."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def parameters(template, *, nts=100, me=True, full_paths=False):
    values = {"ME": ".true." if me else ".false.", "nts": str(nts),
              "steady": ".false.", "target_variable": "1", "printint": "0", "puff": "456",
              "field_path": "'field_0000.h5'" if full_paths else "'field'",
              "jtor_path": "'current_0000.h5'" if full_paths else "'current'"}
    for key, value in values.items():
        template, count = re.subn(rf"(?im)^(\s*{key}\s*=).*?$", rf"\g<1> {value}", template)
        assert count == 1, key
    if full_paths:
        template = re.sub(r"(?im)^\s*(?:field|jtor)_dimensions\([^\n]*\n", "", template)
    else:
        template = template.replace("&INPUT_LST", "&INPUT_LST\n"
            "    field_dimensions = 1000, 1000\n    jtor_dimensions = 1000, 1000", 1)
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
                    f["Br2D"], f["Bz2D"], f["Bphi2D"] = r, z, np.ones_like(r)
                else:
                    f["Jtor"] = 1 + 2*r + 3*z + 10*frame
            if not me and frame == 0:
                (directory / f"{prefix}.h5").symlink_to(f"{prefix}_0000.h5")


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
        ("me_full_filenames", 50, True, 100, True, 50),
        ("non_me_full_filenames", 100, True, 2, False, 0),
        ("field_shape_mismatch", 0, True, 100, True, "Equilibrium dataset shape differs from grid"),
        ("current_shape_mismatch", 0, True, 100, True, "Equilibrium dataset shape differs from grid"),
        ("missing_full_filename", 0, True, 100, True, "Equilibrium input file not found: missing.h5"),
    ]
    with tempfile.TemporaryDirectory(prefix="mhdg-me-restart-") as scratch:
        for name, step, saved_me, nts, me, expected in cases:
            directory = Path(scratch) / name
            directory.mkdir()
            text = parameters(template, nts=nts, me=me, full_paths="full_filename" in name)
            if name == "missing_full_filename":
                text = text.replace("'field_0000.h5'", "'missing.h5'")
            (directory / "param.txt").write_text(text)
            checkpoint(directory / "checkpoint.h5", step, saved_me)
            mesh_file = directory / "inputs/unit_mesh_P1.msh"
            mesh_file.parent.mkdir()
            mesh_file.write_text("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n3\n"
                "1 2 0 0\n2 3 0 0\n3 2 1 0\n$EndNodes\n$Elements\n4\n"
                "1 1 2 1 1 1 2\n2 1 2 2 2 2 3\n3 1 2 3 3 3 1\n"
                "4 2 2 4 4 1 2 3\n$EndElements\n")
            equilibrium_frames(directory, me)
            if name == "missing_full_filename":
                (directory / "missing_0000.h5").symlink_to("field_0000.h5")
            if name in {"field_shape_mismatch", "current_shape_mismatch"}:
                prefix, dataset = ("field", "Bphi2D") if name.startswith("field") else ("current", "Jtor")
                with h5py.File(directory / f"{prefix}_0000.h5", "r+") as f:
                    del f[dataset]
                    f[dataset] = np.ones((4, 2))
            for frame, power in [(0, 0.), (1, 7.)]:
                with h5py.File(directory / f"heat_{frame:04d}.h5", "w") as heat:
                    heat['ip'], heat['jp'] = [2], [2]
                    heat['r2D'] = [[0., 0.], [.5, .5]]
                    heat['z2D'] = [[0., .5], [0., .5]]
                    heat['heating2D_ions'] = np.full((2, 2), power)
                    heat['heating2D_electrons'] = np.full((2, 2), 2 * power)
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
    print(f"ME restart reader checks: {len(cases)} PASS")


if __name__ == "__main__":
    main()
