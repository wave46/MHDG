"""Inspect and select equilibrium trajectories for prescribed regression runs."""

from pathlib import Path

from .support import BundleError


FIELD_NAMES = ("psi", "Bphi", "Br", "Bz", "Jtor")


def trajectory_times(path):
    """Check the grid/time layout without loading discharge-sized field arrays."""
    import h5py
    import numpy as np

    try:
        with h5py.File(path) as source:
            times = np.asarray(source["time"][:], dtype=float)
            r, z = (np.asarray(source[f"grid/{axis}"][:], dtype=float) for axis in ("r", "z"))
            for name, values, minimum in (("time", times, 1), ("R", r, 2), ("Z", z, 2)):
                if (values.ndim != 1 or len(values) < minimum or not np.all(np.isfinite(values))
                        or np.any(np.diff(values) <= 0)):
                    raise BundleError(f"equilibrium {name} must contain finite, increasing coordinates")
            equilibrium = source["equilibrium"]
            if ("Br" in equilibrium) != ("Bz" in equilibrium):
                raise BundleError("equilibrium Br and Bz must both be present or absent")
            if equilibrium["psiSep"].shape != times.shape:
                raise BundleError("equilibrium psiSep shape differs from time")
            for name in FIELD_NAMES:
                if name not in equilibrium and name not in {"psi", "Bphi"}:
                    continue
                if equilibrium[name].shape != (len(times), len(z), len(r)):
                    raise BundleError(f"equilibrium {name} shape differs from time/Z/R")
    except (OSError, KeyError, ValueError, TypeError) as exc:
        raise BundleError(f"cannot inspect equilibrium trajectory {path}: {exc}") from exc
    return times


def select_trajectory_frames(source_path, destination, indices):
    """Copy selected supplied endpoints, preserving primitives and physical times."""
    import h5py

    with h5py.File(source_path) as source, h5py.File(destination, "x") as target:
        target.attrs.update(source.attrs)
        source.copy("grid", target)
        target["time"] = source["time"][indices]
        target["time"].attrs.update(source["time"].attrs)
        equilibrium = target.create_group("equilibrium")
        equilibrium.attrs.update(source["equilibrium"].attrs)
        for name in ("psiSep", "cocos", *FIELD_NAMES):
            if name not in source["equilibrium"]:
                continue
            dataset = source[f"equilibrium/{name}"]
            if name == "cocos":
                source.copy(dataset, equilibrium, name=name)
                continue
            if name == "psiSep":
                selected = equilibrium.create_dataset(name, data=dataset[indices])
            else:
                shape = (len(indices), *dataset.shape[1:])
                selected = equilibrium.create_dataset(name, shape, dtype=dataset.dtype,
                    chunks=(1, *shape[1:]), compression="gzip", compression_opts=1, shuffle=True)
                for output_frame, input_frame in enumerate(indices):
                    selected[output_frame] = dataset[input_frame]
            selected.attrs.update(dataset.attrs)
        target["provenance/source_trajectory"] = Path(source_path).name
        target["provenance/source_frames"] = indices
