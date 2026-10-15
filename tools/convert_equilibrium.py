"""Convert legacy rectangular-grid field/current files to one trajectory HDF5."""
from __future__ import annotations

import argparse
from collections.abc import Sequence
from dataclasses import dataclass
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np


FIELD_NAMES = {"psi": "flux2D", "Bphi": "Bphi2D", "Br": "Br2D", "Bz": "Bz2D", "Jtor": "Jtor"}


@dataclass
class _Frame:
    r: np.ndarray
    z: np.ndarray
    time: float
    psi_sep: float
    poloidal_scale: float | None
    fields: dict[str, np.ndarray]
    attributes: dict[str, dict]
    source_attributes: dict[str, dict]


def convert_equilibria(
    field_files: Sequence[str | Path], output: str | Path, current_files: Sequence[str | Path] | None = None,
) -> Path:
    """Convert ordered files without retiming, regridding or changing normalization.

    Multiple frames require timestamps, including in paired current files. A
    stationary file without a timestamp uses zero. References may be absent;
    Br/Bz must be paired. Separate current files override embedded legacy Jtor.
    An explicit poloidal_scale is preserved and must stay constant across frames.
    Output is published only after success and must not already exist.
    """
    fields = [Path(path) for path in field_files]
    currents = [Path(path) for path in current_files] if current_files is not None else None
    if not fields:
        raise ValueError("At least one magnetic frame is required")
    if currents is not None and len(currents) != len(fields):
        raise ValueError("Supply one current file for each magnetic frame")
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    with tempfile.NamedTemporaryFile(dir=output.parent, prefix=f".{output.name}.", suffix=".tmp", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        _write_trajectory(fields, currents, temporary)
        # Publish without replacing an existing file, including a concurrent creation.
        os.link(temporary, output)
    finally:
        temporary.unlink()
    return output


def _write_trajectory(fields, currents, output):
    previous_time = None
    with h5py.File(output, "w") as destination:
        for index, field_file in enumerate(fields):
            current_file = currents[index] if currents is not None else None
            frame = _read_frame(field_file, current_file, require_time=len(fields) > 1)
            if index == 0:
                _create_layout(destination, frame, fields, currents)
            else:
                _validate_frame(frame, destination, previous_time, field_file)
            destination["time"][index] = frame.time
            destination["equilibrium/psiSep"][index] = frame.psi_sep
            for name, values in frame.fields.items():
                destination[f"equilibrium/{name}"][index] = values
            for kind, attributes in frame.source_attributes.items():
                if attributes:
                    destination.require_group(f"provenance/source_attributes/{index:04d}/{kind}").attrs.update(attributes)
            previous_time = frame.time
            del frame


def _read_frame(field_file, current_file, *, require_time):
    with h5py.File(field_file, "r") as field:
        r, z = _read_grid(field)
        timestamp = _read_scalar(field, "time", required=require_time)
        psi_sep = _read_scalar(field, "psiSep")
        poloidal_scale = _read_scalar(field, "poloidal_scale", required=False)
        if poloidal_scale == 0:
            raise ValueError(f"poloidal_scale must be nonzero: {field_file}")
        attributes = {"grid/r": dict(field["r2D"].attrs), "grid/z": dict(field["z2D"].attrs),
                      "equilibrium/psiSep": dict(field["psiSep"].attrs),
                      "time": dict(field["time"].attrs) if "time" in field else {}}
        if poloidal_scale is not None:
            attributes["equilibrium/poloidal_scale"] = dict(field["poloidal_scale"].attrs)
        values = {}
        for name, legacy in FIELD_NAMES.items():
            if name == "Jtor" and current_file is not None:
                continue
            if legacy not in field and name not in {"psi", "Bphi"}:
                continue
            values[name] = _read_matrix(field, legacy, (len(r), len(z))).T
            attributes[f"equilibrium/{name}"] = dict(field[legacy].attrs)
        if ("Br" in values) != ("Bz" in values):
            raise ValueError(f"Br and Bz must both be present or absent: {field_file}")
        source_attributes = {"magnetic": dict(field.attrs)}
    if current_file is not None:
        with h5py.File(current_file, "r") as current:
            jr, jz = _read_grid(current)
            if not np.array_equal(r, jr) or not np.array_equal(z, jz):
                raise ValueError(f"Magnetic/Jtor R,Z grids differ; no regridding: {current_file}")
            current_time = _read_scalar(current, "time", required=require_time)
            if timestamp is not None and current_time is not None:
                tolerance = 64*np.finfo(float).eps*max(1., abs(timestamp), abs(current_time))
                if abs(timestamp-current_time) > tolerance:
                    raise ValueError(f"Paired magnetic/Jtor timestamps differ: {current_file}")
            values["Jtor"] = _read_matrix(current, "Jtor", (len(r), len(z))).T
            attributes["equilibrium/Jtor"] = dict(current["Jtor"].attrs)
            source_attributes["current"] = dict(current.attrs)
    return _Frame(r, z, timestamp if timestamp is not None else 0., psi_sep, poloidal_scale,
                  values, attributes, source_attributes)


def _read_grid(source):
    r2d = _read_matrix(source, "r2D")
    z2d = _read_matrix(source, "z2D", r2d.shape)
    r, z = r2d[:, 0].copy(), z2d[0, :].copy()
    if len(r) < 2 or len(z) < 2 or np.any(np.diff(r) <= 0) or np.any(np.diff(z) <= 0):
        raise ValueError(f"R,Z axes must contain at least two increasing points: {source.filename}")
    if not np.all(r2d == r[:, None]) or not np.all(z2d == z[None, :]):
        raise ValueError(f"R,Z must form a rectangular grid: {source.filename}")
    return r, z


def _read_matrix(source, name, shape=None):
    if name not in source:
        raise ValueError(f"Missing {name}: {source.filename}")
    dataset = source[name]
    if (dataset.ndim != 2 or any(size < 1 for size in dataset.shape)
            or (shape is not None and dataset.shape != shape)):
        raise ValueError(f"Invalid {name} grid shape: {source.filename}")
    if dataset.dtype.kind not in "fiu":
        raise ValueError(f"{name} must be numeric: {source.filename}")
    values = np.asarray(dataset[...], dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"Nonfinite {name}: {source.filename}")
    return values


def _read_scalar(source, name, *, required=True):
    if name not in source:
        if required:
            raise ValueError(f"Missing {name}: {source.filename}")
        return None
    dataset = source[name]
    if dataset.size != 1 or dataset.dtype.kind not in "fiu":
        raise ValueError(f"{name} must be a numeric scalar or an array containing exactly one value: {source.filename}")
    value = float(np.asarray(dataset[...]).item())
    if not np.isfinite(value):
        raise ValueError(f"Nonfinite {name}: {source.filename}")
    return value


def _create_layout(destination, frame, fields, currents):
    destination["grid/r"], destination["grid/z"] = frame.r, frame.z
    destination.create_dataset("time", (len(fields),), dtype="f8")
    destination.create_dataset("equilibrium/psiSep", (len(fields),), dtype="f8")
    if frame.poloidal_scale is not None:
        destination["equilibrium/poloidal_scale"] = frame.poloidal_scale
    shape = (len(fields), len(frame.z), len(frame.r))
    for name in frame.fields:
        destination.create_dataset(f"equilibrium/{name}", shape, dtype="f8", chunks=(1, *shape[1:]),
                                   compression="gzip", compression_opts=1, shuffle=True)
    for name, attributes in frame.attributes.items():
        destination[name].attrs.update(attributes)
    strings = h5py.string_dtype("utf-8")
    destination.create_dataset("provenance/magnetic_files", data=[path.name for path in fields], dtype=strings)
    if currents is not None:
        destination.create_dataset("provenance/current_files", data=[path.name for path in currents], dtype=strings)


def _validate_frame(frame, destination, previous_time, field_file):
    if frame.time <= previous_time:
        raise ValueError(f"Timestamps must be strictly increasing: {field_file}")
    if not np.array_equal(frame.r, destination["grid/r"][:]) or not np.array_equal(frame.z, destination["grid/z"][:]):
        raise ValueError(f"R,Z grid differs from the first frame; no regridding: {field_file}")
    if set(frame.fields) != set(destination["equilibrium"]) - {"psiSep", "poloidal_scale"}:
        raise ValueError(f"Optional reference fields differ between frames: {field_file}")
    if frame.poloidal_scale != _read_scalar(destination, "equilibrium/poloidal_scale", required=False):
        raise ValueError(f"poloidal_scale differs between frames: {field_file}")
    for name, attributes in frame.attributes.items():
        initial = dict(destination[name].attrs)
        if (initial.keys() != attributes.keys()
                or any(not np.array_equal(initial[key], attributes[key]) for key in initial)):
            raise ValueError(f"Dataset attributes change between frames: {name} in {field_file}")


def _input_paths(path, last_frame):
    path = Path(path)
    if last_frame is None:
        if path.is_file():
            return [path]
        candidate = Path(str(path) + ".h5")
        if candidate.is_file():
            return [candidate]
        raise FileNotFoundError(f"Input file not found: {path}")
    if not 0 <= last_frame <= 9999:
        raise ValueError("last-frame must be in 0..9999")
    stem = str(path)
    if stem.endswith("_0000.h5"):
        stem = stem[:-8]
    elif stem.endswith(".h5"):
        raise ValueError("Numbered input needs a stem or an initial _0000.h5 filename")
    return [Path(f"{stem}_{index:04d}.h5") for index in range(last_frame+1)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field", type=Path, required=True, help="stationary filename or numbered stem/0000 anchor")
    parser.add_argument("--jtor", type=Path, help="optional separate current input, using the same numbering")
    parser.add_argument("--last-frame", type=int, help="convert numbered frames 0000 through this index, inclusive")
    parser.add_argument("--output", type=Path, required=True, help="new trajectory HDF5 filename; must not exist")
    args = parser.parse_args(argv)
    try:
        fields = _input_paths(args.field, args.last_frame)
        currents = _input_paths(args.jtor, args.last_frame) if args.jtor is not None else None
        output = convert_equilibria(fields, args.output, currents)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Conversion failed: {exc}\n")
    print(f"Converted {len(fields)} frame(s) to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
