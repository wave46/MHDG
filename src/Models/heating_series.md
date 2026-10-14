# Heating series input

`heating_series_t` reads legacy stationary heating and a single trajectory file.
It has no dependency on magnetic equilibrium, solver globals or MPI. Solver
integration and legacy numbered-file scheduling follow in the next slice;
`external_heating` currently retains its existing loading path.

The trajectory layout uses these Python/HDF5 shapes:

| Dataset | Shape | Contents |
| --- | --- | --- |
| `/time` | `(nt,)` | Physical timestamps in seconds |
| `/grid/r` | `(nr,)` | R coordinates in metres |
| `/grid/z` | `(nz,)` | Z coordinates in metres |
| `/heating/ions` | `(nt,nz,nr)` | Ion heating power density in W/m³ |
| `/heating/electrons` | `(nt,nz,nr)` | Electron heating power density in W/m³ |

Both power datasets are required. Use an explicit zero array for a channel with
no heating. The shared grid must have finite, strictly increasing R,Z axes with
at least two points each. Times must be finite and strictly increasing. Heating
and magnetic trajectories can have different grids and different timestamps.

`open(path, ierr, message)` resolves a full filename or legacy stationary stem,
reads coordinate/time metadata and checks power dataset shapes. `sample_at(t,
sample, ierr, message)` loads only missing frames into a two-frame cache and
linearly interpolates both powers. Cached primitives remain unchanged when a
caller modifies its returned sample. Samples use Fortran `(nz,nr)` arrays named
`ions` and `electrons`, physical `time`, zero-based `frame_left`/`frame_right`
and `fraction`, with the same endpoint and coverage rules as the
[equilibrium provider](equilibrium_series.md).

Use the same trajectory format with `nt=1` for fixed heating: power datasets
still have shape `(1,nz,nr)`. A single frame supplies heating at any finite
requested time, regardless of its stored timestamp. Legacy stationary files
retain `r2D`, `z2D`, `heating2D_ions` and `heating2D_electrons`, with the existing
Python `(nr,nz)` orientation. Dimensions are inferred; legacy `ip`/`jp` datasets
are not needed. A legacy stationary file's stored time, if present, does not
constrain sampling.

`grid(r,z)` and `time_points(times)` return coordinate copies. `clear()` closes
the read-only input file and releases metadata/cache; finalization also clears
it. Input must remain unchanged while open; concurrent read-only access is
supported. Errors return status and a message rather than stopping the solver.

The provider returns dimensional primitives. It applies no normalization or
spatial interpolation. Solver integration will retain the existing power scale
`refval_time / (refval_specenergydens * refval_mass)` and rebuild nodal heating
when the requested sample or mesh changes. Each provider contains a private `rectangular_series_t` from
`src/InOut/rectangular_series_io.f90`. That shared object owns the file, grid/time
metadata, component shape checks, frame orientation, time brackets and cache-slot
bookkeeping. Heating retains its named power arrays and interpolation; equilibrium
retains its own quantities, optional-reference rules and legacy Jtor handling.
A slot becomes valid only after its provider finishes reading that frame.
`file()` supplies a borrowed HDF5 handle for specialized metadata; the shared
object owns closing it. Create/open independent providers rather than copying
objects that own open handles.
