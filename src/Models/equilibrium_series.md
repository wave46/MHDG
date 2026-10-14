# Equilibrium series input

The `equilibrium_series_t` reader supports legacy stationary and numbered files on a rectangular grid
and a trajectory HDF5. Solver integration follows separately.
It owns its open file and two cached frames. `clear()` releases them; finalization
also closes the file. Create and open each provider independently.

Trajectory datasets use the following **Python/HDF5** shapes:

| Dataset | Shape | Contents |
| --- | --- | --- |
| `/time` | `(nt,)` | Physical timestamps in seconds |
| `/grid/r` | `(nr,)` | R coordinates in metres |
| `/grid/z` | `(nz,)` | Z coordinates in metres |
| `/equilibrium/psi` | `(nt,nz,nr)` | Primitive poloidal flux |
| `/equilibrium/Bphi` | `(nt,nz,nr)` | Toroidal magnetic field in tesla |
| `/equilibrium/psiSep` | `(nt,)` | Separatrix flux, same convention as psi |
| `/equilibrium/Br`, `/equilibrium/Bz` | `(nt,nz,nr)` | Optional paired reference fields in tesla |
| `/equilibrium/Jtor` | `(nt,nz,nr)` | Optional reference current density in A/m² |

The HDF5 Fortran interface reports these field dimensions as `(nr,nz,nt)`.
The frame reader therefore fills an allocated `(nr,nz)` buffer, which the provider
transposes to the existing solver frame orientation `(nz,nr)`. The generic HDF5
helper's `(n1,n2,nframes)` means `(nr,nz,nt)` for this format.

There is one shared R,Z grid; frames cannot carry different grids. Axes must
contain at least two finite, strictly increasing points. Timestamps must be
finite and strictly increasing. No `dt` dataset or configured dimensions are
needed. The reader uses hyperslabs to load individual frames, keeping only
two primitive frames in its cache. It loads the small time and psiSep vectors
when opening the file.
Opening checks field existence and shapes; field values are read on demand by
`sample_at`, only when a required frame is absent from the cache. The stored
dataset-name strings identify which arrays to read and do not themselves read data.

Use the same format for a fixed equilibrium with `nt = 1`: `/time` and
`/equilibrium/psiSep` have shape `(1,)`, and every field has shape `(1,nz,nr)`.
That single frame supplies the equilibrium at any finite requested physical time,
for both steady and transient plasma simulations. Its stored timestamp identifies
the input frame and does not limit simulation-time coverage.

`sample_at(t, sample, ierr, message)` copies an exact frame at a supplied timestamp
or linearly interpolates primitive quantities at an intermediate time. It
returns an error outside a multi-frame trajectory's coverage, allowing only roundoff at the
endpoints. Output arrays use the existing Fortran `(nz,nr)` orientation.
The sample stores physical time, zero-based left/right frame positions and
the interpolation fraction. At an exact frame, both positions identify that
frame and the fraction is zero. Topology and derived geometry are reconstructed
by the solver from the sampled primitives.

Stationary files retain the existing root datasets `r2D`, `z2D`, `flux2D`,
`Bphi2D`, `psiSep`, and optional paired `Br2D`/`Bz2D` and `Jtor`. Their existing
array orientation is preserved; `psiSep` may be a scalar or a one-element array.
They supply the same equilibrium at any finite requested time. Full filenames
and legacy extensionless names are accepted.

For paired legacy inputs, `open(path, ierr, message, jtor_path)` accepts a
separate stationary current file with root `r2D`, `z2D` and `Jtor` datasets.
Its rectangular grid may differ from the magnetic grid, as in the existing reader.
`current_grid(r,z)` supplies the current coordinates, and `sample%jtor` has
shape `(nz_current,nr_current)`. Embedded Jtor uses the magnetic coordinates.
When Jtor is absent, `current_grid` leaves both outputs unallocated. A trajectory
stores Jtor internally on its shared grid and does not take a separate current file.

`open_numbered(path, last_frame, ierr, message, jtor_path)` opens legacy frames
`0000` through `last_frame`, inclusive. Supply an extensionless stem or the full
initial filename `stem_0000.h5`; `stem.h5` identifies a single file and is not
an alias for a numbered sequence. An optional current path follows the same
numbering rule. There is no filename renumbering or intermediate file generation.

Numbered opening reads each field frame's `time` and `psiSep`, checks increasing
finite timestamps, and checks paired current timestamps when supplied. Stored
`dt` values are not interpolation coordinates. Primitive fields are loaded on
cache misses; at that point the reader checks shapes, optional-field consistency
and identical grid coordinates relative to frame 0000. Each magnetic frame must
use the same grid; each separate current frame must use its own unchanged grid.
No regridding is performed. Frame 0000's magnetic file remains open until `clear`;
other numbered files and separate current files are opened briefly for metadata
or a cache miss, then closed. Stationary and trajectory providers retain their
magnetic file handle until `clear` as well. Inputs must remain unchanged while
in use; concurrent read-only inspection is supported.

`grid(r,z)` and `time_points(times)` return copies of the provider coordinates.
Physical inputs are never normalized in the provider. Absence of reference
arrays is represented by unallocated sample components; reconstruction conventions
must be handled explicitly when the provider is connected to the solver.
All reader/sampling failures return `ierr /= 0` and a message to the caller.
