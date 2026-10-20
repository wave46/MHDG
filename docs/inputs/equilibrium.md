# Magnetic equilibrium input

MHDG can evolve the plasma in a fixed magnetic equilibrium or follow a supplied
sequence of equilibria. For rectangular-grid input, each equilibrium is defined
by poloidal flux, toroidal magnetic field and separatrix flux. The solver uses
these data to reconstruct magnetic geometry and topology on the simulation mesh.
Poloidal field and toroidal current can be reconstructed from flux or read from
the input, as selected below.

This guide focuses on rectangular-grid files, with a separate note on
[legacy nodal input](#legacy-nodal-input). Nodal and analytical magnetic
inputs use separate paths; the rectangular-grid reconstruction and COCOS
rules described here do not apply to them.

## Configuration

For an experimental-equilibrium case using a rectangular grid, set the following
entries in the existing `INPUT_LST` block of `param.txt`:

```fortran
&INPUT_LST
  field_from_grid = .true.
  field_path = 'equilibrium.h5'
  jtor_path = ''
  compute_from_flux = .true.
/
```

`equilibrium.h5` is a [single-file equilibrium](#single-file-hdf5-format).
It can contain one fixed frame or a time series. Toroidal-current data, when
supplied, belong in the same file; leave `jtor_path` empty for this format.
Paths are relative to the working directory unless absolute paths are supplied.
Legacy files remain supported, including separate current files, as described
under [legacy inputs](#legacy-inputs).

### Fixed or moving equilibrium

For explicit physical-time runs, `time_mode` in `TIME_LST` selects the numerical
stepping independently of magnetic evolution:

| Equilibrium input | Magnetic evolution | Compatible physical-time modes |
| --- | --- | --- |
| One frame | Fixed equilibrium, including during a plasma transient | `'fixed'`, `'adaptive'` |
| Multiple frames | Moving equilibrium sampled at the requested physical time | `'fixed'`, `'prescribed'`, `'adaptive'` |

A moving input must cover the run's starting time, or restored checkpoint
time, through `t_final`.
A one-frame input supplies the same equilibrium at any finite physical time;
its stored timestamp does not restrict the run.

Explicit physical-time modes require `steady = .false.`, `psdtime = .false.`
and BDF1 (`tis = 1`). Prescribed stepping uses the equilibrium timestamps as
targets; fixed and adaptive steps may fall between frames. Equilibrium-file
timestamps are always in seconds, independently of the namelist's `time_units`.
See [time sampling](#time-sampling) for interpolation and the
[restart guide](../simulation/restart-and-output.md) for physical clocks and
continuation, and [time stepping](../simulation/time-stepping.md) for recovery settings.

Configurations without an explicit `time_mode` retain the
[legacy frame/step relationship](#legacy-time-stepping).

### Stored or flux-derived fields

Both choices require poloidal flux, toroidal field and separatrix flux. They
differ in how the fields applied to the plasma equations are obtained:

| `compute_from_flux` | Poloidal field | Toroidal current, when ohmic heating is enabled |
| --- | --- | --- |
| `.true.` | Reconstructed from poloidal flux | Reconstructed from the Grad–Shafranov operator on flux |
| `.false.` | Supplied `Br` and `Bz` arrays | Supplied `Jtor` array |

In flux-derived mode, supplied `Br`/`Bz` and `Jtor` act as references for the
reconstruction coefficients. Missing reference information can be supplied by
a `cocos` declaration, following the [selection rules](#selecting-the-reconstruction-coefficients).
In stored-field mode, `cocos` cannot replace missing field/current arrays.
Current data and its reconstruction coefficient are not required when ohmic
heating is inactive (`ohmicsrc = .false.`).

Geometry and topology are reconstructed from flux in both modes. With stored
fields, the geometry derived from flux and the supplied field arrays therefore
need to be physically consistent.

## Single-file HDF5 format

The same layout represents fixed and moving equilibria. The shapes below are
those seen from Python/HDF5: `nt` time frames, `nr` radial grid points and `nz`
vertical grid points.

| Dataset | Shape | Physical quantity and units |
| --- | --- | --- |
| `/time` | `(nt,)` | Physical timestamps in seconds |
| `/grid/r` | `(nr,)` | Cylindrical radius $R$ in metres |
| `/grid/z` | `(nz,)` | Vertical coordinate $Z$ in metres |
| `/equilibrium/psi` | `(nt,nz,nr)` | Poloidal flux $\psi$, in the source's flux convention |
| `/equilibrium/Bphi` | `(nt,nz,nr)` | Toroidal magnetic field $B_\phi$ in tesla |
| `/equilibrium/psiSep` | `(nt,)` | Separatrix flux, in the same convention as $\psi$ |
| `/equilibrium/Br`, `/equilibrium/Bz` | `(nt,nz,nr)` | Optional paired poloidal fields $B_R$, $B_Z$ in tesla |
| `/equilibrium/Jtor` | `(nt,nz,nr)` | Optional toroidal current density $J_\phi$ in A/m² |
| `/equilibrium/cocos` | Scalar integer | Optional source COCOS identifier: 1–8 or 11–18 |

The first six datasets are required. Whether the optional field/current arrays
can be omitted depends on the selected field mode and convention information,
as described above. `Br` and `Bz` must either both be present or both be absent.

All frames share one rectangular $R,Z$ grid. Each axis must contain at least
two finite, strictly increasing coordinates. Times must be finite and strictly
increasing, and field values and separatrix flux must be finite. Grid sizes
are inferred from the file; no configured dimensions or `dt` dataset are needed.

For a fixed equilibrium, use `nt = 1`: time and separatrix flux still have
shape `(1,)`, and field arrays have shape `(1,nz,nr)`.

Supply dimensional coordinates, fields and flux. In particular, `psi` is the
source poloidal flux, not a normalized flux label. Its units are Wb/rad for
COCOS 1–8 and Wb for COCOS 11–18; signs and normalization are handled through
the reconstruction coefficients below. Do not normalize the input arrays
using the plasma's solver reference scales.

`cocos` is the supported convention declaration. It may be a scalar or an array
containing exactly one integer value. It describes the source convention and
does not change during time sampling. There is no assumed default COCOS.

## Time sampling

At a supplied timestamp, the solver uses that frame. Between two frames at
$t_i$ and $t_{i+1}$, the input quantities are interpolated linearly. For example,

$$
w = \frac{t-t_i}{t_{i+1}-t_i},
\qquad
\psi(R,Z,t) = (1-w)\psi_i(R,Z) + w\psi_{i+1}(R,Z).
$$

The same temporal interpolation applies to toroidal field, separatrix flux
and any supplied reference fields/current. Geometry and topology are then
reconstructed from the sampled flux; they are not interpolated between
precomputed mesh geometries. A multi-frame input is not extrapolated beyond
its time range, apart from endpoint roundoff tolerance.

In explicit physical-time modes, sampling uses the proposed timestep endpoint
in seconds. Numerical step numbers do not select equilibrium frames, and
sampling does not change the numerical timestep. Checkpoints retain the
physical time and the bracketing frame indices/fraction; continuation locates
the equilibrium from physical time.

Fixed fields are retained between steps. After mesh adaptation, their
representation is rebuilt on the new mesh from the same input data. Moving
fields are also rebuilt on the current mesh. Heating can have its own grid
and timestamps; see the [heating input guide](heating.md).

## Flux reconstruction and conventions

### Physical expressions

For physical cylindrical coordinates $R$ (outward) and $Z$ (upward), the
poloidal magnetic field is reconstructed as

$$
B_R = -\frac{\alpha_B}{R}\frac{\partial\psi}{\partial Z},
\qquad
B_Z = \frac{\alpha_B}{R}\frac{\partial\psi}{\partial R},
\qquad
B_{\mathrm{pol}} = \sqrt{B_R^2+B_Z^2}.
$$

The toroidal current density is

$$
\Delta^*\psi = \frac{\partial^2\psi}{\partial R^2}
- \frac{1}{R}\frac{\partial\psi}{\partial R}
+ \frac{\partial^2\psi}{\partial Z^2},
\qquad
J_\phi = -\frac{\alpha_J}{\mu_0 R}\Delta^*\psi.
$$

Derivatives are with respect to physical coordinates in metres. Fields are in
tesla, and $J_\phi$ is the physical toroidal current density in A/m². The code
uses $\mu_0 = 4\pi\times10^{-7}$ H/m. Toroidal magnetic field is read directly
from the input.

The dimensionless coefficients $\alpha_B$ and $\alpha_J$ account for the flux
normalization and source sign convention. They are resolved internally, rather
than entered as separate parameters. Each coefficient is selected from
$+1$, $-1$, $+1/(2\pi)$ and $-1/(2\pi)$.

The correspondence with input and code names is:

| Input or code name | Physical notation | Meaning |
| --- | --- | --- |
| `psi` (`flux2D` in legacy files) | $\psi$ | Source poloidal flux |
| `psiSep` | $\psi_{\mathrm{sep}}$ | Input separatrix flux |
| `Br`, `Bz`, `Bphi` | $B_R$, $B_Z$, $B_\phi$ | Cylindrical magnetic-field components |
| `Jtor` | $J_\phi$ | Toroidal current density |
| `alpha_B`, `alpha_J` | $\alpha_B$, $\alpha_J$ | Applied field/current reconstruction coefficients |

Flux-derived current is evaluated only in the detected plasma core and is set
to zero at solver nodes outside it. This core restriction belongs to the
reconstructed-current path; stored-field mode instead interpolates the
supplied current.

### COCOS declaration

The convention identifier follows
[Sauter & Medvedev (2013), *Tokamak coordinate conventions: COCOS*](https://doi.org/10.1016/j.cpc.2012.09.010),
Table I and Eqs. 19, 25, 27 and 30. For the component formulas above, the
lookup is:

| COCOS: reduced flux / full flux | $\alpha_B$ | $\alpha_J$ | Handedness $h$ |
| --- | --- | --- | --- |
| 1,5 / 11,15 | $-k$ | $-k$ | $+1$ |
| 2,6 / 12,16 | $+k$ | $-k$ | $-1$ |
| 3,7 / 13,17 | $+k$ | $+k$ | $+1$ |
| 4,8 / 14,18 | $-k$ | $+k$ | $-1$ |

Here $k=1$ for COCOS 1–8 and $k=1/(2\pi)$ for 11–18, with
$\alpha_J=h\alpha_B$. Handedness $h=+1$ describes a right-handed
$(R,\phi,Z)$ basis, and $h=-1$ a right-handed $(R,Z,\phi)$ basis. The indices
paired in each row differ in poloidal-angle conventions that do not change
these component formulas.

### Selecting the reconstruction coefficients

Usable reference fields take precedence over the declaration. A paired
$B_R,B_Z$ reference determines $\alpha_B$; a current reference determines
$\alpha_J$ when ohmic heating is enabled. With one usable reference and a
COCOS declaration, the reference sets the coefficient's magnitude and sign,
while COCOS handedness supplies the relation to the missing coefficient.

The full selection rules, when current is needed for ohmic heating, are:

| Usable reference information | Without `cocos` | With `cocos` |
| --- | --- | --- |
| Poloidal field and current | Infer both coefficients independently | Infer both independently; warn if the declaration conflicts |
| Poloidal field only | Stop: current coefficient unresolved | Infer $\alpha_B$; use $\alpha_J=h\alpha_B$ |
| Current only | Stop: field coefficient unresolved | Infer $\alpha_J$; use $\alpha_B=h\alpha_J$ |
| Neither | Stop | Take both coefficients from COCOS |

When ohmic heating is inactive, current information is not required: a usable
poloidal-field reference alone is sufficient without `cocos`. A reference
whose fit cannot be evaluated supplies no coefficient information. Missing or
unusable references produce warnings when the declaration supplies the missing
information.

Conflicting COCOS metadata produces a warning while retaining reference-derived
coefficients. If field and current references imply different flux
normalizations, each independently selected coefficient is retained with a
warning; that pair does not represent a single consistent COCOS convention.
The separate fit-residual checks still apply: a valid poloidal-field fit with
an applied relative RMS residual above 25% stops a flux-derived run.

This selection preserves source conventions. It does not transform toroidal
field, flows, signed safety factor or drift operators into a common MHDG
coordinate convention. Full coordinate conversion remains necessary before
relying on sign-sensitive drift physics with arbitrary input conventions.

## Legacy inputs

### Stationary rectangular grids

Stationary magnetic files store the following datasets at the HDF5 root.
Shapes use Python/HDF5 notation; the coordinate arrays are also two-dimensional.

| Dataset | Shape | Physical quantity and units |
| --- | --- | --- |
| `r2D` | `(nr,nz)` | Cylindrical radius $R$ in metres |
| `z2D` | `(nr,nz)` | Vertical coordinate $Z$ in metres |
| `flux2D` | `(nr,nz)` | Poloidal flux $\psi$, in the source's flux convention |
| `Bphi2D` | `(nr,nz)` | Toroidal magnetic field $B_\phi$ in tesla |
| `psiSep` | Scalar or one-element array | Separatrix flux, in the same convention as $\psi$ |
| `Br2D`, `Bz2D` | `(nr,nz)` | Optional paired poloidal fields $B_R$, $B_Z$ in tesla |
| `Jtor` | `(nr,nz)` | Optional toroidal current density $J_\phi$ in A/m² |
| `cocos` | Scalar or one-element array | Optional integer source COCOS identifier: 1–8 or 11–18 |

Requirements for the optional fields follow the same
[field-mode and convention rules](#stored-or-flux-derived-fields) as the
single-file format.

If current is supplied separately, select its file with `jtor_path`. Its grid
may differ from the magnetic grid, with datasets at the HDF5 root:

| Dataset | Shape | Physical quantity and units |
| --- | --- | --- |
| `r2D` | `(nr_current,nz_current)` | Cylindrical radius $R$ in metres |
| `z2D` | `(nr_current,nz_current)` | Vertical coordinate $Z$ in metres |
| `Jtor` | `(nr_current,nz_current)` | Toroidal current density $J_\phi$ in A/m² |

Without a separate current path, embedded `Jtor` can be used. Single-file
trajectories share one grid for all quantities and embed current data.

Stationary files supply the same equilibrium at any finite requested time.
Full `.h5` filenames and legacy extensionless names are accepted; full names
make the intended input explicit.

### Numbered files

Legacy moving inputs use files numbered from `0000`, such as
`equilibrium_0000.h5`, `equilibrium_0001.h5`, and so on. Supply the extensionless
stem or the full `_0000.h5` anchor. `stem.h5` selects a single file and is not
an alias for a numbered sequence. A separate numbered current path follows
the same naming rule.

Each frame supplies `time` and `psiSep`. Timestamps must be finite and strictly
increasing; paired current timestamps must match. Stored `dt` datasets do not
determine rectangular-grid interpolation or timestep intervals.

All magnetic frames must retain the same grid, layout, optional fields and
COCOS declaration. If current has a separate grid, that grid must also remain
unchanged. No regridding is performed. In explicit physical-time modes, the
sequence is discovered from frame 0000 through the last contiguous file,
independently of the numerical step count.

### Legacy time stepping

Legacy moving-equilibrium runs without an explicit `time_mode` associate
frames with numerical steps: frame 0000 is the initial equilibrium at time
zero, and frames 1 through `nts` are timestep endpoints. A trajectory must
provide at least `nts + 1` frames. Intervals come from successive physical
timestamps. Restart at step $k$ rebuilds frame $k$ before advancing to $k+1$.
Reaching the time-convergence tolerance does not stop a moving trajectory
early: all remaining steps through `nts` are taken.

Legacy fixed-equilibrium runs use a stationary file or a one-frame trajectory;
the plasma can still be transient. For new physical-time runs, select the
numerical stepping independently of the equilibrium's frame count as described
under [fixed or moving equilibrium](#fixed-or-moving-equilibrium).

### Legacy nodal input

In the 2D experimental-equilibrium path (`testcase` 50–59),
`field_from_grid = .false.` selects global-node input. Each magnetic file
contains:

| Dataset | Contents |
| --- | --- |
| `Br`, `Bz`, `Bt` | Magnetic components in tesla, as vectors in global mesh-node order |
| `flux` | Poloidal-flux vector in the same node order |
| `psiSep` | Separatrix flux, in the same units and convention as `flux` |
| `time`, `dt` | Physical time and interval in seconds for numbered moving-equilibrium files |

MPI partitions read the same global file and select nodes through the mesh
mapping. Separate per-rank equilibrium files are not used. Nodal input
retains its separate legacy current and frame-reading path.

## Converting legacy files

The [converter](../../tools/convert_equilibrium.py) needs Python with NumPy
and h5py. From the repository root, convert a stationary equilibrium/current
pair with:

```bash
python tools/convert_equilibrium.py --field equilibrium.h5 --jtor current_density.h5 --output equilibrium_series.h5
```

For numbered frames 0000 through N, inclusive, supply `--last-frame`:

```bash
python tools/convert_equilibrium.py --field equilibrium_0000.h5 --jtor current_density_0000.h5 --last-frame 28 --output discharge.h5
```

Numbered stems are also accepted. Without `--last-frame`, a full filename
selects one stationary frame even if it ends in `_0000.h5`; an extensionless
name resolves its stationary `.h5` file. Omit `--jtor` to retain embedded
legacy current or leave the reference absent. Poloidal references remain an
optional pair.

Conversion preserves physical values and conventions. It does not normalize,
reconstruct, retime or regrid the equilibrium. All frames, including separate
current files, must share identical rectangular coordinates to fit the
single-file layout. A legacy current file on another grid needs regridding
outside this converter.

Moving frames require increasing finite timestamps and matching paired-current
timestamps. For one stationary field without a timestamp, the stored time is
zero; otherwise its timestamp is retained. Paired stationary timestamps are
compared when both exist.

An explicit `cocos` declaration and its attributes are copied. Its presence and
value must remain constant across frames. If absent, it remains absent; the
converter does not infer a convention from reference fields. Conversion alone
therefore does not resolve an unspecified flux normalization.

The destination must not already exist. It is published only after successful
conversion, so a failed conversion leaves no partial output file. Source files
remain unchanged.

## Checking the applied equilibrium

Startup reports magnetic topology, input/effective separatrix flux, axis
coordinates and minor radius. Rank zero reports global nodal ranges of $B_R$,
$B_Z$ and $B_\phi$ in tesla, and $J_\phi$ in A/m² when active. These are the
fields applied on the current mesh.

Reference fits report unconstrained coefficients, selected coefficients and
residuals. In HDF5, the `/magnetic` group distinguishes a successful reference
comparison from a coefficient selected using convention information:

| Dataset | Interpretation |
| --- | --- |
| `field_fit_valid`, `jtor_comparison_valid` | Whether a usable reference comparison was available |
| `field_fit_alpha`, `jtor_fit_alpha` | Unconstrained reference-fit coefficients |
| `field_applied_alpha`, `jtor_applied_alpha` | Coefficients actually selected |
| `field_alpha_convention_id` | Internal identifier of the selected field coefficient |
| `jtor_applied_sign` | Sign of the applied current coefficient relative to $\alpha_B$ |

The applied-current residual includes the independently selected magnitude;
raw-current comparison diagnostics retain $\alpha_B$ as their reference.
Checkpoint frame indices describe the sampled equilibrium position. The
physical clock remains authoritative on restart, as explained in the
[restart guide](../simulation/restart-and-output.md).

## Implementation notes

### Sampling and file ownership

The [equilibrium provider](../../src/Models/equilibrium_series.f90) returns
dimensional samples. Spatial interpolation, reconstruction, geometry/topology
analysis and solver normalization belong to
[magnetic_field.f90](../../src/Models/magnetic_field.f90) and
[magnetic_geometry_state.f90](../../src/Models/magnetic_geometry_state.f90).

`equilibrium_series_t` contains a private `rectangular_series_t` from
[rectangular_series_io.f90](../../src/InOut/rectangular_series_io.f90), also
used by heating. The shared object owns the file, coordinate/time metadata,
shape checks, time brackets and two-frame cache bookkeeping. Equilibrium owns
its named primitive arrays, separatrix flux, convention/reference validation
and separate legacy-current-grid handling. A cache slot becomes valid only
after its frame has been read successfully.

Opening reads coordinates, time and separatrix-flux vectors and checks dataset
existence/shapes. Field values are loaded on demand through hyperslabs, with
only two primitive frames cached. Numbered input checks timestamps on opening;
grids, optional fields and convention consistency are checked as frames are
loaded. Metadata comparisons of selected coefficients use relative tolerance
`1.e-6`, rather than the unconstrained fit values. Current reference fitting
uses owned-element core quadrature samples.

Providers and their cached raw data survive mesh adaptation. Keep inputs
unchanged while in use; concurrent read-only inspection is supported.
Stationary/trajectory files and the first numbered magnetic file remain open
until `clear()` or finalization. Other numbered and separate-current files are
opened briefly for metadata or a cache miss. The solver closes its provider
on normal completion. Create/open providers independently rather than copying
objects that own file handles.

### Array orientation and interface

For the trajectory's Python/HDF5 shape `(nt,nz,nr)`, the Fortran HDF5 interface
reports `(nr,nz,nt)`. The frame reader fills `(nr,nz)` and transposes it to the
solver sample orientation `(nz,nr)`. The generic frame helper's
`(n1,n2,nframes)` is therefore `(nr,nz,nt)`. Legacy `(nr,nz)` Python matrices
retain their existing orientation on reading.

| Provider method | Purpose |
| --- | --- |
| `open(path, ierr, message, jtor_path)` | Open a single file, with an optional separate legacy-current input |
| `open_numbered(path, last_frame, ierr, message, jtor_path)` | Open legacy frames 0000 through `last_frame`, inclusive |
| `sample_at(t, sample, ierr, message)` | Return an exact/interpolated sample at physical time `t` |
| `grid(r,z)`, `time_points(times)` | Return copies of input coordinates/times |
| `current_grid(r,z)` | Return current coordinates; leave outputs unallocated if current is absent |
| `clear()` | Release metadata/cache and close the owned file; finalization also clears it |

Samples are independent copies of cached primitives. Optional reference arrays
and `sample%cocos` are unallocated when absent. Separate-current samples have
shape `(nz_current,nr_current)`; embedded current uses the magnetic grid.
Samples record physical time, zero-based `frame_left`/`frame_right` and
`fraction`. At an exact frame, both indices equal that frame and the fraction
is zero. Reader/sampling failures return a nonzero `ierr` and message; solver
input failures terminate the run.

### Converter storage and provenance

Conversion processes one frame at a time, writing double precision with
lossless compression and one frame per chunk. Dataset attributes, including
units and convention metadata, are retained and must agree across frames.
Source root attributes are retained under `/provenance/source_attributes`;
source basenames are stored in `/provenance/magnetic_files` and, when supplied,
`/provenance/current_files`. Unrelated legacy datasets and `dt` are not copied.

Python callers can use
`tools.convert_equilibrium.convert_equilibria(field_files, output, current_files)`
with ordered file lists; those lists select frames while retaining timestamps.
