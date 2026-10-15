# Equilibrium series input

The `equilibrium_series_t` reader supports legacy stationary and numbered files on a rectangular grid
and a trajectory HDF5. The rectangular-grid magnetic and Jtor loaders use this provider;
spatial interpolation, normalization, derived fields and topology stay in the solver.
It owns its open file and two cached frames. `clear()` releases them; finalization
also closes the file. Create and open each provider independently.

Equilibrium and heating each contain a private `rectangular_series_t`. It owns
their shared file/grid/time state, common metadata and shape checks, time
bracketing and cache-slot management. Each provider owns its named primitive
arrays and format-specific validation; equilibrium also retains `psiSep`, numbered
input and separate legacy Jtor grids. Successful frame reads mark a reserved
slot valid; a failed read leaves that slot invalid. The public sampling API and
file layout are unchanged by this composition.

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
| `/equilibrium/cocos` | scalar integer | Optional source COCOS identifier, 1–8 or 11–18 |

`cocos` is the only optional convention declaration. It identifies the source
convention using the standard COCOS table; legacy files may store it at the root.
Store it as a scalar or an array containing exactly one integer value. Its presence
and value must stay constant across numbered frames and are unchanged by sampling.
The provider returns it as the optional allocated `sample%cocos`; an absent
identifier remains unallocated. Usable Br/Bz and Jtor references can determine both
reconstruction coefficients without any convention declaration.

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

## Solver input

For the single-file format, set `field_path = 'equilibrium.h5'` and omit
`jtor_path` (or set it to an empty string). Reference Jtor is embedded in the
same file. Existing legacy inputs keep their separate `field_path` and
`jtor_path`; full stationary names and numbered stems or `_0000.h5` anchors
remain supported. Embedded Jtor also works in a legacy stationary file when
no separate current path is supplied.

For legacy configurations without `time_mode`, `ME = .true.` visits frames 0000 through
`nts`, where 0000 represents the initial state at time zero and subsequent
frames are timestep endpoints. A trajectory must therefore supply at least
`nts + 1` frames. Each interval is computed from successive physical timestamps;
the legacy `dt` datasets no longer determine rectangular-grid timesteps.
Restart at step k rebuilds frame k before advancing to k+1. The provider and
its raw-frame cache survive mesh adaptation; rebuilding the representation on
the new mesh does not reopen a trajectory or reread cached frames.

With `ME = .false.`, supply a stationary legacy file or a one-frame trajectory.
This does not require steady plasma: existing fixed-equilibrium transients
remain supported. Independent controller scheduling and adaptive timesteps
follow in later slices. Legacy multi-frame trajectories require ME.

With an explicit `time_mode`, the input determines magnetic evolution: one frame
supplies a fixed equilibrium and multiple frames supply a moving equilibrium.
`ME` is then inferred internally. Legacy numbered sequences are discovered from
frame 0000 through the last contiguous file, independently of the numerical step
count. `prescribed` mode requires multiple frames; fixed numerical steps may use
either fixed or moving equilibrium. A moving input must cover the initial and
final physical times. Sampling uses the proposed solver endpoint in seconds,
including interpolation between supplied frames; numerical step numbers do not
select magnetic frames. Sampling leaves the numerical timestep unchanged.
Fixed fields are retained between steps. Mesh adaptation rebuilds their nodal
representation from the same provider and cache, as for moving fields.
Explicit-mode solver execution remains guarded until controller scheduling and
restart integration are complete.

With `compute_from_flux = .true.`, the poloidal magnetic field and toroidal
current density are reconstructed from the primitive flux `psi`. For physical
cylindrical coordinates R (outward) and Z (upward), the expressions are:

```text
Br     = -alpha_B/R * ∂psi/∂Z
Bz     =  alpha_B/R * ∂psi/∂R
Bpol   = sqrt(Br² + Bz²)

Δ*psi  = ∂²psi/∂R² - (1/R)*∂psi/∂R + ∂²psi/∂Z²
Jtor   = -alpha_J/(mu0*R) * Δ*psi
```

Here R and Z are in metres, Br/Bz/Bpol in tesla, and Jtor in A/m²;
`mu0 = 4*pi*10^-7 H/m` is the permeability used by the code. Derivatives in these
expressions are with respect to physical coordinates. `psi` keeps the input's
flux sign and normalization: COCOS 1–8 use reduced flux (Wb/rad), while 11–18
use full flux (Wb). The dimensionless coefficients `alpha_B` and `alpha_J` account
for that normalization and the source's sign conventions. They are selected
internally, not separate input parameters. Bphi is read directly from the input.

Each reference fit independently selects one of `+1`, `-1`, `+1/(2*pi)` or
`-1/(2*pi)`; unconstrained fits and residuals are reported for comparison.
The current fit uses the existing owned-element core quadrature comparison, and
reconstructed current is set to zero at solver nodes outside the detected plasma core.
`Jtor` denotes the physical toroidal current density, not a covariant component,
parallel current, or integrated current.

The COCOS lookup supplies both coefficients and cylindrical handedness `h`:

| COCOS: reduced flux / full flux | alpha_B | alpha_J | h |
| --- | --- | --- | --- |
| 1,5 / 11,15 | -k | -k | +1 |
| 2,6 / 12,16 | +k | -k | -1 |
| 3,7 / 13,17 | +k | +k | +1 |
| 4,8 / 14,18 | -k | +k | -1 |

Here `k=1` for 1–8 and `k=1/(2*pi)` for 11–18; `alpha_J=h*alpha_B`.
`h=+1` describes a right-handed (R,phi,Z) basis, and `h=-1` a right-handed
(R,Z,phi) basis. Indices paired in each row differ in poloidal-angle conventions
that do not change these component formulas. See
[Sauter & Medvedev (2013), Table I and Eqs. 19, 25, 27, 30](https://doi.org/10.1016/j.cpc.2012.09.010).

Reference fits take precedence over `cocos`. With only one usable reference and
`cocos`, the measured coefficient sets the magnitude and sign; COCOS handedness
supplies the missing coefficient. A metadata mismatch warns and continues. Different reference-derived
magnitudes also warn and continue, retaining each independently selected coefficient;
that pair is not a single consistent COCOS convention. Metadata comparisons use
canonical selected coefficients with relative tolerance `1e-6`, rather than noisy
unconstrained fits. The existing excessive-residual checks remain in effect.

The complete selection rules are below. B means a usable paired Br/Bz fit and J
means a usable Jtor fit. A reference whose fit cannot be evaluated supplies no
coefficient information.

| Reference information | Without `cocos` | With `cocos` |
| --- | --- | --- |
| B and J | Infer both independently | Infer both independently; warn if metadata conflicts |
| B only | Stop if current is needed | Infer alpha_B; alpha_J=h*alpha_B |
| J only | Stop: alpha_B unresolved | Infer alpha_J; alpha_B=h*alpha_J |
| Neither | Stop | Both coefficients from COCOS |

Unresolved current does not block runs with ohmic heating inactive: a usable B
reference alone is sufficient. Missing or unusable references produce warnings
when COCOS supplies the missing information. There is no assumed default COCOS.

This lookup preserves source conventions. It does not transform Bphi, flows,
signed safety factor or drift operators into a global HDG coordinate convention.
A full coordinate conversion, needed before relying on sign-sensitive drift physics
with arbitrary input conventions, remains separate work.

With `compute_from_flux = .false.`, stored Br/Bz are required, and stored Jtor
is required when ohmic heating is active. A declaration cannot replace arrays
in stored-field mode. Nodal/analytical magnetic input retains its existing path,
including its separate legacy current input.

Rank zero reports global minimum/maximum nodal Br, Bz and Bphi in tesla, and
Jtor in A/m² when active. These describe the fields actually applied on the
current mesh. Output continues to distinguish valid reference fits from an
explicitly selected scale through `field_fit_valid`, `field_applied_alpha`,
`field_alpha_convention_id` and `jtor_comparison_valid`. `jtor_fit_alpha` and
`jtor_applied_alpha` record the unconstrained and selected current coefficients.
The legacy `jtor_applied_sign` remains the sign relative to alpha_B; raw current
comparison diagnostics also retain alpha_B as their reference. The applied current
residual includes the independently selected magnitude.

The solver closes the provider on normal completion. Failed input validation
terminates the run; recoverable nonlinear trial statuses are a later stage.

## Legacy conversion

The converter needs Python with NumPy and h5py. From the repository root, convert
a stationary equilibrium/current pair:

```bash
python tools/convert_equilibrium.py --field equilibrium.h5 --jtor current_density.h5 --output equilibrium_series.h5
```

Convert numbered frames 0000 through N, inclusive, by supplying `--last-frame`:

```bash
python tools/convert_equilibrium.py --field equilibrium_0000.h5 --jtor current_density_0000.h5 --last-frame 28 --output discharge.h5
```

Numbered stems are also accepted. Without `--last-frame`, a full filename selects
one stationary frame, even if its name ends in `_0000.h5`; an extensionless name
resolves its stationary `.h5` file. Omit `--jtor` to use embedded legacy Jtor or
leave the reference absent. Br/Bz remain paired optional references.

The output uses the same `(nt,nz,nr)` layout for stationary and moving data.
Conversion reads and writes one frame at a time, uses double precision and
lossless compression with one frame per chunk, and leaves inputs unchanged.
It preserves physical values and conventions without normalization, regridding,
retiming or reconstruction. Moving frames require finite increasing timestamps
and matching paired-current timestamps. For a stationary field without `time`,
its one stored timestamp is zero; otherwise its original timestamp is retained.
Paired stationary timestamps are compared when both are present.

All frames, including separate Jtor files, must share identical rectangular R,Z
coordinates. The legacy solver can still read a separate current grid, but such
input cannot be converted to this shared-grid format without regridding.
Dataset attributes, including provided units/convention metadata, are retained
and must remain consistent across frames. Source root attributes are preserved
per frame under `/provenance/source_attributes`; source basenames are recorded
in `/provenance/magnetic_files` and, when supplied, `/provenance/current_files`.
Unrelated legacy datasets and `dt` are not copied. Do not infer an unspecified
flux normalization from conversion: the primitives and reference fields retain
the source convention.

An explicit legacy `cocos` identifier is copied to the trajectory scalar, including
its dataset attributes. Every source frame must retain its presence and value;
conversion rejects a change of convention. An absent identifier remains absent
rather than being inferred from references.

The destination must not already exist. A failed conversion leaves no partial
published file. Python callers, including regression preparation, can use
`tools.convert_equilibrium.convert_equilibria(field_files, output, current_files)`
with ordered file lists; callers choose the frames, and conversion retains their
timestamps. The regression harness migration to this storage follows separately.
