# External heating input

External heating supplies prescribed ion and electron power-density profiles
to the plasma energy equations. Each channel specifies deposited power per unit
volume in W/m³. Profiles can be fixed or vary in time, with either fixed or
moving magnetic equilibrium.

This guide covers heating supplied on a rectangular cylindrical $R,Z$ grid.
The ion and electron channels are prescribed separately; positive values add
energy to the corresponding equation.

In the 2D NGammaTiTe family, imported heating is added in testcases 50–69,
on top of any testcase-specific sources. Manufactured convergence cases are
excluded. The separate analytical Gaussian-heating option retains its
50–59 scope.

## Configuration

Set these entries in the existing namelist blocks of `param.txt`:

```fortran
&SWITCH_LST
  external_heating = .true.
/
&INPUT_LST
  external_heating_from_grid = .true.
  external_heating_path = 'heating.h5'
/
```

Paths are relative to the working directory unless absolute paths are supplied.
Both power channels must be present in the file. Supply a zero array for a
channel receiving no heating. Nodal heating input is not implemented, so keep
`external_heating_from_grid = .true.` for this input path.

For fixed heating, use a one-frame single file or a legacy stationary file.
Both supply a profile fixed in physical space, regardless of magnetic evolution.

For time-dependent heating, two storage formats are supported:

| Property | Multi-frame single file | Legacy numbered files |
| --- | --- | --- |
| Files | One `heating.h5` containing the whole history | One file per supplied frame: `heating_0000.h5`, `heating_0001.h5`, … |
| Timestamps | The heating file's own `/time` vector | Each file's own scalar `/time`, in seconds |
| Power-array shape in Python/HDF5 | `(nt,nz,nr)` | `(nr,nz)` in each file |
| Magnetic input | Fixed or moving equilibrium | Fixed or moving equilibrium |

Both formats use heating timestamps independently of magnetic-frame times.
For example, heating frames at 0.00, 0.04 and 0.30 s can accompany an
equilibrium history sampled at different times.

Both formats interpolate heating linearly between supplied frames. A frame
index labels input data; it does not prescribe the numerical step number.

## Single-file HDF5 format

The same layout represents fixed and time-dependent heating. Shapes below use
Python/HDF5 notation: `nt` time frames, `nr` radial points and `nz` vertical
points.

| Dataset | Shape | Physical quantity and units |
| --- | --- | --- |
| `/time` | `(nt,)` | Physical timestamps in seconds |
| `/grid/r` | `(nr,)` | Cylindrical radius $R$ in metres |
| `/grid/z` | `(nz,)` | Vertical coordinate $Z$ in metres |
| `/heating/ions` | `(nt,nz,nr)` | Ion heating power density $S_{\mathrm{heating},i}$ in W/m³ |
| `/heating/electrons` | `(nt,nz,nr)` | Electron heating power density $S_{\mathrm{heating},e}$ in W/m³ |

All five datasets are required. Both channels and all frames share one grid.
Each axis must contain at least two finite, strictly increasing coordinates.
Timestamps must be finite and strictly increasing; power-density values must
be finite. Heating and magnetic inputs may have different grids and different
timestamps.

For fixed heating, use `nt = 1`: time has shape `(1,)`, and both power arrays
have shape `(1,nz,nr)`. That frame supplies heating at any finite requested
time, regardless of its stored timestamp.

Supply coordinates and power densities in physical units. File timestamps are
in seconds independently of the namelist's `time_units`. Grid sizes are inferred
from the file.

## Time sampling and restart

With explicit physical `time_mode`, heating is evaluated at the proposed
timestep endpoint. Both channels use the same bracketing heating frames
and linear interpolation, with the same
[interpolation and endpoint-coverage rules](equilibrium.md#time-sampling)
as equilibrium input. A multi-frame heating history must cover the requested
times; there is no extrapolation apart from endpoint roundoff tolerance.
Heating-frame positions are independent of numerical step numbers.

For a steady plasma calculation in fixed equilibrium, the heating history is
sampled at the initial simulation time. During a physical transient, it follows
the physical clock, including when the magnetic field remains fixed.

On physical checkpoint continuation, heating is evaluated at the restored
physical time before advancing. It then follows subsequent timestep endpoints.
Restart classification and completion behavior are described in the
[restart guide](../simulation/restart-and-output.md).

Fixed heating is retained between steps. Mesh adaptation rebuilds the profiles
on the new mesh from the dimensional input data. Time-dependent profiles are
also evaluated on the current mesh.

## Physical units and normalization

Let $m_0$, $n_0$, $u_0$ and $t_0$ be the solver's reference mass, density,
speed and time. The energy-density and power-density scales are

$$
\begin{aligned}
\mathcal E_0 &= m_0 n_0 u_0^2, \\
S_{\mathrm{heating},0} &= \frac{\mathcal E_0}{t_0}, \\
\widehat S_{\mathrm{heating},s}
&= \frac{S_{\mathrm{heating},s}}{S_{\mathrm{heating},0}},
\qquad s\in\{i,e\}.
\end{aligned}
$$

Here $\mathcal E_0$ is in J/m³, $S_{\mathrm{heating},0}$ is in W/m³, and
$\widehat S_{\mathrm{heating},s}$ is the normalized source used in the
corresponding energy equation. The solver performs this conversion after
reading the physical input; no manual normalization of the HDF5 power arrays
is needed.

The correspondence with code names is:

| Code name | Physical notation | Meaning |
| --- | --- | --- |
| `refval_mass` | $m_0$ | Reference mass in kg |
| `refval_density` | $n_0$ | Reference number density in m⁻³ |
| `refval_speed` | $u_0$ | Reference speed in m/s |
| `refval_time` | $t_0$ | Reference time in seconds |
| `refval_specenergydens` | $n_0u_0^2=\mathcal E_0/m_0$ | Energy-density scale per unit reference mass |
| `external_heating_ions`, `external_heating_electrons` | $\widehat S_{\mathrm{heating},i}$, $\widehat S_{\mathrm{heating},e}$ | Normalized heating on simulation nodes |

## Output

When external heating is enabled, solution files save the applied nodal
profiles under `/simulation_parameters/physics`:

| Dataset | Contents |
| --- | --- |
| `external_heating_ions` | Normalized ion heating $\widehat S_{\mathrm{heating},i}$ on the saved mesh |
| `external_heating_electrons` | Normalized electron heating $\widehat S_{\mathrm{heating},e}$ on the saved mesh |

Multiply these values by $S_{\mathrm{heating},0}$ to recover W/m³. The saved
simulation reference scales are under `/simulation_parameters/adimensionalization`;
use the scales belonging to that output when interpreting its normalized heating.

## Legacy inputs

### Stationary files

Legacy stationary files store the following datasets at the HDF5 root.
Shapes use Python/HDF5 notation; coordinate arrays are two-dimensional.

| Dataset | Shape | Physical quantity and units |
| --- | --- | --- |
| `r2D` | `(nr,nz)` | Cylindrical radius $R$ in metres |
| `z2D` | `(nr,nz)` | Vertical coordinate $Z$ in metres |
| `heating2D_ions` | `(nr,nz)` | Ion heating power density $S_{\mathrm{heating},i}$ in W/m³ |
| `heating2D_electrons` | `(nr,nz)` | Electron heating power density $S_{\mathrm{heating},e}$ in W/m³ |

Both channels are required, with zero arrays for inactive channels. Dimensions
are inferred from the arrays; legacy `ip`/`jp` datasets are not needed. A stored
`time`, if present, does not constrain stationary sampling. Full `.h5`
filenames and legacy extensionless names are accepted.

### Numbered files

For a legacy numbered sequence, supply a stem or the full `_0000.h5` anchor,
such as `heating_0000.h5`. Each file uses the stationary layout above and must
retain the same grid and channels.

Files are discovered from `0000` through the last contiguous number,
independently of the magnetic input and numerical step count. A lone
`_0000.h5` file supplies fixed heating and needs no timestamp.

For two or more files, every frame must contain a scalar `/time` in seconds.
These timestamps must be finite and strictly increasing. Magnetic timestamps
are not substituted for missing heating timestamps. Sampling remains linear;
the discovered heating history must cover the requested times.

## Implementation notes

[heating_series.f90](../../src/Models/heating_series.f90) reads dimensional
heating data without depending on magnetic equilibrium, solver globals or MPI.
[external_heating.f90](../../src/Models/external_heating.f90) selects the input,
samples it, applies normalization and interpolates onto the simulation nodes.
Its normalization factor is
`refval_time / (refval_specenergydens * refval_mass)`, equivalent to
$1/S_{\mathrm{heating},0}$.
Only returned sample copies are normalized; cached input frames stay dimensional.

Each `heating_series_t` contains a private `rectangular_series_t` from
[rectangular_series_io.f90](../../src/InOut/rectangular_series_io.f90).
The shared file/grid/time state and two-frame cache bookkeeping follow the
[equilibrium provider ownership rules](equilibrium.md#sampling-and-file-ownership).
Heating owns its two named power arrays and their interpolation. A cache slot
becomes valid only after its provider has successfully read the frame. The
shared series' `file()` accessor returns a borrowed handle; that shared object
owns closing it. Create/open providers independently rather than copying
objects with open handles.

Opening reads coordinate/time metadata and checks channel shapes. Power arrays
are loaded on demand into a two-frame cache. Numbered frames are checked for
matching layout, grid and required channels when first read. Keep inputs
unchanged while open; concurrent read-only access is supported.

| Provider method | Purpose |
| --- | --- |
| `open(path, ierr, message)` | Open a single-file trajectory or legacy stationary input |
| `open_numbered(path, last_frame, ierr, message)` | Open legacy frames 0000 through `last_frame`, reading their own timestamps |
| `sample_at(t, sample, ierr, message)` | Return an exact/interpolated dimensional sample at time `t` in seconds |
| `grid(r,z)`, `time_points(times)` | Return copies of coordinates/timestamps |
| `clear()` | Close the owned file and release metadata/cache; finalization also clears it |

Sample arrays `ions` and `electrons` have Fortran orientation `(nz,nr)`; the
underlying HDF5 dimension ordering is explained under
[array orientation](equilibrium.md#array-orientation-and-interface).
Samples are independent copies and record physical time, zero-based
`frame_left`/`frame_right` and interpolation `fraction`. At an exact frame,
the indices are equal and the fraction is zero.

The solver's heating loader takes an optional physical time in seconds.
Mesh reinitialization invalidates the nodal representation so it is rebuilt
from dimensional cached input at the same requested time.
Normal shutdown closes the provider. Reader/sampling errors return a status
and message; the solver reports input errors and terminates.
