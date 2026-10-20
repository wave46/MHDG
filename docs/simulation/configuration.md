# Configuring and preparing a simulation

This guide covers the existing parameter templates, case inputs and the
high-diffusion initialization workflow for the 2D `NGammaTiTeNeutral` model.
It is a starting point for preparing a case, not a catalogue of every model
or namelist entry. Installation and basic commands remain in the
[root README](../../README.md).

## Choose a starting configuration

The solver reads `param.txt` from its working directory. Copy an appropriate
template from `test/`, then set the paths and physical parameters for your
case before running:

| Template | Intended use | Starting state |
| --- | --- | --- |
| [param_initial.txt](../../test/param_initial.txt) | Short transient at high diffusion to obtain an initial plasma | Mesh, without a saved solution |
| [param_diffred.txt](../../test/param_diffred.txt) | Successive quasi-steady states with reduced diffusion | Saved plasma, using its embedded mesh |
| [param_steady.txt](../../test/param_steady.txt) | Steady solve at specified fixed diffusion | Saved plasma, using its embedded mesh |
| [param_reference.txt](../../test/param_reference.txt) | Steady example with mesh adaptivity and separate particle/heat diffusivities | Saved plasma, using its embedded mesh |

These templates omit `time_mode` and retain legacy stepping. The reference
template's adaptivity concerns the mesh. For a physical transient with a
specified final time, configure [physical time stepping](time-stepping.md);
the reference template also includes a commented example. When enabling an
explicit `time_mode` in one of these templates, also review `target_variable`.
The templates select recycling feedback (`2`), which requires a target-density
file. For constant puff and recycling, select `0` and leave `puff_path` empty;
see [controls](../inputs/controls.md).

The main blocks used in this guide are:

| Namelist | Contents |
| --- | --- |
| `SWITCH_LST` | Model switches, steady/pseudo-time selection, initialization and mesh-reading options |
| `INPUT_LST` | Equilibrium/current/heating paths, auxiliary namelists and output directory |
| `PHYS_LST` | Diffusivities, puff, recycling, effective charge and other physical coefficients |
| `TIME_LST` | Time units, integration scheme, intervals and stopping settings |
| `NUMER_LST` | Newton/time-convergence tolerances, damping and numerical boundary settings |
| `ADAPT_LST` | Mesh adaptivity and geometry path |
| `LSSOLV_LST` | Linear solver and preconditioner |
| `UTILS_LST` | Output frequency and balance diagnostics |

Paths are relative to the working directory unless absolute. Create the
output directory selected by `save_folder`; retain its trailing `/`.

## Prepare geometry and physical inputs

Prepare a Gmsh geometry describing the wall and the boundary regions used
for puffing, pumping and other boundary conditions, then generate the mesh.
Mesh and `.geo` files can be stored anywhere: supply the mesh path on the
command line and set `geometry_path` in `ADAPT_LST` to the corresponding
geometry when using mesh adaptivity.

Create `res/` inside the working directory used for the run:

```bash
mkdir -p res
```

The solver uses this directory for temporary mesh files during preprocessing
and adaptation, independently of the input locations and `save_folder`.
The regression harness creates it when preparing runs; create it yourself
for a standalone launch.

Mesh size is case dependent. Example sizes are roughly 1000 P8 elements as
an adaptive starting mesh for TCV, and 10 thousand P8 elements for WEST
without adaptivity, refined near targets and the far-SOL wall. These are
examples, not resolution criteria; assess the mesh for your geometry and
plasma gradients.

Configure input formats and scientific conventions using the authoritative
guides:

| Input | Configuration and conventions |
| --- | --- |
| Magnetic equilibrium and current | [Equilibrium input](../inputs/equilibrium.md): rectangular grids, fixed/moving histories, flux reconstruction and COCOS |
| Deposited ion/electron heating | [Heating input](../inputs/heating.md): power densities, independent time histories and normalization |
| Puff, effective charge and density feedback | [Controls](../inputs/controls.md): constant/scheduled choices, targets, chords and gains |
| Reduced 1D transport | [Transport model](../models/transport-1d.md): diffusion, pinch and topology settings |

Array sizes are inferred from the HDF5 datasets. Legacy dimension entries
remain accepted but do not set those sizes.

### Auxiliary model files

Select these in `INPUT_LST`:

| Path | Namelist | Used when |
| --- | --- | --- |
| `transport_model_path = 'transport_model.nml'` | `TRANSPORT_MODEL_1D_LST` | `transport_1d = .true.`; see the [transport guide](../models/transport-1d.md) and [sample](../../test/transport_model.nml) |
| `impurity_model_path = 'impurity_model.nml'` | `IMPURITY_RADIATION_LST` | Required with `impurity_radiation = .true.`; see the [sample](../../test/impurity_model.nml) |

Impurity names (`'N'`, `'W'`) and concentrations are arrays, with nonnegative
concentrations expressed as fractions of electron density. The concentration
schedule path also belongs in this auxiliary file; its selection and feedback
rules are in the [controls guide](../inputs/controls.md#configuration).

## Obtain an initial plasma and reduce diffusion

The templates use constant plasma diffusivities, with reduced 1D transport
and imported diffusion profiles off. This workflow uses three launches. Keep
equilibrium, sources and boundary conditions consistent between stages,
and retain each stage's parameter file with its output.

1. **Initialize at high diffusion.** Start from `param_initial.txt`, a mesh
   and no saved plasma. The template supplies 20 m²/s for `diff_n`, `diff_u`,
   `diff_e` and `diff_ee`, with `dt0 = 1e2` in normalized time units and a
   short legacy run. Check that the initial output is finite and physically
   admissible before continuing.

   From `test/`, for example:

   ```bash
   mkdir -p res
   cp param_initial.txt param.txt
   # Edit param.txt for your inputs, physical parameters and output directory.
   ./MHDG-NGammaTiTeNeutral-serial-2D /path/to/mesh.msh
   ```

2. **Reduce diffusion through quasi-steady states.** Use `param_diffred.txt`
   and a saved plasma. With `steady = .false.` and `psdtime = .true.`, the
   solver advances until the time-convergence criterion `tTM` is met, saves
   that stage, then multiplies the plasma diffusivities by `diffred`.
   The large pseudo-time intervals help approach steady states; they do not
   represent elapsed discharge time.

   ```bash
   cp param_diffred.txt param.txt
   # Edit param.txt for the same case and the intended reduction settings.
   ./MHDG-NGammaTiTeNeutral-serial-2D /path/to/mesh.msh /path/to/initial_solution.h5
   ```

3. **Reconverge at the desired fixed diffusion.** Select a converged stage
   near the desired diffusivity, use `param_steady.txt`, and set the exact
   coefficients in `PHYS_LST`. This performs the steady solve used as the
   starting point for subsequent scans or physical transients.

   ```bash
   cp param_steady.txt param.txt
   # Edit param.txt, including the exact desired diffusivities in PHYS_LST.
   ./MHDG-NGammaTiTeNeutral-serial-2D /path/to/mesh.msh /path/to/reduced_stage.h5
   ```

The restart templates use the saved solution's embedded mesh. The mesh
argument remains required; when `readMeshFromSol = .false.`, it must instead
identify the mesh matching the saved solution. See
[restart configuration](restart-and-output.md#starting-from-a-saved-solution).

For example, multiplying 20 m²/s repeatedly by `diffred = 0.6` gives a
stage near 0.56 m²/s, followed by about 0.34 m²/s. To obtain 0.5 m²/s,
reconverge the nearby stage with the desired fixed coefficients.

The legacy loop limit `nts` bounds numerical iterations, not the number of
completed quasi-steady stages. Set the final diffusivities with the fixed
steady solve above.

A converged reduction stage is saved before its coefficients are reduced.
The final save after the loop can follow another reduction without a solve
at the new coefficients. Select a converged stage output and perform the
fixed-diffusion solve above before interpreting it as the desired steady state.
Output names and saved metadata are described in
[restart and output](restart-and-output.md).

### Parameter scans and convergence

From the fixed-diffusion steady solution, vary quantities such as `puff`,
`Re` or the plasma diffusivities:

| Entry in `PHYS_LST` | Coefficient | Input units |
| --- | --- | --- |
| `diff_n` | Plasma-density diffusion | m²/s |
| `diff_u` | Parallel-momentum diffusion | m²/s |
| `diff_e` | Ion-energy diffusion | m²/s |
| `diff_ee` | Electron-energy diffusion | m²/s |

Use smaller parameter increments if Newton fails to reconverge. Damping
controls `dumpnr_min`, `dumpnr_max`, `dumpnr_width` and `dumpnr_n0` belong
in `NUMER_LST`; they set the Newton damping profile. Mesh changes and
physical-time retries have their own configured recovery rules.

## Linear-solver selection

The selected solver must be available in the executable's build. For PETSc
with a direct MUMPS factorization, set entries in `LSSOLV_LST`:

```fortran
  sollib = 3
  kspmethd = 'KSPPREONLY'
  pctype = 'PCLU'
  petsc_factor_solver = 'mumps'
```

An empty `petsc_factor_solver` retains PETSc's default factorization package.
An explicit selection requires `PCLU` and an available backend; an unavailable
backend stops the run. See the [build instructions](../../README.md#building)
for selecting the linked solver libraries.

## Optional neutral limiter and balance reporting

The neutral flux limiter is selected in `PHYS_LST`:

| Entry | Meaning |
| --- | --- |
| `neutral_flux_limiter_mode` | `'off'` (default) or `'lagged_flux_limiter'` |
| `neutral_flux_limiter_tn_source` | `'ti'` (default), using local ion temperature, or `'fixed'` |
| `neutral_flux_limiter_tn_eV` | Positive fixed neutral temperature in eV, required for `'fixed'` |
| `neutral_flux_limiter_eps` | Nonnegative flux regularization in m⁻² s⁻¹ |
| `neutral_flux_limiter_fs_flux_min` | Nonnegative minimum cap flux in m⁻² s⁻¹ |

`neutral_flux_limiter_save_2d = .true.` in `SWITCH_LST` saves nodal limiter
diagnostics and can be used with the limiter off.

`balance_diagnostics_mode` in `UTILS_LST` defaults to `'off'`:

| Mode | Report |
| --- | --- |
| `'summary'` | Compact contents, inputs/outputs and total balances |
| `'equations'` | Individual equation balances |
| `'detailed'` | Physical, discrete and boundary-condition components |

Enabled balance reporting requires a 2D model with plasma temperatures and
neutrals. Physical wall reporting supports Bohm, BohmPump and BohmPuff
boundaries; other physical boundaries produce a warning.
