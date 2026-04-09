# 1D Transport Module

This directory contains the reduced 1D transport-model path used by the
NGammaTiTe solver when `switch%transport_1d = .true.`.

## Purpose

The module does three things:

1. reduce the current 2D solution onto a flux-surface grid,
2. build 1D transport coefficients and pinch on that grid,
3. interpolate and mix those coefficients back into the HDG assembly.

The reduced data and the final model outputs are saved under the top-level
HDF5 group `/transport_1d`. This was chosen to keep them separate from both
the raw solution fields and the static input parameters.

## File Layout

- `flux_surface_transport_data.f90`
  Builds flux-surface-averaged profiles from the current solution.

- `transport_models_1d.f90`
  Parent module. Defines the main transport object, light lifecycle/config
  routines, and the public interfaces implemented in submodules.

- `transport_models_1d_config.f90`
  Scalar model controls grouped in `transport_model_config_t`.

- `transport_models_1d_derived.f90`
  Temporary derived data used only during one update pass.

- `transport_models_1d_update.f90`
  Top-level orchestration of one 1D transport update from flux-surface data.

- `transport_models_1d_diffusion_bohm_gyrobohm.f90`
  Bohm / gyro-Bohm diffusion model and related derived quantities.

- `transport_models_1d_pinch.f90`
  Pinch model definitions and pinch windowing.

- `transport_models_1d_runtime.f90`
  Solver-facing runtime use: interpolation, diffusion mixing, pinch matrix.

- `transport_models_1d_common.f90`
  Shared interpolation and analytic window helpers.

- `transport_models_1d_io.f90`
  HDF5 writing of the reduced profiles, model coefficients, and parameters.

## Update Flow

The 1D model is refreshed once per NR iteration before Jacobian assembly.

High-level sequence:

1. `fs_transport%build_profiles()`
2. `transport_model_1d%update_from_flux_surfaces(fs_transport)`
3. HDG uses:
   - `transport_model_1d%apply_1D_diffusion(...)`
   - `transport_model_1d%compute_1D_pinch_matrix(...)`

Inside `update_from_flux_surfaces(...)`:

1. copy the current flux-surface grid,
2. fill temporary derived profiles from `fs_data`,
3. build projected gradients and `delta_te`,
4. compute Bohm / gyro-Bohm transport,
5. compute pinch.

## Persistent vs Derived State

Persistent on `transport_model_1d_t`:

- `config`
- `rho_grid`
- final coefficients actually used by the solver:
  - `chi_i_fs`
  - `chi_e_fs`
  - `d_fs`
  - `nu_mom_fs`
  - `vpinch_fs`

Temporary on `transport_model_derived_t`:

- thermodynamic profiles,
- geometry / magnetic profiles,
- projected gradients,
- collisionality,
- intermediate Bohm / gyro-Bohm quantities,
- `a_minor`,
- `delta_te`.

This separation is intentional: the main transport object keeps the final
solver-facing state, while the derived helper holds one-update intermediate
data only.

## HDF5 Layout

The solver writes:

- `/transport_1d/profiles`
  - `rho_grid`
  - `shell_weight`
  - `U_fs`
  - `Q_rad_fs`

- `/transport_1d/coefficients`
  - `chi_i_fs`
  - `chi_e_fs`
  - `d_fs`
  - `nu_mom_fs`
  - `vpinch_fs`

- `/transport_1d/params`
  - scalar transport-model controls from `transport_model_config_t`

Important for postprocessing:

- coefficient arrays use `/transport_1d/profiles/rho_grid` as their radial axis
- do not reconstruct a synthetic `linspace(0,1,nrho)` grid

## Extending the Module

To add a new diffusion model:

1. keep the reduced-profile build unchanged,
2. add the model logic in a dedicated submodule,
3. write final outputs into the persistent coefficient arrays on
   `transport_model_1d_t`.

To add a new pinch model:

1. implement it in `transport_models_1d_pinch.f90`,
2. select it in `tm1d_compute_pinch_profile`,
3. keep the runtime use in `transport_models_1d_runtime.f90` unchanged.
