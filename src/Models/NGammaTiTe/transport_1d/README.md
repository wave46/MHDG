# 1D Transport Module

This directory contains the reduced 1D transport-model path used by the
NGammaTiTe solver when `switch%transport_1d = .true.`.

The user-facing settings for this module are read from
`test/transport_model.nml`. That namelist controls:

- the reference radii used by the Bohm / gyro-Bohm model,
- the diffusion replacement window,
- the pinch model and pinch windows,
- the coefficient floors.

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
  Parent module. Defines the main transport object, keeps the light lifecycle
  routines, and declares the public interfaces implemented in submodules.

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

## Model Expressions

The current implementation uses the mixed Bohm / gyro-Bohm structure:

```text
chi_i = c_B,i  * chi_B  + c_gB,i * chi_gB
chi_e = c_B,e  * chi_B  + c_gB,e * chi_gB
```

with

```text
chi_B  = rho_s * c_s * q^2 * a * abs(d p_e / dr) / max(p_e, eps) * Delta_Te
chi_gB = rho_s^2 * c_s * abs(d T_e / dr) / max(T_e, eps)
```

and

```text
Delta_Te = ( T_e(rho_core) - T_e(rho_edge) ) / max( T_e(rho_edge), eps )
```

The particle and momentum transport are then built as

```text
d_fs   = c_B,n * chi_i * chi_e / max(chi_i + chi_e, eps)
nu_mom = Pr * chi_i
```

The pinch models are:

1. Militello-style collisionality suppression

```text
V_pinch = min( 1, exp( 1 - nu_star / max(nu_th, eps) ) )
          * c_pinch * d_fs * r / max(a, eps)^2
```

2. Geometric / Polevoi-style baseline

```text
V_pinch = c_pinch * d_fs * r / max(a, eps)^2
```

3. Constant pinch

```text
V_pinch = V_const
```

In this code the sign convention is radial:

- positive pinch is outward
- negative pinch is inward

The pinch is additionally windowed near the axis and near the outer cutoff
through analytic smoothstep-based ramps in `rho_pol_norm`.

## References

- Mixed Bohm / gyro-Bohm online description from the NTCC/JETTO documentation:
  https://w3.pppl.gov/ntcc/JETTO/mixed_Bohm_gyro_Bohm/

- Final mixed model summary and coefficients:
  https://w3.pppl.gov/ntcc/JETTO/mixed_Bohm_gyro_Bohm/node4.html

- Validation paper for the mixed Bohm / gyro-Bohm model:
  https://scipub.euro-fusion.org/archives/jet-archive/validation-of-a-new-mixed-bohmgyro-bohm-transport-model-on-discharges-of-the-iter-data-base

- Example integrated-model reference used for the Polevoi-style pinch baseline:
  https://doi.org/10.1007/s10894-020-00232-x

- Example ITER/JINTRAC scenario work by E. Militello Asp and collaborators:
  https://nucleus.iaea.org/sites/fusionportal/Shared%20Documents/FEC%202020/fec2020-preprints/preprint1104.pdf

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
