# 1D Transport Module

This directory contains the reduced 1D transport-model path used by the
NGammaTiTe solver when `switch%transport_1d = .true.`.

The user-facing settings for this module are read from
`test/transport_model.nml`. That namelist controls:

- the reference radii used by the Bohm / gyro-Bohm model,
- the topology regions in which reduced transport is applied,
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

1. `fs_transport%build_profiles(region_policy)`
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

## Radial coordinate and topology policy

When magnetic equilibrium geometry is available, reduction and HDG assembly
read `rho_pol_norm`, the increasing-flux normal, and the topology region at the
exact volume or face quadrature point from `magnetic_geometry_cache`.  They do
not interpolate nodal `rho` and do not take a square root of interpolated nodal
psi.  The legacy nodal-psi path remains only as a compatibility fallback when
no equilibrium cache exists.

`transport_region_policy` is independent of the radial coordinate:

- `legacy_all_regions` applies reduced transport everywhere and is the
  compatibility default;
- `core_and_main_sol` includes the axis-connected core and main SOL but leaves
  private flux and undefined regions on parameter-file diffusion with zero 1D
  pinch;
- `core_only` applies reduced transport only in the axis-connected core.

The policy therefore distinguishes equal-rho points in disconnected topology
regions.  The main-SOL profile is not clipped at the LCFS and may extend through
`rho_pol_norm > 1`.

The pinch vector uses the cached increasing-flux normal, so its orientation is
independent of magnetic-field sign.  Its radial scale is `rho_pol_norm/a_minor`;
collisionality uses `a_minor*rho_pol_norm/Rmaj`.

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

$$
\chi_i = c_{\mathrm{B},i}\,\chi_{\mathrm{B}} + c_{\mathrm{gB},i}\,\chi_{\mathrm{gB}},
\qquad
\chi_e = c_{\mathrm{B},e}\,\chi_{\mathrm{B}} + c_{\mathrm{gB},e}\,\chi_{\mathrm{gB}}.
$$

with

$$
\chi_{\mathrm{B}}
= \rho_s c_s q^2 a\,
\frac{\left|\partial_r p_e\right|}{\max(p_e,\varepsilon)}\,
\Delta T_e,
$$

$$
\chi_{\mathrm{gB}}
= \rho_s^2 c_s\,
\frac{\left|\partial_r T_e\right|}{\max(T_e,\varepsilon)},
$$

and

$$
\Delta T_e
= \frac{T_e(\rho_{\mathrm{core}})-T_e(\rho_{\mathrm{edge}})}
        {\max\!\left(T_e(\rho_{\mathrm{edge}}),\varepsilon\right)}.
$$

The particle and momentum transport are then built as

$$
d_{\mathrm{fs}}
= c_{\mathrm{B},n}\,
\frac{\chi_i\chi_e}{\max(\chi_i+\chi_e,\varepsilon)},
\qquad
\nu_{\mathrm{mom}} = \mathrm{Pr}\,\chi_i.
$$

Important normalization note:

- in the HDG formulation used here, the stored ion and electron `\chi`
  coefficients correspond to `2/3` of the more typical heat-conductivity
  convention used in transport literature
- when choosing `\chi`-related coefficients in `transport_model.nml`, this
  factor should be taken into account

The pinch models are:

1. Militello-style collisionality suppression

$$
V_{\mathrm{pinch}}
= \min\!\left(1,\exp\!\left[1-\frac{\nu_*}{\max(\nu_{\mathrm{th}},\varepsilon)}\right]\right)
\; c_{\mathrm{pinch}}\, d_{\mathrm{fs}}\,
\frac{\rho_{\mathrm{pol,n}}}{\max(a_{\mathrm{minor}},\varepsilon)}
$$

2. Geometric / Polevoi-style baseline

$$
V_{\mathrm{pinch}}
= c_{\mathrm{pinch}}\, d_{\mathrm{fs}}\,
\frac{\rho_{\mathrm{pol,n}}}{\max(a_{\mathrm{minor}},\varepsilon)}
$$

3. Constant pinch

$$
V_{\mathrm{pinch}} = V_{\mathrm{const}}.
$$

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
