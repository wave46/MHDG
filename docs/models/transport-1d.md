# Reduced 1D transport

The reduced transport model supplies particle, parallel-momentum, ion-energy
and electron-energy diffusion coefficients, with an optional pinch, for the
NGammaTiTe family of 2D simulations. It reduces the current plasma solution
onto a flux-surface grid, evaluates transport there, and interpolates the
coefficients back into the 2D equations. The plasma is still evolved by the
2D solver.

## Configuration

In the existing namelist blocks of `param.txt`, enable the model and select
its settings file:

```fortran
&SWITCH_LST
  transport_1d = .true.
/
&INPUT_LST
  transport_model_path = 'transport_model.nml'
/
```

Copy [test/transport_model.nml](../../test/transport_model.nml) to your run
directory as a starting point. Its `TRANSPORT_MODEL_1D_LST` namelist controls
diffusion, pinch and the regions where they apply, as described below. Set
`transport_model_path` to the file you intend to use; relative paths are
resolved from the solver's working directory. The default filename is
`transport_model.nml`.

### Diffusion

Ion and electron heat diffusion combine Bohm and gyro-Bohm contributions,
as described in the [NTCC/JETTO mixed-model reference](https://w3.pppl.gov/ntcc/JETTO/mixed_Bohm_gyro_Bohm/node4.html).
Particle and momentum diffusion are then derived from them. The
[diffusion expressions](#diffusion-model) below define these terms and the
electron-temperature drop used in the Bohm contribution.

| Settings | Purpose |
| --- | --- |
| `c_bohm_i`, `c_gyrobohm_i`, `c_bohm_e`, `c_gyrobohm_e` | Ion/electron Bohm and gyro-Bohm weights |
| `rho_core`, `rho_edge` | Reference radii defining the relative electron-temperature drop, $\Delta T_e$ |
| `c_bohm_n`, `c_bohm_n_rho_slope` | Particle-diffusion strength and radial taper |
| `prandtl` | Ratio of momentum diffusivity to modeled ion heat diffusivity |
| `diff_n_min_phys`, `diff_u_min_phys`, `diff_e_min_phys`, `diff_ee_min_phys` | Particle, momentum, ion-energy and electron-energy diffusion floors in m²/s |

The [normalization convention](#coordinates-and-normalization) below matters
when choosing the heat-diffusion weights from literature values.

### Pinch model and equation selection

Pinch adds radial convection to the selected plasma equations. Choose its
radial dependence with `pinch_model`; the [pinch expressions](#pinch-models)
below give the forms used in MHDG.

| `pinch_model` | Model | Strength controls | Reference or expression |
| --- | --- | --- | --- |
| `1` | Collisionality-dependent pinch; built-in default | `c_pinch`, `nu_th` | [Militello Asp et al. (2022), §2][militello] |
| `2` | Geometric / Polevoi-style baseline | `c_pinch` | [Polevoi et al. (2019), related geometric form][polevoi] |
| `3` | Constant speed inside the pinch window | `vpinch_const_phys` in m/s | [Prescribed constant velocity](#pinch-models) |

The supplied sample selects model `2` explicitly. If `pinch_model` is omitted
from the namelist, the solver uses model `1`.

Positive pinch is outward; negative pinch is inward. For models 1 and 2,
choose the sign through `c_pinch`; for model 3, use the sign of
`vpinch_const_phys`.

`pinch_equations` accepts up to four entries:

| ID | Conserved equation |
| --- | --- |
| `1` | Particles |
| `2` | Parallel momentum |
| `3` | Total ion energy |
| `4` | Electron energy |

Omission selects particles only. Examples are `pinch_equations = 1,2,3,4`
for all plasma rows, `3` for ion energy only, and `0` to disable pinch. Zero
entries are unused slots; repeated IDs select an equation once. IDs outside
0–4 and overlong lists are rejected. Each selected equation receives the
same velocity, sign, topology policy and radial window. Neutral rows remain
zero.

### Regions and radial windows

The radial coordinate is the normalized poloidal radius `rho_pol_norm`: zero
at the magnetic axis and one at the last closed flux surface (LCFS). It is
independent of the topology selection, so equal-radius points can belong to
different, disconnected regions.

| `transport_region_policy` | Included regions |
| --- | --- |
| `'legacy_all_regions'` | All regions; compatibility default |
| `'core_and_main_sol'` | Axis-connected core and main scrape-off layer (SOL) |
| `'core_only'` | Axis-connected core |

Excluded regions retain parameter-file diffusion and receive zero 1D pinch.
The two topology-restricted policies require the magnetic geometry cache.
The main-SOL profile may extend beyond `rho_pol_norm = 1`; it is not clipped
at the LCFS. To apply modeled transport there, the corresponding radial cutoff
must also extend beyond one.

Modeled diffusion has full weight below
`rho_diffusion_model_max - rho_blend_width`, then blends smoothly into the
parameter-file coefficients. At and above `rho_diffusion_model_max`, the
parameter-file coefficients are retained. The configured diffusion floors
apply to the modeled coefficients before blending.

Pinch has its own window: it vanishes on axis, rises to full strength over
`rho_pinch_axis_width`, and falls to zero over `rho_pinch_edge_width` before
`rho_pinch_model_max`. It is zero at and beyond that outer cutoff. Both windows
use cubic smoothstep ramps.

## Model expressions

The expressions use flux-surface profiles. Here $p_e$ and $T_e$ are electron
pressure and temperature, $q$ is the safety factor, $\rho_s$ is the sound
gyroradius, $c_s$ is the electron-temperature sound speed, and $a$ is the
minor radius. The normalized poloidal radius is $\rho_{\mathrm{pol,n}}$.
Small denominator bounds are denoted by $\varepsilon$. The
[variable table](#coordinates-and-normalization) maps this notation to input
and code names.

### Diffusion model

The mixed Bohm / gyro-Bohm model is:

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
\max\!\left(1-c_{\mathrm{B},n,\rho}\rho_{\mathrm{pol,n}},0\right)\,
\frac{\chi_i\chi_e}{\max(\chi_i+\chi_e,\varepsilon)},
\qquad
\nu_{\mathrm{mom}} = \mathrm{Pr}\,\chi_i.
$$

The particle-diffusion taper decreases linearly with normalized poloidal
radius and is bounded below by zero. Its default slope is
`c_bohm_n_rho_slope = 0.7`; setting the slope to zero removes the taper.
The same expression applies to selected main-SOL profiles beyond the LCFS.

### Pinch models

The model forms below use MHDG's normalized poloidal radius and outward-positive
velocity convention. The references in the selection table describe related
particle-pinch models; they do not define MHDG's equation selection or sign.

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

## Coordinates and normalization

The pinch expressions use the radial factor $\rho_{\mathrm{pol,n}}/a$.
The inverse aspect ratio entering the electron collisionality is

$$
\epsilon_{\mathrm{asp}}
= \frac{a\,\rho_{\mathrm{pol,n}}}{R_{\mathrm{maj}}}.
$$

The correspondence between code names and physical notation is:

| Input or code name | Symbol | Meaning |
| --- | --- | --- |
| `rho_pol_norm`, `rho_grid` | $\rho_{\mathrm{pol,n}}$ | Normalized poloidal radius; `rho_grid` stores its profile samples |
| `rho_core`, `rho_edge` | $\rho_{\mathrm{core}}$, $\rho_{\mathrm{edge}}$ | Reference radii for the relative electron-temperature drop |
| `a_minor` | $a$, $a_{\mathrm{minor}}$ | Minor radius |
| `Rmaj` | $R_{\mathrm{maj}}$ | Major radius |
| `q_fs` | $q$ | Safety-factor profile |
| `eps_fs` | $\epsilon_{\mathrm{asp}}$ | Inverse-aspect-ratio profile |
| `nuestar_fs`, `nu_th` | $\nu_*$, $\nu_{\mathrm{th}}$ | Normalized electron collisionality and its suppression threshold |
| `c_bohm_n_rho_slope` | $c_{\mathrm{B},n,\rho}$ | Particle-diffusion taper slope |
| `prandtl` | $\mathrm{Pr}$ | Momentum-to-ion-heat diffusivity ratio |
| `c_pinch` | $c_{\mathrm{pinch}}$ | Signed pinch coefficient |
| `vpinch_const_phys` | $V_{\mathrm{const}}$ | Prescribed physical pinch speed in m/s |
| `refval_length`, `refval_time` | $L_0$, $t_0$ | Reference length and time for solver normalization |

With magnetic equilibrium geometry available, transport uses the radial
coordinate, topology region and increasing normalized-flux normal evaluated
at each volume or face quadrature point. The pinch direction follows that
normal independently of magnetic-field sign.

In MHDG's formulation, the stored ion and electron heat-diffusion coefficients
$\chi_i$ and $\chi_e$ correspond to **2/3 of the usual heat-conductivity
convention** used in transport literature. Account for this factor when
choosing the Bohm / gyro-Bohm weights.

The four `diff_*_min_phys` inputs are in m²/s, and `vpinch_const_phys` is in
m/s. They are converted to solver units when applying the model configuration,
using factors $t_0/L_0^2$ for diffusion and $t_0/L_0$ for velocity. Radial
coordinates, window widths and model weights are dimensionless.

The mixed ion and electron diffusion profiles have a numerical lower bound
of `1.e-10` in solver units before particle and momentum diffusion are formed.
The physical diffusion floors specified in the namelist are applied later,
when coefficients are interpolated and blended into the 2D equations. Pinch
models 1 and 2 use the modeled `d_fs` before those assembly floors.

## Output and interpretation

When transport is enabled and its profiles have been initialized, solution
HDF5 files contain the top-level group `/transport_1d`:

| Group | Contents |
| --- | --- |
| `/transport_1d/profiles` | `rho_grid`, `shell_weight`, conservative profiles `U_fs`, and radial-gradient profiles `Q_rad_fs` |
| `/transport_1d/coefficients` | `chi_i_fs`, `chi_e_fs`, `d_fs`, `nu_mom_fs`, `vpinch_fs` |
| `/transport_1d/params` | Model settings, including normalized floors/speed, region policy and equation selection |

Plot the coefficients against the saved radial coordinate
`/transport_1d/profiles/rho_grid`, which can extend beyond the LCFS.

Saved coefficients are the reduced model profiles in solver units. Multiply
diffusion coefficients by $L_0^2/t_0$ to obtain m²/s, and `vpinch_fs` by $L_0/t_0$
to obtain m/s. These profiles precede the local topology masks, diffusion
floors, radial windows and blending used in assembly; they are not a map of
the effective coefficients at every mesh point.

The effective pinch selection is saved as a sorted, zero-padded integer array
at `/transport_1d/params/pinch_equations`.

## References

- Mixed Bohm / gyro-Bohm online description from the NTCC/JETTO documentation:
  https://w3.pppl.gov/ntcc/JETTO/mixed_Bohm_gyro_Bohm/

- Final mixed model summary and coefficients:
  https://w3.pppl.gov/ntcc/JETTO/mixed_Bohm_gyro_Bohm/node4.html

- Validation paper for the mixed Bohm / gyro-Bohm model:
  https://scipub.euro-fusion.org/archives/jet-archive/validation-of-a-new-mixed-bohmgyro-bohm-transport-model-on-discharges-of-the-iter-data-base

- [E. Militello Asp et al., *JINTRAC integrated simulations of ITER scenarios
  including fuelling and divertor power flux control for H, He and DT plasmas*,
  Nuclear Fusion 62 (2022) 126033][militello], §2: collisionality-dependent
  particle pinch.

- [A. R. Polevoi et al., *Impact of Suprathermal Ions on Neutron Yield at
  Pre-DT Phase of ITER Operation* (2019)][polevoi], “1.5D transport simulations”:
  an example of a geometric particle-pinch form.

[militello]: https://scientific-publications.ukaea.uk/wp-content/uploads/MILITELLO_ASP_2022_NUCL-_FUSION_62_126033.PDF
[polevoi]: https://conferences.iaea.org/event/185/contributions/14942/attachments/8096/10672/ID7-POLEVOI.pdf

## Implementation notes

### Update and state

The model is refreshed once per Newton–Raphson iteration before Jacobian
assembly. The sequence in
[Main_utils.f90](../../src/Utils/Main_utils.f90) is:

1. `fs_transport%build_profiles(region_policy)` reduces the current solution.
2. `transport_model_1d%update_from_flux_surfaces(fs_transport)` copies the
   radial grid, builds temporary thermodynamic/geometry profiles, projected
   gradients, electron-temperature drop and collisionality, then computes
   diffusion and pinch.
3. HDG assembly uses `apply_1D_diffusion(...)` and
   `compute_1D_pinch_matrix(...)`.

The transport object retains configuration, `rho_grid` and the five output
coefficient profiles. Thermodynamic profiles, geometry, gradients,
collisionality, `a_minor`, `delta_te` and intermediate Bohm / gyro-Bohm values
belong to a temporary `transport_model_derived_t` used during one update.

Geometry-cache values are evaluated at quadrature points; the model does not
interpolate nodal radius or take the square root of interpolated nodal flux
when that cache exists. The compatibility path without a cache uses nodal
flux. If the supplied normal vanishes, `legacy_all_regions` alone permits
the legacy magnetic-field direction fallback; the restricted policies suppress
pinch at such points.

### Input defaults and checks

The input reader starts from the declared transport defaults, then applies
the entries in `TRANSPORT_MODEL_1D_LST`. Values not listed in that namelist
keep their defaults; the supplied sample makes its choices explicit and
differs from the defaults in several settings.

| Input condition | Reader behavior |
| --- | --- |
| A setting is omitted | Keep that setting's default |
| `transport_model_path` is empty | Use defaults without reading a file |
| The selected file cannot be opened | Print a warning and use defaults |
| The namelist cannot be parsed | Stop the run with an input error |

### Source map

All files below are in
[src/Models/NGammaTiTe/transport_1d/](../../src/Models/NGammaTiTe/transport_1d/):

| File | Responsibility |
| --- | --- |
| `flux_surface_transport_data.f90` | Reduce the current solution to flux-surface profiles |
| `transport_models_1d.f90` | Parent module, public interfaces and lifecycle |
| `transport_models_1d_config.f90` | Configuration, normalization, topology policy and pinch direction |
| `transport_models_1d_derived.f90` | Temporary data for one update |
| `transport_models_1d_update.f90` | Orchestrate the profile/model update |
| `transport_models_1d_diffusion_bohm_gyrobohm.f90` | Bohm / gyro-Bohm diffusion, collisionality and projected gradients |
| `transport_models_1d_pinch.f90` | Pinch models and their radial window |
| `transport_models_1d_runtime.f90` | Interpolation, diffusion mixing and pinch matrix |
| `transport_models_1d_common.f90` | Profile interpolation and analytic window helpers |
| `transport_models_1d_io.f90` | Save coefficients and settings to HDF5 |

### Extending the model

For a new diffusion model, keep the reduced-profile build unchanged, implement
the model in a dedicated submodule, and write its outputs into the persistent
coefficient arrays on `transport_model_1d_t`.

For a new pinch model, implement it in `transport_models_1d_pinch.f90` and
select it in `tm1d_compute_pinch_profile`. Keep its runtime application in
`transport_models_1d_runtime.f90` unchanged.
