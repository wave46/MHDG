# Parallel heat conduction

For ions and electrons, parallel conductive heat flows down the temperature
gradient. Magnetic topology selects the transport law: closed field lines
use classical or corrected conductivity; open field lines
use classical conductivity with optional free-streaming saturation.
Prescribed sheath energy loss is assembled separately.

## Closed field lines: conductivity corrections

`core_heat_conduction_mode` in `PHYS_LST` selects `off` (the default),
`neoclassical` or `mean_free_path`. With `off`, closed surfaces use classical
conductivity. Both active corrections use a singly charged hydrogenic
majority. The closed-surface heat flux is

```math
q_{\parallel,s}=-\kappa_{\parallel,s}\nabla_\parallel T_s,
\qquad \kappa_{\parallel,s}\propto T_s^{5/2}F_s.
```

### Neoclassical correction

The `neoclassical` correction accounts for the change in collisional
parallel transport as the collision length grows relative to the magnetic
connection length.
In a torus, this response sets poloidal temperature variations and the
resulting neoclassical radial heat transport. The fits connect the
Pfirsch–Schlüter, plateau and banana regimes.

[Rozhansky et al. (2009), section 3, equations (36)–(38)](https://pure.mpg.de/pubman/item/item_2141150_3/component/file_2141149/Rozhansky_New.pdf)
describe the ion parallel-conductivity prescription. Its collisionality
interpolation, and the electron interpolation, are based on
[Hinton and Hazeltine (1976), equations (6.125), (6.133) and Table III](https://doi.org/10.1103/RevModPhys.48.239).
The electron factor uses the normalized $K_{11}$ interpolation; the ion
factor uses the $K_2$ interpolation with finite-aspect-ratio dependence.
With the conventions below, the factors used here are:

```math
F_e=\frac{0.58\epsilon^{3/2}}{K_1},\qquad
K_1=\frac{1}{1+2.01\sqrt{\nu_{e}^{\ast}}+1.53\nu_{e}^{\ast}}
 +\frac{0.52\epsilon^3\nu_{e}^{\ast}}{1+0.89\epsilon^{3/2}\nu_{e}^{\ast}},
```

```math
F_i=\frac{1.60\epsilon^{3/2}}{K_2},\qquad
K_2=\frac{0.66+1.88\sqrt{\epsilon}-1.54\epsilon}
 {1+1.03\sqrt{\nu_{i}^{\ast}}+0.31\nu_{i}^{\ast}}
 +\frac{1.17\epsilon^3\nu_{i}^{\ast}}{1+0.74\epsilon^{3/2}\nu_{i}^{\ast}}.
```

Both use the collisionality

```math
\nu_{s}^{\ast}=\frac{L_\parallel}{2\pi\epsilon^{3/2}\tau_s\sqrt{T_{s,J}/m_s}}.
```

With temperatures expressed as energies in joules and singly charged ions,

```math
C=\frac{4\pi\epsilon_0}{e^2},\qquad
\tau_e=\frac{3C^2\sqrt{m_e}T_{e,J}^{3/2}}{4\sqrt{2\pi}n_e\ln\Lambda},\qquad
\tau_i^{\rm Br}=\frac{3C^2\sqrt{m_i}T_{i,J}^{3/2}}{4\sqrt{\pi}n_e\ln\Lambda}.
```

The implementation accepts temperatures in eV, density in m⁻³ and circuit
length in metres. The shared Coulomb logarithm is 15. Electron and ion
reference collision times follow Braginskii's electron–ion and ion–ion
conventions. The fits approach approximately unity at high collisionality.
The fitted values are used across the full collisionality range.

### Mean-free-path correction

`core_heat_conduction_mode='mean_free_path'` selects the algebraic
long-mean-free-path correction, using each species' temperature:

```math
F_s=\frac{1}{1+\lambda_{HC,s}/L_\parallel},\qquad
\lambda_{HC,s}[\mathrm{m}]=7.5\times10^{16}
\frac{T_s[\mathrm{eV}]^2}{n_e[\mathrm{m}^{-3}]}.
```

The factor approaches unity for a short collision length and reduces
conductivity when that length exceeds the magnetic circuit length. At fixed
density and geometry the strongly limited conductivity scales as $T_s^{1/2}$.
The kinetic motivation for heat-flux suppression is developed by
[Luciani, Mora and Virmont (1983), *Nonlocal Heat Transport Due to Steep Temperature Gradients*](https://doi.org/10.1103/PhysRevLett.51.1664),
and [Epperlein and Short (1991), *A practical nonlocal model for electron heat transport in laser plasmas*](https://doi.org/10.1063/1.859789).
Heat-carrying electrons sample temperatures over their collision length,
reducing the local collisional response when that length becomes
comparable to the parallel temperature-variation scale.

Here a local algebraic approximation uses $\lambda_{HC}/L_\parallel$,
with the full magnetic circuit length as its suppression scale. Applying
the same expression to ions with $T_s=T_i$ is the model's ion extension.

### Closed-surface geometry

Each axis-connected closed contour supplies a full poloidal circuit length
and an effective inverse aspect ratio:

```math
L_\parallel=\oint\frac{|B|}{|B_{\rm pol}|}\,dl_{\rm pol},\qquad
\epsilon=\frac{(\oint dl_{\rm pol})^2}{2\pi\oint R\,dl_{\rm pol}}.
```

Circular surfaces give the usual minor-to-major-radius ratio. At large
aspect ratio, the circuit length approaches the familiar expression

```math
\epsilon=\frac{r}{R_0},\qquad L_\parallel\simeq2\pi R_0q.
```

The profiles use the loaded equilibrium's contours and magnetic field.
The current implementation builds them from a rectangular equilibrium
grid, with either flux-derived or supplied poloidal field components.
Temperature-model assembly currently supports this correction in 2D with
`epn=2.5`.

At a regular magnetic axis, poloidal field strength and contour circumference
both vanish linearly with minor radius, leaving a finite circuit length.
The first two small contours determine its smooth continuation in normalized
poloidal flux:

```math
L_\parallel(\psi_N)=L_0+C\psi_N+O(\psi_N^2),\qquad \psi_N=\rho^2.
```

Near an ordinary X-point, the local flux is a quadratic saddle and the
circuit length grows logarithmically. The final two resolved contours
supply the offset and slope, both in metres:

```math
L_\parallel(\psi_N)\sim A+B\ln\frac{1}{1-\psi_N}.
```

These endpoint forms follow from local equilibrium expansions; fitting the
last two contours is the numerical continuation choice. See
[Webster (2009), section VI, equation (21)](https://doi.org/10.1063/1.3046070)
and the [MIT plasma-edge scaling report, p.4](https://library.psfc.mit.edu/catalog/reports/1990/95rr/95rr012/95rr012_full.pdf).
Limiter boundaries use a finite smooth continuation in normalized flux.

## Open field lines: free-streaming limiter

`flux_limiter=.true.` in `SWITCH_LST` enables local saturation on open field
lines, including private flux. The limiter fractions `c_fli` and `c_fle`
control the maximum ion and electron conductive heat flux:

```math
q_{\parallel,s}=
\frac{-\kappa_{\parallel,s}\nabla_\parallel T_s}
 {1+|\kappa_{\parallel,s}\nabla_\parallel T_s|/(c_{{\rm fl},s}q_{{\rm FS},s})},
\qquad q_{{\rm FS},s}\propto\frac{n_sT_s^{3/2}}{\sqrt{m_s}}.
```

This limits the classical conductive flux as it approaches a chosen fraction
of the thermal free-streaming flux. The ion temperature uses thermal energy
after subtracting bulk kinetic energy. The local temperature floor is
10⁻⁵ eV and the minimum free-streaming flux is 10⁻¹⁰ in solver units.

Closed surfaces use a free-streaming multiplier of one, with
either classical or corrected conductivity selected above.
The topology stays fixed during each Newton solve.

## Conductivity ceiling and Newton treatment

`T_fluxlim_maxi` and `T_fluxlim_maxe` specify the conductivity ceilings in eV:

```math
\kappa_{\parallel,s}=\kappa_{s,0}
\min\!\left(T_s^{5/2}F_s,T_{\max,s}^{5/2}\right).
```

The correction factor is one on open field lines and when the core
correction is disabled. MMS testcase 2 uses the uncapped conductivity law.
The specific-thermal-energy floor is 10⁻²⁰ in solver units.

Newton assembly differentiates the conductivity with respect to density
and temperature. On open field lines the current free-streaming multiplier
is updated at each assembly and held fixed during linearization. The
conductive flux uses the element state in volume terms and the trace state
on faces, together with the element gradient.

Thermal stabilization uses the corrected, capped conductivity before local
free-streaming saturation. See the [thermal stabilization description](../numerics/thermal-stabilization.md).
