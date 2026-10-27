# Parallel heat-flux limiter

For ions and electrons, parallel conduction follows the Spitzer–Härm law,
with a limiter that prevents the conductive flux from exceeding a chosen
fraction of the free-streaming flux:

$$
q_{\parallel,s}=\phi_s q_{\mathrm{SH},s},\qquad
q_{\mathrm{SH},s}=-\kappa_{\parallel,s}\nabla_\parallel T_s,\qquad
\phi_s=\frac{1}{1+|q_{\mathrm{SH},s}|/q_{\max,s}}.
$$

Here $s$ denotes ions or electrons, and $q_{\max,s}$ is the allowed flux.
The free-streaming flux scales as
$nT_s^{3/2}$, with the appropriate species mass. The dimensionless inputs
`c_fli` and `c_fle` set its allowed fraction for each species, and
`flux_limiter` enables both limiters. The prescribed sheath energy loss is
assembled separately.

## Temperature limits

The parallel conductivity grows as $T_s^p$, normally with $p=5/2$.
Ion temperature uses thermal energy after subtracting bulk kinetic energy.
The inputs `T_fluxlim_maxi` and `T_fluxlim_maxe`, specified in **eV**, set
separate upper temperatures for the ion and electron conductivities. Above
these temperatures the conductivity stays at its capped value. The
manufactured solution (MMS) testcase 2 bypasses these upper limits.

The conductivity uses a numerical temperature floor barely above zero:

$$
T_{\kappa,\min}\approx2.67\times10^{-20}\ \mathrm{eV}.
$$

With the code's fixed reference scales, the reference ion energy is about
4.00 eV, giving this floor from $\tfrac23\,10^{-20}E_{\rm ref}$.
This is only a numerical positivity check, with no meaningful physical
temperature scale. It keeps the conductivity positive without changing
the energy or temperature stored in the solution.

The **free-streaming cap** uses a separate, fixed temperature floor of
$10^{-5}$ eV. Below it the cap still varies linearly with density. A further
minimum free-streaming flux of $10^{-10}$ in solver flux units prevents a
zero denominator when evaluating the limiter.

## Linearization

Inside the conductivity's temperature limits, its derivative follows the
power law. Below the lower floor or above the upper limit, the conductivity
is constant and its derivative is zero. Exactly at the lower limit the
constant-side derivative is used; exactly at the upper limit the power-law
side is used. Multiplying all conservative state components by the same
factor leaves temperature and conductivity unchanged.

For the currently assembled exact heat-limiter derivative, write
$h=\phi q_{\mathrm{SH}}$. Its response to changes in the classical flux and
the free-streaming flux is

$$
\frac{\partial h}{\partial q_{\mathrm{SH}}}=\phi^2,\qquad
\frac{\partial h}{\partial q_{\mathrm{fs}}}=
\begin{cases}
h(1-\phi)/q_{\mathrm{fs}}, & |q_{\mathrm{fs}}|>10^{-10},\\
0, & |q_{\mathrm{fs}}|\le10^{-10}.
\end{cases}
$$

The linearized expression includes a constant term so that it reproduces
the actual limited flux at the state where it was evaluated. This remains
true when a temperature or flux floor is active. The same coefficient
rules are used in volume and face assembly; face conductivity uses trace
temperature with element gradients.
