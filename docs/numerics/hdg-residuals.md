# Mixed HDG residuals and Newton convergence

In `NUMER_LST`, `newton_convergence` selects `update` (default, existing `tNR`
test), `residual`, or `both`. Residual acceptance requires

$$ r_k\leq\varepsilon_{\mathrm{abs}}+\varepsilon_{\mathrm{rel}}r_0, $$

where `residual_atol` and `residual_rtol` set the two tolerances. Their defaults,
`1.d-6` and `1.d-4`, are provisional. `both` additionally requires the preceding
undamped volume update to satisfy `tNR`, with at least one update on this mesh.
Residual acceptance also requires the existing raw-state admissibility check;
this is separate from the norm and does not certify positivity at every Gauss point.

## Evaluation and output

The [solver-flow overview](../simulation/solver-flow.md) places these checks
within initialization, time stepping and the HDG solve.

Let the current complete state be $X_k=(U_k,Q_k,\widehat U_k)$. The MHDG
Newton loop follows this sequence:

1. Assemble matrices and right-hand sides at $X_k$.
2. Evaluate its nonlinear residual $F(X_k)$ and normalized magnitude $r_k$.
3. If the selected residual criterion and admissibility pass, accept $X_k$.
4. Otherwise solve for a correction and apply the damping factor $\lambda$:

   $$ J_k\Delta X=-F(X_k),\qquad X_{k+1}=X_k+\lambda\Delta X. $$

The next iteration repeats steps 1–2 at $X_{k+1}$ to evaluate $F(X_{k+1})$.
With `update` stopping, the solver instead accepts $X_{k+1}$ when its update
passes `tNR`; the printed $r_k$ still belongs to $X_k$. With `residual` or
`both`, it accepts $X_k$ in step 3, so the printed norm belongs to the saved state.
The reference $r_0$ resets for each nonlinear solve/retry or mesh-matrix rebuild.
The manufactured-solution `Convergence` driver monitors residuals but retains
its existing update-based stopping test.

`residual_diagnostics_mode='summary'` in `UTILS_LST` prints two lines:

```text
HDG residual before update: norm=  3.007E-07; relative=  1.000E+00; target=  1.000E-06
  local peaks (conservation, gradient, trace):   5.096E-05   1.613E-11   3.776E-05
```

`norm` is $r_k$, `relative` is $r_k/r_0$, and `target` is the acceptance bound
above. Each peak is the largest local norm over all equations in that block;
peaks are diagnostic, not additional stopping tests. `equations`
also prints per-equation columns: `U RMS`, `Q RMS`, `trace RMS`, `U peak`,
`Q peak`, `trace peak`. Here `U` means conservation and `Q` means the gradient relation.
`dotiming` adds the maximum-rank residual evaluation time in seconds.

The residual monitor owns `residual_diagnostics_mode`: `off` (default),
`summary`, or `equations`. Balance diagnostics has its own independent mode.
With residual reporting off, `update` skips residual evaluation; `residual` and
`both` still evaluate silently. When $r_0=0$, the absolute tolerance governs; the relative
value is zero if $r_k=0$ and a large sentinel otherwise.

## Residual equations and normalization

For element state $U$, gradient $Q$, and trace $\widehat U$, the assembled blocks give

$$
\begin{aligned}
R_U^e &= A_{uu}U+A_{uq}Q+A_{ul}\widehat U_{\mathrm{local}}+A_{ul,\mathrm{dir}}-S,\\
R_Q^e &= Q+M_Q^{-1}(A_{qu}U+A_{ql}\widehat U_{\mathrm{local}}+A_{ql,\mathrm{dir}}),\\
R_{\widehat U}^e &= A_{lu}U+A_{lq}Q+A_{ll}\widehat U_{\mathrm{global}}-f_h.
\end{aligned}
$$

$S$ and $f_h$ include boundary data and affine linearization corrections, as well
as physical terms. Transient conservation includes the assembled BDF/history
contribution. Reconstruction relies on consistent affine assembly; it does not
independently reevaluate every physical closure. `iAqq` already stores $M_Q^{-1}$.

$M_e,M_f$ are scalar element/face masses and $V_e,S_f$ their measures, using
curved geometry and the cylindrical radius weight in axisymmetry. Equation $j$
uses fixed positive reference $s_j$, set by `residual_scales(j)` (default one).
All references are in solver units; local density or temperature is never a divisor.
For direction $d$, the local squared norms are

$$
\begin{aligned}
\ell_{U,j,e}^2 &= \frac{(R_{U,j}^e)^T M_e^{-1}R_{U,j}^e}{V_e s_j^2},\\
\ell_{Q,j,e}^2 &= \frac{\sum_d(R_{Q,j,d}^e)^T M_eR_{Q,j,d}^e}{N_{\mathrm{dim}}V_e s_j^2},\\
\ell_{\widehat U,j,f}^2 &= \frac{(R_{\widehat U,j}^f)^T M_f^{-1}R_{\widehat U,j}^f}{S_f s_j^2}.
\end{aligned}
$$

Each residual here is a vector over local nodal test functions. Conservation and
trace use inverse mass; the already normalized gradient uses mass. Global block/
equation RMS values are measure-weighted averages of these local squares. The
scalar $r_k$ is their RMS with equal weight per active block and equation; an
empty block is omitted. Norms measure equation mismatch in fixed model reference
units, not solution error. Reference length and time are one in solver coordinates.

Owned elements count once. Neighbouring trace contributions are added before
norm evaluation; ghost elements supply contributions to owned faces. Prescribed
traces enter through the `dir` offsets and have no unknown trace rows. Face
orientations and periodic paired-column mappings follow production assembly.

## Implementation and cost

`hdg_Residual.f90` provides a residual-monitor type owning geometry caches,
norms, convergence history and reporting policy. The pure element evaluator
reads the existing matrix/state containers and an element index; it writes only
that element's three residual vectors. Newton assembly and the
`Convergence` driver call its evaluation/report methods; boundary assembly and
physical balance reporting have no residual dependency.

`begin_solve()` resets convergence history while retaining cached geometry.
`invalidate_mesh()` clears geometry and history after startup or a mesh rebuild.
The MHDG startup, adaptive solver-state rebuild, and standalone Convergence
setup notify the monitor after preparing their mesh matrices. The first-equation
routine itself has no monitoring dependency. These are distinct events;
a magnetic-field change alone does not invalidate the mass matrices.

Evaluation adds local matrix-vector products, trace accumulation and three MPI
reductions; timing adds one reduction. First use includes cache setup. No extra
global solve is needed, and no caches or residual datasets are saved to HDF5.
`computeResidual` remains the volume-change measure used for update/time checks.
