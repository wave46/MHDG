# Intentional exact real comparisons

GNU Fortran's `-Wcompare-reals` diagnoses equality and inequality comparisons
of real values; it is enabled by `-Wextra`. A diagnostic identifies a comparison
for review, not necessarily an error. See the
[GNU Fortran warning documentation](https://gcc.gnu.org/onlinedocs/gcc-9.4.0/gfortran/Error-and-Warning-Options.html).

## Review policy

Keep `-Wcompare-reals` enabled. The table distinguishes required validation and
inactive-state guards from optional fast paths. An exact comparison can be
appropriate without the conditional itself being necessary.

Do not replace these comparisons with arbitrary epsilon tests, equivalent
expressions intended to evade the warning, or broad warning suppression.
Comparisons of independently computed approximations require their own review;
when approximate agreement is intended, any tolerance must have a justified
scale and meaning.

This inventory is not a blanket exemption for a file or procedure. Review new
or changed comparisons and update the corresponding entry when its expression,
location or contract changes. The number of reviewed warnings alone cannot
identify whether a new warning is acceptable.

## Reviewed comparisons

The inventory covers 21 production diagnostics in the selected 2D
NGammaTiTeNeutral and NGammaTiTeNeutralGamma serial/MPI builds. Multiple
operators on one line count separately. Other models, test drivers and external
libraries are outside this inventory. Comparisons link to the current source
lines.

| ID | Source | Comparison | Current purpose | Review disposition |
| --- | --- | --- | --- | --- |
| 1 | [rectangular_series_io.f90:184](../../src/InOut/rectangular_series_io.f90#L184) | `ALL(this%r == other%r)` | Requires identical stored R axes. Approximate matching would relax the fixed-grid input contract. | Keep validation: removing it can apply later frame data on the wrong stored grid. |
| 2 | [rectangular_series_io.f90:184](../../src/InOut/rectangular_series_io.f90#L184) | `ALL(this%z == other%z)` | Requires identical stored Z axes under the same fixed-grid contract. | Keep validation: Z-axis consistency is independently required. |
| 3 | [rectangular_series_io.f90:263](../../src/InOut/rectangular_series_io.f90#L263) | `ANY(r2d /= SPREAD(r,1,SIZE(z)))` | Every R-grid row must reproduce the stored R axis. Accepting a nearly rectangular grid requires an explicit input/interpolation decision. | Keep validation: the reader retains only one R row; unchecked remaining rows could describe an unsupported grid. |
| 4 | [rectangular_series_io.f90:263](../../src/InOut/rectangular_series_io.f90#L263) | `ANY(z2d /= SPREAD(z,2,SIZE(r)))` | Every Z-grid column must reproduce the stored Z axis under the same rectangularity contract. | Keep validation: the reader retains only one Z column; unchecked remaining columns could describe an unsupported grid. |
| 5 | [equilibrium_series.f90:393](../../src/Models/equilibrium_series.f90#L393) | `ALL(left == right)` | Numbered current frames must retain identical stored grid axes. | Keep validation: separate-current frame arrays are interpreted on the initial stored current grid. |
| 6 | [time_stepping.f90:54](../../src/HDG/time_stepping.f90#L54) | `target /= time%t_final` | Target selection clips to the stored final endpoint using MIN. The exact endpoint is exempt from dt_min, permitting a shortened final step. | Keep behavior: removing the endpoint distinction changes acceptance of shortened final timesteps. |
| 7 | [external_heating.f90:52](../../src/Models/external_heating.f90#L52) | `target_time == applied_time` | Heating is reused only at the identical requested time. A mismatch triggers resampling; a tolerance could reuse heating sampled at a different time. | Optional cache optimization, not defensive validation. Removable with extra sampling and mesh interpolation; keep it. |
| 8 | [read_input.f90:321](../../src/InOut/read_input.f90#L321) | `heating_ion_fraction == heating_ion_fraction_unset` | The omitted value is initialized to an exact -HUGE sentinel. Only that value requests the heating_equation default. | Keep default selection: omission and an explicit split have different meanings. |
| 9 | [read_input.f90:353](../../src/InOut/read_input.f90#L353) | `Re_n_pump == Re_n_pump_unset` | The omitted value is initialized to an exact -HUGE sentinel. Only that value requests inheritance from Re_n. | Keep default selection: omission inherits Re_n; an explicit value overrides it. |
| 10 | [read_input.f90:1084](../../src/InOut/read_input.f90#L1084) | `t_final /= time_input_unset` | Tests whether t_final was supplied, using the exact unset sentinel; this is an input-presence check. | Keep input validation: catches a final-time parameter that legacy stepping would otherwise ignore. |
| 11 | [read_input.f90:1089](../../src/InOut/read_input.f90#L1089) | `tfi /= time_input_unset` | Tests whether the deprecated tfi alias was supplied, using the exact unset sentinel. | Keep alias detection: do not overwrite t_final with the unset alias. |
| 12 | [read_input.f90:1091](../../src/InOut/read_input.f90#L1091) | `t_final /= time_input_unset` | Tests whether t_final was also supplied before comparing the two aliases. | Keep presence check: alias-only input must remain accepted. |
| 13 | [read_input.f90:1091](../../src/InOut/read_input.f90#L1091) | `t_final /= tfi` | The current contract requires both supplied aliases to represent exactly the same floating-point value. Approximate agreement would be an explicit compatibility-policy change. | Keep conflict diagnostic: removing it silently gives tfi precedence over an explicitly different t_final. |
| 14 | [read_input.f90:1099](../../src/InOut/read_input.f90#L1099) | `t_start /= 0.d0` | Legacy mode requires the exact default start time. A nonzero start belongs to an explicit physical-time mode. | Keep input validation: catches a start time that legacy stepping would otherwise ignore. |
| 15 | [read_input.f90:1102](../../src/InOut/read_input.f90#L1102) | `t_final == time_input_unset` | Only the exact unset sentinel is replaced by the legacy default. | Keep normalization of omitted input: leaving the sentinel would also leak it into saved Final_time metadata. |
| 16 | [read_input.f90:1109](../../src/InOut/read_input.f90#L1109) | `t_final == time_input_unset` | Only the exact unset sentinel denotes an omitted required final time; finiteness and ordering are validated separately. | Diagnostic-only redundancy: the later final-time ordering check also rejects the unset sentinel for finite t_start. Keep the clearer missing-parameter message. |
| 17 | [physics.f90:871](../../src/Models/NGammaTiTe/physics.f90#L871) | `Re_n == 1.d0` | Full recycling has exactly zero neutral absorption. Values below one retain their nonzero contribution, however small. | Keep functional guard: inactive absorption must return before singular temperature evaluation, including all-zero trial states. |
| 18 | [physics.f90:884](../../src/Models/NGammaTiTe/physics.f90#L884) | `Re_n == 1.d0` | Uses the same full-recycling branch as the flux, preserving a matching zero Jacobian. | Keep functional guard for the same reason; the helper guarantees a fully zero derivative array. |
| 19 | [hdg_ComputeJacobian.f90:1704](../../src/HDG/hdg_ComputeJacobian.f90#L1704) | `puff_wall_integral .EQ. 0.d0` | An exactly zero puff load contributes nothing to relocation. A tolerance could discard a small nonzero source. | Keep the combined fast path: avoid repeated element-volume integration where no wall source contributes. |
| 20 | [hdg_ComputeJacobian.f90:1704](../../src/HDG/hdg_ComputeJacobian.f90#L1704) | `ALL(pump_wall_load .EQ. 0.d0)` | The relocation early return also requires every pump-load component to be exactly zero. | Keep the combined fast path; both loads must be zero before skipping the element. |
| 21 | [hdg_ComputeJacobian.f90:1861](../../src/HDG/hdg_ComputeJacobian.f90#L1861) | `ALL(pump_wall_load .EQ. 0.d0)` | Skips assembly only when every pump-load component is exactly zero; all nonzero loads remain active. | Optional fast path. Removing it forms and applies a zero pump matrix; keep to avoid that extra work. |

## Contract changes needing separate discussion

Approximate grid matching or accepting nearly rectangular coordinates would
change the input contract in entries 1–5. Approximate agreement between `tfi`
and `t_final` would change entry 13. Such changes should define the intended
behavior and relevant tolerance before implementation, with focused input
checks. The final-endpoint and heating-cache decisions in entries 6–7 should
also be revisited if time selection or caching semantics change.

## Replacing equality with an inequality

For a quantity guaranteed to be nonnegative, `x <= 0` and `x == 0` are
equivalent on valid inputs. This is a domain argument, not an approximate
floating-point comparison.

The current static/prescribed puff paths and static cryopump coefficient do
not consistently enforce nonnegativity. Replacing their zero-load checks with
nonpositive-load checks can therefore silently discard negative contributions
that the present code assembles. If negative physical strengths are forbidden,
that needs an explicit input contract; an assembly fast path should not silently
normalize them to inactivity. Assembled nodal load vectors also contain
basis-function contributions, so do not assume componentwise nonnegativity from
the physical coefficient alone.

Keep the current exact-zero tests until the relevant domain and intended
negative-input behavior have been established.
