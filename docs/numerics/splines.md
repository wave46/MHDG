# One-dimensional spline interpolation

`mod_splines` provides a cubic interpolant of sampled data and evaluates its
value and first derivative. This guide covers the retained Fortran utility.
For input interpolation, see the [equilibrium](../inputs/equilibrium.md#time-sampling)
and [heating](../inputs/heating.md#time-sampling-and-restart) guides.

## Data and interpolation conventions

Supply matching one-dimensional double-precision arrays $x_i,y_i$, with
**at least four points** and strictly increasing $x_i$. Use finite values;
the constructor does not validate lengths, ordering or finiteness.

Within each interval, the interpolant has the form

$$
S_i(x)=a_i(x-x_i)^3+b_i(x-x_i)^2+c_i(x-x_i)+d_i,
\qquad
S_i'(x)=3a_i(x-x_i)^2+2b_i(x-x_i)+c_i.
$$

It passes through the supplied values with continuous first and second
derivatives. End conditions are
[not-a-knot](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html):
the first two intervals share one cubic polynomial, as do the last two.
The optional endpoint-slope arguments do not affect this construction.

Outside the supplied range, evaluation continues the first or last **cubic
polynomial**, including its derivative. It does not switch to linear
extrapolation or clamp the value at an endpoint.

The utility performs no normalization or unit conversion. Values retain the
units of $y$; slopes have units of $y/x$.

## Fortran usage

Import `spline` from `mod_splines`, construct an object from the samples,
then evaluate it:

```fortran
program interpolate_profile
  use mod_splines, only: spline
  implicit none
  type(spline) :: sp
  real*8 :: xi(9), yi(9), x
  integer :: i

  xi = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, &
        1.25d0, 1.5d0, 1.75d0, 2.0d0]
  yi = [18.0d0, 18.4921875d0, 18.9375d0, 19.2890625d0, &
        19.5d0, 19.5234375d0, 19.3125d0, 18.8203125d0, 18.0d0]
  sp = spline(xi, yi)

  do i = 0, 10
    x = 0.2d0*i
    print *, x, sp%value(x), sp%slope(x)
  end do
end program interpolate_profile
```

| Call | Result |
| --- | --- |
| `sp = spline(xi, yi)` | Construct the interpolant |
| `sp%value(x)` | Interpolated or extrapolated value |
| `sp%slope(x)` | First derivative with respect to $x$ |
| `sp%indexof(x)` | Interval index used for evaluation; first/last interval outside the data range |

`value`, `slope` and `indexof` are elemental: their evaluation coordinate
can be a scalar or an array. Array evaluation returns matching result shapes.

## Implementation and current limitations

The [module](../../src/Spline/mod_splines.f90) solves a tridiagonal system
for nodal slopes and stores polynomial and derivative coefficients per
interval. The current solver retains the type and build dependency, with
no interpolation-method calls in the tracked solver sources.

| Interface or input | Current limitation |
| --- | --- |
| Optional `y1_slope`, `yn_slope` | Accepted by the constructor but ignored for four or more points; they cannot impose endpoint slopes |
| Three or fewer points | The alternate construction branch leaves evaluation coefficients unallocated; `value` and `slope` cannot evaluate its result |
| `slope2`, `extrema` | Absent from the current type; second-derivative evaluation and extremum search are not provided |

## Attribution

The bundled guide credits **FortranConsoleSpline**, © John Alexiou, and
identifies it as MIT licensed. This checkout has no accompanying license text.
