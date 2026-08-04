MODULE interpolation
  USE printutils
CONTAINS

  FUNCTION binarysearch(length, array, VALUE, delta)
    ! Given an array and a value, returns the index of the element that
    ! is closest to, but less than, the given value.
    ! Uses a binary search algorithm.
    ! "delta" is the tolerance used to determine if two values are equal
    ! if ( abs(x1 - x2) <= delta) then
    !    assume x1 = x2
    ! endif

    IMPLICIT NONE
    INTEGER, INTENT(IN)                 :: length
    REAL, DIMENSION(length), INTENT(IN) :: array
    !f2py depend(length) array
    REAL, INTENT(IN)                    :: VALUE
    REAL, INTENT(IN), OPTIONAL          :: delta

    INTEGER                             :: binarysearch

    INTEGER                             :: left, middle, right
    REAL                                :: d

    IF (PRESENT(delta) .EQV. .TRUE.) THEN
       d = delta
    ELSE
       d = 1e-9
    ENDIF

    left = 1
    right = length
    DO
       IF (left > right) THEN
          EXIT
       ENDIF
       middle = NINT((left + right)/2.0)
       IF (ABS(array(middle) - VALUE) <= d) THEN
          binarySearch = middle
          RETURN
       ELSE IF (array(middle) > VALUE) THEN
          right = middle - 1
       ELSE
          left = middle + 1
       END IF
    END DO
    binarysearch = right

  END FUNCTION binarysearch

  REAL FUNCTION interpolate(x_len, x_array, y_len, y_array, f, x, y, delta)
    ! This function uses bilinear interpolation to estimate the value
    ! of a function f at point (x,y)
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by x_array and the grid y values specified by y_array
    ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
    IMPLICIT NONE
    INTEGER, INTENT(IN)                       :: x_len, y_len
    REAL, DIMENSION(x_len), INTENT(IN)        :: x_array
    REAL, DIMENSION(y_len), INTENT(IN)        :: y_array
    REAL, DIMENSION(x_len, y_len), INTENT(IN) :: f
    REAL, INTENT(IN)                          :: x, y
    REAL, INTENT(IN), OPTIONAL                :: delta
    !f2py depend(x_len) x_array, f
    !f2py depend(y_len) y_array, f

    REAL                                      :: denom, x1, x2, y1, y2
    INTEGER                                   :: i, j

    i = binarysearch(x_len, x_array, x, delta)
    j = binarysearch(y_len, y_array, y, delta)

    IF (i == x_len) THEN
       WRITE (6, *) "Problem IN the binary search"
       WRITE (6, *) "x", x
       WRITE (6, *) "max(x_array)", MAXVAL(x_array)
       WRITE (6, *) "min(x_array)", MINVAL(x_array)
    END IF

    IF (j == y_len) THEN
       WRITE (6, *) "Problem IN the binary search"
       WRITE (6, *) "y", y
       WRITE (6, *) "max(y_array)", MAXVAL(y_array)
       WRITE (6, *) "min(y_array)", MINVAL(y_array)
    END IF

    x1 = x_array(i)
    x2 = x_array(i + 1)

    y1 = y_array(j)
    y2 = y_array(j + 1)

    denom = (x2 - x1)*(y2 - y1)

    interpolate = (f(i, j)*(x2 - x)*(y2 - y) + f(i + 1, j)*(x - x1)*(y2 - y) + &
         f(i, j + 1)*(x2 - x)*(y - y1) + f(i + 1, j + 1)*(x - x1)*(y - y1))/denom

  END FUNCTION interpolate

   !-----------------------------------------------------------------------
   ! Locate the grid cell containing a coordinate and compute local parameter
   ! t in [0,1] for that cell.
   !
   ! If value is outside the grid, we clamp to the first/last cell so that
   ! interpolation routines can still evaluate safely at boundaries.
   !-----------------------------------------------------------------------
   SUBROUTINE find_cell_and_local_coordinate(n, vec, value, idx, t)
      INTEGER, INTENT(IN) :: n
      REAL*8, INTENT(IN) :: vec(n), value
      INTEGER, INTENT(OUT) :: idx
      REAL*8, INTENT(OUT) :: t
      INTEGER :: lo, hi, mid
      REAL*8 :: den

      IF (value <= vec(1)) THEN
         idx = 1
      ELSEIF (value >= vec(n)) THEN
         idx = n - 1
      ELSE
         lo = 1
         hi = n
         DO WHILE (hi - lo > 1)
            mid = (lo + hi)/2
            IF (value >= vec(mid)) THEN
               lo = mid
            ELSE
               hi = mid
            ENDIF
         ENDDO
         idx = lo
      ENDIF

      den = vec(idx + 1) - vec(idx)
      IF (ABS(den) < 1.d-16) THEN
         t = 0.d0
      ELSE
         t = (value - vec(idx))/den
      ENDIF
      t = MAX(0.d0, MIN(1.d0, t))
   END SUBROUTINE find_cell_and_local_coordinate

   !-----------------------------------------------------------------------
   ! Build derivative fields used by bicubic Hermite interpolation.
   !
   ! Method:
   ! - fx, fy: first derivatives from finite differences on a rectilinear grid
   !   (centered in interior, one-sided at boundaries).
   ! - fxy: derivative of fx with respect to y, computed with the same stencil.
   !
   ! Note: this is a local bicubic Hermite patch method (often called
   ! "bicubic spline" in practice), not a global spline solve.
   !-----------------------------------------------------------------------
   SUBROUTINE build_bicubic_derivatives(ny, yvec, nx, xvec, f, fx, fy, fxy)
      INTEGER, INTENT(IN) :: ny, nx
      REAL*8, INTENT(IN) :: yvec(ny), xvec(nx)
      REAL*8, INTENT(IN) :: f(ny, nx)
      REAL*8, INTENT(OUT) :: fx(ny, nx), fy(ny, nx), fxy(ny, nx)
      INTEGER :: i, j
      REAL*8 :: dx, dy

      DO i = 1, ny
          DO j = 1, nx
               IF (j == 1) THEN
                   dx = xvec(2) - xvec(1)
                   fx(i, j) = (f(i, 2) - f(i, 1))/dx
               ELSEIF (j == nx) THEN
                   dx = xvec(nx) - xvec(nx - 1)
                   fx(i, j) = (f(i, nx) - f(i, nx - 1))/dx
               ELSE
                   dx = xvec(j + 1) - xvec(j - 1)
                   fx(i, j) = (f(i, j + 1) - f(i, j - 1))/dx
               ENDIF

               IF (i == 1) THEN
                   dy = yvec(2) - yvec(1)
                   fy(i, j) = (f(2, j) - f(1, j))/dy
               ELSEIF (i == ny) THEN
                   dy = yvec(ny) - yvec(ny - 1)
                   fy(i, j) = (f(ny, j) - f(ny - 1, j))/dy
               ELSE
                   dy = yvec(i + 1) - yvec(i - 1)
                   fy(i, j) = (f(i + 1, j) - f(i - 1, j))/dy
               ENDIF
          ENDDO
      ENDDO

      DO i = 1, ny
          DO j = 1, nx
               IF (i == 1) THEN
                   dy = yvec(2) - yvec(1)
                   fxy(i, j) = (fx(2, j) - fx(1, j))/dy
               ELSEIF (i == ny) THEN
                   dy = yvec(ny) - yvec(ny - 1)
                   fxy(i, j) = (fx(ny, j) - fx(ny - 1, j))/dy
               ELSE
                   dy = yvec(i + 1) - yvec(i - 1)
                   fxy(i, j) = (fx(i + 1, j) - fx(i - 1, j))/dy
               ENDIF
          ENDDO
      ENDDO
   END SUBROUTINE build_bicubic_derivatives

   !-----------------------------------------------------------------------
   ! Evaluate bicubic Hermite interpolant value at (x,y).
   !
   ! Reference formulation:
   ! - Bicubic interpolation as tensor product of 1D cubic Hermite basis.
   ! - Basis functions used here are:
   !     h00(t)= 2t^3-3t^2+1
   !     h10(t)= t^3-2t^2+t
   !     h01(t)=-2t^3+3t^2
   !     h11(t)= t^3-t^2
   !
   ! These correspond to the standard Hermite form described in common
   ! numerical analysis texts (e.g. Numerical Recipes) and the bicubic
   ! interpolation article on Wikipedia.
   !-----------------------------------------------------------------------
   SUBROUTINE eval_bicubic_value(ny, yvec, nx, xvec, f, fx, fy, fxy, y, x, val)
      INTEGER, INTENT(IN) :: ny, nx
      REAL*8, INTENT(IN) :: yvec(ny), xvec(nx)
      REAL*8, INTENT(IN) :: f(ny, nx), fx(ny, nx), fy(ny, nx), fxy(ny, nx)
      REAL*8, INTENT(IN) :: y, x
      REAL*8, INTENT(OUT) :: val
      INTEGER :: iy, ix
      REAL*8 :: ty, tx, dy, dx
      REAL*8 :: h00x, h10x, h01x, h11x, h00y, h10y, h01y, h11y
      REAL*8 :: a0, a1, b0, b1

      CALL find_cell_and_local_coordinate(ny, yvec, y, iy, ty)
      CALL find_cell_and_local_coordinate(nx, xvec, x, ix, tx)

      dy = yvec(iy + 1) - yvec(iy)
      dx = xvec(ix + 1) - xvec(ix)

      h00x = 2.d0*tx**3 - 3.d0*tx**2 + 1.d0
      h10x = tx**3 - 2.d0*tx**2 + tx
      h01x = -2.d0*tx**3 + 3.d0*tx**2
      h11x = tx**3 - tx**2

      h00y = 2.d0*ty**3 - 3.d0*ty**2 + 1.d0
      h10y = ty**3 - 2.d0*ty**2 + ty
      h01y = -2.d0*ty**3 + 3.d0*ty**2
      h11y = ty**3 - ty**2

      a0 = h00x*f(iy,ix) + h10x*dx*fx(iy,ix) + h01x*f(iy,ix+1) + h11x*dx*fx(iy,ix+1)
      a1 = h00x*f(iy+1,ix) + h10x*dx*fx(iy+1,ix) + h01x*f(iy+1,ix+1) + h11x*dx*fx(iy+1,ix+1)
      b0 = h00x*dy*fy(iy,ix) + h10x*dx*dy*fxy(iy,ix) + h01x*dy*fy(iy,ix+1) + h11x*dx*dy*fxy(iy,ix+1)
      b1 = h00x*dy*fy(iy+1,ix) + h10x*dx*dy*fxy(iy+1,ix) + h01x*dy*fy(iy+1,ix+1) + h11x*dx*dy*fxy(iy+1,ix+1)

      val = h00y*a0 + h10y*b0 + h01y*a1 + h11y*b1
   END SUBROUTINE eval_bicubic_value

   !-----------------------------------------------------------------------
   ! Evaluate bicubic Hermite interpolant and first derivatives at (x,y).
   !
   ! dval_dx and dval_dy are obtained by differentiating the Hermite basis
   ! analytically and applying chain rule factors 1/dx and 1/dy.
   !-----------------------------------------------------------------------
   SUBROUTINE eval_bicubic_with_derivatives(ny, yvec, nx, xvec, f, fx, fy, fxy, y, x, val, dval_dy, dval_dx)
      INTEGER, INTENT(IN) :: ny, nx
      REAL*8, INTENT(IN) :: yvec(ny), xvec(nx)
      REAL*8, INTENT(IN) :: f(ny, nx), fx(ny, nx), fy(ny, nx), fxy(ny, nx)
      REAL*8, INTENT(IN) :: y, x
      REAL*8, INTENT(OUT) :: val, dval_dy, dval_dx
      INTEGER :: iy, ix
      REAL*8 :: ty, tx, dy, dx
      REAL*8 :: h00x, h10x, h01x, h11x, h00y, h10y, h01y, h11y
      REAL*8 :: dh00x, dh10x, dh01x, dh11x, dh00y, dh10y, dh01y, dh11y
      REAL*8 :: a0, a1, b0, b1

      CALL find_cell_and_local_coordinate(ny, yvec, y, iy, ty)
      CALL find_cell_and_local_coordinate(nx, xvec, x, ix, tx)

      dy = yvec(iy + 1) - yvec(iy)
      dx = xvec(ix + 1) - xvec(ix)

      h00x = 2.d0*tx**3 - 3.d0*tx**2 + 1.d0
      h10x = tx**3 - 2.d0*tx**2 + tx
      h01x = -2.d0*tx**3 + 3.d0*tx**2
      h11x = tx**3 - tx**2

      h00y = 2.d0*ty**3 - 3.d0*ty**2 + 1.d0
      h10y = ty**3 - 2.d0*ty**2 + ty
      h01y = -2.d0*ty**3 + 3.d0*ty**2
      h11y = ty**3 - ty**2

      dh00x = 6.d0*tx**2 - 6.d0*tx
      dh10x = 3.d0*tx**2 - 4.d0*tx + 1.d0
      dh01x = -6.d0*tx**2 + 6.d0*tx
      dh11x = 3.d0*tx**2 - 2.d0*tx

      dh00y = 6.d0*ty**2 - 6.d0*ty
      dh10y = 3.d0*ty**2 - 4.d0*ty + 1.d0
      dh01y = -6.d0*ty**2 + 6.d0*ty
      dh11y = 3.d0*ty**2 - 2.d0*ty

      a0 = h00x*f(iy,ix) + h10x*dx*fx(iy,ix) + h01x*f(iy,ix+1) + h11x*dx*fx(iy,ix+1)
      a1 = h00x*f(iy+1,ix) + h10x*dx*fx(iy+1,ix) + h01x*f(iy+1,ix+1) + h11x*dx*fx(iy+1,ix+1)
      b0 = h00x*dy*fy(iy,ix) + h10x*dx*dy*fxy(iy,ix) + h01x*dy*fy(iy,ix+1) + h11x*dx*dy*fxy(iy,ix+1)
      b1 = h00x*dy*fy(iy+1,ix) + h10x*dx*dy*fxy(iy+1,ix) + h01x*dy*fy(iy+1,ix+1) + h11x*dx*dy*fxy(iy+1,ix+1)

      val = h00y*a0 + h10y*b0 + h01y*a1 + h11y*b1

      dval_dx = h00y*(dh00x*f(iy,ix)/dx + dh10x*fx(iy,ix) + dh01x*f(iy,ix+1)/dx + dh11x*fx(iy,ix+1)) + &
                     h10y*(dh00x*dy*fy(iy,ix)/dx + dh10x*dy*fxy(iy,ix) + dh01x*dy*fy(iy,ix+1)/dx + dh11x*dy*fxy(iy,ix+1)) + &
                     h01y*(dh00x*f(iy+1,ix)/dx + dh10x*fx(iy+1,ix) + dh01x*f(iy+1,ix+1)/dx + dh11x*fx(iy+1,ix+1)) + &
                     h11y*(dh00x*dy*fy(iy+1,ix)/dx + dh10x*dy*fxy(iy+1,ix) + dh01x*dy*fy(iy+1,ix+1)/dx + dh11x*dy*fxy(iy+1,ix+1))

      dval_dy = dh00y*a0/dy + dh10y*b0/dy + dh01y*a1/dy + dh11y*b1/dy
   END SUBROUTINE eval_bicubic_with_derivatives

   !-----------------------------------------------------------------------
   ! Evaluate bicubic Hermite interpolant, first derivatives, and pure second
   ! derivatives at (x,y).
   !
   ! d2val_dx2 and d2val_dy2 are formed from second derivatives of the
   ! Hermite basis.  The optional mixed derivative is useful when the
   ! interpolant is used to refine magnetic critical points.
   !-----------------------------------------------------------------------
   SUBROUTINE eval_bicubic_with_2nd_derivatives(ny, yvec, nx, xvec, f, fx, fy, fxy, y, x, val, dval_dy, dval_dx, d2val_dy2, d2val_dx2, d2val_dxdy)
      INTEGER, INTENT(IN) :: ny, nx
      REAL*8, INTENT(IN) :: yvec(ny), xvec(nx)
      REAL*8, INTENT(IN) :: f(ny, nx), fx(ny, nx), fy(ny, nx), fxy(ny, nx)
      REAL*8, INTENT(IN) :: y, x
      REAL*8, INTENT(OUT) :: val, dval_dy, dval_dx, d2val_dy2, d2val_dx2
      REAL*8, INTENT(OUT), OPTIONAL :: d2val_dxdy
      INTEGER :: iy, ix
      REAL*8 :: ty, tx, dy, dx
      REAL*8 :: h00x, h10x, h01x, h11x, h00y, h10y, h01y, h11y
      REAL*8 :: dh00x, dh10x, dh01x, dh11x, dh00y, dh10y, dh01y, dh11y
      REAL*8 :: d2h00x, d2h10x, d2h01x, d2h11x, d2h00y, d2h10y, d2h01y, d2h11y
      REAL*8 :: a0, a1, b0, b1
      REAL*8 :: da0_dx, da1_dx, db0_dx, db1_dx

      CALL find_cell_and_local_coordinate(ny, yvec, y, iy, ty)
      CALL find_cell_and_local_coordinate(nx, xvec, x, ix, tx)

      dy = yvec(iy + 1) - yvec(iy)
      dx = xvec(ix + 1) - xvec(ix)

      h00x = 2.d0*tx**3 - 3.d0*tx**2 + 1.d0
      h10x = tx**3 - 2.d0*tx**2 + tx
      h01x = -2.d0*tx**3 + 3.d0*tx**2
      h11x = tx**3 - tx**2

      h00y = 2.d0*ty**3 - 3.d0*ty**2 + 1.d0
      h10y = ty**3 - 2.d0*ty**2 + ty
      h01y = -2.d0*ty**3 + 3.d0*ty**2
      h11y = ty**3 - ty**2

      dh00x = 6.d0*tx**2 - 6.d0*tx
      dh10x = 3.d0*tx**2 - 4.d0*tx + 1.d0
      dh01x = -6.d0*tx**2 + 6.d0*tx
      dh11x = 3.d0*tx**2 - 2.d0*tx

      dh00y = 6.d0*ty**2 - 6.d0*ty
      dh10y = 3.d0*ty**2 - 4.d0*ty + 1.d0
      dh01y = -6.d0*ty**2 + 6.d0*ty
      dh11y = 3.d0*ty**2 - 2.d0*ty

      d2h00x = 12.d0*tx - 6.d0
      d2h10x = 6.d0*tx - 4.d0
      d2h01x = -12.d0*tx + 6.d0
      d2h11x = 6.d0*tx - 2.d0

      d2h00y = 12.d0*ty - 6.d0
      d2h10y = 6.d0*ty - 4.d0
      d2h01y = -12.d0*ty + 6.d0
      d2h11y = 6.d0*ty - 2.d0

      a0 = h00x*f(iy,ix) + h10x*dx*fx(iy,ix) + h01x*f(iy,ix+1) + h11x*dx*fx(iy,ix+1)
      a1 = h00x*f(iy+1,ix) + h10x*dx*fx(iy+1,ix) + h01x*f(iy+1,ix+1) + h11x*dx*fx(iy+1,ix+1)
      b0 = h00x*dy*fy(iy,ix) + h10x*dx*dy*fxy(iy,ix) + h01x*dy*fy(iy,ix+1) + h11x*dx*dy*fxy(iy,ix+1)
      b1 = h00x*dy*fy(iy+1,ix) + h10x*dx*dy*fxy(iy+1,ix) + h01x*dy*fy(iy+1,ix+1) + h11x*dx*dy*fxy(iy+1,ix+1)

      val = h00y*a0 + h10y*b0 + h01y*a1 + h11y*b1

      da0_dx = dh00x*f(iy,ix)/dx + dh10x*fx(iy,ix) + dh01x*f(iy,ix+1)/dx + dh11x*fx(iy,ix+1)
      da1_dx = dh00x*f(iy+1,ix)/dx + dh10x*fx(iy+1,ix) + dh01x*f(iy+1,ix+1)/dx + dh11x*fx(iy+1,ix+1)
      db0_dx = dh00x*dy*fy(iy,ix)/dx + dh10x*dy*fxy(iy,ix) + dh01x*dy*fy(iy,ix+1)/dx + dh11x*dy*fxy(iy,ix+1)
      db1_dx = dh00x*dy*fy(iy+1,ix)/dx + dh10x*dy*fxy(iy+1,ix) + dh01x*dy*fy(iy+1,ix+1)/dx + dh11x*dy*fxy(iy+1,ix+1)

      dval_dx = h00y*da0_dx + h10y*db0_dx + h01y*da1_dx + h11y*db1_dx

      dval_dy = dh00y*a0/dy + dh10y*b0/dy + dh01y*a1/dy + dh11y*b1/dy

      d2val_dx2 = h00y*(d2h00x*f(iy,ix)/dx**2 + d2h10x*fx(iy,ix)/dx + d2h01x*f(iy,ix+1)/dx**2 + d2h11x*fx(iy,ix+1)/dx) + &
                        h10y*(d2h00x*dy*fy(iy,ix)/dx**2 + d2h10x*dy*fxy(iy,ix)/dx + d2h01x*dy*fy(iy,ix+1)/dx**2 + d2h11x*dy*fxy(iy,ix+1)/dx) + &
                        h01y*(d2h00x*f(iy+1,ix)/dx**2 + d2h10x*fx(iy+1,ix)/dx + d2h01x*f(iy+1,ix+1)/dx**2 + d2h11x*fx(iy+1,ix+1)/dx) + &
                        h11y*(d2h00x*dy*fy(iy+1,ix)/dx**2 + d2h10x*dy*fxy(iy+1,ix)/dx + d2h01x*dy*fy(iy+1,ix+1)/dx**2 + d2h11x*dy*fxy(iy+1,ix+1)/dx)

      d2val_dy2 = d2h00y*a0/dy**2 + d2h10y*b0/dy**2 + d2h01y*a1/dy**2 + d2h11y*b1/dy**2

      IF (PRESENT(d2val_dxdy)) THEN
         d2val_dxdy = (dh00y*da0_dx + dh10y*db0_dx + dh01y*da1_dx + dh11y*db1_dx)/dy
      ENDIF
   END SUBROUTINE eval_bicubic_with_2nd_derivatives

  FUNCTION nodesearch(x, y, xy_len, x_array, y_array)
    ! Given a  point (x,y), returns the index of the closest node IN 2D
    IMPLICIT NONE
    INTEGER,INTENT(IN)               :: xy_len
    REAL,INTENT(IN)                  :: x, y
    REAL,DIMENSION(xy_len),INTENT(IN):: x_array, y_array
    REAL*8                           :: d(xy_len)

    INTEGER :: nodesearch

    d = SQRT( (x_array - x)**2 + (y_array - y)**2 )
    nodesearch = MINLOC(d, DIM = 1)

  END FUNCTION nodesearch
#ifndef PARALL
  SUBROUTINE lineintegration(qp_len, x_vec, y_vec, f, Xc, T, nodes2D, nli)
#else
   SUBROUTINE lineintegration(qp_len, x_vec, y_vec, f, Xc, T, nodes2D, ghost_elems, nli)
#endif
   ! Given a set of points (x, y) along a line of sight, returns the line integration of
   ! f. The function f is evaluated IN the closest nodes to points (x, y).
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: qp_len, nodes2D, T(:,:)
   REAL*8, INTENT(IN)  :: x_vec(:), y_vec(:), f(:), Xc(:,:)
   REAL*8, INTENT(OUT) :: nli
   INTEGER             :: i, inp, n_ind, f_ind(2), np_len
   REAL*8              :: x, y, x_prev, y_prev, f_prev, dl, f_cur, d(nodes2D)
   LOGICAL             :: is_inside
   INTEGER             :: element_index
#ifdef PARALL
   INTEGER, INTENT(IN) :: ghost_elems(:)
#endif

   ! Initialize result
   nli = 0.0D0
   x_prev = 0.0D0
   y_prev = 0.0D0
   f_prev = 0.0D0
   dl = 0.0D0
   f_cur = 0.0D0
   d = 0.0D0

   DO i = 1, qp_len
      x = x_vec(i)
      y = y_vec(i)
      np_len = SIZE(Xc, 1)

      ! Check if the point (x, y) is inside the mesh
      CALL point_in_mesh(x, y, Xc, T, is_inside, element_index)

      IF (.NOT. is_inside) THEN
         x_prev = x
         y_prev = y
         f_prev = 0.
         CYCLE  ! Skip this point if it is outside the mesh
      ENDIF
#ifdef PARALL
      ! Check if the element is a ghost element
      IF ((ghost_elems(element_index) .EQ. 1)) THEN
         x_prev = x
         y_prev = y
         f_prev = 0.
         CYCLE
      END IF
#endif

      ! Compute distances to the point (x, y)
      d = SQRT((Xc(T(element_index, :), 1) - x)**2 + (Xc(T(element_index, :), 2) - y)**2)

      ! Find the index of the closest node in the element
      inp = MINLOC(d, DIM = 1)
      f_cur = f((element_index-1)*nodes2D + inp)

      ! Compute the line segment length and update the summation
      IF (i > 1) THEN
         dl = SQRT((x - x_prev)**2 + (y - y_prev)**2) * 1.901e-3
         nli = nli + 0.5D0 * (f_cur + f_prev) * dl
      END IF

      ! Update previous values
      x_prev = x
      y_prev = y
      f_prev = f_cur
   END DO

  END SUBROUTINE lineintegration

  SUBROUTINE point_in_mesh(x, y, Xc, T, is_inside, element_index)
     IMPLICIT NONE
     REAL*8, INTENT(IN) :: x, y       ! Coordinates of the point
     REAL*8, INTENT(IN) :: Xc(:,:)    ! Node coordinates
     INTEGER, INTENT(IN) :: T(:,:)    ! Connectivity matrix
     LOGICAL, INTENT(OUT) :: is_inside
     INTEGER, INTENT(OUT) :: element_index
     INTEGER :: i, j
     REAL*8 :: x_nodes(3), y_nodes(3)
  
     is_inside = .FALSE.
     element_index = -1
  
     ! Loop over all elements in the mesh
     DO i = 1, SIZE(T, 1)
         ! Get the coordinates of the nodes of the current element
         DO j = 1, 3
             x_nodes(j) = Xc(T(i, j), 1)
             y_nodes(j) = Xc(T(i, j), 2)
         END DO
        
         ! Check if the point (x, y) is inside the current element
         IF (point_in_triangle(x, y, x_nodes, y_nodes)) THEN
             is_inside = .TRUE.
             element_index = i
             RETURN
         END IF
     END DO
  END SUBROUTINE point_in_mesh

  LOGICAL FUNCTION point_in_triangle(x, y, x_nodes, y_nodes)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x, y          ! Coordinates of the point
    REAL*8, INTENT(IN) :: x_nodes(3), y_nodes(3)  ! Triangle vertices
    REAL*8 :: area, area1, area2, area3

    ! Compute the area of the triangle
    area = 0.5D0 * ABS(x_nodes(1)*(y_nodes(2)-y_nodes(3)) + &
                       x_nodes(2)*(y_nodes(3)-y_nodes(1)) + &
                       x_nodes(3)*(y_nodes(1)-y_nodes(2)))

    ! Compute the areas of the sub-triangles formed with the point
    area1 = 0.5D0 * ABS(x*(y_nodes(2)-y_nodes(3)) + &
                        x_nodes(2)*(y_nodes(3)-y) + &
                        x_nodes(3)*(y-y_nodes(2)))

    area2 = 0.5D0 * ABS(x_nodes(1)*(y-y_nodes(3)) + &
                        x*(y_nodes(3)-y_nodes(1)) + &
                        x_nodes(3)*(y_nodes(1)-y))

    area3 = 0.5D0 * ABS(x_nodes(1)*(y_nodes(2)-y) + &
                        x_nodes(2)*(y-y_nodes(1)) + &
                        x*(y_nodes(1)-y_nodes(2)))

    ! Check if the sum of the sub-triangle areas equals the total area
    point_in_triangle = ABS(area - (area1 + area2 + area3)) < 1.0D-10
  END FUNCTION point_in_triangle

  !real function lineintegration(qp_len, x_vec, y_vec, np_len, x_array, y_array ,f)
  ! Given a set of points (x, y) along a line of sight, retruns the line integration of
  ! f. The function f is evaluated IN the closest nodes to points (x, y).
  !  implicit none
  !  real,intent(IN) :: qp_len, np_len
  !  real,dimension(qp_len),intent(IN) :: x_vec,  y_vec
  !  real,dimension(np_len),intent(IN) :: x_array,  y_array
  !  real, intent(IN) :: f
  !  integer :: i, iel, inp, n_ind, f_ind(2)
  !  real*8 :: x, y, dl(qp_len-1), x_qp(qp_len), f_qp(qp_len)

  ! Search closest node and evaluate f
  !  DO i = 1, qp_len
  !     x = x_vec(i)
  !     y = y_vec(i)
  !     n_ind = nodesearch(x, y, np_len, x_array, y_array)
  !     f_ind = findloc(Mesh%T(n_ind))
  !     iel = f_ind(1)
  !     ip = f_ind(2)
  !     x_qp(i) = x_array(n_ind)
  !     f_qp(i) = f((iel-1)*refElPol%N2D + ip)
  !  END DO

  ! Integrate
  !  dl = x_qp(2:qp_len) - x_qp(1:qp_len-1)
  !  lineintegration = sum(0.5*(f_qp(2:qp_len) + f_qp(1:qp_len-1))*dl)

  !end function lineintegration

END MODULE interpolation
