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
      END IF

      left = 1
      right = length
      DO
         IF (left > right) THEN
            EXIT
         END IF
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
         WRITE (6, *) "Problem in the binary search"
         WRITE (6, *) "x", x
         WRITE (6, *) "max(x_array)", MAXVAL(x_array)
         WRITE (6, *) "min(x_array)", MINVAL(x_array)
      END IF

      IF (j == y_len) THEN
         WRITE (6, *) "Problem in the binary search"
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

   FUNCTION nodesearch(x, y, xy_len, x_array, y_array)
      ! Given a  point (x,y), returns the index of the closest node in 2D
      IMPLICIT NONE
      INTEGER, INTENT(IN)               :: xy_len
      REAL, INTENT(IN)                  :: x, y
      REAL, DIMENSION(xy_len), INTENT(IN):: x_array, y_array
      REAL*8                           :: d(xy_len)

      INTEGER                     :: nodesearch

      d = SQRT((x_array - x)**2 + (y_array - y)**2)
      nodesearch = MINLOC(d, DIM=1)

   END FUNCTION nodesearch

   SUBROUTINE lineintegration(qp_len, x_vec, y_vec, f, Xc, T, nodes2D, nli)
      ! Given a set of points (x, y) along a line of sight, retruns the line integration of
      ! f. The function f is evaluated in the closest nodes to points (x, y).
      IMPLICIT NONE
      INTEGER, INTENT(IN)          :: qp_len, nodes2D, T(:, :)
      REAL*8, INTENT(IN)           :: x_vec(:), y_vec(:), f(:), Xc(:, :)
      REAL*8, INTENT(OUT)          :: nli
      INTEGER                     :: i, iel, inp, n_ind, f_ind(2), np_len
      REAL*8                      :: x, y, dl(qp_len - 1), x_qp(qp_len), f_qp(qp_len)

      ! Search closest node and evaluate f
      DO i = 1, qp_len
         x = x_vec(i)
         y = y_vec(i)
         np_len = SIZE(Xc, 1)
         n_ind = nodesearch(x, y, np_len, Xc(:, 1), Xc(:, 2))
         f_ind = FINDLOC(T, n_ind)
         iel = f_ind(1)
         inp = f_ind(2)
         x_qp(i) = Xc(n_ind, 1)
         f_qp(i) = f((iel - 1)*nodes2D + inp)
      END DO

      dl = (x_qp(2:qp_len) - x_qp(1:qp_len - 1))*1.901e-3
      nli = SUM(0.5*(f_qp(2:qp_len) + f_qp(1:qp_len - 1))*dl)

   END SUBROUTINE lineintegration

   !real function lineintegration(qp_len, x_vec, y_vec, np_len, x_array, y_array ,f)
   ! Given a set of points (x, y) along a line of sight, retruns the line integration of
   ! f. The function f is evaluated in the closest nodes to points (x, y).
   !  implicit none
   !  real,intent(in)            :: qp_len, np_len
   !  real,dimension(qp_len),intent(in) :: x_vec,  y_vec
   !  real,dimension(np_len),intent(in) :: x_array,  y_array
   !  real, intent(in)           :: f
   !  integer                    :: i, iel, inp, n_ind, f_ind(2)
   !  real*8                     :: x, y, dl(qp_len-1), x_qp(qp_len), f_qp(qp_len)

   ! Search closest node and evaluate f
   !  do i = 1, qp_len
   !     x = x_vec(i)
   !     y = y_vec(i)
   !     n_ind = nodesearch(x, y, np_len, x_array, y_array)
   !     f_ind = findloc(Mesh%T(n_ind))
   !     iel = f_ind(1)
   !     ip = f_ind(2)
   !     x_qp(i) = x_array(n_ind)
   !     f_qp(i) = f((iel-1)*refElPol%N2D + ip)
   !  end do

   ! Integrate
   !  dl = x_qp(2:qp_len) - x_qp(1:qp_len-1)
   !  lineintegration = sum(0.5*(f_qp(2:qp_len) + f_qp(1:qp_len-1))*dl)

   !end function lineintegration

   PURE SUBROUTINE spline_coeff(x, y, n, y2)
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: n
      REAL, DIMENSION(n), INTENT(IN) :: x, y
      REAL, DIMENSION(n), INTENT(OUT) :: y2

      REAL, DIMENSION(n)            :: u
      INTEGER                       :: i, k
      REAL                          :: sig, p

      y2(1) = 0.0
      u(1) = 0.0

      DO i = 2, n - 1
         sig = (x(i) - x(i - 1))/(x(i + 1) - x(i - 1))
         p = sig*y2(i - 1) + 2.0
         y2(i) = (sig - 1.0)/p
         u(i) = (y(i + 1) - y(i))/(x(i + 1) - x(i)) - (y(i) - y(i - 1))/(x(i) - x(i - 1))
         u(i) = (6.0*u(i)/(x(i + 1) - x(i - 1)) - sig*u(i - 1))/p
      END DO

      y2(n) = 0.0

      DO k = n - 1, 1, -1
         y2(k) = y2(k)*y2(k + 1) + u(k)
      END DO
   END SUBROUTINE spline_coeff

   PURE SUBROUTINE spline_eval(xa, ya, y2a, n, x, y)
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: n
      REAL, DIMENSION(n), INTENT(IN) :: xa, ya, y2a
      REAL, INTENT(IN)              :: x
      REAL, INTENT(OUT)             :: y

      INTEGER                       :: klo, khi, k
      REAL                          :: h, a, b

      klo = 1
      khi = n
      DO WHILE (khi - klo > 1)
         k = (khi + klo)/2
         IF (xa(k) > x) THEN
            khi = k
         ELSE
            klo = k
         END IF
      END DO

      h = xa(khi) - xa(klo)
      IF (ABS(h) .LT. 1.e-12) THEN
         y = ya(klo)
         RETURN
      END IF

      a = (xa(khi) - x)/h
      b = (x - xa(klo))/h
      y = a*ya(klo) + b*ya(khi) + ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))*(h**2)/6.0

   END SUBROUTINE spline_eval

   PURE SUBROUTINE spline_eval_vec(xa, ya, y2a, n, xv, yv)
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: n
      REAL, DIMENSION(n), INTENT(IN) :: xa, ya, y2a
      REAL, DIMENSION(:), INTENT(IN) :: xv
      REAL, DIMENSION(SIZE(xv)), INTENT(OUT) :: yv
      INTEGER                       :: i

      DO CONCURRENT(i=1:SIZE(xv))
         CALL spline_eval(xa, ya, y2a, n, xv(i), yv(i))
      END DO

   END SUBROUTINE spline_eval_vec

   ! real pure function spline2d(x_len, x_array, y_len, y_array, f, x, y)
   !    implicit none
   !    integer, intent(in)           :: x_len, y_len
   !    real, dimension(x_len), intent(in) :: x_array
   !    real, dimension(y_len), intent(in) :: y_array
   !    real, dimension(x_len, y_len), intent(in) :: f
   !    real, intent(in)              :: x, y

   !    real                          :: y2(y_len), temp(y_len), x2(x_len)
   !    integer                       :: j
   !    real                          :: val

   !    ! Spline along x for each y_j (contiguous memory access in Fortran)
   !    do j = 1, y_len
   !       call spline_coeff(x_array, f(:, j), x_len, x2)
   !       call spline_eval(x_array, f(:, j), x2, x_len, x, temp(j))
   !    end do

   !    ! Spline along y
   !    call spline_coeff(y_array, temp, y_len, y2)
   !    call spline_eval(y_array, temp, y2, y_len, y, val)
   !    spline2d = val

   ! end function spline2d

   PURE SUBROUTINE spline2d_vec(x_len, x_array, y_len, y_array, f, xv, yv, fv)
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: x_len, y_len
      REAL, DIMENSION(x_len), INTENT(IN) :: x_array
      REAL, DIMENSION(y_len), INTENT(IN) :: y_array
      REAL, DIMENSION(x_len, y_len), INTENT(IN) :: f
      REAL, DIMENSION(:), INTENT(IN) :: xv
      REAL, DIMENSION(SIZE(xv)), INTENT(IN) :: yv
      REAL, DIMENSION(SIZE(xv)), INTENT(OUT) :: fv

      REAL, DIMENSION(x_len, y_len) :: x2mat
      REAL, DIMENSION(SIZE(xv), y_len) :: tmp
      REAL, DIMENSION(y_len)    :: y_work
      REAL, DIMENSION(y_len)    :: y2
      INTEGER                   :: j, p

      ! First stage: for each y_j, evaluate spline in x at all xv(:)
      DO j = 1, y_len
         CALL spline_coeff(x_array, f(:, j), x_len, x2mat(:, j))
         CALL spline_eval_vec(x_array, f(:, j), x2mat(:, j), x_len, xv, tmp(:, j))
      END DO

      ! Second stage: for each query point, spline in y
      DO p = 1, SIZE(xv)
         y_work = tmp(p, :)
         CALL spline_coeff(y_array, y_work, y_len, y2)
         CALL spline_eval(y_array, y_work, y2, y_len, yv(p), fv(p))
      END DO
   END SUBROUTINE spline2d_vec

END MODULE interpolation
