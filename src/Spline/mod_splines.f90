MODULE mod_splines
  IMPLICIT NONE


  TYPE :: spline
     REAL*8, ALLOCATABLE :: x(:), y(:), coeff(:,:), coeff_slope(:,:)
   CONTAINS
     PROCEDURE :: indexof => sp_index_of_x
     PROCEDURE :: VALUE => sp_interpolate_value
     PROCEDURE :: slope => sp_interpolate_slope
  END TYPE spline

  INTERFACE spline
     MODULE PROCEDURE :: sp_calculate_from_data
  END INTERFACE spline

CONTAINS

  ELEMENTAL FUNCTION sp_interpolate_value(sp,x) RESULT(y)
    CLASS(spline), INTENT(in) :: sp
    REAL*8, INTENT(IN)        :: x
    REAL*8                    :: y, t
    INTEGER                   :: k


    k = sp%indexof(x)
    t = (x - sp%x(k))
    y = sp%coeff(k,1)*t**3+sp%coeff(k,2)*t**2+sp%coeff(k,3)*t + sp%coeff(k,4)

  ENDFUNCTION sp_interpolate_value

  ELEMENTAL FUNCTION sp_interpolate_slope(sp,x)  RESULT(yp)
    CLASS(spline), INTENT(in) :: sp
    REAL*8, INTENT(IN)        :: x
    REAL*8                    :: yp, t
    INTEGER                   :: k


    k = sp%indexof(x)
    t = (x - sp%x(k))
    yp = sp%coeff_slope(k,1)*t**2+sp%coeff_slope(k,2)*t+sp%coeff_slope(k,3)

  ENDFUNCTION sp_interpolate_slope

  PURE FUNCTION sp_calculate_from_data(x,y,y1_slope,yn_slope) RESULT(sp)
    USE LinearAlgebra, ONLY: solve_linear_system_sing
    ! =====================================================
    ! Input x and y=f(x), n (dimension of x,y), (Ordered)
    ! y1 and yn are the first derivatives of f in the 1st point and the n-th
    ! Output: array y2(n) containing second derivatives of f(x_i)
    ! =====================================================

    TYPE(spline) :: sp
    REAL*8, INTENT(in) :: x(:), y(:)
    REAL*8 :: y2(SIZE(y)), dx(SIZE(y)-1), divdif(SIZE(y)-1), b(SIZE(y)), s(SIZE(y)), Mat(SIZE(y),3), dzzdx(SIZE(y)-1), dzdxdx(SIZE(y)-1), coeff(SIZE(y)-1,4), coeff_slope(SIZE(y)-1,3)
    REAL*8, OPTIONAL, INTENT(in) :: y1_slope, yn_slope
    REAL*8:: x31, xn
    REAL*8:: p, qn, sig, un, u(SIZE(y))
    INTEGER:: n, i, j

    n = SIZE(x)

    IF(n .GT. 3) THEN

       DO i = 1, n-1
          dx(i) = x(i+1) - x(i)
       ENDDO

       DO i = 1,n-1
          divdif(i) = (y(i+1) - y(i))/dx(i)
       ENDDO

       b = 0
       b(2:n-1) = 3*(dx(2:n-1)*divdif(1:n-2)+dx(1:n-2)*divdif(2:n-1))

       x31 = x(3) - x(1)
       xn  = x(n) - x(n-2)

       b(1) = ((dx(1)+2*x31)*dx(2)*divdif(1)+dx(1)**2*divdif(2))/x31
       b(n) = (dx(n-1)**2*divdif(n-2)+(2*xn+dx(n-1))*dx(n-2)*divdif(n-1))/xn

       ! diagonal -1
       Mat(1,1) = x31
       Mat(2:n-1,1) = dx(1:n-2)
       Mat(n,1) = 0
       ! diagonal 0
       Mat(1,2) = dx(2)
       Mat(2:n-1,2) = 2*(dx(2:n-1)+dx(1:n-2))
       Mat(n,2) = dx(n-2)
       ! diagonal + 1
       Mat(1,3) = 0
       Mat(2:n-1,3) = dx(2:n-1)
       Mat(n,3) = xn

       ! solve linear equation solution for the slopes
       CALL thomas_algorithm(Mat(:,3), Mat(:,2), Mat(:,1), b, s, n)

       ! now coefficients must be calculated

       dzzdx = (divdif-s(1:n-1))/dx
       dzdxdx = (s(2:n)-divdif)/dx

       coeff(:,1) = (dzdxdx-dzzdx)/dx
       coeff(:,2) = 2*dzzdx-dzdxdx
       coeff(:,3) = s(1:n-1)
       coeff(:,4) = y(1:n-1)

       coeff_slope(:,1) = 3.*coeff(:,1)
       coeff_slope(:,2) = 2.*coeff(:,2)
       coeff_slope(:,3) = coeff(:,3)

       sp%x = x
       sp%y = y
       sp%coeff = coeff
       sp%coeff_slope = coeff_slope


    ELSE

       IF (PRESENT(y1_slope)) THEN    ! natural spline conditions
          y2(1) = -0.5
          u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-y1_slope)
       ELSE
          y2(1) = 0
          u(1) = 0
       END IF

       DO i = 2, n-1                            ! tridiag. decomposition
          sig = (x(i)-(i-1))/(x(i+1)-x(i-1))
          p = sig*y2(i-1)+2.
          y2(i) = (sig-1.)/p
          u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
       END DO

       IF (PRESENT(yn_slope)) THEN   ! natural spline conditions
          qn = 0.5
          un=(3./(x(n)-x(n-1)))*(yn_slope-(y(n)-y(n-1))/(x(n)-x(n-1)))
       ELSE
          qn = 0
          un = 0
       END IF

       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

       DO j = n-1, 1, -1          !  backwards substitution tri-diagonale
          y2(j) = y2(j)*y2(j+1)+u(j)
       END DO

       sp%x = x
       sp%y = y
    ENDIF
    RETURN
  END FUNCTION sp_calculate_from_data

  PURE SUBROUTINE thomas_algorithm(a, b, c, d, x, n)
    IMPLICIT NONE

    INTEGER,INTENT(in) :: n
    REAL*8,DIMENSION(n),INTENT(in) :: a,b,c,d
    REAL*8,DIMENSION(n),INTENT(out) :: x
    REAL*8,DIMENSION(n) :: cp,dp
    REAL*8 :: m
    INTEGER i

    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)

    DO i = 2,n
       m = b(i)-cp(i-1)*a(i)
       cp(i) = c(i)/m
       dp(i) = (d(i)-dp(i-1)*a(i))/m
    END DO

    x(n) = dp(n)
    DO i = n-1, 1, -1
       x(i) = dp(i)-cp(i)*x(i+1)
    END DO
  END SUBROUTINE thomas_algorithm

  ELEMENTAL FUNCTION sp_index_of_x(sp,x) RESULT(k_low)
    CLASS(spline), INTENT(in) :: sp
    REAL*8, INTENT(in) :: x
    INTEGER:: n, k, k_low, k_high
    n = SIZE(sp%y)
    k_low = 1
    k_high = n
    IF(x<sp%x(k_low)) THEN
       RETURN
    ELSEIF (x>sp%x(k_high)) THEN
       k_low = k_high-1
       RETURN
    END IF
    DO WHILE(k_high - k_low > 1)
       k = (k_high + k_low) / 2
       IF (sp%x(k) > x) THEN
          k_high = k
       ELSE
          k_low = k
       END IF
    END DO
  END FUNCTION sp_index_of_x

END MODULE mod_splines
