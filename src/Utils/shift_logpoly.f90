MODULE shift_logpoly
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shift_logpoly_1d_5
  PUBLIC :: shift_logpoly_1d_9
  PUBLIC :: shift_logpoly_1d_17
  PUBLIC :: shift_logpoly_2d_9x9

CONTAINS

  SUBROUTINE shift_logpoly_2d_9x9(alpha, shift_x, shift_y, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(9,9)
    REAL*8, INTENT(IN)    :: shift_x, shift_y, log_scale

    CALL shift_logpoly_2d(alpha, shift_x, shift_y, log_scale)
  ENDSUBROUTINE shift_logpoly_2d_9x9

  SUBROUTINE shift_logpoly_1d_9(alpha, shift_x, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(9)
    REAL*8, INTENT(IN)    :: shift_x, log_scale

    CALL shift_logpoly_1d(alpha, shift_x, log_scale)
  ENDSUBROUTINE shift_logpoly_1d_9

  SUBROUTINE shift_logpoly_1d_5(alpha, shift_x, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(5)
    REAL*8, INTENT(IN)    :: shift_x, log_scale

    CALL shift_logpoly_1d(alpha, shift_x, log_scale)
  ENDSUBROUTINE shift_logpoly_1d_5

  SUBROUTINE shift_logpoly_1d_17(alpha, shift_x, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(17)
    REAL*8, INTENT(IN)    :: shift_x, log_scale

    CALL shift_logpoly_1d(alpha, shift_x, log_scale)
  ENDSUBROUTINE shift_logpoly_1d_17

  SUBROUTINE shift_logpoly_2d(alpha, shift_x, shift_y, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(:,:)
    REAL*8, INTENT(IN)    :: shift_x, shift_y, log_scale
    REAL*8                :: shifted(SIZE(alpha,1), SIZE(alpha,2))
    REAL*8                :: x_factor, y_factor
    INTEGER               :: i, j, ip, jp, pow_x, pow_y

    shifted = 0.d0
    DO j = 1, SIZE(alpha,2)
      DO i = 1, SIZE(alpha,1)
        DO jp = 1, j
          pow_y = (j - jp)
          y_factor = 1.d0
          IF (pow_y > 0) y_factor = shift_y**pow_y
          DO ip = 1, i
            pow_x = (i - ip)
            x_factor = 1.d0
            IF (pow_x > 0) x_factor = shift_x**pow_x
            shifted(ip,jp) = shifted(ip,jp) + alpha(i,j)* &
              &binomial_coefficient(i - 1, ip - 1)*binomial_coefficient(j - 1, jp - 1)* &
              &x_factor*y_factor
          END DO
        END DO
      END DO
    END DO

    shifted(1,1) = shifted(1,1) + log_scale
    alpha = shifted
  ENDSUBROUTINE shift_logpoly_2d

  SUBROUTINE shift_logpoly_1d(alpha, shift_x, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(:)
    REAL*8, INTENT(IN)    :: shift_x, log_scale
    REAL*8                :: shifted(SIZE(alpha))
    REAL*8                :: x_factor
    INTEGER               :: i, ip, pow_x

    shifted = 0.d0
    DO i = 1, SIZE(alpha)
      DO ip = 1, i
        pow_x = (i - ip)
        x_factor = 1.d0
        IF (pow_x > 0) x_factor = shift_x**pow_x
        shifted(ip) = shifted(ip) + alpha(i)*binomial_coefficient(i - 1, ip - 1)*x_factor
      END DO
    END DO

    shifted(1) = shifted(1) + log_scale
    alpha = shifted
  ENDSUBROUTINE shift_logpoly_1d

  REAL*8 FUNCTION binomial_coefficient(n, k)
    INTEGER, INTENT(IN) :: n, k
    INTEGER             :: i, kk

    IF (k < 0 .OR. k > n) THEN
      binomial_coefficient = 0.d0
      RETURN
    END IF

    IF (k == 0 .OR. k == n) THEN
      binomial_coefficient = 1.d0
      RETURN
    END IF

    kk = MIN(k, n - k)
    binomial_coefficient = 1.d0
    DO i = 1, kk
      binomial_coefficient = binomial_coefficient*DBLE(n - kk + i)/DBLE(i)
    END DO
  ENDFUNCTION binomial_coefficient

END MODULE shift_logpoly
