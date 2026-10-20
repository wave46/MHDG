PROGRAM test_bdf_coefficients
  USE globals, ONLY: time
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  INTEGER :: order, step, active_order, degree, j
  REAL*8 :: coefficients(7), weights(7), nodes(7), powers(7)
  REAL*8 :: derivative, expected, tolerance

  ! The assembly puts historical coefficients on the RHS:
  ! dt*du/dt = ktis(1)*u_new - SUM(ktis(j+1)*u_previous_j).
  DO order = 1, 6
     time%tis = order
     DO step = 1, order + 1
        time%ik = step
        active_order = MIN(step, order)
        CALL setTimeIntegrationCoefficients(coefficients(1:order+1))
        IF (.NOT. ALL(ieee_is_finite(coefficients(1:order+1)))) ERROR STOP 'nonfinite BDF coefficients'
        IF (ANY(coefficients(active_order+2:order+1) /= 0.d0)) ERROR STOP 'unused startup history'
        IF (active_order == 1) THEN
           IF (ANY(coefficients(1:2) /= 1.d0)) ERROR STOP 'BDF1 compatibility'
        ENDIF

        weights(1) = coefficients(1)
        weights(2:active_order+1) = -coefficients(2:active_order+1)
        DO j = 0, active_order
           nodes(j+1) = -DBLE(j)
        ENDDO
        ! A BDF of order k must differentiate every polynomial of degree <= k
        ! exactly at the endpoint. These identities uniquely determine its weights.
        DO degree = 0, active_order
           powers(1:active_order+1) = 1.d0
           IF (degree > 0) powers(1:active_order+1) = nodes(1:active_order+1)**degree
           derivative = SUM(weights(1:active_order+1)*powers(1:active_order+1))
           expected = 0.d0
           IF (degree == 1) expected = 1.d0
           tolerance = 64.d0*EPSILON(1.d0)*MAX(1.d0, &
                SUM(ABS(weights(1:active_order+1)*powers(1:active_order+1))))
           IF (ABS(derivative-expected) > tolerance) THEN
              WRITE (*,*) 'BDF order, startup step, polynomial degree:', order, step, degree
              WRITE (*,*) 'Derivative, expected:', derivative, expected
              ERROR STOP 'BDF polynomial differentiation'
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  WRITE (*,'(A)') 'BDF1-6 coefficients and startup checks: PASS'

CONTAINS
  ! Extracted verbatim from HDG_computeJacobian; globals supplies its real time state.
  INCLUDE 'bdf_coefficients_routine.F'
END PROGRAM test_bdf_coefficients
