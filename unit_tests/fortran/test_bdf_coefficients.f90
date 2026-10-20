PROGRAM test_bdf_coefficients
  USE globals, ONLY: time, sol
  USE Main_utils, ONLY: initialize_solu0_uiter_uconv, update_solution, update_u0, &
       uiter, uiter_best, qiter_best, global_step => it
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  INTEGER :: order, step, active_order, degree, j
  REAL*8 :: coefficients(7), weights(7), nodes(7), powers(7)
  REAL*8 :: derivative, expected, tolerance

  ! The assembly puts historical coefficients on the RHS:
  ! dt*du/dt = ktis(1)*u_new - SUM(ktis(j+1)*u_previous_j).
  time%bdf_start_step = 0
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
  CALL check_history_startup(0)
  CALL check_history_startup(50)
  WRITE (*,'(A)') 'BDF1-6 coefficients, fresh/restart history and reseeding: PASS'

CONTAINS
  SUBROUTINE check_history_startup(saved_step)
    INTEGER, INTENT(IN) :: saved_step
    INTEGER :: requested, advance, available
    REAL*8 :: weights(7)

    DO requested = 1,6
       time%tis = requested
       time%it = saved_step
       time%ik = saved_step
       time%t = 3.d0
       sol%Nt = saved_step
       ALLOCATE(sol%u(2),sol%q(4),sol%u_tilde(2))
       sol%u = 7.d0
       sol%q = 0.d0
       sol%u_tilde = 7.d0
       CALL initialize_solu0_uiter_uconv()
       IF (time%it /= saved_step .OR. time%ik /= saved_step .OR. &
           sol%Nt /= saved_step .OR. time%t /= 3.d0) ERROR STOP 'history initialization changed clock'
       ! Unavailable columns must never enter the derivative, even if nonzero.
       sol%u0(:,2:) = -123.d0
       DO advance = 1,requested+1
          time%ik = saved_step+advance
          global_step = saved_step+advance
          available = MIN(requested,advance)
          CALL setTimeIntegrationCoefficients(weights(1:requested+1))
          IF (weights(available+1) == 0.d0 .OR. ANY(weights(available+2:requested+1) /= 0.d0)) &
               ERROR STOP 'history startup order'
          IF (MAXVAL(ABS(weights(1)*sol%u-MATMUL(sol%u0,weights(2:requested+1)))) > 1.d-12) &
               ERROR STOP 'restart created a derivative of a constant state'
          CALL update_solution()
       ENDDO
       CALL update_u0(sol%u)
       IF (time%bdf_start_step /= time%ik) ERROR STOP 'single-state reseed retained old startup'
       time%ik = time%ik+1
       CALL setTimeIntegrationCoefficients(weights(1:requested+1))
       IF (ANY(weights(1:2) /= 1.d0) .OR. ANY(weights(3:requested+1) /= 0.d0)) &
            ERROR STOP 'single-state reseed did not return to BDF1'
       DEALLOCATE(sol%u,sol%q,sol%u_tilde,sol%u0,sol%u_conv,sol%q_conv,uiter,uiter_best,qiter_best)
       NULLIFY(sol%u,sol%q,sol%u_tilde,sol%u_conv,sol%q_conv,uiter_best,qiter_best)
    ENDDO
  END SUBROUTINE check_history_startup

  ! Extracted verbatim from HDG_computeJacobian; globals supplies its real time state.
  INCLUDE 'bdf_coefficients_routine.F'
END PROGRAM test_bdf_coefficients
