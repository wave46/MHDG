PROGRAM test_bohm_momentum
  USE extracted_bohm
  IMPLICIT NONE
  REAL*8 :: u(5), uhat(5), trial(5), ds(5), speed, factor, step
  REAL*8 :: tau(5,5), ju(5,5), juhat(5,5), plus(5), minus(5), expected(5), actual(5)
  REAL*8 :: ju_general(5,5), juhat_general(5,5)
  REAL*8 :: velocities(7), velocity, bn, target, worst
  INTEGER :: matrix_case, delta, sign_case, j, i, k, checks, branch_checks, selected, equivalence_checks
  LOGICAL :: full

  u = [1.2d0,.2d0,2.1d0,3.6d0,.3d0]
  uhat = [.8d0,-.15d0,1.6d0,2.2d0,.1d0]
  worst = 0.d0
  checks = 0
  equivalence_checks = 0
  DO matrix_case = 0,1
    full = matrix_case == 1
    tau = 0.d0
    DO i = 1,5
      tau(i,i) = 1.d0+.2d0*i
      IF (full) THEN
        DO j = 1,5
          IF (i /= j) tau(i,j) = .03d0*(i+j)
        END DO
      ENDIF
    END DO
    DO delta = 0,1
      DO sign_case = -1,1,2
        factor = .65d0*sign_case
        CALL compute_cs(uhat,speed)
        CALL compute_dcs_du(uhat,ds)
        speed = factor*speed
        ds = factor*ds
        CALL assemble(u,uhat,ds,speed,delta,tau,full,ju,juhat)
        IF (.NOT. full) THEN
          ! The general assembly must reduce to the diagonal assembly when
          ! supplied the same diagonal tau, despite different loop indices.
          CALL assemble(u,uhat,ds,speed,delta,tau,.TRUE.,ju_general,juhat_general)
          IF (MAXVAL(ABS(ju-ju_general))>1.d-12 .OR. &
              MAXVAL(ABS(juhat-juhat_general))>1.d-12) ERROR STOP 'diagonal/full tau reduction'
          equivalence_checks = equivalence_checks+1
        ENDIF
        CALL residual(u,uhat,speed,delta,tau,expected)
        actual = MATMUL(ju,u)+MATMUL(juhat,uhat)
        IF (MAXVAL(ABS(actual-expected)) > 1.d-12) ERROR STOP 'trace-density residual reconstruction'
        IF (.NOT. full .AND. delta == 1) THEN
          IF (ABS(ju(2,1)) > 1.d-12) ERROR STOP 'active momentum still uses element density'
        ENDIF
        DO j = 1,5
          step = 1.d-6*MAX(ABS(u(j)),.1d0)
          trial = u; trial(j) = trial(j)+step
          CALL residual(trial,uhat,speed,delta,tau,plus)
          trial = u; trial(j) = trial(j)-step
          CALL residual(trial,uhat,speed,delta,tau,minus)
          CALL check((plus-minus)/(2*step),ju(:,j))
          step = 1.d-6*MAX(ABS(uhat(j)),.1d0)
          trial = uhat; trial(j) = trial(j)+step
          CALL trace_residual(u,trial,factor,delta,tau,plus)
          trial = uhat; trial(j) = trial(j)-step
          CALL trace_residual(u,trial,factor,delta,tau,minus)
          CALL check((plus-minus)/(2*step),juhat(:,j))
        END DO
      END DO
    END DO
  END DO

  ! Constant cs=1 at each state, while changing signed parallel velocity.
  velocities = [-3.d0,-1.d0,-.5d0,0.d0,.5d0,1.d0,3.d0]
  branch_checks = 0
  DO sign_case = -1,1,2
    bn = sign_case*.5d0
    DO k = 1,SIZE(velocities)
      velocity = velocities(k)
      trial = [2.d0,2.d0*velocity,2.d0*(.75d0+.5d0*velocity**2),1.5d0,.1d0]
      CALL select_branch(trial,bn,selected,target)
      delta = 1
      IF (sign_case*velocity > 1.d0+1.d-12) delta = 0
      IF (selected /= delta) ERROR STOP 'signed outgoing Bohm branch'
      IF (delta == 1 .AND. ABS(target-sign_case) > 1.d-12) ERROR STOP 'outgoing sonic target'
      branch_checks = branch_checks+1
    END DO
    ! Grazing incidence keeps its original attenuation and active branch.
    CALL select_branch(trial,sign_case*.05d0,selected,target)
    IF (selected /= 1 .OR. ABS(target-sign_case*.5d0)>1.d-12) ERROR STOP 'grazing branch changed'
    branch_checks = branch_checks+1
  END DO
  ! Preserve the existing reciprocal-suppression switch threshold.
  trial = [2.d0,3.d0,3.75d0,1.5d0,.1d0]
  phys%bohm_energy_thresh = trial(3)
  CALL select_branch(trial,.5d0,selected,target)
  IF (selected /= 1 .OR. ABS(target-.5d0)>1.d-12) ERROR STOP 'suppressed branch changed'
  branch_checks = branch_checks+1
  WRITE(*,*) 'Bohm momentum PASS: derivative columns, branch cases, worst error:',checks,branch_checks,worst
  WRITE(*,*) 'Diagonal/full assembly equivalence cases:',equivalence_checks

CONTAINS

  SUBROUTINE residual(u,uhat,speed,delta,tau,result)
    REAL*8, INTENT(IN) :: u(5),uhat(5),speed,tau(5,5)
    INTEGER, INTENT(IN) :: delta
    REAL*8, INTENT(OUT) :: result(5)
    REAL*8 :: target(5)
    target = u
    target(2) = delta*speed*uhat(1)+(1-delta)*u(2)
    result = MATMUL(tau,target-uhat)
  END SUBROUTINE

  SUBROUTINE trace_residual(u,uhat,factor,delta,tau,result)
    REAL*8, INTENT(IN) :: u(5),uhat(5),factor,tau(5,5)
    INTEGER, INTENT(IN) :: delta
    REAL*8, INTENT(OUT) :: result(5)
    REAL*8 :: speed
    CALL compute_cs(uhat,speed)
    speed = factor*speed
    CALL residual(u,uhat,speed,delta,tau,result)
  END SUBROUTINE

  SUBROUTINE check(numerical,analytic)
    REAL*8, INTENT(IN) :: numerical(5),analytic(5)
    worst = MAX(worst,MAXVAL(ABS(numerical-analytic)))
    checks = checks+1
    IF (MAXVAL(ABS(numerical-analytic))>1.d-8) ERROR STOP 'Bohm momentum Jacobian'
  END SUBROUTINE
END PROGRAM
