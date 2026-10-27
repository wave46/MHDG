#ifndef TEMPERATURE
#error This check requires a temperature model.
#endif
PROGRAM test_temperature_exchange
  USE globals, ONLY: phys, switch
  USE physics, ONLY: evaluate_temperature_exchange, temperature_exchange_state_is_evaluable
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite, ieee_is_nan
  IMPLICIT NONE
  REAL*8, PARAMETER :: mu=9.109d-31/3.35d-27
  REAL*8, PARAMETER :: fractions(*)=[-4.d0,-.01d0,-mu/2,-1.d-9,0.d0,1.d-9,mu/2,.1d0,.5d0,1.d0,2.d0,5.d0]
  REAL*8, PARAMETER :: densities(*)=[.001d0,1.d0,10.d0], totals(*)=[.01d0,.4d0,2.d0]
  REAL*8, PARAMETER :: velocities(*)=[.3d0,3.d0], scale_factors(*)=[.1d0,7.d0]
  REAL*8 :: temperature_scale, worst=0.d0
  INTEGER :: checks=0
  CHARACTER(64) :: current_setting

  phys%tie=-2.d0
  phys%Mref=1.3d0
#ifdef KEQUATION
  phys%diff_k_min=.01d0
#endif
  temperature_scale=2.d0/(3.d0*phys%Mref)

  CALL check_source_setting('steady source',.FALSE.,54,.08d0)
  CALL check_source_setting('startup above diffusion threshold',.TRUE.,54,.08d0)
  CALL check_source_setting('MMS bypasses startup multiplier',.TRUE.,2,.08d0)
  CALL check_source_setting('startup below diffusion threshold',.TRUE.,54,.01d0)
  CALL check_zero_temperature_join()
  CALL check_invalid_states()
  WRITE(*,*)'Temperature exchange: PASS; FD checks, worst scaled discrepancy:',checks,worst

CONTAINS

  SUBROUTINE check_source_setting(label,startup,testcase,diffusion)
    CHARACTER(*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: startup
    INTEGER, INTENT(IN) :: testcase
    REAL*8, INTENT(IN) :: diffusion
    REAL*8 :: coefficient, kappa
    INTEGER :: density, energy

    current_setting=label
    switch%psdtime=startup
    switch%testcase=testcase
    phys%diff_ee=diffusion
    coefficient=diffusion
#ifdef KEQUATION
    coefficient=coefficient+phys%diff_k_min
#endif
    kappa=-1.d0/phys%tie
    IF(startup.AND.testcase/=2.AND.coefficient>.038d0)kappa=kappa*coefficient/.038d0
    DO density=1,SIZE(densities)
      DO energy=1,SIZE(totals)
        CALL check_temperature_partitions(densities(density),totals(energy),kappa)
      ENDDO
    ENDDO
  END SUBROUTINE check_source_setting

  SUBROUTINE check_temperature_partitions(density,total,kappa)
    REAL*8, INTENT(IN) :: density,total,kappa
    REAL*8 :: electron
    INTEGER :: speed, fraction

    DO speed=1,SIZE(velocities)
      DO fraction=1,SIZE(fractions)
        electron=total*fractions(fraction)
        CALL check_source_state(density,velocities(speed),total-electron,electron,total,kappa)
      ENDDO
    ENDDO
  END SUBROUTINE check_temperature_partitions

  SUBROUTINE check_source_state(density,velocity,ion,electron,total,kappa)
    REAL*8, INTENT(IN) :: density,velocity,ion,electron,total,kappa
    REAL*8 :: U(5), s, ds(5), expected, scale

    U=state(density,velocity,ion,electron)
    CALL evaluate_temperature_exchange(U,s,ds)
    CALL require(temperature_exchange_state_is_evaluable(U),'single-species undershoot rejected',U)
    CALL require(ieee_is_finite(s).AND.ALL(ieee_is_finite(ds)),'nonfinite exchange',U)
    expected=kappa*density**2*(ion-electron)/reference_collision_temperature(ion,electron,total)**1.5d0
    scale=MAX(ABS(s),SUM(ABS(ds*U)),1.d-100)
    CALL require(ABS(s-expected)<=2.d-11*scale,'collision-temperature value',U)
    CALL require(ABS(DOT_PRODUCT(ds,U)-2.d0*s)<=2.d-11*scale,'degree-two homogeneity',U)
    IF(electron<0.d0)CALL require(s>0.d0.AND.ds(4)-ds(3)<0.d0,'electron restoring response',U)
    IF(ion<0.d0)CALL require(s<0.d0,'ion restoring response',U)
    CALL check_derivatives(U,s,ds,density,total,electron)
    CALL check_source_scaling(U,s,scale)
  END SUBROUTINE check_source_state

  FUNCTION reference_collision_temperature(ion,electron,total) RESULT(collision)
    REAL*8, INTENT(IN) :: ion,electron,total
    REAL*8 :: collision, a, exponent, factor

    IF(electron>=0.d0)THEN
      collision=electron+mu*ion
      RETURN
    ENDIF
    a=mu*total
    exponent=2.d0*(1.d0-mu)*electron/a
    factor=0.d0
    IF(exponent>=-700.d0)factor=EXP(exponent)
    collision=.5d0*a*(1.d0+factor)
  END FUNCTION reference_collision_temperature

  SUBROUTINE check_derivatives(U,s,ds,density,total,electron)
    REAL*8, INTENT(IN) :: U(5),s,ds(5),density,total,electron
    REAL*8 :: base_step, reference, error
    INTEGER :: column

    DO column=1,SIZE(U)
      base_step=2.d-5*MAX(ABS(U(column)),density*total/temperature_scale)
      IF(column==4)base_step=2.d-5*density*MAX(ABS(electron),mu*total)/temperature_scale
      reference=MAX(ABS(ds(column)),ABS(s)/MAX(ABS(U(column)),density*total/temperature_scale),1.d-100)
      error=derivative_error(U,column,ds(column),base_step,reference)
      worst=MAX(worst,error)
      IF(error>3.d-5)WRITE(*,*)'Derivative column, scaled discrepancy:',column,error
      CALL require(error<=3.d-5,'exchange finite difference',U)
      checks=checks+1
    ENDDO
  END SUBROUTINE check_derivatives

  FUNCTION derivative_error(U,column,analytic,base_step,reference) RESULT(error)
    REAL*8, INTENT(IN) :: U(5),analytic,base_step,reference
    INTEGER, INTENT(IN) :: column
    REAL*8 :: error, V(5), step, sp, sm, numerical
    INTEGER :: attempt

    ! Near cancellation, several steps distinguish truncation from roundoff.
    error=HUGE(1.d0)
    DO attempt=0,3
      step=base_step*10.d0**(-attempt)
      V=U
      V(column)=U(column)+step
      CALL evaluate_temperature_exchange(V,sp)
      V(column)=U(column)-step
      CALL evaluate_temperature_exchange(V,sm)
      numerical=(sp-sm)/(2.d0*step)
      error=MIN(error,ABS(numerical-analytic)/reference)
    ENDDO
  END FUNCTION derivative_error

  SUBROUTINE check_source_scaling(U,s,scale)
    REAL*8, INTENT(IN) :: U(5),s,scale
    REAL*8 :: scaled_source, factor
    INTEGER :: i

    DO i=1,SIZE(scale_factors)
      factor=scale_factors(i)
      CALL evaluate_temperature_exchange(factor*U,scaled_source)
      CALL require(ABS(scaled_source-factor**2*s)<=2.d-11*factor**2*scale,'source scaling',U)
    ENDDO
  END SUBROUTINE check_source_scaling

  SUBROUTINE check_zero_temperature_join()
    REAL*8 :: U(5), s, ds(5), left_s, right_s, left_ds(5), right_ds(5)

    current_setting='join at zero electron temperature'
    U=state(1.d0,.3d0,1.d0,0.d0)
    CALL evaluate_temperature_exchange(U,s,ds)
    CALL evaluate_temperature_exchange(state(1.d0,.3d0,1.d0,-mu*1.d-8),left_s,left_ds)
    CALL evaluate_temperature_exchange(state(1.d0,.3d0,1.d0, mu*1.d-8),right_s,right_ds)
    CALL require(MAX(ABS(left_s-s),ABS(right_s-s))<=2.d-7*ABS(s),'value continuity',U)
    CALL require(MAX(MAXVAL(ABS(left_ds-ds)),MAXVAL(ABS(right_ds-ds)))<=2.d-7*MAXVAL(ABS(ds)), &
         'first derivative continuity',U)
  END SUBROUTINE check_zero_temperature_join

  SUBROUTINE check_invalid_states()
    REAL*8 :: U(5)

    current_setting='invalid trial states'
    CALL expect_rejected(state(0.d0,0.d0,1.d0,1.d0),'zero density')
    CALL expect_rejected(state(-1.d0,0.d0,1.d0,1.d0),'negative density')
    U=state(1.d0,0.d0,1.d0,1.d0)
    U(3:4)=-1.d0
    CALL expect_rejected(U,'negative conserved energies')
    CALL expect_rejected(state(1.d0,0.d0,0.d0,0.d0),'zero thermal energy')
    CALL expect_rejected(state(1.d0,3.d0,-.2d0,-.3d0),'both species negative despite positive ion total energy')
    CALL expect_rejected(state(1.d0,0.d0,-1.d0,.2d0),'negative thermal sum from ions')
    CALL expect_rejected(state(1.d0,0.d0,.2d0,-1.d0),'negative thermal sum from electrons')
  END SUBROUTINE check_invalid_states

  SUBROUTINE expect_rejected(U,label)
    REAL*8, INTENT(IN) :: U(5)
    CHARACTER(*), INTENT(IN) :: label
    REAL*8 :: s, ds(5)

    CALL require(.NOT.temperature_exchange_state_is_evaluable(U),label,U)
    CALL evaluate_temperature_exchange(U,s,ds)
    CALL require(ieee_is_nan(s).AND.ALL(ieee_is_nan(ds)),'invalid source must be explicit: '//label,U)
  END SUBROUTINE expect_rejected

  SUBROUTINE require(condition,label,U)
    LOGICAL, INTENT(IN) :: condition
    CHARACTER(*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: U(5)

    IF(condition)RETURN
    WRITE(*,*)TRIM(current_setting),': ',label
    WRITE(*,*)'Conservative state:',U
    ERROR STOP 'temperature-exchange check failed'
  END SUBROUTINE require

  FUNCTION state(density,velocity,ion,electron) RESULT(U)
    REAL*8, INTENT(IN) :: density,velocity,ion,electron
    REAL*8 :: U(5)

    U=[density,density*velocity,density*(ion/temperature_scale+.5d0*velocity**2), &
         density*electron/temperature_scale,.02d0*density]
  END FUNCTION state
END PROGRAM test_temperature_exchange
