PROGRAM test_amjuel_extrapolation
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys
  USE physics, ONLY: neutral_rt, compute_Te, compute_2D_eirene_rate, &
       compute_2D_eirene_rate_du, eirene_rate_te_min_phys, eirene_rate_te_max_phys, &
       eirene_rate_ne_min_phys, eirene_rate_ne_max_phys
  IMPLICIT NONE
  REAL*8 :: alpha(9,9), temperatures(5), densities(3), U(5), V(5), jac(5)
  REAL*8 :: te, rate, plus, minus, numerical, analytic, error, worst, cutoff_rate
  REAL*8, PARAMETER :: temperature_scale = 50.d0, density_scale = 1.d19
  REAL*8, PARAMETER :: step = 2.d-6, tolerance = 5.d-7
  INTEGER :: it, id, component, variable, checks

  phys%Mref = 1.25d0
  neutral_rt%eirene_te_min = eirene_rate_te_min_phys/temperature_scale
  neutral_rt%eirene_te_max = eirene_rate_te_max_phys/temperature_scale
  neutral_rt%eirene_ne_min = eirene_rate_ne_min_phys/density_scale
  neutral_rt%eirene_ne_max = eirene_rate_ne_max_phys/density_scale
  ! A smooth log-rate with a density-dependent temperature slope. Numerical
  ! differentiation of the production value routine is the derivative oracle.
  alpha = 0.d0
  alpha(1,1) = -12.d0
  alpha(2,1) = 0.3d0
  alpha(1,2) = 0.2d0
  alpha(2,2) = 0.08d0
  alpha(3,2) = -0.009d0
  alpha(2,3) = 0.004d0
  alpha(4,3) = 0.0002d0
  temperatures = [0.5d0*neutral_rt%eirene_te_min, neutral_rt%eirene_te_min, &
       1.d0/temperature_scale, neutral_rt%eirene_te_max, 2.d0*neutral_rt%eirene_te_max]
  ! Keep density perturbations away from the nonsmooth clamp endpoints.
  densities = [0.5d0*neutral_rt%eirene_ne_min, &
       SQRT(neutral_rt%eirene_ne_min*neutral_rt%eirene_ne_max), &
       2.d0*neutral_rt%eirene_ne_max]
  worst = 0.d0
  checks = 0
  DO it = 1,SIZE(temperatures)
    DO id = 1,SIZE(densities)
      U = [densities(id), 0.d0, densities(id), &
           1.5d0*phys%Mref*densities(id)*temperatures(it), 0.d0]
      CALL compute_Te(U,te)
      CALL compute_2D_eirene_rate(te,U(1),alpha,rate)
      CALL compute_2D_eirene_rate_du(U(1),U(4),te,U(1),alpha,jac)
      IF (.NOT. ieee_is_finite(rate) .OR. rate <= 0.d0) ERROR STOP 'invalid test rate'
      IF (.NOT. ALL(ieee_is_finite(jac))) ERROR STOP 'nonfinite rate derivative'
      DO component = 1,2
        variable = 1
        IF (component == 2) variable = 4
        V = U
        V(variable) = U(variable)*(1.d0+step)
        CALL evaluate(V,alpha,plus)
        V(variable) = U(variable)*(1.d0-step)
        CALL evaluate(V,alpha,minus)
        numerical = (plus-minus)/(2.d0*step*rate)
        analytic = jac(variable)*U(variable)/rate
        error = ABS(numerical-analytic)/MAX(1.d0,ABS(numerical),ABS(analytic))
        worst = MAX(worst,error)
        checks = checks+1
        IF (error > tolerance) THEN
          WRITE(*,*) 'Te [eV], ne [m^-3], variable:',te*temperature_scale,U(1)*density_scale,variable
          WRITE(*,*) 'analytic, finite difference, scaled error:',analytic,numerical,error
          ERROR STOP 'AMJUEL conservative-variable derivative mismatch'
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  ! The existing log-rate cutoff must still give zero value and sensitivity.
  alpha = 0.d0
  alpha(1,1) = -120.d0
  CALL compute_2D_eirene_rate(te,U(1),alpha,cutoff_rate)
  CALL compute_2D_eirene_rate_du(U(1),U(4),te,U(1),alpha,jac)
  IF (cutoff_rate /= 0.d0 .OR. ANY(jac /= 0.d0)) ERROR STOP 'rate cutoff changed'
  WRITE(*,*) 'AMJUEL extrapolation: PASS; derivative checks, worst scaled error:',checks,worst

CONTAINS

  SUBROUTINE evaluate(state,coefficients,value)
    REAL*8, INTENT(IN) :: state(:), coefficients(:,:)
    REAL*8, INTENT(OUT) :: value
    REAL*8 :: temperature
    CALL compute_Te(state,temperature)
    CALL compute_2D_eirene_rate(temperature,state(1),coefficients,value)
  END SUBROUTINE evaluate

END PROGRAM test_amjuel_extrapolation
