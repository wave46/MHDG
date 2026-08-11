PROGRAM test_neutral_flux_limiter
  USE neutral_flux_limiter, ONLY: neutral_tn_source_invalid, &
       neutral_tn_source_ti, neutral_tn_source_fixed, &
       neutral_flux_limiter_config_t, neutral_flux_limiter_result_t, &
       parse_neutral_tn_source, neutral_free_streaming_speed, &
       compute_neutral_unlimited_flux, &
       apply_neutral_perpendicular_operator, evaluate_neutral_flux_limiter
  IMPLICIT NONE

  CALL test_temperature_sources()
  CALL test_complete_flux()
  CALL test_limiter_algebra()
  WRITE (*, '(A)') 'neutral flux limiter kernel checks: PASS'

CONTAINS

  SUBROUTINE test_temperature_sources()
    TYPE(neutral_flux_limiter_config_t) :: config
    REAL*8 :: speed

    CALL assert_equal(parse_neutral_tn_source('ti'), neutral_tn_source_ti, &
         'Ti source parsing')
    CALL assert_equal(parse_neutral_tn_source('  fixed  '), &
         neutral_tn_source_fixed, 'fixed source parsing')
    CALL assert_equal(parse_neutral_tn_source('electron'), &
         neutral_tn_source_invalid, 'invalid source parsing')

    config%mref = 2.25d0
    config%fixed_tn = 9.d0
    config%tn_source = neutral_tn_source_ti
    speed = neutral_free_streaming_speed(config, 4.d0)
    CALL assert_close(speed, 3.d0, 'Ti cap speed retains Mref')
    config%tn_source = neutral_tn_source_fixed
    speed = neutral_free_streaming_speed(config, 4.d0)
    CALL assert_close(speed, 4.5d0, 'fixed-Tn cap speed retains Mref')
  END SUBROUTINE test_temperature_sources

  SUBROUTINE test_complete_flux()
    REAL*8 :: gradient(2), pressure_flux(2), magnetic_direction(2)
    REAL*8 :: unlimited_flux(2), parallel_flux

    gradient = (/1.d0, -2.d0/)
    pressure_flux = (/0.5d0, 1.5d0/)
    magnetic_direction = (/0.3d0, 0.4d0/)

    CALL compute_neutral_unlimited_flux(dnn=2.d0, &
         neutral_gradient=gradient, pressure_flux=pressure_flux, &
         unlimited_flux=unlimited_flux)
    CALL assert_vector_close(unlimited_flux, (/-2.5d0, 2.5d0/), &
         'complete diffusion and pressure flux')

    parallel_flux = DOT_PRODUCT(unlimited_flux, magnetic_direction)
    CALL apply_neutral_perpendicular_operator(unlimited_flux, &
         magnetic_direction)
    CALL assert_vector_close(unlimited_flux, (/-2.575d0, 2.4d0/), &
         'poloidal components after perpendicular operator')
    CALL assert_close(parallel_flux, 0.25d0, &
         'operator uses the complete unlimited flux')
  END SUBROUTINE test_complete_flux

  SUBROUTINE test_limiter_algebra()
    TYPE(neutral_flux_limiter_config_t) :: config
    TYPE(neutral_flux_limiter_result_t) :: result
    REAL*8 :: unlimited_flux(2)

    config%tn_source = neutral_tn_source_ti
    config%mref = 1.d0
    config%fs_fraction = 0.5d0
    config%fs_flux_min = 0.d0
    config%epsilon = 0.d0
    config%diff_nn_min = 1.d0
    unlimited_flux = (/3.d0, 4.d0/)
    CALL evaluate_neutral_flux_limiter(config=config, ti=100.d0, &
         neutral_density=2.d0, dnn=8.d0, &
         unlimited_flux=unlimited_flux, result=result)
    CALL assert_close(result%flux_cap, 10.d0, &
         'density free-streaming cap')
    CALL assert_close(result%unlimited_flux_norm, 5.d0, &
         'unregularized complete-flux norm')
    CALL assert_close(result%activation_ratio, 0.5d0, &
         'unregularized activation ratio')
    CALL assert_close(result%phi, 2.d0/3.d0, 'gamma-one limiter law')

    config%epsilon = 2.d0
    unlimited_flux = 0.d0
    CALL evaluate_neutral_flux_limiter(config=config, ti=100.d0, &
         neutral_density=2.d0, dnn=8.d0, &
         unlimited_flux=unlimited_flux, result=result)
    CALL assert_close(result%activation_ratio, 0.2d0, &
         'epsilon regularization')
    CALL assert_close(result%phi, 5.d0/6.d0, &
         'regularized limiter factor')

    config%fs_fraction = 0.d0
    config%fs_flux_min = 4.d0
    config%epsilon = 0.d0
    unlimited_flux = (/2.d0, 0.d0/)
    CALL evaluate_neutral_flux_limiter(config=config, ti=100.d0, &
         neutral_density=-3.d0, dnn=8.d0, &
         unlimited_flux=unlimited_flux, result=result)
    CALL assert_close(result%flux_cap, 4.d0, 'minimum cap flux')
    CALL assert_close(result%activation_ratio, 0.5d0, &
         'minimum-cap activation ratio')

    config%fs_fraction = 1.d0
    config%fs_flux_min = 0.d0
    config%diff_nn_min = 2.d0
    unlimited_flux = (/9.d0, 0.d0/)
    CALL evaluate_neutral_flux_limiter(config=config, ti=1.d0, &
         neutral_density=1.d0, dnn=4.d0, &
         unlimited_flux=unlimited_flux, result=result)
    CALL assert_close(result%phi, 0.5d0, 'diffusion floor')
    CALL assert_true(result%phi > 0.d0 .AND. result%phi <= 1.d0, &
         'limiter bounds')
  END SUBROUTINE test_limiter_algebra

  SUBROUTINE assert_equal(value, expected, label)
    INTEGER, INTENT(IN) :: value, expected
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (value /= expected) THEN
       WRITE (*, '(A,1X,A,2(1X,I0))') 'FAIL:', TRIM(label), value, expected
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_equal

  SUBROUTINE assert_close(value, expected, label)
    REAL*8, INTENT(IN) :: value, expected
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (ABS(value - expected) > 1.d-12) THEN
       WRITE (*, '(A,1X,A,2(1X,ES12.4))') 'FAIL:', TRIM(label), value, expected
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_close

  SUBROUTINE assert_vector_close(value, expected, label)
    REAL*8, INTENT(IN) :: value(:), expected(:)
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (MAXVAL(ABS(value - expected)) > 1.d-12) THEN
       WRITE (*, '(A,1X,A)') 'FAIL:', TRIM(label)
       WRITE (*, '(A,*(1X,ES12.4))') 'value:', value
       WRITE (*, '(A,*(1X,ES12.4))') 'expected:', expected
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_vector_close

  SUBROUTINE assert_true(condition, label)
    LOGICAL, INTENT(IN) :: condition
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (.NOT. condition) THEN
       WRITE (*, '(A,1X,A)') 'FAIL:', TRIM(label)
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_true

END PROGRAM test_neutral_flux_limiter
