PROGRAM test_transport_taper
  USE transport_models_1d_config, ONLY: transport_model_config_t, &
       tm1d_config_reset, tm1d_config_apply, tm1d_particle_taper_factor
  IMPLICIT NONE

  TYPE(transport_model_config_t) :: config
  REAL*8 :: untapered, tapered, clipped

  CALL assert_close(config%c_bohm_n_rho_slope, 0.7d0, 'initialized default slope')
  CALL tm1d_config_reset(config)
  CALL assert_close(config%c_bohm_n_rho_slope, 0.7d0, 'reset default slope')

  CALL tm1d_config_apply(config, 1.d0, 1.d0, c_bohm_n_rho_slope=0.d0)
  CALL assert_close(config%c_bohm_n_rho_slope, 0.d0, 'configured zero slope')

  untapered = tm1d_particle_taper_factor(1.2d0, 0.d0)
  tapered = tm1d_particle_taper_factor(1.2d0, 0.7d0)
  clipped = tm1d_particle_taper_factor(2.d0, 0.7d0)
  CALL assert_close(untapered, 1.d0, 'zero slope is exactly untapered')
  CALL assert_close(tapered, 0.16d0, 'default slope at rho above one')
  CALL assert_close(clipped, 0.d0, 'particle taper lower-bound clipping')

  WRITE (*, '(A)') 'transport particle taper checks: PASS'

CONTAINS

  SUBROUTINE assert_close(value, expected, label)
    REAL*8, INTENT(IN) :: value, expected
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (ABS(value - expected) > 1.d-12) THEN
       WRITE (*, '(A,1X,A,2(1X,ES12.4))') 'FAIL:', TRIM(label), value, expected
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_close

END PROGRAM test_transport_taper
