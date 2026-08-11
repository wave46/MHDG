MODULE neutral_flux_limiter
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: neutral_tn_source_invalid = 0
  INTEGER, PARAMETER, PUBLIC :: neutral_tn_source_ti = 1
  INTEGER, PARAMETER, PUBLIC :: neutral_tn_source_fixed = 2

  TYPE, PUBLIC :: neutral_flux_limiter_config_t
     INTEGER :: tn_source = neutral_tn_source_invalid
     REAL*8 :: fixed_tn = 0.d0
     REAL*8 :: mref = 0.d0
     REAL*8 :: epsilon = 0.d0
     REAL*8 :: fs_fraction = 1.d0
     REAL*8 :: fs_flux_min = 0.d0
     REAL*8 :: diff_nn_min = 0.d0
  END TYPE neutral_flux_limiter_config_t

  TYPE, PUBLIC :: neutral_flux_limiter_result_t
     REAL*8 :: phi = 1.d0
     REAL*8 :: flux_cap = 0.d0
     REAL*8 :: unlimited_flux_norm = 0.d0
     REAL*8 :: activation_ratio = 0.d0
  END TYPE neutral_flux_limiter_result_t

  PUBLIC :: parse_neutral_tn_source
  PUBLIC :: neutral_free_streaming_speed
  PUBLIC :: compute_neutral_unlimited_flux
  PUBLIC :: apply_neutral_perpendicular_operator
  PUBLIC :: evaluate_neutral_flux_limiter

CONTAINS

  PURE INTEGER FUNCTION parse_neutral_tn_source(source) RESULT(source_id)
    CHARACTER(LEN=*), INTENT(IN) :: source

    SELECT CASE (TRIM(ADJUSTL(source)))
    CASE ('ti')
       source_id = neutral_tn_source_ti
    CASE ('fixed')
       source_id = neutral_tn_source_fixed
    CASE DEFAULT
       source_id = neutral_tn_source_invalid
    END SELECT
  END FUNCTION parse_neutral_tn_source

  PURE REAL*8 FUNCTION neutral_free_streaming_speed(config, ti) &
       RESULT(speed)
    TYPE(neutral_flux_limiter_config_t), INTENT(IN) :: config
    REAL*8, INTENT(IN) :: ti
    REAL*8 :: cap_temperature

    SELECT CASE (config%tn_source)
    CASE (neutral_tn_source_ti)
       cap_temperature = ti
    CASE (neutral_tn_source_fixed)
       cap_temperature = config%fixed_tn
    CASE DEFAULT
       cap_temperature = 0.d0
    END SELECT
    speed = SQRT(config%mref*cap_temperature)
  END FUNCTION neutral_free_streaming_speed

  PURE SUBROUTINE compute_neutral_unlimited_flux(dnn, neutral_gradient, &
       pressure_flux, unlimited_flux)
    REAL*8, INTENT(IN) :: dnn
    REAL*8, INTENT(IN) :: neutral_gradient(:), pressure_flux(:)
    REAL*8, INTENT(OUT) :: unlimited_flux(:)

    unlimited_flux = -dnn*neutral_gradient - pressure_flux
  END SUBROUTINE compute_neutral_unlimited_flux

  PURE SUBROUTINE apply_neutral_perpendicular_operator(flux, &
       magnetic_direction)
    REAL*8, INTENT(INOUT) :: flux(:)
    REAL*8, INTENT(IN) :: magnetic_direction(:)

    flux = flux - DOT_PRODUCT(flux, magnetic_direction)*magnetic_direction
  END SUBROUTINE apply_neutral_perpendicular_operator

  PURE SUBROUTINE evaluate_neutral_flux_limiter(config, ti, &
       neutral_density, dnn, unlimited_flux, result)
    TYPE(neutral_flux_limiter_config_t), INTENT(IN) :: config
    TYPE(neutral_flux_limiter_result_t), INTENT(OUT) :: result
    REAL*8, INTENT(IN) :: ti, neutral_density, dnn
    REAL*8, INTENT(IN) :: unlimited_flux(:)
    REAL*8 :: free_streaming_speed, gamma_regularized
    REAL*8 :: phi_diffusion_floor

    free_streaming_speed = neutral_free_streaming_speed(config, ti)
    result%flux_cap = MAX(config%fs_fraction*MAX(neutral_density, 0.d0)* &
         free_streaming_speed, config%fs_flux_min)
    result%unlimited_flux_norm = &
         SQRT(DOT_PRODUCT(unlimited_flux, unlimited_flux))
    gamma_regularized = SQRT(result%unlimited_flux_norm**2 + &
         config%epsilon**2)

    IF (result%flux_cap > 0.d0) THEN
       result%activation_ratio = gamma_regularized/result%flux_cap
    ELSEIF (gamma_regularized > 0.d0) THEN
       result%activation_ratio = HUGE(1.d0)
    ELSE
       result%activation_ratio = 0.d0
    ENDIF

    phi_diffusion_floor = MIN(1.d0, config%diff_nn_min/dnn)
    result%phi = MAX(1.d0/(1.d0 + result%activation_ratio), &
         phi_diffusion_floor)
  END SUBROUTINE evaluate_neutral_flux_limiter

END MODULE neutral_flux_limiter
