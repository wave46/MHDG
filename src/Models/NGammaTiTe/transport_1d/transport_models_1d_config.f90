MODULE transport_models_1d_config
  USE magnetic_topology, ONLY: magnetic_region_core, magnetic_region_main_sol
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: transport_model_config_t, tm1d_config_reset, tm1d_config_apply
  PUBLIC :: tm1d_region_policy_from_name, tm1d_region_policy_name
  PUBLIC :: tm1d_region_is_included, tm1d_build_pinch_velocity
  PUBLIC :: transport_region_legacy_all_regions
  PUBLIC :: transport_region_core_and_main_sol, transport_region_core_only

  REAL*8, PARAMETER :: rho_edge_default = 0.99d0
  REAL*8, PARAMETER :: rho_core_default = 0.8d0
  INTEGER, PARAMETER :: transport_region_legacy_all_regions = 0
  INTEGER, PARAMETER :: transport_region_core_and_main_sol = 1
  INTEGER, PARAMETER :: transport_region_core_only = 2
  CHARACTER(LEN=32), PARAMETER :: transport_region_policy_names(0:2) = &
       [CHARACTER(LEN=32) :: 'legacy_all_regions', 'core_and_main_sol', &
       'core_only']

  TYPE :: transport_model_config_t
     REAL*8 :: rho_edge = rho_edge_default
     REAL*8 :: rho_core = rho_core_default
     REAL*8 :: rho_diffusion_model_max = 1.d0
     INTEGER :: transport_region_policy = transport_region_legacy_all_regions
     INTEGER :: pinch_model = 1
     REAL*8 :: c_pinch = 0.5d0
     REAL*8 :: nu_th = 0.04d0
     REAL*8 :: vpinch_const = 0.d0
     REAL*8 :: rho_pinch_axis_width = 0.02d0
     REAL*8 :: rho_pinch_model_max = 0.99d0
     REAL*8 :: rho_pinch_edge_width = 0.03d0
     REAL*8 :: rho_blend_width = 0.02d0
     REAL*8 :: diff_n_min = 0.d0
     REAL*8 :: diff_u_min = 0.d0
     REAL*8 :: diff_e_min = 0.d0
     REAL*8 :: diff_ee_min = 0.d0
     REAL*8 :: c_bohm_i = 1.6d-4
     REAL*8 :: c_gyrobohm_i = 1.75d-2
     REAL*8 :: c_bohm_e = 8.d-5
     REAL*8 :: c_gyrobohm_e = 3.5d-2
     REAL*8 :: c_bohm_n = 1.d0
     REAL*8 :: prandtl = 1.d0
  END TYPE transport_model_config_t

CONTAINS

  SUBROUTINE tm1d_config_reset(config)
    TYPE(transport_model_config_t), INTENT(INOUT) :: config

    config%rho_edge = rho_edge_default
    config%rho_core = rho_core_default
    config%rho_diffusion_model_max = 1.d0
    config%transport_region_policy = transport_region_legacy_all_regions
    config%pinch_model = 1
    config%c_pinch = 0.5d0
    config%nu_th = 0.04d0
    config%vpinch_const = 0.d0
    config%rho_pinch_axis_width = 0.02d0
    config%rho_pinch_model_max = 0.99d0
    config%rho_pinch_edge_width = 0.03d0
    config%rho_blend_width = 0.02d0
    config%diff_n_min = 0.d0
    config%diff_u_min = 0.d0
    config%diff_e_min = 0.d0
    config%diff_ee_min = 0.d0
    config%c_bohm_i = 1.6d-4
    config%c_gyrobohm_i = 1.75d-2
    config%c_bohm_e = 8.d-5
    config%c_gyrobohm_e = 3.5d-2
    config%c_bohm_n = 1.d0
    config%prandtl = 1.d0
  END SUBROUTINE tm1d_config_reset

  SUBROUTINE tm1d_config_apply(config, refval_time, refval_length, rho_edge, &
       rho_core, rho_diffusion_model_max, transport_region_policy, c_bohm_i, &
       c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, prandtl, pinch_model, &
       c_pinch, nu_th, vpinch_const_phys, rho_pinch_axis_width, &
       rho_pinch_model_max, rho_pinch_edge_width, rho_blend_width, &
       diff_n_min_phys, diff_u_min_phys, diff_e_min_phys, diff_ee_min_phys)
    TYPE(transport_model_config_t), INTENT(INOUT) :: config
    REAL*8, INTENT(IN) :: refval_time, refval_length
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core, rho_diffusion_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, prandtl, c_pinch, nu_th, vpinch_const_phys, rho_pinch_axis_width, rho_pinch_model_max, rho_pinch_edge_width, rho_blend_width, diff_n_min_phys, diff_u_min_phys, diff_e_min_phys, diff_ee_min_phys
    INTEGER, INTENT(IN), OPTIONAL :: pinch_model
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: transport_region_policy
    INTEGER :: policy
    LOGICAL :: valid

    IF (PRESENT(rho_edge)) config%rho_edge = rho_edge
    IF (PRESENT(rho_core)) config%rho_core = rho_core
    IF (PRESENT(rho_diffusion_model_max)) config%rho_diffusion_model_max = rho_diffusion_model_max
    IF (PRESENT(transport_region_policy)) THEN
       CALL tm1d_region_policy_from_name(transport_region_policy, policy, valid)
       IF (.NOT. valid) ERROR STOP 'Invalid transport_region_policy'
       config%transport_region_policy = policy
    ENDIF
    IF (PRESENT(c_bohm_i)) config%c_bohm_i = c_bohm_i
    IF (PRESENT(c_gyrobohm_i)) config%c_gyrobohm_i = c_gyrobohm_i
    IF (PRESENT(c_bohm_e)) config%c_bohm_e = c_bohm_e
    IF (PRESENT(c_gyrobohm_e)) config%c_gyrobohm_e = c_gyrobohm_e
    IF (PRESENT(c_bohm_n)) config%c_bohm_n = c_bohm_n
    IF (PRESENT(prandtl)) config%prandtl = prandtl
    IF (PRESENT(pinch_model)) config%pinch_model = pinch_model
    IF (PRESENT(c_pinch)) config%c_pinch = c_pinch
    IF (PRESENT(nu_th)) config%nu_th = nu_th
    IF (PRESENT(vpinch_const_phys)) config%vpinch_const = vpinch_const_phys*refval_time/refval_length
    IF (PRESENT(rho_pinch_axis_width)) config%rho_pinch_axis_width = rho_pinch_axis_width
    IF (PRESENT(rho_pinch_model_max)) config%rho_pinch_model_max = rho_pinch_model_max
    IF (PRESENT(rho_pinch_edge_width)) config%rho_pinch_edge_width = rho_pinch_edge_width
    IF (PRESENT(rho_blend_width)) config%rho_blend_width = rho_blend_width
    IF (PRESENT(diff_n_min_phys)) config%diff_n_min = diff_n_min_phys*refval_time/refval_length**2
    IF (PRESENT(diff_u_min_phys)) config%diff_u_min = diff_u_min_phys*refval_time/refval_length**2
    IF (PRESENT(diff_e_min_phys)) config%diff_e_min = diff_e_min_phys*refval_time/refval_length**2
    IF (PRESENT(diff_ee_min_phys)) config%diff_ee_min = diff_ee_min_phys*refval_time/refval_length**2
  END SUBROUTINE tm1d_config_apply

  SUBROUTINE tm1d_region_policy_from_name(name, policy, valid)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: policy
    LOGICAL, INTENT(OUT) :: valid
    CHARACTER(LEN=32) :: normalized
    INTEGER :: i, code, candidate

    normalized = ADJUSTL(name)
    DO i = 1, LEN_TRIM(normalized)
       code = IACHAR(normalized(i:i))
       IF (code >= IACHAR('A') .AND. code <= IACHAR('Z')) THEN
          normalized(i:i) = ACHAR(code + IACHAR('a') - IACHAR('A'))
       ENDIF
    ENDDO

    policy = transport_region_legacy_all_regions
    valid = .FALSE.
    DO candidate = LBOUND(transport_region_policy_names, 1), &
         UBOUND(transport_region_policy_names, 1)
       IF (TRIM(normalized) /= TRIM(transport_region_policy_names(candidate))) CYCLE
       policy = candidate
       valid = .TRUE.
       RETURN
    ENDDO
  END SUBROUTINE tm1d_region_policy_from_name

  CHARACTER(LEN=32) FUNCTION tm1d_region_policy_name(policy)
    INTEGER, INTENT(IN) :: policy

    IF (policy < LBOUND(transport_region_policy_names, 1) .OR. &
         policy > UBOUND(transport_region_policy_names, 1)) THEN
       tm1d_region_policy_name = 'unknown'
       RETURN
    ENDIF
    tm1d_region_policy_name = transport_region_policy_names(policy)
  END FUNCTION tm1d_region_policy_name

  LOGICAL FUNCTION tm1d_region_is_included(policy, region)
    INTEGER, INTENT(IN) :: policy, region

    SELECT CASE (policy)
    CASE (transport_region_legacy_all_regions)
       tm1d_region_is_included = .TRUE.
    CASE (transport_region_core_and_main_sol)
       tm1d_region_is_included = region == magnetic_region_core .OR. &
            region == magnetic_region_main_sol
    CASE (transport_region_core_only)
       tm1d_region_is_included = region == magnetic_region_core
    CASE DEFAULT
       tm1d_region_is_included = .FALSE.
    END SELECT
  END FUNCTION tm1d_region_is_included

  SUBROUTINE tm1d_build_pinch_velocity(vpinch, policy, b, outward_normal, velocity)
    REAL*8, INTENT(IN) :: vpinch, b(:), outward_normal(:)
    INTEGER, INTENT(IN) :: policy
    REAL*8, INTENT(OUT) :: velocity(2)
    REAL*8 :: direction(2), direction_norm

    velocity = 0.d0
    direction = 0.d0
    IF (SIZE(outward_normal) >= 2) direction = outward_normal(1:2)
    direction_norm = NORM2(direction)
    IF (direction_norm <= 1.d-12) THEN
       IF (policy /= transport_region_legacy_all_regions) RETURN
       IF (SIZE(b) < 2) RETURN
       direction = [b(2), -b(1)]
       direction_norm = NORM2(direction)
    ENDIF
    IF (direction_norm <= 1.d-12) RETURN

    velocity = vpinch*direction/direction_norm
  END SUBROUTINE tm1d_build_pinch_velocity

END MODULE transport_models_1d_config
