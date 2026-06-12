MODULE transport_models_1d_config
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: transport_model_config_t, tm1d_config_reset, tm1d_config_apply

  REAL*8, PARAMETER :: rho_edge_default = 0.99d0
  REAL*8, PARAMETER :: rho_core_default = 0.8d0

  TYPE :: transport_model_config_t
     REAL*8 :: rho_edge = rho_edge_default
     REAL*8 :: rho_core = rho_core_default
     REAL*8 :: rho_diffusion_model_max = 1.d0
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
     REAL*8 :: c_bohm_n_rho_slope = 0.7d0
     REAL*8 :: prandtl = 1.d0
  END TYPE transport_model_config_t

CONTAINS

  SUBROUTINE tm1d_config_reset(config)
    TYPE(transport_model_config_t), INTENT(INOUT) :: config

    config%rho_edge = rho_edge_default
    config%rho_core = rho_core_default
    config%rho_diffusion_model_max = 1.d0
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
    config%c_bohm_n_rho_slope = 0.7d0
    config%prandtl = 1.d0
  END SUBROUTINE tm1d_config_reset

  SUBROUTINE tm1d_config_apply(config, refval_time, refval_length, rho_edge, rho_core, rho_diffusion_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, c_bohm_n_rho_slope, prandtl, pinch_model, c_pinch, nu_th, vpinch_const_phys, rho_pinch_axis_width, rho_pinch_model_max, rho_pinch_edge_width, rho_blend_width, diff_n_min_phys, diff_u_min_phys, diff_e_min_phys, diff_ee_min_phys)
    TYPE(transport_model_config_t), INTENT(INOUT) :: config
    REAL*8, INTENT(IN) :: refval_time, refval_length
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core, rho_diffusion_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, c_bohm_n_rho_slope, prandtl, c_pinch, nu_th, vpinch_const_phys, rho_pinch_axis_width, rho_pinch_model_max, rho_pinch_edge_width, rho_blend_width, diff_n_min_phys, diff_u_min_phys, diff_e_min_phys, diff_ee_min_phys
    INTEGER, INTENT(IN), OPTIONAL :: pinch_model

    IF (PRESENT(rho_edge)) config%rho_edge = rho_edge
    IF (PRESENT(rho_core)) config%rho_core = rho_core
    IF (PRESENT(rho_diffusion_model_max)) config%rho_diffusion_model_max = rho_diffusion_model_max
    IF (PRESENT(c_bohm_i)) config%c_bohm_i = c_bohm_i
    IF (PRESENT(c_gyrobohm_i)) config%c_gyrobohm_i = c_gyrobohm_i
    IF (PRESENT(c_bohm_e)) config%c_bohm_e = c_bohm_e
    IF (PRESENT(c_gyrobohm_e)) config%c_gyrobohm_e = c_gyrobohm_e
    IF (PRESENT(c_bohm_n)) config%c_bohm_n = c_bohm_n
    IF (PRESENT(c_bohm_n_rho_slope)) config%c_bohm_n_rho_slope = c_bohm_n_rho_slope
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

END MODULE transport_models_1d_config
