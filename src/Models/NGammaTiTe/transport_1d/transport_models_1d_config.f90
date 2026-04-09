MODULE transport_models_1d_config
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: transport_model_config_t

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
     REAL*8 :: prandtl = 1.d0
  END TYPE transport_model_config_t

END MODULE transport_models_1d_config
