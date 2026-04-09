SUBMODULE (transport_models_1d) transport_models_1d_diffusion_bohm_gyrobohm
IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE tm1d_compute_collisionality_profile(this, work)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    TYPE(transport_model_workspace_t), INTENT(INOUT) :: work
    REAL*8 :: ne_dim(this%nrho), te_dim(this%nrho), rmaj_dim(this%nrho), lambda_e(this%nrho)

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    ne_dim = MAX(ABS(this%ne_fs)*simpar%refval_density, model_tol)
    te_dim = MAX(ABS(this%te_fs)*simpar%refval_temperature, model_tol)
    rmaj_dim = MAX(ABS(this%Rmaj_fs)*simpar%refval_length, model_tol)

    lambda_e = 31.3d0 - LOG(SQRT(ne_dim)/te_dim)

    work%nuestar_fs = 6.921d-18 * ABS(this%q_fs) * rmaj_dim * ne_dim * MAX(phys%Zeff, 1.d0) * lambda_e / &
         (MAX(this%eps_fs, model_tol)**1.5d0 * te_dim**2)
  END SUBROUTINE tm1d_compute_collisionality_profile

  MODULE SUBROUTINE tm1d_compute_bohm_profile(this, work)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    TYPE(transport_model_workspace_t), INTENT(INOUT) :: work
    REAL*8 :: rho_s_te_fs(this%nrho)

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    rho_s_te_fs = work%cs_te_fs / MAX(this%omega_fs, model_tol)
    work%chi_bohm_fs = rho_s_te_fs * work%cs_te_fs * this%q_fs**2 * this%a_minor * &
         ABS(this%dpe_dr_fs) / MAX(this%pe_fs, model_tol) * this%delta_te
  END SUBROUTINE tm1d_compute_bohm_profile

  MODULE SUBROUTINE tm1d_compute_gyrobohm_profile(this, work)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    TYPE(transport_model_workspace_t), INTENT(INOUT) :: work
    REAL*8 :: rho_s_te_fs(this%nrho)

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    rho_s_te_fs = work%cs_te_fs / MAX(this%omega_fs, model_tol)
    work%chi_gyrobohm_fs = rho_s_te_fs**2 * work%cs_te_fs * ABS(this%dte_dr_fs) / MAX(this%te_fs, model_tol)
  END SUBROUTINE tm1d_compute_gyrobohm_profile

  MODULE SUBROUTINE tm1d_compute_mixed_transport(this, work)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(transport_model_workspace_t), INTENT(INOUT) :: work

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    this%chi_i_fs = MAX(this%c_bohm_i*work%chi_bohm_fs + this%c_gyrobohm_i*work%chi_gyrobohm_fs, 1.d-10)
    this%chi_e_fs = MAX(this%c_bohm_e*work%chi_bohm_fs + this%c_gyrobohm_e*work%chi_gyrobohm_fs, 1.d-10)
    this%d_part_fs = this%c_bohm_n * this%chi_i_fs*this%chi_e_fs / MAX(this%chi_i_fs + this%chi_e_fs, 1.d-10)
    this%nu_mom_fs = this%prandtl * this%chi_i_fs
  END SUBROUTINE tm1d_compute_mixed_transport

  MODULE SUBROUTINE tm1d_build_projected_gradients(this, fs_data)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    REAL*8 :: u1_safe(fs_data%nrho)

    u1_safe = MAX(fs_data%U_fs(1, :), model_tol)

    this%dpe_dr_fs = 2.d0/(3.d0*phys%Mref) * fs_data%Q_rad_fs(4, :)
    this%dte_dr_fs = 2.d0/(3.d0*phys%Mref) * &
         (fs_data%Q_rad_fs(4, :)/u1_safe - fs_data%U_fs(4, :)*fs_data%Q_rad_fs(1, :)/u1_safe**2)
  END SUBROUTINE tm1d_build_projected_gradients

END SUBMODULE transport_models_1d_diffusion_bohm_gyrobohm
