SUBMODULE (transport_models_1d) transport_models_1d_update
IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE tm1d_update_from_flux_surfaces(this, fs_data)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    REAL*8, ALLOCATABLE :: ua(:, :), up(:, :)
    REAL*8, ALLOCATABLE :: chi_bohm_fs(:), chi_gyrobohm_fs(:), nuestar_fs(:), cs_te_fs(:)

    IF (.NOT. fs_data%profiles_built) RETURN

    IF ((.NOT. this%is_initialized) .OR. this%nrho /= fs_data%nrho) THEN
       CALL this%init(fs_data%nrho, this%rho_edge, this%rho_core, this%rho_diffusion_model_max)
    END IF
    IF (.NOT. this%is_initialized) RETURN

    ALLOCATE(ua(fs_data%nrho, fs_data%neq))
    ALLOCATE(up(fs_data%nrho, phys%npv))
    ALLOCATE(chi_bohm_fs(this%nrho))
    ALLOCATE(chi_gyrobohm_fs(this%nrho))
    ALLOCATE(nuestar_fs(this%nrho))
    ALLOCATE(cs_te_fs(this%nrho))

    ua = TRANSPOSE(fs_data%U_fs)
    CALL cons2phys(ua, up)

    this%rho_grid = fs_data%rho_grid

    this%ne_fs = up(:, 1)
    this%pi_fs = up(:, 5)
    this%pe_fs = up(:, 6)
    this%ti_fs = up(:, 7)
    this%te_fs = up(:, 8)
    this%q_fs = fs_data%q_fs
    this%omega_fs = fs_data%omega_fs
    this%Rmaj_fs = fs_data%Rmaj_fs
    this%rmin_fs = fs_data%rmin_fs
    this%eps_fs = fs_data%eps_fs
    this%a_minor = phys%a_minor
    cs_te_fs = SQRT(MAX(this%te_fs*phys%Mref, model_tol))

    CALL tm1d_build_projected_gradients(this, fs_data)
    CALL this%compute_delta_te(fs_data)
    CALL tm1d_compute_collisionality_profile(this, nuestar_fs)
    CALL tm1d_compute_bohm_profile(this, cs_te_fs, chi_bohm_fs)
    CALL tm1d_compute_gyrobohm_profile(this, cs_te_fs, chi_gyrobohm_fs)
    CALL tm1d_compute_mixed_transport(this, chi_bohm_fs, chi_gyrobohm_fs)
    CALL this%compute_pinch_profile(nuestar_fs)

    DEALLOCATE(chi_bohm_fs, chi_gyrobohm_fs, nuestar_fs, cs_te_fs)
    DEALLOCATE(ua, up)
END SUBROUTINE tm1d_update_from_flux_surfaces

  MODULE SUBROUTINE tm1d_compute_delta_te(this, fs_data)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    REAL*8 :: te_core, te_edge

    this%delta_te = 0.d0
    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    CALL fs_data%interp_scalar(this%rho_core, this%te_fs, te_core)
    CALL fs_data%interp_scalar(this%rho_edge, this%te_fs, te_edge)
    te_edge = MAX(te_edge, model_tol)
    this%delta_te = (te_core - te_edge)/te_edge
  END SUBROUTINE tm1d_compute_delta_te

END SUBMODULE transport_models_1d_update
