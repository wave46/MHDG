SUBMODULE (transport_models_1d) transport_models_1d_update
IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE tm1d_update_from_flux_surfaces(this, fs_data)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    TYPE(transport_model_derived_t) :: work

    IF (.NOT. fs_data%profiles_built) RETURN

    IF ((.NOT. this%is_initialized) .OR. this%nrho /= fs_data%nrho) THEN
       CALL this%init(fs_data%nrho, this%config%rho_edge, this%config%rho_core, this%config%rho_diffusion_model_max)
    END IF
    IF (.NOT. this%is_initialized) RETURN

    CALL work%init(this%nrho)

    this%rho_grid = fs_data%rho_grid
    this%a_minor = phys%a_minor
    CALL tm1d_fill_workspace_from_fs(fs_data, work)

    CALL tm1d_build_projected_gradients(this, fs_data, work)
    CALL this%compute_delta_te(fs_data, work)
    CALL tm1d_compute_collisionality_profile(this, work)
    CALL tm1d_compute_bohm_profile(this, work)
    CALL tm1d_compute_gyrobohm_profile(this, work)
    CALL tm1d_compute_mixed_transport(this, work)
    CALL this%compute_pinch_profile(work)

    CALL work%destroy()
  END SUBROUTINE tm1d_update_from_flux_surfaces

  MODULE SUBROUTINE tm1d_compute_delta_te(this, fs_data, work)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    TYPE(transport_model_derived_t), INTENT(IN) :: work
    REAL*8 :: te_core, te_edge

    this%delta_te = 0.d0
    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    CALL fs_data%interp_scalar(this%config%rho_core, work%te_fs, te_core)
    CALL fs_data%interp_scalar(this%config%rho_edge, work%te_fs, te_edge)
    te_edge = MAX(te_edge, model_tol)
    this%delta_te = (te_core - te_edge)/te_edge
  END SUBROUTINE tm1d_compute_delta_te

  SUBROUTINE tm1d_fill_workspace_from_fs(fs_data, work)
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    TYPE(transport_model_derived_t), INTENT(INOUT) :: work
    REAL*8, ALLOCATABLE :: ua(:, :), up(:, :)

    ALLOCATE(ua(fs_data%nrho, fs_data%neq))
    ALLOCATE(up(fs_data%nrho, phys%npv))

    ua = TRANSPOSE(fs_data%U_fs)
    CALL cons2phys(ua, up)

    work%ne_fs = up(:, 1)
    work%pi_fs = up(:, 5)
    work%pe_fs = up(:, 6)
    work%ti_fs = up(:, 7)
    work%te_fs = up(:, 8)
    work%q_fs = fs_data%q_fs
    work%omega_fs = fs_data%omega_fs
    work%Rmaj_fs = fs_data%Rmaj_fs
    work%rmin_fs = fs_data%rmin_fs
    work%eps_fs = fs_data%eps_fs

    DEALLOCATE(ua, up)
  END SUBROUTINE tm1d_fill_workspace_from_fs

END SUBMODULE transport_models_1d_update
