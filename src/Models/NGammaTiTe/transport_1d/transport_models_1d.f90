MODULE transport_models_1d
  USE HDF5
  USE HDF5_io_module
  USE globals
  USE flux_surface_transport_data, ONLY: flux_surface_transport_t
  USE interpolation, ONLY: find_cell_and_local_coordinate
  USE physics, ONLY: cons2phys
  USE transport_models_1d_config, ONLY: transport_model_config_t, &
       tm1d_config_reset, tm1d_config_apply, tm1d_region_policy_name, &
       tm1d_region_is_included, tm1d_build_pinch_velocity, &
       tm1d_particle_taper_factor
  USE transport_models_1d_derived, ONLY: transport_model_derived_t
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: transport_model_1d_t, transport_model_1d

  REAL*8, PARAMETER :: model_tol = 1.d-12

  TYPE :: transport_model_1d_t
     LOGICAL :: is_initialized = .FALSE.
     INTEGER :: nrho = 0
     REAL*8, ALLOCATABLE :: rho_grid(:)
     TYPE(transport_model_config_t) :: config
     REAL*8, ALLOCATABLE :: chi_i_fs(:)
     REAL*8, ALLOCATABLE :: chi_e_fs(:)
     REAL*8, ALLOCATABLE :: d_fs(:)
     REAL*8, ALLOCATABLE :: nu_mom_fs(:)
     REAL*8, ALLOCATABLE :: vpinch_fs(:)
   CONTAINS
     PROCEDURE :: init => tm1d_init
     PROCEDURE :: destroy => tm1d_destroy
     PROCEDURE :: set_config => tm1d_set_config
     PROCEDURE :: update_from_flux_surfaces => tm1d_update_from_flux_surfaces
     PROCEDURE :: compute_pinch_profile => tm1d_compute_pinch_profile
     PROCEDURE :: interp_transport => tm1d_interp_transport
     PROCEDURE :: compute_1D_pinch_matrix => tm1d_compute_1D_pinch_matrix
     PROCEDURE :: apply_1D_diffusion => tm1d_apply_1D_diffusion
     PROCEDURE :: write_hdf5 => tm1d_write_hdf5
     FINAL :: tm1d_finalize
  END TYPE transport_model_1d_t

  TYPE(transport_model_1d_t), SAVE :: transport_model_1d

  INTERFACE
     MODULE SUBROUTINE tm1d_update_from_flux_surfaces(this, fs_data)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
     END SUBROUTINE tm1d_update_from_flux_surfaces
     MODULE SUBROUTINE tm1d_compute_collisionality_profile(this, work)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       TYPE(transport_model_derived_t), INTENT(INOUT) :: work
     END SUBROUTINE tm1d_compute_collisionality_profile
     MODULE SUBROUTINE tm1d_compute_bohm_profile(this, work)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       TYPE(transport_model_derived_t), INTENT(INOUT) :: work
     END SUBROUTINE tm1d_compute_bohm_profile
     MODULE SUBROUTINE tm1d_compute_gyrobohm_profile(this, work)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       TYPE(transport_model_derived_t), INTENT(INOUT) :: work
     END SUBROUTINE tm1d_compute_gyrobohm_profile
     MODULE SUBROUTINE tm1d_compute_mixed_transport(this, work)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       TYPE(transport_model_derived_t), INTENT(INOUT) :: work
     END SUBROUTINE tm1d_compute_mixed_transport
     MODULE SUBROUTINE tm1d_build_projected_gradients(this, fs_data, work)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
       TYPE(transport_model_derived_t), INTENT(INOUT) :: work
     END SUBROUTINE tm1d_build_projected_gradients
     MODULE SUBROUTINE tm1d_compute_pinch_profile(this, work)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       TYPE(transport_model_derived_t), INTENT(IN) :: work
     END SUBROUTINE tm1d_compute_pinch_profile
     MODULE SUBROUTINE tm1d_compute_militello_pinch(this, work)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       TYPE(transport_model_derived_t), INTENT(IN) :: work
     END SUBROUTINE tm1d_compute_militello_pinch
     MODULE SUBROUTINE tm1d_compute_geometric_pinch(this, work)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       TYPE(transport_model_derived_t), INTENT(IN) :: work
     END SUBROUTINE tm1d_compute_geometric_pinch
     MODULE SUBROUTINE tm1d_compute_constant_pinch(this, work)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       TYPE(transport_model_derived_t), INTENT(IN) :: work
     END SUBROUTINE tm1d_compute_constant_pinch
     MODULE SUBROUTINE tm1d_compute_1D_pinch_matrix(this, b, rho, APinch, &
          region, outward_normal)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: b(:), rho
       REAL*8, INTENT(OUT) :: APinch(:,:)
       INTEGER, INTENT(IN) :: region
       REAL*8, INTENT(IN) :: outward_normal(:)
     END SUBROUTINE tm1d_compute_1D_pinch_matrix
     MODULE REAL*8 FUNCTION tm1d_pinch_window(this, rho)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: rho
     END FUNCTION tm1d_pinch_window
     MODULE REAL*8 FUNCTION tm1d_smoothstep01(s)
       REAL*8, INTENT(IN) :: s
     END FUNCTION tm1d_smoothstep01
     MODULE REAL*8 FUNCTION tm1d_axis_ramp(rho, width)
       REAL*8, INTENT(IN) :: rho, width
     END FUNCTION tm1d_axis_ramp
     MODULE REAL*8 FUNCTION tm1d_edge_cutoff(rho, rho_max, width)
       REAL*8, INTENT(IN) :: rho, rho_max, width
     END FUNCTION tm1d_edge_cutoff
     MODULE REAL*8 FUNCTION tm1d_blend_weight(this, rho)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: rho
     END FUNCTION tm1d_blend_weight
     MODULE SUBROUTINE tm1d_interp_profile(this, rho, profile, value)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: rho
       REAL*8, INTENT(IN) :: profile(:)
       REAL*8, INTENT(OUT) :: value
     END SUBROUTINE tm1d_interp_profile
     MODULE SUBROUTINE tm1d_write_hdf5(this, parent_group_id)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       INTEGER(HID_T), INTENT(IN) :: parent_group_id
     END SUBROUTINE tm1d_write_hdf5
     MODULE SUBROUTINE tm1d_interp_transport(this, rho, chi_i, chi_e, d, nu_mom, vpinch)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: rho
       REAL*8, INTENT(OUT) :: chi_i, chi_e, d, nu_mom, vpinch
     END SUBROUTINE tm1d_interp_transport
     MODULE SUBROUTINE tm1d_apply_1D_diffusion(this, rho_pol_norm, diff_iso, &
          diff_ani, region)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: rho_pol_norm(:)
       REAL*8, INTENT(INOUT) :: diff_iso(:, :, :), diff_ani(:, :, :)
       INTEGER, INTENT(IN) :: region(:)
     END SUBROUTINE tm1d_apply_1D_diffusion
  END INTERFACE

CONTAINS

  SUBROUTINE tm1d_init(this, nrho, rho_edge, rho_core, rho_diffusion_model_max)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: nrho
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core, rho_diffusion_model_max

    CALL tm1d_clear_storage(this)

    this%nrho = MAX(0, nrho)
    CALL tm1d_config_apply(this%config, simpar%refval_time, simpar%refval_length, &
         rho_edge=rho_edge, rho_core=rho_core, rho_diffusion_model_max=rho_diffusion_model_max)

    IF (this%nrho <= 0) RETURN

    ALLOCATE(this%chi_i_fs(this%nrho))
    ALLOCATE(this%chi_e_fs(this%nrho))
    ALLOCATE(this%d_fs(this%nrho))
    ALLOCATE(this%nu_mom_fs(this%nrho))
    ALLOCATE(this%vpinch_fs(this%nrho))
    ALLOCATE(this%rho_grid(this%nrho))

    this%chi_i_fs = 0.d0
    this%chi_e_fs = 0.d0
    this%d_fs = 0.d0
    this%nu_mom_fs = 0.d0
    this%vpinch_fs = 0.d0
    this%rho_grid = 0.d0
    this%is_initialized = .TRUE.
  END SUBROUTINE tm1d_init

  SUBROUTINE tm1d_destroy(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this

    CALL tm1d_clear_storage(this)
    CALL tm1d_config_reset(this%config)
  END SUBROUTINE tm1d_destroy

  SUBROUTINE tm1d_clear_storage(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%chi_i_fs)) DEALLOCATE(this%chi_i_fs)
    IF (ALLOCATED(this%chi_e_fs)) DEALLOCATE(this%chi_e_fs)
    IF (ALLOCATED(this%d_fs)) DEALLOCATE(this%d_fs)
    IF (ALLOCATED(this%nu_mom_fs)) DEALLOCATE(this%nu_mom_fs)
    IF (ALLOCATED(this%vpinch_fs)) DEALLOCATE(this%vpinch_fs)
    IF (ALLOCATED(this%rho_grid)) DEALLOCATE(this%rho_grid)

    this%is_initialized = .FALSE.
    this%nrho = 0
  END SUBROUTINE tm1d_clear_storage

  SUBROUTINE tm1d_set_config(this, rho_edge, rho_core, rho_diffusion_model_max, &
       transport_region_policy, c_bohm_i, c_gyrobohm_i, c_bohm_e, &
       c_gyrobohm_e, c_bohm_n, c_bohm_n_rho_slope, prandtl, pinch_model, &
       c_pinch, nu_th, &
       vpinch_const_phys, rho_pinch_axis_width, rho_pinch_model_max, &
       rho_pinch_edge_width, rho_blend_width, diff_n_min_phys, diff_u_min_phys, &
       diff_e_min_phys, diff_ee_min_phys)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core, rho_diffusion_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, c_bohm_n_rho_slope, prandtl, c_pinch, nu_th, vpinch_const_phys, rho_pinch_axis_width, rho_pinch_model_max, rho_pinch_edge_width, rho_blend_width, diff_n_min_phys, diff_u_min_phys, diff_e_min_phys, diff_ee_min_phys
    INTEGER, INTENT(IN), OPTIONAL :: pinch_model
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: transport_region_policy

    CALL tm1d_config_apply(this%config, simpar%refval_time, simpar%refval_length, &
         rho_edge=rho_edge, rho_core=rho_core, rho_diffusion_model_max=rho_diffusion_model_max, &
         transport_region_policy=transport_region_policy, &
         c_bohm_i=c_bohm_i, c_gyrobohm_i=c_gyrobohm_i, c_bohm_e=c_bohm_e, &
         c_gyrobohm_e=c_gyrobohm_e, c_bohm_n=c_bohm_n, &
         c_bohm_n_rho_slope=c_bohm_n_rho_slope, prandtl=prandtl, &
         pinch_model=pinch_model, c_pinch=c_pinch, nu_th=nu_th, &
         vpinch_const_phys=vpinch_const_phys, rho_pinch_axis_width=rho_pinch_axis_width, &
         rho_pinch_model_max=rho_pinch_model_max, rho_pinch_edge_width=rho_pinch_edge_width, &
         rho_blend_width=rho_blend_width, diff_n_min_phys=diff_n_min_phys, &
         diff_u_min_phys=diff_u_min_phys, diff_e_min_phys=diff_e_min_phys, &
         diff_ee_min_phys=diff_ee_min_phys)
  END SUBROUTINE tm1d_set_config

  SUBROUTINE tm1d_finalize(this)
    TYPE(transport_model_1d_t), INTENT(INOUT) :: this

    CALL this%destroy()
  END SUBROUTINE tm1d_finalize

END MODULE transport_models_1d
