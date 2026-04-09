MODULE transport_models_1d
  USE HDF5
  USE HDF5_io_module
  USE globals
  USE flux_surface_transport_data, ONLY: flux_surface_transport_t
  USE interpolation, ONLY: find_cell_and_local_coordinate
  USE physics, ONLY: cons2phys
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: transport_model_1d_t, transport_model_1d

  REAL*8, PARAMETER :: rho_edge_default = 0.99d0
  REAL*8, PARAMETER :: rho_core_default = 0.8d0
  REAL*8, PARAMETER :: model_tol = 1.d-12

  TYPE :: transport_model_1d_t
     LOGICAL :: is_initialized = .FALSE.
     INTEGER :: nrho = 0
     REAL*8, ALLOCATABLE :: rho_grid(:)
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
     REAL*8 :: a_minor = 0.d0
     REAL*8 :: delta_te = 0.d0
     REAL*8 :: c_bohm_i = 1.6d-4
     REAL*8 :: c_gyrobohm_i = 1.75d-2
     REAL*8 :: c_bohm_e = 8.d-5
     REAL*8 :: c_gyrobohm_e = 3.5d-2
     REAL*8 :: c_bohm_n = 1.d0
     REAL*8 :: prandtl = 1.d0
     REAL*8, ALLOCATABLE :: te_fs(:)
     REAL*8, ALLOCATABLE :: ti_fs(:)
     REAL*8, ALLOCATABLE :: ne_fs(:)
     REAL*8, ALLOCATABLE :: pe_fs(:)
     REAL*8, ALLOCATABLE :: pi_fs(:)
     REAL*8, ALLOCATABLE :: q_fs(:)
     REAL*8, ALLOCATABLE :: omega_fs(:)
     REAL*8, ALLOCATABLE :: Rmaj_fs(:)
     REAL*8, ALLOCATABLE :: rmin_fs(:)
     REAL*8, ALLOCATABLE :: eps_fs(:)
     REAL*8, ALLOCATABLE :: dte_dr_fs(:)
     REAL*8, ALLOCATABLE :: dpe_dr_fs(:)
     REAL*8, ALLOCATABLE :: chi_i_fs(:)
     REAL*8, ALLOCATABLE :: chi_e_fs(:)
     REAL*8, ALLOCATABLE :: d_part_fs(:)
     REAL*8, ALLOCATABLE :: nu_mom_fs(:)
     REAL*8, ALLOCATABLE :: vpinch_fs(:)
   CONTAINS
     PROCEDURE :: init => tm1d_init
     PROCEDURE :: destroy => tm1d_destroy
     PROCEDURE :: set_config => tm1d_set_config
     PROCEDURE :: update_from_flux_surfaces => tm1d_update_from_flux_surfaces
     PROCEDURE :: compute_delta_te => tm1d_compute_delta_te
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
     MODULE SUBROUTINE tm1d_compute_delta_te(this, fs_data)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
     END SUBROUTINE tm1d_compute_delta_te
     MODULE SUBROUTINE tm1d_compute_collisionality_profile(this, nuestar_fs)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(OUT) :: nuestar_fs(:)
     END SUBROUTINE tm1d_compute_collisionality_profile
     MODULE SUBROUTINE tm1d_compute_bohm_profile(this, cs_te_fs, chi_bohm_fs)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: cs_te_fs(:)
       REAL*8, INTENT(OUT) :: chi_bohm_fs(:)
     END SUBROUTINE tm1d_compute_bohm_profile
     MODULE SUBROUTINE tm1d_compute_gyrobohm_profile(this, cs_te_fs, chi_gyrobohm_fs)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: cs_te_fs(:)
       REAL*8, INTENT(OUT) :: chi_gyrobohm_fs(:)
     END SUBROUTINE tm1d_compute_gyrobohm_profile
     MODULE SUBROUTINE tm1d_compute_mixed_transport(this, chi_bohm_fs, chi_gyrobohm_fs)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: chi_bohm_fs(:), chi_gyrobohm_fs(:)
     END SUBROUTINE tm1d_compute_mixed_transport
     MODULE SUBROUTINE tm1d_build_projected_gradients(this, fs_data)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
     END SUBROUTINE tm1d_build_projected_gradients
     MODULE SUBROUTINE tm1d_compute_pinch_profile(this, nuestar_fs)
       CLASS(transport_model_1d_t), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: nuestar_fs(:)
     END SUBROUTINE tm1d_compute_pinch_profile
     MODULE SUBROUTINE tm1d_compute_militello_pinch(this, nuestar_fs, vpinch_fs)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: nuestar_fs(:)
       REAL*8, INTENT(OUT) :: vpinch_fs(:)
     END SUBROUTINE tm1d_compute_militello_pinch
     MODULE SUBROUTINE tm1d_compute_geometric_pinch(this, vpinch_fs)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(OUT) :: vpinch_fs(:)
     END SUBROUTINE tm1d_compute_geometric_pinch
     MODULE SUBROUTINE tm1d_compute_constant_pinch(this, vpinch_fs)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(OUT) :: vpinch_fs(:)
     END SUBROUTINE tm1d_compute_constant_pinch
     MODULE SUBROUTINE tm1d_compute_1D_pinch_matrix(this, b, rho, APinch)
       CLASS(transport_model_1d_t), INTENT(IN) :: this
       REAL*8, INTENT(IN) :: b(:), rho
       REAL*8, INTENT(OUT) :: APinch(:,:)
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
  END INTERFACE

CONTAINS

  SUBROUTINE tm1d_init(this, nrho, rho_edge, rho_core, rho_diffusion_model_max)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: nrho
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core, rho_diffusion_model_max

    CALL this%destroy()

    this%nrho = MAX(0, nrho)
    IF (PRESENT(rho_edge)) this%rho_edge = rho_edge
    IF (.NOT. PRESENT(rho_edge)) this%rho_edge = rho_edge_default
    IF (PRESENT(rho_core)) this%rho_core = rho_core
    IF (.NOT. PRESENT(rho_core)) this%rho_core = rho_core_default
    IF (PRESENT(rho_diffusion_model_max)) this%rho_diffusion_model_max = rho_diffusion_model_max
    IF (.NOT. PRESENT(rho_diffusion_model_max)) this%rho_diffusion_model_max = 1.d0

    IF (this%nrho <= 0) RETURN

    ALLOCATE(this%te_fs(this%nrho))
    ALLOCATE(this%ti_fs(this%nrho))
    ALLOCATE(this%ne_fs(this%nrho))
    ALLOCATE(this%pe_fs(this%nrho))
    ALLOCATE(this%pi_fs(this%nrho))
    ALLOCATE(this%q_fs(this%nrho))
    ALLOCATE(this%omega_fs(this%nrho))
    ALLOCATE(this%Rmaj_fs(this%nrho))
    ALLOCATE(this%rmin_fs(this%nrho))
    ALLOCATE(this%eps_fs(this%nrho))
    ALLOCATE(this%dte_dr_fs(this%nrho))
    ALLOCATE(this%dpe_dr_fs(this%nrho))
    ALLOCATE(this%chi_i_fs(this%nrho))
    ALLOCATE(this%chi_e_fs(this%nrho))
    ALLOCATE(this%d_part_fs(this%nrho))
    ALLOCATE(this%nu_mom_fs(this%nrho))
    ALLOCATE(this%vpinch_fs(this%nrho))
    ALLOCATE(this%rho_grid(this%nrho))

    this%te_fs = 0.d0
    this%ti_fs = 0.d0
    this%ne_fs = 0.d0
    this%pe_fs = 0.d0
    this%pi_fs = 0.d0
    this%q_fs = 0.d0
    this%omega_fs = 0.d0
    this%Rmaj_fs = 0.d0
    this%rmin_fs = 0.d0
    this%eps_fs = 0.d0
    this%dte_dr_fs = 0.d0
    this%dpe_dr_fs = 0.d0
    this%chi_i_fs = 0.d0
    this%chi_e_fs = 0.d0
    this%d_part_fs = 0.d0
    this%nu_mom_fs = 0.d0
    this%vpinch_fs = 0.d0
    this%rho_grid = 0.d0
    this%delta_te = 0.d0
    this%is_initialized = .TRUE.
  END SUBROUTINE tm1d_init

  SUBROUTINE tm1d_destroy(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%te_fs)) DEALLOCATE(this%te_fs)
    IF (ALLOCATED(this%ti_fs)) DEALLOCATE(this%ti_fs)
    IF (ALLOCATED(this%ne_fs)) DEALLOCATE(this%ne_fs)
    IF (ALLOCATED(this%pe_fs)) DEALLOCATE(this%pe_fs)
    IF (ALLOCATED(this%pi_fs)) DEALLOCATE(this%pi_fs)
    IF (ALLOCATED(this%q_fs)) DEALLOCATE(this%q_fs)
    IF (ALLOCATED(this%omega_fs)) DEALLOCATE(this%omega_fs)
    IF (ALLOCATED(this%Rmaj_fs)) DEALLOCATE(this%Rmaj_fs)
    IF (ALLOCATED(this%rmin_fs)) DEALLOCATE(this%rmin_fs)
    IF (ALLOCATED(this%eps_fs)) DEALLOCATE(this%eps_fs)
    IF (ALLOCATED(this%dte_dr_fs)) DEALLOCATE(this%dte_dr_fs)
    IF (ALLOCATED(this%dpe_dr_fs)) DEALLOCATE(this%dpe_dr_fs)
    IF (ALLOCATED(this%chi_i_fs)) DEALLOCATE(this%chi_i_fs)
    IF (ALLOCATED(this%chi_e_fs)) DEALLOCATE(this%chi_e_fs)
    IF (ALLOCATED(this%d_part_fs)) DEALLOCATE(this%d_part_fs)
    IF (ALLOCATED(this%nu_mom_fs)) DEALLOCATE(this%nu_mom_fs)
    IF (ALLOCATED(this%vpinch_fs)) DEALLOCATE(this%vpinch_fs)
    IF (ALLOCATED(this%rho_grid)) DEALLOCATE(this%rho_grid)

    this%is_initialized = .FALSE.
    this%nrho = 0
    this%rho_edge = rho_edge_default
    this%rho_core = rho_core_default
    this%rho_diffusion_model_max = 1.d0
    this%pinch_model = 1
    this%c_pinch = 0.5d0
    this%nu_th = 0.04d0
    this%vpinch_const = 0.d0
    this%rho_pinch_axis_width = 0.02d0
    this%rho_pinch_model_max = 0.99d0
    this%rho_pinch_edge_width = 0.03d0
    this%rho_blend_width = 0.02d0
    this%diff_n_min = 0.d0
    this%diff_u_min = 0.d0
    this%diff_e_min = 0.d0
    this%diff_ee_min = 0.d0
    this%a_minor = 0.d0
    this%delta_te = 0.d0
    this%c_bohm_i = 1.6d-4
    this%c_gyrobohm_i = 1.75d-2
    this%c_bohm_e = 8.d-5
    this%c_gyrobohm_e = 3.5d-2
    this%c_bohm_n = 1.d0
    this%prandtl = 1.d0
  END SUBROUTINE tm1d_destroy

  SUBROUTINE tm1d_finalize(this)
    TYPE(transport_model_1d_t), INTENT(INOUT) :: this

    CALL this%destroy()
  END SUBROUTINE tm1d_finalize

  SUBROUTINE tm1d_set_config(this, rho_edge, rho_core, rho_diffusion_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, prandtl, pinch_model, c_pinch, nu_th, vpinch_const_phys, rho_pinch_axis_width, rho_pinch_model_max, rho_pinch_edge_width, rho_blend_width, diff_n_min_phys, diff_u_min_phys, diff_e_min_phys, diff_ee_min_phys)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core, rho_diffusion_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, prandtl, c_pinch, nu_th, vpinch_const_phys, rho_pinch_axis_width, rho_pinch_model_max, rho_pinch_edge_width, rho_blend_width, diff_n_min_phys, diff_u_min_phys, diff_e_min_phys, diff_ee_min_phys
    INTEGER, INTENT(IN), OPTIONAL :: pinch_model

    IF (PRESENT(rho_edge)) this%rho_edge = rho_edge
    IF (PRESENT(rho_core)) this%rho_core = rho_core
    IF (PRESENT(rho_diffusion_model_max)) this%rho_diffusion_model_max = rho_diffusion_model_max
    IF (PRESENT(c_bohm_i)) this%c_bohm_i = c_bohm_i
    IF (PRESENT(c_gyrobohm_i)) this%c_gyrobohm_i = c_gyrobohm_i
    IF (PRESENT(c_bohm_e)) this%c_bohm_e = c_bohm_e
    IF (PRESENT(c_gyrobohm_e)) this%c_gyrobohm_e = c_gyrobohm_e
    IF (PRESENT(c_bohm_n)) this%c_bohm_n = c_bohm_n
    IF (PRESENT(prandtl)) this%prandtl = prandtl
    IF (PRESENT(pinch_model)) this%pinch_model = pinch_model
    IF (PRESENT(c_pinch)) this%c_pinch = c_pinch
    IF (PRESENT(nu_th)) this%nu_th = nu_th
    IF (PRESENT(vpinch_const_phys)) this%vpinch_const = vpinch_const_phys*simpar%refval_time/simpar%refval_length
    IF (PRESENT(rho_pinch_axis_width)) this%rho_pinch_axis_width = rho_pinch_axis_width
    IF (PRESENT(rho_pinch_model_max)) this%rho_pinch_model_max = rho_pinch_model_max
    IF (PRESENT(rho_pinch_edge_width)) this%rho_pinch_edge_width = rho_pinch_edge_width
    IF (PRESENT(rho_blend_width)) this%rho_blend_width = rho_blend_width
    IF (PRESENT(diff_n_min_phys)) this%diff_n_min = diff_n_min_phys*simpar%refval_time/simpar%refval_length**2
    IF (PRESENT(diff_u_min_phys)) this%diff_u_min = diff_u_min_phys*simpar%refval_time/simpar%refval_length**2
    IF (PRESENT(diff_e_min_phys)) this%diff_e_min = diff_e_min_phys*simpar%refval_time/simpar%refval_length**2
    IF (PRESENT(diff_ee_min_phys)) this%diff_ee_min = diff_ee_min_phys*simpar%refval_time/simpar%refval_length**2
  END SUBROUTINE tm1d_set_config


  SUBROUTINE tm1d_interp_transport(this, rho, chi_i, chi_e, d_part, nu_mom, vpinch)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(OUT) :: chi_i, chi_e, d_part, nu_mom, vpinch

    chi_i = 0.d0
    chi_e = 0.d0
    d_part = 0.d0
    nu_mom = 0.d0
    vpinch = 0.d0

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN
    IF (rho <= model_tol) RETURN

    CALL tm1d_interp_profile(this, rho, this%chi_i_fs, chi_i)
    CALL tm1d_interp_profile(this, rho, this%chi_e_fs, chi_e)
    CALL tm1d_interp_profile(this, rho, this%d_part_fs, d_part)
    CALL tm1d_interp_profile(this, rho, this%nu_mom_fs, nu_mom)
    CALL tm1d_interp_profile(this, rho, this%vpinch_fs, vpinch)
  END SUBROUTINE tm1d_interp_transport


  SUBROUTINE tm1d_apply_1D_diffusion(this, rho_pol_norm, diff_iso, diff_ani)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho_pol_norm(:)
    REAL*8, INTENT(INOUT) :: diff_iso(:, :, :), diff_ani(:, :, :)
    INTEGER :: g
    REAL*8 :: rho_g, w, chi_i, chi_e, d_part, nu_mom, vpinch

    IF (.NOT. this%is_initialized) RETURN
    IF (SIZE(diff_iso, 3) /= SIZE(rho_pol_norm)) RETURN
    IF (SIZE(diff_ani, 3) /= SIZE(rho_pol_norm)) RETURN

    DO g = 1, SIZE(rho_pol_norm)
       rho_g = MAX(rho_pol_norm(g), 0.d0)
       w = tm1d_blend_weight(this, rho_g)
       IF (w <= 0.d0) CYCLE

       CALL this%interp_transport(rho_g, chi_i, chi_e, d_part, nu_mom, vpinch)

       d_part = MAX(d_part, this%diff_n_min)
       nu_mom = MAX(nu_mom, this%diff_u_min)
       chi_i = MAX(chi_i, this%diff_e_min)
       chi_e = MAX(chi_e, this%diff_ee_min)

       diff_iso(1, 1, g) = diff_iso(1, 1, g) + w*(d_part - diff_iso(1, 1, g))
       diff_iso(2, 2, g) = diff_iso(2, 2, g) + w*(nu_mom - diff_iso(2, 2, g))
       diff_iso(3, 3, g) = diff_iso(3, 3, g) + w*(chi_i - diff_iso(3, 3, g))
       diff_iso(4, 4, g) = diff_iso(4, 4, g) + w*(chi_e - diff_iso(4, 4, g))

       diff_ani(1, 1, g) = diff_ani(1, 1, g) + w*(d_part - diff_ani(1, 1, g))
       diff_ani(2, 2, g) = diff_ani(2, 2, g) + w*(nu_mom - diff_ani(2, 2, g))
       diff_ani(3, 3, g) = diff_ani(3, 3, g) + w*(chi_i - diff_ani(3, 3, g))
       diff_ani(4, 4, g) = diff_ani(4, 4, g) + w*(chi_e - diff_ani(4, 4, g))
    END DO
  END SUBROUTINE tm1d_apply_1D_diffusion

  SUBROUTINE tm1d_write_hdf5(this, parent_group_id)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: coeffs_group_id, params_group_id
    INTEGER :: ierr

    IF (.NOT. switch%transport_1d) RETURN
    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    CALL HDF5_group_create('coefficients', parent_group_id, coeffs_group_id, ierr)
    CALL HDF5_array1D_saving(coeffs_group_id, this%chi_i_fs, SIZE(this%chi_i_fs), 'chi_i_fs')
    CALL HDF5_array1D_saving(coeffs_group_id, this%chi_e_fs, SIZE(this%chi_e_fs), 'chi_e_fs')
    CALL HDF5_array1D_saving(coeffs_group_id, this%d_part_fs, SIZE(this%d_part_fs), 'd_part_fs')
    CALL HDF5_array1D_saving(coeffs_group_id, this%nu_mom_fs, SIZE(this%nu_mom_fs), 'nu_mom_fs')
    CALL HDF5_array1D_saving(coeffs_group_id, this%vpinch_fs, SIZE(this%vpinch_fs), 'vpinch_fs')
    CALL HDF5_group_close(coeffs_group_id, ierr)

    CALL HDF5_group_create('params', parent_group_id, params_group_id, ierr)
    CALL HDF5_real_saving(params_group_id, this%a_minor, 'a_minor')
    CALL HDF5_real_saving(params_group_id, this%rho_core, 'rho_core')
    CALL HDF5_real_saving(params_group_id, this%rho_edge, 'rho_edge')
    CALL HDF5_real_saving(params_group_id, this%rho_diffusion_model_max, 'rho_diffusion_model_max')
    CALL HDF5_real_saving(params_group_id, this%rho_blend_width, 'rho_blend_width')
    CALL HDF5_real_saving(params_group_id, this%diff_n_min, 'diff_n_min')
    CALL HDF5_real_saving(params_group_id, this%diff_u_min, 'diff_u_min')
    CALL HDF5_real_saving(params_group_id, this%diff_e_min, 'diff_e_min')
    CALL HDF5_real_saving(params_group_id, this%diff_ee_min, 'diff_ee_min')
    CALL HDF5_real_saving(params_group_id, this%c_bohm_i, 'c_bohm_i')
    CALL HDF5_real_saving(params_group_id, this%c_gyrobohm_i, 'c_gyrobohm_i')
    CALL HDF5_real_saving(params_group_id, this%c_bohm_e, 'c_bohm_e')
    CALL HDF5_real_saving(params_group_id, this%c_gyrobohm_e, 'c_gyrobohm_e')
    CALL HDF5_real_saving(params_group_id, this%c_bohm_n, 'c_bohm_n')
    CALL HDF5_real_saving(params_group_id, this%prandtl, 'prandtl')
    CALL HDF5_integer_saving(params_group_id, this%pinch_model, 'pinch_model')
    CALL HDF5_real_saving(params_group_id, this%c_pinch, 'c_pinch')
    CALL HDF5_real_saving(params_group_id, this%nu_th, 'nu_th')
    CALL HDF5_real_saving(params_group_id, this%vpinch_const, 'vpinch_const')
    CALL HDF5_real_saving(params_group_id, this%rho_pinch_axis_width, 'rho_pinch_axis_width')
    CALL HDF5_real_saving(params_group_id, this%rho_pinch_model_max, 'rho_pinch_model_max')
    CALL HDF5_real_saving(params_group_id, this%rho_pinch_edge_width, 'rho_pinch_edge_width')
    CALL HDF5_group_close(params_group_id, ierr)
  END SUBROUTINE tm1d_write_hdf5

END MODULE transport_models_1d
