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
     REAL*8 :: rho_edge = rho_edge_default
     REAL*8 :: rho_core = rho_core_default
     REAL*8 :: rho_model_max = 1.d0
     INTEGER :: pinch_model = 1
     REAL*8 :: c_pinch = 0.5d0
     REAL*8 :: nu_th = 0.04d0
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
     REAL*8, ALLOCATABLE :: cs_te_fs(:)
     REAL*8, ALLOCATABLE :: dte_dr_fs(:)
     REAL*8, ALLOCATABLE :: dpe_dr_fs(:)
     REAL*8, ALLOCATABLE :: nuestar_fs(:)
     REAL*8, ALLOCATABLE :: chi_bohm_fs(:)
     REAL*8, ALLOCATABLE :: chi_gyrobohm_fs(:)
     REAL*8, ALLOCATABLE :: chi_i_fs(:)
     REAL*8, ALLOCATABLE :: chi_e_fs(:)
     REAL*8, ALLOCATABLE :: d_part_fs(:)
     REAL*8, ALLOCATABLE :: nu_mom_fs(:)
     REAL*8, ALLOCATABLE :: pinch_factor_militello_fs(:)
     REAL*8, ALLOCATABLE :: vpinch_militello_fs(:)
     REAL*8, ALLOCATABLE :: vpinch_geometric_fs(:)
     REAL*8, ALLOCATABLE :: vpinch_fs(:)
   CONTAINS
     PROCEDURE :: init => tm1d_init
     PROCEDURE :: destroy => tm1d_destroy
     PROCEDURE :: set_config => tm1d_set_config
     PROCEDURE :: update_from_flux_surfaces => tm1d_update_from_flux_surfaces
     PROCEDURE :: compute_delta_te => tm1d_compute_delta_te
     PROCEDURE :: compute_collisionality_profile => tm1d_compute_collisionality_profile
     PROCEDURE :: compute_bohm_profile => tm1d_compute_bohm_profile
     PROCEDURE :: compute_gyrobohm_profile => tm1d_compute_gyrobohm_profile
     PROCEDURE :: compute_mixed_transport => tm1d_compute_mixed_transport
     PROCEDURE :: compute_pinch_profile => tm1d_compute_pinch_profile
     PROCEDURE :: interp_transport => tm1d_interp_transport
     PROCEDURE :: write_hdf5 => tm1d_write_hdf5
     FINAL :: tm1d_finalize
  END TYPE transport_model_1d_t

  TYPE(transport_model_1d_t), SAVE :: transport_model_1d

CONTAINS

  SUBROUTINE tm1d_init(this, nrho, rho_edge, rho_core, rho_model_max)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: nrho
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core, rho_model_max

    CALL this%destroy()

    this%nrho = MAX(0, nrho)
    IF (PRESENT(rho_edge)) this%rho_edge = rho_edge
    IF (.NOT. PRESENT(rho_edge)) this%rho_edge = rho_edge_default
    IF (PRESENT(rho_core)) this%rho_core = rho_core
    IF (.NOT. PRESENT(rho_core)) this%rho_core = rho_core_default
    IF (PRESENT(rho_model_max)) this%rho_model_max = rho_model_max
    IF (.NOT. PRESENT(rho_model_max)) this%rho_model_max = 1.d0

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
    ALLOCATE(this%cs_te_fs(this%nrho))
    ALLOCATE(this%dte_dr_fs(this%nrho))
    ALLOCATE(this%dpe_dr_fs(this%nrho))
    ALLOCATE(this%nuestar_fs(this%nrho))
    ALLOCATE(this%chi_bohm_fs(this%nrho))
    ALLOCATE(this%chi_gyrobohm_fs(this%nrho))
    ALLOCATE(this%chi_i_fs(this%nrho))
    ALLOCATE(this%chi_e_fs(this%nrho))
    ALLOCATE(this%d_part_fs(this%nrho))
    ALLOCATE(this%nu_mom_fs(this%nrho))
    ALLOCATE(this%pinch_factor_militello_fs(this%nrho))
    ALLOCATE(this%vpinch_militello_fs(this%nrho))
    ALLOCATE(this%vpinch_geometric_fs(this%nrho))
    ALLOCATE(this%vpinch_fs(this%nrho))

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
    this%cs_te_fs = 0.d0
    this%dte_dr_fs = 0.d0
    this%dpe_dr_fs = 0.d0
    this%nuestar_fs = 0.d0
    this%chi_bohm_fs = 0.d0
    this%chi_gyrobohm_fs = 0.d0
    this%chi_i_fs = 0.d0
    this%chi_e_fs = 0.d0
    this%d_part_fs = 0.d0
    this%nu_mom_fs = 0.d0
    this%pinch_factor_militello_fs = 0.d0
    this%vpinch_militello_fs = 0.d0
    this%vpinch_geometric_fs = 0.d0
    this%vpinch_fs = 0.d0
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
    IF (ALLOCATED(this%cs_te_fs)) DEALLOCATE(this%cs_te_fs)
    IF (ALLOCATED(this%dte_dr_fs)) DEALLOCATE(this%dte_dr_fs)
    IF (ALLOCATED(this%dpe_dr_fs)) DEALLOCATE(this%dpe_dr_fs)
    IF (ALLOCATED(this%nuestar_fs)) DEALLOCATE(this%nuestar_fs)
    IF (ALLOCATED(this%chi_bohm_fs)) DEALLOCATE(this%chi_bohm_fs)
    IF (ALLOCATED(this%chi_gyrobohm_fs)) DEALLOCATE(this%chi_gyrobohm_fs)
    IF (ALLOCATED(this%chi_i_fs)) DEALLOCATE(this%chi_i_fs)
    IF (ALLOCATED(this%chi_e_fs)) DEALLOCATE(this%chi_e_fs)
    IF (ALLOCATED(this%d_part_fs)) DEALLOCATE(this%d_part_fs)
    IF (ALLOCATED(this%nu_mom_fs)) DEALLOCATE(this%nu_mom_fs)
    IF (ALLOCATED(this%pinch_factor_militello_fs)) DEALLOCATE(this%pinch_factor_militello_fs)
    IF (ALLOCATED(this%vpinch_militello_fs)) DEALLOCATE(this%vpinch_militello_fs)
    IF (ALLOCATED(this%vpinch_geometric_fs)) DEALLOCATE(this%vpinch_geometric_fs)
    IF (ALLOCATED(this%vpinch_fs)) DEALLOCATE(this%vpinch_fs)

    this%is_initialized = .FALSE.
    this%nrho = 0
    this%rho_edge = rho_edge_default
    this%rho_core = rho_core_default
    this%rho_model_max = 1.d0
    this%pinch_model = 1
    this%c_pinch = 0.5d0
    this%nu_th = 0.04d0
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

  SUBROUTINE tm1d_set_config(this, rho_edge, rho_core, rho_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, prandtl, pinch_model, c_pinch, nu_th)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core, rho_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, prandtl, c_pinch, nu_th
    INTEGER, INTENT(IN), OPTIONAL :: pinch_model

    IF (PRESENT(rho_edge)) this%rho_edge = rho_edge
    IF (PRESENT(rho_core)) this%rho_core = rho_core
    IF (PRESENT(rho_model_max)) this%rho_model_max = rho_model_max
    IF (PRESENT(c_bohm_i)) this%c_bohm_i = c_bohm_i
    IF (PRESENT(c_gyrobohm_i)) this%c_gyrobohm_i = c_gyrobohm_i
    IF (PRESENT(c_bohm_e)) this%c_bohm_e = c_bohm_e
    IF (PRESENT(c_gyrobohm_e)) this%c_gyrobohm_e = c_gyrobohm_e
    IF (PRESENT(c_bohm_n)) this%c_bohm_n = c_bohm_n
    IF (PRESENT(prandtl)) this%prandtl = prandtl
    IF (PRESENT(pinch_model)) this%pinch_model = pinch_model
    IF (PRESENT(c_pinch)) this%c_pinch = c_pinch
    IF (PRESENT(nu_th)) this%nu_th = nu_th
  END SUBROUTINE tm1d_set_config

  SUBROUTINE tm1d_update_from_flux_surfaces(this, fs_data)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    REAL*8, ALLOCATABLE :: ua(:, :), up(:, :)

    IF (.NOT. fs_data%profiles_built) RETURN

    IF ((.NOT. this%is_initialized) .OR. this%nrho /= fs_data%nrho) THEN
       CALL this%init(fs_data%nrho, this%rho_edge, this%rho_core, this%rho_model_max)
    END IF
    IF (.NOT. this%is_initialized) RETURN

    ALLOCATE(ua(fs_data%nrho, fs_data%neq))
    ALLOCATE(up(fs_data%nrho, phys%npv))

    ua = TRANSPOSE(fs_data%U_fs)
    CALL cons2phys(ua, up)

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
    this%cs_te_fs = SQRT(MAX(this%te_fs*phys%Mref, model_tol))

    CALL tm1d_build_projected_gradients(this, fs_data)
    CALL this%compute_delta_te(fs_data)
    CALL this%compute_collisionality_profile()
    CALL this%compute_bohm_profile()
    CALL this%compute_gyrobohm_profile()
    CALL this%compute_mixed_transport()
    CALL this%compute_pinch_profile()

    DEALLOCATE(ua, up)
  END SUBROUTINE tm1d_update_from_flux_surfaces

  SUBROUTINE tm1d_compute_delta_te(this, fs_data)
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

  SUBROUTINE tm1d_compute_collisionality_profile(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8 :: ne_dim(this%nrho), te_dim(this%nrho), rmaj_dim(this%nrho), lambda_e(this%nrho)

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    ne_dim = MAX(ABS(this%ne_fs)*simpar%refval_density, model_tol)
    te_dim = MAX(ABS(this%te_fs)*simpar%refval_temperature, model_tol)
    rmaj_dim = MAX(ABS(this%Rmaj_fs)*simpar%refval_length, model_tol)

    lambda_e = 31.3d0 - LOG(SQRT(ne_dim)/te_dim)

    this%nuestar_fs = 6.921d-18 * ABS(this%q_fs) * rmaj_dim * ne_dim * MAX(phys%Zeff, 1.d0) * lambda_e / &
         (MAX(this%eps_fs, model_tol)**1.5d0 * te_dim**2)
  END SUBROUTINE tm1d_compute_collisionality_profile

  SUBROUTINE tm1d_compute_bohm_profile(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8 :: rho_s_te_fs(this%nrho)

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    rho_s_te_fs = this%cs_te_fs / MAX(this%omega_fs, model_tol)
    this%chi_bohm_fs = rho_s_te_fs * this%cs_te_fs * this%q_fs**2 * this%a_minor * &
         ABS(this%dpe_dr_fs) / MAX(this%pe_fs, model_tol) * this%delta_te
  END SUBROUTINE tm1d_compute_bohm_profile

  SUBROUTINE tm1d_compute_gyrobohm_profile(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8 :: rho_s_te_fs(this%nrho)

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    rho_s_te_fs = this%cs_te_fs / MAX(this%omega_fs, model_tol)
    this%chi_gyrobohm_fs = rho_s_te_fs**2 * this%cs_te_fs * ABS(this%dte_dr_fs) / MAX(this%te_fs, model_tol)
  END SUBROUTINE tm1d_compute_gyrobohm_profile

  SUBROUTINE tm1d_compute_mixed_transport(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    this%chi_i_fs = MAX(this%c_bohm_i*this%chi_bohm_fs + this%c_gyrobohm_i*this%chi_gyrobohm_fs, 1.d-10)
    this%chi_e_fs = MAX(this%c_bohm_e*this%chi_bohm_fs + this%c_gyrobohm_e*this%chi_gyrobohm_fs, 1.d-10)
    this%d_part_fs = this%c_bohm_n * this%chi_i_fs*this%chi_e_fs / MAX(this%chi_i_fs + this%chi_e_fs, 1.d-10)
    this%nu_mom_fs = this%prandtl * this%chi_i_fs
  END SUBROUTINE tm1d_compute_mixed_transport

  SUBROUTINE tm1d_compute_pinch_profile(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    this%pinch_factor_militello_fs = MIN(1.d0, EXP(1.d0 - this%nuestar_fs/MAX(this%nu_th, model_tol)))
    this%vpinch_militello_fs = this%pinch_factor_militello_fs * this%c_pinch * this%d_part_fs * this%rmin_fs / MAX(this%a_minor, model_tol)**2
    this%vpinch_geometric_fs = this%c_pinch * this%d_part_fs * this%rmin_fs / MAX(this%a_minor, model_tol)**2

    SELECT CASE (this%pinch_model)
    CASE (1)
       this%vpinch_fs = this%vpinch_militello_fs
    CASE (2)
       this%vpinch_fs = this%vpinch_geometric_fs
    CASE DEFAULT
       this%vpinch_fs = 0.d0
    END SELECT
  END SUBROUTINE tm1d_compute_pinch_profile


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
    IF (rho > this%rho_model_max) RETURN

    CALL tm1d_interp_profile(this, rho, this%chi_i_fs, chi_i)
    CALL tm1d_interp_profile(this, rho, this%chi_e_fs, chi_e)
    CALL tm1d_interp_profile(this, rho, this%d_part_fs, d_part)
    CALL tm1d_interp_profile(this, rho, this%nu_mom_fs, nu_mom)
    CALL tm1d_interp_profile(this, rho, this%vpinch_fs, vpinch)
  END SUBROUTINE tm1d_interp_transport

  SUBROUTINE tm1d_write_hdf5(this, parent_group_id)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    IF (.NOT. switch%transport_1d) RETURN
    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    CALL HDF5_group_create('transport_model_1d', parent_group_id, group_id, ierr)
    CALL HDF5_array1D_saving(group_id, this%ne_fs, SIZE(this%ne_fs), 'ne_fs')
    CALL HDF5_array1D_saving(group_id, this%ti_fs, SIZE(this%ti_fs), 'ti_fs')
    CALL HDF5_array1D_saving(group_id, this%te_fs, SIZE(this%te_fs), 'te_fs')
    CALL HDF5_array1D_saving(group_id, this%pi_fs, SIZE(this%pi_fs), 'pi_fs')
    CALL HDF5_array1D_saving(group_id, this%pe_fs, SIZE(this%pe_fs), 'pe_fs')
    CALL HDF5_array1D_saving(group_id, this%q_fs, SIZE(this%q_fs), 'q_fs')
    CALL HDF5_array1D_saving(group_id, this%omega_fs, SIZE(this%omega_fs), 'omega_fs')
    CALL HDF5_array1D_saving(group_id, this%Rmaj_fs, SIZE(this%Rmaj_fs), 'Rmaj_fs')
    CALL HDF5_array1D_saving(group_id, this%rmin_fs, SIZE(this%rmin_fs), 'rmin_fs')
    CALL HDF5_array1D_saving(group_id, this%eps_fs, SIZE(this%eps_fs), 'eps_fs')
    CALL HDF5_array1D_saving(group_id, this%cs_te_fs, SIZE(this%cs_te_fs), 'cs_te_fs')
    CALL HDF5_array1D_saving(group_id, this%dte_dr_fs, SIZE(this%dte_dr_fs), 'dte_dr_fs')
    CALL HDF5_array1D_saving(group_id, this%dpe_dr_fs, SIZE(this%dpe_dr_fs), 'dpe_dr_fs')
    CALL HDF5_array1D_saving(group_id, this%nuestar_fs, SIZE(this%nuestar_fs), 'nuestar_fs')
    CALL HDF5_array1D_saving(group_id, this%chi_bohm_fs, SIZE(this%chi_bohm_fs), 'chi_bohm_fs')
    CALL HDF5_array1D_saving(group_id, this%chi_gyrobohm_fs, SIZE(this%chi_gyrobohm_fs), 'chi_gyrobohm_fs')
    CALL HDF5_array1D_saving(group_id, this%chi_i_fs, SIZE(this%chi_i_fs), 'chi_i_fs')
    CALL HDF5_array1D_saving(group_id, this%chi_e_fs, SIZE(this%chi_e_fs), 'chi_e_fs')
    CALL HDF5_array1D_saving(group_id, this%d_part_fs, SIZE(this%d_part_fs), 'd_part_fs')
    CALL HDF5_array1D_saving(group_id, this%nu_mom_fs, SIZE(this%nu_mom_fs), 'nu_mom_fs')
    CALL HDF5_array1D_saving(group_id, this%pinch_factor_militello_fs, SIZE(this%pinch_factor_militello_fs), 'pinch_factor_militello_fs')
    CALL HDF5_array1D_saving(group_id, this%vpinch_militello_fs, SIZE(this%vpinch_militello_fs), 'vpinch_militello_fs')
    CALL HDF5_array1D_saving(group_id, this%vpinch_geometric_fs, SIZE(this%vpinch_geometric_fs), 'vpinch_geometric_fs')
    CALL HDF5_array1D_saving(group_id, this%vpinch_fs, SIZE(this%vpinch_fs), 'vpinch_fs')
    CALL HDF5_real_saving(group_id, this%a_minor, 'a_minor')
    CALL HDF5_real_saving(group_id, this%rho_core, 'rho_core')
    CALL HDF5_real_saving(group_id, this%rho_edge, 'rho_edge')
    CALL HDF5_real_saving(group_id, this%rho_model_max, 'rho_model_max')
    CALL HDF5_real_saving(group_id, this%delta_te, 'delta_te')
    CALL HDF5_real_saving(group_id, this%c_bohm_i, 'c_bohm_i')
    CALL HDF5_real_saving(group_id, this%c_gyrobohm_i, 'c_gyrobohm_i')
    CALL HDF5_real_saving(group_id, this%c_bohm_e, 'c_bohm_e')
    CALL HDF5_real_saving(group_id, this%c_gyrobohm_e, 'c_gyrobohm_e')
    CALL HDF5_real_saving(group_id, this%c_bohm_n, 'c_bohm_n')
    CALL HDF5_real_saving(group_id, this%prandtl, 'prandtl')
    CALL HDF5_integer_saving(group_id, this%pinch_model, 'pinch_model')
    CALL HDF5_real_saving(group_id, this%c_pinch, 'c_pinch')
    CALL HDF5_real_saving(group_id, this%nu_th, 'nu_th')
    CALL HDF5_group_close(group_id, ierr)
  END SUBROUTINE tm1d_write_hdf5

  SUBROUTINE tm1d_build_projected_gradients(this, fs_data)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    REAL*8 :: u1_safe(fs_data%nrho)

    u1_safe = MAX(fs_data%U_fs(1, :), model_tol)

    this%dpe_dr_fs = 2.d0/(3.d0*phys%Mref) * fs_data%Q_rad_fs(4, :)
    this%dte_dr_fs = 2.d0/(3.d0*phys%Mref) * &
         (fs_data%Q_rad_fs(4, :)/u1_safe - fs_data%U_fs(4, :)*fs_data%Q_rad_fs(1, :)/u1_safe**2)
  END SUBROUTINE tm1d_build_projected_gradients


  SUBROUTINE tm1d_interp_profile(this, rho, profile, value)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(IN) :: profile(:)
    REAL*8, INTENT(OUT) :: value
    INTEGER :: ilow, ihigh, i
    REAL*8 :: alpha
    REAL*8 :: rho_grid(this%nrho)

    value = 0.d0

    IF (this%nrho <= 0) RETURN
    IF (SIZE(profile) /= this%nrho) RETURN
    IF (this%nrho == 1) THEN
       value = profile(1)
       RETURN
    END IF

    DO i = 1, this%nrho
       rho_grid(i) = DBLE(i - 1)/DBLE(MAX(this%nrho - 1, 1))
    END DO

    CALL find_cell_and_local_coordinate(this%nrho, rho_grid, MIN(MAX(rho, 0.d0), 1.d0), ilow, alpha)
    ilow = MIN(MAX(ilow, 1), this%nrho)
    ihigh = MIN(ilow + 1, this%nrho)

    IF (ihigh == ilow) THEN
       value = profile(ilow)
    ELSE
       value = (1.d0 - alpha)*profile(ilow) + alpha*profile(ihigh)
    END IF
  END SUBROUTINE tm1d_interp_profile

END MODULE transport_models_1d
