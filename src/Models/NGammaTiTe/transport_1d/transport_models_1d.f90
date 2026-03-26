MODULE transport_models_1d
  USE HDF5
  USE HDF5_io_module
  USE globals
  USE flux_surface_transport_data, ONLY: flux_surface_transport_t
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
     REAL*8 :: a_minor = 0.d0
     REAL*8 :: delta_te = 0.d0
     REAL*8, ALLOCATABLE :: te_fs(:)
     REAL*8, ALLOCATABLE :: ti_fs(:)
     REAL*8, ALLOCATABLE :: ne_fs(:)
     REAL*8, ALLOCATABLE :: pe_fs(:)
     REAL*8, ALLOCATABLE :: pi_fs(:)
     REAL*8, ALLOCATABLE :: q_fs(:)
     REAL*8, ALLOCATABLE :: Rmaj_fs(:)
     REAL*8, ALLOCATABLE :: rmin_fs(:)
     REAL*8, ALLOCATABLE :: eps_fs(:)
     REAL*8, ALLOCATABLE :: cs_te_fs(:)
     REAL*8, ALLOCATABLE :: dte_dr_fs(:)
     REAL*8, ALLOCATABLE :: dpe_dr_fs(:)
   CONTAINS
     PROCEDURE :: init => tm1d_init
     PROCEDURE :: destroy => tm1d_destroy
     PROCEDURE :: set_config => tm1d_set_config
     PROCEDURE :: update_from_flux_surfaces => tm1d_update_from_flux_surfaces
     PROCEDURE :: compute_delta_te => tm1d_compute_delta_te
     PROCEDURE :: write_hdf5 => tm1d_write_hdf5
     FINAL :: tm1d_finalize
  END TYPE transport_model_1d_t

  TYPE(transport_model_1d_t), SAVE :: transport_model_1d

CONTAINS

  SUBROUTINE tm1d_init(this, nrho, rho_edge, rho_core)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: nrho
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core

    CALL this%destroy()

    this%nrho = MAX(0, nrho)
    IF (PRESENT(rho_edge)) this%rho_edge = rho_edge
    IF (.NOT. PRESENT(rho_edge)) this%rho_edge = rho_edge_default
    IF (PRESENT(rho_core)) this%rho_core = rho_core
    IF (.NOT. PRESENT(rho_core)) this%rho_core = rho_core_default

    IF (this%nrho <= 0) RETURN

    ALLOCATE(this%te_fs(this%nrho))
    ALLOCATE(this%ti_fs(this%nrho))
    ALLOCATE(this%ne_fs(this%nrho))
    ALLOCATE(this%pe_fs(this%nrho))
    ALLOCATE(this%pi_fs(this%nrho))
    ALLOCATE(this%q_fs(this%nrho))
    ALLOCATE(this%Rmaj_fs(this%nrho))
    ALLOCATE(this%rmin_fs(this%nrho))
    ALLOCATE(this%eps_fs(this%nrho))
    ALLOCATE(this%cs_te_fs(this%nrho))
    ALLOCATE(this%dte_dr_fs(this%nrho))
    ALLOCATE(this%dpe_dr_fs(this%nrho))

    this%te_fs = 0.d0
    this%ti_fs = 0.d0
    this%ne_fs = 0.d0
    this%pe_fs = 0.d0
    this%pi_fs = 0.d0
    this%q_fs = 0.d0
    this%Rmaj_fs = 0.d0
    this%rmin_fs = 0.d0
    this%eps_fs = 0.d0
    this%cs_te_fs = 0.d0
    this%dte_dr_fs = 0.d0
    this%dpe_dr_fs = 0.d0
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
    IF (ALLOCATED(this%Rmaj_fs)) DEALLOCATE(this%Rmaj_fs)
    IF (ALLOCATED(this%rmin_fs)) DEALLOCATE(this%rmin_fs)
    IF (ALLOCATED(this%eps_fs)) DEALLOCATE(this%eps_fs)
    IF (ALLOCATED(this%cs_te_fs)) DEALLOCATE(this%cs_te_fs)
    IF (ALLOCATED(this%dte_dr_fs)) DEALLOCATE(this%dte_dr_fs)
    IF (ALLOCATED(this%dpe_dr_fs)) DEALLOCATE(this%dpe_dr_fs)

    this%is_initialized = .FALSE.
    this%nrho = 0
    this%rho_edge = rho_edge_default
    this%rho_core = rho_core_default
    this%a_minor = 0.d0
    this%delta_te = 0.d0
  END SUBROUTINE tm1d_destroy

  SUBROUTINE tm1d_finalize(this)
    TYPE(transport_model_1d_t), INTENT(INOUT) :: this

    CALL this%destroy()
  END SUBROUTINE tm1d_finalize

  SUBROUTINE tm1d_set_config(this, rho_edge, rho_core)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge, rho_core

    IF (PRESENT(rho_edge)) this%rho_edge = rho_edge
    IF (PRESENT(rho_core)) this%rho_core = rho_core
  END SUBROUTINE tm1d_set_config

  SUBROUTINE tm1d_update_from_flux_surfaces(this, fs_data)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    REAL*8, ALLOCATABLE :: ua(:, :), up(:, :)

    IF (.NOT. fs_data%profiles_built) RETURN

    IF ((.NOT. this%is_initialized) .OR. this%nrho /= fs_data%nrho) THEN
       CALL this%init(fs_data%nrho, this%rho_edge, this%rho_core)
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
    this%Rmaj_fs = fs_data%Rmaj_fs
    this%rmin_fs = fs_data%rmin_fs
    this%eps_fs = fs_data%eps_fs
    this%a_minor = phys%a_minor
    this%cs_te_fs = SQRT(MAX(this%te_fs*phys%Mref, model_tol))

    CALL tm1d_build_projected_gradients(this, fs_data)
    CALL this%compute_delta_te(fs_data)

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
    CALL HDF5_array1D_saving(group_id, this%Rmaj_fs, SIZE(this%Rmaj_fs), 'Rmaj_fs')
    CALL HDF5_array1D_saving(group_id, this%rmin_fs, SIZE(this%rmin_fs), 'rmin_fs')
    CALL HDF5_array1D_saving(group_id, this%eps_fs, SIZE(this%eps_fs), 'eps_fs')
    CALL HDF5_array1D_saving(group_id, this%cs_te_fs, SIZE(this%cs_te_fs), 'cs_te_fs')
    CALL HDF5_array1D_saving(group_id, this%dte_dr_fs, SIZE(this%dte_dr_fs), 'dte_dr_fs')
    CALL HDF5_array1D_saving(group_id, this%dpe_dr_fs, SIZE(this%dpe_dr_fs), 'dpe_dr_fs')
    CALL HDF5_real_saving(group_id, this%a_minor, 'a_minor')
    CALL HDF5_real_saving(group_id, this%rho_core, 'rho_core')
    CALL HDF5_real_saving(group_id, this%rho_edge, 'rho_edge')
    CALL HDF5_real_saving(group_id, this%delta_te, 'delta_te')
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

END MODULE transport_models_1d
