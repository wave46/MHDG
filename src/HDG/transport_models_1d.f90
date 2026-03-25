MODULE transport_models_1d
  USE globals
  USE flux_surface_transport_data, ONLY: flux_surface_transport_t
  USE physics, ONLY: cons2phys
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: transport_model_1d_t, transport_model_1d

  REAL*8, PARAMETER :: rho_edge_default = 0.99d0
  REAL*8, PARAMETER :: model_tol = 1.d-12

  TYPE :: transport_model_1d_t
     LOGICAL :: is_initialized = .FALSE.
     INTEGER :: nrho = 0
     REAL*8 :: rho_edge = rho_edge_default
     REAL*8, ALLOCATABLE :: te_fs(:)
     REAL*8, ALLOCATABLE :: ti_fs(:)
     REAL*8, ALLOCATABLE :: ne_fs(:)
     REAL*8, ALLOCATABLE :: delta_te(:)
   CONTAINS
     PROCEDURE :: init => tm1d_init
     PROCEDURE :: destroy => tm1d_destroy
     PROCEDURE :: update_from_flux_surfaces => tm1d_update_from_flux_surfaces
     PROCEDURE :: compute_delta_te => tm1d_compute_delta_te
     FINAL :: tm1d_finalize
  END TYPE transport_model_1d_t

  TYPE(transport_model_1d_t), SAVE :: transport_model_1d

CONTAINS

  SUBROUTINE tm1d_init(this, nrho, rho_edge)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: nrho
    REAL*8, INTENT(IN), OPTIONAL :: rho_edge

    CALL this%destroy()

    this%nrho = MAX(0, nrho)
    IF (PRESENT(rho_edge)) this%rho_edge = rho_edge
    IF (.NOT. PRESENT(rho_edge)) this%rho_edge = rho_edge_default

    IF (this%nrho <= 0) RETURN

    ALLOCATE(this%te_fs(this%nrho))
    ALLOCATE(this%ti_fs(this%nrho))
    ALLOCATE(this%ne_fs(this%nrho))
    ALLOCATE(this%delta_te(this%nrho))

    this%te_fs = 0.d0
    this%ti_fs = 0.d0
    this%ne_fs = 0.d0
    this%delta_te = 0.d0
    this%is_initialized = .TRUE.
  END SUBROUTINE tm1d_init

  SUBROUTINE tm1d_destroy(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%te_fs)) DEALLOCATE(this%te_fs)
    IF (ALLOCATED(this%ti_fs)) DEALLOCATE(this%ti_fs)
    IF (ALLOCATED(this%ne_fs)) DEALLOCATE(this%ne_fs)
    IF (ALLOCATED(this%delta_te)) DEALLOCATE(this%delta_te)

    this%is_initialized = .FALSE.
    this%nrho = 0
    this%rho_edge = rho_edge_default
  END SUBROUTINE tm1d_destroy

  SUBROUTINE tm1d_finalize(this)
    TYPE(transport_model_1d_t), INTENT(INOUT) :: this

    CALL this%destroy()
  END SUBROUTINE tm1d_finalize

  SUBROUTINE tm1d_update_from_flux_surfaces(this, fs_data)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(flux_surface_transport_t), INTENT(IN) :: fs_data
    REAL*8, ALLOCATABLE :: ua(:, :), up(:, :)

    IF (.NOT. fs_data%profiles_built) RETURN

    IF ((.NOT. this%is_initialized) .OR. this%nrho /= fs_data%nrho) THEN
       CALL this%init(fs_data%nrho, this%rho_edge)
    END IF
    IF (.NOT. this%is_initialized) RETURN

    ALLOCATE(ua(fs_data%nrho, fs_data%neq))
    ALLOCATE(up(fs_data%nrho, phys%npv))

    ua = TRANSPOSE(fs_data%U_fs)
    CALL cons2phys(ua, up)

    this%ne_fs = up(:, 1)
    this%ti_fs = up(:, 7)
    this%te_fs = up(:, 8)

    CALL this%compute_delta_te(fs_data%rho_grid)

    DEALLOCATE(ua, up)
  END SUBROUTINE tm1d_update_from_flux_surfaces

  SUBROUTINE tm1d_compute_delta_te(this, rho_grid)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: rho_grid(:)
    INTEGER :: iedge
    REAL*8 :: te_edge

    this%delta_te = 0.d0
    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    iedge = tm1d_find_edge_index(this%rho_edge, rho_grid)
    te_edge = MAX(this%te_fs(iedge), model_tol)

    this%delta_te = (this%te_fs - te_edge)/te_edge
  END SUBROUTINE tm1d_compute_delta_te

  INTEGER FUNCTION tm1d_find_edge_index(rho_edge, rho_grid)
    REAL*8, INTENT(IN) :: rho_edge
    REAL*8, INTENT(IN) :: rho_grid(:)
    INTEGER :: i

    tm1d_find_edge_index = SIZE(rho_grid)
    DO i = 1, SIZE(rho_grid)
       IF (rho_grid(i) >= rho_edge) THEN
          tm1d_find_edge_index = i
          RETURN
       END IF
    END DO
  END FUNCTION tm1d_find_edge_index

END MODULE transport_models_1d
