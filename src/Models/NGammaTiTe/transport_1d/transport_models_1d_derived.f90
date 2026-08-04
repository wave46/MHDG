MODULE transport_models_1d_derived
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: transport_model_derived_t

  TYPE :: transport_model_derived_t
     REAL*8 :: a_minor = 0.d0
     REAL*8 :: delta_te = 0.d0
     REAL*8, ALLOCATABLE :: cs_te_fs(:)
     REAL*8, ALLOCATABLE :: ne_fs(:)
     REAL*8, ALLOCATABLE :: pe_fs(:)
     REAL*8, ALLOCATABLE :: te_fs(:)
     REAL*8, ALLOCATABLE :: q_fs(:)
     REAL*8, ALLOCATABLE :: omega_fs(:)
     REAL*8, ALLOCATABLE :: Rmaj_fs(:)
     REAL*8, ALLOCATABLE :: eps_fs(:)
     REAL*8, ALLOCATABLE :: nuestar_fs(:)
     REAL*8, ALLOCATABLE :: dpe_dr_fs(:)
     REAL*8, ALLOCATABLE :: dte_dr_fs(:)
     REAL*8, ALLOCATABLE :: chi_bohm_fs(:)
     REAL*8, ALLOCATABLE :: chi_gyrobohm_fs(:)
   CONTAINS
     PROCEDURE :: init => tm1d_workspace_init
     PROCEDURE :: destroy => tm1d_workspace_destroy
     FINAL :: tm1d_workspace_finalize
  END TYPE transport_model_derived_t

CONTAINS

  SUBROUTINE tm1d_workspace_init(this, nrho)
    CLASS(transport_model_derived_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: nrho

    CALL this%destroy()
    IF (nrho <= 0) RETURN

    this%a_minor = 0.d0
    this%delta_te = 0.d0

    ALLOCATE(this%cs_te_fs(nrho))
    ALLOCATE(this%ne_fs(nrho))
    ALLOCATE(this%pe_fs(nrho))
    ALLOCATE(this%te_fs(nrho))
    ALLOCATE(this%q_fs(nrho))
    ALLOCATE(this%omega_fs(nrho))
    ALLOCATE(this%Rmaj_fs(nrho))
    ALLOCATE(this%eps_fs(nrho))
    ALLOCATE(this%nuestar_fs(nrho))
    ALLOCATE(this%dpe_dr_fs(nrho))
    ALLOCATE(this%dte_dr_fs(nrho))
    ALLOCATE(this%chi_bohm_fs(nrho))
    ALLOCATE(this%chi_gyrobohm_fs(nrho))

    this%cs_te_fs = 0.d0
    this%ne_fs = 0.d0
    this%pe_fs = 0.d0
    this%te_fs = 0.d0
    this%q_fs = 0.d0
    this%omega_fs = 0.d0
    this%Rmaj_fs = 0.d0
    this%eps_fs = 0.d0
    this%nuestar_fs = 0.d0
    this%dpe_dr_fs = 0.d0
    this%dte_dr_fs = 0.d0
    this%chi_bohm_fs = 0.d0
    this%chi_gyrobohm_fs = 0.d0
  END SUBROUTINE tm1d_workspace_init

  SUBROUTINE tm1d_workspace_destroy(this)
    CLASS(transport_model_derived_t), INTENT(INOUT) :: this

    this%a_minor = 0.d0
    this%delta_te = 0.d0

    IF (ALLOCATED(this%cs_te_fs)) DEALLOCATE(this%cs_te_fs)
    IF (ALLOCATED(this%ne_fs)) DEALLOCATE(this%ne_fs)
    IF (ALLOCATED(this%pe_fs)) DEALLOCATE(this%pe_fs)
    IF (ALLOCATED(this%te_fs)) DEALLOCATE(this%te_fs)
    IF (ALLOCATED(this%q_fs)) DEALLOCATE(this%q_fs)
    IF (ALLOCATED(this%omega_fs)) DEALLOCATE(this%omega_fs)
    IF (ALLOCATED(this%Rmaj_fs)) DEALLOCATE(this%Rmaj_fs)
    IF (ALLOCATED(this%eps_fs)) DEALLOCATE(this%eps_fs)
    IF (ALLOCATED(this%nuestar_fs)) DEALLOCATE(this%nuestar_fs)
    IF (ALLOCATED(this%dpe_dr_fs)) DEALLOCATE(this%dpe_dr_fs)
    IF (ALLOCATED(this%dte_dr_fs)) DEALLOCATE(this%dte_dr_fs)
    IF (ALLOCATED(this%chi_bohm_fs)) DEALLOCATE(this%chi_bohm_fs)
    IF (ALLOCATED(this%chi_gyrobohm_fs)) DEALLOCATE(this%chi_gyrobohm_fs)
  END SUBROUTINE tm1d_workspace_destroy

  SUBROUTINE tm1d_workspace_finalize(this)
    TYPE(transport_model_derived_t), INTENT(INOUT) :: this

    CALL this%destroy()
  END SUBROUTINE tm1d_workspace_finalize

END MODULE transport_models_1d_derived
