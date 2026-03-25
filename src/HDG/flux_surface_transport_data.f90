MODULE flux_surface_transport_data
  USE globals
  USE MPI_OMP
  USE HDF5_io_module
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: flux_surface_transport_t, fs_transport

  REAL*8, PARAMETER :: pi_fs = 4.d0*DATAN(1.d0)
  REAL*8, PARAMETER :: rho_step_default = 1.d-2
  REAL*8, PARAMETER :: rho_tol = 1.d-12

  TYPE :: flux_surface_transport_t
     LOGICAL :: is_initialized = .FALSE.
     LOGICAL :: profiles_built = .FALSE.
     INTEGER :: neq = 0
     INTEGER :: nrho = 0
     REAL*8  :: drho = rho_step_default
     REAL*8  :: rho_max = 0.d0
     REAL*8, ALLOCATABLE :: rho_grid(:)
     REAL*8, ALLOCATABLE :: shell_weight(:)
     REAL*8, ALLOCATABLE :: U_sum(:, :)
     REAL*8, ALLOCATABLE :: Q_sum(:, :)
     REAL*8, ALLOCATABLE :: U_fs(:, :)
     REAL*8, ALLOCATABLE :: Q_fs(:, :)
   CONTAINS
     PROCEDURE :: init => fs_init
     PROCEDURE :: destroy => fs_destroy
     PROCEDURE :: reset_accumulators => fs_reset_accumulators
     PROCEDURE :: build_profiles => fs_build_profiles
     PROCEDURE :: write_hdf5 => fs_write_hdf5
  END TYPE flux_surface_transport_t

  TYPE(flux_surface_transport_t), SAVE :: fs_transport

CONTAINS

  SUBROUTINE fs_init(this, neq, rho_max, drho)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: neq
    REAL*8, INTENT(IN) :: rho_max
    REAL*8, INTENT(IN), OPTIONAL :: drho
    INTEGER :: i

    CALL this%destroy()

    this%neq = neq
    this%drho = rho_step_default
    IF (PRESENT(drho)) this%drho = drho
    this%rho_max = MAX(0.d0, rho_max)
    this%nrho = MAX(1, INT(this%rho_max/this%drho + 0.5d0) + 1)

    ALLOCATE(this%rho_grid(this%nrho))
    ALLOCATE(this%shell_weight(this%nrho))
    ALLOCATE(this%U_sum(this%neq, this%nrho))
    ALLOCATE(this%Q_sum(this%neq, this%nrho))
    ALLOCATE(this%U_fs(this%neq, this%nrho))
    ALLOCATE(this%Q_fs(this%neq, this%nrho))

    DO i = 1, this%nrho
       this%rho_grid(i) = (i - 1)*this%drho
    END DO

    CALL this%reset_accumulators()
    this%is_initialized = .TRUE.
  END SUBROUTINE fs_init

  SUBROUTINE fs_destroy(this)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%rho_grid)) DEALLOCATE(this%rho_grid)
    IF (ALLOCATED(this%shell_weight)) DEALLOCATE(this%shell_weight)
    IF (ALLOCATED(this%U_sum)) DEALLOCATE(this%U_sum)
    IF (ALLOCATED(this%Q_sum)) DEALLOCATE(this%Q_sum)
    IF (ALLOCATED(this%U_fs)) DEALLOCATE(this%U_fs)
    IF (ALLOCATED(this%Q_fs)) DEALLOCATE(this%Q_fs)

    this%is_initialized = .FALSE.
    this%profiles_built = .FALSE.
    this%neq = 0
    this%nrho = 0
    this%rho_max = 0.d0
    this%drho = rho_step_default
  END SUBROUTINE fs_destroy

  SUBROUTINE fs_reset_accumulators(this)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this

    IF (.NOT. this%is_initialized) RETURN

    this%shell_weight = 0.d0
    this%U_sum = 0.d0
    this%Q_sum = 0.d0
    this%U_fs = 0.d0
    this%Q_fs = 0.d0
    this%profiles_built = .FALSE.
  END SUBROUTINE fs_reset_accumulators

  SUBROUTINE fs_build_profiles(this)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this
    INTEGER :: iel, g, ieq, irho
    INTEGER :: sizeu
    REAL*8 :: rho_max_local, rho_max_glob
    REAL*8 :: weight_g, rho_g, detJg, dpsi_dxi, dpsi_deta, gradpsi_norm
    REAL*8 :: Xel(refElPol%Nnodes2D, 2)
    REAL*8 :: psiel(refElPol%Nnodes2D)
    REAL*8 :: ue(refElPol%Nnodes2D, phys%Neq)
    REAL*8 :: qe(refElPol%Nnodes2D, phys%Neq*Mesh%Ndim)
    REAL*8 :: xy(refElPol%NGauss2D, 2)
    REAL*8 :: Psig(refElPol%NGauss2D)
    REAL*8 :: ueg(refElPol%NGauss2D, phys%Neq)
    REAL*8 :: qeg(refElPol%NGauss2D, phys%Neq*Mesh%Ndim)
    REAL*8 :: J11(refElPol%NGauss2D), J12(refElPol%NGauss2D), J21(refElPol%NGauss2D), J22(refElPol%NGauss2D)
    REAL*8 :: detJ(refElPol%NGauss2D), iJ11(refElPol%NGauss2D), iJ12(refElPol%NGauss2D), iJ21(refElPol%NGauss2D), iJ22(refElPol%NGauss2D)
    REAL*8 :: gradpsi(2), npsi(2), gradu(2)
    REAL*8, ALLOCATABLE :: ures(:, :), qres(:, :)
#ifdef PARALL
    INTEGER :: ierr
#endif

    sizeu = SIZE(sol%u)
    ALLOCATE(ures(sizeu/phys%Neq, phys%Neq))
    ALLOCATE(qres(sizeu/phys%Neq, phys%Neq*Mesh%Ndim))
    ures = TRANSPOSE(RESHAPE(sol%u, [phys%Neq, sizeu/phys%Neq]))
    qres = TRANSPOSE(RESHAPE(sol%q, [phys%Neq*Mesh%Ndim, sizeu/phys%Neq]))

    rho_max_local = 0.d0
    DO iel = 1, Mesh%Nelems
#ifdef PARALL
       IF (ASSOCIATED(Mesh%ghostElems)) THEN
          IF (Mesh%ghostElems(iel) /= 0) CYCLE
       END IF
#endif
       psiel = phys%magnetic_psi(Mesh%T(iel, :))
       rho_max_local = MAX(rho_max_local, SQRT(MAX(0.d0, MAXVAL(psiel))))
    END DO

    rho_max_glob = rho_max_local
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, rho_max_glob, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif

    IF ((.NOT. this%is_initialized) .OR. ABS(rho_max_glob - this%rho_max) > 0.5d0*this%drho) THEN
       CALL this%init(phys%Neq, rho_max_glob, rho_step_default)
    ELSE
       CALL this%reset_accumulators()
    END IF

    DO iel = 1, Mesh%Nelems
#ifdef PARALL
       IF (ASSOCIATED(Mesh%ghostElems)) THEN
          IF (Mesh%ghostElems(iel) /= 0) CYCLE
       END IF
#endif
       Xel = Mesh%X(Mesh%T(iel, :), :)
       psiel = phys%magnetic_psi(Mesh%T(iel, :))
       ue = ures((iel - 1)*refElPol%Nnodes2D + 1:iel*refElPol%Nnodes2D, :)
       qe = qres((iel - 1)*refElPol%Nnodes2D + 1:iel*refElPol%Nnodes2D, :)

       xy = MATMUL(refElPol%N2D, Xel)
       Psig = MATMUL(refElPol%N2D, psiel)
       ueg = MATMUL(refElPol%N2D, ue)
       qeg = MATMUL(refElPol%N2D, qe)

       J11 = MATMUL(refElPol%Nxi2D, Xel(:, 1))
       J12 = MATMUL(refElPol%Nxi2D, Xel(:, 2))
       J21 = MATMUL(refElPol%Neta2D, Xel(:, 1))
       J22 = MATMUL(refElPol%Neta2D, Xel(:, 2))
       detJ = J11*J22 - J21*J12
       iJ11 = J22/detJ
       iJ12 = -J12/detJ
       iJ21 = -J21/detJ
       iJ22 = J11/detJ

       DO g = 1, refElPol%NGauss2D
          rho_g = SQRT(MAX(0.d0, Psig(g)))
          irho = MIN(MAX(NINT(rho_g/this%drho) + 1, 1), this%nrho)

          detJg = ABS(detJ(g))
          weight_g = refElPol%gauss_weights2D(g)*detJg
          IF (switch%axisym) weight_g = 2.d0*pi_fs*xy(g, 1)*weight_g
          weight_g = weight_g*phys%lscale**3

          this%shell_weight(irho) = this%shell_weight(irho) + weight_g
          this%U_sum(:, irho) = this%U_sum(:, irho) + ueg(g, :)*weight_g

          dpsi_dxi = DOT_PRODUCT(refElPol%Nxi2D(g, :), psiel)
          dpsi_deta = DOT_PRODUCT(refElPol%Neta2D(g, :), psiel)
          gradpsi(1) = iJ11(g)*dpsi_dxi + iJ12(g)*dpsi_deta
          gradpsi(2) = iJ21(g)*dpsi_dxi + iJ22(g)*dpsi_deta
          gradpsi_norm = SQRT(DOT_PRODUCT(gradpsi, gradpsi))

          IF (gradpsi_norm > rho_tol) THEN
             npsi = gradpsi/gradpsi_norm
          ELSE
             npsi = 0.d0
          END IF

          DO ieq = 1, this%neq
             gradu(1) = qeg(g, (ieq - 1)*Mesh%Ndim + 1)
             gradu(2) = qeg(g, (ieq - 1)*Mesh%Ndim + 2)
             this%Q_sum(ieq, irho) = this%Q_sum(ieq, irho) + DOT_PRODUCT(gradu, npsi)*weight_g
          END DO
       END DO
    END DO

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%shell_weight, this%nrho, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%U_sum, this%neq*this%nrho, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%Q_sum, this%neq*this%nrho, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    DO irho = 1, this%nrho
       IF (this%shell_weight(irho) > rho_tol) THEN
          this%U_fs(:, irho) = this%U_sum(:, irho)/this%shell_weight(irho)
          this%Q_fs(:, irho) = this%Q_sum(:, irho)/this%shell_weight(irho)
       END IF
    END DO

    this%profiles_built = .TRUE.

    DEALLOCATE(ures, qres)
  END SUBROUTINE fs_build_profiles

  SUBROUTINE fs_write_hdf5(this, file_or_group_id)
    CLASS(flux_surface_transport_t), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: file_or_group_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    IF (.NOT. this%profiles_built) RETURN

    CALL HDF5_group_create('transport_1d', file_or_group_id, group_id, ierr)
    CALL HDF5_array1D_saving(group_id, this%rho_grid, SIZE(this%rho_grid), 'rho_grid')
    CALL HDF5_array1D_saving(group_id, this%shell_weight, SIZE(this%shell_weight), 'shell_weight')
    CALL HDF5_array2D_saving(group_id, this%U_fs, SIZE(this%U_fs, 1), SIZE(this%U_fs, 2), 'U_fs')
    CALL HDF5_array2D_saving(group_id, this%Q_fs, SIZE(this%Q_fs, 1), SIZE(this%Q_fs, 2), 'Q_fs')
    CALL HDF5_group_close(group_id, ierr)
  END SUBROUTINE fs_write_hdf5

END MODULE flux_surface_transport_data
