MODULE flux_surface_transport_data
  USE globals
  USE MPI_OMP
  USE HDF5_io_module
  USE interpolation, ONLY: find_cell_and_local_coordinate
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: flux_surface_transport_t, fs_transport

  REAL*8, PARAMETER :: pi_fs = 4.d0*DATAN(1.d0)
  REAL*8, PARAMETER :: rho_step_default = 1.d-2
  REAL*8, PARAMETER :: shell_weight_tol = 1.d-12
  REAL*8, PARAMETER :: gradpsi_tol = 1.d-12

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
     REAL*8, ALLOCATABLE :: Q_rad_sum(:, :)
     REAL*8, ALLOCATABLE :: q_sum(:)
     REAL*8, ALLOCATABLE :: Rmaj_sum(:)
     REAL*8, ALLOCATABLE :: rmin_sum(:)
     REAL*8, ALLOCATABLE :: U_fs(:, :)
     REAL*8, ALLOCATABLE :: Q_rad_fs(:, :)
     REAL*8, ALLOCATABLE :: q_fs(:)
     REAL*8, ALLOCATABLE :: Rmaj_fs(:)
     REAL*8, ALLOCATABLE :: rmin_fs(:)
     REAL*8, ALLOCATABLE :: eps_fs(:)
   CONTAINS
     PROCEDURE :: init => fs_init
     PROCEDURE :: destroy => fs_destroy
     PROCEDURE :: reset_accumulators => fs_reset_accumulators
     PROCEDURE :: build_profiles => fs_build_profiles
     PROCEDURE :: reduce_profile_sums => fs_reduce_profile_sums
     PROCEDURE :: finalize_profiles => fs_finalize_profiles
     PROCEDURE :: interp_U => fs_interp_U
     PROCEDURE :: interp_Q_rad => fs_interp_Q_rad
     PROCEDURE :: interp_scalar => fs_interp_scalar
     PROCEDURE :: write_hdf5 => fs_write_hdf5
     FINAL :: fs_finalize
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
    ALLOCATE(this%Q_rad_sum(this%neq, this%nrho))
    ALLOCATE(this%q_sum(this%nrho))
    ALLOCATE(this%Rmaj_sum(this%nrho))
    ALLOCATE(this%rmin_sum(this%nrho))
    ALLOCATE(this%U_fs(this%neq, this%nrho))
    ALLOCATE(this%Q_rad_fs(this%neq, this%nrho))
    ALLOCATE(this%q_fs(this%nrho))
    ALLOCATE(this%Rmaj_fs(this%nrho))
    ALLOCATE(this%rmin_fs(this%nrho))
    ALLOCATE(this%eps_fs(this%nrho))

    DO i = 1, this%nrho
       this%rho_grid(i) = (i - 1)*this%drho
    END DO

    this%is_initialized = .TRUE.
    CALL this%reset_accumulators()
  END SUBROUTINE fs_init

  SUBROUTINE fs_destroy(this)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%rho_grid)) DEALLOCATE(this%rho_grid)
    IF (ALLOCATED(this%shell_weight)) DEALLOCATE(this%shell_weight)
    IF (ALLOCATED(this%U_sum)) DEALLOCATE(this%U_sum)
    IF (ALLOCATED(this%Q_rad_sum)) DEALLOCATE(this%Q_rad_sum)
    IF (ALLOCATED(this%q_sum)) DEALLOCATE(this%q_sum)
    IF (ALLOCATED(this%Rmaj_sum)) DEALLOCATE(this%Rmaj_sum)
    IF (ALLOCATED(this%rmin_sum)) DEALLOCATE(this%rmin_sum)
    IF (ALLOCATED(this%U_fs)) DEALLOCATE(this%U_fs)
    IF (ALLOCATED(this%Q_rad_fs)) DEALLOCATE(this%Q_rad_fs)
    IF (ALLOCATED(this%q_fs)) DEALLOCATE(this%q_fs)
    IF (ALLOCATED(this%Rmaj_fs)) DEALLOCATE(this%Rmaj_fs)
    IF (ALLOCATED(this%rmin_fs)) DEALLOCATE(this%rmin_fs)
    IF (ALLOCATED(this%eps_fs)) DEALLOCATE(this%eps_fs)

    this%is_initialized = .FALSE.
    this%profiles_built = .FALSE.
    this%neq = 0
    this%nrho = 0
    this%rho_max = 0.d0
    this%drho = rho_step_default
  END SUBROUTINE fs_destroy

  SUBROUTINE fs_finalize(this)
    TYPE(flux_surface_transport_t), INTENT(INOUT) :: this

    CALL this%destroy()
  END SUBROUTINE fs_finalize

  SUBROUTINE fs_reset_accumulators(this)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this

    IF (.NOT. this%is_initialized) RETURN

    this%shell_weight = 0.d0
    this%U_sum = 0.d0
    this%Q_rad_sum = 0.d0
    this%q_sum = 0.d0
    this%Rmaj_sum = 0.d0
    this%rmin_sum = 0.d0
    this%U_fs = 0.d0
    this%Q_rad_fs = 0.d0
    this%q_fs = 0.d0
    this%Rmaj_fs = 0.d0
    this%rmin_fs = 0.d0
    this%eps_fs = 0.d0
    this%profiles_built = .FALSE.
  END SUBROUTINE fs_reset_accumulators

  SUBROUTINE fs_build_profiles(this)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this
    REAL*8, ALLOCATABLE :: ures(:, :), qres(:, :)
    REAL*8 :: rho_max_glob

    CALL fs_reshape_solution_fields(ures, qres)
    rho_max_glob = fs_compute_rho_max()
    CALL fs_ensure_grid(this, rho_max_glob)
    CALL fs_accumulate_profiles(this, ures, qres)
    CALL this%reduce_profile_sums()
    CALL this%finalize_profiles()

    DEALLOCATE(ures, qres)
  END SUBROUTINE fs_build_profiles

  SUBROUTINE fs_reduce_profile_sums(this)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this
#ifdef PARALL
    INTEGER :: ierr

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%shell_weight, this%nrho, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%U_sum, this%neq*this%nrho, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%Q_rad_sum, this%neq*this%nrho, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%q_sum, this%nrho, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%Rmaj_sum, this%nrho, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%rmin_sum, this%nrho, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE fs_reduce_profile_sums

  SUBROUTINE fs_finalize_profiles(this)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this
    INTEGER :: irho

    IF (.NOT. this%is_initialized) RETURN

    this%U_fs = 0.d0
    this%Q_rad_fs = 0.d0
    this%q_fs = 0.d0
    this%Rmaj_fs = 0.d0
    this%rmin_fs = 0.d0
    this%eps_fs = 0.d0
    DO irho = 1, this%nrho
       IF (this%shell_weight(irho) > shell_weight_tol) THEN
          this%U_fs(:, irho) = this%U_sum(:, irho)/this%shell_weight(irho)
          this%Q_rad_fs(:, irho) = this%Q_rad_sum(:, irho)/this%shell_weight(irho)
          this%q_fs(irho) = this%q_sum(irho)/this%shell_weight(irho)
          this%Rmaj_fs(irho) = this%Rmaj_sum(irho)/this%shell_weight(irho)
          this%rmin_fs(irho) = this%rmin_sum(irho)/this%shell_weight(irho)
          this%eps_fs(irho) = this%rmin_fs(irho)/MAX(this%Rmaj_fs(irho), shell_weight_tol)
       END IF
    END DO
    this%profiles_built = .TRUE.
  END SUBROUTINE fs_finalize_profiles


  SUBROUTINE fs_interp_U(this, rho, U_rho)
    CLASS(flux_surface_transport_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(OUT) :: U_rho(:)

    CALL fs_interp_profile(this, rho, this%U_fs, U_rho)
  END SUBROUTINE fs_interp_U

  SUBROUTINE fs_interp_Q_rad(this, rho, Q_rad_rho)
    CLASS(flux_surface_transport_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(OUT) :: Q_rad_rho(:)

    CALL fs_interp_profile(this, rho, this%Q_rad_fs, Q_rad_rho)
  END SUBROUTINE fs_interp_Q_rad

  SUBROUTINE fs_interp_scalar(this, rho, profile_fs, value)
    CLASS(flux_surface_transport_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(IN) :: profile_fs(:)
    REAL*8, INTENT(OUT) :: value
    INTEGER :: ilow, ihigh
    REAL*8 :: alpha, rho_eval, rho_low, rho_high

    value = 0.d0

    IF (SIZE(profile_fs) /= this%nrho) RETURN
    IF (.NOT. this%profiles_built) RETURN
    IF (this%nrho <= 0) RETURN

    IF (this%nrho == 1) THEN
       IF (this%shell_weight(1) > shell_weight_tol) value = profile_fs(1)
       RETURN
    END IF

    CALL find_cell_and_local_coordinate(this%nrho, this%rho_grid, rho, ilow, alpha)
    ihigh = MIN(ilow + 1, this%nrho)

    IF (this%shell_weight(ilow) > shell_weight_tol .AND. this%shell_weight(ihigh) > shell_weight_tol) THEN
       value = (1.d0 - alpha)*profile_fs(ilow) + alpha*profile_fs(ihigh)
       RETURN
    END IF

    ilow = fs_find_valid_left(this, ilow)
    ihigh = fs_find_valid_right(this, ihigh)

    IF (ilow <= 0 .AND. ihigh <= 0) RETURN
    IF (ilow <= 0) THEN
       value = profile_fs(ihigh)
       RETURN
    END IF
    IF (ihigh <= 0) THEN
       value = profile_fs(ilow)
       RETURN
    END IF
    IF (ilow == ihigh) THEN
       value = profile_fs(ilow)
       RETURN
    END IF

    rho_low = this%rho_grid(ilow)
    rho_high = this%rho_grid(ihigh)
    rho_eval = MIN(MAX(rho, rho_low), rho_high)
    IF (rho_high <= rho_low + gradpsi_tol) THEN
       value = profile_fs(ilow)
       RETURN
    END IF

    alpha = (rho_eval - rho_low)/(rho_high - rho_low)
    value = (1.d0 - alpha)*profile_fs(ilow) + alpha*profile_fs(ihigh)
  END SUBROUTINE fs_interp_scalar

  SUBROUTINE fs_write_hdf5(this, parent_group_id)
    CLASS(flux_surface_transport_t), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    IF (.NOT. switch%save_reduced_profiles_1D) RETURN
    IF (.NOT. this%profiles_built) RETURN

    CALL HDF5_group_create('transport_1d', parent_group_id, group_id, ierr)
    CALL HDF5_array1D_saving(group_id, this%rho_grid, SIZE(this%rho_grid), 'rho_grid')
    CALL HDF5_array1D_saving(group_id, this%shell_weight, SIZE(this%shell_weight), 'shell_weight')
    CALL HDF5_array2D_saving(group_id, this%U_fs, SIZE(this%U_fs, 1), SIZE(this%U_fs, 2), 'U_fs')
    CALL HDF5_array2D_saving(group_id, this%Q_rad_fs, SIZE(this%Q_rad_fs, 1), SIZE(this%Q_rad_fs, 2), 'Q_rad_fs')
    CALL HDF5_array1D_saving(group_id, this%q_fs, SIZE(this%q_fs), 'q_fs')
    CALL HDF5_array1D_saving(group_id, this%Rmaj_fs, SIZE(this%Rmaj_fs), 'Rmaj_fs')
    CALL HDF5_array1D_saving(group_id, this%rmin_fs, SIZE(this%rmin_fs), 'rmin_fs')
    CALL HDF5_array1D_saving(group_id, this%eps_fs, SIZE(this%eps_fs), 'eps_fs')
    CALL HDF5_group_close(group_id, ierr)

    IF (MPIvar%glob_id == 0 .AND. utils%printint > 0) THEN
       WRITE (6, *) 'Saved reduced 1D flux-surface profiles under /solution/transport_1d'
    END IF
  END SUBROUTINE fs_write_hdf5

  SUBROUTINE fs_reshape_solution_fields(ures, qres)
    REAL*8, ALLOCATABLE, INTENT(OUT) :: ures(:, :), qres(:, :)
    INTEGER :: sizeu

    sizeu = SIZE(sol%u)
    ALLOCATE(ures(sizeu/phys%Neq, phys%Neq))
    ALLOCATE(qres(sizeu/phys%Neq, phys%Neq*Mesh%Ndim))
    ures = TRANSPOSE(RESHAPE(sol%u, [phys%Neq, sizeu/phys%Neq]))
    qres = TRANSPOSE(RESHAPE(sol%q, [phys%Neq*Mesh%Ndim, sizeu/phys%Neq]))
  END SUBROUTINE fs_reshape_solution_fields

  FUNCTION fs_compute_rho_max() RESULT(rho_max_glob)
    REAL*8 :: rho_max_glob
    REAL*8 :: rho_max_local
    INTEGER :: iel
#ifdef PARALL
    INTEGER :: ierr
#endif

    rho_max_local = 0.d0
    DO iel = 1, Mesh%Nelems
       IF (.NOT. fs_is_local_element(iel)) CYCLE
       rho_max_local = MAX(rho_max_local, SQRT(MAX(0.d0, MAXVAL(phys%magnetic_psi(Mesh%T(iel, :))))))
    END DO

    rho_max_glob = rho_max_local
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, rho_max_glob, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
  END FUNCTION fs_compute_rho_max

  SUBROUTINE fs_ensure_grid(this, rho_max_glob)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: rho_max_glob

    IF ((.NOT. this%is_initialized) .OR. ABS(rho_max_glob - this%rho_max) > 0.5d0*this%drho) THEN
       CALL this%init(phys%Neq, rho_max_glob, rho_step_default)
    ELSE
       CALL this%reset_accumulators()
    END IF
  END SUBROUTINE fs_ensure_grid

  SUBROUTINE fs_accumulate_profiles(this, ures, qres)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: ures(:, :), qres(:, :)
    INTEGER :: iel

    DO iel = 1, Mesh%Nelems
       IF (.NOT. fs_is_local_element(iel)) CYCLE
       CALL fs_accumulate_element(this, iel, ures, qres)
    END DO
  END SUBROUTINE fs_accumulate_profiles

  SUBROUTINE fs_accumulate_element(this, iel, ures, qres)
    CLASS(flux_surface_transport_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: iel
    REAL*8, INTENT(IN) :: ures(:, :), qres(:, :)
    INTEGER :: g, ieq, irho
    REAL*8 :: rho_g, weight_g, dpsi_dxi, dpsi_deta, gradpsi_norm
    REAL*8 :: Xel(refElPol%Nnodes2D, 2)
    REAL*8 :: psiel(refElPol%Nnodes2D), q_cylel(refElPol%Nnodes2D)
    REAL*8 :: ue(refElPol%Nnodes2D, phys%Neq)
    REAL*8 :: qe(refElPol%Nnodes2D, phys%Neq*Mesh%Ndim)
    REAL*8 :: xy(refElPol%NGauss2D, 2)
    REAL*8 :: psig(refElPol%NGauss2D), q_cylg(refElPol%NGauss2D)
    REAL*8 :: ueg(refElPol%NGauss2D, phys%Neq)
    REAL*8 :: qeg(refElPol%NGauss2D, phys%Neq*Mesh%Ndim)
    REAL*8 :: J11(refElPol%NGauss2D), J12(refElPol%NGauss2D), J21(refElPol%NGauss2D), J22(refElPol%NGauss2D)
    REAL*8 :: detJ(refElPol%NGauss2D), iJ11(refElPol%NGauss2D), iJ12(refElPol%NGauss2D), iJ21(refElPol%NGauss2D), iJ22(refElPol%NGauss2D)
    REAL*8 :: gradpsi(2), npsi(2), gradu(2)

    Xel = Mesh%X(Mesh%T(iel, :), :)
    psiel = phys%magnetic_psi(Mesh%T(iel, :))
    q_cylel = phys%q_cyl(Mesh%T(iel, :))
    ue = ures((iel - 1)*refElPol%Nnodes2D + 1:iel*refElPol%Nnodes2D, :)
    qe = qres((iel - 1)*refElPol%Nnodes2D + 1:iel*refElPol%Nnodes2D, :)

    xy = MATMUL(refElPol%N2D, Xel)
    psig = MATMUL(refElPol%N2D, psiel)
    q_cylg = MATMUL(refElPol%N2D, q_cylel)
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
       rho_g = SQRT(MAX(0.d0, psig(g)))
       irho = FLOOR(rho_g/this%drho) + 1
       irho = MIN(MAX(irho, 1), this%nrho)
       weight_g = refElPol%gauss_weights2D(g)*ABS(detJ(g))
       IF (switch%axisym) weight_g = 2.d0*pi_fs*xy(g, 1)*weight_g

       this%shell_weight(irho) = this%shell_weight(irho) + weight_g
       this%U_sum(:, irho) = this%U_sum(:, irho) + ueg(g, :)*weight_g
       this%q_sum(irho) = this%q_sum(irho) + q_cylg(g)*weight_g
       this%Rmaj_sum(irho) = this%Rmaj_sum(irho) + xy(g, 1)*weight_g
       this%rmin_sum(irho) = this%rmin_sum(irho) + SQRT((xy(g, 1) - phys%r_axis)**2 + (xy(g, 2) - phys%z_axis)**2)*weight_g

       dpsi_dxi = DOT_PRODUCT(refElPol%Nxi2D(g, :), psiel)
       dpsi_deta = DOT_PRODUCT(refElPol%Neta2D(g, :), psiel)
       gradpsi(1) = iJ11(g)*dpsi_dxi + iJ12(g)*dpsi_deta
       gradpsi(2) = iJ21(g)*dpsi_dxi + iJ22(g)*dpsi_deta
       gradpsi_norm = SQRT(DOT_PRODUCT(gradpsi, gradpsi))

       IF (gradpsi_norm > gradpsi_tol) THEN
          npsi = gradpsi/gradpsi_norm
       ELSE
          npsi = 0.d0
       END IF

       DO ieq = 1, this%neq
          gradu(1) = qeg(g, (ieq - 1)*Mesh%Ndim + 1)
          gradu(2) = qeg(g, (ieq - 1)*Mesh%Ndim + 2)
          this%Q_rad_sum(ieq, irho) = this%Q_rad_sum(ieq, irho) + DOT_PRODUCT(gradu, npsi)*weight_g
       END DO
    END DO
  END SUBROUTINE fs_accumulate_element


  SUBROUTINE fs_interp_profile(this, rho, profile_fs, profile_rho)
    CLASS(flux_surface_transport_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(IN) :: profile_fs(:, :)
    REAL*8, INTENT(OUT) :: profile_rho(:)
    INTEGER :: ilow, ihigh
    REAL*8 :: alpha, rho_eval, rho_low, rho_high

    profile_rho = 0.d0

    IF (SIZE(profile_rho) /= SIZE(profile_fs, 1)) RETURN
    IF (.NOT. this%profiles_built) RETURN
    IF (this%nrho <= 0) RETURN

    IF (this%nrho == 1) THEN
       IF (this%shell_weight(1) > shell_weight_tol) profile_rho = profile_fs(:, 1)
       RETURN
    END IF

    CALL find_cell_and_local_coordinate(this%nrho, this%rho_grid, rho, ilow, alpha)
    ihigh = MIN(ilow + 1, this%nrho)

    IF (this%shell_weight(ilow) > shell_weight_tol .AND. this%shell_weight(ihigh) > shell_weight_tol) THEN
       profile_rho = (1.d0 - alpha)*profile_fs(:, ilow) + alpha*profile_fs(:, ihigh)
       RETURN
    END IF

    ilow = fs_find_valid_left(this, ilow)
    ihigh = fs_find_valid_right(this, ihigh)

    IF (ilow <= 0 .AND. ihigh <= 0) RETURN
    IF (ilow <= 0) THEN
       profile_rho = profile_fs(:, ihigh)
       RETURN
    END IF
    IF (ihigh <= 0) THEN
       profile_rho = profile_fs(:, ilow)
       RETURN
    END IF
    IF (ilow == ihigh) THEN
       profile_rho = profile_fs(:, ilow)
       RETURN
    END IF

    rho_low = this%rho_grid(ilow)
    rho_high = this%rho_grid(ihigh)
    rho_eval = MIN(MAX(rho, rho_low), rho_high)
    IF (rho_high <= rho_low + gradpsi_tol) THEN
       profile_rho = profile_fs(:, ilow)
       RETURN
    END IF

    alpha = (rho_eval - rho_low)/(rho_high - rho_low)
    profile_rho = (1.d0 - alpha)*profile_fs(:, ilow) + alpha*profile_fs(:, ihigh)
  END SUBROUTINE fs_interp_profile

  INTEGER FUNCTION fs_find_valid_left(this, istart)
    CLASS(flux_surface_transport_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: istart
    INTEGER :: i

    fs_find_valid_left = 0
    DO i = MIN(MAX(istart, 1), this%nrho), 1, -1
       IF (this%shell_weight(i) > shell_weight_tol) THEN
          fs_find_valid_left = i
          RETURN
       END IF
    END DO
  END FUNCTION fs_find_valid_left

  INTEGER FUNCTION fs_find_valid_right(this, istart)
    CLASS(flux_surface_transport_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: istart
    INTEGER :: i

    fs_find_valid_right = 0
    DO i = MIN(MAX(istart, 1), this%nrho), this%nrho
       IF (this%shell_weight(i) > shell_weight_tol) THEN
          fs_find_valid_right = i
          RETURN
       END IF
    END DO
  END FUNCTION fs_find_valid_right

  LOGICAL FUNCTION fs_is_local_element(iel)
    INTEGER, INTENT(IN) :: iel

    fs_is_local_element = .TRUE.
#ifdef PARALL
    IF (ASSOCIATED(Mesh%ghostElems)) THEN
       fs_is_local_element = Mesh%ghostElems(iel) == 0
    END IF
#endif
  END FUNCTION fs_is_local_element

END MODULE flux_surface_transport_data
