MODULE solve_petsc
#include "petsc/finclude/petsc.h"
   USE matrices_types
   USE types
   USE globals
   USE Petsc
   USE MPI_OMP, ONLY: MPIvar, OMPvar

   IMPLICIT NONE

   TYPE PETSC_STRUC ! A type to store structures for the PETSC library
      Mat                                :: matK          ! sparse matrix for PETSC
      IS                                 :: is_num, is_perm
      Mat                                :: matKperm      ! permuted sparse matrix for PETSC
      Vec                                :: rhs_vec_perm, solPETSC_vec_perm
      Mat, DIMENSION(:), POINTER         :: matP, matR => NULL()       ! sparse matrices for PETSC (for preconditioner multigrid)
      KSP                                :: ksp           ! solver
      PC                                 :: pc            ! preconditioner
      Vec                                :: solPETSC_vec      ! solution
      Vec                                :: rhs_vec           ! right hand side
      PetscInt                           :: n, nnz        ! matrix n and number of non-zeros
      PetscInt, DIMENSION(:), POINTER    :: rowptr => NULL()
      PetscInt, DIMENSION(:), POINTER    :: loc2glob => NULL()
      PetscInt, DIMENSION(:), POINTER    :: cols => NULL()
      PetscScalar, DIMENSION(:), POINTER :: vals_matK => NULL()
      PetscScalar                        :: residue   ! UNPRECONDTIONED residue norm of the linear system
      PetscInt                           :: its      ! Previous number of iterations
      PetscInt                           :: its_recompute_prec ! Maximum number of iterations reached in the previous solve
      KSPConvergedReason                 :: convergedReason

      IS                        :: is_row_f, is_col_f        ! graph ordering of the fine level (computed once)
      IS, ALLOCATABLE           :: rowperm(:), colperm(:)    ! per-level permutations, 0:nl-1
      Mat, ALLOCATABLE          :: matKt(:)                  ! permuted operators per level, 0:nl-1
      Mat, ALLOCATABLE          :: matPt(:), matRt(:)        ! permuted transfer operators, 1:nl-1
      LOGICAL                   :: mg_perm_built = .FALSE.

   END TYPE PETSC_STRUC

   TYPE(PETSC_STRUC)            :: matPETSc

CONTAINS

   SUBROUTINE InitPETSC
      IMPLICIT NONE

      PetscErrorCode            :: ierr
      PetscBool                 :: initialized

      ! Check if PETSc is already initialized (e.g., by another libraUnable to read file '/home/mcapasso/Documents/MHDG/Work/MHDG_old/test/matK.m'. Input must be a MAT-file or an ASCII file containing numeric data with same number of columns in each row.ry) and if not, initialize it with some options
      CALL PetscInitialized(initialized, ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-mat_mkl_pardiso_65", '20', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-malloc_debug', '', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-malloc_test', '', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-malloc_dump', '', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-info', 'true', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-log_view', '', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-ksp_monitor_true_residual', 'ascii:ksp_residual.log', ierr); CHKERRA(ierr)
      ! Per-level residual logs: use to see whether the coarse solve or a specific level's smoother is stalling.
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mg_coarse_ksp_monitor_true_residual', 'ascii:mg_coarse_residual.log', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mg_levels_1_ksp_monitor_true_residual', 'ascii:mg_level1_residual.log', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mg_levels_2_ksp_monitor_true_residual', 'ascii:mg_level2_residual.log', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mg_levels_3_ksp_monitor_true_residual', 'ascii:mg_level3_residual.log', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-on_error_attach_debugger', 'true', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-start_in_debugger', 'true', ierr); CHKERRA(ierr)
      !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mat_no_inode', '', ierr); CHKERRA(ierr)
      IF (.NOT. initialized) THEN
         CALL PetscInitialize(ierr); CHKERRA(ierr)
         !call petscinitializenoarguments(ierr) ! Default communicator = MPI_COMM_WORLD
      END IF

      IF (ierr .NE. 0) THEN
         PRINT *, 'Unable to initialize PETSc even though it was requested. Aborting...'
         CALL MPI_Abort(MPI_COMM_WORLD, -1, ierr); CHKERRA(ierr)
      END IF

   END SUBROUTINE InitPETSC

   SUBROUTINE FinalizePETSC
      IMPLICIT NONE

      PetscErrorCode            :: ierr
      PetscBool                 :: finalized

      ! Check if PETSc is already finalized (e.g., by another library) and if not, finalize it
      CALL PetscFinalized(finalized, ierr); CHKERRA(ierr)
      IF (.NOT. finalized) THEN
         CALL PetscFinalize(ierr); CHKERRA(ierr)! Default communicator = MPI_COMM_WORLD
      END IF
      IF (ierr .NE. 0) THEN
         PRINT *, 'Unable to finalize PETSc even though it was requested. Aborting...'
         CALL MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      END IF

   END SUBROUTINE FinalizePETSC

   !***********************************************
   ! Initialization of matrix
   ! Part specific to PETSc
   !***********************************************
   SUBROUTINE init_mat_PETSC(matPETSC)
      TYPE(PETSC_STRUC)         :: matPETSC
      PetscErrorCode            :: ierr
      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke

      cks = 0
      cke = 0

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tps)
         CALL system_CLOCK(cks, clock_rate)
      END IF

      ! Link the PETSC instance with the CSR matrix storage
      IF (ASSOCIATED(matPETSC%rowptr)) THEN
         DEALLOCATE (matPETSC%rowptr)
      END IF
      IF (ASSOCIATED(matPETSC%loc2glob)) THEN
         DEALLOCATE (matPETSC%loc2glob)
      END IF
      IF (ASSOCIATED(matPETSC%cols)) THEN
         DEALLOCATE (matPETSC%cols)
      END IF
      IF (ASSOCIATED(matPETSC%vals_matK)) THEN
         DEALLOCATE (matPETSC%vals_matK)
      END IF

      matPETSC%n = matK%n
      matPETSC%nnz = matK%nnz
      ALLOCATE (matPETSC%rowptr(matPETSC%n + 1))
      ALLOCATE (matPETSC%loc2glob(matPETSC%n))
      ALLOCATE (matPETSC%cols(matPETSC%nnz))
      ALLOCATE (matPETSC%vals_matK(matPETSC%nnz))
      matPETSC%rowptr = matK%rowptr
      matPETSC%cols = matK%cols
      matPETSC%vals_matK = matK%vals
      matPETSC%loc2glob = matK%loc2glob

      IF (matK%start) THEN
         matPETSC%residue = 0
         matPETSC%its = 1
         matPETSC%its_recompute_prec = 1

         ! Create matrix and set it up
         CALL MatCreate(MPI_COMM_WORLD, matPETSC%matK, ierr); CHKERRA(ierr)
#ifndef PARALL
         CALL MatSetType(matPETSC%matK, MATSEQAIJ, ierr); CHKERRA(ierr)
         CALL MatSetSizes(matPETSC%matK, matPETSC%n, matPETSC%n, matPETSC%n, matPETSC%n, ierr); CHKERRA(ierr)
#else
         CALL MatSetType(matPETSC%matK, MATMPIAIJ, ierr); CHKERRA(ierr)
         CALL MatSetSizes(matPETSC%matK, matPETSC%n, matPETSC%n, PETSC_DETERMINE, PETSC_DETERMINE, ierr); CHKERRA(ierr)
#endif
         CALL MatSetBlockSize(matPETSC%matK, phys%neq, ierr); CHKERRA(ierr)
         CALL MatSetUp(matPETSC%matK, ierr); CHKERRA(ierr)
         ! Sparsity pattern is fixed across iterations: keep pattern and reject new entries.
         CALL MatSetOption(matPETSC%matK, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr); CHKERRA(ierr)
         CALL MatSetOption(matPETSC%matK, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr); CHKERRA(ierr)

         ! Create rhs and solution vector
         CALL MatCreateVecs(matPETSC%matK, matPETSC%rhs_vec, matPETSC%solPETSC_vec, ierr); CHKERRA(ierr)

         ! if multigrid preconditioner, create the matrices for the levels and the prolongator/restriction operators (setup and preallocation)
         IF (lssolver%pctype .EQ. 'PCMG') THEN
            CALL setup_MG()
         END IF

         ! Create the ksp instance
         CALL KSPCreate(PETSC_COMM_WORLD, matPETSC%ksp, ierr); CHKERRA(ierr)
         ! Build preconditioner based on linear system matrix
         CALL KSPSetOperators(matPETSC%ksp, matPETSC%matK, matPETSC%matK, ierr); CHKERRA(ierr)
      END IF

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tpe)
         CALL system_CLOCK(cke, clock_rate)
         timing%rlstime1 = timing%rlstime1 + (cke - cks)/REAL(clock_rate)
         timing%clstime1 = timing%clstime1 + tpe - tps
      END IF

   END SUBROUTINE init_mat_PETSC

   SUBROUTINE preallocate_matrix_PETSC(matPETSC)

      TYPE(PETSC_STRUC), INTENT(IN)            :: matPETSC
      INTEGER, ALLOCATABLE         :: counter_diag(:), counter_offdiag(:), istart_vec(:), iend_vec(:)
      PetscInt                     :: istart, iend, r, rr, rrr, M, N, I, J
      PetscErrorCode               :: ierr

      ! Get global sizes
      CALL MatGetSize(matPETSC%matK, M, N, ierr); CHKERRA(ierr)

      ! Get local row ownership of the matrix
      CALL MatGetOwnershipRange(matPETSC%matK, istart, iend, ierr); CHKERRA(ierr)

      ! allocate diagonal and offdiagonal non-zero counters and a vector containing the beginning and ends!
      ! of the row ownership, per process. Each process owns X rows of the matrix, a slice of matrix.
      ALLOCATE (counter_diag(N))
      ALLOCATE (counter_offdiag(N))
      ALLOCATE (istart_vec(MPIvar%glob_size))
      ALLOCATE (iend_vec(MPIvar%glob_size))

      counter_diag = 0
      counter_offdiag = 0
      istart_vec = 0
      iend_vec = 0
      istart_vec(MPIvar%glob_id + 1) = istart
      iend_vec(MPIvar%glob_id + 1) = iend

#ifdef PARALL
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, istart_vec, MPIvar%glob_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, iend_vec, MPIvar%glob_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      ! Count non-zero elements
      DO r = 1, matPETSC%n
         DO rr = matPETSC%rowptr(r), (matPETSC%rowptr(r + 1) - 1)
            I = matPETSC%loc2glob(r) - 1;
            J = matPETSC%cols(rr) - 1;
            ! Find out on which process' slice I am (when found, exit the loop, rrr will be kept)
            DO rrr = 1, MPIvar%glob_size
               IF (((I + 1) .GE. istart_vec(rrr)) .AND. ((I + 1) .LE. iend_vec(rrr))) THEN
                  EXIT
               END IF
            END DO

            ! put in or out diagonal (always on diagonal for serial)
            IF ((J .GE. istart_vec(rrr)) .AND. (J .LE. iend_vec(rrr) - 1)) THEN
               counter_diag(I + 1) = counter_diag(I + 1) + 1
            ELSE
               counter_offdiag(I + 1) = counter_offdiag(I + 1) + 1
            END IF
         END DO
      END DO

#ifdef PARALL
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, counter_diag, N, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, counter_offdiag, N, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MatMPIAIJSetPreallocation(matPETSC%matK,PETSC_DECIDE,counter_diag(istart+1:iend),PETSC_DECIDE,counter_offdiag(istart+1:iend), ierr); CHKERRA(ierr)
#else
      CALL MatSeqAIJSetPreallocation(matPETSC%matK, PETSC_DECIDE, counter_diag, ierr); CHKERRA(ierr)
#endif

      DEALLOCATE (counter_diag)
      DEALLOCATE (counter_offdiag)
      DEALLOCATE (istart_vec)
      DEALLOCATE (iend_vec)

   END SUBROUTINE preallocate_matrix_PETSC

   !***********************************************
   ! From the generic CSR storage, fills the PETSC
   ! instance of the matrix
   !***********************************************
   SUBROUTINE build_mat_PETSC(matPETSC)
      TYPE(PETSC_STRUC), INTENT(IN)         :: matPETSC
      PetscErrorCode            :: ierr
      PetscInt                  :: one = 1, r, rr, row_nnz
      PetscInt, ALLOCATABLE     :: loc2glob_petsc(:)
      PetscInt, ALLOCATABLE     :: cols_petsc(:)
      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke

      cks = 0
      cke = 0

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tps)
         CALL system_CLOCK(cks, clock_rate)
      END IF

      ALLOCATE (loc2glob_petsc(matPETSC%n))
      loc2glob_petsc = matPETSC%loc2glob - 1

      ALLOCATE (cols_petsc(SIZE(matPETSC%cols)))
      cols_petsc = matPETSC%cols - 1

      ! Set values one row at a time to reduce PETSc call overhead (this is faster)
      DO r = 1, matPETSC%n
         row_nnz = matPETSC%rowptr(r + 1) - matPETSC%rowptr(r)
         CALL MatSetValues(matPETSC%matK, one, loc2glob_petsc(r), row_nnz, cols_petsc(matPETSC%rowptr(r):matPETSC%rowptr(r + 1) - 1), matPETSC%vals_matK(matPETSC%rowptr(r):matPETSC%rowptr(r + 1) - 1), INSERT_VALUES, ierr); CHKERRA(ierr)
      END DO

      CALL MatAssemblyBegin(matPETSC%matK, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
      CALL MatAssemblyEnd(matPETSC%matK, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)

      DEALLOCATE (loc2glob_petsc)
      DEALLOCATE (cols_petsc)

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tpe)
         CALL system_CLOCK(cke, clock_rate)
         timing%rlstime2 = timing%rlstime2 + (cke - cks)/REAL(clock_rate)
         timing%clstime2 = timing%clstime2 + tpe - tps
      END IF

   END SUBROUTINE build_mat_PETSC

   SUBROUTINE setup_permuted_system(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode                   :: ierr
      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tps)
         CALL system_CLOCK(cks, clock_rate)
      END IF

      CALL hdg_compute_permutation_mpiadj(matPETSC%matK, lssolver%theta, matPETSC%is_num, matPETSC%is_perm, ierr)

      CALL MatCreateSubMatrix(matPETSC%matK, matPETSC%is_perm, matPETSC%is_perm, MAT_INITIAL_MATRIX, matPETSC%matKperm, ierr); CHKERRA(ierr)
      CALL MatSetBlockSize(matPETSC%matKperm, phys%neq, ierr); CHKERRA(ierr)
      CALL MatCreateVecs(matPETSC%matKperm, matPETSC%solPETSC_vec_perm, matPETSC%rhs_vec_perm, ierr); CHKERRA(ierr)
      CALL KSPSetOperators(matPETSC%ksp, matPETSC%matKperm, matPETSC%matKperm, ierr); CHKERRA(ierr)

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tpe)
         CALL system_CLOCK(cke, clock_rate)
         timing%rlstime5 = timing%rlstime5 + (cke - cks)/REAL(clock_rate)
         timing%clstime5 = timing%clstime5 + tpe - tps
      END IF

   END SUBROUTINE setup_permuted_system

   SUBROUTINE setup_ksp_pc_settings(matPETSC, ir)
#include "mkl_service.fi"

      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      INTEGER, OPTIONAL, INTENT(IN)               :: ir
      KSPType                          :: ksp_type
      PCType                           :: pc_type
      PetscViewerAndFormat             :: vf
      PetscErrorCode                   :: ierr
      CHARACTER(LEN=128)        :: prefix

      ! Set solver and preconditioner type from param.txt
      CALL set_ksp_pc_type(ksp_type, pc_type)
      CALL KSPSetType(matPETSC%ksp, ksp_type, ierr); CHKERRA(ierr)

      ! Set the norm type
      IF (lssolver%kspnorm .EQ. 1) THEN
         CALL KSPSetNormType(matPETSC%ksp, KSP_NORM_PRECONDITIONED, ierr); CHKERRA(ierr)
      ELSEIF (lssolver%kspnorm .EQ. 2) THEN
         CALL KSPSetNormType(matPETSC%ksp, KSP_NORM_UNPRECONDITIONED, ierr); CHKERRA(ierr)
      END IF

      ! Set tolerances on convergence
      CALL KSPSetTolerances(matPETSC%ksp, lssolver%rtol, lssolver%atol, PETSC_DEFAULT_REAL, lssolver%kspitmax, ierr); CHKERRA(ierr)

      ! Get the preconditioner and set it
      CALL KSPGetPC(matPETSC%ksp, matPETSC%pc, ierr); CHKERRA(ierr)
      CALL PCSetType(matPETSC%pc, pc_type, ierr); CHKERRA(ierr)

      IF (lssolver%pctype .EQ. 'PCLU') THEN
         !call PCFactorSetMatSolverType(matPETSC%pc, MATSOLVERMUMPS, ierr)
         CALL PCFactorSetMatSolverType(matPETSC%pc,  MATSOLVERMKL_CPARDISO, ierr)
      END IF

#ifdef TEMPERATURE
      IF (lssolver%pctype .EQ. 'PCFIELDSPLIT') THEN
         IF(lssolver%reorder_matrix) THEN
            CALL setup_permuted_system(matPETSC)
         END IF
         CALL setup_fieldsplit_ngamma_energies(matPETSC)
         !call setup_fieldsplit_plasma_neutral(matPETSC)
      END IF
#endif
      SELECT CASE (lssolver%kspmethd)
       CASE ('KSPGMRES')
         IF (lssolver%pctype .EQ. 'PCMG') THEN
            ! For MG, use flexible GMRES to remain robust with nested/nonlinear preconditioning effects.
            CALL KSPSetType(matPETSC%ksp, KSPFGMRES, ierr); CHKERRA(ierr)
         END IF
         ! GMRES + any other preconditioners
         CALL KSPGMRESSetRestart(matPETSC%ksp, lssolver%gmresres, ierr); CHKERRA(ierr)
      END SELECT

      ! Set the settings for the preconditioner, depending on the type of preconditioner selected in param.txt
      SELECT CASE (lssolver%pctype)
       CASE ('PCMG')
         IF(lssolver%reorder_matrix) THEN
            CALL setup_permuted_system(matPETSC)
            matPETSC%is_row_f = matPETSC%is_perm     ! graph ordering, computed once
            matPETSC%is_col_f = matPETSC%is_perm
            CALL set_MG_ksp_pc_settings_reordering(matPETSC, ir)
         ELSE
            CALL set_MG_ksp_pc_settings(matPETSC, ir)
         END IF

       CASE ('PCASM')
         IF(lssolver%reorder_matrix) THEN
            CALL setup_permuted_system(matPETSC)
         END IF
         CALL set_ASM_ksp_pc_settings(matPETSC)
       CASE ('PCAMG')
         CALL set_AMG_ksp_pc_settings(matPETSC)
       CASE ('PCBJACOBI')
         CALL set_BJ_ksp_pc_settings(matPETSC)
       CASE ('PCKSP')
         CALL set_PCKSP_ksp_pc_settings(matPETSC)
      END SELECT

      CALL KSPSetUp(matPETSC%ksp, ierr); CHKERRA(ierr)

      SELECT CASE (lssolver%pctype)
       CASE ('PCMG')
         CALL finalize_MG_subsolvers(matPETSC, ir)
       CASE ('PCASM')
         CALL finalize_ASM_subsolvers(matPETSC)
       CASE ('PCBJACOBI')
         CALL finalize_BJ_subsolvers(matPETSC)
      END SELECT


      ! Print out the residue
      IF (lssolver%kspitrace .EQ. 1) THEN
         CALL PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, vf, ierr); CHKERRA(ierr)
         CALL KSPMonitorSet(matPETSC%ksp, KSPMonitorTrueResidual, vf, PetscViewerAndFormatDestroy, ierr); CHKERRA(ierr)
      END IF

   END SUBROUTINE setup_ksp_pc_settings

   SUBROUTINE fill_rhs_PETSC(matPETSC)
      TYPE(PETSC_STRUC), INTENT(IN)         :: matPETSC
      PetscErrorCode            :: ierr
      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke
      PetscInt, ALLOCATABLE     :: loc2glob_petsc(:)

      cks = 0
      cke = 0

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tps)
         CALL system_CLOCK(cks, clock_rate)
      END IF

      ALLOCATE (loc2glob_petsc(matPETSC%n))
      loc2glob_petsc = matPETSC%loc2glob - 1

      ! Set values in the vector rhs (-1 because PETSc uses 0-based indexing)
      CALL VecSetValues(matPETSC%rhs_vec, matPETSC%n, loc2glob_petsc, rhs%vals, INSERT_VALUES, ierr); CHKERRA(ierr)
      CALL VecAssemblyBegin(matPETSC%rhs_vec, ierr); CHKERRA(ierr)
      CALL VecAssemblyEnd(matPETSC%rhs_vec, ierr); CHKERRA(ierr)

      DEALLOCATE (loc2glob_petsc)

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tpe)
         CALL system_CLOCK(cke, clock_rate)
         timing%rlstime3 = timing%rlstime3 + (cke - cks)/REAL(clock_rate)
         timing%clstime3 = timing%clstime3 + tpe - tps
      END IF
   END SUBROUTINE fill_rhs_PETSC

   !***********************************************
   ! Solve problem with PETSC
   !***********************************************
   SUBROUTINE solve_mat_PETSC(matPETSC, ir)
      TYPE(PETSC_STRUC), INTENT(INOUT)         :: matPETSC
      INTEGER, INTENT(IN)       :: ir
      REAL*8                    :: tps, tpe, tps_total, tpe_total
      INTEGER                   :: cks, clock_rate, cke, cks_total, cke_total, i
      PetscErrorCode            :: ierr
      LOGICAL                   :: reorder, scaling, reused
      PetscBool                 :: reuse_pc
      Vec                       :: c, brcm
      PetscViewer               :: viewer_ksp, viewer_pc
      TYPE(tpetscviewer)             :: viewer_draw
      PetscDraw                      :: draw
      IS                        ::isrow, iscol
      VecScatter                :: scat
      IS                          :: is_num, is_perm
      Vec                         :: bperm, xperm
      VecScatter                  :: sc
      KSP                         :: ksp

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tps_total)
         CALL system_CLOCK(cks_total, clock_rate)
      END IF

      ! Set or not the initial guess to zero
      IF ((.NOT. lssolver%igz) .AND. (lssolver%kspmethd .NE. 'KSPPREONLY')) THEN
         CALL KSPSetInitialGuessNonzero(matPETSC%ksp, PETSC_TRUE, ierr); CHKERRA(ierr)
      END IF

      ! Reuse or not preconditioner
      CALL set_preconditioner_reuse(matPETSC%ksp, ir, verbose=.TRUE.)

      IF (lssolver%pctype .EQ. 'PCASM') THEN
         CALL set_ASM_preconditioner_reuse(matPETSC, ir)
      END IF

      IF (lssolver%pctype .EQ. 'PCMG') THEN
         CALL set_MG_preconditioner_reuse(matPETSC, ir)
      END IF

      IF(lssolver%pctype .EQ. 'PCASM' .OR. lssolver%pctype .EQ. 'PCMG' .OR. lssolver%pctype .EQ. 'PCFIELDSPLIT') THEN

         IF (lssolver%reorder_matrix) THEN

            IF (lssolver%timing) THEN
               CALL cpu_TIME(tps)
               CALL system_CLOCK(cks, clock_rate)
            END IF

            ! matKperm is created in setup_permuted_system on the first iteration
            IF (.NOT. matK%start) THEN
               CALL MatCreateSubMatrix(matPetsc%matk, matPetsc%is_perm, matPetsc%is_perm, MAT_REUSE_MATRIX, matPETSC%matKperm, ierr)

               ! Level operators follow the outer preconditioner reuse policy (lssolver%rprecond)
               IF (lssolver%pctype .EQ. 'PCMG') THEN
                  CALL PCGetReusePreconditioner(matPETSC%pc, reuse_pc, ierr); CHKERRA(ierr)
                  IF (.NOT. reuse_pc) CALL refresh_MG_level_operators(matPETSC)
               END IF
            END IF

            CALL hdg_create_vec_scatter(matPetsc%is_num, matPetsc%rhs_vec, matPETSC%rhs_vec_perm, sc, ierr)

            CALL VecScatterBegin(sc, matPetsc%rhs_vec, matPETSC%rhs_vec_perm, INSERT_VALUES, SCATTER_FORWARD, ierr)
            CALL VecScatterEnd  (sc, matPetsc%rhs_vec, matPETSC%rhs_vec_perm, INSERT_VALUES, SCATTER_FORWARD, ierr)

            !only if you use a nonzero initial guess (-ksp_initial_guess_nonzero / KSPSetInitialGuessNonzero):
            CALL VecScatterBegin(sc, matPetsc%solPETSC_vec, matPETSC%solPETSC_vec_perm, INSERT_VALUES, SCATTER_FORWARD, ierr)
            CALL VecScatterEnd  (sc, matPetsc%solPETSC_vec, matPETSC%solPETSC_vec_perm, INSERT_VALUES, SCATTER_FORWARD, ierr)

            IF (lssolver%timing) THEN
               CALL cpu_TIME(tpe)
               CALL system_CLOCK(cke, clock_rate)
               timing%rlstime5 = timing%rlstime5 + (cke - cks)/REAL(clock_rate)
               timing%clstime5 = timing%clstime5 + tpe - tps
            END IF

            CALL KSPSolve(matPETSC%ksp, matPETSC%rhs_vec_perm, matPETSC%solPETSC_vec_perm, ierr)

            IF (lssolver%timing) THEN
               CALL cpu_TIME(tps)
               CALL system_CLOCK(cks, clock_rate)
            END IF

            CALL VecScatterBegin(sc, matPETSC%solPETSC_vec_perm, matPetsc%solPETSC_vec, INSERT_VALUES, SCATTER_REVERSE, ierr)
            CALL VecScatterEnd  (sc, matPETSC%solPETSC_vec_perm, matPetsc%solPETSC_vec, INSERT_VALUES, SCATTER_REVERSE, ierr)

            IF (lssolver%timing) THEN
               CALL cpu_TIME(tpe)
               CALL system_CLOCK(cke, clock_rate)
               timing%rlstime5 = timing%rlstime5 + (cke - cks)/REAL(clock_rate)
               timing%clstime5 = timing%clstime5 + tpe - tps
            END IF

            CALL VecScatterDestroy(sc, ierr); CHKERRA(ierr)
         ELSE
            CALL KSPSolve(matPETSC%ksp, matPetsc%rhs_vec, matPetsc%solPETSC_vec, ierr)
         END IF
      ELSE
         CALL KSPSolve(matPETSC%ksp, matPetsc%rhs_vec, matPetsc%solPETSC_vec, ierr)
      END IF

      ! call PetscViewerDrawOpen(PETSC_COMM_WORLD, PETSC_NULL_CHARACTER, 'matK pattern', &
      !    PETSC_DECIDE, PETSC_DECIDE, 800, 800, viewer_draw, ierr); CHKERRA(ierr)
      ! call PetscViewerDrawGetDraw(viewer_draw, 0, draw, ierr); CHKERRA(ierr)
      ! call PetscDrawSetPause(draw, 1.0d0, ierr); CHKERRA(ierr) ! -1 = wait for mouse click before closing
      ! call MatView(matPETSC%matKperm, viewer_draw, ierr); CHKERRA(ierr)
      ! call PetscViewerDestroy(viewer_draw, ierr); CHKERRA(ierr)

      ! call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'ksp_view.log', viewer_ksp, ierr); CHKERRA(ierr)
      ! call PetscViewerFileSetMode(viewer_ksp, FILE_MODE_APPEND, ierr); CHKERRA(ierr)
      ! call KSPView(matPETSC%ksp, viewer_ksp, ierr); CHKERRA(ierr)
      ! call PetscViewerDestroy(viewer_ksp, ierr); CHKERRA(ierr)

      ! call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'pc_view.log', viewer_pc, ierr); CHKERRA(ierr)
      ! call PetscViewerFileSetMode(viewer_pc, FILE_MODE_APPEND, ierr); CHKERRA(ierr)
      ! call PCView(matPETSC%pc, viewer_pc, ierr); CHKERRA(ierr)
      ! call PetscViewerDestroy(viewer_pc, ierr); CHKERRA(ierr)


      IF (scaling) THEN
         ! Rescale the solution back to the original system
         CALL VecPointwiseMult(matPETSC%solPETSC_vec, matPETSC%solPETSC_vec, c, ierr); CHKERRA(ierr)
      END IF

      ! Get the number of iterations and the reason for convergence
      CALL reused_preconditioner(ir, reused)

      CALL KSPGetIterationNumber(matPETSC%ksp, matPETSC%its, ierr); CHKERRA(ierr)
      CALL KSPGetConvergedReason(matPETSC%ksp, matPETSC%convergedReason, ierr); CHKERRA(ierr)

      IF (.NOT. reused .OR. ir .EQ. 1) THEN
         matPETSC%its_recompute_prec = matPETSC%its
      END IF

      ! Compute the unpreconditioned relative residue
      CALL compute_residue

      IF (lssolver%timing) THEN
         CALL cpu_TIME(tpe_total)
         CALL system_CLOCK(cke_total, clock_rate)
         timing%rlstime4 = timing%rlstime4 + (cke_total - cks_total)/REAL(clock_rate)
         timing%clstime4 = timing%clstime4 + tpe_total - tps_total
      END IF

   END SUBROUTINE solve_mat_PETSC

   SUBROUTINE set_preconditioner_reuse(ksp, ir, verbose)
      KSP, INTENT(INOUT)              :: ksp
      INTEGER, INTENT(IN)             :: ir
      LOGICAL, INTENT(IN)             :: verbose
      PetscErrorCode                   :: ierr

      IF (lssolver%pctype .NE. 'PCLU') THEN
         SELECT CASE (lssolver%rprecond)
          CASE (1)
            CALL KSPSetReusePreconditioner(ksp, PETSC_FALSE, ierr); CHKERRA(ierr)
          CASE (2)
            IF (MOD(ir - 1, lssolver%Nrprecond) .EQ. 0) THEN
               CALL KSPSetReusePreconditioner(ksp, PETSC_FALSE, ierr); CHKERRA(ierr)
               IF (verbose) THEN
                  IF (MPIvar%glob_id .EQ. 0) THEN
                     IF(ir .EQ. 1) THEN
                        WRITE(*,*) '----- PETSC IS COMPUTING PRECONDITIONER -----'
                     ELSE
                        WRITE (*, *) '----- PETSC IS RECOMPUTING PRECONDITIONER -----'
                     END IF
                  END IF
               END IF
            ELSE
               CALL KSPSetReusePreconditioner(ksp, PETSC_TRUE, ierr); CHKERRA(ierr)
            END IF
          CASE (3)
            IF (matPETSC%its .EQ. lssolver%kspitmax) THEN
               CALL KSPSetReusePreconditioner(ksp, PETSC_FALSE, ierr); CHKERRA(ierr)
               IF (verbose) THEN
                  IF (MPIvar%glob_id .EQ. 0) THEN
                     WRITE (*, *) '----- PETSC IS RECOMPUTING PRECONDITIONER BECAUSE MAX NUMBER OF ITERATIONS WAS REACHED -----'
                  END IF
               END IF
            ELSE
               CALL KSPSetReusePreconditioner(ksp, PETSC_TRUE, ierr); CHKERRA(ierr)
            END IF
          CASE (4)
            IF(ir .EQ. 1) THEN
               IF (verbose) THEN
                  IF (MPIvar%glob_id .EQ. 0) THEN
                     WRITE(*,*) '----- PETSC IS COMPUTING PRECONDITIONER -----'

                  END IF

               END IF
               CALL KSPSetReusePreconditioner(ksp, PETSC_FALSE, ierr); CHKERRA(ierr)
            ELSE
               IF ((REAL(matPETSC%its)/REAL(matPETSC%its_recompute_prec) .GE. 1.25) .OR. matPETSC%its .EQ. lssolver%kspitmax) THEN
                  CALL KSPSetReusePreconditioner(ksp, PETSC_FALSE, ierr); CHKERRA(ierr)
                  IF (verbose) THEN
                     IF (MPIvar%glob_id .EQ. 0) THEN
                        IF(REAL(matPETSC%its)/REAL(matPETSC%its_recompute_prec) .GE. 1.25) THEN
                           WRITE (*, *) '----- PETSC IS RECOMPUTING PRECONDITIONER BECAUSE THE NUMBER OF ITERATIONS INCREASED BY 25% COMPARED TO THE LAST TIME THE PRECONDITIONER WAS RECOMPUTED -----'
                        ELSE
                           WRITE (*, *) '----- PETSC IS RECOMPUTING PRECONDITIONER BECAUSE MAX NUMBER OF ITERATIONS WAS REACHED -----'
                        END IF
                     END IF
                  END IF
               ELSE
                  IF (verbose) THEN
                     IF (MPIvar%glob_id .EQ. 0) THEN
                        WRITE(*,*) '----- PETSC IS REUSING PRECONDITIONER -----'
                     END IF
                  END IF
                  CALL KSPSetReusePreconditioner(ksp, PETSC_TRUE, ierr); CHKERRA(ierr)
               END IF
            END IF
          CASE DEFAULT
            IF (MPIvar%glob_id .EQ. 0) THEN
               WRITE (*, *) 'Wrong choice, value must be 1, 2, 3 or 4, stopping. '
               CALL MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
            END IF
         END SELECT
      ELSEIF (lssolver%pctype .EQ. 'PCMG') THEN
         CALL KSPSetReusePreconditioner(ksp, PETSC_TRUE, ierr); CHKERRA(ierr)
      END IF

   END SUBROUTINE set_preconditioner_reuse

   SUBROUTINE reused_preconditioner(ir, reused)
      INTEGER, INTENT(IN)             :: ir
      LOGICAL, INTENT(OUT)             :: reused
      PetscErrorCode                   :: ierr

      reused = .FALSE.

      IF (lssolver%pctype .NE. 'PCLU') THEN
         SELECT CASE (lssolver%rprecond)
          CASE (2)
            IF (MOD(ir - 1, lssolver%Nrprecond) .EQ. 0) THEN
               reused = .FALSE.
            ELSE
               reused = .TRUE.
            END IF
          CASE (3)
            IF (matPETSC%its .EQ. lssolver%kspitmax) THEN
               reused = .FALSE.
            ELSE
               reused = .TRUE.
            END IF
          CASE (4)
            IF (REAL(matPETSC%its)/REAL(matPETSC%its_recompute_prec) .GE. 1.25) THEN
               reused = .FALSE.
            ELSE
               reused = .TRUE.
            END IF
          CASE DEFAULT
            IF (MPIvar%glob_id .EQ. 0) THEN
               WRITE (*, *) 'Wrong choice, value must be 1, 2, 3 or 4, stopping. '
               CALL MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
            END IF
         END SELECT
      END IF

   END SUBROUTINE reused_preconditioner

   SUBROUTINE set_MG_preconditioner_reuse(matPETSC, ir)
      ! Propagate the preconditioner reuse flag to all MG smoother and coarse sub-KSPs.
      ! Must be called every Newton iteration after the outer set_preconditioner_reuse.
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      INTEGER, INTENT(IN)               :: ir
      PetscErrorCode                   :: ierr
      KSP                              :: ksp_smooth, ksp_coarse
      INTEGER                          :: i

      CALL PCMGGetCoarseSolve(matPETSC%pc, ksp_coarse, ierr); CHKERRA(ierr)
      CALL set_preconditioner_reuse(ksp_coarse, ir, verbose=.FALSE.)

      DO i = 1, lssolver%mglevels - 1
         CALL PCMGGetSmoother(matPETSC%pc, i, ksp_smooth, ierr); CHKERRA(ierr)
         CALL set_preconditioner_reuse(ksp_smooth, ir, verbose=.FALSE.)
      END DO

   END SUBROUTINE set_MG_preconditioner_reuse

   SUBROUTINE set_ASM_preconditioner_reuse(matPETSC, ir)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      INTEGER, INTENT(IN)              :: ir
      PetscErrorCode                   :: ierr
      PetscInt                         :: n_local, first_local
      KSP, POINTER                     :: subksps(:)
      INTEGER                          :: i

      NULLIFY(subksps)
      CALL PCASMGetSubKSP(matPETSC%pc, n_local, first_local, PETSC_NULL_KSP, ierr); CHKERRA(ierr)

      IF (n_local .GT. 0) THEN
         ALLOCATE(subksps(n_local))
         CALL PCASMGetSubKSP(matPETSC%pc, n_local, first_local, subksps, ierr); CHKERRA(ierr)

         DO i = 1, n_local
            CALL set_preconditioner_reuse(subksps(i), ir, verbose=.FALSE.)
         END DO

         DEALLOCATE(subksps)
         NULLIFY(subksps)
      END IF
   END SUBROUTINE set_ASM_preconditioner_reuse

   SUBROUTINE set_PCKSP_ksp_pc_settings(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode                   :: ierr
      KSP                              :: innerksp
      PC                               :: innerpc

      CALL PCKSPGetKSP(matPETSC%pc, innerksp, ierr); CHKERRA(ierr)
      CALL KSPSetType(innerksp, KSPCHEBYSHEV, ierr); CHKERRA(ierr)

      ! Polynomial degree = number of Chebyshev iterations
      CALL KSPSetTolerances(innerksp, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, 15, ierr); CHKERRA(ierr)

      ! --- Inner PC (cheap smoother) ---
      CALL KSPGetPC(innerksp, innerpc, ierr); CHKERRA(ierr)
      CALL PCSetType(innerpc, PCJACOBI, ierr); CHKERRA(ierr)

   END SUBROUTINE set_PCKSP_ksp_pc_settings

   SUBROUTINE set_AMG_ksp_pc_settings(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode                   :: ierr
      ! Tell BoomerAMG the unknowns are interlaced blocks of coupled physical variables (systems AMG).
      CALL MatSetBlockSize(matPETSC%matK, phys%Neq, ierr); CHKERRA(ierr)
      CALL PCHYPRESetType(matPETSC%pc, 'boomeramg', ierr); CHKERRA(ierr)
      CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_relax_type_all', 'symmetric-SOR/Jacobi', ierr); CHKERRA(ierr)
      CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_coarsen_type', 'HMIS', ierr); CHKERRA(ierr)
      CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_interp_type', 'ext+i', ierr); CHKERRA(ierr)
      CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_max_levels', '25', ierr); CHKERRA(ierr)
      CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_num_sweeps', '2', ierr); CHKERRA(ierr)
      CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_agg_nl', '1', ierr); CHKERRA(ierr)
      ! Anisotropic diffusion needs a higher threshold so coarsening follows the dominant (parallel) coupling.
      CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_strong_threshold', '0.5', ierr); CHKERRA(ierr)
   END SUBROUTINE set_AMG_ksp_pc_settings


#ifdef TEMPERATURE
   SUBROUTINE setup_fieldsplit_ngamma_energies(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode                   :: ierr
      PetscInt                         :: istart, iend, ig
      PetscInt                         :: n_local, n1, n2, nsplit
      PetscInt, ALLOCATABLE            :: idx1(:), idx2(:)
      IS                               :: is1, is2
      PC                               :: subpc, subsubpc
      KSP, POINTER                     :: subksps(:), subsubksps(:)
      PetscViewerAndFormat             :: vf
      INTEGER                          :: i

      ! Build index sets for interleaved unknown ordering: [f0, f1, f2, ...]
      CALL MatGetOwnershipRange(matPETSC%matK, istart, iend, ierr); CHKERRA(ierr)

      n_local = iend - istart
      nsplit = 2

      IF (n_local .GT. 0) THEN
         ALLOCATE (idx1(n_local), idx2(n_local))
      ELSE
         ALLOCATE (idx1(0), idx2(0))
      END IF

      idx1 = 0
      idx2 = 0
      n1 = 0
      n2 = 0

      DO ig = istart, iend - 1
         SELECT CASE (MOD(ig, 3))
          CASE (0, 1) ! density or momentum
            n1 = n1 + 1
            idx1(n1) = ig
          CASE (2, 3)    ! energies
            n2 = n2 + 1
            idx2(n2) = ig
         END SELECT
      END DO

      CALL ISCreateGeneral(PETSC_COMM_WORLD, n1, idx1, PETSC_COPY_VALUES, is1, ierr); CHKERRA(ierr)
      CALL ISCreateGeneral(PETSC_COMM_WORLD, n2, idx2, PETSC_COPY_VALUES, is2, ierr); CHKERRA(ierr)

      CALL PCFieldSplitSetIS(matPETSC%pc, 'plasma', is1, ierr); CHKERRA(ierr)
      CALL PCFieldSplitSetIS(matPETSC%pc, 'energies', is2, ierr); CHKERRA(ierr)

      CALL ISDestroy(is1, ierr); CHKERRA(ierr)
      CALL ISDestroy(is2, ierr); CHKERRA(ierr)
      DEALLOCATE (idx1, idx2)

      CALL PCFieldSplitSetType(matPETSC%pc, PC_COMPOSITE_MULTIPLICATIVE, ierr); CHKERRA(ierr)
      CALL PCSetUp(matPETSC%pc, ierr); CHKERRA(ierr)
      CALL KSPSetUp(matPETSC%ksp, ierr); CHKERRA(ierr)

      ALLOCATE (subksps(nsplit))

      CALL PCFieldSplitGetSubKSP(matPETSC%pc, nsplit, subksps, ierr); CHKERRA(ierr)

      ! Plasma split
      CALL KSPSetType(subksps(1), KSPPREONLY, ierr); CHKERRA(ierr)
      CALL KSPGetPC(subksps(1), subpc, ierr); CHKERRA(ierr)
      CALL PCSetType(subpc, PCLU, ierr); CHKERRA(ierr)
      CALL PCFactorSetMatSolverType(subpc,  MATSOLVERMKL_CPARDISO, ierr); CHKERRA(ierr)

      ! Energy split
      CALL KSPSetType(subksps(2), KSPPREONLY, ierr); CHKERRA(ierr)
      CALL KSPGetPC(subksps(2), subpc, ierr); CHKERRA(ierr)
      CALL PCSetType(subpc, PCLU, ierr); CHKERRA(ierr)
      CALL PCFactorSetMatSolverType(subpc,  MATSOLVERMKL_CPARDISO, ierr); CHKERRA(ierr)

      NULLIFY (subksps)

   END SUBROUTINE setup_fieldsplit_ngamma_energies
#endif

   SUBROUTINE set_MG_ksp_pc_settings(matPETSC, ir)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      INTEGER, INTENT(IN)               :: ir
      PetscErrorCode                   :: ierr
      KSP                              :: ksp_smooth, ksp_coarse
      PC                               :: innerpc, pc_coarse
      INTEGER                          :: i
      ! Set the number of levels
      CALL PCMGSetLevels(matPETSC%pc, lssolver%mglevels, PETSC_NULL_MPI_COMM, ierr); CHKERRA(ierr)

      ! Set the form of multigrid cycle
      SELECT CASE (lssolver%mgtypeform)
       CASE (1)
         CALL PCMGSetType(matPETSC%pc, PC_MG_MULTIPLICATIVE, ierr); CHKERRA(ierr)
       CASE (2)
         CALL PCMGSetType(matPETSC%pc, PC_MG_ADDITIVE, ierr); CHKERRA(ierr)
       CASE (3)
         CALL PCMGSetType(matPETSC%pc, PC_MG_FULL, ierr); CHKERRA(ierr)
       CASE (4)
         CALL PCMGSetType(matPETSC%pc, PC_MG_KASKADE, ierr); CHKERRA(ierr)
      END SELECT

      ! Set the type of cycle (V or W)
      SELECT CASE (lssolver%mgcycletype)
       CASE (1)
         CALL PCMGSetCycleType(matPETSC%pc, PC_MG_CYCLE_V, ierr); CHKERRA(ierr)
       CASE (2)
         CALL PCMGSetCycleType(matPETSC%pc, PC_MG_CYCLE_W, ierr); CHKERRA(ierr)
      END SELECT

      ! Set the type of Galerkin coarse grid operator construction (A_coarse = R*A_fine*P)
      CALL PCMGSetGalerkin(matPETSC%pc, PC_MG_GALERKIN_BOTH, ierr); CHKERRA(ierr)

      ! Coarse solve: redundant + ILU(1) offers better global coupling than ASM+ILU(0).
      !call PCMGGetCoarseSolve(matPETSC%pc, ksp_coarse, ierr); CHKERRA(ierr)
      !call KSPSetType(ksp_coarse, KSPPREONLY, ierr); CHKERRA(ierr)
      !call KSPGetPC(ksp_coarse, pc_coarse, ierr); CHKERRA(ierr)
      !call PCSetType(pc_coarse, PCREDUNDANT, ierr); CHKERRA(ierr)

      ! Coarse solve
      CALL PCMGGetCoarseSolve(matPETSC%pc, ksp_coarse, ierr); CHKERRA(ierr)
      CALL KSPSetType(ksp_coarse, KSPPREONLY, ierr); CHKERRA(ierr)
      CALL KSPGetPC(ksp_coarse, pc_coarse, ierr); CHKERRA(ierr)
      CALL PCSetType(pc_coarse, PCLU, ierr); CHKERRA(ierr)

      CALL PCFactorSetMatSolverType(pc_coarse,  MATSOLVERMKL_CPARDISO, ierr); CHKERRA(ierr)
      !call PCFactorSetMatSolverType(pc_coarse, MATSOLVERMUMPS, ierr); CHKERRA(ierr)

      CALL set_preconditioner_reuse(ksp_coarse, ir, verbose=.FALSE.)

      ! Create the multigrid operators (restriction and prolongation)
      DO i = 1, lssolver%mglevels - 1

         ! Set coarse operator and transfer operators for this level
         CALL PCMGSetInterpolation(matPETSC%pc, i, matPETSC%matP(i), ierr); CHKERRA(ierr) ! prolongation PL -> PH
         CALL PCMGSetRestriction(matPETSC%pc, i, matPETSC%matR(i), ierr); CHKERRA(ierr)   ! restriction PH -> PL

         ! Set smoother on all levels
         CALL PCMGGetSmoother(matPETSC%pc, i, ksp_smooth, ierr); CHKERRA(ierr)
         IF (i .EQ. lssolver%mglevels - 1) THEN
            ! Explicitly attach finest-level operators to avoid PCMG fallback behavior.
            CALL KSPSetOperators(ksp_smooth, matPETSC%matK, matPETSC%matK, ierr); CHKERRA(ierr)
         END IF
         ! use plain GMRES (not flexible) since the inner PC (ASM+ILU) never changes
         ! between smoothing steps; FGMRES only adds overhead here.
         CALL KSPSetType(ksp_smooth, KSPGMRES, ierr); CHKERRA(ierr)
         !call KSPGMRESSetRestart(ksp_smooth, 5, ierr); CHKERRA(ierr)

         ! Critical for Krylov smoothers in PCMG: without this the smoother
         ! discards the current MG iterate and restarts from zero every sweep.
         CALL KSPSetInitialGuessNonzero(ksp_smooth, PETSC_TRUE, ierr); CHKERRA(ierr)

         ! Configure the smoother's inner preconditioner
         CALL KSPGetPC(ksp_smooth, innerpc, ierr); CHKERRA(ierr)
         CALL PCSetType(innerpc, PCASM, ierr); CHKERRA(ierr)
         ! Restricted ASM avoids double-counting corrections in the overlap region, sharpening convergence per sweep.
         CALL PCASMSetType(innerpc, PC_ASM_RESTRICT, ierr); CHKERRA(ierr)
         CALL PCASMSetOverlap(innerpc, 3, ierr); CHKERRA(ierr)
         CALL PCASMSetLocalSubdomains(innerpc, 1, PETSC_NULL_IS, PETSC_NULL_IS, ierr); CHKERRA(ierr)

         ! Set the norm type for the smoother, which can be different from the one used for the outer solver
         IF (lssolver%kspnorm .EQ. 1) THEN
            CALL KSPSetNormType(ksp_smooth, KSP_NORM_PRECONDITIONED, ierr); CHKERRA(ierr)
         ELSEIF (lssolver%kspnorm .EQ. 2) THEN
            CALL KSPSetNormType(ksp_smooth, KSP_NORM_UNPRECONDITIONED, ierr); CHKERRA(ierr)
         END IF

         ! A smoother should apply a few cheap relaxation sweeps, not solve accurately.
         ! The finest level's KSP is bound to the full matK (see KSPSetOperators above), so it
         ! needs the most sweeps. Interior levels (e.g. order-2 with orders=[1,2,4]) still need
         ! enough sweeps to damp their own oscillatory modes before restricting further down.
         IF (i .EQ. lssolver%mglevels - 1) THEN
            CALL KSPSetTolerances(ksp_smooth, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, 5, ierr); CHKERRA(ierr)
         ELSE
            CALL KSPSetTolerances(ksp_smooth, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, 5, ierr); CHKERRA(ierr)
         END IF

      END DO

   END SUBROUTINE set_MG_ksp_pc_settings

   SUBROUTINE set_MG_ksp_pc_settings_reordering(matPETSC, ir)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      INTEGER, INTENT(IN)               :: ir
      PetscErrorCode                   :: ierr
      KSP                              :: ksp_smooth, ksp_coarse
      PC                               :: innerpc, pc_coarse
      INTEGER                          :: i
      LOGICAL                          :: perm

      perm = lssolver%reorder_matrix

      CALL PCMGSetLevels(matPETSC%pc, lssolver%mglevels, PETSC_NULL_MPI_COMM, ierr); CHKERRA(ierr)

      SELECT CASE (lssolver%mgtypeform)
       CASE (1); CALL PCMGSetType(matPETSC%pc, PC_MG_MULTIPLICATIVE, ierr); CHKERRA(ierr)
       CASE (2); CALL PCMGSetType(matPETSC%pc, PC_MG_ADDITIVE, ierr); CHKERRA(ierr)
       CASE (3); CALL PCMGSetType(matPETSC%pc, PC_MG_FULL, ierr); CHKERRA(ierr)
       CASE (4); CALL PCMGSetType(matPETSC%pc, PC_MG_KASKADE, ierr); CHKERRA(ierr)
      END SELECT

      SELECT CASE (lssolver%mgcycletype)
       CASE (1); CALL PCMGSetCycleType(matPETSC%pc, PC_MG_CYCLE_V, ierr); CHKERRA(ierr)
       CASE (2); CALL PCMGSetCycleType(matPETSC%pc, PC_MG_CYCLE_W, ierr); CHKERRA(ierr)
      END SELECT

      IF (perm) THEN
         ! Explicit Galerkin on permuted operators. Sets P, R, all level operators
         ! and the operator of the outer KSP.
         CALL PCMGSetGalerkin(matPETSC%pc, PC_MG_GALERKIN_NONE, ierr); CHKERRA(ierr)
         CALL build_explicit_MG(matPETSC, matPETSC%is_row_f, matPETSC%is_col_f)
      ELSE
         CALL PCMGSetGalerkin(matPETSC%pc, PC_MG_GALERKIN_BOTH, ierr); CHKERRA(ierr)
      END IF

      ! Coarse solve
      CALL PCMGGetCoarseSolve(matPETSC%pc, ksp_coarse, ierr); CHKERRA(ierr)
      CALL KSPSetType(ksp_coarse, KSPPREONLY, ierr); CHKERRA(ierr)
      CALL KSPGetPC(ksp_coarse, pc_coarse, ierr); CHKERRA(ierr)
      CALL PCSetType(pc_coarse, PCLU, ierr); CHKERRA(ierr)
      CALL PCFactorSetMatSolverType(pc_coarse, MATSOLVERMUMPS, ierr); CHKERRA(ierr)
      CALL set_preconditioner_reuse(ksp_coarse, ir, verbose=.FALSE.)

      DO i = 1, lssolver%mglevels - 1

         IF (.NOT. perm) THEN
            CALL PCMGSetInterpolation(matPETSC%pc, i, matPETSC%matP(i), ierr); CHKERRA(ierr)
            CALL PCMGSetRestriction(matPETSC%pc, i, matPETSC%matR(i), ierr); CHKERRA(ierr)
         END IF
         ! (with perm, P and R were already set inside build_explicit_MG)

         CALL PCMGGetSmoother(matPETSC%pc, i, ksp_smooth, ierr); CHKERRA(ierr)

         IF (i .EQ. lssolver%mglevels - 1 .AND. .NOT. perm) THEN
            CALL KSPSetOperators(ksp_smooth, matPETSC%matK, matPETSC%matK, ierr); CHKERRA(ierr)
         END IF

         CALL KSPSetType(ksp_smooth, KSPGMRES, ierr); CHKERRA(ierr)
         CALL KSPSetInitialGuessNonzero(ksp_smooth, PETSC_TRUE, ierr); CHKERRA(ierr)

         CALL KSPGetPC(ksp_smooth, innerpc, ierr); CHKERRA(ierr)
         CALL PCSetType(innerpc, PCASM, ierr); CHKERRA(ierr)
         CALL PCASMSetType(innerpc, PC_ASM_RESTRICT, ierr); CHKERRA(ierr)

         IF (i .EQ. lssolver%mglevels - 1) THEN
            CALL PCASMSetOverlap(innerpc, 3, ierr); CHKERRA(ierr)
         ELSE
            CALL PCASMSetOverlap(innerpc, 0, ierr); CHKERRA(ierr)
         END IF
         CALL PCASMSetLocalSubdomains(innerpc, 1, PETSC_NULL_IS, PETSC_NULL_IS, ierr); CHKERRA(ierr)

         IF (lssolver%kspnorm .EQ. 1) THEN
            CALL KSPSetNormType(ksp_smooth, KSP_NORM_PRECONDITIONED, ierr); CHKERRA(ierr)
         ELSEIF (lssolver%kspnorm .EQ. 2) THEN
            CALL KSPSetNormType(ksp_smooth, KSP_NORM_UNPRECONDITIONED, ierr); CHKERRA(ierr)
         END IF

         IF (i .EQ. lssolver%mglevels - 1) THEN
            CALL KSPSetTolerances(ksp_smooth, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, 10, ierr); CHKERRA(ierr)
         ELSE
            CALL KSPSetTolerances(ksp_smooth, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, 10, ierr); CHKERRA(ierr)
         END IF
      END DO
   END SUBROUTINE set_MG_ksp_pc_settings_reordering

   SUBROUTINE build_explicit_MG(matPETSC, is_row_f, is_col_f)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      IS, INTENT(IN)                   :: is_row_f, is_col_f
      PetscErrorCode               :: ierr
      INTEGER                      :: l, nl
      IS, ALLOCATABLE              :: rowperm(:), colperm(:)

      nl = lssolver%mglevels
      ALLOCATE (rowperm(0:nl - 1), colperm(0:nl - 1))
      CALL destroy_MG_level_operators(matPETSC)
      ALLOCATE (matPETSC%matKt(0:nl - 1), matPETSC%matPt(1:nl - 1), matPETSC%matRt(1:nl - 1))
      matPETSC%mg_perm_built = .FALSE.

      ! 1. No PETSc Galerkin: we provide every operator ourselves
      CALL PCMGSetGalerkin(matPETSC%pc, PC_MG_GALERKIN_NONE, ierr); CHKERRA(ierr)

      ! 2. Finest level: only graph ordering, computed once
      rowperm(nl - 1) = is_row_f
      colperm(nl - 1) = is_col_f
      CALL MatCreateSubMatrix(matPETSC%matK, rowperm(nl - 1), colperm(nl - 1), &
         MAT_INITIAL_MATRIX, matPETSC%matKt(nl - 1), ierr); CHKERRA(ierr)


      ! 3. Loop fine -> coarse
      DO l = nl - 1, 1, -1

         ! induced permutations of level l-1 (from the level directly above)
         ! Same permutation for rows and columns, otherwise the diagonal of Kt(l-1) is lost
         CALL induced_perm(matPETSC%matP(l), colperm(l), .TRUE., colperm(l - 1))   ! via columns of P
         rowperm(l - 1) = colperm(l - 1)

         ! permuted transfer operators (rows, cols):
         !   Pt = P(colperm_l , colperm_{l-1}),   Rt = R(rowperm_{l-1} , rowperm_l)
         CALL MatCreateSubMatrix(matPETSC%matP(l), colperm(l), colperm(l - 1), &
            MAT_INITIAL_MATRIX, matPETSC%matPt(l), ierr); CHKERRA(ierr)
         CALL MatCreateSubMatrix(matPETSC%matR(l), rowperm(l - 1), rowperm(l), &
            MAT_INITIAL_MATRIX, matPETSC%matRt(l), ierr); CHKERRA(ierr)

         ! explicit Galerkin product: Kt_{l-1} = Rt * Kt_l * Pt
         CALL MatMatMatMult(matPETSC%matRt(l), matPETSC%matKt(l), matPETSC%matPt(l), MAT_INITIAL_MATRIX, &
            PETSC_DEFAULT_REAL, matPETSC%matKt(l - 1), ierr); CHKERRA(ierr)

         ! give permuted transfer operators to PCMG (replace the natural ones)
         CALL PCMGSetInterpolation(matPETSC%pc, l, matPETSC%matPt(l), ierr); CHKERRA(ierr)
         CALL PCMGSetRestriction(matPETSC%pc, l, matPETSC%matRt(l), ierr); CHKERRA(ierr)
      END DO

      ! 4. Attach operators to every level's smoother and to the coarse solver
      CALL attach_level_operators(matPETSC, matPETSC%matKt)
      matPETSC%mg_perm_built = .TRUE.

   END SUBROUTINE build_explicit_MG

   ! Refresh the numerical values of the permuted level operators. Graph, permutations and
   ! transfer operators are kept: only Kt_nl-1 = K(perm,perm) and Kt_{l-1} = Rt*Kt_l*Pt are recomputed.
   SUBROUTINE refresh_MG_level_operators(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode               :: ierr
      INTEGER                      :: l, nl

      nl = lssolver%mglevels
      CALL MatCreateSubMatrix(matPETSC%matK, matPETSC%is_row_f, matPETSC%is_col_f, &
         MAT_REUSE_MATRIX, matPETSC%matKt(nl - 1), ierr); CHKERRA(ierr)
      DO l = nl - 1, 1, -1
         CALL MatMatMatMult(matPETSC%matRt(l), matPETSC%matKt(l), matPETSC%matPt(l), MAT_REUSE_MATRIX, &
            PETSC_DEFAULT_REAL, matPETSC%matKt(l - 1), ierr); CHKERRA(ierr)
      END DO
   END SUBROUTINE refresh_MG_level_operators

   SUBROUTINE destroy_MG_level_operators(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode               :: ierr
      INTEGER                      :: l

      IF (.NOT. matPETSC%mg_perm_built) RETURN
      DO l = LBOUND(matPETSC%matKt, 1), UBOUND(matPETSC%matKt, 1)
         CALL MatDestroy(matPETSC%matKt(l), ierr); CHKERRA(ierr)
      END DO
      DO l = LBOUND(matPETSC%matPt, 1), UBOUND(matPETSC%matPt, 1)
         CALL MatDestroy(matPETSC%matPt(l), ierr); CHKERRA(ierr)
         CALL MatDestroy(matPETSC%matRt(l), ierr); CHKERRA(ierr)
      END DO
      DEALLOCATE (matPETSC%matKt, matPETSC%matPt, matPETSC%matRt)
      matPETSC%mg_perm_built = .FALSE.
   END SUBROUTINE destroy_MG_level_operators

   SUBROUTINE attach_level_operators(matPETSC, Kt)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      Mat, INTENT(IN)              :: Kt(0:)
      KSP                          :: ksp_l
      PetscErrorCode               :: ierr
      INTEGER                      :: l

      CALL PCMGGetCoarseSolve(matPETSC%pc, ksp_l, ierr); CHKERRA(ierr)
      CALL KSPSetOperators(ksp_l, Kt(0), Kt(0), ierr); CHKERRA(ierr)
      DO l = 1, lssolver%mglevels - 1
         CALL PCMGGetSmoother(matPETSC%pc, l, ksp_l, ierr); CHKERRA(ierr)
         CALL KSPSetOperators(ksp_l, Kt(l), Kt(l), ierr); CHKERRA(ierr)
      END DO
   END SUBROUTINE attach_level_operators


   ! Builds the coarse permutation from the fine one.
   !   T       : transfer operator (P: n_f x n_c  -> by_columns = .true.
   !                                R: n_c x n_f  -> by_columns = .false.)
   !   is_fine : permutation of the fine level (same side as the T index being contracted)
   !   is_coarse (out): permutation of the coarse level
   SUBROUTINE induced_perm(T, is_fine, by_columns, is_coarse)
      Mat, INTENT(IN)              :: T
      IS, INTENT(IN)               :: is_fine
      LOGICAL, INTENT(IN)          :: by_columns
      IS, INTENT(OUT)              :: is_coarse
      PetscErrorCode               :: ierr
      Mat                          :: Tabs
      Vec                          :: pos, ones_f, num, den, key_seq
      VecScatter                   :: ctx
      PetscInt                     :: nf, nc, nloc, k, cs, ce, nf_inc, off
      PetscInt, POINTER            :: idx(:)
      PetscInt, ALLOCATABLE        :: sortidx(:)
      PetscReal, POINTER           :: keyv(:)
      PetscReal, ALLOCATABLE       :: keycopy(:)
      PetscScalar                  :: one

      one = 1.0
      CALL mat_abs_copy(T, Tabs)

      ! vectors in the fine and coarse spaces
      IF (by_columns) THEN
         CALL MatCreateVecs(Tabs, num, pos, ierr)    ! num: coarse (cols), pos: fine (rows)
      ELSE
         CALL MatCreateVecs(Tabs, pos, num, ierr)    ! pos: fine (cols), num: coarse (rows)
      END IF
      CALL VecDuplicate(pos, ones_f, ierr); CALL VecSet(ones_f, one, ierr)
      CALL VecDuplicate(num, den, ierr)

      ! pos(old fine index) = new position (1-based)
      CALL ISGetLocalSize(is_fine, nf, ierr); CHKERRA(ierr)
      CALL MPI_Scan(nf, nf_inc, 1, MPIU_INTEGER, MPI_SUM, PETSC_COMM_WORLD, ierr)
      off = nf_inc - nf   ! IS entries owned by lower ranks

      CALL ISGetIndicesF90(is_fine, idx, ierr); CHKERRA(ierr)
      DO k = 0, nf - 1
         CALL VecSetValue(pos, idx(k + 1), DBLE(off + k + 1), INSERT_VALUES, ierr); CHKERRA(ierr)
      END DO
      CALL ISRestoreIndicesF90(is_fine, idx, ierr); CHKERRA(ierr)
      CALL VecAssemblyBegin(pos, ierr); CALL VecAssemblyEnd(pos, ierr)

      ! weighted mean position of the parents of each coarse DOF
      IF (by_columns) THEN
         CALL MatMultTranspose(Tabs, pos, num, ierr)
         CALL MatMultTranspose(Tabs, ones_f, den, ierr)
      ELSE
         CALL MatMult(Tabs, pos, num, ierr)
         CALL MatMult(Tabs, ones_f, den, ierr)
      END IF
      CALL VecPointwiseDivide(num, num, den, ierr)    ! num now holds the key

      ! global sort of keys (gather to all, sort, keep the local slice)
      CALL VecGetSize(num, nc, ierr)
      CALL VecScatterCreateToAll(num, ctx, key_seq, ierr)
      CALL VecScatterBegin(ctx, num, key_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)
      CALL VecScatterEnd(ctx, num, key_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)
      CALL VecGetArrayReadF90(key_seq, keyv, ierr)
      ALLOCATE (keycopy(nc), sortidx(nc))
      keycopy = keyv(1:nc)
      DO k = 1, nc; sortidx(k) = k - 1; END DO
      CALL PetscSortRealWithPermutation(nc, keycopy, sortidx, ierr)  ! sortidx: 0-based, ascending key
      CALL VecRestoreArrayReadF90(key_seq, keyv, ierr)

      CALL VecGetOwnershipRange(num, cs, ce, ierr)
      nloc = ce - cs
      CALL ISCreateGeneral(PETSC_COMM_WORLD, nloc, sortidx(cs + 1:ce), PETSC_COPY_VALUES, &
         is_coarse, ierr)

      CALL VecScatterDestroy(ctx, ierr); CALL VecDestroy(key_seq, ierr)
      CALL VecDestroy(pos, ierr); CALL VecDestroy(ones_f, ierr)
      CALL VecDestroy(num, ierr); CALL VecDestroy(den, ierr); CALL MatDestroy(Tabs, ierr)
      DEALLOCATE (keycopy, sortidx)
   END SUBROUTINE induced_perm

   SUBROUTINE mat_abs_copy(T, Tabs)
      Mat, INTENT(IN)              :: T
      Mat, INTENT(OUT)             :: Tabs
      PetscErrorCode               :: ierr
      PetscInt                     :: rs, re, i, nc, maxnc
      PetscInt, ALLOCATABLE        :: cols(:)
      PetscScalar, ALLOCATABLE     :: vals(:)
      PetscInt                     :: idxm(1)

      CALL MatDuplicate(T, MAT_COPY_VALUES, Tabs, ierr); CHKERRA(ierr)
      CALL MatGetOwnershipRange(T, rs, re, ierr); CHKERRA(ierr)

      maxnc = 0
      DO i = rs, re - 1
         CALL MatGetRow(T, i, nc, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr); CHKERRA(ierr)
         maxnc = MAX(maxnc, nc)
         CALL MatRestoreRow(T, i, nc, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr); CHKERRA(ierr)
      END DO
      ALLOCATE (cols(MAX(maxnc, 1)), vals(MAX(maxnc, 1)))

      DO i = rs, re - 1
         CALL MatGetRow(T, i, nc, cols, vals, ierr); CHKERRA(ierr)
         idxm(1) = i
         IF (nc .GT. 0) THEN
            CALL MatSetValues(Tabs, 1, idxm, nc, cols(1:nc), ABS(vals(1:nc)), INSERT_VALUES, ierr); CHKERRA(ierr)
         END IF
         CALL MatRestoreRow(T, i, nc, cols, vals, ierr); CHKERRA(ierr)
      END DO

      CALL MatAssemblyBegin(Tabs, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
      CALL MatAssemblyEnd(Tabs, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)

   END SUBROUTINE mat_abs_copy

   SUBROUTINE finalize_MG_subsolvers(matPETSC, ir)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      INTEGER, INTENT(IN)               :: ir
      PetscErrorCode                   :: ierr
      KSP                              :: ksp_smooth, ksp_coarse, ksp_redundant
      KSP, POINTER                     :: subksps(:)
      PC                               :: innerpc, subpc
      PetscInt                         :: n_local, first_local
      INTEGER                          :: i, j

      NULLIFY (subksps)

      ! Finalize coarse-level redundant inner solver.
      !call PCMGGetCoarseSolve(matPETSC%pc, ksp_coarse, ierr); CHKERRA(ierr)
      !call KSPGetPC(ksp_coarse, innerpc, ierr); CHKERRA(ierr)
      !call PCRedundantGetKSP(innerpc, ksp_redundant, ierr); CHKERRA(ierr)
      !call KSPSetType(ksp_redundant, KSPPREONLY, ierr); CHKERRA(ierr)
      !call KSPGetPC(ksp_redundant, subpc, ierr); CHKERRA(ierr)
      !call PCSetType(subpc, PCLU, ierr); CHKERRA(ierr)
      !call PCFactorSetFill(subpc, 20., ierr); CHKERRA(ierr)
      !call PCFactorSetLevels(subpc, 3, ierr); CHKERRA(ierr)

      DO i = 1, lssolver%mglevels - 1

         CALL PCMGGetSmoother(matPETSC%pc, i, ksp_smooth, ierr); CHKERRA(ierr)
         CALL KSPGetPC(ksp_smooth, innerpc, ierr); CHKERRA(ierr)
         CALL KSPSetInitialGuessNonzero(ksp_smooth, PETSC_FALSE, ierr); CHKERRA(ierr)


         CALL set_preconditioner_reuse(ksp_smooth, ir, verbose=.FALSE.)

         CALL PCASMGetSubKSP(innerpc, n_local, first_local, PETSC_NULL_KSP, ierr); CHKERRA(ierr)

         IF (n_local .GT. 0) THEN
            ALLOCATE (subksps(n_local))
            CALL PCASMGetSubKSP(innerpc, n_local, first_local, subksps, ierr); CHKERRA(ierr)

            DO j = 1, n_local
               CALL KSPSetType(subksps(j), KSPPREONLY, ierr); CHKERRA(ierr)
               CALL KSPGetPC(subksps(j), subpc, ierr); CHKERRA(ierr)
               IF (i .EQ. lssolver%mglevels - 1) THEN
                  CALL PCSetType(subpc, PCLU, ierr); CHKERRA(ierr)
                  CALL PCFactorSetMatSolverType(subpc,  MATSOLVERMKL_PARDISO, ierr)
               ELSE
                  CALL PCSetType(subpc, PCILU, ierr); CHKERRA(ierr)
                  CALL PCFactorSetLevels(subpc,1, ierr); CHKERRA(ierr)
               END IF
               !call PCFactorSetMatSolverType(subpc,  MATSOLVERMKL_PARDISO, ierr)
               !call PCFactorSetMatSolverType(subpc,  MATSOLVERMUMPS, ierr)
               !call PCFactorSetMatOrderingType(subpc, 'natural', ierr); CHKERRA(ierr)
               ! ILU(0) was too weak per sweep; higher fill level gives each smoothing pass more punch.
               !call PCFactorSetLevels(subpc, 0, ierr); CHKERRA(ierr)
               !call PCFactorSetFill(subpc, 5., ierr); CHKERRA(ierr)
            END DO

            DEALLOCATE (subksps)
            NULLIFY (subksps)
         END IF
      END DO
   END SUBROUTINE finalize_MG_subsolvers

   SUBROUTINE set_BJ_ksp_pc_settings(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode                   :: ierr
      PetscInt                         :: n_local, blocks
      PetscInt, ALLOCATABLE            :: indices(:)

      CALL MatGetLocalSize(matPETSC%matK, PETSC_NULL_INTEGER, n_local, ierr); CHKERRA(ierr)

      blocks = 1 ! n_local/(phys%neq*mesh%Nnodesperface)
      WRITE (*, *) 'Number of blocks for Block Jacobi preconditioner: ', blocks

      ALLOCATE (indices(blocks + 1))
      indices = n_local/blocks
      CALL PCBJacobiSetLocalBlocks(matPETSC%pc, blocks, indices, ierr); CHKERRA(ierr)
      DEALLOCATE (indices)

   END SUBROUTINE set_BJ_ksp_pc_settings

   SUBROUTINE finalize_BJ_subsolvers(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode                   :: ierr
      PetscInt                         :: n_local, blocks
      KSP, POINTER                     :: subksps(:)
      PC                               :: subpc
      INTEGER                          :: i

      NULLIFY (subksps)

      CALL MatGetLocalSize(matPETSC%matK, PETSC_NULL_INTEGER, n_local, ierr); CHKERRA(ierr)

      blocks = 1 ! n_local/(phys%neq*mesh%Nnodesperface)
      ALLOCATE (subksps(blocks))
      CALL PCBJacobiGetSubKSP(matPETSC%pc, blocks, n_local, subksps, ierr); CHKERRA(ierr)

      WRITE (*, *) 'Setting up Block Jacobi preconditioner...'
      DO i = 1, blocks
         CALL KSPSetType(subksps(i), KSPPREONLY, ierr); CHKERRA(ierr)
         CALL KSPGetPC(subksps(i), subpc, ierr); CHKERRA(ierr)
         CALL PCSetType(subpc, PCLU, ierr); CHKERRA(ierr)
         CALL PCFactorSetMatSolverType(subpc,  MATSOLVERMKL_PARDISO, ierr)
         !call PCFactorSetMatOrderingType(subpc, 'amd', ierr)
      END DO

      DEALLOCATE (subksps)
      NULLIFY (subksps)

   END SUBROUTINE finalize_BJ_subsolvers

   SUBROUTINE set_ASM_ksp_pc_settings(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode                   :: ierr
      !call PCASMSetTotalSubdomains(matPETSC%pc, MPIvar%glob_size*2, PETSC_NULL_INTEGER, PETSC_NULL_VEC, ierr)
      CALL PCASMSetOverlap(matPETSC%pc, 3, ierr); CHKERRA(ierr)
   END SUBROUTINE set_ASM_ksp_pc_settings

   SUBROUTINE finalize_ASM_subsolvers(matPETSC)
      TYPE(PETSC_STRUC), INTENT(INOUT) :: matPETSC
      PetscErrorCode                   :: ierr
      PetscInt                         :: n_local, first_local
      KSP, POINTER                     :: subksps(:)
      PC                               :: subpc
      INTEGER                          :: i

      NULLIFY (subksps)

      CALL PCASMGetSubKSP(matPETSC%pc, n_local, first_local, PETSC_NULL_KSP, ierr); CHKERRA(ierr)

      IF (n_local .GT. 0) THEN
         ALLOCATE (subksps(n_local))
         CALL PCASMGetSubKSP(matPETSC%pc, n_local, first_local, subksps, ierr); CHKERRA(ierr)

         DO i = 1, n_local

            CALL KSPSetType(subksps(i), KSPPREONLY, ierr); CHKERRA(ierr)
            CALL KSPGetPC(subksps(i), subpc, ierr); CHKERRA(ierr)
            !call PCFactorSetMatOrderingType(subpc, 'natural', ierr); CHKERRA(ierr)
            CALL PCSetType(subpc, PCLU, ierr); CHKERRA(ierr)
            !call PCFactorReorderForNonzeroDiagonal(subpc, 1.0e-5, ierr); CHKERRA(ierr)
            !call PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO, ierr); CHKERRA(ierr)
            CALL PCFactorSetMatSolverType(subpc,  MATSOLVERMKL_PARDISO, ierr); CHKERRA(ierr)
            !call PCFactorSetLevels(subpc, 2, ierr); CHKERRA(ierr)
            !call PCFactorSetFill(subpc, 5.0, ierr); CHKERRA(ierr)
            !call PCSetType(subpc, PCHYPRE, ierr); CHKERRA(ierr)
            !call PCHypreSetType(subpc, "euclid", ierr); CHKERRA(ierr)

            !call PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-mat_mkl_pardiso_1", "2", ierr); CHKERRA(ierr)
            !call PCSetFromOptions(subpc, ierr); CHKERRA(ierr)


         END DO

         DEALLOCATE (subksps)
         NULLIFY (subksps)
      END IF

   END SUBROUTINE finalize_ASM_subsolvers

   SUBROUTINE PETSC_retrieve_array(vector, aux_sol, size)
      Vec, INTENT(IN)              :: vector
      INTEGER, INTENT(IN)          :: size
      PetscScalar, INTENT(INOUT)     :: aux_sol(size)
      PetscScalar, POINTER         :: data_array(:)
      PetscErrorCode               :: ierr
      ! copy the vector values into the array
      CALL VecGetArrayReadF90(vector, data_array, ierr); CHKERRA(ierr)
      aux_sol = data_array
      CALL VecRestoreArrayReadF90(vector, data_array, ierr); CHKERRA(ierr)
   END SUBROUTINE PETSC_retrieve_array

   SUBROUTINE PETSC_set_array(vector, aux_sol, size)
      Vec, INTENT(INOUT)           :: vector
      INTEGER, INTENT(IN)          :: size
      PetscScalar, INTENT(IN)      :: aux_sol(size)
      PetscScalar, POINTER         :: data_array(:)
      PetscErrorCode               :: ierr
      ! copy the vector values into the array
      CALL VecGetArrayF90(vector, data_array, ierr); CHKERRA(ierr)
      data_array = aux_sol
      CALL VecRestoreArrayF90(vector, data_array, ierr); CHKERRA(ierr)
   END SUBROUTINE PETSC_set_array

   SUBROUTINE compute_residue
      Vec                          :: true_residue_vec
      PetscReal                    :: alpha = -1
      PetscScalar                  :: true_residue, norm
      PetscErrorCode               :: ierr

      CALL MatCreateVecs(matPETSC%matK, true_residue_vec, PETSC_NULL_VEC, ierr); CHKERRA(ierr)
      CALL VecScale(matPETSC%rhs_vec, alpha, ierr); CHKERRA(ierr)
      CALL MatMultAdd(matPETSC%matK, matPETSC%solPETSC_vec, matPETSC%rhs_vec, true_residue_vec, ierr); CHKERRA(ierr)
      CALL VecScale(matPETSC%rhs_vec, alpha, ierr); CHKERRA(ierr)
      CALL VecNorm(true_residue_vec, NORM_2, true_residue, ierr); CHKERRA(ierr)
      CALL VecNorm(matPETSC%rhs_vec, NORM_2, norm, ierr); CHKERRA(ierr)
      matPETSC%residue = true_residue/norm
      CALL VecDestroy(true_residue_vec, ierr); CHKERRA(ierr)

   END SUBROUTINE compute_residue

   SUBROUTINE set_ksp_pc_type(ksp, pc)
      KSPType, INTENT(OUT)         :: ksp
      PCType, INTENT(OUT)          :: pc
      INTEGER                      :: ierr

      SELECT CASE (lssolver%kspmethd)
       CASE ('KSPCG')
         ksp = KSPCG
       CASE ('KSPBCGS')
         ksp = KSPBCGS
       CASE ('KSPGMRES')
         ksp = KSPGMRES
       CASE ('KSPFGMRES')
         ksp = KSPFGMRES
       CASE ('KSPMINRES')
         ksp = KSPMINRES
       CASE ('KSPPREONLY')
         ksp = KSPPREONLY
       CASE ('KSPRICHARDSON')
         ksp = KSPRICHARDSON
       CASE default
         IF (MPIvar%glob_id .EQ. 0) THEN
            WRITE (*, *) "KSP METHOD NOT in LIST, add it and recompile"
            STOP
         END IF
      END SELECT

      SELECT CASE (lssolver%pctype)
       CASE ('PCNONE')
         pc = PCNONE
       CASE ('PCJACOBI')
         pc = PCJACOBI
       CASE ('PCSOR')
         pc = PCSOR
       CASE ('PCLU')
         pc = PCLU
       CASE ('PCILU')
         pc = PCILU
       CASE ('PCGAMG')
         pc = PCGAMG
       CASE ('PCMG')
         pc = PCMG
       CASE ('PCBJACOBI')
         pc = PCBJACOBI
       CASE ('PCASM')
         pc = PCASM
       CASE ('PCFIELDSPLIT')
         pc = PCFIELDSPLIT
       CASE ('PCKSP')
         pc = PCKSP
       CASE ('PCGASM')
         pc = PCGASM
       CASE ('PCAMG')
         pc = PCHYPRE
       CASE default
         IF (MPIvar%glob_id .EQ. 0) THEN
            WRITE (*, *) "PC TYPE NOT in LIST, add it and recompile"
            CALL MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
         END IF
      END SELECT
   END SUBROUTINE set_ksp_pc_type

   ! Setup and preallocate the matrices for the levels of the multigrid preconditioner, and fill them with the prolongator/restriction operators
   SUBROUTINE setup_MG()
      INTEGER                      :: i
      INTEGER                      :: orders_pol(lssolver%mglevels)
      PetscErrorCode               :: ierr
#ifdef TOR3D
      INTEGER                      :: orders_tor(lssolver%mglevels)
#endif

      ! Populate the arrays containing the poloidal polynomial orders for the levels of the multigrid preconditioner.
      ! The coarsest level is the first one, the finest level is the last one (the one of the original problem)
      orders_pol(lssolver%mglevels) = refElPol%Ndeg
      DO i = 1, lssolver%mglevels - 1
         IF (refElPol%Ndeg/(2*i) .LT. 1) THEN
            orders_pol(lssolver%mglevels - i) = 1
         ELSE
            orders_pol(lssolver%mglevels - i) = refElPol%Ndeg/(2*i)
         END IF
      END DO

#ifdef TOR3D
      ! Populate the arrays containing the toroidal polynomial orders for the levels of the multigrid preconditioner.
      ! The coarsest level is the first one, the finest level is the last one (the one of the original problem)
      orders_tor(lssolver%mglevels) = refElTor%Ndeg
      DO i = 1, lssolver%mglevels - 1
         IF (refElTor%Ndeg/(2*i) .LT. 1) THEN
            orders_tor(lssolver%mglevels - i) = 1
         ELSE
            orders_tor(lssolver%mglevels - i) = refElTor%Ndeg/(2*i)
         END IF
      END DO
#endif

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "Multigrid levels and polynomial orders for the prolongator/restriction operators:"
         DO i = 1, lssolver%mglevels
            WRITE (*, *) "Level ", i, ": poloidal order = ", orders_pol(i)
#ifdef TOR3D
            WRITE (*, *) "Level ", i, ": toroidal order = ", orders_tor(i)
#endif
         END DO
      END IF

      ! allocate the array of matrices R and P for the multigrid preconditioner for each level
      ALLOCATE (matPETSC%matP(lssolver%mglevels - 1))
      ALLOCATE (matPETSC%matR(lssolver%mglevels - 1))

      ! Loop over levels
      DO i = 1, lssolver%mglevels - 1

         ! Create the matrices R and P for the current level
         CALL MatCreate(PETSC_COMM_WORLD, matPETSC%matR(i), ierr); CHKERRA(ierr)
         CALL MatCreate(PETSC_COMM_WORLD, matPETSC%matP(i), ierr); CHKERRA(ierr)

         ! Setup and preallocate the matrices R and P for the current level
#ifndef TOR3D
#ifndef PARALL
         CALL setup_MG_serial_matrices_2D(orders_pol, i)
#else
         CALL setup_MG_parallel_matrices_2D(orders_pol, i)
#endif
#else
#ifndef PARALL

         CALL setup_MG_serial_matrices_3D(orders_pol, orders_tor, i)
#else
         CALL setup_MG_parallel_matrices_3D(orders_pol, orders_tor, i)
#endif
#endif
         ! Set some options for the matrices R and P for the current level
         CALL MatSetOption(matPETSC%matR(i), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr); CHKERRA(ierr)
         CALL MatSetOption(matPETSC%matR(i), MAT_ROW_ORIENTED, PETSC_TRUE, ierr); CHKERRA(ierr)
         CALL MatSetOption(matPETSC%matP(i), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr); CHKERRA(ierr)
         CALL MatSetOption(matPETSC%matP(i), MAT_ROW_ORIENTED, PETSC_TRUE, ierr); CHKERRA(ierr)

         ! Fill the matrices R and P for the current level
#ifndef TOR3D
         CALL create_MG_structures_2D(orders_pol(i), orders_pol(i + 1), orders_pol(lssolver%mglevels), i)
#else
         CALL create_MG_structures_3D(orders_pol(i), orders_pol(i + 1), orders_pol(lssolver%mglevels), orders_tor(i), orders_tor(i + 1), i)
#endif
      END DO

   END SUBROUTINE setup_MG

#ifndef TOR3D
#ifndef PARALL
   SUBROUTINE setup_MG_serial_matrices_2D(orders, i)
      INTEGER, INTENT(IN)          :: orders(:)
      INTEGER, INTENT(IN)          :: i
      INTEGER                      :: nDoF_PH, nDoF_PL, nrows, ncols
      PetscErrorCode               :: ierr

      ! Set type of the matrix
      CALL MatSetType(matPETSC%matR(i), MATSEQAIJ, ierr); CHKERRA(ierr)
      CALL MatSetType(matPETSC%matP(i), MATSEQAIJ, ierr); CHKERRA(ierr)

      ! nDoF for each edge for PL and PH
      nDoF_PL = (orders(i) + 1)*phys%Neq
      nDoF_PH = (orders(i + 1) + 1)*phys%Neq
      nrows = Mesh%Nfaces*nDoF_PL
      ncols = Mesh%Nfaces*nDoF_PH

      ! Set the size and preallocate
      CALL MatSetSizes(matPETSC%matR(i), nrows, ncols, nrows, ncols, ierr); CHKERRA(ierr)
      CALL MatSetSizes(matPETSC%matP(i), ncols, nrows, ncols, nrows, ierr); CHKERRA(ierr)

      ! Preallocation is easy, all rows have the same number of non-zeros
      CALL MatSEQAIJSetPreallocation(matPETSC%matR(i), nDoF_PH, PETSC_NULL_INTEGER, ierr); CHKERRA(ierr)
      CALL MatSEQAIJSetPreallocation(matPETSC%matP(i), nDoF_PL, PETSC_NULL_INTEGER, ierr); CHKERRA(ierr)

   END SUBROUTINE setup_MG_serial_matrices_2D
#else
   SUBROUTINE setup_MG_parallel_matrices_2D(orders, i)
      INTEGER, INTENT(IN)          :: orders(:)
      INTEGER, INTENT(IN)          :: i
      INTEGER                      :: nDoF_PH, nDoF_PL, nrows, ncols
      PetscErrorCode               :: ierr

      ! Set type of the matrix
      CALL MatSetType(matPETSC%matR(i), MATMPIAIJ, ierr); CHKERRA(ierr)
      CALL MatSetType(matPETSC%matP(i), MATMPIAIJ, ierr); CHKERRA(ierr)

      ! nDoF for each edge for PL and PH
      nDoF_PL = (orders(i) + 1)*phys%Neq
      nDoF_PH = (orders(i + 1) + 1)*phys%Neq
      nrows = (Mesh%Nfaces - Mesh%nghostfaces)*nDoF_PL
      ncols = (Mesh%Nfaces - Mesh%nghostfaces)*nDoF_PH

      ! Set the size and preallocate
      CALL MatSetSizes(matPETSC%matR(i), nrows, ncols, PETSC_DETERMINE, PETSC_DETERMINE, ierr); CHKERRA(ierr)
      CALL MatSetSizes(matPETSC%matP(i), ncols, nrows, PETSC_DETERMINE, PETSC_DETERMINE, ierr); CHKERRA(ierr)

      ! Preallocation is easy gladily, all rows have the same number of non-zeros and they are only on the diagonal
      CALL MatMPIAIJSetPreallocation(matPETSC%matR(i), nDoF_PH, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, ierr); CHKERRA(ierr)
      CALL MatMPIAIJSetPreallocation(matPETSC%matP(i), nDoF_PL, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, ierr); CHKERRA(ierr)

   END SUBROUTINE setup_MG_parallel_matrices_2D
#endif
#else
#ifndef PARALL
   SUBROUTINE setup_MG_serial_matrices_3D(orders_pol, orders_tor, i)
      INTEGER, INTENT(IN)          :: orders_pol(:), orders_tor(:)
      INTEGER, INTENT(IN)          :: i
      INTEGER                           :: iel, ifa, k, itor
      INTEGER                           :: nDof_tor_face_PH, nDof_tor_face_PL, nDoF_pol_elem_PH, nDoF_pol_elem_PL, nDof_tor_elem_PH, nDof_tor_elem_PL, nrows, ncols, Ntor, Nelems, Nfaces
      INTEGER, ALLOCATABLE              :: nnz_R(:), nnz_P(:)
      INTEGER                           :: idx_start, idx_end
      PetscErrorCode                    :: ierr

      ! Set type of the matrix
      CALL MatSetType(matPETSC%matR(i), MATSEQAIJ, ierr); CHKERRA(ierr)
      CALL MatSetType(matPETSC%matP(i), MATSEQAIJ, ierr); CHKERRA(ierr)

      ! Setup all variables needed for the preallocation of the matrices R and P for the current level of the multigrid preconditioner
      Nelems = Mesh%Nelems
      Nfaces = Mesh%Nfaces
      Ntor = numer%Ntor
      nDoF_pol_elem_PH = INT((orders_pol(i + 1) + 1)*(orders_pol(i + 1) + 2)*0.5*phys%Neq)
      nDoF_pol_elem_PL = INT((orders_pol(i) + 1)*(orders_pol(i) + 2)*0.5*phys%Neq)
      nDof_tor_face_PH = INT((orders_pol(i + 1) + 1)*(orders_tor(i + 1) + 1)*phys%Neq)
      nDof_tor_face_PL = INT((orders_pol(i) + 1)*(orders_tor(i) + 1)*phys%Neq)
      nDof_tor_elem_PH = Nelems*nDoF_pol_elem_PH + Nfaces*nDof_tor_face_PH
      nDof_tor_elem_PL = Nelems*nDoF_pol_elem_PL + Nfaces*nDof_tor_face_PL

      ! Number of rows for R and P (the sizes are transposed)
      nrows = Ntor*nDof_tor_elem_PL
      ncols = Ntor*nDof_tor_elem_PH

      ! allocate non zero arrays for R and P
      ALLOCATE (nnz_R(nrows))
      ALLOCATE (nnz_P(ncols))

      ! Loop over toroidal elements, fill in the poloidal face contribution to the non-zero pattern of R and P
      DO itor = 1, Ntor
         idx_start = (itor - 1)*nDof_tor_elem_PL + 1
         idx_end = idx_start + Nelems*nDoF_pol_elem_PL - 1

         nnz_R(idx_start:idx_end) = nDoF_pol_elem_PH

         idx_start = (itor - 1)*nDof_tor_elem_PH + 1
         idx_end = idx_start + Nelems*nDoF_pol_elem_PH - 1

         nnz_P(idx_start:idx_end) = nDoF_pol_elem_PL
      END DO

      ! Loop over toroidal elements, fill in the toroidal face contribution to the non-zero pattern of R and P
      DO itor = 1, Ntor

         idx_start = Nelems*nDoF_pol_elem_PL + (itor - 1)*nDof_tor_elem_PL + 1
         idx_end = idx_start + Nfaces*nDof_tor_face_PL - 1

         nnz_R(idx_start:idx_end) = nDof_tor_face_PH

         idx_start = Nelems*nDoF_pol_elem_PH + (itor - 1)*nDof_tor_elem_PH + 1
         idx_end = idx_start + Nfaces*nDof_tor_face_PH - 1

         nnz_P(idx_start:idx_end) = nDof_tor_face_PL
      END DO

      ! Set sizes and preallocate
      CALL MatSetSizes(matPETSC%matR(i), nrows, ncols, PETSC_DECIDE, PETSC_DECIDE, ierr); CHKERRA(ierr)
      CALL MatSetSizes(matPETSC%matP(i), ncols, nrows, PETSC_DECIDE, PETSC_DECIDE, ierr); CHKERRA(ierr)

      ! Preallocate the matrices R and P, nnz_R and nnz_P contain the number of non-zeros per row for R and P, respectively
      CALL MatSEQAIJSetPreallocation(matPETSC%matR(i), PETSC_DEFAULT_INTEGER, nnz_R, ierr); CHKERRA(ierr)
      CALL MatSEQAIJSetPreallocation(matPETSC%matP(i), PETSC_DEFAULT_INTEGER, nnz_P, ierr); CHKERRA(ierr)

      DEALLOCATE (nnz_R, nnz_P)

   END SUBROUTINE setup_MG_serial_matrices_3D
#else
   SUBROUTINE setup_MG_parallel_matrices_3D(orders_pol, orders_tor, i)
      INTEGER, INTENT(IN)               :: orders_pol(lssolver%mglevels), orders_tor(lssolver%mglevels)
      INTEGER, INTENT(IN)               :: i
      INTEGER                           :: iel, ifa, k, itor
      INTEGER                           :: nDof_tor_face_PH, nDof_tor_face_PL, nDoF_pol_elem_PH, nDoF_pol_elem_PL, nDof_tor_elem_PH, nDof_tor_elem_PL, nrows, ncols, Ntor, Nelems_global, Nfaces_global
      INTEGER                           :: nDof_tor_elem_global_PH, nDof_tor_elem_global_PL, nDof_tor_elem_local_PH, nDof_tor_elem_local_PL
      INTEGER, ALLOCATABLE              :: nnz_R(:), nnz_P(:), d_nnz_R(:), o_nnz_R(:), d_nnz_P(:), o_nnz_P(:)
      INTEGER                           :: idx_start, idx_end
      INTEGER                           :: istart, iend, cstart, cend
      PetscErrorCode                    :: ierr

      ! Set type of the matrix
      CALL MatSetType(matPETSC%matR(i), MATMPIAIJ, ierr); CHKERRA(ierr)
      CALL MatSetType(matPETSC%matP(i), MATMPIAIJ, ierr); CHKERRA(ierr)

      ! Setup all variables needed for the preallocation of the matrices R and P for the current level of the multigrid preconditioner
      Nelems_global = Mesh%Nel_glob ! Global number of elements WITHOUT ghost elements
      Nfaces_global = Mesh%Nfa_glob ! Global number of faces WITHOUT ghost faces
      Ntor = numer%Ntor
      nDoF_pol_elem_PH = INT((orders_pol(i + 1) + 1)*(orders_pol(i + 1) + 2)*0.5*phys%Neq)
      nDoF_pol_elem_PL = INT((orders_pol(i) + 1)*(orders_pol(i) + 2)*0.5*phys%Neq)
      nDof_tor_face_PH = INT((orders_pol(i + 1) + 1)*(orders_tor(i + 1) + 1)*phys%Neq)
      nDof_tor_face_PL = INT((orders_pol(i) + 1)*(orders_tor(i) + 1)*phys%Neq)
      nDof_tor_elem_global_PH = Nelems_global*nDoF_pol_elem_PH + Nfaces_global*nDof_tor_face_PH
      nDof_tor_elem_global_PL = Nelems_global*nDoF_pol_elem_PL + Nfaces_global*nDof_tor_face_PL

      ! Number of GLOBAL rows for R and P (the sizes are transposed)
      nrows = Ntor*nDof_tor_elem_global_PL
      ncols = Ntor*nDof_tor_elem_global_PH

      ALLOCATE (nnz_R(nrows))
      ALLOCATE (nnz_P(ncols))

      ! Loop over toroidal elements, fill in the poloidal face contribution to the non-zero pattern of R and P
      DO itor = 1, Ntor
         idx_start = (itor - 1)*nDof_tor_elem_global_PL + 1
         idx_end = idx_start + Nelems_global*nDoF_pol_elem_PL - 1

         nnz_R(idx_start:idx_end) = nDoF_pol_elem_PH

         idx_start = (itor - 1)*nDof_tor_elem_global_PH + 1
         idx_end = idx_start + Nelems_global*nDoF_pol_elem_PH - 1

         nnz_P(idx_start:idx_end) = nDoF_pol_elem_PL
      END DO

      ! Loop over toroidal elements, fill in the toroidal face contribution to the non-zero pattern of R and P
      DO itor = 1, Ntor

         idx_start = Nelems_global*nDoF_pol_elem_PL + (itor - 1)*nDof_tor_elem_global_PL + 1
         idx_end = idx_start + Nfaces_global*nDof_tor_face_PL - 1

         nnz_R(idx_start:idx_end) = nDof_tor_face_PH

         idx_start = Nelems_global*nDoF_pol_elem_PH + (itor - 1)*nDof_tor_elem_global_PH + 1
         idx_end = idx_start + Nfaces_global*nDof_tor_face_PH - 1

         nnz_P(idx_start:idx_end) = nDof_tor_face_PL
      END DO

      ! Set LOCAL sizes and preallocate
      nDof_tor_elem_local_PH = (Mesh%Nelems - Mesh%Nghostelems)*nDoF_pol_elem_PH + (Mesh%Nfaces - Mesh%Nghostfaces)*nDof_tor_face_PH
      nDof_tor_elem_local_PL = (Mesh%Nelems - Mesh%Nghostelems)*nDoF_pol_elem_PL + (Mesh%Nfaces - Mesh%Nghostfaces)*nDof_tor_face_PL

      nrows = Ntor*nDof_tor_elem_local_PL
      ncols = Ntor*nDof_tor_elem_local_PH

      CALL MatSetSizes(matPETSC%matR(i), nrows, ncols, PETSC_DECIDE, PETSC_DECIDE, ierr); CHKERRA(ierr)
      CALL MatSetSizes(matPETSC%matP(i), ncols, nrows, PETSC_DECIDE, PETSC_DECIDE, ierr); CHKERRA(ierr)

      ! Need to call matSetup to be able to use MatGetOwnershipRange
      CALL matSetup(matPETSC%matR(i), ierr); CHKERRA(ierr)
      CALL matSetup(matPETSC%matP(i), ierr); CHKERRA(ierr)

      ! Get local row and columns ranges for R
      CALL MatGetOwnershipRange(matPETSC%matR(i), istart, iend, ierr); CHKERRA(ierr)
      CALL MatGetOwnershipRangeColumn(matPETSC%matR(i), cstart, cend, ierr); CHKERRA(ierr)

      ! Compute the number of diagonal and off diagonal entries for R
      ALLOCATE (d_nnz_R(iend - istart), o_nnz_R(iend - istart))
      CALL compute_diag_offdiag_nnz_3D(Ntor, Nelems_global, Nfaces_global, nDof_tor_elem_global_PL, nDof_tor_elem_global_PH, nDoF_pol_elem_PL, nDoF_pol_elem_PH, nDof_tor_face_PL, nDof_tor_face_PH, &
         istart, iend, cstart, cend, d_nnz_R, o_nnz_R)
      ! allocate R with corresponding number of diagonal entries and off diagonal entries
      CALL MatMPIAIJSetPreallocation(matPETSC%matR(i), PETSC_DEFAULT_INTEGER, d_nnz_R, PETSC_DEFAULT_INTEGER, o_nnz_R, ierr); CHKERRA(ierr)

      ! Get local row and columns ranges for P
      CALL MatGetOwnershipRange(matPETSC%matP(i), istart, iend, ierr); CHKERRA(ierr)
      CALL MatGetOwnershipRangeColumn(matPETSC%matP(i), cstart, cend, ierr); CHKERRA(ierr)

      ! Compute the number of diagonal and off diagonal entries for P
      ALLOCATE (d_nnz_P(iend - istart), o_nnz_P(iend - istart))
      CALL compute_diag_offdiag_nnz_3D(Ntor, Nelems_global, Nfaces_global, nDof_tor_elem_global_PH, nDof_tor_elem_global_PL, nDoF_pol_elem_PH, nDoF_pol_elem_PL, nDof_tor_face_PH, nDof_tor_face_PL, &
         istart, iend, cstart, cend, d_nnz_P, o_nnz_P)
      ! allocate P with corresponding number of diagonal entries and off diagonal entries
      CALL MatMPIAIJSetPreallocation(matPETSC%matP(i), PETSC_DEFAULT_INTEGER, d_nnz_P, PETSC_DEFAULT_INTEGER, o_nnz_P, ierr); CHKERRA(ierr)

      DEALLOCATE (nnz_R, nnz_P, d_nnz_R, o_nnz_R, d_nnz_P, o_nnz_P)

   END SUBROUTINE setup_MG_parallel_matrices_3D

   SUBROUTINE compute_diag_offdiag_nnz_3D(Ntor, Nelems, Nfaces, tor_row_width, tor_col_width, elem_row_width, elem_col_width, face_row_width, face_col_width, &
      istart, iend, cstart, cend, d_nnz, o_nnz)
      INTEGER, INTENT(IN)          :: Ntor, Nelems, Nfaces
      INTEGER, INTENT(IN)          :: tor_row_width, tor_col_width
      INTEGER, INTENT(IN)          :: elem_row_width, elem_col_width
      INTEGER, INTENT(IN)          :: face_row_width, face_col_width
      INTEGER, INTENT(IN)          :: istart, iend, cstart, cend
      INTEGER, INTENT(OUT)         :: d_nnz(iend - istart), o_nnz(iend - istart)

      INTEGER                      :: itor, iel, ifa
      INTEGER                      :: row_start, row_end, local_row_start, local_row_end
      INTEGER                      :: col_start, col_end, diag_nnz

      d_nnz = 0
      o_nnz = 0

      ! Compute the number of diagonal and off diagonal entries for the poloidal faces part
      DO itor = 1, Ntor
         DO iel = 1, Nelems

            ! Compute starting row index
            row_start = (itor - 1)*tor_row_width + (iel - 1)*elem_row_width
            ! Compute ending row index
            row_end = row_start + elem_row_width - 1

            ! Compute the local row range for the current process (intersection between the row range of the current block and the row range owned by the process)
            local_row_start = MAX(row_start, istart)
            local_row_end = MIN(row_end, iend - 1)

            ! if there is an intersection
            IF (local_row_start .LE. local_row_end) THEN
               ! Compute starting colomn index
               col_start = (itor - 1)*tor_col_width + (iel - 1)*elem_col_width
               ! Compute ending column index
               col_end = col_start + elem_col_width - 1
               ! MIN(col_end, cend - 1) gives the last column index of the block that is owned by the process
               ! MAX(col_start, cstart) gives the first column index of the block that is owned by the process
               ! MIN(col_end, cend - 1) - MAX(col_start, cstart) + 1 gives the number of columns of the block that are owned by the process (number of diagonal entries)
               ! if the number of columns of the block that are owned by the process is negative, it means that there is no column of the block that is owned by the process, so we set the number of diagonal entries to 0 with the MAX(0, ...)
               diag_nnz = MAX(0, MIN(col_end, cend - 1) - MAX(col_start, cstart) + 1)
               ! The number of diagonal entries is the number of columns of the block that are owned by the process, the number of off-diagonal entries is the total number of columns of the block minus the number of diagonal entries
               d_nnz(local_row_start - istart + 1:local_row_end - istart + 1) = diag_nnz
               o_nnz(local_row_start - istart + 1:local_row_end - istart + 1) = elem_col_width - diag_nnz
            END IF
         END DO
      END DO

      ! Compute the number of diagonal and off diagonal entries for the toroidal faces part
      DO itor = 1, Ntor
         DO ifa = 1, Nfaces
            ! Compute starting row index
            row_start = Nelems*elem_row_width + (itor - 1)*tor_row_width + (ifa - 1)*face_row_width
            ! Compute ending row index
            row_end = row_start + face_row_width - 1
            ! Find intersection
            local_row_start = MAX(row_start, istart)
            local_row_end = MIN(row_end, iend - 1)
            ! if the intersection exists
            IF (local_row_start .LE. local_row_end) THEN
               ! Compute starting column index
               col_start = Nelems*elem_col_width + (itor - 1)*tor_col_width + (ifa - 1)*face_col_width
               ! Compute ending column index
               col_end = col_start + face_col_width - 1
               ! Compute the number of diagonal entries
               diag_nnz = MAX(0, MIN(col_end, cend - 1) - MAX(col_start, cstart) + 1)
               ! Assign the number of diagonal and off-diagonal entries
               d_nnz(local_row_start - istart + 1:local_row_end - istart + 1) = diag_nnz
               o_nnz(local_row_start - istart + 1:local_row_end - istart + 1) = face_col_width - diag_nnz
            END IF
         END DO
      END DO

   END SUBROUTINE compute_diag_offdiag_nnz_3D

#endif
#endif

#ifndef TOR3D
   SUBROUTINE create_MG_structures_2D(order_PL, order_PH, order_Pog, level_counter)
      USE reference_element, ONLY: create_reference_element, Vandermonde_1D, orthopoly1D, orthopoly1d_deriv, Reference_element_type
      USE LinearAlgebra, ONLY: solve_linear_system, solve_linear_system_sing

      INTEGER, INTENT(IN)                     :: order_PL, order_PH, order_Pog
      INTEGER, INTENT(IN)                     :: level_counter
      TYPE(Reference_element_type)            :: refEl_PL, refEl_PH
      INTEGER                                 :: ifa, ifa_counter, iel, ifl, i, idx, nq_PH, nq_PL, np_Pog, np_PH, np_PL, nDoF_PL, nDoF_PH, Nfaces, Nintfaces
      INTEGER                                 :: nod(order_Pog + 1)
      REAL*8                                  :: xy_Pog(order_Pog + 1, 2), xyg_PH(order_PH + 2, 2), xyg_PL(order_PL + 2, 2)
      REAL*8                                  :: V_Pog(order_Pog + 1, order_Pog + 1),  V_PL(order_PL + 1, order_PL + 1), p_PL(order_PL + 1), p_PH(order_PH + 1), pxi_PH(order_PH + 1), p_Pog(order_Pog + 1), pxi_Pog(order_Pog + 1)
      INTEGER, ALLOCATABLE, DIMENSION(:, :)   :: faces
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: wq_PH, wq_PL
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: xq_PH, xq_PL
      REAL*8, ALLOCATABLE, DIMENSION(:, :)    :: N_PogPH, N_PogPL, N_PH, N_PL, N_PLPH, dNdxi_PogPL, dNdxi_PogPH, M_PH, M_PL, M_PHPL, M_PLPH, R_f_block, P_f_block
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: J_PH, J_PL
      REAL*8, ALLOCATABLE, DIMENSION(:, :, :) :: R_f, P_f

      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! P_og refers to the original polynomial order (highest one), PH refers to the higher order, PL to the lower order
      ! N are the shape functions, M the mass matrices, J the jacobians, R_f and P_f the restriction and prolongation operators at the face level
      ! dNdxi are the derivatives of the shape functions with respect to the reference coordinate xi
      !!!!!!!!!!!!!!!!!!!!!!!!

      nod = 0
      xy_Pog = 0.d0
      xyg_PH = 0.d0
      xyg_PL = 0.d0
      V_Pog = 0.d0
      V_PL = 0.d0
      p_PL = 0.d0
      p_PH = 0.d0
      pxi_PH = 0.d0
      p_Pog = 0.d0
      pxi_Pog = 0.d0

      ! Scale the mesh coordinates to real physical coordinates
      mesh%X = mesh%X*phys%lscale

      ! Get total number of faces and internal faces
      Nfaces = Mesh%Nfaces
      Nintfaces = Mesh%Nintfaces

      ! Fill the faces array with internal and external faces
      ALLOCATE (faces(Nfaces, 2))
      faces(1:Nintfaces, :) = Mesh%intfaces(:, 1:2)
      faces(Nintfaces + 1:Nfaces, :) = Mesh%extfaces(:, 1:2)

      ! Create the reference element PH
      CALL create_reference_element(refEl_PH, 2, order_PH, verbose=0)

      ! Create the reference element PL
      CALL create_reference_element(refEl_PL, 2, order_PL, verbose=0)

      nq_PH = refEl_PH%NGauss1D; ! Number of quadrature points for PH (1D)
      nq_PL = refEl_PL%NGauss1D; ! Number of quadrature points for PL (1D)

      np_Pog = refElPol%Nfacenodes; ! Number of 1D nodes for original element
      np_PH = refEl_PH%Nfacenodes; ! Number of 1D nodes for PH
      np_PL = refEl_PL%Nfacenodes; ! Number of 1D nodes for PL

      ALLOCATE (wq_PH(nq_PH), wq_PL(nq_PL), xq_PH(nq_PH), xq_PL(nq_PL), N_PH(np_PH, nq_PH), N_PL(np_PL, nq_PL))

      wq_PH = refEl_PH%gauss_weights1D ! PH weights
      wq_PL = refEl_PL%gauss_weights1D ! PL weights
      xq_PH = refEl_PH%gauss_points1D ! PH Gauss coordinates
      xq_PL = refEl_PL%gauss_points1D ! PL Gauss coordinates
      N_PH = TRANSPOSE(refEl_PH%N1D) ! PH shape functions
      N_PL = TRANSPOSE(refEl_PL%N1D) ! PL shape functions

      ALLOCATE (N_PLPH(np_PL, nq_PH))
      ALLOCATE (N_PogPH(np_Pog, nq_PH))
      ALLOCATE (dNdxi_PogPH(np_Pog, nq_PH))
      ALLOCATE (N_PogPL(np_Pog, nq_PL))
      ALLOCATE (dNdxi_PogPL(np_Pog, nq_PL))
      N_PLPH = 0.0d0
      N_PogPH = 0.0d0
      dNdxi_PogPH = 0.0d0
      N_PogPL = 0.0d0
      dNdxi_PogPL = 0.0d0

      ! Compute the P_og shape functions and derivatives at PL gauss points
      CALL Vandermonde_1D(V_Pog, refElPol)
      V_Pog = TRANSPOSE(V_Pog)

      DO i = 1, nq_PL
         CALL orthopoly1d_deriv(xq_PL(i), refElPol%Ndeg, p_Pog, pxi_Pog)
         CALL solve_linear_system_sing(V_Pog, p_Pog, N_PogPL(:, i))
         CALL solve_linear_system_sing(V_Pog, pxi_Pog, dNdxi_PogPL(:, i))
      END DO

      ! Compute the P_og shape functions and derivatives at PH gauss points
      DO i = 1, nq_PH
         CALL orthopoly1d_deriv(xq_PH(i), refElPol%Ndeg, p_Pog, pxi_Pog)
         CALL solve_linear_system_sing(V_Pog, p_Pog, N_PogPH(:, i))
         CALL solve_linear_system_sing(V_Pog, pxi_Pog, dNdxi_PogPH(:, i))
      END DO

      ! Compute the PL shape functions at PH gauss points
      CALL Vandermonde_1D(V_PL, refEl_PL)
      V_PL = TRANSPOSE(V_PL)
      DO i = 1, nq_PH
         CALL orthopoly1D(xq_PH(i), refEl_PL%Ndeg, p_PL)
         CALL solve_linear_system_sing(V_PL, p_PL, N_PLPH(:, i))
      END DO

      nDoF_PH = np_PH*phys%neq
      nDoF_PL = np_PL*phys%neq

#ifndef PARALL
      ALLOCATE (R_f(Nfaces, np_PL, np_PH))
      ALLOCATE (P_f(Nfaces, np_PH, np_PL))
#else
      ALLOCATE (R_f(Nfaces - Mesh%Nghostfaces, np_PL, np_PH))
      ALLOCATE (P_f(Nfaces - Mesh%Nghostfaces, np_PH, np_PL))
#endif
      R_f = 0.0d0
      P_f = 0.0d0

      ALLOCATE (R_f_block(np_PL, np_PH))
      ALLOCATE (P_f_block(np_PH, np_PL))
      R_f_block = 0.0d0
      P_f_block = 0.0d0

      ALLOCATE (J_PH(nq_PH))
      ALLOCATE (J_PL(nq_PL))
      J_PH = 0.0d0
      J_PL = 0.0d0

      ! allocate Mass matrices
      ALLOCATE (M_PH(np_PH, np_PH), M_PL(np_PL, np_PL), M_PHPL(np_PH, np_PL), M_PLPH(np_PL, np_PH))
      M_PH = 0.0d0
      M_PL = 0.0d0
      M_PHPL = 0.0d0
      M_PLPH = 0.0d0

      ifa_counter = 1
      ! Loop over faces
      !!$OMP PARALLEL DEFAULT(SHARED) private(ifa, iel, ifl, i, idx, nod, xy_PH, p_PH, J_PH, J_PL, M_PH, M_PL, M_PHPL, M_PLPH, P_f_block, R_f_block)
      !!$OMP do
      DO ifa = 1, Nfaces

#ifdef PARALL
         IF (mesh%ghostfaces(ifa) .EQ. 1) CYCLE ! skip ghost faces
#endif
         ! Get element number
         iel = faces(ifa, 1)

         ! Get local face number in the element
         ifl = faces(ifa, 2)

         ! Get local node number in the elemnt
         nod = refElPol%Face_nodes(ifl, :)

         ! Coordinates of the nodes of the face in the original mesh
         xy_Pog = Mesh%X(Mesh%T(iel, nod), :)

         ! Compute Jacobian at each quadrature point PH
         DO i = 1, nq_PH
            J_PH(i) = NORM2(MATMUL(dNdxi_PogPH(:, i), xy_Pog));
         END DO
         ! Compute Jacobian at each quadrature point PL
         DO i = 1, nq_PL
            J_PL(i) = NORM2(MATMUL(dNdxi_PogPL(:, i), xy_Pog))
         END DO

         ! Apply axisymmetric correction if needed
         IF (switch%axisym) THEN
            xyg_PH = MATMUL(TRANSPOSE(N_PogPH), xy_Pog)
            xyg_PL = MATMUL(TRANSPOSE(N_PogPL), xy_Pog)
            J_PH = J_PH*xyg_PH(:, 1)
            J_PL = J_PL*xyg_PL(:, 1)
         END IF

         ! Assemble mass matrices
         CALL assemble_M_matrix(np_PH, np_PH, J_PH, wq_PH, N_PH, N_PH, M_PH)
         CALL assemble_M_matrix(np_PL, np_PL, J_PL, wq_PL, N_PL, N_PL, M_PL)
         CALL assemble_M_matrix(np_PH, np_PL, J_PH, wq_PH, N_PH, N_PLPH, M_PHPL)
         CALL assemble_M_matrix(np_PL, np_PH, J_PH, wq_PH, N_PLPH, N_PH, M_PLPH)

         ! Compute operators
         CALL solve_linear_system(M_PH, M_PHPL, P_f_block)
         CALL solve_linear_system(M_PL, M_PLPH, R_f_block)

#ifndef PARALL
         P_f(ifa, :, :) = P_f_block
         R_f(ifa, :, :) = R_f_block
#else
         !!$OMP CRITICAL (ifa_counter_lock)
         idx = ifa_counter
         ifa_counter = ifa_counter + 1
         !!$OMP end CRITICAL (ifa_counter_lock)

         P_f(idx, :, :) = P_f_block
         R_f(idx, :, :) = R_f_block
#endif

      END DO
      !!$OMP end do
      !!$OMP end PARALLEL

#ifndef PARALL
      CALL construct_global_operators_2D(R_f, P_f, Nfaces, phys%neq, np_PL, np_PH, nDoF_PL, nDoF_PH, level_counter)
#else
      CALL construct_global_operators_2D(R_f, P_f, Nfaces - Mesh%Nghostfaces, phys%neq, np_PL, np_PH, nDoF_PL, nDoF_PH, level_counter)
#endif
      DEALLOCATE (faces, wq_PH, wq_PL, xq_PH, xq_PL, N_PogPH, N_PogPL, N_PH, N_PL, N_PLPH, dNdxi_PogPL, dNdxi_PogPH, J_PH, J_PL, M_PH, M_PL, M_PHPL, M_PLPH, R_f, P_f, R_f_block, P_f_block)

      CALL free_reference_element_pol(refEl_PH)
      CALL free_reference_element_pol(refEl_PL)

      mesh%X = mesh%X/phys%lscale

   END SUBROUTINE create_MG_structures_2D

   SUBROUTINE construct_global_operators_2D(R_f, P_f, Nfaces, Neq, np_PL, np_PH, nDoF_PL, nDoF_PH, level_counter)
      USE LinearAlgebra, ONLY: col
      INTEGER, INTENT(IN)          :: Nfaces, Neq, np_PL, np_PH, nDoF_PL, nDoF_PH, level_counter
      REAL*8, INTENT(IN)           :: R_f(Nfaces, np_PL, np_PH), P_f(Nfaces, np_PH, np_PL)
      INTEGER                      :: i, j, k, ifa, idx, ifa_counter
      INTEGER                      :: Nfaces_loc
      INTEGER                      :: idx_PL(nDoF_PL), idx_PH(nDoF_PH), idx_PL_local(nDoF_PL), idx_PH_local(nDoF_PH)
      REAL*8                       :: R_ifa(np_PL, np_PH), P_ifa(np_PH, np_PL)
      REAL*8                       :: R_block(nDoF_PL, nDoF_PH), P_block(nDoF_PH, nDoF_PL)
      PetscErrorCode               :: ierr
      Vec                          :: v
      PetscScalar, POINTER         :: data_array(:)
      PetscScalar, ALLOCATABLE     :: R_vals(:), P_vals(:)

#ifndef PARALL
      Nfaces_loc = Nfaces
#else
      Nfaces_loc = Mesh%Nfaces
      ifa_counter = 1
#endif

      ! local indices (-1 because PETSC uses 0-based indexing)
      idx_PH_local = (/(i, i=0, nDoF_PH - 1)/)
      idx_PL_local = (/(i, i=0, nDoF_PL - 1)/)

      ALLOCATE (R_vals(nDoF_PL*nDoF_PH))
      ALLOCATE (P_vals(nDoF_PH*nDoF_PL))

      ! Assemble R and P matrix matrix
      DO ifa = 1, Nfaces_loc

#ifndef PARALL
         ! local operators
         R_ifa = R_f(ifa, :, :)
         P_ifa = P_f(ifa, :, :)

#else
         IF (Mesh%ghostfaces(ifa) .EQ. 1) CYCLE ! skip ghost faces
         ! local operators
         R_ifa = R_f(ifa_counter, :, :)
         P_ifa = P_f(ifa_counter, :, :)
         ifa_counter = ifa_counter + 1
#endif

         ! expand to interleaved variables
         CALL KRON(R_ifa, np_PL, np_PH, Neq, R_block)
         CALL KRON(P_ifa, np_PH, np_PL, Neq, P_block)

         ! Values to be inserted in PETSC format (row-major order)
         R_vals = col(TRANSPOSE(R_block))
         P_vals = col(TRANSPOSE(P_block))

#ifndef PARALL
         idx_PL = (ifa - 1)*nDoF_PL + idx_PL_local
         idx_PH = (ifa - 1)*nDoF_PH + idx_PH_local
#else
         idx_PL = (Mesh%loc2glob_fa(ifa) - 1)*nDoF_PL + idx_PL_local
         idx_PH = (Mesh%loc2glob_fa(ifa) - 1)*nDoF_PH + idx_PH_local
#endif

         CALL MatSetValues(matPETSC%matR(level_counter), nDoF_PL, idx_PL, nDoF_PH, idx_PH, R_vals, INSERT_VALUES, ierr); CHKERRA(ierr)
         CALL MatSetValues(matPETSC%matP(level_counter), nDoF_PH, idx_PH, nDoF_PL, idx_PL, P_vals, INSERT_VALUES, ierr); CHKERRA(ierr)
      END DO

      CALL MatAssemblyBegin(matPETSC%matR(level_counter), MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
      CALL MatAssemblyEnd(matPETSC%matR(level_counter), MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
      CALL MatAssemblyBegin(matPETSC%matP(level_counter), MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
      CALL MatAssemblyEnd(matPETSC%matP(level_counter), MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)

      DEALLOCATE (R_vals)
      DEALLOCATE (P_vals)

   END SUBROUTINE construct_global_operators_2D
#else
   SUBROUTINE create_MG_structures_3D(order_pol_PL, order_pol_PH, order_pol_Pog, order_tor_PL, order_tor_PH, level_counter)
      USE reference_element, ONLY: create_reference_element, Vandermonde_1D, orthopoly1D, orthopoly1d_deriv, Reference_element_type
      USE LinearAlgebra, ONLY: solve_linear_system, solve_linear_system_sing

      INTEGER, INTENT(IN)                     :: order_pol_PL, order_pol_PH, order_pol_Pog
      INTEGER, INTENT(IN)                     :: order_tor_PL, order_tor_PH
      INTEGER, INTENT(IN)                     :: level_counter
      TYPE(Reference_element_type)            :: refEl_pol_PL, refEl_pol_PH, refEl_tor_PL, refEl_tor_PH
      INTEGER                                 :: itor, iel, ifa, ifl, i, idx, Nfaces, Nintfaces, Nelems, Ntor, Neq
      INTEGER                                 :: nq_pol_1D_PH, nq_pol_1D_PL, nq_tor_1D_PH, nq_tor_1D_PL, nq_pol_2D_PH, nq_pol_2D_PL, nq_tor_2D_PH, nq_tor_2D_PL
      INTEGER                                :: np_pol_1D_Pog, np_pol_1D_PH, np_pol_1D_PL, np_tor_1D_PH, np_tor_1D_PL, np_pol_2D_Pog, np_pol_2D_PH, np_pol_2D_PL, np_tor_2D_PH, np_tor_2D_PL, nDoF_pol_PH, nDoF_pol_PL, nDoF_tor_PH, nDoF_tor_PL
      INTEGER                                 :: nod(order_pol_Pog + 1)
      REAL*8                                  :: htor
      REAL*8, ALLOCATABLE, DIMENSION(:, :)       :: xy_Pog
      REAL*8, ALLOCATABLE, DIMENSION(:, :)     :: xyg_PH, xyg_PL
      REAL*8, ALLOCATABLE, DIMENSION(:)        :: xydNorm_PH, xydNorm_PL
      REAL*8                                  :: V_Pog(order_pol_Pog + 1, order_pol_Pog + 1),  V_PL(order_pol_PL + 1, order_pol_PL + 1), p_PL(order_pol_PL + 1), p_PH(order_pol_PH + 1), pxi_PH(order_pol_PH + 1), p_Pog(order_pol_Pog + 1), pxi_Pog(order_pol_Pog + 1)
      INTEGER, ALLOCATABLE, DIMENSION(:, :)   :: faces
      REAL*8, ALLOCATABLE, DIMENSION(:, :)    :: N_pol_1D_PH, N_pol_1D_PL, N_tor_1D_PH, N_tor_1D_PL, N_pol_2D_PH, N_pol_2D_PL, N_pol_2D_PLPH, N_pol_2D_PogPH, dNdxi_pol_2D_PogPH, dNdeta_pol_2D_PogPH, N_pol_2D_PogPL, dNdxi_pol_2D_PogPL, dNdeta_pol_2D_PogPL,&
         N_pol_1D_PogPL, dNdxi_pol_1D_PogPL, N_pol_1D_PogPH, dNdxi_pol_1D_PogPH, N_pol_1D_PLPH, N_tor_1D_PLPH, N_tor_2D_PH, N_tor_2D_PL, N_tor_2D_PLPH,&
         M_PH, M_PL, M_PHPL, M_PLPH
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: wq_tor_2D_PH, wq_tor_2D_PL
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: J_PH, J_PL
      REAL*8, ALLOCATABLE, DIMENSION(:, :)    :: P, R
      REAL*8, ALLOCATABLE, DIMENSION(:, :, :) :: R_el, P_el, R_f, P_f
#ifdef PARALL
      INTEGER                      :: iel_counter, ifa_counter
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! P_og refers to the original polynomial order (highest one), PH refers to the higher order, PL to the lower order
      ! N are the shape functions, M the mass matrices, J the jacobians, R_f and P_f the restriction and prolongation operators at the face level
      ! dNdxi are the derivatives of the shape functions with respect to the reference coordinate xi
      !!!!!!!!!!!!!!!!!!!!!!!!

      nod = 0
      V_Pog = 0.d0
      V_PL = 0.d0
      p_PL = 0.d0
      p_PH = 0.d0
      pxi_PH = 0.d0
      p_Pog = 0.d0
      pxi_Pog = 0.d0

      ! Scale the mesh coordinates to real physical coordinates
      mesh%X = mesh%X*phys%lscale

      ! Get total number of faces and internal faces
      Ntor = numer%ntor
      Nelems = Mesh%Nelems
      Nfaces = Mesh%Nfaces
      Nintfaces = Mesh%Nintfaces

      Neq = phys%neq
      htor = numer%tmax/Ntor

      ! Fill the faces array with internal and external faces
      ALLOCATE (faces(Nfaces, 2))
      faces(1:Nintfaces, :) = Mesh%intfaces(:, 1:2)
      faces(Nintfaces + 1:Nfaces, :) = Mesh%extfaces(:, 1:2)

      ! Create the reference element PH
      CALL create_reference_element(refEl_pol_PH, 2, order_pol_PH, verbose=0)
      CALL create_reference_element(refEl_tor_PH, 1, order_tor_PH, verbose=0)

      ! Create the reference element PL
      CALL create_reference_element(refEl_pol_PL, 2, order_pol_PL, verbose=0)
      CALL create_reference_element(refEl_tor_PL, 1, order_tor_PL, verbose=0)

      ! Number of 1D and 2D Gauss points for poloidal and toroidal
      nq_pol_1D_PH = size(refEl_pol_PH%N1D, 1)
      nq_pol_1D_PL = size(refEl_pol_PL%N1D, 1)
      nq_tor_1D_PH = size(refEl_tor_PH%N1D, 1)
      nq_tor_1D_PL = size(refEl_tor_PL%N1D, 1)
      nq_pol_2D_PH = size(refEl_pol_PH%N2D, 1)
      nq_pol_2D_PL = size(refEl_pol_PL%N2D, 1)
      nq_tor_2D_PL = nq_pol_1D_PL*nq_tor_1D_PL
      nq_tor_2D_PH = nq_pol_1D_PH*nq_tor_1D_PH

      ! Number of 1D and 2D nodes for poloidal and toroidal
      np_pol_1D_Pog = size(refElPol%N1D, 2)
      np_pol_1D_PH = size(refEl_pol_PH%N1D, 2)
      np_pol_1D_PL = size(refEl_pol_PL%N1D, 2)
      np_tor_1D_PH = size(refEl_tor_PH%N1D, 2)
      np_tor_1D_PL = size(refEl_tor_PL%N1D, 2)
      np_pol_2D_PH = size(refEl_pol_PH%N2D, 2)
      np_pol_2D_Pog = size(refElPol%N2D, 2)
      np_pol_2D_PL = size(refEl_pol_PL%N2D, 2)
      np_tor_2D_PH = np_pol_1D_PH*np_tor_1D_PH
      np_tor_2D_PL = np_pol_1D_PL*np_tor_1D_PL

      ! Poloidal and toroidal 1D weights, coordinates and shape functions (at corresponding Gauss points)
      ! Poloidal and toroidal 2D weights, coordinates and shape functions
      ! Toroidal 1D PH shape functions and derivatives (at corresponding Gauss points)
      ! 1D and 2D poloidal shape functions derivatives (at corresponding Gauss points)
      ASSOCIATE (wq_pol_1D_PH => refEl_pol_PH%gauss_weights1D, &
         wq_pol_1D_PL => refEl_pol_PL%gauss_weights1D, &
         xq_pol_1D_PH => refEl_pol_PH%gauss_points1D, &
         xq_pol_1D_PL => refEl_pol_PL%gauss_points1D, &
         wq_tor_1D_PH => refEl_tor_PH%gauss_weights1D, &
         wq_tor_1D_PL => refEl_tor_PL%gauss_weights1D, &
         xq_tor_1D_PH => refEl_tor_PH%gauss_points1D, &
         xq_tor_1D_PL => refEl_tor_PL%gauss_points1D, &
         wq_pol_2D_PH => refEl_pol_PH%gauss_weights2D, &
         wq_pol_2D_PL => refEl_pol_PL%gauss_weights2D, &
         xq_pol_2D_PH => refEl_pol_PH%gauss_points2D, &
         xq_pol_2D_PL => refEl_pol_PL%gauss_points2D)

         CALL allocate_arrays()

         N_pol_1D_PH = TRANSPOSE(refEl_pol_PH%N1D)
         N_pol_1D_PL = TRANSPOSE(refEl_pol_PL%N1D)
         N_tor_1D_PH = TRANSPOSE(refEl_tor_PH%N1D)
         N_tor_1D_PL = TRANSPOSE(refEl_tor_PL%N1D)
         N_pol_2D_PH = TRANSPOSE(refEl_pol_PH%N2D)
         N_pol_2D_PL = TRANSPOSE(refEl_pol_PL%N2D)

         ! Compute the 2D PL shape functions at 2D PH gauss points (poloidal)
         CALL compute_P_sf_at_Gauss_points_2D(refEl_pol_PL, np_pol_2D_PL, nq_pol_2D_PH, xq_pol_2D_PH, N_pol_2D_PLPH)
         ! Compute the 2D P_og shape functions and derivatives at 2D PH Gauss points (poloidal)
         CALL compute_P_sf_at_Gauss_points_2D(refElPol, np_pol_2D_Pog, nq_pol_2D_PH, xq_pol_2D_PH, N_pol_2D_PogPH, dNdxi_pol_2D_PogPH, dNdeta_pol_2D_PogPH)
         ! Compute the 2D P_og shape functions and derivatives at 2D PL Gauss points (poloidal)
         CALL compute_P_sf_at_Gauss_points_2D(refElPol, np_pol_2D_Pog, nq_pol_2D_PL, xq_pol_2D_PL, N_pol_2D_PogPL, dNdxi_pol_2D_PogPL, dNdeta_pol_2D_PogPL)

         ! Compute 1D P_og shape functions at PL Gauss points (poloidal)
         CALL compute_P_sf_at_Gauss_points_1D(refElPol, np_pol_1D_Pog, nq_pol_1D_PL, xq_pol_1D_PL, N_pol_1D_PogPL, dNdxi_pol_1D_PogPL)
         ! Compute 1D _og shape functions at PH Gauss points (poloidal)
         CALL compute_P_sf_at_Gauss_points_1D(refElPol, np_pol_1D_Pog, nq_pol_1D_PH, xq_pol_1D_PH, N_pol_1D_PogPH, dNdxi_pol_1D_PogPH)
         ! % Compute 1D PL shape functions at PH gauss points (poloidal)
         CALL compute_P_sf_at_Gauss_points_1D(refEl_pol_PL, np_pol_1D_PL, nq_pol_1D_PH, xq_pol_1D_PH, N_pol_1D_PLPH)
         ! Compute 1D PL shape functions at PL gauss points (toroidal)
         CALL compute_P_sf_at_Gauss_points_1D(refEl_tor_PL, np_tor_1D_PL, nq_tor_1D_PH, xq_tor_1D_PH, N_tor_1D_PLPH)
         ! Compute 2D PH shape functions and weights at 2D PH gauss points(toroidal)
         CALL compute_toroidal_sf_and_weights_2D(np_pol_1D_PH, nq_pol_1D_PH, np_tor_1D_PH, nq_tor_1D_PH, N_pol_1D_PH, N_tor_1D_PH, wq_pol_1D_PH, wq_tor_1D_PH, N_tor_2D_PH, wq_tor_2D_PH)
         ! Compute 2D PL shape functions and weights at 2D PL gauss points(toroidal)
         CALL compute_toroidal_sf_and_weights_2D(np_pol_1D_PL, nq_pol_1D_PL, np_tor_1D_PL, nq_tor_1D_PL, N_pol_1D_PL, N_tor_1D_PL, wq_pol_1D_PL, wq_tor_1D_PL, N_tor_2D_PL, wq_tor_2D_PL)
         ! Compute 2D PL shape functions at 2D PH gauss points(toroidal)
         CALL compute_toroidal_sf_and_weights_2D(np_pol_1D_PL, nq_pol_1D_PH, np_tor_1D_PL, nq_tor_1D_PH, N_pol_1D_PLPH, N_tor_1D_PLPH, N_tor_2D = N_tor_2D_PLPH)

#ifndef PARALL
         ALLOCATE (R_el(Nelems, np_pol_2D_PL, np_pol_2D_PH))
         ALLOCATE (P_el(Nelems, np_pol_2D_PH, np_pol_2D_PL))
#else
         ALLOCATE (R_el(Nelems - Mesh%nghostelems, np_pol_2D_PL, np_pol_2D_PH))
         ALLOCATE (P_el(Nelems - Mesh%nghostelems, np_pol_2D_PH, np_pol_2D_PL))
         iel_counter = 1
#endif

         ALLOCATE (R(np_pol_2D_PL, np_pol_2D_PH))
         ALLOCATE (P(np_pol_2D_PH, np_pol_2D_PL))
         ALLOCATE (J_PH(nq_pol_2D_PH))
         ALLOCATE (J_PL(nq_pol_2D_PL))
         ALLOCATE (M_PH(np_pol_2D_PH, np_pol_2D_PH), M_PL(np_pol_2D_PL, np_pol_2D_PL), M_PHPL(np_pol_2D_PH, np_pol_2D_PL), M_PLPH(np_pol_2D_PL, np_pol_2D_PH))
         R_el = 0.0d0
         P_el = 0.0d0

         !-------------------------------------------------------------
         ! Construct the operators (R and P) for the poloidal elements
         !-------------------------------------------------------------
         IF (switch%axisym) THEN
            ALLOCATE (xyg_PH(nq_pol_2D_PH, 2), xyg_PL(nq_pol_2D_PL, 2))
         END IF
         ALLOCATE (xy_Pog(Mesh%Nnodesperelem, 2))
         xy_Pog = 0.0d0

         DO iel = 1, Nelems

#ifdef PARALL
            IF (mesh%ghostelems(iel) .EQ. 1) CYCLE ! skip ghost elements
#endif

            ! Coordinates of the nodes of the element in the original mesh
            xy_Pog = mesh%X(mesh%T(iel, :), :)

            ! Compute Jacobian at PH
            CALL compute_jacobian_2D_pol(np_pol_2D_PH, nq_pol_2D_PH, dNdxi_pol_2D_PogPH, dNdeta_pol_2D_PogPH, xy_Pog, J_PH)
            ! Compute Jacobian at PL
            CALL compute_jacobian_2D_pol(np_pol_2D_PH, nq_pol_2D_PL, dNdxi_pol_2D_PogPL, dNdeta_pol_2D_PogPL, xy_Pog, J_PL)

            !-----------------------------
            ! Axisymmetric correction
            !-----------------------------
            IF (switch%axisym) THEN
               ! Get coordinates of Gauss points in physical space
               xyg_PH = MATMUL(transpose(N_pol_2D_PogPH), xy_Pog)
               xyg_PL = MATMUL(transpose(N_pol_2D_PogPL), xy_Pog)
               J_PH = J_PH*xyg_PH(:, 1)
               J_PL = J_PL*xyg_PL(:, 1)
            END IF

            !------------------------------------------------------------
            ! Compute mass matrices and mixed mass matrices
            !------------------------------------------------------------
            CALL assemble_M_matrix(np_pol_2D_PH, np_pol_2D_PH, J_PH, wq_pol_2D_PH, N_pol_2D_PH, N_pol_2D_PH, M_PH)
            CALL assemble_M_matrix(np_pol_2D_PL, np_pol_2D_PL, J_PL, wq_pol_2D_PL, N_pol_2D_PL, N_pol_2D_PL, M_PL)
            CALL assemble_M_matrix(np_pol_2D_PH, np_pol_2D_PL, J_PH, wq_pol_2D_PH, N_pol_2D_PH, N_pol_2D_PLPH, M_PHPL)
            CALL assemble_M_matrix(np_pol_2D_PL, np_pol_2D_PH, J_PH, wq_pol_2D_PH, N_pol_2D_PLPH, N_pol_2D_PH, M_PLPH)

            !------------------------------------------------------------
            ! P = M_PH \ M_PHPL  (solve M_PH * P = M_PHPL)
            ! R = M_PL \ M_PLPH  (solve M_PL * R = M_PLPH)
            !------------------------------------------------------------
            CALL solve_linear_system(M_PH, M_PHPL, P)   ! Compute Prolongation operator
            CALL solve_linear_system(M_PL, M_PLPH, R)   ! Compute Restriction operator

            !--------------------
            ! Store operators
            !--------------------
#ifndef PARALL
            P_el(iel, :, :) = P
            R_el(iel, :, :) = R
#else

            P_el(iel_counter, :, :) = P
            R_el(iel_counter, :, :) = R
            iel_counter = iel_counter + 1
#endif

         END DO

         IF (switch%axisym) THEN
            DEALLOCATE (xyg_PH, xyg_PL)
         END IF
         DEALLOCATE (xy_Pog)
         DEALLOCATE (R, P)
         DEALLOCATE (J_PH, J_PL)
         DEALLOCATE (M_PH, M_PL, M_PHPL, M_PLPH)

         !-----------------------------------------------------------------
         ! Construct the operators (R and P) for each toroidal face
         !-----------------------------------------------------------------
#ifndef PARALL
         ALLOCATE (R_f(Nfaces, np_tor_2D_PL, np_tor_2D_PH))
         ALLOCATE (P_f(Nfaces, np_tor_2D_PH, np_tor_2D_PL))
#else
         ALLOCATE (R_f(Nfaces - mesh%Nghostfaces, np_tor_2D_PL, np_tor_2D_PH))
         ALLOCATE (P_f(Nfaces - mesh%Nghostfaces, np_tor_2D_PH, np_tor_2D_PL))
         ifa_counter = 1
#endif
         ALLOCATE (R(np_tor_2D_PL, np_tor_2D_PH))
         ALLOCATE (P(np_tor_2D_PH, np_tor_2D_PL))
         ALLOCATE (J_PH(nq_tor_2D_PH))
         ALLOCATE (J_PL(nq_tor_2D_PL))
         ALLOCATE(M_PH(np_tor_2D_PH, np_tor_2D_PH), M_PL(np_tor_2D_PL, np_tor_2D_PL), M_PHPL(np_tor_2D_PH, np_tor_2D_PL), M_PLPH(np_tor_2D_PL, np_tor_2D_PH))
         ALLOCATE (xy_Pog(Mesh%Nnodesperface, 2))
         ALLOCATE (xydNorm_PH(nq_pol_1D_PH))
         ALLOCATE (xydNorm_PL(nq_pol_1D_PL))
         xydNorm_PH = 0.0d0
         xydNorm_PL = 0.0d0
         xy_Pog = 0.0d0
         R_f = 0.0d0
         P_f = 0.0d0
         R = 0.0d0
         P = 0.0d0

         IF (switch%axisym) THEN
            ALLOCATE (xyg_PH(nq_pol_1D_PH, 2))
            ALLOCATE (xyg_PL(nq_pol_1D_PL, 2))
         END IF

         DO ifa = 1, Nfaces

#ifdef PARALL
            IF (mesh%ghostfaces(ifa) .EQ. 1) CYCLE ! skip ghost faces
#endif
            ! Element number
            iel = faces(ifa, 1)
            ! Face local numbering
            ifl = faces(ifa, 2)
            ! Face nodes in local numbering
            nod = refElPol%Face_nodes(ifl, :)
            ! Get coordinates of the face
            xy_Pog = mesh%X(mesh%T(iel, nod), :)

            ! Compute the 1D shape function derivative at the physical points of the face PH
            CALL compute_sf_derivatives_at_phys_points(dNdxi_pol_1D_PogPH, xy_Pog, xydNorm_PH)
            ! Compute the 1D shape function derivative at the physical points of the face PL
            CALL compute_sf_derivatives_at_phys_points(dNdxi_pol_1D_PogPL, xy_Pog, xydNorm_PL)

            !------------------------------------------------------------
            ! Build Jacobians for the toroidal face (with axisymmetric correction if needed)
            !------------------------------------------------------------
            IF (switch%axisym) THEN

               ! Compute physical coordinates of the PH Gauss points on the face
               xyg_PH = MATMUL(transpose(N_pol_1D_PogPH), xy_Pog)
               ! Compute physical coordinates of the PL Gauss points on the face
               xyg_PL = MATMUL(transpose(N_pol_1D_PogPL), xy_Pog)

               ! Compute 2D toroidal Jacobian at PH and PL gauss points
               CALL compute_jacobian_2D_tor(nq_pol_1D_PH, nq_tor_1D_PH, xydNorm_PH, htor, xyg_PH, J_PH)
               CALL compute_jacobian_2D_tor(nq_pol_1D_PL, nq_tor_1D_PL, xydNorm_PL, htor, xyg_PL, J_PL)

            ELSE
               ! Compute 2D toroidal Jacobian at PH and PL gauss points
               CALL compute_jacobian_2D_tor(nq_pol_1D_PH, nq_tor_1D_PH, xydNorm_PH, htor, J=J_PH)
               CALL compute_jacobian_2D_tor(nq_pol_1D_PL, nq_tor_1D_PL, xydNorm_PL, htor, J=J_PL)
            END IF

            !------------------------------------------------------------
            ! Compute mass matrices and mixed mass matrices for the toroidal face
            !------------------------------------------------------------
            CALL assemble_M_matrix(np_tor_2D_PH, np_tor_2D_PH, J_PH, wq_tor_2D_PH, N_tor_2D_PH, N_tor_2D_PH, M_PH)
            CALL assemble_M_matrix(np_tor_2D_PL, np_tor_2D_PL, J_PL, wq_tor_2D_PL, N_tor_2D_PL, N_tor_2D_PL, M_PL)
            CALL assemble_M_matrix(np_tor_2D_PH, np_tor_2D_PL, J_PH, wq_tor_2D_PH, N_tor_2D_PH, N_tor_2D_PLPH, M_PHPL)
            CALL assemble_M_matrix(np_tor_2D_PL, np_tor_2D_PH, J_PH, wq_tor_2D_PH, N_tor_2D_PLPH, N_tor_2D_PH, M_PLPH)

            !------------------------------------------------------------
            ! P = M_PH \ M_PHPL  (solve M_PH * P = M_PHPL)
            ! R = M_PL \ M_PLPH  (solve M_PL * R = M_PLPH)
            !------------------------------------------------------------

            CALL solve_linear_system(M_PH, M_PHPL, P) ! Compute Prolongation operator
            CALL solve_linear_system(M_PL, M_PLPH, R) ! Compute Restriction operator

            !------------------------------------------------------------
            ! Store operators:
            !------------------------------------------------------------

#ifndef PARALL
            P_f(ifa, :, :) = P
            R_f(ifa, :, :) = R
#else
            P_f(ifa_counter, :, :) = P
            R_f(ifa_counter, :, :) = R
            ifa_counter = ifa_counter + 1
#endif
         END DO

         IF (switch%axisym) THEN
            DEALLOCATE (xyg_PH)
            DEALLOCATE (xyg_PL)
         END IF

         DEALLOCATE (xy_Pog)
         DEALLOCATE (xydNorm_PH, xydNorm_PL)
         DEALLOCATE (R, P)
         DEALLOCATE (J_PH, J_PL)
         DEALLOCATE (M_PH, M_PL, M_PHPL, M_PLPH)
      END ASSOCIATE

#ifndef PARALL
      CALL construct_global_operators_3D(R_el, P_el, R_f, P_f, Ntor, Nelems, Nfaces, np_pol_2D_PL, np_pol_2D_PH,  np_tor_2D_PL, np_tor_2D_PH, np_pol_1D_PL, np_pol_1D_PH, np_tor_1D_PL, np_tor_1D_PH, Neq, level_counter)
#else
      CALL construct_global_operators_3D(R_el, P_el, R_f, P_f, Ntor, Nelems-Mesh%Nghostelems, Nfaces-Mesh%Nghostfaces, np_pol_2D_PL, np_pol_2D_PH, np_tor_2D_PL, np_tor_2D_PH, np_pol_1D_PL, np_pol_1D_PH, np_tor_1D_PL, np_tor_1D_PH, Neq, level_counter)
#endif
      mesh%X = mesh%X/phys%lscale

      CALL free_arrays()

   CONTAINS

      SUBROUTINE allocate_arrays()
         ALLOCATE (N_pol_1D_PH(np_pol_1D_PH, nq_pol_1D_PH))
         ALLOCATE (N_pol_1D_PL(np_pol_1D_PL, nq_pol_1D_PL))
         ALLOCATE (N_tor_1D_PH(np_tor_1D_PH, nq_tor_1D_PH))
         ALLOCATE (N_tor_1D_PL(np_tor_1D_PL, nq_tor_1D_PL))
         ALLOCATE (N_pol_2D_PH(np_pol_2D_PH, nq_pol_2D_PH))
         ALLOCATE (N_pol_2D_PL(np_pol_2D_PL, nq_pol_2D_PL))
         ALLOCATE (N_pol_2D_PLPH(np_pol_2D_PL, nq_pol_2D_PH))
         ALLOCATE (N_pol_2D_PogPH(np_pol_2D_Pog, nq_pol_2D_PH))
         ALLOCATE (dNdxi_pol_2D_PogPH(np_pol_2D_Pog, nq_pol_2D_PH))
         ALLOCATE (dNdeta_pol_2D_PogPH(np_pol_2D_Pog, nq_pol_2D_PH))
         ALLOCATE (N_pol_2D_PogPL(np_pol_2D_Pog, nq_pol_2D_PL))
         ALLOCATE (dNdxi_pol_2D_PogPL(np_pol_2D_Pog, nq_pol_2D_PL))
         ALLOCATE (dNdeta_pol_2D_PogPL(np_pol_2D_Pog, nq_pol_2D_PL))
         ALLOCATE (N_pol_1D_PogPL(np_pol_1D_Pog, nq_pol_1D_PL))
         ALLOCATE (dNdxi_pol_1D_PogPL(np_pol_1D_Pog, nq_pol_1D_PL))
         ALLOCATE (N_pol_1D_PogPH(np_pol_1D_Pog, nq_pol_1D_PH))
         ALLOCATE (dNdxi_pol_1D_PogPH(np_pol_1D_Pog, nq_pol_1D_PH))
         ALLOCATE (N_pol_1D_PLPH(np_pol_1D_PL, nq_pol_1D_PH))
         ALLOCATE (N_tor_1D_PLPH(np_tor_1D_PL, nq_tor_1D_PH))
         ALLOCATE (N_tor_2D_PH(np_tor_2D_PH, nq_tor_2D_PH))
         ALLOCATE (N_tor_2D_PLPH(np_tor_2D_PL, nq_tor_2D_PH))
         ALLOCATE (wq_tor_2D_PH(nq_tor_2D_PH))
         ALLOCATE (N_tor_2D_PL(np_tor_2D_PL, nq_tor_2D_PL))
         ALLOCATE (wq_tor_2D_PL(nq_tor_2D_PL))

         N_pol_1D_PH = 0.d0
         N_pol_1D_PL = 0.d0
         N_tor_1D_PH = 0.d0
         N_tor_1D_PL = 0.d0
         N_pol_2D_PH = 0.d0
         N_pol_2D_PL = 0.d0
         N_pol_2D_PLPH = 0.d0
         N_pol_2D_PogPH = 0.d0
         dNdxi_pol_2D_PogPH = 0.d0
         dNdeta_pol_2D_PogPH = 0.d0
         N_pol_2D_PogPL = 0.d0
         dNdxi_pol_2D_PogPL = 0.d0
         dNdeta_pol_2D_PogPL = 0.d0
         N_pol_1D_PogPL = 0.d0
         dNdxi_pol_1D_PogPL = 0.d0
         N_pol_1D_PogPH = 0.d0
         dNdxi_pol_1D_PogPH = 0.d0
         N_pol_1D_PLPH = 0.d0
         N_tor_1D_PLPH = 0.d0
         N_tor_2D_PH = 0.d0
         N_tor_2D_PLPH = 0.d0
         wq_tor_2D_PH = 0.d0
         N_tor_2D_PL = 0.d0
         wq_tor_2D_PL = 0.d0

      END SUBROUTINE allocate_arrays

      SUBROUTINE free_arrays()

         DEALLOCATE (N_pol_1D_PH)
         DEALLOCATE (N_pol_1D_PL)
         DEALLOCATE (N_tor_1D_PH)
         DEALLOCATE (N_tor_1D_PL)
         DEALLOCATE (N_pol_2D_PH)
         DEALLOCATE (N_pol_2D_PL)
         DEALLOCATE (N_pol_2D_PLPH)
         DEALLOCATE (N_pol_2D_PogPH)
         DEALLOCATE (dNdxi_pol_2D_PogPH)
         DEALLOCATE (dNdeta_pol_2D_PogPH)
         DEALLOCATE (N_pol_2D_PogPL)
         DEALLOCATE (dNdxi_pol_2D_PogPL)
         DEALLOCATE (dNdeta_pol_2D_PogPL)
         DEALLOCATE (N_pol_1D_PogPL)
         DEALLOCATE (dNdxi_pol_1D_PogPL)
         DEALLOCATE (N_pol_1D_PogPH)
         DEALLOCATE (dNdxi_pol_1D_PogPH)
         DEALLOCATE (N_pol_1D_PLPH)
         DEALLOCATE (N_tor_1D_PLPH)
         DEALLOCATE (N_tor_2D_PH)
         DEALLOCATE (N_tor_2D_PLPH)
         DEALLOCATE (wq_tor_2D_PH)
         DEALLOCATE (N_tor_2D_PL)
         DEALLOCATE (wq_tor_2D_PL)
         DEALLOCATE (faces)
         DEALLOCATE (R_el)
         DEALLOCATE (P_el)
         DEALLOCATE (R_f)
         DEALLOCATE (P_f)

         CALL free_reference_element_pol(refEl_pol_PH)
         CALL free_reference_element_pol(refEl_pol_PL)
         CALL free_reference_element_tor(refEl_tor_PH)
         CALL free_reference_element_tor(refEl_tor_PL)

      END SUBROUTINE free_arrays

   END SUBROUTINE create_MG_structures_3D

   SUBROUTINE construct_global_operators_3D(R_el, P_el, R_f, P_f, Ntor, Nelems, Nfaces, np_pol_2D_PL, np_pol_2D_PH, np_tor_2D_PL, np_tor_2D_PH, np_pol_1D_PL, np_pol_1D_PH, np_tor_1D_PL, np_tor_1D_PH, Neq, level_counter)
      USE LinearAlgebra, ONLY: col
      IMPLICIT NONE
      INTEGER, INTENT(IN)          :: Ntor, Nelems, Nfaces
      INTEGER, INTENT(IN)          :: np_pol_2D_PL, np_pol_2D_PH, np_tor_2D_PL, np_tor_2D_PH
      INTEGER, INTENT(IN)          :: np_pol_1D_PL, np_pol_1D_PH
      INTEGER, INTENT(IN)          :: np_tor_1D_PL, np_tor_1D_PH
      INTEGER, INTENT(IN)          :: Neq, level_counter

      REAL*8, INTENT(IN)           :: R_el(Nelems, np_pol_2D_PL, np_pol_2D_PH)
      REAL*8, INTENT(IN)           :: P_el(Nelems, np_pol_2D_PH, np_pol_2D_PL)
      REAL*8, INTENT(IN)           :: R_f(Nfaces, np_tor_2D_PL, np_tor_2D_PH)
      REAL*8, INTENT(IN)           :: P_f(Nfaces, np_tor_2D_PH, np_tor_2D_PL)

      INTEGER                      :: itor, iel, ifa, i
      INTEGER                      :: Nelems_loc, Nelems_glob, Nfaces_loc, Nfaces_glob, nDoF_pol_1D_PH, nDoF_pol_1D_PL, nDoF_pol_2D_PH, nDoF_pol_2D_PL, nDoF_tor_2D_PH, nDoF_tor_2D_PL, nDoF_tor_3D_PH, nDoF_tor_3D_PL, nDoF_tor_elem_PH, nDoF_tor_elem_PL, nDoF_tor_PH, nDoF_tor_PL
      INTEGER                      :: indices_pol_2D_PH(np_pol_2D_PH*Neq), indices_pol_2D_PL(np_pol_2D_PL*Neq), indices_tor_2D_PH(np_tor_2D_PH*Neq), indices_tor_2D_PL(np_tor_2D_PL*Neq)
      INTEGER, ALLOCATABLE         :: idx_PH(:), idx_PL(:)
      REAL*8                       :: R_iel(np_pol_2D_PL, np_pol_2D_PH), P_iel(np_pol_2D_PH, np_pol_2D_PL)
      REAL*8                       :: R_ifa(np_tor_2D_PL, np_tor_2D_PH), P_ifa(np_tor_2D_PH, np_tor_2D_PL)
      REAL*8, ALLOCATABLE          :: R_block(:, :), P_block(:, :), R_vals(:), P_vals(:)
      PetscErrorCode               :: ierr
      PetscViewer                  :: viewer

#ifdef PARALL
      INTEGER                      :: iel_counter, ifa_counter
#endif

#ifndef PARALL
      ! Global number of elements WITHOUT ghost elements
      Nelems_glob = Nelems
      ! Global number of faces WITHOUT ghost faces
      Nfaces_glob = Nfaces
      ! Local number of elemnets WITH ghost elements
      Nelems_loc = Nelems
      ! Local number of faces WITH ghost faces
      Nfaces_loc = Nfaces
#else
      ! Global number of elements WITHOUT ghost elements
      Nelems_glob = Mesh%Nel_glob
      ! Global number of faces WITHOUT ghost faces
      Nfaces_glob = Mesh%Nfa_glob
      ! Local number of elemnets WITH ghost elements
      Nelems_loc = Mesh%Nelems
      ! Local number of faces WITH ghost faces
      Nfaces_loc = Mesh%Nfaces
#endif
      ! Poloidal edge nDoF
      nDoF_pol_1D_PH = np_pol_1D_PH*Neq
      nDoF_pol_1D_PL = np_pol_1D_PL*Neq
      ! Poloidal face nDoF
      nDoF_pol_2D_PH = np_pol_2D_PH*Neq
      nDoF_pol_2D_PL = np_pol_2D_PL*Neq
      ! Toroidal face nDoF
      nDoF_tor_2D_PH = np_tor_2D_PH*Neq
      nDoF_tor_2D_PL = np_tor_2D_PL*Neq
      ! Toroidal volume nDoF
      nDoF_tor_3D_PH = np_pol_2D_PH*np_tor_1D_PH*Neq
      nDoF_tor_3D_PL = np_pol_2D_PL*np_tor_1D_PL*Neq
      ! Toroidal element nDoF only from elements
      nDoF_tor_elem_PH = Nelems_glob*nDoF_pol_2D_PH
      nDoF_tor_elem_PL = Nelems_glob*nDoF_pol_2D_PL
      ! Total nDoF per toroidal element
      nDoF_tor_PH = Nelems_glob*nDoF_pol_2D_PH + Nfaces_glob*nDoF_tor_2D_PH
      nDoF_tor_PL = Nelems_glob*nDoF_pol_2D_PL + Nfaces_glob*nDoF_tor_2D_PL

      ! Local indices for poloidal face (1, 2, ..., nDoF_pol_2D_P)), starts from 0 because of PETSc's 0-based indexing
      indices_pol_2D_PH = (/(i, i=0, nDoF_pol_2D_PH - 1)/)
      indices_pol_2D_PL = (/(i, i=0, nDoF_pol_2D_PL - 1)/)

      ! Local indices for toroidal face (0, 1, ..., nDoF_tor_2D_P - 1), starts from 0 because of PETSc's 0-based indexing
      indices_tor_2D_PH = (/(i, i=0, nDoF_tor_2D_PH - 1)/)
      indices_tor_2D_PL = (/(i, i=0, nDoF_tor_2D_PL - 1)/)

      ALLOCATE (R_block(nDoF_pol_2D_PL, nDoF_pol_2D_PH))
      ALLOCATE (P_block(nDoF_pol_2D_PH, nDoF_pol_2D_PL))
      ALLOCATE (R_vals(nDoF_pol_2D_PL*nDoF_pol_2D_PH))
      ALLOCATE (P_vals(nDoF_pol_2D_PH*nDoF_pol_2D_PL))
      ALLOCATE (idx_PH(nDoF_pol_2D_PH))
      ALLOCATE (idx_PL(nDoF_pol_2D_PL))

      !--------------------------------------------------------------
      ! Restrict and prolongate over elements
      !--------------------------------------------------------------

      DO itor = 1, Ntor
#ifdef PARALL
         iel_counter = 1
#endif
         DO iel = 1, Nelems_loc

#ifndef PARALL
            R_iel = R_el(iel, :, :)
            P_iel = P_el(iel, :, :)
#else
            IF (Mesh%ghostelems(iel) .EQ. 1) CYCLE ! skip ghost elements
            R_iel = R_el(iel_counter, :, :)
            P_iel = P_el(iel_counter, :, :)
            iel_counter = iel_counter + 1
#endif
            CALL KRON(R_iel, np_pol_2D_PL, np_pol_2D_PH, Neq, R_block)
            CALL KRON(P_iel, np_pol_2D_PH, np_pol_2D_PL, Neq, P_block)

            R_vals = col(TRANSPOSE(R_block))
            P_vals = col(TRANSPOSE(P_block))

#ifndef PARALL
            idx_PH = (itor - 1)*nDoF_tor_PH + (iel - 1)*nDoF_pol_2D_PH + indices_pol_2D_PH
            idx_PL = (itor - 1)*nDoF_tor_PL + (iel - 1)*nDoF_pol_2D_PL + indices_pol_2D_PL
#else
            idx_PH = (itor - 1)*nDoF_tor_PH + (Mesh%loc2glob_el(iel) - 1)*nDoF_pol_2D_PH + indices_pol_2D_PH
            idx_PL = (itor - 1)*nDoF_tor_PL + (Mesh%loc2glob_el(iel) - 1)*nDoF_pol_2D_PL + indices_pol_2D_PL
#endif

            ! The -1 in the indices is because PETSc expects 0-based indexing for the global DOF indices, while Fortran arrays are 1-based.
            CALL MatSetValues(matPETSC%matR(level_counter), nDoF_pol_2D_PL, idx_PL, nDoF_pol_2D_PH, idx_PH, R_vals, INSERT_VALUES, ierr); CHKERRA(ierr)
            CALL MatSetValues(matPETSC%matP(level_counter), nDoF_pol_2D_PH, idx_PH, nDoF_pol_2D_PL, idx_PL, P_vals, INSERT_VALUES, ierr); CHKERRA(ierr)

         END DO
      END DO

      DEALLOCATE (R_block)
      DEALLOCATE (P_block)
      DEALLOCATE (R_vals)
      DEALLOCATE (P_vals)
      DEALLOCATE (idx_PH)
      DEALLOCATE (idx_PL)

      CALL MPI_Barrier(PETSC_COMM_WORLD, ierr); CHKERRA(ierr)

      ALLOCATE (R_block(nDoF_tor_2D_PL, nDoF_tor_2D_PH))
      ALLOCATE (P_block(nDoF_tor_2D_PH, nDoF_tor_2D_PL))
      ALLOCATE (R_vals(nDoF_tor_2D_PL*nDoF_tor_2D_PH))
      ALLOCATE (P_vals(nDoF_tor_2D_PH*nDoF_tor_2D_PL))
      ALLOCATE (idx_PH(nDoF_tor_2D_PH))
      ALLOCATE (idx_PL(nDoF_tor_2D_PL))

      !--------------------------------------------------------------
      ! Restrict and prolongate over faces
      !--------------------------------------------------------------
      DO itor = 1, Ntor
#ifdef PARALL
         ifa_counter = 1
#endif
         DO ifa = 1, Nfaces_loc

#ifndef PARALL
            R_ifa = R_f(ifa, :, :)
            P_ifa = P_f(ifa, :, :)
#else
            IF (Mesh%ghostfaces(ifa) .EQ. 1) CYCLE ! skip ghost faces
            R_ifa = R_f(ifa_counter, :, :)
            P_ifa = P_f(ifa_counter, :, :)
            ifa_counter = ifa_counter + 1
#endif

            CALL KRON(R_ifa, np_tor_2D_PL, np_tor_2D_PH, Neq, R_block)
            CALL KRON(P_ifa, np_tor_2D_PH, np_tor_2D_PL, Neq, P_block)

            R_vals = col(TRANSPOSE(R_block))
            P_vals = col(TRANSPOSE(P_block))
#ifndef PARALL
            idx_PH = nDoF_tor_elem_PH + (itor - 1)*nDoF_tor_PH + (ifa - 1)*nDoF_tor_2D_PH + indices_tor_2D_PH
            idx_PL = nDoF_tor_elem_PL + (itor - 1)*nDoF_tor_PL + (ifa - 1)*nDoF_tor_2D_PL + indices_tor_2D_PL
#else
            idx_PH = nDoF_tor_elem_PH + (itor - 1)*nDoF_tor_PH + (Mesh%loc2glob_fa(ifa) - 1)*nDoF_tor_2D_PH + indices_tor_2D_PH
            idx_PL = nDoF_tor_elem_PL + (itor - 1)*nDoF_tor_PL + (Mesh%loc2glob_fa(ifa) - 1)*nDoF_tor_2D_PL + indices_tor_2D_PL
#endif

            ! The -1 in the indices is because PETSc expects 0-based indexing for the global DOF indices, while Fortran arrays are 1-based.
            CALL MatSetValues(matPETSC%matR(level_counter), nDoF_tor_2D_PL, idx_PL, nDoF_tor_2D_PH, idx_PH, R_vals, INSERT_VALUES, ierr); CHKERRA(ierr)
            CALL MatSetValues(matPETSC%matP(level_counter), nDoF_tor_2D_PH, idx_PH, nDoF_tor_2D_PL, idx_PL, P_vals, INSERT_VALUES, ierr); CHKERRA(ierr)

         END DO
      END DO

      CALL MatAssemblyBegin(matPETSC%matR(level_counter), MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
      CALL MatAssemblyEnd(matPETSC%matR(level_counter), MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
      CALL MatAssemblyBegin(matPETSC%matP(level_counter), MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
      CALL MatAssemblyEnd(matPETSC%matP(level_counter), MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)

      DEALLOCATE (R_block)
      DEALLOCATE (P_block)
      DEALLOCATE (R_vals)
      DEALLOCATE (P_vals)
      DEALLOCATE (idx_PH)
      DEALLOCATE (idx_PL)

   END SUBROUTINE construct_global_operators_3D

   SUBROUTINE compute_jacobian_2D_tor(nq_pol_1D, nq_tor_1D, xydNorm, htor, xyg, J)
      USE LinearAlgebra, ONLY: col, TensorProduct
      INTEGER, INTENT(IN)          :: nq_pol_1D, nq_tor_1D
      REAL*8, INTENT(IN)           :: htor
      REAL*8, INTENT(IN)           :: xydNorm(nq_pol_1D)
      REAL*8, OPTIONAL, INTENT(IN) :: xyg(nq_pol_1D, 2)
      REAL*8, INTENT(OUT)          :: J(nq_pol_1D*nq_tor_1D)
      REAL*8                       :: ones(nq_tor_1D)

      ones = 1.0d0

      IF (PRESENT(xyg)) THEN
         J = col(tensorProduct(xydNorm*xyg(:, 1), ones))*(htor/2)
      ELSE
         J = col(tensorProduct(xydNorm, ones))*(htor/2)
      END IF

   END SUBROUTINE compute_jacobian_2D_tor

   PURE SUBROUTINE compute_jacobian_2D_pol(np, nq, dNdxi, dNdeta, xy, J)
      INTEGER, INTENT(IN)          :: np, nq
      REAL*8, INTENT(IN)           :: dNdxi(np, nq), dNdeta(np, nq), xy(np, 2)
      REAL*8, INTENT(OUT)          :: J(nq)
      REAL*8                       :: dx_dxi(2), dx_deta(2)
      INTEGER                      :: g

      J = 0.0d0

      DO g = 1, nq
         dx_dxi = MATMUL(dNdxi(:, g), xy)
         dx_deta = MATMUL(dNdeta(:, g), xy)
         J(g) = dx_dxi(1)*dx_deta(2) - dx_dxi(2)*dx_deta(1)
      END DO
   END SUBROUTINE compute_jacobian_2D_pol

   SUBROUTINE compute_toroidal_sf_and_weights_2D(np_pol_1D, nq_pol_1D, np_tor_1D, nq_tor_1D, N_pol_1D, N_tor_1D, wq_pol_1D, wq_tor_1D, N_tor_2D, wq_tor_2D)
      USE LinearAlgebra, ONLY: col, TensorProduct
      INTEGER, INTENT(IN)          :: np_pol_1D, nq_pol_1D, np_tor_1D, nq_tor_1D
      REAL*8, INTENT(IN)           :: N_pol_1D(np_pol_1D, nq_pol_1D), N_tor_1D(np_tor_1D, nq_tor_1D)
      REAL*8, OPTIONAL, INTENT(IN)           :: wq_pol_1D(nq_pol_1D), wq_tor_1D(nq_tor_1D)
      REAL*8, INTENT(OUT)          :: N_tor_2D(np_pol_1D*np_tor_1D, nq_pol_1D*nq_tor_1D)
      REAL*8, OPTIONAL, INTENT(OUT)          :: wq_tor_2D(nq_pol_1D*nq_tor_1D)

      INTEGER                      :: i, j, idx

      N_tor_2D = 0.0d0
      IF (PRESENT(wq_tor_2D)) wq_tor_2D = 0.0d0

      DO i = 1, np_tor_1D
         DO j = 1, np_pol_1D
            idx = (i - 1)*np_pol_1D + j
            N_tor_2D(idx, :) = col(tensorProduct(N_pol_1D(j, :), N_tor_1D(i, :)))
         END DO
      END DO

      IF (PRESENT(wq_tor_2D)) THEN
         DO i = 1, nq_tor_1D
            DO j = 1, nq_pol_1D
               idx = (i - 1)*nq_pol_1D + j
               wq_tor_2D(idx) = wq_pol_1D(j)*wq_tor_1D(i)
            END DO
         END DO
      END IF

   END SUBROUTINE compute_toroidal_sf_and_weights_2D

   PURE SUBROUTINE compute_sf_derivatives_at_phys_points(dNdxi, xy, xydNorm)
      IMPLICIT NONE
      REAL*8, INTENT(IN)           :: dNdxi(:, :)
      REAL*8, INTENT(IN)           :: xy(:, :)
      REAL*8, INTENT(OUT)          :: xydNorm(:)

      xydNorm = 0.0d0
      xydNorm = NORM2(MATMUL(TRANSPOSE(dNdxi), xy), 2)

   END SUBROUTINE compute_sf_derivatives_at_phys_points

   SUBROUTINE compute_P_sf_at_Gauss_points_1D(refEl_local, np, nq, xq, N, dNdxi)
      USE LinearAlgebra, ONLY: solve_linear_system_sing
      USE reference_element, ONLY: Vandermonde_1D, orthopoly1d, orthopoly1d_deriv
      TYPE(Reference_element_type), INTENT(IN)           :: refEl_local
      INTEGER, INTENT(IN)           :: np, nq
      REAL*8, INTENT(IN)            :: xq(nq)
      REAL*8, INTENT(OUT)           :: N(np, nq)
      REAL*8, OPTIONAL, INTENT(OUT) :: dNdxi(np, nq)
      REAL*8                        ::V(np, np), p(np), dp(np)
      INTEGER                       :: g

      N = 0.0d0
      IF (PRESENT(dNdxi)) dNdxi = 0.0d0

      ! Compute the P_og shape functions and derivatives at PH gauss points
      CALL Vandermonde_1D(V, refEl_local)
      V = TRANSPOSE(V)

      IF (PRESENT(dNdxi)) THEN
         DO g = 1, nq
            CALL orthopoly1d_deriv(xq(g), refEl_local%Ndeg, p, dp)
            CALL solve_linear_system_sing(V, p, N(:, g))
            CALL solve_linear_system_sing(V, dp, dNdxi(:, g))
         END DO
      ELSE
         DO g = 1, nq
            CALL orthopoly1D(xq(g), refEl_local%Ndeg, p)
            CALL solve_linear_system_sing(V, p, N(:, g))
         END DO
      END IF

   END SUBROUTINE compute_P_sf_at_Gauss_points_1D

   SUBROUTINE compute_P_sf_at_Gauss_points_2D(refEl_local, np_pol_2D_1, nq_pol_2D_2, xq_pol_2D_2, N, dNdxi, dNdeta)
      USE LinearAlgebra, ONLY: solve_linear_system_sing
      USE reference_element, ONLY: Vandermonde_2D, orthopoly2D_deriv
      TYPE(Reference_element_type), INTENT(IN) :: refEl_local
      INTEGER, INTENT(IN)          :: np_pol_2D_1, nq_pol_2D_2
      REAL*8, INTENT(IN)           :: xq_pol_2D_2(nq_pol_2D_2, 2)
      REAL*8, INTENT(OUT)             :: N(np_pol_2D_1, nq_pol_2D_2)
      REAL*8, OPTIONAL, INTENT(OUT)   :: dNdxi(np_pol_2D_1, nq_pol_2D_2), dNdeta(np_pol_2D_1, nq_pol_2D_2)
      REAL*8                       :: V(np_pol_2D_1, np_pol_2D_1), p(np_pol_2D_1), dpxi(np_pol_2D_1), dpeta(np_pol_2D_1)
      INTEGER                      :: g

      N = 0.0d0
      IF (PRESENT(dNdxi)) dNdxi = 0.0d0
      IF (PRESENT(dNdeta)) dNdeta = 0.0d0

      CALL Vandermonde_2D(V, refEl_local)
      V = TRANSPOSE(V)

      DO g = 1, nq_pol_2D_2
         CALL orthopoly2D_deriv(xq_pol_2D_2(g, 1), xq_pol_2D_2(g, 2), refEl_local%Ndeg, np_pol_2D_1, p, dpxi, dpeta)
         CALL solve_linear_system_sing(V, p, N(:, g))
         IF (PRESENT(dNdxi)) CALL solve_linear_system_sing(V, dpxi, dNdxi(:, g))
         IF (PRESENT(dNdeta)) CALL solve_linear_system_sing(V, dpeta, dNdeta(:, g))
      END DO

   END SUBROUTINE compute_P_sf_at_Gauss_points_2D
#endif

   ! Assemble mass matrix
   SUBROUTINE assemble_M_matrix(np1, np2, JJ, wq, N1, N2, M)
      INTEGER, INTENT(IN)          :: np1, np2
      REAL*8, INTENT(IN)           :: JJ(:), wq(:)
      REAL*8, INTENT(IN)           :: N1(:, :), N2(:, :)
      REAL*8, INTENT(OUT)          :: M(np1, np2)
      INTEGER                      :: i, j

      M = 0.0d0

      DO i = 1, np1
         DO j = 1, np2
            M(i, j) = SUM(JJ*wq*N1(i, :)*N2(j, :))
         END DO
      END DO

   END SUBROUTINE assemble_M_matrix

   PURE SUBROUTINE KRON(M, np1, np2, Neq, M_block)
      IMPLICIT NONE
      INTEGER, INTENT(IN)          :: np1, np2, Neq
      REAL*8, INTENT(IN)           :: M(np1, np2)
      REAL*8, INTENT(OUT)          :: M_block(np1*Neq, np2*Neq)
      INTEGER                      :: i, j, k

      M_block = 0.0D0

      DO i = 1, np1
         DO j = 1, np2
            DO k = 1, Neq
               M_block((i - 1)*Neq + k, (j - 1)*Neq + k) = M(i, j)
            END DO
         END DO
      END DO
   END SUBROUTINE KRON

   SUBROUTINE hdg_compute_permutation_mpiadj(A, theta, is_num, is_perm, ierr)
      Mat                          :: A, graph
      PetscReal                    :: theta
      IS                           :: is_num, is_perm, is_part, is_part_rows
      PetscErrorCode               :: ierr
      PetscInt, ALLOCATABLE        :: rowptr(:), adjacency(:), counts(:)
      PetscInt, ALLOCATABLE        :: part_rows(:)
      PetscInt, POINTER            :: part_nodes(:)
      PetscInt                     :: n_local, nn_local, neq, inode, ieq
      INTEGER                      :: rank, nranks, mpierr
      PetscInt, ALLOCATABLE        :: adjwgt(:)

      CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpierr)
      CALL MPI_Comm_size(PETSC_COMM_WORLD, nranks, mpierr)

      IF (rank .EQ. 0) WRITE (*,*) '----- PETSC IS COMPUTING MATRIX PERMUTATION -----'

      neq = phys%neq!
      CALL hdg_build_graph_mpiadj_blk(A, theta, neq, graph, rowptr, adjacency, ierr); CHKERRA(ierr)

      IF (rank .EQ. 0) WRITE (*,*) '----- PETSC IS PARTITIONING THE MATRIX -----'
      CALL hdg_partition_distributed(graph, is_part, ierr); CHKERRA(ierr)
      CALL MatDestroy(graph, ierr); CHKERRA(ierr)
      DEALLOCATE(rowptr, adjacency)!, adjwgt)

      ! ---- expand the node partition (nn_local entries) to a row partition (n_local entries)
      CALL ISGetLocalSize(is_part, nn_local, ierr); CHKERRA(ierr)
      n_local = nn_local*neq
      ALLOCATE(part_rows(n_local))
      CALL ISGetIndicesF90(is_part, part_nodes, ierr); CHKERRA(ierr)
      DO inode = 1, nn_local
         DO ieq = 1, neq
            part_rows((inode - 1)*neq + ieq) = part_nodes(inode)
         END DO
      END DO
      CALL ISRestoreIndicesF90(is_part, part_nodes, ierr); CHKERRA(ierr)
      CALL ISDestroy(is_part, ierr); CHKERRA(ierr)

      CALL ISCreateGeneral(PETSC_COMM_WORLD, n_local, part_rows, PETSC_COPY_VALUES, &
         is_part_rows, ierr); CHKERRA(ierr)
      DEALLOCATE(part_rows)

      ! ---- from here on everything is per row, as before
      ALLOCATE(counts(nranks))
      CALL ISPartitioningCount(is_part_rows, nranks, counts, ierr); CHKERRA(ierr)
      n_local = counts(rank + 1)          ! new number of rows owned by this rank
      DEALLOCATE(counts)

      CALL ISPartitioningToNumbering(is_part_rows, is_num, ierr); CHKERRA(ierr)
      CALL ISInvertPermutation(is_num, n_local, is_perm, ierr); CHKERRA(ierr)
      CALL ISDestroy(is_part_rows, ierr); CHKERRA(ierr)

   END SUBROUTINE hdg_compute_permutation_mpiadj

   SUBROUTINE hdg_build_graph_mpiadj_blk(A, theta, neq, graph, rowptr, adjacency, ierr)
      Mat                          :: A, graph
      PetscReal                    :: theta
      PetscInt                     :: neq
      PetscInt, ALLOCATABLE, INTENT(OUT) :: rowptr(:), adjacency(:)
      PetscErrorCode               :: ierr

      PetscInt                     :: n_global, row_start, row_end, n_local
      PetscInt                     :: nn_global, nn_local, node_start, node_end, inode
      PetscInt                     :: total_send, total_recv, target_row, col_index, owner
      PetscInt                     :: pos, write_pos, row_nnz, tmp_col, k, j, row, nkept, row_start_dummy
      PetscInt, ALLOCATABLE        :: old_rowptr(:), row_counts(:), row_cursor(:)
      PetscInt, ALLOCATABLE        :: send_pairs(:), recv_pairs(:), send_displs(:), recv_displs(:)
      PetscInt, ALLOCATABLE        :: owner_counts(:), owner_displs(:), rank_starts(:)
      PetscInt, ALLOCATABLE        :: send_counts(:), recv_counts(:)
      ! work arrays for the per-node neighbour search (host-associated)
      PetscInt, ALLOCATABLE        :: columns_row(:), blist(:), kept(:)
      PetscScalar, ALLOCATABLE     :: values_row(:)
      PetscReal, ALLOCATABLE       :: slist(:)
      INTEGER                      :: rank, nranks, mpierr, r
      INTEGER, ALLOCATABLE         :: send_counts2(:), recv_counts2(:), send_displs2(:), recv_displs2(:)

      IF (mpivar%glob_id .EQ. 0) WRITE (*,*) '----- PETSC IS COMPUTING THE BLOCK GRAPH -----'

      CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, mpierr)
      CALL MPI_Comm_size(PETSC_COMM_WORLD, nranks, mpierr)
      CALL MatGetSize(A, n_global, n_global, ierr); CHKERRA(ierr)
      CALL MatGetOwnershipRange(A, row_start, row_end, ierr); CHKERRA(ierr)
      n_local = row_end - row_start

      IF (MOD(n_local, neq) .NE. 0 .OR. MOD(n_global, neq) .NE. 0) THEN
         WRITE (*,*) 'hdg_build_graph_mpiadj_blk: ownership range not a multiple of neq'
         ierr = 1
         RETURN
      END IF
      nn_local  = n_local/neq
      nn_global = n_global/neq

      ! node ownership
      ALLOCATE(owner_counts(nranks), owner_displs(nranks), rank_starts(nranks + 1))
      CALL MPI_Allgather(nn_local, 1, MPIU_INTEGER, owner_counts, 1, MPIU_INTEGER, &
         PETSC_COMM_WORLD, mpierr)
      owner_displs(1) = 0
      DO r = 2, nranks
         owner_displs(r) = owner_displs(r - 1) + owner_counts(r - 1)
      END DO
      rank_starts(1:nranks) = owner_displs
      rank_starts(nranks + 1) = nn_global
      node_start = rank_starts(rank + 1)
      node_end   = rank_starts(rank + 2)

      ALLOCATE(row_counts(nn_local), send_counts(nranks))
      row_counts = 0
      send_counts = 0
      ALLOCATE(columns_row(MAX(1, n_local)), values_row(MAX(1, n_local)))
      ALLOCATE(blist(64), slist(64), kept(64))

      ! ---- pass 1: count kept edges per node and reverse edges per owner
      DO inode = node_start, node_end - 1
         CALL node_neighbors(inode, nkept)
         row_counts(inode - node_start + 1) = nkept
         DO k = 1, nkept
            CALL find_owner(kept(k), owner)
            send_counts(owner) = send_counts(owner) + 1
         END DO
      END DO

      ALLOCATE(recv_counts(nranks), send_displs(nranks), recv_displs(nranks))
      CALL MPI_Alltoall(send_counts, 1, MPIU_INTEGER, recv_counts, 1, MPIU_INTEGER, &
         PETSC_COMM_WORLD, mpierr)
      send_displs(1) = 0
      recv_displs(1) = 0
      DO r = 2, nranks
         send_displs(r) = send_displs(r - 1) + send_counts(r - 1)
         recv_displs(r) = recv_displs(r - 1) + recv_counts(r - 1)
      END DO
      total_send = SUM(send_counts)
      total_recv = SUM(recv_counts)
      ALLOCATE(send_pairs(2*total_send), recv_pairs(2*total_recv))
      ALLOCATE(send_counts2(nranks), recv_counts2(nranks), send_displs2(nranks), recv_displs2(nranks))
      send_counts2 = 2*send_counts
      recv_counts2 = 2*recv_counts
      send_displs2 = 2*send_displs
      recv_displs2 = 2*recv_displs

      ! ---- pass 2: pack reversed edges as (destination node, source node)
      send_cursor: BLOCK
         INTEGER, ALLOCATABLE      :: cursor(:)
         ALLOCATE(cursor(nranks))
         cursor = send_displs
         DO inode = node_start, node_end - 1
            CALL node_neighbors(inode, nkept)
            DO k = 1, nkept
               CALL find_owner(kept(k), owner)
               pos = cursor(owner)*2
               send_pairs(pos + 1) = kept(k)
               send_pairs(pos + 2) = inode
               cursor(owner) = cursor(owner) + 1
            END DO
         END DO
         DEALLOCATE(cursor)
      END BLOCK send_cursor

      CALL MPI_Alltoallv(send_pairs, send_counts2, send_displs2, MPIU_INTEGER, &
         recv_pairs, recv_counts2, recv_displs2, MPIU_INTEGER, PETSC_COMM_WORLD, mpierr)

      DO k = 1, total_recv
         target_row = recv_pairs(2*k - 1)
         row_counts(target_row - node_start + 1) = row_counts(target_row - node_start + 1) + 1
      END DO

      ALLOCATE(rowptr(nn_local + 1), old_rowptr(nn_local + 1), row_cursor(nn_local))
      rowptr(1) = 0
      DO k = 1, nn_local
         rowptr(k + 1) = rowptr(k) + row_counts(k)
      END DO
      old_rowptr = rowptr
      ALLOCATE(adjacency(rowptr(nn_local + 1)))
      row_cursor = rowptr(1:nn_local) + 1

      ! ---- pass 3: fill local directed edges
      DO inode = node_start, node_end - 1
         CALL node_neighbors(inode, nkept)
         DO k = 1, nkept
            pos = row_cursor(inode - node_start + 1)
            adjacency(pos) = kept(k)
            row_cursor(inode - node_start + 1) = pos + 1
         END DO
      END DO

      ! append reversed edges
      DO k = 1, total_recv
         target_row = recv_pairs(2*k - 1)
         col_index  = recv_pairs(2*k)
         pos = row_cursor(target_row - node_start + 1)
         adjacency(pos) = col_index
         row_cursor(target_row - node_start + 1) = pos + 1
      END DO

      ! ---- sort + dedupe (unchanged)
      rowptr(1) = 0
      DO row = 1, nn_local
         pos = old_rowptr(row) + 1
         row_nnz = old_rowptr(row + 1) - old_rowptr(row)
         DO j = pos + 1, pos + row_nnz - 1
            tmp_col = adjacency(j)
            k = j - 1
            DO WHILE (k .GE. pos)
               IF (adjacency(k) .LE. tmp_col) EXIT
               adjacency(k + 1) = adjacency(k)
               k = k - 1
            END DO
            adjacency(k + 1) = tmp_col
         END DO

         write_pos = rowptr(row) + 1
         DO j = pos, pos + row_nnz - 1
            IF (j .EQ. pos) THEN
               adjacency(write_pos) = adjacency(j)
               write_pos = write_pos + 1
            ELSEIF (adjacency(j) .NE. adjacency(j - 1)) THEN
               adjacency(write_pos) = adjacency(j)
               write_pos = write_pos + 1
            END IF
         END DO
         rowptr(row + 1) = write_pos - 1
      END DO

      CALL MatCreateMPIAdj(PETSC_COMM_WORLD, nn_local, nn_global, rowptr, adjacency, PETSC_NULL_INTEGER, graph, ierr); CHKERRA(ierr)

      DEALLOCATE(owner_counts, owner_displs, rank_starts, row_counts, row_cursor)
      DEALLOCATE(old_rowptr, send_pairs, recv_pairs)
      DEALLOCATE(send_counts, recv_counts, send_displs, recv_displs)
      DEALLOCATE(send_counts2, recv_counts2, send_displs2, recv_displs2)
      DEALLOCATE(columns_row, values_row, blist, slist, kept)

   CONTAINS

      SUBROUTINE find_owner(node, own)
         PetscInt, INTENT(IN)      :: node
         PetscInt, INTENT(OUT)     :: own
         DO own = 1, nranks
            IF (node .GE. rank_starts(own) .AND. node .LT. rank_starts(own + 1)) EXIT
         END DO
      END SUBROUTINE find_owner

      !Strong off-diagonal block neighbours of node `nd` -> kept(1:nk)
      SUBROUTINE node_neighbors(nd, nk)
         PetscInt, INTENT(IN)      :: nd
         PetscInt, INTENT(OUT)     :: nk
         PetscInt                  :: ii, jj, kk, nl, ncols_l, bc
         PetscReal                 :: aval, smax, thr          ! was: a
         PetscInt, ALLOCATABLE     :: itmp(:)
         PetscReal, ALLOCATABLE    :: rtmp(:)
         LOGICAL                   :: found

         nl = 0
         DO ii = 0, neq - 1
            CALL MatGetRow(A, nd*neq + ii, ncols_l, columns_row, values_row, ierr); CHKERRA(ierr)
            DO jj = 1, ncols_l
               bc = columns_row(jj)/neq
               IF (bc .EQ. nd) CYCLE
               aval = ABS(values_row(jj))        ! was: a = ...
               found = .FALSE.
               DO kk = nl, 1, -1
                  IF (blist(kk) .EQ. bc) THEN
                     slist(kk) = MAX(slist(kk), aval)
                     found = .TRUE.
                     EXIT
                  END IF
               END DO
               IF (.NOT. found) THEN
                  IF (nl .EQ. SIZE(blist)) THEN
                     ALLOCATE(itmp(2*nl), rtmp(2*nl))
                     itmp(1:nl) = blist(1:nl); CALL MOVE_ALLOC(itmp, blist)
                     rtmp(1:nl) = slist(1:nl); CALL MOVE_ALLOC(rtmp, slist)
                     DEALLOCATE(kept); ALLOCATE(kept(2*nl))
                  END IF
                  nl = nl + 1
                  blist(nl) = bc
                  slist(nl) = aval
               END IF
            END DO
            CALL MatRestoreRow(A, nd*neq + ii, ncols_l, columns_row, values_row, ierr); CHKERRA(ierr)
         END DO

         nk = 0
         IF (nl .EQ. 0) RETURN
         smax = MAXVAL(slist(1:nl))
         thr  = theta*smax
         DO kk = 1, nl
            IF (slist(kk) .GT. 0.0 .AND. slist(kk) .GE. thr) THEN
               nk = nk + 1
               kept(nk) = blist(kk)
            END IF
         END DO
      END SUBROUTINE node_neighbors

   END SUBROUTINE hdg_build_graph_mpiadj_blk

   SUBROUTINE hdg_partition_distributed(G, is_part, ierr)
      Mat                          :: G
      IS                           :: is_part
      PetscErrorCode               :: ierr

      MatPartitioning              :: part
      CALL MatPartitioningCreate(PETSC_COMM_WORLD, part, ierr); CHKERRA(ierr)
      CALL MatPartitioningSetAdjacency(part, G, ierr); CHKERRA(ierr)
      !call MatPartitioningSetType(part, 'parmetis', ierr); CHKERRA(ierr)
      CALL MatPartitioningSetType(part, 'ptscotch', ierr); CHKERRA(ierr)
      CALL MatPartitioningSetFromOptions(part, ierr); CHKERRA(ierr)
      CALL MatPartitioningApply(part, is_part, ierr); CHKERRA(ierr)
      CALL MatPartitioningDestroy(part, ierr); CHKERRA(ierr)
   END SUBROUTINE hdg_partition_distributed

   SUBROUTINE hdg_partition_seq(G, is_part, ierr)
      Mat                          :: G
      IS                           :: is_part
      PetscErrorCode               :: ierr

      Mat, POINTER                 :: Gseq(:)
      IS                           :: isrow(1), iscol(1), is_seq
      MatPartitioning              :: part
      PetscMPIInt                  :: rank, nranks, nmpi
      PetscInt                     :: N, M, rstart, rend, ntake, ione, izero, edgecut
      PetscInt, POINTER            :: pp(:)
      PetscInt, ALLOCATABLE        :: partall(:)
      CHARACTER(LEN=16)            :: strn

      ione = 1
      izero = 0
      NULLIFY (Gseq)
      ALLOCATE(Gseq(ione))
      CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
      CALL MPI_Comm_size(PETSC_COMM_WORLD, nranks, ierr)
      CALL MatGetSize(G, N, M, ierr)
      CALL MatGetOwnershipRange(G, rstart, rend, ierr)

      ! rank 0 takes the whole graph, the other ranks take nothing
      ntake = 0
      IF (rank .EQ. 0) ntake = N
      CALL ISCreateStride(PETSC_COMM_SELF, ntake, izero, ione, isrow(1), ierr); CHKERRA(ierr)
      CALL ISCreateStride(PETSC_COMM_SELF, ntake, izero, ione, iscol(1), ierr); CHKERRA(ierr)
      CALL MatCreateSubMatrices(G, ione, isrow, iscol, MAT_INITIAL_MATRIX, Gseq, ierr); CHKERRA(ierr)

      ALLOCATE (partall(N))
      partall = 0

      IF (rank .EQ. 0) THEN
         CALL MatPartitioningCreate(PETSC_COMM_SELF, part, ierr); CHKERRA(ierr)
         CALL MatPartitioningSetAdjacency(part, Gseq(1), ierr); CHKERRA(ierr)
         !call MatPartitioningSetType(part, 'scotch', ierr); CHKERRA(ierr)
         CALL MatPartitioningSetType(part, 'parmetis', ierr); CHKERRA(ierr)
         WRITE (strn, '(I0)') MPIvar%glob_size
         CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mat_partitioning_nparts', TRIM(ADJUSTL(strn)), ierr); CHKERRA(ierr)
         CALL MatPartitioningSetFromOptions(part, ierr); CHKERRA(ierr)
         CALL PetscOptionsClearValue(PETSC_NULL_OPTIONS, '-mat_partitioning_nparts', ierr); CHKERRA(ierr)
         CALL MatPartitioningApply(part, is_seq, ierr); CHKERRA(ierr)
         CALL ISGetIndicesF90(is_seq, pp, ierr); CHKERRA(ierr)
         partall(1:N) = pp(1:N)
         CALL ISRestoreIndicesF90(is_seq, pp, ierr); CHKERRA(ierr)
         CALL ISDestroy(is_seq, ierr); CHKERRA(ierr)
         CALL MatPartitioningDestroy(part, ierr); CHKERRA(ierr)
      END IF

      CALL MatDestroySubMatrices(ione, Gseq, ierr); CHKERRA(ierr)
      CALL ISDestroy(isrow(1), ierr); CHKERRA(ierr)
      CALL ISDestroy(iscol(1), ierr); CHKERRA(ierr)

      nmpi = N
      CALL MPI_Bcast(partall, nmpi, MPIU_INTEGER, 0, PETSC_COMM_WORLD, ierr)

      ! target rank of each locally owned row, same as MatPartitioningApply returns
      CALL ISCreateGeneral(PETSC_COMM_WORLD, rend - rstart, partall(rstart+1:rend), PETSC_COPY_VALUES, is_part, ierr); CHKERRA(ierr)
      DEALLOCATE (partall)
   END SUBROUTINE hdg_partition_seq

   SUBROUTINE hdg_apply_permutation(A, is_perm, Aperm, reuse, ierr)
      Mat                          :: A, Aperm
      IS                           :: is_perm
      PetscErrorCode               :: ierr
      MatReuse                     :: reuse

      CALL MatCreateSubMatrix(A, is_perm, is_perm, reuse, Aperm, ierr); CHKERRA(ierr)

   END SUBROUTINE hdg_apply_permutation

   SUBROUTINE hdg_build_graph_fast(A, theta, G, ierr)
      Mat A, G, P, M
      PetscReal theta
      PetscErrorCode ierr
      PetscInt rstart, rend, nloc, i, j, k, ncols, nkeep, ione, maxrow, ipass
      PetscInt irow(1)
      PetscInt, ALLOCATABLE        :: cols(:), jkeep(:)
      PetscScalar, ALLOCATABLE     :: vals(:), vkeep(:)
      PetscReal amax, cut
      PetscScalar one

      IF (MPIvar%glob_id .EQ. 0) THEN
         PRINT *, 'Building non rescaled graph...'
      END IF

      ione = 1
      one = 1.0
      CALL MatGetOwnershipRange(A, rstart, rend, ierr); CHKERRA(ierr)
      nloc = rend - rstart

      ! ---- pass 1: longest row, only ncols is requested (no copy) ----
      maxrow = 1
      DO i = rstart, rend - 1
         CALL MatGetRow(A, i, ncols, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr); CHKERRA(ierr)
         maxrow = max(maxrow, ncols)
         CALL MatRestoreRow(A, i, ncols, PETSC_NULL_INTEGER, PETSC_NULL_SCALAR, ierr); CHKERRA(ierr)
      END DO
      ALLOCATE (cols(maxrow), vals(maxrow), jkeep(maxrow), vkeep(maxrow))

      ! ---- preallocator matrix ----
      CALL MatCreate(PETSC_COMM_WORLD, P, ierr); CHKERRA(ierr)
      CALL MatSetType(P, MATPREALLOCATOR, ierr); CHKERRA(ierr)
      CALL MatSetSizes(P, nloc, nloc, PETSC_DETERMINE, PETSC_DETERMINE, ierr); CHKERRA(ierr)
      CALL MatSetUp(P, ierr); CHKERRA(ierr)

      ! ---- real matrix: type and sizes only, preallocation comes from P ----
      CALL MatCreate(PETSC_COMM_WORLD, G, ierr); CHKERRA(ierr)
      CALL MatSetType(G, MATAIJ, ierr); CHKERRA(ierr)
      CALL MatSetSizes(G, nloc, nloc, PETSC_DETERMINE, PETSC_DETERMINE, ierr); CHKERRA(ierr)

      ! ---- two passes: ipass=1 fills P, ipass=2 fills G ----
      DO ipass = 1, 2
         IF (ipass .EQ. 1) THEN
            M = P
         ELSE
            M = G
         END IF

         DO i = rstart, rend - 1
            CALL MatGetRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
            amax = 0.0
            DO j = 1, ncols
               IF (cols(j) .NE. i) amax = max(amax, abs(vals(j)))
            END DO
            cut = theta*amax
            nkeep = 0
            DO j = 1, ncols
               IF (cols(j) .NE. i .AND. abs(vals(j)) .GE. cut .AND. abs(vals(j)) .GT. 0.0) THEN
                  nkeep = nkeep + 1
                  jkeep(nkeep) = cols(j)
                  vkeep(nkeep) = one
               END IF
            END DO
            CALL MatRestoreRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)

            IF (nkeep .GT. 0) THEN
               irow(1) = i
               ! entries (i, jkeep(k))
               CALL MatSetValues(M, ione, irow, nkeep, jkeep, vkeep, INSERT_VALUES, ierr); CHKERRA(ierr)
               ! transposed entries (jkeep(k), i): nkeep rows, one column
               CALL MatSetValues(M, nkeep, jkeep, ione, irow, vkeep, INSERT_VALUES, ierr); CHKERRA(ierr)
            END IF
         END DO

         CALL MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
         CALL MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)

         IF (ipass .EQ. 1) THEN
            ! exact diag/off-diag counts from P into G, do not fill
            CALL MatPreallocatorPreallocate(P, PETSC_FALSE, G, ierr); CHKERRA(ierr)
            CALL MatDestroy(P, ierr); CHKERRA(ierr)
         END IF
      END DO

      CALL MatSetOption(G, MAT_SYMMETRIC, PETSC_TRUE, ierr); CHKERRA(ierr)
      DEALLOCATE (cols, vals, jkeep, vkeep)
   END SUBROUTINE hdg_build_graph_fast

   SUBROUTINE hdg_build_graph_scaled(A, theta, G, ierr)
      Mat A, G, P, M
      PetscReal theta
      PetscErrorCode ierr

      PetscInt rstart, rend, nloc, i, j, k, ncols, nkeep, ione, maxrow, ipass
      PetscInt nghost, ntmp, c, pos
      PetscInt irow(1)
      PetscInt, ALLOCATABLE        :: cols(:), jkeep(:), gcols(:), idxto(:)
      PetscScalar, ALLOCATABLE     :: vals(:), vkeep(:)
      PetscReal, ALLOCATABLE       :: w(:)
      PetscReal amax, rmax, diag, wmax, cut, dj, di, tiny_rel
      PetscScalar one
      PetscScalar, POINTER         :: dl(:), dg(:)
      Vec dvec, gvec
      IS isfrom, isto
      VecScatter scat

      ione = 1
      one = 1.0
      tiny_rel = 1.0e-12

      IF (MPIvar%glob_id .EQ. 0) THEN
         PRINT *, 'Building scaled graph...'
      END IF

      CALL MatGetOwnershipRange(A, rstart, rend, ierr); CHKERRA(ierr)
      nloc = rend - rstart

      ! ---------------------------------------------------------------
      ! pass 1: longest row, d_i = sqrt(|a_ii|), count off-process columns
      ! ---------------------------------------------------------------
      CALL VecCreateMPI(PETSC_COMM_WORLD, nloc, PETSC_DETERMINE, dvec, ierr); CHKERRA(ierr)

      maxrow = 1
      nghost = 0
      DO i = rstart, rend - 1
         CALL MatGetRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
         maxrow = max(maxrow, ncols)
         diag = 0.0
         rmax = 0.0
         DO j = 1, ncols
            IF (cols(j) .EQ. i) THEN
               diag = abs(vals(j))
            ELSE
               rmax = max(rmax, abs(vals(j)))
               IF (cols(j) .LT. rstart .OR. cols(j) .GE. rend) nghost = nghost + 1
            END IF
         END DO
         ! fallback for zero/tiny diagonal (saddle-point-like rows, weak diagonal)
         IF (diag .LE. tiny_rel*rmax) diag = rmax
         IF (diag .EQ. 0.0) diag = 1.0
         CALL VecSetValue(dvec, i, sqrt(diag), INSERT_VALUES, ierr); CHKERRA(ierr)
         CALL MatRestoreRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
      END DO
      CALL VecAssemblyBegin(dvec, ierr); CHKERRA(ierr)
      CALL VecAssemblyEnd(dvec, ierr); CHKERRA(ierr)

      ! ---------------------------------------------------------------
      ! collect distinct ghost columns, sorted
      ! ---------------------------------------------------------------
      ALLOCATE (gcols(max(nghost, 1)))
      k = 0
      DO i = rstart, rend - 1
         CALL MatGetRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
         DO j = 1, ncols
            IF (cols(j) .LT. rstart .OR. cols(j) .GE. rend) THEN
               k = k + 1
               gcols(k) = cols(j)
            END IF
         END DO
         CALL MatRestoreRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
      END DO
      ntmp = k
      CALL PetscSortRemoveDupsInt(ntmp, gcols, ierr); CHKERRA(ierr)
      nghost = ntmp

      ! ---------------------------------------------------------------
      ! scatter d_j of ghost columns into a small sequential vector
      ! ---------------------------------------------------------------
      ALLOCATE (idxto(max(nghost, 1)))
      DO k = 1, nghost
         idxto(k) = k - 1
      END DO
      CALL VecCreateSeq(PETSC_COMM_SELF, nghost, gvec, ierr); CHKERRA(ierr)
      CALL ISCreateGeneral(PETSC_COMM_SELF, nghost, gcols, PETSC_COPY_VALUES, isfrom, ierr); CHKERRA(ierr)
      CALL ISCreateGeneral(PETSC_COMM_SELF, nghost, idxto, PETSC_COPY_VALUES, isto, ierr); CHKERRA(ierr)
      CALL VecScatterCreate(dvec, isfrom, gvec, isto, scat, ierr); CHKERRA(ierr)
      CALL VecScatterBegin(scat, dvec, gvec, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRA(ierr)
      CALL VecScatterEnd(scat, dvec, gvec, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRA(ierr)
      CALL VecScatterDestroy(scat, ierr); CHKERRA(ierr)
      CALL ISDestroy(isfrom, ierr); CHKERRA(ierr)
      CALL ISDestroy(isto, ierr); CHKERRA(ierr)
      DEALLOCATE (idxto)

      CALL VecGetArrayRead(dvec, dl, ierr); CHKERRA(ierr)
      CALL VecGetArrayRead(gvec, dg, ierr); CHKERRA(ierr)

      ALLOCATE (cols(maxrow), vals(maxrow), jkeep(maxrow), vkeep(maxrow), w(maxrow))

      ! ---------------------------------------------------------------
      ! preallocator and real matrix
      ! ---------------------------------------------------------------
      CALL MatCreate(PETSC_COMM_WORLD, P, ierr); CHKERRA(ierr)
      CALL MatSetType(P, MATPREALLOCATOR, ierr); CHKERRA(ierr)
      CALL MatSetSizes(P, nloc, nloc, PETSC_DETERMINE, PETSC_DETERMINE, ierr); CHKERRA(ierr)
      CALL MatSetUp(P, ierr); CHKERRA(ierr)

      CALL MatCreate(PETSC_COMM_WORLD, G, ierr); CHKERRA(ierr)
      CALL MatSetType(G, MATAIJ, ierr); CHKERRA(ierr)
      CALL MatSetSizes(G, nloc, nloc, PETSC_DETERMINE, PETSC_DETERMINE, ierr); CHKERRA(ierr)

      DO ipass = 1, 2
         IF (ipass .EQ. 1) THEN
            M = P
         ELSE
            M = G
         END IF

         DO i = rstart, rend - 1
            CALL MatGetRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
            di = dl(i - rstart + 1)

            ! normalized weights w_ij = |a_ij| / (d_i d_j)
            wmax = 0.0
            DO j = 1, ncols
               c = cols(j)
               IF (c .EQ. i) THEN
                  w(j) = 0.0
                  CYCLE
               END IF
               IF (c .GE. rstart .AND. c .LT. rend) THEN
                  dj = dl(c - rstart + 1)
               ELSE
                  pos = bsearch(gcols, nghost, c)
                  dj = dg(pos)
               END IF
               w(j) = abs(vals(j))/(di*dj)
               wmax = max(wmax, w(j))
            END DO
            cut = theta*wmax

            nkeep = 0
            DO j = 1, ncols
               IF (cols(j) .NE. i .AND. w(j) .GE. cut .AND. w(j) .GT. 0.0) THEN
                  nkeep = nkeep + 1
                  jkeep(nkeep) = cols(j)
                  vkeep(nkeep) = one
               END IF
            END DO
            CALL MatRestoreRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)

            IF (nkeep .GT. 0) THEN
               irow(1) = i
               CALL MatSetValues(M, ione, irow, nkeep, jkeep, vkeep, INSERT_VALUES, ierr); CHKERRA(ierr)
               CALL MatSetValues(M, nkeep, jkeep, ione, irow, vkeep, INSERT_VALUES, ierr); CHKERRA(ierr)
            END IF
         END DO

         CALL MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
         CALL MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)

         IF (ipass .EQ. 1) THEN
            CALL MatPreallocatorPreallocate(P, PETSC_FALSE, G, ierr); CHKERRA(ierr)
            CALL MatDestroy(P, ierr); CHKERRA(ierr)
         END IF
      END DO

      CALL MatSetOption(G, MAT_SYMMETRIC, PETSC_TRUE, ierr); CHKERRA(ierr)

      CALL VecRestoreArrayRead(dvec, dl, ierr); CHKERRA(ierr)
      CALL VecRestoreArrayRead(gvec, dg, ierr); CHKERRA(ierr)
      CALL VecDestroy(dvec, ierr); CHKERRA(ierr)
      CALL VecDestroy(gvec, ierr); CHKERRA(ierr)
      DEALLOCATE (cols, vals, jkeep, vkeep, w, gcols)

   CONTAINS

      ! binary search of key in sorted array a(1:n); returns 1-based position
      INTEGER FUNCTION bsearch(a, n, key)
         PetscInt, INTENT(IN)      :: a(:), n, key
         PetscInt lo, hi, mid
         lo = 1
         hi = n
         bsearch = 1
         DO WHILE (lo .LE. hi)
            mid = (lo + hi)/2
            IF (a(mid) .EQ. key) THEN
               bsearch = mid
               RETURN
            ELSE IF (a(mid) .LT. key) THEN
               lo = mid + 1
            ELSE
               hi = mid - 1
            END IF
         END DO
      END FUNCTION bsearch

   END SUBROUTINE hdg_build_graph_scaled

   SUBROUTINE hdg_build_graph(A, theta, G, ierr)
      Mat A, G, Gt
      PetscReal theta
      PetscErrorCode ierr
      PetscInt rstart, rend, nloc, i, j, ncols, nkeep, nzmax, ione, irows
      PetscInt Nglob, Mglob, nbound, maxrow
      PetscInt, ALLOCATABLE        :: cols(:), jkeep(:)
      PetscScalar, ALLOCATABLE     :: vals(:), vkeep(:)
      PetscReal amax, cut
      PetscScalar one
      DOUBLE PRECISION info(MAT_INFO_SIZE)

      ione = 1
      one = 1.0
      CALL MatGetOwnershipRange(A, rstart, rend, ierr); CHKERRA(ierr)
      nloc = rend - rstart

      ! upper bound for the length of any row
      CALL MatGetSize(A, Mglob, Nglob, ierr); CHKERRA(ierr)
      CALL MatGetInfo(A, MAT_LOCAL, info, ierr); CHKERRA(ierr)
      nbound = min(Nglob, nint(info(MAT_INFO_NZ_USED)))
      nbound = max(nbound, 1)
      ALLOCATE (cols(nbound), vals(nbound), jkeep(nbound), vkeep(nbound))

#ifdef TOR3D
      ! pass 1: longest actual row (for preallocation of G)
      maxrow = 1
      DO i = rstart, rend-1
         CALL MatGetRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
         maxrow = max(maxrow, ncols)
         CALL MatRestoreRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
      END DO
#else
      ! the rows for internal faces (the ones with the most nonzeros) give number of nodes per face * 5 coupled traces * number of equations per node
      maxrow = (refElPol%Ndeg + 1)*5*phys%Neq
#endif

      nzmax = maxrow

      CALL MatCreate(PETSC_COMM_WORLD, G, ierr); CHKERRA(ierr)
      CALL MatSetSizes(G, nloc, nloc, PETSC_DETERMINE, PETSC_DETERMINE, ierr); CHKERRA(ierr)
      CALL MatSetType(G, MATAIJ, ierr); CHKERRA(ierr)
      CALL MatSeqAIJSetPreallocation(G, nzmax, PETSC_NULL_INTEGER, ierr); CHKERRA(ierr)
      CALL MatMPIAIJSetPreallocation(G, nzmax, PETSC_NULL_INTEGER, nzmax, PETSC_NULL_INTEGER, ierr); CHKERRA(ierr)
      CALL MatSetOption(G, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRA(ierr)

      ! ! pass 2: keep off-diagonal entries with |a_ij| >= theta*max|a_ik|
      DO i = rstart, rend-1
         CALL MatGetRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
         amax = 0.0
         DO j = 1, ncols
            IF (cols(j) .NE. i) amax = max(amax, abs(vals(j)))
         END DO
         cut = theta*amax
         nkeep = 0
         DO j = 1, ncols
            IF (cols(j) .NE. i .AND. abs(vals(j)) .GE. cut .AND. abs(vals(j)) .GT. 0.0) THEN
               nkeep = nkeep + 1
               jkeep(nkeep) = cols(j)
               vkeep(nkeep) = one
            END IF
         END DO

         IF (nkeep .GT. 0) THEN
            CALL MatSetValues(G, ione, i, nkeep, jkeep(1:nkeep), vkeep(1:nkeep), INSERT_VALUES, ierr); CHKERRA(ierr)
         END IF
         CALL MatRestoreRow(A, i, ncols, cols, vals, ierr); CHKERRA(ierr)
      END DO

      DEALLOCATE (cols, vals, jkeep, vkeep)

      CALL MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)
      CALL MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY, ierr); CHKERRA(ierr)

      CALL MatTranspose(G, MAT_INITIAL_MATRIX, Gt, ierr); CHKERRA(ierr)
      CALL MatAXPY(G, one, Gt, DIFFERENT_NONZERO_PATTERN, ierr); CHKERRA(ierr)
      CALL MatDestroy(Gt, ierr); CHKERRA(ierr)
   END SUBROUTINE hdg_build_graph

   !---------------------------------------------------------------------------
   ! Scatter between the old layout (vold) and the new one (vnew):
   !   forward : vnew(is_num(i)) = vold(i)   (rhs  -> permuted rhs)
   !   reverse : vold(i) = vnew(is_num(i))   (permuted solution -> original)
   ! usage:
   !   call VecScatterBegin(sc, b,    bperm, INSERT_VALUES, SCATTER_FORWARD, ierr)
   !   call VecScatterEnd  (sc, b,    bperm, INSERT_VALUES, SCATTER_FORWARD, ierr)
   !   ... KSPSolve(ksp, bperm, xperm) ...
   !   call VecScatterBegin(sc, x, xperm, INSERT_VALUES, SCATTER_REVERSE, ierr)
   !   call VecScatterEnd  (sc, x, xperm, INSERT_VALUES, SCATTER_REVERSE, ierr)
   !   (reverse: the "to" vector of the call is the one written, i.e. x = old layout;
   !    argument order is always (sc, vold, vnew, ...))
   !---------------------------------------------------------------------------
   SUBROUTINE hdg_create_vec_scatter(is_num, vold, vnew, sc, ierr)
      IS                           :: is_num
      Vec                          :: vold, vnew
      VecScatter                   :: sc
      PetscErrorCode               :: ierr

      MPI_Comm                     :: comm
      IS                           :: is_from
      PetscInt                     :: rstart, rend, ione

      ione  = 1

      comm = PETSC_COMM_WORLD
      CALL VecGetOwnershipRange(vold, rstart, rend, ierr); CHKERRA(ierr)
      CALL ISCreateStride(comm, rend - rstart, rstart, ione, is_from, ierr); CHKERRA(ierr)
      CALL VecScatterCreate(vold, is_from, vnew, is_num, sc, ierr); CHKERRA(ierr)
      CALL ISDestroy(is_from, ierr); CHKERRA(ierr)

   END SUBROUTINE hdg_create_vec_scatter

   SUBROUTINE terminate_PETSC()
      PetscErrorCode               :: ierr
      INTEGER                      :: i
      IF (ASSOCIATED(matPETSC%rowptr)) THEN
         DEALLOCATE (matPETSC%rowptr)
      END IF
      IF (ASSOCIATED(matPETSC%loc2glob)) THEN
         DEALLOCATE (matPETSC%loc2glob)
      END IF
      IF (ASSOCIATED(matPETSC%cols)) THEN
         DEALLOCATE (matPETSC%cols)
      END IF
      IF (ASSOCIATED(matPETSC%vals_matK)) THEN
         DEALLOCATE (matPETSC%vals_matK)
      END IF

      IF (lssolver%reorder_matrix .AND. &
         (lssolver%pctype .EQ. 'PCASM') .OR. lssolver%pctype .EQ. 'PCMG' .OR. lssolver%pctype .EQ. 'PCFIELDSPLIT') THEN
         CALL MatDestroy(matPETSC%matKperm, ierr); CHKERRA(ierr)
         CALL ISDestroy(matPETSC%is_num, ierr); CHKERRA(ierr)
         CALL ISDestroy(matPETSC%is_perm, ierr); CHKERRA(ierr)
         CALL VecDestroy(matPETSC%rhs_vec_perm, ierr); CHKERRA(ierr)
         CALL VecDestroy(matPETSC%solPETSC_vec_perm, ierr); CHKERRA(ierr)
      END IF

      ! destroy everything
      CALL VecDestroy(matPETSC%rhs_vec, ierr); CHKERRA(ierr)
      CALL VecDestroy(matPETSC%solPETSC_vec, ierr); CHKERRA(ierr)
      CALL MatDestroy(matPETSC%matK, ierr); CHKERRA(ierr)
      CALL KSPDestroy(matPETSC%ksp, ierr); CHKERRA(ierr)
      !call PCDestroy(matPETSC%pc, ierr); CHKERRA(ierr)

      IF (lssolver%pctype .EQ. 'PCMG') THEN
         DO i = 1, lssolver%mglevels - 1
            CALL MatDestroy(matPETSC%matP(i), ierr); CHKERRA(ierr)
            CALL MatDestroy(matPETSC%matR(i), ierr); CHKERRA(ierr)
         END DO
         DEALLOCATE (matPETSC%matP, matPETSC%matR)
         CALL destroy_MG_level_operators(matPETSC)
         NULLIFY(matPETSC%matP, matPETSC%matR)

      END IF

   END SUBROUTINE terminate_PETSC

END MODULE solve_petsc
