MODULE solve_petsc
#include "petsc/finclude/petsc.h"
  USE matrices_types
  USE types
  USE globals
  USE MPI_OMP
  USE petsc
  USE initialization
  IMPLICIT NONE


  TYPE PETSC_STRUC ! A type to store structures for the PETSC library
     Mat                               :: matK          ! sparse matrix for PSBLAS
     KSP                               :: ksp           ! solver
     PC                                :: pc            ! preconditioner
     Vec                               :: solPETSC_vec      ! solution
     Vec                               :: rhs_vec           ! right hand side
     PetscInt                          :: n, nnz        ! matrix n and number of non-zeros
     PetscInt, DIMENSION(:), POINTER   :: rowptr => NULL()
     PetscInt, DIMENSION(:), POINTER   :: loc2glob => NULL()
     PetscInt, DIMENSION(:), POINTER   :: cols => NULL()
     PetscScalar, DIMENSION(:), POINTER :: vals_matK => NULL()
     PetscScalar, DIMENSION(:), POINTER :: vals_rhs => NULL()
     PetscScalar                       :: residue   ! UNPRECONDTIONED residue norm of the linear system
     PetscInt                          :: its      ! Previous number of iterations
     KSPConvergedReason                :: convergedReason

  END TYPE PETSC_STRUC

  TYPE(PETSC_STRUC) :: matPETSC


CONTAINS

  !***********************************************
  ! Initialization of matrix
  ! Part specific to PETSc
  !***********************************************
  SUBROUTINE init_mat_PETSC(matPETSC)
    USE petsc
    USE MPI_OMP

    TYPE(PETSC_STRUC)     :: matPETSC
    PetscErrorCode        :: ierr
    REAL*8                :: tps, tpe
    INTEGER               :: cks, clock_rate, cke

    cks = 0
    cke = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! Link the PETSC instance with the CSR matrix storage
    IF (ASSOCIATED(matPETSC%rowptr)) THEN
       DEALLOCATE (matPETSC%rowptr)
    ENDIF
    IF (ASSOCIATED(matPETSC%loc2glob)) THEN
       DEALLOCATE (matPETSC%loc2glob)
    ENDIF
    IF (ASSOCIATED(matPETSC%cols)) THEN
       DEALLOCATE (matPETSC%cols)
    ENDIF
    IF (ASSOCIATED(matPETSC%vals_matK)) THEN
       DEALLOCATE (matPETSC%vals_matK)
    ENDIF
    IF (ASSOCIATED(matPETSC%vals_rhs)) THEN
       DEALLOCATE (matPETSC%vals_rhs)
    ENDIF

    matPETSC%n = matK%n
    matPETSC%nnz = matK%nnz
    ALLOCATE (matPETSC%rowptr(matPETSC%n + 1))
    ALLOCATE (matPETSC%loc2glob(matPETSC%n))
    ALLOCATE (matPETSC%cols(matPETSC%nnz))
    ALLOCATE (matPETSC%vals_matK(matPETSC%nnz))
    ALLOCATE (matPETSC%vals_rhs(matPETSC%n))
    matPETSC%rowptr    = matK%rowptr
    matPETSC%cols      = matK%cols
    matPETSC%vals_matK = matK%vals
    matPETSC%loc2glob  = matK%loc2glob
    matPETSC%vals_rhs  = rhs%vals


    IF (matK%start) THEN
       matPETSC%residue   = 0
       matPETSC%its       = 0
#ifdef PARALL
       ! create matrix and set it up
       CALL MatCreate(MPI_COMM_WORLD,matPETSC%matK , ierr)
       CALL MatSetType(matPETSC%matK, MATMPIAIJ, ierr)
       CALL MatSetSizes(matPETSC%matK, matPETSC%n,matPETSC%n, PETSC_DETERMINE,PETSC_DETERMINE,ierr);
       CALL MatSetUp(matPETSC%matK, ierr)

       ! create rhs and solution vector
       CALL MatCreateVecs(matPETSC%matK, matPETSC%rhs_vec, matPETSC%solPETSC_vec, ierr)
       CALL VecAssemblyBegin(matPETSC%rhs_vec,ierr)
       CALL VecAssemblyEnd(matPETSC%rhs_vec,ierr)
       CALL VecAssemblyBegin(matPETSC%solPETSC_vec,ierr)
       CALL VecAssemblyEnd(matPETSC%solPETSC_vec,ierr)

#else
       ! create matrix and set it up
       CALL MatCreate(MPI_COMM_WORLD,matPETSC%matK , ierr)
       CALL MatSetType(matPETSC%matK, MATSEQAIJ, ierr)
       CALL MatSetSizes(matPETSC%matK, matPETSC%n,matPETSC%n, matPETSC%n,matPETSC%n,ierr);
       CALL MatSetUp(matPETSC%matK, ierr)

       ! Create solution vector / initial guess also
       CALL MatCreateVecs(matPETSC%matK, matPETSC%rhs_vec, matPETSC%solPETSC_vec, ierr)
       CALL VecAssemblyBegin(matPETSC%rhs_vec,ierr)
       CALL VecAssemblyEnd(matPETSC%rhs_vec,ierr)
       CALL VecAssemblyBegin(matPETSC%solPETSC_vec,ierr)
       CALL VecAssemblyEnd(matPETSC%solPETSC_vec,ierr)
#endif
       ! create the ksp instance
       CALL KSPCreate(PETSC_COMM_WORLD, matPETSC%ksp, ierr);
    ENDIF

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime1 = timing%rlstime1 + (cke - cks)/REAL(clock_rate)
       timing%clstime1 = timing%clstime1 + tpe - tps
    END IF
  END SUBROUTINE init_mat_PETSC


  SUBROUTINE preallocate_matrix_PETSC(matPETSC)
    TYPE(PETSC_STRUC)      :: matPETSC
    INTEGER, ALLOCATABLE   :: counter_diag(:), counter_offdiag(:), istart_vec(:), iend_vec(:)
    PetscInt               :: istart, iend, r, rr, rrr,  M, N, I, J
    PetscErrorCode         :: ierr

    ! Get global sizes
    CALL MatGetSize(matPETSC%matK,M,N, ierr)

    ! Get local row ownership of the matrix
    CALL MatGetOwnershipRange(matPETSC%matK, istart, iend, ierr)

    ! allocate diagonal and offdiagonal non-zero counters and a vector containing the beginning and ends!
    ! of the row ownership, per process. Each process owns X rows of the matrix, a slice of matrix.
    ALLOCATE(counter_diag(N))
    ALLOCATE(counter_offdiag(N))
    ALLOCATE(istart_vec(MPIvar%glob_size))
    ALLOCATE(iend_vec(MPIvar%glob_size))

    counter_diag = 0
    counter_offdiag = 0
    istart_vec = 0
    iend_vec = 0
    istart_vec(MPIvar%glob_id+1) = istart
    iend_vec(MPIvar%glob_id+1) = iend

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,istart_vec,MPIvar%glob_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,iend_vec,MPIvar%glob_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    ! count non-zero elements
    DO r = 1, matPETSC%n
       DO rr = matPETSC%rowptr(r), (matPETSC%rowptr(r+1)-1)
          I = matPETSC%loc2glob(r)-1;
          J = matPETSC%cols(rr)-1;

          ! Find out on which process' slice I am (when found, exit the loop, rrr will be kept)
          DO rrr = 1, MPIvar%glob_size
             IF(((I+1) .GE. istart_vec(rrr)) .AND. ((I+1) .LE. iend_vec(rrr))) THEN
                EXIT
             ENDIF
          ENDDO

          ! put in or out diagonal (always on diagonal for serial)
          IF((J .GE. istart_vec(rrr)) .AND. (J .LE. iend_vec(rrr)-1)) THEN
             counter_diag(I+1) = counter_diag(I+1) + 1
          ELSE
             counter_offdiag(I+1) = counter_offdiag(I+1) + 1
          ENDIF
       END DO
    END DO

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,counter_diag,N, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,counter_offdiag,N, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MatMPIAIJSetPreallocation(matPETSC%matK,PETSC_DECIDE,counter_diag(istart+1:iend),PETSC_DECIDE,counter_offdiag(istart+1:iend), ierr); CHKERRA(ierr)
#else
    CALL MatSeqAIJSetPreallocation(matPETSC%matK, PETSC_DECIDE, counter_diag, ierr)
#endif

    DEALLOCATE(counter_diag)
    DEALLOCATE(counter_offdiag)
    DEALLOCATE(istart_vec)
    DEALLOCATE(iend_vec)
  END SUBROUTINE preallocate_matrix_PETSC

  !***********************************************
  ! From the generic CSR storage, fills the PETSC
  ! instance of the matrix
  !***********************************************
  SUBROUTINE build_mat_PETSC(matPETSC)
    TYPE(PETSC_STRUC)     :: matPETSC
    PetscErrorCode        :: ierr
    PetscInt              :: one = 1, r, rr
    REAL*8                :: tps, tpe
    INTEGER               :: cks, clock_rate, cke

    cks = 0
    cke = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! set the values in the matrix
    DO r = 1, matPETSC%n
       DO rr = matPETSC%rowptr(r), (matPETSC%rowptr(r+1)-1)
          CALL MatSetValues(matPETSC%matK,one, matPETSC%loc2glob(r)-1, one, matPETSC%cols(rr)-1, matPETSC%vals_matK(rr), INSERT_VALUES, ierr); CHKERRA(ierr)
       END DO
    END DO

    CALL MatAssemblyBegin(matPETSC%matK,MAT_FINAL_ASSEMBLY, ierr)
    CALL MatAssemblyEnd(matPETSC%matK,MAT_FINAL_ASSEMBLY, ierr)

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime2 = timing%rlstime2 + (cke - cks)/REAL(clock_rate)
       timing%clstime2 = timing%clstime2 + tpe - tps
    END IF
  END SUBROUTINE build_mat_PETSC


  SUBROUTINE fill_vec_PETSC(matPETSC)
    TYPE(PETSC_STRUC)     :: matPETSC
    PetscErrorCode        :: ierr
    REAL*8                :: tps, tpe
    INTEGER               :: cks, clock_rate, cke

    cks = 0
    cke = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! set values
    CALL VecSetValues(matPETSC%rhs_vec,matPETSC%n, matPETSC%loc2glob-1,rhs%vals, INSERT_VALUES, ierr);

    CALL VecAssemblyBegin(matPETSC%rhs_vec,ierr)
    CALL VecAssemblyEnd(matPETSC%rhs_vec,ierr)


    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime3 = timing%rlstime3 + (cke - cks)/REAL(clock_rate)
       timing%clstime3 = timing%clstime3 + tpe - tps
    END IF
  END SUBROUTINE fill_vec_PETSC


  !***********************************************
  ! Solve problem with PETSC
  !***********************************************
  SUBROUTINE solve_mat_PETSC(matPETSC, ir)
    TYPE(PETSC_STRUC)     :: matPETSC
    INTEGER, INTENT(IN)   :: ir
    REAL*8                :: tps, tpe
    INTEGER               :: cks, clock_rate, cke
    !PetscViewer           :: viewer
    !PetscDraw             :: draw
    KSPType               :: ksp_type
    PCType                :: pc_type
    PetscViewerAndFormat  :: vf
    PetscErrorCode        :: ierr

    cks = 0
    cke = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! Build preconditioner based on linear system matrix
    CALL KSPSetOperators(matPETSC%ksp, matPETSC%matK, matPETSC%matK, ierr);

    ! set or not the initial guess to zero
    IF((.NOT. lssolver%igz) .AND. (lssolver%kspmethd .NE. 'KSPPREONLY')) THEN
       CALL KSPSetInitialGuessNonzero(matPETSC%ksp, PETSC_TRUE, ierr)
    ENDIF

    ! set solver and preconditione type from param.txt
    CALL set_ksp_pc_type(ksp_type,pc_type)
    CALL KSPSetType(matPETSC%ksp,ksp_type, ierr);

    ! set the norm type
    IF(lssolver%kspnorm .EQ. 1) THEN
       CALL KSPSetNormType(matPETSC%ksp,KSP_NORM_PRECONDITIONED, ierr)
    ELSEIF(lssolver%kspnorm .EQ. 2) THEN
       CALL KSPSetNormType(matPETSC%ksp,KSP_NORM_UNPRECONDITIONED, ierr)
    ENDIF

    ! set tolerances
    CALL KSPSetTolerances(matPETSC%ksp, lssolver%rtol, lssolver%atol, PETSC_DEFAULT_REAL, lssolver%kspitmax, ierr)

    ! get the preconditioner and set it
    CALL KSPGetPC(matPETSC%ksp, matPETSC%pc, ierr);
    CALL PCSetType(matPETSC%pc, pc_type, ierr);
    ! in case of PCMG set the additional info
    SELECT CASE (lssolver%pctype)
    CASE ('PCMG')
       CALL PCMGSetLevels(matPETSC%pc,lssolver%mglevels,MPI_COMM_WORLD, ierr)
       IF    (lssolver%mgtypeform .EQ. 1) THEN
          CALL PCMGSetType(matPETSC%pc,PC_MG_MULTIPLICATIVE, ierr)
       ELSEIF(lssolver%mgtypeform .EQ. 2) THEN
          CALL PCMGSetType(matPETSC%pc,PC_MG_ADDITIVE, ierr)
       ELSEIF(lssolver%mgtypeform .EQ. 3) THEN
          CALL PCMGSetType(matPETSC%pc,PC_MG_FULL, ierr)
       ELSEIF(lssolver%mgtypeform .EQ. 4) THEN
          CALL PCMGSetType(matPETSC%pc,PC_MG_KASKADE, ierr)
       ENDIF
    END SELECT

    ! set preconditioner up
    CALL PCSetUp(matPETSC%pc, ierr);
    CALL KSPSetUp(matPETSC%ksp, ierr);

    ! print out the residue
    IF(lssolver%kspitrace .EQ. 1) THEN
       CALL PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,vf,ierr)
       CALL KSPMonitorSet(matPETSC%ksp,KSPMonitorTrueResidual,vf,PetscViewerAndFormatDestroy,ierr)
    ENDIF

    ! reuse or not preconditioner
    IF(lssolver%rprecond .EQ. 1) THEN
       CALL KSPSetReusePreconditioner(matPETSC%ksp,PETSC_FALSE,ierr)
    ELSEIF(lssolver%rprecond .EQ. 2) THEN
       IF(MOD(ir-1, lssolver%Nrprecond) .EQ. 0) THEN
          CALL KSPSetReusePreconditioner(matPETSC%ksp,PETSC_FALSE,ierr)
          IF ((MPIvar%glob_id .EQ. 0) .AND. (ir .NE. 1)) THEN
             WRITE(*,*) 'PETSC is re-computing preconditioner'
          ENDIF
       ELSE
          CALL KSPSetReusePreconditioner(matPETSC%ksp,PETSC_TRUE,ierr)
       ENDIF
    ELSEIF(lssolver%rprecond .EQ. 3) THEN
       IF(matPETSC%its .EQ. lssolver%kspitmax) THEN
          CALL KSPSetReusePreconditioner(matPETSC%ksp,PETSC_FALSE,ierr)
          IF (MPIvar%glob_id .EQ. 0) THEN
             WRITE(*,*) 'PETSC is re-computing preconditioner'
          ENDIF
       ELSE
          CALL KSPSetReusePreconditioner(matPETSC%ksp,PETSC_TRUE,ierr)
       ENDIF
    ELSE
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE(*,*) 'Wrong choice, value must be between 1, 2 or 3, stopping. '
          STOP
       ENDIF
    ENDIF

    !call KSPView(matPETSC%ksp,PETSC_VIEWER_STDOUT_WORLD, ierr)

    !call PetscViewerDrawOpen(PETSC_COMM_WORLD, PETSC_NULL_CHARACTER, PETSC_NULL_CHARACTER, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, viewer, ierr);
    !call PetscViewerDrawGetDraw(viewer, 0, draw, ierr);
    !call PetscDrawSetSaveFinalImage(draw, 'matrix_plot.ppm', ierr);
    !call MatView(matPETSC%matK, viewer, ierr)


    CALL KSPSolve(matPETSC%ksp, matPETSC%rhs_vec, matPETSC%solPETSC_vec, ierr); CHKERRA(ierr)

    CALL KSPGetIterationNumber(matPETSC%ksp,matPETSC%its,ierr)
    CALL KSPGetConvergedReason(matPETSC%ksp, matPETSC%convergedReason, ierr)

    ! compute the unpreconditioned relative residue
    CALL compute_residue

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime4 = timing%rlstime4 + (cke - cks)/REAL(clock_rate)
       timing%clstime4 = timing%clstime4 + tpe - tps
    END IF
  END SUBROUTINE solve_mat_PETSC


  SUBROUTINE PETSC_retrieve_array(vector, aux_sol, size)
    Vec, INTENT(IN)           :: vector
    INTEGER, INTENT(IN)       :: size
    PetscScalar, INTENT(INOUT)     :: aux_sol(size)
    PetscScalar, POINTER      :: data_array(:)
    PetscErrorCode            :: ierr
    ! copy the vector values into the array
    CALL VecGetArrayReadF90(vector,data_array,ierr)
    aux_sol = data_array
    CALL VecRestoreArrayReadF90(vector,data_array,ierr)
  END SUBROUTINE PETSC_retrieve_array

  SUBROUTINE PETSC_set_array(vector, aux_sol, size)
    Vec, INTENT(INOUT)           :: vector
    INTEGER, INTENT(IN)       :: size
    PetscScalar, INTENT(IN)     :: aux_sol(size)
    PetscScalar, POINTER      :: data_array(:)
    PetscErrorCode            :: ierr
    ! copy the vector values into the array
    CALL VecGetArrayF90(vector,data_array,ierr)
    data_array = aux_sol
    CALL VecRestoreArrayF90(vector,data_array,ierr)
  END SUBROUTINE PETSC_set_array


  SUBROUTINE compute_residue
    Vec                       :: true_residue_vec
    PetscReal                 :: alpha=-1
    PetscScalar               :: true_residue, norm
    PetscErrorCode            :: ierr

    CALL MatCreateVecs(matPETSC%matK, true_residue_vec, true_residue_vec, ierr)
    CALL VecScale(matPETSC%rhs_vec, alpha, ierr)
    CALL MatMultAdd(matPETSC%matK, matPETSC%solPETSC_vec,matPETSC%rhs_vec, true_residue_vec, ierr)
    CALL VecScale(matPETSC%rhs_vec, alpha, ierr)
    CALL VecNorm(true_residue_vec, NORM_2, true_residue, ierr)
    CALL VecNorm(matPETSC%rhs_vec, NORM_2, norm, ierr)
    matPETSC%residue = true_residue/norm
    CALL VecDestroy(true_residue_vec, ierr)

  END SUBROUTINE compute_residue

  SUBROUTINE set_ksp_pc_type(ksp_type, pc_type)
    KSPType, INTENT(OUT)       :: ksp_type
    PCType, INTENT(OUT)       :: pc_type

    SELECT CASE (lssolver%kspmethd)
    CASE ('KSPCG')
       ksp_type = KSPCG
    CASE ('KSPBCGS')
       ksp_type = KSPBCGS
    CASE ('KSPGMRES')
       ksp_type = KSPGMRES
    CASE ('KSPMINRES')
       ksp_type = KSPMINRES
    CASE ('KSPPREONLY')
       ksp_type = KSPPREONLY
    CASE ('KSPRICHARDSON')
       ksp_type = KSPRICHARDSON
    CASE default
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE(*,*) "KSP METHOD NOT IN LIST, add it and recompile"
          STOP
       ENDIF
    END SELECT

    SELECT CASE (lssolver%pctype)
    CASE ('PCNONE')
       pc_type = PCNONE
    CASE ('PCJACOBI')
       pc_type = PCJACOBI
    CASE ('PCSOR')
       pc_type = PCSOR
    CASE ('PCLU')
       pc_type = PCLU
    CASE ('PCILU')
       pc_type = PCILU
    CASE ('PCGAMG')
       pc_type = PCGAMG
    CASE ('PCMG')
       pc_type = PCMG
    CASE default
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE(*,*) "PC TYPE NOT IN LIST, add it and recompile"
          STOP
       ENDIF
    END SELECT
  END SUBROUTINE set_ksp_pc_type



  SUBROUTINE terminate_PETSC()
    PetscErrorCode        :: ierr
    IF (ASSOCIATED(matPETSC%rowptr)) THEN
       DEALLOCATE (matPETSC%rowptr)
    ENDIF
    IF (ASSOCIATED(matPETSC%loc2glob)) THEN
       DEALLOCATE (matPETSC%loc2glob)
    ENDIF
    IF (ASSOCIATED(matPETSC%cols)) THEN
       DEALLOCATE (matPETSC%cols)
    ENDIF
    IF (ASSOCIATED(matPETSC%vals_rhs)) THEN
       DEALLOCATE (matPETSC%vals_rhs)
    ENDIF
    IF (ASSOCIATED(matPETSC%vals_matK)) THEN
       DEALLOCATE (matPETSC%vals_matK)
    ENDIF

    ! destroy everything
    CALL VecDestroy(matPETSC%rhs_vec, ierr); CHKERRA(ierr)
    CALL VecDestroy(matPETSC%solPETSC_vec, ierr); CHKERRA(ierr)
    CALL MatDestroy(matPETSC%matK, ierr); CHKERRA(ierr)
    CALL KSPDestroy(matPETSC%ksp, ierr); CHKERRA(ierr)
  END SUBROUTINE terminate_PETSC

END MODULE solve_petsc
