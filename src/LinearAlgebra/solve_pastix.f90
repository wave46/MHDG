MODULE solve_pastix
  USE matrices_types
  USE types
  USE globals
  IMPLICIT NONE

#include "pastix_fortran.h"

  TYPE PASTIX_STRUC ! A type to store matrices in CSR format
     ! Variables are defined in a PASTIX specific way
     pastix_int_t :: n    ! Number of (locally owned) columns in the matrix
     pastix_int_t :: nnz  ! Number of (locally owned) non-zeros
     pastix_int_t, DIMENSION(:), POINTER    :: rowptr => NULL()  ! Index of first element of each row in cols and vals
     pastix_int_t, DIMENSION(:), POINTER    :: cols => NULL()  ! Column of each element
     pastix_float_t, DIMENSION(:), POINTER  :: vals => NULL()  ! Value of each element
     pastix_int_t, DIMENSION(:), POINTER    :: loc2glob => NULL()  ! Global column number of the local columns
     ! From here below the variables are PASTIX specific
     pastix_data_ptr_t                     :: pastix_data ! Structure to keep information in PaStiX (0 for first CALL)
     INTEGER                               :: pastix_comm ! MPI communicator used by pastix
     pastix_int_t, DIMENSION(:), POINTER :: perm => NULL() ! permutation tabular
     pastix_int_t, DIMENSION(:), POINTER :: invp => NULL() ! reverse permutation tabular
     pastix_float_t, DIMENSION(:), POINTER :: rhs => NULL() ! Right hand side
     pastix_int_t                          :: nrhs   ! right hand side number
     pastix_int_t, DIMENSION(:), POINTER :: iparm => NULL()! Integer parameters
     REAL(float), DIMENSION(:), POINTER :: dparm => NULL()! Floating point parameters
     INTEGER                               :: nbthread    ! Number of threads in PaStiX
     LOGICAL                               :: start  ! keeps track if it is the first solve
  END TYPE PASTIX_STRUC

  INTEGER, PARAMETER                   :: verbose_level_PASTIX = 0 ! Level of verbose (0 = Silent mode, no messages; 1 = Some messages;
  ! 2 = Many messages; 3 = Like a gossip; 4 = Really talking too much...)
  TYPE(PASTIX_STRUC) :: matPASTIX

CONTAINS

  !***********************************************
  ! Initialization of matrix
  ! Part specific to PASTIX
  !***********************************************
  SUBROUTINE init_mat_PASTIX(matPASTIX)
    USE MPI_OMP

    !   TYPE(MAT_CSR_TYP) :: matCSR
    TYPE(PASTIX_STRUC) :: matPASTIX
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! Link the PASTIX instance with the CSR matrix storage
    IF (ASSOCIATED(matPASTIX%rowptr)) THEN
       DEALLOCATE (matPASTIX%rowptr)
    ENDIF
    IF (ASSOCIATED(matPASTIX%loc2glob)) THEN
       DEALLOCATE (matPASTIX%loc2glob)
    ENDIF
    IF (ASSOCIATED(matPASTIX%cols)) THEN
       DEALLOCATE (matPASTIX%cols)
    ENDIF
    IF (ASSOCIATED(matPASTIX%vals)) THEN
       DEALLOCATE (matPASTIX%vals)
    ENDIF
    IF (ASSOCIATED(matPASTIX%rhs)) THEN
       DEALLOCATE (matPASTIX%rhs)
    ENDIF
    matPASTIX%n = matK%n
    matPASTIX%nnz = matK%nnz
    ALLOCATE (matPASTIX%rowptr(matPASTIX%n + 1))
    ALLOCATE (matPASTIX%loc2glob(matPASTIX%n))
    ALLOCATE (matPASTIX%cols(matPASTIX%nnz))
    ALLOCATE (matPASTIX%vals(matPASTIX%nnz))
    ALLOCATE (matPASTIX%rhs(matPASTIX%n))
    matPASTIX%rowptr = matK%rowptr
    matPASTIX%cols = matK%cols
    matPASTIX%vals = matK%vals
    matPASTIX%loc2glob = matK%loc2glob
    matPASTIX%rhs = rhs%vals

    IF (ASSOCIATED(matPASTIX%iparm)) THEN
       DEALLOCATE (matPASTIX%iparm)
    ENDIF
    IF (ASSOCIATED(matPASTIX%dparm)) THEN
       DEALLOCATE (matPASTIX%dparm)
    ENDIF
    IF (ASSOCIATED(matPASTIX%perm)) THEN
       DEALLOCATE (matPASTIX%perm)
    ENDIF
    IF (ASSOCIATED(matPASTIX%invp)) THEN
       DEALLOCATE (matPASTIX%invp)
    ENDIF
    ! Allocation of arrays specific to PASTIX
    ALLOCATE (matPASTIX%iparm(IPARM_SIZE))
    ALLOCATE (matPASTIX%dparm(DPARM_SIZE))
    ALLOCATE (matPASTIX%perm(matPASTIX%n))
    ALLOCATE (matPASTIX%invp(matPASTIX%n))

    matPASTIX%iparm = 0
    matPASTIX%dparm = 0
    matPASTIX%perm = 0
    matPASTIX%invp = 0

    ! Initializes the pastix instance
    matPASTIX%pastix_comm = MPI_COMM_WORLD
    matPASTIX%pastix_data = 0
    matPASTIX%nrhs = 1
    matPASTIX%iparm(IPARM_MODIFY_PARAMETER) = API_NO
    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_INIT
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_INIT
#ifdef THREAD_FUNNELED
    matPASTIX%iparm(IPARM_THREAD_COMM_MODE) = API_THREAD_FUNNELED
#endif

    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
         matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
         matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
         matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)
    matPASTIX%iparm(IPARM_THREAD_NBR) = OMPvar%Nthreads
    matPASTIX%iparm(IPARM_VERBOSE) = verbose_level_PASTIX
    matPASTIX%iparm(IPARM_SYM) = API_SYM_NO
    matPASTIX%iparm(IPARM_FACTORIZATION) = API_FACT_LU
#ifdef THREAD_FUNNELED
    matPASTIX%iparm(IPARM_THREAD_COMM_MODE) = API_THREAD_FUNNELED
#endif

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime1 = timing%rlstime1 + (cke - cks)/REAL(clock_rate)
       timing%clstime1 = timing%clstime1 + tpe - tps
    END IF
  END SUBROUTINE init_mat_PASTIX

  !                real*8  :: tps, tpe
  !                integer :: cks,clock_rate,cke
  !                if (lssolver%timing) then
  !                                call cpu_time(tps)
  !                                call system_clock(cks,clock_rate)
  !                endif

  !                if (lssolver%timing) then
  !                        call cpu_time(tpe)
  !                        call system_clock(cke,clock_rate)
  !                        timing%rlstime1 = timing%rlstime1+(cke-cks)/real(clock_rate)
  !                        timing%clstime1 = timing%clstime1+tpe-tps
  !                end if

  !***********************************************
  ! From the generic CSR storage, fills the PASTIX
  ! instance of the matrix
  !***********************************************
  SUBROUTINE build_mat_PASTIX(matPASTIX)
    !   TYPE(MAT_CSR_TYP) :: matCSR
    TYPE(PASTIX_STRUC) :: matPASTIX
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! Link the PASTIX instance with the CSR matrix storage
    IF (ASSOCIATED(matPASTIX%rowptr)) THEN
       DEALLOCATE (matPASTIX%rowptr)
    ENDIF
    IF (ASSOCIATED(matPASTIX%loc2glob)) THEN
       DEALLOCATE (matPASTIX%loc2glob)
    ENDIF
    IF (ASSOCIATED(matPASTIX%cols)) THEN
       DEALLOCATE (matPASTIX%cols)
    ENDIF
    IF (ASSOCIATED(matPASTIX%vals)) THEN
       DEALLOCATE (matPASTIX%vals)
    ENDIF
    IF (ASSOCIATED(matPASTIX%rhs)) THEN
       DEALLOCATE (matPASTIX%rhs)
    ENDIF
    matPASTIX%n = matK%n
    matPASTIX%nnz = matK%nnz
    ALLOCATE (matPASTIX%rowptr(matPASTIX%n + 1))
    ALLOCATE (matPASTIX%loc2glob(matPASTIX%n))
    ALLOCATE (matPASTIX%cols(matPASTIX%nnz))
    ALLOCATE (matPASTIX%vals(matPASTIX%nnz))
    ALLOCATE (matPASTIX%rhs(matPASTIX%n))
    matPASTIX%rowptr = matK%rowptr
    matPASTIX%cols = matK%cols
    matPASTIX%vals = matK%vals
    matPASTIX%loc2glob = matK%loc2glob
    matPASTIX%rhs = rhs%vals

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime4 = timing%rlstime4 + (cke - cks)/REAL(clock_rate)
       timing%clstime4 = timing%clstime4 + tpe - tps
    END IF
  END SUBROUTINE build_mat_PASTIX

  !***********************************
  ! PASTIX specific part: consists
  ! only in checking the matrix
  !***********************************
  SUBROUTINE check_mat_PASTIX(matPASTIX)

    TYPE(PASTIX_STRUC) :: matPASTIX

    INTEGER :: nnz
    pastix_data_ptr_t :: data_check
    pastix_int_t :: flagcor
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! Removes multiple entries and checks symmetry
    nnz = matPASTIX%rowptr(matPASTIX%n + 1) - 1
    flagcor = API_YES
    CALL pastix_fortran_checkMatrix(data_check, matPASTIX%pastix_comm, matPASTIX%iparm(IPARM_VERBOSE), &
         matPASTIX%iparm(IPARM_SYM), flagcor, &
         matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
         matPASTIX%loc2glob, matPASTIX%iparm(IPARM_DOF_NBR))
    IF (nnz .NE. (matPASTIX%rowptr(matPASTIX%n + 1) - 1)) THEN
       DEALLOCATE (matPASTIX%cols, matPASTIX%vals)
       matPASTIX%nnz = matPASTIX%rowptr(matPASTIX%n + 1) - 1
       ALLOCATE (matPASTIX%cols(matPASTIX%nnz))
       ALLOCATE (matPASTIX%vals(matPASTIX%nnz))
       CALL pastix_fortran_checkMatrix_END(data_check, matPASTIX%iparm(IPARM_VERBOSE), &
            matPASTIX%cols, matPASTIX%vals, matPASTIX%iparm(IPARM_DOF_NBR))
    ENDIF

    ! The matrix does not need to be checked anymore
    matPASTIX%iparm(IPARM_MATRIX_VERIFICATION) = API_NO

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime2 = timing%rlstime2 + (cke - cks)/REAL(clock_rate)
       timing%clstime2 = timing%clstime2 + tpe - tps
    END IF
  END SUBROUTINE check_mat_PASTIX

  !***********************************************
  ! Analysis of matrix with PASTIX
  !***********************************************
  SUBROUTINE anal_mat_PASTIX(matPASTIX)
    USE matrices_tools, ONLY: dump_CSR

    TYPE(PASTIX_STRUC) :: matPASTIX
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! Analysis phase
    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_ORDERING
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_ANALYSE
    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
         matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
         matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
         matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime3 = timing%rlstime3 + (cke - cks)/REAL(clock_rate)
       timing%clstime3 = timing%clstime3 + tpe - tps
    END IF
  END SUBROUTINE anal_mat_PASTIX

  !***********************************************
  ! LU factorization of matrix with PASTIX
  !***********************************************
  SUBROUTINE LU_mat_pastix(matPASTIX)

    TYPE(PASTIX_STRUC) :: matPASTIX

    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! Factorization phase
    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_NUMFACT
    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
         matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
         matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
         matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime5 = timing%rlstime5 + (cke - cks)/REAL(clock_rate)
       timing%clstime5 = timing%clstime5 + tpe - tps
    END IF

  END SUBROUTINE LU_mat_pastix

  !***********************************************
  ! Solve problem with PASTIX
  !***********************************************
  SUBROUTINE solve_mat_PASTIX(matPASTIX)
    TYPE(PASTIX_STRUC) :: matPASTIX

    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0


    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    ! Solve phase
    matPASTIX%iparm(IPARM_TRANSPOSE_SOLVE) = API_YES
    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_SOLVE
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_SOLVE
    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
         matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
         matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
         matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime6 = timing%rlstime6 + (cke - cks)/REAL(clock_rate)
       timing%clstime6 = timing%clstime6 + tpe - tps
    END IF
  END SUBROUTINE solve_mat_PASTIX

  !SUBROUTINE free_mat_PASTIX(matPASTIX)
  !   TYPE(PASTIX_STRUC) :: matPASTIX
  !
  !   IF (associated(matPASTIX%rowptr)) then
  !      deallocate(matPASTIX%rowptr)
  !   endif
  !   IF (associated(matPASTIX%loc2glob)) then
  !      deallocate(matPASTIX%loc2glob)
  !   endif
  !   IF (associated(matPASTIX%cols)) then
  !      deallocate(matPASTIX%cols)
  !   endif
  !   IF (associated(matPASTIX%vals)) then
  !      deallocate(matPASTIX%vals)
  !   endif
  !   IF (associated(matPASTIX%rhs)) then
  !      deallocate(matPASTIX%rhs)
  !   endif
  !END SUBROUTINE free_mat_PASTIX

  SUBROUTINE terminate_mat_PASTIX()

    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_CLEAN
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_CLEAN

    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
         matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
         matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
         matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)

    IF (ASSOCIATED(matPASTIX%rowptr)) THEN
       DEALLOCATE (matPASTIX%rowptr)
    ENDIF
    IF (ASSOCIATED(matPASTIX%loc2glob)) THEN
       DEALLOCATE (matPASTIX%loc2glob)
    ENDIF
    IF (ASSOCIATED(matPASTIX%cols)) THEN
       DEALLOCATE (matPASTIX%cols)
    ENDIF
    IF (ASSOCIATED(matPASTIX%vals)) THEN
       DEALLOCATE (matPASTIX%vals)
    ENDIF
    IF (ASSOCIATED(matPASTIX%rhs)) THEN
       DEALLOCATE (matPASTIX%rhs)
    ENDIF
    IF (ASSOCIATED(matPASTIX%iparm)) THEN
       DEALLOCATE (matPASTIX%iparm)
    ENDIF
    IF (ASSOCIATED(matPASTIX%dparm)) THEN
       DEALLOCATE (matPASTIX%dparm)
    ENDIF
    IF (ASSOCIATED(matPASTIX%perm)) THEN
       DEALLOCATE (matPASTIX%perm)
    ENDIF
    IF (ASSOCIATED(matPASTIX%invp)) THEN
       DEALLOCATE (matPASTIX%invp)
    ENDIF
  END SUBROUTINE terminate_mat_PASTIX

END MODULE solve_pastix
