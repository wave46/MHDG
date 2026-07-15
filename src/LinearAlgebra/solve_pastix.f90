MODULE solve_pastix
   USE iso_c_binding
   USE spmf
   USE pastixf
   USE matrices_types
   USE types
   USE globals
   USE MPI_OMP, ONLY: MPIvar, OMPvar

   IMPLICIT NONE

   TYPE PASTIX_STRUC ! A type to store matrices in CSR format
      ! Variables are defined in a PASTIX specific way
      INTEGER                                               :: n        ! Number of vertices
      ! integer                                               :: nnz      ! Number of non-zeroes in the graph
      ! integer(kind=pastix_int_t), dimension(:), allocatable :: rows     ! Array of size n+1 of indirections to rows for each vertex
      ! integer(kind=pastix_int_t), dimension(:), allocatable :: cols     ! Array of size nnz that corresponds to the global column numbering
      ! integer(kind=pastix_int_t), dimension(:), allocatable :: loc2glob ! Corresponding numbering from local to global
      ! real(kind=c_double),        dimension(:), allocatable :: vals     ! Values stored in the matrix
      ! From here below the variables are PASTIX specific
      type(pastix_data_t), POINTER                      :: pastix_data => NULL()

      type(spmatrix_t), POINTER :: spm => NULL()      ! Spm format of the matrix
      ! integer                               :: pastix_comm ! MPI communicator used by pastix
      ! pastix_int_t, dimension(:), pointer :: perm => null() ! permutation tabular
      ! pastix_int_t, dimension(:), pointer :: invp => null() ! reverse permutation tabular
      ! pastix_float_t, dimension(:), pointer :: rhs => null() ! Right hand side
      INTEGER(c_int)                                           :: info
      INTEGER(kind=pastix_int_t)                               :: nrhs
      REAL(kind=c_double), DIMENSION(:, :), ALLOCATABLE :: x, b
      REAL(kind=c_double), DIMENSION(:), ALLOCATABLE :: rhs
      INTEGER(kind=pastix_int_t)                       :: iparm(iparm_size)
      REAL(kind=c_double)                              :: dparm(dparm_size)

   END TYPE

   TYPE(PASTIX_STRUC)           :: matPASTIX

CONTAINS

   SUBROUTINE set_nthreads(matPASTIX)
      TYPE(PASTIX_STRUC)        :: matPASTIX
      matPASTIX%iparm(IPARM_THREAD_NBR) = OMPvar%Nthreads
   END SUBROUTINE set_nthreads

   !***********************************************
   ! Initialization of matrix
   ! Part specific to PASTIX
   !***********************************************
   SUBROUTINE init_mat_PASTIX(matPASTIX)

      IMPLICIT NONE

      !   TYPE(MAT_CSR_TYP)     :: matCSR
      TYPE(PASTIX_STRUC)        :: matPASTIX

      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke
      cks = 0.

      IF (lssolver%timing) THEN
         CALL cpu_time(tps)
         CALL system_clock(cks, clock_rate)
      END IF

      matPASTIX%n = matK%n

      matPASTIX%nrhs = 1

      CALL pastixInitParam(matPASTIX%iparm, matPASTIX%dparm)

      CALL set_nthreads(matPASTIX)

      ! matPASTIX%iparm(IPARM_SYM) = API_SYM_NO
      ! Verbose mode - Default: PastixVerboseNo
      ! Possible values : PastixVerboseNot, PastixVerboseNo, PastixVerboseYes
      matPASTIX%iparm(IPARM_VERBOSE) = PastixVerboseNot

      ! Refinement mode - Default: PastixRefineGMRES
      ! Possible values : PastixRefineGMRES, PastixRefineCG, PastixRefineSR, PastixRefineBiCGSTAB
      matPASTIX%iparm(IPARM_REFINEMENT) = PastixRefineGMRES

      ! Factorization mode - Default: PastixFactLU
      ! Possible values : PastixFactLLH, PastixFactLDLT, PastixFactLU, PastixFactLLT, PastixFactLDLH
      matPASTIX%iparm(IPARM_FACTORIZATION) = PastixFactLU

      ! Scheduler mode - Default: PastixSchedDynamic
      ! Possible values : PastixSchedSequential, PastixSchedStatic, PastixSchedParsec, PastixSchedStarpu, PastixSchedDynamic
      matPASTIX%iparm(IPARM_SCHEDULER) = PastixSchedDynamic

      matPASTIX%iparm(IPARM_ITERMAX) = 250 ! default is 250

#ifdef THREAD_FUNNELED
      matPASTIX%iparm(IPARM_MPI_THREAD_LEVEL) = PastixMpiThreadFunneled
#endif

      ! We solve A^T x = b because of the csr/csc format
      matPASTIX%iparm(IPARM_TRANSPOSE_SOLVE) = PastixTrans

      CALL pastixInit(matPASTIX%pastix_data, MPI_COMM_WORLD, matPASTIX%iparm, matPASTIX%dparm)
      ALLOCATE(matPASTIX%spm)
      CALL spmInitDist(matPASTIX%spm, MPI_COMM_WORLD)
      CALL build_mat_PASTIX(matPASTIX)

      IF (lssolver%timing) THEN
         CALL cpu_time(tpe)
         CALL system_clock(cke, clock_rate)
         timing%rlstime1 = timing%rlstime1 + (cke - cks)/REAL(clock_rate)
         timing%clstime1 = timing%clstime1 + tpe - tps
      END IF
   END SUBROUTINE init_mat_PASTIX

   !***********************************************
   ! From the generic CSR storage, fills the PASTIX
   ! instance of the matrix
   !***********************************************
   SUBROUTINE build_mat_PASTIX(matPASTIX)

      TYPE(PASTIX_STRUC)        :: matPASTIX
      REAL(kind=c_double)                                      :: normA
      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke, ierr
      cks =0.

      IF (lssolver%timing) THEN
         CALL cpu_time(tps)
         CALL system_clock(cks, clock_rate)
      END IF

      ! we deallocate and reallocate because the mesh can change
      IF (associated(matPASTIX%spm)) THEN
         matPASTIX%spm%colptr = c_null_ptr
         matPASTIX%spm%rowptr = c_null_ptr
         matPASTIX%spm%values = c_null_ptr
         matPASTIX%spm%loc2glob = c_null_ptr
         CALL spmExit(matPASTIX%spm)
         !deallocate (matPASTIX%spm)
      END IF

      !allocate (matPASTIX%spm)

      !call spmInitDist(matPASTIX%spm, MPI_COMM_WORLD)
      ! call spmInit( matPASTIX%spm )

      matPASTIX%spm%baseval = 1              ! 0 or 1, 1-based because of fortran
      matPASTIX%spm%mtxtype = SpmGeneral   ! PastixGeneral, PastixSymmetric, PastixHermitian
      matPASTIX%spm%flttype = SpmDouble      ! Values are stores in double
      matPASTIX%spm%fmttype = SpmCSC         ! Format in CSC
      matPASTIX%spm%n = matK%n      ! Local number of unknowns
      matPASTIX%spm%nnz = matK%nnz    ! Local number of non zeroes
      matPASTIX%spm%dof = 1              ! Degree of freedom per unknown

#ifdef PARALL
      matPASTIX%spm%replicated = 0
#else
      matPASTIX%spm%replicated = 1
#endif

      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      !call spmAlloc(matPASTIX%spm)

      !call spmGetArray(matPASTIX%spm, colptr=colptr, rowptr=rowptr, dvalues=values)
      !colptr(:) = matK%rowptr(:)  ! swap because csr/csc
      !rowptr(:) = matK%cols(:)
      !values(:) = matK%vals(:)
      matPASTIX%spm%colptr = c_loc(matK%rowptr)
      matPASTIX%spm%rowptr = c_loc(matK%cols)
      matPASTIX%spm%values = c_loc(matK%vals)
      matPASTIX%spm%loc2glob = c_loc(matK%loc2glob)

      CALL spmUpdateComputedFields(matPASTIX%spm)

      ! Scale A for better stability with low-rank computations
      CALL spmNorm(SpmFrobeniusNorm, matPASTIX%spm, normA)
      CALL spmScal(1./normA, matPASTIX%spm)

      IF (allocated(matPASTIX%x)) THEN
         DEALLOCATE (matPASTIX%x)
      END IF
      IF (allocated(matPASTIX%b)) THEN
         DEALLOCATE (matPASTIX%b)
      END IF
      IF (allocated(matPASTIX%rhs)) THEN
         DEALLOCATE (matPASTIX%rhs)
      END IF
      ALLOCATE (matPASTIX%x(matPASTIX%spm%nexp, matPASTIX%nrhs))
      ALLOCATE (matPASTIX%b(matPASTIX%spm%nexp, matPASTIX%nrhs))
      ALLOCATE (matPASTIX%rhs(matPASTIX%spm%nexp))

      matPASTIX%b(:, 1) = rhs%vals/normA
      matPASTIX%x(:, 1) = matPASTIX%b(:, 1)

      IF (lssolver%timing) THEN
         CALL cpu_time(tpe)
         CALL system_clock(cke, clock_rate)
         timing%rlstime4 = timing%rlstime4 + (cke - cks)/REAL(clock_rate)
         timing%clstime4 = timing%clstime4 + tpe - tps
      END IF
   END SUBROUTINE build_mat_PASTIX

   !***********************************
   ! PASTIX specific part: consists
   ! only in checking the matrix
   !***********************************
   SUBROUTINE check_mat_PASTIX(matPASTIX)

      TYPE(PASTIX_STRUC)        :: matPASTIX
      type(spmatrix_t), POINTER                                 :: spm2

      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke
      cks = 0.

      IF (lssolver%timing) THEN
         CALL cpu_time(tps)
         CALL system_clock(cks, clock_rate)
      END IF

      ALLOCATE (spm2)
      CALL spmCheckAndCorrect(matPASTIX%spm, spm2, matPASTIX%info)
      IF (matPASTIX%info .NE. 0) THEN
         CALL spmExit(matPASTIX%spm)
         matPASTIX%spm = spm2
      END IF
      DEALLOCATE (spm2)

      IF (lssolver%timing) THEN
         CALL cpu_time(tpe)
         CALL system_clock(cke, clock_rate)
         timing%rlstime2 = timing%rlstime2 + (cke - cks)/REAL(clock_rate)
         timing%clstime2 = timing%clstime2 + tpe - tps
      END IF
   END SUBROUTINE check_mat_PASTIX

   !***********************************************
   ! Analysis of matrix with PASTIX
   !***********************************************
   SUBROUTINE anal_mat_PASTIX(matPASTIX)

      TYPE(PASTIX_STRUC)        :: matPASTIX
      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke
      cks = 0.

      IF (lssolver%timing) THEN
         CALL cpu_time(tps)
         CALL system_clock(cks, clock_rate)
      END IF

      CALL pastix_task_analyze(matPASTIX%pastix_data, matPASTIX%spm, matPASTIX%info)
      IF (matPASTIX%info .NE. 0) THEN
         error STOP "pastix_task_analyze failed"
      END IF

      IF (lssolver%timing) THEN
         CALL cpu_time(tpe)
         CALL system_clock(cke, clock_rate)
         timing%rlstime3 = timing%rlstime3 + (cke - cks)/REAL(clock_rate)
         timing%clstime3 = timing%clstime3 + tpe - tps
      END IF
   END SUBROUTINE anal_mat_PASTIX

   !***********************************************
   ! LU factorization of matrix with PASTIX
   !***********************************************
   SUBROUTINE LU_mat_pastix(matPASTIX)

      TYPE(PASTIX_STRUC)        :: matPASTIX

      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke
      cks = 0.

      IF (lssolver%timing) THEN
         CALL cpu_time(tps)
         CALL system_clock(cks, clock_rate)
      END IF

      CALL pastix_task_numfact(matPASTIX%pastix_data, matPASTIX%spm, matPASTIX%info)
      IF (matPASTIX%info .NE. 0) THEN
         error STOP "pastix_task_numfact failed"
      END IF

      IF (lssolver%timing) THEN
         CALL cpu_time(tpe)
         CALL system_clock(cke, clock_rate)
         timing%rlstime5 = timing%rlstime5 + (cke - cks)/REAL(clock_rate)
         timing%clstime5 = timing%clstime5 + tpe - tps
      END IF

   END SUBROUTINE LU_mat_pastix

   !***********************************************
   ! Solve problem with PASTIX
   !***********************************************
   SUBROUTINE solve_mat_PASTIX(matPASTIX)
      TYPE(PASTIX_STRUC)        :: matPASTIX

      REAL*8                    :: tps, tpe
      INTEGER                   :: cks, clock_rate, cke
      cks = 0.

      IF (lssolver%timing) THEN
         CALL cpu_time(tps)
         CALL system_clock(cks, clock_rate)
      END IF

      CALL pastix_task_solve_and_refine( matPASTIX%pastix_data, matPASTIX%spm%nexp, matPASTIX%nrhs, matPASTIX%b, matPASTIX%spm%nexp, matPASTIX%x, matPASTIX%spm%nexp, matPASTIX%info )
      ! call pastix_task_solve( matPASTIX%pastix_data, matPASTIX%spm%nexp, matPASTIX%nrhs, matPASTIX%x, matPASTIX%spm%nexp, matPASTIX%info )
      IF (matPASTIX%info .NE. 0) THEN
         error STOP "pastix_task_solve failed"
      END IF
      matPASTIX%rhs = matPASTIX%x(:, 1)

      IF (lssolver%timing) THEN
         CALL cpu_time(tpe)
         CALL system_clock(cke, clock_rate)
         timing%rlstime6 = timing%rlstime6 + (cke - cks)/REAL(clock_rate)
         timing%clstime6 = timing%clstime6 + tpe - tps
      END IF
   END SUBROUTINE solve_mat_PASTIX

   SUBROUTINE terminate_mat_PASTIX()

      matPASTIX%spm%colptr = c_null_ptr
      matPASTIX%spm%rowptr = c_null_ptr
      matPASTIX%spm%values = c_null_ptr
      matPASTIX%spm%loc2glob = c_null_ptr
      CALL spmExit(matPASTIX%spm)
      DEALLOCATE (matPASTIX%spm)

      IF (allocated(matPASTIX%x)) THEN
         DEALLOCATE (matPASTIX%x)
         DEALLOCATE (matPASTIX%b)
         DEALLOCATE (matPASTIX%rhs)
      END IF

      !  Destroy the C data structure (Should be last if used to call MPI_Finalize)
      CALL pastixFinalize(matPASTIX%pastix_data)

   END SUBROUTINE terminate_mat_PASTIX

END MODULE solve_pastix
