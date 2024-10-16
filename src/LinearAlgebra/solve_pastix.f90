MODULE solve_pastix
   use iso_c_binding
   use spmf
   use pastixf
   USE matrices_types
   USE types
   USE globals
   IMPLICIT NONE

   TYPE PASTIX_STRUC ! A type to store matrices in CSR format
      ! Variables are defined in a PASTIX specific way
      integer                                               :: n        ! Number of vertices
      ! integer                                               :: nnz      ! Number of non-zeroes in the graph
      ! integer(kind=pastix_int_t), dimension(:), allocatable :: rows     ! Array of size n+1 of indirections to rows for each vertex
      ! integer(kind=pastix_int_t), dimension(:), allocatable :: cols     ! Array of size nnz that corresponds to the global column numbering
      ! integer(kind=pastix_int_t), dimension(:), allocatable :: loc2glob ! Corresponding numbering from local to global
      ! real(kind=c_double),        dimension(:), allocatable :: vals     ! Values stored in the matrix
      ! From here below the variables are PASTIX specific
      type(pastix_data_t), pointer                      :: pastix_data

      type(spmatrix_t), pointer :: spm      ! Spm format of the matrix
      ! integer                               :: pastix_comm ! MPI communicator used by pastix
      ! pastix_int_t, dimension(:), pointer :: perm => null() ! permutation tabular
      ! pastix_int_t, dimension(:), pointer :: invp => null() ! reverse permutation tabular
      ! pastix_float_t, dimension(:), pointer :: rhs => null() ! Right hand side
      integer(c_int)                                           :: info
      integer(kind=pastix_int_t)                               :: nrhs
      real(kind=c_double), dimension(:, :), allocatable :: x, b
      real(kind=c_double), dimension(:), allocatable :: rhs
      integer(kind=pastix_int_t)                       :: iparm(iparm_size)
      real(kind=c_double)                              :: dparm(dparm_size)

      ! Integer                               :: nbthread = OMPvar%Nthreads    ! Number of threads in PaStiX
   END TYPE

   TYPE(PASTIX_STRUC) :: matPASTIX

CONTAINS

   SUBROUTINE set_nthreads(matPASTIX)
      use MPI_OMP
      TYPE(PASTIX_STRUC) :: matPASTIX
      matPASTIX%iparm(IPARM_THREAD_NBR) = OMPvar%Nthreads
   end SUBROUTINE set_nthreads

   !***********************************************
   ! Initialization of matrix
   ! Part specific to PASTIX
   !***********************************************
   SUBROUTINE init_mat_PASTIX(matPASTIX)
      use pastixf

      IMPLICIT NONE

      !   TYPE(MAT_CSR_TYP) :: matCSR
      TYPE(PASTIX_STRUC) :: matPASTIX
      integer(kind=spm_int_t), dimension(:), pointer :: rowptr
      integer(kind=spm_int_t), dimension(:), pointer :: colptr
      real(kind=c_double), dimension(:), pointer :: values

      real*8  :: tps, tpe
      integer :: cks, clock_rate, cke
      if (lssolver%timing) then
         call cpu_time(tps)
         call system_clock(cks, clock_rate)
      end if

      matPASTIX%n = matK%n

      matPASTIX%nrhs = 1

      call pastixInitParam(matPASTIX%iparm, matPASTIX%dparm)

      call set_nthreads(matPASTIX)

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

      call pastixInit(matPASTIX%pastix_data, MPI_COMM_WORLD, matPASTIX%iparm, matPASTIX%dparm)

      call build_mat_PASTIX(matPASTIX)

      if (lssolver%timing) then
         call cpu_time(tpe)
         call system_clock(cke, clock_rate)
         timing%rlstime1 = timing%rlstime1 + (cke - cks)/real(clock_rate)
         timing%clstime1 = timing%clstime1 + tpe - tps
      end if
   END SUBROUTINE init_mat_PASTIX

   !***********************************************
   ! From the generic CSR storage, fills the PASTIX
   ! instance of the matrix
   !***********************************************
   SUBROUTINE build_mat_PASTIX(matPASTIX)

      use pastixf
      TYPE(PASTIX_STRUC) :: matPASTIX
      integer(kind=spm_int_t), dimension(:), pointer :: rowptr
      integer(kind=spm_int_t), dimension(:), pointer :: colptr
      real(kind=c_double), dimension(:), pointer :: values
      real(kind=c_double)                                      :: normA
      real*8  :: tps, tpe
      integer :: cks, clock_rate, cke

      if (lssolver%timing) then
         call cpu_time(tps)
         call system_clock(cks, clock_rate)
      end if

      ! we deallocate and reallocate because the mesh can change
      if (associated(matPASTIX%spm)) then
         matPASTIX%spm%colptr = c_null_ptr
         matPASTIX%spm%rowptr = c_null_ptr
         matPASTIX%spm%values = c_null_ptr
         matPASTIX%spm%loc2glob = c_null_ptr
         call spmExit(matPASTIX%spm)
         deallocate (matPASTIX%spm)
      end if

      allocate (matPASTIX%spm)

      call spmInitDist(matPASTIX%spm, MPI_COMM_WORLD)
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

      call spmUpdateComputedFields(matPASTIX%spm)
      call spmAlloc(matPASTIX%spm)

      call spmGetArray(matPASTIX%spm, colptr=colptr, rowptr=rowptr, dvalues=values)
      colptr(:) = matK%rowptr(:)  ! swap because csr/csc
      rowptr(:) = matK%cols(:)
      values(:) = matK%vals(:)

      matPASTIX%spm%loc2glob = c_loc(matK%loc2glob)

      ! Scale A for better stability with low-rank computations
      call spmNorm(SpmFrobeniusNorm, matPASTIX%spm, normA)
      call spmScal(1./normA, matPASTIX%spm)

      if (allocated(matPASTIX%x)) then
         deallocate (matPASTIX%x)
         deallocate (matPASTIX%b)
         deallocate (matPASTIX%rhs)
      end if
      allocate (matPASTIX%x(matPASTIX%spm%nexp, matPASTIX%nrhs))
      allocate (matPASTIX%b(matPASTIX%spm%nexp, matPASTIX%nrhs))
      allocate (matPASTIX%rhs(matPASTIX%spm%nexp))

      matPASTIX%b(:, 1) = rhs%vals/normA
      matPASTIX%x(:, 1) = matPASTIX%b(:, 1)

      if (lssolver%timing) then
         call cpu_time(tpe)
         call system_clock(cke, clock_rate)
         timing%rlstime4 = timing%rlstime4 + (cke - cks)/real(clock_rate)
         timing%clstime4 = timing%clstime4 + tpe - tps
      end if
   END SUBROUTINE build_mat_PASTIX

   !***********************************
   ! PASTIX specific part: consists
   ! only in checking the matrix
   !***********************************
   SUBROUTINE check_mat_PASTIX(matPASTIX)

      TYPE(PASTIX_STRUC) :: matPASTIX
      type(spmatrix_t), pointer                                 :: spm2

      real*8  :: tps, tpe
      integer :: cks, clock_rate, cke
      if (lssolver%timing) then
         call cpu_time(tps)
         call system_clock(cks, clock_rate)
      end if

      allocate (spm2)
      call spmCheckAndCorrect(matPASTIX%spm, spm2, matPASTIX%info)
      if (matPASTIX%info .ne. 0) then
         call spmExit(matPASTIX%spm)
         matPASTIX%spm = spm2
      end if
      deallocate (spm2)

      if (lssolver%timing) then
         call cpu_time(tpe)
         call system_clock(cke, clock_rate)
         timing%rlstime2 = timing%rlstime2 + (cke - cks)/real(clock_rate)
         timing%clstime2 = timing%clstime2 + tpe - tps
      end if
   END SUBROUTINE check_mat_PASTIX

   !***********************************************
   ! Analysis of matrix with PASTIX
   !***********************************************
   SUBROUTINE anal_mat_PASTIX(matPASTIX)
      use matrices_tools, only: dump_CSR

      TYPE(PASTIX_STRUC) :: matPASTIX
      integer :: iproc, ierr
      real*8  :: tps, tpe
      integer :: cks, clock_rate, cke
      if (lssolver%timing) then
         call cpu_time(tps)
         call system_clock(cks, clock_rate)
      end if

      call pastix_task_analyze(matPASTIX%pastix_data, matPASTIX%spm, matPASTIX%info)
      if (matPASTIX%info .ne. 0) then
         error stop "pastix_task_analyze failed"
      end if

      if (lssolver%timing) then
         call cpu_time(tpe)
         call system_clock(cke, clock_rate)
         timing%rlstime3 = timing%rlstime3 + (cke - cks)/real(clock_rate)
         timing%clstime3 = timing%clstime3 + tpe - tps
      end if
   END SUBROUTINE anal_mat_PASTIX

   !***********************************************
   ! LU factorization of matrix with PASTIX
   !***********************************************
   SUBROUTINE LU_mat_pastix(matPASTIX)

      TYPE(PASTIX_STRUC) :: matPASTIX

      real*8  :: tps, tpe
      integer :: cks, clock_rate, cke
      if (lssolver%timing) then
         call cpu_time(tps)
         call system_clock(cks, clock_rate)
      end if

      call pastix_task_numfact(matPASTIX%pastix_data, matPASTIX%spm, matPASTIX%info)
      if (matPASTIX%info .ne. 0) then
         error stop "pastix_task_numfact failed"
      end if

      if (lssolver%timing) then
         call cpu_time(tpe)
         call system_clock(cke, clock_rate)
         timing%rlstime5 = timing%rlstime5 + (cke - cks)/real(clock_rate)
         timing%clstime5 = timing%clstime5 + tpe - tps
      end if

   END SUBROUTINE LU_mat_pastix

   !***********************************************
   ! Solve problem with PASTIX
   !***********************************************
   SUBROUTINE solve_mat_PASTIX(matPASTIX)
      TYPE(PASTIX_STRUC) :: matPASTIX

      integer :: Nall, irow, iproc, IERR
      real*8  :: tps, tpe
      integer :: cks, clock_rate, cke
      if (lssolver%timing) then
         call cpu_time(tps)
         call system_clock(cks, clock_rate)
      end if

      call pastix_task_solve_and_refine( matPASTIX%pastix_data, matPASTIX%spm%nexp, matPASTIX%nrhs, matPASTIX%b, matPASTIX%spm%nexp, matPASTIX%x, matPASTIX%spm%nexp, matPASTIX%info )
      ! call pastix_task_solve( matPASTIX%pastix_data, matPASTIX%spm%nexp, matPASTIX%nrhs, matPASTIX%x, matPASTIX%spm%nexp, matPASTIX%info )
      if (matPASTIX%info .ne. 0) then
         error stop "pastix_task_solve failed"
      end if
      matPASTIX%rhs = matPASTIX%x(:, 1)

      if (lssolver%timing) then
         call cpu_time(tpe)
         call system_clock(cke, clock_rate)
         timing%rlstime6 = timing%rlstime6 + (cke - cks)/real(clock_rate)
         timing%clstime6 = timing%clstime6 + tpe - tps
      end if
   END SUBROUTINE solve_mat_PASTIX

   SUBROUTINE terminate_mat_PASTIX()

      matPASTIX%spm%colptr = c_null_ptr
      matPASTIX%spm%rowptr = c_null_ptr
      matPASTIX%spm%values = c_null_ptr
      matPASTIX%spm%loc2glob = c_null_ptr
      call spmExit(matPASTIX%spm)
      deallocate (matPASTIX%spm)

      if (allocated(matPASTIX%x)) then
         deallocate (matPASTIX%x)
         deallocate (matPASTIX%b)
         deallocate (matPASTIX%rhs)
      end if

      !  Destroy the C data structure (Should be last if used to call MPI_Finalize)
      call pastixFinalize(matPASTIX%pastix_data)

   END SUBROUTINE terminate_mat_PASTIX

END MODULE solve_pastix