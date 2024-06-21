!*************************************************
! project: TOKAM3X
! file: MPI_OMP.f90
! date: 02/07/2012
! Initialisation of MPI and OpenMP parallelization
!*************************************************

MODULE MPI_OMP
  USE mpi
  IMPLICIT NONE

  ! Type for MPI parallelization info
  TYPE MPI_typ
     INTEGER :: glob_id, glob_size ! global communicator size and id of the local process in the communicator
#ifdef TOR3D
     INTEGER :: ntor, npol, itor, ipol
#endif
  END TYPE MPI_typ

  ! Type for OMP parallelization info
  TYPE OMP_typ
     INTEGER :: Nthreads ! number of openMP threads
  END TYPE OMP_typ

  TYPE(MPI_typ) :: MPIvar
  TYPE(OMP_typ) :: OMPvar

  PUBLIC :: MPIvar, OMPvar

CONTAINS

  !********************************************
  ! Initialization of MPI/OMP
  !********************************************
  SUBROUTINE init_MPI_OMP()
    INTEGER :: IERR, OMP_GET_MAX_THREADS, MPI_THREAD_provided, MPI_THREAD_required

    ! Initialization of the MPI communicator
    !  MPI_THREAD_required = MPI_THREAD_SINGLE
    MPI_THREAD_required = MPI_THREAD_MULTIPLE
    !#ifdef THREAD_FUNNELED
    !   MPI_THREAD_required = MPI_THREAD_FUNNELED
    !#else
    !   MPI_THREAD_required = MPI_THREAD_MULTIPLE
    !#endif
    CALL MPI_init_thread(MPI_THREAD_required, MPI_THREAD_provided, IERR)
    CALL MPI_Comm_size(MPI_COMM_WORLD, MPIvar%glob_size, IERR)
#ifndef PARALL
    IF (MPIvar%glob_size > 1) THEN
       WRITE (6, *) 'Error: using serial version in a parallel environment'
       STOP
    ENDIF
#endif
    CALL MPI_Comm_rank(MPI_COMM_WORLD, MPIvar%glob_id, IERR)
    IF ((MPI_THREAD_provided .NE. MPI_THREAD_required) .AND. (MPIvar%glob_id .EQ. 0)) THEN
       PRINT *, 'Problem with initialization of multi-threaded MPI.'
       PRINT *, 'Required: ', MPI_THREAD_required
       PRINT *, 'Provided: ', MPI_THREAD_provided
       PRINT *, '(MPI_THREAD_SINGLE = ', MPI_THREAD_SINGLE, ', MPI_THREAD_FUNNELED = ', MPI_THREAD_FUNNELED, &
            ', MPI_THREAD_SERIALIZED = ', MPI_THREAD_SERIALIZED, ', MPI_THREAD_MULTIPLE = ', MPI_THREAD_MULTIPLE, ')'
       PRINT *, 'Exiting...'
       STOP
    ENDIF

    ! Makes sure that all the processes in the communicator have the same number of openMP threads
    IF (MPIvar%glob_id .EQ. 0) THEN
       OMPvar%Nthreads = OMP_GET_MAX_THREADS()
    ENDIF
    CALL MPI_BCast(OMPvar%Nthreads, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
    CALL OMP_SET_NUM_THREADS(OMPvar%Nthreads)
#ifdef TOR3D
    MPIvar%ntor = 1
    MPIvar%npol = 1
    MPIvar%itor = 0
    MPIvar%ipol = 0
#endif
  END SUBROUTINE init_MPI_OMP

#ifdef TOR3D
#ifdef PARALL
  SUBROUTINE set_divisions
    USE globals
    MPIvar%ntor = numer%npartor
    IF (MPIvar%glob_size .LT. MPIvar%ntor) THEN
       WRITE (6, *) "Error: wrong number of MPI toroidal partition in input file or wrong number of MPI processes"
       WRITE (6, *) "Number of processes: ", MPIvar%glob_size
       WRITE (6, *) "Number of MPI toroidal partitions: ", MPIvar%ntor
       STOP
    ENDIF
    IF (MOD(MPIvar%glob_size, MPIvar%ntor) .NE. 0) THEN
       WRITE (6, *) "Error: the number of MPI processes must be a multiple of the number of MPI toroidal divisions (set in input file)"
       WRITE (6, *) "Number of processes: ", MPIvar%glob_size
       WRITE (6, *) "Number of MPI toroidal partitions: ", MPIvar%ntor
       STOP
    ENDIF

    MPIvar%npol = MPIvar%glob_size/MPIvar%ntor
    MPIvar%itor = 1 + MPIvar%glob_id/MPIvar%npol
    MPIvar%ipol = 1 + MOD(MPIvar%glob_id, MPIvar%npol)
  END SUBROUTINE set_divisions
#endif
#endif

END MODULE MPI_OMP
