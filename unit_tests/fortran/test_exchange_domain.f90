PROGRAM test_exchange_domain
  USE globals
  USE Main_utils, ONLY: check_temperature_exchange_domain
  USE solver_status, ONLY: solver_status_t, solve_success, solve_retry
  USE MPI_OMP, ONLY: init_MPI_OMP, MPIvar
  USE mpi, ONLY: MPI_FINALIZE
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan, ieee_positive_inf
  IMPLICIT NONE
  INTEGER :: ierr

  CALL init_MPI_OMP()
  phys%neq=5
  Mesh%Nelems=1
  Mesh%Nnodesperelem=2
  refElPol%Ngauss2d=1
  ALLOCATE(refElPol%N2D(1,2),sol%u(10))
  refElPol%N2D(1,:)=[.5d0,.5d0]

  CALL check_single_species_undershoots()
  CALL check_invalid_states()
  CALL check_nonfinite_states()
  CALL check_gauss_interpolation()
  IF(MPIvar%glob_id==0)WRITE(*,*)'Exchange domain: single-species recovery, both-negative rejection, nonfinite and collective Gauss checks PASS'
  CALL MPI_FINALIZE(ierr)

CONTAINS

  SUBROUTINE check_single_species_undershoots()
    CALL set_uniform_state(state(1.d0,0.d0,2.d0,-.1d0))
    CALL expect_status(solve_success,'negative electron temperature with positive total')
    CALL set_uniform_state(state(1.d0,0.d0,-.1d0,2.d0))
    CALL expect_status(solve_success,'negative ion temperature with positive total')
  END SUBROUTINE check_single_species_undershoots

  SUBROUTINE check_invalid_states()
    CALL expect_collective_rejection(state(0.d0,0.d0,1.d0,1.d0),'zero density')
    CALL expect_collective_rejection(state(-1.d0,0.d0,1.d0,1.d0),'negative density')
    CALL expect_collective_rejection(state(1.d0,0.d0,-.1d0,-.2d0),'both species negative')
    CALL expect_collective_rejection(state(1.d0,0.d0,-1.d0,1.d0),'zero thermal sum')
    CALL expect_collective_rejection(state(1.d0,0.d0,-1.d0,.2d0),'negative thermal sum')
    CALL expect_collective_rejection(state(1.d0,3.d0,-.1d0,-.2d0),'negative thermal sum despite positive ion total energy')
  END SUBROUTINE check_invalid_states

  SUBROUTINE check_nonfinite_states()
    REAL*8 :: U(5)

    U=state(1.d0,0.d0,1.d0,1.d0)
    U(3)=ieee_value(0.d0,ieee_positive_inf)
    CALL expect_collective_rejection(U,'infinite ion energy')
    U=state(1.d0,0.d0,1.d0,1.d0)
    U(4)=ieee_value(0.d0,ieee_quiet_nan)
    CALL expect_collective_rejection(U,'NaN electron energy')
  END SUBROUTINE check_nonfinite_states

  SUBROUTINE expect_collective_rejection(U,label)
    REAL*8, INTENT(IN) :: U(5)
    CHARACTER(*), INTENT(IN) :: label

    CALL set_uniform_state(state(1.d0,0.d0,1.d0,1.d0))
    ! Only one rank is invalid; its rejection must reach every rank.
    IF(MPIvar%glob_id==MPIvar%glob_size-1)CALL set_uniform_state(U)
    CALL expect_status(solve_retry,label)
    CALL set_uniform_state(state(1.d0,0.d0,1.d0,1.d0))
    CALL expect_status(solve_success,'valid state after '//label)
  END SUBROUTINE expect_collective_rejection

  SUBROUTINE check_gauss_interpolation()
    ! Positive nodal thermal energies can extrapolate to a negative Gauss state.
    refElPol%N2D(1,:)=[2.d0,-1.d0]
    sol%u=[1.d0,0.d0,.1d0,.2d0,.1d0,1.d0,0.d0,.1d0,2.d0,.1d0]
    CALL expect_status(solve_retry,'positive nodes with invalid Gauss state')
  END SUBROUTINE check_gauss_interpolation

  FUNCTION state(density,velocity,ion_internal,electron_internal) RESULT(U)
    REAL*8, INTENT(IN) :: density,velocity,ion_internal,electron_internal
    REAL*8 :: U(5)

    U=[density,density*velocity,density*(ion_internal+.5d0*velocity**2),density*electron_internal,.1d0]
  END FUNCTION state

  SUBROUTINE set_uniform_state(U)
    REAL*8, INTENT(IN) :: U(5)
    sol%u=[U,U]
  END SUBROUTINE set_uniform_state

  SUBROUTINE expect_status(expected,label)
    INTEGER, INTENT(IN) :: expected
    CHARACTER(*), INTENT(IN) :: label
    TYPE(solver_status_t) :: status

    CALL check_temperature_exchange_domain(status)
    IF(status%outcome==expected)RETURN
    WRITE(*,*)label,': rank, expected, actual:',MPIvar%glob_id,expected,status%outcome
    ERROR STOP 'exchange source-domain policy'
  END SUBROUTINE expect_status
END PROGRAM test_exchange_domain
