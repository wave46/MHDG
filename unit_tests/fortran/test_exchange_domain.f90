PROGRAM test_exchange_domain
  USE globals
  USE Main_utils, ONLY: check_temperature_exchange_domain, check_newton_convergence, &
       newton_candidate_is_acceptable, handle_newton_result, newton_result_t, solve_nonlinear_system, &
       physical_iterate_valid, seed_divergence_checkpoint, divergence_checkpoint_available, &
       divergence_refinements, uiter, uiter_best, qiter_best, errNR, nb_args
  USE solver_status, ONLY: solver_status_t, solve_success, solve_retry
  USE MPI_OMP, ONLY: init_MPI_OMP, MPIvar
  USE mpi, ONLY: MPI_FINALIZE
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan, ieee_positive_inf
  IMPLICIT NONE
  INTEGER :: ierr

  CALL init_MPI_OMP()
  phys%neq=5
  phys%idx_rhon_eq=5
  Mesh%Nelems=1
  Mesh%Nnodesperelem=2
  refElPol%Ngauss2d=1
  ALLOCATE(refElPol%N2D(1,2),sol%u(10),sol%q(20),sol%u_tilde(10),sol%u0(10,1),uiter(10))
#ifdef PARALL
  ALLOCATE(Mesh%ghostElems(1))
  Mesh%ghostElems=0
#endif
  sol%q=0.d0;sol%u_tilde=1.d0
  adapt%check_neutral_positivity=.TRUE.
  numer%dumpnr=1.d0;numer%dumpnr_min=1.d0;numer%dumpnr_max=1.d0
  numer%dumpnr_n0=0.d0;numer%dumpnr_width=1.d0
  numer%nrp=2;numer%tNR=1.d-4;numer%div=1.d3
  switch%saveNR=.FALSE.;switch%transport_1d=.FALSE.
  refElPol%N2D(1,:)=[.5d0,.5d0]

  CALL check_single_species_undershoots()
  CALL check_invalid_states()
  CALL check_nonfinite_states()
  CALL check_bootstrap_policy()
  CALL check_gauss_interpolation()
  CALL check_driver_domain_failure()
  IF(MPIvar%glob_id==0)WRITE(*,*)'Exchange domain: undershoots, bootstrap policy, Gauss checkpoints and Newton recovery routing PASS'
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
    TYPE(solver_status_t) :: status

    adapt%adaptivity=.TRUE.;adapt%div_adapt=.TRUE.
    refElPol%N2D(1,:)=[.5d0,.5d0]
    CALL set_uniform_state(state(1.d0,0.d0,1.d0,1.d0))
    CALL seed_divergence_checkpoint()
    IF(.NOT.divergence_checkpoint_available)ERROR STOP 'valid Gauss checkpoint rejected'
    ! Positive nodal thermal energies can extrapolate to a negative Gauss state.
    refElPol%N2D(1,:)=[2.d0,-1.d0]
    sol%u=[1.d0,0.d0,.1d0,.2d0,.1d0,1.d0,0.d0,.1d0,2.d0,.1d0]
    CALL expect_status(solve_retry,'positive nodes with invalid Gauss state')
    CALL seed_divergence_checkpoint()
    IF(divergence_checkpoint_available)ERROR STOP 'unevaluable Gauss state seeded checkpoint'
    CALL check_temperature_exchange_domain(status)
    IF(status%outcome/=solve_retry)ERROR STOP 'invalid candidate bypassed source-domain gate'
    refElPol%N2D(1,:)=[.5d0,.5d0]
  END SUBROUTINE check_gauss_interpolation

  SUBROUTINE check_bootstrap_policy()
    TYPE(solver_status_t) :: status
    TYPE(newton_result_t) :: result
    REAL*8 :: retained(10)

    time%legacy_input=.TRUE.;switch%steady=.FALSE.;switch%psdtime=.FALSE.
    switch%ME=.FALSE.;nb_args=1
    adapt%adaptivity=.FALSE.;adapt%div_adapt=.FALSE.
    CALL set_uniform_state(state(1.d0,0.d0,2.d0,-.1d0))
    uiter=sol%u
    CALL check_newton_convergence(status)
    IF(status%outcome/=solve_success.OR.errNR/=0.d0)ERROR STOP 'evaluable zero-update trial rejected'
    IF(physical_iterate_valid)ERROR STOP 'disabled recovery hid negative temperature'
    IF(.NOT.newton_candidate_is_acceptable(physical_iterate_valid))ERROR STOP 'time_init continuation rejected'
    nb_args=2;switch%steady=.TRUE.
    IF(.NOT.newton_candidate_is_acceptable(.FALSE.))ERROR STOP 'diffred continuation rejected'
    switch%steady=.FALSE.;switch%psdtime=.TRUE.
    IF(.NOT.newton_candidate_is_acceptable(.FALSE.))ERROR STOP 'pseudo-time continuation rejected'
    retained=sol%u
    result%status=solver_status_t(solve_retry,0,'Newton iteration limit reached')
    result%iteration_limit_reached=.TRUE.
    CALL handle_newton_result(result)
    IF(ANY(sol%u/=retained))ERROR STOP 'budget output changed raw state'
    IF(result%status%outcome/=solve_retry)ERROR STOP 'budget output relabelled as converged'
    switch%psdtime=.FALSE.;switch%ME=.TRUE.;nb_args=1
    IF(newton_candidate_is_acceptable(.FALSE.))ERROR STOP 'moving equilibrium allowed bootstrap exception'
    time%legacy_input=.FALSE.;switch%steady=.TRUE.
    IF(newton_candidate_is_acceptable(.FALSE.))ERROR STOP 'physical timestep allowed bootstrap exception'
    IF(.NOT.newton_candidate_is_acceptable(.TRUE.))ERROR STOP 'admissible physical state rejected'
  END SUBROUTINE check_bootstrap_policy

  SUBROUTINE check_driver_domain_failure()
    TYPE(newton_result_t) :: result
    REAL*8 :: rejected(10)

    time%legacy_input=.TRUE.;switch%steady=.TRUE.;switch%ME=.FALSE.
    adapt%adaptivity=.TRUE.;adapt%div_adapt=.TRUE.
    adapt%max_divergence_refinements=4
    divergence_refinements=0;divergence_checkpoint_available=.FALSE.
    CALL set_uniform_state(state(1.d0,0.d0,1.d0,1.d0))
    IF(MPIvar%glob_id==MPIvar%glob_size-1)CALL set_uniform_state(state(1.d0,0.d0,-1.d0,.2d0))
    sol%u0(:,1)=sol%u
    rejected=sol%u
    physical_iterate_valid=.TRUE.
    CALL solve_nonlinear_system(result)
    IF(result%status%outcome/=solve_retry.OR.result%iterations/=1)ERROR STOP 'Newton domain rejection lost'
    IF(physical_iterate_valid)ERROR STOP 'assembly failure bypassed domain recovery route'
    IF(divergence_refinements/=0)ERROR STOP 'missing checkpoint consumed a refinement'
    IF(ANY(sol%u/=rejected))ERROR STOP 'failed domain recovery changed current state'
    IF(result%iteration_limit_reached)ERROR STOP 'domain failure mislabeled as budget output'
    physical_iterate_valid=.TRUE.
    CALL solve_nonlinear_system(result,predicted=.TRUE.)
    IF(result%status%outcome/=solve_retry)ERROR STOP 'invalid prediction did not request fallback'
    IF(physical_iterate_valid.OR.divergence_refinements/=0)ERROR STOP 'prediction attempted mesh recovery'
  END SUBROUTINE check_driver_domain_failure

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
