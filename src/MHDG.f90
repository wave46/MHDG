PROGRAM MHDG
  USE Main_utils
  USE MPI_OMP

  IMPLICIT NONE
  LOGICAL :: time_finished
  TYPE(newton_result_t) :: newton_result

  ! Initialize MPI
  CALL init_MPI_OMP()

  ! check initial arguments like mesh, initial solution etc
  CALL check_input_arguments()

  ! Start timing the code
  CALL cpu_TIME(time_start)
  CALL system_CLOCK(clock_start, clock_rate)

  IF (MPIvar%glob_id .EQ. 0) THEN
     WRITE (6, *) "STARTING"
  ENDIF

  ! Number of threads
  Nthreads = OMPvar%Nthreads

  IF (MPIvar%glob_id .EQ. 0) THEN
     WRITE(6,*) "Using ", Nthreads, " threads"
  ENDIF

  ! Read input file param.txt
  CALL read_input()

#ifdef WITH_PETSC
  IF (lssolver%sollib .EQ. 3) THEN
     CALL InitPETSC()
  ENDIF
#endif

  ! Set parallelization division in toroidal and poloidal plane
#ifdef TOR3D
#ifdef PARALL
  CALL set_divisions()
#endif
#endif

  ! Initialization of the simulation parameters !TODO check if I need to pass dt to init_sim
  CALL init_sim(nts, dt)

  CALL load_mesh()


  ! Linear solver: set the start to true
  matK%start = .TRUE.
  IF (lssolver%timing) THEN
     CALL init_solve_timing()
  ENDIF

  ! Initialize marked elements for thresholds
  ALLOCATE (mkelms(Mesh%Nelems))
  mkelms = .FALSE.

  ! Create the reference element based on the mesh type
#ifdef PARALL
  CALL mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
  CALL create_reference_element(refElPol, 2, verbose = 1)

#ifdef TOR3D
  ! create toroidal reference element and toroidal structures
  CALL create_reference_element(refElTor, 1, numer%ptor, verbose = 1)
  CALL create_toroidal_structures(refElTor, refElPol)
  ! Define toroidal discretization
  CALL define_toroidal_discretization()
#endif


  ! Mesh preprocess: create the mesh related structures
  ! used in the HDG scheme
  ierr = 1
#ifdef PARALL
  ! This is still the global mesh; ghost-face data is created during decomposition.
  CALL mesh_preprocess_serial(ierr)
#else
  CALL mesh_preprocess(ierr)
#endif


  IF((ierr .EQ. 0)) THEN
     CALL free_mesh
     IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
        CALL load_gmsh_mesh(mesh_input_name, 1)
     ELSE
        CALL load_gmsh_mesh(mesh_input_name, 0)
     ENDIF

     CALL mesh_preprocess_serial(ierr)

     IF(ierr .EQ. 0) THEN
        WRITE(*,*) "Problem in mesh_preprocess. STOP."
        STOP
     ENDIF
  ENDIF

#ifdef PARALL
  ! if the restart solution is given
  IF (nb_args .EQ. 2) THEN
     ! load the serial solution
     CALL initialize_solution()
     ! decompose the domain over the processes
     CALL domain_decomposition()
     ! decompose the solution over the processes
     CALL solution_decomposition()
  ELSE ! if the restart solution is not given
     ! just decompose the domain
     CALL domain_decomposition()
  ENDIF
#endif

  ! MPI restarts were already loaded before solution decomposition, even with one rank.
  ! Read serial restart data before deciding whether another field is needed.
#ifndef PARALL
  IF (nb_args .EQ. 2) CALL initialize_solution()
#endif

  CALL exit_if_restart_complete()

  CALL initialize_magnetic_field()
  CALL load_magnetic_field_Jtor()
  IF (nb_args .EQ. 1) CALL initialize_solution()

  IF (switch%external_heating) THEN
     ! Initialise external heating
     CALL initialize_external_heating()
     ! Load external heating
     CALL load_external_heating()
  ENDIF
  

  restart_adapt = adapt%rest_adapt
  CALL set_parameters()

#ifdef PARALL
  ! initialise comunication
  CALL init_com()
#endif

  ! Allocation and initialization of the elemental matrices
  CALL init_elmat()

  ! Compute first equation (definition of the gradient)
  CALL HDG_precalculatedfirstequation()

  ! Initialize shock capturing
  IF ((switch%shockcp .GT. 0) .OR. (adapt%shockcp_adapt .GT. 0))  THEN
     CALL initializeShockCapturing()
  ENDIF

  ! add initial perturbation (pertini), magnetic or blob
  CALL add_initial_perturbation()

  ! Initialise puff, only if neutrals are present
  CALL initialize_puff()

  ! Read 1D diffusion profiles if needed
  IF (switch%import_diffusion_1D) THEN
     CALL read_1D_diffusion_profiles()
  ENDIF

  ! Save solution
  IF (switch%transport_1d) THEN
     CALL update_reduced_transport_profiles()
  ENDIF
  CALL setSolName(save_name, mesh_name, 0, .TRUE., .FALSE.)
  CALL HDF5_save_solution(save_name)

  ! Allocate and initialize uiter, uiter_best, qiter_best, u0, u_conv, q_conv
  CALL initialize_solu0_uiter_uconv()

  errNR_adapt = 1e10
  ir_adapt = 0
  ir_check = 0

  it0 = 1
  ! For ME, nts is the final global step, including on restart.
  IF (switch%ME .OR. .NOT. time%legacy_input) it0 = time%it + 1
  dt0 = time%dt

  !*******************************************************
  !                  TIME LOOP
  !*******************************************************
  it = it0
  DO ! ************ TIME LOOP *********************
     CALL advance_time_step(newton_result,time_finished)
     IF (time_finished) EXIT
     CALL handle_newton_result(newton_result)

     !  ! Apply threshold
     !  CALL HDG_applyThreshold()

     ! Check convergence in time advancing and update

     errlstime = computeResidual(sol%u, sol%u0(:, 1), time%dt)
     sol%tres(sol%Nt) = errlstime
     sol%time(sol%Nt) = time%t

     ! Display results
     IF (MOD(time%it, utils%freqdisp) .EQ. 0) THEN
        CALL displayResults()
     END IF

     ! Check for NaN (doesn't work with optimization flags)
     CALL check_for_NaNs()

     IF (.NOT. switch%steady) THEN
        ! Save solution
        IF (MOD(time%it, utils%freqsave) .EQ. 0) THEN
           CALL setSolName(save_name, mesh_name, time%it, .TRUE., .FALSE.)
           CALL HDF5_save_solution(save_name)
        END IF

        !****************************************
        ! Check steady state and update or exit
        !****************************************
        IF (time%legacy_input .AND. .NOT. switch%ME .AND. errlstime .LT. numer%tTM) THEN
           IF (switch%psdtime) THEN
              ! Pseudo-time simulation

              ! Save solution
              CALL setSolName(save_name, mesh_name, it, .TRUE., .TRUE.)
              CALL HDF5_save_solution(save_name)

              ! reduce diffusion
              CALL reduce_diffusion()

              ! update sol%u0
              CALL update_solution()
              time%it = 0

              ! update sol%u_conv = sol%u, sol%q_conv = sol%q
              CALL update_uconv_qconv(sol%u, sol%q)
              CALL free_mesh_loc(Mesh_prec)
              CALL deep_copy_mesh_struct(Mesh, Mesh_prec)

              ! call the adaptive procedure if time refinement is on
              IF((adapt%adaptivity) .AND. (adapt%time_adapt) .AND. (MOD(it,adapt%freq_t_adapt) .EQ. 0)) THEN

                 CALL adaptivity

                 CALL update_u0(sol%u)
                 CALL free_mesh_loc(Mesh_prec)
                 CALL deep_copy_mesh_struct(Mesh, Mesh_prec)
              ENDIF

              ! compute dt
           ELSE
              ! Time advancing simulation
              IF (MPIvar%glob_id .EQ. 0) THEN
                 WRITE (6, *) "**********************"
                 WRITE (6, *) "Time scheme converged!"
                 WRITE (6, *) "**********************"
              END IF
              EXIT ! Here I exit the time advancing scheme if I reach convergence
           END IF
        ELSEIF (errlstime .GT. numer%div) THEN
           WRITE (6, *) 'Problem in the time advancing scheme'
           STOP
        ELSE

           ! update u0
           CALL update_solution()

           ! call the adaptive procedure if time refinement is on
           IF((adapt%adaptivity) .AND. (adapt%time_adapt) .AND. (MOD(it,adapt%freq_t_adapt) .EQ. 0)) THEN

              ! Periodic refinement projects the completed timestep.
              CALL update_uconv_qconv(sol%u, sol%q)
              CALL adaptivity

              CALL update_u0(sol%u)
              CALL free_mesh_loc(Mesh_prec)
              CALL deep_copy_mesh_struct(Mesh, Mesh_prec)
           ENDIF

        END IF
     ELSE
        IF (switch%psdtime) THEN
           ! Save solution
           CALL setSolName(save_name, mesh_name, it, .TRUE., .TRUE.)
           CALL HDF5_save_solution(save_name)

           CALL reduce_diffusion()

           CALL update_solution()
           time%it = 0
        ELSE
           EXIT
        ENDIF
     ENDIF

     IF(checkpoint .GT. 1) THEN
        convergence_counter = convergence_counter + 1
     ENDIF

     it = it + 1
  END DO ! ************ END OF THE TIME LOOP *********************

  ! Save solution
  CALL setSolName(save_name, mesh_name, time%it, .TRUE., .TRUE.)
  CALL HDF5_save_solution(save_name)

  CALL cpu_TIME(time_finish)
  CALL system_CLOCK(clock_end, clock_rate)
  PRINT '("Elapsed cpu-time = ",f10.3," seconds.")', time_finish - time_start
  PRINT '("Elapsed run-time = ",f10.3," seconds.")', (clock_end - clock_start)/REAL(clock_rate)

  ! Print timing infos
  CALL print_timing_infos()

  IF (switch%testcase < 5) THEN
     ALLOCATE (L2err(phys%neq))
     CALL computeL2ErrorAnalyticSol(L2err)
     WRITE (6, *) " "
     DO i = 1, phys%Neq
        WRITE (6, '(A,I1,A,ES16.5)') "L2 error in U(", i, ") = ", L2err(i)
     END DO
     DEALLOCATE (L2err)
  END IF

  ! free all main allocatables and global variables
  CALL close_heating_input()
  CALL close_magnetic_input()
  CALL clear_control_curves()
  CALL free_main()

END PROGRAM MHDG
