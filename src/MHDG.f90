PROGRAM MHDG
  USE Main_utils
  USE MPI_OMP

  IMPLICIT NONE

  WRITE (6, *) "STARTING"

  ! check initial arguments like mesh, initial solution etc
  CALL check_input_arguments()

  ! Start timing the code
  CALL cpu_TIME(time_start)
  CALL system_CLOCK(clock_start, clock_rate)

  ! Initialize MPI
  CALL init_MPI_OMP()
  
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

  CALL read_splines()

  ! Linear solver: set the start to true
  matK%start = .TRUE.
  IF (lssolver%timing) THEN
     CALL init_solve_timing
  ENDIF

  ! Initialize marked elements for thresholds
  ALLOCATE (mkelms(Mesh%Nelems))
  mkelms = .FALSE.

  ! Create the reference element based on the mesh type
  CALL mpi_barrier(MPI_COMM_WORLD,ierr)
  CALL create_reference_element(refElPol, 2, verbose = 1)

  ! create the temp.msh and temp.mesh needed by the adaptivity
  IF((switch%readMeshFromSol) .AND. (adapt%adaptivity)) THEN
     CALL generate_msh_from_solution_mesh('./res/temp.msh')
     CALL convert_msh2mesh('./res/temp')
  ENDIF

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
  CALL mesh_preprocess_serial(ierr)

  IF((ierr .EQ. 0) .AND. (switch%read_gmsh)) THEN
     CALL free_mesh
     IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
        CALL load_gmsh_mesh(mesh_name, 1)
     ELSE
        CALL load_gmsh_mesh(mesh_name, 0)
     ENDIF
     CALL mesh_preprocess_serial(ierr)

     IF(ierr .EQ. 0) THEN
        WRITE(*,*) "Problem in mesh_preprocess. STOP."
        STOP
     ENDIF
  ENDIF

  ! Initialize magnetic field (the Mesh is needed)
  CALL initialize_magnetic_field()

  ! Load magnetic field and, if ohmic src, also Jtor
  CALL load_magnetic_field_Jtor()

  ! set parameters like Nelems, Nfacenodes etc
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

  ! Initialize the solution
  CALL initialize_solution()

  ! add initial perturbation (pertini), magnetic or blob
  CALL add_initial_perturbation()

  ! Initialise puff, only if neutrals are present
  CALL initialize_puff()

  ! Save solution
  CALL setSolName(save_name, mesh_name, 0, .TRUE., .FALSE.)
  CALL HDF5_save_solution(save_name)

  CALL mpi_barrier(mpi_comm_world,ierr)

  ! Allocate and initialize uiter, uiter_best, qiter_best, u0, u_conv, q_conv
  CALL initialize_solu0_uiter_uconv()

  errNR_adapt = 1e10
  ir_adapt = 0
  ir_check = 0

  it0 = 1
  IF (switch%ME .AND. time%it .NE. 0) it0 = time%it
  IF (switch%ME) THEN
     time%dt = time%dt_ME/simpar%refval_time
  ENDIF
  dt0 = time%dt

  !*******************************************************
  !                  TIME LOOP
  !*******************************************************
  DO it = it0, nts ! ************ TIME LOOP *********************

     ! if a new time step starts, it means that the solution converged sol%u_conv = sol%u, sol%q_conv = sol%q
     CALL update_uconv_qconv(sol%u, sol%q)

     ! Actualization of time
     time%t = time%t + time%dt
     time%it = time%it + 1
     time%ik = time%ik + 1
     sol%Nt = sol%Nt + 1

     IF (MPIvar%glob_id .EQ. 0) THEN
        WRITE (6, '(" *", 60("*"), "**")')
        WRITE (6, '(" *", 20X,    "Time iteration   = ", I5, 16X, " *")') time%it
        WRITE (6, '(" *", 60("*"), "**")')
     END IF

     !*******************************************************
     !             Newton-Raphson iterations
     !*******************************************************
     ! uiter = sol%u0
     CALL update_uiter()

     ir = 1

     DO WHILE(ir .LE. numer%nrp) ! ************ NEWTON-RAPHSON LOOP *********************

        ! update nonconstant dumping factor
        CALL update_dumpnr()

        IF (MPIvar%glob_id .EQ. 0) THEN
           WRITE (6, *) "***** NR iteration: ", ir, "*****"
           WRITE (6, *) "NR dumping factor:  ",  numer%dumpnr
        ENDIF

        ! Compute Jacobian
        CALL HDG_computeJacobian()
        ! Set boundary conditions
        CALL hdg_BC()
        ! Compute elemental mapping
        CALL hdg_Mapping()
        ! Assembly the global matrix
        CALL hdg_Assembly()

        !  WRITE (6, *) "Save matrix"
        !  call HDF5_save_CSR_matrix('Mat')
        !  call HDF5_save_CSR_vector('rhs')
        !  stop
        !  call displayMatrixInt(Mesh%F)
        !  call displayMatrixInt(Mesh%extfaces)
        !  call displayVectorInt(Mesh%periodic_faces)
        !  stop
        !  if (ir==10) then
        !   call print_matrices_hdf5
        !   stop
        !  endif

        ! Solve linear system
        CALL solve_global_system(ir)

        ! Compute element-by-element solution
        CALL compute_element_solution()

        ! Check for NaN (should work with optimization flags)
        CALL check_for_NaNs()

        IF (adapt%adaptivity .and. restart_adapt) THEN
          CALL adaptivity
          DEALLOCATE(uiter)
          ALLOCATE(uiter(size(sol%u)))
          uiter = 0.
        ENDIF

        ! Apply threshold
        ! CALL HDG_applyThreshold(mkelms)

        ! Apply filtering
        ! CALL HDG_FilterSolution()

        ! Compute error on oscillations, print max value of oscillation and save solution as check-point if oscillations are lower than threshold
        CALL compute_error_oscillations(error_oscillation, oscillations, min_osc, max_osc, n_osc, ir, ir_check, Mesh_prec)

        ! Save solution
        IF (switch%saveNR) THEN
           CALL setSolName(save_name, mesh_name, ir, .FALSE., .TRUE.)
           CALL HDF5_save_solution(save_name)
        END IF

        ! Check convergence of Newton-Raphson
        errNR = computeResidual(sol%u, uiter, 1.)
        errNR = errNR/numer%dumpnr

        IF (MPIvar%glob_id .EQ. 0) THEN
           WRITE (*, *)   "Error:                  ", errNR
#ifdef WITH_PETSC
           IF(lssolver%sollib .EQ. 3) THEN
              WRITE (*, *) "Relative Residue PETSc: ", matPETSC%residue
              WRITE (*,*)  "Number of Iterations:   ", matPETSC%its
              WRITE (*,*)  "Converged Reason:       ", matPETSC%convergedReason
           END IF
#endif
        ENDIF

        IF (errNR .LT. numer%tNR) THEN
           ! Save check-point solution if NR error is smaller than threshold
           WRITE(*,*) "Solution saved as checkpoint."
           CALL update_uconv_qconv(uiter_best, qiter_best)
           errNR_adapt = 1e10
           ir_adapt = 0
           ir_check = 1
           CALL deep_copy_mesh_struct(Mesh, Mesh_prec)
           EXIT
        ELSEIF (errNR .GT. numer%div) THEN
           WRITE (6, *) 'Problem in the N-R procedure'
           STOP
        ELSE
           uiter = sol%u
           !! ADAPTIVITY
           IF(errNR .LT. errNR_adapt) THEN

              errNR_adapt = errNR
              ir_adapt = ir

              ! if the NR is the lowest reached so far, then save it as best check-point
              WRITE(*,*) "Solution saved as last checkpoint."
              CALL update_uiter_qiter_best(uiter_best, qiter_best, sol%u, sol%q)
              divergence_counter_adapt = 0

           ELSEIF(errNR .GT. errNR_adapt) THEN
              divergence_counter_adapt = divergence_counter_adapt + 1
           ENDIF

           WRITE(*,*) "ir_check: ", ir_check
           ! Call adaptivity if one of the following conditions is respected
           IF ((adapt%adaptivity) .AND. (ir .GT. (ir_check+1)) .AND. ((adapt%osc_adapt .AND. (MAXVAL(oscillations) .GT. adapt%osc_tol)) .OR. ((adapt%NR_adapt) .AND. (MOD(ir,adapt%freq_NR_adapt) .EQ. 0))))THEN !  .or. (flag)) THEN

              WRITE(*,*) "Residue before mapping:", computeResidual(sol%u, uiter, 1.)/numer%dumpnr

              ! call adaptivity precedure
              CALL adaptivity()

              ! u0 also needs to be projected from old mesh to new mesh
              CALL project_u0_newmesh()

              ! update uiter to new mapped solution
              CALL update_uiter()

              IF(ir_check .NE. numer%nrp) THEN
                 ir = ir_check
              ELSE
                 ir = 0
              ENDIF

              WRITE(*,*) "Residue after mapping:", computeResidual(sol%u, uiter, 1.)/numer%dumpnr

              CALL deep_copy_mesh_struct(Mesh, Mesh_prec)

           ENDIF
        END IF
        IF (MPIvar%glob_id .EQ. 0) THEN
           WRITE (6, *) "*********************************"
           WRITE (6, *) " "
           WRITE (6, *) " "
        ENDIF
        ir = ir + 1
     END DO ! ************ END OF NEWTON-RAPHSON LOOP *********************

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
        IF (errlstime .LT. numer%tTM) THEN
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

              CALL deep_copy_mesh_struct(Mesh, Mesh_prec)

              ! call the adaptive procedure if time refinement is on
              IF((adapt%adaptivity) .AND. (adapt%time_adapt) .AND. (MOD(it,adapt%freq_t_adapt) .EQ. 0)) THEN

                 CALL adaptivity

                 CALL update_u0(sol%u)

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

              CALL adaptivity

              CALL update_u0(sol%u)

              CALL deep_copy_mesh_struct(Mesh, Mesh_prec)
           ENDIF

           ! if moving equilibrium case then update the magnetic field, otherwise just continue
           IF(switch%ME) THEN
              ! ReLoad magnetic field and Jtor
              CALL load_magnetic_field()
              CALL loadJtorMap()
              CALL SetPuff()
              time%dt = time%dt_ME/simpar%refval_time
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
  CALL free_main()

END PROGRAM MHDG
