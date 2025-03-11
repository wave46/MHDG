!*****************************************
! project: MHDG
! file: Main_utils.f90
! date: 25/11/2024
! Set of functions and subroutine used in MHDG.f90

MODULE Main_utils

  USE in_out
  USE GMSH_io_module, ONLY: load_gmsh_mesh, HDF5_save_mesh, read_splines, convert_gmsh_to_hdf5, gmsh_mesh2d_write, hdf5_save_mesh_struct
  USE reference_element
  USE preprocess
  USE MPI_OMP
  USE printutils
  USE debug
  USE initialization
#ifdef WITH_PETSC
   USE solve_petsc, only: matPETSC, InitPETSC
#endif
  USE adaptivity_common_module
  USE adaptivity_estimator_module
  USE adaptivity_indicator_module
  USE adaptivity_estimator_indicator_module
  USE Postprocess, only: computeL2ErrorAnalyticSol
#ifdef PARALL
  USE Communications
  USE domain_decomposition_module
#endif
  USE HDG_LimitingTechniques

  IMPLICIT NONE

  INTEGER                      :: Np, Nel, Nfp, Nf, Ndim, Nthreads
  INTEGER                      :: it, ir, ir_check, it0, nts, nu, nut, nb_args, IERR, k, i,is, count_adapt = 0, ir_adapt = 0, checkpoint = 1, divergence_counter_adapt = 0, convergence_counter = 1, order

  LOGICAL, ALLOCATABLE         :: mkelms(:)
  REAL*8                       :: dt, dt0, errNR, errNR_adapt, errlstime
  REAL*8, ALLOCATABLE          :: uiter(:), L2err(:),u0_temp(:,:)
  REAL*8, POINTER              :: uiter_best(:) => NULL(), qiter_best(:) => NULL()
  CHARACTER(LEN=1024)          :: mesh_name,mesh_name_proj, save_name
  CHARACTER ( len = 255 )      :: gmsh_filename
  CHARACTER ( len = 255 )      :: gmsh_filename_mesh, h5_filename
  CHARACTER ( len = 50 )       :: count_adapt_char
  REAL*8                       :: cputtot, runttot
  INTEGER                      :: OMP_GET_MAX_THREADS
  INTEGER                      :: clock_rate,clock_start, clock_end
  REAL*8                       :: time_start, time_finish
  REAL*8,ALLOCATABLE           :: xs(:,:)
  REAL*8, ALLOCATABLE          :: oscillations(:)
  REAL*8                       :: max_osc = -100., min_osc = -100.
  INTEGER,ALLOCATABLE          :: vector_nodes_unique(:,:)
  INTEGER                      :: N_n_vertex, n_osc
  INTEGER*8                    :: file_id
  TYPE(Mesh_type)              :: Mesh_prec
  LOGICAL                      :: restart_adapt
#ifdef PARALL
  TYPE(Mesh_type)              :: Mesh_glob
  REAL*8, POINTER              :: X_glob(:,:) => NULL(), B_glob(:,:) => NULL()
  REAL*8, POINTER              :: u_glob(:) => NULL(), q_glob(:) => NULL(), magnetic_flux_glob(:) => NULL(), u_tilde_glob(:) => NULL()
  INTEGER, POINTER             :: T_glob(:,:) => NULL(), Tb_glob(:,:) => NULL(), F_glob(:,:) => NULL(), N_glob(:,:) => NULL(), flipface_glob(:,:) => NULL(), Tlin_glob(:,:) => NULL(), &
                                  & boundaryFlag_glob(:) => NULL(), extFaces_Glob(:,:) => NULL(), intFaces_Glob(:,:) => NULL()
  CHARACTER(70)                :: npr, nid
#endif
#ifdef TOR3D
  INTEGER                      :: Neq, ntorloc,Nfl, Np1dPol, Np1dTor, Ng1dPol, Ng1dTor,
#endif
  TYPE(Reference_element_type) :: refElPol_prec


CONTAINS

#ifdef PARALL
  SUBROUTINE domain_decomposition()
    IF(MPIvar%glob_size .GT. 1) THEN
       ! split the mesh
       
      CALL free_mesh_loc(Mesh_glob)
      CALL deep_copy_mesh_struct(Mesh, Mesh_glob)
      Mesh_glob%X = Mesh_glob%X*phys%lscale
      CALL split_mesh(MPIvar%glob_size, 3 , .FALSE.)
      CALL mesh_preprocess(ierr)

      WRITE (nid, *) MPIvar%glob_id + 1
      WRITE (npr, *) MPIvar%glob_size
      h5_filename = TRIM(ADJUSTL(mesh_name)) // '_' // TRIM(ADJUSTL(nid)) // '_' // TRIM(ADJUSTL(npr)) // '.h5'
      Mesh%X = Mesh%X*phys%lscale
      CALL HDF5_save_mesh_struct(Mesh, h5_filename)
      Mesh%X = Mesh%X/phys%lscale
      !CALL HDF5_save_mesh(h5_filename, Mesh%Ndim, mesh%Nelems, mesh%Nextfaces, mesh%Nnodes, mesh%Nnodesperelem, mesh%Nnodesperface, mesh%elemType, mesh%T, mesh%X, mesh%Tb, mesh%boundaryFlag)
       
    ELSE
       ALLOCATE(Mesh%ghostelems(Mesh%Nelems))
       ALLOCATE(Mesh%ghostfaces(Mesh%Nfaces))
       ALLOCATE(Mesh%ghostpro(Mesh%Nfaces))
       ALLOCATE(Mesh%ghostloc(Mesh%Nfaces))
       ALLOCATE(Mesh%loc2glob_el(Mesh%Nelems))
       ALLOCATE(Mesh%loc2glob_fa(Mesh%Nfaces))
       ALLOCATE(Mesh%loc2glob_nodes(Mesh%Nnodes))

       Mesh%ghostelems = 0
       Mesh%ghostfaces = 0
       Mesh%ghostpro = 0
       Mesh%ghostLoc = 0
       Mesh%loc2glob_el = [(i, i = 1, Mesh%Nelems)]
       Mesh%loc2glob_fa = [(i, i = 1, Mesh%Nfaces)]
       Mesh%loc2glob_nodes = [(i, i = 1, Mesh%Nnodes)]

       Mesh%Nel_glob = Mesh%Nelems
       Mesh%Nfa_glob = Mesh%Nfaces
       Mesh%Nno_glob = Mesh%Nnodes
    ENDIF

  ENDSUBROUTINE domain_decomposition

  SUBROUTINE solution_decomposition()
    REAL*8, POINTER         :: u_glob(:), u_tilde_glob(:), q_glob(:)
    REAL*8, ALLOCATABLE     :: u_3d(:,:,:), u_tilde_3d(:,:,:), q_4D(:,:,:,:)
    REAL*8, ALLOCATABLE     :: u_3d_local(:,:,:), u_tilde_3d_local(:,:,:), q_4d_local(:,:,:,:)
    INTEGER                 :: i, index

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) "*************************************************"
       WRITE(*,*) "           SOLUTION DECOMPOSITION                "
       WRITE(*,*) "*************************************************"
    ENDIF

    u_glob => sol%u
    u_tilde_glob => sol%u_tilde
    q_glob => sol%q

    ALLOCATE(u_3d(Mesh%Nel_glob, Mesh%Nnodesperelem, phys%neq))
    ALLOCATE(u_tilde_3d(Mesh%Nfa_glob, Mesh%Nnodesperface, phys%neq))
    ALLOCATE(q_4d(Mesh%Nel_glob, Mesh%Nnodesperelem, phys%neq, Mesh%Ndim))
    u_tilde_3d = 0.
    u_3d = 0.
    q_4D = 0.

    CALL reshape_transpose_permute(u_glob, u_3d, phys%neq, Mesh%Nel_glob, Mesh%Nnodesperelem)
    CALL reshape_transpose_permute(u_tilde_glob, u_tilde_3d, phys%neq, Mesh%Nfa_glob, Mesh%Nnodesperface)
    CALL reshape_transpose_permute_4D(q_glob, q_4D, Mesh%Ndim, phys%neq, Mesh%Nel_glob, Mesh%Nnodesperelem)

    ALLOCATE(u_3d_local(Mesh%Nelems, Mesh%Nnodesperelem, phys%neq))
    ALLOCATE(u_tilde_3d_local(Mesh%Nfaces, Mesh%Nnodesperface, phys%neq))
    ALLOCATE(q_4d_local(Mesh%Nelems, Mesh%Nnodesperelem, phys%neq, Mesh%Ndim))
    u_3d_local = 0.
    u_tilde_3d_local = 0.
    q_4d_local = 0.

    DO i = 1, Mesh%Nelems
      index = Mesh%loc2glob_el(i)
      u_3d_local(i,:,:) = u_3d(index,:,:)
      q_4d_local(i,:,:,:) = q_4D(index,:,:,:)
    ENDDO

    DO i = 1, Mesh%Nfaces
      u_tilde_3d_local(i,:,:) = u_tilde_3d(Mesh%loc2glob_fa(i),:,:)
    ENDDO

    DEALLOCATE(sol%u, sol%u_tilde, sol%q)

    ALLOCATE(sol%u(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq))
    ALLOCATE(sol%u_tilde(Mesh%Nfaces*Mesh%Nnodesperface*phys%neq))
    ALLOCATE(sol%q(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq*Mesh%Ndim))
    sol%u = 0
    sol%u_tilde = 0
    sol%q = 0

    CALL flatten_row_major(u_3d_local, sol%u, SIZE(u_3d_local,1), SIZE(u_3d_local,2), SIZE(u_3d_local,3))
    CALL flatten_row_major(u_tilde_3d_local, sol%u_tilde, SIZE(u_tilde_3d_local,1), SIZE(u_tilde_3d_local,2), SIZE(u_tilde_3d_local,3))
    CALL flatten_row_major_4D(q_4d_local, sol%q, SIZE(q_4d_local,1), SIZE(q_4d_local,2), SIZE(q_4d_local,3),SIZE(q_4d_local,4))


    DEALLOCATE(u_3d, u_tilde_3d, q_4d)
    DEALLOCATE(u_3d_local, u_tilde_3d_local, q_4d_local)
    NULLIFY(u_glob, u_tilde_glob, q_glob)

  ENDSUBROUTINE solution_decomposition
#endif

  SUBROUTINE adaptivity()

    ! Start timing
    IF (utils%timing) THEN
       CALL cpu_TIME(timing%tps1)
       CALL system_CLOCK(timing%cks1, timing%clock_rate1)
    END IF

    IF(restart_adapt) THEN
       time%t = time%t + time%dt
       time%it = time%it + 1
       time%ik = time%ik + 1
       sol%Nt = sol%Nt + 1
    ENDIF

    !! Start refining procedure
    count_adapt = count_adapt + 1

    ! Deep copy previous Mesh and reference element
    CALL free_mesh_loc(Mesh_prec)
    CALL deep_copy_mesh_struct(Mesh,Mesh_prec)
    CALL free_reference_element_pol(refElPol_prec)
    CALL deep_copy_refel_struct(refElPol,refElPol_prec)

    ! Rescale to dimensional values
    Mesh%X = Mesh%X*phys%lscale
    Mesh_prec%X = Mesh_prec%X*phys%lscale

    ! Free elmat, mat, magnetic field, Jtor, puff (to save some memory for the adaptivity)
    CALL free_before_adaptivity()
    CALL adaptivity_new(mesh_name,count_adapt,restart_adapt)
    ! Call estimator, estimator_indicator or indicator
    IF (adapt%evaluator .EQ. 2) THEN
       CALL adaptivity_estimator(mesh_name, adapt%param_est, count_adapt, order)
    ELSEIF ((adapt%evaluator .EQ. 1) ) THEN
       CALL adaptivity_indicator(mesh_name, adapt%thr_ind, adapt%param_est, count_adapt, order)
    ELSEIF((adapt%evaluator .EQ. 0) .OR. (restart_adapt)) THEN
       CALL adaptivity_indicator_estimator(mesh_name, adapt%thr_ind, adapt%param_est, count_adapt, order)
    ELSE
       WRITE(*,*) "Choice of adaptivity evaluator not valid. STOP."
    ENDIF


#ifdef PARALL

    Mesh%X = Mesh%X/phys%lscale

    ! Save new mesh and solution
    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE (count_adapt_char, *) count_adapt
       CALL HDF5_save_mesh("./res/new_mesh_n" // TRIM(ADJUSTL(count_adapt_char)) // ".h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
       !CALL HDF5_save_solution("./res/projected_solution_n" // TRIM(ADJUSTL(count_adapt_char)))
    ENDIF

    ! Domain decomposition on the new mesh
    CALL split_mesh(MPIvar%glob_size, 3, .FALSE.)
    CALL mesh_preprocess(ierr)
    Mesh%X = Mesh%X*phys%lscale

    ! Restart communication
    CALL init_com()

    ! Ghost cells are removed from the mesh and the solution and then the result is gathered over the processes
    CALL gather_mesh(Mesh_prec,  T_glob = T_glob, X_glob = X_glob)
    CALL gather_solution(Mesh_in = Mesh_prec, Nnodesperelem = Mesh_prec%Nnodesperelem, u_in = sol%u_conv, q_in = sol%q_conv, u_glob = u_glob, q_glob = q_glob)
    ! Project the check-point solution to the new mesh
    CALL projectSolutionDifferentMeshes_general(T_glob,X_glob,Mesh%T,Mesh%X, u_glob, q_glob, sol%u, sol%q)

    DEALLOCATE(X_glob,T_glob, u_glob, q_glob)
    NULLIFY(X_glob, T_glob, u_glob, q_glob)

    ! Ghost cells are removed from the mesh and the solution and then the result is gathered over the processes
    CALL gather_mesh(Mesh_in = Mesh, T_glob = T_glob, X_glob = X_glob)
    CALL gather_solution(Mesh_in = Mesh, Nnodesperelem = Mesh%Nnodesperelem, u_in = sol%u, q_in = sol%q, u_glob = u_glob, q_glob = q_glob)

    DEALLOCATE(X_glob,T_glob, u_glob, q_glob)
    NULLIFY(X_glob, T_glob, u_glob, q_glob)

#else
    Mesh%X = Mesh%X*phys%lscale
    Mesh%X = Mesh%X/phys%lscale

    ! Project the check-point solution to the new mesh
    CALL projectSolutionDifferentMeshes_general(Mesh_prec%T,Mesh_prec%X,Mesh%T, Mesh%X, sol%u_conv, sol%q_conv, sol%u, sol%q)

#endif
    IF(MPIvar%glob_id .EQ. 0) THEN
#ifndef PARALL
       WRITE(*,*) "Number of elements previous mesh: ", Mesh_prec%Nelems
       WRITE(*,*) "Number of elements current mesh:  ", Mesh%Nelems
#else
       WRITE(*,*) "Number of elements previous mesh: ", Mesh_prec%Nel_glob
       WRITE(*,*) "Number of elements current mesh:  ", Mesh%Nel_glob
#endif
    ENDIF

    ! Extract new trace solution
    CALL update_solution_arrays()

    ! Rescale back to adimensional values
    Mesh%X = Mesh%X/phys%lscale
    Mesh_prec%X = Mesh_prec%X/phys%lscale

    ! reset parameters
    CALL set_parameters()

    ! Re-initialize magnetic field (the Mesh is needed)
    CALL initialize_magnetic_field()

    ! Re-load magnetic field and Jtor
    CALL load_magnetic_field_Jtor()

    ! Re-initialise puff, only if neutrals are present
    CALL initialize_puff()

    ! Re-Allocation and initialization of the elemental matrices
    CALL init_elmat()

    ! restart to first ever iteration
    matK%start = .TRUE.

    ! Re-Initialize shock capturing
    IF ((switch%shockcp .GT. 0) .OR. ((adapt%adaptivity) .AND. (adapt%shockcp_adapt .GT. 0)))  THEN
       CALL initializeShockCapturing()
    ENDIF

    ! Re-Compute first equation (definition of the gradient)
    CALL HDG_precalculatedfirstequation()

    !! UPDATE VARIABLES
    errNR_adapt = 1e10
    IF(restart_adapt) THEN
       time%t = time%t - time%dt
       time%it = time%it - 1
       time%ik = time%ik - 1
       sol%Nt = sol%Nt - 1
       restart_adapt = .FALSE.
    ENDIF

    IF (utils%timing) THEN
       CALL cpu_TIME(timing%tpe1)
       CALL system_CLOCK(timing%cke1, timing%clock_rate1)
       timing%runtadapt = timing%runtadapt + (timing%cke1-timing%cks1)/REAL(timing%clock_rate1)
       timing%cputadapt = timing%cputadapt + timing%tpe1-timing%tps1
    END IF

  ENDSUBROUTINE adaptivity


  SUBROUTINE free_before_adaptivity()
    IF(ASSOCIATED(phys%B)) DEALLOCATE(phys%B)
    IF(ASSOCIATED(phys%magnetic_flux)) DEALLOCATE(phys%magnetic_flux)
    IF(ASSOCIATED(phys%Bperturb)) DEALLOCATE(phys%Bperturb)
    IF(ASSOCIATED(phys%magnetic_psi)) DEALLOCATE(phys%magnetic_psi)
    IF(ASSOCIATED(phys%Jtor))     DEALLOCATE(phys%Jtor)
    IF(ASSOCIATED(phys%puff_exp)) DEALLOCATE(phys%puff_exp)

    NULLIFY(phys%magnetic_flux, phys%Bperturb, phys%magnetic_psi, phys%Jtor, phys%puff_exp)
    CALL free_el_mat()
    CALL free_mat()
  ENDSUBROUTINE free_before_adaptivity


  !************************************************
  ! Display results
  !************************************************
  SUBROUTINE displayResults()
    INTEGER              :: ieq
    REAL*8, ALLOCATABLE   :: uphy(:, :)
    REAL*8               :: Vmax(phys%npv), Vmin(phys%npv)
    nu = SIZE(sol%u)

    ALLOCATE (uphy(nu/phys%Neq, phys%npv))

    ! Compute physical variables
    CALL cons2phys(TRANSPOSE(RESHAPE(sol%u, (/phys%Neq, nu/phys%Neq/))), uphy)
    DO ieq = 1, phys%npv
       Vmax(ieq) = MAXVAL(uphy(:, ieq))
       Vmin(ieq) = MINVAL(uphy(:, ieq))
    END DO

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmax, phys%npv, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmin, phys%npv, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, '(" * ", 60("-"), "*")')
       WRITE (6, '(" * Time (adimensional) = ", E12.5, 27X, " *")') time%t
       WRITE (6, '(" * Dt (adimensional)        = ", E12.5, 27X, " *")') time%dt
       WRITE (6, '(" * ", 45("^"), 14X, " *")')
       WRITE (6, '(" * ", 10("_"), "      Minimum   ", 4X, "     Maximum   ", 14X, " *")')
       DO ieq = 1, phys%npv
          WRITE (6, '(" * ", A7, " -->", ES16.8, 3X, ES16.8, 13X, " *")') &
               & TRIM(phys%phyVarNam(ieq)), Vmin(ieq), Vmax(ieq)
       END DO
       WRITE (6, '(" * ", 60("-"), "*")')
       WRITE (6, '(" * Time residual  = ", 1X, 2(E16.8, 2X), 13X, " *")') sol%tres(it)
       WRITE (6, *) '  '
       WRITE (6, *) '  '
       WRITE (6, *) '  '
    END IF
#ifdef PARALL
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
    DEALLOCATE (uphy)

  END SUBROUTINE displayResults

  !************************************************
  ! Compute the residual
  !************************************************
  FUNCTION computeResidual(u, uref, coeff) RESULT(res)
    REAL*8   :: u(:), uref(:)
    INTEGER  :: nglo
    REAL*8   :: res, sum2, coeff
#ifdef PARALL
    INTEGER  :: ierr
#endif

    sum2 = SUM((u - uref)**2)
    nglo = SIZE(u)

#ifdef PARALL
    CALL mpi_allreduce(MPI_IN_PLACE, sum2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL mpi_allreduce(nu, nglo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    res = SQRT(sum2)/SQRT(REAL(nglo))/coeff
  END FUNCTION computeResidual

  !************************************************
  ! Set the name of the solution
  !************************************************
  SUBROUTINE setSolName(save_name, mesh_name, it, convNR, convT)
    CHARACTER(LEN=1024), INTENT(INOUT):: save_name
    CHARACTER(LEN=1024), INTENT(IN)   :: mesh_name
    INTEGER, INTENT(IN)                  :: it
    LOGICAL, INTENT(IN)                  :: convNR, convT
    CHARACTER(LEN=20)                 :: Num
    INTEGER                             :: l, i

    ! At the beginning, the save name is the mesh name..
    save_name = TRIM(ADJUSTL(mesh_name))

    ! Eliminate path info
    l = LEN(save_name)
    i = INDEX(save_name, '/', .TRUE.)
    save_name = save_name(i + 1:l)
    ! look for the P of the order of the polynomial (P4 for example)
    i = INDEX(save_name, 'P', .TRUE.)
    WRITE (Num, "(i10)") refElPol%nDeg
    save_name(i+1:i+1) = TRIM(ADJUSTL(Num))


#ifdef TOR3D
    ! Add the number of toroidal elements
    WRITE (Num, "(i10)") numer%ntor
    Num = TRIM(ADJUSTL(Num))
    save_name = TRIM(ADJUSTL(save_name))//"_Ntor"//Num

    ! Add the poloidal interpolation in the toroidal direction
    WRITE (Num, "(i10)") numer%ptor
    Num = TRIM(ADJUSTL(Num))
    save_name = TRIM(ADJUSTL(save_name))//"Ptor"//Num

#endif
#ifndef KEQUATION
    ! Diffusion
    WRITE (Num, "(E10.3)") phys%diff_n*simpar%refval_diffusion
#else
    WRITE (Num, "(E10.3)") (phys%diff_n+phys%diff_k_min)*simpar%refval_diffusion
#endif
    save_name = TRIM(ADJUSTL(save_name))//"_DPe"//TRIM(ADJUSTL(Num))
#ifdef TEMPERATURE
    WRITE (Num, "(E10.3)") phys%diff_pari
    save_name = TRIM(ADJUSTL(save_name))//"_DPai"//TRIM(ADJUSTL(Num))

    WRITE (Num, "(E10.3)") phys%diff_pare
    save_name = TRIM(ADJUSTL(save_name))//"_DPae"//TRIM(ADJUSTL(Num))
#endif
    ! Complete the save name: if not converged the NR, I put NR + the iteration number
    IF (.NOT. convNR) THEN
       WRITE (Num, "(i10)") it
       Num = TRIM(ADJUSTL(Num))
       k = INDEX(Num, " ") - 1
       save_name = TRIM(ADJUSTL(save_name))//"_NR"//REPEAT("0", 4 - k)//TRIM(ADJUSTL(Num))
    END IF

    ! Complete the save name: if not converged the time scheme, I put the iteration number
    IF (.NOT. convT) THEN
       WRITE (Num, "(i10)") it/utils%freqsave
       Num = TRIM(ADJUSTL(Num))
       k = INDEX(Num, " ") - 1
       save_name = TRIM(ADJUSTL(save_name))//"_"//REPEAT("0", 4 - k)//TRIM(ADJUSTL(Num))
    END IF

    IF (switch%decoup) THEN
       save_name = TRIM(ADJUSTL(save_name))//'_UNCP'
    ENDIF

    ! Add "Sol_"
#ifdef TOR3D
    save_name = 'Sol3D_'//TRIM(ADJUSTL(save_name))
#else
    save_name = 'Sol2D_'//TRIM(ADJUSTL(save_name))
#endif
    ! Add save_folder
    save_name = TRIM(ADJUSTL(input%save_folder))//save_name
  END SUBROUTINE setSolName

  SUBROUTINE compute_dt(errlstime)
    REAL*8, INTENT(in) :: errlstime

    IF ((errlstime*time%dt) < 1e-3) THEN
       WRITE (6, *) "******** Changing time step ***********"
       time%dt = time%dt*2.
    ENDIF
  END SUBROUTINE compute_dt


#ifdef TOR3D
  !**********************************************
  ! Definition of the toroidal discretization
  !**********************************************
  SUBROUTINE define_toroidal_discretization
    INTEGER :: i, ntorloc, itor, itorg, nnodes_toroidal
    INTEGER :: ind(numer%ptor + 1)
    REAL*8  :: tdiv(numer%ntor + 1), tel(numer%ptor + 1), htor

    ! Toroidal discretization
    htor = numer%tmax/numer%ntor
    tdiv = 0.
    DO i = 1, numer%ntor
       tdiv(i + 1) = i*htor
    END DO
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif

    nnodes_toroidal = ntorloc + 1 + (numer%ptor - 1)*ntorloc
    Mesh%Nnodes_toroidal = nnodes_toroidal
    ALLOCATE (Mesh%toroidal(nnodes_toroidal))

    DO itor = 1, ntorloc
#ifdef PARALL
       itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
       IF (itorg == numer%ntor + 1) itorg = 1
#else
       itorg = itor
#endif
       tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
       ind = (itor - 1)*numer%ptor + (/(i, i=1, numer%ptor + 1)/)
       Mesh%toroidal(ind) = tel
    END DO
  END SUBROUTINE define_toroidal_discretization
#endif

  SUBROUTINE update_solution_arrays()
    DEALLOCATE(sol%u_tilde)
    DEALLOCATE(sol%u_tilde0)
    ALLOCATE(sol%u_tilde(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
    ALLOCATE(sol%u_tilde0(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
    CALL extractFaceSolution()
    sol%u_tilde0 = sol%u_tilde

    IF(restart_adapt) THEN
       DEALLOCATE(sol%u0)
       ALLOCATE (sol%u0(SIZE(sol%u), time%tis))
       sol%u0 = 0.
       sol%u0(:, 1) = sol%u
    ENDIF

    IF(.NOT. ((adapt%NR_adapt) .AND. (MOD(ir,adapt%freq_NR_adapt) .EQ. 0))) THEN
       ! Update check-point solution to the one projected on the new mesh
       DEALLOCATE(sol%u_conv)
       ALLOCATE(sol%u_conv(SIZE(sol%u)))
       sol%u_conv = sol%u
       DEALLOCATE(sol%q_conv)
       ALLOCATE(sol%q_conv(SIZE(sol%q)))
       sol%q_conv = sol%q
    ENDIF

  ENDSUBROUTINE update_solution_arrays

  SUBROUTINE update_solution()
    IF (time%tis .GT. 1 .AND. it .LT. time%tis) THEN
       DO is = it, 1, -1
          sol%u0(:, is + 1) = sol%u0(:, is)
       END DO
    ELSEIF (time%tis .GT. 1) THEN
       DO is = time%tis, 2, -1
          sol%u0(:, is) = sol%u0(:, is - 1)
       END DO
    END IF
    sol%u0(:, 1) = sol%u
  ENDSUBROUTINE update_solution

  SUBROUTINE update_uiter()
    IF(SIZE(uiter) .NE. SIZE(sol%u0(:,1))) THEN
       DEALLOCATE(uiter)
       ALLOCATE(uiter(SIZE(sol%u0(:,1))))
    ENDIF
    uiter = sol%u0(:, 1)
  ENDSUBROUTINE update_uiter

  SUBROUTINE update_uiter_qiter_best(uiter_best, qiter_best, u, q)
    REAL*8, POINTER, INTENT(OUT)  :: uiter_best(:), qiter_best(:)
    REAL*8, INTENT(IN)           :: u(:), q(:)

    IF(ASSOCIATED(uiter_best)) THEN
       IF(SIZE(uiter_best) .NE. SIZE(u)) THEN
          DEALLOCATE(uiter_best)
          DEALLOCATE(qiter_best)
          ALLOCATE(uiter_best(SIZE(u)))
          ALLOCATE(qiter_best(SIZE(q)))
          uiter_best = u
          qiter_best = q
       ENDIF
    ELSE
       ALLOCATE(uiter_best(SIZE(u)))
       ALLOCATE(qiter_best(SIZE(q)))
       uiter_best = u
       qiter_best = q
    ENDIF


  ENDSUBROUTINE update_uiter_qiter_best


  SUBROUTINE set_parameters()
    order = refElPol%nDeg
#ifdef TOR3D
    Ndim = 3                                               ! Number of dimensions
#ifdef PARALL
    ntorloc = numer%ntor/MPIvar%ntor
#else
    ntorloc = numer%ntor
#endif
    Nel = Mesh%Nelems*ntorloc                             ! Number of 3D elements
    Np = refElPol%Nnodes2D*refElTor%Nnodes1D             ! Number of nodes for each 3D element
    Nfl = refElPol%Nnodes1D*refElTor%Nnodes1D             ! Number of nodes in the lateral faces
    Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
    Nf = Mesh%Nfaces                                     ! Number of faces in the 2D mesh
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) + phys%Neq*refElPol%Nnodes2D*Mesh%Nelems ! Size of utilde
    ELSE
       nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
    ENDIF
#else
    nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
#endif
    Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
    nu = phys%Neq*Nel*Np
#else
    Ndim = Mesh%ndim
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfp = refElPol%Nfacenodes*Nf
    nut = Mesh%Nfaces*Mesh%Nnodesperface*phys%Neq
    nu = Mesh%Nelems*Mesh%Nnodesperelem*phys%Neq
#endif
  ENDSUBROUTINE set_parameters

  SUBROUTINE initialize_solu0_uiter_uconv()
    ALLOCATE(uiter(SIZE(sol%u)))
    ALLOCATE(uiter_best(SIZE(sol%u)))
    ALLOCATE(qiter_best(SIZE(sol%q)))
    ALLOCATE(sol%u0(SIZE(sol%u), time%tis))
    ALLOCATE (sol%u_conv(SIZE(sol%u)))
    ALLOCATE (sol%q_conv(SIZE(sol%q)))
    sol%u0 = 0.
    sol%u0(:, 1) = sol%u
    uiter = sol%u
    uiter_best = sol%u
    qiter_best = sol%q
    sol%u_conv = sol%u
    sol%q_conv = sol%q
  ENDSUBROUTINE initialize_solu0_uiter_uconv

  SUBROUTINE update_u0(u)
    REAL*8, INTENT(IN)      :: u(:)

    IF(SIZE(sol%u0,1) .NE. SIZE(u)) THEN
       DEALLOCATE(sol%u0)
       ALLOCATE (sol%u0(SIZE(u), time%tis))
    ENDIF
    sol%u0 = 0.
    sol%u0(:, 1) = u
  ENDSUBROUTINE update_u0

  SUBROUTINE update_uconv_qconv(u_conv, q_conv)
    REAL*8, INTENT(IN)        :: u_conv(:), q_conv(:)
    IF(SIZE(sol%u_conv) .NE. SIZE(u_conv)) THEN
       DEALLOCATE(sol%u_conv)
       ALLOCATE(sol%u_conv(SIZE(u_conv)))
    ENDIF
    IF(SIZE(sol%q_conv) .NE. SIZE(q_conv)) THEN
       DEALLOCATE(sol%q_conv)
       ALLOCATE(sol%q_conv(SIZE(q_conv)))
    ENDIF
    sol%u_conv = u_conv
    sol%q_conv = q_conv
  ENDSUBROUTINE update_uconv_qconv

  SUBROUTINE update_dumpnr()
    numer%dumpnr = numer%dumpnr_min+(numer%dumpnr_max-numer%dumpnr_min)/2.*(1+TANH((ir-numer%dumpnr_n0)/numer%dumpnr_width))
  ENDSUBROUTINE update_dumpnr

  SUBROUTINE initialize_solution()
    IF (nb_args .EQ. 3) THEN
       ALLOCATE(sol%u_tilde(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
       ALLOCATE(sol%u_tilde0(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
       CALL extractFaceSolution()
    ELSE IF (nb_args .EQ. 2) THEN
       ! restart simulation: load solution from file (the name is given in argument)
       CALL HDF5_load_solution(save_name)
       ALLOCATE(sol%u_tilde0(SIZE(sol%u_tilde)))
    ELSE
       CALL init_sol()
    END IF
  ENDSUBROUTINE initialize_solution

  SUBROUTINE load_mesh()

    IF (nb_args .EQ. 3) THEN
       ! find projection points
       CALL load_mesh_h5(mesh_name)
       ALLOCATE(xs(SIZE(Mesh%T,1)*SIZE(Mesh%T,2),2))
       xs = Mesh%X(colint(TRANSPOSE(Mesh%T)),:)
       CALL free_mesh()

       ! Load solution to project from
       CALL load_mesh_h5(mesh_name_proj)
       CALL create_reference_element(refElPol, 2, verbose = 1)
       CALL mesh_preprocess(ierr)
       CALL HDF5_load_solution(save_name)

       ! Project the solution
       !call projectSolutionDifferentMeshes(xs)
       DEALLOCATE(xs)
       CALL free_mesh()
       CALL free_reference_element()
    ENDIF

    IF(switch%readMeshFromSol) THEN
       CALL HDF5_load_mesh_from_solution(save_name)
    ELSE

       ! Load the mesh file from gmsh input or h5
       
      IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
         CALL load_gmsh_mesh(mesh_name, 0)
      ELSE
         CALL load_gmsh_mesh(mesh_name, 1)
      ENDIF
      CALL create_reference_element(refElPol, 2, verbose = 1)

      IF((switch%set_2d_order) .AND. (refElPol%nDeg .NE. switch%order_2d) ) THEN
         CALL mesh_preprocess_serial(ierr)
         CALL set_order_mesh(switch%order_2d)
         CALL free_reference_element_pol(refElPol)
         Mesh%X = Mesh%X*phys%lscale
         Mesh%xmax = Mesh%xmax*phys%lscale
         Mesh%xmin = Mesh%xmin*phys%lscale
         Mesh%ymax = Mesh%ymax*phys%lscale
         Mesh%ymin = Mesh%ymin*phys%lscale
      ENDIF

      IF(switch%gmsh2h5) THEN
         IF(MPIvar%glob_id .EQ. 0) THEN
           h5_filename = TRIM(ADJUSTL(mesh_name)) // '.h5'
           CALL convert_gmsh_to_hdf5(h5_filename, SIZE(Mesh%X,2), SIZE(Mesh%T,1), SIZE(Mesh%Tb,1), SIZE(Mesh%X,1), SIZE(Mesh%T,2), refElPol%nDeg+1, 0, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
         ENDIF
      ENDIF


      IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
         Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
      END IF
      IF (adapt%shockcp_adapt .GT. 0) THEN
         gmsh_filename      = TRIM(ADJUSTL(mesh_name))//'.msh'
         i = 1
         DO
            IF((i+2) .GE. LEN(gmsh_filename)) THEN
               WRITE(*,*) "GMSH file input not found, check input syntax."
               STOP
            ENDIF
            IF((gmsh_filename(i:i) .EQ. 'm') .AND. (gmsh_filename(i+1:i+1) .EQ. 's') .AND. (gmsh_filename(i+2:i+2) .EQ. 'h')) THEN
               gmsh_filename_mesh = TRIM(ADJUSTL(gmsh_filename(1:i-4))) // 'P1.mesh'
               EXIT
            ENDIF
            i = i + 1
         ENDDO
         CALL copy_file(gmsh_filename,"./res/temp.msh")
         CALL copy_file(gmsh_filename_mesh,"./res/temp.mesh")
      ENDIF
       
    ENDIF
  ENDSUBROUTINE load_mesh

  SUBROUTINE check_input_arguments()

    ! Check the number of input arguments
    nb_args = iargc()

    IF (nb_args .LT. 1) THEN
       PRINT *, " Error: the mesh name is needed"
       STOP
    END IF

    CALL getarg(1, mesh_name)
    mesh_name = ADJUSTL(TRIM(mesh_name))

    IF (nb_args .GT. 1) THEN
       CALL getarg(2, save_name)
       save_name = ADJUSTL(TRIM(save_name))
       IF(MPIvar%glob_id .EQ. 0) THEN
         PRINT *, " Restart simulation with solution: ", save_name
       ENDIF
    END IF

    IF (nb_args .GT. 2) THEN
       CALL getarg(3, mesh_name_proj)
       mesh_name_proj = ADJUSTL(TRIM(mesh_name_proj))
       PRINT *, " Projecting solution from: ", mesh_name_proj
    END IF

    IF (nb_args .GT. 3) THEN
       PRINT *, " Too many arguments "
       STOP
    END IF
  ENDSUBROUTINE check_input_arguments

  SUBROUTINE check_for_NaNs()
    DO i = 1, SIZE(sol%u)
       IF(ISNAN(sol%u(i))) THEN
          WRITE (6, *) "NaN detected. STOPPING."
          STOP
       END IF
    END DO
  ENDSUBROUTINE check_for_NaNs

  SUBROUTINE project_u0_newmesh()

    IF(.NOT. ((adapt%NR_adapt) .AND. (MOD(ir,adapt%freq_NR_adapt) .EQ. 0))) THEN
       ALLOCATE(u0_temp(SIZE(sol%u0,1),SIZE(sol%u0,2)))
       u0_temp = sol%u0
    ELSE
       ALLOCATE(u0_temp(SIZE(sol%u_conv),SIZE(sol%u0,2)))
       u0_temp = 0
       u0_temp(:,1) = sol%u_conv
    ENDIF

    DEALLOCATE(sol%u0)
    ALLOCATE(sol%u0(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq,time%tis))
    sol%u0 = 0.

    Mesh%X = Mesh%X*phys%lscale
    Mesh_prec%X = Mesh_prec%X*phys%lscale

#ifndef PARALL
    CALL projectSolutionDifferentMeshes_general_arrays(Mesh_prec%T,Mesh_prec%X, Mesh%T, Mesh%X, u1 = u0_temp(:,1), u2 = sol%u0(:,1))
#else
    ! Ghost cells are removed from the mesh and the solution and then the result is gathered over the processes

    CALL gather_mesh(Mesh_in = Mesh_prec, T_glob = T_glob, X_glob = X_glob)
    CALL gather_solution(Mesh_in = Mesh_prec, Nnodesperelem = Mesh_prec%Nnodesperelem, u_in = u0_temp(:,1), u_glob = u_glob)

    CALL projectSolutionDifferentMeshes_general_arrays(T_glob,X_glob, Mesh%T, Mesh%X, u1 = u_glob, u2 = sol%u0(:,1))
    DEALLOCATE(T_glob, X_glob, u_glob)
    NULLIFY(T_glob, X_glob, u_glob)
#endif



    Mesh%X = Mesh%X/phys%lscale
    Mesh_prec%X = Mesh_prec%X/phys%lscale

    DEALLOCATE(u0_temp)
  ENDSUBROUTINE project_u0_newmesh

  SUBROUTINE reduce_diffusion()
    ! Update the diffusion, the elemental matrices and the solution
    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, *) "************************************************"
       WRITE (6, *) "Reducing diffusion: ", phys%diff_n*switch%diffred*simpar%refval_diffusion

#ifdef KEQUATION
       WRITE (6, *) "K diffusion min: ", phys%diff_k_min*switch%diffred*simpar%refval_diffusion
#endif
       WRITE (6, *) "************************************************"
    END IF
    phys%diff_n = phys%diff_n*switch%diffred
    phys%diff_u = phys%diff_u*switch%diffred
#ifdef TEMPERATURE
    phys%diff_e = phys%diff_e*switch%diffred
    phys%diff_ee = phys%diff_ee*switch%diffred
#endif
#ifdef VORTICITY
    phys%diff_vort = phys%diff_vort*switch%diffred
    phys%diff_pot = phys%diff_pot*switch%diffred
#endif
#ifdef KEQUATION
    phys%diff_k_min = phys%diff_k_min*switch%diffred
    phys%diff_k_max = phys%diff_k_max!*switch%diffred
#endif
  ENDSUBROUTINE reduce_diffusion

  SUBROUTINE print_timing_infos()
    ! Code timing
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%timing) THEN
          cputtot = 1e-8
          runttot = 1e-8
          cputtot = timing%cputpre + timing%cputjac + timing%cputbcd + timing%cputmap + timing%cputass + timing%cputglb + timing%cputsol + timing%cputadapt
          runttot = timing%runtpre + timing%runtjac + timing%runtbcd + timing%runtmap + timing%runtass + timing%runtglb + timing%runtsol + timing%runtadapt
#ifdef PARALL
          cputtot = cputtot + timing%cputcom
          runttot = runttot + timing%runtcom
#endif
          WRITE(6, *) " "
          WRITE(6, *) " "
          WRITE(6, *) " "
          WRITE(6, '(" *", 90("*"), "**")')
          WRITE(6, '(" *", 36X, "CODE TIMING ( Nthreads = ",i2,")", 26X, " *")') Nthreads
          WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")')
          WRITE(6, '(" *", 2X,  "Precal. matr     : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
               &timing%cputpre,timing%cputpre/cputtot*100,timing%runtpre,timing%runtpre/runttot*100,timing%cputpre/timing%runtpre/Nthreads
          WRITE(6, '(" *", 2X,  "Jacobian         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputjac,timing%cputjac/cputtot*100,timing%runtjac,timing%runtjac/runttot*100,timing%cputjac/timing%runtjac/Nthreads
          WRITE(6, '(" *", 2X,  "Mapping          : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputmap,timing%cputmap/cputtot*100,timing%runtmap,timing%runtmap/runttot*100,timing%cputmap/timing%runtmap/Nthreads
          WRITE(6, '(" *", 2X,  "Boundary cond.   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputbcd,timing%cputbcd/cputtot*100,timing%runtbcd,timing%runtbcd/runttot*100,timing%cputbcd/timing%runtbcd/Nthreads
          WRITE(6, '(" *", 2X,  "Assembly         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputass,timing%cputass/cputtot*100,timing%runtass,timing%runtass/runttot*100,timing%cputass/timing%runtass/Nthreads
          WRITE(6, '(" *", 2X,  "Solve glob. syst.: ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputglb,timing%cputglb/cputtot*100,timing%runtglb,timing%runtglb/runttot*100,timing%cputglb/timing%runtglb/Nthreads
          WRITE(6, '(" *", 2X,  "Element solution : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputsol,timing%cputsol/cputtot*100,timing%runtsol,timing%runtsol/runttot*100,timing%cputsol/timing%runtsol/Nthreads
          IF(adapt%adaptivity .EQV. .TRUE.) THEN
             WRITE(6, '(" *", 2X,  "Adaptivity       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
                  &timing%cputadapt,timing%cputadapt/cputtot*100,timing%runtadapt,timing%runtadapt/runttot*100,timing%cputadapt/timing%runtadapt/Nthreads
          ENDIF
#ifdef PARALL
          WRITE(6, '(" *", 2X,  "Communications   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputcom,timing%cputcom/cputtot*100,timing%runtcom,timing%runtcom/runttot*100,timing%cputcom/timing%runtcom/Nthreads
#endif
          WRITE(6, '(" *", 2X, "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")')   &
               cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
          WRITE(6, '(" *", 90("*"), "**")')
          WRITE(6, *) " "
          WRITE(6, *) " "
          WRITE(6, *) " "
       END IF

       IF (lssolver%timing) THEN
          cputtot = timing%clstime1 + timing%clstime2 + timing%clstime3 + timing%clstime4 + timing%clstime5 + timing%clstime6
          runttot = timing%rlstime1 + timing%rlstime2 + timing%rlstime3 + timing%rlstime4 + timing%rlstime5 + timing%rlstime6
          WRITE(6, '(" *", 90("*"), "**")')
          IF (lssolver%sollib .EQ. 1) THEN
             WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PASTIX ( Nthreads = ",i2,")", 14X, " *")') Nthreads
          ELSE IF (lssolver%sollib .EQ. 2) THEN
             WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PSBLAS ( Nthreads = ",i2,")", 14X, " *")') Nthreads
          ELSE IF (lssolver%sollib .EQ. 3) THEN
             WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PETSc ( Nthreads = ",i2,")", 14X, " *")') Nthreads
          ENDIF
          WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")')
          IF (lssolver%sollib .EQ. 1) THEN
             WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
                  &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
             WRITE(6, '(" *", 2X,  "Check mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
             WRITE(6, '(" *", 2X,  "Anal. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
             WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
             WRITE(6, '(" *", 2X,  "LU decomp.       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime5,timing%clstime5/cputtot*100,timing%rlstime5,timing%rlstime5/runttot*100,timing%clstime5/timing%rlstime5/Nthreads
             WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime6,timing%clstime6/cputtot*100,timing%rlstime6,timing%rlstime6/runttot*100,timing%clstime6/timing%rlstime6/Nthreads
             WRITE(6, '(" *", 2X,  "Total time      : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",7X,F5.1 , 10X, " *")') &
                  &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
          ELSEIF (lssolver%sollib .EQ. 2) THEN
             WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
                  &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
             WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
             WRITE(6, '(" *", 2X,  "Build prec       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
             WRITE(6, '(" *", 2X,  "Fill vec         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
             WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime5,timing%clstime5/cputtot*100,timing%rlstime5,timing%rlstime5/runttot*100,timing%clstime5/timing%rlstime5/Nthreads
             WRITE(6, '(" *", 2X,  "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")') &
                  &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
          ELSEIF (lssolver%sollib .EQ. 3) THEN
             WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
                  &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
             WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
             WRITE(6, '(" *", 2X,  "Fill vec         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
             WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
             WRITE(6, '(" *", 2X,  "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")') &
                  &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
          ENDIF
          WRITE (6, '(" *", 90("*"), "**")')
          WRITE (6, *) " "
          WRITE (6, *) " "
          WRITE (6, *) " "
       END IF
    END IF
  ENDSUBROUTINE print_timing_infos

  SUBROUTINE free_main()
#ifdef WITH_PASTIX
    USE solve_pastix, only: terminate_mat_PASTIX
#endif
#ifdef WITH_PETSC
    USE solve_petsc, only: terminate_PETSC, FinalizePETSC
#endif
#ifdef WITH_PSBLAS
    USE solve_psblas, only: terminate_PSBLAS
#endif

    IF(ALLOCATED(uiter)) DEALLOCATE (uiter)
    IF(ASSOCIATED(uiter_best)) DEALLOCATE(uiter_best)
    IF(ASSOCIATED(qiter_best)) DEALLOCATE(qiter_best)
    NULLIFY(uiter_best)
    NULLIFY(qiter_best)
    IF(ALLOCATED(mkelms)) DEALLOCATE (mkelms)
    IF(ALLOCATED(sol%u0)) DEALLOCATE (sol%u0)
    IF(ALLOCATED(oscillations)) DEALLOCATE(oscillations)

    CALL free_splines(splines)

    IF (lssolver%sollib .EQ. 1) THEN
#ifdef WITH_PASTIX
       CALL terminate_mat_PASTIX()
       ! MPI finalization
       CALL MPI_finalize(IERR)
#endif
    ELSEIF (lssolver%sollib .EQ. 2) THEN
#ifdef WITH_PSBLAS
       CALL terminate_PSBLAS()
#endif
    ELSEIF (lssolver%sollib .EQ. 3) THEN
#ifdef WITH_PETSC
       CALL terminate_PETSC()
       CALL FinalizePETSC()
       CALL MPI_finalize(IERR)
#endif
    ENDIF

    CALL free_all
    CALL free_mesh
#ifdef PARALL
    CALL free_mesh_loc(Mesh_glob)
#endif
    CALL free_mesh_loc(Mesh_prec)
    CALL free_reference_element
    CALL free_reference_element_pol(refElPol_prec)

  ENDSUBROUTINE free_main


ENDMODULE Main_utils
