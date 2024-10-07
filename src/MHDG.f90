PROGRAM MHDG
  USE in_out
  USE GMSH_io_module, ONLY: load_gmsh_mesh, HDF5_save_mesh, generate_splines_from_geo_file, convert_gmsh_to_hdf5, gmsh_mesh2d_write
  USE reference_element
  USE preprocess
  USE MPI_OMP
  USE printutils
  USE debug
  USE initialization
#ifdef WITH_PASTIX
  USE solve_pastix
#endif
#ifdef WITH_PSBLAS
  USE solve_psblas
#endif
#ifdef WITH_PETSC
#include "petsc/finclude/petsc.h"
  USE solve_petsc
#endif
  USE adaptivity_common_module
  USE adaptivity_estimator_module
  USE adaptivity_indicator_module
  USE adaptivity_estimator_indicator_module
  USE Postprocess
#ifdef PARALL
  USE Communications
#endif
  USE HDG_LimitingTechniques

  IMPLICIT NONE

  INTEGER                      :: Np, Nel, Nfp, Nf, Ndim, Nthreads
  INTEGER                      :: it, ir, ir_check, it0, nts, nu, nut, nb_args, IERR, k, i,is, order, count_adapt = 0, ir_adapt = 0, checkpoint = 1, divergence_counter_adapt = 0, convergence_counter = 1

  LOGICAL, ALLOCATABLE         :: mkelms(:)
  REAL*8                       :: dt, dt0, errNR, errNR_adapt, errlstime
  REAL*8, ALLOCATABLE          :: uiter(:), uiter_best(:), qiter_best(:), L2err(:),u0_temp(:,:)
  CHARACTER(LEN=1024)          :: mesh_name,mesh_name_proj, save_name
  CHARACTER ( len = 255 )      :: gmsh_filename
  CHARACTER ( len = 255 )      :: gmsh_filename_mesh, h5_filename
  CHARACTER ( len = 50 )       :: count_adapt_char
  REAL*8                       :: cputtot, runttot
  INTEGER                      :: OMP_GET_MAX_THREADS
  INTEGER                      :: clock_rate,clock_start, clock_end
  REAL*8                       :: time_start, time_finish
  REAL*8,ALLOCATABLE           :: xs(:,:)
  INTEGER                      :: switch_save
  REAL*8, ALLOCATABLE          :: error_oscillation(:),oscillations(:)
  REAL*8                       :: max_osc, min_osc
  INTEGER,ALLOCATABLE          :: vector_nodes_unique(:,:)
  INTEGER                      :: N_n_vertex, n_osc
  INTEGER*8                    :: file_id
  TYPE(Mesh_type)              :: Mesh_prec, Mesh_init
#ifdef PARALL
  TYPE(Mesh_type)              :: Mesh_glob
  INTEGER, ALLOCATABLE         :: T_glob(:,:)
  REAL*8, ALLOCATABLE          :: X_glob(:,:)
  REAL*8, POINTER              :: u_glob(:), q_glob(:)
#endif
#ifdef TOR3D
  INTEGER                      :: Neq, ntorloc,Nfl, Np1dPol, Np1dTor, Ng1dPol, Ng1dTor,
#endif
  TYPE(Reference_element_type) :: refElPol_prec


  WRITE (6, *) "STARTING"
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
     PRINT *, " Restart simulation with solution: ", save_name
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

  ! Number of threads
  Nthreads = OMP_GET_MAX_THREADS()

  ! Start timing the code
  CALL cpu_TIME(time_start)
  CALL system_CLOCK(clock_start, clock_rate)

  ! Initialize MPI
  CALL init_MPI_OMP()

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

  IF (nb_args .EQ. 3) THEN
     ! find projection points
     CALL load_mesh(mesh_name)
     ALLOCATE(xs(SIZE(Mesh%T,1)*SIZE(Mesh%T,2),2))
     xs = Mesh%X(colint(TRANSPOSE(Mesh%T)),:)
     CALL free_mesh()

     ! Load solution to project from
     CALL load_mesh(mesh_name_proj)
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
     IF(switch%read_gmsh) THEN
        IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
           CALL load_gmsh_mesh(mesh_name, 0)
        ELSE
           CALL load_gmsh_mesh(mesh_name, 1)
        ENDIF
        CALL create_reference_element(refElPol, 2, verbose = 1)

        IF((switch%set_2d_order) .AND. (refElPol%nDeg .NE. switch%order_2d) ) THEN
           !CALL mesh_preprocess(ierr)
           CALL mesh_preprocess_serial(ierr)
           CALL HDF5_create('./fortran_save_0.h5', file_id, ierr)
           CALL HDF5_array2D_saving_int(file_id, Mesh%T, SIZE(Mesh%T,1), SIZE(Mesh%T,2), 'T_f')
           CALL HDF5_array2D_saving(file_id, Mesh%X, SIZE(Mesh%X,1), SIZE(Mesh%X,2), 'X_f')
           CALL HDF5_close(file_id)
           CALL set_order_mesh(switch%order_2d)
           CALL create_reference_element(refElPol, 2, verbose = 1)
           Mesh%X = Mesh%X*phys%lscale
           Mesh%xmax = Mesh%xmax*phys%lscale
           Mesh%xmin = Mesh%xmin*phys%lscale
           Mesh%ymax = Mesh%ymax*phys%lscale
           Mesh%ymin = Mesh%ymin*phys%lscale
        ENDIF

        IF(switch%gmsh2h5) THEN
           h5_filename = TRIM(ADJUSTL(mesh_name)) // '.h5'
           CALL convert_gmsh_to_hdf5(h5_filename, SIZE(Mesh%X,2), SIZE(Mesh%T,1), SIZE(Mesh%Tb,1), SIZE(Mesh%X,1), SIZE(Mesh%T,2), refElPol%nDeg+1, 0, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
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
     ELSE

#ifndef PARALL
        CALL load_mesh(mesh_name)
#else

        CALL load_mesh_serial(mesh_name)
        ! CALL load_mesh(mesh_name)
#endif

     ENDIF
  ENDIF

  IF(MPIvar%glob_id .EQ. 0) THEN
     WRITE(*,*) "*************************************************"
     WRITE(*,*) "               SPLINE READING                    "
     WRITE(*,*) "*************************************************"
  ENDIF

  IF (switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) THEN
     WRITE(*,*) "Splines read from file ./res/geometries/Circ_InfLim_YesHole_Structured.geo"
     CALL generate_splines_from_geo_file('./res/geometries/Circ_InfLim_YesHole_Structured.geo')
  ELSE
     IF(ANY(Mesh%boundaryFlag .EQ. 5)) THEN
        !WRITE(*,*) "Splines read from file ./res/geometries/West_Mesh_NoHole_farWall.geo"
        !CALL generate_splines_from_geo_file('./res/geometries/West_Mesh_NoHole_farWall.geo')
        WRITE(*,*) "Splines read from file ./res/geometries/TCV_smooth.geo"
        CALL generate_splines_from_geo_file('./res/geometries/TCV_smooth.geo')
     ELSEIF(ANY(Mesh%boundaryFlag .EQ. 8)) THEN
        WRITE(*,*) "Splines read from file ./res/geometries/West_Mesh_YesHole_SmoothCorner.geo"
        CALL generate_splines_from_geo_file('./res/geometries/West_Mesh_YesHole_SmoothCorner.geo')
     ELSE
        WRITE(*,*) "Splines read from file ./res/geometries/West_Mesh_NoHole_SmoothCorner.geo"
        CALL generate_splines_from_geo_file('./res/geometries/West_Mesh_NoHole_SmoothCorner.geo')
     ENDIF
  END IF


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
  order = refElPol%nDeg

  ! create the temp.msh and temp.mesh needed by the adaptivity
  IF((switch%readMeshFromSol) .AND. (adapt%adaptivity)) THEN
     CALL generate_msh_from_solution_mesh('./res/temp.msh')
     CALL convert_msh2mesh('./res/temp')

  ENDIF

#ifdef TOR3D
  !*************************************************
  !                REFERENCE ELEMENT 3D
  !*************************************************
  ! Create the reference element for the toroidal interpolation
  CALL create_reference_element(refElTor, 1, numer%ptor, verbose = 1)
  CALL create_toroidal_structures(refElTor, refElPol)
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
     !CALL mesh_preprocess(ierr)
  ENDIF

  ! #ifdef PARALL
  !
  !   IF(MPIvar%glob_id .eq. 0) THEN
  !     WRITE(*,*) "*************************************************"
  !     WRITE(*,*) "            DOMAIN DECOMPOSITION                 "
  !     WRITE(*,*) "*************************************************"
  !   ENDIF
  !
  !   CALL deep_copy_mesh_struct(Mesh, Mesh_glob)
  !   Mesh_glob%X = Mesh_glob%X*phys%lscale
  !
  !   CALL split_mesh(MPIvar%glob_size,2, .FALSE.)
  !   CALL mpi_barrier(MPI_COMM_WORLD,ierr)
  !   CALL mesh_preprocess(ierr)
  !
  ! #endif



  ! CALL HDF5_save_mesh("./newmesh_pre_init_pre.h5", Mesh%Ndim, mesh%Nelems, mesh%Nextfaces, mesh%Nnodes, mesh%Nnodesperelem, mesh%Nnodesperface, mesh%elemType, mesh%T, mesh%X, mesh%Tb, mesh%boundaryFlag)
  CALL deep_copy_mesh_struct(Mesh, Mesh_init)
  CALL deep_copy_mesh_struct(Mesh, Mesh_init)
#ifdef TOR3D
  ! Define toroidal discretization
  CALL define_toroidal_discretization
#endif

  ! Initialize magnetic field (the Mesh is needed)
  CALL initialize_magnetic_field()

  ! Load magnetic field and, if ohmic src, also Jtor
  IF (switch%ME .EQV. .FALSE.) THEN
     CALL load_magnetic_field()
     IF (switch%ohmicsrc) CALL loadJtorMap()
  END IF


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
     nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) + &
          &phys%Neq*refElPol%Nnodes2D*Mesh%Nelems ! Size of utilde
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

#ifdef PARALL
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
  IF (nb_args .EQ. 3) THEN
     ALLOCATE(sol%u_tilde(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
     ALLOCATE(sol%u_tilde0(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
     CALL extractFaceSolution()
  ELSE IF (nb_args .EQ. 2) THEN
     ! restart simulation: load solution from file (the name is given in argument)
     CALL HDF5_load_solution(save_name)
     ALLOCATE(sol%u_tilde0(SIZE(sol%u_tilde)))
     ! Update magnetic field and, if ohmic src, Jtor also
     ! In case of restart initialize the magnetic configuration to the previous one
     IF (switch%ME .EQV. .TRUE.) THEN
        CALL load_magnetic_field()
        IF (switch%ohmicsrc) CALL loadJtorMap()
     END IF
  ELSE
     CALL init_sol()
  END IF



  ALLOCATE (sol%u_init(SIZE(sol%u)))
  ALLOCATE (sol%q_init(SIZE(sol%q)))
  ALLOCATE (sol%u0_init(SIZE(sol%u)))
  ALLOCATE (sol%q0_init(SIZE(sol%q)))
  ALLOCATE (sol%u_conv(SIZE(sol%u)))
  ALLOCATE (sol%q_conv(SIZE(sol%q)))

  !! without blob!!
  sol%u_init = sol%u
  sol%q_init = sol%q
  sol%u0_init = sol%u
  sol%q0_init = sol%q

  IF (switch%pertini .EQ. 1) THEN
     CALL add_perturbation()
     WRITE(6,*) "Adding perturbation to the initial solution"
  ELSE IF (switch%pertini .EQ. 2) THEN
     CALL add_blob()
     WRITE(6,*) "Adding density blob to initial solution"
  ENDIF

  ! Initialize puff in case of time evolving magnetic configurations
#ifdef NEUTRAL
  IF (switch%ME .EQV. .FALSE.) THEN
     IF (MPIvar%glob_id .EQ. 0) THEN
        WRITE(6,*) 'Puff is analytical'
     ENDIF
  ELSE
     IF (MPIvar%glob_id .EQ. 0) THEN
        WRITE(6,*) 'Puff is experimental'
     ENDIF
     CALL SetPuff()
  END IF
#endif

  ! Save solution
  CALL setSolName(save_name, mesh_name, 0, .TRUE., .FALSE.)
  CALL HDF5_save_solution(save_name)

  CALL mpi_barrier(mpi_comm_world,ierr)


  ! Allocate and initialize uiter and u0
  ALLOCATE(uiter(SIZE(sol%u)))
  ALLOCATE(uiter_best(SIZE(sol%u)))
  ALLOCATE(qiter_best(SIZE(sol%q)))
  errNR_adapt = 1e10
  ir_adapt = 0

  ! Allocate and initialize uiter and u0
  ALLOCATE (sol%u0(SIZE(sol%u), time%tis))

  sol%u0 = 0.
  sol%u0(:, 1) = sol%u

  sol%u_conv = sol%u
  sol%q_conv = sol%q

  switch_save = 1


  IF ((adapt%adaptivity) .AND. (adapt%rest_adapt)) THEN

     time%t = time%t + time%dt
     time%it = time%it + 1
     time%ik = time%ik + 1
     sol%Nt = sol%Nt + 1

     ! Compute Jacobian
     CALL HDG_computeJacobian()
     ! ! Set boundary conditions
     CALL hdg_BC()
     ! ! Compute elemental mapping
     CALL hdg_Mapping()
     ! ! Assembly the global matrix
     CALL hdg_Assembly()
     ! ! Solve linear system
     CALL solve_global_system(ir)

     ! ! Compute element-by-element solution
     CALL compute_element_solution()

     !! START REFINING PROCEDURE
     count_adapt = count_adapt + 1

     Mesh%X = Mesh%X*phys%lscale
     Mesh_init%X = Mesh_init%X*phys%lscale

     !DEALLOCATE ALL BEFORE ADAPTIVITY
     DEALLOCATE(phys%B)
     DEALLOCATE(phys%magnetic_flux)
     IF(ASSOCIATED(phys%Bperturb)) DEALLOCATE(phys%Bperturb)
     IF(ASSOCIATED(phys%Jtor))     DEALLOCATE(phys%Jtor)
     IF(ASSOCIATED(phys%puff_exp)) DEALLOCATE(phys%puff_exp)
     CALL free_el_mat()
     CALL free_mat

     CALL deep_copy_mesh_struct(Mesh,Mesh_prec)
     CALL deep_copy_mesh_struct(Mesh,Mesh_init)
     CALL deep_copy_refel_struct(refElPol,refElPol_prec)
     CALL adaptivity_indicator_estimator(mesh_name, adapt%thr_ind, adapt%param_est, count_adapt, order)
     !CALL adaptivity_estimator(mesh_name, adapt%param_est, count_adapt, order)
     !CALL adaptivity_indicator(mesh_name, adapt%thr_ind, adapt%param_est, count_adapt, order, Mesh_init,Mesh_prec,refElPol_prec)
     !CALL HDF5_save_mesh("./newmesh.h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)

     !CALL HDF5_save_mesh("./old_mesh.h5", Mesh_prec%Ndim, Mesh_prec%Nelems, Mesh_prec%Nextfaces, Mesh_prec%Nnodes, Mesh_prec%Nnodesperelem, Mesh_prec%Nnodesperface, Mesh_prec%elemType, Mesh_prec%T, Mesh_prec%X, Mesh_prec%Tb, Mesh_prec%boundaryFlag)

#ifdef PARALL
     !CALL split_mesh(MPIvar%glob_size, 2, .FALSE.)
     CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
     CALL mesh_preprocess(ierr)

     CALL init_com()

     ! ghost cells are removed from the mesh and the solution and then the result is gathered over the processes
     CALL gather_mesh_solution(Mesh_prec, sol%u_conv, sol%q_conv, T_glob, X_glob, u_glob, q_glob)
     CALL projectSolutionDifferentMeshes_general(Mesh_glob%T,Mesh_glob%X,Mesh%T,Mesh%X, u_glob, q_glob, sol%u, sol%q)
#else
     ! call HDF5_create('./fortran_save_0.h5', file_id, ierr)
     ! call HDF5_array2D_saving_int(file_id, Mesh_prec%T, size(Mesh_prec%T,1), size(Mesh_prec%T,2), 'T_prec_f')
     ! call HDF5_array2D_saving(file_id, Mesh_prec%X, size(Mesh_prec%X,1), size(Mesh_prec%X,2), 'X_prec_f')
     ! call HDF5_array2D_saving_int(file_id, Mesh%T, size(Mesh%T,1), size(Mesh%T,2), 'T_f')
     ! call HDF5_array2D_saving(file_id, Mesh%X, size(Mesh%X,1), size(Mesh%X,2), 'X_f')
     ! call HDF5_array2D_saving_int(file_id, Mesh_prec%Tb, size(Mesh_prec%Tb,1), size(Mesh_prec%Tb,2), 'Tb_prec_f')
     ! call HDF5_array2D_saving_int(file_id, Mesh%Tb, size(Mesh%Tb,1), size(Mesh%Tb,2), 'Tb_f')
     ! call HDF5_close(file_id)

     CALL projectSolutionDifferentMeshes_general(Mesh_prec%T,Mesh_prec%X,Mesh%T, Mesh%X, sol%u_conv, sol%q_conv, sol%u, sol%q)
#endif

     DEALLOCATE(sol%u_tilde)
     DEALLOCATE(sol%u_tilde0)
     ALLOCATE(sol%u_tilde(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
     ALLOCATE(sol%u_tilde0(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))

     CALL extractFaceSolution()

     DEALLOCATE(sol%u0)
     ALLOCATE (sol%u0(SIZE(sol%u), time%tis))
     sol%u0 = 0.
     sol%u0(:, 1) = sol%u

     Mesh%X = Mesh%X/phys%lscale
     Mesh_prec%X = Mesh_prec%X/phys%lscale
     !CALL deep_copy_mesh_struct(Mesh_temp,Mesh)
     WRITE(*,*) "Number of elements previous mesh: ", SIZE(Mesh_prec%T,1)
     WRITE(*,*) "Number of elements current mesh:  ", SIZE(Mesh%T,1)
     CALL deep_copy_mesh_struct(Mesh,Mesh_prec)


     ! Mesh%X = Mesh%X*phys%lscale
     ! CALL projectSolutionDifferentMeshes_general(Mesh_init, Mesh, sol%u0_init, sol%q0_init, sol%u_init, sol%q_init)
     !
     ! call HDF5_create('./fortran_save.h5', file_id, ierr)
     ! call HDF5_array2D_saving_int(file_id, Mesh_init%T, size(Mesh_init%T,1), size(Mesh_init%T,2), 'T_f_init')
     ! call HDF5_array2D_saving(file_id, Mesh_init%X, size(Mesh_init%X,1), size(Mesh_init%X,2), 'X_f_init')
     ! call HDF5_array2D_saving_int(file_id, Mesh%T, size(Mesh%T,1), size(Mesh%T,2), 'T_f')
     ! call HDF5_array2D_saving(file_id, Mesh%X, size(Mesh%X,1), size(Mesh%X,2), 'X_f')
     ! call HDF5_array1D_saving(file_id, sol%u0_init, size(sol%u0_init), 'u0_init_f')
     ! call HDF5_array1D_saving(file_id, sol%q0_init, size(sol%q0_init), 'q0_init_f')
     ! call HDF5_array1D_saving(file_id, sol%u_init, size(sol%u_init), 'u_f')
     ! call HDF5_array1D_saving(file_id, sol%q_init, size(sol%q_init), 'q_f')
     ! call HDF5_close(file_id)
     ! pause
     !
     ! Mesh%X = Mesh%X/phys%lscale



     ! Initialize magnetic field (the Mesh is needed)
     CALL initialize_magnetic_field()

     ! Load magnetic field and, if case is 54 or 59
     CALL load_magnetic_field()
     IF (switch%ohmicsrc) CALL loadJtorMap()


     ! if moving equilibrium case then update the magnetic field, otherwise just continue
     IF(switch%ME) THEN
        CALL loadJtorMap()
        !if (switch%testcase .ge. 80 .and. switch%testcase .le. 89) then
        CALL SetPuff()
        !endif
     ENDIF


     ! Allocation and initialization of the elemental matrices
     CALL init_elmat()

     matK%start = .TRUE.

     ! Initialize shock capturing
     IF ((switch%shockcp .GT. 0) .OR. ((adapt%adaptivity) .AND. (adapt%shockcp_adapt .GT. 0)))  THEN
        CALL initializeShockCapturing()
     ENDIF

     ! Compute first equation (definition of the gradient)
     CALL HDG_precalculatedfirstequation()
     IF(MPIvar%glob_id .EQ. 0) THEN
        WRITE (count_adapt_char, *) count_adapt
        CALL HDF5_save_mesh("./res/new_mesh_n" // TRIM(ADJUSTL(count_adapt_char)) // ".h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
        CALL HDF5_save_solution("./res/projected_solution_n" // TRIM(ADJUSTL(count_adapt_char)))
     ENDIF
     !! UPDATE VARIABLES
     errNR_adapt = 1e10
     time%t = time%t - time%dt
     time%it = time%it - 1
     time%ik = time%ik - 1
     sol%Nt = sol%Nt - 1
  ENDIF

  sol%u0 = 0.
  sol%u0(:, 1) = sol%u
  dt0 = time%dt
  switch_save = 0
  ir_check = 0
  it0 = 1
  IF (switch%ME .AND. time%it .NE. 0) it0 = time%it
  IF (switch%ME) THEN
     time%dt = time%dt_ME/simpar%refval_time
  ENDIF
  dt0 = time%dt
  !time%it = time%it+34
  !*******************************************************
  !                  TIME LOOP
  !*******************************************************
  DO it = it0, nts ! ************ TIME LOOP *********************


     IF(SIZE(sol%u_conv) .NE. SIZE(sol%u)) THEN
        DEALLOCATE(sol%u_conv)
        ALLOCATE(sol%u_conv(SIZE(sol%u)))
     ENDIF
     IF(SIZE(sol%q_conv) .NE. SIZE(sol%q)) THEN
        DEALLOCATE(sol%q_conv)
        ALLOCATE(sol%q_conv(SIZE(sol%q)))
     ENDIF
     sol%u_conv = sol%u
     sol%q_conv = sol%q
     !CALL deep_copy_mesh_struct(Mesh, Mesh_prec)

     ! Compute time step
     !time%dt = dt0*2**(it/100.)

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
     IF(SIZE(uiter) .NE. SIZE(sol%u0(:,1))) THEN
        DEALLOCATE(uiter)
        ALLOCATE(uiter(SIZE(sol%u0(:,1))))
     ENDIF
     uiter = sol%u0(:, 1)
     ir = 1
     DO WHILE(ir .LE. numer%nrp) ! ************ NEWTON-RAPHSON LOOP *********************
        numer%dumpnr = numer%dumpnr_min+(numer%dumpnr_max-numer%dumpnr_min)/2.*(1+TANH((ir-numer%dumpnr_n0)/numer%dumpnr_width))
        IF (MPIvar%glob_id .EQ. 0) THEN
           WRITE (6, *) "***** NR iteration: ", ir, "*****"
           WRITE (6, *) "NR dumping factor: ",  numer%dumpnr
        ENDIF

        ! Compute Jacobian
        CALL HDG_computeJacobian()
        ! Set boundary conditions
        CALL hdg_BC()
        ! Compute elemental mapping
        CALL hdg_Mapping()
        ! Assembly the global matrix
        CALL hdg_Assembly()
        !IF (switch_save.EQ.0) THEN
        !  WRITE (6, *) "Save matrix"
        !        call HDF5_save_CSR_matrix('Mat')
        !        call HDF5_save_CSR_vector('rhs')
        !
        !        stop
        !  switch_save = 1
        !call displayMatrixInt(Mesh%F)
        !call displayMatrixInt(Mesh%extfaces)
        !call displayVectorInt(Mesh%periodic_faces)
        !      stop
        !      if (ir==10) then
        !      call print_matrices_hdf5
        !      stop
        !      endif
        !ENDIF
        ! Solve linear system
        CALL solve_global_system(ir)

        ! Compute element-by-element solution
        CALL compute_element_solution()
        ! Check for NaN (doesn't work with optimization flags)
        DO i = 1, SIZE(sol%u)
           IF(ISNAN(sol%u(i))) THEN
              !IF (sol%u(i) /= sol%u(i)) THEN
              WRITE (6, *) "NaN detected"
              STOP
           END IF
        END DO

        ! Apply threshold
        !      CALL HDG_applyThreshold(mkelms)

        ! Apply filtering
        !      CALL HDG_FilterSolution()

        !!!!! ADAPTIVITY
        CALL unique_2D(Mesh%T(:,1:RefElPol%Nvertices),vector_nodes_unique)
        N_n_vertex = SIZE(vector_nodes_unique)
        IF(.NOT. ALLOCATED(error_oscillation)) THEN
           ALLOCATE(error_oscillation(N_n_vertex))
        ENDIF
        IF(SIZE(error_oscillation) .NE. N_n_vertex) THEN
           DEALLOCATE(error_oscillation)
           ALLOCATE(error_oscillation(N_n_vertex))
        ENDIF
        error_oscillation = 0.

        IF(.NOT. ALLOCATED(oscillations)) THEN
           ALLOCATE(oscillations(Mesh%Nelems))
        ENDIF
        IF(SIZE(oscillations) .NE. Mesh%Nelems) THEN
           DEALLOCATE(oscillations)
           ALLOCATE(oscillations(Mesh%Nelems))
        ENDIF
        oscillations = -100.

        CALL check_oscillations(adapt%thr_ind, error_oscillation, oscillations)
        !CALL relative_difference_solution(sol%u,uiter,flag)

        IF (utils%timing) THEN
           CALL cpu_TIME(timing%tps1)
           CALL system_CLOCK(timing%cks1, timing%clock_rate1)
        END IF


        IF (utils%timing) THEN
           CALL cpu_TIME(timing%tpe1)
           CALL system_CLOCK(timing%cke1, timing%clock_rate1)
           timing%runtadapt = timing%runtadapt + (timing%cke1-timing%cks1)/REAL(timing%clock_rate1)
           timing%cputadapt = timing%cputadapt + timing%tpe1-timing%tps1
        END IF

        max_osc = MAXVAL(oscillations)
        min_osc = MINVAL(oscillations)
        n_osc = COUNT((oscillations .LT. 0.) .AND. (oscillations .GT. -100.))

#ifdef PARALL
        CALL MPI_Allreduce(MPI_IN_PLACE, max_osc, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
        CALL MPI_Allreduce(MPI_IN_PLACE, min_osc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
        CALL MPI_Allreduce(MPI_IN_PLACE, n_osc, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
#ifdef PARALL
        IF(MPIvar%glob_id .EQ. 0) THEN
#endif
           WRITE(*,*) "MAX ERROR OSCILLATION:       ", max_osc
           !WRITE(*,*) "MIN ERROR OSCILLATION:       ", min_osc
           WRITE(*,*) "NUMBER OF OSCILLATIONS:      ", n_osc

           IF(max_osc .LE. adapt%osc_check) THEN
              WRITE(*,*) "Solution saved as checkpoint."
              IF(SIZE(sol%u_conv) .NE. SIZE(sol%u)) THEN
                 DEALLOCATE(sol%u_conv)
                 DEALLOCATE(sol%q_conv)
                 ALLOCATE(sol%u_conv(SIZE(sol%u)))
                 ALLOCATE(sol%q_conv(SIZE(sol%q)))
              ENDIF
              sol%u_conv = sol%u
              sol%q_conv = sol%q
              ir_check = ir
              CALL deep_copy_mesh_struct(Mesh, Mesh_prec)
           ENDIF

#ifdef PARALL
        ENDIF
#endif
        DEALLOCATE(vector_nodes_unique)


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
           errNR_adapt = 1e10
           ir_adapt = 0
           WRITE(*,*) "Solution saved as checkpoint."
           IF(SIZE(sol%u_conv) .NE. SIZE(uiter_best)) THEN
              DEALLOCATE(sol%u_conv)
              DEALLOCATE(sol%q_conv)
              ALLOCATE(sol%u_conv(SIZE(uiter_best)))
              ALLOCATE(sol%q_conv(SIZE(qiter_best)))
           ENDIF
           sol%u_conv = uiter_best
           sol%q_conv = qiter_best
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
              WRITE(*,*) "Solution saved as last checkpoint."
              IF(SIZE(uiter_best) .NE. SIZE(sol%u)) THEN
                 DEALLOCATE(uiter_best)
                 DEALLOCATE(qiter_best)
                 ALLOCATE(uiter_best(SIZE(sol%u)))
                 ALLOCATE(qiter_best(SIZE(sol%q)))
                 uiter_best = sol%u
                 qiter_best = sol%q
              ENDIF
              divergence_counter_adapt = 0
           ELSEIF(errNR .GT. errNR_adapt) THEN
              divergence_counter_adapt = divergence_counter_adapt + 1

           ENDIF

           WRITE(*,*) "ir_check: ", ir_check
           IF ((adapt%adaptivity) .AND. (ir .GT. (ir_check+1)) .AND. ((adapt%osc_adapt .AND. (MAXVAL(oscillations) .GT. adapt%osc_tol)) .OR. ((adapt%NR_adapt) .AND. (MOD(ir,adapt%freq_NR_adapt) .EQ. 0))))THEN !  .or. (flag)) THEN
              IF(.NOT. ((adapt%NR_adapt) .AND. (MOD(ir,adapt%freq_NR_adapt) .EQ. 0))) THEN
                 ALLOCATE(u0_temp(SIZE(sol%u0,1),SIZE(sol%u0,2)))
                 u0_temp = sol%u0
              ELSE
                 ALLOCATE(u0_temp(SIZE(sol%u_conv),SIZE(sol%u0,2)))
                 u0_temp = 0
                 u0_temp(:,1) = sol%u_conv
              ENDIF

              WRITE(*,*) "Residue before mapping:", computeResidual(sol%u, uiter, 1.)/numer%dumpnr

              CALL adaptivity


              DEALLOCATE(sol%u0)
              ALLOCATE(sol%u0(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq,time%tis))
              sol%u0 = 0.

              Mesh%X = Mesh%X*phys%lscale
              Mesh_prec%X = Mesh_prec%X*phys%lscale

              CALL projectSolutionDifferentMeshes_general_arrays(Mesh_prec%T,Mesh_prec%X, Mesh%T, Mesh%X, u1 = u0_temp(:,1), u2 = sol%u0(:,1))

              Mesh%X = Mesh%X/phys%lscale
              Mesh_prec%X = Mesh_prec%X/phys%lscale

              DEALLOCATE(u0_temp)
              DEALLOCATE(uiter)
              ALLOCATE(uiter(SIZE(sol%u)))
              uiter = sol%u0(:,1)

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
     DO i = 1, SIZE(sol%u)
        !IF (sol%u(i) /= sol%u(i)) THEN
        IF(ISNAN(sol%u(i))) THEN
           WRITE (6, *) "NaN detected"
           STOP
        END IF
     END DO

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

              ! Update the diffusion, the elemental matrices and the solution
              IF (MPIvar%glob_id .EQ. 0) THEN
                 WRITE (6, *) "************************************************"
                 WRITE (6, *) "Reducing diffusion: ", phys%diff_n*switch%diffred*simpar%refval_diffusion

#ifndef NEUTRALP
#ifdef NEUTRAL
                 WRITE (6, *) "Neutrals diffusion: ", phys%diff_nn!*switch%diffred*simpar%refval_diffusion
#ifdef KEQUATION
                 WRITE (6, *) "K diffusion min: ", phys%diff_k_min*switch%diffred*simpar%refval_diffusion
#endif
#endif
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
#ifndef NEUTRALP
#ifdef NEUTRAL
              phys%diff_nn = phys%diff_nn!*switch%diffred
#ifdef KEQUATION
              phys%diff_k_min = phys%diff_k_min*switch%diffred
              phys%diff_k_max = phys%diff_k_max!*switch%diffred
#endif
#endif
#endif

              !**********************************
              !           UPDATE SOLUTION
              !**********************************
              ! Update u0
              IF (time%tis .GT. 1 .AND. it .LT. time%tis) THEN
                 DO is = it, 1, -1
                    sol%u0(:, is + 1) = sol%u0(:, is)
                 END DO
              ELSEIF (time%tis > 1) THEN
                 DO is = time%tis, 2, -1
                    sol%u0(:, is) = sol%u0(:, is - 1)
                 END DO
              END IF
              sol%u0(:, 1) = sol%u

              time%it = 0

              IF(SIZE(sol%u_conv) .NE. SIZE(sol%u)) THEN
                 DEALLOCATE(sol%u_conv)
                 ALLOCATE(sol%u_conv(SIZE(sol%u)))
              ENDIF
              IF(SIZE(sol%q_conv) .NE. SIZE(sol%q)) THEN
                 DEALLOCATE(sol%q_conv)
                 ALLOCATE(sol%q_conv(SIZE(sol%q)))
              ENDIF
              sol%u_conv = sol%u
              sol%q_conv = sol%q
              CALL deep_copy_mesh_struct(Mesh, Mesh_prec)

              IF((adapt%adaptivity) .AND. (adapt%time_adapt) .AND. (MOD(it,adapt%freq_t_adapt) .EQ. 0)) THEN
                 !IF((MOD(it,int(50/time%dt)) .eq. 0) .and. (adapt%adaptivity)) THEN
                 CALL adaptivity
                 DEALLOCATE(sol%u0)
                 ALLOCATE (sol%u0(SIZE(sol%u), time%tis))
                 sol%u0 = 0.
                 sol%u0(:, 1) = sol%u
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
           !**********************************
           !           UPDATE SOLUTION
           !**********************************
           ! Update u0
           IF (time%tis .GT. 1 .AND. it .LT. time%tis) THEN
              DO is = it, 1, -1
                 sol%u0(:, is + 1) = sol%u0(:, is)
              END DO
           ELSEIF (time%tis > 1) THEN
              DO is = time%tis, 2, -1
                 sol%u0(:, is) = sol%u0(:, is - 1)
              END DO
           END IF
           sol%u0(:, 1) = sol%u

           IF((adapt%adaptivity) .AND. (adapt%time_adapt) .AND. (MOD(it,adapt%freq_t_adapt) .EQ. 0)) THEN
              CALL adaptivity
              DEALLOCATE(sol%u0)
              ALLOCATE (sol%u0(SIZE(sol%u), time%tis))
              sol%u0 = 0.
              sol%u0(:, 1) = sol%u
              CALL deep_copy_mesh_struct(Mesh, Mesh_prec)
           ENDIF

           ! if moving equilibrium case then update the magnetic field, otherwise just continue
           IF(switch%ME) THEN
              ! ReLoad magnetic field and Jtor
              CALL load_magnetic_field()
              CALL loadJtorMap()
              !if (switch%testcase .ge. 80 .and. switch%testcase .le. 89) then
              CALL SetPuff()
              !endif
              IF (switch%ME) THEN
                 time%dt = time%dt_ME/simpar%refval_time
              ENDIF
           ENDIF

           ! compute dt
           ! CALL compute_dt(errlstime)
        END IF
     ELSE
        IF (switch%psdtime) THEN
           ! Save solution
           CALL setSolName(save_name, mesh_name, it, .TRUE., .TRUE.)
           CALL HDF5_save_solution(save_name)

           ! Update the diffusion, the elemental matrices and the solution
           IF (MPIvar%glob_id .EQ. 0) THEN
              WRITE (6, *) "************************************************"
              WRITE (6, *) "Reducing diffusion: ", phys%diff_n*switch%diffred*simpar%refval_diffusion

#ifndef NEUTRALP
#ifdef NEUTRAL
              WRITE (6, *) "Neutrals diffusion: ", phys%diff_nn!*switch%diffred*simpar%refval_diffusion
#ifdef KEQUATION
              WRITE (6, *) "K diffusion min: ", phys%diff_k_min!*switch%diffred*simpar%refval_diffusion
#endif
#endif
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
#ifndef NEUTRALP
#ifdef NEUTRAL
           phys%diff_nn = phys%diff_nn!*switch%diffred
#ifdef KEQUATION
           phys%diff_k_min = phys%diff_k_min*switch%diffred
           phys%diff_k_max = phys%diff_k_max!*switch%diffred
#endif
#endif
#endif

           !**********************************
           !           UPDATE SOLUTION
           !**********************************
           ! Update u0
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

  IF (switch%testcase < 5) THEN
     ALLOCATE (L2err(phys%neq))
     CALL computeL2ErrorAnalyticSol(L2err)
     WRITE (6, *) " "
     DO i = 1, phys%Neq
        WRITE (6, '(A,I1,A,ES16.5)') "L2 error in U(", i, ") = ", L2err(i)
     END DO
     DEALLOCATE (L2err)
  END IF

  DEALLOCATE (uiter)
  DEALLOCATE(uiter_best)
  DEALLOCATE(qiter_best)
  DEALLOCATE (mkelms)
  DEALLOCATE (sol%u0)
  IF(ALLOCATED(oscillations)) DEALLOCATE(oscillations)
  IF(ALLOCATED(error_oscillation)) DEALLOCATE(error_oscillation)

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

CONTAINS

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
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmax, phys%npv, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmin, phys%npv, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
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
    CALL mpi_allreduce(MPI_IN_PLACE, sum2, 1, mpi_double_precision, mpi_sum, MPI_COMM_WORLD, ierr)
    CALL mpi_allreduce(nu, nglo, 1, mpi_integer, mpi_sum, MPI_COMM_WORLD, ierr)
#endif

    res = SQRT(sum2)/SQRT(DBLE(nglo))/coeff
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
    save_name = 'Sol3D_'//save_name
#else
    save_name = 'Sol2D_'//save_name
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

  SUBROUTINE adaptivity

    IF (utils%timing) THEN
       CALL cpu_TIME(timing%tps1)
       CALL system_CLOCK(timing%cks1, timing%clock_rate1)
    END IF

    !! START REFINING PROCEDURE
    count_adapt = count_adapt + 1

    CALL deep_copy_mesh_struct(Mesh,Mesh_prec)
    CALL deep_copy_refel_struct(refElPol,refElPol_prec)

    Mesh%X = Mesh%X*phys%lscale
    Mesh_prec%X = Mesh_prec%X*phys%lscale
    !DEALLOCATE ALL BEFORE ADAPTIVITY
    DEALLOCATE(phys%B)
    DEALLOCATE(phys%magnetic_flux)
    IF(ASSOCIATED(phys%Bperturb)) DEALLOCATE(phys%Bperturb)
    IF(ASSOCIATED(phys%Jtor))     DEALLOCATE(phys%Jtor)
    IF(ASSOCIATED(phys%puff_exp)) DEALLOCATE(phys%puff_exp)
    CALL free_el_mat()
    CALL free_mat
    IF(adapt%evaluator .EQ. 0) THEN
       CALL adaptivity_indicator_estimator(mesh_name, adapt%thr_ind, adapt%param_est, count_adapt, order)
    ELSEIF(adapt%evaluator .EQ. 1) THEN
       CALL adaptivity_indicator(mesh_name, adapt%thr_ind, adapt%param_est, count_adapt, order)
    ELSEIF(adapt%evaluator .EQ. 2) THEN
       CALL adaptivity_estimator(mesh_name, adapt%param_est, count_adapt, order)
    ELSE
       WRITE(*,*) "Choice of adaptivity evaluator not valid. STOP."
    ENDIF


    ! CALL HDF5_save_mesh("./newmesh.h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)


    CALL HDF5_save_mesh("./old_mesh.h5", Mesh_prec%Ndim, Mesh_prec%Nelems, Mesh_prec%Nextfaces, Mesh_prec%Nnodes, Mesh_prec%Nnodesperelem, Mesh_prec%Nnodesperface, Mesh_prec%elemType, Mesh_prec%T, Mesh_prec%X, Mesh_prec%Tb, Mesh_prec%boundaryFlag)
    ! ALLOCATE(u_proj_1D(SIZE(sol%u)))
    ! u_proj_1D = sol%u
    ! sol%u = sol%u_conv
    ! CALL HDF5_save_solution("./projecting_solution")
    ! sol%u = u_proj_1D
    ! DEALLOCATE(u_proj_1D)

    ! call HDF5_create('./fortran_save_0.h5', file_id, ierr)
    ! call HDF5_array2D_saving_int(file_id, Mesh_prec%T, size(Mesh_prec%T,1), size(Mesh_prec%T,2), 'T_prec_f')
    ! call HDF5_array2D_saving(file_id, Mesh_prec%X, size(Mesh_prec%X,1), size(Mesh_prec%X,2), 'X_prec_f')
    ! call HDF5_array2D_saving_int(file_id, Mesh%T, size(Mesh%T,1), size(Mesh%T,2), 'T_f')
    ! call HDF5_array2D_saving(file_id, Mesh%X, size(Mesh%X,1), size(Mesh%X,2), 'X_f')
    ! call HDF5_array2D_saving_int(file_id, Mesh_prec%Tb, size(Mesh_prec%Tb,1), size(Mesh_prec%Tb,2), 'Tb_prec_f')
    ! call HDF5_array2D_saving_int(file_id, Mesh%Tb, size(Mesh%Tb,1), size(Mesh%Tb,2), 'Tb_f')
    ! call HDF5_close(file_id)

    CALL projectSolutionDifferentMeshes_general(Mesh_prec%T, Mesh_prec%X, Mesh%T, Mesh%X, sol%u_conv, sol%q_conv, sol%u, sol%q)
    DEALLOCATE(sol%u_tilde)
    DEALLOCATE(sol%u_tilde0)
    ALLOCATE(sol%u_tilde(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
    ALLOCATE(sol%u_tilde0(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))

    CALL extractFaceSolution()
    sol%u_tilde0 = sol%u_tilde

    DEALLOCATE(sol%u_conv)
    ALLOCATE(sol%u_conv(SIZE(sol%u)))
    sol%u_conv = sol%u
    DEALLOCATE(sol%q_conv)
    ALLOCATE(sol%q_conv(SIZE(sol%q)))
    sol%q_conv = sol%q

    WRITE(*,*) "Number of elements previous mesh: ", SIZE(Mesh_prec%T,1)
    WRITE(*,*) "Number of elements current mesh:  ", SIZE(Mesh%T,1)

    Mesh%X = Mesh%X/phys%lscale
    Mesh_prec%X = Mesh_prec%X/phys%lscale




    ! Initialize magnetic field (the Mesh is needed)
    CALL initialize_magnetic_field()

    ! Load magnetic field and, if case is 54 or 59
    CALL load_magnetic_field()
    IF (switch%ohmicsrc) CALL loadJtorMap()

    ! if moving equilibrium case then update the magnetic field, otherwise just continue
    IF(switch%ME) THEN
       CALL loadJtorMap()
       !if (switch%testcase .ge. 80 .and. switch%testcase .le. 89) then
       CALL SetPuff()
       !endif
    ENDIF



    ! Allocation and initialization of the elemental matrices
    CALL init_elmat()

    matK%start = .TRUE.

    ! Initialize shock capturing
    IF ((switch%shockcp .GT. 0) .OR. (adapt%shockcp_adapt .GT. 0))  THEN
       CALL initializeShockCapturing()
    ENDIF

    ! Compute first equation (definition of the gradient)
    CALL HDG_precalculatedfirstequation()

    WRITE (count_adapt_char, *) count_adapt
    CALL HDF5_save_mesh("./res/new_mesh_n" // TRIM(ADJUSTL(count_adapt_char)) // ".h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
    CALL HDF5_save_solution("./res/projected_solution_n" // TRIM(ADJUSTL(count_adapt_char)))

    !! UPDATE VARIABLES
    errNR_adapt = 1e10

    IF (utils%timing) THEN
       CALL cpu_TIME(timing%tpe1)
       CALL system_CLOCK(timing%cke1, timing%clock_rate1)
       timing%runtadapt = timing%runtadapt + (timing%cke1-timing%cks1)/REAL(timing%clock_rate1)
       timing%cputadapt = timing%cputadapt + timing%tpe1-timing%tps1
    END IF


  ENDSUBROUTINE adaptivity

END PROGRAM MHDG
