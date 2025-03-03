
MODULE adaptivity_estimator_indicator_module
  USE globals
  USE reference_element
  USE gmsh
  USE adaptivity_common_module
  USE adaptivity_indicator_module
  USE adaptivity_estimator_module
  USE MPI_OMP
  IMPLICIT NONE

CONTAINS

  SUBROUTINE adaptivity_indicator_estimator(mesh_name,thresh, param_adapt, count_adapt, order)
    USE in_out, ONLY: copy_file
    USE gmsh_io_module, ONLY: load_gmsh_mesh, HDF5_save_mesh, convert_gmsh_to_hdf5
    USE preprocess
#ifdef PARALL
    USE Communications, ONLY: gather_1D_vector_int,gather_1D_vector_real,gather_mesh
#endif

    TYPE(gmsh_t)                                :: gmsh
    REAL*8, INTENT(IN)                          :: thresh
    INTEGER, INTENT(IN)                         :: param_adapt, count_adapt, order

    REAL*8,  ALLOCATABLE                        :: u_sol(:,:), u_star_sol(:,:), h(:), h_target(:), h_target_temp(:), error_oscillation(:), error_L2(:), error_L2_vertices(:), error_target(:)
    INTEGER, ALLOCATABLE                        :: vector_nodes_unique(:)
    INTEGER                                     :: i, N_n_vertex, n_el_unstable, ierr
    REAL*8                                      :: eps_plot(Mesh%Nnodes)
    REAL*8                                      :: eg_L2

    REAL*8, POINTER                             :: h_target_on_nodes(:), nodes_glob(:,:)
    INTEGER, POINTER                            :: connectivity_glob(:,:)
#ifdef PARALL
    REAL*8, ALLOCATABLE                         :: h_root(:), error_L2_vertices_root(:), error_oscillation_root(:)
    REAL*8, POINTER                             :: h_glob(:), error_oscillation_glob(:), error_L2_vertices_glob(:)
    INTEGER, POINTER                            :: vector_nodes_unique_glob(:), count_vec_glob(:)
    INTEGER, ALLOCATABLE                        :: count_vec(:), noghost_index(:)
#endif

    CHARACTER(70)                               :: param_adapt_char, count_adapt_char
    ! mesh_name is the mesh path + mesh name + .msh extension ("./Meshes/CircLim.msh")
    ! mesh_name_npne (mesh name no path no extension) is just the name of the mesh ("CircLim")
    ! new_mesh_name_npne (new mesh name no path no extension) is just the name of the mesh + param_adapt + count_adapt ("CircLim_param2_n1")
    ! buffer is a dummy array to store intermediate mesh names
    CHARACTER(1024), INTENT(IN)                 :: mesh_name
    CHARACTER(1024)                             :: mesh_name_npne, new_mesh_name_npne, buffer

    NULLIFY(h_target_on_nodes, nodes_glob, connectivity_glob)
#ifdef PARALL
    NULLIFY(h_glob, error_oscillation_glob, error_L2_vertices_glob, vector_nodes_unique_glob, count_vec_glob)
#endif

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) "*************************************************"
       WRITE(*,*) "        ADAPTIVITY ESTIMATOR-INDICATOR           "
       WRITE(*,*) "*************************************************"
    ENDIF

#ifdef PARALL
    ALLOCATE(noghost_index(Mesh%Nelems-Mesh%nghostelems))
    ! only select the indices of the non-ghost elements
    noghost_index = PACK([(i, i=1, Mesh%Nelems)], Mesh%ghostElems(:) .EQ. 0)
    CALL unique_1D(RESHAPE(Mesh%T(noghost_index,1:refElPol%Nvertices), [SIZE(Mesh%T(noghost_index,1:refElPol%Nvertices),1) * SIZE(Mesh%T(noghost_index,1:refElPol%Nvertices),2)]), vector_nodes_unique)
    CALL gather_1D_vector_int(Mesh%loc2glob_nodes(vector_nodes_unique), vector_nodes_unique_glob, allgather = .FALSE.)
    DEALLOCATE(noghost_index)
#else
    CALL unique_1D(RESHAPE(Mesh%T(:,1:refElPol%Nvertices), [SIZE(Mesh%T(:,1:refElPol%Nvertices),1) * SIZE(Mesh%T(:,1:refElPol%Nvertices),2)]), vector_nodes_unique)
#endif

    N_n_vertex = SIZE(vector_nodes_unique)

    ALLOCATE(h(N_n_vertex))
    ALLOCATE(error_L2_vertices(N_n_vertex))

#ifndef PARALL

    ALLOCATE(error_target(N_n_vertex))
    ALLOCATE(h_target(N_n_vertex))
    error_target = adapt%tol_est
    h_target = 100.
#else

    ALLOCATE(count_vec(N_n_vertex))
    IF(MPIvar%glob_id .EQ. 0) THEN
       N_n_vertex = MAXVAL(vector_nodes_unique_glob)
       ALLOCATE(error_L2_vertices_root(N_n_vertex))
       ALLOCATE(error_target(N_n_vertex))
       ALLOCATE(h_target(N_n_vertex))
       error_L2_vertices_root = 0
       error_target = adapt%tol_est
       h_target = 100.
       N_n_vertex = SIZE(vector_nodes_unique)
    ENDIF

#endif
    ALLOCATE(h_target_temp(N_n_vertex))
    ALLOCATE(error_L2(SIZE(Mesh%T,1)))
    ALLOCATE(error_oscillation(N_n_vertex))
    h = 0.
    error_oscillation = 0.
    error_L2_vertices = 0.
    h_target_temp = 0.
    error_L2 = 0.
    eg_L2 = 0.

    !! COMPUTE h_map !!
#ifndef PARALL
    CALL h_map(N_n_vertex,Mesh%X(vector_nodes_unique,1:2),Mesh%T(:,1:3),vector_nodes_unique, h)
#else
    CALL h_map(N_n_vertex,Mesh%X(vector_nodes_unique,1:2),Mesh%T(:,1:3),vector_nodes_unique, h, count_vec)

    CALL gather_1D_vector_real(h, h_glob, allgather = .FALSE.)
    CALL gather_1D_vector_int(count_vec, count_vec_glob, allgather = .FALSE.)

    IF(MPIvar%glob_id .EQ. 0) THEN
       ALLOCATE(h_root(MAXVAL(vector_nodes_unique_glob)))
       h_root = 0.
       CALL compute_error_on_vertices_root(h_glob, vector_nodes_unique_glob, count_vec_glob, h_root)
    ENDIF
#endif
    !! END OF COMPUTING h_map !!

    !! ESTIMATOR !!
    IF(param_adapt .EQ. 0) THEN
       ! u_sol, u_star_sol are allocated here
       CALL post_process_matrix_solution(Mesh%X,Mesh%T,sol%u,sol%q,u_sol,u_star_sol)

       DO i = 1, phys%npv
          ! error estimation for the mesh and the solution
          CALL calculate_L2_error_two_sols_different_p_scalar_general(Mesh%X,Mesh%T,i, u_sol,u_star_sol, error_L2, eg_L2)
          CALL error_on_vertices(error_L2,Mesh%T,vector_nodes_unique, N_n_vertex, error_L2_vertices)

#ifdef PARALL
          CALL gather_1D_vector_real(error_L2_vertices, error_L2_vertices_glob, allgather = .FALSE.)

          IF(MPIvar%glob_id .EQ. 0) THEN
             CALL compute_error_on_vertices_root(error_L2_vertices_glob, vector_nodes_unique_glob, count_vec_glob, error_L2_vertices_root)

             ! richardson formula only for estimator
             h_target_temp = h_root * ((error_target / error_L2_vertices_root) ** (1./ (order + 1.)))
             h_target = MIN(h_target_temp, h_target)
          ENDIF
          DEALLOCATE(error_L2_vertices_glob)
          NULLIFY(error_L2_vertices_glob)
#else
          ! richardson formula only for estimator
          h_target_temp = EXP( ( LOG(error_target) - LOG( error_L2_vertices ) )/(order+1) + LOG(h) )
          h_target = MIN(h_target_temp, h_target)
#endif
       ENDDO
       DEALLOCATE(u_sol,u_star_sol)
    ELSE
       ! error estimation for the mesh and the solution
       CALL L2_error_estimator_eval(Mesh%X,Mesh%T,sol%u,sol%q,param_adapt,error_L2,eg_L2)
       CALL error_on_vertices(error_L2,Mesh%T,vector_nodes_unique, N_n_vertex, error_L2_vertices)

#ifdef PARALL

       CALL gather_1D_vector_real(error_L2_vertices, error_L2_vertices_glob, allgather = .FALSE.)

       IF(MPIvar%glob_id .EQ. 0) THEN
          CALL compute_error_on_vertices_root(error_L2_vertices_glob, vector_nodes_unique_glob, count_vec_glob, error_L2_vertices_root)

          ! richardson formula only for estimator
          h_target = h_root * ((error_target / error_L2_vertices_root) ** (1./ (order + 1.)))
          h_target = SQRT(h_target*h) !smoothing richardson formula
       ENDIF
#else
       ! richardson formula only for estimator
       h_target = h * ((error_target / error_L2_vertices) ** (1./ (order + 1.)))
       h_target = SQRT(h_target*h) !smoothing richardson formula
#endif
    ENDIF
    !! END OF ESTIMATOR !!

    !! INDICATOR
    CALL hdg_ShockCapturing_adapt(thresh, eps_plot)
    !! error indicator based on the elemental oscillations
#ifndef PARALL
    CALL read_error(eps_plot, error_oscillation)
    WHERE (ABS(error_oscillation) .GT. 1e-10)
       h_target = h * 0.5
    END WHERE
    h_target = MIN(h_target, 0.1)
#else
    CALL read_error(eps_plot, error_oscillation, count_vec)
    CALL gather_1D_vector_real(error_oscillation, error_oscillation_glob, allgather = .FALSE.)

    IF(MPIvar%glob_id .EQ. 0) THEN
       ALLOCATE(error_oscillation_root(MAXVAL(vector_nodes_unique_glob)))
       error_oscillation_root = 0.
       CALL compute_error_on_vertices_root(error_oscillation_glob, vector_nodes_unique_glob, count_vec_glob, error_oscillation_root)

       h_target = h_root  ! Start by assigning h to h_target
       WHERE (ABS(error_oscillation_root) .GT. 1e-10)
          h_target = h_root * 0.5
       END WHERE
       h_target = MIN(h_target, 0.1)

    ENDIF
#endif
    !! END OF INDICATOR !!


#ifndef PARALL
    ALLOCATE(nodes_glob(SIZE(Mesh%X,1),SIZE(Mesh%X,2)))
    ALLOCATE(connectivity_glob(SIZE(Mesh%T,1),SIZE(Mesh%T,2)))
    nodes_glob = Mesh%X
    connectivity_glob = Mesh%T
#else
   CALL gather_mesh(Mesh,connectivity_glob,nodes_glob)
   DEALLOCATE(vector_nodes_unique)
   CALL unique_1D(vector_nodes_unique_glob,vector_nodes_unique)  
   
    IF(MPIvar%glob_id .EQ. 0) THEN
       N_n_vertex = MAXVAL(vector_nodes_unique_glob)
#endif

    ALLOCATE(h_target_on_nodes(SIZE(Mesh%X,1)))
    h_target_on_nodes = 0.5
    DO i=1,SIZE(vector_nodes_unique,1)
       h_target_on_nodes(vector_nodes_unique(i)) = h_target(i)
    ENDDO

    CALL gmsh_create_from_h_target( h_target_on_nodes, Mesh%X, Mesh%T(:,1:3), order)

#ifdef PARALL
    ENDIF
    ! wait for process 0 to finish writing before loading new mesh
    CALL MPI_BARRIER(mpi_comm_world, ierr)
#endif
    CALL free_mesh

    CALL free_reference_element_pol(refElPol)
    CALL create_reference_element(refElPol,2,order, verbose = 0)

    IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
       CALL load_gmsh_mesh("./res/temp",0)
    ELSE
       CALL load_gmsh_mesh("./res/temp",1)
    ENDIF

    CALL free_reference_element_pol(refElPol)
    CALL create_reference_element(refElPol,2,order, verbose = 0)
    CALL mesh_preprocess_serial(ierr)
    CALL read_extended_connectivity('./res/temp.msh')

    Mesh%X = Mesh%X*phys%lscale

    IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
       Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
    END IF

    n_el_unstable = 0
    DO i = 1, SIZE(error_oscillation)
       IF(error_oscillation(i) .GT. 1e-12) THEN
          n_el_unstable = n_el_unstable + 1
       ENDIF
    ENDDO

#ifdef PARALL
    CALL MPI_Allreduce(MPI_IN_PLACE, n_el_unstable, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    N_n_vertex = MAXVAL(vector_nodes_unique_glob)
#endif

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,'(a, F5.2)') "********** Percentage of oscillating elements on previous mesh: ", REAL(n_el_unstable*100)/REAL(N_n_vertex), "%"
    ENDIF

    DEALLOCATE(error_oscillation)
    DEALLOCATE(h)
    IF(ASSOCIATED(h_target_on_nodes)) DEALLOCATE(h_target_on_nodes)    
    NULLIFY(h_target_on_nodes)
    DEALLOCATE(nodes_glob)
    IF(ASSOCIATED(nodes_glob)) DEALLOCATE(nodes_glob) 
    NULLIFY(nodes_glob)
    IF(ASSOCIATED(connectivity_glob)) DEALLOCATE(connectivity_glob)
    NULLIFY(connectivity_glob)
    DEALLOCATE(vector_nodes_unique)
    DEALLOCATE(error_L2_vertices)
    DEALLOCATE(h_target_temp)
    DEALLOCATE(error_L2)

#ifndef PARALL
    DEALLOCATE(h_target)
#else
    IF(ASSOCIATED(error_L2_vertices_glob)) DEALLOCATE(error_L2_vertices_glob)

    IF(MPIvar%glob_id .EQ. 0) THEN
       DEALLOCATE(error_target)
       DEALLOCATE(error_L2_vertices_root)
    ENDIF

    DEALLOCATE(h_glob)
    DEALLOCATE(count_vec_glob)
    DEALLOCATE(error_oscillation_glob)
    DEALLOCATE(vector_nodes_unique_glob)
    DEALLOCATE(count_vec)

    IF(MPIvar%glob_id .EQ. 0) THEN
       DEALLOCATE(h_target)
       DEALLOCATE(h_root)
       DEALLOCATE(error_oscillation_root)
    ENDIF

    NULLIFY(h_glob, error_oscillation_glob, error_L2_vertices_glob, vector_nodes_unique_glob, count_vec_glob)
#endif


  END SUBROUTINE adaptivity_indicator_estimator



END MODULE adaptivity_estimator_indicator_module
