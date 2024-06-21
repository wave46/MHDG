!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! Module for schock capturing adaptivity
!************************************************************

MODULE adaptivity_estimator_indicator_module
  USE globals
  USE reference_element
  USE gmsh
  USE adaptivity_common_module
  USE adaptivity_indicator_module
  USE adaptivity_estimator_module
  IMPLICIT NONE

CONTAINS

  SUBROUTINE adaptivity_indicator_estimator(mesh_name,thresh, param_adapt, count_adapt, order)
    USE in_out, ONLY: copy_file
    USE gmsh_io_module, ONLY: load_gmsh_mesh, HDF5_save_mesh, convert_gmsh_to_hdf5
    USE preprocess

    TYPE(gmsh_t)                                :: gmsh
    REAL*8, INTENT(IN)                          :: thresh
    INTEGER, INTENT(IN)                         :: param_adapt, count_adapt, order

    REAL*8,  ALLOCATABLE                        :: two_d_nodes(:,:)
    INTEGER, ALLOCATABLE                        :: two_d_elements(:,:)
    REAL*8,  ALLOCATABLE                        :: u_sol(:,:), u_star_sol(:,:), h(:), h_target(:), h_target_temp(:), error_oscillation(:), error_L2(:), error_L2_vertices(:), error_L2_init(:), error_target(:), temp_error(:)
    INTEGER, ALLOCATABLE                        :: vector_nodes_unique(:)
    INTEGER                                     :: i, N_n_vertex, n_el_unstable, ierr
    REAL*8                                      :: eps_plot(Mesh%Nnodes)
    REAL*8                                      :: eg_L2

#ifdef PARALL
    INTEGER, ALLOCATABLE                        :: T_nogho(:,:)
    REAL*8, ALLOCATABLE                         :: q_nogho(:), u_nogho(:)
#endif


    CHARACTER(70)                               :: param_adapt_char, count_adapt_char
    ! mesh_name is the mesh path + mesh name + .msh extension ("./Meshes/CircLim.msh")
    ! mesh_name_npne (mesh name no path no extension) is just the name of the mesh ("CircLim")
    ! new_mesh_name_npne (new mesh name no path no extension) is just the name of the mesh + param_adapt + count_adapt ("CircLim_param2_n1")
    ! buffer is a dummy array to store intermediate mesh names
    CHARACTER(1024), INTENT(IN)                 :: mesh_name
    CHARACTER(1024)                             :: mesh_name_npne, new_mesh_name_npne, buffer

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) '*************** Starting refinement procedure with oscillations: ****************'
    ENDIF

    CALL hdg_ShockCapturing_adapt(thresh, eps_plot)

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) '*************** Unique  ****************'
    ENDIF

    CALL unique_1D(RESHAPE(Mesh%T(:,1:RefElPol%Nvertices), [SIZE(Mesh%T(:,1:RefElPol%Nvertices),1) * SIZE(Mesh%T(:,1:RefElPol%Nvertices),2)]), vector_nodes_unique)

    N_n_vertex = SIZE(vector_nodes_unique)

    ALLOCATE(two_d_nodes(N_n_vertex,2))
    ALLOCATE(two_d_elements(SIZE(Mesh%T,1),3))
    ALLOCATE(error_oscillation(N_n_vertex))
    ALLOCATE(h(N_n_vertex))
    ALLOCATE(error_L2_vertices(N_n_vertex))
    ALLOCATE(temp_error(N_n_vertex))
    ALLOCATE(error_target(N_n_vertex))
    ALLOCATE(h_target(N_n_vertex))
    ALLOCATE(h_target_temp(N_n_vertex))
    ALLOCATE(error_L2(SIZE(Mesh%T,1)))
    ALLOCATE(error_L2_init(SIZE(Mesh%T,1)))

    two_d_nodes = 0.
    two_d_elements = 0.
    error_oscillation = 0.
    h = 0.
    error_L2_vertices = 0.
    error_target = adapt%tol_est
    h_target = 100.
    error_L2 = 0.
    error_L2_init = 0.


    two_d_nodes = Mesh%X(vector_nodes_unique,:)
    two_d_elements = Mesh%T(:,1:3)

    !! error indicator based on the elemental oscillations
    CALL read_error(eps_plot, error_oscillation)

    !! use error map to create element size map: h
    CALL h_map(N_n_vertex,two_d_nodes,two_d_elements,vector_nodes_unique, h)

    IF(param_adapt .EQ. 0) THEN
       ! u_sol, u_star_sol are allocated here
       CALL post_process_matrix_solution(Mesh%X,Mesh%T,sol%u,sol%q,u_sol,u_star_sol)

       DO i = 1, phys%npv
          ! error estimation for the mesh and the solution
          CALL calculate_L2_error_two_sols_different_p_scalar_general(Mesh%X,Mesh%T,i, u_sol,u_star_sol, error_L2, eg_L2)
          CALL error_on_vertices(error_L2,Mesh%T,vector_nodes_unique, N_n_vertex, error_L2_vertices)
          ! richardson formula only for estimator
          h_target_temp = EXP( ( LOG(error_target) - LOG( error_L2_vertices ) )/(order+1) + LOG(h) )
          h_target = MIN(h_target_temp, h_target)
       ENDDO

       DEALLOCATE(u_sol,u_star_sol)
    ELSE
       ! error estimation for the mesh and the solution
#ifdef PARALL
       CALL L2_error_estimator_eval(Mesh%X,T_nogho,u_nogho,q_nogho,param_adapt,error_L2,eg_L2)
       CALL error_on_vertices(error_L2,T_nogho,vector_nodes_unique, N_n_vertex, error_L2_vertices)
#else
       CALL L2_error_estimator_eval(Mesh%X,Mesh%T,sol%u,sol%q,param_adapt,error_L2,eg_L2)
       CALL error_on_vertices(error_L2,Mesh%T,vector_nodes_unique, N_n_vertex, error_L2_vertices)
#endif

       ! richardson formula only for estimator
       h_target = EXP( ( LOG(error_target) - LOG( error_L2_vertices ) )/(order+1) + LOG(h) )

    ENDIF


    h_target = SQRT(h_target*h) !smoothing richardson formula

    ! create contribution of estimator and indicator for mmg
    ! if oscillations are detected just use original size/2 otherwise use
    ! estimator calculated with Richardson

    DO i=1, SIZE(h) !here we add the oscillation contribution
       IF (error_oscillation(i) .GT. 0) THEN
          ! refine dividing by 2 the element size when find oscillations
          h_target(i) = h(i)/2
       ENDIF  !if no oscillation don't do anything
    ENDDO

    DO i=1,SIZE(h) !check on coarsening ---> not higher than initial mesh
       IF (h_target(i) .GT. 0.1) THEN
          h_target(i) = 0.1
       ENDIF
    ENDDO


    CALL generate_htarget_sol_file(N_n_vertex, h_target)

    CALL extract_mesh_name_from_fullpath_woext(mesh_name, mesh_name_npne)

    WRITE(param_adapt_char, *) param_adapt
    WRITE(count_adapt_char, *) count_adapt
    new_mesh_name_npne = TRIM(ADJUSTL(mesh_name_npne)) // '_param'// TRIM(ADJUSTL(param_adapt_char)) // '_n' // TRIM(ADJUSTL(count_adapt_char))

    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".mesh"
    CALL mmg_create_mesh_from_h_target(buffer)

    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne))
    CALL convert_mesh2msh(buffer)
    CALL convert_msh2mesh(buffer)

    CALL delete_file("./res/temp.mesh")
    CALL delete_file("./res/temp.msh")
    CALL delete_file("./res/ElSizeMap.sol")


    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".mesh"
    CALL copy_file(buffer, "./res/temp.mesh")
    !CALL delete_file(buffer)

    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".msh"

    CALL open_merge_with_geometry(gmsh, buffer)
    CALL copy_file(buffer, "./res/temp.msh")
    !CALL delete_file(buffer)


    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".sol"
    CALL delete_file(buffer)

    !CALL merge_with_geometry(gmsh)

    WRITE(*,*) "********** Loading mesh P1  **********"
    CALL free_mesh
    IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
       CALL load_gmsh_mesh("./res/temp",0)
    ELSE
       CALL load_gmsh_mesh("./res/temp",1)
    ENDIF

    CALL create_reference_element(refElPol,2,1, verbose = 0)
    CALL mesh_preprocess(ierr)

    Mesh%X = Mesh%X*phys%lscale

    IF(ierr .EQ. 0) THEN
       WRITE(*,*) "Error! Corresponding face in Tb not found. STOP"
       STOP
    ENDIF

    CALL HDF5_save_mesh("./newmesh_pre.h5", Mesh%Ndim, mesh%Nelems, mesh%Nextfaces, mesh%Nnodes, mesh%Nnodesperelem, mesh%Nnodesperface, mesh%elemType, mesh%T, mesh%X, mesh%Tb, mesh%boundaryFlag)

    CALL read_extended_connectivity('./res/temp.msh')

    WRITE(*,*) "********** Increasing order mesh **********"
    CALL set_order_mesh(order)
    CALL free_reference_element
    CALL create_reference_element(refElPol,2,order, verbose = 0)
    CALL mesh_preprocess(ierr)


    Mesh%X = Mesh%X*phys%lscale

    IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
       Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
    END IF

    CALL HDF5_save_mesh("./newmesh_notround.h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)

    WRITE(*,*) "********** Rounding edges **********"
    !CALL round_edges(Mesh)

    ! overwrite the temp.msh file with the new one with rounded edges (still order 1)
    CALL write_msh_file(Mesh%X,Mesh%T)
    ! convert the mesh to .mesh
    CALL convert_msh2mesh('./res/temp')

    CALL HDF5_save_mesh("./newmesh_round.h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)


    n_el_unstable = 0

    DO i = 1, SIZE(error_oscillation)
       IF(error_oscillation(i) .GT. 1e-12) THEN
          n_el_unstable = n_el_unstable + 1
       ENDIF
    ENDDO

    WRITE(*,'(a, F5.2)') "********** Percentage of oscillating elements on previous mesh: ", REAL(n_el_unstable*100)/REAL(N_n_vertex), "%"

    DEALLOCATE(two_d_nodes)
    DEALLOCATE(two_d_elements)
    DEALLOCATE(error_oscillation)
    DEALLOCATE(h)
    DEALLOCATE(error_L2_vertices)
    DEALLOCATE(error_target)
    DEALLOCATE(h_target)
    DEALLOCATE(h_target_temp)
    DEALLOCATE(error_L2)
    DEALLOCATE(error_L2_init)
    DEALLOCATE(temp_error)
    DEALLOCATE(vector_nodes_unique)

  END SUBROUTINE adaptivity_indicator_estimator



END MODULE adaptivity_estimator_indicator_module
