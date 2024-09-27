!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! Module for schock capturing adaptivity
!************************************************************

MODULE adaptivity_estimator_module
  USE globals
  USE reference_element
  USE LinearAlgebra
  USE gmsh
  USE adaptivity_common_module
  USE adaptivity_indicator_module
  IMPLICIT NONE

CONTAINS

  SUBROUTINE adaptivity_estimator(mesh_name, param_adapt, count_adapt, order)
    USE in_out, ONLY: copy_file
    USE gmsh_io_module, ONLY: load_gmsh_mesh, HDF5_save_mesh, convert_gmsh_to_hdf5
    USE preprocess


    TYPE(gmsh_t)                                :: gmsh
    INTEGER, INTENT(IN)                         :: param_adapt, count_adapt, order

    REAL*8,  ALLOCATABLE                        :: two_d_nodes(:,:)
    INTEGER, ALLOCATABLE                        :: two_d_elements(:,:)
    REAL*8,  ALLOCATABLE                        :: u_sol(:,:), u_star_sol(:,:), h(:), h_target(:),h_target_temp(:), error_L2(:), error_L2_vertices(:), error_L2_init(:), error_target(:)
    INTEGER, ALLOCATABLE                        :: vector_nodes_unique(:)
    INTEGER                                     :: i, N_n_vertex,ierr
    REAL*8                                      :: eg_L2, eg_L2_init

#ifdef PARALL
    INTEGER                                     :: N_n_vertex_glob, counter, j, jj, ii
    REAL*8                                      :: dummy_sum
    REAL*8, ALLOCATABLE                         :: h_target_glob_rep(:), h_target_glob(:)
    INTEGER, ALLOCATABLE                        :: T_nogho(:,:)
    REAL*8, ALLOCATABLE                         :: u_nogho(:), q_nogho(:)
    INTEGER, ALLOCATABLE                        :: recvcounts(:), displs(:), vector_nodes_unique_glob(:), temp(:)
    LOGICAL, ALLOCATABLE                        :: already_counted(:)
#endif

    CHARACTER(70)                               :: param_adapt_char, count_adapt_char
    ! mesh_name is the mesh path + mesh name + .msh extension ("./Meshes/CircLim.msh")
    ! mesh_name_npne (mesh name no path no extension) is just the name of the mesh ("CircLim")
    ! new_mesh_name_npne (new mesh name no path no extension) is just the name of the mesh + param_adapt + count_adapt ("CircLim_param2_n1")
    ! buffer is a dummy array to store intermediate mesh names
    CHARACTER(1024), INTENT(IN)                 :: mesh_name
    CHARACTER(1024)                             :: mesh_name_npne, new_mesh_name_npne, buffer

    !!  Refinement part with indicator
    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) '*************** Starting refinement procedure  ****************'
    ENDIF

#ifdef PARALL

    ! remove ghost elements from T, u and q
    ALLOCATE(T_nogho(Mesh%Nelems-Mesh%nghostelems, Mesh%Nnodesperelem))
    ALLOCATE(u_nogho((Mesh%Nelems-Mesh%nghostelems)*phys%Neq*Mesh%Nnodesperelem))
    ALLOCATE(q_nogho((Mesh%Nelems-Mesh%nghostelems)*phys%Neq*Mesh%Ndim*Mesh%Nnodesperelem))

    counter = 1
    DO i = 1, SIZE(Mesh%T,1)
       IF(Mesh%ghostElems(i) .NE. 1) THEN
          T_nogho(counter,:) = Mesh%T(i,:)
          counter = counter + 1
       ENDIF
    ENDDO

    counter = 1
    DO i = 1, SIZE(Mesh%T,1)
       IF(Mesh%ghostElems(i) .NE. 1) THEN
          DO j = 1, Mesh%Nnodesperelem
             DO jj = 1, phys%neq
                u_nogho(counter) = sol%u((i-1)*Mesh%Nnodesperelem*phys%neq+(j-1)*phys%neq+jj)
                counter = counter + 1
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    counter = 1
    DO i = 1, SIZE(Mesh%T,1)
       IF(Mesh%ghostElems(i) .NE. 1) THEN
          DO j = 1, Mesh%Nnodesperelem
             DO ii = 1, phys%neq
                DO jj = 1, Mesh%Ndim
                   q_nogho(counter) = sol%q((i-1)*Mesh%Nnodesperelem*phys%neq*Mesh%Ndim + (j-1)*phys%neq + (ii-1)*Mesh%Ndim+jj)
                   counter = counter + 1
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) '*************** Unique  ****************'
    ENDIF
    CALL unique_1D(RESHAPE(T_nogho(:,1:refElPol%Nvertices), [SIZE(T_nogho(:,1:refElPol%Nvertices),1) * SIZE(T_nogho(:,1:refElPol%Nvertices),2)]), vector_nodes_unique)

#else
    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) '*************** Unique  ****************'
    ENDIF
    CALL unique_1D(RESHAPE(Mesh%T(:,1:refElPol%Nvertices), [SIZE(Mesh%T(:,1:refElPol%Nvertices),1) * SIZE(Mesh%T(:,1:refElPol%Nvertices),2)]), vector_nodes_unique)
#endif

    N_n_vertex = SIZE(vector_nodes_unique)

    ALLOCATE(two_d_nodes(N_n_vertex,2))
    ALLOCATE(h(N_n_vertex))
    ALLOCATE(error_L2_vertices(N_n_vertex))
    ALLOCATE(error_target(N_n_vertex))
    ALLOCATE(h_target(N_n_vertex))
    ALLOCATE(h_target_temp(N_n_vertex))
#ifndef PARALL
    ALLOCATE(two_d_elements(SIZE(Mesh%T,1),3))
    ALLOCATE(error_L2(SIZE(Mesh%T,1)))
    ALLOCATE(error_L2_init(SIZE(Mesh%T,1)))
#else
    ALLOCATE(two_d_elements(SIZE(T_nogho,1),3))
    ALLOCATE(error_L2(SIZE(T_nogho,1)))
    ALLOCATE(error_L2_init(SIZE(T_nogho,1)))

#endif

    two_d_nodes = 0.
    two_d_elements = 0.
    h = 0.
    error_L2_vertices = 0.
    h_target_temp = 0.
    error_target = adapt%tol_est
    h_target = 0.
    error_L2 = 0.
    error_L2_init = 0.
    eg_L2 = 0.
    eg_L2_init = 0.

    two_d_nodes = Mesh%X(vector_nodes_unique,1:2)
#ifndef PARALL
    two_d_elements = Mesh%T(:,1:3)
#else
    two_d_elements = T_nogho(:,1:3)
#endif
    !! use error map to create element size map: h_target
    CALL h_map(N_n_vertex,two_d_nodes,two_d_elements,vector_nodes_unique, h)

    h_target_temp = 0.
    h_target = 100.

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

    ! create contribution of estimator and indicator for mmg
    ! if oscillations are detected just use original size/2 otherwise use
    ! estimator calculated with Richardson

    DO i=1,SIZE(h_target) !check on coarsening ---> not higher than initial mesh
       IF (h_target(i) .GT. 0.1) THEN
          h_target(i) = 0.1
       ENDIF
    ENDDO

#ifndef PARALL
    CALL generate_htarget_sol_file(N_n_vertex, h_target)
#else
    ! reduce the size of global h_target array
    CALL MPI_Reduce(N_n_vertex, N_n_vertex_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD, ierr);

    ALLOCATE(h_target_glob_rep(N_n_vertex_glob))
    ALLOCATE(vector_nodes_unique_glob(N_n_vertex_glob))
    h_target_glob_rep = 0
    vector_nodes_unique_glob = 0

    ! vector containing N_n_vertex from each process
    ALLOCATE(recvcounts(MPIvar%glob_size))
    recvcounts(MPIvar%glob_id+1) = N_n_vertex
    CALL MPI_Allgather(N_n_vertex, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! vector containing the displacemenet
    ALLOCATE(displs(MPIvar%glob_size))
    displs(1) = 0
    DO i = 2, MPIvar%glob_size
       displs(i) = displs(i-1) + recvcounts(i-1)
    ENDDO


    CALL MPI_Gatherv(h_target, N_n_vertex, MPI_REAL8, h_target_glob_rep, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_Gatherv(Mesh%loc2glob_nodes(vector_nodes_unique), N_n_vertex, MPI_INT, vector_nodes_unique_glob, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    DEALLOCATE(displs, recvcounts)

    IF(MPIvar%glob_id .EQ. 0) THEN


       ALLOCATE(temp(SIZE(vector_nodes_unique_glob)))
       temp = vector_nodes_unique_glob

       CALL quicksort_int(temp)

       ! Count total duplicates
       counter = 0
       DO i = 2, SIZE(temp)
          IF (temp(i) .EQ. temp(i - 1)) THEN
             counter = counter + 1
          ENDIF
       END DO

       DEALLOCATE(temp)
       ALLOCATE(already_counted(N_n_vertex_glob))
       ALLOCATE(h_target_glob(N_n_vertex_glob-counter))
       already_counted = .FALSE.
       h_target_glob = 0

       DO i = 1, N_n_vertex_glob
          dummy_sum = 0
          counter = 0
          IF(already_counted(i) .EQV. .FALSE.) THEN
             DO j = 1, N_n_vertex_glob
                IF(vector_nodes_unique_glob(i) .EQ. vector_nodes_unique_glob(j)) THEN
                   dummy_sum = dummy_sum + h_target_glob_rep(j)
                   already_counted(j) = .TRUE.
                   counter = counter + 1
                ENDIF
             ENDDO
             h_target_glob(vector_nodes_unique_glob(i)) = dummy_sum/counter
          ENDIF
       ENDDO

       DEALLOCATE(already_counted)

       CALL generate_htarget_sol_file(SIZE(h_target_glob), h_target_glob)

    ENDIF
    DEALLOCATE(h_target_glob_rep)
    DEALLOCATE(vector_nodes_unique_glob)
#endif

#ifdef PARALL
    IF(MPIvar%glob_id .EQ. 0) THEN
#endif

       CALL extract_mesh_name_from_fullpath_woext(mesh_name, mesh_name_npne)

       WRITE(param_adapt_char, *) param_adapt
       WRITE(count_adapt_char, *) count_adapt
       new_mesh_name_npne = TRIM(ADJUSTL(mesh_name_npne)) // '_param'// TRIM(ADJUSTL(param_adapt_char)) // '_n' // TRIM(ADJUSTL(count_adapt_char))

       buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".mesh"
       IF(MPIvar%glob_id .EQ. 0) THEN
          CALL mmg_create_mesh_from_h_target(buffer)
       ENDIF

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
       
       

#ifdef PARALL
    ENDIF
    ! wait for process 0 to finish writing before loading new meshe
    CALL MPI_BARRIER(mpi_comm_world, ierr)
#endif

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) "********** Loading mesh P1  **********"
    ENDIF
    CALL free_mesh
    IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
       CALL load_gmsh_mesh("./res/temp",0)
    ELSE
       CALL load_gmsh_mesh("./res/temp",1)
    ENDIF

    CALL create_reference_element(refElPol,2,1, verbose = 0)
    CALL mesh_preprocess_serial(ierr)

    Mesh%X = Mesh%X*phys%lscale

    IF(ierr .EQ. 0) THEN
       WRITE(*,*) "Error! Corresponding face in Tb not found. STOP"
       STOP
    ENDIF

    CALL read_extended_connectivity('./res/temp.msh')

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) "********** Increasing order mesh **********"
    ENDIF

    CALL set_order_mesh(order)
    CALL free_reference_element
    CALL create_reference_element(refElPol,2,order, verbose = 0)
    CALL mesh_preprocess_serial(ierr)

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

    DEALLOCATE(two_d_nodes)
    DEALLOCATE(two_d_elements)
    DEALLOCATE(h)
    DEALLOCATE(error_L2_vertices)
    DEALLOCATE(error_target)
    DEALLOCATE(h_target)
    DEALLOCATE(h_target_temp)
    DEALLOCATE(error_L2)
    DEALLOCATE(error_L2_init)
    DEALLOCATE(vector_nodes_unique)

  END SUBROUTINE adaptivity_estimator

  SUBROUTINE adaptivity_estimator_get_error(param_adapt,vector_nodes_unique,error_estimator)

    INTEGER, INTENT(IN)                         :: param_adapt
    INTEGER, INTENT(IN)                         :: vector_nodes_unique(:)
    REAL*8, INTENT(OUT)                         :: error_estimator(:)

    REAL*8,  ALLOCATABLE                        :: error_L2(:), error_L2_vertices(:), error_L2_init(:), error_target(:)
    INTEGER                                     :: N_n_vertex
    REAL*8                                      :: eg_L2

    N_n_vertex = SIZE(error_estimator)
    ALLOCATE(error_L2_vertices(N_n_vertex))
    ALLOCATE(error_target(N_n_vertex))
    ALLOCATE(error_L2(SIZE(Mesh%T,1)))
    ALLOCATE(error_L2_init(SIZE(Mesh%T,1)))

    ! error estimation for the mesh and the solution
    CALL L2_error_estimator_eval(Mesh%X,Mesh%T,sol%u,sol%q,param_adapt,error_L2,eg_L2)
    CALL error_on_vertices(error_L2,Mesh%T,vector_nodes_unique, N_n_vertex, error_L2_vertices)

    error_estimator = error_L2_vertices

    DEALLOCATE(error_L2_vertices)
    DEALLOCATE(error_target)
    DEALLOCATE(error_L2)
    DEALLOCATE(error_L2_init)

  END SUBROUTINE adaptivity_estimator_get_error

  SUBROUTINE calculate_L2_error_two_sols_different_p_scalar_general(X,T,error_param, u_sol,u_star_sol, error_L2, eg_L2)

    REAL*8, INTENT(IN)                :: X(:,:)
    INTEGER, INTENT(IN)               :: T(:,:)
    REAL*8, INTENT(IN)                :: u_sol(:,:), u_star_sol(:,:)
    INTEGER, INTENT(IN)               :: error_param
    REAL*8, INTENT(OUT)               :: error_L2(:)
    REAL*8, INTENT(OUT)               :: eg_L2
    TYPE(Reference_element_type)      :: refElv_star

    CALL create_reference_element(refElv_star,2,refElPol%Ndeg+1, verbose = 0)
    ! elemental error,global error, area elements
    CALL calculate_L2_error_two_sols_different_p_scalar(X,T,u_sol(:,error_param),refElPol,X,T,u_star_sol(:,error_param),refElv_star,error_L2,eg_L2)
    CALL free_reference_element_pol(refElv_star)

  ENDSUBROUTINE calculate_L2_error_two_sols_different_p_scalar_general

  SUBROUTINE post_process_matrix_solution(X,T,u,q,u_sol,u_star_sol)
    USE physics, ONLY: cons2phys

    REAL*8, INTENT(IN)                :: X(:,:)
    INTEGER, INTENT(IN)               :: T(:,:)
    REAL*8, INTENT(IN)                :: u(:), q(:)
    REAL*8, ALLOCATABLE, INTENT(OUT)  :: u_sol(:,:), u_star_sol(:,:)

    TYPE(Reference_element_type)      :: refElv_star
    REAL*8, ALLOCATABLE               :: K(:,:,:), Bt(:,:,:), int_N(:,:,:)
    REAL*8, ALLOCATABLE               :: u_star(:), u_int(:)
    INTEGER                           :: n_elements

#ifndef PARALL
    n_elements = Mesh%Nelems
#else
    n_elements = SIZE(Mesh%T,1) - Mesh%nghostelems
#endif

    ALLOCATE(K((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq,(refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq, n_elements))
    ALLOCATE(Bt((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq,phys%neq*mesh%ndim*(refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2, n_elements))
    ALLOCATE(int_N((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq,phys%neq, n_elements))
    ALLOCATE(u_star((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq*n_elements))
    ALLOCATE(u_int((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq*n_elements))
    ALLOCATE(u_sol(Mesh%Nnodesperelem*n_elements*phys%neq/phys%Neq, phys%npv))
    ALLOCATE(u_star_sol((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq*n_elements/phys%Neq, phys%npv))


    CALL create_reference_element(refElv_star,2,refElPol%Ndeg+1, verbose = 0)
    ! local calculation of the p+1 solution u_star
    CALL hdg_post_process_matrix(X,T,refElv_star,refElPol%Ndeg,refElPol%Ndeg+1, K, Bt, int_N)
    CALL hdg_postprocess_solution(q,u,K,Bt,int_N,refElv_star,refElPol,n_elements, u_star, u_int)

    CALL cons2phys(TRANSPOSE(RESHAPE(u, (/phys%Neq, SIZE(T,1)*SIZE(T,2)/))), u_sol)
    CALL cons2phys(TRANSPOSE(RESHAPE(u_star, (/phys%Neq, SIZE(T,1)*(refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2/))), u_star_sol)

    ! elemental error,global error, area elements
    CALL free_reference_element_pol(refElv_star)

    DEALLOCATE(K)
    DEALLOCATE(Bt)
    DEALLOCATE(int_N)
    DEALLOCATE(u_star)
    DEALLOCATE(u_int)


  ENDSUBROUTINE post_process_matrix_solution

  SUBROUTINE L2_error_estimator_eval(X,T,u,q,error_param,error_L2,eg_L2)
    USE physics, ONLY: cons2phys

    REAL*8, INTENT(IN)                :: X(:,:)
    INTEGER, INTENT(IN)               :: T(:,:)
    REAL*8, INTENT(IN)                :: u(:), q(:)
    INTEGER, INTENT(IN)               :: error_param
    REAL*8, INTENT(OUT)               :: error_L2(:)
    REAL*8, INTENT(OUT)               :: eg_L2

    TYPE(Reference_element_type)      :: refElv_star
    REAL*8, ALLOCATABLE               :: K(:,:,:), Bt(:,:,:), int_N(:,:,:)
    REAL*8, ALLOCATABLE               :: u_star(:), u_int(:), u_sol(:,:), u_star_sol(:,:)
    INTEGER                           :: n_elements

#ifndef PARALL
    n_elements = Mesh%Nelems
#else
    n_elements = SIZE(Mesh%T,1) - Mesh%nghostelems
#endif

    ALLOCATE(K((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq,(refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq, n_elements))
    ALLOCATE(Bt((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq,phys%neq*mesh%ndim*(refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2, n_elements))
    ALLOCATE(int_N((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq,phys%neq, n_elements))
    ALLOCATE(u_star((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq*n_elements))
    ALLOCATE(u_int((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq*n_elements))
    ALLOCATE(u_sol(Mesh%Nnodesperelem*n_elements*phys%neq/phys%Neq, phys%npv))
    ALLOCATE(u_star_sol((refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2*phys%neq*n_elements/phys%Neq, phys%npv))


    CALL create_reference_element(refElv_star,2,refElPol%Ndeg+1, verbose = 0)

    ! local calculation of the p+1 solution u_star
    CALL hdg_post_process_matrix(X,T,refElv_star,refElPol%Ndeg,refElPol%Ndeg+1, K, Bt, int_N)
    CALL hdg_postprocess_solution(q,u,K,Bt,int_N,refElv_star,refElPol,n_elements, u_star, u_int)

    CALL cons2phys(TRANSPOSE(RESHAPE(u, (/phys%Neq, SIZE(T,1)*SIZE(T,2)/))), u_sol)
    CALL cons2phys(TRANSPOSE(RESHAPE(u_star, (/phys%Neq, SIZE(T,1)*(refElPol%Ndeg+2)*(refElPol%Ndeg+3)/2/))), u_star_sol)

    ! elemental error,global error, area elements
    CALL calculate_L2_error_two_sols_different_p_scalar(X,T,u_sol(:,error_param),refElPol,X,T,u_star_sol(:,error_param),refElv_star,error_L2,eg_L2)
    CALL free_reference_element_pol(refElv_star)

    DEALLOCATE(K)
    DEALLOCATE(Bt)
    DEALLOCATE(int_N)
    DEALLOCATE(u_sol)
    DEALLOCATE(u_star_sol)
    DEALLOCATE(u_star)
    DEALLOCATE(u_int)


  ENDSUBROUTINE L2_error_estimator_eval



  SUBROUTINE hdg_post_process_matrix(X, T, refElv, p1, p2, K, Bt, int_N, M)
    TYPE(Reference_element_type), INTENT(IN)                :: refElv
    REAL*8, INTENT(IN)              :: X(:,:)
    INTEGER, INTENT(IN)             :: T(:,:)
    INTEGER, INTENT(IN)             :: p1, p2
    REAL*8, INTENT(OUT)             :: K(:,:,:), Bt(:,:,:), int_N(:,:,:)
    REAL*8, INTENT(OUT), OPTIONAL   :: M(:,:,:)

    ! Declarations
    INTEGER                              :: n_elements, Nv2, iElem
    REAL*8, DIMENSION(:, :), ALLOCATABLE :: shapeFunctions
    REAL*8, DIMENSION(:, :), ALLOCATABLE :: Ke, Bte, int_Ne, Me

    ! Initialization
    n_elements = SIZE(T, 1)

    Nv2 = 0
    IF(refElv%elemType .EQ. 0) THEN
       ! the casting is only necessary for not showing the compiler warning
       Nv2 = INT(0.5*(p2+1)*(p2+2))
    ELSEIF(refElv%elemType .EQ. 1) THEN
       Nv2 = p2*(p2+1)
    ELSE
       WRITE(*,*) "not coded.STOP."
       STOP
    ENDIF

    ALLOCATE(shapeFunctions(Nv2, SIZE(T, 2)))
    ALLOCATE(Ke(SIZE(K,1),SIZE(K,2)))
    ALLOCATE(Bte(SIZE(Bt,1),SIZE(Bt,2)))
    ALLOCATE(int_Ne(SIZE(int_N,1), SIZE(int_N,2)))
    IF(PRESENT(M)) THEN
       ALLOCATE(Me(SIZE(M,1), SIZE(M,2)))
    ENDIF

    shapeFunctions = 0.

    CALL compute_shape_functions_at_interp_points(p1,p2,refElv%elemType,shapeFunctions)

    ! Loop in elements
    IF(PRESENT(M)) THEN
       DO iElem = 1, n_elements

          CALL elemental_matrices(MATMUL(shapeFunctions, X(T(iElem, :), :)), refElv, Ke, Bte, int_Ne, Me)
          K(:, :, iElem) = Ke
          Bt(:, :, iElem) = Bte
          int_N(:, :, iElem) = int_Ne
          M(:, :, iElem) = Me
       ENDDO
    ELSE
       DO iElem = 1, n_elements

          CALL elemental_matrices(MATMUL(shapeFunctions, X(T(iElem, :), :)), refElv, Ke, Bte, int_Ne)
          K(:, :, iElem) = Ke
          Bt(:, :, iElem) = Bte
          int_N(:, :, iElem) = int_Ne
       ENDDO
    ENDIF


    DEALLOCATE(shapeFunctions)
    DEALLOCATE(Ke)
    DEALLOCATE(Bte)
    DEALLOCATE(int_Ne)
    IF(PRESENT(M)) THEN
       DEALLOCATE(Me)
    ENDIF


  ENDSUBROUTINE hdg_post_process_matrix


  SUBROUTINE elemental_matrices(Xe, refElv, Ke, Bte, int_Ne, Me)
    TYPE(Reference_element_type), INTENT(IN)               :: refElv
    REAL*8, DIMENSION(:,:), INTENT(IN)                     :: Xe
    REAL*8, DIMENSION(:, :), INTENT(INOUT)                 :: Ke, Bte, int_Ne
    REAL*8, DIMENSION(:, :), OPTIONAL, INTENT(INOUT)       :: Me
    REAL*8, ALLOCATABLE,DIMENSION(:,:)                     :: K_loc, Bx, By
    REAL*8, ALLOCATABLE,DIMENSION(:)                       :: int_N_loc, Nx_g, Ny_g
    REAL*8, DIMENSION(2, 2)                                :: J, invJ
    REAL*8                                                 :: xg, dvolu, detJ
    INTEGER                                                :: n_gauss, g


    ALLOCATE(K_loc(refElv%Nnodes2D,refElv%Nnodes2D))
    ALLOCATE(Bx(refElv%Nnodes2D,refElv%Nnodes2D))
    ALLOCATE(By(refElv%Nnodes2D,refElv%Nnodes2D))
    ALLOCATE(int_N_loc(refElv%Nnodes2D))
    ALLOCATE(Nx_g(refElv%Nnodes2D))
    ALLOCATE(Ny_g(refElv%Nnodes2D))

    Ke = 0.
    Bte = 0.
    int_Ne = 0.
    IF(PRESENT(Me)) THEN
       Me = 0.
    ENDIF
    K_loc = 0.
    Bx = 0.
    By = 0.
    int_N_loc = 0.
    Nx_g = 0.
    Ny_g = 0.
    xg = 0.
    dvolu = 0.
    detJ = 0.

    !Number of Gauss points in the interior
    n_gauss = refElv%NGauss2D

    ! VOLUME COMPUTATIONS: LOOP IN GAUSS POINTS
    DO g = 1, n_gauss
       !Shape functions and derivatives at the current integration point

       ! Jacobian
       J(1, 1) = dot_PRODUCT(refElv%Nxi2D(g, :), Xe(:, 1))
       J(1, 2) = dot_PRODUCT(refElv%Nxi2D(g, :), Xe(:, 2))
       J(2, 1) = dot_PRODUCT(refElv%Neta2D(g, :), Xe(:, 1))
       J(2, 2) = dot_PRODUCT(refElv%Neta2D(g, :), Xe(:, 2))

       detJ = J(1, 1)*J(2, 2) - J(1, 2)*J(2, 1)

       ! Gauss point position
       xg = dot_PRODUCT(refElv%N2d(g, :), Xe(:, 1))

       ! Integration weight
       dvolu = refElv%gauss_weights2D(g) * detJ

       IF (switch%axisym) THEN
          dvolu = dvolu*xg
       END IF

       ! x and y derivatives
       CALL invert_matrix(J, invJ)

       Nx_g = invJ(1, 1) * refElv%Nxi2D(g, :) + invJ(1, 2) * refElv%Neta2D(g, :)
       Ny_g = invJ(2, 1) * refElv%Nxi2D(g, :) + invJ(2, 2) * refElv%Neta2D(g, :)

       ! Contribution of the current integration point to the elemental matrix
       K_loc = K_loc + (TensorProduct(Nx_g, Nx_g) + TensorProduct(Ny_g,Ny_g)) * dvolu
       IF(PRESENT(Me)) THEN
          Me = Me + TensorProduct(refElv%N2D(g, :), refElv%N2D(g, :)) * dvolu
       ENDIF
       Bx = Bx + TensorProduct(Nx_g, refElv%N2D(g, :)) * dvolu
       By = By + TensorProduct(Ny_g, refElv%N2D(g, :)) * dvolu
       int_N_loc = int_N_loc + refElv%N2D(g, :)*dvolu
    ENDDO

    ! Expand the matrices
    CALL expand_matrix_A(K_loc, phys%neq, Ke)
    CALL expand_matrix_Bt(Bx, By, Bte)
    CALL expand_matrix_A(RESHAPE(int_N_loc, [SIZE(int_N_loc), 1]), phys%neq, int_Ne)

    DEALLOCATE(K_loc)
    DEALLOCATE(Bx)
    DEALLOCATE(By)
    DEALLOCATE(int_N_loc)
    DEALLOCATE(Nx_g)
    DEALLOCATE(Ny_g)

  CONTAINS

    SUBROUTINE expand_matrix_A(A, n, res)
      REAL*8, INTENT(IN)                        :: A(:,:)
      INTEGER, INTENT(IN)                       :: n
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)   :: temp_1, temp_2
      REAL*8, INTENT(OUT)                       :: res(:,:)
      INTEGER                                   :: i, j, k, l


      ALLOCATE(temp_1(SIZE(A,1), SIZE(A,2),n ,n))
      ALLOCATE(temp_2(n, SIZE(A,1), n, SIZE(A,2)))


      temp_1 = 0.
      temp_2 = 0.
      res = 0.

      ! equivalent to matlab repmat res(:,:,1:n+1:n^2) = repmat(A, [1 1 n]);
      DO i = 1, n
         temp_1(:,:,i,i) = A
      ENDDO

      ! [1 2 3 4] --> [3 1 4 2]
      DO i = 1, SIZE(temp_1, 1)
         DO j = 1, SIZE(temp_1, 2)
            DO k = 1, SIZE(temp_1, 3)
               DO l = 1, SIZE(temp_1, 4)
                  temp_2(k, i, l, j) = temp_1(i, j, k, l)
               END DO
            END DO
         END DO
      END DO


      ! reshape
      res = RESHAPE(temp_2, [n*SIZE(A,1), n*SIZE(A,2)])

      DEALLOCATE(temp_1, temp_2)
    ENDSUBROUTINE expand_matrix_A

    SUBROUTINE expand_matrix_Bt(Bx, By, B)
      REAL*8, INTENT(in)    :: Bx(:, :), By(:, :)
      !real*8, INTENT(out)   :: B(phys%Neq*size(Bx, 1), phys%Neq*Mesh%ndim*size(Bx, 2))
      REAL*8, INTENT(out)   :: B(:,:)
      INTEGER*4             :: i, k, j, n, m, neq_dim, neq, ndim


      ! neq_dim = phys%Neq*Mesh%ndim
      ! neq = phys%Neq
      ! ndim = Mesh%ndim

      neq_dim = 2*Mesh%ndim
      neq = 2
      ndim = Mesh%ndim

      n = SIZE(Bx, 1)
      m = SIZE(Bx, 2)
      B = 0.d0

      DO j = 1, m
         DO i = 1, n
            DO k = 1, neq
               B((i - 1)*neq + k, (j - 1)*neq_dim + 1 + (k - 1)*ndim) = Bx(i, j)
               B((i - 1)*neq + k, (j - 1)*neq_dim + 2 + (k - 1)*ndim) = By(i, j)
            END DO
         END DO
      END DO
    ENDSUBROUTINE expand_matrix_Bt


  ENDSUBROUTINE elemental_matrices


  SUBROUTINE hdg_postprocess_solution(qin, uin, K, Bt, int_N, refElv_star, refElv, n_elements, u_star,u_int)
    TYPE(Reference_element_type), INTENT(IN)                 :: refElv
    TYPE(Reference_element_type), INTENT(IN)                 :: refElv_star
    REAL*8, DIMENSION(:), INTENT(IN)                         :: qin, uin
    REAL*8, DIMENSION(:, :, :), INTENT(IN)                   :: K, Bt, int_N
    INTEGER, INTENT(IN)                                      :: n_elements
    REAL*8, DIMENSION(:), INTENT(OUT)                        :: u_star
    REAL*8, DIMENSION(:), INTENT(OUT)                        :: u_int

    REAL*8, DIMENSION(:, :), ALLOCATABLE                     :: Ke, Bte, int_Ne, int_ue_star
    REAL*8, DIMENSION(:, :), ALLOCATABLE                     :: u, q, u_star2D, q_star2D
    REAL*8, DIMENSION(:, :), ALLOCATABLE                     :: shapeFunctions, A
    REAL*8, DIMENSION(:), ALLOCATABLE                        :: BtL, L_int, int_ue, b, sol
    INTEGER, DIMENSION(:), ALLOCATABLE                       :: ind_u_star, ind_L_star
    INTEGER                                                  :: Nv, Nv_star, ielem, i, j, counter
    INTEGER                                                  :: NeqNv, NeqNvStar, NeqNvStar2, NeqNvStar4, NeqNvStar2Elem, NeqNvStar4Elem

    INTEGER                                                  :: SizeK1, SizeK2, SizeBt1, SizeBt2, SizeIntN1, SizeIntN2, SizeIntUeStar1, SizeIntUeStar2, SizeA


    SizeK1 = SIZE(K,1)
    SizeK2 = SIZE(K,2)
    SizeBt1 = SIZE(Bt,1)
    SizeBt2 = SIZE(Bt,2)
    SizeIntN1 = SIZE(int_N,1)
    SizeIntN2 = SIZE(int_N,2)

    Nv_star = refElv_star%Nnodes2D
    Nv = refElv%Nnodes2D
    NeqNv = phys%Neq * Nv
    NeqNvStar = phys%Neq * Nv_star
    NeqNvStar2 = phys%neq * NeqNvStar
    NeqNvStar4 = phys%neq*Mesh%ndim * NeqNvStar
    NeqNvStar2Elem = phys%neq * Nv_star
    NeqNvStar4Elem = phys%neq*Mesh%ndim * Nv_star


    ! Reshape matrices (this might not be needed)
    ALLOCATE(u(Nv, phys%Neq*n_elements))
    ALLOCATE(q(Nv, refElPol%ndim*phys%Neq*n_elements))
    ALLOCATE(u_star2D(Nv_star, phys%Neq*n_elements))
    ALLOCATE(q_star2D(Nv_star, refElPol%ndim*phys%Neq*n_elements))

    ALLOCATE(Ke(SizeK1,SizeK2))
    ALLOCATE(Bte(SizeBt1,SizeBt2))
    ALLOCATE(int_Ne(SizeIntN1,SizeIntN2))
    ALLOCATE(int_ue_star(SizeIntN1,SizeIntN2))

    SizeIntUeStar1 = SizeIntN1
    SizeIntUeStar2 = SizeIntN2

    ! reshape uin to u(Nv, 2*n_elements)
    counter = 1
    DO j = 1, n_elements*Nv
       DO i = 1, phys%Neq
          u(MOD(j - 1, Nv) + 1,(i - 1) * n_elements + (j - 1) / Nv + 1) = uin(counter)
          counter = counter + 1
       ENDDO
    ENDDO

    ! reshape qin to q(Nv, 4*n_elements)
    counter = 1
    DO j = 1, n_elements*Nv
       DO i = 1, phys%Neq*refElPol%ndim
          q(MOD(j - 1, Nv) + 1,(i - 1) * n_elements + (j - 1) / Nv + 1) = qin(counter)
          counter = counter + 1
       ENDDO
    ENDDO


    ALLOCATE(shapeFunctions(Nv_star,Nv))

    CALL compute_shape_functions_at_interp_points(refElv%nDeg, refElv_star%nDeg, refElv%elemType,shapeFunctions)

    ! Interpolate solution at finer mesh
    ALLOCATE(L_int(refElPol%ndim*phys%Neq*n_elements*Nv_star))

    u_star2D = MATMUL(shapeFunctions, u)
    q_star2D = MATMUL(shapeFunctions, q)

    counter = 1
    DO j = 1, n_elements*Nv_star
       DO i = 1, phys%Neq
          u_int(counter) = u_star2D(MOD(j - 1, Nv_star) + 1,(i - 1) * n_elements + (j - 1) / Nv_star + 1)
          counter = counter + 1
       ENDDO
    ENDDO

    ! reshape qin to q(Nv, 4*n_elements)
    counter = 1
    DO j = 1, n_elements*Nv_star
       DO i = 1, phys%Neq*refElPol%ndim
          L_int(counter) = q_star2D(MOD(j - 1, Nv_star) + 1,(i - 1) * n_elements + (j - 1) / Nv_star + 1)
          counter = counter + 1
       ENDDO
    ENDDO

    SizeA = NeqNvStar + phys%neq

    ALLOCATE(A(SizeA, SizeA))
    ALLOCATE(b(SizeA))
    ALLOCATE(sol(SizeA))
    ALLOCATE(ind_u_star(NeqNvStar))
    ALLOCATE(ind_L_star(phys%Neq*refElPol%ndim*Nv_star))
    ALLOCATE(int_ue(phys%Neq))
    ALLOCATE(BtL(NeqNvStar))

    DO ielem = 1, n_elements
       ! Index
       ind_u_star = (ielem - 1) * NeqNvStar2Elem + (/(j, j=1,NeqNvStar2Elem)/)
       ind_L_star = (ielem - 1) * NeqNvStar4Elem + (/(j, j=1,NeqNvStar4Elem)/)

       ! Elemental matrices
       Ke = K(:, :, ielem)
       Bte = Bt(:, :, ielem)
       int_Ne = int_N(:, :, ielem)

       ! Multiplication
       BtL = MATMUL(Bte, L_int(ind_L_star))
       int_ue_star = int_Ne
       int_ue = MATMUL(TRANSPOSE(int_Ne), u_int(ind_u_star))

       ! Lagrange multipliers
       A(1:SizeK1,1:SizeK2) = Ke
       A(SizeK1+1:SizeA,1:SizeIntUeStar1) = TRANSPOSE(int_ue_star)
       A(1:SizeIntUeStar1,SizeK1+1:SizeA) = int_ue_star
       A(SizeIntUeStar1+1:SizeA,SizeIntUeStar1+1:SizeA) = 0.

       b(1:NeqNvStar) = BtL
       b(NeqNvStar+1:SizeA) = int_ue

       ! Elemental solution
       CALL solve_linear_system_sing(A,b,sol)
       ! Postprocessed solution
       u_star(ind_u_star) = sol(1 : SIZE(sol) - phys%neq)
    ENDDO

    DEALLOCATE(u)
    DEALLOCATE(q)
    DEALLOCATE(u_star2D)
    DEALLOCATE(q_star2D)
    DEALLOCATE(A)
    DEALLOCATE(b)
    DEALLOCATE(sol)
    DEALLOCATE(L_int)
    DEALLOCATE(BtL)
    DEALLOCATE(int_ue)
    DEALLOCATE(ind_u_star)
    DEALLOCATE(ind_L_star)
    DEALLOCATE(shapeFunctions)
    DEALLOCATE(Ke)
    DEALLOCATE(Bte)
    DEALLOCATE(int_Ne)
    DEALLOCATE(int_ue_star)

  ENDSUBROUTINE hdg_postprocess_solution

  ! To simplify implementation we assume the order of refEl2 > order of refEl1
  SUBROUTINE calculate_L2_error_two_sols_different_p_scalar(X1, T1, u1, refEl1, X2, T2, u2, refEl2, error, eg)
    TYPE(Reference_element_type), INTENT(IN)    :: refEl1, refEl2
    REAL*8, INTENT(IN)                          :: X1(:,:), X2(:,:)
    INTEGER, INTENT(IN)                         :: T1(:,:), T2(:,:)
    REAL*8, INTENT(IN)                          :: u1(:), u2(:)
    REAL*8, INTENT(OUT)                         :: error(:)
    REAL*8, INTENT(OUT)                         :: eg

    REAL*8                                      :: shapeFunctions_post(refEl2%Nnodes2D,refEl1%Nnodes2D)
    REAL*8                                      :: shapeFunctions_temp(refEl1%Nnodes2D,refEl2%NGauss2D,3), shapeFunctions(refEl2%NGauss2D,refEl1%Nnodes2D,3)
    REAL*8                                      :: error2(SIZE(T2,1)), sol_norm(SIZE(T2,1)), dom_area(SIZE(T2,1))
    REAL*8                                      :: Xe_p1(refEl1%Nnodes2D, SIZE(X1,2)), Xe_p2(refEl2%Nnodes2D, SIZE(X2,2))
    REAL*8                                      :: ue_p1(refEl1%Nnodes2D), ue_p2(refEl2%Nnodes2D)
    INTEGER                                     :: ind_p1(refEl1%Nnodes2D), ind_p2(refEl2%Nnodes2D)
    INTEGER                                     :: nnodes_p1, nnodes_p2, n_elements, iElem, i
    REAL*8                                      :: elem_area, elem_sol_norm, elem_error

    CALL compute_shape_functions_at_interp_points(refEl1%nDeg,refEl2%nDeg,refEl1%elemType,shapeFunctions_post)

    n_elements = SIZE(T2,1)

    nnodes_p1 = refEl1%Nnodes2D
    nnodes_p2 = refEl2%Nnodes2D

    CALL compute_shape_functions_at_points(refEl1, refEl2%gauss_points2D, shapeFunctions_temp)

    shapeFunctions(:,:,1) = TRANSPOSE(shapeFunctions_temp(:,:,1))
    shapeFunctions(:,:,2) = TRANSPOSE(shapeFunctions_temp(:,:,2))
    shapeFunctions(:,:,3) = TRANSPOSE(shapeFunctions_temp(:,:,3))

    DO iElem = 1, n_elements

       Xe_p2 = MATMUL(shapeFunctions_post,X2(T2(ielem,:),:))
       Xe_p1 = X1(T1(ielem,:),:)

       ind_p1 =  (/(i, i=(ielem-1)*nnodes_p1+1,ielem*nnodes_p1)/)
       ind_p2 =  (/(i, i=(ielem-1)*nnodes_p2+1,ielem*nnodes_p2)/)

       ue_p1 = u1(ind_p1)
       ue_p2 = u2(ind_p2)

       CALL compute_elemental_error(ue_p2,ue_p1,refEl2,refEl1,Xe_p2,Xe_p1,shapeFunctions,refEl2%gauss_points2D, elem_error, elem_sol_norm, elem_area)

       dom_area(ielem) = elem_area
       sol_norm(ielem) = elem_sol_norm
       error2(ielem) = elem_error

    ENDDO

    IF(adapt%difference .EQ. 0) THEN 
      !DO i = 1, size(sol_norm)
      !   IF (sol_norm(i)<1.e-3) THEN
      !      sol_norm(i)=sol_norm(i)+1.e-3
      !   ENDIF
      !ENDDO
       ! relative difference
      
       error = SQRT(error2/sol_norm)

    ELSEIF(adapt%difference .EQ. 1) THEN
       ! absolute difference
       error = SQRT(error2/SUM(sol_norm)*SUM(dom_area)/dom_area)
    ELSE
       WRITE(*,*) "Choice of relative/absolute difference for the adaptivity not valid. STOP."
       STOP
    ENDIF

    eg = SQRT(SUM(error2)/SUM(dom_area))

  ENDSUBROUTINE calculate_L2_error_two_sols_different_p_scalar

  SUBROUTINE compute_elemental_error(ue_p2, ue_p1, refEl2, refEl1, Xe_p2, Xe_p1, n, coordGauss, elem_error, sol_norm, elem_area)
    TYPE(Reference_element_type), INTENT(IN)           :: refEl2, refEl1
    REAL*8, DIMENSION(:), INTENT(IN)                   :: ue_p2, ue_p1
    REAL*8, DIMENSION(:,:,:), INTENT(INOUT)            :: n
    REAL*8, DIMENSION(:, :), INTENT(IN)                :: Xe_p2, Xe_p1
    REAL*8, DIMENSION(:, :), INTENT(IN)                :: coordGauss
    REAL*8, INTENT(OUT)                                :: elem_error, sol_norm, elem_area

    REAL*8, DIMENSION(refEl2%NGauss2D,refEl2%Ndim)       :: xyg_xi, xyg_eta, xyg_p2, xyg_p1, xyg_err
    REAL*8, DIMENSION(refEl2%Ndim,refEl2%NGauss2D)       :: gxy0,gxy
    REAL*8, DIMENSION(refEl2%NGauss2D)                  :: detJ, dvolu, ueg_p1, ueg_p2, err_gauss
    LOGICAL                                            :: convergence

    xyg_xi = MATMUL(refEl2%Nxi2D, Xe_p2)
    xyg_eta = MATMUL(refEl2%Neta2D, Xe_p2)

    detJ = xyg_xi(:, 1) * xyg_eta(:, 2) - xyg_xi(:, 2) * xyg_eta(:, 1)
    dvolu = detJ * refEl2%gauss_weights2D

    ! check gauss point position
    xyg_p2 = MATMUL(refEl2%N2D, Xe_p2)
    xyg_p1 = MATMUL(n(:, :, 1), Xe_p1)
    xyg_err = (xyg_p2 - xyg_p1)

    IF (MAXVAL(ABS(xyg_err)) .GT. 1e-10) THEN
       gxy0 = TRANSPOSE(coordGauss)

       CALL recompute_gauss_points(Xe_p1,xyg_p2,refEl1,gxy0,n, gxy,convergence)

       IF (.NOT. convergence) THEN
          WRITE(*,*) "Recomputation of Gauss points failed (not converged). STOP."
          STOP
       ENDIF
    ENDIF

    ueg_p2 = MATMUL(refEl2%N2D,ue_p2)
    ueg_p1 = MATMUL(n(:,:,1),ue_p1)
    err_gauss = ueg_p2-ueg_p1
    elem_error = DOT_PRODUCT(err_gauss,dvolu*err_gauss)
    elem_area = SUM(dvolu)
    sol_norm = DOT_PRODUCT(ueg_p2,dvolu*ueg_p2)

  ENDSUBROUTINE compute_elemental_error

  SUBROUTINE recompute_gauss_points(Xe_M, xyg_g, refEl1, gxy0, n, gxy, convergence)
    TYPE(Reference_element_type), INTENT(IN)    :: refEl1
    REAL*8, DIMENSION(:, :), INTENT(IN)         :: Xe_M, xyg_g
    REAL*8, DIMENSION(:,:)                      :: gxy0
    REAL*8, DIMENSION(:,:), INTENT(OUT)         :: gxy
    REAL*8, DIMENSION(:,:,:), INTENT(INOUT)     :: n
    LOGICAL, INTENT(OUT)                        :: convergence
    ! Declarations for other variables used in the subroutine
    REAL*8                                      :: res_x(SIZE(gxy0,1),SIZE(gxy0,2)), res_fun(SIZE(xyg_g,1), SIZE(xyg_g,2)), ac(SIZE(Xe_M,1),SIZE(Xe_M,2)), bd(SIZE(Xe_M,1),SIZE(Xe_M,2))
    REAL*8                                      :: ac_perm(SIZE(ac,2),1,SIZE(ac,1)), bd_perm(SIZE(bd,2),1,SIZE(bd,1))
    REAL*8                                      :: Vand(SIZE(Xe_M,1),SIZE(Xe_M,1)), invVand(SIZE(Xe_M,1),SIZE(Xe_M,1))
    REAL*8                                      :: invDetJ(1,1,SIZE(n,1)), invJ(2,2,SIZE(n,1)), aux(1,2,SIZE(n,1))
    REAL*8                                      :: p(SIZE(Xe_M,1),SIZE(n,1)), dpxi(SIZE(Xe_M,1),SIZE(n,1)), dpeta(SIZE(Xe_M,1),SIZE(n,1)), shapeFunctions(3*SIZE(n,1),SIZE(Xe_M,1))
    INTEGER                                     :: n_gauss, max_iter, iter,i, j
    REAL*8                                      :: tol


    ! Determine n_gauss, max_iter, tol, etc.
    n_gauss =  SIZE(n,1)
    iter = 1
    convergence = .TRUE.
    res_x = 1.0d0
    res_fun = 1.0d0
    max_iter = 10
    tol = 1e-10
    gxy = gxy0

    CALL vandermonde_2d(Vand,refEl1)
    CALL invert_matrix(Vand,invVand)

    DO WHILE (MAXVAL(ABS(res_fun)) > tol .AND. MAXVAL(ABS(res_x)) > tol .AND. iter < max_iter)

       ac = MATMUL(n(:, :, 2), Xe_M)
       bd = MATMUL(n(:, :, 3), Xe_M)

       invDetJ(1,1,:) = 1.0d0 / (ac(:, 1) * bd(:, 2) - bd(:, 1) * ac(:,2))

       ! permute [1 2] --> [2 3 1]
       DO i = 1, SIZE(ac,1)
          DO j = 1, SIZE(ac,2)
             ac_perm(j, 1, i) = ac(i, j)
             bd_perm(j, 1, i) = bd(i, j)
          END DO
       END DO

       invJ(1, 1, :) = (bd_perm(2, 1, :) - bd_perm(1, 1, :)) * invDetJ(1, 1, :)
       invJ(2, 1, :) = (-ac_perm(2, 1, :) + ac_perm(1, 1, :)) * invDetJ(1, 1, :)

       aux(1,:,:) = TRANSPOSE(xyg_g - MATMUL(n(:, :, 1), Xe_M))
       gxy = gxy0 + RESHAPE(SUM(invJ*aux, 2), [2, n_gauss])

       CALL orthopoly2d_deriv(gxy(i,1), gxy(i,2), refEl1%nDeg, n_gauss, p, dpxi, dpeta)

       shapeFunctions = TRANSPOSE(MATMUL(invVand, RESHAPE([p,dpxi,dpeta],[3*SIZE(Xe_M,1), SIZE(n,1)])))

       n(:, :, 1) = shapeFunctions(:, :)
       n(:, :, 2) = shapeFunctions(n_gauss + 1:2 * n_gauss, :)
       n(:, :, 3) = shapeFunctions(2 * n_gauss + 1:, :)

       res_fun = xyg_g - MATMUL(n(:, :, 1), Xe_M)
       res_x = gxy - gxy0
       gxy0 = gxy
       iter = iter + 1
    END DO

    IF (iter == max_iter) convergence = .FALSE.
  END SUBROUTINE recompute_gauss_points

  SUBROUTINE error_on_vertices(error_L2,T, vector_nodes_unique,N_n_vertex,error_output)

    REAL*8, INTENT(IN)               :: error_L2(:)
    INTEGER, INTENT(IN)              :: vector_nodes_unique(:)
    INTEGER, INTENT(IN)              :: T(:,:)
    REAL*8, INTENT(OUT)              :: error_output(:)
    INTEGER, INTENT(IN)              :: N_n_vertex

    INTEGER                          :: count_vec(N_n_vertex)
    REAL*8                           :: g(N_n_vertex)
    INTEGER                          :: i, j, A, B, C


    g = 0.
    count_vec = 0

    DO i = 1, SIZE(T,1)

       ! Find indexes of the vertex nodes
       A = 0
       B = 0
       C = 0

       DO j = 1, N_n_vertex
          IF ((T(i, 1) + 1) .EQ. (vector_nodes_unique(j) + 1)) THEN
             A = j
          ENDIF
          IF ((T(i, 2) + 1) .EQ. (vector_nodes_unique(j) + 1)) THEN
             B = j
          ENDIF
          IF ((T(i, 3) + 1) .EQ. (vector_nodes_unique(j) + 1)) THEN
             C = j
          ENDIF

          IF(A*B*C .NE. 0) EXIT
       ENDDO

       IF((A .EQ. 0) .OR. (B .EQ. 0) .OR. (C .EQ. 0)) THEN
          WRITE(*,*) "Point not found! error_on_vertices. STOP."
          STOP
       ENDIF

       g(A) = g(A) + error_L2(i)
       count_vec(A) = count_vec(A) + 1

       g(B) = g(B) + error_L2(i)
       count_vec(B) = count_vec(B) + 1

       g(C) = g(C) + error_L2(i)
       count_vec(C) = count_vec(C) + 1
    ENDDO


    ! Calculate h values
    DO i = 1, N_n_vertex
       IF (count_vec(i) .GT. 0) THEN
          error_output(i) = g(i) / count_vec(i)
       ELSE
          WRITE(*,*) "GOT A DIVISION BY 0 IN error_on_vertices ADAPTIVITY, SOMETHING IS WRONG!"
          STOP
       ENDIF
    ENDDO


  ENDSUBROUTINE error_on_vertices

  SUBROUTINE error_on_vertices_target(error_L2,T,vector_nodes_unique,N_n_vertex,error_output)

    REAL*8, INTENT(IN)               :: error_L2(:)
    INTEGER, INTENT(IN)              :: vector_nodes_unique(:)
    INTEGER, INTENT(IN)              :: T(:,:)
    REAL*8, INTENT(OUT)              :: error_output(:)
    INTEGER, INTENT(IN)              :: N_n_vertex

    REAL*8                           :: count_vec(N_n_vertex)
    REAL*8                           :: g(N_n_vertex)
    INTEGER                          :: i, j, A, B, C, n_elements

    n_elements = SIZE(T,1)
    g = 0.
    count_vec = 0.

    DO i = 1, n_elements

       ! Find indices of the vertex nodes
       A = 0
       B = 0
       C = 0

       DO j = 1, N_n_vertex
          IF ((T(i, 1) + 1) .EQ. (vector_nodes_unique(j) + 1)) THEN
             A = j
          ENDIF
          IF ((T(i, 2) + 1) .EQ. (vector_nodes_unique(j) + 1)) THEN
             B = j
          ENDIF
          IF ((T(i, 3) + 1) .EQ. (vector_nodes_unique(j) + 1)) THEN
             C = j
          ENDIF

          IF(A*B*C .NE. 0) EXIT
       ENDDO

       IF((A .EQ. 0) .OR. (B .EQ. 0) .OR. (C .EQ. 0)) THEN
          WRITE(*,*) "Point not found! error_on_vertices. STOP."
          STOP
       ENDIF

       g(A) = g(A) + error_L2(i)
       count_vec(A) = count_vec(A) + 1.

       g(B) = g(B) + error_L2(i)
       count_vec(B) = count_vec(B) + 1.

       g(C) = g(C) + error_L2(i)
       count_vec(C) = count_vec(C) + 1.
    ENDDO

    ! Calculate h values
    DO i = 1, N_n_vertex
       IF (count_vec(i) .GT. 0.) THEN
          error_output(i) = g(i) / count_vec(i)
       ELSE
          WRITE(*,*) "GOT A DIVISION BY 0 IN error_on_vertices ADAPTIVITY, SOMETHING IS WRONG!"
          STOP
       ENDIF
    ENDDO


  ENDSUBROUTINE error_on_vertices_target

END MODULE adaptivity_estimator_module
