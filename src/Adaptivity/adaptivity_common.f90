!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! Module for schock capturing adaptivity
!************************************************************

MODULE adaptivity_common_module
  USE globals
  USE reference_element
  USE gmsh
  USE GMSH_io_module
  IMPLICIT NONE

CONTAINS

   SUBROUTINE adaptivity_console_output()
      CHARACTER(1024) :: buffer

      IF ((adapt%evaluator .EQ. 2)) THEN
         buffer = "            ADAPTIVITY ESTIMATOR                 "
      ELSEIF ((adapt%evaluator .EQ. 1)) THEN
         buffer = "            ADAPTIVITY INDICATOR                 "
      ELSEIF((adapt%evaluator .EQ. 0)) THEN
         buffer = "        ADAPTIVITY ESTIMATOR-INDICATOR           "
      ENDIF

      IF(MPIvar%glob_id .EQ. 0) THEN
         WRITE(*,*) "*************************************************"
         WRITE(*,*) TRIM(buffer)
         WRITE(*,*) "*************************************************"
      ENDIF

   ENDSUBROUTINE adaptivity_console_output

   SUBROUTINE calculate_h_map_elements(nodes,connectivity,h_map)
      REAL*8,INTENT(IN)                :: nodes(:,:)
      INTEGER, INTENT(IN)              :: connectivity(:,:)
      REAL*8, INTENT(OUT)              :: h_map(SIZE(connectivity,1))
      INTEGER                          :: i
      REAL*8, DIMENSION(2,2)           :: J
      REAL*8                           :: detJ

      DO i = 1, SIZE(connectivity,1)      
         CALL jacobian(nodes, connectivity(i,1), connectivity(i,2), connectivity(i,3), J)
         detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
         h_map(i) = SQRT(2.0*detJ/SQRT(3.0))
      ENDDO

   END SUBROUTINE calculate_h_map_elements       

   SUBROUTINE get_h_target_vertices(h_map_elements,h_target_vertices,T)
      REAL*8, INTENT(IN)                              :: h_map_elements(:)
      INTEGER,INTENT(IN)                              :: T(:,:)
      REAL*8, DIMENSION(:), POINTER, INTENT(INOUT)    :: h_target_vertices
      REAL*8, ALLOCATABLE                             :: h_target_nodal(:)
      INTEGER, ALLOCATABLE                            :: nodes_repeats(:)
      INTEGER                                         :: number_of_vertices


      ALLOCATE(h_target_nodal(SIZE(T,1)))
      ALLOCATE(nodes_repeats(SIZE(T,1)))      

      CALL sum_h_target_nodal(T, h_map_elements, h_target_nodal, nodes_repeats)

      number_of_vertices = COUNT(nodes_repeats /= 0)

      ALLOCATE(h_target_vertices(number_of_vertices))

      CALL average_h_target(h_target_nodal, nodes_repeats, h_target_vertices)

      DEALLOCATE(h_target_nodal, nodes_repeats)
   END SUBROUTINE get_h_target_vertices

   SUBROUTINE generate_new_mesh(mesh_name,h_target,count_adapt)
      USE in_out, ONLY: copy_file
      TYPE(gmsh_t)                :: gmsh
      CHARACTER(1024), INTENT(IN) :: mesh_name
      INTEGER, INTENT(IN)         :: count_adapt
      REAL*8, INTENT(IN)          :: h_target(:)
      INTEGER                     :: N_n_vertex
      CHARACTER(1024)             :: mesh_name_npne,new_mesh_name_npne, buffer
      CHARACTER(70)               :: param_adapt_char, count_adapt_char
      



      N_n_vertex = SIZE(h_target)

      CALL generate_htarget_sol_file(N_n_vertex,h_target)

      CALL extract_mesh_name_from_fullpath_woext(mesh_name, mesh_name_npne)

      WRITE(param_adapt_char, *) adapt%param_est
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

      buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".msh"
      CALL open_merge_with_geometry(gmsh, buffer)
      CALL copy_file(buffer, "./res/temp.msh")
      
      buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".sol"
      CALL delete_file(buffer)

   END SUBROUTINE generate_new_mesh

   SUBROUTINE load_new_mesh_gmsh(order)
      USE preprocess
      INTEGER, INTENT(IN) :: order
      INTEGER                     :: ierr

      IF(MPIvar%glob_id .EQ. 0) THEN
         WRITE(*,*) "********** Loading new mesh  **********"
      ENDIF

      CALL free_mesh

      CALL free_reference_element_pol(refElPol)
      CALL create_reference_element(refElPol,2,order, verbose = 0)

      IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
         CALL load_gmsh_mesh("./res/temp",0)
      ELSE
         CALL load_gmsh_mesh("./res/temp",1)
      ENDIF

      CALL mesh_preprocess_serial(ierr)
      CALL read_extended_connectivity('./res/temp.msh')

      Mesh%X = Mesh%X*phys%lscale

      IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
         Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
      END IF

   END SUBROUTINE load_new_mesh_gmsh

   SUBROUTINE load_new_mesh(order)
      USE preprocess
      INTEGER, INTENT(IN) :: order
      INTEGER                     :: ierr

      IF(MPIvar%glob_id .EQ. 0) THEN
         WRITE(*,*) "********** Loading mesh P1  **********"
      ENDIF
      CALL free_mesh

      IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
         CALL load_gmsh_mesh("./res/temp",0)
      ELSE
         CALL load_gmsh_mesh("./res/temp",1)
      ENDIF
      CALL free_reference_element_pol(refElPol)
      CALL create_reference_element(refElPol,2,1, verbose = 0)
      CALL mesh_preprocess_serial(ierr)
  
      Mesh%X = Mesh%X*phys%lscale

      IF(ierr .EQ. 0) THEN
         WRITE(*,*) "Error! Corresponding face in Tb not found. STOP"
         STOP
      ENDIF

      IF(MPIvar%glob_id .EQ. 0) THEN
         CALL HDF5_save_mesh("./newmesh_pre.h5", Mesh%Ndim, mesh%Nelems, mesh%Nextfaces, mesh%Nnodes, mesh%Nnodesperelem, mesh%Nnodesperface, mesh%elemType, mesh%T, mesh%X, mesh%Tb, mesh%boundaryFlag)
      ENDIF
  
      CALL read_extended_connectivity('./res/temp.msh')
  
      CALL set_order_mesh(order)
      CALL free_reference_element_pol(refElPol)
      CALL create_reference_element(refElPol,2,order, verbose = 0)
      CALL mesh_preprocess_serial(ierr)

      Mesh%X = Mesh%X*phys%lscale

      IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
         Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
      END IF

      IF(MPIvar%glob_id .EQ. 0) THEN
         CALL HDF5_save_mesh("./newmesh_notround.h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
      ENDIF

      IF(MPIvar%glob_id .EQ. 0) THEN
         ! overwrite the temp.msh file with the new one with rounded edges (still order 1)
         CALL write_msh_file(Mesh%X,Mesh%T)
         ! convert the mesh to .mesh
         CALL convert_msh2mesh('./res/temp')
      ENDIF

      IF(MPIvar%glob_id .EQ. 0) THEN
         CALL HDF5_save_mesh("./newmesh_round.h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
      ENDIF

   END SUBROUTINE load_new_mesh

   SUBROUTINE jacobian(two_d_nodes, A, B, C, J)
      REAL*8, INTENT(IN)              :: two_d_nodes(:,:)
      REAL*8, INTENT(OUT)             :: J(2,2)
      INTEGER, INTENT(IN)             :: A, B, C
  
      ! Calculate Jacobian matrix
      J(1,1) = two_d_nodes(B,1) - two_d_nodes(A,1)
      J(2,1) = two_d_nodes(B,2) - two_d_nodes(A,2)
      J(1,2) = two_d_nodes(C,1) - two_d_nodes(A,1)
      J(2,2) = two_d_nodes(C,2) - two_d_nodes(A,2)
  
   END SUBROUTINE jacobian

   SUBROUTINE combine_h_target_ind_est(h_map_elements,h_target_elements_est,h_target_elements_ind,h_target_elements)
      REAL*8, INTENT(IN) :: h_map_elements(:)
      REAL*8, INTENT(IN) :: h_target_elements_est(:)
      REAL*8, INTENT(IN) :: h_target_elements_ind(:)
      REAL*8, INTENT(OUT) :: h_target_elements(:)
      REAL*8, PARAMETER :: tol = 1.0E-10
      
      h_target_elements = h_target_elements_est
      WHERE(ABS(h_target_elements_ind-h_map_elements) .LT. tol)
         h_target_elements = h_target_elements_ind
      END WHERE
   ENDSUBROUTINE combine_h_target_ind_est

   SUBROUTINE sum_h_target_nodal(T, h_map_elements, h_target_nodal, nodes_repeats)
      INTEGER, INTENT(IN)               :: T(:,:)
      REAL*8, INTENT(IN)                :: h_map_elements(:)
      REAL*8, INTENT(OUT)               :: h_target_nodal(:)
      INTEGER, INTENT(OUT)              :: nodes_repeats(:)
      INTEGER                           :: i, j

      h_target_nodal = 0.
      nodes_repeats = 0
      DO i=1,SIZE(T,1)
         DO j=1,3
            h_target_nodal(T(i,j)) = h_target_nodal(T(i,j)) + h_map_elements(i)
            nodes_repeats(T(i,j)) = nodes_repeats(T(i,j)) + 1
         ENDDO
      ENDDO
   END SUBROUTINE sum_h_target_nodal

   SUBROUTINE average_h_target(h_target_nodal, nodes_repeats, h_target_vertices)
      REAL*8, INTENT(IN)                              :: h_target_nodal(:)
      INTEGER, INTENT(IN)                             :: nodes_repeats(:)
      REAL*8, DIMENSION(:), POINTER, INTENT(INOUT)      :: h_target_vertices
      INTEGER                                         :: i, j

      j = 1
      DO i=1,SIZE(h_target_nodal)
         IF(nodes_repeats(i) /= 0) THEN
            h_target_vertices(j) = h_target_nodal(i)/REAL(nodes_repeats(i))
            j = j + 1
         ENDIF
      ENDDO
   END SUBROUTINE average_h_target

  SUBROUTINE open_merge_with_geometry(gmsh_l,path2msh)
    TYPE(gmsh_t), INTENT(IN)           :: gmsh_l
    CHARACTER ( len = * ), INTENT(IN) :: path2msh

    CALL gmsh_l%initialize()
    CALL gmsh_l%OPEN(adapt%geometry_path)
    CALL gmsh_l%MERGE(path2msh)
    CALL gmsh_l%option%setNumber("Mesh.MshFileVersion", 2.2)
    CALL gmsh_l%WRITE(path2msh)
    CALL gmsh_l%finalize()

  ENDSUBROUTINE open_merge_with_geometry

  SUBROUTINE set_order_mesh(p)

    INTEGER, INTENT(IN)                      :: p
    TYPE(Reference_element_type)             :: refElLocal
    INTEGER                                  :: n_element_nodes, n_face_nodes, n_int_faces, n_nodes, ini,ind, i, j, elem, iface, counter, element_order, elemType, ielem, meshed_face, &
                                                n_already_meshed_faces, n_boundaries, n_boundary_elements, Nnodesperelem, n_elements,Ndim, Nelems, Nextfaces, Nnodes, &
                                                Nnodesperface, counter1, start, stop_index
    INTEGER                                  :: temp(p-1), elem_pos(2), face_pos(2), face_nodes(3,p-1), face_info(3), elements(2), nodes_face(2, p-1), element_face(2), already_meshed_face_nodes(2, p-1), already_meshed_element_faces(2), &
                                                already_meshed_element(2), ifacenode(p+1)
    INTEGER, ALLOCATABLE, DIMENSION(:,:)     :: int_faces, elem_int_face, Tp, Tb, total_face_info
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: boundaryFlag, unique_boundary_flag
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: aux_coord_logical, local_coord_logical, int_face_meshed
    REAL*8, ALLOCATABLE, DIMENSION(:,:)      :: Xp, Xp_aux, elem_nodes_mod, coord_ref
    INTEGER , ALLOCATABLE, DIMENSION(:,:)    :: Tb_Dirichlet, Tb_LEFT, Tb_RIGHT, Tb_UP, Tb_DOWN, Tb_WALL, Tb_LIM, Tb_IN, Tb_OUT, Tb_ULIM, Tb_PUFF, Tb_PUMP, mesh_info

    IF (utils%printint > 0) THEN
      IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, *) '*************************************************'
       WRITE (6, *) '*              INCREASE ORDER MESH              *'
       WRITE (6, *) '*************************************************'
      ENDIF
    ENDIF

    CALL create_reference_element(refElLocal,2, p, verbose = 1)

    face_nodes = refElLocal%face_nodes(:,2:SIZE(refElLocal%face_nodes,2)-1); ! without vertices
    ALLOCATE(coord_ref(SIZE(refElLocal%coord2d,1), SIZE(refElLocal%coord2d,2)))
    ALLOCATE(int_faces(SIZE(Mesh%intfaces,1), SIZE(Mesh%intfaces,2)))

    coord_ref = refElLocal%coord2d;
    n_element_nodes = SIZE(coord_ref,1);
    n_face_nodes = SIZE(face_nodes,2);


    int_faces = Mesh%intfaces
    elem_int_face = Mesh%F

    n_int_faces = SIZE(int_faces,1);

    DO i = 1, SIZE(elem_int_face,1)
       DO j = 1, SIZE(elem_int_face,2)
          IF(elem_int_face(i,j) .GT. n_int_faces) THEN
             elem_int_face(i,j) = 0
          ENDIF
       ENDDO
    ENDDO

    n_elements = SIZE(Mesh%T,1);
    n_nodes = SIZE(Mesh%X,1);

    ALLOCATE(Xp(n_element_nodes*n_elements,2))
    ALLOCATE(Tp(n_elements,n_element_nodes))
    ALLOCATE(aux_coord_logical(n_element_nodes))
    ALLOCATE(local_coord_logical(n_element_nodes))
    ALLOCATE(int_face_meshed(n_int_faces))
    ALLOCATE(elem_nodes_mod(n_element_nodes,2))

    Xp(1:n_nodes,:) = Mesh%X
    Tp(:,1:3) = Mesh%T

    elem_pos(1) = 1
    elem_pos(2) = 3
    face_pos(1) = 2
    face_pos(2) = 4

    aux_coord_logical = .TRUE.
    aux_coord_logical(1:3) = .FALSE.
    int_face_meshed = .FALSE.

    ini = n_nodes + 1;
    counter1 = 1

    DO elem = 1, n_elements

       local_coord_logical = aux_coord_logical

       ! read the faces infos
       face_info = elem_int_face(elem,:)

       ! count how many faces are valid .neq. 0
       n_already_meshed_faces = 0
       DO iface = 1, 3
          IF(face_info(iface) .NE. 0) THEN
             IF(int_face_meshed(face_info(iface)) .EQV. .TRUE.) THEN
                n_already_meshed_faces = n_already_meshed_faces + 1
             ENDIF
          ENDIF
       ENDDO

       ! if there is at least one face valid
       IF(n_already_meshed_faces .LE. 3) THEN
          ! loop through all faces of the element
          DO iface = 1, 3
             ! if the face is valid
             IF(face_info(iface) .NE. 0) THEN
                ! if it is already remeshed
                IF(int_face_meshed(face_info(iface)) .EQV. .TRUE.) THEN
                   meshed_face = face_info(iface)
                   elements(1) = int_faces(face_info(iface),elem_pos(1))
                   elements(2) = int_faces(face_info(iface),elem_pos(2))

                   IF(elements(1) .NE. elem) THEN
                      already_meshed_element(1) = elements(1)
                      already_meshed_element_faces(1) = int_faces(meshed_face,face_pos(1));
                      already_meshed_face_nodes(1,:) = face_nodes(already_meshed_element_faces(1),:);
                   ELSE
                      element_face(1) = int_faces(meshed_face,face_pos(1));
                      nodes_face(1,:) = face_nodes(element_face(1),:);
                      local_coord_logical(nodes_face(1,:)) = .FALSE.
                   ENDIF

                   IF(elements(2) .NE. elem) THEN
                      already_meshed_element(2) = elements(2)
                      already_meshed_element_faces(2) = int_faces(meshed_face,face_pos(2));
                      already_meshed_face_nodes(2,:) = face_nodes(already_meshed_element_faces(2),:);
                   ELSE
                      element_face(2) = int_faces(meshed_face,face_pos(2));
                      nodes_face(2,:) = face_nodes(element_face(2),:);
                      local_coord_logical(nodes_face(2,:)) = .FALSE.
                   ENDIF

                   IF((elements(1) .NE. elem) .AND. (elements(2) .EQ. elem)) THEN
                      CALL fliplr_int(Tp(already_meshed_element(1),already_meshed_face_nodes(1,:)),temp);
                      Tp(elem,nodes_face(2,:)) = temp
                   ENDIF

                   IF((elements(1) .EQ. elem) .AND. (elements(2) .NE. elem)) THEN
                      CALL fliplr_int(Tp(already_meshed_element(2),already_meshed_face_nodes(2,:)),temp)
                      Tp(elem,nodes_face(1,:)) = temp
                   ENDIF
                ENDIF
                int_face_meshed(face_info(iface)) = .TRUE.
             ENDIF
          ENDDO
       ENDIF

       elem_nodes_mod = 0

       DO i = 1, SIZE(local_coord_logical)
          IF(local_coord_logical(i) .EQV. .TRUE.) THEN
             CALL linear_mapping(Mesh%X(Mesh%T(elem,:),:),coord_ref(i,:), elem_nodes_mod(i,:))
          ENDIF
       ENDDO

       ind = ini

       DO i = 1, SIZE(local_coord_logical)
          IF(local_coord_logical(i) .EQV. .TRUE.) THEN
             Xp(ind,:) = elem_nodes_mod(i,:)
             Tp(elem, i) = ind
             ind = ind + 1
          ENDIF
       ENDDO

       IF((n_element_nodes - n_already_meshed_faces*n_face_nodes - 3 - 1) .GE. 0) THEN
          ini = ini - 1 + n_element_nodes - n_already_meshed_faces*n_face_nodes - 3 + 1
       ENDIF

    ENDDO

    ALLOCATE(Xp_aux(ini-1,2))
    Xp_aux(:,:) = Xp(1:ini-1,:)

    n_face_nodes = SIZE(refElLocal%Face_nodes,2)

    CALL unique_stable(Mesh%boundaryFlag,unique_boundary_flag)

    n_boundaries = SIZE(unique_boundary_flag)

    DO i = 1, n_boundaries
       n_boundary_elements = COUNT(Mesh%boundaryFlag .EQ. unique_boundary_flag(i))
       IF(n_boundary_elements .NE. 0) THEN
          IF(unique_boundary_flag(i) .EQ. 5) THEN
             ALLOCATE(Tb_PUMP(n_boundary_elements, n_face_nodes))
          ELSEIF(unique_boundary_flag(i) .EQ. 6) THEN
             ALLOCATE(Tb_PUFF(n_boundary_elements, n_face_nodes))
          ELSEIF(unique_boundary_flag(i) .EQ. 7) THEN
             ALLOCATE(Tb_LIM(n_boundary_elements, n_face_nodes))
          ELSEIF (unique_boundary_flag(i) .EQ. 8) THEN
             ALLOCATE(Tb_IN(n_boundary_elements, n_face_nodes))
          ELSEIF (unique_boundary_flag(i) .EQ. 9) THEN
             ALLOCATE(Tb_OUT(n_boundary_elements, n_face_nodes))
          ENDIF
       ENDIF
    ENDDO

    start = 0
    stop_index = 0
    DO i = 1, n_boundaries

       start = stop_index + 1
       n_boundary_elements = COUNT(Mesh%boundaryFlag .EQ. unique_boundary_flag(i))
       stop_index = start + n_boundary_elements - 1

       counter = 1
       DO j = start,stop_index
          ielem = Mesh%face_info(j,1)
          iface = Mesh%face_info(j,2)
          ifacenode = refElLocal%face_nodes(iface,:)
          IF(unique_boundary_flag(i) .EQ. 5)  Tb_PUMP(counter,:)  = Tp(ielem,ifacenode)
          IF(unique_boundary_flag(i) .EQ. 6)  Tb_PUFF(counter,:)  = Tp(ielem,ifacenode)
          IF(unique_boundary_flag(i) .EQ. 7)  Tb_LIM(counter,:)   = Tp(ielem,ifacenode)
          IF(unique_boundary_flag(i) .EQ. 8)  Tb_IN(counter,:)    = Tp(ielem,ifacenode)
          IF(unique_boundary_flag(i) .EQ. 9)  Tb_OUT(counter,:)   = Tp(ielem,ifacenode)
          counter = counter + 1
       ENDDO
    ENDDO

    Ndim = SIZE(Xp_aux,2)
    Nelems = SIZE(Tp,1)

    Nextfaces = 0
    IF(ALLOCATED(Tb_IN)) THEN
       Nextfaces = SIZE(Tb_IN,1)
    ENDIF
    IF(ALLOCATED(Tb_LIM)) THEN
       Nextfaces = Nextfaces + SIZE(Tb_LIM,1)
    ENDIF
    IF(ALLOCATED(Tb_OUT)) THEN
       Nextfaces = Nextfaces + SIZE(Tb_OUT,1)
    ENDIF
    IF(ALLOCATED(Tb_PUFF)) THEN
       Nextfaces = Nextfaces + SIZE(Tb_PUFF,1)
    ENDIF
    IF(ALLOCATED(Tb_PUMP)) THEN
       Nextfaces = Nextfaces + SIZE(Tb_PUMP,1)
    ENDIF

    Nnodes = SIZE(Xp_aux,1)
    Nnodesperelem = SIZE(Tp,2)
    Nnodesperface = n_face_nodes
    elemType = Mesh%elemType
    element_order = n_face_nodes

    ALLOCATE(mesh_info(SIZE(Mesh%face_info,1),SIZE(Mesh%face_info,2)))

    mesh_info = Mesh%face_info

    CALL free_mesh

    CALL generate_elemface_info(Tp,Tb_IN, Tb_LIM, Tb_PUFF, Tb_PUMP, Tb_OUT, p+1, total_face_info)
    CALL generate_boundary_names(Tb_Dirichlet, Tb_LEFT, Tb_RIGHT, Tb_UP, Tb_DOWN, Tb_WALL, Tb_LIM, Tb_IN, Tb_OUT, Tb_PUFF, Tb_PUMP, Tb_ULIM, Tb, boundaryFlag, element_order)
    CALL load_mesh2global_var(Ndim, Nelems, Nextfaces, Nnodes, Nnodesperelem, Nnodesperface, elemType, Tp, Xp_aux, Tb, boundaryFlag, total_face_info)

    CALL free_reference_element_pol(refElLocal)
    DEALLOCATE(int_faces, elem_int_face, Tp, Tb)
    DEALLOCATE(total_face_info)
    DEALLOCATE(boundaryFlag, unique_boundary_flag)
    DEALLOCATE(aux_coord_logical, local_coord_logical, int_face_meshed)
    DEALLOCATE(Xp, Xp_aux, elem_nodes_mod, coord_ref)
    IF(ALLOCATED(Tb_Dirichlet)) DEALLOCATE(Tb_Dirichlet)
    IF(ALLOCATED(Tb_PUMP)) DEALLOCATE(Tb_PUMP)
    IF(ALLOCATED(Tb_PUFF)) DEALLOCATE(Tb_PUFF)
    IF(ALLOCATED(Tb_LEFT)) DEALLOCATE(Tb_LEFT)
    IF(ALLOCATED(Tb_RIGHT)) DEALLOCATE(Tb_RIGHT)
    IF(ALLOCATED(Tb_UP)) DEALLOCATE(Tb_UP)
    IF(ALLOCATED(Tb_DOWN)) DEALLOCATE(Tb_DOWN)
    IF(ALLOCATED(Tb_WALL)) DEALLOCATE(Tb_WALL)
    IF(ALLOCATED(Tb_LIM)) DEALLOCATE(Tb_LIM)
    IF(ALLOCATED(Tb_IN)) DEALLOCATE(Tb_IN)
    IF(ALLOCATED(Tb_OUT)) DEALLOCATE(Tb_OUT)
    IF(ALLOCATED(Tb_ULIM)) DEALLOCATE(Tb_ULIM)
  ENDSUBROUTINE set_order_mesh

  SUBROUTINE fliplr_int(arrin, arrout)
    INTEGER, DIMENSION(:), INTENT(IN) :: arrin
    INTEGER, DIMENSION(:), INTENT(INOUT) :: arrout
    INTEGER :: i, n
    INTEGER :: temp

    n = SIZE(arrin)

    DO i = 1, n / 2
       ! Swap elements from the left and right sides
       temp = arrin(i)
       arrout(i) = arrin(n - i + 1)
       arrout(n - i + 1) = temp
    END DO

    ! swap middle element
    IF(MOD(n,2) .EQ. 1) arrout(INT(n/2)+1) = arrin(INT(n/2)+1)
  END SUBROUTINE fliplr_int

  SUBROUTINE linear_mapping(vertCoord, xiVector, X)
    ! Linear mapping between local and cartesian coordinates
    ! Input:
    !   vertCoord: vertexes of the element
    !   xiVector:  point in local coordinates
    ! Output:
    !   X: point in cartesian coordinates

    REAL(8), DIMENSION(:,:), INTENT(IN) :: vertCoord
    REAL(8), DIMENSION(:), INTENT(IN)   :: xiVector
    REAL(8), DIMENSION(:), INTENT(OUT)  :: X

    REAL*8                              :: N(3)


    CALL linear_shape_functions_2D(xiVector, N)

    X(1) = dot_PRODUCT(N, vertCoord(:,1))
    X(2) = dot_PRODUCT(N, vertCoord(:,2))

  CONTAINS

    SUBROUTINE linear_shape_functions_2D(xiVector, N)
      ! Reference triangle is [-1,-1; 1,-1; -1,1]
      REAL(8), DIMENSION(:), INTENT(IN) :: xiVector
      REAL(8), DIMENSION(3), INTENT(OUT) :: N

      REAL(8) :: xi, eta

      xi = xiVector(1)
      eta = xiVector(2)
      N(1) = -xi - eta
      N(2) = 1.0d0 + xi
      N(3) = 1.0d0 + eta
      N = 0.5*N

    ENDSUBROUTINE linear_shape_functions_2D

  ENDSUBROUTINE linear_mapping

  SUBROUTINE inverse_isop_transf(x, Xe, refEl, xieta)
    TYPE(Reference_element_type), INTENT(IN)  :: RefEl
    REAL*8, INTENT(OUT)                       :: xieta(:,:)
    REAL*8, INTENT(IN)                        :: x(:,:), Xe(:,:)
    REAL*8                                    :: x0(SIZE(x,1),SIZE(x,2))
    INTEGER                                   :: maxit, npoints, nnodes, i, j, k, counter, n
    REAL*8, ALLOCATABLE                       :: p(:), dpxi(:), dpeta(:), xind(:,:), x0ind(:,:), xietaind(:,:), rhs(:,:)
    INTEGER, ALLOCATABLE                      :: ind(:)
    REAL*8                                    :: xieta0(SIZE(x,1),SIZE(x,2)), aux_xieta(SIZE(x,1),SIZE(x,2)), Nx(refEl%Nnodes2D), Ny(refEl%Nnodes2D)
    REAL*8                                    :: Vand(refEl%Nnodes2D, refEl%Nnodes2D), invV(refEl%Nnodes2D, refEl%Nnodes2D)
    REAL*8                                    :: Jxx, Jxy, Jyx, Jyy, detJ
    REAL*8                                    :: tol

    maxit   = 5
    tol     = 1e-10
    npoints = SIZE(x,1)
    nnodes  = SIZE(Xe,1)

    CALL inverse_linear_transformation(x, Xe, xieta0)

    ! just fucking brute force it
    DO j = 1, SIZE(xieta0,2)
       DO i = 1, SIZE(xieta0,1)
          IF(ABS(xieta0(i,j)-1.0) .LT. 1e-12) THEN
             xieta0(i,j) = xieta0(i,j) - 1.e-10
          ENDIF
       ENDDO
    ENDDO

    CALL iso_transformation_high_order(xieta0, Xe, refEl, x0)


    n = COUNT(SQRT((x(:,1) - x0(:,1))**2+(x(:,2)-x0(:,2))**2) .GT. (tol*SQRT(x(:,1)**2+x(:,2)**2)+1.e-14))

    IF(n .NE. 0) THEN

       ALLOCATE(ind(n))
       ALLOCATE(rhs(n,SIZE(x,2)))
       ALLOCATE(xind(n,SIZE(x,2)))
       ALLOCATE(x0ind(n, SIZE(x0,2)))
       ALLOCATE(xietaind(n,SIZE(xieta,2)))

       ind = 0.

       counter = 1
       DO i = 1, npoints
          IF((SQRT((x(i,1) - x0(i,1))**2+(x(i,2)-x0(i,2))**2)) .GT. (tol*SQRT(x(i,1)**2+x(i,2)**2)+1.e-14)) THEN
             ind(counter) = counter
             counter = counter + 1
          ENDIF
       ENDDO

       xind     = x(ind,:)
       x0ind    = x0(ind,:)
       xietaind = xieta0(ind,:)

       ALLOCATE (p(refEl%Nnodes2D), dpxi(refEl%Nnodes2D), dpeta(refEl%Nnodes2D))

       DO i = 1, maxit
          IF (ALL((SQRT((xind(:,1)-x0ind(:,1))**2+(xind(:,2)-x0ind(:,2))**2)) .LT. (tol*SQRT(xind(:,1)**2+xind(:,2)**2)+1.e-14))) THEN
             EXIT
          ENDIF

          CALL vandermonde_2d(Vand, refEl)
          CALL invert_matrix(TRANSPOSE(Vand), invV)


          p = 0.d0
          dpxi = 0.d0
          dpeta = 0.d0

          DO j = 1, SIZE(xietaind,1)
             CALL orthopoly2d_deriv(xietaind(j,1), xietaind(j,2), refEl%Ndeg, refEl%Nnodes2D, p, dpxi, dpeta)

             Nx = MATMUL(invV, dpxi)
             Ny = MATMUL(invV, dpeta)

             Jxx = dot_PRODUCT(Nx,Xe(:,1))
             Jxy = dot_PRODUCT(Ny,Xe(:,1))
             Jyx = dot_PRODUCT(Nx,Xe(:,2))
             Jyy = dot_PRODUCT(Ny,Xe(:,2))
             detJ = Jxx*Jyy-Jxy*Jyx
             rhs = xind-x0ind

             xietaind(j,1)=xietaind(j,1)+(rhs(j,1)*Jyy-rhs(j,2)*Jxy)/detJ
             xietaind(j,2)=xietaind(j,2)+(rhs(j,2)*Jxx-rhs(j,1)*Jyx)/detJ

          ENDDO

          ! just fucking brute force it
          DO k = 1, SIZE(xietaind,2)
             DO j = 1, SIZE(xietaind,1)
                IF(ABS(xietaind(j,k)-1.0) .LT. 1e-12) THEN
                   xietaind(j,k) = xietaind(j,k) - 1.e-10
                ENDIF
             ENDDO
          ENDDO

          CALL iso_transformation_high_order(xietaind, Xe, refEl,x0ind)
       ENDDO

       IF(ANY((SQRT((xind(:,1)-x0ind(:,1))**2+(xind(:,2)-x0ind(:,2))**2)) .GT. (tol*SQRT(xind(:,1)**2+xind(:,2)**2)+1.e-14))) THEN
          WRITE(*,*) "inverse_isop_transf non converging."
          STOP
       ENDIF

       aux_xieta = xieta0
       aux_xieta(ind,:) = xietaind
       xieta = aux_xieta

       DEALLOCATE(p)
       DEALLOCATE(dpxi)
       DEALLOCATE(dpeta)
       DEALLOCATE(ind)
       DEALLOCATE(xind)
       DEALLOCATE(x0ind)
       DEALLOCATE(xietaind)
       DEALLOCATE(rhs)

    ELSE
       xieta = xieta0
    ENDIF
  ENDSUBROUTINE inverse_isop_transf

  SUBROUTINE inverse_linear_transformation(x,Xe,xieta)
    USE LinearAlgebra, only: solve_linear_system

    REAL*8, INTENT(IN)      :: x(:,:), Xe(:,:)
    REAL*8, INTENT(OUT)     :: xieta(:,:)
    REAL*8                  :: x1(2),x2(2),x3(2), J(2,2), aux(SIZE(x,1), 2), xieta_temp(SIZE(xieta,2), SIZE(xieta,1))

    ! take vertices
    x1 = Xe(1,:)
    x2 = Xe(2,:)
    x3 = Xe(3,:)


    J(:,1) = (x2-x1)/2
    J(:,2) = (x3-x1)/2

    aux(:,1)  = x(:,1)-(x2(1)+x3(1))/2
    aux(:,2) = x(:,2)-(x2(2)+x3(2))/2

    CALL solve_linear_system(J,TRANSPOSE(aux),xieta_temp)
    xieta = TRANSPOSE(xieta_temp)

  ENDSUBROUTINE inverse_linear_transformation

  SUBROUTINE iso_transformation_high_order(xieta, Xe, refEl, x)
    TYPE(Reference_element_type)      :: refEl
    REAL*8,  INTENT(IN)               :: xieta(:,:)
    REAL*8,  INTENT(IN)               :: Xe(:,:)
    REAL*8,  INTENT(OUT)              :: x(:,:)
    REAL*8                            :: shapeFunctions(refEl%Nnodes2D,SIZE(xieta,1),3)

    shapeFunctions = 0.d0

    CALL compute_shape_functions_at_points(refEl, xieta, shapeFunctions)

    ! take only shape function and not its derivatives in xi and eta
    x(:,1) = MATMUL(TRANSPOSE(shapeFunctions(:,:,1)), Xe(:,1))
    x(:,2) = MATMUL(TRANSPOSE(shapeFunctions(:,:,1)), Xe(:,2))

  ENDSUBROUTINE iso_transformation_high_order

  PURE SUBROUTINE find_matches_int(a, b, indices)
    INTEGER, DIMENSION(:), INTENT(IN)                 :: a
    INTEGER, INTENT(IN)                               :: b
    INTEGER, DIMENSION(:), INTENT(INOUT), ALLOCATABLE :: indices
    INTEGER                                           :: counter
    INTEGER                                           :: i

    counter = COUNT(a .EQ. b)

    ALLOCATE(indices(counter))

    counter = 1
    DO i = 1, SIZE(a)
       IF (a(i) .EQ. b) THEN
          indices(counter) = i
          counter = counter +1
       END IF
    END DO

  END SUBROUTINE find_matches_int

  SUBROUTINE write_msh_file(X, T)

    REAL*8, INTENT(in)                  :: X(:,:)
    INTEGER, INTENT(in)                 :: T(:,:)
    INTEGER, ALLOCATABLE                :: unique_T(:)
    CHARACTER ( len = 255 )             :: buffer
    INTEGER                             :: unit_in,unit_out, i, j, ios

    CALL unique_1D(RESHAPE(T(:,1:3), [SIZE(T(:,1:3),1)*SIZE(T(:,1:3),2)]), unique_T)

    ! rename temp.msh
    CALL rename('./res/temp.msh','./res/temp_origin.msh')
    ! get file_id of temp.msh
    CALL get_unit ( unit_in )
    ! Open the source file for reading
    OPEN(unit=unit_in, file='./res/temp_origin.msh', status='old', action='read', iostat=ios)
    IF (ios /= 0) THEN
       PRINT *, "Error opening source file."
       RETURN
    END IF

    ! get unit and open destination file
    CALL get_unit ( unit_out )
    ! Open the destination file for writing
    OPEN(unit=unit_out, file='./res/temp.msh', status='replace', action='write', iostat=ios)
    IF (ios /= 0) THEN
       PRINT *, "Error opening destination file."
       CLOSE(unit_in)
       RETURN
    END IF

    ! Copy data from source to destination
    DO
       READ (unit_in, '(a)', iostat = ios ) buffer
       IF (ios /= 0) EXIT ! Exit loop when end of file is reached
       WRITE(unit_out, '(a)') buffer(1:100)
       IF(buffer(1:6) .EQ. '$Nodes') EXIT
    END DO

    ! write the number of nodes
    WRITE(unit_out, '(I0)') SIZE(unique_T)

    ! write the new coordinates of the nodes
    DO i = 1, SIZE(unique_T)
       WRITE(unit_out, '(I0, 2F14.10, A)') i, (X(unique_T(i), j), j = 1, 2), ' 0'
    END DO
    WRITE(unit_out, '(A)') "$EndNodes"

    ! copy the connectivity
    ! go down till it encounters $Elements
    DO
       READ (unit_in, '(a)', iostat = ios ) buffer
       IF (ios /= 0) EXIT ! Exit loop if end of file is reached
       IF(buffer(1:9) .EQ. '$Elements') EXIT
    END DO

    ! write $Elements
    WRITE(unit_out, '(A)') "$Elements"

    ! Copy connectivity
    DO
       READ (unit_in, '(a)', iostat = ios ) buffer
       IF (ios /= 0) EXIT ! Exit loop when end of file is reached
       WRITE(unit_out, '(a)') buffer(1:100)
    END DO

    CLOSE(unit_in)
    CLOSE(unit_out)

    CALL delete_file('./res/temp_origin.msh')
    DEALLOCATE(unique_T)

  END SUBROUTINE write_msh_file


  SUBROUTINE read_extended_connectivity(filename)

    CHARACTER * ( * ), INTENT(IN)       :: filename
    INTEGER, ALLOCATABLE                :: T_gmsh(:,:), Tb_gmsh(:,:)
    REAL*8, ALLOCATABLE                 :: X_P1(:,:)
    CHARACTER( LEN = 255 )              :: buffer
    INTEGER                             :: i, ios
    INTEGER ( kind = 4 )                :: unit_gmsh

    ! Open the destination file for writing
    CALL get_unit ( unit_gmsh )
    ! Open the source file for reading
    OPEN(unit=unit_gmsh, file=filename, status='old', action='read', iostat=ios)
    IF (ios /= 0) THEN
       PRINT *, "Error opening destination file."
       CLOSE(unit_gmsh)
       RETURN
    END IF

    ! read till the nodes are found
    DO
       READ (unit_gmsh, '(a)', iostat = ios ) buffer
       IF (ios /= 0) EXIT ! Exit loop if end of file is reached
       IF(buffer(1:6) .EQ. '$Nodes') EXIT
    END DO

    ! read one more line (# of nodes)
    READ (unit_gmsh, '(a)', iostat = ios ) buffer

    ALLOCATE(X_P1(SIZE(Mesh%X,1),2 + Mesh%Ndim))
    ! read coordinates of the nodes
    DO i = 1,SIZE(Mesh%X,1)
       READ(unit_gmsh,*) X_P1(i,:)
    ENDDO

    ! skip $EndNodes, $Elements, #elements
    DO i = 1,3
       ! read one more line
       READ (unit_gmsh, '(a)', iostat = ios ) buffer
    ENDDO

    ALLOCATE(Tb_gmsh(SIZE(Mesh%Tb,1),5 + SIZE(Mesh%Tb,2)))

    DO i = 1,SIZE(Mesh%Tb,1)
       READ(unit_gmsh,*) Tb_gmsh(i,:)
    ENDDO

    ALLOCATE(T_gmsh(SIZE(Mesh%T,1), 5 + SIZE(Mesh%T,2)))
    DO i = 1,SIZE(Mesh%T,1)
       READ(unit_gmsh,*) T_gmsh(i,:)
    ENDDO

    CLOSE(unit_gmsh)

    ALLOCATE(Mesh%Tb_gmsh(SIZE(Tb_gmsh,1),SIZE(Tb_gmsh,2)))
    Mesh%Tb_gmsh = Tb_gmsh
    ALLOCATE(Mesh%T_gmsh(SIZE(T_gmsh,1), SIZE(T_gmsh,2)))
    Mesh%T_gmsh = T_gmsh
    ALLOCATE(Mesh%X_P1(SIZE(X_P1,1), SIZE(X_P1,2)))
    Mesh%X_P1 = X_P1

    DEALLOCATE(Tb_gmsh,T_gmsh, X_P1)

  ENDSUBROUTINE read_extended_connectivity


  SUBROUTINE generate_msh_from_solution_mesh(mesh_name)


    CHARACTER * ( * )                   :: mesh_name
    INTEGER ( kind = 4 )                :: gmsh_unit
    INTEGER                             :: i, n

    ! get unit file and open it
    CALL get_unit ( gmsh_unit )
    OPEN ( unit = gmsh_unit, file = mesh_name, status = 'replace' )

    ! write mesh format
    WRITE ( gmsh_unit, '(a)' ) '$MeshFormat'
    WRITE ( gmsh_unit, '(a)' ) '2.2 0 8'
    WRITE ( gmsh_unit, '(a)' ) '$EndMeshFormat'

    WRITE ( gmsh_unit, '(a)' ) '$PhysicalNames '

    n = 0
    DO i = 1, 10
       IF(COUNT(Mesh%boundaryFlag .EQ. i) .NE. 0) THEN
          n = n + 1
       ENDIF
    ENDDO

    ! + 1 is the domain
    WRITE (gmsh_unit, *) n + 1
    IF(ANY(Mesh%boundaryFlag .EQ. 5)) THEN
       WRITE (gmsh_unit, *) 1,5, '"PUMP"'
    ENDIF
    IF(ANY(Mesh%boundaryFlag .EQ. 6)) THEN
       WRITE (gmsh_unit, *) 1,6, '"PUFF"'
    ENDIF
    IF(ANY(Mesh%boundaryFlag .EQ. 8)) THEN
       WRITE (gmsh_unit, *) 1,1, '"IN"'
    ENDIF
    IF(ANY(Mesh%boundaryFlag .EQ. 9)) THEN
       WRITE (gmsh_unit, *) 1,2, '"OUT"'
    ENDIF
    IF(ANY(Mesh%boundaryFlag .EQ. 7)) THEN
       WRITE (gmsh_unit, *) 1,3, '"LIM"'
    ENDIF

    WRITE (gmsh_unit, *) 2,4, '"DOM"'

    WRITE ( gmsh_unit, '(a)' ) '$EndPhysicalNames '

    ! write nodes
    WRITE ( gmsh_unit, '(a)' ) '$Nodes'
    WRITE ( gmsh_unit, '(i6)' ) SIZE(Mesh%X_P1,1)
    DO i = 1, SIZE(Mesh%X_P1,1)
       WRITE ( gmsh_unit, * ) i, Mesh%X_P1(i,2:)
    END DO
    WRITE ( gmsh_unit, '(a)' ) '$EndNodes'

    ! write elements, Tb extended to gmsh and T extended to gmsh
    WRITE ( gmsh_unit, '(a)' ) '$Elements'
    WRITE ( gmsh_unit, '(i6)' ) SIZE(Mesh%Tb_gmsh,1) + SIZE(Mesh%T_gmsh,1)
    DO i = 1, SIZE(Mesh%Tb_gmsh,1)
       WRITE ( gmsh_unit, *) Mesh%Tb_gmsh(i,:)
    ENDDO
    DO i = 1, SIZE(Mesh%T_gmsh,1)
       WRITE ( gmsh_unit, *) Mesh%T_gmsh(i,:)
    ENDDO

    WRITE ( gmsh_unit, '(a)' ) '$EndElements'

    ! close file
    CLOSE ( unit = gmsh_unit )

  ENDSUBROUTINE generate_msh_from_solution_mesh


  SUBROUTINE convert_msh2mesh(mesh_name)

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)            :: mesh_name
    CHARACTER(LEN = 1024)               :: file_in, file_out
    TYPE(gmsh_t)                        :: gmsh_l

    file_in  = TRIM(ADJUSTL(mesh_name))// '.msh'
    file_out = TRIM(ADJUSTL(mesh_name))// '.mesh'

    CALL gmsh_l%initialize()
    CALL gmsh_l%OPEN(file_in)
    CALL gmsh_l%WRITE(file_out)
    CALL gmsh_l%finalize()

  ENDSUBROUTINE convert_msh2mesh

  SUBROUTINE convert_mesh2msh(mesh_name)

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)            :: mesh_name
    CHARACTER(LEN = 1024)               :: file_in, file_out
    TYPE(gmsh_t)                        :: gmsh_l

    file_in  = TRIM(ADJUSTL(mesh_name))// '.mesh'
    file_out = TRIM(ADJUSTL(mesh_name))// '.msh'

    CALL gmsh_l%initialize()
    CALL gmsh_l%OPEN(file_in)
    CALL gmsh_l%option%setNumber("Mesh.MshFileVersion", 2.2)
    CALL gmsh_l%WRITE(file_out)
    CALL gmsh_l%finalize()

  ENDSUBROUTINE convert_mesh2msh

  SUBROUTINE delete_file(filename)

    CHARACTER(*), INTENT(IN)        :: filename
    INTEGER                         :: fileID, stat

    CALL get_unit ( fileID )   ! get unit file and open it

    OPEN(unit=fileID, iostat=stat, file=filename, status='old')
    IF(stat .NE. 0) THEN
       WRITE(*,*) "Problem opening file to delete."
    ELSE
       CLOSE(unit=fileID, iostat=stat, status='delete')
    ENDIF

    IF(stat .NE. 0) THEN
       WRITE(*,*) "Problem deleting file."
    ENDIF

  ENDSUBROUTINE delete_file

  SUBROUTINE extract_mesh_name_from_fullpath_woext(mesh_name, mesh_name_npne)
    CHARACTER(1024), INTENT(IN)             :: mesh_name
    CHARACTER(1024), INTENT(OUT)            :: mesh_name_npne
    INTEGER                                 :: i, start

    i = LEN(mesh_name)
    start = -1

    DO
       IF(i-2 .EQ. 0) THEN
          WRITE(*,*) "GMSH file input not found, check input sintax."
          STOP
       ENDIF

       ! ! get the index at the dot of .msh
       ! IF((mesh_name(i-2:i-2) .eq. 'm') .and. (mesh_name(i-1:i-1) .eq. 's') .and. (mesh_name(i:i) .eq. 'h')) THEN
       !   end = i-4
       ! ENDIF
       ! get the index at the first slash reading from right to left
       IF(mesh_name(i:i) .EQ. '/') THEN
          start = i+1
          EXIT
       ENDIF
       i = i - 1
    ENDDO

    mesh_name_npne = TRIM(ADJUSTL(mesh_name(start:)))

  ENDSUBROUTINE extract_mesh_name_from_fullpath_woext

  SUBROUTINE mmg_create_mesh_from_h_target(mesh_name)

#include "mmg/mmg2d/libmmg2df.h"

    MMG5_DATA_PTR_T           :: mmgMesh
    MMG5_DATA_PTR_T           :: mmgSol
    INTEGER                   :: ier
    CHARACTER(len=300)        :: filename,filename_sol,fileout
    CHARACTER(*), INTENT(IN)  :: mesh_name


    PRINT*,"  -- Creating new Mesh file from h_target"

    filename       = './res/temp.mesh'
    filename_sol   = './res/ElSizeMap.sol'
    fileout        = mesh_name

    !> ------------------------------ STEP   I --------------------------
    !! 1) Initialisation of mesh and sol structures
    !!   args of InitMesh:
    !! MMG5_ARG_start: we start to give the args of a variadic func
    !! MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
    !! mmgMesh: your MMG5_pMesh (that store your mesh)
    !! MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
    !! mmgSol: your MMG5_pSol (that store your metric) */

    mmgMesh = 0
    mmgSol  = 0

    CALL MMG2D_Init_mesh(MMG5_ARG_start, &
         MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
         MMG5_ARG_end)

    CALL MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_verbose,-1, ier);


    !> 2) Build mesh in MMG5 format
    !! Two solutions: just use the MMG2D_loadMesh function that will read a .mesh(b)
    !! file formatted or manually set your mesh using the MMG2D_Set* functions

    !> with MMG2D_loadMesh function
    CALL MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_3dMedit,2, ier);
    IF ((switch%testcase .LT. 60 .OR. switch%testcase .GT. 80)) THEN
       !CALL MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_nosurf,1, ier);
       CALL MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_xreg,1, ier);

       !call MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_noswap,1, ier);
       !call MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_nomove,1, ier);
       !call MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hausd,1, ier);
       !call MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hausd,0.001, ier);
    END IF

    !CALL MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_xreg,1, ier);
    CALL MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_angle,1, ier);
    CALL MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_angleDetection,0.1, ier);
    !CALL MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hmin,4e-5, ier);
    CALL MMG2D_Set_dparameter(mmgMesh,mmgSol,MMG2D_DPARAM_hgrad,2.4, ier);



    CALL MMG2D_loadMesh(mmgMesh,TRIM(ADJUSTL(filename)),LEN(TRIM(ADJUSTL(filename))),ier)
    IF ( ier /= 1 ) THEN
       WRITE(*,*) "Error loading the .mesh file"
       CALL EXIT(102)
    ENDIF
    !> 3) Build sol in MMG5 format
    !! Two solutions: just use the MMG2D_loadMet function that will read a .sol(b)
    !!    file formatted or manually set your sol using the MMG2D_Set* functions

    !> With MMG2D_loadSol function
    CALL MMG2D_loadSol(mmgMesh,mmgSol,TRIM(ADJUSTL(filename_sol)),LEN(TRIM(ADJUSTL(filename_sol))),ier)
    IF ( ier /= 1 ) THEN
       WRITE(*,*) "Error loading the .sol file"
       CALL EXIT(104)
    ENDIF

    !> 4) (not mandatory): check IF the number of given entities match with mesh size
    CALL MMG2D_Chk_meshData(mmgMesh,mmgSol,ier)
    IF ( ier /= 1 ) THEN
       WRITE(*,*) "Error checking the data for mmg."
       CALL EXIT(107)
    ENDIF

    !> ------------------------------ STEP  II --------------------------
    !! remesh function
    ! NULLIFY(va)
    CALL MMG2D_mmg2dlib(mmgMesh,mmgSol,ier)
    IF ( ier == MMG5_STRONGFAILURE ) THEN
       PRINT*,"BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH"
       STOP MMG5_STRONGFAILURE
    ELSE IF ( ier == MMG5_LOWFAILURE ) THEN
       PRINT*,"BAD ENDING OF MMG2DLIB"
    ELSE
       PRINT*,"MMG2DLIB SUCCEED"
    ENDIF

    !> ------------------------------ STEP III --------------------------
    !! get results
    !! Two solutions: just use the MMG2D_saveMesh/MMG2D_saveSol functions
    !!    that will WRITE .mesh(b)/.sol formatted files or manually get your mesh/sol
    !!    using the MMG2D_getMesh/MMG2D_getSol functions

    !> 1) Automatically save the mesh
    CALL MMG2D_saveMesh(mmgMesh,TRIM(ADJUSTL(fileout)),LEN(TRIM(ADJUSTL(fileout))),ier)
    IF ( ier /= 1 ) THEN
       WRITE(*,*) "Error checking the data for mmg."
       CALL EXIT(106)
    ENDIF

    !> 2) Automatically save the solution
    CALL MMG2D_saveSol(mmgMesh,mmgSol,TRIM(ADJUSTL(fileout)),LEN(TRIM(ADJUSTL(fileout))),ier)
    IF ( ier /= 1 ) THEN
       WRITE(*,*) "Error saving the solution."
       CALL EXIT(107)
    ENDIF

    !> 3) Free the MMG2D5 structures
    CALL MMG2D_Free_all(MMG5_ARG_start, &
         MMG5_ARG_ppMesh,mmgMesh,MMG5_ARG_ppMet,mmgSol, &
         MMG5_ARG_end)

  ENDSUBROUTINE mmg_create_mesh_from_h_target

  SUBROUTINE generate_htarget_sol_file(N_n_vertex, h_target)

    INTEGER, INTENT(IN)         :: N_n_vertex
    REAL*8, INTENT(IN)          :: h_target(:)
    INTEGER                     :: fileID


    CALL get_unit ( fileID )
    ! Open the file for writing
    OPEN(unit=fileID, file='./res/ElSizeMap.sol')

    ! WRITE data to the file
    WRITE(fileID, '(A,I0)') 'MeshVersionFormatted ', 2
    WRITE(fileID, *)
    WRITE(fileID, '(A,I0)') 'Dimension ', 3
    WRITE(fileID, *)
    WRITE(fileID, '(A)') 'SolAtVertices'
    WRITE(fileID, '(I0)') N_n_vertex
    WRITE(fileID, '(I1, 1X, I1)') 1, 1
    WRITE(fileID, *)
    WRITE(fileID, '(F8.6)') h_target
    WRITE(fileID, *)
    WRITE(fileID, '(A)') 'End'

    ! Close the file
    CLOSE(fileID)

  ENDSUBROUTINE generate_htarget_sol_file

  PURE SUBROUTINE unique_1D(list_in, list_out)
    !! From a 1D array of integers list_in extracts the list of unique occurences of values
    !integer, dimension(:), intent(in) :: list_in
    INTEGER, DIMENSION(:), INTENT(in) :: list_in
    !! The list of integers to work on
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: list_out
    !! The list in output, with a single occurence of each value in list in

    INTEGER, DIMENSION(SIZE(list_in)) :: list_in_sorted
    INTEGER :: nlist, n_unique, i

    nlist = SIZE(list_in)

    ! Copy list_in and sort it
    list_in_sorted = list_in
    CALL quicksort_int(list_in_sorted)

    ! The number of jumps in this list gives the number of unique values
    n_unique = COUNT((list_in_sorted(2:nlist)-list_in_sorted(1:nlist-1)) .GT. 0) + 1
    ALLOCATE(list_out(n_unique))

    ! Now loop in the sorted list and each time a jump is found,
    n_unique = 1
    list_out(1) = list_in_sorted(1)
    DO i = 2, nlist
       IF (list_in_sorted(i).NE.list_in_sorted(i-1)) THEN
          n_unique = n_unique + 1
          list_out(n_unique) = list_in_sorted(i)
       ENDIF
    ENDDO

  END SUBROUTINE unique_1D

  PURE RECURSIVE SUBROUTINE quicksort_int(a)
    !! quicksort.f -*-f90-*-
    !! Author: t-nissie, some tweaks by 1AdAstra1
    !! License: GPLv3
    !! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    INTEGER, DIMENSION(:), INTENT(inout) :: a
    INTEGER ::  x, t
    INTEGER :: first, last
    INTEGER i, j

    first = 1

    last = SIZE(a, 1)
    x = a( (first+last) / 2 )
    i = first
    j = last

    DO
       DO WHILE (a(i) < x)
          i=i+1
       ENDDO
       DO WHILE (x < a(j))
          j=j-1
       ENDDO
       IF (i .GE. j) EXIT
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    ENDDO

    IF (first < i - 1) CALL quicksort_int(a(first : i - 1))
    IF (j + 1 < last)  CALL quicksort_int(a(j + 1 : last))

  ENDSUBROUTINE quicksort_int

  PURE RECURSIVE SUBROUTINE quicksort_real(a)
    !! quicksort.f -*-f90-*-
    !! Author: t-nissie, some tweaks by 1AdAstra1
    !! License: GPLv3
    !! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    REAL*8, DIMENSION(:), INTENT(inout) :: a
    !! The array to sort

    REAL*8 ::  x, t
    INTEGER :: first, last
    INTEGER i, j

    first = 1
    last = SIZE(a, 1)
    x = a( (first+last) / 2 )
    i = first
    j = last

    DO
       DO WHILE (a(i) < x)
          i=i+1
       ENDDO
       DO WHILE (x < a(j))
          j=j-1
       ENDDO
       IF (i .GE. j) EXIT
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    ENDDO

    IF (first < i - 1) CALL quicksort_real(a(first : i - 1))
    IF (j + 1 < last)  CALL quicksort_real(a(j + 1 : last))

  ENDSUBROUTINE quicksort_real

  PURE SUBROUTINE unique_stable(arrayin, uniqueArr)
    INTEGER, DIMENSION(:), INTENT(IN) :: arrayin
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: uniqueArr
    INTEGER, DIMENSION(SIZE(arrayin)) :: indices
    LOGICAL, DIMENSION(SIZE(arrayin)) :: isUnique
    INTEGER :: i, j, nUnique


    nUnique = 0
    isUnique = .TRUE.
    indices = 0

    ! Loop through the original array to find unique elements
    DO i = 1, SIZE(arrayin)
       IF (isUnique(i)) THEN
          nUnique = nUnique + 1
          indices(nUnique) = i
          DO j = i + 1, SIZE(arrayin)
             IF (arrayin(j) == arrayin(i)) THEN
                isUnique(j) = .FALSE.
             END IF
          END DO
       END IF
    END DO

    ! Allocate the unique array and populate it in the original order
    ALLOCATE(uniqueArr(nUnique))
    DO i = 1, nUnique
       uniqueArr(i) = arrayin(indices(i))
    END DO

  ENDSUBROUTINE unique_stable


  SUBROUTINE gmsh_create_from_h_target(h_target_on_elements,vertices_coordinates, connectivity, p_order)

      USE, INTRINSIC :: iso_c_binding
      USE gmsh
      USE MPI_OMP, only: OMPvar

      TYPE(gmsh_t) :: gmsh_l
      INTEGER, INTENT(IN) :: p_order
      REAL*8, DIMENSION(:), INTENT(IN) :: h_target_on_elements ! on the elements
      REAL*8, DIMENSION(:,:), INTENT(IN) :: vertices_coordinates ! coordinates of the vertices
      INTEGER, DIMENSION(:,:), INTENT(IN) :: connectivity ! connectivity of the triangles
      REAL*8, ALLOCATABLE :: data_for_gmsh(:)
      INTEGER*8           :: gmsh_dim, number_of_vertices_per_triangle, i, j, start_index
      INTEGER*4           :: size_view,number_of_triangles,ret
      REAL*8              :: sf_index, vertex_coordinates(2),h_target_on_vertex
      !GMSH always have (X,Y,Z) coordinates
      gmsh_dim = 3 
      number_of_vertices_per_triangle = 3 
      number_of_triangles = SIZE(connectivity,1)

      ! initalize gmsh
      CALL gmsh_l%initialize()

      ! Set verbosity level to 2 (Errors and warnings)
      CALL gmsh_l%option%setNumber("General.Verbosity", 2.0)

      !create model
      CALL gmsh_l%model%add("geo")
      ! merge
      CALL gmsh_l%merge(adapt%geometry_path)


      ALLOCATE(data_for_gmsh((number_of_vertices_per_triangle*(number_of_vertices_per_triangle+1))*number_of_triangles))

      ! Prepare data in gmsh format
      DO i = 1, number_of_triangles
         start_index = (i-1)*(number_of_vertices_per_triangle*(number_of_vertices_per_triangle+1))
         DO j = 1, number_of_vertices_per_triangle
            vertex_coordinates = vertices_coordinates(connectivity(i,j),:)
            h_target_on_vertex = h_target_on_elements(i)
            data_for_gmsh(start_index+j) = vertex_coordinates(1)
            data_for_gmsh(start_index+number_of_vertices_per_triangle+j) = vertex_coordinates(2)
            data_for_gmsh(start_index+2*number_of_vertices_per_triangle+j) = 0.0 ! no Z coordinate
            data_for_gmsh(start_index+3*number_of_vertices_per_triangle+j) = h_target_on_vertex
         ENDDO
      ENDDO

      ! Add h_target as a post-processing view
      size_view = gmsh_l%view%add("h_target")
      call gmsh_l%view%addListData(size_view, "ST",number_of_triangles,data_for_gmsh)
      sf_index = gmsh_l%view%getIndex(size_view)

      ! Add the view as a field
      ret = gmsh_l%model%mesh%field%add("PostView")
      call gmsh_l%model%mesh%field%setNumber(ret, "ViewIndex", 0d0)

      ! Apply the view as the current background mesh size field:
      call gmsh_l%model%mesh%field%setAsBackgroundMesh(ret)

      ! ignore characteristic length from geometry
      call gmsh_l%option%setNumber("Mesh.MeshSizeExtendFromBoundary", 0d0)
      call gmsh_l%option%setNumber("Mesh.MeshSizeFromPoints", 0d0)
      call gmsh_l%option%setNumber("Mesh.MeshSizeFromCurvature", 0d0)
      CALL gmsh_l%option%setNumber("Mesh.MeshSizeFactor", 1d0)
      CALL gmsh_l%option%setNumber("General.NumThreads", REAL(OMPvar%Nthreads))

      !Changing the algorithm to Delaunay, the default is Frontal-Delaunay (Don't Know if needed)
      call gmsh_l%option%setNumber("Mesh.Algorithm", 5d0)

      ! Generate the refined mesh
      call gmsh_l%model%mesh%generate(2)
      call gmsh_l%model%mesh%setOrder(p_order)
      call gmsh_l%model%mesh%optimize('HighOrder')
      CALL gmsh_l%option%setNumber("Mesh.MshFileVersion", 2.2)
      call gmsh_l%write('./res/temp.msh')
      CALL gmsh_l%finalize()

      DEALLOCATE(data_for_gmsh)

   END SUBROUTINE gmsh_create_from_h_target

   SUBROUTINE save_copy_new_mesh(mesh_name, count_adapt)
      USE in_out, ONLY: copy_file
      CHARACTER(1024), INTENT(IN) :: mesh_name
      INTEGER, INTENT(IN)         :: count_adapt
      CHARACTER(70)               :: param_adapt_char, count_adapt_char
      CHARACTER(1024)             :: mesh_name_npne,new_mesh_name_npne, buffer

      CALL extract_mesh_name_from_fullpath_woext(mesh_name, mesh_name_npne)

      WRITE(param_adapt_char, *) adapt%param_est
      WRITE(count_adapt_char, *) count_adapt
      new_mesh_name_npne = TRIM(ADJUSTL(mesh_name_npne)) // '_param'// TRIM(ADJUSTL(param_adapt_char)) // '_n' // TRIM(ADJUSTL(count_adapt_char))
      
      buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".msh"
      CALL copy_file("./res/temp.msh",buffer)

   END SUBROUTINE save_copy_new_mesh
END MODULE adaptivity_common_module
