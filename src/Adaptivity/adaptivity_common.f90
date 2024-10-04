!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! Module for schock capturing adaptivity
!************************************************************

MODULE adaptivity_common_module
  USE globals
  USE reference_element
  USE LinearAlgebra
  USE gmsh
  IMPLICIT NONE

CONTAINS

  SUBROUTINE merge_with_geometry(gmsh)
    TYPE(gmsh_t), INTENT(IN)          :: gmsh

    CALL gmsh%initialize()

    IF (switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) THEN
       CALL gmsh%OPEN("./res/geometries/Circ_InfLim_YesHole_Structured.geo")
    ELSE
       !CALL gmsh%open("./res/West_Mesh_farWall_NoHole_SmoothCorner_base.geo")
       IF(ANY(Mesh%boundaryFlag .EQ. 5)) THEN
          !CALL gmsh%OPEN("./res/geometries/West_Mesh_NoHole_farWall.geo")
          CALL gmsh%OPEN("./res/geometries/TCV_smooth.geo")
       ELSEIF(ANY(Mesh%boundaryFlag .EQ. 8)) THEN
          CALL gmsh%OPEN("./res/geometries/West_Mesh_YesHole_SmoothCorner.geo")
       ELSE
          CALL gmsh%OPEN("./res/geometries/West_Mesh_NoHole_SmoothCorner.geo")
       ENDIF
    END IF

    CALL gmsh%MERGE("./res/temp.msh")
    CALL gmsh%option%setNumber("Mesh.MshFileVersion", 2.2)
    CALL gmsh%WRITE("./res/temp.msh")
    CALL gmsh%finalize()

  ENDSUBROUTINE merge_with_geometry

  SUBROUTINE open_merge_with_geometry(gmsh,path2msh)
    TYPE(gmsh_t), INTENT(IN)           :: gmsh
    CHARACTER ( len = * ), INTENT(IN) :: path2msh

    CALL gmsh%initialize()
    IF (switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) THEN
       CALL gmsh%OPEN("./res/geometries/Circ_InfLim_YesHole_Structured.geo")
    ELSE
       !CALL gmsh%open("./res/West_Mesh_farWall_NoHole_SmoothCorner_base.geo")
       IF(ANY(Mesh%boundaryFlag .EQ. 5)) THEN
          !CALL gmsh%OPEN("./res/geometries/West_Mesh_NoHole_farWall.geo")
          CALL gmsh%OPEN("./res/geometries/TCV_smooth.geo")
       ELSEIF(ANY(Mesh%boundaryFlag .EQ. 8)) THEN
          CALL gmsh%OPEN("./res/geometries/West_Mesh_YesHole_SmoothCorner.geo")
       ELSE
          CALL gmsh%OPEN("./res/geometries/West_Mesh_NoHole_SmoothCorner.geo")
       ENDIF
    END IF
    CALL gmsh%MERGE(path2msh)
    CALL gmsh%option%setNumber("Mesh.MshFileVersion", 2.2)
    CALL gmsh%WRITE(path2msh)
    CALL gmsh%finalize()

  ENDSUBROUTINE open_merge_with_geometry

  SUBROUTINE set_order_mesh(p)
    USE gmsh_io_module, ONLY: generate_boundary_names, generate_elemface_info, load_mesh2global_var

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

  SUBROUTINE round_edges(Mesh_loc)
    USE mod_splines
    USE HDF5_io_module

    ! Variables
    TYPE(Mesh_type), INTENT(IN)       :: Mesh_loc

    REAL*8, ALLOCATABLE               :: X(:,:)
    INTEGER, ALLOCATABLE              :: T(:,:), Tb(:,:), extFaces(:)
    INTEGER                           :: n_boundaries, n, i, j, k, l, ibound, isp, counter, counter2, Ndim, Nelems, Nextfaces, Nnodesperface, Nexternalnodes, Nnodesperelem, Ninnerfnodes, Ninnerenodes, Difffenodes, Sumfenodes, tot_n_splines
    REAL*8                            :: RotMat(2,2), elem_nodes_mod(Mesh_loc%Nnodesperelem - refElPol%Nvertices,2)
    REAL*8                            :: xcenter, ycenter, xmin, xmax, ymin, ymax
    INTEGER, ALLOCATABLE              :: Tb_bound(:), Tb_bound_vtx(:), extFaces_bound(:), indices(:)
    REAL*8, ALLOCATABLE               :: nodes_on_spline(:,:),nodes_on_spline0(:,:), nodes_on_spline_proj(:,:), opposite_vertex(:,:), m1(:),  m2(:), q1(:), q2(:), m(:), q(:), displacement(:,:), gamma_par(:,:), norm(:)
    REAL*8, PARAMETER                 :: tol = 1.e-4, tol_iter = 1.e-12
    INTEGER, PARAMETER                :: max_iter = 10
    LOGICAL                           :: check_vtx(2)

    n_boundaries = MAXVAL(Mesh_loc%boundaryFlag)
    Nnodesperface = Mesh_loc%Nnodesperface
    Nnodesperelem = Mesh_loc%Nnodesperelem
    Ninnerfnodes = Nnodesperface - 2
    Ninnerenodes = Mesh_loc%Nnodesperelem - refElPol%Nvertices*(Ninnerfnodes + 1)
    Nexternalnodes = Mesh_loc%Nnodesperelem-Ninnerenodes
    Difffenodes = Ninnerenodes - Ninnerfnodes
    Sumfenodes = Ninnerenodes + Ninnerfnodes
    Ndim = Mesh_loc%Ndim
    Nelems = Mesh_loc%Nelems
    Nextfaces = Mesh_loc%NextFaces

    ALLOCATE(X(Mesh_loc%Nnodes, Ndim))
    ALLOCATE(T(Nelems, Mesh_loc%Nnodesperelem))
    ALLOCATE(Tb(Nextfaces, Nnodesperface))
    ALLOCATE(extFaces(Nextfaces))

    X = Mesh_loc%X
    T = Mesh_loc%T
    Tb = Mesh_loc%Tb
    extFaces = Mesh_loc%extFaces(:,1)

    tot_n_splines = splines(1)%tot_n_splines

    DO k = 1,2
       DO ibound = 1, n_boundaries
          n = COUNT(Mesh_loc%boundaryFlag .EQ. ibound)
          IF(n .EQ. 0) CYCLE

          ! count how many vertices are on the splines that define the current boundary
          n = 0
          DO isp = 1, tot_n_splines
             IF(splines(isp)%boundaryFlag .NE. ibound) CYCLE
             xmin = splines(isp)%xmin - tol
             xmax = splines(isp)%xmax + tol
             ymin = splines(isp)%ymin - tol
             ymax = splines(isp)%ymax + tol

             DO i = 1, NextFaces
                IF (Mesh_loc%boundaryFlag(i) .EQ. ibound) THEN
                   IF(ALL(((X(Tb(i, [1,Nnodesperface]),1) .GE. xmin) .AND. (X(Tb(i, [1,Nnodesperface]),1) .LE. xmax)) .AND. ((X(Tb(i, [1,Nnodesperface]),2) .GE. ymin) .AND. (X(Tb(i, [1,Nnodesperface]),2) .LE. ymax)))) THEN
                      ! as a second check, slower, check if any distance of the points to any point in the spline is less than 0.02
                      check_vtx = .FALSE.
                      counter2 = 1
                      DO l = 1,Nnodesperface,Nnodesperface-1
                         DO j = 1, splines(isp)%n_points
                            IF(NORM2(splines(isp)%points_coord(j,:) - X(Tb(i, l), :)) .LT. 0.02) THEN
                               check_vtx(counter2) = .TRUE.
                               EXIT
                            ENDIF
                         ENDDO
                         counter2 = counter2 + 1
                      ENDDO

                      IF(ALL(check_vtx)) THEN
                         n = n + 1  ! Increment counter for next iteration
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO

          IF(n .EQ. 0) CYCLE

          ALLOCATE(Tb_bound_vtx(2*n))
          IF(k .GT. 1) THEN
             ALLOCATE(Tb_bound(n*(Sumfenodes)))
             ALLOCATE(extFaces_bound(n*(Sumfenodes)))
             ALLOCATE(opposite_vertex(n*(Sumfenodes),2))
             ALLOCATE(displacement(n*(Sumfenodes),2))
             ALLOCATE(gamma_par(n*(Sumfenodes),2))
             ALLOCATE(norm(n*(Sumfenodes)))
             ALLOCATE(m(n))
             ALLOCATE(q(n))
             ALLOCATE(m2(n*(Sumfenodes)))
             ALLOCATE(q2(n*(Sumfenodes)))
          ENDIF

          counter = 1
          DO isp = 1, tot_n_splines
             IF(splines(isp)%boundaryFlag .NE. ibound) CYCLE

             xmin = splines(isp)%xmin - tol
             xmax = splines(isp)%xmax + tol
             ymin = splines(isp)%ymin - tol
             ymax = splines(isp)%ymax + tol

             DO i = 1, NextFaces
                IF (Mesh_loc%boundaryFlag(i) .EQ. ibound) THEN
                   ! as a first check, check if the coordinates of the points are in the range of the spline
                   IF(ALL(((X(Tb(i, [1,Nnodesperface]),1) .GE. xmin) .AND. (X(Tb(i, [1,Nnodesperface]),1) .LE. xmax)) .AND. ((X(Tb(i, [1,Nnodesperface]),2) .GE. ymin) .AND. (X(Tb(i, [1,Nnodesperface]),2) .LE. ymax)))) THEN
                      ! as a second check, slower, check if any distance of the points to any point in the spline is less than 0.02
                      check_vtx = .FALSE.
                      counter2 = 1
                      DO l = 1,Nnodesperface,Nnodesperface-1
                         DO j = 1, splines(isp)%n_points
                            IF(NORM2(splines(isp)%points_coord(j,:) - X(Tb(i, l), :)) .LT. 0.02) THEN
                               check_vtx(counter2) = .TRUE.
                               EXIT
                            ENDIF
                         ENDDO
                         counter2 = counter2 + 1
                      ENDDO

                      IF(ALL(check_vtx))  THEN
                         ! Extract node numbers for face vertices nodes
                         Tb_bound_vtx(counter:counter+1) = Tb(i, [1,Nnodesperface])
                         counter = counter + 2  ! Increment counter for next iteration
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO

          IF(k .GT. 1) THEN
             ! this flattens the array at the same time
             counter = 1
             DO isp = 1, tot_n_splines
                IF(splines(isp)%boundaryFlag .NE. ibound) CYCLE

                xmin = splines(isp)%xmin - tol
                xmax = splines(isp)%xmax + tol
                ymin = splines(isp)%ymin - tol
                ymax = splines(isp)%ymax + tol

                DO i = 1, NextFaces
                   IF (Mesh_loc%boundaryFlag(i) .EQ. ibound) THEN
                      ! as a first check, check if the coordinates of the points are in the range of the spline
                      IF(ALL(((X(Tb(i, 2:Nnodesperface-1),1) .GE. xmin) .AND. (X(Tb(i, 2:Nnodesperface-1),1) .LE. xmax)) .AND. ((X(Tb(i, 2:Nnodesperface-1),2) .GE. ymin) .AND. (X(Tb(i, 2:Nnodesperface-1),2) .LE. ymax)))) THEN
                         ! as a second check, slower, check if any distance of the vertices to any point in the spline is less than 0.02
                         check_vtx = .FALSE.
                         counter2 = 1
                         DO l = 1,Nnodesperface,Nnodesperface-1
                            DO j = 1, splines(isp)%n_points
                               IF(NORM2(splines(isp)%points_coord(j,:) - X(Tb(i, l), :)) .LT. 0.02) THEN
                                  check_vtx(counter2) = .TRUE.
                                  EXIT
                               ENDIF
                            ENDDO
                            counter2 = counter2 + 1
                         ENDDO

                         IF(ALL(check_vtx)) THEN
                            ! Extract node numbers for inner face nodes from Tb
                            Tb_bound(counter:counter+Ninnerfnodes-1) = Tb(i, 2:Nnodesperface-1)
                            extFaces_bound(counter:counter+Ninnerfnodes-1) = extfaces(i)
                            counter = counter + Ninnerfnodes  ! Increment counter for next iteration

                            ! Extract node numbers for inner element nodes from T array
                            Tb_bound(counter:counter+Ninnerenodes-1) = T(extFaces(i), Nexternalnodes+1:)
                            extFaces_bound(counter:counter+Ninnerenodes-1) = extfaces(i)
                            counter = counter + Ninnerenodes  ! Increment counter for next iteration
                         ENDIF
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO

             ! compute displacement inner element nodes
             CALL displacement_innerenodes

             ! move the inner element nodes to the boundary element face
             X(Tb_bound,:) = X(Tb_bound,:) + displacement

          ELSE
             CALL unique_1D(Tb_bound_vtx,Tb_bound)
          ENDIF

          DO isp = 1, tot_n_splines
             IF(splines(isp)%boundaryFlag .NE. ibound) CYCLE

             xmin = splines(isp)%xmin - tol
             xmax = splines(isp)%xmax + tol
             ymin = splines(isp)%ymin - tol
             ymax = splines(isp)%ymax + tol

             ! count number of points in the mesh that are supposed to be on the spline
             !n = COUNT(((X(Tb_bound,1) .ge. xmin) .and. (X(Tb_bound,1) .le. xmax)) .and. ((X(Tb_bound,2) .ge. ymin) .and. (X(Tb_bound,2) .le. ymax)))

             n = 0
             DO i = 1, SIZE(Tb_bound)
                IF(((X(Tb_bound(i),1) .GE. xmin) .AND. (X(Tb_bound(i),1) .LE. xmax)) .AND. ((X(Tb_bound(i),2) .GE. ymin) .AND. (X(Tb_bound(i),2) .LE. ymax))) THEN
                   DO j = 1, splines(isp)%n_points
                      IF(NORM2(splines(isp)%points_coord(j,:) - X(Tb_bound(i),:)) .LT. 0.02) THEN
                         n = n + 1
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO

             IF(n .NE. 0) THEN

                ALLOCATE(nodes_on_spline(n,Ndim))
                ALLOCATE(nodes_on_spline_proj(n,Ndim))
                ALLOCATE(nodes_on_spline0(n,Ndim))
                ALLOCATE(opposite_vertex(n,2))
                ALLOCATE(m1(n))
                ALLOCATE(m2(n))
                ALLOCATE(q1(n))
                ALLOCATE(q2(n))
                ALLOCATE(indices(n))

                ! select the nodes that are on the spline
                counter = 1
                DO i = 1, SIZE(Tb_bound)
                   IF(((X(Tb_bound(i),1) .GE. xmin) .AND. (X(Tb_bound(i),1) .LE. xmax)) .AND. ((X(Tb_bound(i),2) .GE. ymin) .AND. (X(Tb_bound(i),2) .LE. ymax))) THEN
                      DO j = 1, splines(isp)%n_points
                         IF(NORM2(splines(isp)%points_coord(j,:) - X(Tb_bound(i),:)) .LT. 0.02) THEN
                            nodes_on_spline(counter,:) = X(Tb_bound(i),:)
                            indices(counter) = i
                            counter = counter + 1
                            EXIT
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO

                IF(k .EQ. 1) THEN
                   WRITE(*,*) "Cubic Spline Projection Vertices"
                ELSE
                   WRITE(*,*) "Cubic Spline Projection Internals"
                ENDIF

                ! rotate the reference frame
                CALL rotate(nodes_on_spline)

                nodes_on_spline0 = nodes_on_spline

                IF(k .NE. 1) THEN
                   CALL get_line_opposite_vertex
                ENDIF

                counter = 0
                DO

                   IF(k .EQ. 1) THEN
                      CALL project_perp_to_spline
                   ELSE
                      CALL project_from_opposite_vertex
                   ENDIF

                   IF(ALL(NORM2((nodes_on_spline - nodes_on_spline_proj)/nodes_on_spline,2) .LE. tol_iter) .OR. (counter .EQ. max_iter)) THEN
                      nodes_on_spline = nodes_on_spline_proj
                      EXIT
                   ELSE
                      nodes_on_spline = nodes_on_spline_proj
                   ENDIF

                   counter = counter + 1
                ENDDO

                CALL anti_rotate(nodes_on_spline)

                ! move back the inner element nodes
                IF(k .NE. 1) THEN
                   nodes_on_spline_proj = X(Tb_bound(indices),:) - displacement(indices,:)
                   CALL anti_rotate(opposite_vertex)

                   gamma_par(indices,1) = (nodes_on_spline_proj(:,1)-opposite_vertex(:,1)) * NORM2(opposite_vertex-nodes_on_spline,2)/norm(indices)
                   gamma_par(indices,2) = (nodes_on_spline_proj(:,2)-opposite_vertex(:,2)) * NORM2(opposite_vertex-nodes_on_spline,2)/norm(indices)

                   X(Tb_bound(indices),:) =  opposite_vertex + gamma_par(indices,:)
                ELSE
                   X(Tb_bound(indices),:) = nodes_on_spline(:,:)
                ENDIF

                DEALLOCATE(nodes_on_spline)
                DEALLOCATE(nodes_on_spline_proj)
                DEALLOCATE(nodes_on_spline0)
                DEALLOCATE(m1, m2, q1, q2)
                DEALLOCATE(indices)
                DEALLOCATE(opposite_vertex)

             ENDIF
          ENDDO

          DEALLOCATE(Tb_bound)
          DEALLOCATE(Tb_bound_vtx)
          IF(k .GT. 1) THEN
             DEALLOCATE(extFaces_bound)
             DEALLOCATE(displacement)
             DEALLOCATE(gamma_par)
             DEALLOCATE(norm)
             DEALLOCATE(m)
             DEALLOCATE(q)
          ENDIF
       ENDDO

       IF(k .EQ. 1) THEN
          DO i = 1, Nelems
             DO j = 1, Nnodesperelem-3
                CALL linear_mapping(X(T(i,1:3),:),refElPol%coord2d(j+3,:), elem_nodes_mod(j,:))
             ENDDO
             X(T(i,4:),:) = elem_nodes_mod
          ENDDO
       ENDIF

    ENDDO

    CALL generate_fekete_nodes(X,T,refElPol%nDeg,refElPol)

    Mesh_loc%X = X

    DEALLOCATE(X,T,Tb, extFaces)

  CONTAINS

    SUBROUTINE displacement_innerenodes

      ! find coordinates of vertex opposite to each face
      DO i = 1, n*Sumfenodes
         DO j = 1, 3
            IF(COUNT(T(extFaces_bound(i),j) .EQ. Tb_bound_vtx) .EQ. 0) THEN
               opposite_vertex(i,:) = X(T(extFaces_bound(i),j),:)
            ENDIF
         ENDDO
      ENDDO

      ! find slope and intercept of line passing through the vertex and the point to project
      m2 = (opposite_vertex(:,2) - X(Tb_bound,2))/(opposite_vertex(:,1) - X(Tb_bound,1))
      q2 = opposite_vertex(:,2) - m2*opposite_vertex(:,1)

      !  get slope of the line passing by the vertices
      counter = 1
      DO i = 1, 2*n-1, 2
         m(counter) = (X(Tb_bound_vtx(i+1),2) - X(Tb_bound_vtx(i),2))/(X(Tb_bound_vtx(i+1),1) - X(Tb_bound_vtx(i),1))
         counter = counter + 1
      ENDDO

      ! get intercept of the line passing by the vertices
      counter = 1
      DO i = 1, 2*n, 2
         q(counter) = X(Tb_bound_vtx(i),2)-m(counter)*X(Tb_bound_vtx(i),1)
         counter = counter + 1
      ENDDO

      ! find intersection between the two lines (here i use displacement to store the new coordinates)
      counter = 1
      DO i = 1, n
         displacement(counter:counter + Sumfenodes-1 ,1) = (q(i)-q2(counter:counter + Sumfenodes -1))/(m2(counter:counter + Sumfenodes -1)-m(i))
         counter = counter + Sumfenodes
      ENDDO

      displacement(:,2) = m2*displacement(:,1)+q2

      ! a-c = ca vector
      gamma_par = opposite_vertex - displacement

      norm = NORM2(gamma_par,2)

      displacement = displacement - X(Tb_bound,:)

      DO i = 1, n*Sumfenodes
         IF(NORM2(displacement(i,:)) .LT. 1e-12) THEN
            displacement(i,:) = 0
         ENDIF
      ENDDO

      DEALLOCATE(opposite_vertex, m2, q2)

    ENDSUBROUTINE displacement_innerenodes

    SUBROUTINE project_from_opposite_vertex

      ! project to spline and find y1
      nodes_on_spline_proj(:,1) = nodes_on_spline(:,1)
      nodes_on_spline_proj(:,2) = splines(isp)%spline%VALUE(nodes_on_spline(:,1))

      ! get the slope (m1)
      m1 = splines(isp)%spline%slope(nodes_on_spline(:,1))

      ! get the intercept (q1)
      q1 = nodes_on_spline_proj(:,2) - m1*nodes_on_spline_proj(:,1)

      ! find x coordinate of intersection point
      nodes_on_spline(:,1) = (q1-q2)/(m2-m1)

      ! project new point to spline and find y2 and slope m2
      nodes_on_spline_proj(:,1) = nodes_on_spline(:,1)
      nodes_on_spline_proj(:,2) = splines(isp)%spline%VALUE(nodes_on_spline(:,1))
      m1 = splines(isp)%spline%slope(nodes_on_spline(:,1))

    ENDSUBROUTINE project_from_opposite_vertex

    SUBROUTINE get_line_opposite_vertex
      !find coordinates of vertex opposite to each face
      DO i = 1, n
         DO j = 1, 3
            IF(COUNT(T(extFaces_bound(indices(i)),j) .EQ. Tb_bound_vtx) .EQ. 0) THEN
               opposite_vertex(i,:) = X(T(extFaces_bound(indices(i)),j),:)
            ENDIF
         ENDDO
      ENDDO

      CALL rotate(opposite_vertex)

      ! find slope and intercept of line passing through the vertex and the point to project
      m2 = (opposite_vertex(:,2) - nodes_on_spline(:,2))/(opposite_vertex(:,1) - nodes_on_spline(:,1))
      q2 = opposite_vertex(:,2) - m2*opposite_vertex(:,1)

    ENDSUBROUTINE get_line_opposite_vertex

    SUBROUTINE project_perp_to_spline
      ! project to spline and find y1
      nodes_on_spline_proj(:,1) = nodes_on_spline(:,1)
      nodes_on_spline_proj(:,2) = splines(isp)%spline%VALUE(nodes_on_spline(:,1))
      ! get the slope (m1)
      m1 = splines(isp)%spline%slope(nodes_on_spline(:,1))

      ! get the intercept (q1)
      q1 = nodes_on_spline_proj(:,2) - m1*nodes_on_spline_proj(:,1)

      ! get slope of perpendicular line (m2)
      m2 = -1/m1
      ! get intercept of perpendicular line (q2) passing through original point
      q2 = nodes_on_spline0(:,2)-m2*nodes_on_spline0(:,1)

      ! find x coordinate of intersection point
      nodes_on_spline(:,1) = (q2-q1)/(m1-m2)

      ! project new point to spline and find y2 and slope m2
      nodes_on_spline_proj(:,1) = nodes_on_spline(:,1)
      nodes_on_spline_proj(:,2) = splines(isp)%spline%VALUE(nodes_on_spline(:,1))
      m1 = splines(isp)%spline%slope(nodes_on_spline(:,1))

    ENDSUBROUTINE project_perp_to_spline

    SUBROUTINE rotate(nodes)
      REAL*8, INTENT(INOUT)            :: nodes(:,:)

      ! rotate the points so that the straight boundary aligns to the x axis
      xcenter = splines(isp)%xcenter
      ycenter = splines(isp)%ycenter

      nodes(:,1) = nodes(:,1) - xcenter
      nodes(:,2) = nodes(:,2) - ycenter

      RotMat = splines(isp)%RotMat

      nodes = TRANSPOSE(MATMUL(RotMat,TRANSPOSE(nodes)))

      nodes(:,1) = nodes(:,1) + xcenter
      nodes(:,2) = nodes(:,2) + ycenter

    ENDSUBROUTINE rotate

    SUBROUTINE anti_rotate(nodes)
      REAL*8, INTENT(INOUT)            :: nodes(:,:)

      ! rotate the points so that the straight boundary aligns to the x axis
      xcenter = splines(isp)%xcenter
      ycenter = splines(isp)%ycenter

      nodes(:,1) = nodes(:,1) - xcenter
      nodes(:,2) = nodes(:,2) - ycenter

      RotMat = splines(isp)%AntiRotMat

      nodes = TRANSPOSE(MATMUL(RotMat,TRANSPOSE(nodes)))

      nodes(:,1) = nodes(:,1) + xcenter
      nodes(:,2) = nodes(:,2) + ycenter

    ENDSUBROUTINE anti_rotate

  ENDSUBROUTINE round_edges

  SUBROUTINE blending_boundary(Z,m,porder,s,W)
    INTEGER, INTENT(IN)                 :: m, porder
    REAL*8, INTENT(IN)                  :: Z(:,:), s(:)
    REAL*8, INTENT(OUT)                 :: W(:)
    INTEGER                             :: i, ind
    REAL*8                              :: C(porder-1, porder-1), invC(porder-1, porder-1)
    REAL*8                              :: ZC(SIZE(Z,1),2)

    IF (m .EQ. 1) THEN
       ! First vertex ([0,0] in xi-eta)
       W = 1.0D0 - Z(:,1) - Z(:,2)
    ELSEIF (m .EQ. porder + 1) THEN
       ! Second vertex ([1,0] in xi-eta)
       W = Z(:,1)
    ELSEIF (m .EQ. 2*porder + 1) THEN
       ! Third vertex ([0,0] in xi-eta)
       W = Z(:,2)
    ELSE
       ! Edge nodes
       C = 1.0D0
       C(:,1) = s(2:SIZE(s)-1) * (1.0D0 - s(2:SIZE(s)-1))
       DO i = 2, porder-1
          C(:,i) = C(:,i-1) * s(2:SIZE(s)-1)
       ENDDO

       CALL invert_matrix(C,invC)

       IF (m < porder + 1) THEN
          ! First edge
          ZC = Z
          ind = m - 1
       ELSEIF (m .LT. 2*porder + 1) THEN
          ! Second edge
          ZC(:,1) = Z(:,2)
          ZC(:,2) = 1.0D0 - Z(:,1) - Z(:,2)
          ind = m - porder - 1
       ELSE
          ! Third edge
          ZC(:,1) = 1.0D0 - Z(:,1) - Z(:,2)
          ZC(:,2) = Z(:,1)
          ind = m - 2*porder - 1
       ENDIF

       W = 0.0D0
       DO i = 1, porder-1
          W = W + invC(i,ind) * ZC(:,1)**i
       ENDDO
       W = W * (1.0D0 - ZC(:,1) - ZC(:,2))
    ENDIF
  ENDSUBROUTINE blending_boundary

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
  END SUBROUTINE inverse_isop_transf

  SUBROUTINE inverse_linear_transformation(x,Xe,xieta)

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

  SUBROUTINE find_matches_int(a, b, indices)
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
    USE GMSH_io_module, ONLY: get_unit
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
    USE GMSH_io_module, ONLY: get_unit
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


    ! do
    !   read (unit_gmsh, '(a)', iostat = ios ) buffer
    !   if (ios /= 0) exit ! Exit loop if end of file is reached
    !   IF(buffer(1:9) .eq. '$Elements') exit
    ! end do
    !
    ! ! read one more line
    ! read (unit_gmsh, '(a)', iostat = ios ) buffer

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

    DEALLOCATE(Tb_gmsh,T_gmsh)
    WRITE(*,*) "********** Extended connectivity read **********"

  ENDSUBROUTINE read_extended_connectivity


  SUBROUTINE generate_msh_from_solution_mesh(mesh_name)
    USE GMSH_io_module, ONLY: get_unit

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
    USE gmsh

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)            :: mesh_name
    CHARACTER(LEN = 1024)               :: file_in, file_out
    TYPE(gmsh_t)                        :: gmsh

    file_in  = TRIM(ADJUSTL(mesh_name))// '.msh'
    file_out = TRIM(ADJUSTL(mesh_name))// '.mesh'

    CALL gmsh%initialize()
    CALL gmsh%OPEN(file_in)
    CALL gmsh%WRITE(file_out)
    CALL gmsh%finalize()

  ENDSUBROUTINE convert_msh2mesh

  SUBROUTINE convert_mesh2msh(mesh_name)
    USE gmsh

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)            :: mesh_name
    CHARACTER(LEN = 1024)               :: file_in, file_out
    TYPE(gmsh_t)                        :: gmsh

    file_in  = TRIM(ADJUSTL(mesh_name))// '.mesh'
    file_out = TRIM(ADJUSTL(mesh_name))// '.msh'

    CALL gmsh%initialize()
    CALL gmsh%OPEN(file_in)
    CALL gmsh%option%setNumber("Mesh.MshFileVersion", 2.2)
    CALL gmsh%WRITE(file_out)
    CALL gmsh%finalize()

  ENDSUBROUTINE convert_mesh2msh

  SUBROUTINE delete_file(filename)
    USE GMSH_io_module, ONLY: get_unit
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
    USE GMSH_io_module, ONLY: get_unit
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

  SUBROUTINE unique_1D(list_in, list_out)
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

  RECURSIVE SUBROUTINE quicksort_int(a)
    !! quicksort.f -*-f90-*-
    !! Author: t-nissie, some tweaks by 1AdAstra1
    !! License: GPLv3
    !! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    INTEGER, DIMENSION(:), INTENT(inout) :: a
    INTEGER ::  x, t
    INTEGER :: first = 1, last
    INTEGER i, j

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

  SUBROUTINE quicksort_real(a)
    !! quicksort.f -*-f90-*-
    !! Author: t-nissie, some tweaks by 1AdAstra1
    !! License: GPLv3
    !! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    REAL*8, DIMENSION(:), INTENT(inout) :: a
    !! The array to sort

    REAL*8 ::  x, t
    INTEGER :: first = 1, last
    INTEGER i, j

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

  END SUBROUTINE quicksort_real

  SUBROUTINE unique_stable(arrayin, uniqueArr)
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

  SUBROUTINE intersect_stable_int(a,b,c)
    INTEGER, INTENT(IN)                       :: a(:)
    INTEGER, INTENT(IN)                       :: b(:)
    INTEGER, INTENT(OUT), ALLOCATABLE         :: c(:)
    INTEGER, ALLOCATABLE                      :: temp(:)
    INTEGER                                   :: i, j, counter

    counter = 0
    DO i = 1, SIZE(a)
       DO j = 1, SIZE(b)
          IF(a(i) .EQ. b(j)) THEN
             counter = counter + 1
          ENDIF
       ENDDO
    ENDDO

    ALLOCATE(temp(counter))
    counter = 1

    DO i = 1, SIZE(a)
       DO j = 1, SIZE(b)
          IF(a(i) .EQ. b(j)) THEN
             temp(counter) = a(i)
             counter = counter + 1
          ENDIF
       ENDDO
    ENDDO

    ! c is allocated in here
    CALL unique_stable(temp, c)

    DEALLOCATE(temp)

  ENDSUBROUTINE intersect_stable_int

  ! SUBROUTINE move_vertices(mesh_init, refElPol_prec)
  !
  !   TYPE(Reference_element_type), INTENT(IN)  :: refElPol_prec
  !   TYPE(Mesh_type), INTENT(IN)               :: mesh_init
  !   INTEGER                                   :: elemType, nnodes, ndim, nextfaces_wol,nextFaces_prec_wol, i, j, k, iel, counter, correl1,correl2, n_elements_prec
  !   INTEGER, ALLOCATABLE                      :: extFaces(:,:),extFaces_prec(:,:)
  !   INTEGER, ALLOCATABLE                      :: indices(:), index_el_nextto(:), index_el_nextto_unique(:), index_el_nextto_ext(:,:), indices_nodes_correl(:)
  !   INTEGER                                   :: Face_nodes(refElPol_prec%Nfacenodes)
  !   INTEGER                                   :: correl(2,2), iface(2)
  !   REAL*8                                    :: Xe_prec(mesh_init%Nnodesperelem,2), X_vertex_nodes(2,2), bc(3), A(3,3), b(3), xieta_vertex1(1,2), xieta_vertex2(1,2), xieta_node1(1,2)
  !   REAL*8                                    :: shapeFunctions1(refElPol_prec%Nnodes2D,SIZE(xieta_vertex1,1),3),shapeFunctions2(refElPol_prec%Nnodes2D,SIZE(xieta_vertex2,1),3)
  !   REAL*8                                    :: tol=1e-10
  !   REAL*8, ALLOCATABLE                       :: temp(:,:)
  !
  !
  !   elemType        = Mesh%elemType
  !   nnodes          = Mesh%Nnodesperelem
  !   ndim            = Mesh%ndim
  !   n_elements_prec = mesh_init%Nelems
  !   correl          = 0
  !   correl1         = 0
  !   correl2         = 0
  !
  !   ALLOCATE(temp(1,SIZE(X_vertex_nodes,2)))
  !   ! overkill allocation
  !   ALLOCATE(index_el_nextto(n_elements_prec/5))
  !   ALLOCATE(index_el_nextto_ext(n_elements_prec/10,2))
  !
  !   index_el_nextto = 0
  !   index_el_nextto_ext = 0
  !
  !   ! six is twice the number of vertices on a triangle
  !   ALLOCATE(indices_nodes_correl(6))
  !
  !   ! number of ext faces with out limiter
  !   nextfaces_wol = COUNT(Mesh%boundaryFlag .NE. 7) ! count external faces not on limiter
  !
  !   ALLOCATE(extFaces(nextfaces_wol,2))
  !
  !   counter = 1
  !   DO i = 1, Mesh%NextFaces
  !      IF (Mesh%boundaryFlag(i) .NE. 7) THEN
  !         extfaces(counter,:) = Mesh%extfaces(i,:)
  !         counter = counter + 1
  !      ENDIF
  !   ENDDO
  !
  !   nextFaces_prec_wol = COUNT(mesh_init%boundaryFlag .NE. 7) ! count external faces not on limiter
  !   ALLOCATE(extFaces_prec(nextFaces_prec_wol,2))
  !
  !   counter = 1
  !   DO i = 1, mesh_init%NextFaces
  !      IF (mesh_init%boundaryFlag(i) .NE. 7) THEN
  !         extFaces_prec(counter,:) = mesh_init%extfaces(i,:)
  !         counter = counter + 1
  !      ENDIF
  !   ENDDO
  !
  !   DO i = 1, nextfaces_wol
  !      Face_nodes = refElPol_prec%Face_nodes(extfaces(i,2),:)
  !      ! physical coord of first vertex
  !      X_vertex_nodes(1,:) = Mesh%X(Mesh%T(extfaces(i,1),Face_nodes(1)),:)
  !      ! physical coord of last vertex
  !      X_vertex_nodes(2,:) = Mesh%X(Mesh%T(extfaces(i,1),Face_nodes(SIZE(Face_nodes))),:)
  !      tol = 1e-10
  !      DO WHILE(ANY(correl .EQ. 0)) ! find nodes and corresponding elements on old mesh
  !         ! loop on the 2 vertices
  !         DO j  = 1, SIZE(X_vertex_nodes,1)
  !            ! loop on old elements
  !            DO iel = 1, n_elements_prec
  !
  !               Xe_prec = mesh_init%X(mesh_init%T(iel,:),:)
  !               b(1) = X_vertex_nodes(j,1)
  !               b(2) = X_vertex_nodes(j,2)
  !               b(3) = 1.
  !
  !               A(1,:) = Xe_prec(1:3,1)
  !               A(2,:) = Xe_prec(1:3,2)
  !               A(3,:) = 1.
  !
  !               ! find IF the element is inside the triangle
  !               CALL solve_linear_system_sing(A,b,bc)
  !
  !               ! The point is inside
  !               IF(ALL(bc .GE. -tol ) .AND. (ALL(bc .LE. (1 + tol)))) THEN
  !                  correl1 = iel
  !                  correl2 = j
  !               ENDIF
  !
  !               IF((correl1 .NE. 0) .AND. (correl2 .NE. 0)) THEN
  !                  IF(correl(j,1) .NE. 0) THEN
  !                     correl1 = 0
  !                     correl2 = 0
  !                     EXIT
  !                  ELSE
  !                     correl(j,1) = correl2
  !                     correl(j,2) = correl1
  !                  ENDIF
  !               ENDIF
  !               correl1 = 0
  !               correl2 = 0
  !            ENDDO
  !         ENDDO
  !         tol = tol * 10
  !      ENDDO
  !
  !      iface = 0
  !
  !      ! find face (1, 2 or 3) on old mesh corresponding to the two vertices
  !      DO j = 1, SIZE(correl,1)
  !         DO k = 1, SIZE(ExtFaces_prec,1)
  !            IF(extFaces_prec(k,1) .EQ. correl(j,2)) THEN
  !               iface(j) = extFaces_prec(k,2)
  !            ENDIF
  !         ENDDO
  !      ENDDO
  !
  !      ! fix node on element with no external face
  !      IF(ANY(extFaces_prec(:,1) .EQ. correl(1,2))) THEN
  !         correl(1,2) = correl(1,2)
  !      ELSE
  !         ! IF last vertex is ok
  !         IF(ANY(extFaces_prec(:,1) .EQ. correl(2,2))) THEN
  !            correl(1,2) = correl(2,2)
  !            iface(1)  = iface(2)
  !         ENDIF
  !      ENDIF
  !
  !      IF(ANY(extFaces_prec(:,1) .EQ. correl(2,2))) THEN
  !         correl(2,2) = correl(2,2)
  !      ELSE
  !         ! IF first vertex is ok
  !         IF(ANY(extFaces_prec(:,1) .EQ. correl(1,2))) THEN
  !            correl(2,2) = correl(1,2)
  !            iface(2)  = iface(1)
  !         ENDIF
  !      ENDIF
  !
  !      ! IF none of the 2 vertices are on an external element
  !      IF(ANY(iface .EQ. 0)) THEN
  !         ! indices is allocated here
  !         IF(ALLOCATED(indices)) DEALLOCATE(indices)
  !         CALL find_matches_int(iface, 0, indices)
  !
  !         ! indices of the nodes of the two elements found
  !         indices_nodes_correl(1:3) =  mesh_init%T(correl(1,2),1:3)
  !         indices_nodes_correl(4:6) =  mesh_init%T(correl(2,2),1:3)
  !
  !
  !         counter = 1
  !         DO j = 1, n_elements_prec
  !            ! 6 is twice the number of vertices (length of indices_nodes_correl)
  !            DO k = 1, 6
  !               IF(ANY(mesh_init%T(j,:) .EQ. indices_nodes_correl(k))) THEN
  !                  index_el_nextto(counter) = j
  !                  counter = counter + 1
  !               ENDIF
  !            ENDDO
  !         ENDDO
  !
  !         ! list of elements close to the two individuated (index_el_nextto_unique is allocated in here)
  !         CALL unique_1D(index_el_nextto, index_el_nextto_unique)
  !
  !         ! find which one of those are external
  !         counter = 1
  !         ! start from 2 cause the first one is always a 0
  !         DO j = 2, SIZE(index_el_nextto_unique)
  !            DO k = 1, nextFaces_prec_wol
  !               IF(index_el_nextto_unique(j) .EQ. extFaces_prec(k,1)) THEN
  !                  index_el_nextto_ext(counter,1) = index_el_nextto_unique(j)
  !                  index_el_nextto_ext(counter,2) = extFaces_prec(k,2)
  !                  counter = counter + 1
  !               ENDIF
  !            ENDDO
  !         ENDDO
  !         tol = 1e-15
  !         DO j = 1, SIZE(indices)
  !            DO k = 1, SIZE(index_el_nextto_ext,1)
  !               ! you reached the end of the array
  !               IF(ANY(index_el_nextto_ext(k,:) .EQ. 0)) EXIT
  !               ! projection on element
  !               CALL inverse_isop_transf(RESHAPE(X_vertex_nodes(indices(j),:), [1,2]), mesh_init%X(mesh_init%T(index_el_nextto_ext(k,1),:),:), refElPol_prec, xieta_node1)
  !
  !               IF(index_el_nextto_ext(k,2) .EQ. 1) THEN
  !                  IF((xieta_node1(1,1) .GT. -1-tol) .AND. (xieta_node1(1,1) .LT. 1+tol)) THEN
  !                     correl(indices(j),2) = index_el_nextto_ext(k,1)
  !                     iface(indices(j)) = index_el_nextto_ext(k,2)
  !                     EXIT
  !                  ENDIF
  !               ELSEIF(index_el_nextto_ext(k,2) .EQ. 2) THEN
  !                  IF(((xieta_node1(1,1) .GT. -1-tol) .AND. (xieta_node1(1,1) .LT. 1+tol))  .OR. ((xieta_node1(1,2) .GT. -1-tol) .AND. (xieta_node1(1,2) .LT. 1+tol))) THEN
  !                     correl(indices(j),2) = index_el_nextto_ext(k,1)
  !                     iface(indices(j)) = index_el_nextto_ext(k,2)
  !                     EXIT
  !                  ENDIF
  !               ELSEIF(index_el_nextto_ext(k,2) .EQ. 3) THEN
  !                  IF((xieta_node1(1,2) .GT. -1-tol) .AND. (xieta_node1(1,2) .LT. 1+tol)) THEN
  !                     correl(indices(j),2) = index_el_nextto_ext(k,1)
  !                     iface(indices(j)) = index_el_nextto_ext(k,2)
  !                     EXIT
  !                  ENDIF
  !               ENDIF
  !            ENDDO
  !         ENDDO
  !      ENDIF
  !
  !      ! calculate xieta for both vertices
  !      ! find xieta of first and last node of the face that are on the rounded border
  !      CALL inverse_isop_transf(RESHAPE(X_vertex_nodes(1,:),[1, 2]), mesh_init%X(mesh_init%T(correl(1,2),:),:), refElPol_prec, xieta_vertex1)
  !      ! calculate xieta for both vertices
  !      CALL inverse_isop_transf(RESHAPE(X_vertex_nodes(2,:),[1, 2]), mesh_init%X(mesh_init%T(correl(2,2),:),:), refElPol_prec, xieta_vertex2)
  !
  !      ! check on which face it is and impose to stay on the extreme
  !      IF (iface(1) .EQ. 1) THEN
  !         xieta_vertex1(1,2) = -1.0
  !      ELSE IF (iface(1) .EQ. 2) THEN
  !         ! linear projection
  !         xieta_vertex1(1, 1) = -(xieta_vertex1(1, 2) - xieta_vertex1(1, 1)) / 2.0
  !         xieta_vertex1(1, 2) = (xieta_vertex1(1, 2) - xieta_vertex1(1, 1)) / 2.0
  !      ELSE IF (iface(1) .EQ. 3) THEN
  !         xieta_vertex1(1,1) = -1.0
  !      END IF
  !
  !      IF (iface(2) .EQ. 1) THEN
  !         xieta_vertex2(1,2) = -1.0
  !      ELSE IF (iface(2) .EQ. 2) THEN
  !         ! linear projection
  !         xieta_vertex2(1, 1) = -(xieta_vertex2(1, 2) - xieta_vertex2(1,1)) / 2.0
  !         xieta_vertex2(1, 2) =  (xieta_vertex2(1, 2) - xieta_vertex2(1,1)) / 2.0
  !      ELSE IF (iface(2) .EQ. 3) THEN
  !         xieta_vertex2(1,1) = -1.0
  !      END IF
  !
  !      CALL compute_shape_functions_at_points(refElPol_prec, xieta_vertex1, shapeFunctions1)
  !      CALL compute_shape_functions_at_points(refElPol_prec, xieta_vertex2, shapeFunctions2)
  !
  !      Mesh%X(Mesh%T(ExtFaces(i,1),Face_nodes(1)),1) = dot_PRODUCT(shapeFunctions1(:,1,1),mesh_init%X(mesh_init%T(correl(1,2),:),1))
  !      Mesh%X(Mesh%T(ExtFaces(i,1),Face_nodes(1)),2) = dot_PRODUCT(shapeFunctions1(:,1,1),mesh_init%X(mesh_init%T(correl(1,2),:),2))
  !
  !      Mesh%X(Mesh%T(ExtFaces(i,1),Face_nodes(SIZE(Face_nodes))),1) = dot_PRODUCT(shapeFunctions2(:,1,1),mesh_init%X(mesh_init%T(correl(2,2),:),1))
  !      Mesh%X(Mesh%T(ExtFaces(i,1),Face_nodes(SIZE(Face_nodes))),2) = dot_PRODUCT(shapeFunctions2(:,1,1),mesh_init%X(mesh_init%T(correl(2,2),:),2))
  !
  !      correl = 0
  !   ENDDO
  !
  !   DEALLOCATE(temp)
  !   DEALLOCATE(index_el_nextto)
  !   DEALLOCATE(index_el_nextto_ext)
  !   DEALLOCATE(indices_nodes_correl)
  !   IF(ALLOCATED(indices)) DEALLOCATE(indices)
  ! ENDSUBROUTINE move_vertices

  ! SUBROUTINE round_edges_giacomo(mesh_init, refElPol_prec)
  !
  !   TYPE(Reference_element_type), INTENT(IN)  :: refElPol_prec
  !   TYPE(Reference_element_type)              :: refElPol_local
  !   TYPE(Mesh_type), INTENT(IN)               :: mesh_init
  !   INTEGER                                   :: elemType, nnodes, ndim, nextfaces_wol,nextFaces_prec_wol, i, j, k, iel, counter, correl1,correl2, n_elements_prec
  !   INTEGER, ALLOCATABLE                      :: extFaces(:,:),extFaces_prec(:,:)
  !   INTEGER, ALLOCATABLE                      :: indices(:), index_el_nextto(:), index_el_nextto_unique(:), index_el_nextto_ext(:,:), indices_nodes_correl(:)
  !   INTEGER                                   :: Face_nodes(refElPol_prec%Nfacenodes),int_face_nodes(refElPol_prec%Nfacenodes-2)
  !   REAL*8                                    :: X_int_nodes(refElPol_prec%Nfacenodes-2,2)
  !   REAL*8, ALLOCATABLE                       :: X(:,:)
  !   INTEGER                                   :: correl(2,2), iface(2)
  !   REAL*8                                    :: Xe_prec(mesh_init%Nnodesperelem,2), X_vertex_nodes(2,2), bc(3), A(3,3), b(3), xieta_vertex1(1,2), xieta_vertex2(1,2), xieta_node1(1,2)
  !   REAL*8                                    :: shapeFunctions1(refElPol_prec%Nnodes2D,SIZE(xieta_vertex1,1),3),shapeFunctions2(refElPol_prec%Nnodes2D,SIZE(xieta_vertex2,1),3)
  !   REAL*8                                    :: tol=1e-10
  !   REAL*8                                    :: X_new(refElPol_prec%Nfacenodes-2,2)
  !
  !   REAL*8, ALLOCATABLE                       :: temp(:,:)
  !
  !
  !   CALL create_reference_element(refElPol_local,2,4, verbose = 1)
  !   elemType        = Mesh%elemType
  !   nnodes          = Mesh%Nnodesperelem
  !   ndim            = Mesh%ndim
  !   n_elements_prec = mesh_init%Nelems
  !   correl          = 0
  !   correl1         = 0
  !   correl2         = 0
  !
  !   ALLOCATE(temp(1,SIZE(X_vertex_nodes,2)))
  !   ! overkill allocation
  !   ALLOCATE(index_el_nextto(n_elements_prec/2))
  !   ALLOCATE(index_el_nextto_ext(n_elements_prec/2,2))
  !   ! six is twice the number of vertices on a triangle
  !   ALLOCATE(indices_nodes_correl(6))
  !   ALLOCATE(indices(2))
  !
  !   ! number of ext faces with out limiter
  !   nextfaces_wol = COUNT(Mesh%boundaryFlag .ne. 7) ! count external faces not on limiter
  !
  !   ALLOCATE(extFaces(nextfaces_wol,2))
  !   extFaces = 0
  !
  !   counter = 1
  !   DO i = 1, Mesh%NextFaces
  !     IF (Mesh%boundaryFlag(i) .ne. 7) THEN
  !       extFaces(counter,:) = Mesh%extfaces(i,:)
  !       counter = counter + 1
  !     ENDIF
  !   ENDDO
  !
  !
  !   nextFaces_prec_wol = COUNT(mesh_init%boundaryFlag .ne. 7) ! count external faces not on limiter
  !
  !   ALLOCATE(extFaces_prec(nextFaces_prec_wol,2))
  !   extFaces_prec = 0
  !
  !   counter = 1
  !   DO i = 1, mesh_init%NextFaces
  !     IF (mesh_init%boundaryFlag(i) .ne. 7) THEN
  !       extFaces_prec(counter,:) = mesh_init%extfaces(i,:)
  !       counter = counter + 1
  !     ENDIF
  !   ENDDO
  !
  !
  !   ! loop on external faces
  !   DO i = 1, nextfaces_wol
  !     Face_nodes = refElPol_local%Face_nodes(extFaces(i,2),:)! all the nodes on the external face
  !     int_face_nodes = refElPol_local%Face_Nodes(extFaces(i,2),2:SIZE(refElPol_local%Face_Nodes,2)-1) !central nodes on the external face
  !     X_int_nodes = Mesh%X(Mesh%T(ExtFaces(i,1),int_face_nodes),:) !physical coords of central nodes
  !     X_vertex_nodes(1,:) = Mesh%X(Mesh%T(extfaces(i,1),Face_nodes(1)),:) !physical coords of the first vertex
  !     X_vertex_nodes(2,:) = Mesh%X(Mesh%T(extfaces(i,1),Face_nodes(SIZE(Face_nodes))),:) !physical coords of the last vertex
  !     CALL get_coords_rounded_mesh(X_int_nodes,mesh_init%X,mesh_init%T,refElPol_local,X_vertex_nodes,extFaces_prec, X_new) ! find the new coords of the nodes on the rounded borders
  !     Mesh%X(Mesh%T(ExtFaces(i,1),int_face_nodes),:) = X_new(:,:)
  !     CALL new_points_2D(Mesh%X(Mesh%T(ExtFaces(i,1),:),:),refElPol_prec%Ndeg,refElPol_local%coord2d, X)
  !     Mesh%X(Mesh%T(ExtFaces(i,1),:),:) = X
  !   ENDDO
  !
  !   CALL free_reference_element_pol(refElPol_local)
  !   DEALLOCATE(temp)
  !   DEALLOCATE(index_el_nextto)
  !   DEALLOCATE(index_el_nextto_ext)
  !   DEALLOCATE(indices_nodes_correl)
  !   DEALLOCATE(indices)
  !   DEALLOCATE(X)
  ! ENDSUBROUTINE round_edges_giacomo
  !
  ! SUBROUTINE get_coords_rounded_mesh(x_central, X, T, refEl, x_vertex_nodes, extFaces_prec, X_new)
  !   TYPE(Reference_element_type), INTENT(IN)    :: refEl
  !   REAL*8, INTENT(IN)                          :: x_central(:,:), X(:,:), x_vertex_nodes(:,:)
  !   INTEGER, INTENT(IN)                         :: T(:,:), extFaces_prec(:,:)
  !   REAL*8, INTENT(OUT)                         :: X_new(:,:)
  !   REAL*8                                      :: X_all_face(refEl%Nnodes1D,2), A(3,3), b(3), bc(3)
  !   REAL*8                                      :: xieta_temp(refEl%Nnodes1D-2,2),xieta(refEl%Nnodes1D-2,2), xieta_vertex1(1,2), xieta_vertex2(1,2), xieta_node1(1,2), xieta_node2(1,2)
  !   REAL*8                                      :: xieta_check(1,2), xieta_check1(1,2), xieta_check2(1,2), xieta_check_pre(1,2), xieta1(1,2), xieta2(1,2)
  !   INTEGER                                     :: correl(SIZE(x_central,1),1), correl3(size(x_all_face,1),2)
  !   INTEGER, ALLOCATABLE                        :: correl_unique(:), correl3_unique(:), face_nodes(:,:), iface_vec(:,:)
  !   INTEGER, ALLOCATABLE                        :: index_el_nextto(:), index_el_nextto_unique(:), index_el_nextto_ext(:,:), index_nodes_correl(:,:), index_nodes_correl_unique(:), node_check(:), &
  !                                                  indices(:), index_common_nodes(:), index_common_nodes_check(:), index_common_nodes_check_unique(:), index_common_nodes_check1(:), index_common_nodes_check2(:), &
  !                                                  index_common_nodes_check3(:)
  !   REAL*8                                      :: shapeFunctions(refEl%Nnodes2D,refEl%Nnodes1D-2,3), shapeFunctions1(refEl%Nnodes2D,1,3), shapeFunctions2(refEl%Nnodes2D,1,3)
  !   REAL*8, ALLOCATABLE                         :: xieta_vertex_middle1(:,:), xieta_vertex_middle2(:,:), temp(:,:), xieta_vertex_first(:,:), xieta_vertex_last(:,:),&
  !                                                  xieta_vertex_first_outNode(:,:),xieta_vertex_last_outNode(:,:)
  !   INTEGER                                     :: n, n_points, n_elements, correl1, correl2, correl4, ip, iel, iface, counter, index_middleNode_oldMesh,  i, j
  !   REAL*8                                      :: tol
  !
  !
  !   ALLOCATE(xieta_vertex_first(SIZE(correl3,1),2))
  !   ALLOCATE(xieta_vertex_last(SIZE(correl3,1),2))
  !   ALLOCATE(xieta_vertex_first_outNode(SIZE(correl3,1),2))
  !   ALLOCATE(xieta_vertex_last_outNode(SIZE(correl3,1),2))
  !
  !   X_new = 0.
  !   X_all_face = 0.
  !   correl = 0
  !   correl1 = 0
  !   correl2 = 0
  !   correl3 = 0
  !   correl4 = 0
  !   tol = 1e-10
  !   n_elements = SIZE(T,1)
  !   n_points = SIZE(x_central,1)
  !
  !   X_all_face(1,:) = x_vertex_nodes(1,:)
  !   X_all_face(2:2+SIZE(x_central,1)-1,:) = x_central
  !   X_all_face(5,:) = x_vertex_nodes(2,:)
  !
  !
  !   DO WHILE (ANY(correl(:,1) .eq. 0))
  !     DO ip = 1, n_points
  !       DO iel = 1, n_elements
  !
  !         b(1) = x_central(ip,1)
  !         b(2) = x_central(ip,2)
  !         b(3) = 1.
  !
  !         A(1,:) = X(T(iel,1:3),1)
  !         A(2,:) = X(T(iel,1:3),2)
  !         A(3,:) = 1.
  !
  !         ! find IF the element is inside the triangle
  !         CALL solve_linear_system_sing(A,b,bc)
  !
  !
  !         IF (ALL(bc .ge. -tol) .AND. ALL(bc .le. 1.d0 + tol)) THEN
  !           correl4 = iel
  !         ENDIF
  !
  !         IF (correl4 .ne. 0) THEN
  !           IF (correl(ip,1) .ne. 0) THEN
  !             correl4 = 0
  !             EXIT
  !           ELSE
  !             correl(ip,1) = correl4 ! corresponding element
  !           END IF
  !         END IF
  !         correl4 = 0
  !       END DO
  !     END DO
  !     tol = tol * 10.
  !   END DO
  !
  !   CALL unique_stable(correl(:,1), correl_unique)
  !
  !   tol = 1.e-10
  !
  !   ALLOCATE(index_el_nextto(n_elements/2))
  !   ALLOCATE(index_el_nextto_ext(n_elements/2,2))
  !
  !   ! Element nodes are inside only one element of the previous mesh
  !   IF (SIZE(correl_unique) .eq. 1) THEN
  !
  !     iface = 0
  !     counter = 0
  !     !  what if there are two external faces??
  !     DO j = 1, size(ExtFaces_prec,1)
  !       IF (extFaces_prec(j,1) .eq. correl_unique(1)) THEN
  !         iface = extFaces_prec(j,2)
  !         counter = counter + 1
  !       END IF
  !     ENDDO
  !
  !     IF (iface .eq. 0) THEN
  !
  !       ALLOCATE(index_nodes_correl(1,3))
  !
  !       index_el_nextto = 0
  !       index_el_nextto_ext = 0
  !
  !       counter = 1
  !       index_nodes_correl(1,:) = T(correl_unique(1), 1:3)
  !       DO i = 1, n_elements
  !         DO j = 1, SIZE(index_nodes_correl,2)
  !           IF (ANY(T(i,1:3) .eq. index_nodes_correl(1,j))) THEN
  !             index_el_nextTo(counter) = i
  !             counter = counter + 1
  !           END IF
  !         END DO
  !       END DO
  !
  !       ! list of elements close to the two individuated (index_el_nextto_unique is allocated in here)
  !       CALL unique_1D(index_el_nextto, index_el_nextto_unique)
  !
  !       ! find which one of those are external
  !       counter = 1
  !       ! start from 2 cause the first one is always a 0
  !       DO i = 2, SIZE(index_el_nextto_unique)
  !         DO j = 1, SIZE(extFaces_prec,1)
  !           IF(index_el_nextto_unique(i) .eq. extFaces_prec(j,1)) THEN
  !             index_el_nextto_ext(counter,1) = index_el_nextto_unique(i)
  !             index_el_nextto_ext(counter,2) = extFaces_prec(j,2)
  !             counter = counter + 1
  !           ENDIF
  !         ENDDO
  !       ENDDO
  !
  !       DO i = 1, SIZE(index_el_nextTo_ext,1)
  !         CALL inverse_isop_transf(x_central,X(T(index_el_nextTo_ext(i,1),:),:), refEl, xieta_node1)
  !         IF ((xieta_node1(1,1) .le. 1) .AND. (xieta_node1(1,1) .ge. -1)) THEN
  !           correl_unique(1) = index_el_nextTo_ext(i,1)
  !           iface = index_el_nextTo_ext(i,2)
  !           EXIT
  !         END IF
  !       END DO
  !     END IF
  !
  !     ! find xieta of first and last node of the face that are on the rounded border
  !     CALL inverse_isop_transf(RESHAPE(x_vertex_nodes(1,:),[1,2]),X(T(correl_unique(1),:),:)+1e-14,refEl, xieta_vertex1)
  !     CALL inverse_isop_transf(RESHAPE(x_vertex_nodes(2,:),[1,2]),X(T(correl_unique(1),:),:)+1e-14,refEl, xieta_vertex2)
  !
  !     IF(counter .gt. 1) THEN
  !       IF((xieta_vertex1(1,1) .lt. 0) .AND. (xieta_vertex1(1,2) .lt. 0) .AND. (xieta_vertex2(1,1) .gt. 0) .AND. (xieta_vertex2(1,2) .lt. 0)) THEN
  !         iface = 1
  !       ELSEIF((xieta_vertex1(1,1) .gt. 0) .AND. (xieta_vertex1(1,2) .lt. 0) .AND. (xieta_vertex2(1,1) .lt. 0) .AND. (xieta_vertex2(1,2) .gt. 0)) THEN
  !         iface = 2
  !       ELSEIF((xieta_vertex1(1,1) .lt. 0) .AND. (xieta_vertex1(1,2) .gt. 0) .AND. (xieta_vertex2(1,1) .lt. 0) .AND. (xieta_vertex2(1,2) .lt. 0)) THEN
  !         iface = 3
  !       ENDIF
  !     ENDIF
  !
  !     xieta_temp(:,1) = refEl%coord1d(2:refEl%Nnodes1D-1)
  !     xieta_temp(:,2) = refEl%coord1d(2:refEl%Nnodes1D-1)
  !     xieta(:,1) = xieta_vertex1(1,1) + ((xieta_temp(:,1)+1)/2) * (xieta_vertex2(1,1)-xieta_vertex1(1,1))
  !     xieta(:,2) = xieta_vertex1(1,2) + ((xieta_temp(:,2)+1)/2) * (xieta_vertex2(1,2)-xieta_vertex1(1,2))
  !
  !     IF (iface == 1) THEN
  !       xieta(:,2) = -1.d0
  !     ELSEIF (iface == 2) THEN
  !       xieta(:,1) = -(xieta(:,2) - xieta(:,1))/2
  !       xieta(:,2) = (xieta(:,2) - xieta(:,1))/2
  !     ELSEIF (iface == 3) THEN
  !       xieta(:,1) = -1.d0
  !     END IF
  !
  !     CALL compute_shape_functions_at_points(refElPol, xieta, shapeFunctions)
  !     X_new = MATMUL(TRANSPOSE(shapeFunctions(:,:,1)), X(T(correl_unique(1),:),:))
  !
  !   ELSEIF (SIZE(correl_unique) .gt. 1) THEN
  !
  !       !!  external element nodes are inside two elements of the previous mesh
  !       DO WHILE (ANY(correl3(:,1) .eq. 0)) !find nodes and corresponding elements
  !           DO ip = 1, SIZE(x_all_face,1)
  !               DO iel = 1, n_elements
  !
  !                 b(1) = x_all_face(ip,1)
  !                 b(2) = x_all_face(ip,2)
  !                 b(3) = 1.
  !
  !                 A(1,:) = X(T(iel,1:3),1)
  !                 A(2,:) = X(T(iel,1:3),2)
  !                 A(3,:) = 1.
  !
  !                 ! find IF the element is inside the triangle
  !                 CALL solve_linear_system_sing(A,b,bc)
  !
  !                 IF (ALL(bc .ge. -tol) .AND. ALL(bc .le. 1.d0 + tol)) THEN
  !                    correl1 = iel
  !                    correl2 = ip
  !                 ENDIF
  !
  !                 IF((correl1 .ne. 0) .and. (correl2 .ne. 0)) THEN
  !                   IF (correl3(ip,1) .ne. 0 ) THEN
  !                     correl1 = 0
  !                     correl2 = 0
  !                     EXIT
  !                   ELSE
  !                     ! [nodes index, corresponding element]
  !                     correl3(ip,1) = correl2
  !                     correl3(ip,2) = correl1
  !                   ENDIF
  !                 ENDIF
  !                 correl1 = 0
  !                 correl2 = 0
  !               ENDDO
  !           ENDDO
  !           tol = tol*10
  !       ENDDO
  !
  !       ! Fix the first vertex
  !       ! Find to which element the first vertex belongs to
  !       DO i = 2, SIZE(correl3,1)
  !         IF (ANY(ExtFaces_prec(:, 1) .eq. correl3(1, 2))) THEN
  !           EXIT ! exit here
  !         ELSE
  !           IF (ANY(ExtFaces_prec(:, 1) .eq. correl3(i, 2))) THEN
  !             CALL inverse_isop_transf(RESHAPE(x_vertex_nodes(1,:),[1,2]), X(T(correl3(i, 2), :),:), refEl, xieta_check1)
  !             IF ((xieta_check1(1, 1) .ge. -1) .AND. (xieta_check1(1, 1) .le. 1)) THEN
  !               correl3(1, 2) = correl3(i, 2)
  !               EXIT
  !             ELSE IF ((xieta_check1(1, 2) .ge. -1) .AND. (xieta_check1(1, 2) .le. 1)) THEN
  !               correl3(1, 2) = correl3(i, 2)
  !               EXIT
  !             END IF
  !           END IF
  !         END IF
  !       END DO
  !
  !       DO i = 1, SIZE(T,1)
  !         DO j = 1, 3
  !           IF(sqrt((x_vertex_nodes(1,1)-X(T(i,j),1))**2+(x_vertex_nodes(1,2)-X(T(i,j),2))**2) .lt. 1e-12) THEN
  !             correl3(1,2) = i
  !           ENDIF
  !         ENDDO
  !       ENDDO
  !
  !
  !       ! Fix the last vertex
  !       ! Find to which element the last vertex belongs to
  !       DO i = SIZE(correl3, 1) - 1, 1, -1
  !         IF (ANY(ExtFaces_prec(:, 1) .eq. correl3(SIZE(correl3, 1), 2))) THEN
  !           EXIT ! exit here
  !         ELSE
  !           IF (ANY(ExtFaces_prec(:, 1) .eq. correl3(i, 2))) THEN
  !             CALL inverse_isop_transf(RESHAPE(x_vertex_nodes(SIZE(x_vertex_nodes, 1),:),[1,2]), X(T(correl3(i, 2), :),:), refEl, xieta_check2)
  !             IF ((xieta_check2(1, 1) .ge. -1) .AND. (xieta_check2(1, 1) .le. 1)) THEN
  !               correl3(SIZE(correl3, 1), 2) = correl3(i, 2)
  !               EXIT
  !             ELSE IF ((xieta_check2(1, 2) .ge. -1.) .AND. (xieta_check2(1, 2) .le. 1)) THEN
  !               correl3(SIZE(correl3, 1), 2) = correl3(i, 2)
  !               EXIT
  !             END IF
  !           END IF
  !         END IF
  !       END DO
  !
  !       DO i = 1, SIZE(T,1)
  !         DO j = 1, 3
  !           IF(sqrt((x_vertex_nodes(1,1)-X(T(i,j),1))**2+(x_vertex_nodes(1,2)-X(T(i,j),2))**2) .lt. 1e-12) THEN
  !             correl3(SIZE(correl3, 1),2) = i
  !           ENDIF
  !         ENDDO
  !       ENDDO
  !
  !       ! Fix the internal nodes
  !       ! Find to which element the internal nodes belong to
  !       DO i = 2, SIZE(correl3,1)-1
  !         ! if the element is on external surface, make a check IF it is correct or wrong
  !         IF (ANY(ExtFaces_prec(:,1) .eq. correl3(i,2))) THEN
  !           ! check IF the element is completely wrong and far away (no nodes found between consecutive elts)
  !           CALL intersect_stable_int(T(correl3(i,2),1:3),T(correl3(i-1,2),1:3),node_check)
  !
  !           ! The elements are not neighbours
  !           IF(SIZE(node_check) .eq. 0) THEN
  !             ! check on previous element
  !             CALL inverse_isop_transf(RESHAPE(x_all_face(i,:),[1,2]),X(T(correl3(i-1,2),:),:),refEl, xieta_check)
  !             IF((xieta_check(1,1) .ge. -1) .AND. (xieta_check(1,1) .le. 1)) THEN
  !               correl3(i,2) = correl3(i-1,2)
  !             ELSE
  !               ! check on next elements
  !               DO j = i, SIZE(correl3,1)
  !                  CALL inverse_isop_transf(RESHAPE(x_all_face(i,:),[1,2]),X(T(correl3(j,2),:),:),refEl, xieta_check)
  !                  IF((xieta_check(1,1) .ge. -1) .AND. (xieta_check(1,1) .le. 1)) THEN
  !                    correl3(i,2) = correl3(j,2)
  !                    EXIT
  !                  ELSE
  !                    IF((xieta_check(1,2) .ge. -1) .AND. (xieta_check(1,2) .le. 1)) THEN
  !                     correl3(i,2) = correl3(j,2)
  !                     EXIT
  !                    ENDIF
  !                  ENDIF
  !                ENDDO
  !              ENDIF
  !           ELSE
  !              ! node_check is not empty --> elements are neighbours
  !              ! check on same element, IF it is fine keep it
  !               CALL inverse_isop_transf(RESHAPE(x_all_face(i,:),[1,2]),X(T(correl3(i,2),:),:),refEl, xieta_check_pre);
  !               IF ((xieta_check_pre(1,1) .ge. -1) .AND. (xieta_check_pre(1,1) .le. 1)) THEN
  !                   CYCLE
  !               ELSE
  !                   IF ((xieta_check_pre(1,2) .ge. -1) .AND. (xieta_check_pre(1,2) .ge. 1)) THEN
  !                     CYCLE
  !                   ENDIF
  !               ENDIF
  !               ! IF the nodes are not on the same element, start to check all others
  !               DO j=1,SIZE(correl3,1)
  !                   ! check on next elements IF is not empty the node between elts
  !                   CALL inverse_isop_transf(RESHAPE(x_all_face(i,:),[1,2]),X(T(correl3(j,2),:),:),refEl, xieta_check)
  !                   IF ((xieta_check(1,1) .ge. -1) .AND. (xieta_check(1,1) .le. 1)) THEN
  !                       correl3(i,2) = correl3(j,2)
  !                       EXIT
  !                   ELSE
  !                       IF ((xieta_check(1,2) .ge. -1) .AND. (xieta_check(1,2) .le. 1)) THEN
  !                           correl3(i,2) = correl3(j,2)
  !                           EXIT
  !                       ENDIF
  !                   ENDIF
  !               ENDDO
  !           ENDIF
  !
  !         ELSE
  !             ! The element is not on the external surface
  !             ! Check on precedent element, IF it is fine, keep it
  !             CALL inverse_isop_transf(RESHAPE(x_all_face(i,:),[1,2]), X(T(correl3(i - 1, 2), :), :), refEl, xieta_check_pre)
  !             IF ((xieta_check_pre(1, 1) .ge. -1) .AND. (xieta_check_pre(1, 1) .le. 1)) THEN
  !                 correl3(i, 2) = correl3(i - 1, 2)
  !                 CYCLE
  !             ELSE
  !                 IF ((xieta_check_pre(1, 2) .ge. -1) .AND. (xieta_check_pre(1, 2) .le. 1)) THEN
  !                     correl3(i, 2) = correl3(i - 1, 2)
  !                     CYCLE
  !                 END IF
  !             END IF
  !             ! Check on next element, if it is fine, keep it
  !             CALL inverse_isop_transf(RESHAPE(x_all_face(i,:),[1,2]), X(T(correl3(i + 1, 2), :), :), refEl, xieta_check_pre)
  !             IF ((xieta_check_pre(1, 1) .ge. -1) .AND. (xieta_check_pre(1, 1) .le. 1)) THEN
  !                 correl3(i, 2) = correl3(i + 1, 2)
  !                 CYCLE
  !             ELSE
  !                 IF ((xieta_check_pre(1, 2) .ge. -1) .AND. (xieta_check_pre(1, 2) .le. 1)) THEN
  !                   correl3(i, 2) = correl3(i - 1, 2)
  !                   CYCLE
  !                 END IF
  !             END IF
  !
  !             ! Use this only in extreme cases where next and precedent elements are not concerned(!)
  !             DO j = 1, SIZE(correl3, 1)
  !               IF (ANY(ExtFaces_prec(:, 1) .eq. correl3(j, 2))) THEN ! Check from the beginning for an external face that matches up with xieta of the node
  !                   CALL inverse_isop_transf(RESHAPE(x_all_face(i,:),[1,2]), X(T(correl3(j, 2), :), :), refEl, xieta_check)
  !                   IF ((xieta_check(1, 1) .ge. -1) .AND. (xieta_check(1, 1) .le. 1)) THEN
  !                       correl3(i, 2) = correl3(j, 2)
  !                       EXIT
  !                   ELSE
  !                       IF ((xieta_check(1, 2) .ge. -1) .AND. (xieta_check(1, 2) .le. 1)) THEN
  !                         correl3(i, 2) = correl3(j, 2)
  !                         EXIT
  !                       END IF
  !                   END IF
  !               END IF
  !           END DO
  !         ENDIF
  !       ENDDO
  !
  !       ALLOCATE(face_nodes(SIZE(correl3,1), SIZE(refEl%Face_nodes,2)))
  !       ALLOCATE(iface_vec(SIZE(correl3,1),1))
  !
  !       face_nodes = 0
  !       iface_vec = 0
  !       DO i = 1, SIZE(correl3, 1)
  !           DO j = 1, SIZE(ExtFaces_prec, 1)
  !               IF (ExtFaces_prec(j, 1) .eq. correl3(i, 2)) THEN
  !                   face_nodes(i, :) = refEl%Face_nodes(ExtFaces_prec(j, 2), :) ! All the nodes on the external face in the old mesh
  !                   iface_vec(i, 1) = ExtFaces_prec(j, 2) ! All the nodes on the external face in the old mesh
  !               ENDIF
  !           ENDDO
  !       ENDDO
  !
  !
  !       index_el_nextTo = 0
  !       index_el_nextTo_ext = 0
  !
  !       IF (ANY(iface_vec(:,1) .eq. 0)) THEN
  !           n = COUNT(iface_vec(:,1) .eq. 0)
  !           ALLOCATE(indices(n))
  !
  !           counter = 1
  !           DO i = 1, SIZE(iface_vec)
  !             IF(iface_vec(i,1) .eq. 0) THEN
  !               indices(counter) = i
  !               counter = counter + 1
  !             ENDIF
  !           ENDDO
  !
  !           index_nodes_correl = T(correl3(:,2),1:3); ! indices of the nodes of the elements found
  !
  !           CALL unique_1D(RESHAPE(index_nodes_correl,[SIZE(index_nodes_correl)]),index_nodes_correl_unique) !made as an array
  !
  !           counter = 1
  !           DO i = 1,n_elements
  !             DO j = 1,SIZE(index_nodes_correl_unique)
  !                   IF (ANY(T(i,:) .eq. index_nodes_correl_unique(j))) THEN
  !                       index_el_nextTo(counter) = i
  !                       counter = counter + 1
  !                   ENDIF
  !               ENDDO
  !           ENDDO
  !
  !           ! list of elements close to the two individuated
  !           CALL unique_1D(index_el_nextTo,index_el_nextTo_unique)
  !
  !           counter = 1
  !           ! start from 2 cause index_el_nextTo_unique(1) is always equal to 0
  !           DO i=2,SIZE(index_el_nextTo_unique) !find which ones of those are external
  !               DO j=1,SIZE(ExtFaces_prec,1)
  !                   IF (index_el_nextTo_unique(i) .eq. ExtFaces_prec(j,1)) THEN
  !                       index_el_nextTo_ext(counter,1) = index_el_nextTo_unique(i)
  !                       index_el_nextTo_ext(counter,2) = ExtFaces_prec(j,2)
  !                       counter = counter + 1
  !                   ENDIF
  !               ENDDO
  !           ENDDO
  !
  !           DO i = 1,SIZE(indices)
  !               DO j = 1, SIZE(index_el_nextTo_ext,1)
  !                   CALL inverse_isop_transf(RESHAPE(x_all_face(indices(i),:),[1,2]),X(T(index_el_nextTo_ext(j,1),:),:),refEl,xieta_node1) ! %projection on el-y
  !                   IF (index_el_nextTo_ext(j,2) .eq. 1) THEN
  !                       IF ((xieta_node1(1,1) .gt. -1) .AND. (xieta_node1(1,1) .lt. 1)) THEN
  !                           correl3(indices(i),2) = index_el_nextTo_ext(j,1);
  !                           iface_vec(indices(i),1)  = index_el_nextTo_ext(j,2);
  !                           face_nodes(indices(i),:) = refEl%Face_nodes(index_el_nextTo_ext(j,2),:);
  !                           EXIT
  !                       ENDIF
  !                   ELSEIF (index_el_nextTo_ext(j,2) .eq. 2) THEN
  !                       IF (((xieta_node1(1,1) .gt. -1) .AND. (xieta_node1(1,1) .lt. 1))  .OR. ((xieta_node1(1,2) .gt. -1) .AND. (xieta_node1(1,2) .lt. 1))) THEN
  !                           correl3(indices(i),2) = index_el_nextTo_ext(j,1);
  !                           iface_vec(indices(i),1) = index_el_nextTo_ext(j,2);
  !                           face_nodes(indices(i),:) = refEl%Face_nodes(index_el_nextTo_ext(j,2),:)
  !                           EXIT
  !                       ENDIF
  !                   ELSEIF (index_el_nextTo_ext(j,2) .eq. 3) THEN
  !                       IF ((xieta_node1(1,2) .gt. -1) .AND. (xieta_node1(1,2) .lt. 1)) THEN
  !                           correl3(indices(i),2) = index_el_nextTo_ext(j,1);
  !                           iface_vec(indices(i),1) = index_el_nextTo_ext(j,2);
  !                           face_nodes(indices(i),:) = refEl%Face_nodes(index_el_nextTo_ext(j,2),:);
  !                           EXIT
  !                       ENDIF
  !                   ENDIF
  !               ENDDO
  !           ENDDO
  !       ENDIF
  !
  !       ! find common external vertex in old mesh
  !       CALL unique_stable(correl3(:,2),correl3_unique) !vertex elements
  !
  !       IF(ALLOCATED(index_common_nodes_check_unique)) DEALLOCATE(index_common_nodes_check_unique)
  !       IF (SIZE(correl3_unique) .eq. 2) THEN
  !            CALL intersect_stable_int(T(correl3_unique(1),1:3),T(correl3_unique(2),1:3),index_common_nodes_check);
  !            ALLOCATE(index_common_nodes_check_unique(SIZE(index_common_nodes_check)))
  !            index_common_nodes_check_unique = index_common_nodes_check
  !       ELSEIF (SIZE(correl3_unique) .eq. 3) THEN
  !            CALL intersect_stable_int(T(correl3_unique(1),1:3),T(correl3_unique(2),1:3),index_common_nodes_check1);
  !            CALL intersect_stable_int(T(correl3_unique(2),1:3),T(correl3_unique(3),1:3),index_common_nodes_check2);
  !            CALL intersect_stable_int(T(correl3_unique(3),1:3),T(correl3_unique(1),1:3),index_common_nodes_check3);
  !            ALLOCATE(index_common_nodes_check(SIZE(index_common_nodes_check1)+SIZE(index_common_nodes_check2) + SIZE(index_common_nodes_check3)))
  !            !index_common_nodes_check = RESHAPE([index_common_nodes_check1,index_common_nodes_check2,index_common_nodes_check3], [SIZE(index_common_nodes_check1)+SIZE(index_common_nodes_check2) + SIZE(index_common_nodes_check3)]);
  !            index_common_nodes_check = [index_common_nodes_check1, index_common_nodes_check2, index_common_nodes_check3]
  !            CALL unique_stable(index_common_nodes_check,index_common_nodes_check_unique);
  !       ENDIF
  !
  !       index_middleNode_oldMesh = -1 ! check condition
  !
  !       DO i = 1, SIZE(correl3,1)-1
  !         CALL intersect_stable_int(T(correl3(i,2),1:3),T(correl3(i+1,2),1:3),index_common_nodes)
  !         ! two vertices in common
  !         IF(SIZE(index_common_nodes) .eq. 2) THEN
  !           IF (ANY(T(correl3(i,2),face_nodes(i,:)) .eq. index_common_nodes(1))) THEN
  !                   index_middleNode_oldMesh = index_common_nodes(1)
  !               ELSEIF (ANY(T(correl3(i,2),face_nodes(i,:)) .eq. index_common_nodes(2))) THEN
  !                   index_middleNode_oldMesh = index_common_nodes(2)
  !               ELSEIF (ANY(T(correl3(i+1,2),face_nodes(i+1,:)) .eq. index_common_nodes(1))) THEN
  !                   index_middleNode_oldMesh = index_common_nodes(1)
  !               ELSEIF (ANY(T(correl3(i+1,2),face_nodes(i+1,:)) .eq. index_common_nodes(2))) THEN
  !                   index_middleNode_oldMesh = index_common_nodes(2)
  !               ELSE
  !                   WRITE(*,*) "Error in get_coords_rounded_mesh. STOP."
  !                   STOP
  !               ENDIF
  !           ELSEIF (SIZE(index_common_nodes) .eq. 3) THEN   ! three vertices in common
  !               IF(ALLOCATED(indices)) DEALLOCATE(indices)
  !               IF (iface_vec(i+1,1) .eq. 1) THEN
  !                   index_middleNode_oldMesh = index_common_nodes(2);
  !                   CALL find_matches_int(index_common_nodes_check_unique, index_middleNode_oldMesh, indices);
  !                   IF (SIZE(indices) .eq. 0) THEN
  !                     index_middleNode_oldMesh = index_common_nodes(1);
  !                   ENDIF
  !               ELSEIF (iface_vec(i+1,1) .eq. 2) THEN
  !                   index_middleNode_oldMesh = index_common_nodes(3);
  !                   CALL find_matches_int(index_common_nodes_check_unique, index_middleNode_oldMesh, indices);
  !                   IF (SIZE(indices) .eq. 0) THEN
  !                       index_middleNode_oldMesh = index_common_nodes(2);
  !                   ENDIF
  !               ELSEIF (iface_vec(i+1,1) .eq. 3) THEN
  !                   index_middleNode_oldMesh = index_common_nodes(1);
  !                   CALL find_matches_int(index_common_nodes_check_unique, index_middleNode_oldMesh, indices)
  !                   IF (SIZE(indices) .eq. 0) THEN
  !                       index_middleNode_oldMesh = index_common_nodes(3);
  !                   ENDIF
  !               ENDIF
  !           ELSE  ! one vertices in common
  !               index_middleNode_oldMesh = index_common_nodes(1)
  !               IF (index_middleNode_oldMesh .eq. -1) THEN !middle node not found
  !                   WRITE(*,*) "Error in get_coords_rounded_mesh. Middle node not found. STOP."
  !                   STOP
  !               ENDIF
  !           ENDIF
  !
  !           ! this is useless?
  !           ! ALLOCATE(xieta_vertex_middle1(SIZE(correl3_unique)-1,2))
  !           ! ALLOCATE(xieta_vertex_middle2(SIZE(correl3_unique)-1,2))
  !           !
  !           !
  !           ! xieta_vertex_middle1 = 0.
  !           ! xieta_vertex_middle2 = 0.
  !           ! temp = 0.
  !           ! WRITE(*,*) "init13"
  !           ! CALL inverse_isop_transf(RESHAPE(X(index_middleNode_oldMesh,:), [1,2]),X(T(correl3_unique(1),:),:),refEl,temp) ! find xieta of common node in refEl1
  !           ! xieta_vertex_middle1(i,:) = temp(1,:)
  !           ! temp = 0.
  !           ! WRITE(*,*) "init15"
  !           ! CALL inverse_isop_transf(RESHAPE(X(index_middleNode_oldMesh,:), [1,2]),X(T(correl3_unique(SIZE(correl3_unique)),:),:),refEl, temp); ! find xieta of common node in refEl2
  !           ! xieta_vertex_middle2(i,:) = temp(1,:)
  !           ! DEALLOCATE(temp, xieta_vertex_middle1, )
  !         ENDDO
  !
  !         ALLOCATE(temp(1,2))
  !         ! find xieta of first node of the face that are on the rounded border
  !         CALL inverse_isop_transf(RESHAPE(x_all_face(1,:), [1, 2]),X(T(correl3_unique(1),:),:),refEl, temp)
  !         DO i = 1, SIZE(xieta_vertex_first,1)
  !           xieta_vertex_first(i,:) = temp(1,:)
  !         ENDDO
  !         ! find xieta of last node of the face that are on the rounded border
  !         CALL inverse_isop_transf(RESHAPE(x_all_face(SIZE(x_all_face,1),:), [1,2]),X(T(correl3_unique(SIZE(correl3_unique)),:),:),refEl, temp)
  !         DO i = 1, SIZE(xieta_vertex_last,1)
  !           xieta_vertex_last(i,:) = temp(1,:)
  !         ENDDO
  !         ! find xieta of outside node in refEl1
  !         CALL inverse_isop_transf(RESHAPE(x_vertex_nodes(SIZE(x_vertex_nodes,1),:), [1,2]),X(T(correl3_unique(1),:),:),refEl, temp)
  !         DO i = 1, SIZE(xieta_vertex_first_outNode,1)
  !           xieta_vertex_first_outNode(i,:) = temp(1,:)
  !         ENDDO
  !
  !         ! find xieta of outside node in refEl2
  !         CALL inverse_isop_transf(RESHAPE(x_vertex_nodes(1,:), [1,2]),X(T(correl3_unique(SIZE(correl3_unique)),:),:),refEl, temp)
  !         DO i = 1, SIZE(xieta_vertex_last_outNode,1)
  !           xieta_vertex_last_outNode(i,:) = temp(1,:)
  !         ENDDO
  !
  !         DEALLOCATE(temp)
  !
  !         ! loop on internal nodes --> find new coords for each
  !         DO i = 2, SIZE(correl3,1)-1
  !           ! nodes in first el
  !           IF(correl3(i, 2) .eq. correl3_unique(1)) THEN
  !             IF (iface_vec(i, 1) .eq. 1) THEN
  !               xieta1(1, 1) = xieta_vertex_first(1, 1) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex_first_outNode(1, 1) - xieta_vertex_first(1, 1))
  !               xieta1(1, 2) = -1
  !             ELSE IF (iface_vec(i, 1) .eq. 2) THEN
  !               xieta1(1, 1) = xieta_vertex_first(1, 1) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex_first_outNode(1, 1) - xieta_vertex_first(1, 1))
  !               xieta1(1, 2) = xieta_vertex_first(1, 2) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex_first_outNode(1, 2) - xieta_vertex_first(1, 2))
  !             ELSE IF (iface_vec(i, 1) .eq. 3) THEN
  !               xieta1(1, 1) = -1
  !               xieta1(1, 2) = xieta_vertex_first(1, 2) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex_first_outNode(1, 2) - xieta_vertex_first(1, 2))
  !             END IF
  !
  !             CALL compute_shape_functions_at_points(refEl, xieta1, shapeFunctions1)
  !             X_new(i-1, :) = RESHAPE(MATMUL(TRANSPOSE(shapeFunctions1(:,:,1)), X(T(correl3_unique(1), :), :)), [SIZE(X_new,2)])
  !
  !           ELSE IF (correl3(i, 2) .eq. correl3_unique(SIZE(correl3_unique))) THEN
  !             IF (iface_vec(i, 1) == 1) THEN
  !               xieta1(1, 1) = xieta_vertex_last_outNode(SIZE(xieta_vertex_last_outNode,1), 1) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex_last(SIZE(xieta_vertex_last,1), 1) - xieta_vertex_last_outNode(SIZE(xieta_vertex_last_outNode,1), 1))
  !               xieta1(1, 2) = -1
  !             ELSE IF (iface_vec(i, 1) == 2) THEN
  !               xieta1(1, 1) = xieta_vertex_last_outNode(SIZE(xieta_vertex_last_outNode,1), 1) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex_last(SIZE(xieta_vertex_last,1), 1) - xieta_vertex_last_outNode(SIZE(xieta_vertex_last_outNode,1), 1))
  !               xieta1(1, 2) = xieta_vertex_last_outNode(SIZE(xieta_vertex_last_outNode,1), 2) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex_last(SIZE(xieta_vertex_last,1), 2) - xieta_vertex_last_outNode(SIZE(xieta_vertex_last_outNode,1), 2))
  !             ELSE IF (iface_vec(i, 1) == 3) THEN
  !               xieta1(1, 1) = -1
  !               xieta1(1, 2) = xieta_vertex_last_outNode(SIZE(xieta_vertex_last_outNode,1), 2) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex_last(SIZE(xieta_vertex_last,1), 2) - xieta_vertex_last_outNode(SIZE(xieta_vertex_last_outNode,1), 2))
  !             END IF
  !
  !             ! Assuming shapeFunctions2, N_vertex2, and X are declared and defined elsewhere
  !             CALL compute_shape_functions_at_points(refEl, xieta1, shapeFunctions1)
  !             X_new(i - 1, :) = RESHAPE(MATMUL(TRANSPOSE(shapeFunctions1(:,:,1)), X(T(correl3_unique(SIZE(correl3_unique)), :), :)), [SIZE(X_new,2)])
  !
  !           ELSE
  !             CALL inverse_isop_transf(RESHAPE(x_all_face(i, :),[1,2]), X(T(correl3(i, 2), :), :), refEl, xieta_vertex1)
  !             CALL inverse_isop_transf(RESHAPE(x_all_face(i, :),[1,2]), X(T(correl3(i, 2), :), :), refEl, xieta_vertex2)
  !
  !             ! Case where element is on the external face for Ne > 2
  !             IF (iface_vec(i, 1) == 1) THEN
  !               xieta1(1, 1) = xieta_vertex1(1, 1) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex2(1,1) - xieta_vertex1(1, 1))
  !               xieta1(1, 2) = -1
  !             ELSE IF (iface_vec(i, 1) == 2) THEN
  !               xieta1(1, 1) = xieta_vertex1(1, 1) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex2(1,1) - xieta_vertex1(1, 1))
  !               xieta1(1, 2) = xieta_vertex1(1, 2) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex2(1,2) - xieta_vertex1(1, 2))
  !             ELSE IF (iface_vec(i, 1) == 3) THEN
  !               xieta1(1, 1) = -1
  !               xieta1(1, 2) = xieta_vertex1(1, 2) + ((refEl%coord1d(i) + 1) / 2) * (xieta_vertex2(1,2) - xieta_vertex1(1, 2))
  !             END IF
  !
  !             CALL  compute_shape_functions_at_points(refEl, xieta1, shapeFunctions2)
  !             X_new(i-1, :) = RESHAPE(MATMUL(TRANSPOSE(shapeFunctions2(:,:,1)), X(T(correl3(i, 2), :), :)), [SIZE(X_new,2)])
  !           ENDIF
  !       ENDDO
  !
  !       CALL nodes_correct_order(X_new,x_vertex_nodes) ! function for ordering nodes based on the distance from a vertex (in case they are not ordered)
  !     ENDIF
  ! END SUBROUTINE get_coords_rounded_mesh
  !
  ! SUBROUTINE nodes_correct_order(X, x_vertex_nodes)
  !     REAL*8, INTENT(INOUT)         :: X(:,:)
  !     REAL*8, INTENT(IN)            :: x_vertex_nodes(:,:)
  !     REAL*8                        :: distance(size(X,1),2), X_temp(size(X,1),2)
  !     REAL*8                        :: correct_order_nodes(size(X,1),2)
  !     INTEGER                       :: num_int_nodes, i
  !
  !     num_int_nodes = size(X,1)
  !     distance = 0
  !
  !     DO i = 1, num_int_nodes
  !       distance(i,2) = i
  !     END DO
  !
  !     distance(:,1) = SQRT((X(:,1) - x_vertex_nodes(1,1))**2 + (X(:,2) - x_vertex_nodes(1,2))**2)
  !
  !     CALL sort_rows(distance, correct_order_nodes)
  !
  !     X_temp = 0
  !
  !     DO i=1,num_int_nodes
  !         X_temp(i,:) = X(INT(correct_order_nodes(i,2)),:)
  !     ENDDO
  !
  !     X = X_temp
  ! ENDSUBROUTINE nodes_correct_order
  !
  ! SUBROUTINE sort_rows(matrix, matrix_sorted)
  !     REAL*8, INTENT(IN)                :: matrix(:,:)
  !     REAL*8, INTENT(OUT)               :: matrix_sorted(:,:)
  !     REAL*8                            :: buffer(SIZE(matrix,2))
  !     INTEGER                           :: nsize, irow, krow
  !
  !     nsize = SIZE(matrix,1)
  !
  !     DO irow = 1, nsize
  !         krow = minloc( matrix( irow:nsize, 1 ), dim=1) + irow - 1
  !         buffer = matrix(irow,:)
  !         matrix_sorted(irow,:) = matrix( krow,:)
  !         matrix_sorted(krow,:) = buffer(:)
  !     ENDDO
  ! ENDSUBROUTINE sort_rows

  ! SUBROUTINE new_points_2D(p, porder, sin, pnew)
  !
  !   ! p: nodal coordinates of the element (x-y)
  !   ! porder: degree of interpolation
  !   ! s: nodal coordinates of the element (xi-eta)
  !   ! note: Reference triangle is [0,0; 1,0; 0,1];
  !
  !   IMPLICIT NONE
  !   INTEGER, INTENT(IN) :: porder
  !   REAL*8, DIMENSION(:,:), INTENT(IN)               :: sin
  !   REAL*8, DIMENSION(:,:), INTENT(IN)               :: p
  !   REAL*8                                           :: s(SIZE(sin,1),SIZE(sin,2))
  !   REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: pnew
  !   REAL*8, DIMENSION(:), ALLOCATABLE                :: sc, sc_des
  !   REAL*8, DIMENSION(:, :), ALLOCATABLE             :: sb, xp, V
  !   REAL*8, DIMENSION(:,:), ALLOCATABLE              :: alpha
  !   INTEGER                                          :: i, j, n, counter, start, end
  !
  !   s = sin
  !
  !   s = 0.5D0 * (s + 1.0D0)
  !
  !   n = COUNT(s(:,2) .lt. 1e-12)
  !
  !   ALLOCATE(sc(n))
  !   ALLOCATE(sc_des(n))
  !
  !   counter = 1
  !   DO i = 1, SIZE(s,1)
  !     IF(s(i,2) .lt. 1e-12) THEN
  !       sc(counter) = s(i,1)
  !       counter = counter + 1
  !     ENDIF
  !   ENDDO
  !
  !
  !   CALL quicksort_real(sc)
  !
  !   ! invert the order (descending)
  !   DO i = 0, SIZE(sc)-1
  !     sc_des(SIZE(sc)-i) = sc(i+1)
  !   ENDDO
  !
  !   ! Boundary nodes on local coordinates (xi-eta)
  !   ALLOCATE(sb(SIZE(sc)+SIZE(sc_des)-1 + SIZE(sc_des)-2, 2))
  !   sb = 0.0
  !   start = 1
  !   end   = SIZE(sc)
  !   sb(start:end, 1) = sc
  !   start = end + 1
  !   end = start + SIZE(sc_des(2:SIZE(sc_des))) - 1
  !   sb(start:end, 1) = sc_des(2:SIZE(sc_des))
  !   sb(start:end, 2) = sc(2:SIZE(sc_des))
  !   start = end + 1
  !   end = start  + SIZE(sc_des(2:SIZE(sc_des)-1)) -1
  !   sb(start:end,2) = sc_des(2:SIZE(sc_des)-1)
  !
  !
  !   ! Boundary nodes on cartesian coordinates (x-y)
  !   ALLOCATE(xp(SIZE(sb,1),SIZE(sb,2)))
  !   xp = 0.0
  !
  !   DO i = 1, SIZE(sb, 1)
  !     DO j = 1, SIZE(s, 1)
  !       IF (NORM2(sb(i, :) - s(j, 1:2)) .lt. 1.0E-8) THEN
  !         xp(i, :) = p(j, :)
  !       END IF
  !     END DO
  !   END DO
  !
  !   ! Compute V matrix
  !   ALLOCATE(V(SIZE(sb, 1), SIZE(sb, 1)))
  !
  !   V = 0.
  !   DO i = 1, SIZE(sb, 1)
  !     CALL blending_boundary(sb, i, porder, sc, V(:, i))
  !   END DO
  !
  !   ! Solve for alpha
  !   ALLOCATE(alpha(SIZE(xp, 1),SIZE(xp, 2)))
  !   CALL solve_linear_system_sing(V, xp(:,1), alpha(:,1))
  !   CALL solve_linear_system_sing(V, xp(:,2), alpha(:,2))
  !
  !   DEALLOCATE(V)
  !
  !   ALLOCATE(V(SIZE(s,1), SIZE(sb,1)))
  !
  !   V = 0.
  !   ! Compute pnew
  !
  !   DO i = 1, SIZE(sb, 1)
  !     CALL blending_boundary(s(:, 1:2), i, porder, sc, V(:, i))
  !   END DO
  !
  !   ALLOCATE(pnew(SIZE(V, 1), SIZE(alpha, 2)))
  !
  !   pnew(:,1) = MATMUL(V, alpha(:,1))
  !   pnew(:,2) = MATMUL(V, alpha(:,2))
  !
  !   DEALLOCATE(sc)
  !   DEALLOCATE(sc_des)
  !   DEALLOCATE(sb)
  !   DEALLOCATE(xp)
  !   DEALLOCATE(V)
  !   DEALLOCATE(alpha)
  ! ENDSUBROUTINE new_points_2D



END MODULE adaptivity_common_module
