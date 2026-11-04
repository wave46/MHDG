!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! module for schock capturing adaptivity
!************************************************************

MODULE adaptivity_common_module
   USE, INTRINSIC               :: iso_c_binding, ONLY: c_double, c_int, c_size_t
   USE globals
   USE reference_element
   USE gmsh
   USE GMSH_io_module
   IMPLICIT NONE

   TYPE                         :: gmsh_entity_mesh_data
      INTEGER                   :: DIMENSION = -1
      INTEGER                   :: tag = -1
      INTEGER                   :: element_type = -1
      INTEGER(c_size_t), ALLOCATABLE :: element_tags(:)
      INTEGER(c_size_t), ALLOCATABLE :: node_tags(:)
   END TYPE gmsh_entity_mesh_data

CONTAINS

   SUBROUTINE adaptivity_console_output()
      CHARACTER(1024)           :: buffer

      IF ((adapt%evaluator .EQ. 2)) THEN
         buffer = "            ADAPTIVITY ESTIMATOR                 "
      ELSEIF ((adapt%evaluator .EQ. 1)) THEN
         buffer = "            ADAPTIVITY INDICATOR                 "
      ELSEIF ((adapt%evaluator .EQ. 0)) THEN
         buffer = "        ADAPTIVITY ESTIMATOR-INDICATOR           "
      END IF

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "*************************************************"
         WRITE (*, *) TRIM(buffer)
         WRITE (*, *) "*************************************************"
      END IF

   END SUBROUTINE adaptivity_console_output

   SUBROUTINE calculate_h_map_elements(nodes, connectivity, h_map)
      REAL*8, INTENT(IN)                :: nodes(:, :)
      INTEGER, INTENT(IN)               :: connectivity(:, :)
      REAL*8, INTENT(OUT)               :: h_map(SIZE(connectivity, 1))
      INTEGER                           :: i
      REAL*8                            :: side1, side2, side3, s, area

      DO i = 1, SIZE(connectivity, 1)
         side1 = SQRT((nodes(connectivity(i,1),1) - nodes(connectivity(i,2),1))**2 + (nodes(connectivity(i,1),2) - nodes(connectivity(i,2),2))**2)
         side2 = SQRT((nodes(connectivity(i,2),1) - nodes(connectivity(i,3),1))**2 + (nodes(connectivity(i,2),2) - nodes(connectivity(i,3),2))**2)
         side3 = SQRT((nodes(connectivity(i,3),1) - nodes(connectivity(i,1),1))**2 + (nodes(connectivity(i,3),2) - nodes(connectivity(i,1),2))**2)

         s = (side1 + side2 + side3)/2.0
         area = SQRT(s*(s - side1)*(s - side2)*(s - side3))
         !double circumradius
         h_map(i) = (side1*side2*side3)/(2.0*area)
      END DO
   END SUBROUTINE calculate_h_map_elements

   SUBROUTINE get_h_target_vertices(h_map_elements, h_target_vertices, T)
      REAL*8, INTENT(IN)                              :: h_map_elements(:)
      INTEGER, INTENT(IN)                              :: T(:, :)
      REAL*8, DIMENSION(:), POINTER, INTENT(INOUT)    :: h_target_vertices
      REAL*8, ALLOCATABLE                             :: h_target_nodal(:)
      INTEGER, ALLOCATABLE                            :: nodes_repeats(:)
      INTEGER                                         :: number_of_vertices

      ALLOCATE (h_target_nodal(MAXVAL(T)))
      ALLOCATE (nodes_repeats(MAXVAL(T)))

      CALL sum_h_target_nodal(T, h_map_elements, h_target_nodal, nodes_repeats)

      number_of_vertices = COUNT(nodes_repeats /= 0)

      ALLOCATE (h_target_vertices(number_of_vertices))

      CALL average_h_target(h_target_nodal, nodes_repeats, h_target_vertices)

      DEALLOCATE (h_target_nodal, nodes_repeats)
   END SUBROUTINE get_h_target_vertices

   ! subroutine compute_nodal_characteristic_lengths_tri(num_nodes, num_elements, num_dims, &
   !      X, T, estimated_lc)
   !   !
   !   ! Computes the average length of all edges connected to each node in a mesh
   !   ! composed only of triangles.
   !   !
   !   ! Arguments:
   !   integer, intent(in)      :: num_nodes              ! Total number of nodes
   !   integer, intent(in)      :: num_elements           ! Total number of elements
   !   integer, intent(in)      :: num_dims               ! Number of spatial dimensions (2 or 3)
   !   real(kind=8), intent(in) :: X(num_nodes, num_dims) ! Nodal coordinates (real(kind=8) for double precision)
   !   integer, intent(in)      :: T(num_elements, 3)     ! Connectivity matrix (1-based indexing assumed for Fortran)
   !   ! Fixed to 3 columns for triangles
   !   real(kind=8), intent(out) :: estimated_lc(num_nodes) ! Output: Estimated characteristic lengths for each node
   !
   !   ! Local variables
   !   integer                  :: i, j, k, elem_idx, node_idx_1, node_idx_2
   !   real(kind=8)             :: edge_length
   !   real(kind=8), allocatable :: node_lc_sum(:)
   !   integer, allocatable     :: node_edge_counts(:)
   !
   !   ! --- Define edges for triangles (local indices within the element's row in T) ---
   !   ! For a triangle with nodes (T_elem_row(1), T_elem_row(2), T_elem_row(3))
   !   ! Edges are: (node 1, node 2), (node 2, node 3), (node 3, node 1)
   !   integer, parameter       :: local_triangle_edges(3, 2) = RESHAPE([1, 2, 2, 3, 3, 1], shape=[3, 2])
   !   integer, parameter       :: num_nodes_per_elem = 3 ! Fixed for triangles
   !   integer, parameter       :: num_edges_per_elem = 3 ! Fixed for triangles
   !
   !
   !   ! --- Initialize arrays ---
   !   allocate(node_lc_sum(num_nodes))
   !   allocate(node_edge_counts(num_nodes))
   !   node_lc_sum = 0.0_8
   !   node_edge_counts = 0
   !
   !   ! --- Iterate through each element ---
   !   do elem_idx = 1, num_elements
   !      ! --- Iterate through each edge of the current triangle ---
   !      do i = 1, num_edges_per_elem
   !         ! Get the global node indices (Fortran uses 1-based indexing, T is assumed 1-based)
   !         node_idx_1 = T(elem_idx, local_triangle_edges(i, 1))
   !         node_idx_2 = T(elem_idx, local_triangle_edges(i, 2))
   !
   !         ! Calculate the Euclidean distance (length) of the edge
   !         edge_length = 0.0_8
   !         do j = 1, num_dims
   !            edge_length = edge_length + (X(node_idx_1, j) - X(node_idx_2, j))**2
   !         end do
   !         edge_length = SQRT(edge_length)
   !
   !         ! Accumulate the sum of edge lengths and count for both nodes on this edge
   !         node_lc_sum(node_idx_1) = node_lc_sum(node_idx_1) + edge_length
   !         node_edge_counts(node_idx_1) = node_edge_counts(node_idx_1) + 1
   !
   !         ! Make sure we don't double-count if the edge is self-referencing (shouldn't happen)
   !         if (node_idx_1 /= node_idx_2) then
   !            node_lc_sum(node_idx_2) = node_lc_sum(node_idx_2) + edge_length
   !            node_edge_counts(node_idx_2) = node_edge_counts(node_idx_2) + 1
   !         end if
   !      end do ! i (edges per element)
   !   end do ! elem_idx
   !
   !   ! --- Calculate the average for each node ---
   !   do i = 1, num_nodes
   !      if (node_edge_counts(i) > 0) then
   !         estimated_lc(i) = node_lc_sum(i) / real(node_edge_counts(i), kind=8)
   !      else
   !         ! Handle isolated nodes (nodes with no connected edges).
   !         ! This shouldn't happen in a valid mesh.
   !         WRITE(*,*) 'Warning: Node', i, 'has no connected edges. Assigning default LC = 1.0_8'
   !         estimated_lc(i) = 1.0_8 ! Default value, adjust as needed
   !      end if
   !   end do
   !
   !   deallocate(node_lc_sum)
   !   deallocate(node_edge_counts)
   !
   ! end subroutine compute_nodal_characteristic_lengths_tri

   SUBROUTINE load_new_mesh_gmsh(order)
      USE preprocess
      INTEGER, INTENT(IN)       :: order
      INTEGER                     :: ierr

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "********** Loading new mesh  **********"
      END IF

      CALL free_mesh

      CALL free_reference_element_pol(refElPol)
      CALL create_reference_element(refElPol, 2, order, verbose=0)

      IF ((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
         CALL load_gmsh_mesh("./res/temp", 1)
      ELSE
         CALL load_gmsh_mesh("./res/temp", 1)
      END IF

      CALL mesh_preprocess_serial(ierr)

      Mesh%X = Mesh%X*phys%lscale

   END SUBROUTINE load_new_mesh_gmsh

   SUBROUTINE combine_h_target_ind_est(h_map_elements, h_target_elements_est, h_target_elements_ind, h_target_elements)
      REAL*8, INTENT(IN)        :: h_map_elements(:)
      REAL*8, INTENT(IN)        :: h_target_elements_est(:)
      REAL*8, INTENT(IN)        :: h_target_elements_ind(:)
      REAL*8, INTENT(OUT)       :: h_target_elements(:)
      REAL*8, PARAMETER         :: tol = 1.0E-10

      h_target_elements = h_target_elements_est
      WHERE (ABS(h_target_elements_ind - h_map_elements) .GT. tol)
         h_target_elements = h_target_elements_ind
      END WHERE
   END SUBROUTINE combine_h_target_ind_est

   SUBROUTINE sum_h_target_nodal(T, h_map_elements, h_target_nodal, nodes_repeats)
      INTEGER, INTENT(IN)               :: T(:, :)
      REAL*8, INTENT(IN)                :: h_map_elements(:)
      REAL*8, INTENT(OUT)               :: h_target_nodal(:)
      INTEGER, INTENT(OUT)              :: nodes_repeats(:)
      INTEGER                           :: i, j

      h_target_nodal = 0.
      nodes_repeats = 0
      DO i = 1, SIZE(T, 1)
         DO j = 1, 3
            h_target_nodal(T(i, j)) = h_target_nodal(T(i, j)) + h_map_elements(i)
            nodes_repeats(T(i, j)) = nodes_repeats(T(i, j)) + 1
         END DO
      END DO
   END SUBROUTINE sum_h_target_nodal

   SUBROUTINE average_h_target(h_target_nodal, nodes_repeats, h_target_vertices)
      REAL*8, INTENT(IN)                              :: h_target_nodal(:)
      INTEGER, INTENT(IN)                             :: nodes_repeats(:)
      REAL*8, DIMENSION(:), POINTER, INTENT(INOUT)      :: h_target_vertices
      INTEGER                                         :: i, j

      j = 1
      DO i = 1, SIZE(h_target_nodal)
         IF (nodes_repeats(i) /= 0) THEN
            h_target_vertices(j) = h_target_nodal(i)/REAL(nodes_repeats(i))
            j = j + 1
         END IF
      END DO
   END SUBROUTINE average_h_target

   ! subroutine set_order_mesh(p)

   !    integer, intent(in)                      :: p
   !    TYPE(Reference_element_type)             :: refElLocal
   !    integer                                  :: n_element_nodes, n_face_nodes, n_int_faces, n_nodes, ini,ind, i, j, elem, iface, counter, element_order, elemType, ielem, meshed_face, &
   !       n_already_meshed_faces, n_boundaries, n_boundary_elements, Nnodesperelem, n_elements, Ndim, Nelems, Nextfaces, Nnodes, &
   !       Nnodesperface, counter1, start, stop_index
   !    integer                                  :: temp(p-1), elem_pos(2), face_pos(2), face_nodes(3,p-1), face_info(3), elements(2), nodes_face(2, p-1), element_face(2), already_meshed_face_nodes(2, p-1), already_meshed_element_faces(2), &
   !       already_meshed_element(2), ifacenode(p + 1)
   !    integer, allocatable, dimension(:, :)     :: int_faces, elem_int_face, Tp, Tb
   !    integer, allocatable, dimension(:)       :: boundaryFlag, unique_boundary_flag
   !    logical, allocatable, dimension(:)       :: aux_coord_logical, local_coord_logical, int_face_meshed
   !    real*8, allocatable, dimension(:, :)      :: Xp, Xp_aux, elem_nodes_mod, coord_ref
   !    integer , allocatable, dimension(:,:)    :: Tb_Dirichlet, Tb_LEFT, Tb_RIGHT, Tb_UP, Tb_DOWN, Tb_WALL, Tb_LIM, Tb_IN, Tb_OUT, Tb_ULIM, Tb_PUFF, Tb_PUMP, mesh_info

   !    if (utils%printint > 0) then
   !       if (MPIvar%glob_id .eq. 0) then
   !          WRITE (6, *) '*************************************************'
   !          WRITE (6, *) '*              INCREASE ORDER MESH              *'
   !          WRITE (6, *) '*************************************************'
   !       end if
   !    end if

   !    call create_reference_element(refElLocal, 2, p, verbose=0)

   !    face_nodes = refElLocal%face_nodes(:, 2:SIZE(refElLocal%face_nodes, 2) - 1); ! without vertices
   !    allocate (coord_ref(SIZE(refElLocal%coord2d, 1), SIZE(refElLocal%coord2d, 2)))
   !    allocate (int_faces(SIZE(Mesh%intfaces, 1), SIZE(Mesh%intfaces, 2)))

   !    coord_ref = refElLocal%coord2d;
   !    n_element_nodes = SIZE(coord_ref, 1);
   !    n_face_nodes = SIZE(face_nodes, 2);
   !    int_faces = Mesh%intfaces
   !    elem_int_face = Mesh%F

   !    n_int_faces = SIZE(int_faces, 1);
   !    do i = 1, SIZE(elem_int_face, 1)
   !       do j = 1, SIZE(elem_int_face, 2)
   !          if (elem_int_face(i, j) .gt. n_int_faces) then
   !             elem_int_face(i, j) = 0
   !          end if
   !       end do
   !    end do

   !    n_elements = SIZE(Mesh%T, 1);
   !    n_nodes = SIZE(Mesh%X, 1);
   !    allocate (Xp(n_element_nodes*n_elements, 2))
   !    allocate (Tp(n_elements, n_element_nodes))
   !    allocate (aux_coord_logical(n_element_nodes))
   !    allocate (local_coord_logical(n_element_nodes))
   !    allocate (int_face_meshed(n_int_faces))
   !    allocate (elem_nodes_mod(n_element_nodes, 2))

   !    Xp(1:n_nodes, :) = Mesh%X
   !    Tp(:, 1:3) = Mesh%T

   !    elem_pos(1) = 1
   !    elem_pos(2) = 3
   !    face_pos(1) = 2
   !    face_pos(2) = 4

   !    aux_coord_logical = .true.
   !    aux_coord_logical(1:3) = .false.
   !    int_face_meshed = .false.

   !    ini = n_nodes + 1;
   !    counter1 = 1

   !    do elem = 1, n_elements

   !       local_coord_logical = aux_coord_logical

   !       ! read the faces infos
   !       face_info = elem_int_face(elem, :)

   !       ! count how many faces are valid .neq. 0
   !       n_already_meshed_faces = 0
   !       do iface = 1, 3
   !          if (face_info(iface) .ne. 0) then
   !             if (int_face_meshed(face_info(iface)) .eqv. .true.) then
   !                n_already_meshed_faces = n_already_meshed_faces + 1
   !             end if
   !          end if
   !       end do

   !       ! if there is at least one face valid
   !       if (n_already_meshed_faces .le. 3) then
   !          ! loop through all faces of the element
   !          do iface = 1, 3
   !             ! if the face is valid
   !             if (face_info(iface) .ne. 0) then
   !                ! if it is already remeshed
   !                if (int_face_meshed(face_info(iface)) .eqv. .true.) then
   !                   meshed_face = face_info(iface)
   !                   elements(1) = int_faces(face_info(iface), elem_pos(1))
   !                   elements(2) = int_faces(face_info(iface), elem_pos(2))

   !                   if (elements(1) .ne. elem) then
   !                      already_meshed_element(1) = elements(1)
   !                      already_meshed_element_faces(1) = int_faces(meshed_face, face_pos(1));
   !                      already_meshed_face_nodes(1, :) = face_nodes(already_meshed_element_faces(1), :);
   !                   else
   !                      element_face(1) = int_faces(meshed_face, face_pos(1));
   !                      nodes_face(1, :) = face_nodes(element_face(1), :);
   !                      local_coord_logical(nodes_face(1, :)) = .false.
   !                   end if

   !                   if (elements(2) .ne. elem) then
   !                      already_meshed_element(2) = elements(2)
   !                      already_meshed_element_faces(2) = int_faces(meshed_face, face_pos(2));
   !                      already_meshed_face_nodes(2, :) = face_nodes(already_meshed_element_faces(2), :);
   !                   else
   !                      element_face(2) = int_faces(meshed_face, face_pos(2));
   !                      nodes_face(2, :) = face_nodes(element_face(2), :);
   !                      local_coord_logical(nodes_face(2, :)) = .false.
   !                   end if

   !                   if ((elements(1) .ne. elem) .and. (elements(2) .eq. elem)) then
   !                      call fliplr_int(Tp(already_meshed_element(1), already_meshed_face_nodes(1, :)), temp);
   !                      Tp(elem, nodes_face(2, :)) = temp
   !                   end if

   !                   if ((elements(1) .eq. elem) .and. (elements(2) .ne. elem)) then
   !                      call fliplr_int(Tp(already_meshed_element(2), already_meshed_face_nodes(2, :)), temp)
   !                      Tp(elem, nodes_face(1, :)) = temp
   !                   end if
   !                end if
   !                int_face_meshed(face_info(iface)) = .true.
   !             end if
   !          end do
   !       end if

   !       elem_nodes_mod = 0

   !       do i = 1, SIZE(local_coord_logical)
   !          if (local_coord_logical(i) .eqv. .true.) then
   !             call linear_mapping(Mesh%X(Mesh%T(elem, :), :), coord_ref(i, :), elem_nodes_mod(i, :))
   !          end if
   !       end do

   !       ind = ini

   !       do i = 1, SIZE(local_coord_logical)
   !          if (local_coord_logical(i) .eqv. .true.) then
   !             Xp(ind, :) = elem_nodes_mod(i, :)
   !             Tp(elem, i) = ind
   !             ind = ind + 1
   !          end if
   !       end do

   !       if ((n_element_nodes - n_already_meshed_faces*n_face_nodes - 3 - 1) .ge. 0) then
   !          ini = ini - 1 + n_element_nodes - n_already_meshed_faces*n_face_nodes - 3 + 1
   !       end if

   !    end do

   !    allocate (Xp_aux(ini - 1, 2))
   !    Xp_aux(:, :) = Xp(1:ini - 1, :)

   !    n_face_nodes = SIZE(refElLocal%Face_nodes, 2)

   !    call unique_stable(Mesh%boundaryFlag, unique_boundary_flag)

   !    n_boundaries = SIZE(unique_boundary_flag)

   !    do i = 1, n_boundaries
   !       n_boundary_elements = COUNT(Mesh%boundaryFlag .eq. unique_boundary_flag(i))
   !       if (n_boundary_elements .ne. 0) then
   !          if (unique_boundary_flag(i) .eq. 5) then
   !             allocate (Tb_PUMP(n_boundary_elements, n_face_nodes))
   !          elseif (unique_boundary_flag(i) .eq. 6) then
   !             allocate (Tb_PUFF(n_boundary_elements, n_face_nodes))
   !          elseif (unique_boundary_flag(i) .eq. 7) then
   !             allocate (Tb_LIM(n_boundary_elements, n_face_nodes))
   !          elseif (unique_boundary_flag(i) .eq. 8) then
   !             allocate (Tb_IN(n_boundary_elements, n_face_nodes))
   !          elseif (unique_boundary_flag(i) .eq. 9) then
   !             allocate (Tb_OUT(n_boundary_elements, n_face_nodes))
   !          end if
   !       end if
   !    end do

   !    start = 0
   !    stop_index = 0
   !    do i = 1, n_boundaries

   !       start = stop_index + 1
   !       n_boundary_elements = COUNT(Mesh%boundaryFlag .eq. unique_boundary_flag(i))
   !       stop_index = start + n_boundary_elements - 1

   !       counter = 1
   !       do j = start, stop_index
   !          ielem = Mesh%extfaces(j, 1)
   !          iface = Mesh%extfaces(j, 2)
   !          ifacenode = refElLocal%face_nodes(iface, :)
   !          if (unique_boundary_flag(i) .eq. 5) Tb_PUMP(counter, :) = Tp(ielem, ifacenode)
   !          if (unique_boundary_flag(i) .eq. 6) Tb_PUFF(counter, :) = Tp(ielem, ifacenode)
   !          if (unique_boundary_flag(i) .eq. 7) Tb_LIM(counter, :) = Tp(ielem, ifacenode)
   !          if (unique_boundary_flag(i) .eq. 8) Tb_IN(counter, :) = Tp(ielem, ifacenode)
   !          if (unique_boundary_flag(i) .eq. 9) Tb_OUT(counter, :) = Tp(ielem, ifacenode)
   !          counter = counter + 1
   !       end do
   !    end do

   !    Ndim = SIZE(Xp_aux, 2)
   !    Nelems = SIZE(Tp, 1)

   !    Nextfaces = 0
   !    if (ALLOCATED(Tb_IN)) then
   !       Nextfaces = SIZE(Tb_IN, 1)
   !    end if
   !    if (ALLOCATED(Tb_LIM)) then
   !       Nextfaces = Nextfaces + SIZE(Tb_LIM, 1)
   !    end if
   !    if (ALLOCATED(Tb_OUT)) then
   !       Nextfaces = Nextfaces + SIZE(Tb_OUT, 1)
   !    end if
   !    if (ALLOCATED(Tb_PUFF)) then
   !       Nextfaces = Nextfaces + SIZE(Tb_PUFF, 1)
   !    end if
   !    if (ALLOCATED(Tb_PUMP)) then
   !       Nextfaces = Nextfaces + SIZE(Tb_PUMP, 1)
   !    end if

   !    Nnodes = SIZE(Xp_aux, 1)
   !    Nnodesperelem = SIZE(Tp, 2)
   !    Nnodesperface = n_face_nodes
   !    elemType = Mesh%elemType
   !    element_order = n_face_nodes

   !    allocate (mesh_info(SIZE(Mesh%extfaces, 1), SIZE(Mesh%extfaces, 2)))

   !    mesh_info = Mesh%extfaces

   !    call free_mesh

   !    !call generate_elemface_info(Tp,Tb_IN, Tb_LIM, Tb_PUFF, Tb_PUMP, Tb_OUT, p+1, total_face_info)
   !    call generate_boundary_names(Tb_Dirichlet, Tb_LEFT, Tb_RIGHT, Tb_UP, Tb_DOWN, Tb_WALL, Tb_LIM, Tb_IN, Tb_OUT, Tb_PUFF, Tb_PUMP, Tb_ULIM, Tb, boundaryFlag, element_order)
   !    call load_mesh2global_var(Ndim, Nelems, Nextfaces, Nnodes, Nnodesperelem, Nnodesperface, elemType, Tp, Xp_aux, Tb, boundaryFlag)

   !    call free_reference_element_pol(refElLocal)
   !    deallocate (int_faces, elem_int_face, Tp, Tb)
   !    deallocate (boundaryFlag, unique_boundary_flag)
   !    deallocate (aux_coord_logical, local_coord_logical, int_face_meshed)
   !    deallocate (Xp, Xp_aux, elem_nodes_mod, coord_ref)
   !    if (ALLOCATED(Tb_Dirichlet)) deallocate (Tb_Dirichlet)
   !    if (ALLOCATED(Tb_PUMP)) deallocate (Tb_PUMP)
   !    if (ALLOCATED(Tb_PUFF)) deallocate (Tb_PUFF)
   !    if (ALLOCATED(Tb_LEFT)) deallocate (Tb_LEFT)
   !    if (ALLOCATED(Tb_RIGHT)) deallocate (Tb_RIGHT)
   !    if (ALLOCATED(Tb_UP)) deallocate (Tb_UP)
   !    if (ALLOCATED(Tb_DOWN)) deallocate (Tb_DOWN)
   !    if (ALLOCATED(Tb_WALL)) deallocate (Tb_WALL)
   !    if (ALLOCATED(Tb_LIM)) deallocate (Tb_LIM)
   !    if (ALLOCATED(Tb_IN)) deallocate (Tb_IN)
   !    if (ALLOCATED(Tb_OUT)) deallocate (Tb_OUT)
   !    if (ALLOCATED(Tb_ULIM)) deallocate (Tb_ULIM)
   ! end subroutine set_order_mesh

   SUBROUTINE set_order_mesh(mesh_filename, geometry_filename, p, ordered_mesh_name)
      USE mpi, ONLY: MPI_BARRIER, MPI_COMM_WORLD

      CHARACTER(*), INTENT(IN)  :: mesh_filename, geometry_filename
      INTEGER, INTENT(IN)       :: p
      CHARACTER(*), INTENT(OUT) :: ordered_mesh_name
      INTEGER                   :: barrier_error

      IF (LEN_TRIM(geometry_filename) .EQ. 0) &
         CALL mesh_order_error('set_2d_order requires geometry_path')

      ordered_mesh_name = './res/mesh_ordered_input'
      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE(*,*) 'Increasing mesh order on the supplied CAD geometry.'
         CALL elevate_mesh_order_with_gmsh(TRIM(mesh_filename), TRIM(geometry_filename), &
            TRIM(ordered_mesh_name)//'.msh', p)
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, barrier_error)
   END SUBROUTINE set_order_mesh

   SUBROUTINE elevate_mesh_order_with_gmsh(mesh_filename, geometry_filename, output_filename, p)
      CHARACTER(*), INTENT(IN)                   :: mesh_filename, geometry_filename, output_filename
      INTEGER, INTENT(IN)                        :: p
      TYPE(gmsh_t)                               :: gmsh_l
      TYPE(gmsh_entity_mesh_data), ALLOCATABLE   :: source_entities(:)
      INTEGER(c_size_t), ALLOCATABLE             :: source_node_tags(:)
      REAL(c_double), ALLOCATABLE                :: source_coordinates(:)
      INTEGER, ALLOCATABLE                       :: source_node_dimensions(:), source_node_entities(:)
      INTEGER                                    :: source_entity_count

      CALL gmsh_l%initialize()
      CALL gmsh_l%option%setNumber('General.Verbosity', 2.d0)
      CALL gmsh_l%open(mesh_filename)
      CALL capture_source_mesh(gmsh_l, source_node_tags, source_coordinates, source_node_dimensions, &
         source_node_entities, source_entities, source_entity_count)

      CALL gmsh_l%clear()
      CALL gmsh_l%open(geometry_filename)
      CALL classify_cad_point_nodes(gmsh_l, source_coordinates, source_node_dimensions, source_node_entities, &
         source_entities, source_entity_count)
      CALL add_source_nodes_to_cad(gmsh_l, source_node_tags, source_coordinates, source_node_dimensions, &
         source_node_entities)
      CALL add_source_elements_to_cad(gmsh_l, source_entities, source_entity_count)

      CALL gmsh_l%model%mesh%setOrder(p)
      CALL gmsh_l%model%mesh%optimize('HighOrder')
      CALL gmsh_l%option%setNumber('Mesh.MshFileVersion', 2.2d0)
      CALL gmsh_l%write(output_filename)
      CALL gmsh_l%finalize()
   END SUBROUTINE elevate_mesh_order_with_gmsh

   SUBROUTINE capture_source_mesh(gmsh_l, node_tags, coordinates, node_dimensions, node_entities, &
      entity_data, entity_count)
      TYPE(gmsh_t), INTENT(INOUT)                         :: gmsh_l
      INTEGER(c_size_t), ALLOCATABLE, INTENT(OUT)         :: node_tags(:)
      REAL(c_double), ALLOCATABLE, INTENT(OUT)            :: coordinates(:)
      INTEGER, ALLOCATABLE, INTENT(OUT)                   :: node_dimensions(:), node_entities(:)
      TYPE(gmsh_entity_mesh_data), ALLOCATABLE, INTENT(OUT) :: entity_data(:)
      INTEGER, INTENT(OUT)                                :: entity_count
      INTEGER(c_int), ALLOCATABLE                         :: dimension_tags(:,:), element_types(:)
      INTEGER(c_size_t), ALLOCATABLE                      :: entity_node_tags(:)
      REAL(c_double), ALLOCATABLE                         :: unused_coordinates(:), unused_parameters(:)
      INTEGER, ALLOCATABLE                                :: node_index(:)
      INTEGER                                             :: entity, node, position, DIMENSION, tag, source_element_type

      CALL gmsh_l%model%mesh%getNodes(node_tags, coordinates, unused_parameters, returnParametricCoord=.FALSE.)
      IF (SIZE(node_tags) .EQ. 0) CALL mesh_order_error('source mesh has no nodes')
      ALLOCATE(node_index(INT(MAXVAL(node_tags))))
      node_index = 0
      DO node = 1, SIZE(node_tags)
         node_index(INT(node_tags(node))) = node
      END DO

      ALLOCATE(node_dimensions(SIZE(node_tags)), node_entities(SIZE(node_tags)))
      node_dimensions = -1
      node_entities = -1
      CALL gmsh_l%model%getEntities(dimension_tags)
      ALLOCATE(entity_data(SIZE(dimension_tags,2)))
      entity_count = 0

      DO entity = 1, SIZE(dimension_tags,2)
         DIMENSION = INT(dimension_tags(1,entity))
         tag = INT(dimension_tags(2,entity))
         IF (DIMENSION .LT. 0 .OR. DIMENSION .GT. 2) CYCLE
         CALL gmsh_l%model%mesh%getNodes(entity_node_tags, unused_coordinates, unused_parameters, &
            dim=DIMENSION, tag=tag, includeBoundary=.FALSE., returnParametricCoord=.FALSE.)
         DO node = 1, SIZE(entity_node_tags)
            position = node_index(INT(entity_node_tags(node)))
            node_dimensions(position) = DIMENSION
            node_entities(position) = tag
         END DO

         IF (DIMENSION .EQ. 0) CYCLE
         source_element_type = MERGE(1, 2, DIMENSION .EQ. 1)
         CALL gmsh_l%model%mesh%getElementTypes(element_types, dim=DIMENSION, tag=tag)
         IF (ANY(element_types .NE. source_element_type)) &
            CALL mesh_order_error('set_2d_order requires a first-order line/triangle mesh')
         IF (SIZE(element_types) .EQ. 0) CYCLE

         entity_count = entity_count + 1
         entity_data(entity_count)%DIMENSION = DIMENSION
         entity_data(entity_count)%tag = tag
         entity_data(entity_count)%element_type = source_element_type
         CALL gmsh_l%model%mesh%getElementsByType(source_element_type, &
            entity_data(entity_count)%element_tags, entity_data(entity_count)%node_tags, tag=tag)
      END DO

      IF (entity_count .EQ. 0) CALL mesh_order_error('source mesh has no line or triangle elements')
      IF (ANY(node_dimensions .EQ. -1)) CALL mesh_order_error('source mesh contains unclassified nodes')
   END SUBROUTINE capture_source_mesh

   SUBROUTINE classify_cad_point_nodes(gmsh_l, coordinates, node_dimensions, node_entities, &
      entity_data, entity_count)
      TYPE(gmsh_t), INTENT(INOUT)              :: gmsh_l
      REAL(c_double), INTENT(IN)                :: coordinates(:)
      INTEGER, INTENT(INOUT)                    :: node_dimensions(:), node_entities(:)
      TYPE(gmsh_entity_mesh_data), INTENT(IN)   :: entity_data(:)
      INTEGER, INTENT(IN)                       :: entity_count
      INTEGER(c_int), ALLOCATABLE               :: curves(:,:), curve_boundary(:,:)
      INTEGER, ALLOCATABLE                      :: point_tags(:)
      REAL(c_double), ALLOCATABLE               :: point_coordinates(:), no_parameters(:)
      REAL(c_double)                            :: tolerance, distance, nearest_distance, coordinate_scale
      INTEGER                                   :: curve, boundary, point_count, point, node, nearest_node, point_tag

      CALL gmsh_l%model%getEntities(curves, dim=1)
      ALLOCATE(point_tags(MAX(1,2*SIZE(curves,2))), no_parameters(0))
      point_count = 0
      DO curve = 1, SIZE(curves,2)
         IF (.NOT. entity_has_elements(entity_data, entity_count, 1, INT(curves(2,curve)))) CYCLE
         CALL gmsh_l%model%getBoundary(curves(:,curve:curve), curve_boundary, combined=.FALSE., oriented=.FALSE.)
         DO boundary = 1, SIZE(curve_boundary,2)
            IF (curve_boundary(1,boundary) .NE. 0_c_int) CYCLE
            point_tag = ABS(INT(curve_boundary(2,boundary)))
            IF (point_count .EQ. 0 .OR. .NOT. ANY(point_tags(1:point_count) .EQ. point_tag)) THEN
               point_count = point_count + 1
               point_tags(point_count) = point_tag
            END IF
         END DO
      END DO

      coordinate_scale = MAX(1.d0, MAXVAL(coordinates(1::3))-MINVAL(coordinates(1::3)), &
         MAXVAL(coordinates(2   ::3))-MINVAL(coordinates(2::3)))
      tolerance = 1.d-8*coordinate_scale
      DO point = 1, point_count
         CALL gmsh_l%model%getValue(0, point_tags(point), no_parameters, point_coordinates)
         nearest_node = 0
         nearest_distance = HUGE(1.d0)
         DO node = 1, SIZE(node_dimensions)
            IF (node_dimensions(node) .GT. 1) CYCLE
            distance = SQRT(SUM((coordinates(3*node-2:3*node)-point_coordinates(1:3))**2))
            IF (distance .LT. nearest_distance) THEN
               nearest_distance = distance
               nearest_node = node
            END IF
         END DO
         IF (nearest_node .EQ. 0 .OR. nearest_distance .GT. tolerance) &
            CALL mesh_order_error('CAD endpoint does not match a source mesh vertex')
         node_dimensions(nearest_node) = 0
         node_entities(nearest_node) = point_tags(point)
      END DO
   END SUBROUTINE classify_cad_point_nodes

   LOGICAL FUNCTION entity_has_elements(entity_data, entity_count, DIMENSION, tag)
      TYPE(gmsh_entity_mesh_data), INTENT(IN) :: entity_data(:)
      INTEGER, INTENT(IN)                     :: entity_count, DIMENSION, tag
      INTEGER                                 :: entity

      entity_has_elements = .FALSE.
      DO entity = 1, entity_count
         IF (entity_data(entity)%DIMENSION .EQ. DIMENSION .AND. entity_data(entity)%tag .EQ. tag) THEN
            entity_has_elements = .TRUE.
            RETURN
         END IF
      END DO
   END FUNCTION entity_has_elements

   SUBROUTINE add_source_nodes_to_cad(gmsh_l, node_tags, coordinates, node_dimensions, node_entities)
      TYPE(gmsh_t), INTENT(INOUT)          :: gmsh_l
      INTEGER(c_size_t), INTENT(IN)        :: node_tags(:)
      REAL(c_double), INTENT(IN)           :: coordinates(:)
      INTEGER, INTENT(IN)                  :: node_dimensions(:), node_entities(:)
      INTEGER(c_int), ALLOCATABLE          :: cad_entities(:,:)
      INTEGER(c_size_t), ALLOCATABLE       :: entity_node_tags(:)
      REAL(c_double), ALLOCATABLE          :: entity_coordinates(:), parameters(:)
      INTEGER                              :: DIMENSION, entity, tag, count_nodes, node, local_node, added_nodes

      added_nodes = 0
      DO DIMENSION = 0, 2
         CALL gmsh_l%model%getEntities(cad_entities, dim=DIMENSION)
         DO entity = 1, SIZE(cad_entities,2)
            tag = INT(cad_entities(2,entity))
            count_nodes = COUNT(node_dimensions .EQ. DIMENSION .AND. node_entities .EQ. tag)
            IF (count_nodes .EQ. 0) CYCLE
            ALLOCATE(entity_node_tags(count_nodes), entity_coordinates(3*count_nodes))
            local_node = 0
            DO node = 1, SIZE(node_tags)
               IF (node_dimensions(node) .NE. DIMENSION .OR. node_entities(node) .NE. tag) CYCLE
               local_node = local_node + 1
               entity_node_tags(local_node) = node_tags(node)
               entity_coordinates(3*local_node-2:3*local_node) = coordinates(3*node-2:3*node)
            END DO

            IF (DIMENSION .EQ. 0) THEN
               CALL gmsh_l%model%mesh%addNodes(DIMENSION, tag, entity_node_tags, entity_coordinates)
            ELSE
               CALL gmsh_l%model%getParametrization(DIMENSION, tag, entity_coordinates, parameters)
               CALL gmsh_l%model%mesh%addNodes(DIMENSION, tag, entity_node_tags, entity_coordinates, &
                  parametricCoord=parameters)
               DEALLOCATE(parameters)
            END IF
            added_nodes = added_nodes + count_nodes
            DEALLOCATE(entity_node_tags, entity_coordinates)
         END DO
      END DO
      IF (added_nodes .NE. SIZE(node_tags)) CALL mesh_order_error('source mesh and CAD entity tags do not match')
   END SUBROUTINE add_source_nodes_to_cad

   SUBROUTINE add_source_elements_to_cad(gmsh_l, entity_data, entity_count)
      TYPE(gmsh_t), INTENT(INOUT)             :: gmsh_l
      TYPE(gmsh_entity_mesh_data), INTENT(IN) :: entity_data(:)
      INTEGER, INTENT(IN)                     :: entity_count
      INTEGER                                 :: entity

      DO entity = 1, entity_count
         CALL gmsh_l%model%mesh%addElementsByType(entity_data(entity)%tag, entity_data(entity)%element_type, &
            entity_data(entity)%element_tags, entity_data(entity)%node_tags)
      END DO
   END SUBROUTINE add_source_elements_to_cad

   SUBROUTINE mesh_order_error(message)
      CHARACTER(*), INTENT(IN)  :: message

      WRITE(*,*) 'Unable to increase mesh order on the CAD geometry: ', TRIM(message)
      ERROR STOP 1
   END SUBROUTINE mesh_order_error

   SUBROUTINE fliplr_int(arrin, arrout)
      INTEGER, DIMENSION(:), INTENT(IN) :: arrin
      INTEGER, DIMENSION(:), INTENT(INOUT) :: arrout
      INTEGER                   :: i, n
      INTEGER                   :: temp

      n = SIZE(arrin)

      DO i = 1, n/2
         ! Swap elements from the left and right sides
         temp = arrin(i)
         arrout(i) = arrin(n - i + 1)
         arrout(n - i + 1) = temp
      END DO

      ! swap middle element
      IF (MOD(n, 2) .EQ. 1) arrout(INT(n/2) + 1) = arrin(INT(n/2) + 1)
   END SUBROUTINE fliplr_int

   SUBROUTINE linear_mapping(vertCoord, xiVector, X)
      ! Linear mapping between local and cartesian coordinates
      ! Input:
      !   vertCoord: vertexes of the element
      !   xiVector:  point in local coordinates
      ! Output:
      !   X: point in cartesian coordinates

      REAL(8), DIMENSION(:, :), INTENT(IN) :: vertCoord
      REAL(8), DIMENSION(:), INTENT(IN)   :: xiVector
      REAL(8), DIMENSION(:), INTENT(OUT)  :: X

      REAL*8                              :: N(3)

      CALL linear_shape_functions_2D(xiVector, N)

      X(1) = dot_PRODUCT(N, vertCoord(:, 1))
      X(2) = dot_PRODUCT(N, vertCoord(:, 2))

   CONTAINS

      SUBROUTINE linear_shape_functions_2D(xiVector, N)
         ! Reference triangle is [-1,-1; 1,-1; -1,1]
         REAL(8), DIMENSION(:), INTENT(IN) :: xiVector
         REAL(8), DIMENSION(3), INTENT(OUT) :: N

         REAL(8)                :: xi, eta

         xi = xiVector(1)
         eta = xiVector(2)
         N(1) = -xi - eta
         N(2) = 1.0d0 + xi
         N(3) = 1.0d0 + eta
         N = 0.5*N

      END SUBROUTINE linear_shape_functions_2D

   END SUBROUTINE linear_mapping

   SUBROUTINE inverse_isop_transf(x, Xe, refEl, xieta)
      TYPE(Reference_element_type), INTENT(IN)  :: RefEl
      REAL*8, INTENT(OUT)                       :: xieta(:, :)
      REAL*8, INTENT(IN)                        :: x(:, :), Xe(:, :)
      REAL*8                                    :: x0(SIZE(x, 1), SIZE(x, 2))
      INTEGER                                   :: maxit, npoints, nnodes, i, j, k, counter, n
      REAL*8, ALLOCATABLE                       :: p(:), dpxi(:), dpeta(:), xind(:, :), x0ind(:, :), xietaind(:, :), rhs(:, :)
      INTEGER, ALLOCATABLE                      :: ind(:)
      REAL*8                                    :: xieta0(SIZE(x,1),SIZE(x,2)), aux_xieta(SIZE(x,1),SIZE(x,2)), Nx(refEl%Nnodes2D), Ny(refEl%Nnodes2D)
      REAL*8                                    :: Vand(refEl%Nnodes2D, refEl%Nnodes2D), invV(refEl%Nnodes2D, refEl%Nnodes2D)
      REAL*8                                    :: Jxx, Jxy, Jyx, Jyy, detJ
      REAL*8                                    :: tol

      maxit = 5
      tol = 1e-10
      npoints = SIZE(x, 1)
      nnodes = SIZE(Xe, 1)

      CALL inverse_linear_transformation(x, Xe, xieta0)

      ! just fucking brute force it
      DO j = 1, SIZE(xieta0, 2)
         DO i = 1, SIZE(xieta0, 1)
            IF (ABS(xieta0(i, j) - 1.0) .LT. 1e-12) THEN
               xieta0(i, j) = xieta0(i, j) - 1.e-10
            END IF
         END DO
      END DO

      CALL iso_transformation_high_order(xieta0, Xe, refEl, x0)

      n = COUNT(SQRT((x(:, 1) - x0(:, 1))**2 + (x(:, 2) - x0(:, 2))**2) .GT. (tol*SQRT(x(:, 1)**2 + x(:, 2)**2) + 1.e-14))

      IF (n .NE. 0) THEN

         ALLOCATE (ind(n))
         ALLOCATE (rhs(n, SIZE(x, 2)))
         ALLOCATE (xind(n, SIZE(x, 2)))
         ALLOCATE (x0ind(n, SIZE(x0, 2)))
         ALLOCATE (xietaind(n, SIZE(xieta, 2)))

         ind = 0.

         counter = 1
         DO i = 1, npoints
            IF ((SQRT((x(i, 1) - x0(i, 1))**2 + (x(i, 2) - x0(i, 2))**2)) .GT. (tol*SQRT(x(i, 1)**2 + x(i, 2)**2) + 1.e-14)) THEN
               ind(counter) = counter
               counter = counter + 1
            END IF
         END DO

         xind = x(ind, :)
         x0ind = x0(ind, :)
         xietaind = xieta0(ind, :)

         ALLOCATE (p(refEl%Nnodes2D), dpxi(refEl%Nnodes2D), dpeta(refEl%Nnodes2D))

         DO i = 1, maxit
            IF (ALL((SQRT((xind(:,1)-x0ind(:,1))**2+(xind(:,2)-x0ind(:,2))**2)) .LT. (tol*SQRT(xind(:,1)**2+xind(:,2)**2)+1.e-14))) THEN
               EXIT
            END IF

            CALL vandermonde_2d(Vand, refEl)
            CALL invert_matrix(TRANSPOSE(Vand), invV)

            p = 0.d0
            dpxi = 0.d0
            dpeta = 0.d0

            DO j = 1, SIZE(xietaind, 1)
               CALL orthopoly2d_deriv(xietaind(j, 1), xietaind(j, 2), refEl%Ndeg, refEl%Nnodes2D, p, dpxi, dpeta)

               Nx = MATMUL(invV, dpxi)
               Ny = MATMUL(invV, dpeta)

               Jxx = dot_PRODUCT(Nx, Xe(:, 1))
               Jxy = dot_PRODUCT(Ny, Xe(:, 1))
               Jyx = dot_PRODUCT(Nx, Xe(:, 2))
               Jyy = dot_PRODUCT(Ny, Xe(:, 2))
               detJ = Jxx*Jyy - Jxy*Jyx
               rhs = xind - x0ind

               xietaind(j, 1) = xietaind(j, 1) + (rhs(j, 1)*Jyy - rhs(j, 2)*Jxy)/detJ
               xietaind(j, 2) = xietaind(j, 2) + (rhs(j, 2)*Jxx - rhs(j, 1)*Jyx)/detJ

            END DO

            ! just fucking brute force it
            DO k = 1, SIZE(xietaind, 2)
               DO j = 1, SIZE(xietaind, 1)
                  IF (ABS(xietaind(j, k) - 1.0) .LT. 1e-12) THEN
                     xietaind(j, k) = xietaind(j, k) - 1.e-10
                  END IF
               END DO
            END DO

            CALL iso_transformation_high_order(xietaind, Xe, refEl, x0ind)
         END DO

         IF(ANY((SQRT((xind(:,1)-x0ind(:,1))**2+(xind(:,2)-x0ind(:,2))**2)) .GT. (tol*SQRT(xind(:,1)**2+xind(:,2)**2)+1.e-14))) THEN
            WRITE (*, *) "inverse_isop_transf non converging."
            STOP
         END IF

         aux_xieta = xieta0
         aux_xieta(ind, :) = xietaind
         xieta = aux_xieta

         DEALLOCATE (p)
         DEALLOCATE (dpxi)
         DEALLOCATE (dpeta)
         DEALLOCATE (ind)
         DEALLOCATE (xind)
         DEALLOCATE (x0ind)
         DEALLOCATE (xietaind)
         DEALLOCATE (rhs)

      ELSE
         xieta = xieta0
      END IF
   END SUBROUTINE inverse_isop_transf

   SUBROUTINE inverse_linear_transformation(x, Xe, xieta)
      USE LinearAlgebra, ONLY: solve_linear_system

      REAL*8, INTENT(IN)        :: x(:, :), Xe(:, :)
      REAL*8, INTENT(OUT)       :: xieta(:, :)
      REAL*8                    :: x1(2), x2(2), x3(2), J(2, 2), aux(SIZE(x, 1), 2), xieta_temp(SIZE(xieta, 2), SIZE(xieta, 1))

      ! take vertices
      x1 = Xe(1, :)
      x2 = Xe(2, :)
      x3 = Xe(3, :)

      J(:, 1) = (x2 - x1)/2
      J(:, 2) = (x3 - x1)/2

      aux(:, 1) = x(:, 1) - (x2(1) + x3(1))/2
      aux(:, 2) = x(:, 2) - (x2(2) + x3(2))/2

      CALL solve_linear_system(J, TRANSPOSE(aux), xieta_temp)
      xieta = TRANSPOSE(xieta_temp)

   END SUBROUTINE inverse_linear_transformation

   SUBROUTINE iso_transformation_high_order(xieta, Xe, refEl, x)
      TYPE(Reference_element_type)      :: refEl
      REAL*8, INTENT(IN)               :: xieta(:, :)
      REAL*8, INTENT(IN)               :: Xe(:, :)
      REAL*8, INTENT(OUT)              :: x(:, :)
      REAL*8                            :: shapeFunctions(refEl%Nnodes2D, SIZE(xieta, 1), 3)

      shapeFunctions = 0.d0

      CALL compute_shape_functions_at_points(refEl, xieta, shapeFunctions)

      ! take only shape function and not its derivatives in xi and eta
      x(:, 1) = MATMUL(TRANSPOSE(shapeFunctions(:, :, 1)), Xe(:, 1))
      x(:, 2) = MATMUL(TRANSPOSE(shapeFunctions(:, :, 1)), Xe(:, 2))

   END SUBROUTINE iso_transformation_high_order

   PURE SUBROUTINE find_matches_int(a, b, indices)
      INTEGER, DIMENSION(:), INTENT(IN)                 :: a
      INTEGER, INTENT(IN)                               :: b
      INTEGER, DIMENSION(:), INTENT(INOUT), ALLOCATABLE :: indices
      INTEGER                                           :: counter
      INTEGER                                           :: i

      counter = COUNT(a .EQ. b)

      ALLOCATE (indices(counter))

      counter = 1
      DO i = 1, SIZE(a)
         IF (a(i) .EQ. b) THEN
            indices(counter) = i
            counter = counter + 1
         END IF
      END DO

   END SUBROUTINE find_matches_int

   SUBROUTINE delete_file(filename)

      CHARACTER(*), INTENT(IN)        :: filename
      INTEGER                         :: fileID, stat

      CALL get_unit(fileID)   ! get unit file and open it

      OPEN (unit=fileID, iostat=stat, file=filename, status='old')
      IF (stat .NE. 0) THEN
         WRITE (*, *) "Problem opening file to delete."
      ELSE
         CLOSE (unit=fileID, iostat=stat, status='delete')
      END IF

      IF (stat .NE. 0) THEN
         WRITE (*, *) "Problem deleting file."
      END IF

   END SUBROUTINE delete_file

   SUBROUTINE extract_mesh_name_from_fullpath_woext(mesh_name, mesh_name_npne)
      CHARACTER(1024), INTENT(IN)             :: mesh_name
      CHARACTER(1024), INTENT(OUT)            :: mesh_name_npne
      INTEGER                                 :: i, start

      i = LEN(mesh_name)
      start = -1

      DO
         IF (i - 2 .EQ. 0) THEN
            WRITE (*, *) "GMSH file input not found, check input sintax."
            STOP
         END IF

         ! ! get the index at the dot of .msh
         ! if((mesh_name(i-2:i-2) .eq. 'm') .and. (mesh_name(i-1:i-1) .eq. 's') .and. (mesh_name(i:i) .eq. 'h')) then
         !   end = i-4
         ! ENDIF
         ! get the index at the first slash reading from right to left
         IF (mesh_name(i:i) .EQ. '/') THEN
            start = i + 1
            EXIT
         END IF
         i = i - 1
      END DO

      mesh_name_npne = TRIM(ADJUSTL(mesh_name(start:)))

   END SUBROUTINE extract_mesh_name_from_fullpath_woext

   PURE SUBROUTINE unique_1D(list_in, list_out)
      !! From a 1D array of integers list_in extracts the list of unique occurences of values
      !integer, dimension(:), intent(in) :: list_in
      INTEGER, DIMENSION(:), INTENT(IN) :: list_in
      !! The list of integers to work on
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: list_out
      !! The list in output, with a single occurence of each value in list in

      INTEGER, DIMENSION(SIZE(list_in)) :: list_in_sorted
      INTEGER                   :: nlist, n_unique, i

      nlist = SIZE(list_in)

      ! Copy list_in and sort it
      list_in_sorted = list_in
      CALL quicksort_int(list_in_sorted)

      ! The number of jumps in this list gives the number of unique values
      n_unique = COUNT((list_in_sorted(2:nlist) - list_in_sorted(1:nlist - 1)) .GT. 0) + 1
      ALLOCATE (list_out(n_unique))

      ! Now loop in the sorted list and each time a jump is found,
      n_unique = 1
      list_out(1) = list_in_sorted(1)
      DO i = 2, nlist
         IF (list_in_sorted(i) .NE. list_in_sorted(i - 1)) THEN
            n_unique = n_unique + 1
            list_out(n_unique) = list_in_sorted(i)
         END IF
      END DO

   END SUBROUTINE unique_1D

   PURE RECURSIVE SUBROUTINE quicksort_int(a)
      !! quicksort.f -*-f90-*-
      !! Author: t-nissie, some tweaks by 1AdAstra1
      !! License: GPLv3
      !! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
      INTEGER, DIMENSION(:), INTENT(INOUT) :: a
      INTEGER                   ::  x, t
      INTEGER                   :: first, last
      INTEGER i, j

      first = 1

      last = SIZE(a, 1)
      x = a((first + last)/2)
      i = first
      j = last

      DO
         DO WHILE (a(i) < x)
            i = i + 1
         END DO
         DO WHILE (x < a(j))
            j = j - 1
         END DO
         IF (i .GE. j) EXIT
         t = a(i); a(i) = a(j); a(j) = t
         i = i + 1
         j = j - 1
      END DO

      IF (first < i - 1) CALL quicksort_int(a(first:i - 1))
      IF (j + 1 < last) CALL quicksort_int(a(j + 1:last))

   END SUBROUTINE quicksort_int

   PURE RECURSIVE SUBROUTINE quicksort_real(a)
      !! quicksort.f -*-f90-*-
      !! Author: t-nissie, some tweaks by 1AdAstra1
      !! License: GPLv3
      !! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
      REAL*8, DIMENSION(:), INTENT(INOUT) :: a
      !! The array to sort

      REAL*8                    ::  x, t
      INTEGER                   :: first, last
      INTEGER i, j

      first = 1
      last = SIZE(a, 1)
      x = a((first + last)/2)
      i = first
      j = last

      DO
         DO WHILE (a(i) < x)
            i = i + 1
         END DO
         DO WHILE (x < a(j))
            j = j - 1
         END DO
         IF (i .GE. j) EXIT
         t = a(i); a(i) = a(j); a(j) = t
         i = i + 1
         j = j - 1
      END DO

      IF (first < i - 1) CALL quicksort_real(a(first:i - 1))
      IF (j + 1 < last) CALL quicksort_real(a(j + 1:last))

   END SUBROUTINE quicksort_real

   PURE SUBROUTINE unique_stable(arrayin, uniqueArr)
      INTEGER, DIMENSION(:), INTENT(IN) :: arrayin
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: uniqueArr
      INTEGER, DIMENSION(SIZE(arrayin)) :: indices
      LOGICAL, DIMENSION(SIZE(arrayin)) :: isUnique
      INTEGER                   :: i, j, nUnique

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

      ! allocate the unique array and populate it in the original order
      ALLOCATE (uniqueArr(nUnique))
      DO i = 1, nUnique
         uniqueArr(i) = arrayin(indices(i))
      END DO

   END SUBROUTINE unique_stable

   SUBROUTINE gmsh_create_from_h_target(h_target_on_elements, vertices_coordinates, connectivity, p_order)

      USE, INTRINSIC            :: iso_c_binding
      USE gmsh
      USE MPI_OMP, ONLY: OMPvar

      TYPE(gmsh_t)              :: gmsh_l
      INTEGER, INTENT(IN)       :: p_order
      REAL*8, DIMENSION(:), INTENT(IN) :: h_target_on_elements ! on the elements
      REAL*8, DIMENSION(:, :), INTENT(IN) :: vertices_coordinates ! coordinates of the vertices
      INTEGER, DIMENSION(:, :), INTENT(IN) :: connectivity ! connectivity of the triangles
      REAL*8, ALLOCATABLE       :: data_for_gmsh(:)
      INTEGER*8                 :: gmsh_dim, number_of_vertices_per_triangle, i, j, start_index
      INTEGER*4                 :: size_view, number_of_triangles, ret
      REAL*8                    :: sf_index, vertex_coordinates(2), h_target_on_vertex
      REAL*8, PARAMETER         :: minimum_mesh_size = 0.5d-4
      REAL*8, PARAMETER         :: target_size_quantum = 2.d-3*minimum_mesh_size
      REAL*8, PARAMETER         :: view_coordinate_quantum = 1.d-11

      !GMSH always have (X,Y,Z) coordinates
      gmsh_dim = 3
      number_of_vertices_per_triangle = 3
      number_of_triangles = SIZE(connectivity, 1)

      ! initalize gmsh
      CALL gmsh_l%initialize()

      ! Set verbosity level to 2 (Errors and warnings)
      CALL gmsh_l%option%setNumber("General.Verbosity", 2d0)

      !create model
      CALL gmsh_l%model%add("geo")
      ! merge
      CALL gmsh_l%merge(adapt%geometry_path)

      ALLOCATE (data_for_gmsh((number_of_vertices_per_triangle*(number_of_vertices_per_triangle + 1))*number_of_triangles))

      ! Prepare data in gmsh format
      DO i = 1, number_of_triangles
         start_index = (i - 1)*(number_of_vertices_per_triangle*(number_of_vertices_per_triangle + 1))
         DO j = 1, number_of_vertices_per_triangle
            vertex_coordinates = view_coordinate_quantum*ANINT(vertices_coordinates(connectivity(i,j),:)/view_coordinate_quantum)
            h_target_on_vertex = target_size_quantum*ANINT(h_target_on_elements(i)/target_size_quantum)
            data_for_gmsh(start_index + j) = vertex_coordinates(1)
            data_for_gmsh(start_index + number_of_vertices_per_triangle + j) = vertex_coordinates(2)
            data_for_gmsh(start_index + 2*number_of_vertices_per_triangle + j) = 0.0 ! no Z coordinate
            data_for_gmsh(start_index + 3*number_of_vertices_per_triangle + j) = h_target_on_vertex
         END DO
      END DO

      ! Add h_target as a post-processing view
      size_view = gmsh_l%view%add("h_target")
      CALL gmsh_l%view%addListData(size_view, "ST", number_of_triangles, data_for_gmsh)
      sf_index = gmsh_l%view%getIndex(size_view)

      ! Add the view as a field
      ret = gmsh_l%model%mesh%field%add("PostView")
      CALL gmsh_l%model%mesh%field%setNumber(ret, "ViewIndex", 0d0)

      ! Apply the view as the current background mesh size field:
      CALL gmsh_l%model%mesh%field%setAsBackgroundMesh(ret)

      ! ignore characteristic length from geometry
      CALL gmsh_l%option%setNumber("Mesh.RandomSeed", 1d0)
      CALL gmsh_l%option%setNumber("Mesh.MeshSizeExtendFromBoundary", 0d0)
      CALL gmsh_l%option%setNumber("Mesh.MeshSizeFromPoints", 0d0)
      CALL gmsh_l%option%setNumber("Mesh.MeshSizeFromCurvature", 0d0)
      CALL gmsh_l%option%setNumber("Mesh.MeshSizeFactor", 1d0)
      CALL gmsh_l%option%setNumber("General.NumThreads", 1d0) !real(OMPvar%Nthreads) should be use what specified by OMP_NUM_THREADS
      CALL gmsh_l%option%setNumber("Mesh.MeshSizeMin", 0.5e-4)
      !Changing the algorithm to Delaunay, the default is Frontal-Delaunay (Don't Know if needed)
      CALL gmsh_l%option%setNumber("Mesh.Algorithm", 5d0)
      CALL gmsh_l%option%setNumber("Geometry.OCCBooleanPreserveNumbering", 1d0)

      ! Generate the refined mesh
      CALL gmsh_l%model%mesh%generate(2)
      CALL gmsh_l%model%mesh%setOrder(p_order)
      CALL gmsh_l%model%mesh%optimize('HighOrder')
      CALL gmsh_l%option%setNumber("Mesh.MshFileVersion", 2.2)
      CALL gmsh_l%write('./res/temp.msh')
      CALL gmsh_l%finalize()

      DEALLOCATE (data_for_gmsh)

   END SUBROUTINE gmsh_create_from_h_target

   SUBROUTINE save_copy_new_mesh(mesh_name, count_adapt)
      USE in_out, ONLY: copy_file
      CHARACTER(1024), INTENT(IN) :: mesh_name
      INTEGER, INTENT(IN)         :: count_adapt
      CHARACTER(70)               :: param_adapt_char, count_adapt_char
      CHARACTER(1024)             :: mesh_name_npne, new_mesh_name_npne, buffer

      CALL extract_mesh_name_from_fullpath_woext(mesh_name, mesh_name_npne)

      WRITE (param_adapt_char, *) adapt%param_est
      WRITE (count_adapt_char, *) count_adapt

      new_mesh_name_npne = TRIM(ADJUSTL(mesh_name_npne))//'_param'//TRIM(ADJUSTL(param_adapt_char))//'_n'//TRIM(ADJUSTL(count_adapt_char))

      buffer = "./res/"//TRIM(ADJUSTL(new_mesh_name_npne))//".msh"

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "Mesh saved as: ", TRIM(ADJUSTL(buffer))
         WRITE (*,*)
      END IF

      CALL copy_file("./res/temp.msh", buffer)

   END SUBROUTINE save_copy_new_mesh
END MODULE adaptivity_common_module
