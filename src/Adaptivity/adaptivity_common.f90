!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! Module for schock capturing adaptivity
!************************************************************

MODULE adaptivity_common_module
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int, c_size_t
  USE globals, ONLY: adapt, free_mesh, geom, Mesh, phys, refElPol, switch
  USE MPI_OMP, ONLY: MPIvar
  USE reference_element, ONLY: create_reference_element, free_reference_element_pol
  USE gmsh, ONLY: gmsh_t
  USE GMSH_io_module, ONLY: load_gmsh_mesh
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: adaptivity_console_output
  PUBLIC :: calculate_h_map_elements
  PUBLIC :: combine_h_target_ind_est
  PUBLIC :: gmsh_create_from_h_target
  PUBLIC :: load_new_mesh_gmsh
  PUBLIC :: quicksort_int
  PUBLIC :: save_copy_new_mesh
  PUBLIC :: set_order_mesh
  PUBLIC :: unique_1D
  PUBLIC :: unique_stable

  TYPE :: gmsh_entity_mesh_data
     INTEGER :: dimension = -1
     INTEGER :: tag = -1
     INTEGER :: element_type = -1
     INTEGER(c_size_t), ALLOCATABLE :: element_tags(:)
     INTEGER(c_size_t), ALLOCATABLE :: node_tags(:)
  END TYPE gmsh_entity_mesh_data

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

   SUBROUTINE calculate_h_map_elements(nodes, connectivity, h_map)
      REAL*8, INTENT(IN)                :: nodes(:,:)
      INTEGER, INTENT(IN)               :: connectivity(:,:)
      REAL*8, INTENT(OUT)               :: h_map(SIZE(connectivity,1))
      INTEGER                           :: i

      DO i = 1, SIZE(connectivity,1)
         ! Gmsh's isotropic size field represents a desired linear edge length.
         h_map(i) = triangle_max_edge(nodes, connectivity(i,:))
      ENDDO
   END SUBROUTINE calculate_h_map_elements

   PURE FUNCTION triangle_max_edge(nodes, triangle) RESULT(max_edge)
      REAL*8, INTENT(IN) :: nodes(:,:)
      INTEGER, INTENT(IN) :: triangle(:)
      REAL*8 :: max_edge, edge_lengths(3)

      edge_lengths(1) = NORM2(nodes(triangle(1),:) - nodes(triangle(2),:))
      edge_lengths(2) = NORM2(nodes(triangle(2),:) - nodes(triangle(3),:))
      edge_lengths(3) = NORM2(nodes(triangle(3),:) - nodes(triangle(1),:))
      max_edge = MAXVAL(edge_lengths)
   END FUNCTION triangle_max_edge

   SUBROUTINE average_element_targets_on_nodes(element_targets, connectivity, node_targets)
      REAL*8, INTENT(IN)  :: element_targets(:)
      INTEGER, INTENT(IN) :: connectivity(:,:)
      REAL*8, INTENT(OUT) :: node_targets(:)
      INTEGER, ALLOCATABLE :: node_counts(:)
      INTEGER              :: i, j, node

      ALLOCATE(node_counts(SIZE(node_targets)))
      node_targets = 0.d0
      node_counts = 0
      DO i = 1, SIZE(connectivity,1)
         DO j = 1, SIZE(connectivity,2)
            node = connectivity(i,j)
            node_targets(node) = node_targets(node) + element_targets(i)
            node_counts(node) = node_counts(node) + 1
         ENDDO
      ENDDO

      DO node = 1, SIZE(node_targets)
         IF (node_counts(node) .GT. 0) &
            node_targets(node) = node_targets(node)/REAL(node_counts(node), KIND=8)
      ENDDO
      DEALLOCATE(node_counts)
   END SUBROUTINE average_element_targets_on_nodes


   SUBROUTINE load_new_mesh_gmsh(order)
      USE preprocess, ONLY: mesh_preprocess_serial
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

      Mesh%X = Mesh%X*phys%lscale

      IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
         Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
      END IF

   END SUBROUTINE load_new_mesh_gmsh

   SUBROUTINE combine_h_target_ind_est(h_map_elements,h_target_elements_est,h_target_elements_ind,h_target_elements)
      REAL*8, INTENT(IN) :: h_map_elements(:)
      REAL*8, INTENT(IN) :: h_target_elements_est(:)
      REAL*8, INTENT(IN) :: h_target_elements_ind(:)
      REAL*8, INTENT(OUT) :: h_target_elements(:)
      REAL*8, PARAMETER :: tol = 1.0E-10
      
      h_target_elements = h_target_elements_est
      WHERE(ABS(h_target_elements_ind-h_map_elements) .GT. tol)
         h_target_elements = h_target_elements_ind
      END WHERE
   ENDSUBROUTINE combine_h_target_ind_est

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
    ENDIF
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
    INTEGER                                             :: entity, node, position, dimension, tag, source_element_type

    CALL gmsh_l%model%mesh%getNodes(node_tags, coordinates, unused_parameters, returnParametricCoord=.FALSE.)
    IF (SIZE(node_tags) .EQ. 0) CALL mesh_order_error('source mesh has no nodes')
    ALLOCATE(node_index(INT(MAXVAL(node_tags))))
    node_index = 0
    DO node = 1, SIZE(node_tags)
       node_index(INT(node_tags(node))) = node
    ENDDO

    ALLOCATE(node_dimensions(SIZE(node_tags)), node_entities(SIZE(node_tags)))
    node_dimensions = -1
    node_entities = -1
    CALL gmsh_l%model%getEntities(dimension_tags)
    ALLOCATE(entity_data(SIZE(dimension_tags,2)))
    entity_count = 0

    DO entity = 1, SIZE(dimension_tags,2)
       dimension = INT(dimension_tags(1,entity))
       tag = INT(dimension_tags(2,entity))
       IF (dimension .LT. 0 .OR. dimension .GT. 2) CYCLE
       CALL gmsh_l%model%mesh%getNodes(entity_node_tags, unused_coordinates, unused_parameters, &
            dim=dimension, tag=tag, includeBoundary=.FALSE., returnParametricCoord=.FALSE.)
       DO node = 1, SIZE(entity_node_tags)
          position = node_index(INT(entity_node_tags(node)))
          node_dimensions(position) = dimension
          node_entities(position) = tag
       ENDDO

       IF (dimension .EQ. 0) CYCLE
       source_element_type = MERGE(1, 2, dimension .EQ. 1)
       CALL gmsh_l%model%mesh%getElementTypes(element_types, dim=dimension, tag=tag)
       IF (ANY(element_types .NE. source_element_type)) &
          CALL mesh_order_error('set_2d_order requires a first-order line/triangle mesh')
       IF (SIZE(element_types) .EQ. 0) CYCLE

       entity_count = entity_count + 1
       entity_data(entity_count)%dimension = dimension
       entity_data(entity_count)%tag = tag
       entity_data(entity_count)%element_type = source_element_type
       CALL gmsh_l%model%mesh%getElementsByType(source_element_type, &
            entity_data(entity_count)%element_tags, entity_data(entity_count)%node_tags, tag=tag)
    ENDDO

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
          ENDIF
       ENDDO
    ENDDO

    coordinate_scale = MAX(1.d0, MAXVAL(coordinates(1::3))-MINVAL(coordinates(1::3)), &
         MAXVAL(coordinates(2::3))-MINVAL(coordinates(2::3)))
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
          ENDIF
       ENDDO
       IF (nearest_node .EQ. 0 .OR. nearest_distance .GT. tolerance) &
          CALL mesh_order_error('CAD endpoint does not match a source mesh vertex')
       node_dimensions(nearest_node) = 0
       node_entities(nearest_node) = point_tags(point)
    ENDDO
  END SUBROUTINE classify_cad_point_nodes

  LOGICAL FUNCTION entity_has_elements(entity_data, entity_count, dimension, tag)
    TYPE(gmsh_entity_mesh_data), INTENT(IN) :: entity_data(:)
    INTEGER, INTENT(IN)                     :: entity_count, dimension, tag
    INTEGER                                 :: entity

    entity_has_elements = .FALSE.
    DO entity = 1, entity_count
       IF (entity_data(entity)%dimension .EQ. dimension .AND. entity_data(entity)%tag .EQ. tag) THEN
          entity_has_elements = .TRUE.
          RETURN
       ENDIF
    ENDDO
  END FUNCTION entity_has_elements

  SUBROUTINE add_source_nodes_to_cad(gmsh_l, node_tags, coordinates, node_dimensions, node_entities)
    TYPE(gmsh_t), INTENT(INOUT)          :: gmsh_l
    INTEGER(c_size_t), INTENT(IN)        :: node_tags(:)
    REAL(c_double), INTENT(IN)           :: coordinates(:)
    INTEGER, INTENT(IN)                  :: node_dimensions(:), node_entities(:)
    INTEGER(c_int), ALLOCATABLE          :: cad_entities(:,:)
    INTEGER(c_size_t), ALLOCATABLE       :: entity_node_tags(:)
    REAL(c_double), ALLOCATABLE          :: entity_coordinates(:), parameters(:)
    INTEGER                              :: dimension, entity, tag, count_nodes, node, local_node, added_nodes

    added_nodes = 0
    DO dimension = 0, 2
       CALL gmsh_l%model%getEntities(cad_entities, dim=dimension)
       DO entity = 1, SIZE(cad_entities,2)
          tag = INT(cad_entities(2,entity))
          count_nodes = COUNT(node_dimensions .EQ. dimension .AND. node_entities .EQ. tag)
          IF (count_nodes .EQ. 0) CYCLE
          ALLOCATE(entity_node_tags(count_nodes), entity_coordinates(3*count_nodes))
          local_node = 0
          DO node = 1, SIZE(node_tags)
             IF (node_dimensions(node) .NE. dimension .OR. node_entities(node) .NE. tag) CYCLE
             local_node = local_node + 1
             entity_node_tags(local_node) = node_tags(node)
             entity_coordinates(3*local_node-2:3*local_node) = coordinates(3*node-2:3*node)
          ENDDO

          IF (dimension .EQ. 0) THEN
             CALL gmsh_l%model%mesh%addNodes(dimension, tag, entity_node_tags, entity_coordinates)
          ELSE
             CALL gmsh_l%model%getParametrization(dimension, tag, entity_coordinates, parameters)
             CALL gmsh_l%model%mesh%addNodes(dimension, tag, entity_node_tags, entity_coordinates, &
                  parametricCoord=parameters)
             DEALLOCATE(parameters)
          ENDIF
          added_nodes = added_nodes + count_nodes
          DEALLOCATE(entity_node_tags, entity_coordinates)
       ENDDO
    ENDDO
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
    ENDDO
  END SUBROUTINE add_source_elements_to_cad

  SUBROUTINE mesh_order_error(message)
    CHARACTER(*), INTENT(IN) :: message

    WRITE(*,*) 'Unable to increase mesh order on the CAD geometry: ', TRIM(message)
    ERROR STOP 1
  END SUBROUTINE mesh_order_error

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
      TYPE(gmsh_t) :: gmsh_l
      INTEGER, INTENT(IN) :: p_order
      REAL*8, DIMENSION(:), INTENT(IN) :: h_target_on_elements ! on the elements
      REAL*8, DIMENSION(:,:), INTENT(IN) :: vertices_coordinates ! coordinates of the vertices
      INTEGER, DIMENSION(:,:), INTENT(IN) :: connectivity ! connectivity of the triangles
      REAL*8, ALLOCATABLE :: data_for_gmsh(:), h_target_on_nodes(:)
      INTEGER*8           :: number_of_vertices_per_triangle, i, j, start_index
      INTEGER*4           :: size_view,number_of_triangles,ret
      REAL*8              :: sf_index, vertex_coordinates(2)
      REAL*8, PARAMETER   :: minimum_mesh_size = 0.5d-4
      REAL*8, PARAMETER   :: target_size_quantum = 2.d-3*minimum_mesh_size
      REAL*8, PARAMETER   :: view_coordinate_quantum = 1.d-11
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

      !Set minimal size for mesh elements
      CALL gmsh_l%option%setNumber("Mesh.MeshSizeMin", minimum_mesh_size)


      ALLOCATE(data_for_gmsh((number_of_vertices_per_triangle*(number_of_vertices_per_triangle+1))*number_of_triangles))
      ALLOCATE(h_target_on_nodes(SIZE(vertices_coordinates,1)))
      ! Shared view nodes must carry one value for continuous linear interpolation.
      CALL average_element_targets_on_nodes(h_target_on_elements, connectivity, h_target_on_nodes)
      ! Solver roundoff must not select a different Gmsh topology.
      WHERE (h_target_on_nodes .GT. 0.d0)
         h_target_on_nodes = target_size_quantum*ANINT(h_target_on_nodes/target_size_quantum)
      END WHERE

      ! Prepare data in gmsh format
      DO i = 1, number_of_triangles
         start_index = (i-1)*(number_of_vertices_per_triangle*(number_of_vertices_per_triangle+1))
         DO j = 1, number_of_vertices_per_triangle
            vertex_coordinates = view_coordinate_quantum*ANINT(&
               vertices_coordinates(connectivity(i,j),:)/view_coordinate_quantum)
            data_for_gmsh(start_index+j) = vertex_coordinates(1)
            data_for_gmsh(start_index+number_of_vertices_per_triangle+j) = vertex_coordinates(2)
            data_for_gmsh(start_index+2*number_of_vertices_per_triangle+j) = 0.0 ! no Z coordinate
            data_for_gmsh(start_index+3*number_of_vertices_per_triangle+j) = &
               h_target_on_nodes(connectivity(i,j))
         ENDDO
      ENDDO

      ! Add h_target as a post-processing view
      size_view = gmsh_l%view%add("h_target")
      call gmsh_l%view%addListData(size_view, "ST",number_of_triangles,data_for_gmsh)
      sf_index = gmsh_l%view%getIndex(size_view)

      ! Add the view as a field
      ret = gmsh_l%model%mesh%field%add("PostView")
      call gmsh_l%model%mesh%field%setNumber(ret, "ViewIndex", sf_index)

      ! Apply the view as the current background mesh size field:
      call gmsh_l%model%mesh%field%setAsBackgroundMesh(ret)

      ! ignore characteristic length from geometry
      call gmsh_l%option%setNumber("Mesh.MeshSizeExtendFromBoundary", 0d0)
      call gmsh_l%option%setNumber("Mesh.MeshSizeFromPoints", 0d0)
      call gmsh_l%option%setNumber("Mesh.MeshSizeFromCurvature", 0d0)
      CALL gmsh_l%option%setNumber("Mesh.MeshSizeFactor", 1d0)
      ! Stabilize adaptive topology and avoid Gmsh's multi-threaded 2D meshing stall.
      CALL gmsh_l%option%setNumber("General.NumThreads", 1.d0)

      !Changing the algorithm to Delaunay, the default is Frontal-Delaunay (Don't Know if needed)
      call gmsh_l%option%setNumber("Mesh.Algorithm", 5d0)

      ! Generate the refined mesh
      call gmsh_l%model%mesh%generate(2)
      call gmsh_l%model%mesh%setOrder(p_order)
      call gmsh_l%model%mesh%optimize('HighOrder')
      CALL gmsh_l%option%setNumber("Mesh.MshFileVersion", 2.2)
      call gmsh_l%write('./res/temp.msh')
      CALL gmsh_l%finalize()

      DEALLOCATE(data_for_gmsh, h_target_on_nodes)

   END SUBROUTINE gmsh_create_from_h_target

   SUBROUTINE save_copy_new_mesh(mesh_name, count_adapt)
      USE in_out, ONLY: copy_file
      CHARACTER(1024), INTENT(IN) :: mesh_name
      INTEGER, INTENT(IN)         :: count_adapt
      CHARACTER(70)               :: count_adapt_char
      CHARACTER(1024)             :: mesh_name_npne,new_mesh_name_npne, buffer

      CALL extract_mesh_name_from_fullpath_woext(mesh_name, mesh_name_npne)

      WRITE(count_adapt_char, *) count_adapt
      new_mesh_name_npne = TRIM(ADJUSTL(mesh_name_npne)) // '_n' // TRIM(ADJUSTL(count_adapt_char))
      
      buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".msh"
      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*,*) "Mesh saved as: ", TRIM(ADJUSTL(new_mesh_name_npne))
      ENDIF
      CALL copy_file("./res/temp.msh",buffer)

   END SUBROUTINE save_copy_new_mesh
END MODULE adaptivity_common_module
