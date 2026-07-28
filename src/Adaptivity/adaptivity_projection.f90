MODULE adaptivity_projection_module
  USE globals, ONLY: phys, refElPol, utils
  USE MPI_OMP
  USE element_mapping_module, ONLY: map_physical_to_nearest_reference, map_physical_to_reference
  USE reference_element, ONLY: compute_shape_functions_at_points

  IMPLICIT NONE
  PRIVATE

  REAL*8, PARAMETER :: linear_candidate_tolerance = 1.d-10
  REAL*8, PARAMETER :: curved_owner_tolerance = 1.d-8
  REAL*8, PARAMETER :: recovery_distance_ratio = 1.d-2
  REAL*8, PARAMETER :: recovery_box_padding_ratio = 0.5d0

  PUBLIC :: projectSolutionDifferentMeshes_general
  PUBLIC :: projectSolutionDifferentMeshes_general_arrays

CONTAINS

  SUBROUTINE projectSolutionDifferentMeshes_general(T1, X1, T2, X2, u1, q1, u2, q2)
    INTEGER, INTENT(IN)                                    :: T1(:,:), T2(:,:)
    REAL*8, INTENT(IN)                                     :: X1(:,:), X2(:,:)
    REAL*8, POINTER, INTENT(IN)                            :: u1(:)
    REAL*8, POINTER, OPTIONAL, INTENT(IN)                  :: q1(:)
    REAL*8, POINTER, INTENT(INOUT)                         :: u2(:)
    REAL*8, POINTER, OPTIONAL, INTENT(INOUT)               :: q2(:)

    REAL*8, ALLOCATABLE                                    :: u_old(:,:), u_new(:,:)
    REAL*8, ALLOCATABLE                                    :: q_old(:,:,:), q_new(:,:,:)
    INTEGER                                                :: old_points, new_points, gradient_dimensions

    CALL validate_optional_gradient_pair(PRESENT(q1), PRESENT(q2))
    CALL validate_mesh_geometry(T1, X1, T2, X2)
    IF(.NOT. ASSOCIATED(u1)) ERROR STOP "Old projection solution is not associated"

    old_points = SIZE(T1,1)*SIZE(T1,2)
    new_points = SIZE(T2,1)*SIZE(T2,2)
    gradient_dimensions = SIZE(X1,2)
    IF(SIZE(u1) .NE. old_points*phys%neq) ERROR STOP "Old projection solution has an invalid size"

    CALL resize_real_pointer(u2, new_points*phys%neq)
    ALLOCATE(u_old(old_points,phys%neq), u_new(new_points,phys%neq))
    CALL unpack_flat_solution(u1, u_old)

    IF(PRESENT(q1)) THEN
       IF(.NOT. ASSOCIATED(q1)) ERROR STOP "Old projection gradient is not associated"
       IF(SIZE(q1) .NE. old_points*phys%neq*gradient_dimensions) THEN
          ERROR STOP "Old projection gradient has an invalid size"
       ENDIF

       CALL resize_real_pointer(q2, new_points*phys%neq*gradient_dimensions)
       ALLOCATE(q_old(old_points,phys%neq,gradient_dimensions))
       ALLOCATE(q_new(new_points,phys%neq,gradient_dimensions))
       CALL unpack_flat_gradient(q1, q_old)
       CALL project_solution_between_meshes(T1, X1, T2, X2, &
            u_old, q_old, u_new, q_new)
       CALL pack_flat_gradient(q_new, q2)
       DEALLOCATE(q_old, q_new)
    ELSE
       CALL project_solution_between_meshes(T1, X1, T2, X2, &
            u_old=u_old, u_new=u_new)
    ENDIF

    CALL pack_flat_solution(u_new, u2)
    DEALLOCATE(u_old, u_new)
  ENDSUBROUTINE projectSolutionDifferentMeshes_general

  SUBROUTINE projectSolutionDifferentMeshes_general_arrays(T1, X1, T2, X2, u1, q1, u2, q2)
    INTEGER, INTENT(IN)                    :: T1(:,:), T2(:,:)
    REAL*8, INTENT(IN)                     :: X1(:,:), X2(:,:)
    REAL*8, INTENT(IN)                     :: u1(:)
    REAL*8, OPTIONAL, INTENT(IN)           :: q1(:)
    REAL*8, INTENT(INOUT)                  :: u2(:)
    REAL*8, OPTIONAL, INTENT(INOUT)        :: q2(:)

    REAL*8, ALLOCATABLE                    :: u_old(:,:), u_new(:,:)
    REAL*8, ALLOCATABLE                    :: q_old(:,:,:), q_new(:,:,:)
    INTEGER                                :: old_points, new_points, gradient_dimensions

    CALL validate_optional_gradient_pair(PRESENT(q1), PRESENT(q2))
    CALL validate_mesh_geometry(T1, X1, T2, X2)

    old_points = SIZE(T1,1)*SIZE(T1,2)
    new_points = SIZE(T2,1)*SIZE(T2,2)
    gradient_dimensions = SIZE(X1,2)
    IF(SIZE(u1) .NE. old_points*phys%neq .OR. SIZE(u2) .NE. new_points*phys%neq) THEN
       ERROR STOP "Projection solution arrays have invalid sizes"
    ENDIF

    ALLOCATE(u_old(old_points,phys%neq), u_new(new_points,phys%neq))
    CALL unpack_flat_solution(u1, u_old)

    IF(PRESENT(q1)) THEN
       IF(SIZE(q1) .NE. old_points*phys%neq*gradient_dimensions .OR. &
            SIZE(q2) .NE. new_points*phys%neq*gradient_dimensions) THEN
          ERROR STOP "Projection gradient arrays have invalid sizes"
       ENDIF

       ALLOCATE(q_old(old_points,phys%neq,gradient_dimensions))
       ALLOCATE(q_new(new_points,phys%neq,gradient_dimensions))
       CALL unpack_flat_gradient(q1, q_old)
       CALL project_solution_between_meshes(T1, X1, T2, X2, &
            u_old, q_old, u_new, q_new)
       CALL pack_flat_gradient(q_new, q2)
       DEALLOCATE(q_old, q_new)
    ELSE
       CALL project_solution_between_meshes(T1, X1, T2, X2, &
            u_old=u_old, u_new=u_new)
    ENDIF

    CALL pack_flat_solution(u_new, u2)
    DEALLOCATE(u_old, u_new)
  ENDSUBROUTINE projectSolutionDifferentMeshes_general_arrays

  SUBROUTINE project_solution_between_meshes(old_connectivity, old_coordinates, new_connectivity, new_coordinates, &
       u_old, q_old, u_new, q_new)
    INTEGER, INTENT(IN)                    :: old_connectivity(:,:), new_connectivity(:,:)
    REAL*8, INTENT(IN)                     :: old_coordinates(:,:), new_coordinates(:,:)
    REAL*8, INTENT(IN)                     :: u_old(:,:)
    REAL*8, OPTIONAL, INTENT(IN)           :: q_old(:,:,:)
    REAL*8, INTENT(OUT)                    :: u_new(:,:)
    REAL*8, OPTIONAL, INTENT(OUT)          :: q_new(:,:,:)

    REAL*8                                 :: target_points(SIZE(new_connectivity,1)*SIZE(new_connectivity,2),2)
    REAL*8                                 :: reference_points(SIZE(target_points,1),2)
    INTEGER                                :: point_elements(SIZE(target_points,1))
    INTEGER                                :: nearest_elements(SIZE(target_points,1))
    REAL*8                                 :: nearest_distances(SIZE(target_points,1))
    REAL*8                                 :: nearest_sizes(SIZE(target_points,1))
    INTEGER                                :: global_missing

    CALL validate_projection_shapes(old_connectivity, new_connectivity, u_old, q_old, u_new, q_new)
    CALL collect_target_points(new_connectivity, new_coordinates, target_points)
    CALL print_projection_banner()

    point_elements = 0
    reference_points = HUGE(1.d0)
    CALL assign_linear_projection_candidates(target_points, old_connectivity, old_coordinates, point_elements, &
         reference_points)
    IF(ANY(point_elements .EQ. 0)) THEN
       CALL assign_curved_projection_candidates(target_points, old_connectivity, old_coordinates, point_elements, &
            reference_points)
    ENDIF

    nearest_elements = 0
    nearest_distances = HUGE(1.d0)
    nearest_sizes = 0.d0
    IF(ANY(point_elements .EQ. 0)) THEN
       CALL recover_nearest_projection_points(target_points, old_connectivity, old_coordinates, point_elements, &
            reference_points, nearest_elements, nearest_distances, nearest_sizes)
    ENDIF

    CALL report_unmatched_projection_points(target_points, point_elements, nearest_elements, nearest_distances, &
         nearest_sizes, global_missing)
    IF(global_missing .NE. 0) ERROR STOP "Adaptive projection has unmatched points"

    CALL interpolate_projected_solution(reference_points, point_elements, old_connectivity, u_old, q_old, u_new, q_new)
  ENDSUBROUTINE project_solution_between_meshes

  SUBROUTINE assign_linear_projection_candidates(target_points, old_connectivity, old_coordinates, point_elements, &
       reference_points)
    REAL*8, INTENT(IN)                     :: target_points(:,:), old_coordinates(:,:)
    INTEGER, INTENT(IN)                    :: old_connectivity(:,:)
    INTEGER, INTENT(INOUT)                 :: point_elements(:)
    REAL*8, INTENT(INOUT)                  :: reference_points(:,:)

    REAL*8                                 :: element_coordinates(SIZE(old_connectivity,2),2)
    REAL*8                                 :: triangle_vertices(3,2), candidate_reference(2)
    LOGICAL                                :: owned
    INTEGER                                :: element, point

    DO element = 1, SIZE(old_connectivity,1)
       element_coordinates = old_coordinates(old_connectivity(element,:),:)
       triangle_vertices = element_coordinates(1:3,:)
       DO point = 1, SIZE(target_points,1)
          IF(point_elements(point) .NE. 0) CYCLE
          IF(.NOT. point_is_in_linear_triangle(target_points(point,:), triangle_vertices, &
               linear_candidate_tolerance)) CYCLE

          CALL confirm_curved_element_ownership(target_points(point,:), element_coordinates, candidate_reference, owned)
          IF(owned) THEN
             point_elements(point) = element
             reference_points(point,:) = candidate_reference
          ENDIF
       ENDDO
    ENDDO
  ENDSUBROUTINE assign_linear_projection_candidates

  SUBROUTINE assign_curved_projection_candidates(target_points, old_connectivity, old_coordinates, point_elements, &
       reference_points)
    REAL*8, INTENT(IN)                     :: target_points(:,:), old_coordinates(:,:)
    INTEGER, INTENT(IN)                    :: old_connectivity(:,:)
    INTEGER, INTENT(INOUT)                 :: point_elements(:)
    REAL*8, INTENT(INOUT)                  :: reference_points(:,:)

    REAL*8                                 :: element_coordinates(SIZE(old_connectivity,2),2)
    REAL*8                                 :: candidate_reference(2), xmin, xmax, ymin, ymax, bounding_scale, padding
    LOGICAL                                :: owned
    INTEGER                                :: element, point

    DO point = 1, SIZE(target_points,1)
       IF(point_elements(point) .NE. 0) CYCLE
       DO element = 1, SIZE(old_connectivity,1)
          element_coordinates = old_coordinates(old_connectivity(element,:),:)
          CALL curved_element_bounding_box(element_coordinates, xmin, xmax, ymin, ymax, bounding_scale)
          padding = MAX(1.d-12, curved_owner_tolerance*MAX(1.d0,bounding_scale))
          IF(.NOT. point_is_in_padded_box(target_points(point,:), xmin, xmax, ymin, ymax, padding)) CYCLE

          CALL confirm_curved_element_ownership(target_points(point,:), element_coordinates, candidate_reference, owned)
          IF(owned) THEN
             point_elements(point) = element
             reference_points(point,:) = candidate_reference
             EXIT
          ENDIF
       ENDDO
    ENDDO
  ENDSUBROUTINE assign_curved_projection_candidates

  SUBROUTINE confirm_curved_element_ownership(target_point, element_coordinates, reference_point, owned)
    REAL*8, INTENT(IN)                     :: target_point(2), element_coordinates(:,:)
    REAL*8, INTENT(OUT)                    :: reference_point(2)
    LOGICAL, INTENT(OUT)                   :: owned

    REAL*8                                 :: target_points(1,2), reference_points(1,2)
    LOGICAL                                :: converged(1)

    target_points(1,:) = target_point
    CALL map_physical_to_reference(target_points, element_coordinates, refElPol, reference_points, converged)
    reference_point = reference_points(1,:)
    owned = converged(1) .AND. reference_point_is_in_triangle(reference_point, curved_owner_tolerance)
  ENDSUBROUTINE confirm_curved_element_ownership

  SUBROUTINE recover_nearest_projection_points(target_points, old_connectivity, old_coordinates, point_elements, &
       reference_points, nearest_elements, nearest_distances, nearest_sizes)
    REAL*8, INTENT(IN)                     :: target_points(:,:), old_coordinates(:,:)
    INTEGER, INTENT(IN)                    :: old_connectivity(:,:)
    INTEGER, INTENT(INOUT)                 :: point_elements(:), nearest_elements(:)
    REAL*8, INTENT(INOUT)                  :: reference_points(:,:), nearest_distances(:), nearest_sizes(:)

    REAL*8                                 :: nearest_reference(2), nearest_mapped(2), nearest_distance, nearest_size
    LOGICAL                                :: nearest_valid
    INTEGER                                :: nearest_element, point

    DO point = 1, SIZE(target_points,1)
       IF(point_elements(point) .NE. 0) CYCLE

       CALL find_nearest_old_element(target_points(point,:), old_connectivity, old_coordinates, .TRUE., &
            nearest_element, nearest_reference, nearest_mapped, nearest_distance, nearest_size, nearest_valid)
       IF(.NOT. nearest_valid) THEN
          CALL find_nearest_old_element(target_points(point,:), old_connectivity, old_coordinates, .FALSE., &
               nearest_element, nearest_reference, nearest_mapped, nearest_distance, nearest_size, nearest_valid)
       ENDIF

       IF(nearest_valid) THEN
          nearest_elements(point) = nearest_element
          nearest_distances(point) = nearest_distance
          nearest_sizes(point) = nearest_size
       ENDIF
       IF(nearest_valid .AND. nearest_distance .LE. recovery_distance_ratio*MAX(1.d-12,nearest_size)) THEN
          point_elements(point) = nearest_element
          reference_points(point,:) = nearest_reference
       ENDIF
    ENDDO
  ENDSUBROUTINE recover_nearest_projection_points

  SUBROUTINE find_nearest_old_element(target_point, old_connectivity, old_coordinates, padded_candidates_only, &
       nearest_element, nearest_reference, nearest_mapped, nearest_distance, nearest_size, valid)
    REAL*8, INTENT(IN)                     :: target_point(2), old_coordinates(:,:)
    INTEGER, INTENT(IN)                    :: old_connectivity(:,:)
    LOGICAL, INTENT(IN)                    :: padded_candidates_only
    INTEGER, INTENT(OUT)                   :: nearest_element
    REAL*8, INTENT(OUT)                    :: nearest_reference(2), nearest_mapped(2), nearest_distance, nearest_size
    LOGICAL, INTENT(OUT)                   :: valid

    REAL*8                                 :: element_coordinates(SIZE(old_connectivity,2),2)
    REAL*8                                 :: candidate_reference(2), candidate_mapped(2), candidate_distance
    REAL*8                                 :: candidate_size, xmin, xmax, ymin, ymax, bounding_scale, padding
    LOGICAL                                :: candidate_valid
    INTEGER                                :: element

    valid = .FALSE.
    nearest_element = 0
    nearest_reference = 0.d0
    nearest_mapped = 0.d0
    nearest_distance = HUGE(1.d0)
    nearest_size = 0.d0

    DO element = 1, SIZE(old_connectivity,1)
       element_coordinates = old_coordinates(old_connectivity(element,:),:)
       CALL curved_element_bounding_box(element_coordinates, xmin, xmax, ymin, ymax, bounding_scale)
       candidate_size = linear_triangle_size(element_coordinates(1:3,:))
       padding = recovery_box_padding_ratio*MAX(bounding_scale,candidate_size)
       IF(padded_candidates_only .AND. &
            .NOT. point_is_in_padded_box(target_point, xmin, xmax, ymin, ymax, padding)) CYCLE

       CALL map_physical_to_nearest_reference(target_point, element_coordinates, refElPol, candidate_reference, &
            candidate_mapped, candidate_distance, candidate_valid)
       CALL retain_nearest_element(element, candidate_reference, candidate_mapped, candidate_distance, candidate_size, &
            candidate_valid, nearest_element, nearest_reference, nearest_mapped, nearest_distance, nearest_size, valid)
    ENDDO
  ENDSUBROUTINE find_nearest_old_element

  SUBROUTINE retain_nearest_element(candidate_element, candidate_reference, candidate_mapped, candidate_distance, &
       candidate_size, candidate_valid, nearest_element, nearest_reference, nearest_mapped, nearest_distance, &
       nearest_size, valid)
    INTEGER, INTENT(IN)                    :: candidate_element
    REAL*8, INTENT(IN)                     :: candidate_reference(2), candidate_mapped(2), candidate_distance
    REAL*8, INTENT(IN)                     :: candidate_size
    LOGICAL, INTENT(IN)                    :: candidate_valid
    INTEGER, INTENT(INOUT)                 :: nearest_element
    REAL*8, INTENT(INOUT)                  :: nearest_reference(2), nearest_mapped(2), nearest_distance, nearest_size
    LOGICAL, INTENT(INOUT)                 :: valid

    REAL*8                                 :: tie_tolerance

    IF(.NOT. candidate_valid) RETURN
    tie_tolerance = 100.d0*EPSILON(1.d0)*MAX(1.d0,candidate_distance,nearest_distance)
    IF(valid .AND. candidate_distance .GT. nearest_distance-tie_tolerance) THEN
       IF(ABS(candidate_distance-nearest_distance) .GT. tie_tolerance .OR. candidate_element .GE. nearest_element) RETURN
    ENDIF

    valid = .TRUE.
    nearest_element = candidate_element
    nearest_reference = candidate_reference
    nearest_mapped = candidate_mapped
    nearest_distance = candidate_distance
    nearest_size = candidate_size
  ENDSUBROUTINE retain_nearest_element

  SUBROUTINE report_unmatched_projection_points(target_points, point_elements, nearest_elements, nearest_distances, &
       nearest_sizes, global_missing)
    REAL*8, INTENT(IN)                     :: target_points(:,:), nearest_distances(:), nearest_sizes(:)
    INTEGER, INTENT(IN)                    :: point_elements(:), nearest_elements(:)
    INTEGER, INTENT(OUT)                   :: global_missing

    INTEGER                                :: local_missing, point, ierr

    local_missing = COUNT(point_elements .EQ. 0)
    global_missing = local_missing
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, global_missing, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF(global_missing .EQ. 0) RETURN

    IF(MPIvar%glob_id .EQ. 0) WRITE(*,*) "Projection unmatched points after recovery: ", global_missing
    IF(local_missing .EQ. 0) RETURN

    DO point = 1, SIZE(target_points,1)
       IF(point_elements(point) .NE. 0) CYCLE
       IF(nearest_elements(point) .NE. 0) THEN
          WRITE(*,*) "Unmatched projection point rank/index/xy/nearest/dist/local h/ratio: ", &
               MPIvar%glob_id, point, target_points(point,1), target_points(point,2), nearest_elements(point), &
               nearest_distances(point), nearest_sizes(point), &
               nearest_distances(point)/MAX(1.d-12,nearest_sizes(point))
       ELSE
          WRITE(*,*) "Unmatched projection point rank/index/xy/no nearest element: ", &
               MPIvar%glob_id, point, target_points(point,1), target_points(point,2)
       ENDIF
    ENDDO
  ENDSUBROUTINE report_unmatched_projection_points

  SUBROUTINE interpolate_projected_solution(reference_points, point_elements, old_connectivity, u_old, q_old, &
       u_new, q_new)
    REAL*8, INTENT(IN)                     :: reference_points(:,:), u_old(:,:)
    INTEGER, INTENT(IN)                    :: point_elements(:), old_connectivity(:,:)
    REAL*8, OPTIONAL, INTENT(IN)           :: q_old(:,:,:)
    REAL*8, INTENT(OUT)                    :: u_new(:,:)
    REAL*8, OPTIONAL, INTENT(OUT)          :: q_new(:,:,:)

    INTEGER                                :: old_dofs(SIZE(old_connectivity,2))
    REAL*8, ALLOCATABLE                    :: element_reference_points(:,:), shape_functions(:,:,:)
    INTEGER, ALLOCATABLE                   :: point_indices(:)
    INTEGER                                :: element, point, direction

    u_new = 0.d0
    IF(PRESENT(q_new)) q_new = 0.d0
    DO element = 1, SIZE(old_connectivity,1)
       CALL find_matching_indices(point_elements, element, point_indices)
       IF(SIZE(point_indices) .EQ. 0) THEN
          DEALLOCATE(point_indices)
          CYCLE
       ENDIF

       ALLOCATE(element_reference_points(SIZE(point_indices),2))
       ALLOCATE(shape_functions(SIZE(old_connectivity,2),SIZE(point_indices),3))
       element_reference_points = reference_points(point_indices,:)
       shape_functions = 0.d0
       CALL compute_shape_functions_at_points(refElPol, element_reference_points, shape_functions)

       old_dofs = (element-1)*SIZE(old_connectivity,2) + (/ (point, point=1,SIZE(old_connectivity,2)) /)
       u_new(point_indices,:) = MATMUL(TRANSPOSE(shape_functions(:,:,1)),u_old(old_dofs,:))
       IF(PRESENT(q_old)) THEN
          DO direction = 1, SIZE(q_old,3)
             q_new(point_indices,:,direction) = &
                  MATMUL(TRANSPOSE(shape_functions(:,:,1)),q_old(old_dofs,:,direction))
          ENDDO
       ENDIF

       DEALLOCATE(element_reference_points, shape_functions, point_indices)
    ENDDO
  ENDSUBROUTINE interpolate_projected_solution

  SUBROUTINE collect_target_points(connectivity, coordinates, target_points)
    INTEGER, INTENT(IN)                    :: connectivity(:,:)
    REAL*8, INTENT(IN)                     :: coordinates(:,:)
    REAL*8, INTENT(OUT)                    :: target_points(:,:)

    INTEGER                                :: element, first_point, last_point

    DO element = 1, SIZE(connectivity,1)
       first_point = (element-1)*SIZE(connectivity,2)+1
       last_point = element*SIZE(connectivity,2)
       target_points(first_point:last_point,:) = coordinates(connectivity(element,:),:)
    ENDDO
  ENDSUBROUTINE collect_target_points

  PURE LOGICAL FUNCTION point_is_in_linear_triangle(point, vertices, tolerance)
    REAL*8, INTENT(IN)                     :: point(2), vertices(3,2), tolerance

    REAL*8                                 :: edge1(2), edge2(2), offset(2), determinant, scale
    REAL*8                                 :: barycentric(3)

    edge1 = vertices(2,:)-vertices(1,:)
    edge2 = vertices(3,:)-vertices(1,:)
    offset = point-vertices(1,:)
    determinant = edge1(1)*edge2(2)-edge1(2)*edge2(1)
    scale = MAX(MAXVAL(ABS(edge1)),MAXVAL(ABS(edge2)))
    IF(scale .LE. TINY(1.d0) .OR. ABS(determinant) .LE. 100.d0*EPSILON(1.d0)*scale**2) THEN
       point_is_in_linear_triangle = .FALSE.
       RETURN
    ENDIF

    barycentric(2) = (offset(1)*edge2(2)-offset(2)*edge2(1))/determinant
    barycentric(3) = (edge1(1)*offset(2)-edge1(2)*offset(1))/determinant
    barycentric(1) = 1.d0-barycentric(2)-barycentric(3)
    point_is_in_linear_triangle = ALL(barycentric .GE. -tolerance) .AND. &
         ALL(barycentric .LE. 1.d0+tolerance)
  ENDFUNCTION point_is_in_linear_triangle

  PURE LOGICAL FUNCTION reference_point_is_in_triangle(point, tolerance)
    REAL*8, INTENT(IN)                     :: point(2), tolerance

    reference_point_is_in_triangle = point(1) .GE. -1.d0-tolerance .AND. &
         point(2) .GE. -1.d0-tolerance .AND. point(1)+point(2) .LE. tolerance
  ENDFUNCTION reference_point_is_in_triangle

  PURE SUBROUTINE curved_element_bounding_box(element_coordinates, xmin, xmax, ymin, ymax, scale)
    REAL*8, INTENT(IN)                     :: element_coordinates(:,:)
    REAL*8, INTENT(OUT)                    :: xmin, xmax, ymin, ymax, scale

    xmin = MINVAL(element_coordinates(:,1))
    xmax = MAXVAL(element_coordinates(:,1))
    ymin = MINVAL(element_coordinates(:,2))
    ymax = MAXVAL(element_coordinates(:,2))
    scale = MAX(xmax-xmin,ymax-ymin)
  ENDSUBROUTINE curved_element_bounding_box

  PURE LOGICAL FUNCTION point_is_in_padded_box(point, xmin, xmax, ymin, ymax, padding)
    REAL*8, INTENT(IN)                     :: point(2), xmin, xmax, ymin, ymax, padding

    point_is_in_padded_box = point(1) .GE. xmin-padding .AND. point(1) .LE. xmax+padding .AND. &
         point(2) .GE. ymin-padding .AND. point(2) .LE. ymax+padding
  ENDFUNCTION point_is_in_padded_box

  PURE REAL*8 FUNCTION linear_triangle_size(vertices)
    REAL*8, INTENT(IN)                     :: vertices(3,2)

    linear_triangle_size = MIN(SQRT(SUM((vertices(1,:)-vertices(2,:))**2)), &
         SQRT(SUM((vertices(1,:)-vertices(3,:))**2)), SQRT(SUM((vertices(2,:)-vertices(3,:))**2)))
  ENDFUNCTION linear_triangle_size

  PURE SUBROUTINE find_matching_indices(values, target_value, indices)
    INTEGER, INTENT(IN)                    :: values(:), target_value
    INTEGER, ALLOCATABLE, INTENT(OUT)      :: indices(:)
    INTEGER                                :: point, match

    ALLOCATE(indices(COUNT(values .EQ. target_value)))
    match = 1
    DO point = 1, SIZE(values)
       IF(values(point) .EQ. target_value) THEN
          indices(match) = point
          match = match+1
       ENDIF
    ENDDO
  ENDSUBROUTINE find_matching_indices

  SUBROUTINE unpack_flat_solution(flat_values, values)
    REAL*8, INTENT(IN)                     :: flat_values(:)
    REAL*8, INTENT(OUT)                    :: values(:,:)
    INTEGER                                :: point, first_value, last_value

    DO point = 1, SIZE(values,1)
       first_value = (point-1)*SIZE(values,2)+1
       last_value = point*SIZE(values,2)
       values(point,:) = flat_values(first_value:last_value)
    ENDDO
  ENDSUBROUTINE unpack_flat_solution

  SUBROUTINE pack_flat_solution(values, flat_values)
    REAL*8, INTENT(IN)                     :: values(:,:)
    REAL*8, INTENT(OUT)                    :: flat_values(:)
    INTEGER                                :: point, first_value, last_value

    DO point = 1, SIZE(values,1)
       first_value = (point-1)*SIZE(values,2)+1
       last_value = point*SIZE(values,2)
       flat_values(first_value:last_value) = values(point,:)
    ENDDO
  ENDSUBROUTINE pack_flat_solution

  SUBROUTINE unpack_flat_gradient(flat_values, values)
    REAL*8, INTENT(IN)                     :: flat_values(:)
    REAL*8, INTENT(OUT)                    :: values(:,:,:)
    INTEGER                                :: point, equation, direction, counter

    counter = 1
    DO point = 1, SIZE(values,1)
       DO equation = 1, SIZE(values,2)
          DO direction = 1, SIZE(values,3)
             values(point,equation,direction) = flat_values(counter)
             counter = counter+1
          ENDDO
       ENDDO
    ENDDO
  ENDSUBROUTINE unpack_flat_gradient

  SUBROUTINE pack_flat_gradient(values, flat_values)
    REAL*8, INTENT(IN)                     :: values(:,:,:)
    REAL*8, INTENT(OUT)                    :: flat_values(:)
    INTEGER                                :: point, equation, direction, counter

    counter = 1
    DO point = 1, SIZE(values,1)
       DO equation = 1, SIZE(values,2)
          DO direction = 1, SIZE(values,3)
             flat_values(counter) = values(point,equation,direction)
             counter = counter+1
          ENDDO
       ENDDO
    ENDDO
  ENDSUBROUTINE pack_flat_gradient

  SUBROUTINE resize_real_pointer(values, required_size)
    REAL*8, POINTER, INTENT(INOUT)          :: values(:)
    INTEGER, INTENT(IN)                    :: required_size

    IF(ASSOCIATED(values)) THEN
       IF(SIZE(values) .EQ. required_size) RETURN
       DEALLOCATE(values)
    ENDIF
    ALLOCATE(values(required_size))
  ENDSUBROUTINE resize_real_pointer

  SUBROUTINE validate_optional_gradient_pair(old_present, new_present)
    LOGICAL, INTENT(IN)                    :: old_present, new_present

    IF(old_present .NEQV. new_present) ERROR STOP "Projection requires q_old and q_new together"
  ENDSUBROUTINE validate_optional_gradient_pair

  SUBROUTINE validate_mesh_geometry(old_connectivity, old_coordinates, new_connectivity, new_coordinates)
    INTEGER, INTENT(IN)                    :: old_connectivity(:,:), new_connectivity(:,:)
    REAL*8, INTENT(IN)                     :: old_coordinates(:,:), new_coordinates(:,:)

    IF(SIZE(old_coordinates,2) .NE. 2 .OR. SIZE(new_coordinates,2) .NE. 2) THEN
       ERROR STOP "Adaptive projection supports only two-dimensional meshes"
    ENDIF
    IF(SIZE(old_connectivity,1) .EQ. 0 .OR. SIZE(new_connectivity,1) .EQ. 0) THEN
       ERROR STOP "Adaptive projection requires nonempty meshes"
    ENDIF
    IF(SIZE(old_connectivity,2) .NE. SIZE(new_connectivity,2)) THEN
       ERROR STOP "Projection between differing polynomial orders is unsupported"
    ENDIF
    IF(SIZE(old_connectivity,2) .NE. refElPol%Nnodes2D) THEN
       ERROR STOP "Projection connectivity does not match the active polynomial order"
    ENDIF
    IF(ANY(old_connectivity .LT. 1) .OR. ANY(old_connectivity .GT. SIZE(old_coordinates,1)) .OR. &
         ANY(new_connectivity .LT. 1) .OR. ANY(new_connectivity .GT. SIZE(new_coordinates,1))) THEN
       ERROR STOP "Projection connectivity contains invalid node indices"
    ENDIF
  ENDSUBROUTINE validate_mesh_geometry

  SUBROUTINE validate_projection_shapes(old_connectivity, new_connectivity, u_old, q_old, u_new, q_new)
    INTEGER, INTENT(IN)                    :: old_connectivity(:,:), new_connectivity(:,:)
    REAL*8, INTENT(IN)                     :: u_old(:,:)
    REAL*8, OPTIONAL, INTENT(IN)           :: q_old(:,:,:)
    REAL*8, INTENT(OUT)                    :: u_new(:,:)
    REAL*8, OPTIONAL, INTENT(OUT)          :: q_new(:,:,:)

    INTEGER                                :: old_points, new_points

    CALL validate_optional_gradient_pair(PRESENT(q_old),PRESENT(q_new))
    old_points = SIZE(old_connectivity,1)*SIZE(old_connectivity,2)
    new_points = SIZE(new_connectivity,1)*SIZE(new_connectivity,2)
    IF(SIZE(u_old,1) .NE. old_points .OR. SIZE(u_new,1) .NE. new_points .OR. &
         SIZE(u_old,2) .NE. SIZE(u_new,2)) ERROR STOP "Projection solution shapes are inconsistent"
    IF(PRESENT(q_old)) THEN
       IF(SIZE(q_old,1) .NE. old_points .OR. SIZE(q_new,1) .NE. new_points .OR. &
            SIZE(q_old,2) .NE. SIZE(u_old,2) .OR. SIZE(q_new,2) .NE. SIZE(u_new,2) .OR. &
            SIZE(q_old,3) .NE. SIZE(q_new,3)) ERROR STOP "Projection gradient shapes are inconsistent"
    ENDIF
  ENDSUBROUTINE validate_projection_shapes

  SUBROUTINE print_projection_banner()
    IF(MPIvar%glob_id .NE. 0 .OR. utils%printint .LE. 0) RETURN
    WRITE(6,*) '*************************************************'
    WRITE(6,*) '*      PROJECTING SOLUTION TO NEW MESH          *'
    WRITE(6,*) '*************************************************'
  ENDSUBROUTINE print_projection_banner

END MODULE adaptivity_projection_module
