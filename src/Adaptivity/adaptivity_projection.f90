MODULE adaptivity_projection_module
  USE globals
  USE MPI_OMP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: projectSolutionDifferentMeshes_general
  PUBLIC :: projectSolutionDifferentMeshes_general_arrays

CONTAINS

  SUBROUTINE projectSolutionDifferentMeshes_general(T1, X1, T2, X2, u1, q1, u2, q2)

    INTEGER, INTENT(IN)                                    :: T1(:,:), T2(:,:)
    REAL*8, INTENT(IN)                                     :: X1(:,:), X2(:,:)
    REAL*8, POINTER, DIMENSION(:), INTENT(IN)              :: u1(:)
    REAL*8, POINTER, DIMENSION(:), OPTIONAL, INTENT(IN)    :: q1(:)
    REAL*8, POINTER, DIMENSION(:), INTENT(INOUT)           :: u2(:)
    REAL*8, POINTER, DIMENSION(:), OPTIONAL, INTENT(INOUT) :: q2(:)
    REAL*8, ALLOCATABLE, DIMENSION(:,:)                    :: u1_2D, u2_2D
    REAL*8, ALLOCATABLE, DIMENSION(:,:,:)                  :: q1_3D, q2_3D
    INTEGER                                                :: i,j,k, counter


    ALLOCATE(u1_2D(SIZE(T1,1)*SIZE(T1,2), phys%neq))
    ALLOCATE(u2_2D(SIZE(T2,1)*SIZE(T2,2), phys%neq))

    IF(PRESENT(q1)) THEN
       ALLOCATE(q1_3D(SIZE(T1,1)*SIZE(T1,2), phys%neq, refElPol%ndim))
       ALLOCATE(q2_3D(SIZE(T2,1)*SIZE(T2,2), phys%neq, refElPol%ndim))
       q1_3D = 0.
       q2_3D = 0.
    ENDIF

    u1_2D = 0.
    u2_2D = 0.

    counter = 1
    ! reshape u sol
    DO i = 1, SIZE(u1_2D,1)
       DO j = 1, SIZE(u1_2D,2)
          u1_2D(i,j) = u1(counter)
          counter = counter + 1
       ENDDO
    ENDDO

    IF(PRESENT(q1)) THEN
       counter = 1
       ! q1_3D not as easy
       DO i = 1, SIZE(q1_3D,1)
          DO j = 1, SIZE(q1_3D,2)
             DO k = 1, SIZE(q1_3D,3)
                q1_3D(i,j,k) = q1(counter)
                counter = counter + 1
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF(PRESENT(q1)) THEN
       CALL projectSolutionDifferentMeshes_Mod(T1, X1, T2, X2, u1_2D, q1_3D, u2_2D, q2_3D)
    ELSE
       CALL projectSolutionDifferentMeshes_Mod(T1, X1, T2, X2, u_old = u1_2D, u_new=u2_2D)
    ENDIF

    IF(SIZE(u2) .NE. SIZE(u2_2D)) THEN
       DEALLOCATE(u2)
       ALLOCATE(u2(SIZE(u2_2D)))
    ENDIF

    IF(PRESENT(q1)) THEN
       IF(SIZE(q2) .NE. SIZE(q2_3D)) THEN
          DEALLOCATE(q2)
          ALLOCATE(q2(SIZE(q2_3D)))
       ENDIF
    ENDIF
    ! solu is easy to reshape
    counter = 1
    DO i = 1, SIZE(u2_2D,1)
       DO j = 1, SIZE(u2_2D,2)
          u2(counter)   = u2_2D(i,j)
          counter = counter + 1
       ENDDO
    ENDDO

    IF(PRESENT(q1)) THEN
       ! solq is not as easy
       counter = 1
       DO i = 1, SIZE(q2_3D,1)
          DO j = 1, SIZE(q2_3D,2)
             DO k = 1, SIZE(q2_3D,3)
                q2(counter)   = q2_3D(i,j,k)
                counter = counter + 1
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    DEALLOCATE(u2_2D)
    DEALLOCATE(u1_2D)

    IF(PRESENT(q1)) THEN
       DEALLOCATE(q2_3D)
       DEALLOCATE(q1_3D)
    ENDIF

  ENDSUBROUTINE projectSolutionDifferentMeshes_general

  SUBROUTINE projectSolutionDifferentMeshes_general_arrays(T1, X1, T2, X2, u1, q1, u2, q2)

    INTEGER, INTENT(IN)                                :: T1(:,:), T2(:,:)
    REAL*8, INTENT(IN)                                 :: X1(:,:), X2(:,:)
    REAL*8, DIMENSION(:), INTENT(IN)                   :: u1(:)
    REAL*8, DIMENSION(:), OPTIONAL,INTENT(IN)          :: q1(:)
    REAL*8, DIMENSION(:), INTENT(INOUT)                :: u2
    REAL*8, DIMENSION(:), OPTIONAL,INTENT(INOUT)       :: q2
    REAL*8, ALLOCATABLE, DIMENSION(:,:)                :: u1_2D, u2_2D
    REAL*8, ALLOCATABLE, DIMENSION(:,:,:)              :: q1_3D, q2_3D
    INTEGER                                            :: i,j,k, counter

    ALLOCATE(u1_2D(SIZE(T1,1)*SIZE(T1,2), phys%neq))
    ALLOCATE(u2_2D(SIZE(T2,1)*SIZE(T2,2), phys%neq))


    IF(PRESENT(q1)) THEN
       ALLOCATE(q1_3D(SIZE(T1,1)*SIZE(T1,2), phys%neq, refElPol%ndim))
       ALLOCATE(q2_3D(SIZE(T2,1)*SIZE(T2,2), phys%neq, refElPol%ndim))
       q1_3D = 0.
       q2_3D = 0.
    ENDIF

    u1_2D = 0.
    u2_2D = 0.


    counter = 1
    ! reshape u sol
    DO i = 1, SIZE(u1_2D,1)
       DO j = 1, SIZE(u1_2D,2)
          u1_2D(i,j) = u1(counter)
          counter = counter + 1
       ENDDO
    ENDDO


    IF(PRESENT(q1)) THEN
       counter = 1
       DO i = 1, SIZE(q1_3D,1)
          DO j = 1, SIZE(q1_3D,2)
             DO k = 1, SIZE(q1_3D,3)
                q1_3D(i,j,k) = q1(counter)
                counter = counter + 1
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF(PRESENT(q1)) THEN
       CALL projectSolutionDifferentMeshes_mod(T1, X1, T2, X2, u1_2D, q1_3D, u2_2D, q2_3D)
    ELSE
       CALL projectSolutionDifferentMeshes_mod(T1, X1, T2, X2, u_old=u1_2D, u_new = u2_2D)
    ENDIF

    ! solu is easy to reshape
    counter = 1
    DO i = 1, SIZE(u2_2D,1)
       DO j = 1, SIZE(u2_2D,2)
          u2(counter)   = u2_2D(i,j)
          counter = counter + 1
       ENDDO
    ENDDO

    IF(PRESENT(q1)) THEN
       ! solq is not as easy
       counter = 1
       DO i = 1, SIZE(q2_3D,1)
          DO j = 1, SIZE(q2_3D,2)
             DO k = 1, SIZE(q2_3D,3)
                q2(counter)   = q2_3D(i,j,k)
                counter = counter + 1
             ENDDO
          ENDDO
       ENDDO
       DEALLOCATE(q2_3D)
       DEALLOCATE(q1_3D)
    ENDIF

    DEALLOCATE(u2_2D)
    DEALLOCATE(u1_2D)


  ENDSUBROUTINE projectSolutionDifferentMeshes_general_arrays

  PURE SUBROUTINE curved_element_bounding_box(element_coordinates, xmin, xmax, ymin, ymax, element_size)
    REAL*8, INTENT(IN)          :: element_coordinates(:,:)
    REAL*8, INTENT(OUT)         :: xmin, xmax, ymin, ymax, element_size

    xmin = MINVAL(element_coordinates(:,1))
    xmax = MAXVAL(element_coordinates(:,1))
    ymin = MINVAL(element_coordinates(:,2))
    ymax = MAXVAL(element_coordinates(:,2))
    element_size = MAX(xmax-xmin, ymax-ymin)
  ENDSUBROUTINE curved_element_bounding_box

  PURE LOGICAL FUNCTION point_in_padded_bounding_box(point, xmin, xmax, ymin, ymax, padding)
    REAL*8, INTENT(IN)          :: point(1,2)
    REAL*8, INTENT(IN)          :: xmin, xmax, ymin, ymax, padding

    point_in_padded_bounding_box = point(1,1) .GE. xmin-padding .AND. point(1,1) .LE. xmax+padding .AND. &
         point(1,2) .GE. ymin-padding .AND. point(1,2) .LE. ymax+padding
  ENDFUNCTION point_in_padded_bounding_box

  SUBROUTINE find_nearest_curved_element(target_point, old_connectivity, old_coordinates, candidate_box_padding_ratio, &
       nearest_element, nearest_valid_point, nearest_distance, nearest_element_size)
    USE adaptivity_common_module, ONLY: inverse_isop_transf, clamp_to_curved_triangle

    REAL*8, INTENT(IN)          :: target_point(1,2)
    INTEGER, INTENT(IN)         :: old_connectivity(:,:)
    REAL*8, INTENT(IN)          :: old_coordinates(:,:)
    REAL*8, INTENT(IN)          :: candidate_box_padding_ratio
    INTEGER, INTENT(OUT)        :: nearest_element
    REAL*8, INTENT(OUT)         :: nearest_valid_point(1,2), nearest_distance, nearest_element_size

    REAL*8                      :: element_coordinates(Mesh%Nnodesperelem, refElPol%Ndim)
    REAL*8                      :: reference_point(1,2), clamped_point(1,2)
    REAL*8                      :: xmin, xmax, ymin, ymax, element_size, bbox_pad, distance
    INTEGER                     :: element
    LOGICAL                     :: inverse_converged, clamp_valid

    nearest_element = 0
    nearest_valid_point = target_point
    nearest_distance = HUGE(1.d0)
    nearest_element_size = 0.d0

    DO element = 1, SIZE(old_connectivity,1)
       element_coordinates = old_coordinates(old_connectivity(element,:),:)
       CALL curved_element_bounding_box(element_coordinates, xmin, xmax, ymin, ymax, element_size)
       bbox_pad = MAX(1.d-10, candidate_box_padding_ratio*element_size)

       IF(.NOT. point_in_padded_bounding_box(target_point, xmin, xmax, ymin, ymax, bbox_pad)) CYCLE

       CALL inverse_isop_transf(target_point, element_coordinates, refElPol, reference_point, inverse_converged)
       IF(.NOT. inverse_converged) CYCLE

       CALL clamp_to_curved_triangle(reference_point, element_coordinates, refElPol, clamped_point, clamp_valid)
       IF(.NOT. clamp_valid) CYCLE

       distance = SQRT((target_point(1,1)-clamped_point(1,1))**2 + (target_point(1,2)-clamped_point(1,2))**2)
       IF(distance .LT. nearest_distance) THEN
          nearest_element = element
          nearest_valid_point = clamped_point
          nearest_distance = distance
          nearest_element_size = element_size
       ENDIF
    ENDDO
  ENDSUBROUTINE find_nearest_curved_element

  SUBROUTINE report_unmatched_projection_points(target_points, point_elements, nearest_candidate_elements, &
       nearest_candidate_distances, nearest_candidate_sizes)
    REAL*8, INTENT(IN)          :: target_points(:,:)
    INTEGER, INTENT(IN)         :: point_elements(:), nearest_candidate_elements(:)
    REAL*8, INTENT(IN)          :: nearest_candidate_distances(:), nearest_candidate_sizes(:)

    INTEGER                     :: point, local_missing, global_missing, ierr

    local_missing = COUNT(point_elements .EQ. 0)
    global_missing = local_missing
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, global_missing, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF(MPIvar%glob_id .EQ. 0 .AND. global_missing .NE. 0) THEN
       WRITE(*,*) "Projection unmatched points after fallbacks: ", global_missing
    ENDIF

    IF(local_missing .EQ. 0) RETURN

    WRITE(*,*) "Projection unmatched points on rank: ", MPIvar%glob_id, local_missing
    DO point = 1, SIZE(target_points,1)
       IF(point_elements(point) .NE. 0) CYCLE

       IF(nearest_candidate_elements(point) .NE. 0) THEN
          WRITE(*,*) "Unmatched projection point rank/index/xy/nearest/dist/local h/ratio: ", &
               MPIvar%glob_id, point, target_points(point,1), target_points(point,2), nearest_candidate_elements(point), &
               nearest_candidate_distances(point), nearest_candidate_sizes(point), &
               nearest_candidate_distances(point)/MAX(1.d-12, nearest_candidate_sizes(point))
       ELSE
          WRITE(*,*) "Unmatched projection point rank/index/xy/no nearest candidate: ", &
               MPIvar%glob_id, point, target_points(point,1), target_points(point,2)
       ENDIF
    ENDDO
    WRITE(*,*) "WARNING: projection points remain unmatched; their projected values remain zero."
  ENDSUBROUTINE report_unmatched_projection_points

  SUBROUTINE recover_nearest_projection_points(target_points, old_connectivity, old_coordinates, point_elements, &
       interpolation_points)
    REAL*8, INTENT(IN)          :: target_points(:,:), old_coordinates(:,:)
    INTEGER, INTENT(IN)         :: old_connectivity(:,:)
    INTEGER, INTENT(INOUT)      :: point_elements(:)
    REAL*8, INTENT(INOUT)       :: interpolation_points(:,:)

    INTEGER                     :: nearest_candidate_elements(SIZE(target_points,1))
    REAL*8                      :: nearest_candidate_distances(SIZE(target_points,1))
    REAL*8                      :: nearest_candidate_sizes(SIZE(target_points,1))
    REAL*8                      :: target_point(1,2), nearest_valid_point(1,2)
    REAL*8                      :: nearest_tolerance, nearest_distance, nearest_element_size
    INTEGER                     :: point, nearest_element

    nearest_candidate_elements = 0
    nearest_candidate_distances = HUGE(1.d0)
    nearest_candidate_sizes = 0.d0
    nearest_tolerance = 2.d-4

    DO point = 1, SIZE(target_points,1)
       IF(point_elements(point) .NE. 0) CYCLE
       target_point(1,:) = target_points(point,:)
       CALL find_nearest_curved_element(target_point, old_connectivity, old_coordinates, 0.5d0, nearest_element, &
            nearest_valid_point, nearest_distance, nearest_element_size)
       nearest_candidate_elements(point) = nearest_element
       nearest_candidate_distances(point) = nearest_distance
       nearest_candidate_sizes(point) = nearest_element_size

       IF(nearest_element .NE. 0 .AND. nearest_distance .LE. nearest_tolerance*MAX(1.d-12, nearest_element_size)) THEN
          point_elements(point) = nearest_element
          interpolation_points(point,:) = nearest_valid_point(1,:)
       ENDIF
    ENDDO

    CALL report_unmatched_projection_points(target_points, point_elements, nearest_candidate_elements, &
         nearest_candidate_distances, nearest_candidate_sizes)
  ENDSUBROUTINE recover_nearest_projection_points

  PURE SUBROUTINE invert_2x2_matrix(matrix, inverse_matrix)
    REAL*8, INTENT(IN)          :: matrix(2,2)
    REAL*8, INTENT(OUT)         :: inverse_matrix(2,2)
    REAL*8                      :: determinant

    determinant = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)

    inverse_matrix(1,1) =  matrix(2,2)/determinant
    inverse_matrix(1,2) = -matrix(1,2)/determinant
    inverse_matrix(2,1) = -matrix(2,1)/determinant
    inverse_matrix(2,2) =  matrix(1,1)/determinant
  ENDSUBROUTINE invert_2x2_matrix

  SUBROUTINE find_points_in_linear_elements(target_points, old_connectivity, old_coordinates, point_elements)
    REAL*8, INTENT(IN)          :: target_points(:,:), old_coordinates(:,:)
    INTEGER, INTENT(IN)         :: old_connectivity(:,:)
    INTEGER, INTENT(INOUT)      :: point_elements(:)

    REAL*8                      :: triangle_vertices(3, refElPol%Ndim)
    REAL*8                      :: edge_matrix(2,2), inverse_edge_matrix(2,2)
    REAL*8                      :: barycentric_weights(3), point_offset(2)
    REAL*8                      :: barycentric_tolerance
    INTEGER                     :: point, element

    barycentric_tolerance = 1.d-10

    DO element = 1, SIZE(old_connectivity,1)
       triangle_vertices = old_coordinates(old_connectivity(element,1:3),:)
       edge_matrix(:,1) = triangle_vertices(2,:) - triangle_vertices(1,:)
       edge_matrix(:,2) = triangle_vertices(3,:) - triangle_vertices(1,:)
       CALL invert_2x2_matrix(edge_matrix, inverse_edge_matrix)

       DO point = 1, SIZE(target_points,1)
          IF(point_elements(point) .NE. 0) CYCLE
          point_offset = target_points(point,:) - triangle_vertices(1,:)
          barycentric_weights(2:3) = MATMUL(inverse_edge_matrix, point_offset)
          barycentric_weights(1) = 1.d0 - SUM(barycentric_weights(2:3))

          IF(barycentric_weights(1) .GE. -barycentric_tolerance .AND. &
               barycentric_weights(2) .GE. -barycentric_tolerance .AND. &
               barycentric_weights(3) .GE. -barycentric_tolerance .AND. &
               barycentric_weights(1) .LE. 1.d0+barycentric_tolerance .AND. &
               barycentric_weights(2) .LE. 1.d0+barycentric_tolerance .AND. &
               barycentric_weights(3) .LE. 1.d0+barycentric_tolerance) THEN
             point_elements(point) = element
          ENDIF
       ENDDO
    ENDDO
  ENDSUBROUTINE find_points_in_linear_elements

  SUBROUTINE find_points_in_curved_elements(target_points, old_connectivity, old_coordinates, point_elements)
    USE adaptivity_common_module, ONLY: inverse_isop_transf

    REAL*8, INTENT(IN)          :: target_points(:,:), old_coordinates(:,:)
    INTEGER, INTENT(IN)         :: old_connectivity(:,:)
    INTEGER, INTENT(INOUT)      :: point_elements(:)

    REAL*8                      :: element_coordinates(Mesh%Nnodesperelem, refElPol%Ndim)
    REAL*8                      :: target_point(1,2), reference_point(1,2)
    REAL*8                      :: curved_tolerance, bounding_box_padding
    REAL*8                      :: xmin, xmax, ymin, ymax, element_size
    INTEGER                     :: point, element
    LOGICAL                     :: inverse_converged

    curved_tolerance = 1.d-8

    DO point = 1, SIZE(target_points,1)
       IF(point_elements(point) .NE. 0) CYCLE
       target_point(1,:) = target_points(point,:)

       DO element = 1, SIZE(old_connectivity,1)
          element_coordinates = old_coordinates(old_connectivity(element,:),:)
          CALL curved_element_bounding_box(element_coordinates, xmin, xmax, ymin, ymax, element_size)
          bounding_box_padding = MAX(1.d-12, curved_tolerance*MAX(1.d0, element_size))

          IF(.NOT. point_in_padded_bounding_box(target_point, xmin, xmax, ymin, ymax, bounding_box_padding)) CYCLE

          CALL inverse_isop_transf(target_point, element_coordinates, refElPol, reference_point, inverse_converged)
          IF(.NOT. inverse_converged) CYCLE

          IF(reference_point(1,1) .GE. -1.d0-curved_tolerance .AND. &
               reference_point(1,2) .GE. -1.d0-curved_tolerance .AND. &
               reference_point(1,1) .LE.  1.d0+curved_tolerance .AND. &
               reference_point(1,2) .LE.  1.d0+curved_tolerance .AND. &
               reference_point(1,1)+reference_point(1,2) .LE. curved_tolerance) THEN
             point_elements(point) = element
             EXIT
          ENDIF
       ENDDO
    ENDDO
  ENDSUBROUTINE find_points_in_curved_elements

  SUBROUTINE interpolate_solution_at_projection_points(interpolation_points, point_elements, old_connectivity, &
       old_coordinates, u_old, q_old, u_new, q_new)
    USE adaptivity_common_module, ONLY: find_matches_int, inverse_isop_transf
    USE reference_element, ONLY: compute_shape_functions_at_points

    REAL*8, INTENT(IN)          :: interpolation_points(:,:), old_coordinates(:,:)
    INTEGER, INTENT(IN)         :: point_elements(:), old_connectivity(:,:)
    REAL*8, INTENT(IN)          :: u_old(:,:)
    REAL*8, OPTIONAL,INTENT(IN) :: q_old(:,:,:)
    REAL*8, INTENT(OUT)         :: u_new(:,:)
    REAL*8, OPTIONAL,INTENT(OUT):: q_new(:,:,:)

    REAL*8                      :: element_coordinates(Mesh%Nnodesperelem, refElPol%Ndim)
    INTEGER                     :: old_element_dofs(SIZE(old_connectivity,2))
    INTEGER                     :: point, element, nodes_per_element
    REAL*8, ALLOCATABLE         :: shape_functions(:,:,:)
    REAL*8, ALLOCATABLE         :: element_points(:,:), reference_points(:,:)
    REAL*8, ALLOCATABLE         :: old_element_u(:,:), old_element_q(:,:,:)
    INTEGER, ALLOCATABLE        :: point_indices(:)

    u_new = 0.d0
    IF(PRESENT(q_new)) q_new = 0.d0

    nodes_per_element = SIZE(old_connectivity,2)
    ALLOCATE(old_element_u(nodes_per_element, SIZE(u_old,2)))
    IF(PRESENT(q_old)) ALLOCATE(old_element_q(nodes_per_element, SIZE(q_old,2),2))

    DO element = 1, SIZE(old_connectivity,1)
       CALL find_matches_int(point_elements, element, point_indices)

       ALLOCATE(reference_points(SIZE(point_indices), SIZE(interpolation_points,2)))
       ALLOCATE(element_points(SIZE(point_indices), SIZE(interpolation_points,2)))
       ALLOCATE(shape_functions(nodes_per_element, SIZE(point_indices), 3))

       reference_points = 0.d0
       element_points = interpolation_points(point_indices,:)
       shape_functions = 0.d0
       element_coordinates = old_coordinates(old_connectivity(element,:),:)

       CALL inverse_isop_transf(element_points, element_coordinates, refElPol, reference_points)

       CALL compute_shape_functions_at_points(refElPol, reference_points, shape_functions)

       old_element_dofs = (element-1)*nodes_per_element + (/ (point, point=1, nodes_per_element) /)
       old_element_u = u_old(old_element_dofs, :)
       u_new(point_indices,:) = MATMUL(TRANSPOSE(shape_functions(:,:,1)), old_element_u)

       IF(PRESENT(q_old)) THEN
          old_element_q = q_old(old_element_dofs, :, :)
          q_new(point_indices,:,1) = MATMUL(TRANSPOSE(shape_functions(:,:,1)), old_element_q(:,:,1))
          q_new(point_indices,:,2) = MATMUL(TRANSPOSE(shape_functions(:,:,1)), old_element_q(:,:,2))
       ENDIF

       DEALLOCATE(reference_points)
       DEALLOCATE(element_points)
       DEALLOCATE(point_indices)
       DEALLOCATE(shape_functions)
    ENDDO

    DEALLOCATE(old_element_u)
    IF(PRESENT(q_old)) DEALLOCATE(old_element_q)
  ENDSUBROUTINE interpolate_solution_at_projection_points

  SUBROUTINE projectSolutionDifferentMeshes_Mod(old_connectivity, old_coordinates, new_connectivity, new_coordinates, &
       u_old, q_old, u_new, q_new)
    USE linearAlgebra, ONLY: colint

    INTEGER, INTENT(IN)         :: old_connectivity(:,:), new_connectivity(:,:)
    REAL*8, INTENT(IN)          :: old_coordinates(:,:), new_coordinates(:,:)
    REAL*8                      :: target_points(SIZE(new_connectivity,1)*SIZE(new_connectivity,2), 2)
    REAL*8                      :: interpolation_points(SIZE(new_connectivity,1)*SIZE(new_connectivity,2), 2)

    REAL*8, INTENT(IN)          :: u_old(:,:)
    REAL*8, OPTIONAL,INTENT(IN) :: q_old(:,:,:)

    REAL*8, INTENT(OUT)         :: u_new(:,:)
    REAL*8, OPTIONAL,INTENT(OUT):: q_new(:,:,:)

    INTEGER                     :: point_elements(SIZE(target_points,1))
    target_points = new_coordinates(colint(TRANSPOSE(new_connectivity)),:)
    interpolation_points = target_points


    point_elements = 0

    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*      PROJECTING SOLUTION TO NEW MESH          *'
          WRITE (6, *) '*************************************************'
       END IF
    ENDIF

    CALL find_points_in_linear_elements(target_points, old_connectivity, old_coordinates, point_elements)

    IF(ANY(point_elements .EQ. 0)) THEN
       CALL find_points_in_curved_elements(target_points, old_connectivity, old_coordinates, point_elements)
    ENDIF

    CALL recover_nearest_projection_points(target_points, old_connectivity, old_coordinates, point_elements, &
         interpolation_points)

    CALL interpolate_solution_at_projection_points(interpolation_points, point_elements, old_connectivity, &
         old_coordinates, u_old, q_old, u_new, q_new)

  END SUBROUTINE projectSolutionDifferentMeshes_Mod

END MODULE adaptivity_projection_module
