MODULE element_mapping_module
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: Reference_element_type
  USE LinearAlgebra, ONLY: invert_matrix
  USE reference_element, ONLY: compute_shape_functions_at_points, orthopoly2d_deriv, vandermonde_2d

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: inverse_max_iterations = 50
  INTEGER, PARAMETER :: closest_max_iterations = 50
  INTEGER, PARAMETER :: max_backtracking_steps = 10
  REAL*8, PARAMETER  :: residual_tolerance = 1.d-10
  REAL*8, PARAMETER  :: reference_step_tolerance = 1.d-12

  PUBLIC :: map_physical_to_reference
  PUBLIC :: map_physical_to_nearest_reference
  PUBLIC :: map_reference_to_physical

CONTAINS

  SUBROUTINE map_physical_to_reference(physical_points, element_coordinates, refEl, reference_points, converged)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: physical_points(:,:), element_coordinates(:,:)
    REAL*8, INTENT(OUT)                     :: reference_points(:,:)
    LOGICAL, INTENT(OUT)                    :: converged(:)

    REAL*8                                  :: inverse_vandermonde(refEl%Nnodes2D,refEl%Nnodes2D)
    REAL*8                                  :: tolerance
    LOGICAL                                 :: linear_map_valid
    INTEGER                                 :: point

    CALL validate_mapping_shapes(physical_points, element_coordinates, refEl, reference_points, converged)

    converged = .FALSE.
    CALL map_physical_to_linear_reference(physical_points, element_coordinates, reference_points, linear_map_valid)
    IF(.NOT. linear_map_valid) RETURN

    CALL prepare_mapping_operator(refEl, inverse_vandermonde)
    tolerance = residual_tolerance*MAX(1.d0, mapping_element_scale(element_coordinates))

    DO point = 1, SIZE(physical_points,1)
       CALL solve_inverse_mapping_point(physical_points(point,:), element_coordinates, refEl, inverse_vandermonde, &
            tolerance, reference_points(point,:), converged(point))
    ENDDO
  ENDSUBROUTINE map_physical_to_reference

  SUBROUTINE map_physical_to_nearest_reference(physical_point, element_coordinates, refEl, reference_point, &
       mapped_point, distance, valid)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: physical_point(2), element_coordinates(:,:)
    REAL*8, INTENT(OUT)                     :: reference_point(2), mapped_point(2), distance
    LOGICAL, INTENT(OUT)                    :: valid

    REAL*8                                  :: inverse_vandermonde(refEl%Nnodes2D,refEl%Nnodes2D)
    REAL*8                                  :: seeds(8,2), candidate_reference(2), candidate_mapped(2)
    REAL*8                                  :: candidate_distance
    LOGICAL                                 :: candidate_valid
    INTEGER                                 :: seed

    CALL validate_element_coordinates(element_coordinates, refEl)
    CALL prepare_mapping_operator(refEl, inverse_vandermonde)
    CALL build_closest_point_seeds(physical_point, element_coordinates, seeds)

    valid = .FALSE.
    reference_point = 0.d0
    mapped_point = 0.d0
    distance = HUGE(1.d0)

    DO seed = 1, SIZE(seeds,1)
       candidate_reference = seeds(seed,:)
       CALL solve_closest_mapping_point(physical_point, element_coordinates, refEl, inverse_vandermonde, &
            candidate_reference, candidate_mapped, candidate_distance, candidate_valid)
       CALL retain_closer_candidate(candidate_reference, candidate_mapped, candidate_distance, candidate_valid, &
            reference_point, mapped_point, distance, valid)
    ENDDO
  ENDSUBROUTINE map_physical_to_nearest_reference

  SUBROUTINE solve_inverse_mapping_point(physical_point, element_coordinates, refEl, inverse_vandermonde, &
       tolerance, reference_point, converged)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: physical_point(2), element_coordinates(:,:), inverse_vandermonde(:,:)
    REAL*8, INTENT(IN)                      :: tolerance
    REAL*8, INTENT(INOUT)                   :: reference_point(2)
    LOGICAL, INTENT(OUT)                    :: converged

    REAL*8                                  :: mapped_point(2), trial_reference(2), trial_mapped(2)
    REAL*8                                  :: residual(2), update(2), objective, trial_objective, step_size
    LOGICAL                                 :: update_valid, step_accepted
    INTEGER                                 :: iteration

    converged = .FALSE.
    CALL regularize_reference_triangle_apex(reference_point)
    CALL map_single_reference_to_physical(reference_point, element_coordinates, refEl, mapped_point)
    IF(.NOT. ALL(ieee_is_finite(mapped_point))) RETURN

    DO iteration = 1, inverse_max_iterations
       residual = physical_point-mapped_point
       objective = SUM(residual**2)
       IF(SQRT(objective) .LE. tolerance) THEN
          converged = .TRUE.
          RETURN
       ENDIF

       CALL compute_mapping_update(reference_point, residual, element_coordinates, refEl, inverse_vandermonde, &
            update, update_valid)
       IF(.NOT. update_valid) RETURN

       CALL find_decreasing_mapping_step(physical_point, reference_point, update, objective, element_coordinates, &
            refEl, .FALSE., trial_reference, trial_mapped, trial_objective, step_size, step_accepted)
       IF(.NOT. step_accepted) RETURN

       reference_point = trial_reference
       mapped_point = trial_mapped
       IF(step_size .LE. reference_step_tolerance) THEN
          converged = SQRT(trial_objective) .LE. tolerance
          RETURN
       ENDIF
    ENDDO

    converged = SQRT(SUM((physical_point-mapped_point)**2)) .LE. tolerance
  ENDSUBROUTINE solve_inverse_mapping_point

  SUBROUTINE solve_closest_mapping_point(physical_point, element_coordinates, refEl, inverse_vandermonde, &
       reference_point, mapped_point, distance, valid)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: physical_point(2), element_coordinates(:,:), inverse_vandermonde(:,:)
    REAL*8, INTENT(INOUT)                   :: reference_point(2)
    REAL*8, INTENT(OUT)                     :: mapped_point(2), distance
    LOGICAL, INTENT(OUT)                    :: valid

    REAL*8                                  :: projected_reference(2), trial_reference(2), trial_mapped(2)
    REAL*8                                  :: residual(2), update(2), objective, trial_objective, step_size
    REAL*8                                  :: improvement_tolerance
    LOGICAL                                 :: update_valid, step_accepted
    INTEGER                                 :: iteration

    CALL project_to_reference_triangle(reference_point, projected_reference)
    reference_point = projected_reference
    CALL regularize_reference_triangle_apex(reference_point)
    CALL map_single_reference_to_physical(reference_point, element_coordinates, refEl, mapped_point)
    valid = ALL(ieee_is_finite(mapped_point))
    IF(.NOT. valid) THEN
       distance = HUGE(1.d0)
       RETURN
    ENDIF

    objective = SUM((physical_point-mapped_point)**2)
    DO iteration = 1, closest_max_iterations
       residual = physical_point-mapped_point
       CALL compute_mapping_update(reference_point, residual, element_coordinates, refEl, inverse_vandermonde, &
            update, update_valid)
       IF(.NOT. update_valid) EXIT

       CALL find_decreasing_mapping_step(physical_point, reference_point, update, objective, element_coordinates, &
            refEl, .TRUE., trial_reference, trial_mapped, trial_objective, step_size, step_accepted)
       IF(.NOT. step_accepted) EXIT

       improvement_tolerance = 100.d0*EPSILON(1.d0)*MAX(1.d0,objective)
       reference_point = trial_reference
       mapped_point = trial_mapped
       IF(step_size .LE. reference_step_tolerance .OR. objective-trial_objective .LE. improvement_tolerance) THEN
          objective = trial_objective
          EXIT
       ENDIF
       objective = trial_objective
    ENDDO

    distance = SQRT(MAX(0.d0,objective))
    valid = ALL(ieee_is_finite(reference_point)) .AND. ALL(ieee_is_finite(mapped_point)) .AND. &
         ieee_is_finite(distance) .AND. reference_point_is_in_triangle(reference_point, 10.d0*EPSILON(1.d0))
  ENDSUBROUTINE solve_closest_mapping_point

  SUBROUTINE find_decreasing_mapping_step(physical_point, reference_point, update, objective, element_coordinates, &
       refEl, constrain_to_triangle, trial_reference, trial_mapped, trial_objective, step_size, accepted)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: physical_point(2), reference_point(2), update(2), objective
    REAL*8, INTENT(IN)                      :: element_coordinates(:,:)
    LOGICAL, INTENT(IN)                     :: constrain_to_triangle
    REAL*8, INTENT(OUT)                     :: trial_reference(2), trial_mapped(2), trial_objective, step_size
    LOGICAL, INTENT(OUT)                    :: accepted

    REAL*8                                  :: unconstrained_trial(2), damping
    INTEGER                                 :: backtrack

    accepted = .FALSE.
    damping = 1.d0
    DO backtrack = 1, max_backtracking_steps
       unconstrained_trial = reference_point+damping*update
       IF(constrain_to_triangle) THEN
          CALL project_to_reference_triangle(unconstrained_trial, trial_reference)
       ELSE
          trial_reference = unconstrained_trial
       ENDIF
       CALL regularize_reference_triangle_apex(trial_reference)
       CALL map_single_reference_to_physical(trial_reference, element_coordinates, refEl, trial_mapped)
       trial_objective = SUM((physical_point-trial_mapped)**2)

       IF(ALL(ieee_is_finite(trial_mapped)) .AND. trial_objective .LT. objective) THEN
          accepted = .TRUE.
          step_size = SQRT(SUM((trial_reference-reference_point)**2))
          RETURN
       ENDIF
       damping = 0.5d0*damping
    ENDDO

    trial_reference = reference_point
    CALL map_single_reference_to_physical(reference_point, element_coordinates, refEl, trial_mapped)
    trial_objective = objective
    step_size = 0.d0
  ENDSUBROUTINE find_decreasing_mapping_step

  SUBROUTINE compute_mapping_update(reference_point, residual, element_coordinates, refEl, inverse_vandermonde, &
       update, valid)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: reference_point(2), residual(2), element_coordinates(:,:)
    REAL*8, INTENT(IN)                      :: inverse_vandermonde(:,:)
    REAL*8, INTENT(OUT)                     :: update(2)
    LOGICAL, INTENT(OUT)                    :: valid

    REAL*8                                  :: jacobian(2,2), inverse_jacobian(2,2)

    CALL compute_mapping_jacobian(reference_point, element_coordinates, refEl, inverse_vandermonde, jacobian)
    CALL invert_2x2_matrix(jacobian, inverse_jacobian, valid)
    IF(valid) THEN
       update = MATMUL(inverse_jacobian,residual)
       valid = ALL(ieee_is_finite(update))
    ELSE
       update = 0.d0
    ENDIF
  ENDSUBROUTINE compute_mapping_update

  SUBROUTINE compute_mapping_jacobian(reference_point, element_coordinates, refEl, inverse_vandermonde, jacobian)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: reference_point(2), element_coordinates(:,:), inverse_vandermonde(:,:)
    REAL*8, INTENT(OUT)                     :: jacobian(2,2)

    REAL*8                                  :: polynomial(refEl%Nnodes2D)
    REAL*8                                  :: derivative_xi(refEl%Nnodes2D), derivative_eta(refEl%Nnodes2D)
    REAL*8                                  :: shape_derivative_xi(refEl%Nnodes2D), shape_derivative_eta(refEl%Nnodes2D)

    polynomial = 0.d0
    derivative_xi = 0.d0
    derivative_eta = 0.d0
    CALL orthopoly2d_deriv(reference_point(1), reference_point(2), refEl%Ndeg, refEl%Nnodes2D, polynomial, &
         derivative_xi, derivative_eta)

    shape_derivative_xi = MATMUL(inverse_vandermonde, derivative_xi)
    shape_derivative_eta = MATMUL(inverse_vandermonde, derivative_eta)
    jacobian(:,1) = MATMUL(TRANSPOSE(element_coordinates),shape_derivative_xi)
    jacobian(:,2) = MATMUL(TRANSPOSE(element_coordinates),shape_derivative_eta)
  ENDSUBROUTINE compute_mapping_jacobian

  SUBROUTINE invert_2x2_matrix(matrix, inverse_matrix, valid)
    REAL*8, INTENT(IN)  :: matrix(2,2)
    REAL*8, INTENT(OUT) :: inverse_matrix(2,2)
    LOGICAL, INTENT(OUT):: valid

    REAL*8              :: determinant, matrix_scale, determinant_tolerance

    determinant = matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
    matrix_scale = MAXVAL(ABS(matrix))
    determinant_tolerance = 100.d0*EPSILON(1.d0)*matrix_scale**2
    valid = ALL(ieee_is_finite(matrix)) .AND. ieee_is_finite(determinant) .AND. &
         matrix_scale .GT. 0.d0 .AND. ABS(determinant) .GT. determinant_tolerance
    IF(.NOT. valid) THEN
       inverse_matrix = 0.d0
       RETURN
    ENDIF

    inverse_matrix = RESHAPE((/matrix(2,2), -matrix(2,1), -matrix(1,2), matrix(1,1)/),(/2,2/))/determinant
  ENDSUBROUTINE invert_2x2_matrix

  SUBROUTINE prepare_mapping_operator(refEl, inverse_vandermonde)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(OUT)                     :: inverse_vandermonde(:,:)
    REAL*8                                  :: vandermonde(refEl%Nnodes2D,refEl%Nnodes2D)

    CALL vandermonde_2d(vandermonde, refEl)
    CALL invert_matrix(TRANSPOSE(vandermonde), inverse_vandermonde)
  ENDSUBROUTINE prepare_mapping_operator

  SUBROUTINE map_physical_to_linear_reference(physical_points, element_coordinates, reference_points, valid)
    REAL*8, INTENT(IN)  :: physical_points(:,:), element_coordinates(:,:)
    REAL*8, INTENT(OUT) :: reference_points(:,:)
    LOGICAL, INTENT(OUT):: valid

    REAL*8              :: jacobian(2,2), inverse_jacobian(2,2), center(2)
    REAL*8              :: offsets(SIZE(physical_points,1),2)

    jacobian(:,1) = 0.5d0*(element_coordinates(2,:)-element_coordinates(1,:))
    jacobian(:,2) = 0.5d0*(element_coordinates(3,:)-element_coordinates(1,:))
    CALL invert_2x2_matrix(jacobian, inverse_jacobian, valid)
    IF(.NOT. valid) THEN
       reference_points = 0.d0
       RETURN
    ENDIF

    center = 0.5d0*(element_coordinates(2,:)+element_coordinates(3,:))
    offsets = physical_points-SPREAD(center,DIM=1,NCOPIES=SIZE(physical_points,1))
    reference_points = MATMUL(offsets,TRANSPOSE(inverse_jacobian))
  ENDSUBROUTINE map_physical_to_linear_reference

  SUBROUTINE build_closest_point_seeds(physical_point, element_coordinates, seeds)
    REAL*8, INTENT(IN)  :: physical_point(2), element_coordinates(:,:)
    REAL*8, INTENT(OUT) :: seeds(8,2)

    REAL*8              :: physical_points(1,2), initial_reference(1,2)
    LOGICAL             :: linear_map_valid

    physical_points(1,:) = physical_point
    CALL map_physical_to_linear_reference(physical_points, element_coordinates, initial_reference, linear_map_valid)
    IF(linear_map_valid) THEN
       CALL project_to_reference_triangle(initial_reference(1,:), seeds(1,:))
    ELSE
       seeds(1,:) = (/-1.d0/3.d0, -1.d0/3.d0/)
    ENDIF

    seeds(2,:) = (/-1.d0/3.d0, -1.d0/3.d0/)
    seeds(3,:) = (/-1.d0, -1.d0/)
    seeds(4,:) = (/ 1.d0, -1.d0/)
    seeds(5,:) = (/-1.d0,  1.d0/)
    seeds(6,:) = (/ 0.d0, -1.d0/)
    seeds(7,:) = (/ 0.d0,  0.d0/)
    seeds(8,:) = (/-1.d0,  0.d0/)
  ENDSUBROUTINE build_closest_point_seeds

  SUBROUTINE retain_closer_candidate(candidate_reference, candidate_mapped, candidate_distance, candidate_valid, &
       reference_point, mapped_point, distance, valid)
    REAL*8, INTENT(IN)    :: candidate_reference(2), candidate_mapped(2), candidate_distance
    LOGICAL, INTENT(IN)   :: candidate_valid
    REAL*8, INTENT(INOUT) :: reference_point(2), mapped_point(2), distance
    LOGICAL, INTENT(INOUT):: valid

    IF(.NOT. candidate_valid) RETURN
    IF(valid .AND. candidate_distance .GE. distance) RETURN

    valid = .TRUE.
    reference_point = candidate_reference
    mapped_point = candidate_mapped
    distance = candidate_distance
  ENDSUBROUTINE retain_closer_candidate

  SUBROUTINE map_reference_to_physical(reference_points, element_coordinates, refEl, physical_points)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: reference_points(:,:), element_coordinates(:,:)
    REAL*8, INTENT(OUT)                     :: physical_points(:,:)
    REAL*8                                  :: shape_functions(refEl%Nnodes2D,SIZE(reference_points,1),3)

    CALL validate_forward_mapping_shapes(reference_points, element_coordinates, refEl, physical_points)
    shape_functions = 0.d0
    CALL compute_shape_functions_at_points(refEl, reference_points, shape_functions)
    physical_points = MATMUL(TRANSPOSE(shape_functions(:,:,1)),element_coordinates)
  ENDSUBROUTINE map_reference_to_physical

  SUBROUTINE map_single_reference_to_physical(reference_point, element_coordinates, refEl, physical_point)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: reference_point(2), element_coordinates(:,:)
    REAL*8, INTENT(OUT)                     :: physical_point(2)
    REAL*8                                  :: reference_points(1,2), physical_points(1,2)

    reference_points(1,:) = reference_point
    CALL map_reference_to_physical(reference_points, element_coordinates, refEl, physical_points)
    physical_point = physical_points(1,:)
  ENDSUBROUTINE map_single_reference_to_physical

  PURE SUBROUTINE project_to_reference_triangle(point, projected_point)
    REAL*8, INTENT(IN)  :: point(2)
    REAL*8, INTENT(OUT) :: projected_point(2)

    REAL*8, PARAMETER   :: vertices(2,3) = RESHAPE((/-1.d0,-1.d0, 1.d0,-1.d0, -1.d0,1.d0/),(/2,3/))
    REAL*8              :: edge(2), edge_point(2), parameter, distance_squared, best_distance
    INTEGER             :: face, next_vertex

    IF(reference_point_is_in_triangle(point,0.d0)) THEN
       projected_point = point
       RETURN
    ENDIF

    best_distance = HUGE(1.d0)
    projected_point = vertices(:,1)
    DO face = 1, 3
       next_vertex = MOD(face,3)+1
       edge = vertices(:,next_vertex)-vertices(:,face)
       parameter = DOT_PRODUCT(point-vertices(:,face),edge)/SUM(edge**2)
       edge_point = vertices(:,face)+MAX(0.d0,MIN(1.d0,parameter))*edge
       distance_squared = SUM((point-edge_point)**2)
       IF(distance_squared .LT. best_distance) THEN
          best_distance = distance_squared
          projected_point = edge_point
       ENDIF
    ENDDO
  ENDSUBROUTINE project_to_reference_triangle

  PURE LOGICAL FUNCTION reference_point_is_in_triangle(point, tolerance)
    REAL*8, INTENT(IN) :: point(2), tolerance

    reference_point_is_in_triangle = point(1) .GE. -1.d0-tolerance .AND. &
         point(2) .GE. -1.d0-tolerance .AND. point(1)+point(2) .LE. tolerance
  ENDFUNCTION reference_point_is_in_triangle

  PURE SUBROUTINE regularize_reference_triangle_apex(reference_point)
    REAL*8, INTENT(INOUT) :: reference_point(2)

    IF(ABS(reference_point(2)-1.d0) .LT. 1.d-12) reference_point(2) = reference_point(2)-1.d-10
  ENDSUBROUTINE regularize_reference_triangle_apex

  PURE REAL*8 FUNCTION mapping_element_scale(element_coordinates)
    REAL*8, INTENT(IN) :: element_coordinates(:,:)

    mapping_element_scale = MAX(MAXVAL(element_coordinates(:,1))-MINVAL(element_coordinates(:,1)), &
         MAXVAL(element_coordinates(:,2))-MINVAL(element_coordinates(:,2)))
  ENDFUNCTION mapping_element_scale

  SUBROUTINE validate_mapping_shapes(physical_points, element_coordinates, refEl, reference_points, converged)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: physical_points(:,:), element_coordinates(:,:)
    REAL*8, INTENT(IN)                      :: reference_points(:,:)
    LOGICAL, INTENT(IN)                     :: converged(:)

    CALL validate_element_coordinates(element_coordinates, refEl)
    IF(SIZE(physical_points,2) .NE. 2 .OR. SIZE(reference_points,1) .NE. SIZE(physical_points,1) .OR. &
         SIZE(reference_points,2) .NE. 2 .OR. SIZE(converged) .NE. SIZE(physical_points,1)) THEN
       ERROR STOP "Invalid point shapes in map_physical_to_reference"
    ENDIF
  ENDSUBROUTINE validate_mapping_shapes

  SUBROUTINE validate_forward_mapping_shapes(reference_points, element_coordinates, refEl, physical_points)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: reference_points(:,:), element_coordinates(:,:), physical_points(:,:)

    CALL validate_element_coordinates(element_coordinates, refEl)
    IF(SIZE(reference_points,2) .NE. 2 .OR. SIZE(physical_points,1) .NE. SIZE(reference_points,1) .OR. &
         SIZE(physical_points,2) .NE. 2) ERROR STOP "Invalid point shapes in map_reference_to_physical"
  ENDSUBROUTINE validate_forward_mapping_shapes

  SUBROUTINE validate_element_coordinates(element_coordinates, refEl)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: element_coordinates(:,:)

    IF(SIZE(element_coordinates,1) .NE. refEl%Nnodes2D .OR. SIZE(element_coordinates,2) .NE. 2) THEN
       ERROR STOP "Element coordinates do not match the reference element"
    ENDIF
  ENDSUBROUTINE validate_element_coordinates

END MODULE element_mapping_module
