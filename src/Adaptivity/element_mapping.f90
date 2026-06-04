MODULE element_mapping_module
  USE globals, ONLY: Reference_element_type
  USE LinearAlgebra, ONLY: invert_matrix, solve_linear_system
  USE reference_element, ONLY: compute_shape_functions_at_points, orthopoly2d_deriv, vandermonde_2d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: map_physical_to_reference
  PUBLIC :: map_reference_to_physical

CONTAINS

  SUBROUTINE map_physical_to_reference(physical_points, element_coordinates, refEl, reference_points, converged)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: physical_points(:,:), element_coordinates(:,:)
    REAL*8, INTENT(OUT)                     :: reference_points(:,:)
    LOGICAL, OPTIONAL, INTENT(OUT)          :: converged
    REAL*8                                  :: mapped_points(SIZE(physical_points,1),SIZE(physical_points,2))
    REAL*8, ALLOCATABLE                     :: unconverged_physical(:,:), unconverged_mapped(:,:)
    REAL*8, ALLOCATABLE                     :: unconverged_reference(:,:)
    INTEGER, ALLOCATABLE                    :: unconverged_indices(:)
    LOGICAL                                 :: point_unconverged(SIZE(physical_points,1)), newton_converged
    REAL*8                                  :: initial_reference(SIZE(physical_points,1),SIZE(physical_points,2))
    REAL*8                                  :: updated_reference(SIZE(physical_points,1),SIZE(physical_points,2))
    REAL*8                                  :: tolerance

    tolerance = 1e-10
    IF(PRESENT(converged)) converged = .TRUE.

    CALL map_physical_to_linear_reference(physical_points, element_coordinates, initial_reference)
    CALL regularize_reference_triangle_apex(initial_reference)
    CALL map_reference_to_physical(initial_reference, element_coordinates, refEl, mapped_points)
    CALL find_points_outside_mapping_tolerance(physical_points, mapped_points, tolerance, point_unconverged)

    IF(ANY(point_unconverged)) THEN
       CALL find_true_indices(point_unconverged, unconverged_indices)

       ALLOCATE(unconverged_physical(SIZE(unconverged_indices),SIZE(physical_points,2)))
       ALLOCATE(unconverged_mapped(SIZE(unconverged_indices),SIZE(mapped_points,2)))
       ALLOCATE(unconverged_reference(SIZE(unconverged_indices),SIZE(reference_points,2)))
       unconverged_physical = physical_points(unconverged_indices,:)
       unconverged_mapped = mapped_points(unconverged_indices,:)
       unconverged_reference = initial_reference(unconverged_indices,:)

       CALL refine_reference_points_newton(unconverged_physical, element_coordinates, refEl, tolerance, &
            unconverged_reference, unconverged_mapped, newton_converged)

       IF(.NOT. newton_converged) THEN
          IF(PRESENT(converged)) THEN
             converged = .FALSE.
             reference_points = initial_reference
             DEALLOCATE(unconverged_indices, unconverged_physical, unconverged_mapped, unconverged_reference)
             RETURN
          ELSE
             WRITE(*,*) "map_physical_to_reference non converging."
             STOP
          ENDIF
       ENDIF

       updated_reference = initial_reference
       updated_reference(unconverged_indices,:) = unconverged_reference
       reference_points = updated_reference

       DEALLOCATE(unconverged_indices, unconverged_physical, unconverged_mapped, unconverged_reference)
    ELSE
       reference_points = initial_reference
    ENDIF
  ENDSUBROUTINE map_physical_to_reference

  PURE SUBROUTINE find_points_outside_mapping_tolerance(physical_points, mapped_points, tolerance, point_unconverged)
    REAL*8, INTENT(IN)  :: physical_points(:,:), mapped_points(:,:), tolerance
    LOGICAL, INTENT(OUT):: point_unconverged(:)

    point_unconverged = SQRT((physical_points(:,1)-mapped_points(:,1))**2 + &
         (physical_points(:,2)-mapped_points(:,2))**2) .GT. &
         (tolerance*SQRT(physical_points(:,1)**2+physical_points(:,2)**2)+1.e-14)
  ENDSUBROUTINE find_points_outside_mapping_tolerance

  PURE SUBROUTINE find_true_indices(mask, indices)
    LOGICAL, INTENT(IN)              :: mask(:)
    INTEGER, ALLOCATABLE, INTENT(OUT):: indices(:)
    INTEGER                          :: point, match

    ALLOCATE(indices(COUNT(mask)))
    match = 1
    DO point = 1, SIZE(mask)
       IF(mask(point)) THEN
          indices(match) = point
          match = match + 1
       ENDIF
    ENDDO
  ENDSUBROUTINE find_true_indices

  SUBROUTINE refine_reference_points_newton(physical_points, element_coordinates, refEl, tolerance, reference_points, &
       mapped_points, converged)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: physical_points(:,:), element_coordinates(:,:), tolerance
    REAL*8, INTENT(INOUT)                   :: reference_points(:,:), mapped_points(:,:)
    LOGICAL, INTENT(OUT)                    :: converged
    INTEGER, PARAMETER                      :: max_iterations = 5
    REAL*8                                  :: residual(SIZE(physical_points,1),SIZE(physical_points,2))
    REAL*8                                  :: vandermonde(refEl%Nnodes2D,refEl%Nnodes2D)
    REAL*8                                  :: inverse_vandermonde(refEl%Nnodes2D,refEl%Nnodes2D)
    REAL*8                                  :: jacobian_x_xi, jacobian_x_eta, jacobian_y_xi, jacobian_y_eta, determinant
    LOGICAL                                 :: point_unconverged(SIZE(physical_points,1)), singular_jacobian
    INTEGER                                 :: iteration, point

    CALL vandermonde_2d(vandermonde, refEl)
    CALL invert_matrix(TRANSPOSE(vandermonde), inverse_vandermonde)

    singular_jacobian = .FALSE.
    DO iteration = 1, max_iterations
       CALL find_points_outside_mapping_tolerance(physical_points, mapped_points, tolerance, point_unconverged)
       IF(.NOT. ANY(point_unconverged)) EXIT

       residual = physical_points-mapped_points
       DO point = 1, SIZE(reference_points,1)
          CALL compute_mapping_jacobian(reference_points(point,:), element_coordinates, refEl, inverse_vandermonde, &
               jacobian_x_xi, jacobian_x_eta, jacobian_y_xi, jacobian_y_eta, determinant)
          IF(mapping_jacobian_is_nearly_singular(jacobian_x_xi, jacobian_x_eta, jacobian_y_xi, &
               jacobian_y_eta, determinant)) THEN
             singular_jacobian = .TRUE.
             EXIT
          ENDIF

          reference_points(point,1) = reference_points(point,1) + &
               (residual(point,1)*jacobian_y_eta-residual(point,2)*jacobian_x_eta)/determinant
          reference_points(point,2) = reference_points(point,2) + &
               (residual(point,2)*jacobian_x_xi-residual(point,1)*jacobian_y_xi)/determinant
       ENDDO

       IF(singular_jacobian) EXIT

       CALL regularize_reference_triangle_apex(reference_points)
       CALL map_reference_to_physical(reference_points, element_coordinates, refEl, mapped_points)
    ENDDO

    CALL find_points_outside_mapping_tolerance(physical_points, mapped_points, tolerance, point_unconverged)
    converged = .NOT. singular_jacobian .AND. .NOT. ANY(point_unconverged)
  ENDSUBROUTINE refine_reference_points_newton

  SUBROUTINE compute_mapping_jacobian(reference_point, element_coordinates, refEl, inverse_vandermonde, &
       jacobian_x_xi, jacobian_x_eta, jacobian_y_xi, jacobian_y_eta, determinant)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: reference_point(:), element_coordinates(:,:), inverse_vandermonde(:,:)
    REAL*8, INTENT(OUT)                     :: jacobian_x_xi, jacobian_x_eta, jacobian_y_xi, jacobian_y_eta, determinant
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

    jacobian_x_xi = DOT_PRODUCT(shape_derivative_xi,element_coordinates(:,1))
    jacobian_x_eta = DOT_PRODUCT(shape_derivative_eta,element_coordinates(:,1))
    jacobian_y_xi = DOT_PRODUCT(shape_derivative_xi,element_coordinates(:,2))
    jacobian_y_eta = DOT_PRODUCT(shape_derivative_eta,element_coordinates(:,2))
    determinant = jacobian_x_xi*jacobian_y_eta-jacobian_x_eta*jacobian_y_xi
  ENDSUBROUTINE compute_mapping_jacobian

  PURE LOGICAL FUNCTION mapping_jacobian_is_nearly_singular(jacobian_x_xi, jacobian_x_eta, jacobian_y_xi, &
       jacobian_y_eta, determinant)
    REAL*8, INTENT(IN) :: jacobian_x_xi, jacobian_x_eta, jacobian_y_xi, jacobian_y_eta, determinant
    REAL*8             :: jacobian_scale, determinant_tolerance

    jacobian_scale = MAX(ABS(jacobian_x_xi), ABS(jacobian_x_eta), ABS(jacobian_y_xi), ABS(jacobian_y_eta))
    determinant_tolerance = 100.d0*EPSILON(1.d0)*jacobian_scale**2
    mapping_jacobian_is_nearly_singular = ABS(determinant) .LE. determinant_tolerance
  ENDFUNCTION mapping_jacobian_is_nearly_singular

  PURE SUBROUTINE regularize_reference_triangle_apex(reference_points)
    REAL*8, INTENT(INOUT) :: reference_points(:,:)

    WHERE(ABS(reference_points(:,2)-1.d0) .LT. 1.d-12)
       reference_points(:,2) = reference_points(:,2) - 1.d-10
    ENDWHERE
  ENDSUBROUTINE regularize_reference_triangle_apex

  SUBROUTINE map_physical_to_linear_reference(physical_points, element_coordinates, reference_points)
    REAL*8, INTENT(IN)  :: physical_points(:,:), element_coordinates(:,:)
    REAL*8, INTENT(OUT) :: reference_points(:,:)
    REAL*8              :: vertex1(2), vertex2(2), vertex3(2), jacobian(2,2)
    REAL*8              :: offsets(SIZE(physical_points,1),2)
    REAL*8              :: transposed_reference_points(SIZE(reference_points,2),SIZE(reference_points,1))

    vertex1 = element_coordinates(1,:)
    vertex2 = element_coordinates(2,:)
    vertex3 = element_coordinates(3,:)

    jacobian(:,1) = (vertex2-vertex1)/2
    jacobian(:,2) = (vertex3-vertex1)/2

    offsets(:,1) = physical_points(:,1)-(vertex2(1)+vertex3(1))/2
    offsets(:,2) = physical_points(:,2)-(vertex2(2)+vertex3(2))/2

    CALL solve_linear_system(jacobian, TRANSPOSE(offsets), transposed_reference_points)
    reference_points = TRANSPOSE(transposed_reference_points)
  ENDSUBROUTINE map_physical_to_linear_reference

  SUBROUTINE map_reference_to_physical(reference_points, element_coordinates, refEl, physical_points)
    TYPE(Reference_element_type), INTENT(IN) :: refEl
    REAL*8, INTENT(IN)                      :: reference_points(:,:), element_coordinates(:,:)
    REAL*8, INTENT(OUT)                     :: physical_points(:,:)
    REAL*8                                  :: shape_functions(refEl%Nnodes2D,SIZE(reference_points,1),3)

    shape_functions = 0.d0
    CALL compute_shape_functions_at_points(refEl, reference_points, shape_functions)

    physical_points(:,1) = MATMUL(TRANSPOSE(shape_functions(:,:,1)), element_coordinates(:,1))
    physical_points(:,2) = MATMUL(TRANSPOSE(shape_functions(:,:,1)), element_coordinates(:,2))
  ENDSUBROUTINE map_reference_to_physical

END MODULE element_mapping_module
