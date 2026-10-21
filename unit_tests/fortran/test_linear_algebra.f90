PROGRAM test_linear_algebra
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE LinearAlgebra, ONLY: mymatmul, matmul_dgemm, solve_linear_system, &
       solve_linear_system_sing, invert_matrix, compute_lup_decomp, eig, cross_product
  IMPLICIT NONE
  REAL*8 :: a(2,2), inverse(2,2), identity(2,2), x(2,2), rhs(2,2), solution(2,2), single(2)
  REAL*8 :: left(2,3), right(3,2), product(2,2), lower(2,2), upper(2,2)
  REAL*8 :: eigenvectors(2,2), eigenvalues(2,2), e1(3), e2(3), cross(3), reverse(3)
  INTEGER :: permutation(2,2)

  ! Rectangular products exercise DGEMM's dimensions and leading dimensions.
  left = RESHAPE([1.d0,2.d0,3.d0,4.d0,5.d0,6.d0], SHAPE(left))
  right = RESHAPE([2.d0,-1.d0,3.d0,4.d0,0.d0,-2.d0], SHAPE(right))
  CALL mymatmul(left,right,product)
  CALL check(product,MATMUL(left,right),'DGEMM subroutine')
  CALL check(matmul_dgemm(left,right),MATMUL(left,right),'DGEMM function')

  ! A zero leading entry requires pivoting. Test both DGESV RHS layouts.
  a = RESHAPE([0.d0,1.d0,2.d0,3.d0],SHAPE(a))
  x = RESHAPE([1.d0,-2.d0,3.d0,4.d0],SHAPE(x))
  rhs = MATMUL(a,x)
  CALL solve_linear_system(a,rhs,solution)
  CALL check(solution,x,'DGESV multiple RHS')
  CALL solve_linear_system_sing(a,rhs(:,1),single)
  CALL check(RESHAPE(single,[2,1]),x(:,:1),'DGESV single RHS')

  identity = 0.d0
  identity(1,1) = 1.d0
  identity(2,2) = 1.d0
  CALL invert_matrix(a,inverse)
  CALL check(MATMUL(a,inverse),identity,'DGETRF/DGETRI inverse')
  CALL compute_lup_decomp(a,2,lower,upper,permutation)
  CALL check(MATMUL(REAL(permutation,8),a),MATMUL(lower,upper),'DGETRF factors')

  ! Check real right-eigenvector residuals without fixing eigenvector signs.
  a = RESHAPE([2.d0,0.d0,1.d0,3.d0],SHAPE(a))
  CALL eig(a,eigenvectors,eigenvalues)
  CALL check(MATMUL(a,eigenvectors),MATMUL(eigenvectors,eigenvalues),'DGEEV residual')
  IF (ABS(eigenvalues(1,1)+eigenvalues(2,2)-5.d0) > 1.d-12) &
       ERROR STOP 'DGEEV eigenvalues'
  IF (MINVAL(SUM(eigenvectors**2,DIM=1)) < 0.5d0) ERROR STOP 'DGEEV eigenvectors'

  e1 = [1.d0,0.d0,0.d0]
  e2 = [0.d0,1.d0,0.d0]
  CALL cross_product(e1,e2,cross)
  CALL cross_product(e2,e1,reverse)
  CALL check(RESHAPE(cross,[3,1]),RESHAPE([0.d0,0.d0,1.d0],[3,1]),'cross product orientation')
  CALL check(RESHAPE(cross+reverse,[3,1]),RESHAPE([0.d0,0.d0,0.d0],[3,1]),'cross product antisymmetry')
  WRITE(*,*) 'Linear algebra: PASS; products, solves, inverse, LU, eigenvectors and cross products'
CONTAINS
  SUBROUTINE check(actual,expected,label)
    REAL*8, INTENT(IN) :: actual(:,:), expected(:,:)
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. ALL(ieee_is_finite(actual))) THEN
       WRITE(*,*) label
       ERROR STOP 'nonfinite result'
    ENDIF
    IF (MAXVAL(ABS(actual-expected)) > 1.d-12*MAX(1.d0,MAXVAL(ABS(expected)))) THEN
       WRITE(*,*) label
       ERROR STOP 'incorrect result'
    ENDIF
  END SUBROUTINE check
END PROGRAM test_linear_algebra
