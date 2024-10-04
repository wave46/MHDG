!*****************************************
! project: MHDG
! file: LinearAlgebra.f90
! date: 04/01/2016
! Module collecting different routines for
! linear algebra operations
!*****************************************

MODULE LinearAlgebra

CONTAINS

  SUBROUTINE compute_lup_decomp(A, n, L, U, P)
    INTEGER*4, INTENT(IN)    :: n
    REAL*8, INTENT(IN)       :: A(:, :)
    REAL*8, INTENT(INOUT)    :: L(:, :), U(:, :)
    INTEGER*4, INTENT(INOUT) :: P(:, :)
    INTEGER*4                :: i, j, itemp(1:n)
    REAL*8                   :: temp(1:n, 1:n)

    temp = A
    CALL dgetrf(n, n, temp, n, itemp, i)
    L = temp
    U = temp
    DO i = 1, n
       L(i, i) = 1.d0
       DO j = 1, n
          IF (i .GT. j) THEN
             U(i, j) = 0.d0
          ELSEIF (i .LT. j) THEN
             L(i, j) = 0.d0
          END IF
       END DO
    END DO
    CALL find_permutation(n, P, itemp)

  END SUBROUTINE compute_lup_decomp

  SUBROUTINE find_permutation(n, P, per)
    INTEGER*4, INTENT(IN)    :: n
    INTEGER*4, INTENT(INOUT) :: per(1:n)
    INTEGER*4, INTENT(OUT)   :: P(1:n, 1:n)
    INTEGER*4                :: i, temp, auxper(1:n)

    P = 0
    auxper = (/(i, i=1, n)/)

    DO i = 1, n
       temp = per(i)
       per(i) = auxper(per(i))
       auxper(temp) = auxper(i)
    END DO

    DO i = 1, n
       P(i, per(i)) = 1
    END DO
    !  CALL check_permutation(P,A,L,U,n)
  END SUBROUTINE find_permutation

  SUBROUTINE check_permutation(P, A, L, U, n)
    INTEGER*4, INTENT(IN) :: n, P(1:n, 1:n)
    REAL*8, INTENT(IN)    :: A(1:n, 1:n), L(1:n, 1:n), U(1:n, 1:n)
    REAL*8                :: temp1(1:n, 1:n), temp2(1:n, 1:n)
    INTEGER*4             :: i, j

    temp1 = MATMUL(DBLE(P), A)
    temp2 = MATMUL(L, U)

    DO i = 1, n
       DO j = 1, n
          IF (temp1(i, j) - temp2(i, j) .GT. 1.d-10) THEN
             WRITE (*, *) 'temp1-temp2 = ', temp1(i, j) - temp2(i, j), ' at i,j = ', i, j
             WRITE (*, *) 'P'
             WRITE (*, *) P(1, :)
             WRITE (*, *) P(2, :)
             WRITE (*, *) P(3, :)
             WRITE (*, *) 'P*A'
             WRITE (*, *) temp1(1, :)
             WRITE (*, *) temp1(2, :)
             WRITE (*, *) temp1(3, :)
             WRITE (*, *) 'L*U'
             WRITE (*, *) temp2(1, :)
             WRITE (*, *) temp2(2, :)
             WRITE (*, *) temp2(3, :)
             error STOP "P*A != L*U"
          END IF
       END DO
    END DO

  END SUBROUTINE check_permutation

  !***************************
  ! Invert a square matrix
  !***************************
  SUBROUTINE invert_matrix(M, I)
    REAL*8, INTENT(IN)   :: M(:, :)
    REAL*8, INTENT(OUT)  :: I(:, :)
    INTEGER*4            :: n
    REAL*8               :: work(1:SIZE(M, 1)) ! work array for LAPACK
    INTEGER*4            :: ipiv(1:SIZE(M, 1)) ! pivot indices
    INTEGER*4            :: info

    EXTERNAL DGETRF
    EXTERNAL DGETRI

    n = SIZE(M, 1)

    ! store the matrix IN inverse to prevent it from being overwritten by LAPACK
    I = M

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    CALL DGETRF(n, n, I, n, ipiv, info)
    IF (info .NE. 0) THEN
       STOP 'Matrix is numerically singular!'
    END IF

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF
    CALL DGETRI(n, I, n, ipiv, work, n, info)
    IF (info .NE. 0) THEN
       STOP 'Matrix inversion failed!'
    END IF

  END SUBROUTINE invert_matrix

  !***************************
  ! Solve a linear system with
  ! multiple rhs
  !***************************
  SUBROUTINE solve_linear_system(A, b, x)
    REAL*8, INTENT(IN) :: A(:, :)
    REAL*8, INTENT(IN) :: b(:, :)
    REAL*8, INTENT(OUT):: x(:, :)
    INTEGER            :: n, m, info
    INTEGER            :: p(1:SIZE(b, 1))
    REAL*8,ALLOCATABLE :: temp(:,:)

    EXTERNAL DGESV

    ALLOCATE(temp(SIZE(A, 1), SIZE(A, 2)))
    n = SIZE(b, 1) ! size of the matrix (n x n)
    m = SIZE(b, 2) ! number of RHS
    temp = A
    x = b
    CALL DGESV(n, m, temp, n, p, x, n, info)

    IF (INFO .GT. 0) THEN
       WRITE (*, *) 'The diagonal element of the triangular factor of A,'
       WRITE (*, *) 'U(', INFO, ',', INFO, ') is zero, so that'
       WRITE (*, *) 'A is singular; the solution could not be computed.'
       STOP
    END IF
    DEALLOCATE(temp)

  END SUBROUTINE solve_linear_system

  !***************************
  ! Solve a linear system with
  ! single rhs
  !***************************
  SUBROUTINE solve_linear_system_sing(A, b, x)
    REAL*8, INTENT(IN) :: A(:, :)
    REAL*8, INTENT(IN) :: b(:)
    REAL*8, INTENT(OUT):: x(:)
    INTEGER            :: n, m, info
    INTEGER            :: p(1:SIZE(b, 1))
    REAL*8             :: temp(SIZE(A, 1), SIZE(A, 2))
    EXTERNAL DGESV

    n = SIZE(b, 1) ! size of the matrix (n x n)
    m = 1 ! number of RHS
    temp = A
    x = b
    CALL DGESV(n, m, temp, n, p, x, n, info)

    IF (INFO .GT. 0) THEN
       WRITE (*, *) 'The diagonal element of the triangular factor of A,'
       WRITE (*, *) 'U(', INFO, ',', INFO, ') is zero, so that'
       WRITE (*, *) 'A is singular; the solution could not be computed.'
       STOP
    END IF

  END SUBROUTINE solve_linear_system_sing

  !***************************
  ! Matrix multiplication using
  ! dgemm (avoid stack overflow)
  !***************************
  SUBROUTINE mymatmul(A,B,C)
    REAL*8, INTENT(IN) :: A(:, :)
    REAL*8, INTENT(IN) :: B(:,:)
    REAL*8, INTENT(OUT):: C(:,:)
    INTEGER            :: a1,a2,b1,b2,c1,c2
    EXTERNAL DGEMM

    a1 = SIZE(A,1)
    a2 = SIZE(A,2)
    b1 = SIZE(B,1)
    b2 = SIZE(B,2)
    c1 = SIZE(C,1)
    c2 = SIZE(C,2)

    ! Little check
    IF (a2.NE.b1) THEN
       WRITE(6,*) "Error: number of colums of A different from number of rows of B"
       STOP
    ENDIF
    IF (a1.NE.c1) THEN
       WRITE(6,*) "Error: number of rows of A different from number of rows of C"
       STOP
    ENDIF
    IF (b2.NE.c2) THEN
       WRITE(6,*) "Error: number of colums of B different from number of columns of C"
       STOP
    ENDIF

    CALL DGEMM('N','N',a1,b2,a2,1.,A,a1,B,a2,0.,C,a1)


  END SUBROUTINE mymatmul



  !***************************
  ! Compute eigenvalues and
  ! right eigenvectors
  !***************************
  SUBROUTINE eig(A, V, D)
    REAL*8, INTENT(IN)  :: A(:, :)
    REAL*8, INTENT(OUT) :: V(:, :)
    REAL*8, INTENT(OUT) :: D(:, :)
    INTEGER             :: m, lwork, i, info
    INTEGER, PARAMETER  :: nb = 64
    REAL*8,ALLOCATABLE  :: work(:)
    REAL*8, ALLOCATABLE :: Vmat(:, :), wr(:), wi(:)
    REAL*8              :: dummy(SIZE(A,1), SIZE(A,1))

    EXTERNAL dgeev

    m = SIZE(A, 1)

    ALLOCATE (wr(1:m), wi(1:m))
    ALLOCATE (Vmat(1:m, 1:m))
    D = 0.d0
    Vmat = A

    ! Query the optimal workspace
    lwork = -1
    CALL dgeev('N', 'V', m, Vmat, m, wr, wi, dummy, 1, V, m, dummy, lwork, info)

    ! Compute the eigenvalues and right eigenvectors of A
    lwork = MAX((nb+2)*m, NINT(dummy(1,1)))
    ALLOCATE (work(lwork))
    CALL dgeev('N', 'V', m, Vmat, m, wr, wi, dummy, 1, V, m, work, lwork, info)

    ! Only real eigenvalues are supposed
    DO i = 1, m
       D(i, i) = wr(i)
    END DO
    DEALLOCATE (wr, wi, Vmat,work)
  END SUBROUTINE eig

  !*******************************************
  ! Tensor product
  !*******************************************
  FUNCTION TensorProduct(Vl, Vr) RESULT(Tens)
    REAL, DIMENSION(:)                  :: Vl, Vr
    REAL, DIMENSION(SIZE(Vl), SIZE(Vr)) :: Tens

    INTEGER :: i, j, Nl, Nr

    Nl = SIZE(Vl); Nr = SIZE(Vr)

    DO j = 1, Nr
       DO i = 1, Nl
          Tens(i, j) = Vl(i)*Vr(j)
       END DO
    END DO
  END FUNCTION TensorProduct

  !*******************************************
  ! Tensor product cumulative
  !*******************************************
  SUBROUTINE TensorProductCumul(M, Vl, Vr)
    REAL, DIMENSION(:, :) :: M
    REAL, DIMENSION(:)    :: Vl, Vr

    INTEGER :: i, j, Nl, Nr

    Nl = SIZE(Vl); Nr = SIZE(Vr)

    DO j = 1, Nr
       DO i = 1, Nl
          M(i, j) = M(i, j) + Vl(i)*Vr(j)
       END DO
    END DO
  END SUBROUTINE TensorProductCumul

  !*******************************************
  ! Tensor sum for integers
  !*******************************************
  FUNCTION TensorSumInt(Vl, Vr) RESULT(Tens)
    INTEGER, DIMENSION(:)                  :: Vl, Vr
    INTEGER, DIMENSION(SIZE(Vl), SIZE(Vr)) :: Tens

    INTEGER :: i, j, Nl, Nr

    Nl = SIZE(Vl); Nr = SIZE(Vr)

    DO j = 1, Nr
       DO i = 1, Nl
          Tens(i, j) = Vl(i) + Vr(j)
       END DO
    END DO
  END FUNCTION TensorSumInt

  !****************************************************
  ! Transform matrices IN vectors for integers: colint
  !****************************************************
  FUNCTION colint(M) RESULT(V)
    INTEGER, DIMENSION(:, :)                  :: M
    INTEGER, DIMENSION(SIZE(M, 1)*SIZE(M, 2)) :: V

    INTEGER :: i, j

    DO j = 1, SIZE(M, 2)
       DO i = 1, SIZE(M, 1)
          V((j - 1)*SIZE(M, 1) + i) = M(i, j)
       END DO
    END DO
  END FUNCTION colint

  !****************************************************
  ! Transform matrices IN vectors for reals: col
  !****************************************************
  FUNCTION col(M) RESULT(V)
    REAL, DIMENSION(:, :)                  :: M
    REAL, DIMENSION(SIZE(M, 1)*SIZE(M, 2)) :: V

    INTEGER :: i, j

    DO j = 1, SIZE(M, 2)
       DO i = 1, SIZE(M, 1)
          V((j - 1)*SIZE(M, 1) + i) = M(i, j)
       END DO
    END DO
  END FUNCTION col

END MODULE LinearAlgebra

SUBROUTINE cross_product(a, b, c)
  REAL*8, INTENT(IN) :: a(3), b(3)
  REAL*8, INTENT(OUT):: c(3)

  c = 0.
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)

END SUBROUTINE cross_product

SUBROUTINE ijk_cross_product(i, j, k)
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(OUT):: j, k

  SELECT CASE (i)
  CASE (1)
     j = 3
     k = 2
  CASE (2)
     j = 1
     k = 3
  CASE (3)
     j = 2
     k = 1
  CASE DEFAULT
     WRITE (*, *) "Error, wrong case IN cross-product"
  END SELECT

END SUBROUTINE ijk_cross_product
