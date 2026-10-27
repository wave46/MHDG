MODULE hdg_stabilization
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: add_tau_derivative
CONTAINS

  SUBROUTINE add_tau_derivative(matrix,rhs,row_dofs,trace_dofs,offsets,NN,N,trace,jump,dtau_du)
    REAL*8, INTENT(INOUT) :: matrix(:,:),rhs(:)
    INTEGER, INTENT(IN) :: row_dofs(:),trace_dofs(:),offsets(:)
    REAL*8, INTENT(IN) :: NN(:,:),N(:),trace(:),jump(:),dtau_du(:,:)
    INTEGER :: i,j,rows(SIZE(offsets)),cols(SIZE(offsets))
    REAL*8 :: coefficient(SIZE(trace))

    DO i = 1,SIZE(jump)
      IF (ALL(dtau_du(i,:) == 0.d0)) CYCLE
      rows = row_dofs(offsets+i)
      coefficient = jump(i)*dtau_du(i,:)
      DO j = 1,SIZE(trace)
        cols = trace_dofs(offsets+j)
        matrix(rows,cols) = matrix(rows,cols)+coefficient(j)*NN
      ENDDO
      rhs(rows) = rhs(rows)+DOT_PRODUCT(coefficient,trace)*N
    ENDDO
  END SUBROUTINE add_tau_derivative
END MODULE hdg_stabilization
