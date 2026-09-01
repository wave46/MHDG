!*****************************************
! project: MHDG
! file: residual_norm.f90
! Compute global residual norms over owned
! volume unknowns.
!*****************************************

MODULE residual_norm_module

  USE globals, ONLY: Mesh, phys, numer, refElTor
  USE MPI_OMP, ONLY: MPIvar
#ifdef PARALL
  USE mpi
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: computeResidual

CONTAINS

  FUNCTION computeResidual(u, uref, coeff) RESULT(res)
    REAL*8, INTENT(IN) :: u(:), uref(:), coeff
    REAL*8             :: res, sum2
    INTEGER            :: nglo
#ifdef PARALL
    INTEGER            :: ierr, iel, ind_start, ind_end, nunk_elem, nloc
#ifdef TOR3D
    INTEGER            :: itor, ntor_owned
#endif
#endif

#ifdef PARALL
    sum2 = 0.0D0
    nloc = 0

#ifdef TOR3D
    nunk_elem = phys%Neq*refElTor%Nnodes3D
    ntor_owned = numer%ntor
    IF (MPIvar%ntor .GT. 1) ntor_owned = numer%ntor/MPIvar%ntor

    DO itor = 1, ntor_owned
       DO iel = 1, Mesh%Nelems
          IF (Mesh%ghostElems(iel) .NE. 0) CYCLE
          ind_start = ((itor - 1)*Mesh%Nelems + iel - 1)*nunk_elem + 1
          ind_end = ind_start + nunk_elem - 1
          sum2 = sum2 + SUM((u(ind_start:ind_end) - uref(ind_start:ind_end))**2)
          nloc = nloc + nunk_elem
       END DO
    END DO
#else
    nunk_elem = phys%Neq*Mesh%Nnodesperelem
    DO iel = 1, Mesh%Nelems
       IF (Mesh%ghostElems(iel) .NE. 0) CYCLE
       ind_start = (iel - 1)*nunk_elem + 1
       ind_end = ind_start + nunk_elem - 1
       sum2 = sum2 + SUM((u(ind_start:ind_end) - uref(ind_start:ind_end))**2)
       nloc = nloc + nunk_elem
    END DO
#endif

    nglo = nloc
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, sum2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, nglo, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    sum2 = SUM((u - uref)**2)
    nglo = SIZE(u)
#endif

    res = SQRT(sum2)/SQRT(DBLE(nglo))/coeff
  END FUNCTION computeResidual

END MODULE residual_norm_module
