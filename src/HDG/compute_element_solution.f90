!*****************************************
! project: MHDG
! file: compute_solution.f90
! date: 24/01/2017
! Solve the global system, compute face
! and element solution
!*****************************************

SUBROUTINE compute_element_solution
  USE globals
  USE LinearAlgebra
  USE printUtils
  USE solve_pastix
  USE MPI_OMP
  IMPLICIT NONE

  !**********************************
  ! Retrieve the element solution
  !**********************************
#ifdef TOR3D
  !***********************************************************************
  !                            VERSION 3D TOROIDAL
  !***********************************************************************

  INTEGER               :: itor, iel, iel3d, ifa, i
  INTEGER               :: Ndim, Neq, N2D, Np2D, Np, Nfl, Ntorloc, Np1Dpol, Np1Dtor, Nf, Nfdir
  INTEGER               :: ind_ue(refElTor%Nnodes3D*phys%neq), ind_uf(refElTor%Nft*phys%neq)
  INTEGER*4             :: ind_ug(refElTor%Nnodes3D*phys%Neq*3)
  INTEGER               :: dd, delta(refElPol%Nfaces + 2)
  INTEGER               :: ind_dim(refElPol%Nfaces + 2), ind_sta(refElPol%Nfaces + 2), aux
  INTEGER*4             :: Fe(refElPol%Nfaces)

  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tps1)
     CALL system_CLOCK(timing%cks1, timing%clock_rate1)
  END IF


  Ndim = 3                                               ! Number of dimensions
  Neq = phys%neq                                        ! Number of equations
  N2D = Mesh%Nelems                                     ! Number elements in the 2D mesh
  Np2D = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
  Np1Dpol = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Np = refElTor%Nnodes3D                               ! Number of nodes for each 3D element
  Nfl = refElTor%Nfl                                    ! Number of nodes in the lateral faces
  Nf = refElPol%Nfaces
  Nfdir = Mesh%Ndir

  ind_uf = 0
  ind_dim(1) = Np2D*Neq
  ind_dim(2:refElPol%Nfaces + 1) = Nfl*Neq
  ind_dim(refElPol%Nfaces + 2) = Np2D*Neq
  ind_sta(1) = 1
  aux = 0
  DO i = 2, refElPol%Nfaces + 2
    ind_sta(i) = 1 + ind_dim(i - 1) + aux
    aux = aux + ind_dim(i - 1)
  END DO

#ifdef PARALL
  IF (MPIvar%ntor .GT. 1) THEN
    ntorloc = numer%ntor/MPIvar%ntor + 1
  ELSE
    ntorloc = numer%ntor
  ENDIF
#else
  ntorloc = numer%ntor
#endif

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(itor,iel,iel3d,ind_ue,ind_uf,Fe,dd,delta,ifa)
  !$OMP DO
  DO itor = 1, ntorloc
    DO iel = 1, N2d
      iel3d = (itor - 1)*N2d+iel  !3d numbering of the element
      ind_ue = (iel3d-1)*Np*Neq + (/(i, i=1, Np*Neq)/)

      ! Index for the face solution
      Fe = Mesh%F(iel, :)
      dd = 1 + (itor - 1)*Neq*(N2d*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)
      delta(1) = dd + (iel - 1)*Np2D*Neq
      delta(2:refElPol%Nfaces + 1) = dd + N2d*Np2D*Neq + (Fe - 1)*Nfl*Neq
      delta(refElPol%Nfaces + 2) = dd + N2d*Np2D*Neq + (Mesh%Nfaces - Nfdir)*Nfl*Neq + (iel - 1)*Np2D*Neq
        IF (MPIvar%ntor .EQ. 1) THEN
           IF (itor == numer%ntor) THEN
          delta(refElPol%Nfaces + 2) = 1 + (iel - 1)*Np2D*Neq
           END IF
        ENDIF

        DO ifa = 1, Nf + 2
           DO i = 0, ind_dim(ifa) - 1
          ind_uf(i + ind_sta(ifa)) = delta(ifa) + i
           END DO
        END DO

      ! elemental solutions
        sol%u(ind_ue) = MATMUL(elMat%UU(:, :, iel3d), sol%u_tilde(ind_uf)) + elMat%U0(:, iel3d)

      ! Index for the element gradient
      ind_ug = (iel3d-1)*neq*Np*Ndim + (/(i, i=1, Neq*Np*Ndim)/)

      ! elemental solutions
        sol%q(ind_ug) = MATMUL(elMat%LL(:, :, iel), sol%u_tilde(ind_uf)) + elMat%L0(:, iel)
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

#else
  !***********************************************************************
  !                            VERSION 2D
  !***********************************************************************
  INTEGER*4             :: i, iel, N2d, Np, Nfp, Neq, Nf, Ndim
  INTEGER*4             :: ind_ue(Mesh%Nnodesperelem*phys%Neq)
  INTEGER*4             :: ind_uf(refElPol%Nfacenodes*phys%Neq*refElPol%Nfaces)
  INTEGER*4             :: ind_ug(Mesh%Nnodesperelem*phys%Neq*Mesh%Ndim)
  INTEGER*4             :: Fe(refElPol%Nfaces)

  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tps1)
     CALL system_CLOCK(timing%cks1, timing%clock_rate1)
  END IF

  Ndim = 2
  N2D = Mesh%Nelems
  Np = Mesh%Nnodesperelem
  neq = phys%neq
  Nfp = Mesh%Nnodesperface
  Nf = refElPol%Nfaces

  DO iel = 1, N2d
    ! Index for the face solution
    Fe = Mesh%F(iel, :)
     ind_uf = RESHAPE(tensorSumInt((/(i, i=1, Neq*Nfp)/), (Fe - 1)*Neq*Nfp), (/Neq*Nfp*Nf/))

    ! Index for the element solution
    ind_ue = (iel - 1)*neq*Np + (/(i, i=1, Neq*Np)/)

    ! elemental solutions
     sol%u(ind_ue) = MATMUL(elMat%UU(:, :, iel), sol%u_tilde(ind_uf)) + elMat%U0(:, iel)

    ! Index for the element gradient
    ind_ug = (iel - 1)*neq*Np*Ndim + (/(i, i=1, Neq*Np*Ndim)/)

    ! elemental solutions
     sol%q(ind_ug) = MATMUL(elMat%LL(:, :, iel), sol%u_tilde(ind_uf)) + elMat%L0(:, iel)
  END DO

#endif


  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tpe1)
     CALL system_CLOCK(timing%cke1, timing%clock_rate1)
     timing%runtsol = timing%runtsol + (timing%cke1 - timing%cks1)/REAL(timing%clock_rate1)
    timing%cputsol = timing%cputsol + timing%tpe1 - timing%tps1
  END IF

END SUBROUTINE compute_element_solution

