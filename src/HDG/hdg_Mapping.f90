SUBROUTINE HDG_mapping()
  USE globals
  USE LinearAlgebra
  USE printUtils
  USE MPI_OMP
  USE in_out

  IMPLICIT NONE
  INTEGER*4             :: Ndim, N2d, neq, Np, Nfp, Nfe, i, ifa, iel, iel3, Nfg
#ifdef TOR3D
  INTEGER*4             :: ind(refElPol%Nfaces, phys%neq*refElTor%Nfl), perm(refElTor%Nfl*phys%Neq)
  INTEGER*4             :: j, Nfl, itor, Np1Dtor, Np1Dpol, Np2D, ntorloc
#else
  INTEGER*4             :: ind(refElPol%Nfaces, phys%neq*Mesh%Nnodesperface)
  INTEGER*4             :: perm(refElPol%Nfacenodes*phys%Neq)
#endif
  REAL*8, ALLOCATABLE    :: LL(:, :), UU(:, :), L0(:), U0(:)
  REAL*8, ALLOCATABLE    :: Auq_iAqq(:, :), Mu(:, :), Ml(:, :), Md(:), aux(:, :), auxAqq(:, :), auxAquUU(:, :)
  REAL*8, POINTER        :: iAqq(:, :), Aqu(:, :), Aql(:, :)
  REAL*8, POINTER        :: Auq(:, :), Auu(:, :), Aul(:, :)
  REAL*8, POINTER        :: Aql_dir(:), Aul_dir(:), f(:)

  IF (MPIvar%glob_id .EQ. 0) THEN
     IF (utils%printint > 1) THEN
        WRITE (6, *) '*************************************************'
        WRITE (6, *) '*            MAPPING                            *'
        WRITE (6, *) '*************************************************'
     END IF
  END IF

  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tps1)
     CALL system_CLOCK(timing%cks1, timing%clock_rate1)
  END IF

  Neq = phys%neq
#ifdef TOR3D
  Ndim = 3                                               ! Number of dimensions
  N2D = Mesh%Nelems                                     ! Number elements in the 2D mesh
  Np2D = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
  Np1Dpol = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Np = Np2D*refElTor%Nnodes1D                          ! Number of nodes for each 3D element
  Nfl = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
  Nfg = Np2D*2 + refElPol%Nfaces*Nfl                      ! Number of nodes in all the faces of a 3D element
  Nfe = refElPol%Nfaces                                 ! Number of faces in the 2D element
#else
  Ndim = 2
  N2D = Mesh%Nelems
  Np = Mesh%Nnodesperelem
  Nfe = refElPol%Nfaces
  Nfg = Mesh%Nnodesperface*Nfe
  Nfp = Mesh%Nnodesperface
#endif

#ifdef TOR3D
  ind = 0
  DO i = 1, Nfe
     DO j = 1, Nfl*Neq
        ind(i, j) = refElPol%Nnodes2D*Neq + Neq*Nfl*(i - 1) + j
     END DO
  END DO

  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Np1Dpol, Np1Dtor, Neq, perm)
#else
  ! local assembly indexes
  ind = 0
  DO ifa = 1, Nfe
     ind(ifa, :) = neq*Nfp*(ifa - 1) + ((/(i, i=1, neq*Nfp)/))
  END DO

  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Neq*Nfp, Neq, perm)
#endif

#ifdef TOR3D
#ifdef PARALL
  IF (MPIvar%ntor .GT. 1) THEN
     ntorloc = numer%ntor/MPIvar%ntor + 1
  ELSE
     ntorloc = numer%ntor
  ENDIF
#else
  ntorloc = numer%ntor
#endif
#endif

  !$OMP PARALLEL DEFAULT(SHARED) &
#ifdef TOR3D
  !$OMP PRIVATE(iel,itor,iel3,LL,UU,L0,U0,Auq_iAqq,Auu,Mu,Ml,Md,f,aux,auxAqq,auxAquUU,iAqq,Aqu,Aql,Auq,Aul,Aql_dir,Aul_dir,ifa)
#else
  !$OMP PRIVATE(iel,iel3,LL,UU,L0,U0,Auq_iAqq,Auu,Mu,Ml,Md,f,aux,auxAqq,auxAquUU,iAqq,Aqu,Aql,Auq,Aul,Aql_dir,Aul_dir,ifa)
#endif

  ALLOCATE (LL(Neq*Np*Ndim, Neq*Nfg))
  ALLOCATE (L0(Neq*Np*Ndim))
  ALLOCATE (UU(Neq*Np, Neq*Nfg))
  ALLOCATE (U0(Neq*Np))
  ALLOCATE (Auq_iAqq(Neq*Np, Neq*Ndim*Np))
  ALLOCATE (Ml(Neq*Np, Neq*Nfg))
  ALLOCATE (Mu(Neq*Np, Neq*Np))
  ALLOCATE (Md(Neq*Np))
  ALLOCATE (aux(Neq*Ndim*Np, Neq*Np))
  ALLOCATE (auxAqq(Neq*Ndim*Np, Neq*Ndim*Np))
  ALLOCATE (auxAquUU(Neq*Ndim*Np, Neq*Nfg))


  !*****************
  ! Loop in elements
  !*****************
#ifdef TOR3D
  !$OMP DO SCHEDULE(STATIC) &
#else
  !$OMP DO SCHEDULE(STATIC)
#endif

#ifdef TOR3D
  !$OMP COLLAPSE(2)
  DO itor = 1, ntorloc
#endif
     DO iel = 1, N2d

#ifdef TOR3D
        iel3 = (itor - 1)*N2d+iel
#else
        iel3 = iel
#endif
        iAqq => elmat%iAqq(:, :, iel3)
        Aqu => elmat%Aqu(:, :, iel3)
        Aql => elmat%Aql(:, :, iel3)
        Auq => elmat%Auq(:, :, iel3)
        Auu => elmat%Auu(:, :, iel3)
        Aul => elmat%Aul(:, :, iel3)
        Aql_dir => elmat%Aql_dir(:, iel3)
        Aul_dir => elmat%Aul_dir(:, iel3)
        f => elmat%S(:, iel3)

        ! Auq_iAqq=matmul(Auq,iAqq)
        CALL mymatmul(Auq,iAqq,Auq_iAqq)
        !CALL DGEMM('N','N',Neq*Np,Neq*Ndim*Np,Neq*Ndim*Np,1.,Auq,Neq*Np,Auq,Neq*Ndim*Np,1.,Auq_iAqq,Neq*Np)

        !Mu = Auu - matmul(Auq_iAqq, Aqu)
        CALL mymatmul(Auq_iAqq,Aqu,Mu)
        Mu = Auu - Mu

        !Ml = Aul - matmul(Auq_iAqq, Aql)
        CALL  mymatmul(Auq_iAqq,Aql,Ml)
        Ml = Aul - Ml

        Md = f - (Aul_dir - MATMUL(Auq_iAqq, Aql_dir))
        CALL solve_linear_system(Mu, -Ml, UU)
        CALL solve_linear_system_sing(Mu, Md, U0)

        !auxAquUU = matmul(Aqu, UU) ! Avoid stack overflow
        CALL mymatmul(Aqu,UU,auxAquUU) ! Avoid stack overflow
        auxAquUU = auxAquUU + Aql

        !LL = -matmul(iAqq, auxAquUU)
        CALL mymatmul(iAqq,auxAquUU,LL)
        LL = -LL

        L0 = -MATMUL(iAqq, Aql_dir + (MATMUL(Aqu, U0)))

        ! Flip the faces
        DO ifa = 1, refElPol%Nfaces
           IF (Mesh%flipface(iel, ifa)) THEN
              LL(:, ind(ifa, :)) = LL(:, ind(ifa, perm(:)))
              UU(:, ind(ifa, :)) = UU(:, ind(ifa, perm(:)))
           END IF
        END DO

        ! Store the mapping
        elMat%LL(:, :, iel3) = LL
        elMat%UU(:, :, iel3) = UU
        elMat%L0(:, iel3) = L0
        elMat%U0(:, iel3) = U0
        NULLIFY (iAqq, Aqu, Aql, Auq, Aul, Aql_dir, Aul_dir, f)

     END DO
#ifdef TOR3D
  END DO
#endif

  !$OMP END DO
  DEALLOCATE (LL, UU, L0, U0, Auq_iAqq, Mu, Ml, Md, aux, auxAqq, auxAquUU)
  !$OMP END PARALLEL


  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tpe1)
     CALL system_CLOCK(timing%cke1, timing%clock_rate1)
     timing%runtmap = timing%runtmap + (timing%cke1 - timing%cks1)/REAL(timing%clock_rate1)
     timing%cputmap = timing%cputmap + timing%tpe1 - timing%tps1
  END IF


CONTAINS

#ifdef TOR3D
  !*****************************************
  ! Set permutations for flipping faces
  !****************************************
  SUBROUTINE set_permutations(Np1Dpol, Np1Dtor, Neq, perm)
    INTEGER, INTENT(IN)  :: Np1Dpol, Np1Dtor, Neq
    INTEGER, INTENT(OUT) :: perm(:)
    INTEGER              :: i, j, k
    INTEGER              :: temp_pol(Np1Dpol), temp_tor(Np1Dtor), aux(Np1Dpol, Np1Dtor)

    temp_pol = (/(i, i=1, Np1Dpol)/)
    temp_tor = Np1Dpol*((/(i, i=1, Np1Dtor)/) - 1)
    aux = TensorSumInt(temp_pol, temp_tor)
    DO j = 1, Np1Dtor
       DO i = 1, Np1Dpol
          DO k = 1, Neq
             perm((j - 1)*Np1Dpol*Neq + (i - 1)*Neq + k) = (aux(Np1Dpol - i + 1, j) - 1)*Neq + k
          END DO
       END DO
    END DO

  END SUBROUTINE set_permutations

#else

  !*****************************************
  ! Set permutations for flipping faces
  !****************************************
  SUBROUTINE set_permutations(n, m, perm)
    INTEGER*4, INTENT(IN)  :: n, m
    INTEGER*4, INTENT(OUT) :: perm(:)
    INTEGER*4              :: i
    INTEGER*4              :: temp(m, n/m), templr(m, n/m)

    IF (MOD(n, m) .NE. 0) THEN
       WRITE (6, *) 'Error! n must be a multiple of m'
       STOP
    END IF

    templr = 0
    temp = RESHAPE((/(i, i=1, n)/), (/m, n/m/))
    DO i = 1, n/m
       templr(:, i) = temp(:, n/m - i + 1)
    END DO
    perm = RESHAPE(templr, (/n/))

  END SUBROUTINE set_permutations
#endif
END SUBROUTINE HDG_mapping
