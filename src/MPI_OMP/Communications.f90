!*****************************************
! project: MHDG
! file: Communications.f90
! date: 10/03/2017
! Exchange solution between ghost faces
!*****************************************
MODULE Communications
  USE MPI_OMP
  USE globals
  USE PrintUtils

CONTAINS

  SUBROUTINE init_Com()
    INTEGER, PARAMETER                   :: etq = 100
    INTEGER, DIMENSION(MPI_STATUS_SIZE)  :: stat
    INTEGER                             :: code
    INTEGER                             :: npro, rbuf, lpro, i, nf2sd, ct
    INTEGER                             :: psd, prv, ifa, fcount
    INTEGER, DIMENSION(MPIvar%glob_size) :: proext, proint
    INTEGER, ALLOCATABLE                :: faces_ext(:), faces_int(:)!,connpro(:)

#ifdef TOR3D
    INTEGER                             :: ecount, ne2sd, req
    INTEGER, ALLOCATABLE                :: elems_ext(:), elems_int(:)
#endif


    ! Number of processes
    npro = MPIvar%glob_size

    ! Local process
    lpro = MPIvar%glob_id

    ! faces and processes to receive
    ALLOCATE (Mesh%pr2rv(Mesh%nghostfaces))
    ALLOCATE (Mesh%fc2rv(Mesh%nghostfaces))
    Mesh%pr2rv = Mesh%ghostpro
    Mesh%fc2rv = 0
    ct = 1
    DO i = 1, Mesh%Nfaces
       IF (Mesh%ghostfaces(i) .EQ. 1) THEN
          Mesh%fc2rv(ct) = i
          ct = ct + 1
       END IF
    END DO
    ! ***********************************************************
    ! proext contains the number of ghost faces that each process
    ! compute for the local process, proint contains the number of
    ! ghost faces that the local process compute for the other
    ! processes
    ! ***********************************************************
    proext = 0
    proint = 0
    DO i = 1, Mesh%nghostfaces
       IF (Mesh%ghostpro(i).EQ.0) CYCLE ! Add by Benjamin
       proext(Mesh%ghostpro(i)) = proext(Mesh%ghostpro(i)) + 1
    END DO
    DO i = 1, npro
       rbuf = 0
       CALL MPI_Scatter(proext, 1, MPI_INTEGER, rbuf, 1, MPI_INTEGER, i - 1, MPI_COMM_WORLD, code)
       proint(i) = rbuf
    END DO

    ! Number of faces to send
    nf2sd = 0
    DO i = 1, npro
       nf2sd = nf2sd + proint(i)
    END DO
    ALLOCATE (Mesh%fc2sd(nf2sd))
    ALLOCATE (Mesh%pr2sd(nf2sd))

    ! Communicate to each process the local number of the faces he needs to send,
    ! and to whom
    fcount = 1
    DO i = 1, npro

       IF (MPIvar%glob_id .EQ. (i - 1)) CYCLE
       ! Process to send
       IF (proext(i) .EQ. 0) THEN
          psd = MPI_PROC_NULL
       ELSE
          psd = i - 1
       END IF
       ! Process to receive
       IF (proint(i) .EQ. 0) THEN
          prv = MPI_PROC_NULL
       ELSE
          prv = i - 1
       END IF
       ALLOCATE (faces_ext(proext(i)))
       ALLOCATE (faces_int(proint(i)))
       faces_ext = 0
       faces_int = 0

       ! prepare the vector with the index of the faces the local
       ! process need to receive from a connected process
       ct = 1
       DO ifa = 1, Mesh%nghostfaces
          IF (Mesh%ghostpro(ifa) .EQ. i) THEN
             faces_ext(ct) = Mesh%ghostLoc(ifa)
             ct = ct + 1
          END IF
       END DO

       CALL MPI_SEND(faces_ext, proext(i), MPI_INTEGER, psd, etq, MPI_COMM_WORLD, code)
       CALL MPI_RECV(faces_int, proint(i), MPI_INTEGER, prv, etq, MPI_COMM_WORLD, stat, code)

       Mesh%fc2sd(fcount:fcount + proint(i) - 1) = faces_int
       Mesh%pr2sd(fcount:fcount + proint(i) - 1) = i
       fcount = fcount + proint(i)

       DEALLOCATE (faces_ext, faces_int)
    END DO

#ifdef TOR3D
    ! elements and processes to receive
    IF (Mesh%nghostelems .GT. 0) THEN
       ALLOCATE (Mesh%pe2rv(Mesh%nghostelems))
       ALLOCATE (Mesh%el2rv(Mesh%nghostelems))
       Mesh%pe2rv = Mesh%ghelspro
       Mesh%el2rv = 0
       ct = 1
       DO i = 1, Mesh%Nelems
          IF (Mesh%ghostelems(i) .EQ. 1) THEN
             Mesh%el2rv(ct) = i
             ct = ct + 1
          END IF
       END DO
    ENDIF
    ! ***********************************************************
    ! proext contains the number of ghost elements that each process
    ! compute for the local process, proint contains the number of
    ! ghost elements that the local process compute for the other
    ! processes
    ! ***********************************************************
    proext = 0
    proint = 0
    DO i = 1, Mesh%nghostelems
       proext(Mesh%ghelspro(i)) = proext(Mesh%ghelspro(i)) + 1
    END DO

    DO i = 1, npro
       rbuf = 0
       CALL MPI_Scatter(proext, 1, MPI_INTEGER, rbuf, 1, MPI_INTEGER, i - 1, MPI_COMM_WORLD, code)
       proint(i) = rbuf
    END DO

    ! Number of elements to send
    ne2sd = 0
    DO i = 1, npro
       ne2sd = ne2sd + proint(i)
    END DO
    ALLOCATE (Mesh%el2sd(ne2sd))
    ALLOCATE (Mesh%pe2sd(ne2sd))

    ! Communicate to each process the local number of the faces he needs to send,
    ! and to whom
    ecount = 1
    DO i = 1, npro

       IF (MPIvar%glob_id .EQ. (i - 1)) CYCLE
       ! Process to send
       IF (proext(i) .EQ. 0) THEN
          psd = MPI_PROC_NULL
       ELSE
          psd = i - 1
       END IF
       ! Process to receive
       IF (proint(i) .EQ. 0) THEN
          prv = MPI_PROC_NULL
       ELSE
          prv = i - 1
       END IF
       ALLOCATE (elems_ext(proext(i)))
       ALLOCATE (elems_int(proint(i)))
       elems_ext = 0
       elems_int = 0

       ! prepare the vector with the index of the faces the local
       ! process need to receive from a connected process
       ct = 1
       DO ifa = 1, Mesh%nghostelems
          IF (Mesh%ghelspro(ifa) .EQ. i) THEN
             elems_ext(ct) = Mesh%ghelsLoc(ifa)
             ct = ct + 1
          END IF
       END DO

       CALL MPI_SEND(elems_ext, proext(i), MPI_INTEGER, psd, etq, MPI_COMM_WORLD, code)
       CALL MPI_RECV(elems_int, proint(i), MPI_INTEGER, prv, etq, MPI_COMM_WORLD, stat, code)

       Mesh%el2sd(ecount:ecount + proint(i) - 1) = elems_int
       Mesh%pe2sd(ecount:ecount + proint(i) - 1) = i
       ecount = ecount + proint(i)

       DEALLOCATE (elems_ext, elems_int)
    END DO
#endif
    ! call syncroprint_vector_int(Mesh%fc2sd)
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! call syncroprint_vector_int(Mesh%pr2sd)
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! call syncroprint_vector_int(Mesh%fc2rv)
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! stop
  END SUBROUTINE init_com

#ifdef TOR3D
  SUBROUTINE exchangeSol()
    INTEGER, PARAMETER  :: etq = 100
    INTEGER            :: Neq, Nfl, Np2d, i, j, Fi, itor, N2d, Nfdir
    INTEGER            :: indt(refElTor%Nfl*phys%Neq)
    INTEGER            :: indp(refElPol%Nnodes2D*phys%Neq)
    INTEGER            :: nf2sd, ne2sd
    REAL*8, ALLOCATABLE :: buffrv(:, :), buffsd(:, :)
    REAL*8             :: fbufsd(refElTor%Nfl*phys%Neq)
    REAL*8             :: fbufrv(refElTor%Nfl*phys%Neq)
    INTEGER            :: code
    INTEGER, ALLOCATABLE:: req(:), stat(:, :)
    INTEGER            :: dd, delta
    INTEGER            :: prtorsd, prtorrv, iel, ifa
    INTEGER :: Ntorloc,auxel(refElPol%Nnodes2D*phys%neq),auxfl(phys%Neq*refElTor%Nfl)
    ! integer            :: perm(1:refElPol%Nfacenodes*phys%Neq)

#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       ntorloc = numer%ntor/MPIvar%ntor + 1
       !       ntorass = ntorloc-1
    ELSE
       ntorloc = numer%ntor
       !       ntorass = ntorloc
    ENDIF
#else
    ntorloc = numer%ntor
    ntorass = ntorloc
#endif
    Neq = phys%Neq
    Nfl = refElTor%Nfl
    Np2D = refElPol%Nnodes2D
    N2d = Mesh%Nelems
    Nfdir = Mesh%Ndir
    nf2sd = SIZE(Mesh%fc2sd)
    ne2sd = SIZE(Mesh%el2sd)

    !************************************************************************************************
    !
    !       COMMUNICATIONS FOR POLOIDAL DISTRIBUTION OF THE MESH
    !
    !************************************************************************************************
    !*************************************************
    !  Communication for poloidal faces
    !*************************************************
    ALLOCATE (req(ne2sd + Mesh%nghostelems))
    ALLOCATE (stat(MPI_STATUS_SIZE, ne2sd + Mesh%nghostelems))
    req = mpi_request_null

    ! Allocate buffers
    ALLOCATE (buffrv(Np2D*phys%Neq, Mesh%nghostelems))
    ALLOCATE (buffsd(Np2D*phys%Neq, ne2sd))

    auxel = (/(j, j=0, Np2D*Neq - 1)/)
    DO itor = 1, ntorloc
       buffrv = 0.
       buffsd = 0.
       dd = 1 + (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq

       ! Filling the send buffer
       DO i = 1, ne2sd
          Fi = Mesh%el2sd(i)
          delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
          indp = delta + auxel
          buffsd(:, i) = sol%u_tilde(indp)
       END DO

       ! Receiving
       DO i = 1, Mesh%nghostelems
          CALL MPI_IRECV(buffrv(:, i), Neq*Np2D, MPI_REAL8, Mesh%pe2rv(i) - 1, etq, MPI_COMM_WORLD, req(i), code)
       END DO

       ! Sending
       DO i = 1, ne2sd
          CALL MPI_ISEND(buffsd(:,i),Neq*Np2D,MPI_REAL8,Mesh%pe2sd(i)-1,etq,&
               &MPI_COMM_WORLD,req(Mesh%nghostelems+i),code)
       END DO

       CALL MPI_WAITALL(SIZE(req), req, stat, code)

       ! Storing at the right place
       DO i = 1, Mesh%nghostelems
          Fi = Mesh%el2rv(i)!<--here Fi is the poloidal element
          delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
          indp = delta + auxel
          sol%u_tilde(indp) = buffrv(:, i)
       END DO

    END DO
    DEALLOCATE (buffrv, buffsd, req, stat)

    !*************************************************
    !  Communication for toroidal faces
    !*************************************************
    ALLOCATE (req(nf2sd + Mesh%nghostfaces))
    ALLOCATE (stat(MPI_STATUS_SIZE, nf2sd + Mesh%nghostfaces))
    req = mpi_request_null

    ! Allocate buffers
    ALLOCATE (buffrv(Nfl*phys%Neq, Mesh%nghostfaces))
    ALLOCATE (buffsd(Nfl*phys%Neq, nf2sd))

    auxfl = (/(j, j=0, Nfl*Neq - 1)/)
    DO itor = 1, ntorloc
       buffrv = 0.
       buffsd = 0.
       dd = 1 + (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq

       ! Filling the send buffer
       DO i = 1, nf2sd
          Fi = Mesh%fc2sd(i)
          delta = dd + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
          indt = delta + auxfl
          buffsd(:, i) = sol%u_tilde(indt)
       END DO
       ! Receiving
       DO i = 1, Mesh%nghostfaces
          CALL MPI_IRECV(buffrv(:, i), Neq*Nfl, MPI_REAL8, Mesh%pr2rv(i) - 1, etq, MPI_COMM_WORLD, req(i), code)
       END DO

       ! Sending
       DO i = 1, nf2sd
          CALL MPI_ISEND(buffsd(:, i), Neq*Nfl, MPI_REAL8, Mesh%pr2sd(i) - 1, &
               &etq, MPI_COMM_WORLD, req(Mesh%nghostfaces + i), code)
       END DO

       CALL MPI_WAITALL(SIZE(req), req, stat, code)
       ! Storing at the right place
       DO i = 1, Mesh%nghostfaces
          Fi = Mesh%fc2rv(i)
          delta = dd + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
          indt = delta + auxfl
          sol%u_tilde(indt) = buffrv(:, i)
       END DO

    END DO
    DEALLOCATE (buffrv, buffsd, req, stat)

    !************************************************************************************************
    !
    !       COMMUNICATIONS FOR TOROIDAL DISTRIBUTION OF THE MESH
    !
    !************************************************************************************************
    IF (MPIvar%ntor .GT. 1) THEN

       !******************************************
       !
       !************ Sending backward ************
       !
       !******************************************
       ! Process to send to
       IF (MPIvar%itor == 1) THEN
          prtorsd = (MPIvar%ntor - 1)*MPIvar%npol + MPIvar%ipol - 1
       ELSE
          prtorsd = (MPIvar%itor - 2)*MPIvar%npol + MPIvar%ipol - 1
       END IF

       ! Process to receive from
       IF (MPIvar%itor == MPIvar%ntor) THEN
          prtorrv = MPIvar%ipol - 1
       ELSE
          prtorrv = MPIvar%itor*MPIvar%npol + MPIvar%ipol - 1
       END IF

       !******************************************
       ! Sending backward toroidal faces
       !******************************************
       ALLOCATE (req(2*(Mesh%Nfaces - Nfdir)))
       ALLOCATE (stat(MPI_STATUS_SIZE, 2*(Mesh%Nfaces - Nfdir)))
       req = mpi_request_null
       ! Allocate buffers
       ALLOCATE (buffrv(Nfl*phys%Neq, (Mesh%Nfaces - Nfdir)))
       ALLOCATE (buffsd(Nfl*phys%Neq, (Mesh%Nfaces - Nfdir)))
       buffrv = 0.
       buffsd = 0.

       dd = 1
       auxfl = (/(j, j=0, Nfl*Neq - 1)/)
       ! Filling the send buffer
       DO Fi = 1, Mesh%Nfaces
          IF (Fi .GT. Mesh%Nintfaces) THEN
             iel = Mesh%extfaces(Fi - Mesh%Nintfaces, 1)
             ifa = Mesh%extfaces(Fi - Mesh%Nintfaces, 2)
             IF (Mesh%Fdir(iel, ifa)) CYCLE
          END IF
          delta = dd + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
          indt = delta + auxfl
          buffsd(:, Fi) = sol%u_tilde(indt)
       END DO

       ! Receiving
       DO Fi = 1, Mesh%Nfaces
          IF (Fi .GT. Mesh%Nintfaces) THEN
             iel = Mesh%extfaces(Fi - Mesh%Nintfaces, 1)
             ifa = Mesh%extfaces(Fi - Mesh%Nintfaces, 2)
             IF (Mesh%Fdir(iel, ifa)) CYCLE
          END IF
          CALL MPI_IRECV(buffrv(:, Fi), Neq*Nfl, MPI_REAL8, prtorrv, etq, MPI_COMM_WORLD, req(Fi), code)
       END DO

       ! Sending
       DO Fi = 1, Mesh%Nfaces
          IF (Fi .GT. Mesh%Nintfaces) THEN
             iel = Mesh%extfaces(Fi - Mesh%Nintfaces, 1)
             ifa = Mesh%extfaces(Fi - Mesh%Nintfaces, 2)
             IF (Mesh%Fdir(iel, ifa)) CYCLE
          END IF
          CALL MPI_ISEND(buffsd(:, Fi), Neq*Nfl, MPI_REAL8, prtorsd, &
               &etq, MPI_COMM_WORLD, req(Mesh%Nfaces - Nfdir + i), code)
       END DO

       CALL MPI_WAITALL(SIZE(req), req, stat, code)

       ! Storing at the right place
       dd = 1 + (ntorloc - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq
       DO Fi = 1, Mesh%Nfaces
          IF (Fi .GT. Mesh%Nintfaces) THEN
             iel = Mesh%extfaces(Fi - Mesh%Nintfaces, 1)
             ifa = Mesh%extfaces(Fi - Mesh%Nintfaces, 2)
             IF (Mesh%Fdir(iel, ifa)) CYCLE
          END IF
          delta = dd + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
          indt = delta + auxfl
          sol%u_tilde(indt) = buffrv(:, Fi)
       END DO

       DEALLOCATE (buffrv, buffsd, req, stat)

       !******************************************
       ! Sending backward poloidal faces
       !******************************************
       ALLOCATE (req(2*Mesh%Nelems))
       ALLOCATE (stat(MPI_STATUS_SIZE, 2*Mesh%Nelems))
       req = mpi_request_null

       ! Allocate buffers
       ALLOCATE (buffrv(Np2D*phys%Neq, Mesh%Nelems))
       ALLOCATE (buffsd(Np2D*phys%Neq, Mesh%Nelems))
       buffsd = 0.
       buffrv = 0.

       ! Filling the send buffer
       auxel =  (/(j, j=0, Np2D*Neq - 1)/)
       dd = 1 + (N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq
       DO Fi = 1, Mesh%Nelems
          delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
          indp = delta + auxel
          buffsd(:, Fi) = sol%u_tilde(indp)
       END DO

       ! Receiving
       DO Fi = 1, Mesh%Nelems
          CALL MPI_IRECV(buffrv(:, Fi), Neq*Np2D, MPI_REAL8, prtorrv, etq, MPI_COMM_WORLD, req(Fi), code)
       END DO

       ! Sending
       DO Fi = 1, Mesh%Nelems
          CALL MPI_ISEND(buffsd(:, Fi), Neq*Np2D, MPI_REAL8, prtorsd, etq, MPI_COMM_WORLD, req(Mesh%Nelems + Fi), code)
       END DO

       CALL MPI_WAITALL(SIZE(req), req, stat, code)

       ! Storing at the right place
       dd = 1 + ntorloc*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq
       DO Fi = 1, Mesh%Nelems
          delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
          indp = delta + auxel
          sol%u_tilde(indp) = buffrv(:, Fi)
       END DO

       !******************************************
       !
       !************ Sending forward ************
       !
       !******************************************

       ! Process to receive from
       IF (MPIvar%itor == 1) THEN
          prtorrv = (MPIvar%ntor - 1)*MPIvar%npol + MPIvar%ipol - 1
       ELSE
          prtorrv = (MPIvar%itor - 2)*MPIvar%npol + MPIvar%ipol - 1
       END IF

       ! Process to send to
       IF (MPIvar%itor == MPIvar%ntor) THEN
          prtorsd = MPIvar%ipol - 1
       ELSE
          prtorsd = MPIvar%itor*MPIvar%npol + MPIvar%ipol - 1
       END IF
       !******************************************
       ! Sending forward poloidal faces
       !******************************************
       buffsd = 0.
       buffrv = 0.

       ! Filling the send buffer
       dd = 1 + (ntorloc - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq
       DO Fi = 1, Mesh%Nelems
          delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
          indp = delta + auxel
          buffsd(:, Fi) = sol%u_tilde(indp)
       END DO

       ! Receiving
       DO Fi = 1, Mesh%Nelems
          CALL MPI_IRECV(buffrv(:, Fi), Neq*Np2D, MPI_REAL8, prtorrv, etq, MPI_COMM_WORLD, req(Fi), code)
       END DO

       ! Sending
       DO Fi = 1, Mesh%Nelems
          CALL MPI_ISEND(buffsd(:, Fi), Neq*Np2D, MPI_REAL8, prtorsd, etq, MPI_COMM_WORLD, req(Mesh%Nelems + Fi), code)
       END DO

       CALL MPI_WAITALL(SIZE(req), req, stat, code)

       ! Storing at the right place
       dd = 1
       DO Fi = 1, Mesh%Nelems
          delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
          indp = delta + auxel
          sol%u_tilde(indp) = buffrv(:, Fi)
       END DO

       DEALLOCATE (buffrv, buffsd, req, stat)
    END IF

  END SUBROUTINE exchangeSol
#else
  SUBROUTINE exchangeSol()
    INTEGER, PARAMETER  :: etq = 100
    INTEGER            :: Neq, Nfp, i, j, Fi
    INTEGER            :: ind(Mesh%Nnodesperface*phys%Neq)
    INTEGER            :: nf2sd
    REAL*8, ALLOCATABLE :: buffrv(:, :), buffsd(:, :)
    INTEGER            :: code
    INTEGER, ALLOCATABLE:: req(:), stat(:, :)
    INTEGER            :: perm(refElPol%Nfacenodes*phys%Neq),aux(phys%Neq*Mesh%Nnodesperface)

    Neq = phys%Neq
    Nfp = Mesh%Nnodesperface
    nf2sd = SIZE(Mesh%fc2sd)

    ! Set permutations for ghostfaces that need to be flipped
    CALL set_permutations(Neq*Nfp, Neq, perm)

    ALLOCATE (req(nf2sd + Mesh%nghostfaces))
    ALLOCATE (stat(MPI_STATUS_SIZE, nf2sd + Mesh%nghostfaces))
    req = mpi_request_null

    ! Allocate buffers
    ALLOCATE (buffrv(Mesh%Nnodesperface*phys%Neq, Mesh%nghostfaces))
    ALLOCATE (buffsd(Mesh%Nnodesperface*phys%Neq, nf2sd))
    buffrv = 0.
    buffsd = 0.
    ! Filling the send buffer
    aux =  (/(j, j=1, Neq*Nfp)/)
    DO i = 1, nf2sd
       Fi = Mesh%fc2sd(i)
       ind = (Fi - 1)*Neq*Nfp + aux
       buffsd(:, i) = sol%u_tilde(ind)
    END DO

    ! Receiving
    DO i = 1, Mesh%nghostfaces
       CALL MPI_IRECV(buffrv(:, i), Neq*Nfp, MPI_REAL8, Mesh%pr2rv(i) - 1, etq, MPI_COMM_WORLD, req(i), code)
    END DO

    ! Sending
    DO i = 1, nf2sd
       CALL MPI_ISEND(buffsd(:, i), Neq*Nfp, MPI_REAL8, Mesh%pr2sd(i) - 1, &
            &etq, MPI_COMM_WORLD, req(Mesh%nghostfaces + i), code)
    END DO

    CALL MPI_WAITALL(SIZE(req), req, stat, code)

    ! Storing at the right place
    DO i = 1, Mesh%nghostfaces
       Fi = Mesh%fc2rv(i)
       ind = (Fi - 1)*Neq*Nfp + aux
       sol%u_tilde(ind) = buffrv(:, i)
    END DO

    DEALLOCATE (buffrv, buffsd, req, stat)

  CONTAINS
    !*****************************************
    ! Set permutations for flipping faces
    !****************************************
    SUBROUTINE set_permutations(n, m, perm)
      INTEGER, INTENT(IN)  :: n, m
      INTEGER, INTENT(OUT) :: perm(:)
      INTEGER              :: i
      INTEGER              :: temp(m, n/m), templr(m, n/m)

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

  END SUBROUTINE exchangeSol
#endif

#ifdef PARALL

  SUBROUTINE gather_1D_vector_real(vector_local, vector_global, allgather)
    USE MPI_OMP
    USE GLOBALS
    REAL*8, INTENT(IN)                  :: vector_local(:)
    REAL*8, POINTER, INTENT(OUT)        :: vector_global(:)
    LOGICAL, INTENT(IN)                 :: allgather
    INTEGER                             :: recvcounts(MPIvar%glob_size), displs(MPIvar%glob_size)
    INTEGER                             :: i, ierr

    recvcounts(MPIvar%glob_id+1) = SIZE(vector_local)
    CALL MPI_Allgather(recvcounts(MPIvar%glob_id+1), 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ALLOCATE(vector_global(SUM(recvcounts)))
    vector_global = 0.

    ! vector containing the displacement
    displs(1) = 0
    DO i = 2, MPIvar%glob_size
       displs(i) = displs(i-1) + recvcounts(i-1)
    ENDDO

    IF(allgather) THEN
       CALL MPI_Allgatherv(vector_local, recvcounts(MPIvar%glob_id+1), MPI_REAL8, vector_global, recvcounts, displs, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ELSE
       CALL MPI_Gatherv(vector_local, recvcounts(MPIvar%glob_id+1), MPI_REAL8, vector_global, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    ENDIF


  ENDSUBROUTINE gather_1D_vector_real

  SUBROUTINE gather_2D_matrix_int(matrix_local, matrix_global, allgather)
    USE MPI_OMP
    USE GLOBALS
    INTEGER, INTENT(IN)                 :: matrix_local(:,:)
    INTEGER, POINTER, INTENT(OUT)       :: matrix_global(:,:)

    INTEGER, ALLOCATABLE                :: matrix_local_transpose(:,:)
    INTEGER, ALLOCATABLE                :: matrix_global_transpose(:,:)
    LOGICAL, INTENT(IN)                 :: allgather
    INTEGER                             :: recvcounts(MPIvar%glob_size), displs(MPIvar%glob_size)
    INTEGER                             :: i, ierr

    ALLOCATE(matrix_local_transpose(SIZE(matrix_local,2),SIZE(matrix_local,1)))
    matrix_local_transpose = TRANSPOSE(matrix_local)

    recvcounts(MPIvar%glob_id+1) = SIZE(matrix_local,1)*SIZE(matrix_local,2)
    CALL MPI_Allgather(recvcounts(MPIvar%glob_id+1), 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    !ALLOCATE(matrix_global(SUM(recvcounts)/SIZE(matrix_local,2),SIZE(matrix_local,2)))
    ALLOCATE(matrix_global_transpose(SIZE(matrix_local,2),SUM(recvcounts)/SIZE(matrix_local,2)))
    ALLOCATE(matrix_global(SUM(recvcounts)/SIZE(matrix_local,2),SIZE(matrix_local,2)))
    matrix_global_transpose = 0
    matrix_global = 0
    ! vector containing the displacement
    displs(1) = 0
    DO i = 2, MPIvar%glob_size
       displs(i) = displs(i-1) + recvcounts(i-1)
    ENDDO

    IF(allgather) THEN
       CALL MPI_Allgatherv(matrix_local_transpose, SIZE(matrix_local,1)*SIZE(matrix_local,2), MPI_INTEGER, matrix_global_transpose, recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    ELSE
       CALL MPI_Gatherv(matrix_local, SIZE(matrix_local,1)*SIZE(matrix_local,2), MPI_INTEGER, matrix_global_transpose, recvcounts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ENDIF

    matrix_global = TRANSPOSE(matrix_global_transpose)
    DEALLOCATE(matrix_local_transpose)
    DEALLOCATE(matrix_global_transpose)

  ENDSUBROUTINE gather_2D_matrix_int

  SUBROUTINE gather_2D_matrix_real(matrix_local, matrix_global, allgather)
    USE MPI_OMP
    USE GLOBALS
    REAL*8, INTENT(IN)                  :: matrix_local(:,:)
    REAL*8, POINTER, INTENT(OUT)        :: matrix_global(:,:)

    REAL*8, ALLOCATABLE                 :: matrix_local_transpose(:,:)
    REAL*8, ALLOCATABLE                 :: matrix_global_transpose(:,:)
    LOGICAL, INTENT(IN)                 :: allgather
    INTEGER                             :: recvcounts(MPIvar%glob_size), displs(MPIvar%glob_size)
    INTEGER                             :: i, ierr

    ALLOCATE(matrix_local_transpose(SIZE(matrix_local,2),SIZE(matrix_local,1)))
    matrix_local_transpose = TRANSPOSE(matrix_local)

    matrix_local_transpose = 0.

    recvcounts(MPIvar%glob_id+1) = SIZE(matrix_local,1)*SIZE(matrix_local,2)
    CALL MPI_Allgather(recvcounts(MPIvar%glob_id+1), 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    !ALLOCATE(matrix_global(SUM(recvcounts)/SIZE(matrix_local,2),SIZE(matrix_local,2)))
    ALLOCATE(matrix_global_transpose(SIZE(matrix_local,2),SUM(recvcounts)/SIZE(matrix_local,2)))
    ALLOCATE(matrix_global(SUM(recvcounts)/SIZE(matrix_local,2),SIZE(matrix_local,2)))
    matrix_global_transpose = 0.
    matrix_global = 0.
    ! vector containing the displacement
    displs(1) = 0
    DO i = 2, MPIvar%glob_size
       displs(i) = displs(i-1) + recvcounts(i-1)
    ENDDO

    IF(allgather) THEN
       CALL MPI_Allgatherv(matrix_local_transpose, SIZE(matrix_local,1)*SIZE(matrix_local,2), MPI_REAL8, matrix_global_transpose, recvcounts, displs, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ELSE
       CALL MPI_Gatherv(matrix_local, SIZE(matrix_local,1)*SIZE(matrix_local,2), MPI_REAL8, matrix_global_transpose, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    ENDIF

    matrix_global = TRANSPOSE(matrix_global_transpose)
    DEALLOCATE(matrix_local_transpose)
    DEALLOCATE(matrix_global_transpose)

  ENDSUBROUTINE gather_2D_matrix_real
#endif

#ifdef PARALL

  SUBROUTINE gather_mesh_solution(Mesh_in, u_in, q_in, T_glob, F_glob, flipface_glob, X_glob, u_glob, q_glob)
    TYPE(Mesh_type)                         :: Mesh_in
    REAL*8, INTENT(IN)                      :: u_in(:), q_in(:)
    INTEGER, POINTER, INTENT(OUT)           :: T_glob(:,:), F_glob(:,:), flipface_glob(:,:)!, Tb_glob(:,:)
    REAL*8, POINTER, INTENT(OUT)            :: X_glob(:,:)
    REAL*8, POINTER, INTENT(OUT)            :: u_glob(:), q_glob(:)
    REAL*8, ALLOCATABLE                     :: u_2d(:,:)
    REAL*8, ALLOCATABLE                     :: u_3d(:,:,:), temp3d(:,:,:), u_nogho3d(:,:,:)
    REAL*8, ALLOCATABLE                     :: q_4d(:,:,:,:), q_nogho4d(:,:,:,:)
    REAL*8, ALLOCATABLE                     :: u_nogho(:), q_nogho(:)
    INTEGER                                 :: i, j, ii, jj, counter, ierr


    ! reshape u_in and q_in as 3d and 4d arrays of shape [Nelems,Nnodesperelem, nphys] and [Nelems,Nnodesperelem, nphys, ndim]
    ALLOCATE(u_3d(SIZE(Mesh_in%T,1),SIZE(Mesh_in%T,2),phys%neq))
    ALLOCATE(q_4d(SIZE(Mesh_in%T,1),SIZE(Mesh_in%T,2),phys%neq, Mesh_in%Ndim))
    u_3d = 0.
    q_4D = 0.

    ! equivalent to:
    ! u_2d = TRANSPOSE(RESHAPE(u_in,[phys%neq, SIZE(u_in)/phys%neq]))
    ! temp3d = RESHAPE(u_2d, [SIZE(Mesh_in%T,2), SIZE(Mesh_in%T,1), phys%neq])
    ! permute(temp3d,u_3d)
    CALL reshape_transpose_permute(u_in, u_3d, phys%neq, SIZE(Mesh_in%T,1), SIZE(Mesh_in%T,2))
    CALL reshape_transpose_permute_4D(q_in, q_4D, Mesh_in%Ndim, phys%neq, SIZE(Mesh_in%T,1), SIZE(Mesh_in%T,2))

    !remove ghost elements from T, X, u and q
    ALLOCATE(T_glob(Mesh_in%Nel_glob, Mesh_in%Nnodesperelem))
    ALLOCATE(flipface_glob(Mesh_in%Nel_glob, refElPol%Nfaces))
    !ALLOCATE(Tb_glob(126, Mesh_in%Nnodesperface))
    ALLOCATE(X_glob(Mesh_in%Nno_glob, Mesh_in%Ndim))
    ALLOCATE(F_glob(Mesh_in%Nel_glob, refElPol%Nfaces))
    ALLOCATE(u_nogho3d(Mesh_in%Nel_glob,Mesh_in%Nnodesperelem, phys%neq))
    ALLOCATE(q_nogho4d(Mesh_in%Nel_glob,Mesh_in%Nnodesperelem, phys%neq, Mesh_in%Ndim))

    T_glob = 0
    F_glob = 0
    flipface_glob = 0
    !Tb_glob = 0
    X_glob = -1.e30
    u_nogho3d = 0.
    q_nogho4d = 0.


    DO i = 1, SIZE(Mesh_in%T,1)
      IF(Mesh_in%ghostElems(i) .NE. 1) THEN
        T_glob(Mesh_in%loc2glob_el(i),:) = Mesh_in%loc2glob_nodes(Mesh_in%T(i,:))
        flipface_glob(Mesh_in%loc2glob_el(i),:) = merge(1, 0, Mesh_in%flipface(i,:))
        F_glob(Mesh_in%loc2glob_el(i),:) = Mesh_in%loc2glob_fa(Mesh_in%F(i,:))
        X_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:)),:) = Mesh_in%X(Mesh_in%T(i,:),:)
        u_nogho3d(Mesh_in%loc2glob_el(i),:,:) = u_3d(i,:,:)
        q_nogho4d(Mesh_in%loc2glob_el(i),:,:,:) = q_4d(i,:,:,:)
      ENDIF
    ENDDO

    ! DO i = 1, SIZE(Mesh_in%Tb,1)
    !   IF(Mesh_in%ghostFaces(Mesh_in%Nfaces-Mesh_in%Nextfaces- 1 + i) .NE. 1) THEN
    !     Tb_glob(Mesh_in%loc2glob_fa(Mesh_in%Nfaces-Mesh_in%Nextfaces -1 +i),:) = Mesh_in%loc2glob_nodes(Mesh_in%Tb(i,:))
    !   ENDIF
    ! ENDDO

    ! reduce results over processes
    CALL MPI_Allreduce(MPI_IN_PLACE, T_glob, SIZE(T_glob,1)*SIZE(T_glob,2), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, F_glob, SIZE(F_glob,1)*SIZE(F_glob,2), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, flipface_glob, SIZE(flipface_glob,1)*SIZE(flipface_glob,2), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    !CALL MPI_Allreduce(MPI_IN_PLACE, Tb_glob, SIZE(Tb_glob,1)*SIZE(Tb_glob,2), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, X_glob, SIZE(X_glob,1)*SIZE(X_glob,2), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, u_nogho3d, SIZE(u_nogho3d,1)*SIZE(u_nogho3d,2)*SIZE(u_nogho3d,3), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, q_nogho4d, SIZE(q_nogho4d,1)*SIZE(q_nogho4d,2)*SIZE(q_nogho4d,3)*SIZE(q_nogho4d,4), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! reflatten back the arrays
    ALLOCATE(u_glob(SIZE(u_nogho3d,1)*SIZE(u_nogho3d,2)*SIZE(u_nogho3d,3)))
    ALLOCATE(q_glob(SIZE(q_nogho4d,1)*SIZE(q_nogho4d,2)*SIZE(q_nogho4d,3)*SIZE(q_nogho4d,4)))
    u_glob = 0.
    q_glob = 0.
    CALL flatten_row_major(u_nogho3d, u_glob, SIZE(u_nogho3d,1), SIZE(u_nogho3d,2), SIZE(u_nogho3d,3))
    CALL flatten_row_major_4D(q_nogho4d, q_glob, SIZE(q_nogho4d,1), SIZE(q_nogho4d,2), SIZE(q_nogho4d,3),SIZE(q_nogho4d,4))

    DEALLOCATE(u_3d)
    DEALLOCATE(q_4d)
    DEALLOCATE(u_nogho3d)
    DEALLOCATE(q_nogho4d)

  ENDSUBROUTINE gather_mesh_solution

  SUBROUTINE gather_mesh(Mesh_in, T_glob,  X_glob)
    TYPE(Mesh_type)                         :: Mesh_in
    INTEGER, POINTER, INTENT(OUT)           :: T_glob(:,:)
    REAL*8, POINTER, INTENT(OUT)            :: X_glob(:,:)
    INTEGER                                 :: i, j, ii, jj, counter, ierr


    !remove ghost elements from T, X, u and q
    ALLOCATE(T_glob(Mesh_in%Nel_glob, Mesh_in%Nnodesperelem))
    ALLOCATE(X_glob(Mesh_in%Nno_glob, Mesh_in%Ndim))

    T_glob = 0
    X_glob = -1.e30

    DO i = 1, SIZE(Mesh_in%T,1)
      IF(Mesh_in%ghostElems(i) .NE. 1) THEN
        T_glob(Mesh_in%loc2glob_el(i),:) = Mesh_in%loc2glob_nodes(Mesh_in%T(i,:))
        X_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:)),:) = Mesh_in%X(Mesh_in%T(i,:),:)
      ENDIF
    ENDDO

    ! reduce results over processes
    CALL MPI_Allreduce(MPI_IN_PLACE, T_glob, SIZE(T_glob,1)*SIZE(T_glob,2), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, X_glob, SIZE(X_glob,1)*SIZE(X_glob,2), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

  ENDSUBROUTINE gather_mesh

  SUBROUTINE gather_magnetic_field(Mesh_in, B_glob, B_flux_glob)
    TYPE(Mesh_type)                         :: Mesh_in
    REAL*8, POINTER, INTENT(OUT)            :: B_glob(:,:), B_flux_glob(:)
    INTEGER                                 :: i, j, ii, jj, counter, ierr


    !remove ghost elements from T, X, u and q
    ALLOCATE(B_glob(Mesh_in%Nno_glob, 3))
    ALLOCATE(B_flux_glob(Mesh_in%Nno_glob))

    B_glob = -1.e30
    B_flux_glob = -1.e30

    DO i = 1, SIZE(Mesh_in%T,1)
      IF(Mesh_in%ghostElems(i) .NE. 1) THEN
        B_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:)),:) = phys%B(Mesh_in%T(i,:),:)
        B_flux_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:))) = phys%magnetic_flux(Mesh_in%T(i,:))
      ENDIF
    ENDDO

    ! reduce results over processes
    CALL MPI_Allreduce(MPI_IN_PLACE, B_glob, SIZE(B_glob,1)*SIZE(B_glob,2), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, B_flux_glob, SIZE(B_flux_glob,1), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
  ENDSUBROUTINE gather_magnetic_field



  SUBROUTINE reshape_transpose_permute(u_in, u_3d, neq, Nelems, Nnodesperelem)
    IMPLICIT NONE
    REAL, INTENT(IN)    :: u_in(:)
    REAL, INTENT(OUT)   :: u_3d(Nelems, Nnodesperelem, neq)
    INTEGER, INTENT(IN) :: neq, Nelems, Nnodesperelem
    INTEGER             :: i, j, k

    ! Iterate through the elements to directly reshape, transpose, and permute
    DO k = 1, neq
        DO i = 1, Nelems
            DO j = 1, Nnodesperelem
                ! Calculate the index in the flattened input array u_in
                ! Reshape, transpose, and permute in one go:
                u_3d(i, j, k) = u_in(k + (j - 1) * neq + (i - 1) * neq * Nnodesperelem)
            ENDDO
        ENDDO
    ENDDO
  END SUBROUTINE reshape_transpose_permute

  SUBROUTINE reshape_transpose_permute_4D(q_in, q_4d, ndim, neq, Nelems, Nnodesperelem)
    IMPLICIT NONE
    REAL, INTENT(IN)    :: q_in(:)
    REAL, INTENT(OUT)   :: q_4d(Nelems, Nnodesperelem, neq, ndim)
    INTEGER, INTENT(IN) :: ndim, neq, Nelems, Nnodesperelem
    INTEGER             :: i, j, k, l

    ! Iterate through the elements to directly reshape, transpose, and permute
    DO l = 1, ndim
      DO k = 1, neq
          DO i = 1, Nelems
              DO j = 1, Nnodesperelem
                  ! Calculate the index in the flattened input array u_in
                  ! Reshape, transpose, and permute in one go:
                  q_4d(i, j, k, l) = q_in(l + (k-1)*ndim + (j - 1) * neq * ndim + (i - 1) * neq * Nnodesperelem * ndim)
              ENDDO
          ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE reshape_transpose_permute_4D

  SUBROUTINE permute(input, output)
    IMPLICIT NONE
    REAL, INTENT(IN)  :: input(:,:,:)
    REAL, INTENT(OUT) :: output(SIZE(input,2), SIZE(input,1), SIZE(input,3))
    INTEGER           :: i, j, k

    ! Swap the first two dimensions of the input array into the output array
    DO k = 1, SIZE(input, 3)
       DO j = 1, SIZE(input, 2)
          DO i = 1, SIZE(input, 1)
             output(j, i, k) = input(i, j, k)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE permute

  SUBROUTINE flatten_row_major(input, output, dim1, dim2, dim3)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim1, dim2, dim3
    REAL, INTENT(IN)    :: input(dim1, dim2, dim3)
    REAL, INTENT(OUT)   :: output(dim1 * dim2 * dim3)
    INTEGER             :: i, j, k, index

    index = 1  ! Initialize the index for the 1D output array

    ! Iterate through the 3D array in row-major order
    DO k = 1, dim1
        DO j = 1, dim2
            DO i = 1, dim3
                output(index) = input(k, j, i)
                index = index + 1
            ENDDO
        ENDDO
    ENDDO
  END SUBROUTINE flatten_row_major

  SUBROUTINE flatten_row_major_4D(input, output, dim1, dim2, dim3, dim4)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim1, dim2, dim3, dim4
    REAL, INTENT(IN)    :: input(dim1, dim2, dim3, dim4)
    REAL, INTENT(OUT)   :: output(dim1 * dim2 * dim3 * dim4)
    INTEGER             :: i, j, k, l, index

    index = 1  ! Initialize the index for the 1D output array

    ! Iterate through the 3D array in row-major order
    DO l = 1, dim1
      DO k = 1, dim2
          DO j = 1, dim3
              DO i = 1, dim4
                  output(index) = input(l, k, j, i)
                  index = index + 1
              ENDDO
          ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE flatten_row_major_4D

  SUBROUTINE reshape_permute_flatten(u_in, u_glob, neq, dim1, dim2)
    IMPLICIT NONE
    REAL, INTENT(IN) :: u_in(:)
    REAL, INTENT(OUT) :: u_glob(:)
    INTEGER, INTENT(IN) :: neq, dim1, dim2
    INTEGER :: n, i, j, k
    INTEGER :: num_elements

    ! Calculate the total number of elements in u_in and u_glob
    num_elements = SIZE(u_in)

    ! Check if the sizes match
    IF (SIZE(u_glob) /= num_elements) THEN
        PRINT *, "Error: Mismatch in size of input and output arrays."
        RETURN
    END IF

    ! Combined loop for reshaping, transposing, permuting, and flattening
    DO n = 1, num_elements
        ! Calculate the original indices in terms of the 3D array
        k = (n - 1) / (dim1 * dim2) + 1             ! Third dimension (phys%neq)
        j = MOD((n - 1) / dim1, dim2) + 1            ! Second dimension (size of Mesh_in%T,2)
        i = MOD(n - 1, dim1) + 1                    ! First dimension (size of Mesh_in%T,1)

        ! Map directly from the input 1D array to the output 1D array with the new ordering
        u_glob(n) = u_in((k - 1) * (dim1 * dim2) + (j - 1) * dim1 + i)
    ENDDO
END SUBROUTINE reshape_permute_flatten

#endif

END MODULE Communications
