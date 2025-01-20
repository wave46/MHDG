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

#ifdef PARALL
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
#endif

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
#ifdef PARALL
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

  SUBROUTINE gather_1D_vector_int(vector_local, vector_global, allgather)
    USE MPI_OMP
    USE GLOBALS
    INTEGER, INTENT(IN)                  :: vector_local(:)
    INTEGER, POINTER, INTENT(OUT)        :: vector_global(:)
    LOGICAL, INTENT(IN)                  :: allgather
    INTEGER                              :: recvcounts(MPIvar%glob_size), displs(MPIvar%glob_size)
    INTEGER                              :: i, ierr

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
       CALL MPI_Allgatherv(vector_local, recvcounts(MPIvar%glob_id+1), MPI_INT, vector_global, recvcounts, displs, MPI_INT, MPI_COMM_WORLD, ierr)
    ELSE
       CALL MPI_Gatherv(vector_local, recvcounts(MPIvar%glob_id+1), MPI_INT, vector_global, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    ENDIF


  ENDSUBROUTINE gather_1D_vector_int

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
       CALL MPI_Gatherv(matrix_local_transpose, SIZE(matrix_local,1)*SIZE(matrix_local,2), MPI_INTEGER, matrix_global_transpose, recvcounts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
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

  SUBROUTINE gather_mesh(Mesh_in, T_glob, X_glob, Tb_glob, F_glob, N_glob, intfaces_glob, extfaces_glob, boundaryFlag_glob, Tlin_glob, periodic_faces_glob, elemSize_glob, flag_elems_sc_glob, scdiff_nodes_glob)
    USE preprocess, only: createNodalConnectivity, computeElementSize
    TYPE(Mesh_type)                         :: Mesh_in
    INTEGER, POINTER, INTENT(OUT)           :: T_glob(:,:)
    REAL*8, POINTER, INTENT(OUT)            :: X_glob(:,:)
    REAL*8, OPTIONAL, POINTER, INTENT(OUT)  :: elemSize_glob(:), scdiff_nodes_glob(:,:)
    INTEGER, OPTIONAL, POINTER, INTENT(OUT) :: Tb_glob(:,:), F_glob(:,:), N_glob(:,:), intfaces_glob(:,:), extfaces_glob(:,:), Tlin_glob(:,:)
    INTEGER, OPTIONAL, POINTER, INTENT(OUT) :: boundaryFlag_glob(:), periodic_faces_glob(:), flag_elems_sc_glob(:)
    INTEGER, ALLOCATABLE                    :: intfaces(:,:), extfaces(:,:), Tb(:,:), boundaryFlag(:), periodic_faces(:)
    INTEGER                                 :: i, counter, index, ierr

    !remove ghost elements from T, X
    ALLOCATE(T_glob(Mesh_in%Nel_glob, Mesh_in%Nnodesperelem))
    ALLOCATE(X_glob(Mesh_in%Nno_glob, Mesh_in%Ndim))

    T_glob = 0
    X_glob = -1.e30

    IF(PRESENT(F_glob)) THEN
       ALLOCATE(F_glob(Mesh_in%Nel_glob, refElPol%Nfaces))
       F_glob = 0
    ENDIF
    IF(PRESENT(intfaces_glob)) THEN
       ALLOCATE(intfaces(Mesh_in%Nintfaces_nogho, 5))
       intfaces = 0
    ENDIF
    IF(PRESENT(extfaces_glob)) THEN
       ALLOCATE(extfaces(Mesh_in%Nextfaces_nogho, 2))
       extfaces = 0
    ENDIF
    IF(PRESENT(Tlin_glob)) THEN
       ALLOCATE(Tlin_glob(Mesh_in%Nel_glob, refElPol%Nvertices))
       Tlin_glob = 0
    ENDIF
    IF(PRESENT(flag_elems_sc_glob)) THEN
      ALLOCATE(flag_elems_sc_glob(Mesh_in%Nel_glob))
      flag_elems_sc_glob = 0
    ENDIF
    IF(PRESENT(scdiff_nodes_glob)) THEN
      ALLOCATE(scdiff_nodes_glob(Mesh_in%Nel_glob, Mesh_in%Nnodesperelem))
      scdiff_nodes_glob = 0
    ENDIF

    DO i = 1, Mesh_in%Nelems
       IF(Mesh_in%ghostElems(i) .EQ. 0) THEN
          X_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:)),:) = Mesh_in%X(Mesh_in%T(i,:),:)
          index = Mesh_in%loc2glob_el(i)
          T_glob(index,:) = Mesh_in%loc2glob_nodes(Mesh_in%T(i,:))
          IF(PRESENT(F_glob)) F_glob(index,:) = Mesh_in%loc2glob_fa(Mesh_in%F(i,:))
          IF(PRESENT(Tlin_glob)) Tlin_glob(index,:) = Mesh_in%loc2glob_nodes(Mesh_in%Tlin(i,:))
          IF(PRESENT(flag_elems_sc_glob)) flag_elems_sc_glob(index) = Mesh_in%flag_elems_sc(i)
          IF(PRESENT(scdiff_nodes_glob)) scdiff_nodes_glob(index,:) = Mesh_in%scdiff_nodes(i,:)
       ENDIF
    ENDDO

    IF(PRESENT(Tb_glob)) THEN
       ALLOCATE(Tb(Mesh_in%Nextfaces_nogho,Mesh_in%Nnodesperface))
       counter = 1
       DO i = 1, Mesh_in%NextFaces
          IF(Mesh_in%ghostFaces(Mesh_in%Nintfaces+i) .EQ. 0) THEN
             Tb(counter,:) = Mesh_in%loc2glob_nodes(Mesh_in%Tb(i,:))
             counter = counter + 1
          ENDIF
       ENDDO
    ENDIF

    IF(PRESENT(intfaces_glob)) THEN
      counter = 1
       DO i = 1, Mesh_in%Nintfaces
          IF(Mesh_in%ghostFaces(i) .EQ. 0) THEN
             intfaces(counter,1)   = Mesh_in%loc2glob_el(Mesh_in%intfaces(i,1)) ! number of the triangle
             intfaces(counter,2)   = Mesh_in%intfaces(i,2) ! number of the face
             intfaces(counter,3)   = Mesh_in%loc2glob_el(Mesh_in%intfaces(i,3)) ! number of the neighbour triangle
             intfaces(counter,4:5) = Mesh_in%intfaces(i,4:5) ! number of the face of the neighbour triangle and number of the node sharing the first knot
             counter = counter + 1
          ENDIF
       ENDDO
    ENDIF

    IF(PRESENT(extfaces_glob)) THEN
      counter = 1
      DO i = 1, Mesh_in%NextFaces
          IF (Mesh_in%ghostFaces(Mesh_in%Nintfaces+i) .EQ. 0) THEN
             extfaces(counter,1) = Mesh_in%loc2glob_el(Mesh_in%extfaces(i,1)) ! number of the triangle
             extfaces(counter,2) = Mesh_in%extfaces(i,2) ! number of the face
             counter = counter + 1
          ENDIF
       ENDDO
    ENDIF


    IF(PRESENT(boundaryFlag_glob)) THEN
       ALLOCATE(boundaryFlag(Mesh_in%Nextfaces_nogho))
       boundaryFlag = 0

       counter = 1
       DO i = 1, Mesh_in%Nextfaces
          IF((Mesh_in%boundaryFlag(i) .NE. 0) .AND. (Mesh_in%ghostFaces(Mesh_in%Nintfaces+i) .EQ. 0)) THEN
             boundaryFlag(counter) = Mesh_in%boundaryFlag(i)
             counter = counter + 1
          ENDIF
       ENDDO
    ENDIF

    IF(PRESENT(periodic_faces_glob)) THEN
       ALLOCATE(periodic_faces(Mesh_in%Nextfaces_nogho))
       periodic_faces = 0

       counter = 1
       DO i = 1, Mesh_in%Nextfaces
          IF((Mesh_in%boundaryFlag(i) .NE. 0) .AND. (Mesh_in%ghostFaces(Mesh_in%Nintfaces+i) .EQ. 0)) THEN
             periodic_faces(counter) = Mesh_in%periodic_faces(i)
             counter = counter + 1
          ENDIF
       ENDDO
    ENDIF

    ! stack arrays
    IF(PRESENT(intfaces_glob))        CALL gather_2D_matrix_int(intfaces, intfaces_glob, .FALSE.)
    IF(PRESENT(extfaces_glob))        CALL gather_2D_matrix_int(extfaces, extfaces_glob, .FALSE.)
    IF(PRESENT(Tb_glob))              CALL gather_2D_matrix_int(Tb, Tb_glob, .FALSE.)
    IF(PRESENT(boundaryFlag_glob))    CALL gather_1D_vector_int(boundaryFlag, boundaryFlag_glob, .FALSE.)
    IF(PRESENT(periodic_faces_glob))  CALL gather_1D_vector_int(periodic_faces, periodic_faces_glob, .FALSE.)

    ! reduce results over processes
    CALL MPI_Allreduce(MPI_IN_PLACE, T_glob, SIZE(T_glob,1)*SIZE(T_glob,2), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, X_glob, SIZE(X_glob,1)*SIZE(X_glob,2), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! I believe it is much easier to reconstruct the nodal connectivity rather than trying to reassemble it from processes (lots of communication)
    IF(PRESENT(N_glob)) THEN
      CALL createNodalConnectivity(T_glob, Mesh%Nel_glob, Mesh%Nno_glob, Mesh%Nnodesperelem, N_glob)
    ENDIF

    IF(PRESENT(F_glob))             CALL MPI_Allreduce(MPI_IN_PLACE, F_glob, SIZE(F_glob,1)*SIZE(F_glob,2), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF(PRESENT(Tlin_glob))          CALL MPI_Allreduce(MPI_IN_PLACE, Tlin_glob, SIZE(Tlin_glob,1)*SIZE(Tlin_glob,2), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF(PRESENT(flag_elems_sc_glob))    CALL MPI_Allreduce(MPI_IN_PLACE, flag_elems_sc_glob, SIZE(flag_elems_sc_glob), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF(PRESENT(scdiff_nodes_glob))  CALL MPI_Allreduce(MPI_IN_PLACE, scdiff_nodes_glob, SIZE(scdiff_nodes_glob,1)*SIZE(scdiff_nodes_glob,2), MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! elemsize cannot be gathered as there is no information on the order of elements in the global mesh
    IF(PRESENT(elemSize_glob)) THEN
      CALL computeElementSize(T_glob, X_glob, Tlin_glob, elemSize_glob)
    ENDIF

    IF(PRESENT(intfaces_glob)) DEALLOCATE(intfaces)
    IF(PRESENT(extfaces_glob)) DEALLOCATE(extfaces)
    IF(PRESENT(Tb_glob)) DEALLOCATE(Tb)
    IF(PRESENT(boundaryFlag_glob)) DEALLOCATE(boundaryFlag)
    IF(PRESENT(periodic_faces_glob)) DEALLOCATE(periodic_faces)
  ENDSUBROUTINE gather_mesh

  SUBROUTINE gather_solution(Mesh_in, Nnodesperelem, Nnodesperface, u_tilde_in, u_in, q_in, u_glob, u_tilde_glob, q_glob)
    TYPE(Mesh_type)                         :: Mesh_in
    INTEGER, INTENT(IN)                     :: Nnodesperelem
    INTEGER, INTENT(IN), OPTIONAL           :: Nnodesperface
    REAL*8, INTENT(IN)                      :: u_in(:)
    REAL*8, INTENT(IN), OPTIONAL            :: u_tilde_in(:)
    REAL*8, INTENT(IN), OPTIONAL            :: q_in(:)
    REAL*8, POINTER, INTENT(OUT)            :: u_glob(:)
    REAL*8, POINTER, INTENT(OUT), OPTIONAL  :: q_glob(:)
    REAL*8, POINTER, INTENT(OUT), OPTIONAL  :: u_tilde_glob(:)
    REAL*8, ALLOCATABLE                     :: u_3d(:,:,:), u_nogho_3d(:,:,:), u_tilde_3d(:,:,:), u_tilde_nogho_3d(:,:,:)
    REAL*8, ALLOCATABLE                     :: q_4d(:,:,:,:), q_nogho_4d(:,:,:,:)
    INTEGER                                 :: i, ierr

    ! reshape u_in and q_in as 3d and 4d arrays of shape [Nelems,Nnodesperelem, nphys] and [Nelems,Nnodesperelem, nphys, ndim]
    ALLOCATE(u_3d(Mesh_in%Nelems, Nnodesperelem,phys%neq))
    u_3d = 0.

    IF(PRESENT(u_tilde_in)) THEN
      ALLOCATE(u_tilde_3d(Mesh_in%Nfaces, Nnodesperface, phys%neq))
      u_tilde_3d = 0.
    ENDIF
    IF(PRESENT(q_in)) THEN
       ALLOCATE(q_4d(Mesh_in%Nelems, Nnodesperelem,phys%neq, Mesh_in%Ndim))
       q_4D = 0.
    ENDIF

    ! equivalent to:
    ! u_2d = TRANSPOSE(RESHAPE(u_in,[phys%neq, SIZE(u_in)/phys%neq]))
    ! temp3d = RESHAPE(u_2d, [Nnodesperelem, Mesh_in%Nelems, phys%neq])
    ! permute(temp3d,u_3d)
                            CALL reshape_transpose_permute(u_in, u_3d, phys%neq, Mesh_in%Nelems, Nnodesperelem)
    IF(PRESENT(u_tilde_in)) CALL reshape_transpose_permute(u_tilde_in, u_tilde_3d, phys%neq, Mesh_in%Nfaces, Nnodesperface)
    IF(PRESENT(q_in))       CALL reshape_transpose_permute_4D(q_in, q_4D, Mesh_in%Ndim, phys%neq, Mesh_in%Nelems, Nnodesperelem)

    ! allocation
    ALLOCATE(u_nogho_3d(Mesh_in%Nel_glob, Nnodesperelem, phys%neq))
    u_nogho_3d = 0.
    IF(PRESENT(u_tilde_in)) THEN
      ALLOCATE(u_tilde_nogho_3d(Mesh_in%Nfa_glob, Nnodesperface, phys%neq))
      u_tilde_nogho_3d = 0.
    ENDIF
    IF(PRESENT(q_in)) THEN
      ALLOCATE(q_nogho_4d(Mesh_in%Nel_glob, Nnodesperelem, phys%neq, Mesh_in%Ndim))
      q_nogho_4d = 0.
    ENDIF

    ! filtering out ghost cells
    DO i = 1, Mesh_in%Nelems
      IF(Mesh_in%ghostElems(i) .EQ. 0) THEN
        u_nogho_3d(Mesh_in%loc2glob_el(i),:,:) = u_3d(i,:,:)
      ENDIF
    ENDDO

    IF(PRESENT(q_in)) THEN
      ! filtering out ghost cells
      DO i = 1, Mesh_in%Nelems
        IF(Mesh_in%ghostElems(i) .EQ. 0) THEN
          q_nogho_4d(Mesh_in%loc2glob_el(i),:,:,:) = q_4d(i,:,:,:)
        ENDIF
      ENDDO
    ENDIF

    ! filtering out ghost faces
    IF(PRESENT(u_tilde_in)) THEN
      DO i = 1, Mesh_in%Nfaces
        IF(Mesh_in%ghostFaces(i) .EQ. 0) THEN
          u_tilde_nogho_3d(Mesh_in%loc2glob_fa(i),:,:) = u_tilde_3d(i,:,:)
        ENDIF
      ENDDO
    ENDIF

    ! reduce results over processes
                            CALL MPI_Allreduce(MPI_IN_PLACE, u_nogho_3d, SIZE(u_nogho_3d,1)*SIZE(u_nogho_3d,2)*SIZE(u_nogho_3d,3), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF(PRESENT(u_tilde_in)) CALL MPI_Allreduce(MPI_IN_PLACE, u_tilde_nogho_3d, SIZE(u_tilde_nogho_3d,1)*SIZE(u_tilde_nogho_3d,2)*SIZE(u_tilde_nogho_3d,3), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF(PRESENT(q_in))       CALL MPI_Allreduce(MPI_IN_PLACE, q_nogho_4d, SIZE(q_nogho_4d,1)*SIZE(q_nogho_4d,2)*SIZE(q_nogho_4d,3)*SIZE(q_nogho_4d,4), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! reflatten back the arrays
    ALLOCATE(u_glob(SIZE(u_nogho_3d,1)*SIZE(u_nogho_3d,2)*SIZE(u_nogho_3d,3)))
    u_glob = 0.

    IF(PRESENT(u_tilde_in)) THEN
       ALLOCATE(u_tilde_glob(SIZE(u_tilde_nogho_3d,1)*SIZE(u_tilde_nogho_3d,2)*SIZE(u_tilde_nogho_3d,3)))
       u_tilde_glob = 0.
    ENDIF

    IF(PRESENT(q_in)) THEN
       ALLOCATE(q_glob(SIZE(q_nogho_4d,1)*SIZE(q_nogho_4d,2)*SIZE(q_nogho_4d,3)*SIZE(q_nogho_4d,4)))
       q_glob = 0.
    ENDIF

                            CALL flatten_row_major(u_nogho_3d, u_glob, SIZE(u_nogho_3d,1), SIZE(u_nogho_3d,2), SIZE(u_nogho_3d,3))
    IF(PRESENT(u_tilde_in)) CALL flatten_row_major(u_tilde_nogho_3d, u_tilde_glob, SIZE(u_tilde_nogho_3d,1), SIZE(u_tilde_nogho_3d,2), SIZE(u_tilde_nogho_3d,3))
    IF(PRESENT(q_in))       CALL flatten_row_major_4D(q_nogho_4d, q_glob, SIZE(q_nogho_4d,1), SIZE(q_nogho_4d,2), SIZE(q_nogho_4d,3),SIZE(q_nogho_4d,4))

                            DEALLOCATE(u_3d, u_nogho_3d)
    IF(PRESENT(q_in))       DEALLOCATE(q_4d, q_nogho_4d)
    IF(PRESENT(u_tilde_in)) DEALLOCATE(u_tilde_3d, u_tilde_nogho_3d)

  ENDSUBROUTINE gather_solution

  SUBROUTINE gather_magnetic_field(Mesh_in, B_glob, magnetic_flux_glob, magnetic_psi_glob, Bperturb_glob, Jtor_glob)

    TYPE(Mesh_type)                         :: Mesh_in
    REAL*8, POINTER, INTENT(OUT)            :: B_glob(:,:), magnetic_flux_glob(:)
    REAL*8, OPTIONAL, POINTER, INTENT(OUT)  :: magnetic_psi_glob(:), Bperturb_glob(:,:), Jtor_glob(:)
    INTEGER                                 :: i, ierr


    !remove ghost elements from T, X, u and q
    ALLOCATE(B_glob(Mesh_in%Nno_glob, 3))
    ALLOCATE(magnetic_flux_glob(Mesh_in%Nno_glob))
    B_glob = -1.e30
    magnetic_flux_glob = -1.e30

    IF(PRESENT(magnetic_psi_glob)) THEN
       ALLOCATE(magnetic_psi_glob(Mesh_in%Nno_glob))
       magnetic_psi_glob = -1.e30
    ENDIF
    IF(PRESENT(Jtor_glob)) THEN
       ALLOCATE(Jtor_glob(Mesh_in%Nno_glob))
       Jtor_glob = -1.e30
    ENDIF
    IF(PRESENT(Bperturb_glob)) THEN
       ALLOCATE(Bperturb_glob(Mesh_in%Nno_glob,3))
       Bperturb_glob = -1.e30
    ENDIF

    DO i = 1, SIZE(Mesh_in%T,1)
       IF(Mesh_in%ghostElems(i) .NE. 1) THEN
          B_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:)),:) = phys%B(Mesh_in%T(i,:),:)
          magnetic_flux_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:))) = phys%magnetic_flux(Mesh_in%T(i,:))
          IF(PRESENT(magnetic_psi_glob)) magnetic_psi_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:))) = phys%magnetic_psi(Mesh_in%T(i,:))
          IF(PRESENT(Bperturb_glob))     Bperturb_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:)),:) = phys%Bperturb(Mesh_in%T(i,:),:)
          IF(PRESENT(Jtor_glob))         Jtor_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:))) = phys%Jtor(Mesh_in%T(i,:))
       ENDIF
    ENDDO

    ! reduce results over processes
    CALL MPI_Allreduce(MPI_IN_PLACE, B_glob, SIZE(B_glob,1)*SIZE(B_glob,2), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, magnetic_flux_glob, SIZE(magnetic_flux_glob,1), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF(PRESENT(magnetic_psi_glob)) CALL MPI_Allreduce(MPI_IN_PLACE, magnetic_psi_glob, SIZE(magnetic_psi_glob,1), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF(PRESENT(Bperturb_glob)) CALL MPI_Allreduce(MPI_IN_PLACE, Bperturb_glob, SIZE(Bperturb_glob,1), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF(PRESENT(Jtor_glob)) CALL MPI_Allreduce(MPI_IN_PLACE, Jtor_glob, SIZE(Jtor_glob,1), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

  ENDSUBROUTINE gather_magnetic_field

  SUBROUTINE gather_additional(Mesh_in, scdiff_nodes_glob)

    TYPE(Mesh_type)                         :: Mesh_in
    REAL*8, POINTER, INTENT(OUT)            :: scdiff_nodes_glob(:,:)
    INTEGER                                 :: i, ierr


    !remove ghost elements from T, X, u and q
    ALLOCATE(scdiff_nodes_glob(Mesh_in%Nno_glob, 3))

    scdiff_nodes_glob = -1.e30

    DO i = 1, SIZE(Mesh_in%T,1)
       IF(Mesh_in%ghostElems(i) .NE. 1) THEN
          scdiff_nodes_glob(Mesh_in%loc2glob_nodes(Mesh_in%T(i,:)),:) = phys%B(Mesh_in%T(i,:),:)
       ENDIF
    ENDDO

    ! reduce results over processes
    CALL MPI_Allreduce(MPI_IN_PLACE, scdiff_nodes_glob, SIZE(scdiff_nodes_glob,1)*SIZE(scdiff_nodes_glob,2), MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

  ENDSUBROUTINE gather_additional

  SUBROUTINE reshape_transpose_permute(u_in, u_3d, neq, Nelems, Nnodesperelem)
    IMPLICIT NONE
    REAL*8, INTENT(IN)  :: u_in(:)
    REAL*8, INTENT(OUT) :: u_3d(Nelems, Nnodesperelem, neq)
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
    REAL*8, INTENT(IN)    :: q_in(:)
    REAL*8, INTENT(OUT)   :: q_4d(Nelems, Nnodesperelem, neq, ndim)
    INTEGER, INTENT(IN)   :: ndim, neq, Nelems, Nnodesperelem
    INTEGER               :: i, j, k, l

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
    REAL*8, INTENT(IN)  :: input(:,:,:)
    REAL*8, INTENT(OUT) :: output(SIZE(input,2), SIZE(input,1), SIZE(input,3))
    INTEGER             :: i, j, k

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
    REAL*8, INTENT(IN)  :: input(dim1, dim2, dim3)
    REAL*8, INTENT(OUT) :: output(dim1 * dim2 * dim3)
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
    REAL*8, INTENT(IN)    :: input(dim1, dim2, dim3, dim4)
    REAL*8, INTENT(OUT)   :: output(dim1 * dim2 * dim3 * dim4)
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

  SUBROUTINE reshape_permute_flatten(u_in, u_glob, dim1, dim2)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: u_in(:)
    REAL*8, INTENT(OUT) :: u_glob(:)
    INTEGER, INTENT(IN) :: dim1, dim2
    INTEGER :: n, i, j, k, dim3
    INTEGER :: Nelems

    ! Calculate the total number of elements in u_in and u_glob
    Nelems = SIZE(u_in)
    dim3 = dim1*dim2

    ! Check if the sizes match
    IF (SIZE(u_glob) .NE. Nelems) THEN
       PRINT *, "Error: Mismatch in size of input and output arrays."
       RETURN
    END IF

    ! Combined loop for reshaping, transposing, permuting, and flattening
    DO n = 1, Nelems
       ! Calculate the original indices in terms of the 3D array
       k = (n - 1) / dim3 + 1             ! Third dimension (phys%neq)
       j = MOD((n - 1) / dim1, dim2) + 1            ! Second dimension (size of Mesh_in%T,2)
       i = MOD(n - 1, dim1) + 1                    ! First dimension (size of Mesh_in%T,1)

       ! Map directly from the input 1D array to the output 1D array with the new ordering
       u_glob(n) = u_in((k - 1) * (dim1 * dim2) + (j - 1) * dim1 + i)
    ENDDO
  END SUBROUTINE reshape_permute_flatten

#endif

END MODULE Communications
