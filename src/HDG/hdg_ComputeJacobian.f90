!*****************************************
! project: MHDG
! file: hdg_ConvectionMatrices.f90
! date: 20/12/2016
! Generate the matrices that change
! during the NR iterative process
!*****************************************

SUBROUTINE HDG_computeJacobian()
  USE globals
  USE LinearAlgebra
  USE analytical
  USE physics
  USE printUtils
  USE MPI_OMP
  USE Debug
  USE hdg_limitingtechniques, ONLY:HDG_ShockCapturing

  IMPLICIT NONE

  !***********************************************************************
  !
  !              COMPUTATION OF THE JACOBIANS
  !
  !***********************************************************************
  INTEGER*4             :: Ndim,Neq,N2D,Npel,Npfl,Nfre,Ng1d,Ng2d
  INTEGER*4             :: iel,ifa,iface,i,j
  INTEGER*4             :: sizeu,sizel
  REAL*8,ALLOCATABLE    :: ures(:,:),lres(:,:),u0res(:,:,:)
  REAL*8,ALLOCATABLE    :: Xel(:,:),Xfl(:,:)
  REAL*8,ALLOCATABLE    :: tau_save(:,:)
  REAL*8,ALLOCATABLE    :: xy_g_save(:,:)
  LOGICAL               :: isdir
#ifdef TEMPERATURE
  REAL*8                :: coefi,coefe
#endif
#ifdef TOR3D
  ! Definitions in 3D
  INTEGER*4             :: itor,iel3,ntorloc,itorg,Np1Dpol,Np1Dtor,Np2D,Npfp,Ng1dpol,Ng1dtor,Ngfp,Ngfl,Nfdir,Ngvo
  INTEGER*4             :: inde(refElTor%Nnodes3D)
  INTEGER*4             :: indft(refElTor%Nfl),indfp(Mesh%Nnodesperelem),indl(Mesh%Nnodesperelem)
  INTEGER*4             :: ind_loc(refElPol%Nfaces,refElTor%Nfl*phys%Neq),perm(refElTor%Nfl*phys%Neq)
  INTEGER               :: dd
  INTEGER               :: ind_dim(refElPol%Nfaces + 2),ind_sta(refElPol%Nfaces + 2),aux
  REAL*8                :: tdiv(numer%ntor + 1),tel(refElTor%Nnodes1d),tg(1),htor
  REAL*8,ALLOCATABLE    :: ue(:,:),u0e(:,:,:)
  REAL*8,ALLOCATABLE    :: ufp(:,:),uft(:,:)
  REAL*8,ALLOCATABLE    :: uefp(:,:),ueft(:,:)
  REAL*8,ALLOCATABLE    :: qe(:,:)
  REAL*8,ALLOCATABLE    :: qefp(:,:),qeft(:,:)
  REAL*8,ALLOCATABLE    :: qres(:,:)
  REAL*8,ALLOCATABLE    :: Bel(:,:),fluxel(:),Bfl(:,:),Bfp(:,:)
  INTEGER               :: indbe(refElTor%Nnodes3d),indbp(Mesh%Nnodesperelem),indbt(refElTor%Nfl)
  REAL*8                :: Jtorel(refElTor%Nnodes3d)
#else
  ! Definitions in 2D
  LOGICAL               :: save_tau
  INTEGER*4             :: inde(Mesh%Nnodesperelem)
  INTEGER*4             :: indf(refElPol%Nfacenodes)
  INTEGER*4             :: ind_loc(refElPol%Nfaces,refElPol%Nfacenodes*phys%Neq),perm(refElPol%Nfacenodes*phys%Neq)
  REAL*8                :: ue(Mesh%Nnodesperelem,phys%Neq),u0e(Mesh%Nnodesperelem,phys%Neq,time%tis)
  REAL*8                :: uf(refElPol%Nfacenodes,phys%Neq),uef(refElPol%Nfacenodes,phys%Neq)
  INTEGER               :: indtausave(refElPol%Nfaces*refElPol%Ngauss1d)
  INTEGER               :: inddiff_nn_Vol(refElPol%NGauss2D)
  REAL*8                :: tau_save_el(refElPol%Nfaces*refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Nfaces*refElPol%Ngauss1d,2)
  REAL*8                :: qe(Mesh%Nnodesperelem,phys%Neq*2),qef(refElPol%Nfacenodes,phys%Neq*2)
  REAL*8,ALLOCATABLE    :: qres(:,:)
  REAL*8                :: Bel(refElPol%Nnodes2d,3),fluxel(refElPol%Nnodes2d),psiel(refElPol%Nnodes2d),Bfl(refElPol%Nfacenodes,3),psifl(refElPol%Nfacenodes)
#ifdef KEQUATION
  real*8                :: omegael(refElPol%Nnodes2d),q_cylel(refElPol%Nnodes2d),q_cylfl(refElPol%Nfacenodes)
#endif
  REAL*8                :: Jtorel(refElPol%Nnodes2d)
  REAL*8                :: n,El_n,nn,El_nn,totaln
  REAL*8                :: diff_nn_Vol_el(refElPol%NGauss2D),v_nn_Vol_el(refElPol%NGauss2D,Mesh%Ndim),Xg_el(refElPol%NGauss2D,Mesh%Ndim)
  REAL*8                :: diff_nn_Fac_el(refElPol%Nfaces*refElPol%NGauss1D),v_nn_Fac_el(refElPol%Nfaces*refElPol%NGauss1D,Mesh%Ndim)
#endif
#ifdef PARALL
  INTEGER               :: ierr
#endif

  IF (utils%printint > 1) THEN
    WRITE (6,*) '*************************************************'
    WRITE (6,*) '*          COMPUTING JACOBIAN                   *'
    WRITE (6,*) '*************************************************'
  END IF

  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tps1)
     CALL system_CLOCK(timing%cks1,timing%clock_rate1)
  END IF

  ! Reset matrices
  elMat%Auq = 0.
  elMat%Auu = 0.
  elMat%Aul = 0.
  elMat%Alq = 0.
  elMat%Alu = 0.
  elMat%All = 0.
  elMat%S = 0.
  elMat%fh = 0.

#ifdef TEMPERATURE
  coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1 + phys%epn)
  coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1 + phys%epn)
#endif

#ifdef TOR3D
  !*************************************************************
  !               3D stuff
  !*************************************************************
  Ndim = 3                         ! Number of dimensions
  Neq = phys%Neq                  ! Number of equations
  N2D = Mesh%Nelems               ! Number of 2D elements
  Np2D = refElPol%Nnodes2D         ! Number of nodes in the 2D elements
  Np1Dpol = refElPol%Nnodes1D         ! Number of nodes in the 1D poloidal segments
  Np1Dtor = refElTor%Nnodes1D         ! Number of nodes in the 1D toroidal segments
  Npel = Np2D*refElTor%Nnodes1D    ! Number of nodes of each element
  Npfl = Np1Dpol*Np1Dtor           ! Number of nodes of each lateral face
  Npfp = Np2D                      ! Number of nodes of each poloidal face
  Ng2D = refElPol%Ngauss2d         ! Number of Gauss points in the 2D element
  Ng1Dpol = refElPol%Ngauss1d         ! Number of Gauss points in the 1D poloidal segments
  Ng1Dtor = refEltor%Ngauss1d         ! Number of Gauss points in the 1D toroidal segments
  Ngvo = Ng2D*Ng1dtor              ! Number of Gauss points for volume computations
  Ngfl = Ng1dtor*Ng1dpol           ! Number of Gauss points for toroidal faces computations
  Ngfp = Ng2D                      ! Number of Gauss points for poloidal faces computations
  Nfre = refElPol%Nfaces           ! Number of faces in the reference element
  Nfdir = Mesh%Ndir
  ! Toroidal discretization
  htor = numer%tmax/numer%ntor
  tdiv = 0.
  DO i = 1,numer%ntor
    tdiv(i + 1) = i*htor
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

  ! Local indices
  ind_loc = 0
  DO i = 1,Nfre
    DO j = 1,Npfl*Neq
      ind_loc(i,j) = refElPol%Nnodes2D*Neq + Neq*Npfl*(i - 1) + j
    END DO
  END DO

  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Np1Dpol,Np1Dtor,Neq,perm)

#else
  !*************************************************************
  !               2D stuff
  !*************************************************************
  Ndim = 2
  Neq = phys%Neq
  N2D = Mesh%Nelems
  Npel = refElPol%Nnodes2D
  Npfl = refElPol%Nfacenodes
  Nfre = refElPol%Nfaces
  Ng1d = refElPol%Ngauss1d
  Ng2d = refElPol%Ngauss2d

  ind_loc = 0
  DO i = 1,Nfre
    DO j = 1,Neq*Npfl
      ind_loc(i,j) = Neq*Npfl*(i - 1) + j
    END DO
  END DO

  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Neq*Npfl,Neq,perm)

  save_tau = switch%saveTau
  IF (save_tau) THEN
     ALLOCATE (tau_save(refElPol%Nfaces*Mesh%Nelems*refElPol%Ngauss1d,phys%neq))
     ALLOCATE (xy_g_save(refElPol%Nfaces*Mesh%Nelems*refElPol%Ngauss1d,2))
    tau_save = 0.
    xy_g_save = 0.
  ENDIF

  ! Compute shock capturing diffusion
  IF (switch%shockcp.GT.0) THEN
     CALL HDG_ShockCapturing()
  END IF

#endif
  ! reshape
  sizeu = SIZE(sol%u)
  sizel = SIZE(sol%u_tilde)
  ALLOCATE (ures(sizeu/Neq,Neq))
  ALLOCATE (lres(sizel/Neq,Neq))
  ALLOCATE (u0res(sizeu/Neq,Neq,time%tis))
  ures = TRANSPOSE(RESHAPE(sol%u,[Neq,sizeu/Neq]))
  lres = TRANSPOSE(RESHAPE(sol%u_tilde,[Neq,sizel/Neq]))
  DO i = 1,time%tis
     u0res(:,:,i) = TRANSPOSE(RESHAPE(sol%u0(:,i),[Neq,sizeu/Neq]))
  END DO
  ALLOCATE (qres(sizeu/Neq,Neq*Ndim))
  qres = TRANSPOSE(RESHAPE(sol%q,[Neq*Ndim,sizeu/Neq]))

#ifdef TOR3D
  !********************************************
  !
  !                 3D routines
  !
  !********************************************
  !************************************
  !   Loop in elements
  !************************************
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(itor,iel,itorg,iel3,tel,Xel,indbe,Bel,fluxel,inde,qe,ue,u0e,ifa,indbp,Bfp,dd,indfp,ufp,indl)&
  !$OMP PRIVATE(uefp,qefp,iface,Xfl,indbt,Bfl,isdir,indft,uft,ueft,qeft,i,Jtorel)
  ALLOCATE(Xel(Mesh%Nnodesperelem,2))
  ALLOCATE(Xfl(refElPol%Nfacenodes,2))
  ALLOCATE(Bel(refElTor%Nnodes3d,3),fluxel(refElTor%Nnodes3d),Bfl(refElTor%Nfl,3),Bfp(Mesh%Nnodesperelem,3))
  ALLOCATE(ue(refElTor%Nnodes3D,phys%Neq),u0e(refElTor%Nnodes3D,phys%Neq,time%tis))
  ALLOCATE(ufp(Mesh%Nnodesperelem,phys%Neq),uft(refElTor%Nfl,phys%Neq))
  ALLOCATE(uefp(Mesh%Nnodesperelem,phys%Neq),ueft(refElTor%Nfl,phys%Neq))
  ALLOCATE(qe(refElTor%Nnodes3D,phys%Neq*3))
  ALLOCATE(qefp(Mesh%Nnodesperelem,phys%Neq*3),qeft(refElTor%Nfl,phys%Neq*3))
  !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  DO itor = 1,ntorloc
    DO iel = 1,N2D
      ! I made a perfectly nested loop to enable omp parallelization
#ifdef PARALL
      itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
        IF (itorg == numer%ntor + 1) itorg = 1
#else
      itorg = itor
#endif
      tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))


      !call mpi_barrier(mpi_comm_world,ierr)
      !if (mpivar%glob_id.eq.0) then
      !write(6,*) "TEST 3"
      !endif
      !flush(6)
      !call mpi_barrier(mpi_comm_world,ierr)

      ! Index of 3D element
      iel3 = (itor - 1)*N2d+iel

      ! Coordinates of the nodes of the element
      Xel = Mesh%X(Mesh%T(iel,:),:)

      ! Magnetic field of the nodes of the element
      indbe = colint(tensorSumInt(Mesh%T(iel,:),(itor - 1)*(Np1Dtor - 1)*Mesh%Nnodes + &
        &Mesh%Nnodes*((/(i,i=1,Np1Dtor)/) - 1)))

      Bel = phys%B(indbe,:)
      fluxel = phys%magnetic_flux(indbe)

      ! Ohmic heating (toroidal current)
      IF (switch%ohmicsrc) THEN
        Jtorel = phys%Jtor(indbe)
      ELSE
        Jtorel = 0.
      END IF

      ! Indices to extract the elemental solution
      inde = (iel3 - 1)*Npel + (/(i,i=1,Npel)/)

      qe = qres(inde,:)
      ue = ures(inde,:)
      u0e = u0res(inde,:,:)

      ! Compute the matrices for the element
      CALL elemental_matrices_volume(iel3,Xel,tel,Bel,fluxel,qe,ue,u0e,Jtorel)

      !------------- First poloidal face-----------------
      ifa = 1

      ! Magnetic field of the nodes of the face
      indbp = (itor - 1)*Mesh%Nnodes*(Np1Dtor - 1) + Mesh%T(iel,:)
      Bfp = phys%B(indbp,:)

      ! Face solution
      dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Npfl)
      indfp = dd + (iel - 1)*Np2D+(/(i,i=1,Np2d)/)
      ufp = lres(indfp,:)

      ! Elements solution
      indl = (/(i,i=1,Np2d)/)
      uefp = ures(inde(indl),:)
      qefp = qres(inde(indl),:)


      ! Compute the matrices for the element
      CALL elemental_matrices_faces_pol(iel3,ifa,iel,Xel,tdiv(itorg),Bfp,qefp,uefp,ufp)

      !-------------- Toroidal faces ---------------------
      DO ifa=1,refElPol%Nfaces
        iface = Mesh%F(iel,ifa)
        Xfl = Mesh%X(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)

        ! Magnetic field of the nodes of the face
        indbt = colint(tensorSumInt(Mesh%T(iel,refElPol%face_nodes(ifa,:)),(itor-1)*&
          &(Np1Dtor-1)*Mesh%Nnodes+Mesh%Nnodes*((/(i,i=1,Np1Dtor) /)-1)))
        Bfl = phys%B(indbt,:)

        ! Face solution
        isdir = Mesh%Fdir(iel,ifa)
        IF (isdir) THEN
          uft = 0.
        ELSE
          dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Npfl) + N2D*Np2D
          indft = dd + (iface - 1)*Npfl + (/(i,i=1,Npfl)/)
          uft = lres(indft,:)
        ENDIF

        ! Elements solution
        ueft = ures(inde(refElTor%faceNodes3(ifa,:)),:)
        qeft = qres(inde(refElTor%faceNodes3(ifa,:)),:)

           IF (iface.LE.Mesh%Nintfaces) THEN
              CALL elemental_matrices_faces_int(iel3,ifa+1,iel,Xfl,tel,Bfl,qeft,ueft,uft)
           ELSE
              CALL elemental_matrices_faces_ext(iel3,ifa+1,iel,Xfl,tel,Bfl,qeft,ueft,uft,isdir)
           ENDIF
        ! Flip faces
           IF (Mesh%flipface(iel,ifa)) THEN
          elMat%Alq(ind_loc(ifa,:),:,iel3) = elMat%Alq(ind_loc(ifa,perm),:,iel3)
          elMat%Alu(ind_loc(ifa,:),:,iel3) = elMat%Alu(ind_loc(ifa,perm),:,iel3)
              elMat%ALL(ind_loc(ifa,:),:,iel3) = elMat%ALL(ind_loc(ifa,perm),:,iel3)
              elMat%ALL(:,ind_loc(ifa,:),iel3) = elMat%ALL(:,ind_loc(ifa,perm),iel3)
          elMat%fh(ind_loc(ifa,:),iel3) = elMat%fh(ind_loc(ifa,perm),iel3)
           END IF
      END DO

      !------------- Second poloidal face-----------------
      ifa = refElPol%Nfaces + 2

      ! Magnetic field of the nodes of the face
      indbp = itor*Mesh%Nnodes*(Np1Dtor - 1) + Mesh%T(iel,:)
      Bfp = phys%B(indbp,:)

      ! Face solution

        IF (itor==numer%ntor) THEN
        dd = 0
        ELSE
        dd = itor*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Npfl)
        ENDIF

      indfp = dd + (iel - 1)*Np2D+(/(i,i=1,Np2d)/)
      ufp = lres(indfp,:)

      ! Elements solution
      indl = (/(i,i=Npel - Np2d+1,Npel)/)
      uefp = ures(inde(indl),:)
      qefp = qres(inde(indl),:)

      ! Compute the matrices for the element
      CALL elemental_matrices_faces_pol(iel3,ifa,iel,Xel,tdiv(itorg+1),Bfp,qefp,uefp,ufp)
    END DO
  END DO
  !$OMP END DO
  DEALLOCATE(Xel,Xfl,Bel,fluxel,Bfl,Bfp)
  DEALLOCATE(ue,u0e,ufp,uft,uefp,ueft,qe,qefp,qeft)
  !$OMP END PARALLEL

  DEALLOCATE (ures,lres,u0res)
  DEALLOCATE (qres)

  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tpe1)
     CALL system_CLOCK(timing%cke1,timing%clock_rate1)
     timing%runtjac = timing%runtjac + (timing%cke1 - timing%cks1)/REAL(timing%clock_rate1)
    timing%cputjac = timing%cputjac + timing%tpe1 - timing%tps1
  END IF


CONTAINS

  !***************************************************
  ! Volume computation in 3D
  !***************************************************
  SUBROUTINE elemental_matrices_volume(iel,Xel,tel,Bel,fluxel,qe,ue,u0e,Jtorel)
    INTEGER,INTENT(IN)          :: iel
    REAL*8,INTENT(IN)           :: Xel(:,:),tel(:)
    REAL*8,INTENT(IN)           :: Bel(:,:),fluxel(:),Jtorel(:)
    REAL*8,INTENT(IN)           :: qe(:,:)
    REAL*8,INTENT(IN)           :: ue(:,:)
    REAL*8,INTENT(IN)           :: u0e(:,:,:)
    INTEGER*4                   :: g,NGaussPol,NGaussTor,igtor,igpol,i,j,k,iord
    REAL*8                      :: dvolu,dvolu1d,htor
    REAL*8                      :: J11(Ng2d),J12(Ng2d)
    REAL*8                      :: J21(Ng2d),J22(Ng2d)
    REAL*8                      :: detJ(Ng2d)
    REAL*8                      :: iJ11(Ng2d),iJ12(Ng2d)
    REAL*8                      :: iJ21(Ng2d),iJ22(Ng2d)
    REAL*8                      :: fluxg(Ng2d)
    REAL*8                      :: xy(Ng2d,2),teg(Ng1dtor)
    REAL*8                      :: ueg(Ngvo,neq),upg(Ngvo,phys%npv),u0eg(Ngvo,neq,time%tis)
    REAL*8                      :: qeg(Ngvo,neq*Ndim)
    REAL*8                      :: force(Ngvo,Neq)
    INTEGER*4,DIMENSION(Npel)   :: ind_ass,ind_asq
    REAL*8                      :: ktis(time%tis + 1)
    REAL*8                      :: Nxg(Np2D),Nyg(Np2D),Nx_ax(Np2D)
    REAL*8,POINTER              :: N1g(:),N2g(:),N3g(:)
    REAL*8                      :: N1xg(Np1Dtor),N1xg_cart(Np1Dtor)
    REAL*8,DIMENSION(Npel)      :: Nrr,Nr,Nz,Nt,Ni,NNbb
    REAL*8,DIMENSION(Npel,Npel) :: NNi,NxNi,NyNi,NtNi,NNxy
    REAL*8                      :: NxyzNi(Npel,Npel,Ndim),Nxyzg(Npel,Ndim)

    REAL*8                      :: kmult(Npfl,Npfl)
    REAL*8,PARAMETER            :: tol = 1e-12
    INTEGER*4                   :: ind2(Ng2d)
    REAL*8                      :: Bmod_nod(Npel),b_nod(Npel,3),b(Ngvo,3),Bmod(Ngvo),divbg,driftg(3),gradbmod(3)
    REAL*8                      :: bg(3),Jtor(Ngvo)
    REAL*8                      :: diff_iso_vol(Neq,Neq,Ngvo),diff_ani_vol(Neq,Neq,Ngvo)
    REAL*8                      :: sigma,sigmax,sigmay,x0,y0,A

    REAL*8,ALLOCATABLE  :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)

    !index i from 2nd to 3rd term with a 4th term as step
    ind_ass = (/(i,i=0,Neq*(Npel - 1),Neq)/)
    ind_asq = (/(i,i=0,Neq*(Npel - 1)*Ndim,Neq*Ndim)/)

    !***********************************
    !    Volume computation
    !***********************************

    ! Gauss points position
    xy = MATMUL(refElPol%N2D,Xel)

    ! Gauss points position in the toroidal direction
    teg = MATMUL(refElTor%N1d,tel)

    ! toroidal element size
    htor = tel(Np1Dtor) - tel(1)

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = SQRT(Bel(:,1)**2 + Bel(:,2)**2 + Bel(:,3)**2)
    b_nod(:,1) = Bel(:,1)/Bmod_nod
    b_nod(:,2) = Bel(:,2)/Bmod_nod
    b_nod(:,3) = Bel(:,3)/Bmod_nod

    ! Magnetic field norm and direction at Gauss points
    Bmod = MATMUL(refElTor%N3D,Bmod_nod)
    b = MATMUL(refElTor%N3D,b_nod)

    ! toroidal current at Gauss points
    IF (switch%ohmicsrc) THEN
       Jtor = MATMUL(refElPol%N3D,Jtorel)
    ELSE
      Jtor = 0.
    END IF

    ! Solution at Gauss points
    ueg = MATMUL(refElTor%N3D,ue)
    qeg = MATMUL(refElTor%N3D,qe)

    ! Compute diffusion at Gauss points
    CALL setLocalDiff(xy,ueg,diff_iso_vol,diff_ani_vol,Bmod)

    ! Solution at previous time steps,at Gauss points
    DO i = 1,time%tis
       u0eg(:,:,i) = MATMUL(refElTor%N3D,u0e(:,:,i))
    END DO

    ! Physical variables at Gauss points
    CALL cons2phys(ueg,upg)

    ! Constant sources
    ! Body force at the integration points
    CALL body_force(xy(:,1),xy(:,2),teg,force)

#ifdef NEUTRAL
!    ! Some neutral for WEST
!    IF (switch%testcase .ge. 50 .and. switch%testcase .le. 59) THEN
!      DO g = 1,Ng2d
!        DO igtor =1,Ng1Dtor
!          !IF (xy(g,1)*phys%lscale .gt. 2.36 .and. xy(g,2)*phys%lscale .lt. -0.69 ) THEN
!          IF (xy(g,1)*phys%lscale .gt. 2.446 .and. xy(g,1)*phys%lscale .lt. 2.59 .and. xy(g,2)*phys%lscale .gt. -0.7964 .and. xy(g,2)*phys%lscale .lt. -0.7304 ) THEN
!            ! Case moving equilibrium
!            ! force(g,5) = phys%puff_exp(time%it+1)
!            i = (igtor-1)*Ng2d+g
!            force(i,5) = phys%puff
!          ENDIF
!          !x0 = 2.213
!          !y0 = -0.6968
!          !sigmax = 0.02
!          !sigmay = 0.01
!          !A = phys%lscale**2/(pi*sigmax*sigmay)
!          !force(g,5) = phys%puff*A*exp(-((xy(g,1)*phys%lscale - x0)**2)/(2*sigmax**2) - ((xy(g,2)*phys%lscale - y0)**2)/(2*sigmay**2))
!        ENDDO
!      ENDDO
!    ENDIF
#endif
    ! Some sources for West cases
    IF (switch%testcase .GE. 51 .AND. switch%testcase .LE. 55) THEN
       fluxg = MATMUL(refElPol%N2D,fluxel)
      DO g = 1,Ngvo
        IF (switch%testcase == 51) THEN
             IF (fluxg(g) .LE. -0.88 .AND. fluxg(g) .GE. -0.90) THEN
            !force(g,1) = 3.20119388718018e-05
            force(g,1) = 4.782676673609557e-05
          END IF
        ELSE IF (switch%testcase == 52) THEN
          sigma = 0.1
          x0 = 0.
          A = (phys%lscale**2)/(2*PI*sigma**2)
#ifndef NEUTRAL
             force(g,1) = 0.4*A*EXP(-((fluxg(g) - x0)**2)/(2*sigma**2))
#endif
#ifdef TEMPERATURE
          force(g,3) = 0.
             force(g,4) = 30.*A*EXP(-((fluxg(g) - x0)**2)/(2*sigma**2))
#endif
          !IF (fluxg(g) .le. -0.90 .and. fluxg(g) .ge. -1.) THEN
          !  force(g,1) = 9.45155008295538e-06
          !END IF
        ELSE IF (switch%testcase == 53) THEN
             IF (fluxg(g) .LE. -0.90) THEN
            force(g,1) = 7.24032211339971e-06
          END IF
#ifdef TEMPERATURE
          force(g,3) = 18*force(g,1)
          force(g,4) = force(g,3)
#endif
        ELSE IF (switch%testcase == 54) THEN
#ifndef NEUTRAL
          sigma = 0.22
          x0 = 0.
          A = (phys%lscale**2)/(2*PI*sigma**2)
             IF (fluxg(g) .LE. 0.35) THEN
                force(g,1) = 1.*A*EXP(-((fluxg(g) - x0)**2)/(2*sigma**2))
          ENDIF
          !IF (fluxg(g) .le. -1.03) THEN
          !  force(g,1) = 0.000115575293741846
          !END IF
#endif
        ELSE IF (switch%testcase == 55) THEN
             IF (fluxg(g) .LE. -0.88 .AND. fluxg(g) .GE. -0.90) THEN
            force(g,1) = 10
          END IF
#ifdef TEMPERATURE
          force(g,3) = 18*force(g,1)
          force(g,4) = force(g,3)
#endif
        END IF
      END DO
    END IF
    !! end sources

    ! Loop in 2D Gauss points
    J11 = MATMUL(refElPol%Nxi2D,Xel(:,1))                           ! ng x 1
    J12 = MATMUL(refElPol%Nxi2D,Xel(:,2))                           ! ng x 1
    J21 = MATMUL(refElPol%Neta2D,Xel(:,1))                          ! ng x 1
    J22 = MATMUL(refElPol%Neta2D,Xel(:,2))                          ! ng x 1
    detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
    iJ11 = J22/detJ
    iJ12 = -J12/detJ
    iJ21 = -J21/detJ
    iJ22 = J11/detJ

    ! Coefficient time integration scheme
    CALL setTimeIntegrationCoefficients(ktis)

    NgaussPol = refElPol%NGauss2D
    NgaussTor = refElTor%NGauss1D


    ! Allocate temporary matrices
    ALLOCATE(Auq(Npel,Npel, neq*neq*ndim  ))
    ALLOCATE(Auu(Npel,Npel, neq*neq  ))
    ALLOCATE(rhs(Npel,Neq))
    Auq = 0.
    Auu = 0.
    rhs = 0.
    DO igtor = 1,NGaussTor
      N1g => refElTor%N1D(igtor,:)         ! Toroidal shape function
      N1xg_cart = refElTor%Nxi1D(igtor,:)*2/htor       ! Toroidal shape function derivative
      dvolu1d = 0.5*refElTor%gauss_weights1D(igtor)*htor ! Toroidal contribution to the elemental volume

      DO igpol = 1,NGaussPol
        g = (igtor - 1)*NGaussPol + igpol

        ! Poloidal shape functions and derivatives
        N2g => refElPol%N2D(igpol,:)
        Nxg = iJ11(igpol)*refElPol%Nxi2D(igpol,:) + iJ12(igpol)*refElPol%Neta2D(igpol,:)
        Nyg = iJ21(igpol)*refElPol%Nxi2D(igpol,:) + iJ22(igpol)*refElPol%Neta2D(igpol,:)

        ! 3D integration weight
        dvolu = refElPol%gauss_weights2D(igpol)*detJ(igpol)*dvolu1d

        IF (switch%axisym) THEN
          dvolu = dvolu*xy(igpol,1)
          N1xg = N1xg_cart/xy(igpol,1)
          Nx_ax = Nxg + 1./xy(igpol,1)*N2g
        ELSE
          N1xg = N1xg_cart
          Nx_ax = Nxg
        END IF

        ! 3D shape functions
        N3g => refElTor%N3D(g,:)               ! 3D shape function
        Nrr = col(TensorProduct(Nx_ax,N1g))    ! 3D shape function,derivative in r for computing the divergence
        Nr = col(TensorProduct(Nxg,N1g))      ! 3D shape function,derivative in r
        Nz = col(TensorProduct(Nyg,N1g))      ! 3D shape function,derivative in z
        Nt = col(TensorProduct(N2g,N1xg))     ! 3D shape function,derivative in t

        ! Shape functions products
        Ni = N3g*dvolu                                                ! Npel x 1
        NNi = tensorProduct(N3g,Ni)                                    ! Npel x Npel
        NxNi = tensorProduct(Nr,Ni)                                     ! Npel x Npel
        NyNi = tensorProduct(Nz,Ni)                                     ! Npel x Npel
        NtNi = tensorProduct(Nt,Ni)                                     ! Npel x Npel
        NNxy = b(g,1)*NxNi + b(g,2)*NyNi + b(g,3)*NtNi                      ! Npel x Npel
        NxyzNi(:,:,1) = NxNi
        NxyzNi(:,:,2) = NyNi
        NxyzNi(:,:,3) = NtNi                                           ! Npel x Npel x 3
        NNbb = (Nr*b(g,1) + Nz*b(g,2) + Nt*b(g,3))*dvolu                   ! Npel x 1
        Nxyzg(:,1) = Nr*dvolu
        Nxyzg(:,2) = Nz*dvolu
        Nxyzg(:,3) = Nt*dvolu

          divbg = dot_PRODUCT(Nrr,b_nod(:,1)) + dot_PRODUCT(Nz,b_nod(:,2)) + dot_PRODUCT(Nt,b_nod(:,3))

        ! Diamagnetic drift !TODO: verify drift intensity in isothermal and non-isothermal cases
        driftg = 0.
        gradbmod = 0.
          gradbmod(1) = dot_PRODUCT(Nr,Bmod_nod)
          gradbmod(2) = dot_PRODUCT(Nz,Bmod_nod)
          gradbmod(3) = dot_PRODUCT(Nt,Bmod_nod)
        bg = b(g,:)
          CALL cross_product(bg,gradbmod,driftg)
        driftg = phys%dfcoef*driftg/Bmod(g)

        CALL assemblyVolumeContribution(Auq,Auu,rhs,b(g,:),divbg,driftg,Bmod(g),force(g,:),&
          &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upg(g,:),&
          &ueg(g,:),qeg(g,:),u0eg(g,:,:),xy(g,:),Jtor(g))
      END DO ! END loop in volume Gauss points
    END DO
    CALL do_assembly(Auq,Auu,rhs,ind_ass,ind_asq,iel)
    DEALLOCATE(Auq,Auu,rhs)

  END SUBROUTINE elemental_matrices_volume

  !*****************************************
  ! Poloidal faces computations in 3D
  !*****************************************
  SUBROUTINE elemental_matrices_faces_pol(iel,ifa,iel2,Xfp,tg,Bfp,qef,uef,uf)
    INTEGER*4,INTENT(IN)  :: iel,ifa,iel2
    REAL*8,INTENT(IN)     :: Xfp(:,:),tg(1)
    REAL*8,INTENT(IN)     :: Bfp(:,:)
    REAL*8,INTENT(IN)     :: qef(:,:)
    REAL*8,INTENT(IN)     :: uef(:,:),uf(:,:)
    REAL*8                :: ufg(Ng2d,Neq),uefg(Ng2d,Neq),upgf(Ng2d,phys%npv)
    REAL*8                :: qfg(Ng2d,Neq*Ndim)
    INTEGER*4             :: ind_asf(Np2D),ind_ash(Np2D)
    INTEGER*4             :: g,NGauss,i,j,lel
    REAL*8                :: dsurf(Ng2d),bn
    REAL*8                :: xyf(Ng2d,2)
    REAL*8                :: J11(Ng2d),J12(Ng2d)
    REAL*8                :: J21(Ng2d),J22(Ng2d)
    REAL*8                :: detJ(Ng2d)
    REAL*8                :: iJ11(Ng2d),iJ12(Ng2d)
    REAL*8                :: iJ21(Ng2d),iJ22(Ng2d)
    REAL*8,POINTER        :: Nfg(:)
    REAL*8                :: NNif(Np2D,Np2D),Nif(Np2D),Nfbn(Np2D)
    REAL*8                :: n_g(Ng2d,3)
    INTEGER*4             :: ind_ff(Neq*Np2D),ind_fe(Neq*Np2D),ind_fg(Neq*Ndim*Np2D)
    REAL*8,PARAMETER     :: tol = 1e-12
    REAL*8                :: Bmod_nod(Np2D),b_nod(Np2D,3),b(Ng2d,3),Bmod(Ng2d)
    REAL*8                :: tau(Neq,Neq)
    REAL*8                :: diff_iso_fac(Neq,Neq,Ng2D),diff_ani_fac(Neq,Neq,Ng2D)

    ind_asf = (/(i,i=0,Neq*(Np2D-1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Np2D-1)*Ndim,Neq*Ndim)/)

    ! Gauss points position
    xyf = MATMUL(refElPol%N2D,Xfp)

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = SQRT(Bfp(:,1)**2 + Bfp(:,2)**2 + Bfp(:,3)**2)
    b_nod(:,1) = Bfp(:,1)/Bmod_nod
    b_nod(:,2) = Bfp(:,2)/Bmod_nod
    b_nod(:,3) = Bfp(:,3)/Bmod_nod
    ! Magnetic field norm and direction at Gauss points
    Bmod = MATMUL(refElPol%N2D,Bmod_nod)
    b = MATMUL(refElPol%N2D,b_nod)

    ! Trace solution at face Gauss points
    ufg = MATMUL(refElPol%N2D,uf)

    ! Compute diffusion at Gauss points
    CALL setLocalDiff(xyf,ufg,diff_iso_fac,diff_ani_fac,Bmod)

    ! Loop in 2D Gauss points
    Ngauss = Ng2d

    ! Physical variables related to the trace solution
    CALL cons2phys(ufg,upgf)

    J11 = MATMUL(refElPol%Nxi2D,Xfp(:,1))                           ! ng x 1
    J12 = MATMUL(refElPol%Nxi2D,Xfp(:,2))                           ! ng x 1
    J21 = MATMUL(refElPol%Neta2D,Xfp(:,1))                          ! ng x 1
    J22 = MATMUL(refElPol%Neta2D,Xfp(:,2))                          ! ng x 1
    detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
    iJ11 = J22/detJ
    iJ12 = -J12/detJ
    iJ21 = -J21/detJ
    iJ22 = J11/detJ
    dsurf = refElPol%gauss_weights2D*detJ

    IF (ifa == 1) THEN
      ind_ff = (/(i,i=1,Np2D*Neq)/)
      ind_fe = (/(i,i=1,Np2D*Neq)/)
      ind_fg = (/(i,i=1,Np2D*Ndim*Neq)/)
      ! Exterior normal
      n_g = 0.; n_g(:,3) = -1
    ELSE
      ind_ff = Np2D*Neq + refElPol%Nfaces*Npfl*Neq + (/(i,i=1,Np2D*Neq)/)
      ind_fe = Np2d*(Np1dTor - 1)*Neq + (/(i,i=1,Np2D*Neq)/)
      ind_fg = Np2d*(Np1dTor - 1)*Ndim*Neq + (/(i,i=1,Np2D*Ndim*Neq)/)
      ! Exterior normal
      n_g = 0.; n_g(:,3) = 1
    ENDIF

    ! Element solution at face Gauss points
    uefg = MATMUL(refElPol%N2D,uef)
    ! Gradient solution at face gauss points
    qfg = MATMUL(refElPol%N2D,qef)

    DO g = 1,NGauss

      ! Shape functions product
      Nfg => refElPol%N2D(g,:)
       bn = dot_PRODUCT(b(g,:),n_g(g,:))
      NNif = tensorProduct(Nfg,Nfg)*dsurf(g)
      Nif = Nfg*dsurf(g)
      Nfbn = bn*Nfg*dsurf(g)

      ! Compute the stabilization term
      tau = 0.
      IF (numer%stab == 1) THEN
        ! Constant stabilization
        DO i = 1,Neq
          tau(i,i) = numer%tau(i)
        END DO
      ELSE
        ! Non constant stabilization
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b(g,:),Bmod(g),n_g(g,:),iel2,ifa,0.,xyf(g,:),tau)
        ELSE
          CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),0.,iel2,tau)
        ENDIF
      END IF

      ! Assembly local contributions
      CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
        n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)
    END DO
    !      END DO

  END SUBROUTINE elemental_matrices_faces_pol

  !***************************************************
  ! Interior faces computation in 3D
  !***************************************************
  SUBROUTINE elemental_matrices_faces_int(iel,ifa,iel2,Xfl,tel,Bfl,qef,uef,uf)
    INTEGER,INTENT(IN)        :: iel,ifa,iel2
    REAL*8,INTENT(IN)         :: Xfl(:,:),tel(:)
    REAL*8,INTENT(IN)         :: Bfl(:,:)
    REAL*8,INTENT(IN)         :: qef(:,:)
    REAL*8,INTENT(IN)         :: uef(:,:),uf(:,:)

    INTEGER*4                 :: g,i,j,k,igtor,igpol
    REAL*8                    :: xyf(Ng1Dpol,2),teg(Ng1dtor)
    REAL*8                    :: xyDer(Ng1Dpol,2),xydNorm_g(Ng1Dpol)
    REAL*8                    :: ufg(Ngfl,neq),uefg(Ngfl,neq),upgf(Ngfl,phys%npv)
    REAL*8                    :: dsurf(Ngfl),dsurfg
    REAL*8                    :: qfg(Ngfl,neq*Ndim)
    INTEGER*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
    INTEGER*4,DIMENSION(Npfl) :: ind_asf,ind_ash
    INTEGER*4,DIMENSION(Npfl) :: indf,ind_if,ind_jf,ind_kf
    INTEGER*4                 :: ind(Ng1Dpol)
    INTEGER*4                 :: permsing(Npfl)
    REAL*8                    :: t_g(Ng1dpol,2),n_g(Ngfl,3),bn
    REAL*8                    :: NNif(Npfl,Npfl),Nif(Npfl),Nfbn(Npfl)
    REAL*8,POINTER            :: Nfg(:)
    REAL*8                    :: tau(Neq,Neq)
    REAL*8,PARAMETER          :: tol = 1e-12
    REAL*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ngfl,3),Bmod(Ngfl)
    REAL*8                    :: diff_iso_fac(Neq,Neq,Ngfl),diff_ani_fac(Neq,Neq,Ngfl)

    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    ! toroidal element size
    htor = tel(Np1Dtor) - tel(1)

    ! Gauss points position in the toroidal direction
    teg = MATMUL(refElTor%N1d,tel)

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = SQRT(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
    b_nod(:,1) = Bfl(:,1)/Bmod_nod
    b_nod(:,2) = Bfl(:,2)/Bmod_nod
    b_nod(:,3) = Bfl(:,3)/Bmod_nod

    ! Indices
    ind_ff = Np2d*Neq + (ifa - 2)*Npfl*Neq + (/(i,i=1,Npfl*Neq)/)
    ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*Neq))
    ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*3*Neq))

    ! Coordinates,derivatives and trace solution at face Gauss points
    IF (Mesh%flipFace(iel2,ifa - 1)) THEN
       CALL set_permutations(Np1Dpol,Np1Dtor,1,permsing)
       ufg = MATMUL(refElTor%sFTF,uf(permsing,:))
    ELSE
       ufg = MATMUL(refElTor%sFTF,uf)
    END IF

    xyf = MATMUL(refElPol%N1D,Xfl)
    ! Shape function derivatives at Gauss points
    xyDer = MATMUL(refElPol%Nxi1D,Xfl)
    ! Magnetic field norm and direction at Gauss points
    Bmod = MATMUL(refElTor%sFTF,Bmod_nod)
    b = MATMUL(refElTor%sFTF,b_nod)

    ! Element solution at face Gauss points
    uefg = MATMUL(refElTor%sFTF,uef)
    ! Gradient solution at face gauss points
    qfg = MATMUL(refElTor%sFTF,qef)
    ! Compute diffusion at faces Gauss points
    CALL setLocalDiff(xyf,uefg,diff_iso_fac,diff_ani_fac,Bmod)

    ! Physical variables at face Gauss points
    CALL cons2phys(ufg,upgf)

    ! Compute dsurf
    xydNorm_g = SQRT(xyDer(:,1)**2 + xyDer(:,2)**2)
    dsurf = col(tensorProduct(refElPol%gauss_weights1D*xydNorm_g,refElTor%gauss_weights1D*0.5*htor))

    ! Compute exterior normal
    t_g(:,1) = xyDer(:,1)/xydNorm_g
    t_g(:,2) = xyDer(:,2)/xydNorm_g
    n_g = 0.
    DO i = 1,Ng1dTor
      ind = (i - 1)*Ng1dPol + (/(j,j=1,Ng1dPol)/)
      n_g(ind,1) = t_g(:,2)
      n_g(ind,2) = -t_g(:,1)
    END DO

    !*****************************
    ! Loop in face Gauss points
    !*****************************
    DO igtor = 1,Ng1dTor
      DO igpol = 1,Ng1dPol

        g = (igtor - 1)*Ng1dPol + igpol

        ! Face shape functions
        Nfg => refElTor%sFTF(g,:)

        IF (switch%axisym) THEN
          dsurfg = dsurf(g)*xyf(igpol,1)
        ELSE
          dsurfg = dsurf(g)
        END IF

        ! Shape functions product
          bn = dot_PRODUCT(b(g,:),n_g(g,:))
        NNif = tensorProduct(Nfg,Nfg)*dsurfg
        Nif = Nfg*dsurfg
        Nfbn = bn*Nfg*dsurfg

        ! Compute the stabilization term
        tau = 0.
        IF (numer%stab == 1) THEN
          ! Constant stabilization
          DO i = 1,Neq
            tau(i,i) = numer%tau(i)
          END DO
        ELSE
          ! Non constant stabilization
          ! Compute tau in the Gauss points
          IF (numer%stab < 6) THEN
            CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b(g,:),Bmod(g),n_g(g,:),iel2,ifa,0.,xyf(g,:),tau)
          ELSE
            CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),0.,iel2,tau)
          ENDIF
        END IF

        ! Assembly local contributions
        CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
          n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)

      END DO ! Gauss points
    END DO
    !      END DO ! 2 elements

  END SUBROUTINE elemental_matrices_faces_int

  !***************************************************
  ! Exterior faces computation in 3D
  !***************************************************
  SUBROUTINE elemental_matrices_faces_ext(iel,ifa,iel2,Xfl,tel,Bfl,qef,uef,uf,isdir)
    INTEGER,INTENT(IN)                             :: iel,ifa,iel2
    REAL*8,INTENT(IN)         :: Xfl(:,:),tel(:)
    REAL*8,INTENT(IN)         :: Bfl(:,:)
    LOGICAL,INTENT(IN)        :: isdir
    REAL*8,INTENT(IN)         :: qef(:,:)
    REAL*8,INTENT(INOUT)      :: uef(:,:),uf(:,:)

    INTEGER*4                 :: i,j,k,g,igtor,igpol
    REAL*8                    :: xyf(Ng1Dpol,2),teg(Ng1dtor)
    REAL*8                    :: xyDer(Ng1Dpol,2),xydNorm_g(Ng1Dpol)
    REAL*8                    :: ufg(Ngfl,neq),uefg(Ngfl,neq),upgf(Ngfl,phys%npv)
    REAL*8                    :: dsurf(Ngfl),dsurfg
    REAL*8                    :: qfg(Ngfl,neq*Ndim)
    INTEGER*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
    INTEGER*4,DIMENSION(Npfl) :: ind_asf,ind_ash
    INTEGER*4,DIMENSION(Npfl) :: indf,ind_if,ind_jf,ind_kf
    INTEGER*4                 :: ind(Ng1Dpol)
    INTEGER*4                 :: permsing(Npfl)
    REAL                      :: isext
    REAL*8                    :: t_g(Ng1Dpol,2),n_g(Ngfl,3),bn
    REAL*8                    :: NNif(Npfl,Npfl),Nif(Npfl),Nfbn(Npfl)
    REAL*8,POINTER            :: Nfg(:)
    REAL*8                    :: tau(Neq,Neq)

    REAL*8,PARAMETER         :: tol = 1e-12
    REAL*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ngfl,3),Bmod(Ngfl)
    REAL*8                    :: diff_iso_fac(Neq,Neq,Ngfl),diff_ani_fac(Neq,Neq,Ngfl)

    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    ! toroidal element size
    htor = tel(Np1Dtor) - tel(1)

    ! Gauss points position in the toroidal direction
    teg = MATMUL(refElTor%N1d,tel)

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = SQRT(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
    b_nod(:,1) = Bfl(:,1)/Bmod_nod
    b_nod(:,2) = Bfl(:,2)/Bmod_nod
    b_nod(:,3) = Bfl(:,3)/Bmod_nod

    ! Magnetic field norm and direction at Gauss points
    Bmod = MATMUL(refElTor%sFTF,Bmod_nod)
    b = MATMUL(refElTor%sFTF,b_nod)

    ! Indices
    ind_ff = Np2d*Neq + (ifa - 2)*Npfl*Neq + (/(i,i=1,Npfl*Neq)/)
    ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*Neq))
    ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*3*Neq))

    ! Trace solution at face Gauss points
    xyf = MATMUL(refElPol%N1D,Xfl)
    xyDer = MATMUL(refElPol%Nxi1D,Xfl)

    IF (isdir) THEN
       CALL analytical_solution(xyf(:,1),xyf(:,2),teg,ufg)
    ELSE
#ifdef PARALL
      IF (Mesh%flipface(iel2,ifa - 1)) THEN
          CALL set_permutations(Np1Dpol,Np1Dtor,1,permsing)
        uf = uf(permsing,:)
      END IF
      ! TODO: VERIFY IF I NEED TO FLIP ALSO xyf,b and Bmod in this case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
       ufg = MATMUL(refElTor%sFTF,uf)
    END IF

    ! Element solution at face Gauss points
    uefg = MATMUL(refElTor%sFTF,uef)
    ! Gradient solution at face gauss points
    qfg = MATMUL(refElTor%sFTF,qef)

    ! Compute diffusion at faces Gauss points
    CALL setLocalDiff(xyf,uefg,diff_iso_fac,diff_ani_fac)

    ! Physical variables at face Gauss points
    CALL cons2phys(ufg,upgf)

    ! Compute dsurf
    xydNorm_g = SQRT(xyDer(:,1)**2 + xyDer(:,2)**2)
    dsurf = col(tensorProduct(refElPol%gauss_weights1D*xydNorm_g,refElTor%gauss_weights1D*0.5*htor))

    ! Compute exterior normal
    t_g(:,1) = xyDer(:,1)/xydNorm_g
    t_g(:,2) = xyDer(:,2)/xydNorm_g
    n_g = 0.
    DO i = 1,Ng1dTor
      ind = (i - 1)*Ng1dPol + (/(j,j=1,Ng1dPol)/)
      n_g(ind,1) = t_g(:,2)
      n_g(ind,2) = -t_g(:,1)
    END DO

    !*****************************
    ! Loop in face Gauss points
    !*****************************
    DO igtor = 1,Ng1dTor
      DO igpol = 1,Ng1dPol

        g = (igtor - 1)*Ng1dPol + igpol

        ! Face shape functions
        Nfg => refElTor%sFTF(g,:)

        IF (switch%axisym) THEN
          dsurfg = dsurf(g)*xyf(igpol,1)
        ELSE
          dsurfg = dsurf(g)
        END IF

        ! Shape functions product
          bn = dot_PRODUCT(b(g,:),n_g(g,:))
        NNif = tensorProduct(Nfg,Nfg)*dsurfg
        Nif = Nfg*dsurfg
        Nfbn = bn*Nfg*dsurfg

        ! Compute the stabilization term
        isext = 1.
#ifdef PARALL
          IF (Mesh%boundaryFlag(Mesh%F(iel2,ifa - 1) - Mesh%Nintfaces) .EQ. 0) THEN
          isext = 0.
        END IF
#endif
        tau = 0.
        IF (numer%stab == 1) THEN
          ! Constant stabilization
          DO i = 1,Neq
            tau(i,i) = numer%tau(i)
          END DO
        ELSE
          ! Non constant stabilization
          ! Compute tau in the Gauss points
          IF (numer%stab < 6) THEN
                CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),qfg(g,:), b(g,:),n_g(g,:),iel2,isext,xyf(g,:),tau)

          ELSE
            CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),isext,iel2,tau)
          ENDIF
        END IF

        ! Assembly local contributions
#ifdef PARALL
          IF (Mesh%boundaryFlag(Mesh%F(iel2,ifa - 1) - Mesh%Nintfaces) .EQ. 0) THEN
          ! Ghost face: assembly it as interior
          CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
            n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)
        ELSE
          CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
            n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)
        ENDIF
#else
        CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
          n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)

#endif
      END DO ! Gauss points
    END DO

  END SUBROUTINE elemental_matrices_faces_ext

  !*****************************************
  ! Set permutations for flipping faces
  !****************************************
  SUBROUTINE set_permutations(Np1Dpol,Np1Dtor,Neq,perm)
    INTEGER,INTENT(IN)  :: Np1Dpol,Np1Dtor,Neq
    INTEGER,INTENT(OUT) :: perm(:)
    INTEGER              :: i,j,k
    INTEGER              :: temp_pol(Np1Dpol),temp_tor(Np1Dtor),aux(Np1Dpol,Np1Dtor)

    temp_pol = (/(i,i=1,Np1Dpol)/)
    temp_tor = Np1Dpol*((/(i,i=1,Np1Dtor)/) - 1)
    aux = TensorSumInt(temp_pol,temp_tor)
    DO j = 1,Np1Dtor
      DO i = 1,Np1Dpol
        DO k = 1,Neq
          perm((j - 1)*Np1Dpol*Neq + (i - 1)*Neq + k) = (aux(Np1Dpol - i + 1,j) - 1)*Neq + k
        END DO
      END DO
    END DO
  END SUBROUTINE set_permutations

#else
!TOR3D
  phys%heating_amplitude = phys%heating_power/2./PI**2/phys%heating_sigmar/phys%heating_sigmaz/phys%r_axis
  !********************************************
  !
  !                 2D routines
  !
  !********************************************

  !************************************
  !   Loop in elements in 2D
  !************************************
#ifdef KEQUATION
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(iel,ifa,iface,inde,indf,Xel,Xfl,i,qe,qef,ue,uef,uf,u0e,Bel,Bfl,fluxel,omegael,q_cylel,psiel,psifl,q_cylfl,isdir,Jtorel,El_n,El_nn) &
  !$OMP PRIVATE(Xg_el,diff_nn_Vol_el,diff_nn_Fac_el,v_nn_Vol_el,v_nn_Fac_el,xy_g_save,xy_g_save_el,tau_save,tau_save_el)&
  !$OMP FIRSTPRIVATE(phys,Mesh)
#else
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(iel,ifa,iface,inde,indf,Xel,Xfl,i,qe,qef,ue,uef,uf,u0e,Bel,Bfl,fluxel,psiel,psifl,isdir,Jtorel,El_n,El_nn) &
  !$OMP PRIVATE(Xg_el,diff_nn_Vol_el,diff_nn_Fac_el,v_nn_Vol_el,v_nn_Fac_el,xy_g_save,xy_g_save_el,tau_save,tau_save_el) &
  !$OMP FIRSTPRIVATE(phys,Mesh)
#endif  
  ALLOCATE(Xel(Mesh%Nnodesperelem,2))
  ALLOCATE(Xfl(refElPol%Nfacenodes,2))

  n = 0.
  nn = 0.
  !$OMP DO SCHEDULE(STATIC) REDUCTION(+:n,nn)
  DO iel = 1,N2D

    ! Coordinates of the nodes of the element
    Xel = Mesh%X(Mesh%T(iel,:),:)

    ! Magnetic field of the nodes of the element
    Bel = phys%B(Mesh%T(iel,:),:)
    fluxel = phys%magnetic_flux(Mesh%T(iel,:))
    
    ! Normalized magnetic flux of the nodes of the element: PSI el
    psiel = phys%magnetic_psi(Mesh%T(iel,:))

#ifdef KEQUATION
    !omega and q_cyl on nodes of the element
    
     IF (switch%testcase == 60) THEN
      q_cylel = geom%q
      ! to finish this
        omegael = SQRT(Bel(:,1)**2+Bel(:,2)**2+Bel(:,3)**2)*simpar%refval_charge/simpar%refval_mass*simpar%refval_time
     ELSE
      q_cylel = phys%q_cyl(Mesh%T(iel,:))
      omegael = phys%omega(Mesh%T(iel,:))
     ENDIF
#endif
    
    ! Ohmic heating (toroidal current)
    IF (switch%ohmicsrc) THEN
      Jtorel = phys%Jtor(Mesh%T(iel,:))
    ELSE
      Jtorel = 0.
    END IF

    ! Indices to extract the elemental and face solution
    inde = (iel - 1)*Npel + (/(i,i=1,Npel)/)

    qe = qres(inde,:)
    ue = ures(inde,:)
    u0e = u0res(inde,:,:)

    ! Compute the matrices for the element
#ifndef KEQUATION
    CALL elemental_matrices_volume(iel,Xel,Bel,fluxel,psiel,qe,ue,u0e,Jtorel,El_n,El_nn,diff_nn_Vol_el,v_nn_Vol_el,Xg_el)
#else
    CALL elemental_matrices_volume(iel,Xel,Bel,fluxel,omegael,q_cylel,psiel,qe,ue,u0e,Jtorel,El_n,El_nn,diff_nn_Vol_el,v_nn_Vol_el,Xg_el)
#endif
     IF (save_tau) THEN
       inddiff_nn_Vol = (iel - 1)*refElPol%NGauss2D+(/(i,i=1,refElPol%NGauss2D)/)
       phys%diff_nn_Vol(inddiff_nn_Vol) = diff_nn_Vol_el
       phys%v_nn_Vol(inddiff_nn_Vol,:) = v_nn_Vol_el
       Mesh%Xg(inddiff_nn_Vol,:) = Xg_el
     ENDIF

    ! Compute total plasma and neutral density (don't add contribution of ghost elements)
#ifdef PARALL
     IF (Mesh%ghostElems(iel) .EQ. 0) THEN
#endif
      n  = n + El_n
      nn = nn + El_nn
#ifdef PARALL
    ENDIF
#endif
    ! Loop in local faces
     IF (save_tau) THEN
       diff_nn_Fac_el = 0.
       v_nn_Fac_el = 0.
       tau_save_el = 0.
       xy_g_save_el = 0;
     ENDIF
    
    DO ifa=1,refElPol%Nfaces
      iface = Mesh%F(iel,ifa)
      isdir = Mesh%Fdir(iel,ifa)

      ! Coordinates of the nodes of the face
      Xfl = Mesh%X(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)

      ! Magnetic field of the nodes of the face
      Bfl = phys%B(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)
   
      ! Normalized magnetic flux of the nodes of the face: PSI fl
      psifl = phys%magnetic_psi(Mesh%T(iel,refElPol%face_nodes(ifa,:)))
            
      ! Face solution
      indf = (iface-1)*Npfl + (/(i,i=1,Npfl)/)
      uf = lres(indf,:)

      ! Elements solution
      inde = (iel - 1)*Npel + (/(i,i=1,Npel)/)
      uef = ures(inde(refElPol%face_nodes(ifa,:)),:)
      qef = qres(inde(refElPol%face_nodes(ifa,:)),:)
#ifndef KEQUATION
        IF (iface.LE.Mesh%Nintfaces) THEN
        CALL elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)		 
        ELSE
           IF (Mesh%periodic_faces(iface-Mesh%Nintfaces).EQ.0) THEN
          CALL elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
           ELSE
          ! periodic face
          CALL elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
        endif
      endif
#else
      if (switch%testcase == 60) then
        q_cylfl(:) = geom%q
      else
        q_cylfl(:) = phys%q_cyl(Mesh%T(iel,refElPol%face_nodes(ifa,:)))
      endif

      if (iface.le.Mesh%Nintfaces) then
        CALL elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)		 
      else
        if (Mesh%periodic_faces(iface-Mesh%Nintfaces).eq.0) then
          CALL elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,Bfl,psifl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
        else
          ! periodic face
          CALL elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
        endif
      endif
#endif
	 
      ! Flip faces
        IF (Mesh%flipface(iel,ifa)) THEN
        elMat%Alq(ind_loc(ifa,:),:,iel) = elMat%Alq(ind_loc(ifa,perm),:,iel)
        elMat%Alu(ind_loc(ifa,:),:,iel) = elMat%Alu(ind_loc(ifa,perm),:,iel)
           elMat%ALL(ind_loc(ifa,:),:,iel) = elMat%ALL(ind_loc(ifa,perm),:,iel)
           elMat%ALL(:,ind_loc(ifa,:),iel) = elMat%ALL(:,ind_loc(ifa,perm),iel)
        elMat%fh(ind_loc(ifa,:),iel) = elMat%fh(ind_loc(ifa,perm),iel)
        END IF
    END DO
	
     IF (save_tau) THEN
       indtausave = (iel - 1)*refElPol%Nfaces*refElPol%Ngauss1d+(/(i,i=1,refElPol%Nfaces*refElPol%Ngauss1d)/)
	   phys%diff_nn_Fac(indtausave) = diff_nn_Fac_el
	   phys%v_nn_Fac(indtausave,:) = v_nn_Fac_el
	   tau_save(indtausave,:) = tau_save_el
	   Mesh%Xgf(indtausave,:) = xy_g_save_el
       xy_g_save(indtausave,:) = xy_g_save_el
     END IF


  END DO
  !$OMP END DO
  DEALLOCATE(Xel,Xfl)
  !$OMP END PARALLEL

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, n, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, nn, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

  IF (MPIvar%glob_id.EQ.0) THEN
     IF((switch%ME .EQV. .TRUE.) .AND. (switch%testcase .GE. 80)) THEN
        WRITE(6,*) 'D_n = ', phys%ME_diff_n*simpar%refval_length**2/simpar%refval_time
     ENDIF
     totaln = n + nn
     WRITE(6,*) 'n = ',n
     WRITE(6,*) 'nn = ',nn
     WRITE(6,*) 'total n = ',totaln
   ENDIF

  DEALLOCATE (ures,lres,u0res)
  DEALLOCATE (qres)

  IF (save_tau) THEN
     WRITE (6,*) "Saving Dnn in 2D gauss points"
     WRITE (6,*) "Saving tau in the faces"
     CALL saveMatrix(tau_save,'tau_save')
     CALL saveMatrix(xy_g_save,'xy_g_save')
     DEALLOCATE (tau_save,xy_g_save)
     WRITE (6,*) "Done saving tau!"
  ENDIF

  IF (utils%printint > 1) THEN
    WRITE (6,*) "Done!"
  END IF


  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tpe1)
     CALL system_CLOCK(timing%cke1,timing%clock_rate1)
     timing%runtjac = timing%runtjac + (timing%cke1 - timing%cks1)/REAL(timing%clock_rate1)
    timing%cputjac = timing%cputjac + timing%tpe1 - timing%tps1
  END IF


CONTAINS

  !***************************************************
  ! Volume computation in 2D
  !***************************************************
#ifndef KEQUATION
  SUBROUTINE elemental_matrices_volume(iel,Xel,Bel,fluxel,psiel,qe,ue,u0e,Jtorel,El_n,El_nn,diff_nn_Vol_el,v_nn_Vol_el,Xg_el)
#else
  SUBROUTINE elemental_matrices_volume(iel,Xel,Bel,fluxel,omegael,q_cylel,psiel,qe,ue,u0e,Jtorel,El_n,El_nn,diff_nn_Vol_el,v_nn_Vol_el,Xg_el)
#endif
      INTEGER,INTENT(IN)            :: iel
      REAL*8,INTENT(IN)             :: Xel(:,:)
      REAL*8,INTENT(IN)             :: Bel(:,:),fluxel(:),psiel(:),Jtorel(:)
#ifdef KEQUATION
      REAL*8,INTENT(IN)             :: omegael(:),q_cylel(:)
#endif
      REAL*8,INTENT(IN)             :: qe(:,:)
      REAL*8,INTENT(IN)             :: ue(:,:),u0e(:,:,:)
      REAL*8,INTENT(OUT)            :: El_n,El_nn
      REAL*8,INTENT(OUT)            :: diff_nn_Vol_el(Ng2D),v_nn_Vol_el(Ng2D,ndim),Xg_el(Ng2D,ndim)
      INTEGER*4                     :: g,NGauss,i
      REAL*8                        :: dvolu
      REAL*8                        :: xy(Ng2d,ndim),ueg(Ng2d,neq),u0eg(Ng2d,neq,time%tis)
      REAL*8                        :: force(Ng2d,Neq)
      REAL*8                        :: qeg(Ng2d,neq*Ndim)
      REAL*8                        :: J11(Ng2d),J12(Ng2d)
      REAL*8                        :: J21(Ng2d),J22(Ng2d)
      REAL*8                        :: detJ(Ng2d)
      REAL*8                        :: iJ11(Ng2d),iJ12(Ng2d)
      REAL*8                        :: iJ21(Ng2d),iJ22(Ng2d)
      REAL*8                        :: fluxg(Ng2d),max_flux2D,min_flux2D, Psig(Ng2d)
      INTEGER*4,DIMENSION(Npel)     :: ind_ass,ind_asq
      REAL*8                        :: ktis(time%tis + 1)
      REAL*8,DIMENSION(Npel)        :: Ni,Nxg,Nyg,NNbb,Nx_ax
      REAL*8,DIMENSION(Npel,Npel)   :: NxNi,NyNi,NNxy,NNi
      REAL*8                        :: NxyzNi(Npel,Npel,3),Nxyzg(Npel,3)
      REAL*8                        :: upg(Ng2d,phys%npv)
      REAL*8                        :: Bmod_nod(Npel),b_nod(Npel,3),b(Ng2d,3),Bmod(Ng2d),divbg,driftg(3),gradbmod(3)
#ifdef KEQUATION
      REAL*8                        :: b_tor_nod(Npel),b_tor(Ng2d),gradbtor(3)
      REAL*8                        :: omega(Ng2d),q_cyl(Ng2d)
#endif
    real*8                        :: bg(3), Jtor(Ng2d)
    real*8                        :: diff_iso_vol(Neq,Neq,Ng2d),diff_ani_vol(Neq,Neq,Ng2d)
    real*8,allocatable            :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)
    real*8                        :: auxdiffsc(Ng2d)
    real*8                        :: Pi,sigma,sigmax,sigmay,x0,y0,A,r
    real*8                        :: th_n = 1.e-14
    real*8                        :: Vnng(Ndim)
    
      IF (save_tau) THEN
       Xg_el = 0.
       diff_nn_Vol_el = 0.
       v_nn_Vol_el = 0.
      ENDIF

    ind_ass = (/(i,i=0,Neq*(Npel - 1),Neq)/)
    ind_asq = (/(i,i=0,Neq*(Npel - 1)*Ndim,Neq*Ndim)/)

      g = 0
    force = 0.
    El_n  = 0.
    El_nn  = 0.
    Pi = 3.1415926535
    !***********************************
    !    Volume computation
    !***********************************

    ! Gauss points position
      xy = MATMUL(refElPol%N2D,Xel)
      IF (save_tau) THEN
       Xg_el = xy;
      ENDIF

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
      Bmod_nod = SQRT(Bel(:,1)**2 + Bel(:,2)**2 + Bel(:,3)**2)
    b_nod(:,1) = Bel(:,1)/Bmod_nod
    b_nod(:,2) = Bel(:,2)/Bmod_nod
    b_nod(:,3) = Bel(:,3)/Bmod_nod
#ifdef KEQUATION
    ! Toroidal magnetic field absolute value at element nodes
      b_tor_nod = ABS(Bel(:,3))
#endif

    ! Magnetic field norm and direction at Gauss points
      Bmod = MATMUL(refElPol%N2D,Bmod_nod)
      b = MATMUL(refElPol%N2D,b_nod)

#ifdef KEQUATION
    ! Toroidal magnetic field absolute value at Gauss points
      b_tor = MATMUL(refElPol%N2D,b_tor_nod)

    ! omega and q_cyl at Gauss points
      omega = MATMUL(refElPol%N2D,omegael)
      q_cyl = MATMUL(refElPol%N2D,q_cylel)
#endif
   
    ! Normalized magnetic flux at Gauss points: PSI
      Psig = MATMUL(refElPol%N2D,psiel)
  
    ! toroidal current at Gauss points
    IF (switch%ohmicsrc) THEN
         Jtor = MATMUL(refElPol%N2D,Jtorel)
    ELSE
      Jtor = 0.
    END IF

    ! Solution at Gauss points
      ueg = MATMUL(refElPol%N2D,ue)
      qeg = MATMUL(refElPol%N2D,qe)

    ! Compute diffusion at Gauss points
#ifndef KEQUATION
    CALL setLocalDiff(xy,ueg,qeg,diff_iso_vol,diff_ani_vol)
#else
    CALL setLocalDiff(xy,ueg,qeg,diff_iso_vol,diff_ani_vol,q_cyl)
#endif

    
    if (save_tau) then
       diff_nn_Vol_el = diff_iso_vol(5,5,:)
      ENDIF
    
      IF (switch%shockcp.GT.0) THEN
         auxdiffsc = MATMUL(refElPol%N2D,Mesh%scdiff_nodes(iel,:))
         DO i=1,Neq
        diff_iso_vol(i,i,:) = diff_iso_vol(i,i,:)+auxdiffsc
         END DO
      ENDIF

    ! Solution at previous time steps,at Gauss points
      DO i = 1,time%tis
         u0eg(:,:,i) = MATMUL(refElPol%N2D,u0e(:,:,i))
      END DO

    ! Physical variables at Gauss points
    CALL cons2phys(ueg,upg)

    ! Constant sources
    ! Body force at the integration points
    CALL body_force(xy(:,1),xy(:,2),force)
    
    ! Some sources to limit low density and temeprauture values
    !DO g=1, Ng2d
    !   IF (ueg(g,1) .lt. 1.e-7) force(g,1) = 1.e-7 - ueg(g,1) !Th at n = 1.00E+12 [m^(-3)]
    !   IF (upg(g,7) .lt. 6.e-4) force(g,3) = 3./2.*ueg(g,1)*(1.e-3 - upg(g,7)) !Th at Ti 0.03 eV 
    !   IF (upg(g,8) .lt. 6.e-4) force(g,4) = 3./2.*ueg(g,1)*(1.e-3 - upg(g,8)) !Th at Te 0.03 eV
    !   IF (ueg(g,5) .gt. 1.e+0) force(g,5) = 1.e+0 - ueg(g,5) !Th at nEe
    !    IF (ueg(g,3) .lt. 2.e-4) force(g,3) = 2.e-5
    !    IF (ueg(g,4) .lt. 2.e-4) force(g,4) = 2.e-5
    !END DO
   
    ! Some sources for West cases
      IF (switch%testcase .GE. 50 .AND. switch%testcase .LE. 59) THEN
      ! Compute flux surfaces and normalise them
         fluxg = MATMUL(refElPol%N2D,fluxel)
         max_flux2D = MAXVAL(phys%magnetic_flux)
         min_flux2D = MINVAL(phys%magnetic_flux)
      fluxg = (fluxg - min_flux2D)/(max_flux2D - min_flux2D)
      DO g = 1,Ng2d
            IF (ueg(g,1) .LT. th_n) THEN
          force(g,1) = th_n - ueg(g,1)
        ENDIF
        ! WEST CASE with analytical Gaussian sources on density and energies, no puff.
        IF (switch%testcase == 52) THEN
          sigma = phys%sigma_source
          x0 = 0.
               A = (phys%lscale**2)/SQRT((2*Pi*sigma**2))
               IF (fluxg(g) .LE. phys%fluxg_trunc) THEN
                  force(g,1) = phys%density_source*A*EXP(-((fluxg(g) - x0)**2)/(2*sigma**2))
#ifdef TEMPERATURE
                  force(g,3) = phys%ener_source_e*A*EXP(-((fluxg(g) - x0)**2)/(2*sigma**2))
                  force(g,4) = phys%ener_source_ee*A*EXP(-((fluxg(g) - x0)**2)/(2*sigma**2))
#endif
          ENDIF
!#ifdef NEUTRAL
!        ! WEST CASE only with analytical puff. Energy source is given by JTOR.
!        ELSE IF (switch%testcase == 54) THEN
!          ! location of the buffer differs for the mesh 3895_P4 and for the mesh 26830_P4
!          !IF (xy(g,1)*phys%lscale .gt. 2.36 .and. xy(g,2)*phys%lscale .lt. -0.69 ) THEN
!          IF (xy(g,1)*phys%lscale .gt. 2.446 .and. xy(g,1)*phys%lscale .lt. 2.59 .and. xy(g,2)*phys%lscale .gt. -0.7964 .and. xy(g,2)*phys%lscale .lt. -0.7304 ) THEN
!            force(g,5) = phys%puff
!          ENDIF
!        ELSE IF(switch%testcase .eq. 59) THEN
!          ! location of the buffer differs for the mesh 3895_P4 and for the mesh 26830_P4
!          !IF (xy(g,1)*phys%lscale .gt. 2.36 .and. xy(g,2)*phys%lscale .lt. -0.69 ) THEN
!          IF (xy(g,1)*phys%lscale .gt. 2.446 .and. xy(g,1)*phys%lscale .lt. 2.59 .and. xy(g,2)*phys%lscale .gt. -0.7964 .and. xy(g,2)*phys%lscale .lt. -0.7304 ) THEN
!           ! If it is a time initialization simulation, the puff is analytical (from param.txt) otherwise experimental
!            IF(switch%time_init) THEN
!              force(g,5) = phys%puff
!            ELSE
!              force(g,5) = phys%puff_exp(time%it+1)
!            ENDIF
!          ENDIF
!#endif
        ENDIF
#ifdef TEMPERATURE
      ! additional heating

      if (phys%heating_amplitude>1e-10) then  
        if (abs(xy(g,1)-(phys%r_axis+phys%heating_dr))<3.*abs(phys%heating_sigmar)) then
          if (abs(xy(g,2)-(phys%z_axis+phys%heating_dz))<3.*abs(phys%heating_sigmaz)) then
            force(g,phys%heating_equation) = force(g,phys%heating_equation)+phys%heating_amplitude*exp(-((xy(g,1)-(phys%r_axis+phys%heating_dr))**2)/(phys%heating_sigmar**2)) &
                                                    *exp(-((xy(g,2)-(phys%z_axis+phys%heating_dz))**2)/(phys%heating_sigmaz**2))
          endif
        endif
      endif
#endif
      END DO
    END IF
    
    ! Some sources for ITER cases    
      IF (switch%testcase .GE. 80) THEN
       ! Compute flux surfaces and normalise them
         fluxg = MATMUL(refElPol%N2D,fluxel)
         max_flux2D = MAXVAL(phys%magnetic_flux)
         min_flux2D = MINVAL(phys%magnetic_flux)
       fluxg = (fluxg - min_flux2D)/(max_flux2D - min_flux2D)
       IF (switch%testcase == 81) THEN
          sigma = phys%sigma_source
          x0 = 0.
            A = (phys%lscale**2)/SQRT((2*Pi*sigma**2))
          ! Only energy sources: density from neutral model
            IF (fluxg(g) .LE. phys%fluxg_trunc) THEN
#ifdef NEUTRAL   
#ifdef TEMPERATURE
               force(g,3) = phys%ener_source_e*A*EXP(-((fluxg(g)-x0)**2)/(2*sigma**2))
               force(g,4) = phys%ener_source_ee*A*EXP(-((fluxg(g)-x0)**2)/(2*sigma**2))
#endif
#endif
          ENDIF
       ELSE IF (switch%testcase == 82) THEN
          sigma = phys%sigma_source
            A = (phys%lscale**2)/SQRT((2*Pi*sigma**2))
          ! Only energy sources: density from neutral model
#ifdef NEUTRAL
#ifdef TEMPERATURE
            force(g,3) = phys%ener_source_e*A*EXP(-((fluxg(g) - x0)**2)/(2*sigma**2))
            force(g,4) = phys%ener_source_ee*A*EXP(-((fluxg(g) - x0)**2)/(2*sigma**2))
#endif
#endif
       ENDIF
    ENDIF 
    ! end sources

    ! Some sources for Circular cases
      IF (switch%testcase .GE. 60) THEN
      DO g=1,Ng2D
			IF (switch%testcase==61) THEN
               r =   SQRT ( (xy(g,1)*phys%lscale-geom%R0)**2 +(xy(g,2)*phys%lscale)**2 )
               IF (r .LE. 0.4) THEN
					        force(g,1) = 0. !2.838272668283863e-05
#ifdef TEMPERATURE
                  force(g,3) = 8*2.838272668283863e-05 !force(g,1)
                  force(g,4) = 8*2.838272668283863e-05
#endif
					 END IF
      ELSE IF (switch%testcase==62) THEN
#ifdef TEMPERATURE
            Pi = 3.1415926535
            sigma = 0.3
            x0 = geom%R0
            A = (phys%lscale**2)/(2*Pi*sigma**2)
            force(g,1) = 0.
               force(g,3) = 18.*A*EXP(-((xy(g,1)*phys%lscale - x0)**2 + (xy(g,2)*phys%lscale)**2)/(2*sigma**2))
            force(g,4) = force(g,3)
#endif
			END IF
			END DO
		END IF
		!! end sources


    ! Loop in 2D Gauss points
    Ngauss = Ng2d
      J11 = MATMUL(refElPol%Nxi2D,Xel(:,1))                           ! ng x 1
      J12 = MATMUL(refElPol%Nxi2D,Xel(:,2))                           ! ng x 1
      J21 = MATMUL(refElPol%Neta2D,Xel(:,1))                          ! ng x 1
      J22 = MATMUL(refElPol%Neta2D,Xel(:,2))                          ! ng x 1
    detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
    iJ11 = J22/detJ
    iJ12 = -J12/detJ
    iJ21 = -J21/detJ
    iJ22 = J11/detJ

    ! Coefficient time integration scheme
      CALL setTimeIntegrationCoefficients(ktis)
    ! Allocate temporary matrices
      ALLOCATE(Auq(Npel,Npel, neq*neq*ndim  ))
      ALLOCATE(Auu(Npel,Npel, neq*neq  ))
      ALLOCATE(rhs(Npel,Neq))
    Auq = 0.
    Auu = 0.
    rhs = 0.
    ! Loop in 2D Gauss points
    DO g = 1,NGauss

      ! Integration weight
      dvolu = refElPol%gauss_weights2D(g)*detJ(g)
      IF (switch%axisym) THEN
      	dvolu = dvolu*xy(g,1)
      END IF

      ! Check if total density is costant
      El_n  = El_n  + ueg(g,1)*2*3.1416*dvolu*phys%lscale**3
      El_nn = El_nn + ueg(g,5)*2*3.1416*dvolu*phys%lscale**3

      ! x and y derivatives of the shape functions
      Nxg = iJ11(g)*refElPol%Nxi2D(g,:) + iJ12(g)*refElPol%Neta2D(g,:)
      Nyg = iJ21(g)*refElPol%Nxi2D(g,:) + iJ22(g)*refElPol%Neta2D(g,:)

      ! Shape functions products
      Ni = refElPol%N2D(g,:)*dvolu
      NNi = tensorProduct(Ni,refElPol%N2D(g,:))                        ! Npel x Npel
      NxNi = tensorProduct(Nxg,Ni)                                     ! Npel x Npel
      NyNi = tensorProduct(Nyg,Ni)                                     ! Npel x Npel
      NNxy = b(g,1)*NxNi + b(g,2)*NyNi                                                        ! Npel x Npel
      NxyzNi = 0.
      NxyzNi(:,:,1) = NxNi
      NxyzNi(:,:,2) = NyNi                                            ! Npel x Npel x 2
      NNbb = (Nxg*b(g,1) + Nyg*b(g,2))*dvolu                             ! Npel x 1
      Nxyzg = 0.
      Nxyzg(:,1) = Nxg*dvolu
      Nxyzg(:,2) = Nyg*dvolu

      ! Divergence of b at the Gauss points
      IF (switch%axisym) THEN
        Nx_ax = Nxg + 1./xy(g,1)*refElPol%N2D(g,:)
      ELSE
        Nx_ax = Nxg
      END IF
         divbg = dot_PRODUCT(Nx_ax,b_nod(:,1)) + dot_PRODUCT(Nyg,b_nod(:,2))

      ! Diamagnetic drift !TODO: verify drift intensity in isothermal and non-isothermal cases
      driftg = 0.
      gradbmod = 0.
         gradbmod(1) = dot_PRODUCT(Nxg,Bmod_nod)
         gradbmod(2) = dot_PRODUCT(Nyg,Bmod_nod)
      bg = b(g,:)
         CALL cross_product(bg,gradbmod,driftg)
      driftg = phys%dfcoef*driftg/Bmod(g)

#ifdef KEQUATION
      ! Gradient of toroidal magnetic field on Gauss point
      gradbtor = 0.
         gradbtor(1) = dot_PRODUCT(Nxg,b_tor_nod)
         gradbtor(2) = dot_PRODUCT(Nyg,b_tor_nod)
#endif
#ifndef KEQUATION
      CALL assemblyVolumeContribution(Auq,Auu,rhs,b(g,:),Psig(g),divbg,driftg,Bmod(g),force(g,:),&
        &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upg(g,:),&
        &ueg(g,:),qeg(g,:),u0eg(g,:,:),xy(g,:),Jtor(g),Vnng)
#else
      CALL assemblyVolumeContribution(Auq,Auu,rhs,b(g,:),Psig(g),divbg,driftg,Bmod(g),b_tor(g),gradbtor,omega(g),q_cyl(g),force(g,:),&
        &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upg(g,:),&
        &ueg(g,:),qeg(g,:),u0eg(g,:,:),xy(g,:),Jtor(g),Vnng)
#endif
        
         IF (save_tau) THEN
         v_nn_Vol_el(g,:) = Vnng
         ENDIF
      
    END DO ! END loop in volume Gauss points
      CALL do_assembly(Auq,Auu,rhs,ind_ass,ind_asq,iel)
      DEALLOCATE(Auq,Auu,rhs)

  END SUBROUTINE elemental_matrices_volume

  !***************************************************
  ! Interior faces computation in 2D
  !***************************************************
#ifndef KEQUATION
  SUBROUTINE elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
#else
  SUBROUTINE elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)    
#endif    
    integer,intent(IN)        :: iel,ifa
    real*8,intent(IN)         :: Xfl(:,:)
    real*8,intent(IN)         :: Bfl(:,:), psifl(:)
    real*8,intent(IN)         :: qef(:,:)
    real*8,intent(IN)         :: uef(:,:),uf(:,:)
#ifdef KEQUATION
    real*8,intent(IN)             :: q_cylfl(:)
#endif
    real*8,intent(out)        :: diff_nn_Fac_el(:),v_nn_Fac_el(:,:),tau_save_el(:,:),xy_g_save_el(:,:)
    integer*4                 :: g,NGauss,i,indsave(Ng1d)
    real*8                    :: dline,xyDerNorm_g
    real*8                    :: ufg(Ng1d,neq),uefg(Ng1d,neq)
    real*8                    :: xyf(Ng1d,ndim)
    real*8                    :: xyDer(Ng1d,ndim)
    real*8                    :: qfg(Ng1d,neq*Ndim)
    integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
    integer*4,dimension(Npfl)  :: ind_asf,ind_ash
    real*8                    :: t_g(ndim),n_g(ndim),bn
    real*8                    :: NNif(Npfl,Npfl),Nif(Npfl),Nfbn(Npfl)
    real*8                    :: upgf(Ng1d,phys%npv)
    real*8                    :: tau(Neq,Neq),Vnng(Ndim)
    real*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ng1d,3),Bmod(Ng1d),Psig(Ng1d)
    real*8                    :: diff_iso_fac(Neq,Neq,Ng1d),diff_ani_fac(Neq,Neq,Ng1d)
    real*8                    :: auxdiffsc(Ng1d)
#ifdef KEQUATION
    real*8                    :: q_cyl(Ng1d)
#endif
    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    !***********************************
    ! Faces computations
    !***********************************
    NGauss = Ng1d

    ! Indices
      ind_fe = RESHAPE(tensorSumInt((/(i,i=1,neq)/),neq*(refElPol%face_nodes(ifa,:) - 1)),(/neq*Npfl/))
    ind_ff = (ifa - 1)*neq*Npfl + (/(i,i=1,neq*Npfl)/)
      ind_fg = RESHAPE(tensorSumInt((/(i,i=1,neq*ndim)/),neq*ndim*(refElPol%face_nodes(ifa,:) - 1)),(/neq*Npfl*ndim/))

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
      Bmod_nod = SQRT(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
    b_nod(:,1) = Bfl(:,1)/Bmod_nod
    b_nod(:,2) = Bfl(:,2)/Bmod_nod
    b_nod(:,3) = Bfl(:,3)/Bmod_nod

    ! Trace solution at face Gauss points
    IF (Mesh%flipFace(iel,ifa)) THEN
         ufg = MATMUL(refElPol%N1D,uf((/(i,i=Npfl,1,-1)/),:))
    ELSE
         ufg = MATMUL(refElPol%N1D,uf)
    END IF


    ! Gauss points position and derivatives
      xyf = MATMUL(refElPol%N1D,Xfl)
      xyDer = MATMUL(refElPol%Nxi1D,Xfl)

    ! Magnetic field norm and direction at Gauss points
    Bmod = matmul(refElPol%N1D,Bmod_nod)
    b = matmul(refElPol%N1D,b_nod)

#ifdef KEQUATION
    ! q_cyl at Gauss points
    q_cyl = matmul(refElPol%N1D,q_cylfl)
#endif
   
    ! Normalaized magnetic flux at Gauss points: PSI
      Psig = MATMUL(refElPol%N1d,psifl)

    ! Element solution at face Gauss points
      uefg = MATMUL(refElPol%N1D,uef)
    ! Gradient solution at face gauss points
      qfg = MATMUL(refElPol%N1D,qef)

    ! Compute diffusion at faces Gauss points
#ifndef KEQUATION
    CALL setLocalDiff(xyf,uefg,qfg,diff_iso_fac,diff_ani_fac)
#else
    CALL setLocalDiff(xyf,uefg,qfg,diff_iso_fac,diff_ani_fac,q_cyl)
#endif
    if (save_tau) then
       indsave = (ifa - 1)*Ngauss + (/(i,i=1,Ngauss)/)
       diff_nn_Fac_el(indsave) = diff_iso_fac(5,5,:)
      END IF

      IF (switch%shockcp.GT.0) THEN
         auxdiffsc = MATMUL(refElPol%N1D,Mesh%scdiff_nodes(iel,refElPol%face_nodes(ifa,:)))
         DO i=1,Neq
        diff_iso_fac(i,i,:) = diff_iso_fac(i,i,:)+auxdiffsc
         END DO
      ENDIF

    ! Physical variables at face Gauss points
    CALL cons2phys(ufg,upgf)

    ! Loop in 1D Gauss points
    DO g = 1,NGauss

      ! Calculate the integration weight
         xyDerNorm_g = NORM2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
      IF (switch%axisym) THEN
        dline = dline*xyf(g,1)
      END IF

      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]

      ! Shape functions products
         bn = dot_PRODUCT(b(g,1:2),n_g)
      NNif = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Nif = refElPol%N1D(g,:)*dline
      Nfbn = bn*refElPol%N1D(g,:)*dline

      ! Compute the stabilization term
      tau = 0.
      IF (numer%stab == 1) THEN
        ! Constant stabilization
        DO i = 1,Neq
          tau(i,i) = numer%tau(i)
        END DO
      ELSE
        ! Non constant stabilization
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
#ifndef KEQUATION
          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),qfg(g,:),b(g,:),n_g,iel,ifa,0.,xyf(g,:),tau)
#else
          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),qfg(g,:),b(g,:),n_g,iel,ifa,0.,xyf(g,:),q_cyl(g),tau)
#endif
        ELSE
          CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g,xyf(g,:),0.,iel,tau)
        ENDIF
      END IF

! Assembly local contributions
#ifdef DKLINEARIZED
      CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),q_cyl(g),xyf(g,:),&
      n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,Vnng)
#else
      CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,Vnng)
#endif

         IF (save_tau) THEN
        DO i = 1,Neq
          tau_save_el((ifa - 1)*Ngauss + g,i) = tau(i,i)
        END DO
        v_nn_Fac_el((ifa -1)*Ngauss + g,:) = Vnng
        xy_g_save_el((ifa - 1)*Ngauss + g,:) = xyf(g,:)
         ENDIF
      

    END DO ! Gauss points
!stop


  END SUBROUTINE elemental_matrices_faces_int

  !***************************************************
  ! Exterior faces computation in 2D
  !***************************************************
#ifndef KEQUATION
  SUBROUTINE elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
#else
  SUBROUTINE elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,Bfl,psifl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)    
#endif    
    integer,intent(IN)        :: iel,ifa
    real*8,intent(IN)         :: Xfl(:,:)
    real*8,intent(IN)         :: Bfl(:,:), psifl(:)
    logical,intent(IN)        :: isdir
    real*8,intent(IN)         :: qef(:,:)
    real*8,intent(INOUT)      :: uef(:,:),uf(:,:)
#ifdef KEQUATION
    real*8,intent(IN)             :: q_cylfl(:)
#endif
    real*8,intent(out)        :: diff_nn_Fac_el(:),v_nn_Fac_el(:,:),tau_save_el(:,:),xy_g_save_el(:,:)
    integer*4                 :: g,NGauss,i,indsave(Ng1d)
    real*8                    :: dline,xyDerNorm_g
    real*8                    :: ufg(Ng1d,neq),uefg(Ng1d,neq)
    real*8                    :: xyf(Ng1d,ndim)
    real*8                    :: xyDer(Ng1d,ndim)
    real*8                    :: qfg(Ng1d,neq*Ndim)
    integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
    integer*4,dimension(Npfl)  :: ind_asf,ind_ash
    real                      :: isext
    real*8                    :: t_g(ndim),n_g(ndim),bn
    real*8                    :: NNif(Npfl,Npfl),Nif(Npfl),Nfbn(Npfl)
    real*8                    :: tau(Neq,Neq)
    real*8                    :: upgf(Ng1d,phys%npv)
    real*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ng1d,3),Bmod(Ng1d), Psig(Ng1d)
    real*8                    :: diff_iso_fac(Neq,Neq,Ng1d),diff_ani_fac(Neq,Neq,Ng1d)
    real*8                    :: auxdiffsc(Ng1d)
    real*8                    :: Vnng(Ndim)
#ifdef KEQUATION
    real*8                    :: q_cyl(Ng1d)
#endif	
    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    !***********************************
    ! Faces computations
    !***********************************
    NGauss = Ng1d

    ! Indices
      ind_fe = RESHAPE(tensorSumInt((/(i,i=1,neq)/),neq*(refElPol%face_nodes(ifa,:) - 1)),(/neq*Npfl/))
    ind_ff = (ifa - 1)*neq*Npfl + (/(i,i=1,neq*Npfl)/)
      ind_fg = RESHAPE(tensorSumInt((/(i,i=1,neq*ndim)/),neq*ndim*(refElPol%face_nodes(ifa,:) - 1)),(/neq*Npfl*ndim/))

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
      Bmod_nod = SQRT(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
    b_nod(:,1) = Bfl(:,1)/Bmod_nod
    b_nod(:,2) = Bfl(:,2)/Bmod_nod
    b_nod(:,3) = Bfl(:,3)/Bmod_nod
    ! Magnetic field norm and direction at Gauss points
      Bmod = MATMUL(refElPol%N1D,Bmod_nod)
      b = MATMUL(refElPol%N1D,b_nod)

#ifdef KEQUATION
    ! q_cyl at Gauss points
    q_cyl = matmul(refElPol%N1D,q_cylfl)
#endif
    ! Normalaized magnetic flux at Gauss points: PSI
      Psig = MATMUL(refElPol%N1D,psifl)

    ! Trace solution at face Gauss points
      xyf = MATMUL(refElPol%N1D,Xfl)
    IF (isdir) THEN
      CALL analytical_solution(iel,xyf(:,1),xyf(:,2),ufg)
    ELSE
#ifdef PARALL
      IF (Mesh%flipFace(iel,ifa)) THEN
        uf = uf((/(i,i=Npfl,1,-1)/),:)
      ENDIF
      ! TODO: VERIFY IF I NEED TO FLIP ALSO xyf,b and Bmod in this case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#endif
         ufg = MATMUL(refElPol%N1D,uf)
    END IF

    ! Element solution at face Gauss points
      uefg = MATMUL(refElPol%N1D,uef)
    ! Gradient solution at face gauss points
      qfg = MATMUL(refElPol%N1D,qef)

    ! Compute diffusion at faces Gauss points
#ifndef KEQUATION
    CALL setLocalDiff(xyf,uefg,qfg,diff_iso_fac,diff_ani_fac)
#else
    CALL setLocalDiff(xyf,uefg,qfg,diff_iso_fac,diff_ani_fac,q_cyl)
#endif
    if (save_tau) then
       indsave = (ifa -1)*Ngauss + (/(i,i=1,Ngauss)/)
       diff_nn_Fac_el(indsave) = diff_iso_fac(5,5,:)
      END IF
	
      IF (switch%shockcp.GT.0) THEN
         auxdiffsc = MATMUL(refElPol%N1D,Mesh%scdiff_nodes(iel,refElPol%face_nodes(ifa,:)))
         DO i=1,Neq
        diff_iso_fac(i,i,:) = diff_iso_fac(i,i,:)+auxdiffsc
         END DO
      ENDIF

    ! Physical variables at face Gauss points
    CALL cons2phys(ufg,upgf)

    ! Shape function derivatives at Gauss points
      xyDer = MATMUL(refElPol%Nxi1D,Xfl)

    ! Loop in 1D Gauss points
    DO g = 1,NGauss

      ! Calculate the integration weight
         xyDerNorm_g = NORM2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
      IF (switch%axisym) THEN
        dline = dline*xyf(g,1)
      END IF
      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]

      ! Compute the stabilization term
      isext = 1.
#ifdef PARALL
         IF (Mesh%boundaryFlag(Mesh%F(iel,ifa) - Mesh%Nintfaces) .EQ. 0) THEN
        isext = 0.
      END IF
#endif
      tau = 0.
      IF (numer%stab == 1) THEN
        ! Constant stabilization
        DO i = 1,Neq
          tau(i,i) = numer%tau(i)
        END DO
      ELSE
        ! Non constant stabilization
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
#ifndef KEQUATION
          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),qfg(g,:),b(g,:),n_g,iel,ifa,isext,xyf(g,:),tau)
#else
          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),qfg(g,:),b(g,:),n_g,iel,ifa,isext,xyf(g,:),q_cyl(g),tau)
#endif
        ELSE
          CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g,xyf(g,:),isext,iel,tau)
        ENDIF
      END IF

      ! Shape functions products
         bn = dot_PRODUCT(b(g,1:2),n_g)
      NNif = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Nfbn = bn*refElPol%N1D(g,:)*dline
      Nif = refElPol%N1D(g,:)*dline

      !call displayMatrix(NNif)
      !call displayVector(Nif)
      !call displayMatrix(uf)
      !call displayMatrix(ufg)
      !stop

      ! Assembly local contributions
#ifdef PARALL
         IF (Mesh%boundaryFlag(Mesh%F(iel,ifa) - Mesh%Nintfaces) .EQ. 0) THEN
        ! Ghost face: assembly it as interior
#ifndef DKLINEARIZED
        CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),&
          n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau)

      ELSE
        CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),&
          n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau)
      ENDIF
#else

        CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),q_cyl(g),xyf(g,:),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau)
      ELSE
        CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),q_cyl(g),xyf(g,:),&
          n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau)
      ENDIF
#endif

#else
#ifndef DKLINEARIZED
      CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,Vnng)
#else
      CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),q_cyl(g),xyf(g,:),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,Vnng)
#endif

#endif
      if (save_tau) then
        DO i = 1,Neq
          tau_save_el((ifa - 1)*Ngauss + g,i) = tau(i,i)
        END DO
        v_nn_Fac_el((ifa - 1)*Ngauss + g,:) = Vnng
        xy_g_save_el((ifa - 1)*Ngauss + g,:) = xyf(g,:)
         ENDIF

    END DO ! Gauss points

  END SUBROUTINE elemental_matrices_faces_ext

  !*******************************************
  !           AUXILIARY ROUTINES
  !*******************************************

  !*****************************************
  ! Set permutations for flipping faces
  !****************************************
  SUBROUTINE set_permutations(n,m,perm)
      INTEGER,INTENT(IN)  :: n,m
      INTEGER,INTENT(OUT) :: perm(:)
      INTEGER              :: i
      INTEGER              :: temp(m,n/m),templr(m,n/m)

      IF (MOD(n,m) .NE. 0) THEN
      WRITE (6,*) 'Error! n must be a multiple of m'
      STOP
    END IF

    templr = 0
      temp = RESHAPE((/(i,i=1,n)/),(/m,n/m/))
    DO i = 1,n/m
      templr(:,i) = temp(:,n/m - i + 1)
    END DO
      perm = RESHAPE(templr,(/n/))
  END SUBROUTINE set_permutations

#endif
!TOR3D

  !********************************************************************************************
  !
  !
  !                             ROUTINES FOR 2D AND 3D COMPUTATIONS
  !
  !
  !********************************************************************************************
  SUBROUTINE setTimeIntegrationCoefficients(ktis)
      REAL*8,INTENT(out) :: ktis(:)
      INTEGER :: it

    ktis = 0.

      IF (time%ik .LT. time%tis) THEN
      it = time%ik
      ELSE
      it = time%tis
      END IF

    SELECT CASE (it)
      CASE (1)
      ktis(1) = 1.
      ktis(2) = 1.
      CASE (2)
      ktis(1) = 1.5
      ktis(2) = 2
      ktis(3) = -0.5
      CASE (3)
      ktis(1) = 11./6.
      ktis(2) = 3.
      ktis(3) = -1.5
      ktis(4) = 1./3.
      CASE (4)
      ktis(1) = 25./12.
      ktis(2) = 4.
      ktis(3) = -3.
      ktis(4) = 4./3.
      ktis(4) = -0.25
      CASE (5)
      ktis(1) = 137./60.
      ktis(2) = 5.
      ktis(3) = -5.
      ktis(4) = 10./3.
      ktis(5) = -1.25
      ktis(6) = 0.2
      CASE (6)
      ktis(1) = 147./60.
      ktis(2) = 6.
      ktis(3) = -7.5
      ktis(4) = 20./3.
      ktis(5) = -3.75
      ktis(6) = 1.2
      ktis(7) = -1.6
      CASE default
         WRITE (6,*) 'Formula not available'
         STOP
    END SELECT
  END SUBROUTINE setTimeIntegrationCoefficients

  !********************************************************************
  !
  !         ASSEMBLY VOLUME CONTRIBUTION
  !
  !********************************************************************
#ifndef KEQUATION
  SUBROUTINE assemblyVolumeContribution(Auq,Auu,rhs,b3,psi,divb,drift,Bmod,f,&
      &ktis,diffiso,diffani,Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upe,ue,qe,u0e,xy,Jtor,Vnng)
#else
  SUBROUTINE assemblyVolumeContribution(Auq,Auu,rhs,b3,psi,divb,drift,Bmod,btor,gradBtor,omega,q_cyl,f,&
    &ktis,diffiso,diffani,Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upe,ue,qe,u0e,xy,Jtor,Vnng)
#endif
        REAL*8,INTENT(inout)      :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)
        REAL*8,INTENT(IN)         :: b3(:),psi,divb,drift(:),f(:),ktis(:),Bmod
#ifdef KEQUATION
    real*8,intent(IN)         :: btor,gradBtor(:), omega, q_cyl
#ifdef DKLINEARIZED
    real*8                    :: ddk_dU(Neq), ddk_dU_U
    real*8                    :: gradddk(Ndim)
#endif
#endif
    real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
    real*8,intent(IN)         :: Ni(:),NNi(:,:),Nxyzg(:,:),NNxy(:,:),NxyzNi(:,:,:),NNbb(:)
    real*8,intent(IN)         :: upe(:),ue(:),xy(:),Jtor
    real*8,intent(INOUT)      :: u0e(:,:),Vnng(:)
    real*8,intent(IN)         :: qe(:)
    integer*4                 :: i,j,k,iord,ii,alpha,beta,z
    integer*4,dimension(Npel) :: ind_i,ind_j,ind_k
    real*8,dimension(neq,neq) :: A
    real*8,dimension(neq,Ndim):: APinch
    real*8                    :: kcoeff
    real*8                    :: Qpr(Ndim,Neq),exb(3),bb(3)
    real*8                    :: W2(Neq),dW2_dU(Neq,Neq),QdW2(Ndim,Neq)
    real*8                    :: qq(3,Neq),b(Ndim)
    real*8                    :: grad_n(3),gradpar_n
    real*8                    :: auxvec(Neq)
#ifdef TEMPERATURE
    real*8,dimension(neq,neq) :: GG
    real*8                    :: Telect
    real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
    real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
    real*8                    :: W,dW_dU(Neq,Neq),s,ds_dU(Neq),Zet(ndim,Neq)
    real*8                    :: Sohmic,dSohmic_dU(Neq) ! Ohmic heating
    real*8                    :: W3(Neq),dW3_dU(Neq,Neq),QdW3(Ndim,Neq)
    real*8                    :: W4(Neq),dW4_dU(Neq,Neq),QdW4(Ndim,Neq)
#endif
#ifdef NEUTRAL
#ifdef KEQUATION
        REAL*8                    :: gamma_I,ce, dissip,r
        REAL*8                    :: ddissip_du(Neq)
#endif
    real*8                    :: Ax(Neq,Neq),Ay(Neq,Neq)
    real*8                    :: niz,nrec,fGammacx,fGammarec
    real*8                    :: dniz_dU(Neq),dnrec_dU(Neq),dfGammacx_dU(Neq),dfGammarec_dU(Neq)
#ifdef TEMPERATURE
        REAL*8                    :: sigmaviz,sigmavrec,sigmavcx,Tloss,Tlossrec,fEiiz,fEirec,fEicx
    !amjuel radiation losses 
    real*8                    :: sigmavEiz,sigmavErec
    real*8                    :: dsigmavEiz_dU(Neq),dsigmavErec_dU(Neq)
        REAL*8                    :: dsigmaviz_dU(Neq),dsigmavrec_dU(Neq),dsigmavcx_dU(Neq),dTloss_dU(Neq),dTlossrec_dU(Neq)
        REAL*8                    :: dfEiiz_dU(Neq),dfEirec_dU(Neq),dfEicx_dU(Neq)
#ifdef DNNLINEARIZED
    real*8                    :: Dnn_dU(Neq), Dnn_dU_U
    real*8                    :: gradDnn(Ndim)
#endif
#ifdef NEUTRALP
        REAL*8                    :: Dnn,Dpn,Alphanp,Betanp,GammaLim,Gammaredpn,Tmin
        REAL*8                    :: Anp(Neq),Vpn(Neq),dVpn_dU(Neq,Neq),dDpn_dU(Neq),gmpn(Ndim),gmipn(Ndim),Taupn(Ndim,Neq)
#endif
#endif
    real*8                    :: Sn(Neq,Neq),Sn0(Neq)
#endif




        REAL*8 :: kmult(SIZE(Auq,1),SIZE(Auq,2))



    b = b3(1:Ndim)

    bb = 0.
    bb = b3

    ! Jacobian for convection term
    CALL jacobianMatrices(ue,A)
    
    ! Jacobian for pinch term
    CALL computePinch(b,psi,APinch)

    ! Compute Q^T^(k-1)
        Qpr = RESHAPE(qe,(/Ndim,Neq/))

    ! Split diffusion matrices/vectors for the momentum equation
    CALL compute_W2(ue,W2,diffiso(1,1),diffiso(2,2))
    CALL compute_dW2_dU(ue,dW2_dU,diffiso(1,1),diffiso(2,2))
        QdW2 = MATMUL(Qpr,dW2_dU)

    qq = 0.
    qq(1:Ndim,:) = Qpr

    ! Perpendicular gradient of density
    grad_n=qq(:,1)
        gradpar_n = dot_PRODUCT(grad_n,b3)

#ifdef TEMPERATURE
    ! Jacobian for the curvature term
    CALL GimpMatrix(ue,divb,GG)

    ! Compute V(U^(k-1))
        CALL computeVi(ue,Vveci)
        CALL computeVe(ue,Vvece)

    ! Compute dV_dU (k-1)
        CALL compute_dV_dUi(ue,dV_dUi)
        CALL compute_dV_dUe(ue,dV_dUe)

    ! Compute Alpha(U^(k-1))
    Alphai = computeAlphai(ue)
    Alphae = computeAlphae(ue)

    ! Compute dAlpha/dU^(k-1)
        CALL compute_dAlpha_dUi(ue,dAlpha_dUi)
        CALL compute_dAlpha_dUe(ue,dAlpha_dUe)

        gmi = dot_PRODUCT(MATMUL(Qpr,Vveci),b)    ! scalar
        gme = dot_PRODUCT(MATMUL(Qpr,Vvece),b)    ! scalar
        Taui = MATMUL(Qpr,dV_dUi)                 ! Ndim x Neq
        Taue = MATMUL(Qpr,dV_dUe)                 ! Ndim x Neq

    ! Parallel current term
    ! Compute W(U^(k-1))
        CALL compute_W(ue,W)
    ! Compute dW_dU(U^(k-1))
        CALL compute_dW_dU(ue,dW_dU)

    ! Split diffusion matrices/vectors for the energies equations
    CALL compute_W3(ue,W3,diffiso(1,1),diffiso(2,2),diffiso(3,3))
    CALL compute_dW3_dU(ue,dW3_dU,diffiso(1,1),diffiso(2,2),diffiso(3,3))
        QdW3 = MATMUL(Qpr,dW3_dU)

    CALL compute_W4(ue,W4,diffiso(1,1),diffiso(4,4))
    CALL compute_dW4_dU(ue,dW4_dU,diffiso(1,1),diffiso(4,4))
        QdW4 = MATMUL(Qpr,dW4_dU)

    ! Temperature exchange terms
    ! s(U^(k-1))
        CALL compute_S(ue,s)
    ! ds_du(U^(k-1))
        CALL compute_dS_dU(ue,ds_dU)

    !Ohmic Heating
    IF (switch%ohmicsrc) THEN
      !Compute Sohmic(U^(k-1))
           CALL compute_Sohmic(ue,Sohmic)
      !Compute dSohmic_dU(U^(k-1))
           CALL compute_dSohmic_dU(ue,dSohmic_dU)
    ENDIF

        Zet = MATMUL(Qpr,dW_dU)       ! Ndim x Neq
    
#ifdef NEUTRALP
    ! Compute Vpn(U^(k-1))
    CALL computeVpn(ue,Vpn)
	   ! Compute dVpn_dU(U^(k-1))
	   CALL compute_dVpn_dU(ue,dVpn_dU)
        gmpn = MATMUL(Qpr,Vpn)                      ! Ndim x 1
        Taupn = MATMUL(Qpr,dVpn_dU)                 ! Ndim x Neq
	   ! Compute Dpn(U^(k-1))
    CALL computeDpn(ue,Qpr,Vpn,Dpn)
    ! Compute dDpn_dU(U^(k-1))
    CALL compute_dDpn_dU(ue,Qpr,Vpn,dDpn_dU) 
    ! Reduce Grad Pn for low collision regime 
    ! Threshold set at 0.5xGradPn for Ti = 0.2 eV 
    Gammaredpn = 1.
    Tmin = 0.2/simpar%refval_temperature
        IF (Tmin/upe(7) .LE. 1.) Gammaredpn = Gammaredpn*Tmin/upe(7)
    Dnn = Gammaredpn*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upe(7)*Dpn 
    ! Compute Gammaredpn(U^(k-1))
    !CALL computeGammared(ue,Gammaredpn) 
    !gmipn = matmul(Qpr,Vveci) 
    !CALL computeGammaLim(ue,Qpr,Vpn,GammaLim)
    ! Set Grad Ti = 0. for low collision regime 
    ! (back to diffusion equation for neutral density)
    !CALL computeAlphaCoeff(ue,Qpr,Vpn,Alphanp)  
    !CALL computeBetaCoeff(ue,Qpr,Vpn,Betanp)  
    !Dpn = Alphanp*Dpn
    !dDpn_dU = Alphanp*dDpn_dU
    !Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upe(7)*Dpn 
    !IF (Dnn .gt. phys%diff_nn) Dnn = phys%diff_nn
    !IF (Dpn .gt. phys%diff_nn) THEN
    !   Dpn = Alphanp*Dpn !0.
    !   dDpn_dU = Alphanp*dDpn_dU !0.
    !   Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upe(7)*phys%diff_nn
    !   END IF
    ! Set Gamma Convective = cs_n*n_n for low collision regime 
    !IF (Dpn .gt. phys%diff_nn) THEN
    !   Dpn = 0.
    !   dDpn_dU = 0.
    !   CALL jacobianMatricesNP(ue,Anp)
    !ELSE
    !   Anp = 0.
    !END IF
#endif
#endif
#ifdef KEQUATION
        IF ((switch%testcase .GE. 50) .AND.(switch%testcase .LE. 59)) THEN
      r = xy(1)
        ELSEIF ((switch%testcase .GE. 60) .AND.(switch%testcase .LE. 69)) THEN
      r = xy(1) + geom%R0/simpar%refval_length
    endif
    call compute_gamma_I(ue,qq,btor,gradBtor,r,gamma_I)
    call compute_ce(ue,qq,btor,gradBtor,r,omega,q_cyl,ce)
    call compute_dissip(ue,dissip)
    call compute_ddissip_du(ue,ddissip_du)
    if ((ue(6)<1.e-20) .or. (ue(1)<1e-20) .or.(ue(3)<1.e-20) .or. (ue(4)<1e-20)) then
      dissip =  abs(gamma_I)*dissip/phys%k_max
      ddissip_du = abs(gamma_I)
      gamma_I = 0.
        ELSEIF (ue(6)>phys%k_max) THEN
           dissip =  -1.*ABS(gamma_I)*dissip/phys%k_max
           ddissip_du = -1.*ABS(gamma_I)
      gamma_I = 0.
        ELSE
           IF (gamma_I>0) THEN
        dissip = ce*dissip
        ddissip_du = ce*ddissip_du
      else
        dissip = abs(gamma_I)*dissip/phys%k_max
        ddissip_du = abs(gamma_I)
        gamma_I = 0.
      endif
    endif
#ifdef DKLINEARIZED
    call compute_ddk_dU(ue,xy,q_cyl,ddk_dU)

    ddk_dU_u = dot_product(ddk_dU,ue)
#endif
#endif



#ifdef NEUTRAL
    !Neutral Source Terms needed in the plasma and neutral density equations
        CALL compute_niz(ue,niz)
        CALL compute_nrec(ue,nrec)
        CALL compute_dniz_dU(ue,dniz_dU)
        CALL compute_dnrec_dU(ue,dnrec_dU)
#ifdef TEMPERATURE
        CALL compute_sigmaviz(ue,sigmaviz)
        CALL compute_sigmavrec(ue,sigmavrec)
        CALL compute_dsigmaviz_dU(ue,dsigmaviz_dU)
        CALL compute_dsigmavrec_dU(ue,dsigmavrec_dU)
#endif
    !Neutral Source Terms needed in the plasma momentum equation
        CALL compute_fGammacx(ue,fGammacx)
        CALL compute_dfGammacx_dU(ue,dfGammacx_dU)
        CALL compute_fGammarec(ue,fGammarec)
        CALL compute_dfGammarec_dU(ue,dfGammarec_dU)
#ifdef TEMPERATURE
        CALL compute_sigmavcx(ue,sigmavcx)
        CALL compute_dsigmavcx_dU(ue,dsigmavcx_dU)
    !Neutral Source Terms needed in the ion energy equation
        CALL compute_fEiiz(ue,fEiiz)
        CALL compute_dfEiiz_dU(ue,dfEiiz_dU)
        CALL compute_fEirec(ue,fEirec)
        CALL compute_dfEirec_dU(ue,dfEirec_dU)
        CALL compute_fEicx(ue,fEicx)
        CALL compute_dfEicx_dU(ue,dfEicx_dU)
    !Neutral Source Terms needed in the electron energy equation
        CALL compute_Tloss(ue,Tloss)
        CALL compute_dTloss_dU(ue,dTloss_dU)
        CALL compute_Tlossrec(ue,Tlossrec)
        CALL compute_dTlossrec_dU(ue,dTlossrec_dU)
    !Amjuel energy losses
    call compute_sigmavEiz(ue,sigmavEiz)
    call compute_sigmavErec(ue,sigmavErec)
    call compute_dsigmavEiz_dU(ue,dsigmavEiz_dU)
    call compute_dsigmavErec_dU(ue,dsigmavErec_dU)
#ifdef DNNLINEARIZED
        CALL compute_Dnn_dU(ue,Dnn_dU)

        Dnn_dU_u = dot_PRODUCT(Dnn_dU,Ue)
#endif

#endif

    !Assembly the matrix for neutral sources
#ifdef TEMPERATURE
    call assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
      &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
      &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,Tloss,dTloss_dU,Tlossrec,dTlossrec_dU,sigmavEiz,dsigmavEiz_dU,sigmavErec,dsigmavErec_dU,Sn,Sn0)
#else
        CALL assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,Sn,Sn0)
#endif
#endif
!NEUTRAL

    ! Assembly local matrix
    ! Loop in equations
    DO i = 1,Neq
      ! ind_i = i + ind_ass
           IF (.NOT. switch%steady) THEN
        ! Time derivative contribution
#ifdef VORTICITY
        ! In the vorticity model I don't assemble the mass matrix for the potential equation
              IF (i .NE. 4) THEN
#endif
          z = i+(i-1)*Neq
          Auu(:,:,z)= Auu(:,:,z)+ ktis(1)*NNi/time%dt
#ifdef VORTICITY
        END IF
#endif
      END IF
#ifndef TEMPERATURE
      IF (i == 2) THEN
        ! Curvature contribution (isothermal)
        z = i+(i-2)*Neq
              IF (switch%logrho) THEN
          Auu(:,:,z) = Auu(:,:,z) - phys%a*divb*upe(1)*NNi
          rhs(:,i)=rhs(:,i)+Ni*phys%a*divb*upe(1)*(1-ue(1))
              ELSE
          Auu(:,:,z) = Auu(:,:,z) - phys%a*divb*NNi
              ENDIF

        ! split diffusion momentum equation (LU) (isothermal)
        DO j = 1,Neq
           z = i+(j-1)*Neq
           DO k = 1,Ndim
              Auu(:,:,z) =Auu(:,:,z) + (NxyzNi(:,:,k)*QdW2(k,j))
           END DO
                 Auu(:,:,z) = Auu(:,:,z) - (dot_PRODUCT(QdW2(:,j),b))*NNxy
        END DO
      END IF
#ifdef VORTICITY
           IF (switch%driftdia .AND. i.NE.4) THEN
#else
        IF (switch%driftdia) THEN
#endif
          ! B x GradB drift (isothermal)
          DO k = 1,Ndim
            z = i+(i-1)*Neq
                    Auu(:,:,z)= Auu(:,:,z) +TRANSPOSE(NxyzNi(:,:,k))*drift(k)
          END DO
        END IF
#else
        IF (i == 2) THEN
          DO j = 1,Neq
            z = i+(j-1)*Neq
            ! Curvature contribution (non-isothermal)
            Auu(:,:,z) = Auu(:,:,z) - GG(i,j)*NNi

            ! split diffusion momentum equation (LU) (non-isothermal)
            DO k = 1,Ndim
               Auu(:,:,z) =Auu(:,:,z) + (NxyzNi(:,:,k)*QdW2(k,j))
            END DO
                    Auu(:,:,z) = Auu(:,:,z) - (dot_PRODUCT(QdW2(:,j),b))*NNxy
          END DO
        END IF

        IF (switch%driftdia) THEN
          Telect = upe(8)
          ! B x GradB drift (non-isothermal)
          DO k = 1,Ndim
            z = i+(i-1)*Neq
                    Auu(:,:,z)= Auu(:,:,z) + Telect*TRANSPOSE(NxyzNi(:,:,k))*drift(k)
          END DO
        END IF

        ! Parallel diffusion for the temperature
        IF (i == 3) THEN
          DO j = 1,4
                    Auu(:,:,i+(j-1)*Neq) = Auu(:,:,i+(j-1)*Neq)+coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_PRODUCT(Taui(:,j),b)))*NNxy + &
                         &(dot_PRODUCT(Zet(:,j),b) + ds_dU(j))*NNi
            DO k = 1,Ndim
              z = i+(k-1)*Neq+(j-1)*Neq*Ndim
              Auq(:,:,z) = Auq(:,:,z)+ coefi*Alphai*Vveci(j)*b(k)*NNxy
              IF (j == 4) THEN
                Auq(:,:,z) = Auq(:,:,z)+W*NNi*b(k)
              END IF
            END DO
            ! split diffusion electron energy equation (LU)
            z = i+(j-1)*Neq
            DO k = 1,Ndim
              Auu(:,:,z) = Auu(:,:,z) + (NxyzNi(:,:,k)*QdW3(k,j))
            END DO
                    Auu(:,:,z) = Auu(:,:,z) - (dot_PRODUCT(QdW3(:,j),b))*NNxy
          END DO
                 rhs(:,i) = rhs(:,i) + coefi*Alphai*(dot_PRODUCT(MATMUL(TRANSPOSE(Taui),b),ue))*NNbb + s*Ni
        ELSEIF (i == 4) THEN
          DO j = 1,4
            z = i+(j-1)*Neq
                    Auu(:,:,z)=Auu(:,:,z)+coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_PRODUCT(Taue(:,j),b)))*NNxy - (dot_PRODUCT(Zet(:,j),b) + ds_dU(j))*NNi
            IF (switch%ohmicsrc) THEN
              Auu(:,:,z) = Auu(:,:,z) - dSohmic_dU(j)*(Jtor**2)*NNi
            ENDIF
            DO k = 1,Ndim
              z = i+(k-1)*Neq+(j-1)*Neq*Ndim
              Auq(:,:,z)=Auq(:,:,z)+coefe*Alphae*Vvece(j)*b(k)*NNxy
              IF (j == 4) THEN
                Auq(:,:,z)=Auq(:,:,z)- W*NNi*b(k)
              END IF
            END DO
            ! split diffusion electron energy equation (LU)
            z = i+(j-1)*Neq
            DO k = 1,Ndim
              Auu(:,:,z) = Auu(:,:,z) + (NxyzNi(:,:,k)*QdW4(k,j))
            END DO
                    Auu(:,:,z) = Auu(:,:,z) - (dot_PRODUCT(QdW4(:,j),b))*NNxy
          END DO
                 rhs(:,i) = rhs(:,i)+coefe*Alphae*(dot_PRODUCT(MATMUL(TRANSPOSE(Taue),b),ue))*NNbb - s*Ni
          IF (switch%ohmicsrc) THEN
            rhs(:,i) = rhs(:,i) + Sohmic*(Jtor**2)*Ni
          ENDIF
#ifdef DNNLINEARIZED
          ELSEIF (i == 5) THEN
                 DO j = 1,5
              z = i+(j-1)*Neq
                    DO k = 1,Ndim
                !z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                Auu(:,:,z) =Auu(:,:,z) + (NxyzNi(:,:,k)*Dnn_dU(j)*Qpr(k,i))
                    ENDDO
                 ENDDO
              
            DO k = 1, Ndim
              rhs(:,i) = rhs(:,i)+Dnn_dU_U*Qpr(k,i)*Nxyzg(:,k)
                 ENDDO
#endif
#ifdef KEQUATION
        ELSEIF (i==6) THEN
          DO j=1,6
            z = i+(j-1)*Neq
                    IF (j==6) THEN
              Auu(:,:,z) = Auu(:,:, z) - (gamma_I-ddissip_du(j))*NNi 
            ENDIF
          END DO
          rhs(:,i) = rhs(:,i) + dissip*Ni
#endif
#ifdef NEUTRALP		  
		      ELSEIF (i == 5) THEN
		         DO j = 1,5
			           DO k = 1,Ndim
			              z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                    Auu(:,:,i+(j-1)*Neq) = Auu(:,:,i+(j-1)*Neq) + (Dpn*Taupn(k,j) + dDpn_dU(j)*gmpn(k))*NxyzNi(:,:,k) 
                    !Auu(:,:,i+(j-1)*Neq) = Auu(:,:,i+(j-1)*Neq) - Gammaredpn*(Dpn*Taui(k,j) + dDpn_dU(j)*gmipn(k))*NxyzNi(:,:,k) 
                    Auq(:,:,z) = Auq(:,:,z) + Dpn*Vpn(j)*NxyzNi(:,:,k)
                    !Auq(:,:,z) = Auq(:,:,z) - Gammaredpn*Dpn*Vveci(j)*NxyzNi(:,:,k)
                    IF (j == 5) THEN
                       Auq(:,:,z) = Auq(:,:,z) + Dnn*NxyzNi(:,:,k)
                    END IF
				         END DO
           END DO
           DO k = 1,Ndim
                    rhs(:,i) = rhs(:,i) + dot_PRODUCT(dDpn_dU,ue)*gmpn(k)*Nxyzg(:,k)
              !rhs(:,i) = rhs(:,i) - Gammaredpn*(Dpn*dot_product(Taui(k,:),ue) + dot_product(dDpn_dU,ue)*gmipn(k))*Nxyzg(:,k)
           END DO
#endif
		END IF
#endif
#ifdef KEQUATION
#ifdef DKLINEARIZED
    ! Contribution from linearized dk term assuming so far that Dk is the same in all plasma equations
        if (i .ne. 5) then
          DO j = 1,6
            z = i+(j-1)*Neq
            do k = 1,Ndim
              Auu(:,:,z) =Auu(:,:,z) + ddk_dU(j)*Qpr(k,i)*(NxyzNi(:,:,k)-b(k)*NNxy)
            enddo
          enddo
          DO k = 1, Ndim
            rhs(:,i) = rhs(:,i)+ddk_dU_U*Qpr(k,i)*(Nxyzg(:,k)-b(k)*NNbb)
          enddo
        endif
#endif
#endif
       
	! Convection contribution
        DO j = 1,Neq
           z = i+(j-1)*Neq
           Auu(:,:,z)= Auu(:,:,z) - A(i,j)*NNxy
#ifdef NEUTRAL
!#ifdef NEUTRALP
!          IF (i == 5) Auu(:,:,z) = Auu(:,:,z) - Anp(j)*NNxy
!#endif
          !Sources
          Auu(:,:,z) = Auu(:,:,z) + Sn(i,j)*NNi
#endif
        END DO

        ! Pinch contribution
        z = i+(i-1)*Neq
        Auu(:,:,z) =  Auu(:,:,z) - (APinch(i,1)*NxyzNi(:,:,1) + APinch(i,2)*NxyzNi(:,:,2))

#ifndef TEMPERATURE
        ! Added term for n=exp(x) change of variable
              IF (switch%logrho) THEN
                 CALL logrhojacobianVector(ue,upe,auxvec)
          rhs(:,i)=rhs(:,i)+NNbb*auxvec(i)

                 DO j = 1,Neq
                    IF (i==1 .AND. j==1) THEN
              z = i+(j-1)*Neq
                       Auu(:,:,z)=Auu(:,:,z)+upe(2)*TRANSPOSE(NNxy) ! TODO: this is a first order linearization!!!
                    ENDIF
                 END DO
              ENDIF
#endif

        DO k = 1,Ndim
        ! split diffusion contributions (LQ)
	        IF (i==2) THEN
            DO j = 1,Neq
                z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                Auq(:,:,z) = Auq(:,:,z) + W2(j)*(NxyzNi(:,:,k) -NNxy*b(k))
            END DO
#ifdef TEMPERATURE
          ELSEIF(i==3) THEN
            DO j = 1,Neq
	         z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                 Auq(:,:,z) = Auq(:,:,z) + W3(j)*(NxyzNi(:,:,k) -NNxy*b(k))
           END DO
          ELSEIF(i==4) THEN
            DO j = 1,Neq
                z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                Auq(:,:,z) = Auq(:,:,z) + W4(j)*(NxyzNi(:,:,k) -NNxy*b(k))
            END DO
#endif
           ENDIF

          ! Diagonal terms for perpendicular diffusion
          z = i+(k-1)*Neq+(i-1)*Neq*Ndim
          kmult = diffiso(i,i)*NxyzNi(:,:,k) - diffani(i,i)*NNxy*b(k)
#ifdef VORTICITY
                 IF (i==4) THEN
            kmult = kmult/ue(1)
                 ENDIF
#endif

          Auq(:,:,z)=Auq(:,:,z)+kmult

#ifndef TEMPERATURE
          ! Added term for n=exp(x) change of variable \Grad \chi **2
                 IF (switch%logrho) THEN
                    IF (i==1) THEN
              Auq(:,:,z)=Auq(:,:,z)-2*(diffiso(i,i)*grad_n(k) - diffani(i,i)*b(k)*gradpar_n)*NNi
              rhs(:,i)=rhs(:,i)-Ni*(diffiso(i,i)*grad_n(k)*grad_n(k) - diffani(i,i)*gradpar_n**2/Ndim )
                    ENDIF
                 ENDIF
#endif

#ifdef VORTICITY
                 IF (switch%bxgradb) THEN
            ! B x GradB current in the vorticity equation
            IF (i==3) THEN
              ii =1
              z = i+(ii-1)*Neq
                       IF (switch%logrho) THEN
                          Auu(:,:,z)= Auu(:,:,z) +2*upe(1)*TRANSPOSE(NxyzNi(:,:,k))*drift(k) ! TODO: first order linearization
                       ELSE
                          Auu(:,:,z)= Auu(:,:,z) +2*TRANSPOSE(NxyzNi(:,:,k))*drift(k)
            ENDIF
                    ENDIF
                 ENDIF
                 IF (switch%driftexb .AND. i .NE. 4) THEN
            ! ExB terms
            kcoeff = phys%dfcoef*numer%exbdump/Bmod
                    CALL ijk_cross_product(k,alpha,beta)
            ii = 4
            z = i+(k-1)*Neq+(ii-1)*Neq*Ndim
            Auq(:,:,z)=Auq(:,:,z)+kcoeff*(NxyzNi(:,:,alpha)*b3(beta) - NxyzNi(:,:,beta)*b3(alpha))*ue(i)
                    CALL cross_product(qq(:,ii),bb,exb)
            z = i+(i-1)*Neq
            Auu(:,:,z)=Auu(:,:,z)+kcoeff*exb(k)*NxyzNi(:,:,k)
            rhs(:,i)=rhs(:,i)+kcoeff*exb(k)*ue(i)*Nxyzg(:,k)
          ENDIF

          ! Non-diagonal terms for perpendicular diffusion
          DO ii = 1,Neq
            IF (ii == i) CYCLE ! diagonal already assembled
                    IF (ABS(diffiso(i,ii)) < 1e-12 .AND. ABS(diffani(i,ii)) < 1e-12) CYCLE
            kcoeff = 1.
            ! Non-linear correction for non-linear diffusive terms.
            ! TODO: find a smarter way to include it,avoiding if statements and model dependencies (i==3,ii==1 only holds for Isothermal+Vorticity model)

            !IF ((i == 3 .or. i == 4) .and. ii == 1) then
                    IF ((i == 3) .AND. ii == 1) THEN
              z = i+(ii-1)*Neq
              kcoeff = 1./ue(1)
              Auu(:,:,z)=Auu(:,:,z)-kcoeff**2*(diffiso(i,ii)*Qpr(k,ii)*NxyzNi(:,:,k) - diffani(i,ii)*Qpr(k,ii)*b(k)*NNxy)
              rhs(:,i)=rhs(:,i)-kcoeff*(diffiso(i,ii)*Qpr(k,ii)*Nxyzg(:,k) - diffani(i,ii)*Qpr(k,ii)*b(k)*NNbb)
            ENDIF
            z=i+(k-1)*Neq+(ii-1)*Ndim*Neq
            Auq(:,:,z)=Auq(:,:,z)+kcoeff*diffiso(i,ii)*NxyzNi(:,:,k) - kcoeff*diffani(i,ii)*NNxy*b(k)

            !write(6,*) "kcoeff",kcoeff
            !write(6,*) "i:",i,"ii:",ii, "diffiso(i,ii)",diffiso(i,ii)
            !write(6,*) "i:",i,"ii:",ii, "diffani(i,ii)",diffani(i,ii)
            !call HDF5_save_matrix(kcoeff*diffiso(i,ii)*NxyzNi(:,:,k) - kcoeff*diffani(i,ii)*NNxy*b(k),'fava')
            !stop

          END DO
#endif
        END DO ! loop in k: 1-Ndim

#ifdef VORTICITY
        ! The vorticity is the source term in the potential equation
        IF (i == 4) THEN
          j=3
          z = i+(j-1)*Neq
          Auu(:,:,z)=Auu(:,:,z) +NNi
          !rhs(:,i)=rhs(:,i)-ue(j)*Ni ! doens't work very well like this
        ENDIF
#endif

#ifdef VORTICITY
        !if (switch%testcase.eq.7 .and. switch%logrho .and. i.eq.1 .and. upe(1).gt.1) then
        !  rhs(:,i)=rhs(:,i)-100*(upe(1)-1.)*Ni
        !endif
              IF (switch%testcase.EQ.7) THEN
                 IF ( (xy(1)-geom%R0)/phys%lscale .GT. 0.4 ) THEN
            ! Implicit sources to take into account parallel losses
                    IF (switch%logrho) THEN
                       IF (i==1) THEN
                rhs(:,i)=rhs(:,i)-phys%diagsource(i)*Ni
                       ELSE IF (i==3) THEN
                j=4
                z = i+(j-1)*Neq
                Auu(:,:,z)=Auu(:,:,z) +phys%diagsource(i)*NNi
                       ENDIF
                    ELSE
                       IF (i==1) THEN
                z = i+(i-1)*Neq
                Auu(:,:,z)=Auu(:,:,z) +phys%diagsource(i)*NNi
                       ELSE IF (i==3) THEN
                j=4
                z = i+(j-1)*Neq
                Auu(:,:,z)=Auu(:,:,z) +phys%diagsource(i)*NNi
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
#endif
      END DO ! Loop in equations

      ! Assembly RHS
           IF (.NOT. switch%steady) THEN
#ifdef VORTICITY
        u0e(4,:) = 0.
#endif
        DO iord = 1,time%tis
          ! Time derivative contribution
          rhs=rhs+ktis(iord + 1)*tensorProduct(Ni,u0e(:,iord))/time%dt
        END DO
      END IF

      ! Linear body force contribution
      rhs = rhs+tensorProduct(Ni,f)
#ifdef NEUTRAL
      rhs = rhs-tensorProduct(Ni,Sn0)
#endif
!#ifdef NEUTRALP
!      rhs = rhs+tensorProduct(Ni,fth)
!#endif
    END SUBROUTINE assemblyVolumeContribution

    !********************************************************************
    !
    !         ASSEMBLY INTERIOR FACES CONTRIBUTION
    !
    !********************************************************************

#ifdef DKLINEARIZED
  SUBROUTINE assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,&
    &ind_fg,b3,Bmod,psi,q_cyl,xyf,n,diffiso,diffani,NNif,Nif,Nfbn,uf,upf,qf,tau,Vnng,ifa)
#else
    SUBROUTINE assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,&
        &ind_fg,b3,Bmod,psi,n,diffiso,diffani,NNif,Nif,Nfbn,uf,upf,qf,tau,Vnng,ifa)
#endif
      integer*4,intent(IN)      :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
      real*8,intent(IN)         :: b3(:),n(:),Bmod, psi
      real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
      real*8,intent(IN)         :: NNif(:,:),Nif(:),Nfbn(:)
      real*8,intent(IN)         :: uf(:),upf(:)
      real*8,intent(IN)         :: qf(:)
#ifdef KEQUATION
#ifdef DKLINEARIZED
      real*8,intent(IN)         :: q_cyl, xyf(:)
#endif
#endif
      real*8,optional,intent(INOUT) :: tau(:,:),Vnng(:)
      real*8                     :: kcoeff
      real*8                     :: b(Ndim)
      integer*4,optional         :: ifa
      integer*4                  :: i,j,k,ii,alpha,beta
      integer*4,dimension(size(ind_asf))  :: ind_if,ind_jf,ind_kf
      real*8,dimension(neq,neq) :: A
      real*8,dimension(neq,Ndim):: APinch
      real*8                    :: auxvec(Neq)
      real*8                    :: nn(3),qq(3,Neq),bb(3)
      real*8                    :: bn,kmult(size(ind_asf),size(ind_asf)),kmultf(size(ind_asf))
      real*8                    :: Qpr(Ndim,Neq),exb(3)
      real*8                    :: W2(Neq),dW2_dU(Neq,Neq),QdW2(Ndim,Neq)
#ifdef TEMPERATURE
      real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
      real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
      real*8                    :: W3(Neq),dW3_dU(Neq,Neq),QdW3(Ndim,Neq)
      real*8                    :: W4(Neq),dW4_dU(Neq,Neq),QdW4(Ndim,Neq)
#ifdef DNNLINEARIZED
      real*8                    :: Dnn_dU(Neq), Dnn_dU_U
      real*8                    :: gradDnn(Ndim)
#endif
#ifdef KEQUATION
#ifdef DKLINEARIZED
      real*8                    :: ddk_dU(Neq), ddk_dU_U
      real*8                    :: gradddk(Ndim)
#endif
#endif
#ifdef NEUTRALP
      real*8                    :: Dnn,Dpn,GammaLim,Alphanp,Betanp,Gammaredpn,Tmin
   	  real*8                    :: Anp(Neq),Vpn(Neq),dVpn_dU(Neq,Neq),gmpn(Ndim),gmipn(Ndim),Taupn(Ndim,Neq),dDpn_dU(Neq)
#endif
#endif

      b = b3(1:Ndim)
      bb = b3
      ! Jacobian matrices
           bn = dot_PRODUCT(b,n)
      CALL jacobianMatrices(uf,A)

      ! Jacobian for pinch term
      CALL computePinch(b,psi,APinch)  
   
      ! Compute Q^T^(k-1)
           Qpr = RESHAPE(qf,(/Ndim,Neq/))

      ! Split diffusion vector/matrix for momentum equation
      CALL compute_W2(uf,W2,diffiso(1,1),diffiso(2,2))
      CALL compute_dW2_dU(uf,dW2_dU,diffiso(1,1),diffiso(2,2))
           QdW2 = MATMUL(Qpr,dW2_dU)

      nn = 0.
      qq = 0.
      nn(1:Ndim) = n
      qq(1:Ndim,:) = Qpr
      
#ifdef TEMPERATURE
      ! Compute V(U^(k-1))
           CALL computeVi(uf,Vveci)
           CALL computeVe(uf,Vvece)

      ! Compute dV_dU (k-1)
           CALL compute_dV_dUi(uf,dV_dUi)
           CALL compute_dV_dUe(uf,dV_dUe)

      ! Split diffusion vector/matrix for the energies equations
      CALL compute_W3(uf,W3,diffiso(1,1),diffiso(2,2),diffiso(3,3))
      CALL compute_dW3_dU(uf,dW3_dU,diffiso(1,1),diffiso(2,2),diffiso(3,3))
           QdW3 = MATMUL(Qpr,dW3_dU)

      CALL compute_W4(uf,W4,diffiso(1,1),diffiso(4,4))
      CALL compute_dW4_dU(uf,dW4_dU,diffiso(1,1),diffiso(4,4))
           QdW4 = MATMUL(Qpr,dW4_dU)

      ! Compute Alpha(U^(k-1))
      Alphai = computeAlphai(uf)
      Alphae = computeAlphae(uf)

      ! Compute dAlpha/dU^(k-1)
           CALL compute_dAlpha_dUi(uf,dAlpha_dUi)
           CALL compute_dAlpha_dUe(uf,dAlpha_dUe)

           gmi = dot_PRODUCT(MATMUL(Qpr,Vveci),b)  ! scalar
           gme = dot_PRODUCT(MATMUL(Qpr,Vvece),b)
           Taui = MATMUL(Qpr,dV_dUi)      ! 2x3
           Taue = MATMUL(Qpr,dV_dUe)      ! 2x3
#ifdef DNNLINEARIZED
           CALL compute_Dnn_dU(uf,Dnn_dU)

      Dnn_dU_u = dot_product(Dnn_dU,uf)
#endif
#ifdef KEQUATION
#ifdef DKLINEARIZED
      call compute_ddk_dU(uf,xyf,q_cyl,ddk_dU)

      ddk_dU_u = dot_product(ddk_dU,uf)
#endif
#endif
#ifdef NEUTRALP
    ! Compute Vpn(U^(k-1))
    CALL computeVpn(uf,Vpn)
           gmpn = MATMUL(Qpr,Vpn)                     ! Ndim x 1
	  ! Compute dVpn-dU(U^(k-1))
	  CALL compute_dVpn_dU(uf,dVpn_dU)
           Taupn = MATMUL(Qpr,dVpn_dU)                 ! Ndim x Neq
	  ! Compute Dpn(U^(k-1))
    CALL computeDpn(uf,Qpr,Vpn,Dpn)
    ! Compute dDpn_dU(U^(k-1))
    CALL compute_dDpn_dU(uf,Qpr,Vpn,dDpn_dU)
    ! Reduce Grad Pn for low collision regime 
    ! Threshold set at 0.5xGradPn for Ti = 0.2 eV 
    Gammaredpn = 1.
    Tmin = 0.2/simpar%refval_temperature
           IF (Tmin/upf(7) .LE. 1.) Gammaredpn = Gammaredpn*Tmin/upf(7)
    Dnn = Gammaredpn*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*Dpn 
    ! Comput Gammaredpn(U^(k-1))
    !CALL computeGammared(uf,Gammaredpn)
    !gmipn = matmul(Qpr,Vveci)
    !CALL computeGammaLim(ue,Qpr,Vpn,GammaLim)
    ! Set Grad Ti = 0. for low collision regime 
    ! (back to diffusion equation for neutral density)
    !CALL computeAlphaCoeff(uf,Qpr,Vpn,Alphanp)  
    !CALL computeBetaCoeff(uf,Qpr,Vpn,Betanp)  
    !Dpn = Alphanp*Dpn
    !dDpn_dU = Alphanp*dDpn_dU
    !Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*Dpn 
    !IF (Dnn .gt. phys%diff_nn) Dnn = phys%diff_nn
    !IF (Dpn .gt. phys%diff_nn) THEN
    !   Dpn = Alphanp*Dpn !0.
    !   dDpn_dU = Alphanp*dDpn_dU !0.
    !   Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*phys%diff_nn
    !   END IF
    ! Set Gamma Convective = cs_n*n_n for low collision regime 
    !IF (Dpn .gt. phys%diff_nn) THEN
    !   Dpn = 0.
    !   dDpn_dU = 0.
    !   CALL jacobianMatricesNP(uf,Anp)
    !ELSE
    !   Anp = 0.
    !END IF
#endif
#endif

      ! Assembly local matrix
      DO i = 1,Neq
        ind_if = ind_asf + i
        DO k = 1,Ndim
          ind_kf = ind_ash + k + (i - 1)*Ndim
          kmult = NNif*(n(k)*diffiso(i,i) - bn*b(k)*diffani(i,i))

#ifdef VORTICITY
                 IF (i==4) THEN
            kmult=kmult/uf(1)
                 ENDIF
#endif
          ! Diagonal terms for perpendicular diffusion
          elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
          elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult

          ! Split diffusion contribution
          IF(i == 2) THEN
            ! Assembly LQ momentum equatuation
            j = 1
            ind_kf = ind_ash + k + (j - 1)*Ndim
             kmult = NNif*W2(j)*(n(k) - bn*b(k))
             elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
             elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
           !Assembly LU momentum equation
             DO j=1,Neq
               ind_jf = ind_asf+j
               kmult = QdW2(k,j)*NNif*(n(k)-bn*b(k))
               elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                       elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel)  = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
             END DO
#ifdef TEMPERATURE
         ELSEIF(i==3) THEN
            ! Assembly LQ ion energy
            DO j=1,Neq
              ind_kf = ind_ash + k + (j - 1)*Ndim
               kmult = NNif*W3(j)*(n(k) - bn*b(k))
               elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
               elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
            END DO
           ! Assembly LU ion energy
             DO j=1,Neq
               ind_jf = ind_asf+j
               kmult = QdW3(k,j)*NNif*(n(k)-bn*b(k))
               elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                       elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel)  = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
             END DO
          ELSEIF(i == 4) THEN
             ! Assembly LQ electron energy
             j = 1
             ind_kf = ind_ash + k + (j - 1)*Ndim
              kmult = NNif*W4(j)*(n(k) - bn*b(k))
              elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
              elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
             !Assembly LU electron energy
              DO j=1,Neq
                ind_jf = ind_asf+j
                kmult = QdW4(k,j)*NNif*(n(k)-bn*b(k))
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                       elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel)  = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
              END DO
#endif
          ENDIF
#ifdef VORTICITY

                 IF (switch%driftexb .AND. i .NE. 4) THEN
            ! ExB terms
            kcoeff = phys%dfcoef*numer%exbdump/Bmod
            ii = 4
                    CALL ijk_cross_product(k,alpha,beta)
            ind_kf = ind_ash + k + (ii - 1)*Ndim
            kmult = kcoeff*NNif*(nn(alpha)*b3(beta) - nn(beta)*b3(alpha))*uf(i)
            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
            elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
                    CALL cross_product(qq(:,ii),bb,exb)
            kmult = kcoeff*exb(k)*NNif*b(k)
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
                    elMat%ALL(ind_ff(ind_if),ind_ff(ind_if),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_if),iel) - kmult
            kmultf = kcoeff*exb(k)*uf(i)*Nif*b(k)
            elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
            elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf

            !tau(i,i) = tau(i,i) + abs(kcoeff*exb(k)*b(k))
            tau(i,i) = 100
          ENDIF
          DO ii = 1,Neq
            IF (ii == i) CYCLE ! diagonal alredy assembled
                    IF (ABS(diffiso(i,ii)) < 1e-12 .AND. ABS(diffani(i,ii)) < 1e-12) CYCLE
            ind_kf = ind_ash + k + (ii - 1)*Ndim
            kcoeff = 1.
            ! Non-linear correction for non-linear diffusive terms.
            ! TODO: find a smarter way to include it,avoiding if statements and model dependencies (i==3,ii==1 only holds for Isothermal+Vorticity model)

            !IF ((i == 3 .or. i == 4) .and. ii == 1) then
                    IF ((i == 3 ) .AND. ii == 1) THEN
              ! Non-linear term in the vorticity equation (\Grad// n/n b)
              ind_jf = ind_asf + ii
              kcoeff = 1./uf(1)
              kmult = kcoeff**2*(diffiso(i,ii)*Qpr(k,1)*n(k)*NNif - diffani(i,ii)*Qpr(k,1)*b(k)*NNif*bn)
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult
                       elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) + kmult
              kmultf = kcoeff*(diffiso(i,ii)*Qpr(k,1)*n(k)*Nif - diffani(i,ii)*Qpr(k,1)*b(k)*Nfbn)
              elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) + kmultf
              elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) + kmultf
            ENDIF
            kmult = NNif*kcoeff*(n(k)*diffiso(i,ii) - bn*b(k)*diffani(i,ii))
            elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
          END DO
#endif
        END DO ! k-loop

        ! Convection contribution
        DO j = 1,Neq
          ind_jf = ind_asf + j
          kmult = bn*A(i,j)*NNif
!#ifdef NEUTRALP
!          IF (i == 5) kmult = kmult + bn*Anp(j)*NNif
!#endif
          elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult
                 elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) + kmult
!#ifdef TEMPERATURE
!#ifdef NEUTRAL
!          !X component neutral convective velocity
!          kmult = n(1)*Ax(i,j)*NNif
!          elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
!          elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
!          !Y component neutral convective velocity
!          kmult = n(2)*Ay(i,j)*NNif
!          elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
!          elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
!#endif
!#endif
        END DO ! j-loop
        ! Pinch contribution
        elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) + (APinch(i,1)*n(1) + APinch(i,2)*n(2))*NNif
              elMat%ALL(ind_ff(ind_if),ind_ff(ind_if),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_if),iel) + (APinch(i,1)*n(1) + APinch(i,2)*n(2))*NNif

#ifndef TEMPERATURE
        ! Added term for n=exp(x) change of variable
              IF (switch%logrho) THEN
                 CALL logrhojacobianVector(uf,upf,auxvec)
          kmultf = Nfbn*auxvec(i)
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel)-kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) -kmultf
              ENDIF
#endif

#ifdef TEMPERATURE
        ! Parallel diffusion for the temperature
        IF (i == 3) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
                    kmult = coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_PRODUCT(Taui(:,j),b)))*NNif*bn
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                    elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = coefi*Alphai*Vveci(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
              elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
                 kmultf = coefi*Alphai*(dot_PRODUCT(MATMUL(TRANSPOSE(Taui),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
        ELSEIF (i == 4) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
                    kmult = coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_PRODUCT(Taue(:,j),b)))*NNif*bn
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                    elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = coefe*Alphae*Vvece(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
              elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
                 kmultf = coefe*Alphae*(dot_PRODUCT(MATMUL(TRANSPOSE(Taue),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
#ifdef DNNLINEARIZED
      ELSEIF (i == 5) THEN
        DO j=1,Neq
          ind_jf = ind_asf+j
                    DO k=1,Ndim

            kmult = Dnn_dU(j)*Qpr(k,i)*n(k)*NNif
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                       elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel)  = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
                    ENDDO
        END DO
        kmultf = Dnn_dU_U*(Qpr(1,i)*n(1)+Qpr(2,i)*n(2))*Nif
        elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
        elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
#endif
#ifdef NEUTRALP
       ELSEIF (i == 5) THEN
          DO j = 1,Neq
             ind_jf = ind_asf + j
             DO k = 1,Ndim
                ind_kf = k + (j - 1)*Ndim + ind_ash
                kmult = (Dpn*Taupn(k,j) + dDpn_dU(j)*gmpn(k))*n(k)*NNif 
                !kmult = kmult - Gammaredpn*(Dpn*Taui(k,j) + dDpn_dU(j)*gmipn(k))*n(k)*NNif 
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                       elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
                kmult = Dpn*Vpn(j)*n(k)*NNif
                !kmult = kmult - Gammaredpn*Dpn*Vveci(j)*n(k)*NNif
                IF (j == 5) THEN
                   kmult = kmult + Dnn*n(k)*NNif
                END IF 
                elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
                elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult
             END DO
          END DO
                 kmultf = dot_PRODUCT(dDpn_dU,uf)*(gmpn(1)*n(1) + gmpn(2)*n(2))*Nif
          !kmultf = kmultf - Gammaredpn*(Dpn*(dot_product(Taui(1,:),uf)*n(1) + dot_product(Taui(2,:),uf)*n(2)) + dot_product(dDpn_dU,uf)*(gmipn(1)*n(1) + gmipn(2)*n(2)))*Nif             
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf 
#endif                                                                           
       END IF
#ifdef KEQUATION
#ifdef DKLINEARIZED
       if (i .ne. 5) then
        DO j=1,Neq
          ind_jf = ind_asf+j
          do k=1,Ndim
            kmult = ddk_dU(j)*Qpr(k,i)*(n(k)-b(k)*bn)*NNif
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel)  = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
          enddo
        enddo
        kmultf = ddk_dU_U*((Qpr(1,i)*n(1)+Qpr(2,i)*n(2))*Nif-(Qpr(1,i)*b(1)+Qpr(2,i)*b(2))*Nfbn)
        elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
        elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
      endif
#endif
#endif

#endif
      END DO  ! i-Loop

      ! Assembly stabilization terms
      IF (numer%stab < 6) THEN
        DO i = 1,Neq
          ind_if = i + ind_asf
          kmult = tau(i,i)*NNif
          elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) + kmult
          elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
                 elMat%ALL(ind_ff(ind_if),ind_ff(ind_if),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_if),iel) - kmult
          elMat%Alu(ind_ff(ind_if),ind_fe(ind_if),iel) = elMat%Alu(ind_ff(ind_if),ind_fe(ind_if),iel) + kmult
        END DO
      ELSE
        DO i = 1,Neq
          ind_if = i + ind_asf
          DO j = 1,Neq
            ind_jf = j + ind_asf
            kmult = tau(i,j)*NNif
            elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) + kmult
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                    elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
            elMat%Alu(ind_ff(ind_if),ind_fe(ind_jf),iel) = elMat%Alu(ind_ff(ind_if),ind_fe(ind_jf),iel) + kmult
          END DO
        END DO
      ENDIF
      !************* End stabilization terms************************

    END SUBROUTINE assemblyIntFacesContribution

    !********************************************************************
    !
    !         ASSEMBLY EXTERIOR FACES CONTRIBUTION
    !
    !********************************************************************

#ifdef DKLINEARIZED
    SUBROUTINE assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,&
      &ind_fg,b3,Bmod,psi,q_cyl,xyf,n,diffiso,diffani,NNif,Nif,Nfbn,uf,upf,qf,tau,Vnng,ifa)
#else
    SUBROUTINE assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,&
        &ind_fg,b3,Bmod,psi,n,diffiso,diffani,NNif,Nif,Nfbn,uf,upf,qf,tau,Vnng,ifa)
#endif
      integer*4,intent(IN)      :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
      logical                   :: isdir
      real*8,intent(IN)         :: b3(:),n(:),Bmod, psi
      real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
      real*8,intent(IN)         :: NNif(:,:),Nif(:),Nfbn(:)
      real*8,intent(IN)         :: uf(:),upf(:)
      real*8,intent(IN)         :: qf(:)
#ifdef KEQUATION
#ifdef DKLINEARIZED
      real*8,intent(IN)         :: q_cyl, xyf(:)
#endif
#endif
      real*8,optional,intent(INOUT) :: tau(:,:),Vnng(:)
      integer*4,optional         :: ifa
      real*8                    :: kcoeff
      integer*4                 :: i,j,k,ii,alpha,beta
      integer*4,dimension(Npfl)  :: ind_if,ind_jf,ind_kf
      real*8,dimension(neq,neq) :: A
      real*8,dimension(neq,Ndim):: APinch
      real*8                    :: auxvec(neq)
      real*8                    :: bn,kmult(Npfl,Npfl),kmultf(Npfl)
      real*8                    :: Qpr(Ndim,Neq),exb(3)
      real*8                    :: nn(3),qq(3,Neq),b(Ndim),bb(3)
      real*8                    :: W2(Neq), dW2_dU(Neq,Neq), QdW2(Ndim,Neq)
#ifdef TEMPERATURE
      real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
      real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
      real*8                    :: W3(Neq), dW3_dU(Neq,Neq), QdW3(Ndim,Neq)
      real*8                    :: W4(Neq), dW4_dU(Neq,Neq), QdW4(Ndim,Neq)
#ifdef DNNLINEARIZED
      real*8                    :: Dnn_dU(Neq), Dnn_dU_U
      real*8                    :: gradDnn(Ndim)
#endif
#ifdef KEQUATION
#ifdef DKLINEARIZED
      real*8                    :: ddk_dU(Neq), ddk_dU_U
      real*8                    :: gradddk(Ndim)
#endif
#endif
#ifdef NEUTRALP
      real*8                    :: Dnn,Dpn,GammaLim,Alphanp,Betanp,Gammaredpn,Tmin
      real*8                    :: Anp(Neq),Vpn(Neq),dVpn_dU(Neq,Neq),gmpn(Ndim),gmipn(Ndim),Taupn(Ndim,Neq),dDpn_dU(Neq)
#endif
#endif

      b = b3(1:Ndim)
      bb = b3
      ! Jacobian matrices
           bn = dot_PRODUCT(b,n)
      CALL jacobianMatrices(uf,A)

      ! Jacobian matrices Pinch
      CALL computePinch(b,psi,APinch)

      ! Compute Q^T^(k-1)
           Qpr = RESHAPE(qf,(/Ndim,Neq/))

      ! Split diffusion vector/matrix for the momentum equation
      CALL compute_W2(uf,W2,diffiso(1,1),diffiso(2,2))
      CALL compute_dW2_dU(uf,dW2_dU,diffiso(1,1),diffiso(2,2))
           QdW2 = MATMUL(Qpr,dW2_dU)

      nn = 0.
      qq = 0.
      nn(1:Ndim) = n
      qq(1:Ndim,:) = Qpr

#ifdef TEMPERATURE

      ! Compute V(U^(k-1))
           CALL computeVi(uf,Vveci)
           CALL computeVe(uf,Vvece)

      ! Compute dV_dU (k-1)
           CALL compute_dV_dUi(uf,dV_dUi)
           CALL compute_dV_dUe(uf,dV_dUe)

      ! Split diffusion vector/matrix for the energies equations
      CALL compute_W3(uf,W3,diffiso(1,1),diffiso(2,2),diffiso(3,3))
      CALL compute_dW3_dU(uf,dW3_dU,diffiso(1,1),diffiso(2,2),diffiso(3,3))
           QdW3 = MATMUL(Qpr,dW3_dU)

      CALL compute_W4(uf,W4,diffiso(1,1),diffiso(4,4))
      CALL compute_dW4_dU(uf,dW4_dU,diffiso(1,1),diffiso(4,4))
           QdW4 = MATMUL(Qpr,dW4_dU)

      ! Compute Alpha(U^(k-1))
      Alphai = computeAlphai(uf)
      Alphae = computeAlphae(uf)

      ! Compute dAlpha/dU^(k-1)
           CALL compute_dAlpha_dUi(uf,dAlpha_dUi)
           CALL compute_dAlpha_dUe(uf,dAlpha_dUe)

           gmi = dot_PRODUCT(MATMUL(Qpr,Vveci),b)  ! scalar
           gme = dot_PRODUCT(MATMUL(Qpr,Vvece),b)
           Taui = MATMUL(Qpr,dV_dUi)               ! 2x3
           Taue = MATMUL(Qpr,dV_dUe)               ! 2x3

#ifdef DNNLINEARIZED
           CALL compute_Dnn_dU(uf,Dnn_dU)

      Dnn_dU_u = dot_product(Dnn_dU,uf)
#endif

#ifdef KEQUATION
#ifdef DKLINEARIZED
      call compute_ddk_dU(uf,xyf,q_cyl,ddk_dU)

      ddk_dU_u = dot_product(ddk_dU,uf)
#endif
#endif
#ifdef NEUTRALP
      ! Compute Vpn(U^(k-1))
      CALL computeVpn(uf,Vpn)
           gmpn = MATMUL(Qpr,Vpn)                       ! Ndim x 1
      ! Compute dVpn-dU(U^(k-1))
	    CALL compute_dVpn_dU(uf,dVpn_dU)
           Taupn = MATMUL(Qpr,dVpn_dU)                 ! Ndim x Neq
	    ! Compute Dpn(U^(k-1))
      CALL computeDpn(uf,Qpr,Vpn,Dpn)
      ! Compute dDpn_dU(U^(k-1))
      CALL compute_dDpn_dU(uf,Qpr,Vpn,dDpn_dU)
      ! Reduce Grad Pn for low collision regime 
      ! Threshold set at 0.5xGradPn for Ti = 0.2 eV 
      Gammaredpn = 1.
      Tmin = 0.2/simpar%refval_temperature
           IF (Tmin/upf(7) .LE. 1.) Gammaredpn = Gammaredpn*Tmin/upf(7)
      Dnn = Gammaredpn*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*Dpn 
      ! Comput Gammaredpn(U^(k-1))
      !CALL computeGammared(uf,Gammaredpn)
      !gmipn = matmul(Qpr,Vveci)
      !CALL computeGammaLim(ue,Qpr,Vpn,GammaLim)
      ! Set Grad Ti = 0. for low collision regime 
      ! (back to diffusion equation for neutral density)
      !CALL computeAlphaCoeff(uf,Qpr,Vpn,Alphanp)  
      !CALL computeBetaCoeff(uf,Qpr,Vpn,Betanp)  
      !Dpn = Alphanp*Dpn
      !dDpn_dU = Alphanp*dDpn_dU
      !Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*Dpn      
      !IF (Dnn .gt. phys%diff_nn) Dnn = phys%diff_nn    
      !IF (Dpn .gt. phys%diff_nn) THEN
      ! Dpn = Alphanp*Dpn !0.
      ! dDpn_dU = Alphanp*dDpn_dU !0.
      ! Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*phys%diff_nn
      ! END IF  
      ! Set Gamma Convective = cs_n*n_n for low collision regime 
      !IF (Dpn .gt. phys%diff_nn) THEN
      !   Dpn = 0.
      !   dDpn_dU = 0.
      !   CALL jacobianMatricesNP(uf,Anp)
      !ELSE
      !   Anp = 0.
      !END IF   
#endif
#endif

      ! Assembly local matrix
      DO i = 1,Neq
        ind_if = ind_asf + i
        DO k = 1,Ndim
          ind_kf = ind_ash + k + (i - 1)*Ndim
          kmult = NNif*(n(k)*diffiso(i,i) - bn*b(k)*diffani(i,i))
#ifdef VORTICITY
                 IF (i==4) THEN
            kmult=kmult/uf(1)
                 ENDIF
#endif
          ! Diagonal terms for perpendicular diffusion
          elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
          IF(i == 2) THEN
            ! Assembly LQ
            j = 1 ! other terms are 0 anyway in vector W2
            ind_kf = ind_ash + k + (j - 1)*Ndim
             kmult = NNif*W2(j)*(n(k) - bn*b(k))
             elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
          ! Assembly LU
             DO j=1,Neq
               ind_jf = ind_asf+j
                       IF (.NOT. isdir) THEN
                kmult = QdW2(k,j)*NNif*(n(k)-bn*b(k))
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
               ENDIF
             END DO
#ifdef TEMPERATURE
         ELSEIF (i ==3) THEN
              ! Assembly LQ
              DO j=1,Neq ! here there are 2 non-zero elements in vector W3
               ind_kf = ind_ash + k + (j - 1)*Ndim
               kmult = NNif*W3(j)*(n(k) - bn*b(k))
               elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
             END DO
            ! Assembly LU
               DO j=1,Neq
                 ind_jf = ind_asf+j
                       IF (.NOT. isdir) THEN
                   kmult = QdW3(k,j)*NNif*(n(k)-bn*b(k))
                  elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                 ENDIF
              END DO
          ELSEIF (i ==4) THEN
              ! Assembly LQ
              j = 1 !other terms are 0 anyway in vector W4
              ind_kf = ind_ash + k + (j - 1)*Ndim
              kmult = NNif*W4(j)*(n(k) - bn*b(k))
               elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
            ! Assembly LU
               DO j=1,Neq
                 ind_jf = ind_asf+j
                       IF (.NOT. isdir) THEN
                   kmult = QdW4(k,j)*NNif*(n(k)-bn*b(k))
                  elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                  ENDIF
                 !ind_jf = ind_asf + j
                 !IF (.not. isdir) THEN
                !    kmult = coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_product(Taue(:,j),b)))*NNif*bn
                !    elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                ! END IF
               END DO

#endif
END IF



#ifdef VORTICITY
                 IF (switch%driftexb .AND. i .NE. 4) THEN
            ! ExB terms
            ii = 4
            kcoeff = phys%dfcoef*numer%exbdump/Bmod
                    CALL ijk_cross_product(k,alpha,beta)
            ind_kf = ind_ash + k + (ii - 1)*Ndim
            kmult = kcoeff*NNif*(nn(alpha)*b3(beta) - nn(beta)*b3(alpha))*uf(i)
            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
                    IF (.NOT. isdir) THEN
                       CALL cross_product(qq(:,ii),bb,exb)
              kmult = kcoeff*exb(k)*NNif*b(k)
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
              kmultf = kcoeff*exb(k)*uf(i)*Nif*b(k)
              elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
            ENDIF
                    tau(i,i) = tau(i,i) + ABS(kcoeff*exb(k)*b(k))
            tau(i,i) = 100
          ENDIF
          DO ii = 1,Neq
            IF (ii == i) CYCLE ! diagonal alredy assembled
                    IF (ABS(diffiso(i,ii)) < 1e-12 .AND. ABS(diffani(i,ii)) < 1e-12) CYCLE
            ind_kf = ind_ash + k + (ii - 1)*Ndim
            kcoeff = 1.
            ! Non-linear correction for non-linear diffusive terms.
            ! TODO: find a smarter way to include it,avoiding if statements and model dependencies (i==3,ii==1 only holds for Isothermal+Vorticity model)



            !               IF ((i == 3 .or. i == 4) .and. ii == 1) then
                    IF ((i == 3 ) .AND. ii == 1) THEN


              kcoeff = 1./uf(1)
                       IF (.NOT. isdir) THEN
                ind_jf = ind_asf + ii
                kmult = kcoeff**2*(diffiso(i,ii)*Qpr(k,1)*n(k)*NNif - diffani(i,ii)*Qpr(k,1)*b(k)*NNif*bn)
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult
                elMat%S(ind_fe(ind_if),iel)  = elMat%S(ind_fe(ind_if),iel) + kcoeff*(diffiso(i,ii)*Qpr(k,1)*n(k)*Nif - &
                  &diffani(i,ii)*Qpr(k,1)*b(k)*Nfbn)
              END IF
            ENDIF
            kmult = NNif*kcoeff*(n(k)*diffiso(i,ii) - bn*b(k)*diffani(i,ii))
            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
          END DO
#endif
        END DO ! k-loop

        ! Convection contribution
              IF (.NOT. isdir) THEN
          DO j = 1,Neq
            ind_jf = ind_asf + j
            kmult = bn*A(i,j)*NNif
!#ifdef NEUTRALP
!            IF (i == 5) kmult = kmult + bn*Anp(j)*NNif
!#endif
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult
!#ifdef TEMPERATURE
!#ifdef NEUTRAL
!            !X component neutral convective velocity
!            kmult = n(1)*Ax(i,j)*NNif
!            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
!            !Y component neutral convective velocity
!            kmult = n(2)*Ay(i,j)*NNif
!            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
!#endif
!#endif
          END DO ! j-loop
          ! Pinch Contribution
          elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) + (APinch(i,1)*n(1) + APinch(i,2)*n(2))*NNif
        ENDIF


#ifndef TEMPERATURE
        ! Added term for n=exp(x) change of variable
              IF (.NOT. isdir) THEN
                 IF (switch%logrho) THEN
                    CALL logrhojacobianVector(uf,upf,auxvec)
            kmultf = Nfbn*auxvec(i)
            elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel)-kmultf
                 ENDIF
              ENDIF
#endif



#ifdef TEMPERATURE
        ! Parallel diffusion for the temperature
        IF (i == 3) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
                    IF (.NOT. isdir) THEN
                       kmult = coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_PRODUCT(Taui(:,j),b)))*NNif*bn
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            END IF
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = coefi*Alphai*Vveci(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
                 kmultf = coefi*Alphai*(dot_PRODUCT(MATMUL(TRANSPOSE(Taui),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
        ELSEIF (i == 4) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
                    IF (.NOT. isdir) THEN
                       kmult = coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_PRODUCT(Taue(:,j),b)))*NNif*bn
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            END IF
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = coefe*Alphae*Vvece(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
                 kmultf = coefe*Alphae*(dot_PRODUCT(MATMUL(TRANSPOSE(Taue),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
#ifdef DNNLINEARIZED
        ELSEIF (i == 5) THEN
            DO j=1,Neq
              ind_jf = ind_asf+j
              DO k = 1,Ndim

                kmult = Dnn_dU(j)*Qpr(k,i)*n(k)*NNif
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                    ENDDO
                 ENDDO
            kmultf = Dnn_dU_U*(Qpr(1,i)*n(1)+Qpr(2,i)*n(2))*Nif
            elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
#endif           

#ifdef NEUTRALP
       ELSEIF (i == 5) THEN
          DO j = 1,Neq
             ind_jf = ind_asf + j
             DO k = 1,Ndim
                ind_kf = k + (j - 1)*Ndim + ind_ash
                kmult = (Dpn*Taupn(k,j) + dDpn_dU(j)*gmpn(k))*n(k)*NNif 
                !kmult = kmult - Gammaredpn*(Dpn*Taui(k,j) + dDpn_dU(j)*gmipn(k))*n(k)*NNif 
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                kmult = Dpn*Vpn(j)*n(k)*NNif
                !kmult = kmult - Gammaredpn*Dpn*Vveci(j)*n(k)*NNif
                IF (j == 5) THEN
                   kmult = kmult + Dnn*n(k)*NNif
                END IF
                elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
             END DO
          END DO
                 kmultf = dot_PRODUCT(dDpn_dU,uf)*(gmpn(1)*n(1) + gmpn(2)*n(2))*Nif
          !kmultf = kmultf - Gammaredpn*(Dpn*(dot_product(Taui(1,:),uf)*n(1) + dot_product(Taui(2,:),uf)*n(2)) + dot_product(dDpn_dU,uf)*(gmipn(1)*n(1) + gmipn(2)*n(2)))*Nif                             
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
#endif
      END IF

#ifdef KEQUATION
#ifdef DKLINEARIZED
        if (i .ne. 5) then
          DO j = 1,Neq
            ind_jf = ind_asf + j
            DO k = 1,Ndim
              kmult = ddk_dU(j)*Qpr(k,i)*(n(k)-bn*b(k))*NNif
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            enddo
          enddo
          kmultf = ddk_dU_U*((Qpr(1,i)*n(1)+Qpr(2,i)*n(2))*Nif-(Qpr(1,i)*b(1)+Qpr(2,i)*b(2))*Nfbn)
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
        endif
#endif
#endif

!if below for TEMPERATURE FLAG
#endif 
      END DO  ! i-Loop

      ! Assembly stabilization terms
      IF (numer%stab < 6) THEN
        DO i = 1,Neq
          ind_if = i + ind_asf
          kmult = tau(i,i)*NNif
          elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) + kmult
                 IF (.NOT. isdir) THEN
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
          ENDIF
        END DO
      ELSE
        DO i = 1,Neq
          ind_if = i + ind_asf
          DO j = 1,Neq
            ind_jf = j + ind_asf
            kmult = tau(i,j)*NNif
            elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) + kmult
                    IF (.NOT. isdir) THEN
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            ENDIF
          END DO
        END DO
      ENDIF
      !************* End stabilization terms************************
    END SUBROUTINE assemblyExtFacesContribution

    SUBROUTINE do_assembly(Auq,Auu,rhs,ind_ass,ind_asq,iel)
           REAL*8,INTENT(in)    :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)
           INTEGER*4,INTENT(in) :: iel,ind_ass(:),ind_asq(:)
           INTEGER :: i,j,k,z
           INTEGER*4,DIMENSION(Npel) :: ind_i,ind_j,ind_k

      DO i = 1,Neq
        ind_i = i + ind_ass
        elMat%S(ind_i,iel)=elMat%S(ind_i,iel)+rhs(:,i)
        DO j = 1,Neq
          ind_j = j + ind_ass
          z = i+(j-1)*Neq
          elMat%Auu(ind_i,ind_j,iel)=elMat%Auu(ind_i,ind_j,iel)+Auu(:,:,z)
          DO k = 1,Ndim
            z = i+(k-1)*Neq+(j-1)*Ndim*Neq
            ind_k = ind_asq + k + (j - 1)*Ndim
            elMat%Auq(ind_i,ind_k,iel)=elMat%Auq(ind_i,ind_k,iel)+Auq(:,:,z)
          END DO
        END DO
      END DO

    END SUBROUTINE do_assembly

#ifdef NEUTRAL
  !********************************************************************
  !
  !         ASSEMBLY NEUTRAL SOURCE MATRIX
  !
  !********************************************************************
#ifdef TEMPERATURE
  SUBROUTINE assemblyNeutral(U,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
      &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
      &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,Tloss,dTloss_dU,Tlossrec,dTlossrec_dU,sigmavEiz,dsigmavEiz_dU,sigmavErec,dsigmavErec_dU,Sn,Sn0)
#else
    SUBROUTINE assemblyNeutral(U,niz,dniz_dU,nrec,dnrec_dU,fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,Sn,Sn0)
#endif
             REAL*8, INTENT(IN) :: niz,nrec,fGammacx,fGammarec
             REAL*8, INTENT(IN) :: U(:),dniz_dU(:),dnrec_dU(:),dfGammacx_dU(:),dfGammarec_dU(:)
#ifndef TEMPERATURE
             REAL*8             :: sigmaviz,sigmavrec,sigmavcx
#else
             REAL*8, INTENT(IN) :: sigmaviz,sigmavrec,sigmavcx,fEiiz,fEirec,fEicx,Tloss,Tlossrec
             REAL*8, INTENT(IN) :: dsigmaviz_dU(:),dsigmavrec_dU(:),dsigmavcx_dU(:),dTloss_dU(:),dTlossrec_dU(:)
             REAL*8, INTENT(IN) :: dfEiiz_dU(:),dfEirec_dU(:),dfEicx_dU(:)
      !Amjuel energy losses
      real*8, intent(IN) :: sigmavEiz,sigmavErec
      real*8, intent(IN) :: dsigmavEiz_dU(:),dsigmavErec_dU(:)
#endif
             REAL*8             :: ad,ad4,RE,Sn(:,:),Sn0(:), Ti,Te

      Sn   = 0.
      Sn0  = 0.
      RE   = 0.2
      !RE   = 1.
      ad   = 1e19*1.374e-07 !n0*t0 !1e19*1.374e-07 !old 1.3737e12
      ad4  =  1e19*1.374e-07**3/1.901e-3**2*1.60217662e-19/3.35e-27 ! n0*t0/u0^2/m_i*e = 1e19*1.374e-07**3/1.901e-3**2*1.60217662e-19/3.35e-27 !old (ad*1.6e-19)/((1.3839e4**2)*3.35e-27)

#ifndef TEMPERATURE
      sigmaviz   = 3.01e-14
      sigmavrec  = 1.3638e-20
      sigmavcx   = 4.0808e-15
#endif


      !Assembly Source Terms in plasma density equation

      Sn(1,:)   = ad*(-dniz_dU(:)*sigmaviz + dnrec_dU(:)*sigmavrec)
#ifdef TEMPERATURE

      Sn(1,:)   = Sn(1,:) + ad*(-niz*dsigmaviz_dU(:) + nrec*dsigmavrec_dU(:))
#endif
      !Assembly Source Terms in plasma momentum equation

      Sn(2,:) = ad*(dfGammacx_dU(:)*sigmavcx + dfGammarec_dU(:)*sigmavrec)
#ifdef TEMPERATURE


      Sn(2,:)   = Sn(2,:) + ad*( fGammacx*dsigmavcx_dU(:) + fGammarec*dsigmavrec_dU(:))
      
      !Assembly Source Terms in ion energy equation

      Sn(3,:) = ad*(-RE*dfEiiz_dU(:)*sigmaviz + dfEirec_dU(:)*sigmavrec+dfEicx_dU(:)*sigmavcx)
      Sn(3,:) = Sn(3,:) + ad*(-RE*fEiiz*dsigmaviz_dU(:) + fEirec*dsigmavrec_dU(:) + fEicx*dsigmavcx_dU(:))
      !Assembly Source Terms in electron energy equation


      !Sn(4,:) =  ad4*(dniz_dU(:)*sigmaviz*Tloss + niz*dsigmaviz_dU(:)*Tloss + niz*sigmaviz*dTloss_dU(:) +&
            !   &dnrec_dU(:)*sigmavrec*Tlossrec + nrec*dsigmavrec_dU(:)*Tlossrec + nrec*sigmavrec*dTlossrec_dU(:))
      !AMJUEL rates
      Sn(4,:) =  ad4*(dniz_dU(:)*sigmavEiz + niz*dsigmavEiz_dU(:) +&
        &dnrec_dU(:)*sigmavErec + nrec*dsigmavErec_dU(:))

      !modification with recombination gain
      Sn(4,:) =  Sn(4,:)+ ad4*(-1.*dnrec_dU(:)*sigmavrec*13.6 - nrec*dsigmavrec_dU(:)*13.6)


#endif
      !Assembly Source Terms in neutral density equation
      Sn(5,:) = -Sn(1,:)

      !Assembly RHS Neutral Source Terms
      Sn0(1)    = ad*(niz*sigmaviz - nrec*sigmavrec)
      Sn0(2)    = ad*(-fGammacx*sigmavcx - fGammarec*sigmavrec)
#ifdef AMJUELSPLINES
             Sn0(1)    = Sn0(1) + ad*(niz*dot_PRODUCT(dsigmaviz_dU,U) - nrec*dot_PRODUCT(dsigmavrec_dU,U))
             Sn0(2)    = Sn0(2) + ad*(- fGammarec*dot_PRODUCT(dsigmavrec_dU,U))
#endif
#ifdef TEMPERATURE
      Sn0(3)    = ad*(RE*fEiiz*sigmaviz - fEirec*sigmavrec - fEicx*sigmavcx)
      !Sn0(4)    = ad4*(-niz*sigmaviz*Tloss - nrec*sigmavrec*Tlossrec)
      !modification with recombination gain
      Sn0(4)    = ad4*(nrec*sigmavrec*13.6)
      Sn0(4)    = Sn0(4) + ad4*(-niz*sigmavEiz- nrec*sigmavErec)
#ifdef AMJUELSPLINES
             Sn0(3)    = Sn0(3) + ad*(RE*fEiiz*dot_PRODUCT(dsigmaviz_dU,U) - fEirec*dot_PRODUCT(dsigmavrec_dU,U))
             !Sn0(4)    = Sn0(4) + ad4*(-niz*dot_PRODUCT(dsigmaviz_dU,U)*Tloss - nrec*dot_PRODUCT(dsigmavrec_dU,U)*Tlossrec)
      Sn0(4)    = Sn0(4) + ad4*(-niz*dot_product(dsigmavEiz_dU,U) - nrec*dot_product(dsigmavErec_dU,U))
      Sn0(4)    = Sn0(4) +  ad4*( nrec*dot_PRODUCT(dsigmavrec_dU,U)*13.6)
#endif
#endif
      Sn0(5)  = -Sn0(1)

      !Thresholds: 
#ifdef TEMPERATURE
                   Ti = 2./(3.*phys%Mref)*(U(3)/U(1) - 1./2.*(U(2)/U(1))**2)
                   Te = 2./(3.*phys%Mref)*U(4)/U(1)
                   !if ((Ti .le. 2.e-5) .or. (Te .le. 2.e-5)) then
             !         Sn(1,:) = - abs(Sn(1,:))
             !         Sn0(1) = - abs(Sn0(1))
             !         Sn(2,:) = - abs(Sn(2,:))
             !         Sn0(2) = - abs(Sn0(2))
                      !Sn(3,1) = ad*( - RE*fEiiz*dsigmaviz_dU(1))
                      !Sn(3,2) = 0.
                      !Sn(3,3) = ad*(-RE*dfEiiz_dU(3)*sigmaviz)
                      !Sn(3,4) = ad*(-RE*fEiiz*dsigmaviz_dU(4))
                      !Sn(3,5) = ad*(-RE*dfEiiz_dU(Neq)*sigmaviz)
                      !Sn0(3) = ad*(RE*fEiiz*sigmaviz)
                     !Sn(3,:) = - abs(Sn(3,:))
                     !Sn0(3)  = abs(Sn0(3))
                     !Sn(4,:) = - abs(Sn(4,:))
                     !Sn0(4)  = abs(Sn0(4))
#ifdef AMJUELSPLINES
                      !Sn0(3)    = Sn0(3) + ad*(RE*fEiiz*dot_product(dsigmaviz_dU,U))
#endif
                      !Sn(4,1) = 0. !3./2.*6.e-10
                      !Sn(4,2) = 0.
                      !Sn(4,3) = 0.
                      !Sn(4,4) = 0. !-1.
                      !Sn(4,5) = 0.
                      !Sn0(4) = 0.
             !         Sn(5,:) = - abs(Sn(5,:))
             !         Sn0(5) = - abs(Sn0(5))
             !         Sn(:,:) = 0.
             !         Sn0(:) = 0.
                   !endif
#endif

    END SUBROUTINE assemblyNeutral
#endif


         END SUBROUTINE hdg_ComputeJacobian
