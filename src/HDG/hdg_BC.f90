!*****************************************
! project: MHDG
! file: hdg_BC.f90
! date: 25/01/2017
! Set boundary conditions in the global
! problem
!*****************************************

SUBROUTINE HDG_BC()
  USE globals
  USE LinearAlgebra
  USE physics
  USE analytical
  USE printUtils
  USE in_out
  USE MPI_OMP

  IMPLICIT NONE
#ifdef TOR3D
  !***********************************************************************
  !
  !                            VERSION 3D TOROIDAL
  !
  !***********************************************************************
  INTEGER                   :: itor,iel3,ifa,ifl,iel,fl,bc,Neq,Npel,i,j,Ndim,Fi,N2D,Ngfl
  INTEGER*4                 :: Np2D,Np1Dpol,Np1Dtor,Npfl,ntorloc,itorg,Nfdir
  INTEGER                   :: nod(refElPol%Nfacenodes)
  INTEGER                   :: NGaussPol,NGaussTor,dd,delta
  INTEGER*4                 :: ind_uf(refElTor%Nfl*phys%neq),ind_fe(refElTor%Nfl*phys%neq),ind_ff(refElTor%Nfl*phys%neq)
  INTEGER*4                 :: ind_fg(refElTor%Nfl*3*phys%neq)
  INTEGER*4                 :: ind_loc(refElPol%Nfaces,refElTor%Nfl*phys%Neq),ind(refElPol%Ngauss1d)
  REAL*8                    :: tel(refElTor%Nnodes1d)
  REAL*8                    :: htor,Xf(Mesh%Nnodesperface,Mesh%Ndim)
  REAL*8                    :: uf(refElTor%Nfl*phys%neq)
  REAL*8                    :: thetafg(refElTor%Ngauss1d),dsurf(refElTor%Ngl),n_g(refElTor%Ngl,3),t_g(refElPol%Ngauss1d,2)
  INTEGER*4                 :: ind_qf(refElTor%Nfl*3*phys%Neq)
  REAL*8                    :: qqf(refElTor%Nfl,3*phys%Neq),qfg(refElTor%Ngl,phys%neq*3)
  REAL*8                    :: coefi,coefe
  REAL*8                    :: tdiv(numer%ntor + 1)
  REAL*8                    :: tau_save_el(refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Ngauss1d,2)
  REAL*8                    :: xyf(refElPol%Ngauss1d,2),xyd_g(refElPol%Ngauss1d,2),xydNorm_g(refElPol%Ngauss1d)
  INTEGER                   :: ind_asf(refElTor%Nfl),ind_ash(refElTor%Nfl)
  REAL*8                    :: Bfl(refElTor%Nfl,3)
  INTEGER                   :: indbt(refElTor%Nfl)
  REAL*8                    :: Bmod_nod(refElTor%Nfl),b_nod(refElTor%Nfl,3)
  REAL*8                    :: uex(refElTor%Ngl,phys%neq),uexpg(refElTor%Ngl,phys%npv)
  REAL*8                    :: upg(refElTor%Ngl,phys%npv)
  REAL*8,POINTER            :: Nfi(:,:)
  REAL*8                    :: ug(refElTor%Ngl),ufg(refElTor%Ngl,phys%neq)
  REAL*8                    :: diff_iso_fac(phys%neq,phys%neq,refElTor%Ngl),diff_ani_fac(phys%neq,phys%neq,refElTor%Ngl)
  REAL*8                    :: Bmod(refElTor%Ngl),b(refElTor%Ngl,3)

  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tps1)
     CALL system_CLOCK(timing%cks1,timing%clock_rate1)
  END IF

  Ndim = 3                                               ! Number of dimensions
  Neq = phys%neq
  N2D = Mesh%Nelems                                     ! Number elements in the 2D mesh
  Np2D = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
  Np1Dpol = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Npel = Np2D*refElTor%Nnodes1D                          ! Number of nodes for each 3D element
  Npfl = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
  Nfdir = Mesh%Ndir
  NGaussPol = refElPol%Ngauss1D
  NGaussTor = refElTor%Ngauss1D
  Ngfl = refElTor%Ngl

  ! Some indices
  ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
  ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

  ! Toroidal discretization
  htor = numer%tmax/numer%ntor
  tdiv = 0.
  DO i = 1,numer%ntor
    tdiv(i + 1) = i*htor
  END DO

  coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1 + phys%epn)
  coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1 + phys%epn)

#ifdef PARALL
  IF (MPIvar%ntor .GT. 1) THEN
    ntorloc = numer%ntor/MPIvar%ntor + 1
  ELSE
    ntorloc = numer%ntor
  ENDIF
#else
  ntorloc = numer%ntor
#endif

  ! Initialize matrices
  elMat%Aul_dir = 0.
  elMat%Aql_dir = 0.

  !************************
  ! Loop in external faces
  !************************
  DO itor = 1,ntorloc
#ifdef PARALL
    itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
     IF (itorg == numer%ntor + 1) itorg = 1
#else
    itorg = itor
#endif
    tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
    DO ifa = 1,Mesh%Nextfaces

      ! Global numbering of this face
      Fi = ifa + Mesh%Nintfaces

      ! Element to which this face belongs
      iel = Mesh%extfaces(ifa,1)
      iel3 = (itor - 1)*N2d+iel

      ! Face local numbering
      ifl = Mesh%extfaces(ifa,2)

      ! Nodes in local numbering
      nod = refElPol%Face_nodes(ifl,:)

      ! Coordinates of the nodes of the face
      Xf = Mesh%X(Mesh%T(iel,nod),:)

      !****************************************************
      !                      Magnetic field
      !****************************************************
      ! Magnetic field of the nodes of the face
      indbt = colint(tensorSumInt(Mesh%T(iel,nod),(itor - 1)*(Np1Dtor - 1)*Mesh%Nnodes +&
        &Mesh%Nnodes*((/(i,i=1,Np1Dtor)/) - 1)))
      Bfl = phys%B(indbt,:)
      ! Magnetic field norm and direction at element nodes
        Bmod_nod = SQRT(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
      b_nod(:,1) = Bfl(:,1)/Bmod_nod
      b_nod(:,2) = Bfl(:,2)/Bmod_nod
      b_nod(:,3) = Bfl(:,3)/Bmod_nod

      ! Gauss points position
        xyf = MATMUL(refElPol%N1D,Xf)
        thetafg = MATMUL(refElTor%N1D,tel)

      ! Set the face solution
      CALL analytical_solution(iel,xyf(:,1),xyf(:,2),thetafg,uex)

      ! Physical variables at Gauss points with analytical sol
      CALL cons2phys(uex,uexpg)

      ! Exterior normal & elemental surface
        xyd_g = MATMUL(refelPol%Nxi1D,Xf)
        xydNorm_g = SQRT(xyd_g(:,1)**2 + xyd_g(:,2)**2)
      dsurf = col(tensorProduct(refElPol%gauss_weights1D*xydNorm_g,refElTor%gauss_weights1D*0.5*htor))
      t_g(:,1) = xyd_g(:,1)/xydNorm_g
      t_g(:,2) = xyd_g(:,2)/xydNorm_g
      n_g = 0.
      DO i = 1,NGaussTor
        ind = (i - 1)*NGaussPol + (/(j,j=1,NGaussPol)/)
        n_g(ind,1) = t_g(:,2)
        n_g(ind,2) = -t_g(:,1)
      END DO

      ! Assembly indices
      dd = 1 + (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Npfl)*Neq
      delta = dd + (N2D*Np2D+(Fi - 1)*Npfl)*Neq
      ind_uf = delta + (/(i,i=0,Npfl*Neq - 1)/)
      ind_ff = Np2d*Neq + (ifl - 1)*Npfl*Neq + (/(i,i=1,Npfl*Neq)/)
      ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifl,:) - 1)*Neq))
      ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifl,:) - 1)*3*Neq))
      ind_qf = (iel3 - 1)*Ndim*Neq*Npel + ind_fg
        qqf = TRANSPOSE(RESHAPE(sol%q(ind_qf),(/Ndim*Neq,Npfl/)))

      ! Solution in this face
      uf = sol%u_tilde(ind_uf)

      ! Magnetic field norm and direction at Gauss points
        Bmod = MATMUL(refElTor%sFTF,Bmod_nod)
        b = MATMUL(refElTor%sFTF,b_nod)

      ! Face shape functions
      Nfi => refElTor%sFTF

      ! Gradient solution at face gauss points
        qfg = MATMUL(Nfi,qqf)

      ! Set the face solution
        ufg = MATMUL(Nfi,TRANSPOSE(RESHAPE(uf,[neq,Npfl])))

      ! Physical variables at Gauss points
      CALL cons2phys(ufg,upg)

      ! Compute diffusion at faces Gauss points
      CALL setLocalDiff(xyf,ufg,qfg,diff_iso_fac,diff_ani_fac)

      ! Type of boundary
      fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL
        IF (fl .EQ. 0) CYCLE
#endif
      bc = phys%bcflags(fl)

      SELECT CASE (bc)

      CASE (bc_dirichlet)
        CALL set_dirichletstr_bc()
      CASE (bc_dirichlet_weak_form)
        CALL set_dirichletwf_bc()
      CASE (bc_dirichlet_and_Neumann)
        CALL set_diric_neum_bc()
      CASE(bc_Transmission)
        CALL set_Transmission_bc()
      CASE (bc_Bohm)
        CALL set_Bohm_bc()
      CASE DEFAULT
        WRITE (6,*) "Error: wrong boundary type"
        STOP
      END SELECT

    END DO
  END DO


  IF (utils%timing) THEN
     CALL cpu_TIME(timing%tpe1)
     CALL system_CLOCK(timing%cke1,timing%clock_rate1)
     timing%runtbcd = timing%runtbcd + (timing%cke1 - timing%cks1)/REAL(timing%clock_rate1)
    timing%cputbcd = timing%cputbcd + timing%tpe1 - timing%tps1
  END IF


CONTAINS

  !****************************
  ! Dirichlet in weak form
  !****************************
  SUBROUTINE set_dirichletwf_bc
    INTEGER                   :: g,igpol,igtor
    REAL*8                    :: f(Neq,Npfl),dsurfg
    REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)

    ! Face shape functions
    Nfi => refElTor%sFTF

    !*****************************
    ! Loop in face Gauss points
    !*****************************
    DO igtor = 1,NGaussTor
      DO igpol = 1,NGaussPol

        g = (igtor - 1)*NGaussPol + igpol

        ! Calculate the integration weight
        IF (switch%axisym) THEN
          dsurfg = dsurf(g)*xyf(igpol,1)
        ELSE
          dsurfg = dsurf(g)
        END IF
        NiNi = dsurfg*tensorProduct(Nfi(g,:),Nfi(g,:))
        Ni = Nfi(g,:)*dsurfg

        ! Assembly contribution
        CALL assembly_dirichletwf_bc(iel3,ind_asf,ind_ff,uex(g,:),NiNi,Ni)

      END DO
    END DO

  END SUBROUTINE set_dirichletwf_bc

  !********************************
  ! Transmission boundary condition
  !********************************
  SUBROUTINE set_Transmission_bc
    INTEGER                   :: g,igpol,igtor
    REAL*8                    :: dsurfg
    REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)

    ! Face shape functions
    Nfi => refElTor%sFTF

    !*****************************
    ! Loop in face Gauss points
    !*****************************
    DO igtor = 1,NGaussTor
      DO igpol = 1,NGaussPol

        g = (igtor - 1)*NGaussPol + igpol

        ! Calculate the integration weight
        IF (switch%axisym) THEN
          dsurfg = dsurf(g)*xyf(igpol,1)
        ELSE
          dsurfg = dsurf(g)
        END IF
        NiNi = dsurfg*tensorProduct(Nfi(g,:),Nfi(g,:))
        Ni = Nfi(g,:)*dsurfg

        ! Assembly contribution
        CALL assembly_trasm_bc(iel3,ind_asf,ind_ff,NiNi,Ni,ufg(g,:))

      END DO
    END DO

  END SUBROUTINE set_Transmission_bc

  !***************************************************************
  ! Dirichlet in weak form in some variables and Neumann in others
  !***************************************************************
  SUBROUTINE set_diric_neum_bc
    INTEGER                   :: g,igpol,igtor
    REAL*8                    :: f(Neq,Npfl),dsurfg
    REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
    REAL                      :: tau_stab(Neq,Neq)


    ! Face shape functions
    Nfi => refElTor%sFTF

    !*****************************
    ! Loop in face Gauss points
    !*****************************
    DO igtor = 1,NGaussTor
      DO igpol = 1,NGaussPol

        g = (igtor - 1)*NGaussPol + igpol

        ! Calculate the integration weight
        IF (switch%axisym) THEN
          dsurfg = dsurf(g)*xyf(igpol,1)
        ELSE
          dsurfg = dsurf(g)
        END IF
        NiNi = dsurfg*tensorProduct(Nfi(g,:),Nfi(g,:))
        Ni = Nfi(g,:)*dsurfg


        !****************** Stabilization part******************************
        IF (numer%stab > 1) THEN
          ! Compute tau in the Gauss points
          IF (numer%stab < 6) THEN
            CALL computeTauGaussPoints(upg(g,:),ufg(g,:),b(g,:),n_g(g,:),iel,ifa,1.,xyf(g,:),tau_stab)
          ELSE
            CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),1.,iel,tau_stab)
          ENDIF
          !****************** End stabilization part******************************
        ELSE
          tau_stab=0.
          DO i = 1,Neq
            tau_stab(i,i) = numer%tau(i)
          END DO
        END IF

        ! Assembly contribution
        CALL assembly_diric_neum_bc(iel3,ind_asf,ind_ash,ind_ff,ind_fg,uex(g,:),NiNi,Ni,b(g,:),n_g(g,:),tau_stab,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g))
      END DO
    END DO

  END SUBROUTINE set_diric_neum_bc

  !****************************
  ! Dirichlet in strong form
  !****************************
  SUBROUTINE set_dirichletstr_bc
    INTEGER                   :: g,igpol,igtor
    REAL*8                    :: Ni(Npfl)
    REAL*8                    :: dsurfg

    ! Magnetic field norm and direction at Gauss points
    b = MATMUL(refElTor%sFTF,b_nod)

    ! Face shape functions
    Nfi => refElTor%sFTF

    ! Gradient solution at face gauss points
    qfg = MATMUL(Nfi,qqf)

    !*****************************
    ! Loop in face Gauss points
    !*****************************
    DO igtor = 1,NGaussTor
      DO igpol = 1,NGaussPol

        g = (igtor - 1)*NGaussPol + igpol

        ! Calculate the integration weight
        IF (switch%axisym) THEN
          dsurfg = dsurf(g)*xyf(igpol,1)
        ELSE
          dsurfg = dsurf(g)
        END IF

        ! Shaspe functions product
        Ni = Nfi(g,:)*dsurfg

        CALL assembly_dirichletstr_bc(iel3,ind_asf,ind_ash,ind_fe,ind_fg,Ni,qfg(g,:),uex(g,:),uexpg(g,:),b(g,:),n_g(g,:))

      END DO
    END DO
  END SUBROUTINE set_dirichletstr_bc

  !****************************
  ! Bohm
  !****************************
  SUBROUTINE set_Bohm_bc()
    INTEGER                   :: i,j,k
    INTEGER                   :: g,igpol,igtor
    REAL*8                    :: dsurfg
    INTEGER                    :: setval,delta
    REAL*8                    :: inc,sn
    REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
    REAL*8                    :: soundSpeed
    REAL*8                    :: bn
    LOGICAL                   :: ntang
    REAL                      :: tau_stab(Neq,Neq)

    !*****************************
    ! Loop in face Gauss points
    !*****************************
    DO igtor = 1,NGaussTor
      DO igpol = 1,NGaussPol

        g = (igtor - 1)*NGaussPol + igpol

        ! Calculate the integration weight
        IF (switch%axisym) THEN
          dsurfg = dsurf(g)*xyf(igpol,1)
        ELSE
          dsurfg = dsurf(g)
        END IF

          bn = dot_PRODUCT(b(g,:),n_g(g,:))

        ! Sound speed
#ifndef TEMPERATURE
          SoundSpeed = SQRT(phys%a)
#else
        SoundSpeed = upg(g,9)
          IF (switch%decoup) THEN
             SoundSpeed = SQRT(phys%Mref)
          ENDIF
#endif
        ! tangency
        ntang = .TRUE.
          inc = bn/NORM2(b(g,1:2))

        setval = ufg(g,2)

          IF (ABS(inc) .LE. phys%bohmth) THEN
          sn = 1.
          delta = 0
          ntang = .FALSE.
        ELSE
             sn = SIGN(1.,inc)
          delta = 1
             IF ((ABS(upg(g,2))) .LE. SoundSpeed) THEN
            setval = sn*SoundSpeed
            delta = 1
             ELSE IF ((ABS(upg(g,2))) .GT. SoundSpeed) THEN
            delta = 0
            setval = upg(g,2)
          END IF
             IF (numer%bohmtypebc==1) THEN
            delta=1
             ENDIF
        END IF

        ! Impose the velocity at the sound speed in case of subsonic
        ! velocity.
        NiNi = tensorProduct(Nfi(g,:),Nfi(g,:))*dsurfg
        Ni = Nfi(g,:)*dsurfg

        !****************** Stabilization part******************************
        IF (numer%stab > 1) THEN
          ! Compute tau in the Gauss points
          IF (numer%stab < 6) THEN
            CALL computeTauGaussPoints(upg(g,:),ufg(g,:),b(g,:),n_g(g,:),iel,ifa,1.,xyf(g,:),tau_stab)
          ELSE
            CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),1.,iel,tau_stab)
          ENDIF
          !****************** End stabilization part******************************
        ELSE
          tau_stab=0.
          DO i = 1,Neq
            tau_stab(i,i) = numer%tau(i)
          END DO
        END IF

        ! Assembly Bohm contribution
          IF (numer%bohmtypebc.EQ.0) THEN
          CALL assembly_bohm_bc(iel3,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
            &ufg(g,:),b(g,:),n_g(g,:),tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang)
          ELSE
          CALL assembly_bohm_bc_new(iel3,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
            &ufg(g,:),upg(g,:),b(g,:),n_g(g,:),tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang)
          ENDIF

      END DO
    END DO

  END SUBROUTINE set_Bohm_bc

#else
  !***********************************************************************
  !
  !                            VERSION 2D
  !
  !***********************************************************************
  INTEGER                   :: ifa,ifl,iel,fl,bc,Npfl,Neq,Npel,i,Ndim,Fi,Ng1d
  INTEGER                   :: nod(refElPol%Nfacenodes)
  INTEGER                   :: ind_uf(refElPol%Nfacenodes*phys%neq),ind_ff(refElPol%Nfacenodes*phys%neq)
  INTEGER                   :: ind_fe(refElPol%Nfacenodes*phys%neq)
  REAL*8                    :: Xf(Mesh%Nnodesperface,2)
  INTEGER                   :: ind_fG(refElPol%Nfacenodes*phys%neq*2)
  INTEGER                   :: ind_qf(refElPol%Nfacenodes*phys%neq*2)
  REAL*8                    :: uf(refElPol%Nfacenodes*phys%neq),qf(refElPol%Nfacenodes,phys%neq*2)
  REAL*8                    :: ue(refElPol%Nfacenodes,phys%neq)
  REAL*8                    :: coefi,coefe
  REAL*8                    :: tau_save_el(refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Ngauss1d,2)
  REAL*8                    :: v_nn_Bou_el(refElPol%Ngauss1D,Mesh%Ndim)
  LOGICAL                   :: save_tau
  REAL*8,ALLOCATABLE        :: tau_save(:,:)
  REAL*8,ALLOCATABLE        :: xy_g_save(:,:)
  REAL*8,ALLOCATABLE        :: ures(:,:)
  INTEGER*4                 :: inde(Mesh%Nnodesperelem)
  INTEGER                   :: sizeu
  INTEGER                   :: indtausave(refElPol%Ngauss1d)
  INTEGER                   :: ind_asf(refElPol%Nfacenodes),ind_ash(refElPol%Nfacenodes)
  REAL*8                    :: Bfl(refElPol%Nfacenodes,3),psifl(refElPol%Nfacenodes)
  REAL*8                    :: Bmod_nod(refElPol%Nfacenodes),b_nod(refElPol%Nfacenodes,3)
  REAL*8                    :: xyg(refElPol%Ngauss1d,2),xyder(refElPol%Ngauss1d,2)
  REAL*8                    :: Bmod(refElPol%Ngauss1d),b(refElPol%Ngauss1d,3), psig(refElPol%Ngauss1d)
  REAL*8                    :: ufg(refElPol%Ngauss1d,phys%neq),qfg(refElPol%Ngauss1d,phys%neq*2)
  REAL*8                    :: ueg(refElPol%Ngauss1d,phys%neq)
  REAL*8                    :: upg(refElPol%Ngauss1d,phys%npv),uexpg(refElPol%Ngauss1d,phys%npv)
  REAL*8                    :: uex(refElPol%Ngauss1d,phys%neq)
  REAL*8                    :: diff_iso_fac(phys%neq,phys%neq,refElPol%Ngauss1d)
  REAL*8                    :: diff_ani_fac(phys%neq,phys%neq,refElPol%Ngauss1d)
#ifdef KEQUATION
  real*8                    :: q_cylfl(refElPol%Nfacenodes),q_cyl(refElPol%Ngauss1d)
#endif
#ifdef PARALL
#ifdef SAVEFLUX
  INTEGER                   :: ierr
#endif
#endif
#ifdef SAVEFLUX
  real*8                    :: totalflux_pump, totalflux_puff, totalflux_parallel, totalflux_perpendicular,totalflux_neutral,totalflux_numerical
  real*8                    :: faceflux_pump, faceflux_puff, faceflux_parallel, faceflux_perpendicular,faceflux_neutral,faceflux_numerical

  totalflux_pump = 0.
  totalflux_puff = 0.
  totalflux_parallel = 0.
  totalflux_perpendicular = 0.
  totalflux_neutral = 0.
  totalflux_numerical = 0.
#endif

  save_tau = switch%saveTau
  Ndim = 2
  Npel = refElPol%Nnodes2D
  Npfl = refElPol%Nfacenodes
  Neq = phys%neq
  Ng1d = refElPol%Ngauss1d

  coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1 + phys%epn)
  coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1 + phys%epn)

  ! Some indices
  ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
  ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

  IF (save_tau) THEN
     ALLOCATE (tau_save(Mesh%Nextfaces*refElPol%Ngauss1d,phys%neq))
     ALLOCATE (xy_g_save(Mesh%Nextfaces*refElPol%Ngauss1d,2))
    tau_save = 0.
    xy_g_save = 0.
    v_nn_Bou_el = 0.
    tau_save_el = 0.
    xy_g_save_el = 0.
  ENDIF
  sizeu = SIZE(sol%u)

  ALLOCATE (ures(sizeu/Neq,Neq))
  ures = TRANSPOSE(RESHAPE(sol%u,[Neq,sizeu/Neq]))

  ! Initialize matrices
  elMat%Aul_dir = 0.
  elMat%Aql_dir = 0.
  !************************
  ! Loop in external faces
  !************************
  DO ifa = 1,Mesh%Nextfaces

    ! Global numbering of this face
    Fi = ifa + Mesh%Nintfaces

    ! Element to which this face belongs
    iel = Mesh%extfaces(ifa,1)

    ! Face local numbering
    ifl = Mesh%extfaces(ifa,2)

    ! Nodes in local numbering
    nod = refElPol%Face_nodes(Mesh%extfaces(ifa,2),:)

    ! Coordinates of the nodes of the face
    Xf = Mesh%X(Mesh%T(iel,nod),:)

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field of the nodes of the face
    Bfl = phys%B(Mesh%T(iel,nod),:)

    ! Magnetic field norm and direction at element nodes
     Bmod_nod = SQRT(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
    b_nod(:,1) = Bfl(:,1)/Bmod_nod
    b_nod(:,2) = Bfl(:,2)/Bmod_nod
    b_nod(:,3) = Bfl(:,3)/Bmod_nod

#ifdef KEQUATION
    if (switch%testcase == 60) then
      q_cylfl = geom%q
    else
      q_cylfl = phys%q_cyl(Mesh%T(iel,nod))
    endif
#endif

    ! Normalized magnetic flux of the nodes of the face: PSI
    psifl = phys%magnetic_psi(Mesh%T(iel,nod))

    ! Assembly indices
    ind_uf = (Fi - 1)*Neq*Npfl + (/(i,i=1,Neq*Npfl)/)
    ind_ff = (ifl - 1)*Neq*Npfl + (/(i,i=1,Neq*Npfl)/)
     ind_fe = RESHAPE(tensorSumInt((/(i,i=1,neq)/),neq*(nod - 1)),(/neq*Npfl/))
     ind_fG = RESHAPE(tensorSumInt((/(i,i=1,neq*Ndim)/),neq*Ndim*(refElPol%face_nodes(ifl,:) - 1)),(/neq*Ndim*Npfl/))
    ind_qf = (iel - 1)*Ndim*Neq*Npel + ind_fG


    ! Gauss points position
     xyg = MATMUL(refElPol%N1D,Xf)

    ! Shape function derivatives at Gauss points
     xyDer = MATMUL(refElPol%Nxi1D,Xf)

    ! Solution at nodes
    uf = sol%u_tilde(ind_uf)
     ! Elemental solution at nodes PROBABLY WRONG
     inde = (iel - 1)*Npel + (/(i,i=1,Npel)/)
     !ue = ures(Mesh%T(iel,nod),:)
     ue = ures(inde(nod),:)

    ! Solution gradient at nodes
     qf = TRANSPOSE(RESHAPE(sol%q(ind_qf),(/Ndim*Neq,Npfl/)))
     IF (switch%testcase .NE. 54) THEN !this we shouldn't actually call for WEST case
        ! Analytical solution at face Gauss points
        CALL analytical_solution(iel,xyg(:,1),xyg(:,2),uex)
     ENDIF

    ! Solution at face Gauss points
     ufg = MATMUL(refElPol%N1D,TRANSPOSE(RESHAPE(uf,[neq,Npfl])))

     ! Elemental solution at face Gauss points
     ueg = MATMUL(refElPol%N1D,ue)

    ! Gradient solution at face gauss points
     qfg = MATMUL(refElPol%N1D,qf)

    ! Magnetic field norm and direction at Gauss points
     Bmod = MATMUL(refElPol%N1d,Bmod_nod)
     b = MATMUL(refElPol%N1d,b_nod)

    ! Normalized magnetic flux at Gauss points: PSI
     psig = MATMUL(refElPol%N1d,psifl)

    ! Compute diffusion at faces Gauss points
#ifndef KEQUATION
    CALL setLocalDiff(xyg,ufg,qfg,diff_iso_fac,diff_ani_fac)
#else
    CALL setLocalDiff(xyg,ufg,qfg,diff_iso_fac,diff_ani_fac,q_cyl)
#endif
    if (save_tau) then
       indtausave = (ifa - 1)*refElPol%Ngauss1d+(/(i,i=1,refElPol%Ngauss1d)/)
       phys%diff_nn_Bou(indtausave) = diff_iso_fac(5,5,:)
     ENDIF

    ! Physical variables at Gauss points (division by 0 during convergence test!!!)
     IF ( switch%testcase .NE. 2 ) THEN
    CALL cons2phys(ufg,upg)
    end if

    if (switch%testcase .ne. 54) then !this we shouldn't actually call for WEST case
      ! Physical variables at Gauss points with analytical sol
      CALL cons2phys(uex,uexpg)
     ENDIF

#ifdef SAVEFLUX
    !Initialization of variables for flux control to avoid NaN if not Bohm boundary
    faceflux_pump = 0.
    faceflux_puff = 0.
    faceflux_parallel = 0.
    faceflux_perpendicular = 0.
    faceflux_neutral = 0.
    faceflux_numerical = 0.
#endif
    ! Type of boundary
    fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL
     IF (fl .EQ. 0) CYCLE
#endif
    bc = phys%bcflags(fl)

    SELECT CASE (bc)

    CASE (bc_dirichlet)
      CALL set_dirichletstr_bc()
    CASE (bc_dirichlet_weak_form)
      CALL set_dirichletwf_bc()
        IF (save_tau) THEN
         DO i = 1, Neq
            tau_save_el(:,i) = numer%tau(i)
            xy_g_save_el = xyg
         END DO
        ENDIF
    CASE (bc_dirichlet_and_Neumann)
      CALL set_diric_neum_bc()
    CASE(bc_NeumannH)
      CALL set_NeumannH_bc()
    CASE(bc_Transmission)
      CALL set_Transmission_bc()
    CASE(bc_periodic)
      CYCLE ! done in assembly
      !         CALL set_periodic_bc()
#ifndef SAVEFLUX
    CASE (bc_Bohm)
      CALL set_Bohm_bc(v_nn_Bou_el,tau_save_el,xy_g_save_el)
    CASE (bc_BohmPump)

      CALL set_Bohm_bc(v_nn_Bou_el,tau_save_el,xy_g_save_el)
    CASE (bc_BohmPuff)
      CALL set_Bohm_bc(v_nn_Bou_el,tau_save_el,xy_g_save_el)
#else
    CASE (bc_Bohm)
      CALL set_Bohm_bc(v_nn_Bou_el,tau_save_el,xy_g_save_el,faceflux_pump,faceflux_puff,faceflux_parallel,faceflux_perpendicular,faceflux_neutral,faceflux_numerical)
    CASE (bc_BohmPump)
      CALL set_Bohm_bc(v_nn_Bou_el,tau_save_el,xy_g_save_el,faceflux_pump,faceflux_puff,faceflux_parallel,faceflux_perpendicular,faceflux_neutral,faceflux_numerical)
    CASE (bc_BohmPuff)
      CALL set_Bohm_bc(v_nn_Bou_el,tau_save_el,xy_g_save_el,faceflux_pump,faceflux_puff,faceflux_parallel,faceflux_perpendicular,faceflux_neutral,faceflux_numerical)
#endif
    CASE (bc_iter_core)
      CALL set_itercore_bc()
    CASE DEFAULT
      WRITE (6,*) "Error: wrong boundary type"
      STOP
    END SELECT
#ifdef SAVEFLUX
#ifdef PARALL
     IF (Mesh%ghostFaces(Fi) .EQ. 0) THEN
#endif
    totalflux_pump = totalflux_pump + faceflux_pump
    totalflux_puff = totalflux_puff + faceflux_puff
    totalflux_parallel = totalflux_parallel + faceflux_parallel
    totalflux_perpendicular = totalflux_perpendicular + faceflux_perpendicular
    totalflux_neutral = totalflux_neutral + faceflux_neutral
    totalflux_numerical = totalflux_numerical + faceflux_numerical
#ifdef PARALL
  ENDIF
#endif
#endif
     IF (save_tau) THEN
       indtausave = (ifa - 1)*refElPol%Ngauss1d+(/(i,i=1,refElPol%Ngauss1d)/)
       phys%v_nn_Bou(indtausave,:) = v_nn_Bou_el
       tau_save(indtausave,:) = tau_save_el
       Mesh%Xgb(indtausave,:) = xy_g_save_el
       xy_g_save(indtausave,:) = xy_g_save_el
     ENDIF

  END DO

#ifdef SAVEFLUX
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, totalflux_pump, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, totalflux_puff, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, totalflux_parallel, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, totalflux_perpendicular, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, totalflux_neutral, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, totalflux_numerical, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  IF (MPIvar%glob_id.EQ.0) THEN
     WRITE(6,*) 'pumf = ',totalflux_pump
     WRITE(6,*) 'puff = ',totalflux_puff
     WRITE(6,*) 'plasma parallel = ',totalflux_parallel
     WRITE(6,*) 'plasma gradient = ',totalflux_perpendicular
     WRITE(6,*) 'neutral flux = ',totalflux_neutral
     WRITE(6,*) 'numerical flux = ',totalflux_numerical
     WRITE(6,*) 'net flux = ',totalflux_parallel-totalflux_perpendicular-totalflux_neutral+totalflux_puff-totalflux_pump
  ENDIF
#endif
  IF (save_tau) THEN
     WRITE (6,*) "Saving tau in the boundary faces"
     CALL saveMatrix(tau_save,'tau_save_bound')
     CALL saveMatrix(xy_g_save,'xy_g_save_bound')
     DEALLOCATE (tau_save,xy_g_save)
     WRITE (6,*) "Done saving tau!"
  ENDIF
  DEALLOCATE(ures)

CONTAINS

  !****************************
  ! Dirichlet in weak form
  !****************************
  SUBROUTINE set_dirichletwf_bc
    INTEGER                   :: g
    REAL*8                    :: dline,xyDerNorm_g
    REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)

    ! Loop in 1D Gauss points
    DO g = 1,Ng1d

      ! Calculate the integration weight
       xyDerNorm_g = NORM2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
      IF (switch%axisym) THEN
        dline = dline*xyg(g,1)
      END IF
      NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Ni = refElPol%N1D(g,:)*dline

      CALL assembly_dirichletwf_bc(iel,ind_asf,ind_ff,uex(g,:),NiNi,Ni)

    END DO

  END SUBROUTINE set_dirichletwf_bc


  !   !****************************
  !   ! Periodic bc
  !   !****************************
  !   SUBROUTINE set_periodic_bc
  !      integer                   :: g,i,ifan,ifln,ieln,Fin,ind_ufn(refElPol%Nfacenodes*phys%neq)
  !      real*8                    :: dline,xyDerNorm_g
  !      real*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
  !      real*8                    :: ufn(refElPol%Nfacenodes*phys%neq)
  !      real*8                    :: ufgn(refElPol%Ngauss1d,phys%neq),uaux(refElPol%Ngauss1d,phys%neq)

  !      ! Solution in the corresponding face
  !      ifan = Mesh%periodic_faces(ifa)
  !      ieln = Mesh%extfaces(ifan,1)
  !      ifln = Mesh%extfaces(ifan,2)
  !						Fin = ifan + Mesh%Nintfaces
  !						ind_ufn = (Fin - 1)*Neq*Npfl + (/(i,i=1,Neq*Npfl)/)
  !						ufn = sol%u_tilde(ind_ufn)
  !						uaux = matmul(refElPol%N1D,transpose(reshape(ufn,[neq,Npfl])))

  !      ! Flip solution in the corresponding face
  !      do i=1,refElPol%Ngauss1d
  !         ufgn(i,:) = uaux(refElPol%Ngauss1d+1-i,:)
  !      end do

  !      ! Loop in 1D Gauss points
  !      DO g = 1,Ng1d

  !         ! Calculate the integration weight
  !         xyDerNorm_g = norm2(xyDer(g,:))
  !         dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
  !         IF (switch%axisym) THEN
  !            dline = dline*xyg(g,1)
  !         END IF
  !         NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
  !         Ni = refElPol%N1D(g,:)*dline
  !         CALL assembly_dirichletwf_bc(iel,ind_asf,ind_ff,ufgn(g,:),NiNi,Ni)
  !      END DO

  !   END SUBROUTINE set_periodic_bc

  !***************************************************************
  ! Neumann in all variables
  !***************************************************************
  SUBROUTINE set_NeumannH_bc
    INTEGER                   :: g
    REAL*8                    :: dline,xyDerNorm_g
    REAL*8                    :: t_g(Ndim),n_g(Ndim)
    REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
    REAL                      :: tau_stab(Neq,Neq)
    ! Loop in 1D Gauss points
    DO g = 1,Ng1d

      ! Calculate the integration weight
       xyDerNorm_g = NORM2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
      IF (switch%axisym) THEN
        dline = dline*xyg(g,1)
      END IF

      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]


      !****************** Stabilization part******************************
      IF (numer%stab > 1) THEN
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
#ifndef KEQUATION
          CALL computeTauGaussPoints(upg(g,:),ufg(g,:),qfg(g,:),b(g,1:2),n_g,iel,ifa,1.,xyg(g,:),tau_stab)
#else
          CALL computeTauGaussPoints(upg(g,:),ufg(g,:),qfg(g,:),b(g,1:2),n_g,iel,ifa,1.,xyg(g,:),q_cyl(g),tau_stab)
#endif
        ELSE
          CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b(g,1:2),n_g,xyg(g,:),1.,iel,tau_stab)
        ENDIF
      ELSE
        tau_stab=0.
        DO i = 1,Neq
          tau_stab(i,i) = numer%tau(i)
        END DO
      END IF
      !****************** End stabilization part******************************


      NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Ni = refElPol%N1D(g,:)*dline

      CALL assembly_neum_bc(iel,ind_asf,ind_ash,ind_ff,ind_fg,NiNi,n_g,tau_stab)

    END DO

  END SUBROUTINE set_NeumannH_bc


  !***************************************************************
  ! Dirichlet in weak form in some variables and Neumann in others
  !***************************************************************
  SUBROUTINE set_diric_neum_bc
    INTEGER                   :: g
    REAL*8                    :: dline,xyDerNorm_g
    REAL*8                    :: t_g(Ndim),n_g(Ndim)
    REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
    REAL                      :: tau_stab(Neq,Neq)
    ! Loop in 1D Gauss points
    DO g = 1,Ng1d

      ! Calculate the integration weight
       xyDerNorm_g = NORM2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
      IF (switch%axisym) THEN
        dline = dline*xyg(g,1)
      END IF

      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]


      !****************** Stabilization part******************************
      IF (numer%stab > 1) THEN
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
#ifndef KEQUATION
          CALL computeTauGaussPoints(upg(g,:),ufg(g,:),qfg(g,:),b(g,1:2),n_g,iel,ifa,1.,xyg(g,:),tau_stab)
#else
          CALL computeTauGaussPoints(upg(g,:),ufg(g,:),qfg(g,:),b(g,1:2),n_g,iel,ifa,1.,xyg(g,:),q_cyl(g),tau_stab)
#endif
        ELSE
          CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b(g,1:2),n_g,xyg(g,:),1.,iel,tau_stab)
        ENDIF
      ELSE
        tau_stab=0.
        DO i = 1,Neq
          tau_stab(i,i) = numer%tau(i)
        END DO
      END IF
      !****************** End stabilization part******************************


      NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Ni = refElPol%N1D(g,:)*dline
      CALL assembly_diric_neum_bc(iel,ind_asf,ind_ash,ind_ff,ind_fg,uex(g,:),NiNi,Ni,b(g,1:2),n_g,tau_stab,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g))
    END DO

  END SUBROUTINE set_diric_neum_bc

  !****************************
  ! Dirichlet in strong form
  !****************************
  SUBROUTINE set_dirichletstr_bc
    INTEGER                   :: g
    REAL*8                    :: dline,xyDerNorm_g
    REAL*8                    :: t_g(Ndim),n_g(Ndim)
    REAL*8                    :: Ni(Npfl)

    ! Loop in 1D Gauss points
    DO g = 1,Ng1d

      ! Calculate the integration weight
       xyDerNorm_g = NORM2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
      IF (switch%axisym) THEN
        dline = dline*xyg(g,1)
      END IF

      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]

      ! Shape functions product
      Ni = refElPol%N1D(g,:)*dline

      CALL assembly_dirichletstr_bc(iel,ind_asf,ind_ash,ind_fe,ind_fg,Ni,qfg(g,:),uex(g,:),uexpg(g,:),b(g,1:2),n_g)

    END DO

  END SUBROUTINE set_dirichletstr_bc


  !***************************************************************
  ! Boundary condition to model the ITER core
  !***************************************************************
  SUBROUTINE set_itercore_bc
    INTEGER                   :: g
    REAL*8                    :: dline,xyDerNorm_g
    REAL*8                    :: t_g(Ndim),n_g(Ndim)
    REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
    REAL                      :: tau_stab(Neq,Neq)
    ! Loop in 1D Gauss points
    DO g = 1,Ng1d

      ! Calculate the integration weight
       xyDerNorm_g = NORM2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
      IF (switch%axisym) THEN
        dline = dline*xyg(g,1)
      END IF

      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]


      !****************** Stabilization part******************************
      IF (numer%stab > 1) THEN
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
#ifndef KEQUATION
          CALL computeTauGaussPoints(upg(g,:),ufg(g,:),qfg(g,:),b(g,1:2),n_g,iel,ifa,1.,xyg(g,:),tau_stab)
#else
          CALL computeTauGaussPoints(upg(g,:),ufg(g,:),qfg(g,:),b(g,1:2),n_g,iel,ifa,1.,xyg(g,:),q_cyl(g),tau_stab)
#endif
        ELSE
          CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b(g,1:2),n_g,xyg(g,:),1.,iel,tau_stab)
        ENDIF
      ELSE
        tau_stab=0.
        DO i = 1,Neq
          tau_stab(i,i) = numer%tau(i)
        END DO
      END IF
      !****************** End stabilization part******************************


      NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Ni = refElPol%N1D(g,:)*dline
      CALL assembly_itercore_bc(iel,ind_asf,ind_ash,ind_ff,ind_fg,uex(g,:),NiNi,Ni,b(g,1:2),n_g,tau_stab,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g))
    END DO

  END SUBROUTINE set_itercore_bc



  !********************************
  ! Transmission boundary condition
  !********************************
  SUBROUTINE set_Transmission_bc
    INTEGER                   :: g
    REAL*8                    :: dline,xyDerNorm_g
    REAL*8                    :: t_g(Ndim),n_g(Ndim)
    REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
    ! Loop in 1D Gauss points
    DO g = 1,Ng1d

      ! Calculate the integration weight
       xyDerNorm_g = NORM2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
      IF (switch%axisym) THEN
        dline = dline*xyg(g,1)
      END IF
      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]


      ! Shaspe functions product
      NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Ni = refElPol%N1D(g,:)*dline
      CALL assembly_trasm_bc(iel,ind_asf,ind_ff,NiNi,Ni,ufg(g,:))
    END DO

  END SUBROUTINE set_Transmission_bc


  !****************************
  ! Bohm
  !****************************
#ifndef SAVEFLUX
  SUBROUTINE set_Bohm_bc(v_nn_Bou_el,tau_save_el,xy_g_save_el)
#else
  SUBROUTINE set_Bohm_bc(v_nn_Bou_el,tau_save_el,xy_g_save_el,faceflux_pump,faceflux_puff,faceflux_parallel,faceflux_perpendicular,faceflux_neutral,faceflux_numerical)
#endif
      INTEGER                   :: g,i
      REAL*8                    :: dline,xyDerNorm_g
      REAL*8                    :: setval
      REAL*8                    :: t_g(Ndim),n_g(Ndim)
      REAL*8                    :: inc,sn
      INTEGER                   :: delta
      REAL*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
      REAL*8                    :: soundSpeed
    REAL*8                    :: dcs_du(Neq)
      REAL*8                    :: bn
      LOGICAL                   :: ntang
      REAL*8,INTENT(out)        :: v_nn_Bou_el(:,:),tau_save_el(:,:),xy_g_save_el(:,:)
      REAL*8                    :: Vnng(Ndim)
      REAL                      :: tau_stab(Neq,Neq)
#ifdef SAVEFLUX
    real*8,intent(out)        :: faceflux_pump, faceflux_puff,faceflux_parallel,faceflux_perpendicular,faceflux_neutral,faceflux_numerical
    real*8                    :: flgflux_pump, flgflux_puff,flgflux_parallel,flgflux_perpendicular,flgflux_neutral,flgflux_numerical

    faceflux_pump = 0.
    faceflux_puff = 0.
    faceflux_parallel = 0.
    faceflux_perpendicular = 0.
    faceflux_neutral = 0.
    faceflux_numerical = 0.
#endif
    ! Loop in 1D Gauss points
    DO g = 1,Ng1d

      ! Calculate the integration weight
         xyDerNorm_g = NORM2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
      IF (switch%axisym) THEN
        dline = dline*xyg(g,1)
      END IF

      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]
         bn = dot_PRODUCT(b(g,1:2),n_g)
      ! Sound speed
#ifndef TEMPERATURE
         SoundSpeed = SQRT(phys%a)
#else

         !SoundSpeed = upg(g,9)
         CALL compute_cs(ufg(g,:),SoundSpeed)
         call compute_dcs_du(ufg(g,:),dcs_du)
         IF (switch%decoup) THEN
            SoundSpeed = SQRT(phys%Mref)
         ENDIF
#endif
      ! tangency
      ntang = .TRUE.
#ifdef BOHMLIMIT
      !use this flag as a proxy for turning off/on the bohm limit technique
      if (any(ufg(:,4)<2.e-7)) then
        ntang = .false.
      endif
#endif
      inc = bn/norm2(b(g,1:2))

#ifdef NGAMMA
      setval = ufg(g,2)

      !IF (abs(inc) .le. phys%bohmth) THEN
      !   !sn = 1.
      !   !delta = 0 ! TODO try to module delta in the WEST case
      !   !Modification BC
      !   sn = sign(1.,bn)
      !   delta = 0.
      !   setval = sn
      !   !End modification BC
      !   ntang = .FALSE.
      !ELSE
      !  sn = sign(1.,bn)
      !  delta = 1
      !  IF ((abs(upg(g,2))) .le. SoundSpeed) THEN
      !     setval = sn*SoundSpeed
      !     delta = 1
      !  ELSE IF ((abs(upg(g,2))) .gt. SoundSpeed) THEN
      !    !delta = 0
      !    !setval = upg(g,2)
      !    !Modification BC
      !    delta = 0.
      !    setval = sn
      !    !End modification BC
      !  END IF

      ! NEW BOHM BC
      ! Above angle threshold: Impose outgoing supersonic velocity everywhere
      ! Below angle thresgold: Impose linear decrease of supersonic velocity
      !                        to 0 as a function of the angle of incidence
         setval = ABS(upg(g,2))
         sn = SIGN(1.,bn)
      delta = 1.
         IF (ABS(inc) .LE. phys%bohmth) THEN
            setval = sn*SoundSpeed/phys%bohmth*ABS(inc)
            dcs_du = sn*dcs_du/phys%bohmth*ABS(inc)
      ELSE
            IF (ABS(upg(g,2)) .LE. SoundSpeed) THEN
            setval = sn*SoundSpeed
               dcs_du = sn*dcs_du
            ELSE IF (ABS(upg(g,2)) .GT. SoundSpeed) THEN
            delta = 0.
            !setval = sn*setval
         END IF
         !if (numer%bohmtypebc.eq.1) then
         !  delta=1
         !endif
      END IF
#endif

      ! Shape functions
      NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Ni = refElPol%N1D(g,:)*dline
      !****************** Stabilization part******************************
      IF (numer%stab > 1) THEN
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
#ifndef KEQUATION
          CALL computeTauGaussPoints(upg(g,:),ufg(g,:),qfg(g,:),b(g,1:2),n_g,iel,ifa,1.,xyg(g,:),tau_stab)
#else
          CALL computeTauGaussPoints(upg(g,:),ufg(g,:),qfg(g,:),b(g,1:2),n_g,iel,ifa,1.,xyg(g,:),q_cyl(g),tau_stab)
#endif
        ELSE
          CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b(g,1:2),n_g,xyg(g,:),1.,iel,tau_stab)
        ENDIF

            IF (save_tau) THEN
          DO i = 1,Neq
            tau_save_el(g,i) = tau_stab(i,i)
          END DO
          xy_g_save_el(g,:) = xyg(g,:)
            ENDIF
      ELSE
        tau_stab=0.
        DO i = 1,Neq
          tau_stab(i,i) = numer%tau(i)
        END DO
      END IF
      !****************** End stabilization part******************************


      ! Assembly Bohm contribution
         IF (numer%bohmtypebc.EQ.0) THEN
#ifndef SAVEFLUX
#ifndef DKLINEARIZED
        CALL assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
          &ufg(g,:),upg(g,:),ueg(g,:),b(g,1:2),psig(g),n_g,tau_stab,setval,dcs_du,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang,Vnng)
#else
        CALL assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
          &ufg(g,:),upg(g,:),ueg(g,:),b(g,1:2),psig(g),q_cyl(g),xyg(g,:),n_g,tau_stab,setval,dcs_du,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang,Vnng)
#endif
#else
#ifndef DKLINEARIZED
        CALL assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
          &ufg(g,:),upg(g,:),ueg(g,:),b(g,1:2),psig(g),n_g,tau_stab,setval,dcs_du,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),dline,ntang,Vnng,flgflux_pump,flgflux_puff,flgflux_parallel,flgflux_perpendicular,flgflux_neutral,flgflux_numerical)
#else
        CALL assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
        &ufg(g,:),upg(g,:),ueg(g,:),b(g,1:2),psig(g),q_cyl(g),xyg(g,:),n_g,tau_stab,setval,dcs_du,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),dline,ntang,Vnng,flgflux_pump,flgflux_puff,flgflux_parallel,flgflux_perpendicular,flgflux_neutral,flgflux_numerical)
#endif
        !summing conribution from each part of the face
        faceflux_pump = faceflux_pump+flgflux_pump
        faceflux_puff = faceflux_puff+flgflux_puff
        faceflux_parallel = faceflux_parallel+flgflux_parallel
        faceflux_perpendicular = faceflux_perpendicular+flgflux_perpendicular
        faceflux_neutral = faceflux_neutral+flgflux_neutral
            faceflux_numerical = faceflux_numerical+flgflux_numerical
#endif
         ELSE
        CALL assembly_bohm_bc_new(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
          &ufg(g,:),upg(g,:),b(g,1:2),n_g,tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang)
         ENDIF

         IF (save_tau) THEN
         v_nn_Bou_el(g,:) = Vnng
         ENDIF

    END DO ! END loop in Gauss points

  END SUBROUTINE set_Bohm_bc

#endif
!TOR3D

  !*******************************************
  !   AUXILIARY ROUTINES FOR 2D AND 3D
  !*******************************************

  !*********************************
  ! Assembly Dirichlet in weak form
  !*********************************
!#ifdef NEUTRAL
!  SUBROUTINE assembly_dirichletwf_bc(iel,ind_asf,ind_ff,ufg,NiNi,Ni,xyfg)
!    real*8       :: xyfg(:)
!    real*8       :: kmult_neutral(Npfl)
!    real*8       :: p1,p2,p3,p4,p5,p6,p7,p8,p9,a1,a2,a3,b1,b2,b3,c1,c2,c3,fit
!    integer*4    :: iel,ind_asf(:),ind_ff(:)
!    real*8       :: ufg(:)
!    real*8       :: NiNi(:,:),Ni(:)
!    real*8       :: kmult(Neq*Npfl)
!    integer*4    :: ind(Npfl),i
!
!    !xyfg(:) = xyfg(:)*phys%lscale
!    ! Dirichlet for neutrals
!    !i=Neq
!    !ind = i+ind_asf
!    if ((xyfg(1)*phys%lscale) .le. 3.4 ) then
!      a1 = 0.002492
!      b1 = - 0.7678
!      c1 = 0.03091
!      a2 = 0.004729
!      b2 = - 0.7951
!      c2 = 0.07134
!      a3 = 0.009261
!      b3 = - 0.8927
!      c3 = 0.3785
!    else
!      a1 = 0.002138;
!      b1 = - 0.7698;
!      c1 = 0.02522;
!      a2 = 0.004928;
!      b2 = -0.7942;
!      c2 = 0.06691;
!      a3 = 0.009552;
!      b3 = -0.8594;
!      c3 = 0.3858;
!    end if
!    fit = a1*exp(-((xyfg(2)-b1)/c1)**2) + a2*exp(-((xyfg(2)-b2)/c2)**2) + a3*exp(-((xyfg(2)-b3)/c3)**2)
!    !fit = p1*xyfg(2)**8 + p2*xyfg(2)**7 + p3*xyfg(2)**6 + p4*xyfg(2)**5 + p5*xyfg(2)**4&
!    !       &+ p6*xyfg(2)**3 + p7*xyfg(2)**2 + p8*xyfg(2) + p9
!    kmult_neutral = (fit*Ni)
!    !elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult_neutral
!    kmult = col(tensorProduct(ufg(:),Ni))
!    DO i = 1,Neq
!      ind = i + ind_asf
!      elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
!      if (i.eq.Neq) then
!        elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult_neutral
!      else
!        elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
!      endif
!      !         elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - NiNi
!      !         elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - kmult(ind)
!    END DO
!
!  END SUBROUTINE assembly_dirichletwf_bc
!
!#else

  SUBROUTINE assembly_dirichletwf_bc(iel,ind_asf,ind_ff,ufg,NiNi,Ni)
      INTEGER*4    :: iel,ind_asf(:),ind_ff(:)
      REAL*8       :: ufg(:)
      REAL*8       :: NiNi(:,:),Ni(:)
      REAL*8       :: kmult(Neq*Npfl)
      INTEGER*4    :: ind(Npfl),i

    kmult = col(tensorProduct(ufg(:),Ni))
    DO i = 1,Neq
      ind = i + ind_asf
         elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
      elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
      !         elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - NiNi
      !         elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - kmult(ind)
    END DO

  END SUBROUTINE assembly_dirichletwf_bc

!#endif
!NEUTRAL

  !   !*********************************
  !   ! Assembly periodic bc
  !   !*********************************
  !   SUBROUTINE assembly_periodic_bc(iel,ind_ff,ind_ff,NiNi)
  !      integer*4    :: iel,ind_ff(:)
  !      real*8       :: NiNi(:,:)
  !      real*8       :: kmult(Neq*Npfl)
  !      integer*4    :: ind(Npfl),i

  !      kmult = col(tensorProduct(ufg(:),Ni))
  !      DO i = 1,Neq
  !         ind = i + ind_asf
  !         elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
  !         elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
  !      END DO

  !   END SUBROUTINE assembly_periodic_bc



  !*********************************
  ! Trasmissive conditions
  !*********************************
  SUBROUTINE assembly_trasm_bc(iel,ind_asf,ind_ff,NiNi,Ni,ug)
      INTEGER*4    :: iel,ind_asf(:),ind_ff(:)
      REAL*8       :: NiNi(:,:),Ni(:),ug(:)
      REAL*8       :: kmult(Neq*Npfl)
      INTEGER*4    :: i
      INTEGER*4    :: ind(Npfl)

    kmult = col(tensorProduct(ug(:),Ni))
    DO i = 1,Neq
      ind = i + ind_asf
         elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
      elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
    END DO

  END SUBROUTINE assembly_trasm_bc




  !**************************
  ! Assembly  Neumann for all
  !**************************
  SUBROUTINE assembly_neum_bc(iel,ind_asf,ind_ash,ind_ff,ind_fg,NiNi,ng,tau)
      INTEGER*4    :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fg(:)
      REAL*8       :: NiNi(:,:),ng(:)
      REAL*8       :: tau(:,:)
      INTEGER*4    :: i,idm
      INTEGER*4    :: ind(Npfl),indi(Npfl),indj(Npfl)

    DO i = 1,Neq
      ind = i + ind_asf
         elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - tau(i,i)*NiNi
      elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*NiNi
      DO idm = 1,Ndim
        indi = ind_asf + i
        indj = ind_ash + idm + (i - 1)*Ndim
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + NiNi*ng(idm)
      END DO
    END DO

  END SUBROUTINE assembly_neum_bc


  !*********************************
  ! Assembly Dirichlet and Neumann
  !*********************************
  SUBROUTINE assembly_diric_neum_bc(iel,ind_asf,ind_ash,ind_ff,ind_fg,ufg,NiNi,Ni,bg,ng,tau,diffiso,diffani)
      INTEGER*4    :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fg(:)
      REAL*8       :: NiNi(:,:),Ni(:),ufg(:),bg(:),ng(:)
      REAL*8       :: tau(:,:)
      REAL*8       :: diffiso(:,:),diffani(:,:)
      REAL*8       :: kmult(Neq*Npfl),bn
      INTEGER*4    :: i,j, idm
      INTEGER*4    :: ind(Npfl),indi(Npfl),indj(Npfl)

    !      ndir = 2
    !      ! Dirichlet part
    !      kmult= col(tensorProduct(uexg(:),Ni))
    !      DO i=1,Neq
    !         ind = i+ind_asf
    !         elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
    !         if (i.le.ndir) then
    !            elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
    !         else
    !            elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + numer%tau(i)*NiNi
    !         endif
    !      END DO
    !
    !      ! Neumann part
    !                           DO k = ndir+1,Neq
    !                              DO idm = 1,Ndim
    !                                 indi = ind_asf+k
    !            indj = ind_ash+idm+(k-1)*Ndim
    !            elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)+NiNi*n(idm)
    !                              END DO
    !                                                END DO

#ifndef VORTICITY
      WRITE(6,*) "Wrong bc! Should be used only for vorticity"
      STOP
#endif
      bn = dot_PRODUCT(bg,ng)

    ! Dirichlet values
    kmult = col(tensorProduct(ufg(:),Ni))

    ! Dirichlet weak form for the first equation
    i = 1
    ind = i + ind_asf
      elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
    elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)

    ! Dirichlet weak form for the second equation
    i = 2
    ind = i + ind_asf
      elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
    elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)

      IF (switch%dirivortcore) THEN
      ! Dirichlet weak form for the vorticity equation
      i = 3
      ind = i + ind_asf
         elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
      elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
      ELSE
      ! Neumann homogeneous for the vorticity
      i = 3
      ind = i + ind_asf
         elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - tau(i,i)*NiNi
      elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*NiNi
      !write(6,*) "tau(i,i)",tau(i,i)
      DO idm = 1,Ndim
        indi = ind_asf + i
        indj = ind_ash + idm + (i - 1)*Ndim
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + &
          &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i))

        !elMat%Alq(ind_ff(indi), ind_fG(indj), iel) = elMat%Alq(ind_ff(indi), ind_fG(indj), iel) + NiNi*ng(idm)
      END DO
      ENDIF

    !if (iel==1) then
    !call hdf5_save_matrix( elMat%Alq(:,:,iel), 'Alq')
    !stop
    !endif
    ! Neumann homogeneous for the potential
    i = 4
    ind = i + ind_asf
      elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - tau(i,i)*NiNi
    elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*NiNi

    !elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - 0.1*NiNi
    !elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + 0.1*NiNi

    DO idm = 1,Ndim
      indi = ind_asf + i
      indj = ind_ash + idm + (i - 1)*Ndim
      !elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + &
      !  &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) )

      elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + &
        &NiNi*ng(idm)

    END DO
    j = 1
    DO idm = 1,Ndim
      indi = ind_asf + i
      indj = ind_ash + idm + (j - 1)*Ndim
      elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + NiNi*ng(idm)/ufg(1)*phys%Mref
    END DO

    !if (iel==1) then
    !call hdf5_save_matrix( elMat%Alq(:,:,iel), 'Alq')
    !stop
    !endif

  END SUBROUTINE assembly_diric_neum_bc



  !*********************************
  ! ITER core boundary condition
  !*********************************
  SUBROUTINE assembly_itercore_bc(iel,ind_asf,ind_ash,ind_ff,ind_fg,ufg,NiNi,Ni,bg,ng,tau,diffiso,diffani)
      INTEGER*4    :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fg(:)
      REAL*8       :: NiNi(:,:),Ni(:),ufg(:),bg(:),ng(:)
      REAL*8       :: tau(:,:)
      REAL*8       :: diffiso(:,:),diffani(:,:)
      REAL*8       :: bn
      INTEGER*4    :: i,j, idm
      INTEGER*4    :: indi(Npfl),indj(Npfl)
      REAL*8       :: source_dens_coeff,source_ener_coeff


#ifndef NEUTRAL
      WRITE(6,*) "Wrong bc! Should be used only with neutrals"
      STOP
#endif


      bn = dot_PRODUCT(bg,ng)

    source_dens_coeff = phys%part_source/simpar%refval_density/(Mesh%core_area*phys%lscale**2)/(simpar%refval_diffusion)*phys%lscale
    source_ener_coeff = 0.5*phys%ener_source/simpar%refval_specenergydens/(Mesh%core_area*phys%lscale**2)/(simpar%refval_diffusion)*phys%lscale/simpar%refval_mass

    ! *********** Equation 1 --> Flux plasma equal to flux neutrals + source
    ! stabilization part
    i = 1
    indi = i + ind_asf
      elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
    elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
    ! flux diagonal part
    j = 1
    DO idm = 1,Ndim
      indj = ind_ash + idm + (j - 1)*Ndim
      !elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - &
      !  &NiNi*ng(idm)*phys%diff_n
      elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) -  &
        NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i))
    END DO
   ! flux non-diagonal part
    j = phys%Neq
    DO idm = 1,Ndim
      indj = ind_ash + idm + (j - 1)*Ndim
      elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - &
        &NiNi*ng(idm)*phys%diff_nn
    END DO
    ! flux rhs-part
    elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - source_dens_coeff*Ni



    ! *********** Equation 2 --> zero momentum
    ! stabilization part
    i = 2
    indi = i + ind_asf
      elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - numer%tau(i)*NiNi
     ! flux diagonal part
    j = 2
    DO idm = 1,Ndim
      indj = ind_ash + idm + (j - 1)*Ndim
      !elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*ng(idm)*phys%diff_u
      elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i))
    END DO


#ifdef TEMPERATURE
    ! *********** Equation 3 --> Ions energy flux equal 10 MW
    ! stabilization part
    i = 3
    indi = i + ind_asf
      elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
    elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
    ! flux diagonal part
    j = 3
    DO idm = 1,Ndim
      indj = ind_ash + idm + (j - 1)*Ndim
      !elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - &
      !  &NiNi*ng(idm)*phys%diff_e
      elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i))
    END DO
    ! flux rhs-part
    elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - source_ener_coeff*Ni


    ! *********** Equation 4 --> Electrons energy flux equal 10 MW
    ! stabilization part
    i = 4
    indi = i + ind_asf
      elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
    elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
    ! flux diagonal part
    j = 4
    DO idm = 1,Ndim
      indj = ind_ash + idm + (j - 1)*Ndim
      !elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - &
      !  &NiNi*ng(idm)*phys%diff_ee
      elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i))
    END DO
    ! flux rhs-part
    elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - source_ener_coeff*Ni
#endif


    ! *********** Equation Neq --> Flux neutrals equal to flux plasma - source
    ! stabilization part
    i = phys%Neq
    indi = i + ind_asf
      elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
    elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
    ! flux diagonal part
    j = phys%Neq
    DO idm = 1,Ndim
      indj = ind_ash + idm + (j - 1)*Ndim
      elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - &
        &NiNi*ng(idm)*phys%diff_nn
    END DO
!   ! flux non-diagonal part
!    j = 1
!    DO idm = 1,Ndim
!      indj = ind_ash + idm + (j - 1)*Ndim
!      elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - &
!        &NiNi*ng(idm)*phys%diff_n
!    END DO
!    ! flux rhs-part
!    elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - source_dens_coeff*Ni




  END SUBROUTINE assembly_itercore_bc



  !**********************************
  ! Assembly Dirichlet in strong form
  !**********************************
  SUBROUTINE assembly_dirichletstr_bc(iel,ind_asf,ind_ash,ind_fe,ind_fg,Ni,qg,ufg,upg,bg,ng)
      INTEGER*4        :: iel,ind_asf(:),ind_ash(:),ind_fe(:),ind_fg(:)
      REAL*8           :: Ni(:),ufg(:),upg(:),bg(:),ng(:)
      REAL*8           :: qg(:)
      REAL*8           :: bn,An(Neq,Neq)
      INTEGER          :: i,j,idm
      REAL*8           :: kmultstab(Neq*Npfl),kmultconv(Neq*Npfl),kmultfseq(Ndim*Neq*Npfl)
      REAL*8           :: uexn(Ndim*neq)
      INTEGER*4        :: ind(Npfl),ind_1f(Npfl)
      REAL*8           :: Qpr(Ndim,Neq)
#ifndef TEMPERATURE
      REAL*8           :: auxvec(Neq)
#else
      REAL*8           :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
      REAL*8           :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)
#endif

      bn = dot_PRODUCT(bg,ng)

    CALL jacobianMatricesFace(ufg,bn,An)
    ! Contribution of the current integration point to the elemental matrix
    kmultstab = col(tensorProduct(ufg,Ni))
      kmultconv = col(tensorProduct(MATMUL(An,ufg),Ni))
#ifndef TEMPERATURE
      IF (switch%logrho) THEN
         CALL logrhojacobianVector(ufg,upg,auxvec)
      kmultconv = kmultconv +col(tensorProduct(auxvec,bn*Ni))
      ENDIF
#endif
      uexn = RESHAPE(tensorProduct(ng,ufg),(/Ndim*neq/)) ! Ndim*neq x 1
    kmultfseq = col(tensorProduct(uexn,Ni))  ! Ndim*neq*npf
    DO i = 1,Neq
      ind = i + ind_asf
      elMat%Aul_dir(ind_fe(ind),iel) = elMat%Aul_dir(ind_fe(ind),iel) - numer%tau(i)*kmultstab(ind) + kmultconv(ind)
      DO idm = 1,Ndim
        ind_1f = ind_ash + idm + (i - 1)*Ndim
        elMat%Aql_dir(ind_fG(ind_1f),iel) = elMat%Aql_dir(ind_fG(ind_1f),iel) - kmultfseq(ind_1f)
      END DO
    END DO

    ! Compute Q^T^(k-1)
      Qpr = RESHAPE(qg,(/Ndim,Neq/))

#ifdef TEMPERATURE
    ! Compute V(U^(k-1))
      CALL computeVi(ufg,Vveci)
      CALL computeVe(ufg,Vvece)

    ! Compute dV_dU (k-1)
      CALL compute_dV_dUi(ufg,dV_dUi)
      CALL compute_dV_dUe(ufg,dV_dUe)

    ! Compute Gamma and tau
      gmi = dot_PRODUCT(MATMUL(Qpr,Vveci),bg)   ! scalar
      gme = dot_PRODUCT(MATMUL(Qpr,Vvece),bg)
      Taui = MATMUL(Qpr,dV_dUi)                       ! 2x3
      Taue = MATMUL(Qpr,dV_dUe)

    ! Compute Alpha(U^(k-1))
    Alphai = computeAlphai(ufg)
    Alphae = computeAlphae(ufg)

    ! Compute dAlpha/dU^(k-1)
      CALL compute_dAlpha_dUi(ufg,dAlpha_dUi)
      CALL compute_dAlpha_dUe(ufg,dAlpha_dUe)

    DO i = 3,4
      DO j = 1,4
        IF (i == 3) THEN
          elMat%Aul_dir(ind_fe(i + ind_asf),iel) = elMat%Aul_dir(ind_fe(i + ind_asf),iel) - &
                    &coefi*Ni*bn*((gmi*dAlpha_dUi(j) + Alphai*(dot_PRODUCT(Taui(:,j),bg)))*ufg(j))
        ELSE
          elMat%Aul_dir(ind_fe(i + ind_asf),iel) = elMat%Aul_dir(ind_fe(i + ind_asf),iel) - &
                    &coefe*Ni*bn*((gme*dAlpha_dUe(j) + Alphae*(dot_PRODUCT(Taue(:,j),bg)))*ufg(j))
        END IF
      END DO
    END DO
#endif

  END SUBROUTINE assembly_dirichletstr_bc


  !*********************************
  ! Assembly Bohm
  !*********************************
#ifndef SAVEFLUX
#ifndef DKLINEARIZED
    SUBROUTINE assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg,ufg,upfg,uefg,bg,psig,ng,tau,setval,dcs_du,delta,diffiso,diffani,ntang,Vnng)
#else
    SUBROUTINE assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg,ufg,upfg,uefg,bg,psig,q_cyl,xyf,ng,tau,setval,dcs_du,delta,diffiso,diffani,ntang,Vnng)
#endif
#else
#ifndef DKLINEARIZED
    SUBROUTINE assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg,ufg,upfg,uefg,bg,psig,ng,tau,setval,dcs_du,delta,diffiso,diffani,dline,ntang,Vnng,flgflux_pump,flgflux_puff,flgflux_parallel,flgflux_perpendicular,flgflux_neutral,flgflux_numerical)
#else
    SUBROUTINE assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg,ufg,upfg,uefg,bg,psig,q_cyl,xyf,ng,tau,setval,dcs_du,delta,diffiso,diffani,dline,ntang,Vnng,flgflux_pump,flgflux_puff,flgflux_parallel,flgflux_perpendicular,flgflux_neutral,flgflux_numerical)
#endif
#endif
    integer*4        :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:),bc,delta
    real*8           :: NiNi(:,:),Ni(:),ufg(:),upfg(:),uefg(:),bg(:),psig,ng(:),tau(:,:),setval,dcs_du(:)
    real*8           :: diffiso(:,:),diffani(:,:)
    logical          :: ntang
    real*8           :: qfg(:),Vnng(:)
    real*8           :: bn,Abohm(Neq,Neq),APinch(Neq,Ndim)
    integer          :: i,j,k,idm,Neqstab,Neqgrad
    integer*4        :: ind(Npfl),indi(Npfl),indj(Npfl),indk(Npfl),ind_jf(Npfl),ind_kf(Npfl)
    real*8           :: Qpr(Ndim,Neq),kcoeff, recycling_coeff,  cryopump_coeff,puff_coeff
    real*8           :: W2(Neq), dW2_dU(Neq,Neq), QdW2(Ndim,Neq)
    real*8           :: kmult(Npfl,Npfl),kmultf(Npfl)
#ifdef TEMPERATURE
        REAL*8           :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
        REAL*8           :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)
        REAL*8           :: W3(Neq), dW3_dU(Neq,Neq), QdW3(Ndim,Neq)
        REAL*8           :: W4(Neq), dW4_dU(Neq,Neq), QdW4(Ndim,Neq)
#ifdef DNNLINEARIZED
    real*8                    :: Dnn_dU(Neq), Dnn_dU_U
    real*8                    :: gradDnn(Ndim)
#endif
#ifdef DKLINEARIZED
    real*8                 ::       q_cyl, xyf(:), ddk_dU(Neq), ddk_dU_u
#endif
#ifdef NEUTRALP
        REAL*8           :: Dnn,Dpn,GammaLim,Alphanp,Betanp,Gammaredpn,Tmin
        REAL*8           :: AbohmNP(Neq),Vpn(Neq),dVpn_dU(Neq,Neq),gmpn(Ndim),gmipn(Ndim),Taupn(Ndim,Neq),dDpn_dU(Neq)
#endif
#endif
#ifdef SAVEFLUX
    real*8           :: dline
    real*8,intent(out)::  flgflux_pump,flgflux_puff,flgflux_parallel,flgflux_perpendicular,flgflux_neutral,flgflux_numerical
#endif
#ifdef BOHMLIMIT
    real*8           :: U3_min = 2.e-8 ! 1e16[m^-3]*0.01^[eV]/n0/T0
    logical          :: bohm_limit = .true.    ! false if U3<U3_min

    if (ufg(3)<U3_min) then
      bohm_limit = .false.
    endif

#endif

    Neqstab = Neq
    Neqgrad = Neq
#ifdef VORTICITY
    Neqgrad = 2
    IF (switch%fixdPotLim) THEN
      Neqstab = 3
    ENDIF
        IF (.NOT. (ntang)) THEN
      tau = 1.
    ENDIF
#endif

    !if (iel==195) then
    !call hdf5_save_matrix( elMat%Alu(:,:,iel),"Alu")
    !stop
    !endif

    !tau=1.
        bn = dot_PRODUCT(bg,ng)

    ! Compute Q^T^(k-1)
        Qpr = RESHAPE(qfg,(/Ndim,Neq/))

    ! Split diffusion matrices/vectors for the momentum equation
    CALL compute_W2(uf,W2,diffiso(1,1),diffiso(2,2))
    CALL compute_dW2_dU(uf,dW2_dU,diffiso(1,1),diffiso(2,2))
        QdW2 = MATMUL(Qpr,dW2_dU)

    ! Case diagonal matrix stabilization
    IF (numer%stab < 6) THEN
      DO i = 1,Neqstab
        ind = i + ind_asf
              elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - tau(i,i)*NiNi
              !dirichlet boundary on minimal allowed values
              !IF ((i .ne. 5) .and. ((ufg(1)<1e-5) .or. (ufg(3)<1.e-8*3/2*phys%Mref) .or. (ufg(4)<1.e-8*3/2*phys%Mref))) THEN
              !IF (.false.) THEN
              ! if (i==1) then
              !    elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel)  - tau(i,i)*1e-5*Ni
              ! elseif (i == 3) then
              !    elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel)  - tau(i,i)*1.e-8*3/2*phys%Mref*Ni
              ! elseif (i==4) then
              !    elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel)  - tau(i,i)*1.e-8*3/2*phys%Mref*Ni
              ! endif
        IF (i == 2) THEN
          !elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) + tau(i,i)*(delta*setval)*NiNi
          !if (iel==195) then
          !write(6,*) "delta",delta
          !write(6,*) "setval",setval
          !write(6,*) "tau(i,i)",tau(i,i)
          !endif
          !elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*(1 - delta)*NiNi
          !Modification BC
          !elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*(1 - delta)*setval*NiNi
          !End modification BC
          !NEW BOHM BC: Impose everywhere a value for velocity
          elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) + tau(i,i)*delta*setval*NiNi
          elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*(1. - delta)*NiNi
                 !elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*NiNi
                 !weird investigations
                 elMat%ALL(ind_ff(ind),ind_ff(ind),iel) = elMat%ALL(ind_ff(ind),ind_ff(ind),iel) - tau(i,i)*NiNi
                 elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*NiNi
                 !linearization of sound speed
                 DO j = 1,Neqstab
                  indj = j + ind_asf
                  elMat%Alu(ind_ff(ind),ind_fe(indj),iel) = elMat%Alu(ind_ff(ind),ind_fe(indj),iel) + tau(i,i)*delta*dcs_du(j)*ufg(1)*NiNi
                 ENDDO
        ELSE
          elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*NiNi
        END IF
      END DO
    ELSE
      ! Case full matrix stabilization
      DO i = 1,Neqstab
        indi = i + ind_asf
        DO j = 1,Neqstab
          indj = j + ind_asf
                 elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) - tau(i,j)*NiNi
          IF (j == 1) THEN
            elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) +&
              &tau(i,j)*NiNi + tau(i,j + 1)*(delta*setval)*NiNi
          ELSE IF (j == 2) THEN
            elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) + tau(i,j)*(1 - delta)*NiNi
          ELSE
            elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) + tau(i,j)*NiNi
          ENDIF
        END DO
      END DO
    ENDIF


#ifdef TEMPERATURE

    !IF (ntang) THEN

      ! Compute V(U^(k-1))
        CALL computeVi(ufg,Vveci)
        CALL computeVe(ufg,Vvece)

      ! Compute dV_dU (k-1)
        CALL compute_dV_dUi(ufg,dV_dUi)
        CALL compute_dV_dUe(ufg,dV_dUe)

      ! Split diffusion matrices/vectors for the energies equations
      CALL compute_W3(uf,W3,diffiso(1,1),diffiso(2,2),diffiso(3,3))
      CALL compute_dW3_dU(uf,dW3_dU,diffiso(1,1),diffiso(2,2),diffiso(3,3))
        QdW3 = MATMUL(Qpr,dW3_dU)

      CALL compute_W4(uf,W4,diffiso(1,1),diffiso(4,4))
      CALL compute_dW4_dU(uf,dW4_dU,diffiso(1,1),diffiso(4,4))
        QdW4 = MATMUL(Qpr,dW4_dU)

      ! Compute Alpha(U^(k-1))
      Alphai = computeAlphai(ufg)
      Alphae = computeAlphae(ufg)

      ! Compute dAlpha/dU^(k-1)
        CALL compute_dAlpha_dUi(ufg,dAlpha_dUi)
        CALL compute_dAlpha_dUe(ufg,dAlpha_dUe)

      ! Jacobian matrix for convection part
      CALL jacobianMatricesBohm(ufg,Abohm)

      ! Jacobian matrix for pinch part
      CALL computePinch(bg,psig,APinch)

        gmi = dot_PRODUCT(MATMUL(Qpr,Vveci),bg)  ! scalar
        gme = dot_PRODUCT(MATMUL(Qpr,Vvece),bg)             ! scalar
        Taui = MATMUL(Qpr,dV_dUi)                      ! 2x3
        Taue = MATMUL(Qpr,dV_dUe)      ! 2x3
#ifdef DNNLINEARIZED
      call compute_Dnn_dU(ufg,Dnn_dU)
      Dnn_dU_u = dot_product(Dnn_dU,ufg)
#endif
#ifdef KEQUATION
#ifdef DKLINEARIZED
      call compute_ddk_dU(ufg,xyf,q_cyl,ddk_dU)
      ddk_dU_u = dot_product(ddk_dU,ufg)
#endif
#endif
#ifdef BOHMLIMIT
    IF (ntang) then
#endif
      ! Parallel diffusion for temperature
      !if (.not. ((ufg(1)<1e-5) .or. (ufg(3)<1.e-8*3/2*phys%Mref) .or. (ufg(4)<1.e-8*3/2*phys%Mref))) then
        !IF (.true.) THEN
      DO i = 1,2
        indi = ind_asf + i
        IF (i == 1) THEN
          DO j = 1,Neq
            indj = ind_asf + j
                 elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) + Abohm(i,j)*NiNi*bn
                 elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) -&
                      &coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_PRODUCT(Taui(:,j),bg)))*NiNi*bn
            DO k = 1,Ndim
              elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                &coefi*Alphai*Vveci(j)*bg(k)*NiNi*bn
            END DO
          END DO
              elMat%fh(ind_ff(indi+2),iel) = elMat%fh(ind_ff(indi+2),iel) - coefi*Alphai*( dot_PRODUCT (MATMUL(TRANSPOSE(Taui),bg),ufg)  )*Ni*bn
        ELSE
          DO j = 1,Neq
            indj = ind_asf + j
                 elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) + Abohm(i,j)*NiNi*bn
                 elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) -&
                      &coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_PRODUCT(Taue(:,j),bg)))*NiNi*bn
            DO k = 1,Ndim
              elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                &coefe*Alphae*Vvece(j)*bg(k)*NiNi*bn
            END DO
          END DO
              elMat%fh(ind_ff(indi+2),iel) = elMat%fh(ind_ff(indi+2),iel) - coefe*Alphae*( dot_PRODUCT (MATMUL(TRANSPOSE(Taue),bg),ufg)  )*Ni*bn
        ENDIF
      END DO
      !endif
    !END IF ! tangency
#endif

    ! Perpendicular diffusion
#ifdef BOHMLIMIT
    IF (ntang)  then
#else
    IF (ntang) THEN
#endif

      DO k = 1,Neqgrad
#ifdef NEUTRAL
              IF (k .EQ. 5) CYCLE
#endif
        DO idm = 1,Ndim
          indi = ind_asf + k
          indj = ind_ash + idm + (k - 1)*Ndim
          ! Non-tangent case
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
            &NiNi*(ng(idm)*diffiso(k,k)-bn*bg(idm)*diffani(k,k) )
            IF(k == 2) THEN
              ! Assembly LQ
              j = 1 ! all other terms are 0 anyway in vector W2
              ind_kf = ind_ash + idm + (j - 1)*Ndim
               kmult = NiNi*W2(j)*(ng(idm) - bn*bg(idm))
               elMat%Alq(ind_ff(indi),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(indi),ind_fG(ind_kf),iel) - kmult
            ! Assembly LU
               DO j=1,Neq
                 ind_jf = ind_asf+j
                   kmult = QdW2(idm,j)*NiNi*(ng(idm)-bn*bg(idm))
                       elMat%ALL(ind_ff(indi),ind_ff(ind_jf),iel)  = elMat%ALL(ind_ff(indi),ind_ff(ind_jf),iel) - kmult
               END DO
#ifdef TEMPERATURE
            ELSEIF (k ==3) THEN
                ! Assembly LQ
                DO j=1,Neq ! here there are non-zero values in vector W3
                ind_kf = ind_ash + idm + (j - 1)*Ndim
                 kmult = NiNi*W3(j)*(ng(idm) - bn*bg(idm))
                 elMat%Alq(ind_ff(indi),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(indi),ind_fG(ind_kf),iel) - kmult
               END DO
              ! Assembly LU
                 DO j=1,Neq
                   ind_jf = ind_asf+j
                     kmult = QdW3(idm,j)*NiNi*(ng(idm)-bn*bg(idm))
                       elMat%ALL(ind_ff(indi),ind_ff(ind_jf),iel)  = elMat%ALL(ind_ff(indi),ind_ff(ind_jf),iel) - kmult
                 END DO
            ELSEIF (k ==4) THEN
                ! Assembly LQ
                j = 1 ! all other terms are zero anyway in vector W4
                ind_kf = ind_ash + idm + (j - 1)*Ndim
                 kmult = NiNi*W4(j)*(ng(idm) - bn*bg(idm))
                 elMat%Alq(ind_ff(indi),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(indi),ind_fG(ind_kf),iel) - kmult
              ! Assembly LU
                 DO j=1,Neq
                   ind_jf = ind_asf+j
                     kmult = QdW4(idm,j)*NiNi*(ng(idm)-bn*bg(idm))
                       elMat%ALL(ind_ff(indi),ind_ff(ind_jf),iel)  = elMat%ALL(ind_ff(indi),ind_ff(ind_jf),iel) - kmult
                 END DO
#endif
              END IF
        END DO
#ifdef KEQUATION
#ifdef DKLINEARIZED
            DO j=1,Neq
             ind_jf = ind_asf + j
             kmult = ddk_dU(j)*Qpr(idm,k)*(ng(idm)-bn*bg(idm))*NiNi
             elMat%All(ind_ff(indi),ind_ff(ind_jf),iel) = elMat%All(ind_ff(indi),ind_ff(ind_jf),iel) - kmult
            enddo
            kmultf = ddk_dU_U*(Qpr(idm,k)*ng(idm)-bn*bg(idm))*Ni
            elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - kmultf

#endif
#endif
      END DO

    ELSE

      DO k = 1,Neq
#ifdef NEUTRAL
              IF (k .EQ. 5) CYCLE
#endif
        DO idm = 1,Ndim
          indi = ind_asf + k
          indj = ind_ash + idm + (k - 1)*Ndim
          ! Tangent case
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*ng(idm)
        END DO
      END DO
    ENDIF

#ifdef VORTICITY
    IF (ntang) THEN
      k = 3
      ! Vorticity equation
           IF (switch%dirivortlim) THEN
        ! Dirichlet weak form for the vorticity equation: set vorticity to 0!!!
        indi = k + ind_asf
              elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - numer%tau(k)*NiNi
        !elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(k)*kmult(indi)
           ELSE
        DO idm = 1,Ndim
          indi = ind_asf + k
          indj = ind_ash + idm + (k - 1)*Ndim
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*ng(idm)
        END DO
           END IF
    END IF

    ! Potential equation
    IF (ntang) THEN
      IF (switch%fixdPotLim) THEN
        !*****************************************
        ! Fixed potential on the limiter
        !*****************************************
        i = 4
        j = 4
        indi = ind_asf + i
        indj = ind_asf + j
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) - numer%tau(4)*NiNi
        elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(4)*Ni*phys%Mref*phys%Potfloat
      ELSE
        !*****************************************
        ! Bohm condition for the potential (Robin)
        !*****************************************
        ! Diagonal part (\Gradpar Phi)
        i = 4
        indi = ind_asf + i
        DO idm = 1,Ndim
          indj = ind_ash + idm + (i - 1)*Ndim
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + NiNi*bg(idm)
        END DO

        ! Non diagonal part --> non linear part (\Gradpar n / n)
        i = 4
        j = 1
        indi = ind_asf + i
        indj = ind_asf + j
        DO idm = 1,Ndim
          indk = ind_ash + idm + (j - 1)*Ndim
          kcoeff = phys%Mref*bg(idm)
          ! Linear part
          elMat%Alq(ind_ff(indi),ind_fG(indk),iel) = elMat%Alq(ind_ff(indi),ind_fG(indk),iel) - kcoeff*NiNi/ufg(1)
          ! Non linear correction
                 elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) + &
            &kcoeff*Qpr(idm,1)*NiNi/ufg(1)**2
          elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) + kcoeff*Qpr(idm,1)*Ni/ufg(1)
        END DO

        ! Non diagonal part --> linear part (\Gammma * \Lambda)
        i = 4
        j = 2
        indi = ind_asf + i
        indj = ind_asf + j
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) + &
          &phys%etapar/phys%c2*NiNi*phys%Mref*phys%Potfloat

        ! Diagonal part --> non linear part (\Gammma * \phi)
        ! \phi^(k+1)*\Gamma^k
        i = 4
        j = 4
        indi = ind_asf + i
        indj = ind_asf + j
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) - &
          &phys%etapar/phys%c2*NiNi*ufg(2)
        ! \phi^k*\Gamma^(k+1)
        i = 4
        j = 2
        indi = ind_asf + i
        indj = ind_asf + j
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) - &
          &phys%etapar/phys%c2*NiNi*ufg(4)
        ! \phi^k*\Gamma^k
        elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - phys%etapar/phys%c2*Ni*ufg(4)*ufg(2)

      END IF
    END IF

#endif


#ifdef NEUTRAL
#ifdef NEUTRALP
      ! Compute Vpn(U^(k-1))
      CALL computeVpn(ufg,Vpn)
        gmpn = MATMUL(Qpr,Vpn)                       ! Ndim x 1
	     ! Compute dVpn-dU(U^(k-1))
	     CALL compute_dVpn_dU(ufg,dVpn_dU)
        Taupn = MATMUL(Qpr,dVpn_dU)                 ! Ndim x Neq
	     ! Compute Dpn(U^(k-1))
      CALL computeDpn(ufg,Qpr,Vpn,Dpn)
      ! Compute dDpn_dU(U^(k-1))
      CALL compute_dDpn_dU(ufg,Qpr,Vpn,dDpn_dU)
      ! Reduce Grad Pn for low collision regime
      ! Threshold set at 0.5xGradPn for Ti = 0.2 eV
      Gammaredpn = 1.
      Tmin = 0.2/simpar%refval_temperature
        IF (Tmin/upfg(7) .LE. 1.) Gammaredpn = Gammaredpn*Tmin/upfg(7)
      Dnn = Gammaredpn*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upfg(7)*Dpn
      ! Comput Gammaredpn(U^(k-1))
      !CALL computeGammared(ufg,Gammaredpn)
      !gmipn = matmul(Qpr,Vveci)
      !CALL computeGammaLim(ue,Qpr,Vpn,GammaLim)
      ! Set Grad Ti = 0. for low collision regime
      ! (back to diffusion equation for neutral density)
      !CALL computeAlphaCoeff(ufg,Qpr,Vpn,Alphanp)
      !CALL computeBetaCoeff(ufg,Qpr,Vpn,Betanp)
      !Dpn = Alphanp*Dpn
      !dDpn_dU = Alphanp*dDpn_dU
      !Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upfg(7)*Dpn
      !IF (Dnn .gt. phys%diff_nn) Dnn = phys%diff_nn
      !IF (Dpn .gt. phys%diff_nn) THEN
      !   Dpn = Alphanp*Dpn !0.
      !   dDpn_dU = Alphanp*dDpn_dU !0.
      !   Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upfg(7)*phys%diff_nn
       !END IF
      ! Set Gamma Convective = cs_n*n_n for low collision regime
      !IF (Dpn .gt. phys%diff_nn) THEN
      !   Dpn = 0.
      !   dDpn_dU = 0.
      !   CALL jacobianMatricesBohmNP(ufg,AbohmNP)
      !ELSE
      !   AbohmNP = 0.
      !END IF
#endif

    bc = phys%bcflags(fl)

    SELECT CASE (bc)

    CASE (bc_Bohm)
       recycling_coeff =  phys%Re
       cryopump_coeff = 0.
       puff_coeff = 0.
    CASE (bc_BohmPump)
       recycling_coeff =  phys%Re_pump

       cryopump_coeff = phys%cryopump_power/(Mesh%pump_area*phys%lscale**2)/(simpar%refval_diffusion)*phys%lscale
       if (switch%testcase .ge. 50 .and. switch%testcase .le. 59) recycling_coeff = phys%Re_pump
       puff_coeff = 0.
    CASE (bc_BohmPuff)
       recycling_coeff =  phys%Re
       cryopump_coeff = 0.
       puff_coeff = phys%puff/simpar%refval_density/(Mesh%puff_area*phys%lscale**2)/(simpar%refval_diffusion)*phys%lscale
    CASE DEFAULT
      WRITE (6,*) "Error: wrong boundary type"
      STOP
    END SELECT

#ifdef SAVEFLUX
    !***************** flux control part ****************************

    !contribution from pump
    flgflux_pump = cryopump_coeff*ufg(5)	!cryopump modification
    !dimensionalizing and multiplying by the surface under this gauss point
    flgflux_pump = flgflux_pump*2.*PI*dline*simpar%refval_density*simpar%refval_speed*simpar%refval_length**2

    !contribution from puff
    flgflux_puff = puff_coeff
    !dimensionalizing and multiplying by the surface under this gauss point
    flgflux_puff = flgflux_puff*2.*PI*dline*simpar%refval_density*simpar%refval_speed*simpar%refval_length**2

    !contribution from parallel flux onto the wall
        flgflux_parallel = uefg(2)*bn
    !dimensionalizing and multiplying by the surface under this gauss point (multiplied by the local recycling)
    flgflux_parallel = recycling_coeff*flgflux_parallel*2.*PI*dline*simpar%refval_density*simpar%refval_speed*simpar%refval_length**2

    !Contribution from perpendicular plasma flux
    flgflux_perpendicular = recycling_coeff*(diffiso(1,1)*(Qpr(1,1)*ng(1) + Qpr(2,1)*ng(2))-diffani(1,1)*(Qpr(1,1)*bn*bg(1)+Qpr(2,1)*bn*bg(2)))*2.*PI*dline*simpar%refval_density*simpar%refval_speed*simpar%refval_length**2!-diffani(1,1)*(Qpr(1,1)*bn*bg(1)-Qpr(1,2)*bn*bg(2))

    !Neutral flux
    flgflux_neutral = (diffiso(5,5)*(Qpr(1,5)*ng(1) + Qpr(2,5)*ng(2)))*2.*PI*dline*simpar%refval_density*simpar%refval_speed*simpar%refval_length**2

        !flux neutral numerical
        flgflux_numerical = tau(5,5)* (uefg(5)-ufg(5))*2.*PI*dline*simpar%refval_density*simpar%refval_speed*simpar%refval_length**2


    !***************** end of flux control part *********************
#endif
#ifndef RHSBC
    ! Convective part
    k = 5
    ! Plasma flux
    indi = k+ind_asf
    indj = 2+ind_asf
    !if (ntang) then
      elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) + bn*NiNi*recycling_coeff
    !endif
    ! Pinch
    indj = 1+ind_asf
    elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) + (APinch(1,1)*ng(1) + APinch(1,2)*ng(2))*NiNi*recycling_coeff
    !Neutrals flux
!#ifdef NEUTRALP
    DO j=1,Neq
       indj = ind_asf + j
           elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) + Abohm(k,j)*NiNi*bn
       !elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) + (Abohm(k,j) + AbohmNP(j))*NiNi*bn
    END DO
!#endif
    !cryopump modification Should it be ALU?
    indj = k+ind_asf !5th equation and 5th conservative variable: pump_power*U5
    elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) - cryopump_coeff*NiNi

    ! Neutrals velocity
    !indj = Neq+ind_asf
    !elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) - (Ax(Neq,Neq)*ng(1) + Ay(Neq,Neq)*ng(2))*NiNi

    ! diffusive diagonal part
    DO idm = 1,Ndim
       k = 5
#ifndef NEUTRALP
       indi = ind_asf+k
       indj = ind_ash+idm+(k-1)*Ndim
       !if (ntang) then
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*ng(idm)*diffiso(k,k)
       !else
       !  elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*ng(idm)*diffiso(k,k)
       !endif
#ifdef DNNLINEARIZED
       DO j=1,Neq
        indj = ind_asf + j
        kmult = Dnn_dU(j)*Qpr(idm,k)*ng(idm)*NiNi
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) - kmult
           ENDDO
       kmultf = Dnn_dU_U*(Qpr(idm,k)*ng(idm))*Ni
       elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - kmultf
#endif
#else
       indi = ind_asf + k
       DO j=1,Neq
          indj = ind_asf + j
          indk = ind_ash + idm + (j-1)*Ndim
          kmult = (Dpn*Taupn(idm,j) + dDpn_dU(j)*gmpn(idm))*ng(idm)*NiNi
          !kmult = kmult - Gammaredpn*(Dpn*Taui(idm,j) + dDpn_dU(j)*gmipn(idm))*ng(idm)*NiNi
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) - kmult
          kmult = Dpn*Vpn(j)*ng(idm)*NiNi
          !kmult = kmult - Gammaredpn*Dpn*Vveci(j)*ng(idm)*NiNi
          IF (j == 5) THEN
             kmult = kmult + Dnn*ng(idm)*NiNi
          END IF
          elMat%Alq(ind_ff(indi),ind_fG(indk),iel) = elMat%Alq(ind_ff(indi),ind_fG(indk),iel) - kmult
       END DO
           kmultf = dot_PRODUCT(dDpn_dU,ufg)*gmpn(idm)*ng(idm)*Ni
       !kmultf = kmultf - Gammaredpn*(Dpn*(dot_product(Taui(1,:),ufg)*ng(1) + dot_product(Taui(2,:),ufg)*ng(2)) + dot_product(dDpn_dU,ufg)*(gmipn(1)*ng(1) + gmipn(2)*ng(2)))*Ni
       elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - kmultf
#endif
    END DO
    elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - puff_coeff*Ni

    ! diffusive non-diagonal part
    DO idm = 1,Ndim
      k = 5
      j=1
      indi = ind_asf+k
      indj = ind_ash+idm+(j-1)*Ndim
      !if (ntang) then
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*(ng(idm)*diffiso(j,j)-bn*bg(idm)*diffani(j,j))*recycling_coeff
      !else
      !  elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*ng(idm)*phys%diff_n*recycling_coeff
      !endif
    END DO
#else
  ! diffusive diagonal part
       k  = 5
       indi = ind_asf+k
       DO idm = 1,Ndim
#ifndef NEUTRALP
          indj = ind_ash+idm+(k-1)*Ndim
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*ng(idm)*diffiso(k,k)
#else
          DO j=1,Neq
              indj = ind_asf + j
              indk = ind_ash + idm + (j-1)*Ndim
              kmult = (Dpn*Taupn(idm,j) + dDpn_dU(j)*gmpn(idm))*NiNi*(ng(idm) - bn*bg(idm))
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) - kmult
              kmult = Dpn*Vpn(j)*NiNi*(ng(idm) - bn*bg(idm))
              elMat%Alq(ind_ff(indi),ind_fG(indk),iel) = elMat%Alq(ind_ff(indi),ind_fG(indk),iel) - kmult
          END DO
           kmultf = dot_PRODUCT(dDpn_dU,ufg)*gmpn(idm)*Ni*(ng(idm) - bn*bg(idm))
          elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - kmultf
#endif
      END DO
      ! Plasma outflux: now in the RHS
      elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - recycling_coeff*(ufg(2)*bn - diffiso(1,1)*(Qpr(1,1)*(ng(1) - bn*bg(1)) + Qpr(2,1)*(ng(2)- bn*bg(2))))*Ni
      ! Puff
      elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - puff_coeff*Ni
#endif
#endif


  END SUBROUTINE assembly_bohm_bc


  !*********************************
  ! Assembly Bohm
  !*********************************
  SUBROUTINE assembly_bohm_bc_new(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg,&
      &ufg,upg,bg,ng,tau,setval,delta,diffiso,diffani,ntang)
        INTEGER*4, INTENT(IN)        :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:),delta
        REAL*8, INTENT(IN)           :: NiNi(:,:),Ni(:),ufg(:),upg(:),bg(:),ng(:),tau(:,:),setval
        REAL*8, INTENT(IN)           :: diffiso(:,:),diffani(:,:), qfg(:)
        LOGICAL, INTENT(IN)          :: ntang

        REAL*8                       :: bn,Abohm(Neq,Neq)
        INTEGER                      :: i,j,k,idm
        INTEGER*4                    :: indi(Npfl),indj(Npfl)
        REAL*8                       :: Qpr(Ndim,Neq)
#ifdef TEMPERATURE
        REAL*8                       :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
        REAL*8                       :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)
#endif
#ifdef NEUTRAL
        REAL*8                       :: recycling_coeff,  puff_coeff
#endif
#ifdef VORTICITY
        REAL*8                       :: kcoeff
#endif


        bn = dot_PRODUCT(bg,ng)

    ! Compute Q^T^(k-1)
        Qpr = RESHAPE(qfg,(/Ndim,Neq/))


    !*****************************************************
    ! First equation--> Neumann homog.
    !*****************************************************
    i = 1
    indi = i + ind_asf
    ! Stabilization part
        elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
    elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
    ! Gradient part
    DO idm = 1,Ndim
      indj = ind_ash + idm + (i - 1)*Ndim
      elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
        &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) )
    END DO



    !*****************************************************
    ! Second equation--> Dirich. or Neumann
    !*****************************************************
    i = 2
    indi = i + ind_asf
        IF (numer%bohmtypebc.EQ.1) THEN
      ! Dirichlet type always if non-tangent
      IF (ntang) THEN
        ! >>>>>>>> Non tangent case (delta=1) <<<<<<<<<<
              elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - numer%tau(i)*NiNi
              IF (switch%logrho) THEN
          elMat%Alu(ind_ff(indi),ind_fe(indi - 1),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi - 1),iel) + numer%tau(i)*setval*upg(1)*NiNi
          elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(i)*setval*upg(1)*(1-ufg(1))*Ni
              ELSE
          elMat%Alu(ind_ff(indi),ind_fe(indi - 1),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi - 1),iel) + numer%tau(i)*setval*NiNi
              ENDIF
        !write(6,*) "All:"
        !call displayMatrix(elMat%All(ind_ff(indi),ind_ff(indi),iel))
        !write(6,*) "Alu:"
        !call displayMatrix(elMat%Alu(ind_ff(indi),ind_fe(indi - 1),iel))
        !write(6,*) "fh: "
        !call displayVector(elMat%fh(ind_ff(indi),iel))
        !stop
      ELSE
        ! >>>>>>>> Tangent case (delta=0) <<<<<<<<<<
        ! Stabilization part
              elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
        elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
        ! Gradient part
        DO idm = 1,Ndim
          indj = ind_ash + idm + (i - 1)*Ndim
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
            &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) )
        END DO
      END IF
        ELSEIF (numer%bohmtypebc.EQ.2) THEN
      ! Dirichlet if u<=soundspeed, Neumann if u>soundspeed if non-tangent
      IF (ntang) THEN
        ! >>>>>>>> Non tangent case  <<<<<<<<<<
              IF (delta==1) THEN
          ! Velocity is subsonic--> Dirichlet
                 elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - numer%tau(i)*NiNi
                 IF (switch%logrho) THEN
            elMat%Alu(ind_ff(indi),ind_fe(indi - 1),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi - 1),iel) + numer%tau(i)*setval*upg(1)*NiNi
            elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(i)*setval*upg(1)*(1-ufg(1))*Ni
                 ELSE
            elMat%Alu(ind_ff(indi),ind_fe(indi - 1),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi - 1),iel) + numer%tau(i)*setval*NiNi
                 ENDIF
              ELSE
          ! Velocity is supersonic--> Neumann (delta=0)
          ! Stabilization part
                 elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
          elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
          ! Gradient part
          DO idm = 1,Ndim
            indj = ind_ash + idm + (i - 1)*Ndim
            elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
              &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) )
          END DO
              ENDIF
      ELSE
        ! >>>>>>>> Tangent case (delta=0) <<<<<<<<<<
        ! Stabilization part
              elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
        elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
        ! Gradient part
        DO idm = 1,Ndim
          indj = ind_ash + idm + (i - 1)*Ndim
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
            &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) )
        END DO
      END IF
        ELSE
           WRITE(6,*) "Wrong type in Bohm bc"
           STOP
        ENDIF



#ifdef TEMPERATURE
    !*****************************************************
    ! Third and fourth equation - TEMPERATURE MODEL
    !*****************************************************
    IF (ntang) THEN

      ! Compute V(U^(k-1))
           CALL computeVi(ufg,Vveci)
           CALL computeVe(ufg,Vvece)

      ! Compute dV_dU (k-1)
           CALL compute_dV_dUi(ufg,dV_dUi)
           CALL compute_dV_dUe(ufg,dV_dUe)

      ! Compute Alpha(U^(k-1))
      Alphai = computeAlphai(ufg)
      Alphae = computeAlphae(ufg)

      ! Compute dAlpha/dU^(k-1)
           CALL compute_dAlpha_dUi(ufg,dAlpha_dUi)
           CALL compute_dAlpha_dUe(ufg,dAlpha_dUe)

      ! Jacobian matrix for convection part
      CALL jacobianMatricesBohm(ufg,Abohm)

           gmi = dot_PRODUCT(MATMUL(Qpr,Vveci),bg)  ! scalar
           gme = dot_PRODUCT(MATMUL(Qpr,Vvece),bg)             ! scalar
           Taui = MATMUL(Qpr,dV_dUi)                      ! 2x3
           Taue = MATMUL(Qpr,dV_dUe)      ! 2x3

      ! Parallel diffusion for temperature
      DO i = 1,2
        indi = ind_asf + i
        IF (i == 1) THEN
          DO j = 1,Neq
            indj = ind_asf + j
                    elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) + Abohm(i,j)*NiNi*bn
                    elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) -&
                         &coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_PRODUCT(Taui(:,j),bg)))*NiNi*bn
            DO k = 1,Ndim
              elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                &coefi*Alphai*Vveci(j)*bg(k)*NiNi*bn
            END DO
          END DO
                 elMat%fh(ind_ff(indi+2),iel) = elMat%fh(ind_ff(indi+2),iel) - coefi*Alphai*( dot_PRODUCT (MATMUL(TRANSPOSE(Taui),bg),ufg)  )*Ni*bn
        ELSE
          DO j = 1,Neq
            indj = ind_asf + j
                    elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) + Abohm(i,j)*NiNi*bn
                    elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi + 2),ind_ff(indj),iel) -&
                         &coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_PRODUCT(Taue(:,j),bg)))*NiNi*bn
            DO k = 1,Ndim
              elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                &coefe*Alphae*Vvece(j)*bg(k)*NiNi*bn
            END DO
          END DO
                 elMat%fh(ind_ff(indi+2),iel) = elMat%fh(ind_ff(indi+2),iel) - coefe*Alphae*( dot_PRODUCT (MATMUL(TRANSPOSE(Taue),bg),ufg)  )*Ni*bn
        ENDIF
      END DO


      DO i=3,4
        ! Stabilization part--> use tau for parallel diffusion
        indi = i + ind_asf
              elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
        elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
      END DO
    ELSE
      DO i=3,4
        ! Stabilization part--> use tau for perpendicular diffusion
        indi = i + ind_asf
              elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(1,1)*NiNi ! TODO: I assume the diffusion is the same in all eqs.
        elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(1,1)*NiNi ! TODO: I assume the diffusion is the same in all eqs.
      END DO

    END IF ! tangency

    ! Perpendicular diffusion for energies
    DO k = 3,4
      DO idm = 1,Ndim
        indi = ind_asf + k
        indj = ind_ash + idm + (k - 1)*Ndim
        ! Non-tangent case
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
          &NiNi*(ng(idm)*diffiso(k,k)-bn*bg(idm)*diffani(k,k) )
      END DO
    END DO
#endif

    !! Perpendicular diffusion
    !IF (ntang) THEN
    !   DO k = 1,Neqgrad
    !      DO idm = 1,Ndim
    !         indi = ind_asf + k
    !         indj = ind_ash + idm + (k - 1)*Ndim
    !         ! Non-tangent case
    !         elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
    !                                                &NiNi*(ng(idm)*diffiso(k,k)-bn*bg(idm)*diffani(k,k) )
    !      END DO
    !   END DO
    !ELSE
    !   DO k = 1,Neq
    !      DO idm = 1,Ndim
    !         indi = ind_asf + k
    !         indj = ind_ash + idm + (k - 1)*Ndim
    !         ! Tangent case
    !         elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*ng(idm)
    !      END DO
    !   END DO
    !ENDIF

#ifdef VORTICITY

    !*****************************************************
    ! Third equation - VORTICITY MODEL
    !*****************************************************
    ! Vorticity equation
    IF (ntang) THEN
      !! NON-TANGENT CASE FOR VORTICITY
      i = 3
      indi = i + ind_asf
      ! Vorticity equation
           IF (switch%dirivortlim) THEN
        ! Dirichlet weak form for the vorticity equation: set vorticity to 0!!!
              elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - numer%tau(i)*NiNi
        !elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(k)*kmult(indi)
           ELSE
        ! Neumann for the vorticity equation
        ! Stabilization part
              elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
        elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
        ! Gradient part
        DO idm = 1,Ndim
          indj = ind_ash + idm + (i - 1)*Ndim
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - &
            &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) )
        END DO
           END IF
    ELSE
      !! TANGENT CASE FOR VORTICITY
      ! Neumann for the vorticity equation
      ! Stabilization part
      i = 3
      indi = i + ind_asf
           elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
      elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
      ! Gradient part
      DO idm = 1,Ndim
        indj = ind_ash + idm + (i - 1)*Ndim
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - &
          &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) )
      END DO

    ENDIF
    !*****************************************************
    ! Fourth equation - VORTICITY MODEL
    !*****************************************************
    ! Potential equation
    IF (ntang) THEN
      !! NON-TANGENT CASE FOR POTENTIAL
      IF (switch%fixdPotLim) THEN
        !*****************************************
        ! Fixed potential on the limiter
        !*****************************************
        i = 4
        indi = ind_asf + i
              elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - numer%tau(i)*NiNi
        elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(i)*Ni*phys%Mref*phys%Potfloat
      ELSE
        !*****************************************
        ! Bohm condition for the potential (Robin)
        !*****************************************
        ! Diagonal part (\Gradpar Phi)
        i = 4
        indi = ind_asf + i
        DO idm = 1,Ndim
          indj = ind_ash + idm + (i - 1)*Ndim
          elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + NiNi*bg(idm)
        END DO

        ! Non diagonal part --> non linear part (\Gradpar n / n)
        i = 4
        j = 1
        indi = ind_asf + i
        indj = ind_asf + j
        DO idm = 1,Ndim
          indk = ind_ash + idm + (j - 1)*Ndim
          kcoeff = phys%Mref*bg(idm)
          ! Linear part
          elMat%Alq(ind_ff(indi),ind_fG(indk),iel) = elMat%Alq(ind_ff(indi),ind_fG(indk),iel) - kcoeff*NiNi/ufg(1)
          ! Non linear correction
                 elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) + &
            &kcoeff*Qpr(idm,1)*NiNi/ufg(1)**2
          elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) + kcoeff*Qpr(idm,1)*Ni/ufg(1)
        END DO

        ! Non diagonal part --> linear part (\Gammma * \Lambda)
        i = 4
        j = 2
        indi = ind_asf + i
        indj = ind_asf + j
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) + &
          &phys%etapar/phys%c2*NiNi*phys%Mref*phys%Potfloat

        ! Diagonal part --> non linear part (\Gammma * \phi)
        ! \phi^(k+1)*\Gamma^k
        i = 4
        j = 4
        indi = ind_asf + i
        indj = ind_asf + j
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) - &
          &phys%etapar/phys%c2*NiNi*ufg(2)
        ! \phi^k*\Gamma^(k+1)
        i = 4
        j = 2
        indi = ind_asf + i
        indj = ind_asf + j
              elMat%ALL(ind_ff(indi),ind_ff(indj),iel) = elMat%ALL(ind_ff(indi),ind_ff(indj),iel) - &
          &phys%etapar/phys%c2*NiNi*ufg(4)
        ! \phi^k*\Gamma^k
        elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - phys%etapar/phys%c2*Ni*ufg(4)*ufg(2)


        ! TODO: check if I need stabilization part
      END IF
    ELSE
      !! TANGENT CASE FOR POTENTIAL
      ! Neumann for the potential equation
      ! Stabilization part
      i = 4
      indi = i + ind_asf
           elMat%ALL(ind_ff(indi),ind_ff(indi),iel) = elMat%ALL(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
      elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
      ! Gradient part
      DO idm = 1,Ndim
        indj = ind_ash + idm + (i - 1)*Ndim
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - &
          &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) )
      END DO
    END IF

#endif

#ifdef NEUTRAL



    bc = phys%bcflags(fl)

    SELECT CASE (bc)

    CASE (bc_Bohm)
       recycling_coeff =  phys%Re
       puff_coeff = 0.
    CASE (bc_BohmPump)
       recycling_coeff =  phys%Re_pump
       puff_coeff = 0.
    CASE (bc_BohmPuff)
       recycling_coeff =  0.
       puff_coeff = phys%puff/simpar%refval_density/(Mesh%puff_area*phys%lscale**2)/(simpar%refval_diffusion)*phys%lscale
!       write(6,*) "puff coeff " , puff_coeff
!       stop
    CASE DEFAULT
      WRITE (6,*) "Error: wrong boundary type"
      STOP
    END SELECT

    ! convective part
    k = Neq
    indi = k+ind_asf
    indj = 2+ind_asf
        IF (ntang) THEN
      elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) + bn*NiNi*recycling_coeff
        ENDIF

    ! diffusive diagonal part
    DO idm = 1,Ndim
      k = Neq
      indi = ind_asf+k
      indj = ind_ash+idm+(k-1)*Ndim
           IF (ntang) THEN
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*ng(idm)*phys%diff_nn
           ELSE
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*ng(idm)*phys%diff_nn
           ENDIF
    END DO
    elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - puff_coeff*Ni

    ! diffusive non-diagonal part
    DO idm = 1,Ndim
      k = phys%Neq
      j=1
      indi = ind_asf+k
      indj = ind_ash+idm+(j-1)*Ndim
           IF (ntang) THEN
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*ng(idm)*phys%diff_n*recycling_coeff
           ELSE
        elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*ng(idm)*phys%diff_n*recycling_coeff
           ENDIF
    END DO
#endif

  END SUBROUTINE assembly_bohm_bc_new
END SUBROUTINE HDG_BC
