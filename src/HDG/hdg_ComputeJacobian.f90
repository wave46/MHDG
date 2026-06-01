!*****************************************
! project: MHDG
! file: hdg_ConvectionMatrices.f90
! date: 20/12/2016
! Generate the matrices that change
! during the NR iterative process
!*****************************************

SUBROUTINE HDG_computeJacobian()
  USE globals
  USE LinearAlgebra, only: tensorProduct, tensorSumInt
  USE analytical, only: body_force, analytical_solution
  USE physics
  USE transport_models_1d, ONLY: transport_model_1d
  USE hdg_limitingtechniques, ONLY:HDG_ShockCapturing
  USE diagnostics

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
  LOGICAL               :: save_tau, limiter_diagnostics, limiter_active
  INTEGER*4             :: inde(Mesh%Nnodesperelem),expected_diag_size
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
  REAL*8                :: external_heating_ions_el(refElPol%Nnodes2d),external_heating_electrons_el(refElPol%Nnodes2d)
  real*8                :: omegael(refElPol%Nnodes2d),q_cylel(refElPol%Nnodes2d),q_cylfl(refElPol%Nfacenodes),omegafl(refElPol%Nfacenodes)
  REAL*8                :: Jtorel(refElPol%Nnodes2d)
  REAL*8                :: n,El_n,nn,El_nn
  REAL*8                :: diff_nn_Vol_el(refElPol%NGauss2D),v_nn_Vol_el(refElPol%NGauss2D,Mesh%Ndim),Xg_el(refElPol%NGauss2D,Mesh%Ndim)
  REAL*8                :: diff_nn_Fac_el(refElPol%Nfaces*refElPol%NGauss1D),v_nn_Fac_el(refElPol%Nfaces*refElPol%NGauss1D,Mesh%Ndim)
#endif
#ifdef PARALL
  INTEGER               :: ierr
#endif

  IF (utils%printint .GT. 1) THEN
    IF(MPIvar%glob_id .EQ. 0) THEN
      WRITE (6,*) '*************************************************'
      WRITE (6,*) '*          COMPUTING JACOBIAN                   *'
      WRITE (6,*) '*************************************************'
    ENDIF
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
  limiter_diagnostics = TRIM(ADJUSTL(phys%neutral_flux_limiter_mode)) .NE. 'off'
  limiter_active = TRIM(ADJUSTL(phys%neutral_flux_limiter_mode)) .EQ. 'lagged_flux_limiter'
  IF (save_tau) THEN
     ALLOCATE (tau_save(refElPol%Nfaces*Mesh%Nelems*refElPol%Ngauss1d,phys%neq))
     ALLOCATE (xy_g_save(refElPol%Nfaces*Mesh%Nelems*refElPol%Ngauss1d,2))
    tau_save = 0.
    xy_g_save = 0.
  ENDIF
  expected_diag_size = Mesh%Nelems*Mesh%Nnodesperelem
  IF (limiter_diagnostics) THEN
    IF (ALLOCATED(phys%neutral_flux_limiter_Dnn_Nod)) THEN
      IF (SIZE(phys%neutral_flux_limiter_Dnn_Nod) .NE. expected_diag_size) THEN
       DEALLOCATE(phys%neutral_flux_limiter_Dnn_Nod,phys%neutral_flux_limiter_phi_Nod,phys%neutral_flux_limiter_Deff_Nod, &
          &phys%neutral_flux_limiter_Gamma_unlim_Nod,phys%neutral_flux_limiter_Gamma_max_Nod,phys%neutral_flux_limiter_activation_ratio_Nod, &
          &phys%neutral_flux_limiter_Gamma_lim_Nod)
      ENDIF
    ENDIF

    IF (.NOT. ALLOCATED(phys%neutral_flux_limiter_Dnn_Nod)) THEN
      ALLOCATE(phys%neutral_flux_limiter_Dnn_Nod(expected_diag_size))
      ALLOCATE(phys%neutral_flux_limiter_phi_Nod(expected_diag_size))
      ALLOCATE(phys%neutral_flux_limiter_Deff_Nod(expected_diag_size))
      ALLOCATE(phys%neutral_flux_limiter_Gamma_unlim_Nod(expected_diag_size))
      ALLOCATE(phys%neutral_flux_limiter_Gamma_max_Nod(expected_diag_size))
      ALLOCATE(phys%neutral_flux_limiter_activation_ratio_Nod(expected_diag_size))
      ALLOCATE(phys%neutral_flux_limiter_Gamma_lim_Nod(expected_diag_size))
    ENDIF

    phys%neutral_flux_limiter_Dnn_Nod = 0.d0
    phys%neutral_flux_limiter_phi_Nod = 0.d0
    phys%neutral_flux_limiter_Deff_Nod = 0.d0
    phys%neutral_flux_limiter_Gamma_unlim_Nod = 0.d0
    phys%neutral_flux_limiter_Gamma_max_Nod = 0.d0
    phys%neutral_flux_limiter_activation_ratio_Nod = 0.d0
    phys%neutral_flux_limiter_Gamma_lim_Nod = 0.d0
  ENDIF
#ifdef NEUTRAL
  phys%neutral_wall_source_puff_total = 0.d0
  phys%neutral_wall_source_pump_total = 0.d0
  IF (switch%neutral_wall_sources_in_elements) THEN
    IF (ALLOCATED(phys%neutral_wall_source_puff_Nod)) THEN
      IF (SIZE(phys%neutral_wall_source_puff_Nod) .NE. expected_diag_size) THEN
       DEALLOCATE(phys%neutral_wall_source_puff_Nod,phys%neutral_wall_source_pump_Nod,phys%neutral_wall_source_net_Nod)
      ENDIF
    ENDIF
    IF (.NOT. ALLOCATED(phys%neutral_wall_source_puff_Nod)) THEN
      ALLOCATE(phys%neutral_wall_source_puff_Nod(expected_diag_size))
      ALLOCATE(phys%neutral_wall_source_pump_Nod(expected_diag_size))
      ALLOCATE(phys%neutral_wall_source_net_Nod(expected_diag_size))
    ENDIF
    phys%neutral_wall_source_puff_Nod = 0.d0
    phys%neutral_wall_source_pump_Nod = 0.d0
    phys%neutral_wall_source_net_Nod = 0.d0
  ENDIF
#endif

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

#ifndef TOR3D
#ifdef NEUTRAL
  CALL diag%reset_particles_content()
  CALL account_neutral_reaction_source_totals()
  IF (switch%neutral_wall_sources_in_elements) THEN
    CALL account_neutral_wall_source_element_totals()
  ENDIF
#endif
#endif

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

        CALL assemblyVolumeContribution(Auq,Auu,rhs,b(g,:),divbg,driftg,force(g,:),&
          &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upg(g,:),&
          &ueg(g,:),qeg(g,:),u0eg(g,:,:),Jtor(g))
      END DO ! END loop in volume Gauss points
    END DO
    CALL do_assembly(Auq,Auu,rhs,ind_ass,ind_asq,iel)
    DEALLOCATE(Auq,Auu,rhs)

  ENDSUBROUTINE elemental_matrices_volume

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
          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b(g,:),Bmod(g),n_g(g,:),iel2,0.,xyf(g,:),tau)
        ELSE
          CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),0.,iel2,tau)
        ENDIF
      END IF

      ! Assembly local contributions
      CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),&
        n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
    END DO
    !      END DO

  ENDSUBROUTINE elemental_matrices_faces_pol

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
            CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b(g,:),Bmod(g),n_g(g,:),iel2,0.,xyf(g,:),tau)
          ELSE
            CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),0.,iel2,tau)
          ENDIF
        END IF

        ! Assembly local contributions
        CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),&
          n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)

      END DO ! Gauss points
    END DO
    !      END DO ! 2 elements

  ENDSUBROUTINE elemental_matrices_faces_int

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
          CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),&
            n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
        ELSE
          CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),&
            n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
        ENDIF
#else
        CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),&
          n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)

#endif
      END DO ! Gauss points
    END DO

  ENDSUBROUTINE elemental_matrices_faces_ext

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
  ENDSUBROUTINE set_permutations

#else
!TOR3D
  IF(switch%testcase .NE. 60) THEN
    ! large aspect ratio assumption and low width of source
     phys%heating_amplitude = phys%heating_power/2./PI**2/phys%heating_sigmar/phys%heating_sigmaz/(phys%r_axis+phys%heating_dr)
  ENDIF
  !********************************************
  !
  !                 2D routines
  !
  !********************************************

  !************************************
  !   Loop in elements in 2D
  !************************************

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(iel,ifa,iface,inde,indf,Xel,Xfl,i,qe,qef,ue,uef,uf,u0e,Bel,Bfl,fluxel,omegael,q_cylel,psiel,external_heating_ions_el,external_heating_electrons_el,psifl,q_cylfl,omegafl,isdir,Jtorel,El_n,El_nn) &
  !$OMP PRIVATE(Xg_el,diff_nn_Vol_el,diff_nn_Fac_el,v_nn_Vol_el,v_nn_Fac_el,xy_g_save,xy_g_save_el,tau_save,tau_save_el)&
  !$OMP FIRSTPRIVATE(phys)

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
    
    !external heating
    IF (switch%external_heating) THEN
      external_heating_ions_el = phys%external_heating_ions(Mesh%T(iel,:))
      external_heating_electrons_el = phys%external_heating_electrons(Mesh%T(iel,:))
    ELSE
      external_heating_ions_el = 0.
      external_heating_electrons_el = 0.
    ENDIF

    ! Normalized magnetic flux of the nodes of the element: PSI el
    psiel = phys%magnetic_psi(Mesh%T(iel,:))


    !omega and q_cyl on nodes of the element

     IF (switch%testcase == 60) THEN
      q_cylel = geom%q
      ! to finish this
        omegael = SQRT(Bel(:,1)**2+Bel(:,2)**2+Bel(:,3)**2)*simpar%refval_charge/simpar%refval_mass*simpar%refval_time
     ELSE
      q_cylel = phys%q_cyl(Mesh%T(iel,:))
      omegael = phys%omega(Mesh%T(iel,:))
     ENDIF


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
    CALL elemental_matrices_volume(iel,Xel,Bel,fluxel,omegael,q_cylel,psiel,external_heating_ions_el,external_heating_electrons_el,qe,ue,u0e,Jtorel,El_n,El_nn,diff_nn_Vol_el,v_nn_Vol_el,Xg_el)

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
#ifdef NEUTRAL
      nn = nn + El_nn
#endif
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

      if (switch%testcase == 60) then
        q_cylfl(:) = geom%q
        omegafl(:) = SQRT(Bfl(:,1)**2+Bfl(:,2)**2+Bfl(:,3)**2)*simpar%refval_charge/simpar%refval_mass*simpar%refval_time
      else
        q_cylfl(:) = phys%q_cyl(Mesh%T(iel,refElPol%face_nodes(ifa,:)))
        omegafl(:) = phys%omega(Mesh%T(iel,refElPol%face_nodes(ifa,:)))
      endif

      if (iface.le.Mesh%Nintfaces) then
        CALL elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,omegafl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
      else

        if (Mesh%periodic_faces(iface-Mesh%Nintfaces).eq.0) then
          CALL elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,Bfl,psifl,omegafl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
        else
          ! periodic face
          CALL elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,omegafl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
        endif
      endif


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

#ifdef NEUTRAL
  CALL diag%account_particle_content(n, nn)
#endif
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, n, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, nn, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

  IF (MPIvar%glob_id.EQ.0) THEN
     IF((switch%ME .EQV. .TRUE.) .AND. (switch%testcase .GE. 80)) THEN
        WRITE(6,*) 'D_n = ', phys%ME_diff_n*simpar%refval_length**2/simpar%refval_time
     ENDIF
     IF ((switch%ME .EQV. .TRUE.) .AND. (switch%diff_reverse_Ip .EQV. .TRUE.)) THEN
        WRITE(6,*) 'D_n = ', phys%ME_diff_n*simpar%refval_length**2/simpar%refval_time
        WRITE(6,*) 'mu = ', phys%ME_diff_u*simpar%refval_length**2/simpar%refval_time
        WRITE(6,*) 'chi_i', phys%ME_diff_e*simpar%refval_length**2/simpar%refval_time
        WRITE(6,*) 'chi_e', phys%ME_diff_ee*simpar%refval_length**2/simpar%refval_time
     ENDIF
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

#ifndef TOR3D
#ifdef NEUTRAL
  SUBROUTINE account_neutral_reaction_source_totals()

    INTEGER :: iel,g
    INTEGER :: ind_nodes(Mesh%Nnodesperelem)
    REAL*8  :: Xel(Mesh%Nnodesperelem,2),ue(Mesh%Nnodesperelem,phys%Neq)
    REAL*8  :: xy(refElPol%Ngauss2d,2),ueg(refElPol%Ngauss2d,phys%Neq)
    REAL*8  :: J11(refElPol%Ngauss2d),J12(refElPol%Ngauss2d)
    REAL*8  :: J21(refElPol%Ngauss2d),J22(refElPol%Ngauss2d)
    REAL*8  :: detJ(refElPol%Ngauss2d),dvolu,dim_factor
    REAL*8  :: niz,nrec,sigmaviz,sigmavrec,sigmavcx
    REAL*8  :: total_ionization,total_recombination,total_charge_exchange

    total_ionization = 0.d0
    total_recombination = 0.d0
    total_charge_exchange = 0.d0
    dim_factor = 2.d0*PI*simpar%refval_density*simpar%refval_speed*simpar%refval_length**2

    DO iel = 1,Mesh%Nelems
#ifdef PARALL
      IF (Mesh%ghostElems(iel) .NE. 0) CYCLE
#endif
      ind_nodes = (iel - 1)*Mesh%Nnodesperelem + (/(g,g=1,Mesh%Nnodesperelem)/)
      Xel = Mesh%X(Mesh%T(iel,:),:)
      ue = ures(ind_nodes,:)
      xy = MATMUL(refElPol%N2D,Xel)
      ueg = MATMUL(refElPol%N2D,ue)
      J11 = MATMUL(refElPol%Nxi2D,Xel(:,1))
      J12 = MATMUL(refElPol%Nxi2D,Xel(:,2))
      J21 = MATMUL(refElPol%Neta2D,Xel(:,1))
      J22 = MATMUL(refElPol%Neta2D,Xel(:,2))
      detJ = J11*J22 - J21*J12

      DO g = 1,refElPol%Ngauss2d
        dvolu = refElPol%gauss_weights2D(g)*detJ(g)
        IF (switch%axisym) dvolu = dvolu*xy(g,1)
        CALL compute_niz(ueg(g,:),niz)
        CALL compute_nrec(ueg(g,:),nrec)
        CALL compute_sigmaviz(ueg(g,:),sigmaviz)
        CALL compute_sigmavrec(ueg(g,:),sigmavrec)
        CALL compute_sigmavcx(ueg(g,:),sigmavcx)
        total_ionization = total_ionization + niz*sigmaviz*dvolu*dim_factor
        total_recombination = total_recombination + nrec*sigmavrec*dvolu*dim_factor
        total_charge_exchange = total_charge_exchange + niz*sigmavcx*dvolu*dim_factor
      ENDDO
    ENDDO

    CALL diag%account_volume_particle_reactions(total_ionization, total_recombination, total_charge_exchange)
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, total_ionization, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, total_recombination, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, total_charge_exchange, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  ENDSUBROUTINE account_neutral_reaction_source_totals

  SUBROUTINE account_neutral_wall_source_element_totals()

    INTEGER :: ifac,el,fa,fl,bc,g,inn
    INTEGER :: ind_nodes(refElPol%Nfacenodes)
    REAL*8  :: Xf(refElPol%Nfacenodes,2),uef(refElPol%Nfacenodes,phys%Neq)
    REAL*8  :: xyg(refElPol%NGauss1D,2),xyder(refElPol%NGauss1D,2)
    REAL*8  :: dline,nng,puff_coeff,pump_coeff
    REAL*8  :: total_puff,total_pump

    inn = phys%idx_rhon_eq
    total_puff = 0.d0
    total_pump = 0.d0

    IF (inn <= 0) RETURN

    DO ifac = 1,Mesh%Nextfaces
      fl = Mesh%boundaryFlag(ifac)
#ifdef PARALL
      IF (fl .EQ. 0) CYCLE
      IF (Mesh%ghostFaces(Mesh%Nintfaces + ifac) .NE. 0) CYCLE
#endif

      bc = phys%bcflags(fl)
      IF ((bc .NE. bc_BohmPuff) .AND. (bc .NE. bc_BohmPump)) CYCLE

      el = Mesh%extfaces(ifac,1)
      fa = Mesh%extfaces(ifac,2)
      ind_nodes = (el - 1)*Npel + refElPol%face_nodes(fa,:)
      Xf = Mesh%X(Mesh%T(el,refElPol%face_nodes(fa,:)),:)
      uef = ures(ind_nodes,:)
      xyg = MATMUL(refElPol%N1D,Xf)
      xyder = MATMUL(refElPol%Nxi1D,Xf)

      puff_coeff = 0.d0
      pump_coeff = 0.d0
      SELECT CASE (bc)
      CASE (bc_BohmPuff)
        puff_coeff = phys%puff/simpar%refval_density/(Mesh%puff_area*phys%lscale**2)/(simpar%refval_diffusion)*phys%lscale
      CASE (bc_BohmPump)
        pump_coeff = phys%cryopump_power/(Mesh%pump_area*phys%lscale**2)/(simpar%refval_diffusion)*phys%lscale
      END SELECT

      DO g = 1,refElPol%NGauss1D
        dline = refElPol%gauss_weights1D(g)*NORM2(xyder(g,:))
        IF (switch%axisym) dline = dline*xyg(g,1)
        nng = DOT_PRODUCT(refElPol%N1D(g,:),uef(:,inn))
        total_puff = total_puff + puff_coeff*dline*2.d0*PI*simpar%refval_density*simpar%refval_speed*simpar%refval_length**2
        total_pump = total_pump + pump_coeff*nng*dline*2.d0*PI*simpar%refval_density*simpar%refval_speed*simpar%refval_length**2
      ENDDO
    ENDDO

    CALL diag%account_wall_particle_sources(total_puff, -total_pump)
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, total_puff, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, total_pump, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    phys%neutral_wall_source_puff_total = total_puff
    phys%neutral_wall_source_pump_total = total_pump

  ENDSUBROUTINE account_neutral_wall_source_element_totals
#endif
#endif

  !***************************************************
  ! Volume computation in 2D
  !***************************************************
  SUBROUTINE elemental_matrices_volume(iel,Xel,Bel,fluxel,omegael,q_cylel,psiel,external_heating_ions_el,external_heating_electrons_el,qe,ue,u0e,Jtorel,El_n,El_nn,diff_nn_Vol_el,v_nn_Vol_el,Xg_el)

      INTEGER,INTENT(IN)            :: iel
      REAL*8,INTENT(IN)             :: Xel(:,:)
      REAL*8,INTENT(IN)             :: Bel(:,:),fluxel(:),psiel(:),Jtorel(:)
      REAL*8,INTENT(IN)             :: external_heating_ions_el(:),external_heating_electrons_el(:)
      REAL*8,INTENT(IN)             :: omegael(:),q_cylel(:)
      REAL*8,INTENT(IN)             :: qe(:,:)
      REAL*8,INTENT(IN)             :: ue(:,:),u0e(:,:,:)
      REAL*8,INTENT(OUT)            :: El_n,El_nn
      REAL*8,INTENT(OUT)            :: diff_nn_Vol_el(Ng2D),v_nn_Vol_el(Ng2D,ndim),Xg_el(Ng2D,ndim)
      INTEGER*4                     :: g,NGauss,i,inn,ind_limiter
      REAL*8                        :: dvolu
      REAL*8                        :: xy(Ng2d,ndim),ueg(Ng2d,neq),u0eg(Ng2d,neq,time%tis)
      REAL*8                        :: force(Ng2d,Neq)
      REAL*8                        :: qeg(Ng2d,neq*Ndim)
      REAL*8                        :: J11(Ng2d),J12(Ng2d)
      REAL*8                        :: J21(Ng2d),J22(Ng2d)
      REAL*8                        :: detJ(Ng2d)
      REAL*8                        :: iJ11(Ng2d),iJ12(Ng2d)
      REAL*8                        :: iJ21(Ng2d),iJ22(Ng2d)
      REAL*8                        :: fluxg(Ng2d),max_flux2D,min_flux2D, Psig(Ng2d),rho_pol_norm(Ng2d)
      INTEGER*4,DIMENSION(Npel)     :: ind_ass,ind_asq
      REAL*8                        :: ktis(time%tis + 1)
      REAL*8,DIMENSION(Npel)        :: Ni,Nxg,Nyg,NNbb,Nx_ax
      REAL*8,DIMENSION(Npel,Npel)   :: NxNi,NyNi,NNxy,NNi
      REAL*8                        :: NxyzNi(Npel,Npel,3),Nxyzg(Npel,3)
      REAL*8                        :: upg(Ng2d,phys%npv)
      REAL*8                        :: Bmod_nod(Npel),b_nod(Npel,3),b(Ng2d,3),Bmod(Ng2d),divbg,driftg(3),gradbmod(3)
#ifdef KEQUATION      
      REAL*8                        :: b_tor_nod(Npel),b_tor(Ng2d),gradbtor(3)
#endif
      REAL*8                        :: omega(Ng2d),q_cyl(Ng2d)
    real*8                        :: bg(3), Jtor(Ng2d)
    real*8                        :: diff_iso_vol(Neq,Neq,Ng2d),diff_ani_vol(Neq,Neq,Ng2d)
    real*8,allocatable            :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)
    real*8                        :: auxdiffsc(Ng2d)
    real*8                        :: Pi,sigma,x0,A,r
    real*8                        :: th_n = 1.e-14
    real*8                        :: Vnng(Ndim)
    REAL*8                        :: limiter_phi,limiter_Dnn,limiter_Gamma_max,limiter_ratio
    REAL*8                        :: limiter_Gamma_unlim(Ndim),limiter_Gamma_unlim_abs
    REAL*8                        :: external_heating_ions_gauss(Ng2d), external_heating_electrons_gauss(Ng2d)

      inn = phys%idx_rhon_eq

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
      ! Toroidal magnetic field absolute value at Gauss points
      b_tor = MATMUL(refElPol%N2D,b_tor_nod)
#endif

    ! Magnetic field norm and direction at Gauss points
      Bmod = MATMUL(refElPol%N2D,Bmod_nod)
      b = MATMUL(refElPol%N2D,b_nod)




    ! omega and q_cyl at Gauss points
      omega = MATMUL(refElPol%N2D,omegael)
      q_cyl = MATMUL(refElPol%N2D,q_cylel)


    ! Normalized magnetic flux at Gauss points: PSI
      Psig = MATMUL(refElPol%N2D,psiel)
      rho_pol_norm = SQRT(MAX(Psig,1.e-10))

    ! toroidal current at Gauss points
    IF (switch%ohmicsrc) THEN
         Jtor = MATMUL(refElPol%N2D,Jtorel)
    ELSE
      Jtor = 0.
    END IF

    ! External heating at Gauss points
    IF (switch%external_heating) THEN
      external_heating_ions_gauss = MATMUL(refElPol%N2D,external_heating_ions_el)
      external_heating_electrons_gauss = MATMUL(refElPol%N2D,external_heating_electrons_el)
    ELSE
      external_heating_ions_gauss = 0.
      external_heating_electrons_gauss = 0.
    END IF

    ! Solution at Gauss points
      ueg = MATMUL(refElPol%N2D,ue)
      qeg = MATMUL(refElPol%N2D,qe)

    ! Compute diffusion at Gauss points
#ifndef KEQUATION
    CALL setLocalDiff(xy,ueg,diff_iso_vol,diff_ani_vol)
#else
    CALL setLocalDiff(xy,ueg,diff_iso_vol,diff_ani_vol,q_cyl)
#endif

    IF (switch%import_diffusion_1D) THEN
      CALL add_1D_diff(rho_pol_norm,diff_iso_vol,diff_ani_vol)
    ENDIF

    IF (switch%transport_1d) THEN
      CALL transport_model_1d%apply_1D_diffusion(rho_pol_norm,diff_iso_vol,diff_ani_vol)
    ENDIF


    if (save_tau) then
       diff_nn_Vol_el = diff_iso_vol(inn,inn,:)
      ENDIF

    IF (limiter_diagnostics) THEN
      DO i = 1,Npel
        ind_limiter = (iel - 1)*Mesh%Nnodesperelem + i
        CALL compute_Dnn(ue(i,:), limiter_Dnn)
        CALL compute_neutral_flux_limiter(ue(i,:), qe(i,:), limiter_phi, limiter_Gamma_unlim, &
          &limiter_Gamma_max, limiter_ratio)
        limiter_Gamma_unlim_abs = SQRT(DOT_PRODUCT(limiter_Gamma_unlim, limiter_Gamma_unlim))
        phys%neutral_flux_limiter_Dnn_Nod(ind_limiter) = limiter_Dnn
        phys%neutral_flux_limiter_phi_Nod(ind_limiter) = limiter_phi
        phys%neutral_flux_limiter_Deff_Nod(ind_limiter) = limiter_phi*limiter_Dnn
        phys%neutral_flux_limiter_Gamma_unlim_Nod(ind_limiter) = limiter_Gamma_unlim_abs
        phys%neutral_flux_limiter_Gamma_max_Nod(ind_limiter) = limiter_Gamma_max
        phys%neutral_flux_limiter_activation_ratio_Nod(ind_limiter) = limiter_ratio
        phys%neutral_flux_limiter_Gamma_lim_Nod(ind_limiter) = limiter_phi*limiter_Gamma_unlim_abs
      END DO
    ENDIF

    IF (limiter_active) THEN
      DO g = 1,Ng2D
        CALL compute_neutral_flux_limiter(ueg(g,:), qeg(g,:), limiter_phi, limiter_Gamma_unlim, &
          &limiter_Gamma_max, limiter_ratio)
        diff_iso_vol(inn,inn,g) = limiter_phi*diff_iso_vol(inn,inn,g)
      END DO
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
      ! external heating
      IF (switch%external_heating) THEN
        ! ion heating
        force(g,3) = force(g,3)+external_heating_ions_gauss(g)
        ! electron heating
        force(g,4) = force(g,4)+external_heating_electrons_gauss(g)
      ENDIF

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
#ifdef NEUTRAL
      El_nn = El_nn + ueg(g,inn)*2*3.1416*dvolu*phys%lscale**3
#endif
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
      CALL assemblyVolumeContribution(Auq,Auu,rhs,b(g,:),Psig(g),divbg,driftg,force(g,:),&
        &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upg(g,:),&
        &ueg(g,:),qeg(g,:),u0eg(g,:,:),Jtor(g))
#else
      CALL assemblyVolumeContribution(Auq,Auu,rhs,b(g,:),Psig(g),divbg,driftg,b_tor(g),gradbtor,omega(g),q_cyl(g),force(g,:),&
        &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upg(g,:),&
        &ueg(g,:),qeg(g,:),u0eg(g,:,:),xy(g,:),Jtor(g))
#endif

         IF (save_tau) THEN
         v_nn_Vol_el(g,:) = Vnng
         ENDIF

    END DO ! END loop in volume Gauss points
#ifdef NEUTRAL
    IF (switch%neutral_wall_sources_in_elements) THEN
      CALL add_neutral_wall_sources_to_element(iel,Xel,ue,Auu,rhs)
    ENDIF
    IF (switch%neutral_recycling_in_elements) THEN
      CALL add_neutral_recycling_to_element(iel,Xel,ue,qe,Auu,Auq,rhs)
    ENDIF
#endif
      CALL do_assembly(Auq,Auu,rhs,ind_ass,ind_asq,iel)
      DEALLOCATE(Auq,Auu,rhs)

  ENDSUBROUTINE elemental_matrices_volume

#ifdef NEUTRAL
  SUBROUTINE add_neutral_wall_sources_to_element(iel,Xel,ue,Auu,rhs)

    INTEGER,INTENT(IN)       :: iel
    REAL*8,INTENT(IN)        :: Xel(:,:),ue(:,:)
    REAL*8,INTENT(INOUT)     :: Auu(:,:,:),rhs(:,:)
    INTEGER                  :: ifa,iface,ibf,fl,bc
    INTEGER                  :: g,a,inn,z,ind_source
    REAL*8                   :: Xfl(refElPol%Nfacenodes,2)
    REAL*8                   :: xyg(refElPol%Ngauss1d,2)
    REAL*8                   :: xyder(refElPol%Ngauss1d,2)
    REAL*8                   :: xyv(refElPol%Ngauss2d,2)
    REAL*8                   :: J11(refElPol%Ngauss2d),J12(refElPol%Ngauss2d)
    REAL*8                   :: J21(refElPol%Ngauss2d),J22(refElPol%Ngauss2d)
    REAL*8                   :: detJ(refElPol%Ngauss2d)
    REAL*8                   :: dline,dvolu,element_volume
    REAL*8                   :: puff_coeff,pump_coeff,puff_flux,pump_flux
    REAL*8                   :: puff_integral,pump_integral,puff_source,pump_source
    REAL*8,DIMENSION(Npel)   :: Ni
    REAL*8,DIMENSION(Npel,Npel) :: NNi

    inn = phys%idx_rhon_eq
    IF (inn <= 0) RETURN

#ifdef PARALL
    IF (Mesh%ghostElems(iel) .NE. 0) RETURN
#endif

    z = inn + (inn - 1)*Neq
    puff_integral = 0.d0
    pump_integral = 0.d0

    DO ifa = 1,refElPol%Nfaces
      iface = Mesh%F(iel,ifa)
      IF (iface <= Mesh%Nintfaces) CYCLE

      ibf = iface - Mesh%Nintfaces
      IF (Mesh%periodic_faces(ibf) .NE. 0) CYCLE

      fl = Mesh%boundaryFlag(ibf)
#ifdef PARALL
      IF (fl .EQ. 0) CYCLE
#endif
      bc = phys%bcflags(fl)

      IF ((bc .NE. bc_BohmPuff) .AND. (bc .NE. bc_BohmPump)) CYCLE

      puff_coeff = 0.d0
      pump_coeff = 0.d0
      SELECT CASE (bc)
      CASE (bc_BohmPuff)
        puff_coeff = phys%puff/simpar%refval_density/(Mesh%puff_area*phys%lscale**2)/(simpar%refval_diffusion)*phys%lscale
      CASE (bc_BohmPump)
        pump_coeff = phys%cryopump_power/(Mesh%pump_area*phys%lscale**2)/(simpar%refval_diffusion)*phys%lscale
      END SELECT

      Xfl = Xel(refElPol%face_nodes(ifa,:),:)
      xyg = MATMUL(refElPol%N1D,Xfl)
      xyder = MATMUL(refElPol%Nxi1D,Xfl)

      DO g = 1,refElPol%Ngauss1d
        dline = refElPol%gauss_weights1D(g)*NORM2(xyder(g,:))
        IF (switch%axisym) dline = dline*xyg(g,1)

        puff_integral = puff_integral + puff_coeff*dline
        pump_integral = pump_integral + pump_coeff*dline
      ENDDO
    ENDDO

    IF ((puff_integral .NE. 0.d0) .OR. (pump_integral .NE. 0.d0)) THEN
      xyv = MATMUL(refElPol%N2D,Xel)
      J11 = MATMUL(refElPol%Nxi2D,Xel(:,1))
      J12 = MATMUL(refElPol%Nxi2D,Xel(:,2))
      J21 = MATMUL(refElPol%Neta2D,Xel(:,1))
      J22 = MATMUL(refElPol%Neta2D,Xel(:,2))
      detJ = J11*J22 - J21*J12

      element_volume = 0.d0
      DO g = 1,refElPol%Ngauss2d
        dvolu = refElPol%gauss_weights2D(g)*detJ(g)
        IF (switch%axisym) dvolu = dvolu*xyv(g,1)
        element_volume = element_volume + dvolu
      ENDDO

      IF (element_volume <= 0.d0) THEN
        WRITE(6,*) 'Negative or zero element volume while spreading neutral wall sources in element ',iel
        STOP
      ENDIF

      puff_source = puff_integral/element_volume
      pump_source = pump_integral/element_volume

      DO g = 1,refElPol%Ngauss2d
        dvolu = refElPol%gauss_weights2D(g)*detJ(g)
        IF (switch%axisym) dvolu = dvolu*xyv(g,1)
        Ni = refElPol%N2D(g,:)*dvolu
        NNi = tensorProduct(Ni,refElPol%N2D(g,:))
        rhs(:,inn) = rhs(:,inn) + puff_source*Ni
        Auu(:,:,z) = Auu(:,:,z) + pump_source*NNi
      ENDDO

      IF (ALLOCATED(phys%neutral_wall_source_puff_Nod)) THEN
        DO a = 1,Npel
          ind_source = (iel - 1)*Mesh%Nnodesperelem + a
          puff_flux = puff_source*simpar%refval_density*simpar%refval_speed/simpar%refval_length
          pump_flux = pump_source*ue(a,inn)*simpar%refval_density*simpar%refval_speed/simpar%refval_length
          phys%neutral_wall_source_puff_Nod(ind_source) = puff_flux
          phys%neutral_wall_source_pump_Nod(ind_source) = pump_flux
          phys%neutral_wall_source_net_Nod(ind_source) = puff_flux - pump_flux
        ENDDO
      ENDIF
    ENDIF

  ENDSUBROUTINE add_neutral_wall_sources_to_element

  SUBROUTINE add_neutral_recycling_to_element(iel,Xel,ue,qe,Auu,Auq,rhs)

    INTEGER,INTENT(IN)       :: iel
    REAL*8,INTENT(IN)        :: Xel(:,:),ue(:,:),qe(:,:)
    REAL*8,INTENT(INOUT)     :: Auu(:,:,:),Auq(:,:,:),rhs(:,:)

    WRITE(6,*) 'neutral_recycling_in_elements is not implemented yet. Element = ',iel
    STOP

  ENDSUBROUTINE add_neutral_recycling_to_element
#endif

  !***************************************************
  ! Interior faces computation in 2D
  !***************************************************

  SUBROUTINE elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,omegafl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)

    integer,intent(IN)        :: iel,ifa
    real*8,intent(IN)         :: Xfl(:,:)
    real*8,intent(IN)         :: Bfl(:,:), psifl(:)
    real*8,intent(IN)         :: qef(:,:)
    real*8,intent(IN)         :: uef(:,:),uf(:,:)
    real*8,intent(IN)             :: q_cylfl(:)
    real*8,intent(in)         :: omegafl(:)
    real*8,intent(out)        :: diff_nn_Fac_el(:),v_nn_Fac_el(:,:),tau_save_el(:,:),xy_g_save_el(:,:)
    integer*4                 :: g,NGauss,i,indsave(Ng1d),inn
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
    real*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ng1d,3),Bmod(Ng1d),Psig(Ng1d),rho_pol_norm(Ng1d)
    real*8                    :: diff_iso_fac(Neq,Neq,Ng1d),diff_ani_fac(Neq,Neq,Ng1d)
    real*8                    :: auxdiffsc(Ng1d)
    real*8                    :: q_cyl(Ng1d)
    real*8                    :: omega(Ng1d)
    REAL*8                    :: limiter_phi,limiter_Gamma_max,limiter_ratio
    REAL*8                    :: limiter_Gamma_unlim(Ndim)
    inn = phys%idx_rhon_eq

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


    ! q_cyl and omega at Gauss points
    q_cyl = matmul(refElPol%N1D,q_cylfl)
    omega = matmul(refElPol%N1D,omegafl)

    ! Normalaized magnetic flux at Gauss points: PSI
      Psig = MATMUL(refElPol%N1d,psifl)
      rho_pol_norm = SQRT(MAX(Psig,1.e-10))

    ! Element solution at face Gauss points
      uefg = MATMUL(refElPol%N1D,uef)
    ! Gradient solution at face gauss points
      qfg = MATMUL(refElPol%N1D,qef)

    ! Compute diffusion at faces Gauss points
#ifndef KEQUATION
    CALL setLocalDiff(xyf,uefg,diff_iso_fac,diff_ani_fac)
#else
    CALL setLocalDiff(xyf,uefg,diff_iso_fac,diff_ani_fac,q_cyl)
#endif

    IF (switch%import_diffusion_1D) THEN
      CALL add_1D_diff(rho_pol_norm,diff_iso_fac,diff_ani_fac)
    ENDIF

    IF (switch%transport_1d) THEN
      CALL transport_model_1d%apply_1D_diffusion(rho_pol_norm,diff_iso_fac,diff_ani_fac)
    ENDIF

    IF (limiter_active) THEN
      DO g = 1,Ng1d
        CALL compute_neutral_flux_limiter(uefg(g,:), qfg(g,:), limiter_phi, limiter_Gamma_unlim, &
          &limiter_Gamma_max, limiter_ratio)
        diff_iso_fac(inn,inn,g) = limiter_phi*diff_iso_fac(inn,inn,g)
      END DO
    ENDIF

    if (save_tau) then
       indsave = (ifa - 1)*Ngauss + (/(i,i=1,Ngauss)/)
       diff_nn_Fac_el(indsave) = diff_iso_fac(inn,inn,:)
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
            CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),qfg(g,:),b(g,:),n_g,iel,0.,xyf(g,:),tau,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g))
        ELSE
          CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g,xyf(g,:),0.,iel,tau)
        ENDIF
      END IF

! Assembly local contributions
#ifdef DKLINEARIZED
      CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Psig(g),q_cyl(g),xyf(g,:),&
      n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
#else
      CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Psig(g),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
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


  ENDSUBROUTINE elemental_matrices_faces_int

  !***************************************************
  ! Exterior faces computation in 2D
  !***************************************************

  SUBROUTINE elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,Bfl,psifl,omegafl,q_cylfl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)

    integer,intent(IN)        :: iel,ifa
    real*8,intent(IN)         :: Xfl(:,:)
    real*8,intent(IN)         :: Bfl(:,:), psifl(:)
    logical,intent(IN)        :: isdir
    real*8,intent(IN)         :: qef(:,:)
    real*8,intent(INOUT)      :: uef(:,:),uf(:,:)
    real*8,intent(IN)             :: q_cylfl(:)
    real*8,intent(in)         :: omegafl(:)
    real*8,intent(out)        :: diff_nn_Fac_el(:),v_nn_Fac_el(:,:),tau_save_el(:,:),xy_g_save_el(:,:)
    integer*4                 :: g,NGauss,i,indsave(Ng1d),inn
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
    real*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ng1d,3),Bmod(Ng1d), Psig(Ng1d), rho_pol_norm(Ng1d)
    real*8                    :: diff_iso_fac(Neq,Neq,Ng1d),diff_ani_fac(Neq,Neq,Ng1d)
    real*8                    :: auxdiffsc(Ng1d)
    real*8                    :: Vnng(Ndim)
    real*8                    :: q_cyl(Ng1d)
    real*8                    :: omega(Ng1d)
    REAL*8                    :: limiter_phi,limiter_Gamma_max,limiter_ratio
    REAL*8                    :: limiter_Gamma_unlim(Ndim)
    inn = phys%idx_rhon_eq

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


    ! q_cyl and omega at Gauss points
    q_cyl = matmul(refElPol%N1D,q_cylfl)
    omega = matmul(refElPol%N1D,omegafl)
    ! Normalaized magnetic flux at Gauss points: PSI
    Psig = MATMUL(refElPol%N1D,psifl)
    rho_pol_norm = SQRT(MAX(Psig,1.e-10))

    ! Trace solution at face Gauss points
    xyf = MATMUL(refElPol%N1D,Xfl)
    IF (isdir) THEN
      CALL analytical_solution(iel,xyf(:,1),xyf(:,2),ufg)
    ELSE
#ifdef PARALL
      ! FOR SOME UNKNOWN REASON EXTERNAL FACES DO NOT NEED TO BE FLIPPED ONLY AT THE VERY FIRST ITERATION.
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
    CALL setLocalDiff(xyf,uefg,diff_iso_fac,diff_ani_fac)
#else
    CALL setLocalDiff(xyf,uefg,diff_iso_fac,diff_ani_fac,q_cyl)
#endif


    IF (switch%import_diffusion_1D) THEN
      CALL add_1D_diff(rho_pol_norm,diff_iso_fac,diff_ani_fac)
    ENDIF

    IF (switch%transport_1d) THEN
      CALL transport_model_1d%apply_1D_diffusion(rho_pol_norm,diff_iso_fac,diff_ani_fac)
    ENDIF

    IF (limiter_active) THEN
      DO g = 1,Ng1d
        CALL compute_neutral_flux_limiter(uefg(g,:), qfg(g,:), limiter_phi, limiter_Gamma_unlim, &
          &limiter_Gamma_max, limiter_ratio)
        diff_iso_fac(inn,inn,g) = limiter_phi*diff_iso_fac(inn,inn,g)
      END DO
    ENDIF

    if (save_tau) then
       indsave = (ifa -1)*Ngauss + (/(i,i=1,Ngauss)/)
       diff_nn_Fac_el(indsave) = diff_iso_fac(inn,inn,:)
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
            CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),qfg(g,:),b(g,:),n_g,iel,isext,xyf(g,:),tau,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g))
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
        CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Psig(g),&
          n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)

      ELSE
        CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Psig(g),&
          n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
      ENDIF
#else

        CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),,Psig(g),q_cyl(g),xyf(g,:),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
      ELSE
        CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Psig(g),q_cyl(g),xyf(g,:),&
          n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
      ENDIF
#endif

#else
#ifndef DKLINEARIZED
      CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Psig(g),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
#else
      CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Psig(g),q_cyl(g),xyf(g,:),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),qfg(g,:),tau)
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

  ENDSUBROUTINE elemental_matrices_faces_ext

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
  ENDSUBROUTINE set_permutations

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
  ENDSUBROUTINE setTimeIntegrationCoefficients

  !********************************************************************
  !
  !         ASSEMBLY VOLUME CONTRIBUTION
  !
  !********************************************************************
#ifndef KEQUATION
  SUBROUTINE assemblyVolumeContribution(Auq,Auu,rhs,b3,psi,divb,drift,f,&
      &ktis,diffiso,diffani,Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upe,ue,qe,u0e,Jtor)
#else
  SUBROUTINE assemblyVolumeContribution(Auq,Auu,rhs,b3,psi,divb,drift,btor,gradBtor,omega,q_cyl,f,&
    &ktis,diffiso,diffani,Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upe,ue,qe,u0e,xy,Jtor)
#endif
        REAL*8,INTENT(inout)      :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)
        REAL*8,INTENT(IN)         :: b3(:),psi,divb,drift(:),f(:),ktis(:)
#ifdef KEQUATION
    real*8,intent(IN)         :: btor,gradBtor(:), omega, q_cyl,xy(:)
#ifdef DKLINEARIZED
    real*8                    :: ddk_dU(Neq), ddk_dU_U
    real*8                    :: gradddk(Ndim)
#endif
#endif
    real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
    real*8,intent(IN)         :: Ni(:),NNi(:,:),Nxyzg(:,:),NNxy(:,:),NxyzNi(:,:,:),NNbb(:)
    real*8,intent(IN)         :: upe(:),ue(:),Jtor
    real*8,intent(INOUT)      :: u0e(:,:)
    real*8,intent(IN)         :: qe(:)
#ifdef VORTICITY
    real*8                     :: kcoeff,exb(3)
    integer*4                  :: alpha,beta,ii
#endif
    integer*4                 :: i,j,k,iord,z,inn,ign,ik
    real*8,dimension(neq,neq) :: A
    real*8,dimension(neq,Ndim):: APinch
    real*8                    :: Qpr(Ndim,Neq),bb(3)
    real*8                    :: W2(Neq),dW2_dU(Neq,Neq),QdW2(Ndim,Neq)
    real*8                    :: qq(3,Neq),b(Ndim)
    real*8                    :: grad_n(3),gradpar_n

#ifdef TEMPERATURE
    real*8,dimension(neq,neq) :: GG
#ifdef NEUTRALGAMMA
    real*8,dimension(neq,neq) :: GGn
    real*8                    :: Etan
    real*8                    :: dEtan_dU(Neq),gmGamman(Ndim)
    real*8                    :: Vun(Neq),dVun_dU(Neq,Neq),TauGamman(Ndim,Neq)
#endif
    real*8                    :: Telect
    real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
    real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
    real*8                    :: q_fs_i,q_fs_e,flux_limiter_i,flux_limiter_e, q_sh_i, q_sh_e 
    real*8                    :: flux_limiter_i_ratio, flux_limiter_e_ratio, fl_deriv_i, fl_deriv_e
    real*8                    :: dq_fs_i_dU(Neq), dq_fs_e_dU(Neq)
    real*8                    :: W,dW_dU(Neq,Neq),s,ds_dU(Neq),Zet(ndim,Neq)
    real*8                    :: Sohmic,dSohmic_dU(Neq) ! Ohmic heating
    real*8                    :: W3(Neq),dW3_dU(Neq,Neq),QdW3(Ndim,Neq)
    real*8                    :: W4(Neq),dW4_dU(Neq,Neq),QdW4(Ndim,Neq)
#ifdef NEUTRALPNEW
    real*8                    :: W5p(Neq),dW5p_dU(Neq,Neq),QdW5p(Ndim,Neq),dW5p_dU_u(Neq)
#endif
#else
    real*8                    :: auxvec(Neq)
#endif
#ifdef NEUTRAL
#ifdef KEQUATION
        REAL*8                    :: gamma_I,ce, dissip,r
        REAL*8                    :: ddissip_du(Neq)
#endif
    real*8                    :: niz,nrec,fGammacx,fGammarec
    real*8                    :: dniz_dU(Neq),dnrec_dU(Neq),dfGammacx_dU(Neq),dfGammarec_dU(Neq)
#ifdef NEUTRALGAMMA
    real*8                    :: fGammaN
    real*8                    :: dfGammaN_dU(Neq)
#endif
#ifdef TEMPERATURE
        REAL*8                    :: sigmaviz,sigmavrec,sigmavcx,fEiiz,fEirec,fEicx
#ifdef NEUTRALGAMMA
    real*8                    :: fEiN
    real*8                    :: dfEiN_dU(Neq)
#endif
    !amjuel radiation losses
    real*8                    :: sigmavEiz,sigmavErec
    real*8                    :: dsigmavEiz_dU(Neq),dsigmavErec_dU(Neq)
    real*8                    :: cooling_factor
    real*8                    :: dcooling_factor_dU(Neq)
        REAL*8                    :: dsigmaviz_dU(Neq),dsigmavrec_dU(Neq),dsigmavcx_dU(Neq)
        REAL*8                    :: dfEiiz_dU(Neq),dfEirec_dU(Neq),dfEicx_dU(Neq)
    real*8                    :: Dnn_dU(Neq), Dnn_dU_U
    REAL*8                    :: neutral_limiter_phi,neutral_limiter_Gamma_max,neutral_limiter_ratio
    REAL*8                    :: neutral_limiter_Gamma_unlim(Ndim)
#endif
    real*8                    :: Sn(Neq,Neq),Sn0(Neq)
#endif




        REAL*8 :: kmult(SIZE(Auq,1),SIZE(Auq,2))



    inn = phys%idx_rhon_eq
    ign = phys%idx_gamman_eq
    ik = phys%idx_k_eq

    b = b3(1:Ndim)

    bb = 0.
    bb = b3

    ! Jacobian for convection term
    CALL jacobianMatrices(ue,A)

    ! Jacobian for pinch term
    APinch = 0.d0
    IF (switch%transport_1d) THEN
      CALL transport_model_1d%compute_1D_pinch_matrix(b,SQRT(MAX(psi,0.d0)),APinch)
    ENDIF

    ! Compute Q^T^(k-1)
        Qpr = RESHAPE(qe,(/Ndim,Neq/))

#ifdef TEMPERATURE
    neutral_limiter_phi = 1.d0
    IF (limiter_active) THEN
      CALL compute_neutral_flux_limiter(ue, qe, neutral_limiter_phi, neutral_limiter_Gamma_unlim, &
        &neutral_limiter_Gamma_max, neutral_limiter_ratio)
    ENDIF
#endif

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
#ifdef NEUTRALGAMMA
    CALL GimpMatrixN(ue,divb,GGn)
    CALL computeEtan(ue,Etan)
    CALL compute_dEtan_dU(ue,dEtan_dU)
    CALL computeVun(ue,Vun)
    CALL compute_dVun_dU(ue,dVun_dU)
    gmGamman = MATMUL(Qpr,Vun)
    TauGamman = MATMUL(Qpr,dVun_dU)
#endif

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

    ! Compute flux limiters

    IF (switch%flux_limiter) THEN
    
      call compute_free_streaming_heat_flux_electrons(ue,q_fs_e)
      call compute_free_streaming_heat_flux_ions(ue,q_fs_i)
      ! Compute spitzer-harm heat fluxes
      q_sh_e = coefe*Alphae*gme
      q_sh_i = coefi*Alphai*gmi

      call compute_flux_limiter(q_fs_e,q_sh_e,phys%c_fle,flux_limiter_e)
      call compute_flux_limiter(q_fs_i,q_sh_i,phys%c_fli,flux_limiter_i)

  
    
      ! Derivatives of flux limiters
      flux_limiter_e_ratio = ABS(q_sh_e)/(q_fs_e)/phys%c_fle
      flux_limiter_i_ratio = ABS(q_sh_i)/(q_fs_i)/phys%c_fli

      fl_deriv_e = q_sh_e*flux_limiter_e_ratio/q_fs_e
      fl_deriv_i = q_sh_i*flux_limiter_i_ratio/q_fs_i

      call compute_dfree_streaming_heat_flux_electrons_dU(ue,dq_fs_e_dU)
      call compute_dfree_streaming_heat_flux_ions_dU(ue,dq_fs_i_dU)
    ELSE 
      flux_limiter_e = 1.
      flux_limiter_i = 1.
      fl_deriv_e = 0.
      fl_deriv_i = 0.
      dq_fs_e_dU = 0.
      dq_fs_i_dU = 0.
    END IF

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
#ifdef NEUTRALPNEW
    CALL compute_W5p(ue,W5p)
    CALL compute_dW5p_dU(ue,dW5p_dU)
        QdW5p = MATMUL(Qpr,dW5p_dU)
        dW5p_dU_u = MATMUL(dW5p_dU,ue)
    IF (limiter_active) THEN
      W5p = neutral_limiter_phi*W5p
      dW5p_dU = neutral_limiter_phi*dW5p_dU
      QdW5p = neutral_limiter_phi*QdW5p
      dW5p_dU_u = neutral_limiter_phi*dW5p_dU_u
    ENDIF
#endif

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
    IF ((ue(ik)<1.e-20) .or. (ue(1)<1.e-20) .or.(ue(3)<1.e-20) .or. (ue(4)<1.e-20)) THEN
      dissip =  abs(gamma_I)*dissip/phys%k_max
      ddissip_du = abs(gamma_I)
      gamma_I = 0.
        ELSEIF (ue(ik)>phys%k_max) THEN
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
#ifdef NEUTRALGAMMA
        CALL compute_fGammaN(ue,fGammaN)
        CALL compute_dfGammaN_dU(ue,dfGammaN_dU)
#endif
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
#ifdef NEUTRALGAMMA
        CALL compute_fEiN(ue,fEiN)
        CALL compute_dfEiN_dU(ue,dfEiN_dU)
#endif
    !Neutral Source Terms needed in the electron energy equation
    IF (switch%impurity_radiation) THEN
        CALL compute_cooling_factor(ue,cooling_factor)
        CALL compute_dcooling_factor_dU(ue,dcooling_factor_dU)
    ELSE
        cooling_factor = 0.
        dcooling_factor_dU = 0.
    ENDIF

    call compute_sigmavEiz(ue,sigmavEiz)
    call compute_sigmavErec(ue,sigmavErec)
    call compute_dsigmavEiz_dU(ue,dsigmavEiz_dU)
    call compute_dsigmavErec_dU(ue,dsigmavErec_dU)

        CALL compute_Dnn_dU(ue,Dnn_dU)
        Dnn_dU_u = dot_PRODUCT(Dnn_dU,Ue)
        IF (limiter_active) THEN
          Dnn_dU = neutral_limiter_phi*Dnn_dU
          Dnn_dU_u = neutral_limiter_phi*Dnn_dU_u
        ENDIF

#endif

    !Assembly the matrix for neutral sources
#ifdef TEMPERATURE
IF (switch%impurity_radiation) THEN
#ifdef NEUTRALGAMMA
  call assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
    &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,fGammaN,dfGammaN_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
    &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,fEiN,dfEiN_dU,Sn,Sn0, &
    sigmavEiz=sigmavEiz,dsigmavEiz_dU=dsigmavEiz_dU,sigmavErec=sigmavErec,dsigmavErec_dU=dsigmavErec_dU,&
    cooling_factor=cooling_factor,dcooling_factor_dU=dcooling_factor_dU)
#else
  call assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
    &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
    &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,Sn,Sn0, &
    sigmavEiz=sigmavEiz,dsigmavEiz_dU=dsigmavEiz_dU,sigmavErec=sigmavErec,dsigmavErec_dU=dsigmavErec_dU,&
    cooling_factor=cooling_factor,dcooling_factor_dU=dcooling_factor_dU)
#endif
ELSE
#ifdef NEUTRALGAMMA
  call assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
    &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,fGammaN,dfGammaN_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
    &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,fEiN,dfEiN_dU,Sn,Sn0, &
    sigmavEiz=sigmavEiz,dsigmavEiz_dU=dsigmavEiz_dU,sigmavErec=sigmavErec,dsigmavErec_dU=dsigmavErec_dU)
#else
  call assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
    &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
    &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,Sn,Sn0, &
    sigmavEiz=sigmavEiz,dsigmavEiz_dU=dsigmavEiz_dU,sigmavErec=sigmavErec,dsigmavErec_dU=dsigmavErec_dU)
#endif
ENDIF
#else
#ifdef NEUTRALGAMMA
        CALL assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,fGammaN,dfGammaN_dU,Sn,Sn0)
#else
        CALL assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,Sn,Sn0)
#endif
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
                    Auu(:,:,i+(j-1)*Neq) = Auu(:,:,i+(j-1)*Neq)+flux_limiter_i**2*(coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_PRODUCT(Taui(:,j),b)))+fl_deriv_i*dq_fs_i_dU(j))*NNxy + &
                         &(dot_PRODUCT(Zet(:,j),b) + ds_dU(j))*NNi
            DO k = 1,Ndim
              z = i+(k-1)*Neq+(j-1)*Neq*Ndim
              Auq(:,:,z) = Auq(:,:,z)+ flux_limiter_i**2*coefi*Alphai*Vveci(j)*b(k)*NNxy
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
                 rhs(:,i) = rhs(:,i) + flux_limiter_i**2*coefi*Alphai*(dot_PRODUCT(MATMUL(TRANSPOSE(Taui),b),ue))*NNbb + s*Ni
        ELSEIF (i == 4) THEN
          DO j = 1,4
            z = i+(j-1)*Neq
                    Auu(:,:,z)=Auu(:,:,z)+flux_limiter_e**2*(coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_PRODUCT(Taue(:,j),b))) + fl_deriv_e*dq_fs_e_dU(j))*NNxy - &
                    & (dot_PRODUCT(Zet(:,j),b) + ds_dU(j))*NNi 
            IF (switch%ohmicsrc) THEN
              Auu(:,:,z) = Auu(:,:,z) - dSohmic_dU(j)*(Jtor**2)*NNi
            ENDIF
            DO k = 1,Ndim
              z = i+(k-1)*Neq+(j-1)*Neq*Ndim
              Auq(:,:,z)=Auq(:,:,z)+flux_limiter_e**2*coefe*Alphae*Vvece(j)*b(k)*NNxy
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
                 rhs(:,i) = rhs(:,i)+flux_limiter_e**2*coefe*Alphae*(dot_PRODUCT(MATMUL(TRANSPOSE(Taue),b),ue))*NNbb - s*Ni
          IF (switch%ohmicsrc) THEN
            rhs(:,i) = rhs(:,i) + Sohmic*(Jtor**2)*Ni
          ENDIF
          ELSEIF (i == inn) THEN
                 DO j = 1,Neq
              z = i+(j-1)*Neq
                    DO k = 1,Ndim
                Auu(:,:,z) =Auu(:,:,z) + (NxyzNi(:,:,k)*Dnn_dU(j)*Qpr(k,i))
#ifdef NEUTRALPNEW
                Auu(:,:,z) =Auu(:,:,z) + (NxyzNi(:,:,k)*QdW5p(k,j))
#endif
                    ENDDO
                 ENDDO

            DO k = 1, Ndim
              rhs(:,i) = rhs(:,i)+Dnn_dU_U*Qpr(k,i)*Nxyzg(:,k)
#ifdef NEUTRALPNEW
              rhs(:,i) = rhs(:,i)+dot_PRODUCT(Qpr(k,:),dW5p_dU_u)*Nxyzg(:,k)
#endif
                 ENDDO
#ifdef KEQUATION
        ELSEIF (i==ik) THEN
          DO j=1,Neq
            z = i+(j-1)*Neq
                    IF (j==ik) THEN
              Auu(:,:,z) = Auu(:,:, z) - (gamma_I-ddissip_du(j))*NNi
            ENDIF
          END DO
          rhs(:,i) = rhs(:,i) + dissip*Ni
#endif
#ifdef NEUTRALGAMMA
        ELSEIF (i == ign) THEN
          DO j = 1,Neq
            z = i+(j-1)*Neq
            Auu(:,:,z) = Auu(:,:,z) - GGn(i,j)*NNi
            DO k = 1,Ndim
              z = i+(k-1)*Neq+(j-1)*Neq*Ndim
              Auu(:,:,i+(j-1)*Neq) = Auu(:,:,i+(j-1)*Neq) + Etan*TauGamman(k,j)*NxyzNi(:,:,k)
              Auq(:,:,z) = Auq(:,:,z) + Etan*Vun(j)*NxyzNi(:,:,k)
            END DO
          END DO
          DO k = 1,Ndim
            rhs(:,i) = rhs(:,i) + Etan*dot_PRODUCT(TauGamman(k,:),ue)*Nxyzg(:,k)
          END DO
#endif
		END IF
#endif
#ifdef KEQUATION
#ifdef DKLINEARIZED
    ! Contribution from linearized dk term assuming so far that Dk is the same in all plasma equations
        if (i .ne. inn) then
          DO j = 1,Neq
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
#ifdef NEUTRALPNEW
          ELSEIF(i==inn) THEN
            DO j = 1,Neq
                z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                Auq(:,:,z) = Auq(:,:,z) + W5p(j)*NxyzNi(:,:,k)
            END DO
#endif
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
    ENDSUBROUTINE assemblyVolumeContribution

    !********************************************************************
    !
    !         ASSEMBLY INTERIOR FACES CONTRIBUTION
    !
    !********************************************************************

#ifdef DKLINEARIZED
  SUBROUTINE assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,&
    &ind_fg,b3,psi,q_cyl,xyf,n,diffiso,diffani,NNif,Nif,Nfbn,uf,qf,tau)
#else
    SUBROUTINE assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,&
        &ind_fg,b3,psi,n,diffiso,diffani,NNif,Nif,Nfbn,uf,qf,tau)
#endif
      integer*4,intent(IN)      :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
      real*8,intent(IN)         :: b3(:),n(:), psi
      real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
      real*8,intent(IN)         :: NNif(:,:),Nif(:),Nfbn(:)
      real*8,intent(IN)         :: uf(:)
      real*8,intent(IN)         :: qf(:)
#ifdef KEQUATION
#ifdef DKLINEARIZED
      real*8,intent(IN)         :: q_cyl, xyf(:)
#endif
#endif
      real*8,optional,intent(INOUT) :: tau(:,:)
      real*8                     :: b(Ndim)
#ifdef VORTICITY
      real*8                     :: kcoeff,exb(3)
      integer*4                  :: alpha,beta,ii
#endif
      integer*4                  :: i,j,k,inn,ign,ik
      integer*4,dimension(size(ind_asf))  :: ind_if,ind_jf,ind_kf
      real*8,dimension(neq,neq) :: A
      real*8,dimension(neq,Ndim):: APinch
#ifndef TEMPERATURE
      real*8                    :: auxvec(neq)
#endif
      real*8                    :: nn(3),qq(3,Neq),bb(3)
      real*8                    :: bn,kmult(size(ind_asf),size(ind_asf)),kmultf(size(ind_asf))
      real*8                    :: Qpr(Ndim,Neq)
      real*8                    :: W2(Neq),dW2_dU(Neq,Neq),QdW2(Ndim,Neq)
#ifdef TEMPERATURE
      real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
      real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
      real*8                    :: q_fs_i,q_fs_e,flux_limiter_i,flux_limiter_e,q_sh_i,q_sh_e
      real*8                    :: flux_limiter_i_ratio,flux_limiter_e_ratio, fl_deriv_i, fl_deriv_e
      real*8                    :: dq_fs_i_dU(Neq), dq_fs_e_dU(Neq)
      real*8                    :: W3(Neq),dW3_dU(Neq,Neq),QdW3(Ndim,Neq)
      real*8                    :: W4(Neq),dW4_dU(Neq,Neq),QdW4(Ndim,Neq)
#ifdef NEUTRALPNEW
      real*8                    :: W5p(Neq),dW5p_dU(Neq,Neq),QdW5p(Ndim,Neq)
#endif
      real*8                    :: Dnn_dU(Neq), Dnn_dU_U
      REAL*8                    :: neutral_limiter_phi,neutral_limiter_Gamma_max,neutral_limiter_ratio
      REAL*8                    :: neutral_limiter_Gamma_unlim(Ndim)
#ifdef KEQUATION
#ifdef DKLINEARIZED
      real*8                    :: ddk_dU(Neq), ddk_dU_U
      real*8                    :: gradddk(Ndim)
#endif
#endif
#ifdef NEUTRALGAMMA
      real*8                    :: Etan
      real*8                    :: Vun(Neq),dEtan_dU(Neq),gmGamman(Ndim)
      real*8                    :: dVun_dU(Neq,Neq),TauGamman(Ndim,Neq)
#endif
#endif

      inn = phys%idx_rhon_eq
      ign = phys%idx_gamman_eq
      ik = phys%idx_k_eq

      b = b3(1:Ndim)
      bb = b3
      ! Jacobian matrices
           bn = dot_PRODUCT(b,n)
      CALL jacobianMatrices(uf,A)

      ! Jacobian for pinch term
      APinch = 0.d0
      IF (switch%transport_1d) THEN
      CALL transport_model_1d%compute_1D_pinch_matrix(b,SQRT(MAX(psi,0.d0)),APinch)
    ENDIF

      ! Compute Q^T^(k-1)
           Qpr = RESHAPE(qf,(/Ndim,Neq/))

#ifdef TEMPERATURE
      neutral_limiter_phi = 1.d0
      IF (limiter_active) THEN
        CALL compute_neutral_flux_limiter(uf, qf, neutral_limiter_phi, neutral_limiter_Gamma_unlim, &
          &neutral_limiter_Gamma_max, neutral_limiter_ratio)
      ENDIF
#endif

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
#ifdef NEUTRALPNEW
      CALL compute_W5p(uf,W5p)
      CALL compute_dW5p_dU(uf,dW5p_dU)
           QdW5p = MATMUL(Qpr,dW5p_dU)
      IF (limiter_active) THEN
        W5p = neutral_limiter_phi*W5p
        dW5p_dU = neutral_limiter_phi*dW5p_dU
        QdW5p = neutral_limiter_phi*QdW5p
      ENDIF
#endif

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

      ! Compute flux limiters

      IF (switch%flux_limiter) THEN

        CALL compute_free_streaming_heat_flux_electrons(uf,q_fs_e)
        CALL compute_free_streaming_heat_flux_ions(uf,q_fs_i)

        ! Compute spitzer-härm heat flux
        q_sh_e = coefe*Alphae*gme
        q_sh_i = coefi*Alphai*gmi

        CALL compute_flux_limiter(q_fs_e,q_sh_e,phys%c_fle,flux_limiter_e)
        CALL compute_flux_limiter(q_fs_i,q_sh_i,phys%c_fli,flux_limiter_i)

        ! Derivative of flux limiter
        flux_limiter_e_ratio = ABS(q_sh_e)/q_fs_e/phys%c_fle
        flux_limiter_i_ratio = ABS(q_sh_i)/q_fs_i/phys%c_fli

        fl_deriv_e = q_sh_e*flux_limiter_e_ratio/q_fs_e
        fl_deriv_i = q_sh_i*flux_limiter_i_ratio/q_fs_i

        call compute_dfree_streaming_heat_flux_electrons_dU(uf,dq_fs_e_dU)
        call compute_dfree_streaming_heat_flux_ions_dU(uf,dq_fs_i_dU)

      ELSE
        flux_limiter_e = 1.
        flux_limiter_i = 1.
        fl_deriv_e = 0.
        fl_deriv_i = 0.
        dq_fs_e_dU = 0.
        dq_fs_i_dU = 0.
      END IF

           CALL compute_Dnn_dU(uf,Dnn_dU)
      Dnn_dU_u = dot_product(Dnn_dU,uf)
      IF (limiter_active) THEN
        Dnn_dU = neutral_limiter_phi*Dnn_dU
        Dnn_dU_u = neutral_limiter_phi*Dnn_dU_u
      ENDIF
#ifdef KEQUATION
#ifdef DKLINEARIZED
      call compute_ddk_dU(uf,xyf,q_cyl,ddk_dU)

      ddk_dU_u = dot_product(ddk_dU,uf)
#endif
#endif
#ifdef NEUTRALGAMMA
      CALL computeEtan(uf,Etan)
      CALL compute_dEtan_dU(uf,dEtan_dU)
      CALL computeVun(uf,Vun)
      CALL compute_dVun_dU(uf,dVun_dU)

      gmGamman = MATMUL(Qpr,Vun)
      TauGamman = MATMUL(Qpr,dVun_dU)
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
                    kmult = flux_limiter_i**2*(coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_PRODUCT(Taui(:,j),b))) + &
                    & fl_deriv_i*dq_fs_i_dU(j))*NNif*bn
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                    elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = flux_limiter_i**2*coefi*Alphai*Vveci(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
              elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
                 kmultf = flux_limiter_i**2*coefi*Alphai*(dot_PRODUCT(MATMUL(TRANSPOSE(Taui),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
        ELSEIF (i == 4) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
                    kmult = flux_limiter_e**2*(coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_PRODUCT(Taue(:,j),b)))+&
                    fl_deriv_e*dq_fs_e_dU(j))*NNif*bn
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                    elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = flux_limiter_e**2*coefe*Alphae*Vvece(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
              elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
                 kmultf = flux_limiter_e**2*coefe*Alphae*(dot_PRODUCT(MATMUL(TRANSPOSE(Taue),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
      ELSEIF (i == inn) THEN
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
#ifdef NEUTRALPNEW
        DO j = 1,Neq
          ind_jf = ind_asf + j
          DO k = 1,Ndim
            ind_kf = k + (j - 1)*Ndim + ind_ash
            kmult = W5p(j)*n(k)*NNif
            elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
            elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult

            kmult = QdW5p(k,j)*n(k)*NNif
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
          END DO
        END DO
        kmultf = dot_product(matmul(transpose(QdW5p),n),uf)*Nif
        elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
        elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
#endif
#ifdef NEUTRALGAMMA
       ELSEIF (i == ign) THEN
          DO j = 1,Neq
             ind_jf = ind_asf + j
             DO k = 1,Ndim
                ind_kf = k + (j - 1)*Ndim + ind_ash
                kmult = Etan*TauGamman(k,j)*NNif*n(k)
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%ALL(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
                kmult = Etan*Vun(j)*NNif*n(k)
                elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
                elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult
             END DO
          END DO
          kmultf = Etan*(dot_PRODUCT(TauGamman(1,:),uf)*n(1) + dot_PRODUCT(TauGamman(2,:),uf)*n(2))*Nif
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
#endif
#endif
       END IF
#ifdef KEQUATION
#ifdef DKLINEARIZED
       if (i .ne. inn) then
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

    ENDSUBROUTINE assemblyIntFacesContribution

    !********************************************************************
    !
    !         ASSEMBLY EXTERIOR FACES CONTRIBUTION
    !
    !********************************************************************

#ifdef DKLINEARIZED
    SUBROUTINE assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,&
      &ind_fg,b3,psi,q_cyl,xyf,n,diffiso,diffani,NNif,Nif,Nfbn,uf,qf,tau)
#else
    SUBROUTINE assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,&
        &ind_fg,b3,psi,n,diffiso,diffani,NNif,Nif,Nfbn,uf,qf,tau)
#endif
      integer*4,intent(IN)      :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
      logical                   :: isdir
      real*8,intent(IN)         :: b3(:),n(:), psi
      real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
      real*8,intent(IN)         :: NNif(:,:),Nif(:),Nfbn(:)
      real*8,intent(IN)         :: uf(:)
      real*8,intent(IN)         :: qf(:)
#ifdef KEQUATION
#ifdef DKLINEARIZED
      real*8,intent(IN)         :: q_cyl, xyf(:)
#endif
#endif
      real*8,optional,intent(INOUT) :: tau(:,:)
#ifdef VORTICITY
      integer*4                 :: alpha,beta,ii
      real*8                    :: exb(3),kcoeff
#endif
      integer*4                 :: i,j,k,inn,ign,ik
      integer*4,dimension(Npfl)  :: ind_if,ind_jf,ind_kf
      real*8,dimension(neq,neq) :: A
      real*8,dimension(neq,Ndim):: APinch
#ifndef TEMPERATURE
      real*8                    :: auxvec(neq)
#endif
      real*8                    :: bn,kmult(Npfl,Npfl),kmultf(Npfl)
      real*8                    :: Qpr(Ndim,Neq)
      real*8                    :: nn(3),qq(3,Neq),b(Ndim),bb(3)
      real*8                    :: W2(Neq), dW2_dU(Neq,Neq), QdW2(Ndim,Neq)
#ifdef TEMPERATURE
      real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
      real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
      real*8                    :: q_fs_i,q_fs_e,flux_limiter_i,flux_limiter_e,q_sh_e,q_sh_i
      real*8                    :: flux_limiter_i_ratio,flux_limiter_e_ratio, fl_deriv_i, fl_deriv_e
      real*8                    :: dq_fs_i_dU(Neq), dq_fs_e_dU(Neq)  
      real*8                    :: W3(Neq), dW3_dU(Neq,Neq), QdW3(Ndim,Neq)
      real*8                    :: W4(Neq), dW4_dU(Neq,Neq), QdW4(Ndim,Neq)
#ifdef NEUTRALPNEW
      real*8                    :: W5p(Neq), dW5p_dU(Neq,Neq), QdW5p(Ndim,Neq)
#endif
      real*8                    :: Dnn_dU(Neq), Dnn_dU_U
      REAL*8                    :: neutral_limiter_phi,neutral_limiter_Gamma_max,neutral_limiter_ratio
      REAL*8                    :: neutral_limiter_Gamma_unlim(Ndim)
#ifdef KEQUATION
#ifdef DKLINEARIZED
      real*8                    :: ddk_dU(Neq), ddk_dU_U
      real*8                    :: gradddk(Ndim)
#endif
#endif
#ifdef NEUTRALGAMMA
      real*8                    :: Etan
      real*8                    :: Vun(Neq),dEtan_dU(Neq),gmGamman(Ndim)
      real*8                    :: dVun_dU(Neq,Neq),TauGamman(Ndim,Neq)
#endif
#endif

      inn = phys%idx_rhon_eq
      ign = phys%idx_gamman_eq
      ik = phys%idx_k_eq

      b = b3(1:Ndim)
      bb = b3
      ! Jacobian matrices
           bn = dot_PRODUCT(b,n)
      CALL jacobianMatrices(uf,A)

      ! Jacobian matrices Pinch
      APinch = 0.d0
      IF (switch%transport_1d) THEN
      CALL transport_model_1d%compute_1D_pinch_matrix(b,SQRT(MAX(psi,0.d0)),APinch)
    ENDIF

      ! Compute Q^T^(k-1)
           Qpr = RESHAPE(qf,(/Ndim,Neq/))

#ifdef TEMPERATURE
      neutral_limiter_phi = 1.d0
      IF (limiter_active) THEN
        CALL compute_neutral_flux_limiter(uf, qf, neutral_limiter_phi, neutral_limiter_Gamma_unlim, &
          &neutral_limiter_Gamma_max, neutral_limiter_ratio)
      ENDIF
#endif

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
#ifdef NEUTRALPNEW
      CALL compute_W5p(uf,W5p)
      CALL compute_dW5p_dU(uf,dW5p_dU)
           QdW5p = MATMUL(Qpr,dW5p_dU)
      IF (limiter_active) THEN
        W5p = neutral_limiter_phi*W5p
        dW5p_dU = neutral_limiter_phi*dW5p_dU
        QdW5p = neutral_limiter_phi*QdW5p
      ENDIF
#endif

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

      ! Compute flux limiters
      IF (switch%flux_limiter) THEN
        CALL compute_free_streaming_heat_flux_electrons(uf,q_fs_e)
        CALL compute_free_streaming_heat_flux_ions(uf,q_fs_i)

        ! compute spitzer-harm heat flux 
        q_sh_e = coefe*Alphae*gme
        q_sh_i = coefi*Alphai*gmi

        CALL compute_flux_limiter(q_fs_e,q_sh_e,phys%c_fle, flux_limiter_e)
        CALL compute_flux_limiter(q_fs_i,q_sh_i,phys%c_fli, flux_limiter_i)

        ! Derivative of flux limiter
        flux_limiter_e_ratio = ABS(q_sh_e)/(q_fs_e)/phys%c_fle
        flux_limiter_i_ratio = ABS(q_sh_i)/(q_fs_i)/phys%c_fli

        fl_deriv_e = q_sh_e*flux_limiter_e_ratio/q_fs_e
        fl_deriv_i = q_sh_i*flux_limiter_i_ratio/q_fs_i

        call compute_dfree_streaming_heat_flux_electrons_dU(uf,dq_fs_e_dU)
        call compute_dfree_streaming_heat_flux_ions_dU(uf,dq_fs_i_dU)
      
      ELSE
        flux_limiter_e = 1.
        flux_limiter_i = 1.
        fl_deriv_e = 0.
        fl_deriv_i = 0.
        dq_fs_e_dU = 0.
        dq_fs_i_dU = 0.
      END IF
      
           CALL compute_Dnn_dU(uf,Dnn_dU)
      Dnn_dU_u = dot_product(Dnn_dU,uf)
      IF (limiter_active) THEN
        Dnn_dU = neutral_limiter_phi*Dnn_dU
        Dnn_dU_u = neutral_limiter_phi*Dnn_dU_u
      ENDIF

#ifdef KEQUATION
#ifdef DKLINEARIZED
      call compute_ddk_dU(uf,xyf,q_cyl,ddk_dU)

      ddk_dU_u = dot_product(ddk_dU,uf)
#endif
#endif
#ifdef NEUTRALGAMMA
      CALL computeEtan(uf,Etan)
      CALL compute_dEtan_dU(uf,dEtan_dU)
      CALL computeVun(uf,Vun)
      CALL compute_dVun_dU(uf,dVun_dU)

      gmGamman = MATMUL(Qpr,Vun)
      TauGamman = MATMUL(Qpr,dVun_dU)
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
                       kmult = flux_limiter_i**2*(coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_PRODUCT(Taui(:,j),b)))+&
                       &fl_deriv_i*dq_fs_i_dU(j))*NNif*bn
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            END IF
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = flux_limiter_i**2*coefi*Alphai*Vveci(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
                 kmultf = flux_limiter_i**2*coefi*Alphai*(dot_PRODUCT(MATMUL(TRANSPOSE(Taui),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
        ELSEIF (i == 4) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
                    IF (.NOT. isdir) THEN
                       kmult = flux_limiter_e**2*(coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_PRODUCT(Taue(:,j),b)))+&
                       &fl_deriv_e*dq_fs_e_dU(j))*NNif*bn
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            END IF
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = flux_limiter_e**2*coefe*Alphae*Vvece(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
                 kmultf = flux_limiter_e**2*coefe*Alphae*(dot_PRODUCT(MATMUL(TRANSPOSE(Taue),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
        ELSEIF (i == inn) THEN
            DO j=1,Neq
              ind_jf = ind_asf+j
              DO k = 1,Ndim

                kmult = Dnn_dU(j)*Qpr(k,i)*n(k)*NNif
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                    ENDDO
                 ENDDO
            kmultf = Dnn_dU_U*(Qpr(1,i)*n(1)+Qpr(2,i)*n(2))*Nif
            elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
#ifdef NEUTRALPNEW
            DO j = 1,Neq
              ind_jf = ind_asf + j
              DO k = 1,Ndim
                ind_kf = k + (j - 1)*Ndim + ind_ash
                kmult = W5p(j)*n(k)*NNif
                elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
                IF (.NOT. isdir) THEN
                  elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - QdW5p(k,j)*n(k)*NNif
                END IF
              END DO
            END DO
            kmultf = dot_product(matmul(transpose(QdW5p),n),uf)*Nif
            elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
#endif
#ifdef NEUTRALGAMMA
       ELSEIF (i == ign) THEN
          DO j = 1,Neq
             ind_jf = ind_asf + j
             DO k = 1,Ndim
                ind_kf = k + (j - 1)*Ndim + ind_ash
                kmult = Etan*TauGamman(k,j)*NNif*n(k)
                IF (.not. isdir) THEN
                  elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                END IF
                kmult = Etan*Vun(j)*NNif*n(k)
                elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
             END DO
          END DO
          kmultf = Etan*(dot_PRODUCT(TauGamman(1,:),uf)*n(1) + dot_PRODUCT(TauGamman(2,:),uf)*n(2))*Nif
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
#endif
#endif
      END IF

#ifdef KEQUATION
#ifdef DKLINEARIZED
        if (i .ne. inn) then
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
    ENDSUBROUTINE assemblyExtFacesContribution

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

    ENDSUBROUTINE do_assembly

#ifdef NEUTRAL
  !********************************************************************
  !
  !         ASSEMBLY NEUTRAL SOURCE MATRIX
  !
  !********************************************************************
#ifdef TEMPERATURE
#ifdef NEUTRALGAMMA
  SUBROUTINE assemblyNeutral(U,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
      &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,fGammaN,dfGammaN_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
      &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,fEiN,dfEiN_dU,Sn,Sn0,&
      sigmavEiz,dsigmavEiz_dU,sigmavErec,dsigmavErec_dU,cooling_factor,dcooling_factor_dU)
#else
  SUBROUTINE assemblyNeutral(U,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
      &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
      &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,Sn,Sn0,&
      sigmavEiz,dsigmavEiz_dU,sigmavErec,dsigmavErec_dU,cooling_factor,dcooling_factor_dU)
#endif
#else
#ifdef NEUTRALGAMMA
    SUBROUTINE assemblyNeutral(U,niz,dniz_dU,nrec,dnrec_dU,fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,fGammaN,dfGammaN_dU,Sn,Sn0)
#else
    SUBROUTINE assemblyNeutral(U,niz,dniz_dU,nrec,dnrec_dU,fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,Sn,Sn0)
#endif
#endif
             REAL*8, INTENT(IN) :: niz,nrec,fGammacx,fGammarec
             REAL*8, INTENT(IN) :: U(:),dniz_dU(:),dnrec_dU(:),dfGammacx_dU(:),dfGammarec_dU(:)
#ifdef NEUTRALGAMMA
             REAL*8, INTENT(IN) :: fGammaN
             REAL*8, INTENT(IN) :: dfGammaN_dU(:)
#endif
#ifndef TEMPERATURE
             REAL*8             :: sigmaviz,sigmavrec,sigmavcx
#else
      REAL*8, INTENT(IN)        :: sigmaviz,sigmavrec,sigmavcx,fEiiz,fEirec,fEicx      
      REAL*8, INTENT(IN)        :: dsigmaviz_dU(:),dsigmavrec_dU(:),dsigmavcx_dU(:)
      REAL*8, INTENT(IN)        :: dfEiiz_dU(:),dfEirec_dU(:),dfEicx_dU(:)
#ifdef NEUTRALGAMMA
      REAL*8, INTENT(IN)        :: fEiN
      REAL*8, INTENT(IN)        :: dfEiN_dU(:)
#endif
      REAL*8, INTENT(IN), OPTIONAL :: sigmavEiz,sigmavErec,dsigmavEiz_dU(:),dsigmavErec_dU(:)
      REAL*8, INTENT(IN), OPTIONAL :: cooling_factor,dcooling_factor_dU(:)
#endif
             REAL*8             :: RE,Sn(:,:),Sn0(:)
#ifdef TEMPERATURE
             REAL*8             :: recombination_energy
#endif
             INTEGER*4          :: ign, inn

      Sn   = 0.
      Sn0  = 0.
      RE   = 0.
      inn  = phys%idx_rhon_eq
      ign  = phys%idx_gamman_eq
#ifdef TEMPERATURE
      recombination_energy = neutral_rt%recombination_energy
#endif

#ifndef TEMPERATURE
      sigmaviz   = 3.01e-14*simpar%refval_density*simpar%refval_time
      sigmavrec  = 1.3638e-20*simpar%refval_density*simpar%refval_time
      sigmavcx   = 4.0808e-15*simpar%refval_density*simpar%refval_time
#endif


      !Assembly Source Terms in plasma density equation
      Sn(1,:)   = -dniz_dU(:)*sigmaviz + dnrec_dU(:)*sigmavrec
#ifdef TEMPERATURE

      Sn(1,:)   = Sn(1,:) - niz*dsigmaviz_dU(:) + nrec*dsigmavrec_dU(:)
#endif
      !Assembly Source Terms in plasma momentum equation

      Sn(2,:) = dfGammacx_dU(:)*sigmavcx + dfGammarec_dU(:)*sigmavrec
#ifdef NEUTRALGAMMA
      Sn(2,:) = Sn(2,:) - (dfGammaN_dU(:)*sigmaviz + dfGammaN_dU(:)*sigmavcx)
#endif
#ifdef TEMPERATURE
      Sn(2,:)   = Sn(2,:) + fGammacx*dsigmavcx_dU(:) + fGammarec*dsigmavrec_dU(:)
#ifdef NEUTRALGAMMA
      Sn(2,:)   = Sn(2,:) - (fGammaN*dsigmaviz_dU(:) + fGammaN*dsigmavcx_dU(:))
#endif

      !Assembly Source Terms in ion energy equation

      Sn(3,:) = -RE*dfEiiz_dU(:)*sigmaviz + dfEirec_dU(:)*sigmavrec + dfEicx_dU(:)*sigmavcx
      Sn(3,:) = Sn(3,:) - RE*fEiiz*dsigmaviz_dU(:) + fEirec*dsigmavrec_dU(:) + fEicx*dsigmavcx_dU(:)
#ifdef NEUTRALGAMMA
      Sn(3,:) = Sn(3,:) - (dfEiN_dU(:)*sigmaviz + dfEiN_dU(:)*sigmavcx + &
           &fEiN*dsigmaviz_dU(:) + fEiN*dsigmavcx_dU(:))
#endif
      !Assembly Source Terms in electron energy equation
      Sn(4,:) = dniz_dU(:)*sigmavEiz + niz*dsigmavEiz_dU(:) + &
        &dnrec_dU(:)*sigmavErec + nrec*dsigmavErec_dU(:)
      IF (PRESENT(cooling_factor)) THEN
        Sn(4,:) = Sn(4,:) + phys%impurity_concentration*(nrec*dcooling_factor_dU(:) + dnrec_dU(:)*cooling_factor)
      endif

      Sn(4,:) = Sn(4,:) - recombination_energy*(dnrec_dU(:)*sigmavrec + nrec*dsigmavrec_dU(:))


#endif
      !Assembly Source Terms in neutral density equation
      Sn(inn,:) = -Sn(1,:)
#ifdef NEUTRALGAMMA
      Sn(ign,:) = -Sn(2,:)
#endif

      !Assembly RHS Neutral Source Terms
      Sn0(1)    = niz*sigmaviz - nrec*sigmavrec
      Sn0(2)    = -fGammacx*sigmavcx - fGammarec*sigmavrec
#ifdef NEUTRALGAMMA
      Sn0(2)    = Sn0(2) + fGammaN*sigmaviz + fGammaN*sigmavcx
#endif
#ifdef TEMPERATURE
      Sn0(1)    = Sn0(1) + niz*dot_PRODUCT(dsigmaviz_dU,U) - nrec*dot_PRODUCT(dsigmavrec_dU,U)
      Sn0(2)    = Sn0(2) - fGammarec*dot_PRODUCT(dsigmavrec_dU,U)
#ifdef NEUTRALGAMMA
      Sn0(2)    = Sn0(2) + fGammaN*dot_PRODUCT(dsigmaviz_dU,U)
#endif
      Sn0(3)    = RE*fEiiz*sigmaviz - fEirec*sigmavrec - fEicx*sigmavcx
#ifdef NEUTRALGAMMA
      Sn0(3)    = Sn0(3) + fEiN*sigmaviz + fEiN*sigmavcx
#endif
      Sn0(4)    = nrec*sigmavrec*recombination_energy
      Sn0(4)    = Sn0(4) - niz*sigmavEiz - nrec*sigmavErec
      Sn0(3)    = Sn0(3) + RE*fEiiz*dot_PRODUCT(dsigmaviz_dU,U) - fEirec*dot_PRODUCT(dsigmavrec_dU,U)
#ifdef NEUTRALGAMMA
      Sn0(3)    = Sn0(3) + fEiN*dot_PRODUCT(dsigmaviz_dU,U)
#endif
      Sn0(4)    = Sn0(4) - niz*dot_product(dsigmavEiz_dU,U) - nrec*dot_product(dsigmavErec_dU,U)
      Sn0(4)    = Sn0(4) + nrec*dot_PRODUCT(dsigmavrec_dU,U)*recombination_energy
      IF (PRESENT(cooling_factor)) THEN
        Sn0(4)    = Sn0(4) - phys%impurity_concentration*nrec*cooling_factor
      ENDIF
#endif
      Sn0(inn)  = -Sn0(1)
#ifdef NEUTRALGAMMA
      Sn0(ign)  = -Sn0(2)
#endif

    ENDSUBROUTINE assemblyNeutral
#endif


         ENDSUBROUTINE hdg_ComputeJacobian
