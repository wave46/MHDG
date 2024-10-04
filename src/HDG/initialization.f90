!************************************************************
! project: MHDG
! file: initialization.f90
! date: 03/01/2017
! Initialization of the physics depending on the model,of
! the elemental matrices,of the solution
!************************************************************
MODULE initialization
  USE globals
  USE printutils
  USE analytical
  USE physics
  USE MPI_OMP
  USE LinearAlgebra, ONLY: col,tensorProduct,solve_linear_system


  IMPLICIT NONE
CONTAINS

  !*********************************************
  ! Initialization of the simulation parameters
  !*********************************************
  SUBROUTINE init_sim(nts,dt)

    INTEGER,INTENT(out) :: nts
    REAL,INTENT(out) :: dt

    CALL initPhys()

    ! the time is set to zero here but in case of a restart this value is
    ! overwritten by the loaded solution final time
    time%t = 0.

    ! time step is initialized to dt0.
    dt = time%dt0
    time%dt = dt

    ! time count
    time%it = 0
    time%ik = 0

    ! number of time step of the simulation: if it is a steady state simulation
    ! ndt is set to 1,otherwise the input value is used
    IF (switch%steady) THEN
       nts = 1
    ELSE
       nts = time%nts
    END IF

    ! Allocate and initialize time residual
    ALLOCATE (sol%tres(nts))
    ALLOCATE (sol%time(nts))
    sol%tres = 0.
    sol%time = 0.
    sol%Nt = 0
  END SUBROUTINE init_sim

  !*******************************************
  ! Initialization of the elemental matrices
  !*******************************************
  SUBROUTINE init_elmat
    INTEGER :: Neq,Ndim,Nel,Np,Nfg,Nf

#ifdef TOR3D
    INTEGER :: Ntorloc,N2d,Np1d,Np2d,Nfl
#endif

    Neq = phys%Neq                      ! N. of equations
#ifdef TOR3D
    Ndim = 3                             ! N. of dimensions
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
    N2d = Mesh%Nelems                   ! N. of 2D elements
    Nel = N2d*ntorloc                   ! N. of 3D elements
    Np1d = refElTor%Nnodes1D             ! N. of nodes for each toroidal 1d element
    Np2d = refElPol%Nnodes2D             ! N. of nodes for each poloidal 2D element
    Np = Np2d*Np1d                     ! N. of nodes for each 3D element
    Nfl = refElPol%Nnodes1D*Np1d        ! N. of nodes in the lateral faces
    Nfg = Np2d*2 + refElPol%Nfaces*Nfl    ! N. of nodes in all the faces of a 3D element
    Nf = Mesh%Nfaces                   ! N. of faces in the 2D mesh
#else
    Ndim = 2
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfg = refElPol%Nfacenodes*Nf
#endif

    ALLOCATE (elmat%iAqq(Neq*Ndim*Np,Neq*Ndim*Np,Nel))
    ALLOCATE (elmat%Aqu(Neq*Ndim*Np,Neq*Np,Nel))
    ALLOCATE (elmat%Aql(Neq*Ndim*Np,Neq*Nfg,Nel))
    ALLOCATE (elmat%Auq(Neq*Np,Ndim*Neq*Np,Nel))
    ALLOCATE (elmat%Auu(Neq*Np,Neq*Np,Nel))
    ALLOCATE (elmat%Aul(Neq*Np,Neq*Nfg,Nel))
    ALLOCATE (elmat%Alq(Neq*Nfg,Neq*Ndim*Np,Nel))
    ALLOCATE (elmat%Alu(Neq*Nfg,Neq*Np,Nel))
    ALLOCATE (elmat%ALL(Neq*Nfg,Neq*Nfg,Nel))
    ALLOCATE (elmat%Aql_dir(Neq*Np*Ndim,Nel))
    ALLOCATE (elmat%Aul_dir(Neq*Np,Nel))
    ALLOCATE (elmat%S(Neq*Np,Nel))
    ALLOCATE (elmat%fH(Neq*Nfg,Nel))

    IF (switch%ME) THEN
       ALLOCATE(phys%puff_exp(time%nts))
    END IF

    IF (switch%saveTau) THEN
       ALLOCATE(phys%diff_nn_Vol(Mesh%Nelems*refElPol%NGauss2D))
       ALLOCATE(phys%diff_nn_Fac(Mesh%Nelems*refElPol%Nfaces*refElPol%NGauss1D))
       ALLOCATE(phys%diff_nn_Bou(Mesh%Nextfaces*refElPol%NGauss1D))
       ALLOCATE(phys%v_nn_Vol(Mesh%Nelems*refElPol%NGauss2D,Mesh%Ndim))
       ALLOCATE(phys%v_nn_Fac(Mesh%Nelems*refElPol%Nfaces*refElPol%NGauss1D,Mesh%Ndim))
       ALLOCATE(phys%v_nn_Bou(Mesh%Nextfaces*refElPol%NGauss1D,Mesh%Ndim))
       ALLOCATE(Mesh%Xg(Mesh%Nelems*refElPol%NGauss2D,Mesh%Ndim))
       ALLOCATE(Mesh%Xgf(Mesh%Nelems*refElPol%Nfaces*refElPol%NGauss1D,Mesh%Ndim))
       ALLOCATE(Mesh%Xgb(Mesh%Nextfaces*refElPol%NGauss1D,Mesh%Ndim))
    END IF

    elmat%iAqq = 0.
    elmat%Aqu = 0.
    elmat%Aql = 0.
    elmat%Auq = 0.
    elmat%Auu = 0.
    elmat%Aul = 0.
    elmat%Alq = 0.
    elmat%Alu = 0.
    elmat%All = 0.
    elmat%Aql_dir = 0.
    elmat%Aul_dir = 0.
    elmat%S = 0.
    elmat%fh = 0.
    ALLOCATE (elmat%UU(Neq*Np,Neq*Nfg,Nel))
    ALLOCATE (elmat%U0(Neq*Np,Nel))
    ALLOCATE (elmat%LL(Neq*Np*Ndim,Neq*Nfg,Nel))
    ALLOCATE (elmat%L0(Neq*Np*Ndim,Nel))
    elmat%UU = 0.d0
    elmat%U0 = 0.d0
    elmat%LL = 0.d0
    elmat%L0 = 0.d0
    IF (switch%ME) THEN
       phys%puff_exp = 0.
    END IF
    IF (switch%saveTau) THEN
       phys%diff_nn_Vol = 0.
       phys%diff_nn_Fac = 0.
       phys%diff_nn_Bou = 0.
       phys%v_nn_Vol = 0.
       phys%v_nn_Fac = 0.
       phys%v_nn_Bou = 0.
       Mesh%Xg = 0.
       Mesh%Xgf = 0.
       Mesh%Xgb = 0.
    END IF

  ENDSUBROUTINE init_elmat

  SUBROUTINE init_solve_timing
    USE globals
    timing%cputpre   =1.e-8
    timing%cputmap   =1.e-8
    timing%cputass   =1.e-8
    timing%cputbcd   =1.e-8
    timing%cputsol   =1.e-8
    timing%cputjac   =1.e-8
    timing%cputglb   =1.e-8
    timing%cputcom   =1.e-8
    timing%cputadapt =1.e-8

    timing%runtpre=1.e-8
    timing%runtmap=1.e-8
    timing%runtass=1.e-8
    timing%runtbcd=1.e-8
    timing%runtsol=1.e-8
    timing%runtjac=1.e-8
    timing%runtglb=1.e-8
    timing%runtcom=1.e-8
    timing%runtadapt=1.e-8

    timing%clstime1=1.e-8
    timing%clstime2=1.e-8
    timing%clstime3=1.e-8
    timing%clstime4=1.e-8
    timing%clstime5=1.e-8
    timing%clstime6=1.e-8

    timing%rlstime1=1.e-8
    timing%rlstime2=1.e-8
    timing%rlstime3=1.e-8
    timing%rlstime4=1.e-8
    timing%rlstime5=1.e-8
    timing%rlstime6=1.e-8
  END SUBROUTINE init_solve_timing

  !********************************
  ! Initialization of the solution
  !********************************
  SUBROUTINE init_sol
    INTEGER :: Neq,Nel,Ndim,Np,Nf,Nfg,sizeutilde,sizeu
    INTEGER :: Ngvo

#ifdef TOR3D
    INTEGER :: N2d, Nfl,Ntorloc,Np1d,Np2d
#endif

    Neq = phys%Neq
#ifdef TOR3D
    Ndim = 3                             ! N. of dimensions
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
    N2d = Mesh%Nelems                   ! N. of 2D elements
    Nel = N2d*ntorloc                   ! N. of 3D elements
    Np1d = refElTor%Nnodes1D             ! N. of nodes for each toroidal 1d element
    Np2d = refElPol%Nnodes2D             ! N. of nodes for each poloidal 2D element
    Np = Np2d*Np1d                     ! N. of nodes for each 3D element
    Nfl = refElPol%Nnodes1D*Np1d        ! N. of nodes in the lateral faces
    Nfg = Np2d*2 + refElPol%Nfaces*Nfl    ! N. of nodes in all the faces of a 3D element
    Nf = Mesh%Nfaces                   ! N. of faces in the 2D mesh
    sizeu = Neq*Nel*Np                    ! Size of u
    Ngvo = refElPol%Ngauss2d*refEltor%Ngauss1d
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d) + Neq*Np2d*N2d! Size of utilde
    ELSE
       sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d)! Size of utilde
    ENDIF
#else
    sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d)! Size of utilde
#endif
#else
    Ndim = 2
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfg = refElPol%Nfacenodes*Nf
    sizeu = Neq*Nel*Np
    sizeutilde = Neq*Mesh%Nfaces*Mesh%Nnodesperface
    Ngvo = refElPol%Ngauss2d
#endif

    ! Allocation of the solution vector
    ALLOCATE (sol%u(sizeu))
    ALLOCATE (sol%u_tilde(sizeutilde))
    ALLOCATE (sol%u_tilde0(sizeutilde))
    ALLOCATE (sol%q(sizeu*Ndim))
    sol%u = 0.
    sol%u_tilde = 0.
    sol%q = 0.
    ! Initialize the solution
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6,*) "*** Initializing the solution"
       END IF
    ENDIF
    IF (switch%init.EQ.1) THEN
       ! The solution is intialized in each node to the analytical solution
       IF (MPIvar%glob_id .EQ. 0) THEN
          IF (utils%printint > 0) THEN
             WRITE (6,*) "******* Initializing the solution to the analytic solution"
          END IF
       ENDIF
       CALL init_sol_analytic()
    ELSEIF (switch%init.EQ.2) THEN
       ! The solution is intialized in each node to the analytical solution
       IF (MPIvar%glob_id .EQ. 0) THEN
          IF (utils%printint > 0) THEN
             WRITE (6,*) "******* Initializing the solution with L2 projection"
          END IF
       ENDIF
       CALL init_sol_l2proj()
    ELSE
       WRITE(6,*) "Wrong initialization type"
       STOP
    ENDIF
    ! Extract the face solution from the elemental one
    CALL extractFaceSolution()

    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6,*) "Done!"
       END IF
    ENDIF
  CONTAINS
    !***********************************************************
    ! Initialization of the solution using the analytic solution
    !***********************************************************
    SUBROUTINE init_sol_analytic
      INTEGER             :: iel,i
      INTEGER             :: ind(Np)
      REAL*8              :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
      REAL*8              :: ue(Np,Neq)
      REAL*8,ALLOCATABLE  :: u(:,:)
      REAL*8,ALLOCATABLE  :: qx(:,:),qy(:,:),auxq(:,:)
      REAL*8              :: uex(Np,Neq),uey(Np,Neq)

#ifdef TOR3D
      INTEGER             :: iel3, itor,itorg
      REAL*8              :: htor,tel(refElTor%Nnodes1d), tdiv(numer%ntor + 1), uet(Np,Neq)
      REAL*8,ALLOCATABLE  :: qt(:,:)
#endif


      ALLOCATE (u(Nel*Np,phys%Neq))
      u = 0.
      ALLOCATE (qx(Nel*Np,phys%Neq))
      ALLOCATE (qy(Nel*Np,phys%Neq))
      ALLOCATE (auxq(Nel*Np*phys%Neq,Ndim))
      qx = 0.; qy = 0.
#ifdef TOR3D
      ALLOCATE (qt(Nel*Np,phys%Neq))
      qt = 0.
#endif

#ifdef TOR3D

      !****************************************
      !          3D
      !****************************************
      htor = numer%tmax/numer%ntor
      tdiv = 0.
      DO i = 1,numer%ntor
         tdiv(i + 1) = i*htor
      END DO

      DO itor = 1,ntorloc
#ifdef PARALL
         itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
         IF (itorg == numer%ntor + 1) itorg = 1
#else
         itorg = itor
#endif
         tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
         DO iel = 1,Mesh%Nelems
            iel3 = (itor - 1)*N2d+iel
            ind = (iel3 - 1)*Np + (/(i,i=1,Np)/)
            Xe = Mesh%X(Mesh%T(iel,:),:)
            CALL analytical_solution(Xe(:,1),Xe(:,2),tel,ue)
            CALL analytical_gradient(Xe(:,1),Xe(:,2),tel,ue,uex,uey,uet)
            qx(ind,:) = uex
            qy(ind,:) = uey
            qt(ind,:) = uet
            u(ind,:) = ue
         END DO
      END DO
#else
      !****************************************
      !          2D
      !****************************************
      DO iel = 1,Mesh%Nelems

         ind = (iel - 1)*Np + (/(i,i=1,Np)/)
         Xe = Mesh%X(Mesh%T(iel,:),:)
         CALL analytical_solution(iel,Xe(:,1),Xe(:,2),ue)
         CALL analytical_gradient(Xe(:,1),Xe(:,2),ue,uex,uey)
         qx(ind,:) = uex
         qy(ind,:) = uey
         u(ind,:) = ue
      END DO
#endif

      !****************************************
      !          common
      !****************************************
      sol%u = RESHAPE(TRANSPOSE(u),(/Nel*Np*phys%Neq/))

      auxq(:,1) = RESHAPE(TRANSPOSE(qx),(/Nel*Np*phys%Neq/))
      auxq(:,2) = RESHAPE(TRANSPOSE(qy),(/Nel*Np*phys%Neq/))
#ifdef TOR3D
      auxq(:,3) = RESHAPE(TRANSPOSE(qt),(/Nel*Np*phys%Neq/))
#endif
      sol%q = RESHAPE(TRANSPOSE(auxq),(/Nel*Np*phys%Neq*Ndim/))
      DEALLOCATE (qx,qy,auxq)
#ifdef TOR3D
      DEALLOCATE (qt)
#endif
      DEALLOCATE (u)
    END SUBROUTINE init_sol_analytic

    !***********************************************************
    ! Initialization of the solution using an L2 projection
    !***********************************************************
    SUBROUTINE init_sol_l2proj
      REAL*8,ALLOCATABLE    :: u(:,:)
      REAL*8,ALLOCATABLE    :: qx(:,:),qy(:,:),auxq(:,:)

#ifdef TOR3D
      REAL*8, ALLOCATABLE   :: qt(:,:)
#else
      INTEGER               :: iel,i, g
      INTEGER               :: ind(Np)
      REAL*8                :: dvolu
      REAL*8                :: Xe(Mesh%Nnodesperelem,Mesh%Ndim), xyg(Ngvo,2)
      REAL*8                :: detJ(Ngvo), J11(Ngvo),J12(Ngvo), J21(Ngvo),J22(Ngvo), iJ21(Ngvo),iJ22(Ngvo), iJ11(Ngvo),iJ12(Ngvo)
      REAL*8                :: M(Np,Np),rhs_u(Np,Neq),rhs_ux(Np,Neq),rhs_uy(Np,Neq), ue(Np,Neq), uex(Np,Neq),uey(Np,Neq),ug(Ngvo,Neq),ugx(Ngvo,Neq),ugy(Ngvo,Neq)
      REAL*8,PARAMETER      :: tol = 1e-12
#endif
      ALLOCATE (u(Nel*Np,phys%Neq))
      u = 0.
      ALLOCATE (qx(Nel*Np,phys%Neq))
      ALLOCATE (qy(Nel*Np,phys%Neq))
      ALLOCATE (auxq(Nel*Np*phys%Neq,Ndim))
      qx = 0.; qy = 0.
#ifdef TOR3D
      ALLOCATE (qt(Nel*Np,phys%Neq))
      qt = 0.
#endif

#ifdef TOR3D
      WRITE(6,*) "Not coded yet"
      STOP
      !         !****************************************
      !         !          3D
      !         !****************************************
      !         htor = numer%tmax/numer%ntor
      !         tdiv = 0.
      !         DO i = 1,numer%ntor
      !            tdiv(i + 1) = i*htor
      !         END DO

      !         DO itor = 1,ntorloc
      !#ifdef PARALL
      !            itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
      !            if (itorg == numer%ntor + 1) itorg = 1
      !#else
      !            itorg = itor
      !#endif
      !            tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
      !            DO iel = 1,Mesh%Nelems
      !               iel3 = (itor - 1)*N2d+iel
      !               ind = (iel3 - 1)*Np + (/(i,i=1,Np)/)
      !               Xe = Mesh%X(Mesh%T(iel,:),:)
      !               CALL analytical_solution(Xe(:,1),Xe(:,2),tel,ue)
      !               CALL analytical_gradient(Xe(:,1),Xe(:,2),tel,ue,uex,uey,uet)
      !               qx(ind,:) = uex
      !               qy(ind,:) = uey
      !               qt(ind,:) = uet
      !               u(ind,:) = ue
      !            END DO
      !         END DO
#else
      !****************************************
      !          2D
      !****************************************
      DO iel = 1,Mesh%Nelems

         ind = (iel - 1)*Np + (/(i,i=1,Np)/)
         Xe = Mesh%X(Mesh%T(iel,:),:)
         J11 = MATMUL(refElPol%Nxi2D,Xe(:,1))                           ! ng x 1
         J12 = MATMUL(refElPol%Nxi2D,Xe(:,2))                           ! ng x 1
         J21 = MATMUL(refElPol%Neta2D,Xe(:,1))                          ! ng x 1
         J22 = MATMUL(refElPol%Neta2D,Xe(:,2))                          ! ng x 1
         detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
         iJ11 = J22/detJ
         iJ12 = -J12/detJ
         iJ21 = -J21/detJ
         iJ22 = J11/detJ

         ! Solution at Gauss points
         xyg = MATMUL(refElPol%N2D,Xe)
         CALL analytical_solution(iel,xyg(:,1),xyg(:,2),ug)
         CALL analytical_gradient(xyg(:,1),xyg(:,2),ug,ugx,ugy)


         ! Initialize mass matrix and rhs
         M=0.
         rhs_u=0.
         rhs_ux=0.
         rhs_uy=0.
         DO g=1,Ngvo


            IF (detJ(g) < tol) THEN
               error STOP "Negative jacobian"
            END IF

            ! Integration weight
            dvolu = refElPol%gauss_weights2D(g)*detJ(g)
            IF (switch%axisym) THEN
               dvolu = dvolu*xyg(g,1)
            END IF


            M = M + TensorProduct(refElPol%N2D(g,:),refElPol%N2D(g,:))*dvolu
            rhs_u = rhs_u+tensorProduct(refElPol%N2D(g,:),ug(g,:))*dvolu
            rhs_ux = rhs_ux+tensorProduct(refElPol%N2D(g,:),ugx(g,:))*dvolu
            rhs_uy = rhs_uy+tensorProduct(refElPol%N2D(g,:),ugy(g,:))*dvolu
         END DO
         CALL solve_linear_system(M,rhs_u,ue)
         CALL solve_linear_system(M,rhs_ux,uex)
         CALL solve_linear_system(M,rhs_uy,uey)
         qx(ind,:) = uex
         qy(ind,:) = uey
         u(ind,:) = ue
      END DO

#endif

      !****************************************
      !          common
      !****************************************
      sol%u = RESHAPE(TRANSPOSE(u),(/Nel*Np*phys%Neq/))

      auxq(:,1) = RESHAPE(TRANSPOSE(qx),(/Nel*Np*phys%Neq/))
      auxq(:,2) = RESHAPE(TRANSPOSE(qy),(/Nel*Np*phys%Neq/))
#ifdef TOR3D
      auxq(:,3) = RESHAPE(TRANSPOSE(qt),(/Nel*Np*phys%Neq/))
#endif
      sol%q = RESHAPE(TRANSPOSE(auxq),(/Nel*Np*phys%Neq*Ndim/))
      DEALLOCATE (qx,qy,auxq)
#ifdef TOR3D
      DEALLOCATE (qt)
#endif
      DEALLOCATE (u)
    ENDSUBROUTINE init_sol_l2proj


  ENDSUBROUTINE init_sol



  !     subroutine reset_variables()
  !         integer      :: iel,i,iface,ifa,ieq
  !         integer      :: Np,Neq,Nf,Nel,Nut
  !         integer      :: ind(Mesh%Nnodesperelem)
  !         real*8       :: Xe(Mesh%Nnodesperelem,Mesh%Ndim),Xf(Mesh%Nnodesperface,Mesh%Ndim)
  !         real*8       :: ue(refEl%Nnodes2D,phys%Neq)
  !         integer      :: ind_uf(Mesh%Nnodesperface),faceNodes(Mesh%Nnodesperface)
  !         real*8,allocatable :: u(:,:)
  !         real*8,allocatable :: qx(:,:),qy(:,:),auxq(:,:)
  !         real*8,allocatable :: u_tilde(:,:)




  !!         integer :: neq,Ne,Nf,Nfe,unkF,Np,Nfp,iElem,ifa,iFace,i
  !!         integer :: ind_uf(1:Mesh%Nnodesperface),faceNodes(1:Mesh%Nnodesperface)
  !!         integer :: ind_ue(1:Mesh%Nnodesperelem)
  !!         logical :: alreadydone(1:Mesh%ukf)
  !
  !!         real*8,allocatable :: u(:,:)
  !!         real*8              :: tdiv(numer%ntor + 1)
  !!#ifdef TOR3D
  !!         real*8              :: htor,tel(refElTor%Nnodes1d)
  !!#endif
  !!         real*8,allocatable :: qx(:,:),qy(:,:),auxq(:,:)
  !!         real*8              :: uex(Np,Neq),uey(Np,Neq)
  !!#ifdef TOR3D
  !!         real*8,allocatable :: qt(:,:)
  !!         real*8              :: uet(Np,Neq)
  !!#endif
  !!         ALLOCATE (u(Nel*Np,phys%Neq))
  !!         u = 0.
  !!         ALLOCATE (qx(Nel*Np,phys%Neq))
  !!         ALLOCATE (qy(Nel*Np,phys%Neq))
  !!         ALLOCATE (auxq(Nel*Np*phys%Neq,Ndim))
  !!         qx = 0.; qy = 0.
  !!#ifdef TOR3D
  !!         ALLOCATE (qt(Nel*Np,phys%Neq))
  !!         qt = 0.
  !!#endif
  !         neq = phys%Neq
  !         Nel = Mesh%Nelems
  !         Np = Mesh%Nnodesperelem
  !         Nfp = Mesh%Nnodesperface
  !         Nut = size(sol%u_tilde)
  !
  !         ALLOCATE (u(Nel*Np,Neq))
  !         ALLOCATE (qx(Nel*Np,Neq))
  !         ALLOCATE (qy(Nel*Np,Neq))
  !         ALLOCATE (auxq(Nel*Np*Neq,Ndim))
  !         ALLOCATE (u_tilde(Nut,neq))
  !
  !         u = transpose(reshape(sol%u,[Neq,Nel*Np]))
  !         u_tilde = transpose(reshape(sol%u_tilde,[Neq,Nut]))
  !         q = transpose(reshape(sol%u,[Neq*Ndim,Nel*Np]))
  !
  !
  !
  !#ifdef TOR3D
  ! write(6,*) "Reset variables error in 3D: Not coded yet"
  ! stop
  !#else
  !         !****************************************
  !         !          2D
  !         !****************************************
  !         DO iel = 1,Mesh%Nelems

  !            ind = (iel - 1)*Np + (/(i,i=1,Np)/)
  !            Xe = Mesh%X(Mesh%T(iel,:),:)
  !            CALL analytical_solution(Xe(:,1),Xe(:,2),ue)
  !            CALL analytical_gradient(Xe(:,1),Xe(:,2),ue,uex,uey)
  !            do ieq = 1,phys%neq
  !               if (switch%reset_eqs(ieq).ne.0) then
  !                  qx(ind,ieq) = uex
  !                  qy(ind,ieq) = uey
  !                  u(ind,ieq) = ue
  !
  !               endif
  !            end do
  !         END DO
  !
  !
  !         DO iFace = 1,Mesh%Nintfaces
  !            iel = Mesh%intfaces(iFace,1)
  !            ifa = Mesh%intfaces(iFace,2)
  !            faceNodes = refElPol%Face_nodes(ifa,:)
  !            Xf = Mesh%X(Mesh%T(iel,faceNodes),:)
  !            CALL analytical_solution(Xf(:,1),Xf(:,2),uf)
  !            ind_uf = (iFace - 1)*Nfp + (/(i,i=1,Nfp)/)
  !            do ieq = 1,phys%neq
  !               if (switch%reset_eqs(ieq).ne.0) then
  !                  u_tilde(ind_uf,ieq) = uf
  !               endif
  !            end do
  !         END DO

  !         DO iFace = 1,Mesh%Nextfaces
  !            iel = Mesh%extfaces(iFace,1)
  !            ifa = Mesh%extfaces(iFace,2)
  !            faceNodes = refElPol%Face_nodes(ifa,:)
  !            Xf = Mesh%X(Mesh%T(iel,faceNodes),:)
  !            CALL analytical_solution(Xf(:,1),Xf(:,2),uf)
  !            IF (.not. Mesh%Fdir(iel,ifa)) THEN
  !               ind_uf = Mesh%Nintfaces*Nfp + (iFace - 1)*Nfp + (/(i,i=1,Nfp)/)
  !               do ieq = 1,phys%neq
  !                  if (switch%reset_eqs(ieq).ne.0) then
  !                     u_tilde(ind_uf,:) = u(ind_ue(faceNodes),:)
  !                  endif
  !               end do
  !            END IF
  !         END DO
  !
  !#endif

  !!         !****************************************
  !!         !          common
  !!         !****************************************
  !!         sol%u = reshape(transpose(u),(/Nel*Np*phys%Neq/))

  !!         auxq(:,1) = reshape(transpose(qx),(/Nel*Np*phys%Neq/))
  !!         auxq(:,2) = reshape(transpose(qy),(/Nel*Np*phys%Neq/))
  !!#ifdef TOR3D
  !!         auxq(:,3) = reshape(transpose(qt),(/Nel*Np*phys%Neq/))
  !!#endif
  !!         sol%q = reshape(transpose(auxq),(/Nel*Np*phys%Neq*Ndim/))
  !!         DEALLOCATE (qx,qy,auxq)
  !!#ifdef TOR3D
  !!         DEALLOCATE (qt)
  !!#endif
  !!         DEALLOCATE (u)
  !
  !     end subroutine reset_variables

  !***************************************************************
  ! Extract face solution: routine to define a nodal face solution
  ! equal to the elemental solution at face nodes
  !***************************************************************
#ifdef TOR3D
  SUBROUTINE extractFaceSolution
    INTEGER :: neq,Ne,Nf,Nfe,unkF,Np,N2d,Np2d,Nfl,iel,iElem,ifa,iFace,i,itor,ntorloc,nut
    INTEGER :: c,Np1Dpol,Np1Dtor,Nfdir,sh
    INTEGER :: ind_ue(refElTor%Nnodes3D),ind2(refElPol%Nnodes2D)
    INTEGER :: ind3(refElPol%Nnodes2D*refElTor%Nnodes1D),indl(refElPol%Nnodes1D*refElTor%Nnodes1D)
    REAL*8,ALLOCATABLE :: u(:,:),u_tilde(:,:)

    sol%u_tilde = 0.d0
    neq = phys%Neq
    N2D = Mesh%Nelems                  ! Number of 2D elements
    Np2D = refElPol%Nnodes2D            ! Number of nodes for each 2D element
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
    Ne = N2D*ntorloc                  ! Number of 3D elements
    Nf = Mesh%Nfaces
    Np1Dpol = refElPol%Nnodes1D         ! Number of nodes in the 1D poloidal segments
    Np1Dtor = refElTor%Nnodes1D         ! Number of nodes in the 1D toroidal segments
    Nfl = Np1Dpol*Np1Dtor              ! Number of nodes in the lateral faces
    Nfe = refElPol%Nfaces
    unkF = Mesh%ukf
    Np = Np2D*refElTor%Nnodes1D       ! Number of nodes for each 3D element
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       nut = ntorloc*(Nfl*Nf + Np2d*N2d) + Np2d*N2d ! Size of utilde per equation
    ELSE
       nut = ntorloc*(Nfl*Nf + Np2d*N2d) ! Size of utilde per equation
    ENDIF
#else
    nut = ntorloc*(Nfl*Nf + Np2d*N2d) ! Size of utilde per equation
#endif
    !                                                nut  = ntorloc*(Nfl*Nf + Np2D*N2D)  ! Size of utilde per equation
    Nfdir = Mesh%Ndir

    ! Indices
    indl = (/(i,i=1,Nfl)/)
    ind2 = (/(i,i=1,Np2D)/)
    ind3 = (/(i,i=1,Np)/)

    ALLOCATE (u(Ne*Np,neq))
    ALLOCATE (u_tilde(nut,neq))
    u_tilde = 0.d0
    u = TRANSPOSE(RESHAPE(sol%u,(/neq,Ne*Np/)))

    ! Loop in elements
    c = 0
    DO itor = 1,ntorloc
       ! Poloidal faces
       DO iel = 1,N2D
          iElem = (itor - 1)*N2D+iel
          ind_ue = (iElem - 1)*Np + ind3
          u_tilde(c + ind2,:) = u(ind_ue(ind2),:)
          c = c + Np2D
       END DO

       ! Toroidal interior faces
       DO iFace = 1,Mesh%Nintfaces
          iel = Mesh%intfaces(iFace,1)
          ifa = Mesh%intfaces(iFace,2)
          iElem = (itor - 1)*N2D+iel
          ind_ue = (iElem - 1)*Np + ind3
          u_tilde(c + indl,:) = u(ind_ue(refElTor%faceNodes3(ifa,:)),:)
          c = c + Nfl
       END DO

       ! Toroidal exterior faces
       DO iFace = 1,Mesh%Nextfaces
          iel = Mesh%extfaces(iFace,1)
          ifa = Mesh%extfaces(iFace,2)
          IF (.NOT. Mesh%Fdir(iel,ifa)) THEN
             iElem = (itor - 1)*N2D+iel
             ind_ue = (iElem - 1)*Np + ind3
             u_tilde(c + indl,:) = u(ind_ue(refElTor%faceNodes3(ifa,:)),:)
             c = c + Nfl
          END IF
       END DO
    END DO

#ifdef PARALL
    ! Add solution on toroidal ghost faces
    IF (MPIvar%ntor .GT. 1) THEN
       sh = (Np1Dtor - 1)*Np2d
       DO iel = 1,N2D
          iElem = (ntorloc - 1)*N2D+iel
          ind_ue = (iElem - 1)*Np + ind3
          u_tilde(c + ind2,:) = u(ind_ue(ind2 + sh),:)
          c = c + Np2D
       END DO
    ENDIF
#endif
    sol%u_tilde = RESHAPE(TRANSPOSE(u_tilde),(/nut*neq/))

    DEALLOCATE (u,u_tilde)
  END SUBROUTINE extractFaceSolution
#else
  SUBROUTINE extractFaceSolution
    INTEGER :: neq,Ne,Nf,Nfe,unkF,Np,Nfp,iElem,ifa,iFace,i
    INTEGER :: ind_uf(1:Mesh%Nnodesperface),faceNodes(1:Mesh%Nnodesperface)
    INTEGER :: ind_ue(1:Mesh%Nnodesperelem)
    REAL*8,ALLOCATABLE :: u(:,:),u_tilde(:,:)

    sol%u_tilde = 0.d0
    neq = phys%Neq
    Ne = Mesh%Nelems
    Nf = Mesh%Nfaces
    Nfe = refElPol%Nfaces
    unkF = Mesh%ukf
    Np = Mesh%Nnodesperelem
    Nfp = Mesh%Nnodesperface

    ALLOCATE (u(1:Ne*Np,1:neq))
    ALLOCATE (u_tilde(1:Nf*Nfp,1:neq))
    u_tilde = 0.d0
    u = TRANSPOSE(RESHAPE(sol%u,(/neq,Ne*Np/)))

    DO iFace = 1,Mesh%Nintfaces
       iElem = Mesh%intfaces(iFace,1)
       ifa = Mesh%intfaces(iFace,2)
       ind_ue = (iElem - 1)*Np + (/(i,i=1,Np)/)
       ind_uf = (iFace - 1)*Nfp + (/(i,i=1,Nfp)/)
       faceNodes = refElPol%Face_nodes(ifa,:)
       u_tilde(ind_uf,:) = u(ind_ue(faceNodes),:)
    END DO

    DO iFace = 1,Mesh%Nextfaces
       iElem = Mesh%extfaces(iFace,1)
       ifa = Mesh%extfaces(iFace,2)
       IF (.NOT. Mesh%Fdir(iElem,ifa)) THEN
          ind_ue = (iElem - 1)*Np + (/(i,i=1,Np)/)
          ind_uf = Mesh%Nintfaces*Nfp + (iFace - 1)*Nfp + (/(i,i=1,Nfp)/)
          faceNodes = refElPol%Face_nodes(ifa,:)
          u_tilde(ind_uf,:) = u(ind_ue(faceNodes),:)
       END IF
    END DO
    sol%u_tilde = RESHAPE(TRANSPOSE(u_tilde),(/Nf*Nfp*neq/))

    DEALLOCATE (u,u_tilde)
  END SUBROUTINE extractFaceSolution
#endif

  SUBROUTINE add_perturbation()
    INTEGER             :: Np
    INTEGER             :: iel,i,imod,nmod,ieq,indl, iel2, iphi
    INTEGER,ALLOCATABLE :: ind(:)
    REAL*8              :: Xe(Mesh%Nnodesperelem,Mesh%Ndim),pertphi,perttheta
    REAL*8              :: phi,amp
    REAL*8,ALLOCATABLE  :: u(:,:)

#ifdef TOR3D
    INTEGER             :: itor, itorg, itheta, Np1d, Np2d, ntorloc
    REAL*8              :: tdiv(numer%ntor + 1), tel(refElTor%Nnodes1d)
    REAL*8              :: htor, theta
#endif

    ALLOCATE(u(SIZE(sol%u)/phys%neq,phys%neq))
    u = TRANSPOSE(RESHAPE(sol%u,[phys%neq,SIZE(sol%u)/phys%neq]))

    amp = 1e-3
    nmod = 10

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
    Np1d = refElTor%Nnodes1D             ! N. of nodes for each toroidal 1d element
    Np2d = refElPol%Nnodes2D             ! N. of nodes for each poloidal 2D element
    Np = Np2d*Np1d                     ! N. of nodes for each 3D element
#else
    Np = refElPol%Nnodes2D
#endif

    ALLOCATE(ind(Np))

#ifdef TOR3D
    htor = numer%tmax/numer%ntor
    tdiv = 0.
    DO i = 1,numer%ntor
       tdiv(i + 1) = i*htor
    END DO

    DO itor = 1,ntorloc
#ifdef PARALL
       itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
       IF (itorg == numer%ntor + 1) itorg = 1
#else
       itorg = itor
       tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
#endif
#endif
       DO iel2 = 1,Mesh%Nelems
          iel = iel2
#ifdef TOR3D
          iel = (itor - 1)*Mesh%Nelems+iel2
#endif
          ind = (iel - 1)*Np + (/(i,i=1,Np)/)
          Xe = Mesh%X(Mesh%T(iel2,:),:)*simpar%refval_length
          perttheta = 1.

          DO imod=1,Nmod
#ifdef TOR3D
             DO itheta =1,refElTor%Nnodes1d
                theta = tel(itheta)
                perttheta = 1+amp*(COS(imod*theta))
#endif
                DO iphi = 1,refElPol%Nnodes2D
                   phi = ATAN2(Xe(iphi,2),Xe(iphi,1)-geom%R0)
                   pertphi = (1+amp*(COS(imod*phi)))
                   indl = iphi
#ifdef TOR3D
                   indl = (itheta-1)*refElPol%Nnodes2D+iphi
#endif
                   DO ieq = 1,phys%neq
                      u(ind(indl),ieq) = u(ind(indl),ieq)*pertphi*perttheta
                   END DO ! ieq
                END DO ! iphi
#ifdef TOR3D
             END DO ! itheta
#endif
          END DO ! modes
       END DO ! elements 2d
#ifdef TOR3D
    END DO ! toroidal loop
#endif
    sol%u = col(TRANSPOSE(u))
    DEALLOCATE(ind,u)

  END SUBROUTINE add_perturbation

  SUBROUTINE add_blob()
    INTEGER             :: iel,i
    INTEGER             :: ind(refElPol%Nnodes2D)
    REAL*8              :: Xe(refElPol%Nnodes2D,Mesh%Ndim)
    REAL*8              :: xmax,xmin,ymax,ymin,xm,ym,smod,rs,xsource,ysource
    REAL*8              :: dsource(refElPol%Nnodes2D),aux(refElPol%Nnodes2D)
    REAL*8,ALLOCATABLE  :: u(:,:)

#ifdef TOR3D
    WRITE(6,*) "Blob perturbation not implemented yet"
    STOP
#endif

    ALLOCATE(u(SIZE(sol%u)/phys%neq,phys%neq))
    u = TRANSPOSE(RESHAPE(sol%u,[phys%neq,SIZE(sol%u)/phys%neq]))

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax+xmin)
    ym = 0.5*(ymax+ymin)

    DO iel = 1,Mesh%Nelems
       ind = (iel - 1)*refElPol%Nnodes2D + (/(i,i=1,refElPol%Nnodes2D)/)
       Xe = Mesh%X(Mesh%T(iel,:),:)
       smod = 0.1
       rs = 0.04/simpar%refval_length
       xsource = xm+0.85*(xmax-xm)
       ysource = ym
       dsource   = SQRT((Xe(:,1)-xsource)**2+(Xe(:,2)-ysource)**2)
       aux = -dsource**2/rs**2
       DO i=1,refElPol%Nnodes2D
          IF (aux(i).GT.-30) THEN
             u(ind(i),1) =  u(ind(i),1)+smod*EXP(aux(i))
          ENDIF
       END DO
    END DO ! elements 2d

    sol%u = col(TRANSPOSE(u))

    DEALLOCATE(u)
  END SUBROUTINE add_blob

  SUBROUTINE projectSolutionDifferentMeshes_general(T1, X1, T2, X2, u1, q1, u2, q2)

    INTEGER, INTENT(IN)                                    :: T1(:,:), T2(:,:)
    REAL*8, INTENT(IN)                                     :: X1(:,:), X2(:,:)
    REAL*8, POINTER, DIMENSION(:), INTENT(IN)              :: u1(:)
    REAL*8, POINTER, DIMENSION(:), OPTIONAL, INTENT(IN)    :: q1(:)
    REAL*8, POINTER, DIMENSION(:), INTENT(INOUT)           :: u2(:)
    REAL*8, POINTER, DIMENSION(:), OPTIONAL, INTENT(INOUT) :: q2(:)
    REAL*8, ALLOCATABLE, DIMENSION(:,:)                    :: u1_2D, u2_2D
    REAL*8, ALLOCATABLE, DIMENSION(:,:,:)                  :: q1_3D, q2_3D
    INTEGER                                                :: i,j,k, counter


    ALLOCATE(u1_2D(SIZE(T1,1)*SIZE(T1,2), phys%neq))
    ALLOCATE(u2_2D(SIZE(T2,1)*SIZE(T2,2), phys%neq))

    IF(PRESENT(q1)) THEN
       ALLOCATE(q1_3D(SIZE(T1,1)*SIZE(T1,2), phys%neq, refElPol%ndim))
       ALLOCATE(q2_3D(SIZE(T2,1)*SIZE(T2,2), phys%neq, refElPol%ndim))
       q1_3D = 0.
       q2_3D = 0.
    ENDIF

    u1_2D = 0.
    u2_2D = 0.

    counter = 1
    ! reshape u sol
    DO i = 1, SIZE(u1_2D,1)
       DO j = 1, SIZE(u1_2D,2)
          u1_2D(i,j) = u1(counter)
          counter = counter + 1
       ENDDO
    ENDDO

    IF(PRESENT(q1)) THEN
       counter = 1
       ! q1_3D not as easy
       DO i = 1, SIZE(q1_3D,1)
          DO j = 1, SIZE(q1_3D,2)
             DO k = 1, SIZE(q1_3D,3)
                q1_3D(i,j,k) = q1(counter)
                counter = counter + 1
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF(PRESENT(q1)) THEN
       CALL projectSolutionDifferentMeshes_Mod(T1, X1, T2, X2, u1_2D, q1_3D, u2_2D, q2_3D)
    ELSE
       CALL projectSolutionDifferentMeshes_Mod(T1, X1, T2, X2, u_old = u1_2D, u_new=u2_2D)
    ENDIF

    IF(SIZE(u2) .NE. SIZE(u2_2D)) THEN
       DEALLOCATE(u2)
       ALLOCATE(u2(SIZE(u2_2D)))
    ENDIF

    IF(PRESENT(q1)) THEN
       IF(SIZE(q2) .NE. SIZE(q2_3D)) THEN
          DEALLOCATE(q2)
          ALLOCATE(q2(SIZE(q2_3D)))
       ENDIF
    ENDIF
    ! solu is easy to reshape
    counter = 1
    DO i = 1, SIZE(u2_2D,1)
       DO j = 1, SIZE(u2_2D,2)
          u2(counter)   = u2_2D(i,j)
          counter = counter + 1
       ENDDO
    ENDDO

    IF(PRESENT(q1)) THEN
       ! solq is not as easy
       counter = 1
       DO i = 1, SIZE(q2_3D,1)
          DO j = 1, SIZE(q2_3D,2)
             DO k = 1, SIZE(q2_3D,3)
                q2(counter)   = q2_3D(i,j,k)
                counter = counter + 1
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    DEALLOCATE(u2_2D)
    DEALLOCATE(u1_2D)

    IF(PRESENT(q1)) THEN
       DEALLOCATE(q2_3D)
       DEALLOCATE(q1_3D)
    ENDIF

  ENDSUBROUTINE projectSolutionDifferentMeshes_general

  SUBROUTINE projectSolutionDifferentMeshes_general_arrays(T1, X1, T2, X2, u1, q1, u2, q2)

    INTEGER, INTENT(IN)                                :: T1(:,:), T2(:,:)
    REAL*8, INTENT(IN)                                 :: X1(:,:), X2(:,:)
    REAL*8, DIMENSION(:), INTENT(IN)                   :: u1(:)
    REAL*8, DIMENSION(:), OPTIONAL,INTENT(IN)          :: q1(:)
    REAL*8, DIMENSION(:), INTENT(INOUT)                :: u2
    REAL*8, DIMENSION(:), OPTIONAL,INTENT(INOUT)       :: q2
    REAL*8, ALLOCATABLE, DIMENSION(:,:)                :: u1_2D, u2_2D
    REAL*8, ALLOCATABLE, DIMENSION(:,:,:)              :: q1_3D, q2_3D
    INTEGER                                            :: i,j,k, counter

    ALLOCATE(u1_2D(SIZE(T1,1)*SIZE(T1,2), phys%neq))
    ALLOCATE(u2_2D(SIZE(T2,1)*SIZE(T2,2), phys%neq))


    IF(PRESENT(q1)) THEN
       ALLOCATE(q1_3D(SIZE(T1,1)*SIZE(T1,2), phys%neq, refElPol%ndim))
       ALLOCATE(q2_3D(SIZE(T2,1)*SIZE(T2,2), phys%neq, refElPol%ndim))
       q1_3D = 0.
       q2_3D = 0.
    ENDIF

    u1_2D = 0.
    u2_2D = 0.


    counter = 1
    ! reshape u sol
    DO i = 1, SIZE(u1_2D,1)
       DO j = 1, SIZE(u1_2D,2)
          u1_2D(i,j) = u1(counter)
          counter = counter + 1
       ENDDO
    ENDDO


    IF(PRESENT(q1)) THEN
       counter = 1
       DO i = 1, SIZE(q1_3D,1)
          DO j = 1, SIZE(q1_3D,2)
             DO k = 1, SIZE(q1_3D,3)
                q1_3D(i,j,k) = q1(counter)
                counter = counter + 1
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF(PRESENT(q1)) THEN
       CALL projectSolutionDifferentMeshes_mod(T1, X1, T2, X2, u1_2D, q1_3D, u2_2D, q2_3D)
    ELSE
       CALL projectSolutionDifferentMeshes_mod(T1, X1, T2, X2, u_old=u1_2D, u_new = u2_2D)
    ENDIF

    ! IF(SIZE(u) .ne. SIZE(u2_2D)) THEN
    !   DEALLOCATE(u)
    !   ALLOCATE(u(SIZE(u2_2D)))
    ! ENDIF
    !
    ! IF(SIZE(q) .ne. SIZE(q2_3D)) THEN
    !   DEALLOCATE(q)
    !   ALLOCATE(q(SIZE(q2_3D)))
    ! ENDIF

    ! solu is easy to reshape
    counter = 1
    DO i = 1, SIZE(u2_2D,1)
       DO j = 1, SIZE(u2_2D,2)
          u2(counter)   = u2_2D(i,j)
          counter = counter + 1
       ENDDO
    ENDDO

    IF(PRESENT(q1)) THEN
       ! solq is not as easy
       counter = 1
       DO i = 1, SIZE(q2_3D,1)
          DO j = 1, SIZE(q2_3D,2)
             DO k = 1, SIZE(q2_3D,3)
                q2(counter)   = q2_3D(i,j,k)
                counter = counter + 1
             ENDDO
          ENDDO
       ENDDO
       DEALLOCATE(q2_3D)
       DEALLOCATE(q1_3D)
    ENDIF

    DEALLOCATE(u2_2D)
    DEALLOCATE(u1_2D)


  ENDSUBROUTINE projectSolutionDifferentMeshes_general_arrays

  SUBROUTINE projectSolutionDifferentMeshes_Mod(T1, X1, T2, X2, u_old, q_old, u_new, q_new)
    USE linearAlgebra, ONLY: solve_linear_system_sing, colint, invert_matrix
    USE adaptivity_common_module, ONLY: unique_1D, find_matches_int, inverse_isop_transf
    USE reference_element, ONLY: compute_shape_functions_at_points

    INTEGER, INTENT(IN)         :: T1(:,:), T2(:,:)
    REAL*8, INTENT(IN)          :: X1(:,:), X2(:,:)
    REAL*8                      :: xs(SIZE(T2,1)*SIZE(T2,2), 2)

    REAL*8, INTENT(IN)          :: u_old(:,:)
    REAL*8, OPTIONAL,INTENT(IN) :: q_old(:,:,:)

    REAL*8, INTENT(OUT)         :: u_new(:,:)
    REAL*8, OPTIONAL,INTENT(OUT):: q_new(:,:,:)

    REAL*8                      :: X_old(SIZE(X1,1), SIZE(X1,2)), Xe_elem(Mesh%Nnodesperelem, refElPol%Ndim)
    REAL*8                      :: A(3,3), b(3, SIZE(xs,1)), bcc(3)
    REAL*8                      :: tol, detA, a11,a12,a13,a21,a22,a23,a31,a32,a33, b1, b2, b3
    INTEGER                     :: T_old(SIZE(T1,1), SIZE(T1,2)), ind(SIZE(T1,2)), correl(SIZE(xs,1)), indcheck(SIZE(xs,1))
    INTEGER                     :: i, j, counter, n_elements, np_perelem, iel
    REAL*8,  ALLOCATABLE        :: shapeFunctions(:,:,:)
    REAL*8,  ALLOCATABLE        :: x(:,:), xieta(:,:)
    REAL*8, ALLOCATABLE         :: u_old_ind(:,:), q_old_ind(:,:,:)
    REAL*8, ALLOCATABLE         :: u_prov(:,:), q_prov(:,:,:)
    INTEGER, ALLOCATABLE        :: correl_unique(:), indices(:)
    INTEGER                     :: elem_full(SIZE(T1,1))

    u_new = 0
    IF(PRESENT(q_new)) THEN
       q_new = 0
    ENDIF

    ALLOCATE(u_old_ind(SIZE(T1,2), SIZE(u_old,2)))
    u_old_ind = 0
    ALLOCATE(u_prov(SIZE(xs,1), SIZE(u_old,2)))
    u_prov = 0

    IF(PRESENT(q_old)) THEN
       ALLOCATE(q_prov(SIZE(xs,1), SIZE(q_old,2), SIZE(q_old,3)))
       ALLOCATE(q_old_ind(SIZE(T1,2), SIZE(q_old,2),2))
       q_prov = 0
       q_old_ind = 0
    ENDIF


    n_elements = SIZE(T_old,1)
    np_perelem = SIZE(T_old,2)
    xs = X2(colint(TRANSPOSE(T2)),:)
    X_old = X1
    T_old = T1
    elem_full = 0
    u_prov = 0


    correl = 0
    indcheck = 0

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,*) "*** Projecting the solution: find corresponding elements"
    ENDIF

    !$OMP parallel private(iel, counter, i, tol, Xe_elem, A, b, bcc, detA, a11,a12,a13,a21,a22,a23,a31,a32,a33, b1, b2, b3) shared( xs, X_old, T_old, correl)
    tol = 1e-11
    A(3,:) = 1.
    b3 = 1.

    DO
       !$OMP DO SCHEDULE(STATIC)
       DO iel = 1, n_elements
          Xe_elem = X_old(T_old(iel,1:3),:)
          A(1,:) = Xe_elem(1:3,1)
          A(2,:) = Xe_elem(1:3,2)

          detA = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) + &
               &A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

          a11 =   A(2,2)*A(3,3) - A(2,3)*A(3,2)
          a12 = - A(1,2)*A(3,3) + A(1,3)*A(3,2)
          a13 =   A(1,2)*A(2,3) - A(1,3)*A(2,2)
          a21 = - A(2,1)*A(3,3) + A(2,3)*A(3,1)
          a22 =   A(1,1)*A(3,3) - A(1,3)*A(3,1)
          a23 = - A(1,1)*A(2,3) + A(1,3)*A(2,1)
          a31 =   A(2,1)*A(3,2) - A(2,2)*A(3,1)
          a32 = - A(1,1)*A(3,2) + A(1,2)*A(3,1)
          a33 =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

          DO i=1,SIZE(xs,1)
             IF(correl(i) .NE. 0) CYCLE
             b1 = xs(i,1)
             b2 = xs(i,2)

             bcc(1) = (a11*b1+a12*b2+a13*b3)/detA
             bcc(2) = (a21*b1+a22*b2+a23*b3)/detA
             bcc(3) = (a31*b1+a32*b2+a33*b3)/detA

             IF ( bcc(1)>=-tol .AND. bcc(2)>=-tol .AND. bcc(3)>=-tol .AND. bcc(1)<=1+tol .AND. bcc(2)<=1+tol .AND. bcc(3)<=1+tol) THEN
                correl(i) = iel
             ENDIF
          ENDDO
       ENDDO
       !$OMP END DO

       !$OMP BARRIER
       IF(ALL(correl .NE. 0)) THEN
          EXIT
       ENDIF

       tol = tol * 10

    ENDDO
    !$OMP END PARALLEL


    ! IF(ANY(correl .eq. 0)) THEN
    !   WRITE(*,*) "Couldn't find a point in projection. STOP."
    !
    !
    !   DO i = 1, SIZE(correl)
    !     IF(correl(i) .eq. 0) THEN
    !       WRITE(*,*) i, MOD(i, Mesh%Nnodesperelem)
    !     ENDIF
    !   ENDDO
    !
    !   STOP
    ! ENDIF

    WRITE(*,*) "TWO"

    !$OMP parallel private(iel, indices, xieta, shapeFunctions, Xe_elem, x, ind, u_old_ind, q_old_ind) shared(xs, X_old, T_old, u_new, q_new, correl, n_elements, refElPol, np_perelem, indcheck)
    !$OMP DO SCHEDULE(STATIC)
    DO iel = 1, n_elements
       IF(iel .EQ. 0) CYCLE

       CALL find_matches_int(correl, iel, indices)

       ALLOCATE(xieta(SIZE(indices), SIZE(xs,2)))
       ALLOCATE(x(SIZE(indices), SIZE(xs,2)))
       ALLOCATE(shapeFunctions(np_perelem, SIZE(indices), 3))

       xieta = 0
       x = 0
       shapeFunctions = 0

       !Xe_elem = X_old(T_old(iel,1:3),:)
       Xe_elem = X_old(T_old(iel,:),:)
       x = xs(indices,:)

       ! Find the corresponding point in the reference element (linear approach so far)

       ! a11 =   (Xe_elem(3,2)-Xe_elem(1,2))
       ! a12 =   0.5*(Xe_elem(2,1)+Xe_elem(3,1))
       ! a13 =   (Xe_elem(3,1)-Xe_elem(1,1))
       ! a21 =   0.5*(Xe_elem(2,2)+Xe_elem(3,2))
       ! a22 =   (Xe_elem(2,1)-Xe_elem(1,1))
       ! a31 =   (Xe_elem(2,2)-Xe_elem(1,2))
       ! a32 =   0.5*(Xe_elem(2,1)+Xe_elem(3,1))
       !
       ! d = 0.5*(a22*a11-a13*a31)
       ! d = 1./d
       !
       ! xieta(:,1)= d*(a11*(x(:,1)-a12) - a13*(x(:,2)-a21))
       ! xieta(:,2)= d*(a22*(x(:,2)-a21) - a31*(x(:,1)-a12))

       ! this goddamn function always gives problems
       CALL inverse_isop_transf(x, Xe_elem, refElPol, xieta)

       ! just fucking brute force it
       DO j = 1, SIZE(xieta,2)
          DO i = 1, SIZE(xieta,1)
             IF(ABS(xieta(i,j)-1.0) .LT. 1e-12) THEN
                xieta(i,j) = xieta(i,j) - 1.e-10
             ENDIF
          ENDDO
       ENDDO

       CALL compute_shape_functions_at_points(refElPol, xieta, shapeFunctions)

       ind = (iel-1)*np_perelem + (/ (j, j=1, np_perelem) /)

       u_old_ind = u_old(ind, :)
       u_new(indices,:) = MATMUL(TRANSPOSE(shapeFunctions(:,:,1)), u_old_ind)

       IF(PRESENT(q_old)) THEN
          q_old_ind = q_old(ind, :, :)
          q_new(indices,:,1) = MATMUL(TRANSPOSE(shapeFunctions(:,:,1)), q_old_ind(:,:,1))
          q_new(indices,:,2) = MATMUL(TRANSPOSE(shapeFunctions(:,:,1)), q_old_ind(:,:,2))
       ENDIF

       DEALLOCATE(xieta)
       DEALLOCATE(x)
       DEALLOCATE(indices)
       DEALLOCATE(shapeFunctions)
    END DO
    !$OMP END DO
    !$OMP end parallel

    WRITE(*,*) "THREE"


    ! DO WHILE (ANY(indcheck .eq. 0))
    !   tol = tol * 10
    !
    !   CALL find_matches_int(indcheck, 0, unknownnodes)
    !
    !   !$OMP parallel private(iel, i, Xe_elem, A, b, bcc, detA, a11,a12,a13,a21,a22,a23,a31,a32,a33, b1, b2, b3, d) shared(xs, X_old, T_old, correl, tol)
    !    A(3,:) = 1.
    !    b3 = 1.
    !   !$OMP DO SCHEDULE(STATIC)
    !    DO iel = 1, n_elements
    !        Xe_elem = X_old(T_old(iel,1:3),:)
    !        A(1,:) = Xe_elem(1:3,1)
    !        A(2,:) = Xe_elem(1:3,2)
    !
    !        detA = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) + &
    !                  &A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
    !
    !        a11 =   A(2,2)*A(3,3) - A(2,3)*A(3,2)
    !        a12 = - A(1,2)*A(3,3) + A(1,3)*A(3,2)
    !        a13 =   A(1,2)*A(2,3) - A(1,3)*A(2,2)
    !        a21 = - A(2,1)*A(3,3) + A(2,3)*A(3,1)
    !        a22 =   A(1,1)*A(3,3) - A(1,3)*A(3,1)
    !        a23 = - A(1,1)*A(2,3) + A(1,3)*A(2,1)
    !        a31 =   A(2,1)*A(3,2) - A(2,2)*A(3,1)
    !        a32 = - A(1,1)*A(3,2) + A(1,2)*A(3,1)
    !        a33 =   A(1,1)*A(2,2) - A(1,2)*A(2,1)
    !
    !        DO i=1,SIZE(xs,1)
    !
    !              b1 = xs(i,1)
    !              b2 = xs(i,2)
    !
    !              bcc(1) = (a11*b1+a12*b2+a13*b3)/detA
    !              bcc(2) = (a21*b1+a22*b2+a23*b3)/detA
    !              bcc(3) = (a31*b1+a32*b2+a33*b3)/detA
    !
    !              IF ( bcc(1)>=-tol .and. bcc(2)>=-tol .and. bcc(3)>=-tol .and. bcc(1)<=1+tol .and. bcc(2)<=1+tol .and. bcc(3)<=1+tol) then
    !                correl(i) = iel
    !              ENDIF
    !        ENDDO
    !    ENDDO
    !    !$OMP END DO
    !    !$OMP END PARALLEL
    !
    !   CALL unique_1D(correl, correl_unique)
    !
    !   WRITE(*,*) "FOUR"
    !
    !   !!$OMP parallel do schedule(static) private(i, iel, indices, xieta, x, shapeFunctions, Xe_elem, ind, u_old_ind, q_old_ind) shared(indcheck, correl, correl_unique, xs, X_old, T_old, refElPol, u_old, q_old, u_prov, q_prov)
    !   DO i = 1, SIZE(correl_unique)
    !     iel = correl_unique(i)
    !
    !     IF(iel .eq. 0) CYCLE
    !
    !     CALL find_matches_int(correl, iel, indices)
    !     indcheck(indices) = 1
    !
    !     ALLOCATE(xieta(size(indices), size(xs,2)))
    !     ALLOCATE(x(size(indices), size(xs,2)))
    !     ALLOCATE(shapeFunctions(np_perelem, size(indices), 3))
    !
    !     Xe_elem = X_old(T_old(iel,:),:)
    !     x = xs(indices,:)
    !
    !     ! Find the corresponding point in the reference element (linear approach so far)
    !     a11 =   (Xe_elem(3,2)-Xe_elem(1,2))
    !     a12 =   0.5*(Xe_elem(2,1)+Xe_elem(3,1))
    !     a13 =   (Xe_elem(3,1)-Xe_elem(1,1))
    !     a21 =   0.5*(Xe_elem(2,2)+Xe_elem(3,2))
    !     a22 =   (Xe_elem(2,1)-Xe_elem(1,1))
    !     a31 =   (Xe_elem(2,2)-Xe_elem(1,2))
    !     a32 =   0.5*(Xe_elem(2,1)+Xe_elem(3,1))
    !
    !     d = 0.5*(a22*a11-a13*a31)
    !     d = 1./d
    !
    !     xieta(:,1)= d*(a11*(x(:,1)-a12) - a13*(x(:,2)-a21))
    !     xieta(:,2)= d*(a22*(x(:,2)-a21) - a31*(x(:,1)-a12))
    !
    !     ! this goddamn function always gives problems
    !     !CALL inverse_isop_transf(x, Xe_elem, refElPol, xieta)
    !
    !     ! call HDF5_create('fortran_save_01.h5', file_id, ierr)
    !     ! call HDF5_array2D_saving_int(file_id, T_old, size(T_old,1), size(T_old,2), 'T_old')
    !     ! call HDF5_array2D_saving(file_id, X_old, size(X_old,1), size(X_old,2), 'X_old')
    !     ! call HDF5_array2D_saving(file_id, Xe_elem, size(Xe_elem,1), size(Xe_elem,2), 'Xe_f')
    !     ! call HDF5_array2D_saving(file_id, x, size(x,1),size(x,2), 'x_f')
    !     ! call HDF5_integer_saving(file_id, iel, 'iel_f')
    !     ! call HDF5_close(file_id)
    !
    !     ! just fucking brute force it
    !     DO j = 1, SIZE(xieta,1)
    !       IF(ANY(abs(xieta(j,:)-1.0) .lt. 1e-12)) THEN
    !         xieta(j,:) = xieta(j,:) - 1.e-11
    !       ENDIF
    !     ENDDO
    !
    !     CALL compute_shape_functions_at_points(RefElPol, xieta, shapeFunctions)
    !
    !     ind = (iel-1)*np_perelem + (/ (j, j=1, np_perelem) /)
    !
    !     u_old_ind = u_old(ind, :)
    !     u_prov(indices,:) = MATMUL(TRANSPOSE(shapeFunctions(:,:,1)), u_old_ind)
    !
    !     IF(PRESENT(q_old)) THEN
    !       q_old_ind = q_old(ind, :, :)
    !       q_prov(indices,:,1) = MATMUL(TRANSPOSE(shapeFunctions(:,:,1)), q_old_ind(:,:,1))
    !       q_prov(indices,:,2) = MATMUL(TRANSPOSE(shapeFunctions(:,:,1)), q_old_ind(:,:,2))
    !     ENDIF
    !
    !     DEALLOCATE(xieta)
    !     DEALLOCATE(x)
    !     DEALLOCATE(indices)
    !     DEALLOCATE(shapeFunctions)
    !   END DO
    !   !!$OMP end parallel do
    !
    !   u_new(unknownnodes,:) = u_prov(unknownnodes,:)
    !   IF(PRESENT(q_old)) THEN
    !     q_new(unknownnodes,:,:) = q_prov(unknownnodes,:,:)
    !   ENDIF
    !   DEALLOCATE(unknownnodes)
    ! ENDDO

    ! IF(ANY(indcheck .eq. 0)) THEN
    !   WRITE(*,*) "indcheck equal to 0. STOPPING."
    !   STOP
    ! ENDIF

    IF(ALLOCATED(correl_unique)) DEALLOCATE(correl_unique)

    DEALLOCATE(u_old_ind)
    DEALLOCATE(u_prov)
    IF(PRESENT(q_old)) THEN
       DEALLOCATE(q_old_ind)
       DEALLOCATE(q_prov)
    ENDIF

  END SUBROUTINE projectSolutionDifferentMeshes_Mod


#ifdef WITH_PETSC
  SUBROUTINE InitPETSC
#include "petsc/finclude/petsc.h"
    USE petsc, ONLY: PetscInitialize
    !use petscsys
    USE MPI_OMP
    IMPLICIT NONE

    PetscErrorCode :: ierr
    PetscBool      :: initialized

    ! Init PETSC
    CALL PetscInitialized(initialized, ierr)
    IF (.NOT. initialized) THEN
       CALL petscinitializenoarguments(ierr) ! Default communicator = MPI_COMM_WORLD
    ENDIF

    IF (ierr .NE. 0) THEN
       PRINT*,'Unable to initialize PETSc even though it was requested. Aborting...'
       CALL MPI_Abort(MPI_COMM_WORLD,-1,ierr)
    ENDIF

  END SUBROUTINE InitPETSC

  SUBROUTINE FinalizePETSC
#include "petsc/finclude/petsc.h"
    USE petsc, ONLY: PetscFinalize
    USE MPI_OMP
    IMPLICIT NONE

    PetscErrorCode :: ierr
    PetscBool      :: finalized

    ! Init PETSC

    CALL PetscFinalized(finalized, ierr)
    IF (.NOT. finalized) THEN
       CALL PetscFinalize(ierr) ! Default communicator = MPI_COMM_WORLD
    ENDIF
    IF (ierr .NE. 0) THEN
       PRINT*,'Unable to finalize PETSc even though it was requested. Aborting...'
       CALL MPI_Abort(MPI_COMM_WORLD,-1,ierr)
    ENDIF

  END SUBROUTINE FinalizePETSC
#endif
END MODULE initialization
