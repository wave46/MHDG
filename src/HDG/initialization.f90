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
  USE analytical, only: analytical_solution, analytical_gradient
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

    IF (utils%printint > 0) THEN
       IF(MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*            INITIALIZING SOLUTION              *'
          WRITE (6, *) '*************************************************'
       ENDIF
    END IF


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
    ALLOCATE(sol%u(sizeu))
    ALLOCATE(sol%u_tilde(sizeutilde))
    ALLOCATE(sol%u_tilde0(sizeutilde))
    ALLOCATE(sol%q(sizeu*Ndim))
    sol%u = 0.
    sol%u_tilde = 0.
    sol%q = 0.
    sol%u_tilde0 = 0.

    IF (switch%init.EQ.1) THEN
       ! The solution is intialized in each node to the analytical solution
       IF (MPIvar%glob_id .EQ. 0) THEN
          IF (utils%printint > 0) THEN
             WRITE (6, *) '*      Initializing analytical solution         *'
          END IF
       ENDIF
       CALL init_sol_analytic()
    ELSEIF (switch%init.EQ.2) THEN
       ! The solution is intialized in each node to the analytical solution
       IF (MPIvar%glob_id .EQ. 0) THEN
          IF (utils%printint > 0) THEN
             WRITE (6, *) '* Initializing the solution with L2 projection  *'
          END IF
       ENDIF
       CALL init_sol_l2proj()
    ELSE
       WRITE(6,*) "Wrong initialization type. STOP."
       STOP
    ENDIF
    ! Extract the face solution from the elemental one
    CALL extractFaceSolution()

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
      auxq = 0.
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
      qx = 0.; qy = 0.;auxq = 0.
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
    u = 0.d0
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
    u = 0.d0
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
          IF (Mesh%flipFace(iElem,ifa)) THEN
             ind_uf = Mesh%Nintfaces*Nfp + (iFace - 1)*Nfp + (/(i,i=Nfp,1,-1)/)
          ELSE
             ind_uf = Mesh%Nintfaces*Nfp + (iFace - 1)*Nfp + (/(i,i=1,Nfp)/)
          ENDIF
          faceNodes = refElPol%Face_nodes(ifa,:)
          u_tilde(ind_uf,:) = u(ind_ue(faceNodes),:)
       END IF
    END DO
    sol%u_tilde = RESHAPE(TRANSPOSE(u_tilde),(/Nf*Nfp*neq/))

    DEALLOCATE (u,u_tilde)
  END SUBROUTINE extractFaceSolution
#endif

  SUBROUTINE add_initial_perturbation()
    IF (switch%pertini .EQ. 1) THEN
       CALL add_perturbation()
       WRITE(6,*) "Adding perturbation to the initial solution"
    ELSE IF (switch%pertini .EQ. 2) THEN
       CALL add_blob()
       WRITE(6,*) "Adding density blob to initial solution"
    ENDIF
  ENDSUBROUTINE add_initial_perturbation

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
    ind = 0

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
       smod = 0.2
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

  PURE SUBROUTINE curved_element_bounding_box(element_coordinates, xmin, xmax, ymin, ymax, element_size)
    REAL*8, INTENT(IN)          :: element_coordinates(:,:)
    REAL*8, INTENT(OUT)         :: xmin, xmax, ymin, ymax, element_size

    xmin = MINVAL(element_coordinates(:,1))
    xmax = MAXVAL(element_coordinates(:,1))
    ymin = MINVAL(element_coordinates(:,2))
    ymax = MAXVAL(element_coordinates(:,2))
    element_size = MAX(xmax-xmin, ymax-ymin)
  ENDSUBROUTINE curved_element_bounding_box

  PURE LOGICAL FUNCTION point_in_padded_bounding_box(point, xmin, xmax, ymin, ymax, padding)
    REAL*8, INTENT(IN)          :: point(1,2)
    REAL*8, INTENT(IN)          :: xmin, xmax, ymin, ymax, padding

    point_in_padded_bounding_box = point(1,1) .GE. xmin-padding .AND. point(1,1) .LE. xmax+padding .AND. &
         point(1,2) .GE. ymin-padding .AND. point(1,2) .LE. ymax+padding
  ENDFUNCTION point_in_padded_bounding_box

  SUBROUTINE find_nearest_curved_element(target_point, old_connectivity, old_coordinates, candidate_box_padding_ratio, &
       nearest_element, nearest_valid_point, nearest_distance, nearest_element_size)
    USE adaptivity_common_module, ONLY: inverse_isop_transf, clamp_to_curved_triangle

    REAL*8, INTENT(IN)          :: target_point(1,2)
    INTEGER, INTENT(IN)         :: old_connectivity(:,:)
    REAL*8, INTENT(IN)          :: old_coordinates(:,:)
    REAL*8, INTENT(IN)          :: candidate_box_padding_ratio
    INTEGER, INTENT(OUT)        :: nearest_element
    REAL*8, INTENT(OUT)         :: nearest_valid_point(1,2), nearest_distance, nearest_element_size

    REAL*8                      :: element_coordinates(Mesh%Nnodesperelem, refElPol%Ndim)
    REAL*8                      :: reference_point(1,2), clamped_point(1,2)
    REAL*8                      :: xmin, xmax, ymin, ymax, element_size, bbox_pad, distance
    INTEGER                     :: element
    LOGICAL                     :: inverse_converged, clamp_valid

    nearest_element = 0
    nearest_valid_point = target_point
    nearest_distance = HUGE(1.d0)
    nearest_element_size = 0.d0

    DO element = 1, SIZE(old_connectivity,1)
       element_coordinates = old_coordinates(old_connectivity(element,:),:)
       CALL curved_element_bounding_box(element_coordinates, xmin, xmax, ymin, ymax, element_size)
       bbox_pad = MAX(1.d-10, candidate_box_padding_ratio*element_size)

       IF(.NOT. point_in_padded_bounding_box(target_point, xmin, xmax, ymin, ymax, bbox_pad)) CYCLE

       CALL inverse_isop_transf(target_point, element_coordinates, refElPol, reference_point, inverse_converged)
       IF(.NOT. inverse_converged) CYCLE

       CALL clamp_to_curved_triangle(reference_point, element_coordinates, refElPol, clamped_point, clamp_valid)
       IF(.NOT. clamp_valid) CYCLE

       distance = SQRT((target_point(1,1)-clamped_point(1,1))**2 + (target_point(1,2)-clamped_point(1,2))**2)
       IF(distance .LT. nearest_distance) THEN
          nearest_element = element
          nearest_valid_point = clamped_point
          nearest_distance = distance
          nearest_element_size = element_size
       ENDIF
    ENDDO
  ENDSUBROUTINE find_nearest_curved_element

  PURE SUBROUTINE invert_2x2_matrix(matrix, inverse_matrix)
    REAL*8, INTENT(IN)          :: matrix(2,2)
    REAL*8, INTENT(OUT)         :: inverse_matrix(2,2)
    REAL*8                      :: determinant

    determinant = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)

    inverse_matrix(1,1) =  matrix(2,2)/determinant
    inverse_matrix(1,2) = -matrix(1,2)/determinant
    inverse_matrix(2,1) = -matrix(2,1)/determinant
    inverse_matrix(2,2) =  matrix(1,1)/determinant
  ENDSUBROUTINE invert_2x2_matrix

  SUBROUTINE find_points_in_linear_elements(target_points, old_connectivity, old_coordinates, point_elements)
    REAL*8, INTENT(IN)          :: target_points(:,:), old_coordinates(:,:)
    INTEGER, INTENT(IN)         :: old_connectivity(:,:)
    INTEGER, INTENT(INOUT)      :: point_elements(:)

    REAL*8                      :: triangle_vertices(3, refElPol%Ndim)
    REAL*8                      :: edge_matrix(2,2), inverse_edge_matrix(2,2)
    REAL*8                      :: barycentric_weights(3), point_offset(2)
    REAL*8                      :: barycentric_tolerance
    INTEGER                     :: point, element

!!$OMP parallel private(element, point, barycentric_tolerance, triangle_vertices, edge_matrix, inverse_edge_matrix, barycentric_weights, point_offset) shared(target_points, old_coordinates, old_connectivity, point_elements)
    barycentric_tolerance = 1.d-10

!!$OMP DO SCHEDULE(STATIC)
    DO element = 1, SIZE(old_connectivity,1)
       triangle_vertices = old_coordinates(old_connectivity(element,1:3),:)
       edge_matrix(:,1) = triangle_vertices(2,:) - triangle_vertices(1,:)
       edge_matrix(:,2) = triangle_vertices(3,:) - triangle_vertices(1,:)
       CALL invert_2x2_matrix(edge_matrix, inverse_edge_matrix)

       DO point = 1, SIZE(target_points,1)
          IF(point_elements(point) .NE. 0) CYCLE
          point_offset = target_points(point,:) - triangle_vertices(1,:)
          barycentric_weights(2:3) = MATMUL(inverse_edge_matrix, point_offset)
          barycentric_weights(1) = 1.d0 - SUM(barycentric_weights(2:3))

          IF(barycentric_weights(1) .GE. -barycentric_tolerance .AND. &
               barycentric_weights(2) .GE. -barycentric_tolerance .AND. &
               barycentric_weights(3) .GE. -barycentric_tolerance .AND. &
               barycentric_weights(1) .LE. 1.d0+barycentric_tolerance .AND. &
               barycentric_weights(2) .LE. 1.d0+barycentric_tolerance .AND. &
               barycentric_weights(3) .LE. 1.d0+barycentric_tolerance) THEN
             point_elements(point) = element
          ENDIF
       ENDDO
    ENDDO
!!$OMP END DO
!!$OMP END PARALLEL
  ENDSUBROUTINE find_points_in_linear_elements

  SUBROUTINE find_points_in_curved_elements(target_points, old_connectivity, old_coordinates, point_elements)
    USE adaptivity_common_module, ONLY: inverse_isop_transf

    REAL*8, INTENT(IN)          :: target_points(:,:), old_coordinates(:,:)
    INTEGER, INTENT(IN)         :: old_connectivity(:,:)
    INTEGER, INTENT(INOUT)      :: point_elements(:)

    REAL*8                      :: element_coordinates(Mesh%Nnodesperelem, refElPol%Ndim)
    REAL*8                      :: target_point(1,2), reference_point(1,2)
    REAL*8                      :: curved_tolerance, bounding_box_padding
    REAL*8                      :: xmin, xmax, ymin, ymax, element_size
    INTEGER                     :: point, element
    LOGICAL                     :: inverse_converged

    curved_tolerance = 1.d-8

    DO point = 1, SIZE(target_points,1)
       IF(point_elements(point) .NE. 0) CYCLE
       target_point(1,:) = target_points(point,:)

       DO element = 1, SIZE(old_connectivity,1)
          element_coordinates = old_coordinates(old_connectivity(element,:),:)
          CALL curved_element_bounding_box(element_coordinates, xmin, xmax, ymin, ymax, element_size)
          bounding_box_padding = MAX(1.d-12, curved_tolerance*MAX(1.d0, element_size))

          IF(.NOT. point_in_padded_bounding_box(target_point, xmin, xmax, ymin, ymax, bounding_box_padding)) CYCLE

          CALL inverse_isop_transf(target_point, element_coordinates, refElPol, reference_point, inverse_converged)
          IF(.NOT. inverse_converged) CYCLE

          IF(reference_point(1,1) .GE. -1.d0-curved_tolerance .AND. &
               reference_point(1,2) .GE. -1.d0-curved_tolerance .AND. &
               reference_point(1,1) .LE.  1.d0+curved_tolerance .AND. &
               reference_point(1,2) .LE.  1.d0+curved_tolerance .AND. &
               reference_point(1,1)+reference_point(1,2) .LE. curved_tolerance) THEN
             point_elements(point) = element
             EXIT
          ENDIF
       ENDDO
    ENDDO
  ENDSUBROUTINE find_points_in_curved_elements

  SUBROUTINE projectSolutionDifferentMeshes_Mod(old_connectivity, old_coordinates, new_connectivity, new_coordinates, &
       u_old, q_old, u_new, q_new)
    USE linearAlgebra, ONLY: colint
    USE adaptivity_common_module, ONLY: find_matches_int, inverse_isop_transf
    USE reference_element, ONLY: compute_shape_functions_at_points

    INTEGER, INTENT(IN)         :: old_connectivity(:,:), new_connectivity(:,:)
    REAL*8, INTENT(IN)          :: old_coordinates(:,:), new_coordinates(:,:)
    REAL*8                      :: target_points(SIZE(new_connectivity,1)*SIZE(new_connectivity,2), 2)
    REAL*8                      :: interpolation_points(SIZE(new_connectivity,1)*SIZE(new_connectivity,2), 2)

    REAL*8, INTENT(IN)          :: u_old(:,:)
    REAL*8, OPTIONAL,INTENT(IN) :: q_old(:,:,:)

    REAL*8, INTENT(OUT)         :: u_new(:,:)
    REAL*8, OPTIONAL,INTENT(OUT):: q_new(:,:,:)

    REAL*8                      :: element_coordinates(Mesh%Nnodesperelem, refElPol%Ndim)
    REAL*8                      :: nearest_tol, nearest_distance, nearest_element_size
    REAL*8                      :: target_point(1,2), nearest_valid_point(1,2)
    INTEGER                     :: old_element_dofs(SIZE(old_connectivity,2))
    INTEGER                     :: point_elements(SIZE(target_points,1))
    INTEGER                     :: nearest_candidate_elements(SIZE(target_points,1))
    REAL*8                      :: nearest_candidate_distances(SIZE(target_points,1))
    REAL*8                      :: nearest_candidate_sizes(SIZE(target_points,1))
    INTEGER                     :: i, j, n_elements, np_perelem, iel
    INTEGER                     :: nearest_element
    INTEGER                     :: local_missing, global_missing, ierr
    REAL*8,  ALLOCATABLE        :: shape_functions(:,:,:)
    REAL*8,  ALLOCATABLE        :: element_points(:,:), reference_points(:,:)
    REAL*8, ALLOCATABLE         :: old_element_u(:,:), old_element_q(:,:,:)
    INTEGER, ALLOCATABLE        :: point_indices(:)
    u_new = 0
    IF(PRESENT(q_new)) THEN
       q_new = 0
    ENDIF

    ALLOCATE(old_element_u(SIZE(old_connectivity,2), SIZE(u_old,2)))
    old_element_u = 0

    IF(PRESENT(q_old)) THEN
       ALLOCATE(old_element_q(SIZE(old_connectivity,2), SIZE(q_old,2),2))
       old_element_q = 0
    ENDIF


    n_elements = SIZE(old_connectivity,1)
    np_perelem = SIZE(old_connectivity,2)
    target_points = new_coordinates(colint(TRANSPOSE(new_connectivity)),:)
    interpolation_points = target_points


    point_elements = 0
    nearest_candidate_elements = 0
    nearest_candidate_distances = HUGE(1.d0)
    nearest_candidate_sizes = 0.d0

    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*      PROJECTING SOLUTION TO NEW MESH          *'
          WRITE (6, *) '*************************************************'
       END IF
    ENDIF

    CALL find_points_in_linear_elements(target_points, old_connectivity, old_coordinates, point_elements)

    IF(ANY(point_elements .EQ. 0)) THEN
       CALL find_points_in_curved_elements(target_points, old_connectivity, old_coordinates, point_elements)
    ENDIF

    IF(ANY(point_elements .EQ. 0)) THEN
       nearest_tol = 2.d-4

       DO i=1,SIZE(target_points,1)
          IF(point_elements(i) .NE. 0) CYCLE
          target_point(1,:) = target_points(i,:)
          CALL find_nearest_curved_element(target_point, old_connectivity, old_coordinates, 0.5d0, nearest_element, &
               nearest_valid_point, nearest_distance, nearest_element_size)
          nearest_candidate_elements(i) = nearest_element
          nearest_candidate_distances(i) = nearest_distance
          nearest_candidate_sizes(i) = nearest_element_size

          IF(nearest_element .NE. 0 .AND. nearest_distance .LE. nearest_tol*MAX(1.d-12, nearest_element_size)) THEN
             point_elements(i) = nearest_element
             interpolation_points(i,:) = nearest_valid_point(1,:)
          ENDIF
       ENDDO

    ENDIF

    local_missing = COUNT(point_elements .EQ. 0)
    global_missing = local_missing
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, global_missing, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF(MPIvar%glob_id .EQ. 0 .AND. global_missing .NE. 0) THEN
       WRITE(*,*) "Projection unmatched points after fallbacks: ", global_missing
    ENDIF

    IF(local_missing .NE. 0) THEN
      WRITE(*,*) "Projection unmatched points on rank: ", MPIvar%glob_id, local_missing
      DO i=1,SIZE(target_points,1)
         IF(point_elements(i) .NE. 0) CYCLE

         IF(nearest_candidate_elements(i) .NE. 0) THEN
            WRITE(*,*) "Unmatched projection point rank/index/xy/nearest/dist/local h/ratio: ", &
                 MPIvar%glob_id, i, target_points(i,1), target_points(i,2), nearest_candidate_elements(i), &
                 nearest_candidate_distances(i), nearest_candidate_sizes(i), &
                 nearest_candidate_distances(i)/MAX(1.d-12, nearest_candidate_sizes(i))
         ELSE
            WRITE(*,*) "Unmatched projection point rank/index/xy/no nearest candidate: ", &
                 MPIvar%glob_id, i, target_points(i,1), target_points(i,2)
         ENDIF
      ENDDO
      WRITE(*,*) "Couldn't find a point in projection. STOP."
    ENDIF


!!$OMP parallel private(iel, point_indices, reference_points, shape_functions, element_coordinates, element_points, old_element_dofs, old_element_u, old_element_q) shared(interpolation_points, old_coordinates, old_connectivity, u_new, q_new, point_elements, n_elements, refElPol, np_perelem)
!!$OMP DO SCHEDULE(STATIC)
    DO iel = 1, n_elements
       IF(iel .EQ. 0) CYCLE

       CALL find_matches_int(point_elements, iel, point_indices)

       ALLOCATE(reference_points(SIZE(point_indices), SIZE(interpolation_points,2)))
       ALLOCATE(element_points(SIZE(point_indices), SIZE(interpolation_points,2)))
       ALLOCATE(shape_functions(np_perelem, SIZE(point_indices), 3))

       reference_points = 0
       element_points = 0
       shape_functions = 0

       element_coordinates = old_coordinates(old_connectivity(iel,:),:)
       element_points = interpolation_points(point_indices,:)

       ! this goddamn function always gives problems
       CALL inverse_isop_transf(element_points, element_coordinates, refElPol, reference_points)

       ! just fucking brute force it
       DO j = 1, SIZE(reference_points,2)
          DO i = 1, SIZE(reference_points,1)
             IF(ABS(reference_points(i,j)-1.0) .LT. 1e-12) THEN
                reference_points(i,j) = reference_points(i,j) - 1.e-10
             ENDIF
          ENDDO
       ENDDO

       CALL compute_shape_functions_at_points(refElPol, reference_points, shape_functions)

       old_element_dofs = (iel-1)*np_perelem + (/ (j, j=1, np_perelem) /)

       old_element_u = u_old(old_element_dofs, :)
       u_new(point_indices,:) = MATMUL(TRANSPOSE(shape_functions(:,:,1)), old_element_u)

       IF(PRESENT(q_old)) THEN
          old_element_q = q_old(old_element_dofs, :, :)
          q_new(point_indices,:,1) = MATMUL(TRANSPOSE(shape_functions(:,:,1)), old_element_q(:,:,1))
          q_new(point_indices,:,2) = MATMUL(TRANSPOSE(shape_functions(:,:,1)), old_element_q(:,:,2))
       ENDIF

       DEALLOCATE(reference_points)
       DEALLOCATE(element_points)
       DEALLOCATE(point_indices)
       DEALLOCATE(shape_functions)
    END DO
!!$OMP END DO
!!$OMP end parallel


    DEALLOCATE(old_element_u)
    IF(PRESENT(q_old)) THEN
       DEALLOCATE(old_element_q)
    ENDIF

  END SUBROUTINE projectSolutionDifferentMeshes_Mod


#ifdef WITH_PETSC
  SUBROUTINE InitPETSC
#include "petsc/finclude/petsc.h"
    USE petsc, ONLY: PetscInitialize, PetscInitialized
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
    USE petsc, ONLY: PetscFinalize, PetscFinalized
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
