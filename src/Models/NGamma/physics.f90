!*****************************************
! project: MHDG
! file: physics.f90
! date: 15/02/2017
! Define the physics of the model
!  ******** N-Gamma system Isothermal ****
!*****************************************
MODULE physics

  USE globals
  USE magnetic_field
  IMPLICIT NONE

CONTAINS

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE initPhys()

    ! number of equation of the problem
    phys%Neq = 2
#ifdef NEUTRAL
    phys%Neq = 3
#endif

    ! number of physical variables
    phys%npv = 3
#ifdef NEUTRAL
    phys%npv = 4
#endif

    ALLOCATE (phys%phyVarNam(phys%npv))
    ALLOCATE (phys%conVarNam(phys%Neq))

    ! Set the name of the physical variables
    phys%phyVarNam(1) = "rho"
    phys%phyVarNam(2) = "u"
    phys%phyVarNam(3) = "Mach"
#ifdef NEUTRAL
    phys%phyVarNam(4)= "rhon"   ! density neutral
#endif

    ! Set the name of the conservative variables
    IF (switch%logrho) THEN
       phys%conVarNam(1) = "log(rho)"
    ELSE
       phys%conVarNam(1) = "rho"
    ENDIF
    phys%conVarNam(2) = "Gamma"
#ifdef NEUTRAL
    phys%conVarNam(3) = "rhon"  ! U3 = rhon
#endif

    simpar%model = 'N-Gamma'
#ifdef NEUTRAL
    simpar%model = 'N-Gamma-Neutral'
#endif
    simpar%Ndim = 2
#ifdef TOR3D
    simpar%Ndim = 3
#endif
    simpar%Neq = phys%Neq
    ALLOCATE (simpar%physvar_refval(phys%npv))
    ALLOCATE (simpar%consvar_refval(phys%Neq))
    simpar%physvar_refval(1) = simpar%refval_density
    simpar%physvar_refval(2) = simpar%refval_speed
    simpar%physvar_refval(3) = 1.
#ifdef NEUTRAL
    simpar%physvar_refval(4) = simpar%refval_neutral
#endif
    IF (switch%logrho) THEN
       simpar%consvar_refval(1) = LOG(simpar%refval_density)
    ELSE
       simpar%consvar_refval(1) = simpar%refval_density
    ENDIF

    simpar%consvar_refval(2) = simpar%refval_momentum
#ifdef NEUTRAL
    simpar%consvar_refval(3) = simpar%refval_neutral
#endif

  END SUBROUTINE initPhys

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE phys2cons(up, ua)
    REAL*8, DIMENSION(:, :), INTENT(in)  :: up
    REAL*8, DIMENSION(:, :), INTENT(out) :: ua

    IF (switch%logrho) THEN
       ua(:, 1) = LOG(up(:, 1))
    ELSE
       ua(:, 1) = up(:, 1)
    ENDIF
    ua(:, 2) = up(:, 1)*up(:, 2)
#ifdef NEUTRAL
    ua(:,3) = up(:,4)
#endif
  END SUBROUTINE phys2cons

  !*******************************************
  ! Convert conservative variable to physical
  ! variables
  !*******************************************
  SUBROUTINE cons2phys(ua, up)
    REAL*8, DIMENSION(:, :), INTENT(in)  :: ua
    REAL*8, DIMENSION(:, :), INTENT(out) :: up
    REAL*8 :: dens(SIZE(up,1))
    INTEGER :: i

    IF (switch%logrho) THEN
       up(:, 1) = EXP(ua(:, 1))                ! density
    ELSE
       up(:, 1) = ua(:, 1)                     ! density
    ENDIF
    dens=up(:, 1)
    IF (switch%thresh.NE.0) THEN
       DO i=1,SIZE(up,1)
          dens(i) = MAX(dens(i),numer%thr)
       END DO
    ENDIF
    up(:, 2) = ua(:, 2)/dens              ! u parallel
    up(:, 3) = ua(:, 2)/dens/SQRT(phys%a) ! Mach
#ifdef NEUTRAL
    up(:,4) = ABS(ua(:,3))                    ! density neutral
#endif
  END SUBROUTINE cons2phys

  ! ******************************
  ! Split diffusion terms
  ! ******************************
  SUBROUTINE compute_W2(U,W2)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: W2(:)
    W2 = 0.
    W2(1) = (phys%diff_n-phys%diff_u)*U(2)/U(1)
  END SUBROUTINE compute_W2

  SUBROUTINE compute_dW2_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:,:)
    res = 0.
    res(1,1) = -U(2)*(phys%diff_n-phys%diff_u)/(U(1)**2)
    res(1,2) = 1.*(phys%diff_n-phys%diff_u)/U(1)
  END SUBROUTINE compute_dW2_dU

  SUBROUTINE jacobianMatrices(U, A)
    REAL*8, INTENT(in)  :: U(:)
    REAL*8, INTENT(out) :: A(:, :)
    REAL*8 :: dens

    A = 0.d0
    IF (switch%logrho) THEN
       dens=EXP(U(1))
       A(1,1) = -U(2)/dens
       A(1,2) = 1/dens
       A(2,1) = -U(2)**2/dens+phys%a*dens
       A(2,2) = 2*U(2)/dens
    ELSE
       dens=U(1)
       !      if (switch%thresh.ne.0) then
       !         dens = max(dens,numer%thr)
       !      endif
       A(1, 2) = 1.
       A(2, 1) = (-1*U(2)**2/dens**2 + phys%a)
       A(2, 2) = 2*U(2)/dens
    ENDIF

  END SUBROUTINE jacobianMatrices

  !*****************************************
  ! Jacobian matrix for face computations
  !****************************************
  SUBROUTINE jacobianMatricesFace(U,bn,An)
    REAL*8,INTENT(in)  :: U(:),bn
    REAL*8,INTENT(out) :: An(:,:)

    An = 0.
    CALL jacobianMatrices(U,An)
    An = bn*An
  END SUBROUTINE jacobianMatricesFace


  SUBROUTINE logrhojacobianVector(U,Up,V)
    REAL*8, INTENT(in)  :: U(:),Up(:)
    REAL*8, INTENT(out) :: V(:)

    V = 0.d0
    V(1) = U(2)*U(1)/Up(1)
    V(2) = (U(1)-1)*(U(2)**2/Up(1) - phys%a*Up(1))

  END SUBROUTINE logrhojacobianVector


  !*****************************************
  ! Set the perpendicular diffusion
  !****************************************
  SUBROUTINE setLocalDiff(xy, u, d_iso, d_ani, Bmod)
    REAL*8, INTENT(in)  :: xy(:, :)
    REAL*8, INTENT(in)  :: Bmod(:)
    REAL*8, INTENT(in)  :: u(:,:)
    REAL*8, INTENT(out) :: d_iso(:, :, :), d_ani(:, :, :)
    REAL*8              :: iperdiff(SIZE(xy, 1))

    ! d_iso(Neq,Neq,Ngauss),d_ani(Neq,Neq,Ngauss)
    ! first index corresponds to the equation
    ! second index corresponds to the unknown
    ! third index correspond to the gauss point
    ! example d_iso(2,3,ig) is the diffusion in the non-diagonal
    ! diffusion term in the second equation/third variable
    ! d_iso>0,d_ani=0 is an isotropic diffusion
    ! d_iso>0,d_ani=d_iso is a perpendicular diffusion
    ! d_iso=0,d_ani<0 is a (positive) parallel diffusion

    d_iso = 0.
    d_ani = 0.
    !*****************************
    ! Diagonal terms
    !*****************************
    d_iso(1, 1, :) = phys%diff_n
    d_iso(2, 2, :) = phys%diff_u
#ifdef NEUTRAL
    d_iso(3,3,:) = phys%diff_nn
#endif

    d_ani(1, 1, :) = phys%diff_n
    d_ani(2, 2, :) = phys%diff_u
#ifdef NEUTRAL
    d_ani(3,3,:) = 0.
#endif


    !write(6,*) "inside d_iso", d_iso(1,1,1),d_iso(2,2,1),d_iso(3,3,1)


    !*****************************
    ! Non diagonal terms
    !*****************************
    ! No non-diagonal terms defined for this model

    !    call computeIperDiffusion(xy, u, iperdiff)
    !    d_iso(1, 1, :) = d_iso(1, 1, :)*iperdiff
    !    d_iso(2, 2, :) = d_iso(2, 2, :)*iperdiff
    !#ifdef NEUTRAL
    !    d_iso(3, 3, :) = d_iso(3, 3, :)*iperdiff
    !#endif

    CALL computeIperDiffusion(xy, u, iperdiff)
    d_iso(1, 1, :) = d_iso(1, 1, :) + iperdiff
    d_iso(2, 2, :) = d_iso(2, 2, :) + iperdiff

#ifdef NEUTRAL
    !    d_iso(3, 3, :) = d_iso(3, 3, :) + iperdiff
#endif

    !    if (maxval(iperdiff)>1e-12) then
    !    write(6,*) "iperdiff: ", iperdiff
    !    endif
    !write(6,*) "u: ", u, "iperdiff ", iperdiff

  END SUBROUTINE setLocalDiff

  !*******************************************
  ! Compute local diffusion in points
  !*******************************************
  SUBROUTINE computeIperDiffusion(X, u, ipdiff)
    REAL*8, INTENT(IN)  :: X(:, :), u(:,:)
    REAL*8, INTENT(OUT) :: ipdiff(:)
    REAL*8              :: xcorn, ycorn,d,dref,maxamp
    REAL*8              :: rad(SIZE(X, 1))
    REAL*8              :: h,rhog
    INTEGER             :: g, opt

    ipdiff = 0.

    IF (switch%difcor .GT. 0) THEN
       SELECT CASE (switch%difcor)
       CASE (1)
          ! Circular case with infinitely small limiter
          xcorn = geom%R0
          ycorn = -0.75
       CASE (2)
          ! Circular case with infinitely small limiter
          xcorn = geom%R0
          ycorn = -0.287
       CASE (3)
          ! West
          xcorn = 2.7977
          ycorn = -0.5128
       CASE DEFAULT
          WRITE (6, *) "Case not valid"
          STOP
       END SELECT

       !!**********************************************************
       !! Gaussian around the corner
       !!**********************************************************
       h = 10e-3
       rad = SQRT((X(:, 1)*phys%lscale - xcorn)**2 + (X(:, 2)*phys%lscale - ycorn)**2)
       ipdiff = 1 + numer%dc_coe*EXP(-(2*rad/h)**2)
    END IF

    maxamp = 4.
    opt = 2
    IF (switch%limrho .EQ. 2 .AND. MINVAL(u(:,1)) .LT. numer%minrho) THEN
       DO g = 1,SIZE(u,1)
          rhog = u(g,1)
          IF (rhog<0.) rhog=0.
          IF (rhog<numer%minrho) THEN
             d = numer%minrho-rhog ! 0 < d < minrho
             IF (opt.EQ.1) THEN
                dref = maxamp * d/numer%minrho  ! 0 < dref < maxamp
                !			            ipdiff(g) = exp(dref) ! 1 < ipdiff(g) < exp(maxamp)
                ipdiff(g) = ipdiff(g) + EXP(dref) - 1 ! 0 < ipdiff(g) < exp(maxamp)-1
             ELSE IF (opt==2) THEN
                ipdiff(g) = ipdiff(g) +  1./((1.-d/numer%minrho)*2+1./50.) - 1.
             ENDIF
          ENDIF
       END DO
    ENDIF
  END SUBROUTINE computeIperDiffusion

  ! ******************************
  ! Neutral Source terms
  ! ******************************
#ifdef NEUTRAL

  SUBROUTINE compute_niz(U,niz)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: niz,U1,U3
    REAL,PARAMETER :: tol = 1e-10
    U1 = U(1)
    U3 = U(3)
    IF (U1<tol) U1=tol
    IF (U3<tol) U3=tol
    niz = U1*U3
  END SUBROUTINE compute_niz


  SUBROUTINE compute_dniz_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U3
    REAL,PARAMETER :: tol = 1e-10
    U1 = U(1)
    U3 = U(3)
    IF (U1<tol) U1=tol
    IF (U3<tol) U3=tol
    res = 0.
    res(1) = U3
    res(3) = U1
  END SUBROUTINE compute_dniz_dU


  SUBROUTINE compute_nrec(U,nrec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: nrec,U1
    REAL,PARAMETER :: tol = 1e-10
    U1 = U(1)
    IF (U1<tol) U1=tol
    nrec = U1**2
  END SUBROUTINE compute_nrec


  SUBROUTINE compute_dnrec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1
    REAL,PARAMETER :: tol = 1e-10
    U1 = U(1)
    IF (U1<tol) U1=tol
    res = 0.
    res(1) = 2.*U1
  END SUBROUTINE compute_dnrec_dU


  SUBROUTINE compute_fGammacx(U,fGammacx)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: fGammacx,U2,U3
    REAL,PARAMETER :: tol = 1e-10
    U2 = U(2)
    U3 = U(3)
    IF (U3<tol) U3=tol
    fGammacx = U2*U3
  END SUBROUTINE compute_fGammacx


  SUBROUTINE compute_dfGammacx_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U2,U3
    REAL,PARAMETER :: tol = 1e-10
    U2 = U(2)
    U3 = U(3)
    IF (U3<tol) U3=tol
    res = 0.
    res(2) = U3
    res(3) = U2
  END SUBROUTINE compute_dfGammacx_dU


  SUBROUTINE compute_fGammarec(U,fGammarec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: fGammarec,U1,U2
    REAL,PARAMETER :: tol = 1e-10
    U1 = U(1)
    U2 = U(2)
    IF (U1<tol) U1=tol
    fGammarec = U1*U2
  END SUBROUTINE compute_fGammarec


  SUBROUTINE compute_dfGammarec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U2
    REAL,PARAMETER :: tol = 1e-10
    U1 = U(1)
    U2 = U(2)
    IF (U1<tol) U1=tol
    res = 0.
    res(1) = U2
    res(2) = U1
  END SUBROUTINE compute_dfGammarec_dU

#endif
  !NEUTRAL




  !*******************************************
  ! Compute the stabilization tensor tau
  !*******************************************
  SUBROUTINE computeTauGaussPoints(up, uc, b, bmod, n, iel, ifa, isext, xy, tau)
    REAL*8, INTENT(in)  :: up(:), uc(:), b(:), bmod,n(:), xy(:)
    REAL, INTENT(in)    :: isext
    INTEGER, INTENT(in) :: ifa, iel
    REAL*8, INTENT(out) :: tau(:, :)
    INTEGER             :: ndim
#ifdef NEUTRAL
    REAL*8              :: tau_aux(3),diff_iso(3,3,1),diff_ani(3,3,1)
#else
    REAL*8              :: tau_aux(2),diff_iso(2,2,1),diff_ani(2,2,1)
#endif
    REAL*8              :: bn, bnorm,xyd(1,SIZE(xy)),uu(1,SIZE(uc)),bmod_d(1)


    ndim = SIZE(n)
    bn = dot_PRODUCT(b(1:ndim), n)
    bnorm = NORM2(b(1:ndim))
    xyd(1,:) = xy(:)
    uu(1,:) = uc(:)
    bmod_d = bmod

    CALL setLocalDiff(xyd, uu, diff_iso, diff_ani, bmod_d)


    !    write(6,*) "diff_iso", diff_iso(1,1,1),diff_iso(2,2,1),diff_iso(3,3,1), " phys%diff_n ", phys%diff_n, " phys%diff_u ",phys%diff_u,  " phys%diff_nn ",phys%diff_nn

    IF (numer%stab == 2) THEN
       !      tau_aux = abs(up(2)*bn)
       !      tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
       !      tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
       !#ifdef NEUTRAL
       !      tau_aux(3) =              phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
       !!      tau_aux(3) = tau_aux(3) + phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
       !#endif
       tau_aux = ABS(up(2)*bn)
       tau_aux(1) = tau_aux(1) + diff_iso(1,1,1)
       tau_aux(2) = tau_aux(2) + diff_iso(2,2,1)
#ifdef NEUTRAL
       tau_aux(3) =              diff_iso(3,3,1)
#endif

    ELSEIF (numer%stab == 3) THEN
       tau_aux = MAX(ABS((uc(2) + SQRT(phys%a))*bn/up(1)), ABS((uc(2) - SQRT(phys%a))*bn/up(1)))
       tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
       tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifdef NEUTRAL
       tau_aux(3) = phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#endif

    ELSEIF (numer%stab == 4) THEN
       tau_aux = ABS((up(2)*bn))
       tau_aux(1) = tau_aux(1) + phys%diff_n
       tau_aux(2) = tau_aux(2) + phys%diff_u
#ifdef NEUTRAL
       tau_aux(3) =  phys%diff_nn
#endif

    ELSEIF (numer%stab == 5) THEN
       tau_aux = MAX(ABS((uc(2) + SQRT(phys%a))*bn/up(1)), ABS((uc(2) - SQRT(phys%a))*bn/up(1)))
       tau_aux(1) = tau_aux(1) + phys%diff_n
       tau_aux(2) = tau_aux(2) + phys%diff_u
#ifdef NEUTRAL
       tau_aux(3) =  phys%diff_nn
#endif

    ELSE
       WRITE (6, *) "Wrong stabilization type: ", numer%stab
       STOP
    ENDIF
    tau(1, 1) = tau_aux(1)
    tau(2, 2) = tau_aux(2)
#ifdef NEUTRAL
    tau(3,3) = tau_aux(3)
#endif
  END SUBROUTINE computeTauGaussPoints

  SUBROUTINE computeTauGaussPoints_matrix(up, uc, b, n, xy, isext, iel, tau)
    REAL*8, INTENT(in)  :: up(:), uc(:), b(:), n(:), xy(:), isext
    REAL*8, INTENT(out) :: tau(:, :)
    INTEGER, INTENT(in) :: iel
    REAL*8              :: bn, bnorm
    REAL*8              :: U1, U2, U3, U4
    REAL*8              :: x, y

    x = xy(1)
    y = xy(2)

    U1 = uc(1)
    U2 = uc(2)
    U3 = uc(3)
    U4 = uc(4)

    bn = dot_PRODUCT(b, n)
    bnorm = NORM2(b)
    !************************************
    !
    ! *****     CONVECTIVE PART  ********
    !
    !************************************
    tau(1, 1) = ABS((uc(2)*bn)/uc(1))
    tau(2, 2) = ABS((uc(2)*bn)/uc(1))
    tau(3, 3) = ABS((uc(2)*bn)/uc(1))
    tau(4, 4) = ABS((uc(2)*bn)/uc(1))

    !************************************
    !
    ! *****     DIFFUSIVE PART  ********
    !
    !************************************
    tau(1, 1) = tau(1, 1) + phys%diff_n
    tau(2, 2) = tau(2, 2) + phys%diff_u
  END SUBROUTINE computeTauGaussPoints_matrix


END MODULE physics
