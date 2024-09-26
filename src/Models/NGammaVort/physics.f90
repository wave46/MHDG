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
    phys%Neq = 4

    ! number of physical variables
    phys%npv = 5

    ALLOCATE (phys%phyVarNam(phys%npv))
    ALLOCATE (phys%conVarNam(phys%Neq))

    ! Set the name of the physical variables
    phys%phyVarNam(1) = "rho"
    phys%phyVarNam(2) = "u"
    phys%phyVarNam(3) = "Vorticity"
    phys%phyVarNam(4) = "Electric potential"
    phys%phyVarNam(5) = "Mach"

    ! Set the name of the conservative variables
    IF (switch%logrho) THEN
       phys%conVarNam(1) = "log(rho)"
    ELSE
       phys%conVarNam(1) = "rho"
    ENDIF
    phys%conVarNam(2) = "Gamma"
    phys%conVarNam(3) = "Vorticity"
    phys%conVarNam(4) = "Electric potential"

    simpar%model = 'N-Gamma-Vorticity'
    simpar%Ndim = 2
#ifdef TOR3D
    simpar%Ndim = 3
#endif
    simpar%Neq = phys%Neq
    ALLOCATE (simpar%physvar_refval(phys%npv))
    ALLOCATE (simpar%consvar_refval(phys%Neq))
    simpar%physvar_refval(1) = simpar%refval_density
    simpar%physvar_refval(2) = simpar%refval_speed
    simpar%physvar_refval(3) = simpar%refval_vorticity
    simpar%physvar_refval(4) = simpar%refval_potential
    simpar%physvar_refval(5) = 1.
    IF (switch%logrho) THEN
       simpar%consvar_refval(1) = LOG(simpar%refval_density)
    ELSE
       simpar%consvar_refval(1) = simpar%refval_density
    ENDIF
    simpar%consvar_refval(2) = simpar%refval_momentum
    simpar%consvar_refval(3) = simpar%refval_vorticity
    simpar%consvar_refval(4) = simpar%refval_potential
  END SUBROUTINE initPhys

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE phys2cons(up,ua)
    REAL*8,DIMENSION(:,:),INTENT(in)  :: up
    REAL*8,DIMENSION(:,:),INTENT(out) :: ua

    IF (switch%logrho) THEN
       ua(:, 1) = LOG(up(:, 1))
    ELSE
       ua(:, 1) = up(:, 1)
    ENDIF
    ua(:,2) = up(:,1)*up(:,2)
    ua(:,3) = up(:,3)
    ua(:,4) = up(:,4)
  END SUBROUTINE phys2cons

  !*******************************************
  ! Convert conservative variable to physical
  ! variables
  !*******************************************
  SUBROUTINE cons2phys(ua,up)
    REAL*8,DIMENSION(:,:),INTENT(in)  :: ua
    REAL*8,DIMENSION(:,:),INTENT(out) :: up

    IF (switch%logrho) THEN
       up(:, 1) = EXP(ua(:, 1))
    ELSE
       up(:, 1) = ua(:, 1)
    ENDIF
    up(:,2) = ua(:,2)/up(:, 1)
    up(:,3) = ua(:,3)
    up(:,4) = ua(:,4)
    up(:,5) = ua(:,2)/up(:, 1)/SQRT(phys%a)
  END SUBROUTINE cons2phys

  SUBROUTINE jacobianMatrices(U,A)
    REAL*8,INTENT(in)  :: U(:)
    REAL*8,INTENT(out) :: A(:,:)
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
       A(1, 2) = 1.
       A(2, 1) = (-1*U(2)**2/dens**2 + phys%a)
       A(2, 2) = 2*U(2)/dens
    ENDIF

    IF (switch%convvort) THEN
       A(3,1) = -U(3)*U(2)/dens**2
       A(3,2) = U(3)/dens
       A(3,3) = U(2)/dens
    END IF
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
    REAL*8 :: dens

    V = 0.d0
    V(1) = U(2)*U(1)/Up(1)
    V(2) = (U(1)-1)*(U(2)**2/Up(1) - phys%a*Up(1))

  END SUBROUTINE logrhojacobianVector


  !*****************************************
  ! Set the perpendicular diffusion
  !****************************************
  SUBROUTINE setLocalDiff(xy,d_iso,d_ani,Bmod)
    REAL*8,INTENT(in)  :: xy(:,:)
    REAL*8,INTENT(out) :: d_iso(:,:,:),d_ani(:,:,:)
    REAL*8              :: Bmod(:)
    REAL*8              :: iperdiff(SIZE(xy,1))
    REAL*8              :: coeff(SIZE(Bmod))

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
    coeff = 1./Bmod**2

    !*****************************
    ! Diagonal terms
    !*****************************
    d_iso(1,1,:) = phys%diff_n
    d_iso(2,2,:) = phys%diff_u
    d_iso(3,3,:) = phys%diff_vort
    d_iso(4,4,:) = coeff

    d_ani(1,1,:) = phys%diff_n
    d_ani(2,2,:) = phys%diff_u
    d_ani(3,3,:) = phys%diff_vort
    d_ani(4,4,:) = coeff

    !*****************************
    ! Non diagonal terms
    !*****************************
    !   !! term +\Div(1/eta \Grad_par \phi) in the vorticity equation
    d_iso(3,4,:) = 0.
    d_ani(3,4,:) = 1./phys%etapar*phys%c1 ! Negative parallel diffusion
    !   !! term -\Div(1/n\Grad_par \n) in the vorticity equation
    d_iso(3,1,:) = 0.
    d_ani(3,1,:) = -1./phys%etapar*phys%c2
    !   !! term -\Div(1/B^2\Grad_perp n) in the potential equation
    d_iso(4,1,:) = coeff*phys%Mref
    d_ani(4,1,:) = coeff*phys%Mref





    !d_iso(4,1,:) = 0.
    !d_ani(4,1,:) = 0.

    !d_iso(3,1,:) = 0.
    !d_ani(3,1,:) = 0.
    !
    !d_iso(3,4,:) = 0.
    !d_ani(3,4,:) = 0.


    IF (switch%testcase.GE.5 .AND. switch%testcase.LE.7) THEN
       d_iso(4,1,:) = 0.
       d_ani(4,1,:) = 0.
    ENDIF

    !write(6,*) "diff_pot: ",1./Bmod**2*phys%Mref*0.0232

    IF (switch%testcase .EQ. 1) THEN
       d_iso(4,4,:) = phys%diff_pot
       d_ani(4,4,:) = phys%diff_pot
       d_iso(3,4,:) = 0.
       d_ani(3,4,:) = phys%diff_pari ! Negative parallel diffusion
       !! term -\Div(1/n\Grad_par \n) in the vorticity equation
       d_iso(3,1,:) = 0.
       d_ani(3,1,:) = -phys%diff_pare
       !! term -\Div(1/n 1/B^2\Grad_perp n) in the potential equation
       d_iso(4,1,:) = phys%diff_ee
       d_ani(4,1,:) = phys%diff_ee
    ENDIF
    IF (switch%difcor .GT. 0) THEN
       CALL computeIperDiffusion(xy,iperdiff)
       d_iso(1,1,:) = d_iso(1,1,:)*iperdiff
       d_iso(2,2,:) = d_iso(2,2,:)*iperdiff
       d_iso(3,3,:) = d_iso(3,3,:)*iperdiff
       d_iso(4,4,:) = d_iso(4,4,:)*iperdiff
    ENDIF

  END SUBROUTINE setLocalDiff

  !*******************************************
  ! Compute local diffusion in points
  !*******************************************
  SUBROUTINE computeIperDiffusion(X,ipdiff)
    REAL*8,INTENT(IN)  :: X(:,:)
    REAL*8,INTENT(OUT) :: ipdiff(:)
    REAL*8             :: xcorn,ycorn
    REAL*8             :: rad(SIZE(X,1))
    REAL*8             :: h
    REAL*8,PARAMETER   :: tol = 1.e-6
    INTEGER            :: i

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
       WRITE (6,*) "Case not valid"
       STOP
    END SELECT

    !!**********************************************************
    !! Gaussian around the corner
    !!**********************************************************
    h = 10e-3
    rad = SQRT((X(:,1)*phys%lscale - xcorn)**2 + (X(:,2)*phys%lscale - ycorn)**2)
    ipdiff = 1 + numer%dc_coe*EXP(-(2*rad/h)**2)

  END SUBROUTINE computeIperDiffusion

  !***********************************************************************
  !
  !    COMPUTATION OF THE STABILIZATION PARAMETER
  !
  !***********************************************************************

  !*******************************************
  ! Compute the stabilization tensor tau
  !*******************************************
  SUBROUTINE computeTauGaussPoints(up,uc,b,bmod,n,iel,ifa,isext,xy,tau)
    REAL*8,INTENT(in)  :: up(:),uc(:),b(:),bmod,n(:),xy(:)
    REAL,INTENT(in) :: isext
    INTEGER,INTENT(in)  :: ifa,iel
    REAL*8,INTENT(out) :: tau(:,:)
    REAL*8              :: tau_aux(4),bort(2)
    REAL*8 :: xc,yc,rad,h,aux,bn,bnorm,bortn,coef_scale
    REAL*8 :: U1,U2,U3,U4
    INTEGER :: ndim
    U1 = uc(1)
    U2 = uc(2)
    U3 = uc(3)
    U4 = uc(4)

    tau = 0.
    ndim = SIZE(n)
    bn = dot_PRODUCT(b(1:ndim),n)
    bnorm = NORM2(b(1:ndim))
    bort(1) = b(1)
    bort(2) = -b(2)
    bortn = dot_PRODUCT(bort,n)




    bortn=1.


    IF (numer%stab == 2 .OR. numer%stab == 3) THEN
       tau_aux = ABS(uc(2)*bn/uc(1))
#ifdef TOR3D
       IF (ABS(n(3)) > 0.1) THEN
          ! Poloidal face
          coef_scale = refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
          tau_aux(1) = tau_aux(1) + phys%diff_n*ABS(bortn)*coef_scale
          tau_aux(2) = tau_aux(2) + phys%diff_u*ABS(bortn)*coef_scale
          tau_aux(3) = tau_aux(3) + (phys%Mref/phys%etapar*ABS(bn)+phys%diff_vort*ABS(bortn))*coef_scale
          tau_aux(4) = tau_aux(4) + phys%Mref/bmod**2*ABS(bortn)*coef_scale
       ELSE
#endif


          ! Toroidal face
          coef_scale = refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          coef_scale = 1.
          tau_aux(1) = tau_aux(1) + phys%diff_n*ABS(bortn)*coef_scale
          tau_aux(2) = tau_aux(2) + phys%diff_u*ABS(bortn)*coef_scale
          tau_aux(3) = tau_aux(3) + (1./phys%etapar*ABS(bn)+phys%diff_vort*ABS(bortn))*coef_scale
          tau_aux(4) = tau_aux(4) + 1./bmod**2*ABS(bortn)*coef_scale
          IF (ABS(isext - 1.) < 1e-8) THEN
             tau_aux(3) = phys%diff_vort*ABS(bortn)*coef_scale
          ENDIF
#ifdef TOR3D
       ENDIF
#endif
    ELSE
       WRITE (6,*) "Wrong stabilization type: ",numer%stab
       STOP
    ENDIF
    tau(1,1) = tau_aux(1)
    tau(2,2) = tau_aux(2)
    tau(3,3) = tau_aux(3)
    tau(4,4) = tau_aux(4)

    !tau = tau/(refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale)
  END SUBROUTINE computeTauGaussPoints

  SUBROUTINE computeTauGaussPoints_matrix(up,uc,b,n,xy,isext,iel,tau)
    REAL*8,INTENT(in)  :: up(:),uc(:),b(:),n(:),xy(:),isext
    REAL*8,INTENT(out) :: tau(:,:)
    INTEGER,INTENT(in) :: iel
    REAL*8              :: bn,bnorm
    REAL*8,PARAMETER :: eps = 1e-12
    REAL*8 :: U1,U2,U3,U4
    REAL*8 :: t2,t3,t4,t5,t6,t7,t8,t9
    REAL*8 :: t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
    REAL*8 :: t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
    REAL*8 :: t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
    REAL*8 :: t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
    REAL*8 :: t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
    REAL*8 :: t60,t61,t62,t63,t64,t65,t66,t67,t68,t69
    REAL*8 :: t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80
    REAL*8 :: x,y,r,h,coef,r0,rc

    x = xy(1)
    y = xy(2)

    U1 = uc(1)
    U2 = uc(2)
    U3 = uc(3)
    U4 = uc(4)

    bn = dot_PRODUCT(b,n)
    bnorm = NORM2(b)
    !************************************
    !
    ! *****     CONVECTIVE PART  ********
    !
    !************************************
    tau(1,1) = ABS((uc(2)*bn)/uc(1))
    tau(2,2) = ABS((uc(2)*bn)/uc(1))
    tau(3,3) = ABS((uc(2)*bn)/uc(1))
    tau(4,4) = ABS((uc(2)*bn)/uc(1))

    !************************************
    !
    ! *****     DIFFUSIVE PART  ********
    !
    !************************************
    tau(1,1) = tau(1,1) + phys%diff_n
    tau(2,2) = tau(2,2) + phys%diff_u
  END SUBROUTINE computeTauGaussPoints_matrix


END MODULE physics
