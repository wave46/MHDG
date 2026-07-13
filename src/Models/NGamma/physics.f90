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
      phys%phyVarNam(4) = "rhon"   ! density neutral
#endif

      ! Set the name of the conservative variables
      IF (switch%logrho) THEN
         phys%conVarNam(1) = "log(rho)"
      ELSE
         phys%conVarNam(1) = "rho"
      END IF
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
      END IF

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
      REAL*8, DIMENSION(:, :), INTENT(IN)  :: up
      REAL*8, DIMENSION(:, :), INTENT(OUT) :: ua

      IF (switch%logrho) THEN
         ua(:, 1) = LOG(up(:, 1))
      ELSE
         ua(:, 1) = up(:, 1)
      END IF
      ua(:, 2) = up(:, 1)*up(:, 2)
#ifdef NEUTRAL
      ua(:, 3) = up(:, 4)
#endif
   END SUBROUTINE phys2cons

   !*******************************************
   ! Convert conservative variable to physical
   ! variables
   !*******************************************
   PURE SUBROUTINE cons2phys(ua, up)
      REAL*8, DIMENSION(:, :), INTENT(IN)  :: ua
      REAL*8, DIMENSION(:, :), INTENT(OUT) :: up
      REAL*8, DIMENSION(SIZE(up, 1))       :: dens


      IF (switch%logrho) THEN
         up(:, 1) = EXP(ua(:, 1))                ! density
      ELSE
         up(:, 1) = ua(:, 1)
      END IF

      IF (switch%thresh .NE. 0) THEN
         dens = MAX(up(:, 1),numer%thr)
      ELSE
         dens = up(:, 1)
      END IF

      up(:, 2) = ua(:, 2)/dens              ! u parallel
      up(:, 3) = ua(:, 2)/dens/SQRT(phys%a) ! Mach

#ifdef NEUTRAL
      up(:, 4) = ABS(ua(:, 3))                                                 ! density neutral
#endif

   END SUBROUTINE cons2phys

   ! ******************************
   ! Split diffusion terms
   ! ******************************
   PURE SUBROUTINE compute_W2(U, W2, diff_n, diff_u)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8, INTENT(IN)        :: diff_n, diff_u
      REAL*8, INTENT(OUT)       :: W2(:)
      W2 = 0.
      W2(1) = (diff_n - diff_u)*U(2)/U(1)
   END SUBROUTINE compute_W2

   PURE SUBROUTINE compute_dW2_dU(U, dW2_dU, diff_n, diff_u)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8, INTENT(IN)        :: diff_n, diff_u
      REAL*8, INTENT(OUT)       :: dW2_dU(:, :)
      dW2_dU = 0.
      dW2_dU(1, 1) = -U(2)/(U(1)**2)
      dW2_dU(1, 2) = 1./U(1)

      dW2_dU = (diff_n - diff_u)*dW2_dU
   END SUBROUTINE compute_dW2_dU

   PURE SUBROUTINE jacobianMatrices(U, A)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8, INTENT(OUT)       :: A(:, :)
      REAL*8                    :: dens

      A = 0.d0
      IF (switch%logrho) THEN
         dens = EXP(U(1))
         A(1, 1) = -U(2)/dens
         A(1, 2) = 1/dens
         A(2, 1) = -U(2)**2/dens + phys%a*dens
         A(2, 2) = 2*U(2)/dens
      ELSE
         dens = U(1)
         A(1, 2) = 1.
         A(2, 1) = (-1*U(2)**2/dens**2 + phys%a)
         A(2, 2) = 2*U(2)/dens
      END IF

   END SUBROUTINE jacobianMatrices

   !*****************************************
   ! Jacobian matrix for face computations
   !****************************************
   SUBROUTINE jacobianMatricesFace(U, bn, An)
      REAL*8, INTENT(IN)        :: U(:), bn
      REAL*8, INTENT(OUT)       :: An(:, :)

      An = 0.
      CALL jacobianMatrices(U, An)
      An = bn*An
   END SUBROUTINE jacobianMatricesFace

   PURE SUBROUTINE logrhojacobianVector(U,Up,V)
      REAL*8, INTENT(IN)        :: U(:),Up(:)
      REAL*8, INTENT(OUT)       :: V(:)

      V = 0.d0
      V(1) = U(2)*U(1)/Up(1)
      V(2) = (U(1)-1)*(U(2)**2/Up(1) - phys%a*Up(1))

   END SUBROUTINE logrhojacobianVector


   !*****************************************
   ! Set the perpendicular diffusion
   !****************************************
   SUBROUTINE setLocalDiff(xy, u, d_iso, d_ani)
      REAL*8, INTENT(IN)        :: xy(:, :)
      REAL*8, INTENT(IN)        :: u(:, :)
      REAL*8, INTENT(OUT)       :: d_iso(:, :, :), d_ani(:, :, :)
      REAL*8                    :: iperdiff(SIZE(d_iso, 3))

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
      d_iso(3, 3, :) = phys%diff_nn
#endif

      d_ani(1, 1, :) = phys%diff_n
      d_ani(2, 2, :) = phys%diff_u
#ifdef NEUTRAL
      d_ani(3, 3, :) = 0.
#endif

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
      REAL*8, INTENT(IN)        :: X(:, :), u(:, :)
      REAL*8, INTENT(OUT)       :: ipdiff(:)
      REAL*8                    :: xcorn, ycorn, d, dref, maxamp
      REAL*8                    :: rad(SIZE(X, 1))
      REAL*8                    :: h, rhog
      INTEGER                   :: g, opt

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
            WRITE (6, *) "case not valid"
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
      IF (switch%limrho .EQ. 2 .AND. MINVAL(u(:, 1)) .LT. numer%minrho) THEN
         DO g = 1, SIZE(u, 1)
            rhog = u(g, 1)
            IF (rhog < 0.) rhog = 0.
            IF (rhog < numer%minrho) THEN
               d = numer%minrho - rhog ! 0 < d < minrho
               IF (opt .EQ. 1) THEN
                  dref = maxamp*d/numer%minrho  ! 0 < dref < maxamp
                  !                                    ipdiff(g) = exp(dref) ! 1 < ipdiff(g) < exp(maxamp)
                  ipdiff(g) = ipdiff(g) + EXP(dref) - 1 ! 0 < ipdiff(g) < exp(maxamp)-1
               ELSE IF (opt == 2) THEN
                  ipdiff(g) = ipdiff(g) + 1./((1.-d/numer%minrho)*2 + 1./50.) - 1.
               END IF
            END IF
         END DO
      END IF
   END SUBROUTINE computeIperDiffusion

   !*****************************************
   ! Pinch term
   !***************************************
   PURE SUBROUTINE computePinch(b, psi, APinch)
      REAL*8, INTENT(IN)        :: b(:), psi
      REAL*8, INTENT(OUT)       :: APinch(:, :)
      REAL*8                    :: v_p, bnorm(2)

      APinch = 0.

      bnorm = b(:)/NORM2(b)
      v_p = phys%v_p*(psi**2 + psi**2*TANH((0.95 - psi)/0.02))
      !if (v_p .lt. 1.e-4/simpar%refval_speed) v_p = 0.

      APinch(1, 1) = v_p*bnorm(2)
      APinch(1, 2) = v_p*(-bnorm(1))

   END SUBROUTINE computePinch

   ! ******************************
   ! Neutral Source terms
   ! ******************************
#ifdef NEUTRAL
   SUBROUTINE compute_RN(E, theta, RN)
      ! Compute the recycling coefficeint RN(E,theta) interpolating the TRIM data
      USE interpolation
      REAL*8, INTENT(IN)        :: E, theta
      REAL*8, INTENT(OUT)       :: RN
      REAL*8                    :: E_clipped, theta_clipped
      INTEGER                   :: ip, jp

      RN = 1.

      ip = SIZE(phys%E)
      jp = SIZE(phys%theta)

      E_clipped = MAX(1e-20, MIN(1e3 - 1e-20, E))
      theta_clipped = MAX(1e-20, MIN(90 - 1e-20, theta))

      RN = interpolate(ip, phys%E, jp, phys%theta, phys%RN_DW, E_clipped, theta_clipped, 1e-12)

   END SUBROUTINE compute_RN

   SUBROUTINE compute_niz(U, niz)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8                    :: niz, U1, U3
      REAL, PARAMETER           :: tol = 1e-10
      U1 = U(1)
      U3 = U(3)
      IF (U1 < tol) U1 = tol
      IF (U3 < tol) U3 = tol
      niz = U1*U3
   END SUBROUTINE compute_niz

   SUBROUTINE compute_dniz_dU(U, res)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8                    :: res(:), U1, U3
      REAL, PARAMETER           :: tol = 1e-10
      U1 = U(1)
      U3 = U(3)
      IF (U1 < tol) U1 = tol
      IF (U3 < tol) U3 = tol
      res = 0.
      res(1) = U3
      res(3) = U1
   END SUBROUTINE compute_dniz_dU

   SUBROUTINE compute_nrec(U, nrec)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8                    :: nrec, U1
      REAL, PARAMETER           :: tol = 1e-10
      U1 = U(1)
      IF (U1 < tol) U1 = tol
      nrec = U1**2
   END SUBROUTINE compute_nrec

   SUBROUTINE compute_dnrec_dU(U, res)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8                    :: res(:), U1
      REAL, PARAMETER           :: tol = 1e-10
      U1 = U(1)
      IF (U1 < tol) U1 = tol
      res = 0.
      res(1) = 2.*U1
   END SUBROUTINE compute_dnrec_dU

   SUBROUTINE compute_fGammacx(U, fGammacx)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8                    :: fGammacx, U2, U3
      REAL, PARAMETER           :: tol = 1e-10
      U2 = U(2)
      U3 = U(3)
      IF (U3 < tol) U3 = tol
      fGammacx = U2*U3
   END SUBROUTINE compute_fGammacx

   SUBROUTINE compute_dfGammacx_dU(U, res)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8                    :: res(:), U2, U3
      REAL, PARAMETER           :: tol = 1e-10
      U2 = U(2)
      U3 = U(3)
      IF (U3 < tol) U3 = tol
      res = 0.
      res(2) = U3
      res(3) = U2
   END SUBROUTINE compute_dfGammacx_dU

   SUBROUTINE compute_fGammarec(U, fGammarec)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8                    :: fGammarec, U1, U2
      REAL, PARAMETER           :: tol = 1e-10
      U1 = U(1)
      U2 = U(2)
      IF (U1 < tol) U1 = tol
      fGammarec = U1*U2
   END SUBROUTINE compute_fGammarec

   SUBROUTINE compute_dfGammarec_dU(U, res)
      REAL*8, INTENT(IN)        :: U(:)
      REAL*8                    :: res(:), U1, U2
      REAL, PARAMETER           :: tol = 1e-10
      U1 = U(1)
      U2 = U(2)
      IF (U1 < tol) U1 = tol
      res = 0.
      res(1) = U2
      res(2) = U1
   END SUBROUTINE compute_dfGammarec_dU

#endif
   !NEUTRAL

   !*******************************************
   ! Compute the stabilization tensor tau
   !*******************************************
   SUBROUTINE computeTauGaussPoints(up, uc, b, n, iel, isext, xy, edge_length, tau)
      REAL*8, INTENT(IN)        :: up(:), uc(:), b(:), n(:), xy(:)
      REAL, INTENT(IN)          :: isext
      REAL, INTENT(IN)          :: edge_length
      INTEGER, INTENT(IN)       :: iel
      REAL*8, INTENT(OUT)       :: tau(:, :)
      INTEGER                   :: ndim
#ifdef NEUTRAL
      REAL*8                    :: tau_aux(3), diff_iso(3, 3, 1), diff_ani(3, 3, 1)
#else
      REAL*8                    :: tau_aux(2), diff_iso(2, 2, 1), diff_ani(2, 2, 1)
#endif
      REAL*8                    :: bn, bn_perp, xyd(1, SIZE(xy)), uu(1, SIZE(uc))

      tau = 0.
      ndim = SIZE(n)
      bn = dot_PRODUCT(b(1:ndim), n)

      bn_perp = dot_PRODUCT([b(1) - b(2)], n)

      xyd(1, :) = xy(:)
      uu(1, :) = uc(:)

      CALL setLocalDiff(xyd, uu, diff_iso, diff_ani)

      IF (numer%stab == 2) THEN

         tau_aux = ABS(up(2)*bn)
         tau_aux(1) = tau_aux(1) + diff_iso(1, 1, 1)
         tau_aux(2) = tau_aux(2) + diff_iso(2, 2, 1)
#ifdef NEUTRAL
         tau_aux(3) = diff_iso(3, 3, 1)
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
         tau_aux(3) = phys%diff_nn
#endif

      ELSEIF (numer%stab == 5) THEN

#ifdef TOR3D
         IF (ABS(n(3)) > 0.1) THEN
            tau_aux = MAX(ABS((uc(2)/uc(1) + SQRT(phys%a))*bn), ABS((uc(2)/uc(1) - SQRT(phys%a))*bn))
            ! Poloidal face
            tau_aux(1) = tau_aux(1) + phys%diff_n/refElPol%nDeg*mesh%elemSize(iel)
            tau_aux(2) = tau_aux(2) + phys%diff_u/refElPol%nDeg*mesh%elemSize(iel)
#ifdef NEUTRAL
            tau_aux(3) = phys%diff_nn/refElPol%nDeg*mesh%elemSize(iel)
#endif
         ELSE
#endif
            bn = dot_PRODUCT(b(1:2), n(1:2))
            tau_aux = MAX(ABS((uc(2)/uc(1) + SQRT(phys%a))*bn), ABS((uc(2)/uc(1) - SQRT(phys%a))*bn))

            tau_aux(1) = tau_aux(1) + phys%diff_n!/edge_length*abs(bn_perp)!/refElPol%nDeg*mesh%elemSize(iel)
            tau_aux(2) = tau_aux(2) + phys%diff_u!/edge_length*abs(bn_perp)!/refElPol%nDeg*mesh%elemSize(iel)
#ifdef NEUTRAL
            tau_aux(3) = phys%diff_nn/refElPol%nDeg*mesh%elemSize(iel)
#endif
#ifdef TOR3D
         END IF
#endif
      ELSE
         WRITE (6, *) "Wrong stabilization type: ", numer%stab
         STOP
      END IF

      tau(1, 1) = tau_aux(1)
      tau(2, 2) = tau_aux(2)
#ifdef NEUTRAL
      tau(3, 3) = tau_aux(3)
#endif

   END SUBROUTINE computeTauGaussPoints

   SUBROUTINE computeTauGaussPoints_matrix(up, uc, b, n, xy, isext, iel, tau)
      REAL*8, INTENT(IN)        :: up(:), uc(:), b(:), n(:), xy(:), isext
      REAL*8, INTENT(OUT)       :: tau(:, :)
      INTEGER, INTENT(IN)       :: iel
      REAL*8                    :: bn, bnorm
      REAL*8                    :: U1, U2, U3, U4
      REAL*8                    :: x, y

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
