!*****************************************
! project: MHDG
! file: physics.f90
! date: 15/02/2017
! Define the physics of the model
!  ******** N-Gamma system Isothermal ****
!*****************************************
MODULE physics

  USE globals
  IMPLICIT NONE

CONTAINS

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE initPhys()

    ! number of equation of the problem
    phys%Neq = 3

    ! number of physical variables
    phys%npv = 3

    ALLOCATE (phys%phyVarNam(phys%npv))
    ALLOCATE (phys%conVarNam(phys%Neq))

    ! Set the name of the physical variables
    phys%phyVarNam(1) = "rho"
    phys%phyVarNam(2) = "Mach"
    phys%phyVarNam(3) = "u"

    ! Set the name of the conservative variables
    phys%conVarNam(1) = "rho"
    phys%conVarNam(2) = "Gamma"
    phys%conVarNam(3) = "u"

  END SUBROUTINE initPhys

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE phys2cons(up, ua)
    REAL*8, DIMENSION(:, :), INTENT(in)  :: up
    REAL*8, DIMENSION(:, :), INTENT(out) :: ua

    ua(:, 1) = up(:, 1)
    ua(:, 2) = up(:, 1)*up(:, 2)
    ua(:, 3) = up(:, 3)
  END SUBROUTINE phys2cons

  !*******************************************
  ! Convert conservative variable to physical
  ! variables
  !*******************************************
  SUBROUTINE cons2phys(ua, up)
    REAL*8, DIMENSION(:, :), INTENT(in)  :: ua
    REAL*8, DIMENSION(:, :), INTENT(out) :: up

    up(:, 1) = ua(:, 1)
    up(:, 2) = ua(:, 2)/ua(:, 1)/SQRT(phys%a)
    up(:, 3) = ua(:, 3)

  END SUBROUTINE cons2phys

  !*****************************************
  ! Jacobian matrices
  !****************************************
  !                 SUBROUTINE jacobianMatrices(U,b,Ax,Ay)
  !                 real*8, intent(in)  :: U(:),b(:)
  !                 real*8, intent(out) :: Ax(:,:),Ay(:,:)

  !                        Ax = 0.d0
  !                        Ax(1,2) = b(1)
  !                        Ax(2,1) = (-1*U(2)**2/U(1)**2+phys%a)*b(1)
  !                        Ax(2,2) = 2*U(2)/U(1)*b(1)

  !                        Ay = 0.d0
  !                        Ay(1,2) = b(2)
  !                        Ay(2,1) = (-1*U(2)**2/U(1)**2+phys%a)*b(2)
  !                        Ay(2,2) = 2*U(2)/U(1)*b(2)

  !                 END SUBROUTINE jacobianMatrices

  SUBROUTINE jacobianMatrices(U, A)
    REAL*8, INTENT(in)  :: U(:)
    REAL*8, INTENT(out) :: A(:, :)

    A = 0.d0
    A(1, 2) = 1.
    A(2, 1) = (-1*U(2)**2/U(1)**2 + phys%a)
    A(2, 2) = 2*U(2)/U(1)

  END SUBROUTINE jacobianMatrices

  !*****************************************
  ! Jacobian matrix for face computations
  !****************************************
  SUBROUTINE jacobianMatricesFace(U, bn, An)
    REAL*8, INTENT(in)  :: U(:), bn
    REAL*8, INTENT(out) :: An(:, :)

    An = 0.d0
    An(1, 2) = bn
    An(2, 1) = (-1*U(2)**2/U(1)**2 + phys%a)*bn
    An(2, 2) = 2*U(2)/U(1)*bn

  END SUBROUTINE jacobianMatricesFace

  !*****************************************
  ! Set the perpendicular diffusion
  !****************************************
  SUBROUTINE setLocalDiff(xy, d_iso, d_ani)
    REAL*8, INTENT(in)  :: xy(:, :)
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
    d_iso(3, 3, :) = phys%diff_n

    d_ani(1, 1, :) = phys%diff_n
    d_ani(2, 2, :) = phys%diff_u
    d_ani(3, 3, :) = phys%diff_n

    !*****************************
    ! Non diagonal terms
    !*****************************
    ! No non-diagonal terms defined for this model

    IF (switch%difcor .GT. 0) THEN
       CALL computeIperDiffusion(xy, iperdiff)
       d_iso(1, 1, :) = d_iso(1, 1, :)*iperdiff
       d_iso(2, 2, :) = d_iso(2, 2, :)*iperdiff
       d_iso(3, 3, :) = d_iso(3, 3, :)*iperdiff
    ENDIF

  END SUBROUTINE setLocalDiff

  !*******************************************
  ! Compute local diffusion in points
  !*******************************************
  SUBROUTINE computeIperDiffusion(X, ipdiff)
    REAL*8, INTENT(IN)  :: X(:, :)
    REAL*8, INTENT(OUT) :: ipdiff(:)
    REAL*8             :: xcorn, ycorn
    REAL*8             :: rad(SIZE(X, 1))
    REAL*8             :: h
    REAL*8, PARAMETER   :: tol = 1.e-6
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
       WRITE (6, *) "Case not valid"
       STOP
    END SELECT

    !!**********************************************************
    !! Gaussian around the corner
    !!**********************************************************
    h = 10e-3
    rad = SQRT((X(:, 1)*phys%lscale - xcorn)**2 + (X(:, 2)*phys%lscale - ycorn)**2)
    ipdiff = 1 + numer%dc_coe*EXP(-(2*rad/h)**2)

  END SUBROUTINE computeIperDiffusion

  !*******************************************
  ! Compute the stabilization tensor tau
  !*******************************************
  SUBROUTINE computeTauGaussPoints(up, uc, b, bmod,n, iel, ifa, isext, xy, tau)
    REAL*8, INTENT(in)  :: up(:), uc(:), b(:), bmod, n(:), xy(:)
    REAL, INTENT(in) :: isext
    INTEGER, INTENT(in)  :: ifa, iel
    REAL*8, INTENT(out) :: tau(:, :)
    REAL*8              :: tau_aux(phys%Neq)
    REAL*8 :: xc, yc, rad, h, aux, bn, bnorm

    bn = dot_PRODUCT(b, n)
    bnorm = NORM2(b)

    IF (numer%stab == 2) THEN
       tau_aux = ABS(uc(2)*bn/uc(1))
       tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
       tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
       tau_aux(3) = tau_aux(3) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale

    ELSEIF (numer%stab == 3) THEN
       tau_aux = MAX(ABS((uc(2) + SQRT(phys%a))*bn/uc(1)), ABS((uc(2) - SQRT(phys%a))*bn/uc(1)))
       tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
       tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
       tau_aux(3) = tau_aux(3) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale

    ELSEIF (numer%stab == 4) THEN
       tau_aux = ABS((uc(2)*bn)/uc(1))
       tau_aux(1) = tau_aux(1) + phys%diff_n
       tau_aux(2) = tau_aux(2) + phys%diff_u
       tau_aux(3) = tau_aux(3) + phys%diff_n
    ELSE
       WRITE (6, *) "Wrong stabilization type: ", numer%stab
       STOP
    ENDIF
    tau(1, 1) = tau_aux(1)
    tau(2, 2) = tau_aux(2)
    tau(3, 3) = tau_aux(3)
  END SUBROUTINE computeTauGaussPoints

  SUBROUTINE computeTauGaussPoints_matrix(up, uc, b, n, xy, isext, iel, tau)
    REAL*8, INTENT(in)  :: up(:), uc(:), b(:), n(:), xy(:), isext
    REAL*8, INTENT(out) :: tau(:, :)
    INTEGER, INTENT(in) :: iel
    REAL*8              :: bn, bnorm
    REAL*8, PARAMETER :: eps = 1e-12
    REAL*8 :: U1, U2, U3, U4
    REAL*8 :: t2, t3, t4, t5, t6, t7, t8, t9
    REAL*8 :: t10, t11, t12, t13, t14, t15, t16, t17, t18, t19
    REAL*8 :: t20, t21, t22, t23, t24, t25, t26, t27, t28, t29
    REAL*8 :: t30, t31, t32, t33, t34, t35, t36, t37, t38, t39
    REAL*8 :: t40, t41, t42, t43, t44, t45, t46, t47, t48, t49
    REAL*8 :: t50, t51, t52, t53, t54, t55, t56, t57, t58, t59
    REAL*8 :: t60, t61, t62, t63, t64, t65, t66, t67, t68, t69
    REAL*8 :: t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80
    REAL*8 :: x, y, r, h, coef, r0, rc

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

  !*****************************************
  ! Jacobian matrices
  !****************************************
  !                 SUBROUTINE jacobianMatricesBound(u,v,nx,ny,Anp,Anm)
  !                 USE LinearAlgebra
  !                 real*8, intent(in)  :: u,v,nx,ny
  !                 real*8, intent(out) :: Anp(1:3,1:3),Anm(1:3,1:3)
  !   real*8              :: An(1:3,1:3),Vm(3,3),Dm(3,3),invVm(3,3)
  !
  !                        An = 0.d0
  !                        An(1,2) = nx
  !                        An(1,3) = ny
  !                        An(2,1) = phys%a - u**2.d0*nx - u*v*ny
  !                        An(2,2) = u*2.d0*nx + v*ny
  !                        An(2,3) = u*ny
  !                        An(3,1) = (phys%a-v**2)*ny - u*v*nx
  !                        An(3,2) = v*nx
  !                        An(3,3) = 2*v*ny + u*nx
  !
  !                        CALL eig(An,Vm,Dm)
  !                        CALL invert_matrix(Vm,invVm)
  !                        Anp = 0.5*(An+matmul(Vm,matmul(abs(Dm),invVm)) )
  !                        Anm = 0.5*(An-matmul(Vm,matmul(abs(Dm),invVm)))

  !                 END SUBROUTINE jacobianMatricesBound

  !***********************************************************************
  !
  !                          MAGNETIC FIELD
  !
  !***********************************************************************
#ifdef TOR3D
  SUBROUTINE defineMagneticField(x, y, t, b, divb, drift)
    REAL*8, INTENT(in)      :: x(:), y(:), t(:)
    REAL*8, INTENT(out)     :: b(:, :)
    REAL*8, INTENT(out), OPTIONAL     ::divb(:), drift(:, :)
    REAL*8                  :: xc, yc, R0, q, r
    REAL*8                  :: xmax, xmin, ymax, ymin, xm, ym, p, divbp
    REAL*8                  :: xx, yy, tt, Bp, Bt, Br, Bz, BB, dmax, B0, xr, yr
    INTEGER*4               :: i, j, ind, N2D, N1D

    N2d = SIZE(X, 1)
    N1d = SIZE(t, 1)
    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    xc = -0.5
    yc = -0.5
    ! Initialization
    b = 0.
    divbp = 0.
    IF (PRESENT(divb)) THEN
       divb = 0.
       drift = 0.
    ENDIF

    DO i = 1, N2d
       DO j = 1, N1d
          xx = x(i)
          yy = y(i)
          tt = t(j)
          ind = (j - 1)*N2d+i

          SELECT CASE (switch%testcase)
          CASE (1)
             IF (switch%axisym) THEN
                WRITE (6, *) "This is NOT an axisymmetric test case!"
                STOP
             END IF
             ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
             Br = (yy - yc)
             Bz = (-xx + xc)
             Bt = 1.

          CASE (2)
             IF (.NOT. switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             ! Axysimmetric case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
             Br = (yy - ym)/xx
             Bz = (-xx + xm)/xx
             Bt = 1.
             IF (PRESENT(divb)) THEN
                divbp = -(xm**2*ym+xx**2*ym-xm**2*yy-xx**2*yy+3*ym*yy**2-3*ym**2*yy+ym**3-yy**3-2*xm*xx*ym+2*xm*xx*yy)/(xx**4*((xm-xx)**2/xx**2+(ym-yy)**2/xx**2+1)**(1.5))

             ENDIF
          CASE (50:59)
             WRITE (6, *) "Error in defineMagneticField: you should not be here!"
             STOP
          CASE (60:69)

             ! Circular case with limiter
             R0 = geom%R0
             q = geom%q
             B0 = 2 ! not influential
             xr = xx*phys%lscale
             yr = yy*phys%lscale

             r = SQRT((xr - R0)**2 + yr**2)
             Br = -B0*yr/(xr*q*SQRT(1 - (r/R0)**2))
             Bz = B0*(xr - R0)/(xr*q*SQRT(1 - (r/R0)**2))
             Bt = B0*R0/xr

             IF (PRESENT(divb)) THEN
                IF (switch%axisym) THEN
                   divbp = -yy/xx/SQRT(R0**2*q**2 + (1 - q**2)*r**2)*phys%lscale
                ELSE
                   WRITE (6, *) "Not coded: usually here you should have an axisym simulation"
                   STOP
                END IF
             ENDIF
             IF (PRESENT(divb)) THEN
                IF (switch%driftdia) THEN
                   drift(:, 2) = -1./R0*phys%lscale
                END IF
             ENDIF
          CASE DEFAULT
             WRITE (6, *) "Error! Test case not valid"
             STOP
          END SELECT

          Bp = SQRT(Br**2 + Bz**2)
          BB = SQRT(Bp**2 + Bt**2)
          b(ind, 1) = Br/BB
          b(ind, 2) = Bz/BB
          b(ind, 3) = Bt/BB
          IF (PRESENT(divb)) THEN
             divb(ind) = divbp
          ENDIF
       END DO
    END DO
  END SUBROUTINE defineMagneticField
#else
  SUBROUTINE defineMagneticField(x, y, b, divb, drift)
    REAL*8, INTENT(in)      :: x(:), y(:)
    REAL*8, INTENT(out)     :: b(:, :)
    REAL*8, INTENT(out), OPTIONAL     ::divb(:), drift(:, :)
    REAL*8                  :: xc, yc, R0, q, r(SIZE(x))

    ! Initialization
    b = 0.
    IF (PRESENT(divb)) THEN
       divb = 0.
       drift = 0.
    ENDIF
    SELECT CASE (switch%testcase)
    CASE (1)
       ! Circular field centered in [xc, yc], n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y)
       ! Case 9 of the Matlab version: for convergence purpose
       xc = 0.
       yc = 0.
       b(:, 1) = 0.1*(y - yc)/SQRT(y**2 - 2*y*yc + yc**2 + x**2 - 2*x*xc + xc**2)
       b(:, 2) = 0.1*(-x + xc)/SQRT(y**2 - 2*y*yc + yc**2 + x**2 - 2*x*xc + xc**2)
       IF (PRESENT(divb)) THEN
          divb(:) = 0.
       ENDIF
    CASE (2)
       ! Axisymmetric case with div(b)~=0
       xc = 0.
       yc = 0.
       b(:, 1) = 1./30.*(x - y**2 + 2)
       b(:, 2) = 1./30.*(x*y + y)
       IF (PRESENT(divb)) THEN
          divb(:) = 1./30.+((1./30.)*x - (1./30.)*y**2 + 1./15.)/x + 1./30.*(x + 1)
       ENDIF
    CASE (50:59)
       WRITE (6, *) "Error in defineMagneticField: you should not be here!"
       STOP
    CASE (60:69)

       ! Circular case with limiter
       R0 = geom%R0
       q = geom%q
       r = phys%lscale*SQRT((x - R0/phys%lscale)**2 + y**2)
       b(:, 1) = -phys%lscale*y/SQRT(R0**2*q**2 + (1 - q**2)*r**2)
       b(:, 2) = phys%lscale*(x - R0/phys%lscale)/SQRT(R0**2*q**2 + (1 - q**2)*r**2)

       IF (PRESENT(divb)) THEN
          IF (switch%axisym) THEN
             divb(:) = -y/x/SQRT(R0**2*q**2 + (1 - q**2)*r**2)*phys%lscale
          ELSE
             WRITE (6, *) "Not coded: usually here you should have an axisym simulation"
             STOP
          END IF
       ENDIF

       IF (PRESENT(divb)) THEN
          IF (switch%driftdia) THEN
             drift(:, 2) = -1./R0*phys%lscale
          END IF
       ENDIF
    CASE DEFAULT
       WRITE (6, *) "Error! Test case not valid"
       STOP
    END SELECT
  END SUBROUTINE defineMagneticField
#endif

  SUBROUTINE loadMagneticField()
    USE interpolation
    USE HDF5
    USE HDF5_io_module
    INTEGER        :: i, ierr, ip, jp
    INTEGER(HID_T) :: file_id
    REAL*8, POINTER, DIMENSION(:, :) :: r2D, z2D, flux2D, Br2D, Bz2D, Bphi2D
    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: bx, by, bmod, divb, bmodx, bmody, driftx, drifty
    REAL*8, ALLOCATABLE, DIMENSION(:)   :: xvec, yvec
    REAL*8                            :: x, y

    WRITE (6, *) "******* Loading magnetic field *******"
    ! Allocate storing space in phys
    ALLOCATE (phys%b(Mesh%Nnodes, Mesh%Ndim))
    ALLOCATE (phys%divb(Mesh%Nnodes))
    ALLOCATE (phys%drift(Mesh%Nnodes, Mesh%Ndim))
    ALLOCATE (phys%flux2d(Mesh%Nnodes))

    ! Dimensions of the file storing the magnetic field for West
    ip = 541
    jp = 391
    ALLOCATE (r2D(ip, jp))
    ALLOCATE (z2D(ip, jp))
    ALLOCATE (flux2D(ip, jp))
    ALLOCATE (Br2D(ip, jp))
    ALLOCATE (Bz2D(ip, jp))
    ALLOCATE (Bphi2D(ip, jp))
    ALLOCATE (bx(ip, jp))
    ALLOCATE (by(ip, jp))
    ALLOCATE (bmod(ip, jp))
    ALLOCATE (bmodx(ip, jp))
    ALLOCATE (bmody(ip, jp))
    ALLOCATE (divb(ip, jp))
    ALLOCATE (driftx(ip, jp))
    ALLOCATE (drifty(ip, jp))

    ! Read file
    CALL HDF5_open('WEST_far_465.h5', file_id, IERR)
    CALL HDF5_array2D_reading(file_id, r2D, 'r2D')
    CALL HDF5_array2D_reading(file_id, z2D, 'z2D')
    CALL HDF5_array2D_reading(file_id, flux2D, 'flux2D')
    CALL HDF5_array2D_reading(file_id, Br2D, 'Br2D')
    CALL HDF5_array2D_reading(file_id, Bz2D, 'Bz2D')
    CALL HDF5_array2D_reading(file_id, Bphi2D, 'Bphi2D')
    CALL HDF5_close(file_id)

    ! Apply length scale
    r2D = r2D/phys%lscale
    z2D = z2D/phys%lscale

    ! Compute b
    bmod = SQRT(Br2D**2 + Bz2D**2 + Bphi2D**2)
    bx = -Br2D/bmod
    by = -Bz2D/bmod

    ! Compute divergence of b
    divb = 0.
    IF (switch%axisym) THEN
       ! 1/r*(d r*br/dr)+dbz/dz
       divb(2:ip - 1, 2:jp - 1) = 1./r2D(2:ip - 1, 2:jp - 1)*(r2D(2:ip - 1, 3:jp)*bx(2:ip - 1, 3:jp) - &
            r2D(2:ip - 1, 1:jp - 2)*bx(2:ip - 1, 1:jp - 2))/(r2D(2:ip - 1, 3:jp) - &
            r2D(2:ip - 1, 1:jp - 2)) + (by(3:ip, 2:jp - 1) - by(1:ip - 2, 2:jp - 1))/(z2D(3:ip, 2:jp - 1) - z2D(1:ip - 2, 2:jp - 1))

    ELSE
       ! dbr/dr+dbz/dz
       divb(2:ip - 1, 2:jp - 1) = (bx(2:ip - 1, 3:jp - 1) - bx(2:ip - 1, 1:jp - 2))/(r2D(2:ip - 1, 3:jp) - r2D(2:ip - 1, 1:jp - 2)) + &
            (by(3:ip, 2:jp - 1) - by(1:ip - 2, 2:jp - 1))/(z2D(3:ip, 2:jp - 1) - z2D(1:ip - 2, 2:jp - 1))
    END IF

    ! Compute drift velocity
    driftx = 0.
    drifty = 0.
    IF (switch%driftdia) THEN
       bmodx = (bmod(2:ip - 1, 3:jp) - bmod(2:ip - 1, 1:jp - 2))/(r2D(2:ip - 1, 3:jp) - r2D(2:ip - 1, 1:jp - 2))
       bmody = (bmod(3:ip, 2:jp - 1) - bmod(1:ip - 2, 2:jp - 1))/(z2D(3:ip, 2:jp - 1) - z2D(1:ip - 2, 2:jp - 1))
       driftx(2:ip - 1, 2:jp - 1) = -Bphi2D(2:ip - 1, 2:jp - 1)*bmody/bmod(2:ip - 1, 2:jp - 1)**3
       drifty(2:ip - 1, 2:jp - 1) = Bphi2D(2:ip - 1, 2:jp - 1)*bmodx/bmod(2:ip - 1, 2:jp - 1)**3
    END IF

    ! Interpolate
    ALLOCATE (xvec(jp))
    ALLOCATE (yvec(ip))
    xvec = r2D(1, :)
    yvec = z2D(:, 1)
    DO i = 1, Mesh%Nnodes
       x = Mesh%X(i, 1)
       y = Mesh%X(i, 2)
       phys%b(i, 1) = interpolate(ip, yvec, jp, xvec, bx, y, x, 1e-12)
       phys%b(i, 2) = interpolate(ip, yvec, jp, xvec, by, y, x, 1e-12)
       phys%divb(i) = interpolate(ip, yvec, jp, xvec, divb, y, x, 1e-12)
       phys%drift(i, 1) = interpolate(ip, yvec, jp, xvec, driftx, y, x, 1e-12)
       phys%drift(i, 2) = interpolate(ip, yvec, jp, xvec, drifty, y, x, 1e-12)
       phys%flux2D(i) = interpolate(ip, yvec, jp, xvec, flux2D, y, x, 1e-12)
    END DO

    ! Free memory
    DEALLOCATE (Br2D, Bz2D, Bphi2D, xvec, yvec)
    DEALLOCATE (r2D, z2D, flux2D, bx, by, bmod, bmodx, bmody, divb, driftx, drifty)

  END SUBROUTINE loadMagneticField

  SUBROUTINE loadMagneticFieldTemporalEvolution()
    USE HDF5
    USE HDF5_io_module
    USE MPI_OMP
    INTEGER        :: i, ierr, k
    CHARACTER(LEN=20) :: fname = 'Evolving_equilibrium'
    CHARACTER(10)  :: npr, nid, nit
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id

    WRITE (6, *) "******* Loading magnetic field *******"

    ! Allocate storing space in phys
    ALLOCATE (phys%Br(Mesh%Nnodes))
    ALLOCATE (phys%Bz(Mesh%Nnodes))
    ALLOCATE (phys%Bt(Mesh%Nnodes))
    ALLOCATE (phys%flux2d(Mesh%Nnodes))

    ! File name
    WRITE (nit, "(i10)") time%it
    nit = TRIM(ADJUSTL(nit))
    k = INDEX(nit, " ") - 1

    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))//'.h5'
    END IF

    WRITE (6, *) 'Magnetic field loaded from file: ', TRIM(ADJUSTL(fname_complete))

    ! Read file
    CALL HDF5_open(fname_complete, file_id, IERR)
    CALL HDF5_array1D_reading(file_id, phys%Br, 'Br')
    CALL HDF5_array1D_reading(file_id, phys%Bz, 'Bz')
    CALL HDF5_array1D_reading(file_id, phys%Bt, 'Bt')
    CALL HDF5_array1D_reading(file_id, phys%flux2D, 'flux')
    CALL HDF5_close(file_id)

  END SUBROUTINE loadMagneticFieldTemporalEvolution

END MODULE physics
