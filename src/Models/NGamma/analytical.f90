!*****************************************
! project: MHDG
! file: analytical.f90
! date: 15/02/2017
!  ******** N-Gamma system Isothermal ****
! Define the analytical solution and the
! body forces for different cases
!*****************************************
MODULE analytical
  USE prec_const
  USE globals, ONLY: switch
  USE physics
  USE PrintUtils
  IMPLICIT NONE

CONTAINS

#ifdef TOR3D

  !***********************************************************************
  !
  !                            VERSION 3D TOROIDAL
  !
  !***********************************************************************

  !*****************************************
  ! Analytical solution
  !****************************************
  SUBROUTINE analytical_solution(x, y, t, u)
    REAL*8, DIMENSION(:), INTENT(IN)        :: x, y, t
    REAL*8, DIMENSION(:, :), INTENT(OUT)     :: u
    REAL*8, DIMENSION(SIZE(u, 1), phys%npv)  :: up
    INTEGER:: i, j, ind, N2D, N1D
    REAL*8 :: a, b, r, xx, yy, tt, xmax, xmin, ymax, ymin, xm, ym
    REAL*8 :: aux

    u = 0.
    a = 2*pi
    N2D = SIZE(x, 1)
    N1D = SIZE(t, 1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    DO i = 1, N2d
       DO j = 1, N1d  ! TODO: check if I need to invert i and j
          xx = x(i)
          yy = y(i)
          tt = t(j)
          r = SQRT((xx - xm)**2 + (yy - ym)**2)
          ind = (j - 1)*N2d+i
          SELECT CASE (switch%testcase)
          CASE (1)
             IF (switch%axisym) THEN
                WRITE (6, *) "This is NOT an axisymmetric test case!"
                STOP
             END IF
             ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
             up(ind, 1) = 2 + SIN(a*xx)*SIN(a*yy)
             up(ind, 2) = COS(a*xx)*COS(a*yy)
          CASE (2)
             ! Axisimmetric case with div(b)~=0
             IF (.NOT. switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             up(ind, 1) = 2 + SIN(a*xx)*SIN(a*yy)
             up(ind, 2) = COS(a*xx)*COS(a*yy)
          CASE (50:64)
             up(ind, 1) = 1.
             up(ind, 2) = 0.
          CASE (65)
             up(ind, 1) = 1.
             up(ind, 2) = 0.
             r = SQRT((xx*phys%lscale - geom%R0)**2 + (yy*phys%lscale - 0.75)**2)
             IF (r .LE. 0.05) THEN
                up(ind, 2) = 1.
             END IF
          CASE (80:89)
             up(ind, 1) = 1.
             up(ind, 2) = 0.
          CASE DEFAULT
             WRITE (6, *) "Error! Test case not valid"
             STOP
          END SELECT
       END DO
    END DO
    ! Convert physical variables to conservative variables
    CALL phys2cons(up, u)
  END SUBROUTINE analytical_solution

  !*****************************************
  ! Analytical gradient
  !****************************************
  SUBROUTINE analytical_gradient(x, y, t, u, ux, uy, ut)
    REAL*8, DIMENSION(:), INTENT(IN)        :: x, y, t
    REAL*8, DIMENSION(:, :), INTENT(IN)      :: u
    REAL*8, DIMENSION(:, :), INTENT(OUT)     :: ux, uy, ut
    REAL*8, DIMENSION(SIZE(u, 1), phys%npv)  :: upx, upy, upt
    REAL*8, DIMENSION(SIZE(u, 1), phys%npv)   :: up
    INTEGER:: i, j, ind, N1D, N2D
    REAL*8 :: a, b, r, xx, yy, tt, xmax, xmin, ymax, ymin, xm, ym
    REAL*8 :: aux

    N2D = SIZE(x, 1)
    N1D = SIZE(t, 1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)
    upx = 0.
    upy = 0.
    upt = 0.
    CALL cons2phys(u, up)
    a = 2*pi

    DO i = 1, N2d
       DO j = 1, N1d
          xx = x(i)
          yy = y(i)
          tt = t(j)
          r = SQRT((xx - xm)**2 + (yy - ym)**2)
          ind = (j - 1)*N2d+i
          SELECT CASE (switch%testcase)
          CASE (1)
             IF (switch%axisym) THEN
                WRITE (6, *) "This is NOT an axisymmetric test case!"
                STOP
             END IF
             ! Circular field centered in [xc, yc], n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
             upx(ind, 1) = a*COS(a*xx)*SIN(a*yy)
             upx(ind, 2) = -a*SIN(a*xx)*COS(a*yy)

             upy(ind, 1) = a*SIN(a*xx)*COS(a*yy)
             upy(ind, 2) = -a*COS(a*xx)*SIN(a*yy)

          CASE (2)
             ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
             IF (.NOT. switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             upx(ind, 1) = a*COS(a*xx)*SIN(a*yy)
             upx(ind, 2) = -a*SIN(a*xx)*COS(a*yy)

             upy(ind, 1) = a*SIN(a*xx)*COS(a*yy)
             upy(ind, 2) = -a*COS(a*xx)*SIN(a*yy)
          CASE (3)
             ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
             IF (.NOT. switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             upt(ind, 1) = +a*COS(a*tt)
             upt(ind, 2) = -a*SIN(a*tt)
          CASE (5)
             ! Do nothing
          CASE (6)
             ! Do nothing
          CASE (50:)
             ! Do nothing
          CASE DEFAULT
             WRITE (6, *) "Error! Test case not valid"
             STOP
          END SELECT
       END DO
    END DO

    ! Convert physical variables to conservative variables
    ux(:, 1) = upx(:, 1)
    uy(:, 1) = upy(:, 1)
    ut(:, 1) = upt(:, 1)
    ux(:, 2) = (upx(:, 1)*up(:, 2) + up(:, 1)*upx(:, 2))
    uy(:, 2) = (upy(:, 1)*up(:, 2) + up(:, 1)*upy(:, 2))
    ut(:, 2) = (upt(:, 1)*up(:, 2) + up(:, 1)*upt(:, 2))

  END SUBROUTINE analytical_gradient

  !*****************************************
  ! Body forces
  !****************************************
  SUBROUTINE body_force(x, y, t, f)
    REAL*8, DIMENSION(:), INTENT(IN) :: x, y, t
    REAL*8, DIMENSION(:, :), INTENT(OUT) :: f
    REAL*8  :: xc, yc, D, mu, k
    INTEGER:: i, j, N2D, N1D, ind
    REAL*8 :: a, b, r, r2, xx, yy, tt, xmax, xmin, ymax, ymin, xm, ym
    REAL*8  :: t1, t2, t3, t4, t5, t6, t7, t8, t9, ta, tb, cx, cy, sx, sy, cx2, cy2
    f = 0.
    a = 2*pi
    N2D = SIZE(x, 1)
    N1D = SIZE(t, 1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    f = 0.
    DO i = 1, N2d
       DO j = 1, N1d
          xx = x(i)
          yy = y(i)
          tt = t(j)
          r = SQRT((xx - xm)**2 + (yy - ym)**2)
          r2 = ((xm - xx)**2 + (ym - yy)**2)**(1.5)
          ind = (j - 1)*N2d+i
          SELECT CASE (switch%testcase)
          CASE (1)
             ! Circular field centered in [xc, yc], n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y)
             IF (switch%axisym) THEN
                WRITE (6, *) "This is NOT an axisymmetric test case!"
                STOP
             END IF
             a = 2*pi
             b = 2*pi
             !                                                                                                                xc = 0.
             !                                                                                                                yc = 0.
             !                                                                                                                xm = -0.5
             !                                                                                                                ym = -0.5
             D = phys%diff_n
             mu = phys%diff_u
             k = phys%a

             t1 = COS(a*xx)*COS(a*yy)
             t2 = SIN(a*xx)*SIN(a*yy)
             t3 = COS(a*xx)*SIN(a*yy)
             t4 = COS(a*yy)*SIN(a*xx)
             t5 = COS(a*yy)*SIN(a*yy)
             t8 = COS(a*xx)**2*COS(a*yy)**2
             ta = (xm**2 - 2*ym*yy - 2*xm*xx + xx**2 + ym**2 + yy**2 + 1)

             f(ind,1) = (2*D*a**2*SIN(a*xx)*SIN(a*yy)-2*a*xm*COS(a*xx)*SIN(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+2*a*xx*COS(a*xx)*SIN(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+2*a*ym*COS(a*yy)*SIN(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-2*a*yy*COS(a*yy)*SIN(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+D*a**2*xm**2*SIN(a*xx)*SIN(a*yy)+D*a**2*xx**2*SIN(a*xx)*SIN(a*yy)+D*a**2*ym**2*SIN(a*xx)*SIN(a*yy)+D*a**2*yy**2*SIN(a*xx)*SIN(a*yy)+D*a*xm*COS(a*xx)*SIN(a*yy)-D*a*xx*COS(a*xx)*SIN(a*yy)+D*a*ym*COS(a*yy)*SIN(a*xx)-D*a*yy*COS(a*yy)*SIN(a*xx)+a*xm*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*xx*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*ym*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*yy*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*xm*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*xx*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*ym*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*yy*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-2*D*a**2*xm*ym*COS(a*xx)*COS(a*yy)+2*D*a**2*xx*ym*COS(a*xx)*COS(a*yy)+2*D*a**2*xm*yy*COS(a*xx)*COS(a*yy)-2*D*a**2*xx*yy*COS(a*xx)*COS(a*yy)-2*D*a**2*xm*xx*SIN(a*xx)*SIN(a*yy)-2*D*a**2*ym*yy*SIN(a*xx)*SIN(a*yy))/((xm-xx)**2+(ym-yy)**2+1);
             f(ind,2) = ((a*mu*xx*SIN(2*a*yy))/2-(a*mu*xm*SIN(2*a*yy))/2-(a*mu*ym*SIN(2*a*xx))/2+(a*mu*yy*SIN(2*a*xx))/2-2*a**2*mu*xm*ym+2*a**2*mu*xx*ym+2*a**2*mu*xm*yy-2*a**2*mu*xx*yy+4*a**2*mu*COS(a*xx)*COS(a*yy)+4*a**2*mu*xm*ym*COS(a*xx)**2-4*a**2*mu*xx*ym*COS(a*xx)**2-4*a**2*mu*xm*yy*COS(a*xx)**2+4*a**2*mu*xx*yy*COS(a*xx)**2+4*a**2*mu*xm*ym*COS(a*yy)**2-4*a**2*mu*xx*ym*COS(a*yy)**2-4*a**2*mu*xm*yy*COS(a*yy)**2+4*a**2*mu*xx*yy*COS(a*yy)**2+2*a**2*mu*xm**2*COS(a*xx)*COS(a*yy)+2*a**2*mu*xx**2*COS(a*xx)*COS(a*yy)+2*a**2*mu*ym**2*COS(a*xx)*COS(a*yy)+2*a**2*mu*yy**2*COS(a*xx)*COS(a*yy)-2*a*mu*xm*COS(a*yy)*SIN(a*xx)+2*a*mu*xx*COS(a*yy)*SIN(a*xx)-2*a*mu*ym*COS(a*xx)*SIN(a*yy)+2*a*mu*yy*COS(a*xx)*SIN(a*yy)+8*a**2*mu*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)-2*a*xm*COS(a*xx)**2*COS(a*yy)*SIN(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+2*a*xx*COS(a*xx)**2*COS(a*yy)*SIN(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+4*a*ym*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-4*a*xm*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+4*a*xx*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-4*a*yy*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+2*a*ym*COS(a*xx)*COS(a*yy)**2*SIN(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-2*a*yy*COS(a*xx)*COS(a*yy)**2*SIN(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-8*a**2*mu*xm*ym*COS(a*xx)**2*COS(a*yy)**2+8*a**2*mu*xx*ym*COS(a*xx)**2*COS(a*yy)**2+8*a**2*mu*xm*yy*COS(a*xx)**2*COS(a*yy)**2-8*a**2*mu*xx*yy*COS(a*xx)**2*COS(a*yy)**2+3*a*xm*COS(a*xx)**2*COS(a*yy)**3*SIN(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-3*a*xx*COS(a*xx)**2*COS(a*yy)**3*SIN(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-3*a*ym*COS(a*xx)**3*COS(a*yy)**2*SIN(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+3*a*yy*COS(a*xx)**3*COS(a*yy)**2*SIN(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+a*k*xm*COS(a*yy)*SIN(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-a*k*xx*COS(a*yy)*SIN(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-a*k*ym*COS(a*xx)*SIN(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+a*k*yy*COS(a*xx)*SIN(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-4*a**2*mu*xm*xx*COS(a*xx)*COS(a*yy)-4*a**2*mu*ym*yy*COS(a*xx)*COS(a*yy)-4*a**2*mu*xm*ym*SIN(a*xx)*SIN(a*yy)+4*a**2*mu*xx*ym*SIN(a*xx)*SIN(a*yy)+4*a**2*mu*xm*yy*SIN(a*xx)*SIN(a*yy)-4*a**2*mu*xx*yy*SIN(a*xx)*SIN(a*yy)+2*a*mu*ym*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)+2*a*mu*xm*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)-2*a*mu*xx*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)-2*a*mu*yy*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)+4*a**2*mu*xm**2*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)+4*a**2*mu*xx**2*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)+4*a**2*mu*ym**2*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)+4*a**2*mu*yy**2*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)-8*a**2*mu*xm*xx*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)-8*a**2*mu*ym*yy*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy))/(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1);
             !                                                                                                                f(ind,1) = (2*D*a**2*t2-2*a*xm*t3*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+2*a*xx*t3*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+2*a*ym*t4*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-2*a*yy*t4*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+D*a**2*xm**2*t2+D*a**2*xx**2*t2+D*a**2*ym**2*t2+D*a**2*yy**2*t2+D*a*xm*t3-D*a*xx*t3+D*a*ym*t4-D*a*yy*t4+a*xm*t1**2*sin(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*xx*t1**2*sin(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*ym*cos(a*xx)**2*t5*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*yy*cos(a*xx)**2*t5*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*xm*cos(a*xx)*t2**2*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*xx*cos(a*xx)*t2**2*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*ym*t4**2*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*yy*t4**2*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-2*D*a**2*xm*ym*t1+2*D*a**2*xx*ym*t1+2*D*a**2*xm*yy*t1-2*D*a**2*xx*yy*t1-2*D*a**2*xm*xx*t2-2*D*a**2*ym*yy*t2)/((xm-xx)**2+(ym-yy)**2+1)
             !                                                                                                                f(ind,2) = ((a*mu*xx*sin(2*a*yy))/2-(a*mu*xm*sin(2*a*yy))/2-(a*mu*ym*sin(2*a*xx))/2+(a*mu*yy*sin(2*a*xx))/2-2*a**2*mu*xm*ym+2*a**2*mu*xx*ym+2*a**2*mu*xm*yy-2*a**2*mu*xx*yy+4*a**2*mu*t1+4*a**2*mu*xm*ym*cos(a*xx)**2-4*a**2*mu*xx*ym*cos(a*xx)**2-4*a**2*mu*xm*yy*cos(a*xx)**2+4*a**2*mu*xx*yy*cos(a*xx)**2+4*a**2*mu*xm*ym*cos(a*yy)**2-4*a**2*mu*xx*ym*cos(a*yy)**2-4*a**2*mu*xm*yy*cos(a*yy)**2+4*a**2*mu*xx*yy*cos(a*yy)**2+2*a**2*mu*xm**2*t1+2*a**2*mu*xx**2*t1+2*a**2*mu*ym**2*t1+2*a**2*mu*yy**2*t1-2*a*mu*xm*t4+2*a*mu*xx*t4-2*a*mu*ym*t3+2*a*mu*yy*t3+8*a**2*mu*t1*t2-2*a*xm*cos(a*xx)**2*t4*ta**(0.5)+2*a*xx*cos(a*xx)**2*t4*ta**(0.5)+4*a*ym*t1**2*sin(a*xx)*ta**(0.5)-4*a*xm*cos(a*xx)**2*t5*ta**(0.5)+4*a*xx*cos(a*xx)**2*t5*ta**(0.5)-4*a*yy*t1**2*sin(a*xx)*ta**(0.5)+2*a*ym*t1**2*sin(a*yy)*ta**(0.5)-2*a*yy*t1**2*sin(a*yy)*ta**(0.5)-8*a**2*mu*xm*ym*t8+8*a**2*mu*xx*ym*t8+8*a**2*mu*xm*yy*t8-8*a**2*mu*xx*yy*t8+3*a*xm*cos(a*xx)**2*cos(a*yy)**3*sin(a*xx)*ta**(0.5)-3*a*xx*cos(a*xx)**2*cos(a*yy)**3*sin(a*xx)*ta**(0.5)-3*a*ym*cos(a*xx)**3*cos(a*yy)**2*sin(a*yy)*ta**(0.5)+3*a*yy*cos(a*xx)**3*cos(a*yy)**2*sin(a*yy)*ta**(0.5)+a*k*xm*t4*ta**(0.5)-a*k*xx*t4*ta**(0.5)-a*k*ym*t3*ta**(0.5)+a*k*yy*t3*ta**(0.5)-4*a**2*mu*xm*xx*t1-4*a**2*mu*ym*yy*t1-4*a**2*mu*xm*ym*t2+4*a**2*mu*xx*ym*t2+4*a**2*mu*xm*yy*t2-4*a**2*mu*xx*yy*t2+2*a*mu*ym*t1**2*sin(a*xx)+2*a*mu*xm*cos(a*xx)**2*t5-2*a*mu*xx*cos(a*xx)**2*t5-2*a*mu*yy*t1**2*sin(a*xx)+4*a**2*mu*xm**2*t1*t2+4*a**2*mu*xx**2*t1*t2+4*a**2*mu*ym**2*t1*t2+4*a**2*mu*yy**2*t1*t2-8*a**2*mu*xm*xx*t1*t2-8*a**2*mu*ym*yy*t1*t2)/ta
          CASE (2)
             ! Axisimmetric case with div(b)~=0
             IF (.NOT. switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             a = 2*pi
             b = 2*pi
             !                                                                                                                xc = 0.
             !                                                                                                                yc = 0.
             D = phys%diff_n
             mu = phys%diff_u
             k = phys%a

             cx = COS(a*xx)
             cy = COS(a*yy)
             sx = SIN(a*xx)
             sy = SIN(a*yy)
             cx2 = cx**2
             cy2 = cy**2

             t1 = COS(a*xx)*COS(a*yy)
             t2 = SIN(a*xx)*SIN(a*yy)
             t3 = COS(a*xx)*SIN(a*yy)
             t4 = COS(a*yy)*SIN(a*xx)
             t5 = ((xm**2 - 2*ym*yy - 2*xm*xx + 2*xx**2 + ym**2 + yy**2)/xx**2)**(0.5)
             t6 = ((xm**2 - 2*ym*yy - 2*xm*xx + 2*xx**2 + ym**2 + yy**2)/xx**2)**(1.5)

             f(ind, 1) = (1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy&
                  &**2)**2*(-D*a*t3*xm**4 - D*a*t3*xx**4*6.0D0 + a*t3*xx**5*(1.0D0/xx**2*&
                  &(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0)*2.0D0 - t1*x&
                  &x*ym**3*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*&
                  &2.0D0 - t1*xx**3*ym*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)*&
                  &*2 + 1.0D0)*2.0D0 + t1*xx*yy**3*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**&
                  &2*(ym - yy)**2 + 1.0D0)*2.0D0 + t1*xx**3*yy*SQRT(1.0D0/xx**2*(xm - xx)**2 +&
                  &1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*2.0D0 + D*a**2*t2*xx**5*6.0D0 + a*cx*t2*&
                  &*2*xx**5*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3&
                  &.0D0/2.0D0) - D*a*t3*xm**2*xx**2*1.1D+1 - D*a**2*t2*xm*xx**4*1.0D+1 + D*&
                  &a**2*t2*xm**4*xx - D*a*t3*xm**2*ym**2 - D*a*t3*xx**2*ym**2*5.0D0 + D*a**&
                  &2*t1*xx**4*ym*4.0D0 + D*a**2*t2*xx*ym**4 - D*a*t3*xm**2*yy**2 - D*a*t3*x&
                  &x**2*yy**2*5.0D0 - D*a**2*t1*xx**4*yy*4.0D0 + D*a**2*t2*xx*yy**4 - a*sx*&
                  &t1**2*xx**5*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*&
                  &*(3.0D0/2.0D0) + D*a**2*t2*xm**2*xx**3*9.0D0 - D*a**2*t2*xm**3*xx**2*4&
                  &.0D0 + D*a**2*t1*xx**2*ym**3*2.0D0 + D*a**2*t2*xx**3*ym**2*5.0D0 - D*a**&
                  &2*t1*xx**2*yy**3*2.0D0 + D*a**2*t2*xx**3*yy**2*5.0D0 + D*a*t3*xm*xx**3&
                  &*1.2D+1 + D*a*t3*xm**3*xx*5.0D0 - D*a*t4*xm*ym**3 - D*a*t4*xm**3*ym + D*a*&
                  &t4*xx*ym**3*2.0D0 + D*a*t4*xx**3*ym*2.0D0 + D*a*t4*xm*yy**3 + D*a*t4*xm*&
                  &*3*yy - D*a*t4*xx*yy**3*2.0D0 - D*a*t4*xx**3*yy*2.0D0 - a*t3*xm*xx**4*(1&
                  &.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0)&
                  &*2.0D0 + a*t4*xx**4*ym*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**&
                  &2 + 1.0D0)**(3.0D0/2.0D0)*2.0D0 - a*t4*xx**4*yy*(1.0D0/xx**2*(xm - xx)**&
                  &2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0)*2.0D0 - t1*t2*xx*ym**&
                  &3*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0) - t1*t2*&
                  &xx**3*ym*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)&
                  &+ t1*t2*xx*yy**3*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2&
                  &+ 1.0D0) + t1*t2*xx**3*yy*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym&
                  &- yy)**2 + 1.0D0) + t1*xm*xx**2*ym*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx&
                  &**2*(ym - yy)**2 + 1.0D0)*4.0D0 - t1*xm**2*xx*ym*SQRT(1.0D0/xx**2*(xm - xx&
                  &)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*2.0D0 - t1*xm*xx**2*yy*SQRT(1.0D0&
                  &/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*4.0D0 + t1*xm**2*xx*&
                  &yy*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*2.0D0&
                  &- t1*xx*ym*yy**2*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2&
                  &+ 1.0D0)*6.0D0 + t1*xx*ym**2*yy*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx*&
                  &*2*(ym - yy)**2 + 1.0D0)*6.0D0 + D*a**2*t1*xm**2*xx**2*ym*6.0D0 - D*a**2*t&
                  &2*xm*xx**2*ym**2*4.0D0 + D*a**2*t2*xm**2*xx*ym**2*2.0D0 - D*a**2*t1*xm&
                  &**2*xx**2*yy*6.0D0 - D*a**2*t2*xm*xx**2*yy**2*4.0D0 + D*a**2*t2*xm**2*&
                  &xx*yy**2*2.0D0 + D*a**2*t1*xx**2*ym*yy**2*6.0D0 - D*a**2*t1*xx**2*ym**&
                  &2*yy*6.0D0 + D*a**2*t2*xx*ym**2*yy**2*6.0D0 + D*a*t3*xm*xx*ym**2*3.0D0&
                  &- D*a*t4*xm*xx**2*ym*4.0D0 + D*a*t4*xm**2*xx*ym*4.0D0 + D*a*t3*xm*xx*yy&
                  &**2*3.0D0 + D*a*t4*xm*xx**2*yy*4.0D0 - D*a*t4*xm**2*xx*yy*4.0D0 + D*a*t3&
                  &*xm**2*ym*yy*2.0D0 - D*a*t4*xm*ym*yy**2*3.0D0 + D*a*t4*xm*ym**2*yy*3.0&
                  &D0 + D*a*t3*xx**2*ym*yy*1.0D+1 + D*a*t4*xx*ym*yy**2*6.0D0 - D*a*t4*xx*ym&
                  &**2*yy*6.0D0 + t1*t2*xm*xx**2*ym*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/x&
                  &x**2*(ym - yy)**2 + 1.0D0)*2.0D0 - t1*t2*xm**2*xx*ym*SQRT(1.0D0/xx**2*(x&
                  &m - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0) - t1*t2*xm*xx**2*yy*SQRT(1.0D&
                  &0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*2.0D0 + t1*t2*xm**2&
                  &*xx*yy*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0) - t&
                  &1*t2*xx*ym*yy**2*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**&
                  &2 + 1.0D0)*3.0D0 + t1*t2*xx*ym**2*yy*SQRT(1.0D0/xx**2*(xm - xx)**2 + 1.0D0&
                  &/xx**2*(ym - yy)**2 + 1.0D0)*3.0D0 - a*cx*t2**2*xm*xx**4*(1.0D0/xx**2*(x&
                  &m - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0) - D*a**2*t1*xm&
                  &*xx*ym**3*2.0D0 - D*a**2*t1*xm*xx**3*ym*8.0D0 - D*a**2*t1*xm**3*xx*ym*&
                  &2.0D0 + D*a**2*t1*xm*xx*yy**3*2.0D0 + D*a**2*t1*xm*xx**3*yy*8.0D0 + D*a*&
                  &*2*t1*xm**3*xx*yy*2.0D0 - D*a**2*t2*xx*ym*yy**3*4.0D0 - D*a**2*t2*xx*y&
                  &m**3*yy*4.0D0 - D*a**2*t2*xx**3*ym*yy*1.0D+1 + a*sx*t1**2*xm*xx**4*(1.&
                  &0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0) +&
                  &a*sy*t4**2*xx**4*ym*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2&
                  &+ 1.0D0)**(3.0D0/2.0D0) - a*sy*t4**2*xx**4*yy*(1.0D0/xx**2*(xm - xx)**2&
                  &+ 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0) - D*a**2*t1*xm*xx*ym*y&
                  &y**2*6.0D0 + D*a**2*t1*xm*xx*ym**2*yy*6.0D0 + D*a**2*t2*xm*xx**2*ym*yy&
                  &*8.0D0 - D*a**2*t2*xm**2*xx*ym*yy*4.0D0 - D*a*t3*xm*xx*ym*yy*6.0D0 - a*c&
                  &x2*cy*sy*xx**4*ym*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1&
                  &.0D0)**(3.0D0/2.0D0) + a*cx2*cy*sy*xx**4*yy*(1.0D0/xx**2*(xm - xx)**2 +&
                  &1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0)))/xx

             f(ind, 2) = (1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy&
                  &**2)**2*(a**2*mu*xx**2*ym**3*2.0D0 - a**2*mu*xx**2*yy**3*2.0D0 + a*mu*&
                  &t4*xm**4*2.0D0 + a*mu*t4*xx**4*1.2D+1 - a*t4*t6*xx**5*1.0D+1 + (a*mu*xm*&
                  &*4*SIN(a*yy*2.0D0))/2.0D0 + a*mu*xx**4*SIN(a*yy*2.0D0)*3.0D0 + a**2*mu&
                  &*t1*xx**5*1.2D+1 + a**2*mu*xx**4*ym*4.0D0 - a**2*mu*xx**4*yy*4.0D0 + (a*&
                  &mu*xm*ym**3*SIN(a*xx*2.0D0))/2.0D0 + (a*mu*xm**3*ym*SIN(a*xx*2.0D0))&
                  &/2.0D0 - a*mu*xx*ym**3*SIN(a*xx*2.0D0) - a*mu*xx**3*ym*SIN(a*xx*2.0D0)&
                  &- a*mu*xm*xx**3*SIN(a*yy*2.0D0)*6.0D0 - (a*mu*xm*yy**3*SIN(a*xx*2.0D0&
                  &))/2.0D0 - a*mu*xm**3*xx*SIN(a*yy*2.0D0)*(5.0D0/2.0D0) - (a*mu*xm**3*y&
                  &y*SIN(a*xx*2.0D0))/2.0D0 + a*mu*xx*yy**3*SIN(a*xx*2.0D0) + a*mu*xx**3*&
                  &yy*SIN(a*xx*2.0D0) + a**2*mu*xx**2*yy**3*COS(a*xx)**2*4.0D0 - a**2*cx2&
                  &*mu*xx**4*ym*8.0D0 - a**2*cy2*mu*xx**4*ym*8.0D0 + a**2*cx2*mu*xx**4*yy&
                  &*8.0D0 + a**2*cy2*mu*xx**4*yy*8.0D0 + a**2*mu*t1*t2*xx**5*2.4D+1 - a**2*&
                  &mu*t1*xm*xx**4*2.0D+1 + a**2*mu*t1*xm**4*xx*2.0D0 + a*mu*t4*xm**2*xx**&
                  &2*2.2D+1 + a*mu*t4*xm**2*ym**2*2.0D0 + a**2*mu*t1*xx*ym**4*2.0D0 + a*mu*&
                  &t4*xx**2*ym**2*1.0D+1 + a**2*mu*t2*xx**4*ym*8.0D0 + a*mu*t4*xm**2*yy**&
                  &2*2.0D0 + a**2*mu*t1*xx*yy**4*2.0D0 + a*mu*t4*xx**2*yy**2*1.0D+1 - a**2*&
                  &mu*t2*xx**4*yy*8.0D0 - a**2*mu*xm*xx*ym**3*2.0D0 - a**2*mu*xm*xx**3*ym&
                  &*8.0D0 - a**2*mu*xm**3*xx*ym*2.0D0 + a**2*mu*xm*xx*yy**3*2.0D0 + a**2*mu&
                  &*xm*xx**3*yy*8.0D0 + a**2*mu*xm**3*xx*yy*2.0D0 + cx2*t5*xx*yy**3*COS(a&
                  &*yy)**2*2.0D0 + a*mu*xm**2*xx**2*SIN(a*yy*2.0D0)*(1.1D+1/2.0D0) + (a*m&
                  &u*xm**2*ym**2*SIN(a*yy*2.0D0))/2.0D0 + a*mu*xx**2*ym**2*SIN(a*yy*2.0&
                  &D0)*(5.0D0/2.0D0) + (a*mu*xm**2*yy**2*SIN(a*yy*2.0D0))/2.0D0 + a*mu*xx&
                  &**2*yy**2*SIN(a*yy*2.0D0)*(5.0D0/2.0D0) - a**2*cx2*mu*xx**2*ym**3*4.&
                  &0D0 - a**2*cy2*mu*xx**2*ym**3*4.0D0 + a**2*cy2*mu*xx**2*yy**3*4.0D0 + a*&
                  &*2*mu*t1*xm**2*xx**3*1.8D+1 - a**2*mu*t1*xm**3*xx**2*8.0D0 + a**2*mu*t&
                  &1*xx**3*ym**2*1.0D+1 + a**2*mu*t2*xx**2*ym**3*4.0D0 + a**2*mu*t1*xx**3&
                  &*yy**2*1.0D+1 - a**2*mu*t2*xx**2*yy**3*4.0D0 + a**2*mu*xm**2*xx**2*ym*&
                  &6.0D0 - a**2*mu*xm**2*xx**2*yy*6.0D0 + a**2*mu*xx**2*ym*yy**2*6.0D0 - a*&
                  &*2*mu*xx**2*ym**2*yy*6.0D0 + a*cx2*t4*xx**5*(1.0D0/xx**2*(xm*xx*(-2.&
                  &0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*2.&
                  &0D0 - cx2*cy2*xx**3*ym*SQRT(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 +&
                  &xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*2.0D0 + cx2*cy2*xx**3*yy*SQRT(1.0D0/&
                  &xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*&
                  &2.0D0 - a*t3*xx**4*ym*(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2&
                  &+ xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*1.0D+1 + a*t3*xx**4*yy*(1.&
                  &0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**&
                  &2))**(3.0D0/2.0D0)*1.0D+1 - cx2*cy2*t5*xx*ym**3*2.0D0 - a*mu*t4*xm*xx*&
                  &*3*2.4D+1 - a*mu*t4*xm**3*xx*1.0D+1 + a*mu*t3*xm*ym**3*2.0D0 + a*mu*t3*x&
                  &m**3*ym*2.0D0 - a*mu*t3*xx*ym**3*4.0D0 - a*mu*t3*xx**3*ym*4.0D0 - a*mu*t&
                  &3*xm*yy**3*2.0D0 - a*mu*t3*xm**3*yy*2.0D0 + a*mu*t3*xx*yy**3*4.0D0 + a*m&
                  &u*t3*xx**3*yy*4.0D0 + a*t4*t6*xm*xx**4*1.0D+1 + a**2*cx**2*cy2*mu*xx**&
                  &4*ym*1.6D+1 + a**2*cx2*cy2*mu*xx**2*ym**3*8.0D0 - a**2*cx2*cy2*mu*xx**&
                  &2*yy**3*8.0D0 - a*t1**2*t6*xx**4*yy*SIN(a*xx)*4.0D0 - a**2*cx2*mu*xm**&
                  &2*xx**2*ym*1.2D+1 - a**2*cy2*mu*xm**2*xx**2*ym*1.2D+1 + a**2*cx2*mu*xm&
                  &**2*xx**2*yy*1.2D+1 + a**2*cy2*mu*xm**2*xx**2*yy*1.2D+1 - a**2*cx2*mu*&
                  &xx**2*ym*yy**2*1.2D+1 + a**2*cx2*mu*xx**2*ym**2*yy*1.2D+1 - a**2*cy2*m&
                  &u*xx**2*ym*yy**2*1.2D+1 + a**2*cy2*mu*xx**2*ym**2*yy*1.2D+1 + a**2*mu*&
                  &t1*t2*xm**2*xx**3*3.6D+1 - a**2*mu*t1*t2*xm**3*xx**2*1.6D+1 + a**2*mu*&
                  &t1*t2*xx**3*ym**2*2.0D+1 + a**2*mu*t1*t2*xx**3*yy**2*2.0D+1 - a**2*mu*&
                  &t1*xm*xx**2*ym**2*8.0D0 + a**2*mu*t1*xm**2*xx*ym**2*4.0D0 + a**2*mu*t2&
                  &*xm**2*xx**2*ym*1.2D+1 - a**2*mu*t1*xm*xx**2*yy**2*8.0D0 + a**2*mu*t1*&
                  &xm**2*xx*yy**2*4.0D0 - a**2*mu*t2*xm**2*xx**2*yy*1.2D+1 + a**2*mu*t1*x&
                  &x*ym**2*yy**2*1.2D+1 + a**2*mu*t2*xx**2*ym*yy**2*1.2D+1 - a**2*mu*t2*x&
                  &x**2*ym**2*yy*1.2D+1 + a*cx2*cy*sy*xx**5*(1.0D0/xx**2*(xm*xx*(-2.0D0&
                  &) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*4.0D0&
                  &- a*cx2*cy*mu*sy*xx**4*1.2D+1 - cx2*cy2*t2*xx**3*ym*SQRT(1.0D0/xx**2*&
                  &(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)) + cx2*cy&
                  &2*t2*xx*yy**3*SQRT(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + x&
                  &x**2*2.0D0 + ym**2 + yy**2)) + a**2*mu*xm*xx*ym**3*COS(a*yy)**2*4.0D0 + a*&
                  &*2*mu*xm**3*xx*ym*COS(a*yy)**2*4.0D0 - a*cx2*t4*t6*xm*xx**4*2.0D0 - cx&
                  &2*cy2*t5*xm**2*xx*ym*2.0D0 - cx2*cy2*t5*xm*xx**2*yy*4.0D0 + cx2*cy2*t5&
                  &*xm**2*xx*yy*2.0D0 - cx2*cy2*t5*xx*ym*yy**2*6.0D0 + a*mu*t3*xm*xx**2*y&
                  &m*8.0D0 - a*mu*t3*xm**2*xx*ym*8.0D0 - a*mu*t4*xm*xx*ym**2*6.0D0 - a*mu*t&
                  &3*xm*xx**2*yy*8.0D0 + a*mu*t3*xm**2*xx*yy*8.0D0 - a*mu*t4*xm*xx*yy**2*&
                  &6.0D0 + a*mu*t3*xm*ym*yy**2*6.0D0 - a*mu*t3*xm*ym**2*yy*6.0D0 - a*mu*t4*&
                  &xm**2*ym*yy*4.0D0 - a*mu*t3*xx*ym*yy**2*1.2D+1 + a*mu*t3*xx*ym**2*yy*1&
                  &.2D+1 - a*mu*t4*xx**2*ym*yy*2.0D+1 + a*mu*xm*xx**2*ym*SIN(a*xx*2.0D0)*&
                  &2.0D0 - a*mu*xm**2*xx*ym*SIN(a*xx*2.0D0)*2.0D0 - a*mu*xm*xx**2*yy*SIN(&
                  &a*xx*2.0D0)*2.0D0 + a*mu*xm**2*xx*yy*SIN(a*xx*2.0D0)*2.0D0 - a*mu*xm*x&
                  &x*ym**2*SIN(a*yy*2.0D0)*(3.0D0/2.0D0) + a*mu*xm*ym*yy**2*SIN(a*xx*2.&
                  &0D0)*(3.0D0/2.0D0) - a*mu*xm*ym**2*yy*SIN(a*xx*2.0D0)*(3.0D0/2.0D0) -&
                  &a*mu*xx*ym*yy**2*SIN(a*xx*2.0D0)*3.0D0 + a*mu*xx*ym**2*yy*SIN(a*xx*2&
                  &.0D0)*3.0D0 - a*mu*xm*xx*yy**2*SIN(a*yy*2.0D0)*(3.0D0/2.0D0) - a**2*cx&
                  &2*cy2*mu*xx**4*yy*1.6D+1 - a*mu*xm**2*ym*yy*SIN(a*yy*2.0D0) - a*mu*xx*&
                  &*2*ym*yy*SIN(a*yy*2.0D0)*5.0D0 - a*cx2*cy**3*sx*t6*xx**5*3.0D0 + a**2*&
                  &cx2*mu*xm*xx*ym**3*4.0D0 + a**2*cx2*mu*xm*xx**3*ym*1.6D+1 + a**2*cx2*m&
                  &u*xm**3*xx*ym*4.0D0 + a**2*cy2*mu*xm*xx**3*ym*1.6D+1 - a**2*cx2*mu*xm*&
                  &xx*yy**3*4.0D0 - a**2*cx2*mu*xm*xx**3*yy*1.6D+1 - a**2*cx2*mu*xm**3*xx&
                  &*yy*4.0D0 - a**2*cy2*mu*xm*xx*yy**3*4.0D0 - a**2*cy2*mu*xm*xx**3*yy*1.&
                  &6D+1 - a**2*cy2*mu*xm**3*xx*yy*4.0D0 - a**2*mu*t1*t2*xm*xx**4*4.0D+1 + a&
                  &**2*mu*t1*t2*xm**4*xx*4.0D0 - a*mu*sx*t1**2*xm*ym**3*2.0D0 - a*mu*sx*t&
                  &1**2*xm**3*ym*2.0D0 + a*mu*sx*t1**2*xx**3*ym*4.0D0 + a*mu*sx*t1**2*xm*&
                  &yy**3*2.0D0 + a*mu*sx*t1**2*xm**3*yy*2.0D0 - a*mu*sx*t1**2*xx**3*yy*4.&
                  &0D0 + a**2*mu*t1*t2*xx*ym**4*4.0D0 + a**2*mu*t1*t2*xx*yy**4*4.0D0 - a**2&
                  &*mu*t2*xm*xx*ym**3*4.0D0 - a**2*mu*t2*xm*xx**3*ym*1.6D+1 - a**2*mu*t2*&
                  &xm**3*xx*ym*4.0D0 + a**2*mu*t2*xm*xx*yy**3*4.0D0 + a**2*mu*t2*xm*xx**3&
                  &*yy*1.6D+1 + a**2*mu*t2*xm**3*xx*yy*4.0D0 - a**2*mu*t1*xx*ym*yy**3*8.0&
                  &D0 - a**2*mu*t1*xx*ym**3*yy*8.0D0 - a**2*mu*t1*xx**3*ym*yy*2.0D+1 + a*sx&
                  &*t1**2*t6*xx**4*ym*4.0D0 + a*sy*t1**2*t6*xx**4*ym*2.0D0 - a*sy*t1**2*t&
                  &6*xx**4*yy*2.0D0 - a**2*mu*xm*xx*ym*yy**2*6.0D0 + a**2*mu*xm*xx*ym**2*&
                  &yy*6.0D0 - a*cy*mu*sy*xm**4*COS(a*xx)**2*2.0D0 - cy2*t2*xx*ym**3*COS(a&
                  &*xx)**2*SQRT(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2&
                  &.0D0 + ym**2 + yy**2)) + cy2*t2*xx**3*yy*COS(a*xx)**2*SQRT(1.0D0/xx**2*(&
                  &xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)) + cy2*t5*&
                  &xm*xx**2*ym*COS(a*xx)**2*4.0D0 + cx2*t5*xx*ym**2*yy*COS(a*yy)**2*6.0&
                  &D0 + a*mu*t1**2*xx*ym**3*SIN(a*xx)*4.0D0 - a*mu*t1**2*xx*yy**3*SIN(a*x&
                  &x)*4.0D0 + a**2*cy2*mu*xx**2*ym*yy**2*COS(a*xx)**2*2.4D+1 - a*cx2*cy*m&
                  &u*sy*xm**2*xx**2*2.2D+1 - a*cx2*cy*mu*sy*xm**2*ym**2*2.0D0 - a*cx2*cy*&
                  &mu*sy*xx**2*ym**2*1.0D+1 - a*cx2*cy*mu*sy*xm**2*yy**2*2.0D0 - a**2*cx2&
                  &*cy2*mu*xm*xx*ym**3*8.0D0 - a**2*cx2*cy2*mu*xm**3*xx*ym*8.0D0 + a**2*c&
                  &x2*cy2*mu*xm*xx*yy**3*8.0D0 + a**2*cx2*cy2*mu*xm*xx**3*yy*3.2D+1 + a**&
                  &2*cx2*cy2*mu*xm**3*xx*yy*8.0D0 + a*cx**3*cy2*sy*t6*xx**4*yy*3.0D0 - a*&
                  &*2*cx2*mu*xm*xx*ym**2*yy*1.2D+1 - a**2*cy2*mu*xm*xx*ym**2*yy*1.2D+1 -&
                  &a*cy2*sy*xx**4*ym*COS(a*xx)**3*(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*&
                  &2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*3.0D0 - a*mu*sx&
                  &*t1**2*xm*xx**2*ym*8.0D0 + a*mu*sx*t1**2*xm*xx**2*yy*8.0D0 - a*mu*sx*t&
                  &1**2*xm*ym*yy**2*6.0D0 + a*mu*sx*t1**2*xm*ym**2*yy*6.0D0 + a*mu*sx*t1*&
                  &*2*xx*ym*yy**2*1.2D+1 - a*mu*sx*t1**2*xx*ym**2*yy*1.2D+1 - a**2*mu*t1*&
                  &t2*xx*ym*yy**3*1.6D+1 - a**2*mu*t1*t2*xx*ym**3*yy*1.6D+1 - a**2*mu*t1*&
                  &t2*xx**3*ym*yy*4.0D+1 + a**2*mu*t1*xm*xx**2*ym*yy*1.6D+1 - a**2*mu*t1*&
                  &xm**2*xx*ym*yy*8.0D0 - a**2*mu*t2*xm*xx*ym*yy**2*1.2D+1 + a**2*mu*t2*x&
                  &m*xx*ym**2*yy*1.2D+1 - cx2*t2*xm**2*xx*ym*COS(a*yy)**2*SQRT(1.0D0/xx&
                  &**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)) - cy&
                  &2*t2*xx*ym*yy**2*COS(a*xx)**2*SQRT(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*&
                  &yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*3.0D0 + a*mu*t1**2*xm**2*xx&
                  &*ym*SIN(a*xx)*8.0D0 - a*mu*t1**2*xm**2*xx*yy*SIN(a*xx)*8.0D0 + a**2*cx&
                  &2*cy2*mu*xm**2*xx**2*ym*2.4D+1 - a**2*cx2*cy2*mu*xm**2*xx**2*yy*2.4D&
                  &+1 - a**2*cx2*cy2*mu*xx**2*ym**2*yy*2.4D+1 + a*mu*t4*xm*xx*ym*yy*1.2D+&
                  &1 + a**2*cx**2*mu*xm*xx*ym*yy**2*1.2D+1 + a**2*cy**2*mu*xm*xx*ym*yy**2&
                  &*1.2D+1 - a**2*mu*t1*t2*xm*xx**2*ym**2*1.6D+1 + a**2*mu*t1*t2*xm**2*xx&
                  &*ym**2*8.0D0 - a**2*mu*t1*t2*xm*xx**2*yy**2*1.6D+1 + a**2*mu*t1*t2*xm*&
                  &*2*xx*yy**2*8.0D0 + a**2*mu*t1*t2*xx*ym**2*yy**2*2.4D+1 - a*cy*mu*sy*x&
                  &x**2*yy**2*COS(a*xx)**2*1.0D+1 - a**2*cy2*mu*xm*xx**3*ym*COS(a*xx)**&
                  &2*3.2D+1 + a*cy**3*sx*t6*xm*xx**4*COS(a*xx)**2*3.0D0 + a*cx2*cy*mu*sy*&
                  &xm*xx**3*2.4D+1 + cx2*cy2*t2*xm*xx**2*ym*SQRT(1.0D0/xx**2*(xm*xx*(-2&
                  &.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*2.0D0 - cx2*cy2*t2&
                  &*xm*xx**2*yy*SQRT(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx&
                  &**2*2.0D0 + ym**2 + yy**2))*2.0D0 + cx2*cy2*t2*xm**2*xx*yy*SQRT(1.0D0/xx&
                  &**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)) + a*&
                  &mu*xm*xx*ym*yy*SIN(a*yy*2.0D0)*3.0D0 + cx2*cy2*t2*xx*ym**2*yy*SQRT(1&
                  &.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy*&
                  &*2))*3.0D0 - a*cx2*cy*sy*t6*xm*xx**4*4.0D0 + a*cx2*mu*sy*xm**3*xx*COS(&
                  &a*yy)*1.0D+1 - a**2*cx2*mu*xm*xx*ym*yy**2*COS(a*yy)**2*2.4D+1 + a*cx2*&
                  &cy*mu*sy*xm*xx*ym**2*6.0D0 + a*cx2*cy*mu*sy*xx**2*ym*yy*2.0D+1 + a*cx2&
                  &*cy*mu*xm**2*ym*yy*SIN(a*yy)*4.0D0 + a**2*cx2*cy2*mu*xm*xx*ym**2*yy*&
                  &2.4D+1 + a**2*mu*t1*t2*xm*xx**2*ym*yy*3.2D+1 - a**2*mu*t1*t2*xm**2*xx*&
                  &ym*yy*1.6D+1 + a*cy*mu*sy*xm*xx*yy**2*COS(a*xx)**2*6.0D0 - a*cy*mu*sy*&
                  &xm*xx*ym*yy*COS(a*xx)**2*1.2D+1))/xx

          CASE (50:)
             !Do nothing
          CASE DEFAULT
             WRITE (6, *) "Error! Test case not valid"
             STOP

          END SELECT
       END DO
    END DO
  END SUBROUTINE body_force
#else
  !***********************************************************************
  !
  !                            VERSION 2D
  !
  !***********************************************************************
  !*****************************************
  ! Analytical solution
  !****************************************
  SUBROUTINE analytical_solution(x, y, u)
    REAL*8, DIMENSION(:), INTENT(IN)        :: x, y
    REAL*8, DIMENSION(:, :), INTENT(OUT)     :: u
    REAL*8, DIMENSION(SIZE(u, 1), phys%npv)  :: up
    REAL*8 :: a,b,xx,yy,tt,xmax,xmin,ymax,ymin,xm,ym
    REAL*8 :: dsource(SIZE(x)),aux(SIZE(x)),xsource,ysource,smod
    INTEGER:: np,i
    REAL*8 :: r(SIZE(x)),rs

    np = SIZE(x)
    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    up = 0.
    u = 0.
    a = 2*pi


    SELECT CASE (switch%testcase)
    CASE (1)
       IF (switch%axisym) THEN
          WRITE (6, *) "This is NOT an axisymmetric test case!"
          STOP
       END IF
       ! Circular field centered in [xc, yc], n = 2+sin(a*x)*sin(a*y),  u = cos(a*x)*cos(a*y)
       ! Case 9 of the Matlab version: for convergence purpose
       up(:, 1) = 2 + SIN(a*x)*SIN(a*y)
       up(:, 2) = COS(a*x)*COS(a*y)
    CASE (2)
       ! Axisimmetric case with div(b)~=0
       IF (.NOT. switch%axisym) THEN
          WRITE (6, *) "This is an axisymmetric test case!"
          STOP
       END IF
       up(:, 1) = 2 + SIN(a*x)*SIN(a*y)
       up(:, 2) = COS(a*x)*COS(a*y)
    CASE (5)
       IF (switch%axisym) THEN
          WRITE (6,*) "This is NOT an axisymmetric test case!"
          STOP
       END IF
       smod = 1.
       rs = 0.04/simpar%refval_length
       xsource = xm-0.5*(xmax-xm)
       ysource = ym
       dsource   = SQRT((x-xsource)**2+(y-ysource)**2)
       aux = -dsource**2/rs**2
       up(:,1) = 1e-3
       DO i=1,np
          IF (aux(i).GT.-30) THEN
             up(i,1) =  up(i,1)+smod*EXP(aux(i))
          ENDIF
       END DO
    CASE (6)
       IF (.NOT.switch%axisym) THEN
          WRITE (6,*) "This is an axisymmetric test case!"
          STOP
       END IF
       !
       smod = 1.
       rs = 0.04/simpar%refval_length
       xsource = xm-0.5*(xmax-xm)
       ysource = ym
       dsource   = SQRT((x-xsource)**2+(y-ysource)**2)
       aux = -dsource**2/rs**2
       up(:,1) = 1e-6
       DO i=1,np
          IF (aux(i).GT.-30) THEN
             up(i,1) =  up(i,1)+smod*EXP(aux(i))
          ENDIF
       END DO

    CASE (50:64)


       up(:, 1) = 1.
       up(:, 2) = 0.
#ifdef NEUTRAL
       up(:,4)  = 0.
#endif
    CASE (65)
       up(:, 1) = 1.
       up(:, 2) = 0.
       r = SQRT((x*phys%lscale - geom%R0)**2 + (y*phys%lscale - 0.75)**2)
       DO i = 1, SIZE(x)
          IF (r(i) .LE. 0.05) THEN
             up(i, 2) = 1.
          END IF
       END DO
    CASE (80:89)
       up(:, 1) = 1.
       up(:, 2) = 0.
#ifdef NEUTRAL
       up(:,4)  = 0.
#endif

    CASE DEFAULT
       WRITE (6, *) "Error! Test case not valid"
       STOP
    END SELECT
    ! Convert physical variables to conservative variables
    CALL phys2cons(up, u)



  END SUBROUTINE analytical_solution

  !*****************************************
  ! Analytical gradient
  !****************************************
  SUBROUTINE analytical_gradient(x, y, u, ux, uy)
    REAL*8, DIMENSION(:), INTENT(IN)        :: x, y
    REAL*8, DIMENSION(:, :), INTENT(IN)      :: u
    REAL*8, DIMENSION(:, :), INTENT(OUT)     :: ux, uy
    REAL*8, DIMENSION(SIZE(u, 1), phys%npv)  :: upx, upy
    REAL*8, DIMENSION(SIZE(u, 1), phys%npv)  :: up
    REAL*8 :: a

    upx = 0.
    upy = 0.
    ux = 0.
    uy = 0.
    CALL cons2phys(u, up)
    a = 2*pi
    SELECT CASE (switch%testcase)
    CASE (1)
       IF (switch%axisym) THEN
          WRITE (6, *) "This is NOT an axisymmetric test case!"
          STOP
       END IF
       ! Circular field centered in [xc, yc], n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
       upx(:, 1) = a*COS(a*x)*SIN(a*y)
       upx(:, 2) = -a*SIN(a*x)*COS(a*y)

       upy(:, 1) = a*SIN(a*x)*COS(a*y)
       upy(:, 2) = -a*COS(a*x)*SIN(a*y)

    CASE (2)
       ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
       IF (.NOT. switch%axisym) THEN
          WRITE (6, *) "This is an axisymmetric test case!"
          STOP
       END IF
       upx(:, 1) = a*COS(a*x)*SIN(a*y)
       upx(:, 2) = -a*SIN(a*x)*COS(a*y)

       upy(:, 1) = a*SIN(a*x)*COS(a*y)
       upy(:, 2) = -a*COS(a*x)*SIN(a*y)
    CASE (5)
       ! Do nothing
    CASE (6)
       ! Do nothing
    CASE (50:)
       ! Do nothing
    CASE DEFAULT
       WRITE (6, *) "Error! Test case not valid"
       STOP
    END SELECT
    ! Convert physical variables to conservative variables
    ux(:, 2) = (upx(:, 1)*up(:, 2) + up(:, 1)*upx(:, 2))
    uy(:, 2) = (upy(:, 1)*up(:, 2) + up(:, 1)*upy(:, 2))

    IF (switch%logrho) THEN
       ux(:,1) = upx(:,1)/up(:,1)
       uy(:,1) = upy(:,1)/up(:,1)
    ELSE
       ux(:, 1) = upx(:, 1)
       uy(:, 1) = upy(:, 1)
    END IF

  END SUBROUTINE analytical_gradient

  !*****************************************
  ! Body forces
  !****************************************
  SUBROUTINE body_force(x, y, f)
    REAL*8, DIMENSION(:), INTENT(IN) :: x, y
    REAL*8, DIMENSION(:, :), INTENT(OUT) :: f
    INTEGER                            :: i, n
    REAL*8  :: a, b, xc, yc, D, mu, k
    REAL*8  :: csii, csie, kpari, kpare, Mref, epn, tie, pcourr
    INTEGER:: j, ind, N1D, N2D
    REAL*8 :: r, xx, yy, tt, xmax, xmin, ymax, ymin, xm, ym
    REAL*8 :: aux
    REAL*8 :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
    REAL*8 :: t11, t12, t13, t14, t15, t16, t17, t18, t19, t20
    REAL*8 :: t21, t22, t23, t24, t25, t26, t27, t28, t29, t30
    REAL*8 :: t31, t32, t33, t34, t35, t36, t37, t38, t39, t40
    REAL*8 :: t41, t42, t43, t44, t45, t46, t47, t48, t49, t50
    REAL*8 :: t51, t52, t53, t54, t55, t56, t57, t58, t59, t60
    REAL*8 :: t61, t62, t63, t64, t65, t66, t67, t68, t69, t70
    REAL*8 :: cx, cy, sx, sy, cx2, sx2, cy2, sy2, r2, cc, ss, ct, st

    n = SIZE(x)

    f = 0.
    a = 2*pi

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)


    SELECT CASE (switch%testcase)
    CASE (1)
       ! Circular field centered in [xc, yc], n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y)
       ! Case 9 of the Matlab version: for convergence purpose
       a = 2*pi
       b = 2*pi
       xc = 0.
       yc = 0.
       D = phys%diff_n
       mu = phys%diff_u
       k = phys%a
       IF (switch%logrho) THEN
          f(:,1) = (a*(2*D*a*SIN(a*x)*SIN(a*y)+D*xc*COS(a*x)*SIN(a*y)-D*x*COS(a*x)*SIN(a*y)+D*yc*COS(a*y)*SIN(a*x)-D*y*COS(a*y)*SIN(a*x)+D*xc**3*COS(a*x)*SIN(a*y)-D*x**3*COS(a*x)*SIN(a*y)+D*yc**3*COS(a*y)*SIN(a*x)-D*y**3*COS(a*y)*SIN(a*x)+3*D*a*xc**2*SIN(a*x)*SIN(a*y)+D*a*xc**4*SIN(a*x)*SIN(a*y)+3*D*a*x**2*SIN(a*x)*SIN(a*y)+D*a*x**4*SIN(a*x)*SIN(a*y)+3*D*a*yc**2*SIN(a*x)*SIN(a*y)+D*a*yc**4*SIN(a*x)*SIN(a*y)+3*D*a*y**2*SIN(a*x)*SIN(a*y)+D*a*y**4*SIN(a*x)*SIN(a*y)+3*D*xc*x**2*COS(a*x)*SIN(a*y)-3*D*xc**2*x*COS(a*x)*SIN(a*y)+D*xc*yc**2*COS(a*x)*SIN(a*y)+D*xc**2*yc*COS(a*y)*SIN(a*x)-D*x*yc**2*COS(a*x)*SIN(a*y)+D*x**2*yc*COS(a*y)*SIN(a*x)+D*xc*y**2*COS(a*x)*SIN(a*y)-D*xc**2*y*COS(a*y)*SIN(a*x)-D*x*y**2*COS(a*x)*SIN(a*y)-D*x**2*y*COS(a*y)*SIN(a*x)+3*D*yc*y**2*COS(a*y)*SIN(a*x)-3*D*yc**2*y*COS(a*y)*SIN(a*x)-2*xc*COS(a*x)*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**1.5+2*x*COS(a*x)*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**1.5+2*yc*COS(a*y)*SIN(a*x)*((xc-x)**2+(yc-y)**2+1)**1.5-2*y*COS(a*y)*SIN(a*x)*((xc-x)**2+(yc-y)**2+1)**1.5+6*D*a*xc**2*x**2*SIN(a*x)*SIN(a*y)+2*D*a*xc**2*yc**2*SIN(a*x)*SIN(a*y)+2*D*a*x**2*yc**2*SIN(a*x)*SIN(a*y)+2*D*a*xc**2*y**2*SIN(a*x)*SIN(a*y)+2*D*a*x**2*y**2*SIN(a*x)*SIN(a*y)+6*D*a*yc**2*y**2*SIN(a*x)*SIN(a*y)-2*D*a*xc*yc*COS(a*x)*COS(a*y)+2*D*a*x*yc*COS(a*x)*COS(a*y)+2*D*a*xc*y*COS(a*x)*COS(a*y)-2*D*a*x*y*COS(a*x)*COS(a*y)-6*D*a*xc*x*SIN(a*x)*SIN(a*y)-6*D*a*yc*y*SIN(a*x)*SIN(a*y)-2*D*xc*x*yc*COS(a*y)*SIN(a*x)+2*D*xc*x*y*COS(a*y)*SIN(a*x)-2*D*xc*yc*y*COS(a*x)*SIN(a*y)+2*D*x*yc*y*COS(a*x)*SIN(a*y)+xc*COS(a*x)*COS(a*y)**2*SIN(a*x)*((xc-x)**2+(yc-y)**2+1)**1.5-x*COS(a*x)*COS(a*y)**2*SIN(a*x)*((xc-x)**2+(yc-y)**2+1)**1.5-yc*COS(a*x)**2*COS(a*y)*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**1.5+y*COS(a*x)**2*COS(a*y)*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**1.5-xc*COS(a*x)*SIN(a*x)*SIN(a*y)**2*((xc-x)**2+(yc-y)**2+1)**1.5+x*COS(a*x)*SIN(a*x)*SIN(a*y)**2*((xc-x)**2+(yc-y)**2+1)**1.5+yc*COS(a*y)*SIN(a*x)**2*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**1.5-y*COS(a*y)*SIN(a*x)**2*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**1.5-2*D*a*xc*yc**3*COS(a*x)*COS(a*y)-2*D*a*xc**3*yc*COS(a*x)*COS(a*y)+2*D*a*x*yc**3*COS(a*x)*COS(a*y)+2*D*a*x**3*yc*COS(a*x)*COS(a*y)+2*D*a*xc*y**3*COS(a*x)*COS(a*y)+2*D*a*xc**3*y*COS(a*x)*COS(a*y)-2*D*a*x*y**3*COS(a*x)*COS(a*y)-2*D*a*x**3*y*COS(a*x)*COS(a*y)-4*D*a*xc*x**3*SIN(a*x)*SIN(a*y)-4*D*a*xc**3*x*SIN(a*x)*SIN(a*y)-4*D*a*yc*y**3*SIN(a*x)*SIN(a*y)-4*D*a*yc**3*y*SIN(a*x)*SIN(a*y)-6*D*a*xc*x**2*yc*COS(a*x)*COS(a*y)+6*D*a*xc**2*x*yc*COS(a*x)*COS(a*y)+6*D*a*xc*x**2*y*COS(a*x)*COS(a*y)-6*D*a*xc**2*x*y*COS(a*x)*COS(a*y)-6*D*a*xc*yc*y**2*COS(a*x)*COS(a*y)+6*D*a*xc*yc**2*y*COS(a*x)*COS(a*y)+6*D*a*x*yc*y**2*COS(a*x)*COS(a*y)-6*D*a*x*yc**2*y*COS(a*x)*COS(a*y)-4*D*a*xc*x*yc**2*SIN(a*x)*SIN(a*y)-4*D*a*xc*x*y**2*SIN(a*x)*SIN(a*y)-4*D*a*xc**2*yc*y*SIN(a*x)*SIN(a*y)-4*D*a*x**2*yc*y*SIN(a*x)*SIN(a*y)+8*D*a*xc*x*yc*y*SIN(a*x)*SIN(a*y)))/((SIN(a*x)*SIN(a*y)+2)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**2)
       ELSE
          f(:,1) = (2*D*a**2*SIN(a*x)*SIN(a*y)-2*a*xc*COS(a*x)*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**0.5+2*a*x*COS(a*x)*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**0.5+2*a*yc*COS(a*y)*SIN(a*x)*((xc-x)**2+(yc-y)**2+1)**0.5-2*a*y*COS(a*y)*SIN(a*x)*((xc-x)**2+(yc-y)**2+1)**0.5+D*a**2*xc**2*SIN(a*x)*SIN(a*y)+D*a**2*x**2*SIN(a*x)*SIN(a*y)+D*a**2*yc**2*SIN(a*x)*SIN(a*y)+D*a**2*y**2*SIN(a*x)*SIN(a*y)+D*a*xc*COS(a*x)*SIN(a*y)-D*a*x*COS(a*x)*SIN(a*y)+D*a*yc*COS(a*y)*SIN(a*x)-D*a*y*COS(a*y)*SIN(a*x)+a*xc*COS(a*x)*COS(a*y)**2*SIN(a*x)*((xc-x)**2+(yc-y)**2+1)**0.5-a*x*COS(a*x)*COS(a*y)**2*SIN(a*x)*((xc-x)**2+(yc-y)**2+1)**0.5-a*yc*COS(a*x)**2*COS(a*y)*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**0.5+a*y*COS(a*x)**2*COS(a*y)*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**0.5-a*xc*COS(a*x)*SIN(a*x)*SIN(a*y)**2*((xc-x)**2+(yc-y)**2+1)**0.5+a*x*COS(a*x)*SIN(a*x)*SIN(a*y)**2*((xc-x)**2+(yc-y)**2+1)**0.5+a*yc*COS(a*y)*SIN(a*x)**2*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**0.5-a*y*COS(a*y)*SIN(a*x)**2*SIN(a*y)*((xc-x)**2+(yc-y)**2+1)**0.5-2*D*a**2*xc*yc*COS(a*x)*COS(a*y)+2*D*a**2*x*yc*COS(a*x)*COS(a*y)+2*D*a**2*xc*y*COS(a*x)*COS(a*y)-2*D*a**2*x*y*COS(a*x)*COS(a*y)-2*D*a**2*xc*x*SIN(a*x)*SIN(a*y)-2*D*a**2*yc*y*SIN(a*x)*SIN(a*y))/((xc-x)**2+(yc-y)**2+1)
       ENDIF


       f(:,2) = ((a*mu*x*SIN(2*a*y))/2-(a*mu*xc*SIN(2*a*y))/2-(a*mu*yc*SIN(2*a*x))/2+(a*mu*y*SIN(2*a*x))/2-2*a**2*mu*xc*yc+2*a**2*mu*x*yc+2*a**2*mu*xc*y-2*a**2*mu*x*y+4*a**2*mu*COS(a*x)*COS(a*y)+4*a**2*mu*xc*yc*COS(a*x)**2-4*a**2*mu*x*yc*COS(a*x)**2-4*a**2*mu*xc*y*COS(a*x)**2+4*a**2*mu*x*y*COS(a*x)**2+4*a**2*mu*xc*yc*COS(a*y)**2-4*a**2*mu*x*yc*COS(a*y)**2-4*a**2*mu*xc*y*COS(a*y)**2+4*a**2*mu*x*y*COS(a*y)**2+2*a**2*mu*xc**2*COS(a*x)*COS(a*y)+2*a**2*mu*x**2*COS(a*x)*COS(a*y)+2*a**2*mu*yc**2*COS(a*x)*COS(a*y)+2*a**2*mu*y**2*COS(a*x)*COS(a*y)-2*a*mu*xc*COS(a*y)*SIN(a*x)+2*a*mu*x*COS(a*y)*SIN(a*x)-2*a*mu*yc*COS(a*x)*SIN(a*y)+2*a*mu*y*COS(a*x)*SIN(a*y)+8*a**2*mu*COS(a*x)*COS(a*y)*SIN(a*x)*SIN(a*y)-2*a*xc*COS(a*x)**2*COS(a*y)*SIN(a*x)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5+2*a*x*COS(a*x)**2*COS(a*y)*SIN(a*x)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5+4*a*yc*COS(a*x)*COS(a*y)**2*SIN(a*x)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5-4*a*xc*COS(a*x)**2*COS(a*y)*SIN(a*y)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5+4*a*x*COS(a*x)**2*COS(a*y)*SIN(a*y)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5-4*a*y*COS(a*x)*COS(a*y)**2*SIN(a*x)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5+2*a*yc*COS(a*x)*COS(a*y)**2*SIN(a*y)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5-2*a*y*COS(a*x)*COS(a*y)**2*SIN(a*y)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5-8*a**2*mu*xc*yc*COS(a*x)**2*COS(a*y)**2+8*a**2*mu*x*yc*COS(a*x)**2*COS(a*y)**2+8*a**2*mu*xc*y*COS(a*x)**2*COS(a*y)**2-8*a**2*mu*x*y*COS(a*x)**2*COS(a*y)**2+3*a*xc*COS(a*x)**2*COS(a*y)**3*SIN(a*x)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5-3*a*x*COS(a*x)**2*COS(a*y)**3*SIN(a*x)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5-3*a*yc*COS(a*x)**3*COS(a*y)**2*SIN(a*y)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5+3*a*y*COS(a*x)**3*COS(a*y)**2*SIN(a*y)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5+a*k*xc*COS(a*y)*SIN(a*x)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5-a*k*x*COS(a*y)*SIN(a*x)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5-a*k*yc*COS(a*x)*SIN(a*y)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5+a*k*y*COS(a*x)*SIN(a*y)*(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)**0.5-4*a**2*mu*xc*x*COS(a*x)*COS(a*y)-4*a**2*mu*yc*y*COS(a*x)*COS(a*y)-4*a**2*mu*xc*yc*SIN(a*x)*SIN(a*y)+4*a**2*mu*x*yc*SIN(a*x)*SIN(a*y)+4*a**2*mu*xc*y*SIN(a*x)*SIN(a*y)-4*a**2*mu*x*y*SIN(a*x)*SIN(a*y)+2*a*mu*yc*COS(a*x)*COS(a*y)**2*SIN(a*x)+2*a*mu*xc*COS(a*x)**2*COS(a*y)*SIN(a*y)-2*a*mu*x*COS(a*x)**2*COS(a*y)*SIN(a*y)-2*a*mu*y*COS(a*x)*COS(a*y)**2*SIN(a*x)+4*a**2*mu*xc**2*COS(a*x)*COS(a*y)*SIN(a*x)*SIN(a*y)+4*a**2*mu*x**2*COS(a*x)*COS(a*y)*SIN(a*x)*SIN(a*y)+4*a**2*mu*yc**2*COS(a*x)*COS(a*y)*SIN(a*x)*SIN(a*y)+4*a**2*mu*y**2*COS(a*x)*COS(a*y)*SIN(a*x)*SIN(a*y)-8*a**2*mu*xc*x*COS(a*x)*COS(a*y)*SIN(a*x)*SIN(a*y)-8*a**2*mu*yc*y*COS(a*x)*COS(a*y)*SIN(a*x)*SIN(a*y))/(xc**2-2*yc*y-2*xc*x+x**2+yc**2+y**2+1)
       !               f(:,1) = cos(a * x)**2 * a  * sin(b * y) * cos(b * y) * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y) * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * x-0.2e1 * xc) / 0.20e2+sin(a * x) * cos(b * y)**2 * b * cos(a * x) * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * y-0.2e1 * yc) / 0.20e2-D * (-sin(a * x) * a**2 * sin(b * y)+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * cos(b * y) * b / 0.10e2) * (0.2e1 * x-0.2e1 * xc) / 0.20e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) * (0.2e1 * x-0.2e1 * xc) / 0.20e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * a**2 * sin(b * y) / 0.10e2-(y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * cos(b * y) * b / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * sin(a * x) * cos(b * y) * b * (0.2e1 * x-0.2e1 * xc) / 0.20e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * cos(b * y) * b / 0.10e2) / 0.10e2-sin(a * x) * sin(b * y) * b**2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * cos(b * y) * b / 0.10e2) * (0.2e1 * y-0.2e1 * yc) / 0.20e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * ((y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) / 0.10e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) * (0.2e1 * y-0.2e1 * yc) / 0.20e2+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * cos(b * y) * b / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * sin(a * x) * cos(b * y) * b * (0.2e1 * y-0.2e1 * yc) / 0.20e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * sin(b * y) * b**2 / 0.10e2) / 0.10e2)
       !            f(:,2) = cos(a * x)**3 * a  * sin(b * y) * cos(b * y)**2 * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y)**2 * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * a  / 0.5e1-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x)**2 * cos(b * y)**2 * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * x-0.2e1 * xc) / 0.20e2+sin(a * x) * cos(b * y)**3 * b * cos(a * x)**2 * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x)**2 * cos(b * y) * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(b * y) * b / 0.5e1-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x)**2 * cos(b * y)**2 * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * y-0.2e1 * yc) / 0.20e2+k * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) / 0.10e2+k * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * cos(b * y) * b / 0.10e2-mu * (-0.3e1 * cos(a * x) * a**2 * sin(b * y) * cos(b * y) * sin(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * a**2 * cos(b * y)+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) / 0.10e2) * (0.2e1 * x-0.2e1 * xc) / 0.20e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) * (0.2e1 * x-0.2e1 * xc) / 0.20e2+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (-0.3e1 * cos(a * x) * a**2 * sin(b * y) * cos(b * y) * sin(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * a**2 * cos(b * y)) / 0.10e2-(y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) * (0.2e1 * x-0.2e1 * xc) / 0.20e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * cos(b * y)**2 * b-sin(a * x)**2 * cos(b * y)**2 * b * a -cos(a * x)**2 * a  * sin(b * y)**2 * b+(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * sin(b * y) * b) / 0.10e2) / 0.10e2-0.3e1 * sin(a * x) * cos(b * y) * b**2 * cos(a * x) * sin(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * b**2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) / 0.10e2) * (0.2e1 * y-0.2e1 * yc) / 0.20e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * ((y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) / 0.10e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) * (0.2e1 * y-0.2e1 * yc) / 0.20e2+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * cos(b * y)**2 * b-sin(a * x)**2 * cos(b * y)**2 * b * a -cos(a * x)**2 * a  * sin(b * y)**2 * b+(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * sin(b * y) * b) / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) * (0.2e1 * y-0.2e1 * yc) / 0.20e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (-0.3e1 * sin(a * x) * cos(b * y) * b**2 * cos(a * x) * sin(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * b**2) / 0.10e2) / 0.10e2)
    CASE (2)
       ! Axisimmetric case with div(b)~=0
       IF (.NOT. switch%axisym) THEN
          WRITE (6, *) "This is an axisymmetric test case!"
          STOP
       END IF
       a = 2*pi
       b = 2*pi
       xc = 0.
       yc = 0.
       D = phys%diff_n
       mu = phys%diff_u
       k = phys%a
       Mref = phys%Mref
       DO i=1,SIZE(x)
          xx = x(i)
          yy = y(i)
          f(i, 1) = (1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy&
               &**2)**2*(-D*a*xm**4*COS(a*xx)*SIN(a*yy)-D*a*xx**4*COS(a*xx)*SIN(a*&
               &yy)*6.0D0+a*xx**5*COS(a*xx)*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D&
               &0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0-xx*ym**3*COS(a*xx)*&
               &cos(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0&
               &)*2.0D0-xx**3*ym*COS(a*xx)*COS(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1&
               &.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0+xx*yy**3*COS(a*xx)*COS(a*yy)*sq&
               &rt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0+xx**&
               &3*yy*COS(a*xx)*COS(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(&
               &ym-yy)**2+1.0D0)*2.0D0+D*a**2*xx**5*SIN(a*xx)*SIN(a*yy)*6.0D0-D*a*&
               &*2*xm*xx**4*SIN(a*xx)*SIN(a*yy)*1.0D+1+D*a**2*xm**4*xx*SIN(a*xx)*s&
               &in(a*yy)+D*a**2*xx*ym**4*SIN(a*xx)*SIN(a*yy)+D*a**2*xx*yy**4*SIN(a&
               &*xx)*SIN(a*yy)-a*xx**5*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*(1.0D0/xx*&
               &*2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+a*xx**5&
               &*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx&
               &**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+D*a**2*xx**2*ym**3*COS(a*xx)*&
               &cos(a*yy)*2.0D0-D*a**2*xx**2*yy**3*COS(a*xx)*COS(a*yy)*2.0D0+D*a**&
               &2*xm**2*xx**3*SIN(a*xx)*SIN(a*yy)*9.0D0-D*a**2*xm**3*xx**2*SIN(a*x&
               &x)*SIN(a*yy)*4.0D0+D*a**2*xx**3*ym**2*SIN(a*xx)*SIN(a*yy)*5.0D0+D*&
               &a**2*xx**3*yy**2*SIN(a*xx)*SIN(a*yy)*5.0D0+D*a*xm*xx**3*COS(a*xx)*&
               &sin(a*yy)*1.2D+1+D*a*xm**3*xx*COS(a*xx)*SIN(a*yy)*5.0D0-D*a*xm*ym*&
               &*3*COS(a*yy)*SIN(a*xx)-D*a*xm**3*ym*COS(a*yy)*SIN(a*xx)+D*a*xx*ym*&
               &*3*COS(a*yy)*SIN(a*xx)*2.0D0+D*a*xx**3*ym*COS(a*yy)*SIN(a*xx)*2.0D&
               &0+D*a*xm*yy**3*COS(a*yy)*SIN(a*xx)+D*a*xm**3*yy*COS(a*yy)*SIN(a*xx&
               &)-D*a*xx*yy**3*COS(a*yy)*SIN(a*xx)*2.0D0-D*a*xx**3*yy*COS(a*yy)*si&
               &n(a*xx)*2.0D0-a*xm*xx**4*COS(a*xx)*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)*&
               &*2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0+a*xx**4*ym*c&
               &os(a*yy)*SIN(a*xx)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+&
               &1.0D0)**(3.0D0/2.0D0)*2.0D0-a*xx**4*yy*COS(a*yy)*SIN(a*xx)*(1.0D0/&
               &xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D&
               &0+xm*xx**2*ym*COS(a*xx)*COS(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D&
               &0/xx**2*(ym-yy)**2+1.0D0)*4.0D0-xm**2*xx*ym*COS(a*xx)*COS(a*yy)*sq&
               &rt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0-xm*x&
               &x**2*yy*COS(a*xx)*COS(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**&
               &2*(ym-yy)**2+1.0D0)*4.0D0+xm**2*xx*yy*COS(a*xx)*COS(a*yy)*SQRT(1.0&
               &D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0-xx*ym*yy**&
               &2*COS(a*xx)*COS(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-&
               &yy)**2+1.0D0)*6.0D0+xx*ym**2*yy*COS(a*xx)*COS(a*yy)*SQRT(1.0D0/xx*&
               &*2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*6.0D0+D*a**2*xx**4*ym*&
               &cos(a*xx)*COS(a*yy)*4.0D0-D*a**2*xx**4*yy*COS(a*xx)*COS(a*yy)*4.0D&
               &0-D*a*xm**2*xx**2*COS(a*xx)*SIN(a*yy)*1.1D+1-D*a*xm**2*ym**2*COS(a&
               &*xx)*SIN(a*yy)-D*a*xx**2*ym**2*COS(a*xx)*SIN(a*yy)*5.0D0-D*a*xm**2&
               &*yy**2*COS(a*xx)*SIN(a*yy)-D*a*xx**2*yy**2*COS(a*xx)*SIN(a*yy)*5.0&
               &D0+D*a*xm*xx*ym**2*COS(a*xx)*SIN(a*yy)*3.0D0-D*a*xm*xx**2*ym*COS(a&
               &*yy)*SIN(a*xx)*4.0D0+D*a*xm**2*xx*ym*COS(a*yy)*SIN(a*xx)*4.0D0+D*a&
               &*xm*xx*yy**2*COS(a*xx)*SIN(a*yy)*3.0D0+D*a*xm*xx**2*yy*COS(a*yy)*s&
               &in(a*xx)*4.0D0-D*a*xm**2*xx*yy*COS(a*yy)*SIN(a*xx)*4.0D0-D*a*xm*ym&
               &*yy**2*COS(a*yy)*SIN(a*xx)*3.0D0+D*a*xm*ym**2*yy*COS(a*yy)*SIN(a*x&
               &x)*3.0D0+D*a*xm**2*ym*yy*COS(a*xx)*SIN(a*yy)*2.0D0+D*a*xx*ym*yy**2&
               &*COS(a*yy)*SIN(a*xx)*6.0D0-D*a*xx*ym**2*yy*COS(a*yy)*SIN(a*xx)*6.0&
               &D0+D*a*xx**2*ym*yy*COS(a*xx)*SIN(a*yy)*1.0D+1-D*a**2*xm*xx*ym**3*c&
               &os(a*xx)*COS(a*yy)*2.0D0-D*a**2*xm*xx**3*ym*COS(a*xx)*COS(a*yy)*8.&
               &0D0-D*a**2*xm**3*xx*ym*COS(a*xx)*COS(a*yy)*2.0D0+D*a**2*xm*xx*yy**&
               &3*COS(a*xx)*COS(a*yy)*2.0D0+D*a**2*xm*xx**3*yy*COS(a*xx)*COS(a*yy)&
               &*8.0D0+D*a**2*xm**3*xx*yy*COS(a*xx)*COS(a*yy)*2.0D0-D*a**2*xx*ym*y&
               &y**3*SIN(a*xx)*SIN(a*yy)*4.0D0-D*a**2*xx*ym**3*yy*SIN(a*xx)*SIN(a*&
               &yy)*4.0D0-D*a**2*xx**3*ym*yy*SIN(a*xx)*SIN(a*yy)*1.0D+1-xx*ym**3*c&
               &os(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2&
               &+1.0D0/xx**2*(ym-yy)**2+1.0D0)-xx**3*ym*COS(a*xx)*COS(a*yy)*SIN(a*&
               &xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1&
               &.0D0)+xx*yy**3*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/&
               &xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)+xx**3*yy*COS(a*xx)*&
               &cos(a*yy)*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx&
               &**2*(ym-yy)**2+1.0D0)+a*xm*xx**4*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*&
               &(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D&
               &0)-a*xx**4*ym*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*(1.0D0/xx**2*(xm-xx&
               &)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+a*xx**4*yy*COS(a&
               &*xx)**2*COS(a*yy)*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(y&
               &m-yy)**2+1.0D0)**(3.0D0/2.0D0)-a*xm*xx**4*COS(a*xx)*SIN(a*xx)*SIN(&
               &a*yy)**2*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3&
               &.0D0/2.0D0)+a*xx**4*ym*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy)*(1.0D0/xx*&
               &*2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)-a*xx**4&
               &*yy*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0&
               &/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+D*a**2*xm**2*xx**2*ym*COS(&
               &a*xx)*COS(a*yy)*6.0D0-D*a**2*xm**2*xx**2*yy*COS(a*xx)*COS(a*yy)*6.&
               &0D0+D*a**2*xx**2*ym*yy**2*COS(a*xx)*COS(a*yy)*6.0D0-D*a**2*xx**2*y&
               &m**2*yy*COS(a*xx)*COS(a*yy)*6.0D0-D*a**2*xm*xx**2*ym**2*SIN(a*xx)*&
               &sin(a*yy)*4.0D0+D*a**2*xm**2*xx*ym**2*SIN(a*xx)*SIN(a*yy)*2.0D0-D*&
               &a**2*xm*xx**2*yy**2*SIN(a*xx)*SIN(a*yy)*4.0D0+D*a**2*xm**2*xx*yy**&
               &2*SIN(a*xx)*SIN(a*yy)*2.0D0+D*a**2*xx*ym**2*yy**2*SIN(a*xx)*SIN(a*&
               &yy)*6.0D0+D*a**2*xm*xx**2*ym*yy*SIN(a*xx)*SIN(a*yy)*8.0D0-D*a**2*x&
               &m**2*xx*ym*yy*SIN(a*xx)*SIN(a*yy)*4.0D0+xm*xx**2*ym*COS(a*xx)*COS(&
               &a*yy)*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*&
               &(ym-yy)**2+1.0D0)*2.0D0-xm**2*xx*ym*COS(a*xx)*COS(a*yy)*SIN(a*xx)*&
               &sin(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0&
               &)-xm*xx**2*yy*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/x&
               &x**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0+xm**2*xx*yy*co&
               &s(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+&
               &1.0D0/xx**2*(ym-yy)**2+1.0D0)-xx*ym*yy**2*COS(a*xx)*COS(a*yy)*SIN(&
               &a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2&
               &+1.0D0)*3.0D0+xx*ym**2*yy*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*&
               &sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*3.0D0-D*&
               &a*xm*xx*ym*yy*COS(a*xx)*SIN(a*yy)*6.0D0-D*a**2*xm*xx*ym*yy**2*COS(&
               &a*xx)*COS(a*yy)*6.0D0+D*a**2*xm*xx*ym**2*yy*COS(a*xx)*COS(a*yy)*6.&
               &0D0))/xx

          f(i,2) = (1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy&
               &**2)**2*(a*mu*xm**4*COS(a*yy)*SIN(a*xx)*2.0D0+a*mu*xx**4*COS(a*yy)&
               &*SIN(a*xx)*1.2D+1+a**2*mu*xx**5*COS(a*xx)*COS(a*yy)*1.2D+1-xx*ym**&
               &3*COS(a*xx)**2*COS(a*yy)**2*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**&
               &2*(ym-yy)**2+1.0D0)*2.0D0-xx**3*ym*COS(a*xx)**2*COS(a*yy)**2*SQRT(&
               &1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0+xx*yy**&
               &3*COS(a*xx)**2*COS(a*yy)**2*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**&
               &2*(ym-yy)**2+1.0D0)*2.0D0+xx**3*yy*COS(a*xx)**2*COS(a*yy)**2*SQRT(&
               &1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0+D*a**2*&
               &xx**2*ym**3*COS(a*xx)**2*COS(a*yy)**2*2.0D0-D*a**2*xx**2*yy**3*cos&
               &(a*xx)**2*COS(a*yy)**2*2.0D0-a**2*mu*xm*xx**4*COS(a*xx)*COS(a*yy)*&
               &2.0D+1+a**2*mu*xm**4*xx*COS(a*xx)*COS(a*yy)*2.0D0+a**2*mu*xx*ym**4&
               &*COS(a*xx)*COS(a*yy)*2.0D0+a**2*mu*xx*yy**4*COS(a*xx)*COS(a*yy)*2.&
               &0D0-D*a**2*xx**2*ym**3*COS(a*xx)**2*SIN(a*yy)**2-D*a**2*xx**2*ym**&
               &3*COS(a*yy)**2*SIN(a*xx)**2+D*a**2*xx**2*yy**3*COS(a*xx)**2*SIN(a*&
               &yy)**2+D*a**2*xx**2*yy**3*COS(a*yy)**2*SIN(a*xx)**2-D*a*xm**4*COS(&
               &a*xx)**2*COS(a*yy)*SIN(a*yy)-D*a*xx**4*COS(a*xx)**2*COS(a*yy)*SIN(&
               &a*yy)*6.0D0+a*mu*xm**2*xx**2*COS(a*yy)*SIN(a*xx)*2.2D+1+a*mu*xm**2&
               &*ym**2*COS(a*yy)*SIN(a*xx)*2.0D0+a*mu*xx**2*ym**2*COS(a*yy)*SIN(a*&
               &xx)*1.0D+1+a*mu*xm**2*yy**2*COS(a*yy)*SIN(a*xx)*2.0D0+a*mu*xx**2*y&
               &y**2*COS(a*yy)*SIN(a*xx)*1.0D+1+a*xx**5*COS(a*xx)**2*COS(a*yy)*sin&
               &(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0&
               &D0/2.0D0)*4.0D0+a**2*mu*xx**4*ym*SIN(a*xx)*SIN(a*yy)*8.0D0-a**2*mu&
               &*xx**4*yy*SIN(a*xx)*SIN(a*yy)*8.0D0-a**2*mu*xx**2*ym**3*COS(a*xx)*&
               &*2*SIN(a*yy)**2-a**2*mu*xx**2*ym**3*COS(a*yy)**2*SIN(a*xx)**2+a**2&
               &*mu*xx**2*yy**3*COS(a*xx)**2*SIN(a*yy)**2+a**2*mu*xx**2*yy**3*COS(&
               &a*yy)**2*SIN(a*xx)**2+a**2*mu*xx**2*ym**3*SIN(a*xx)**2*SIN(a*yy)**&
               &2*2.0D0-a**2*mu*xx**2*yy**3*SIN(a*xx)**2*SIN(a*yy)**2*2.0D0+a*mu*x&
               &m**4*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy)+a*mu*xx**4*COS(a*yy)*SIN(a*x&
               &x)**2*SIN(a*yy)*6.0D0+xm*xx**2*ym*COS(a*xx)**2*COS(a*yy)**2*SQRT(1&
               &.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*4.0D0-xm**2*xx&
               &*ym*COS(a*xx)**2*COS(a*yy)**2*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx&
               &**2*(ym-yy)**2+1.0D0)*2.0D0-xm*xx**2*yy*COS(a*xx)**2*COS(a*yy)**2*&
               &sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*4.0D0+xm&
               &**2*xx*yy*COS(a*xx)**2*COS(a*yy)**2*SQRT(1.0D0/xx**2*(xm-xx)**2+1.&
               &0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0-xx*ym*yy**2*COS(a*xx)**2*COS(a*y&
               &y)**2*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*6.&
               &0D0+xx*ym**2*yy*COS(a*xx)**2*COS(a*yy)**2*SQRT(1.0D0/xx**2*(xm-xx)&
               &**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*6.0D0+a**2*mu*xm**2*xx**3*COS(a*&
               &xx)*COS(a*yy)*1.8D+1-a**2*mu*xm**3*xx**2*COS(a*xx)*COS(a*yy)*8.0D0&
               &+a**2*mu*xx**3*ym**2*COS(a*xx)*COS(a*yy)*1.0D+1+a**2*mu*xx**3*yy**&
               &2*COS(a*xx)*COS(a*yy)*1.0D+1-Mref*a*xx**5*COS(a*yy)*SIN(a*xx)*(1.0&
               &D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*2&
               &.0D0-a*xx**5*COS(a*xx)**2*COS(a*yy)**3*SIN(a*xx)*(1.0D0/xx**2*(xm-&
               &xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+a**2*mu*xx**2*&
               &ym**3*SIN(a*xx)*SIN(a*yy)*4.0D0-a**2*mu*xx**2*yy**3*SIN(a*xx)*SIN(&
               &a*yy)*4.0D0+D*a**2*xx**4*ym*COS(a*xx)**2*COS(a*yy)**2*4.0D0-D*a**2&
               &*xx**4*yy*COS(a*xx)**2*COS(a*yy)**2*4.0D0-D*a**2*xx**4*ym*COS(a*xx&
               &)**2*SIN(a*yy)**2*2.0D0-D*a**2*xx**4*ym*COS(a*yy)**2*SIN(a*xx)**2*&
               &2.0D0+D*a**2*xx**4*yy*COS(a*xx)**2*SIN(a*yy)**2*2.0D0+D*a**2*xx**4&
               &*yy*COS(a*yy)**2*SIN(a*xx)**2*2.0D0-a*mu*xm*xx**3*COS(a*yy)*SIN(a*&
               &xx)*2.4D+1-a*mu*xm**3*xx*COS(a*yy)*SIN(a*xx)*1.0D+1+a*mu*xm*ym**3*&
               &cos(a*xx)*SIN(a*yy)*2.0D0+a*mu*xm**3*ym*COS(a*xx)*SIN(a*yy)*2.0D0-&
               &a*mu*xx*ym**3*COS(a*xx)*SIN(a*yy)*4.0D0-a*mu*xx**3*ym*COS(a*xx)*si&
               &n(a*yy)*4.0D0-a*mu*xm*yy**3*COS(a*xx)*SIN(a*yy)*2.0D0-a*mu*xm**3*y&
               &y*COS(a*xx)*SIN(a*yy)*2.0D0+a*mu*xx*yy**3*COS(a*xx)*SIN(a*yy)*4.0D&
               &0+a*mu*xx**3*yy*COS(a*xx)*SIN(a*yy)*4.0D0-a**2*mu*xx**4*ym*COS(a*x&
               &x)**2*SIN(a*yy)**2*2.0D0-a**2*mu*xx**4*ym*COS(a*yy)**2*SIN(a*xx)**&
               &2*2.0D0+a**2*mu*xx**4*yy*COS(a*xx)**2*SIN(a*yy)**2*2.0D0+a**2*mu*x&
               &x**4*yy*COS(a*yy)**2*SIN(a*xx)**2*2.0D0+a**2*mu*xx**4*ym*SIN(a*xx)&
               &**2*SIN(a*yy)**2*4.0D0-a**2*mu*xx**4*yy*SIN(a*xx)**2*SIN(a*yy)**2*&
               &4.0D0+a*xm*xx**4*COS(a*xx)**2*COS(a*yy)**3*SIN(a*xx)*(1.0D0/xx**2*&
               &(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)-a*xx**4*ym&
               &*COS(a*xx)**3*COS(a*yy)**2*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0&
               &/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+a*xx**4*yy*COS(a*xx)**3*co&
               &s(a*yy)**2*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**&
               &2+1.0D0)**(3.0D0/2.0D0)+a**2*mu*xm**2*xx**2*ym*SIN(a*xx)*SIN(a*yy)&
               &*1.2D+1-a**2*mu*xm**2*xx**2*yy*SIN(a*xx)*SIN(a*yy)*1.2D+1+a**2*mu*&
               &xx**2*ym*yy**2*SIN(a*xx)*SIN(a*yy)*1.2D+1-a**2*mu*xx**2*ym**2*yy*s&
               &in(a*xx)*SIN(a*yy)*1.2D+1-D*a**2*xm*xx*ym**3*COS(a*xx)**2*COS(a*yy&
               &)**2*2.0D0-D*a**2*xm*xx**3*ym*COS(a*xx)**2*COS(a*yy)**2*8.0D0-D*a*&
               &*2*xm**3*xx*ym*COS(a*xx)**2*COS(a*yy)**2*2.0D0+D*a**2*xm*xx*yy**3*&
               &cos(a*xx)**2*COS(a*yy)**2*2.0D0+D*a**2*xm*xx**3*yy*COS(a*xx)**2*co&
               &s(a*yy)**2*8.0D0+D*a**2*xm**3*xx*yy*COS(a*xx)**2*COS(a*yy)**2*2.0D&
               &0+a*mu*xm**2*xx**2*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy)*1.1D+1+a*mu*xm&
               &**2*ym**2*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy)+a*mu*xx**2*ym**2*COS(a*&
               &yy)*SIN(a*xx)**2*SIN(a*yy)*5.0D0+a*mu*xm**2*yy**2*COS(a*yy)*SIN(a*&
               &xx)**2*SIN(a*yy)+a*mu*xx**2*yy**2*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy)&
               &*5.0D0+D*a**2*xm*xx*ym**3*COS(a*xx)**2*SIN(a*yy)**2+D*a**2*xm*xx*y&
               &m**3*COS(a*yy)**2*SIN(a*xx)**2+D*a**2*xm*xx**3*ym*COS(a*xx)**2*sin&
               &(a*yy)**2*4.0D0+D*a**2*xm*xx**3*ym*COS(a*yy)**2*SIN(a*xx)**2*4.0D0&
               &+D*a**2*xm**3*xx*ym*COS(a*xx)**2*SIN(a*yy)**2+D*a**2*xm**3*xx*ym*c&
               &os(a*yy)**2*SIN(a*xx)**2-D*a**2*xm*xx*yy**3*COS(a*xx)**2*SIN(a*yy)&
               &**2-D*a**2*xm*xx*yy**3*COS(a*yy)**2*SIN(a*xx)**2-D*a**2*xm*xx**3*y&
               &y*COS(a*xx)**2*SIN(a*yy)**2*4.0D0-D*a**2*xm*xx**3*yy*COS(a*yy)**2*&
               &sin(a*xx)**2*4.0D0-D*a**2*xm**3*xx*yy*COS(a*xx)**2*SIN(a*yy)**2-D*&
               &a**2*xm**3*xx*yy*COS(a*yy)**2*SIN(a*xx)**2+a*xx**5*COS(a*xx)**2*co&
               &s(a*yy)*SIN(a*xx)*SIN(a*yy)**2*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2&
               &*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0-a*mu*xm*xx*ym**2*COS(a*yy)&
               &*SIN(a*xx)*6.0D0+a*mu*xm*xx**2*ym*COS(a*xx)*SIN(a*yy)*8.0D0-a*mu*x&
               &m**2*xx*ym*COS(a*xx)*SIN(a*yy)*8.0D0-a*mu*xm*xx*yy**2*COS(a*yy)*si&
               &n(a*xx)*6.0D0-a*mu*xm*xx**2*yy*COS(a*xx)*SIN(a*yy)*8.0D0+a*mu*xm**&
               &2*xx*yy*COS(a*xx)*SIN(a*yy)*8.0D0+a*mu*xm*ym*yy**2*COS(a*xx)*SIN(a&
               &*yy)*6.0D0-a*mu*xm*ym**2*yy*COS(a*xx)*SIN(a*yy)*6.0D0-a*mu*xm**2*y&
               &m*yy*COS(a*yy)*SIN(a*xx)*4.0D0-a*mu*xx*ym*yy**2*COS(a*xx)*SIN(a*yy&
               &)*1.2D+1+a*mu*xx*ym**2*yy*COS(a*xx)*SIN(a*yy)*1.2D+1-a*mu*xx**2*ym&
               &*yy*COS(a*yy)*SIN(a*xx)*2.0D+1-xx*ym**3*COS(a*xx)**2*COS(a*yy)**2*&
               &sin(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy&
               &)**2+1.0D0)-xx**3*ym*COS(a*xx)**2*COS(a*yy)**2*SIN(a*xx)*SIN(a*yy)&
               &*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)+xx*yy**&
               &3*COS(a*xx)**2*COS(a*yy)**2*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(&
               &xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)+xx**3*yy*COS(a*xx)**2*COS(&
               &a*yy)**2*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx*&
               &*2*(ym-yy)**2+1.0D0)+a**2*mu*xm*xx*ym**3*COS(a*xx)**2*SIN(a*yy)**2&
               &+a**2*mu*xm*xx*ym**3*COS(a*yy)**2*SIN(a*xx)**2+a**2*mu*xm*xx**3*ym&
               &*COS(a*xx)**2*SIN(a*yy)**2*4.0D0+a**2*mu*xm*xx**3*ym*COS(a*yy)**2*&
               &sin(a*xx)**2*4.0D0+a**2*mu*xm**3*xx*ym*COS(a*xx)**2*SIN(a*yy)**2+a&
               &**2*mu*xm**3*xx*ym*COS(a*yy)**2*SIN(a*xx)**2-a**2*mu*xm*xx*yy**3*c&
               &os(a*xx)**2*SIN(a*yy)**2-a**2*mu*xm*xx*yy**3*COS(a*yy)**2*SIN(a*xx&
               &)**2-a**2*mu*xm*xx**3*yy*COS(a*xx)**2*SIN(a*yy)**2*4.0D0-a**2*mu*x&
               &m*xx**3*yy*COS(a*yy)**2*SIN(a*xx)**2*4.0D0-a**2*mu*xm**3*xx*yy*cos&
               &(a*xx)**2*SIN(a*yy)**2-a**2*mu*xm**3*xx*yy*COS(a*yy)**2*SIN(a*xx)*&
               &*2-a**2*mu*xm*xx*ym**3*SIN(a*xx)**2*SIN(a*yy)**2*2.0D0-a**2*mu*xm*&
               &xx**3*ym*SIN(a*xx)**2*SIN(a*yy)**2*8.0D0-a**2*mu*xm**3*xx*ym*SIN(a&
               &*xx)**2*SIN(a*yy)**2*2.0D0+a**2*mu*xm*xx*yy**3*SIN(a*xx)**2*SIN(a*&
               &yy)**2*2.0D0+a**2*mu*xm*xx**3*yy*SIN(a*xx)**2*SIN(a*yy)**2*8.0D0+a&
               &**2*mu*xm**3*xx*yy*SIN(a*xx)**2*SIN(a*yy)**2*2.0D0+D*a**2*xm**2*xx&
               &**2*ym*COS(a*xx)**2*COS(a*yy)**2*6.0D0-D*a**2*xm**2*xx**2*yy*COS(a&
               &*xx)**2*COS(a*yy)**2*6.0D0+D*a**2*xx**2*ym*yy**2*COS(a*xx)**2*COS(&
               &a*yy)**2*6.0D0-D*a**2*xx**2*ym**2*yy*COS(a*xx)**2*COS(a*yy)**2*6.0&
               &D0-a**2*mu*xx*ym*yy**3*COS(a*xx)*COS(a*yy)*8.0D0-a**2*mu*xx*ym**3*&
               &yy*COS(a*xx)*COS(a*yy)*8.0D0-a**2*mu*xx**3*ym*yy*COS(a*xx)*COS(a*y&
               &y)*2.0D+1-D*a**2*xm**2*xx**2*ym*COS(a*xx)**2*SIN(a*yy)**2*3.0D0-D*&
               &a**2*xm**2*xx**2*ym*COS(a*yy)**2*SIN(a*xx)**2*3.0D0+D*a**2*xm**2*x&
               &x**2*yy*COS(a*xx)**2*SIN(a*yy)**2*3.0D0+D*a**2*xm**2*xx**2*yy*COS(&
               &a*yy)**2*SIN(a*xx)**2*3.0D0-D*a**2*xx**2*ym*yy**2*COS(a*xx)**2*sin&
               &(a*yy)**2*3.0D0-D*a**2*xx**2*ym*yy**2*COS(a*yy)**2*SIN(a*xx)**2*3.&
               &0D0+D*a**2*xx**2*ym**2*yy*COS(a*xx)**2*SIN(a*yy)**2*3.0D0+D*a**2*x&
               &x**2*ym**2*yy*COS(a*yy)**2*SIN(a*xx)**2*3.0D0-D*a*xm*ym**3*COS(a*x&
               &x)*COS(a*yy)**2*SIN(a*xx)-D*a*xm**3*ym*COS(a*xx)*COS(a*yy)**2*SIN(&
               &a*xx)+D*a*xx*ym**3*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*2.0D0+D*a*xx**&
               &3*ym*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*2.0D0+D*a*xm*xx**3*COS(a*xx)&
               &**2*COS(a*yy)*SIN(a*yy)*1.2D+1+D*a*xm*yy**3*COS(a*xx)*COS(a*yy)**2&
               &*SIN(a*xx)+D*a*xm**3*xx*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*5.0D0+D*a&
               &*xm**3*yy*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)-D*a*xx*yy**3*COS(a*xx)*&
               &cos(a*yy)**2*SIN(a*xx)*2.0D0-D*a*xx**3*yy*COS(a*xx)*COS(a*yy)**2*s&
               &in(a*xx)*2.0D0+a*xx**4*ym*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*(1.0D0/&
               &xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*4.0D&
               &0-a*xm*xx**4*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)&
               &**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*4.0D0-a*xx**4*yy*&
               &cos(a*xx)*COS(a*yy)**2*SIN(a*xx)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx*&
               &*2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*4.0D0-a**2*mu*xm*xx*ym**3*SIN(&
               &a*xx)*SIN(a*yy)*4.0D0-a**2*mu*xm*xx**3*ym*SIN(a*xx)*SIN(a*yy)*1.6D&
               &+1-a**2*mu*xm**3*xx*ym*SIN(a*xx)*SIN(a*yy)*4.0D0+a**2*mu*xm*xx*yy*&
               &*3*SIN(a*xx)*SIN(a*yy)*4.0D0+a**2*mu*xm*xx**3*yy*SIN(a*xx)*SIN(a*y&
               &y)*1.6D+1+a**2*mu*xm**3*xx*yy*SIN(a*xx)*SIN(a*yy)*4.0D0-a**2*mu*xm&
               &**2*xx**2*ym*COS(a*xx)**2*SIN(a*yy)**2*3.0D0-a**2*mu*xm**2*xx**2*y&
               &m*COS(a*yy)**2*SIN(a*xx)**2*3.0D0+a**2*mu*xm**2*xx**2*yy*COS(a*xx)&
               &**2*SIN(a*yy)**2*3.0D0+a**2*mu*xm**2*xx**2*yy*COS(a*yy)**2*SIN(a*x&
               &x)**2*3.0D0-a**2*mu*xx**2*ym*yy**2*COS(a*xx)**2*SIN(a*yy)**2*3.0D0&
               &-a**2*mu*xx**2*ym*yy**2*COS(a*yy)**2*SIN(a*xx)**2*3.0D0+a**2*mu*xx&
               &**2*ym**2*yy*COS(a*xx)**2*SIN(a*yy)**2*3.0D0+a**2*mu*xx**2*ym**2*y&
               &y*COS(a*yy)**2*SIN(a*xx)**2*3.0D0+a**2*mu*xm**2*xx**2*ym*SIN(a*xx)&
               &**2*SIN(a*yy)**2*6.0D0-a**2*mu*xm**2*xx**2*yy*SIN(a*xx)**2*SIN(a*y&
               &y)**2*6.0D0+D*a**2*xx**5*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*1&
               &.2D+1+a**2*mu*xx**2*ym*yy**2*SIN(a*xx)**2*SIN(a*yy)**2*6.0D0-a**2*&
               &mu*xx**2*ym**2*yy*SIN(a*xx)**2*SIN(a*yy)**2*6.0D0+a*mu*xm*ym**3*co&
               &s(a*xx)*SIN(a*xx)*SIN(a*yy)**2+a*mu*xm**3*ym*COS(a*xx)*SIN(a*xx)*s&
               &in(a*yy)**2-a*mu*xx*ym**3*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*2.0D0-a&
               &*mu*xx**3*ym*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*2.0D0-a*mu*xm*xx**3*&
               &cos(a*yy)*SIN(a*xx)**2*SIN(a*yy)*1.2D+1-a*mu*xm*yy**3*COS(a*xx)*si&
               &n(a*xx)*SIN(a*yy)**2-a*mu*xm**3*xx*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy&
               &)*5.0D0-a*mu*xm**3*yy*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2+a*mu*xx*yy*&
               &*3*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*2.0D0+a*mu*xx**3*yy*COS(a*xx)*&
               &sin(a*xx)*SIN(a*yy)**2*2.0D0-a**2*mu*xm*xx**2*ym**2*COS(a*xx)*COS(&
               &a*yy)*8.0D0+a**2*mu*xm**2*xx*ym**2*COS(a*xx)*COS(a*yy)*4.0D0-a**2*&
               &mu*xm*xx**2*yy**2*COS(a*xx)*COS(a*yy)*8.0D0+a**2*mu*xm**2*xx*yy**2&
               &*COS(a*xx)*COS(a*yy)*4.0D0+a**2*mu*xx*ym**2*yy**2*COS(a*xx)*COS(a*&
               &yy)*1.2D+1+Mref*a*xm*xx**4*COS(a*yy)*SIN(a*xx)*(1.0D0/xx**2*(xm-xx&
               &)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0-Mref*a*xx*&
               &*4*ym*COS(a*xx)*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-&
               &yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0+Mref*a*xx**4*yy*COS(a*xx)*SIN(a&
               &*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0&
               &/2.0D0)*2.0D0+a**2*mu*xx**5*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy&
               &)*1.2D+1-D*a*xm**2*xx**2*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*1.1D+1-D&
               &*a*xm**2*ym**2*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)-D*a*xx**2*ym**2*co&
               &s(a*xx)**2*COS(a*yy)*SIN(a*yy)*5.0D0-D*a*xm**2*yy**2*COS(a*xx)**2*&
               &cos(a*yy)*SIN(a*yy)-D*a*xx**2*yy**2*COS(a*xx)**2*COS(a*yy)*SIN(a*y&
               &y)*5.0D0+a**2*mu*xm*xx**2*ym*yy*COS(a*xx)*COS(a*yy)*1.6D+1-a**2*mu&
               &*xm**2*xx*ym*yy*COS(a*xx)*COS(a*yy)*8.0D0-D*a*xm*xx**2*ym*COS(a*xx&
               &)*COS(a*yy)**2*SIN(a*xx)*4.0D0+D*a*xm**2*xx*ym*COS(a*xx)*COS(a*yy)&
               &**2*SIN(a*xx)*4.0D0+D*a*xm*xx**2*yy*COS(a*xx)*COS(a*yy)**2*SIN(a*x&
               &x)*4.0D0-D*a*xm**2*xx*yy*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*4.0D0+D*&
               &a*xm*xx*ym**2*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*3.0D0-D*a*xm*ym*yy*&
               &*2*COS(a*xx)*COS(a*yy)**2*SIN(a*xx)*3.0D0+D*a*xm*ym**2*yy*COS(a*xx&
               &)*COS(a*yy)**2*SIN(a*xx)*3.0D0+D*a*xx*ym*yy**2*COS(a*xx)*COS(a*yy)&
               &**2*SIN(a*xx)*6.0D0-D*a*xx*ym**2*yy*COS(a*xx)*COS(a*yy)**2*SIN(a*x&
               &x)*6.0D0+D*a*xm*xx*yy**2*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*3.0D0+D*&
               &a*xm**2*ym*yy*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*2.0D0+D*a*xx**2*ym*&
               &yy*COS(a*xx)**2*COS(a*yy)*SIN(a*yy)*1.0D+1-a**2*mu*xm*xx*ym*yy**2*&
               &sin(a*xx)*SIN(a*yy)*1.2D+1+a**2*mu*xm*xx*ym**2*yy*SIN(a*xx)*SIN(a*&
               &yy)*1.2D+1-D*a**2*xm*xx**4*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)&
               &*2.0D+1+D*a**2*xm**4*xx*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*2.&
               &0D0+D*a**2*xx*ym**4*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*2.0D0+&
               &D*a**2*xx*yy**4*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*2.0D0+a*mu&
               &*xm*xx**2*ym*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*4.0D0-a*mu*xm**2*xx*&
               &ym*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*4.0D0-a*mu*xm*xx**2*yy*COS(a*x&
               &x)*SIN(a*xx)*SIN(a*yy)**2*4.0D0+a*mu*xm**2*xx*yy*COS(a*xx)*SIN(a*x&
               &x)*SIN(a*yy)**2*4.0D0-a*mu*xm*xx*ym**2*COS(a*yy)*SIN(a*xx)**2*SIN(&
               &a*yy)*3.0D0+a*mu*xm*ym*yy**2*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*3.0D&
               &0-a*mu*xm*ym**2*yy*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*3.0D0-a*mu*xx*&
               &ym*yy**2*COS(a*xx)*SIN(a*xx)*SIN(a*yy)**2*6.0D0+a*mu*xx*ym**2*yy*c&
               &os(a*xx)*SIN(a*xx)*SIN(a*yy)**2*6.0D0-a*mu*xm*xx*yy**2*COS(a*yy)*s&
               &in(a*xx)**2*SIN(a*yy)*3.0D0-a*mu*xm**2*ym*yy*COS(a*yy)*SIN(a*xx)**&
               &2*SIN(a*yy)*2.0D0-a*mu*xx**2*ym*yy*COS(a*yy)*SIN(a*xx)**2*SIN(a*yy&
               &)*1.0D+1+a*mu*xm*xx*ym*yy*COS(a*yy)*SIN(a*xx)*1.2D+1-a**2*mu*xm*xx&
               &**4*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*2.0D+1+a**2*mu*xm**4*x&
               &x*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*2.0D0+a**2*mu*xx*ym**4*c&
               &os(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*2.0D0+a**2*mu*xx*yy**4*COS(&
               &a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*2.0D0-D*a**2*xm*xx*ym*yy**2*co&
               &s(a*xx)**2*COS(a*yy)**2*6.0D0+D*a**2*xm*xx*ym**2*yy*COS(a*xx)**2*c&
               &os(a*yy)**2*6.0D0+D*a**2*xm**2*xx**3*COS(a*xx)*COS(a*yy)*SIN(a*xx)&
               &*SIN(a*yy)*1.8D+1-D*a**2*xm**3*xx**2*COS(a*xx)*COS(a*yy)*SIN(a*xx)&
               &*SIN(a*yy)*8.0D0+D*a**2*xx**3*ym**2*COS(a*xx)*COS(a*yy)*SIN(a*xx)*&
               &sin(a*yy)*1.0D+1+D*a**2*xx**3*yy**2*COS(a*xx)*COS(a*yy)*SIN(a*xx)*&
               &sin(a*yy)*1.0D+1+D*a**2*xm*xx*ym*yy**2*COS(a*xx)**2*SIN(a*yy)**2*3&
               &.0D0+D*a**2*xm*xx*ym*yy**2*COS(a*yy)**2*SIN(a*xx)**2*3.0D0-D*a**2*&
               &xm*xx*ym**2*yy*COS(a*xx)**2*SIN(a*yy)**2*3.0D0-D*a**2*xm*xx*ym**2*&
               &yy*COS(a*yy)**2*SIN(a*xx)**2*3.0D0-a*xm*xx**4*COS(a*xx)**2*COS(a*y&
               &y)*SIN(a*xx)*SIN(a*yy)**2*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-&
               &yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0+a*xx**4*ym*COS(a*xx)*COS(a*yy)*&
               &*2*SIN(a*xx)**2*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-&
               &yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0-a*xx**4*yy*COS(a*xx)*COS(a*yy)*&
               &*2*SIN(a*xx)**2*SIN(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-&
               &yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0+xm*xx**2*ym*COS(a*xx)**2*COS(a*&
               &yy)**2*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2&
               &*(ym-yy)**2+1.0D0)*2.0D0-xm**2*xx*ym*COS(a*xx)**2*COS(a*yy)**2*sin&
               &(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**&
               &2+1.0D0)-xm*xx**2*yy*COS(a*xx)**2*COS(a*yy)**2*SIN(a*xx)*SIN(a*yy)&
               &*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0+x&
               &m**2*xx*yy*COS(a*xx)**2*COS(a*yy)**2*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D&
               &0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)-xx*ym*yy**2*COS(a&
               &*xx)**2*COS(a*yy)**2*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)*&
               &*2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*3.0D0+xx*ym**2*yy*COS(a*xx)**2*co&
               &s(a*yy)**2*SIN(a*xx)*SIN(a*yy)*SQRT(1.0D0/xx**2*(xm-xx)**2+1.0D0/x&
               &x**2*(ym-yy)**2+1.0D0)*3.0D0+a**2*mu*xm**2*xx**3*COS(a*xx)*COS(a*y&
               &y)*SIN(a*xx)*SIN(a*yy)*1.8D+1-a**2*mu*xm**3*xx**2*COS(a*xx)*COS(a*&
               &yy)*SIN(a*xx)*SIN(a*yy)*8.0D0+a**2*mu*xx**3*ym**2*COS(a*xx)*COS(a*&
               &yy)*SIN(a*xx)*SIN(a*yy)*1.0D+1+a**2*mu*xx**3*yy**2*COS(a*xx)*COS(a&
               &*yy)*SIN(a*xx)*SIN(a*yy)*1.0D+1+a**2*mu*xm*xx*ym*yy**2*COS(a*xx)**&
               &2*SIN(a*yy)**2*3.0D0+a**2*mu*xm*xx*ym*yy**2*COS(a*yy)**2*SIN(a*xx)&
               &**2*3.0D0-a**2*mu*xm*xx*ym**2*yy*COS(a*xx)**2*SIN(a*yy)**2*3.0D0-a&
               &**2*mu*xm*xx*ym**2*yy*COS(a*yy)**2*SIN(a*xx)**2*3.0D0-a**2*mu*xm*x&
               &x*ym*yy**2*SIN(a*xx)**2*SIN(a*yy)**2*6.0D0+a**2*mu*xm*xx*ym**2*yy*&
               &sin(a*xx)**2*SIN(a*yy)**2*6.0D0-D*a**2*xm*xx**2*ym**2*COS(a*xx)*co&
               &s(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0+D*a**2*xm**2*xx*ym**2*COS(a*xx)*&
               &cos(a*yy)*SIN(a*xx)*SIN(a*yy)*4.0D0-D*a**2*xm*xx**2*yy**2*COS(a*xx&
               &)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0+D*a**2*xm**2*xx*yy**2*COS(a*&
               &xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*4.0D0+D*a**2*xx*ym**2*yy**2*COS(&
               &a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*1.2D+1-D*a*xm*xx*ym*yy*COS(a*x&
               &x)**2*COS(a*yy)*SIN(a*yy)*6.0D0-a**2*mu*xm*xx**2*ym**2*COS(a*xx)*c&
               &os(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0+a**2*mu*xm**2*xx*ym**2*COS(a*xx&
               &)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*4.0D0-a**2*mu*xm*xx**2*yy**2*COS(a&
               &*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0+a**2*mu*xm**2*xx*yy**2*co&
               &s(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*4.0D0+a**2*mu*xx*ym**2*yy**2&
               &*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*1.2D+1+a*mu*xm*xx*ym*yy*c&
               &os(a*yy)*SIN(a*xx)**2*SIN(a*yy)*6.0D0-D*a**2*xx*ym*yy**3*COS(a*xx)&
               &*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0-D*a**2*xx*ym**3*yy*COS(a*xx)*&
               &cos(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0-D*a**2*xx**3*ym*yy*COS(a*xx)*c&
               &os(a*yy)*SIN(a*xx)*SIN(a*yy)*2.0D+1-a**2*mu*xx*ym*yy**3*COS(a*xx)*&
               &cos(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0-a**2*mu*xx*ym**3*yy*COS(a*xx)*&
               &cos(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0-a**2*mu*xx**3*ym*yy*COS(a*xx)*&
               &cos(a*yy)*SIN(a*xx)*SIN(a*yy)*2.0D+1+D*a**2*xm*xx**2*ym*yy*COS(a*x&
               &x)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*1.6D+1-D*a**2*xm**2*xx*ym*yy*COS(&
               &a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0+a**2*mu*xm*xx**2*ym*yy*c&
               &os(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*1.6D+1-a**2*mu*xm**2*xx*ym*&
               &yy*COS(a*xx)*COS(a*yy)*SIN(a*xx)*SIN(a*yy)*8.0D0))/xx

       END DO
    CASE(5:6)
       ! Do nothing
    CASE (50:)
       !Do nothing
    CASE DEFAULT
       WRITE (6, *) "Error! Test case not valid"
       STOP

    END SELECT
  END SUBROUTINE body_force
#endif

END MODULE analytical
