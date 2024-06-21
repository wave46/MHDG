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
          CASE (2)
             ! Axisimmetric case with div(b)~=0
             IF (.NOT. switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             up(ind, 1) = 2 + SIN(a*xx)*SIN(a*yy)
          CASE (50:64)
             up(ind, 1) = 1.
          CASE (65)
             up(ind, 1) = 1.
             r = SQRT((xx*phys%lscale - geom%R0)**2 + (yy*phys%lscale - 0.75)**2)
             IF (r .LE. 0.05) THEN
                up(ind, 2) = 1.
             END IF
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

             upy(ind, 1) = a*SIN(a*xx)*COS(a*yy)

          CASE (2)
             ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
             IF (.NOT. switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             upx(ind, 1) = a*COS(a*xx)*SIN(a*yy)

             upy(ind, 1) = a*SIN(a*xx)*COS(a*yy)
          CASE (3)
             ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
             IF (.NOT. switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             upt(ind, 1) = +a*COS(a*tt)
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
    INTEGER:: i
    REAL*8 :: a, r(SIZE(x))

    up = 0.
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
    CASE (2)
       ! Axisimmetric case with div(b)~=0
       IF (.NOT. switch%axisym) THEN
          WRITE (6, *) "This is an axisymmetric test case!"
          STOP
       END IF
       up(:, 1) = 2 + SIN(a*x)*SIN(a*y)
    CASE (50:64)
       up(:, 1) = 1.
    CASE (65)
       up(:, 1) = 1.
       r = SQRT((x*phys%lscale - geom%R0)**2 + (y*phys%lscale - 0.75)**2)
       DO i = 1, SIZE(x)
          IF (r(i) .LE. 0.05) THEN
             up(i, 2) = 1.
          END IF
       END DO
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

       upy(:, 1) = a*SIN(a*x)*COS(a*y)

    CASE (2)
       ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
       IF (.NOT. switch%axisym) THEN
          WRITE (6, *) "This is an axisymmetric test case!"
          STOP
       END IF
       upx(:, 1) = a*COS(a*x)*SIN(a*y)

       upy(:, 1) = a*SIN(a*x)*COS(a*y)
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
    ux(:, 1) = upx(:, 1)
    uy(:, 1) = upy(:, 1)

  END SUBROUTINE analytical_gradient

  !*****************************************
  ! Body forces
  !****************************************
  SUBROUTINE body_force(x, y, f)
    REAL*8, DIMENSION(:), INTENT(IN) :: x, y
    REAL*8, DIMENSION(:, :), INTENT(OUT) :: f
    INTEGER                            :: i, n
    REAL*8  :: a, b, xc, yc, D, mu, k

    n = SIZE(x)

    f = 0.
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
       f(:,1) = COS(a * x)**2 * a  * SIN(b * y) * COS(b * y) * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+SIN(a * x) * SIN(b * y)) * SIN(a * x) * a  * COS(b * y) * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+SIN(a * x) * SIN(b * y)) * COS(a * x) * COS(b * y) * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * x-0.2e1 * xc) / 0.20e2+SIN(a * x) * COS(b * y)**2 * b * COS(a * x) * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+SIN(a * x) * SIN(b * y)) * COS(a * x) * SIN(b * y) * b * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+SIN(a * x) * SIN(b * y)) * COS(a * x) * COS(b * y) * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * y-0.2e1 * yc) / 0.20e2-D * (-SIN(a * x) * a**2 * SIN(b * y)+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * COS(a * x) * a  * SIN(b * y) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * SIN(a * x) * COS(b * y) * b / 0.10e2) * (0.2e1 * x-0.2e1 * xc) / 0.20e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * COS(a * x) * a  * SIN(b * y) * (0.2e1 * x-0.2e1 * xc) / 0.20e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * SIN(a * x) * a**2 * SIN(b * y) / 0.10e2-(y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * SIN(a * x) * COS(b * y) * b / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * SIN(a * x) * COS(b * y) * b * (0.2e1 * x-0.2e1 * xc) / 0.20e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * COS(a * x) * a  * COS(b * y) * b / 0.10e2) / 0.10e2-SIN(a * x) * SIN(b * y) * b**2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * COS(a * x) * a  * SIN(b * y) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * SIN(a * x) * COS(b * y) * b / 0.10e2) * (0.2e1 * y-0.2e1 * yc) / 0.20e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * ((y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * COS(a * x) * a  * SIN(b * y) / 0.10e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * COS(a * x) * a  * SIN(b * y) * (0.2e1 * y-0.2e1 * yc) / 0.20e2+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * COS(a * x) * a  * COS(b * y) * b / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * SIN(a * x) * COS(b * y) * b * (0.2e1 * y-0.2e1 * yc) / 0.20e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * SIN(a * x) * SIN(b * y) * b**2 / 0.10e2) / 0.10e2)
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
       f(:, 1) = 0.
    CASE (50:)
       !Do nothing
    CASE DEFAULT
       WRITE (6, *) "Error! Test case not valid"
       STOP

    END SELECT
  END SUBROUTINE body_force
#endif

END MODULE analytical
