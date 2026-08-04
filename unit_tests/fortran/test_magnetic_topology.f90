PROGRAM test_magnetic_topology
  USE magnetic_topology
  IMPLICIT NONE

  CALL test_limited_wall_contact()
  CALL test_limited_sharp_corner_contact()
  CALL test_high_order_wall_and_grid_resolution()
  CALL test_lower_single_null_regions()
  WRITE (*, '(A)') 'magnetic topology analytical checks: PASS'

CONTAINS

  SUBROUTINE test_limited_wall_contact()
    TYPE(equilibrium_geometry_t) :: geometry
    INTEGER, PARAMETER :: nr = 81, nz = 81
    REAL*8, PARAMETER :: r0 = 2.d0, a = 0.5d0
    REAL*8 :: r(nr), z(nz), psi(nz, nr)
    REAL*8 :: wall(4, 2), ref_nodes(2)
    INTEGER :: faces(4, 2), ierr, ir, iz, region
    REAL*8 :: psi_n, rho, normal(2)
    CHARACTER(LEN=256) :: message

    DO ir = 1, nr
       r(ir) = r0 - 0.75d0 + 1.5d0*REAL(ir - 1)/REAL(nr - 1)
    ENDDO
    DO iz = 1, nz
       z(iz) = -0.75d0 + 1.5d0*REAL(iz - 1)/REAL(nz - 1)
    ENDDO
    DO iz = 1, nz
       DO ir = 1, nr
          psi(iz, ir) = (r(ir) - r0)**2 + z(iz)**2
       ENDDO
    ENDDO

    wall(1, :) = (/r0 - a, -0.6d0/)
    wall(2, :) = (/r0 + 0.65d0, -0.6d0/)
    wall(3, :) = (/r0 + 0.65d0, 0.6d0/)
    wall(4, :) = (/r0 - a, 0.6d0/)
    faces(1, :) = (/1, 2/)
    faces(2, :) = (/2, 3/)
    faces(3, :) = (/3, 4/)
    faces(4, :) = (/4, 1/)
    ref_nodes = (/-1.d0, 1.d0/)

    CALL geometry%init(r, z, psi, ierr, message)
    CALL assert_ok(ierr, message)
    CALL geometry%analyze(1.03d0*a*a, wall, faces, ref_nodes, ierr, message)
    CALL assert_ok(ierr, message)
    CALL assert_true(geometry%topology_kind == topology_limited, &
         'circular wall-contact case was not classified as limited')
    CALL assert_close(geometry%r_axis, r0, 2.d-6, 'limited magnetic-axis R')
    CALL assert_close(geometry%z_axis, 0.d0, 2.d-6, 'limited magnetic-axis Z')
    CALL assert_close(geometry%psi_lcfs, a*a, 2.d-6, 'wall-derived LCFS flux')
    CALL assert_close(geometry%a_minor, a, 2.d-5, 'wall-derived minor radius')
    CALL assert_true(geometry%lcfs_extrema_refined, &
         'limited LCFS extrema were not refined')

    CALL geometry%evaluate_topology(r0 + 0.25d0, 0.d0, psi_n, rho, normal, region)
    CALL assert_true(region == magnetic_region_core, 'limited core classification')
    CALL assert_close(rho, 0.5d0, 2.d-5, 'limited normalized poloidal radius')
    CALL assert_true(normal(1) > 0.999d0, 'limited outward normal')
    CALL geometry%evaluate_topology(r0 + 0.6d0, 0.d0, psi_n, rho, normal, region)
    CALL assert_true(region == magnetic_region_main_sol, 'limited SOL classification')

    psi = -psi
    CALL geometry%init(r, z, psi, ierr, message)
    CALL assert_ok(ierr, message)
    CALL geometry%analyze(-1.03d0*a*a, wall, faces, ref_nodes, ierr, message)
    CALL assert_ok(ierr, message)
    CALL geometry%evaluate_topology(r0 + 0.25d0, 0.d0, psi_n, rho, normal, region)
    CALL assert_close(rho, 0.5d0, 2.d-5, 'reversed-flux normalized poloidal radius')
    CALL assert_true(normal(1) > 0.999d0, 'reversed-flux outward normal')
  END SUBROUTINE test_limited_wall_contact

  SUBROUTINE test_limited_sharp_corner_contact()
    TYPE(equilibrium_geometry_t) :: geometry
    INTEGER, PARAMETER :: nr = 41, nz = 41
    REAL*8, PARAMETER :: r0 = 2.d0, a = 0.5d0
    REAL*8 :: r(nr), z(nz), psi(nz, nr), wall(7, 2), ref_nodes(2)
    INTEGER :: faces(7, 2), ierr, ir, iz
    CHARACTER(LEN=256) :: message

    CALL fill_circular_flux(r0, nr, nz, r, z, psi)
    wall(1, :) = (/r0 - 0.8d0, -0.7d0/)
    wall(2, :) = (/r0 + 0.8d0, -0.7d0/)
    wall(3, :) = (/r0 + 0.8d0,  0.7d0/)
    wall(4, :) = (/r0 - 0.8d0,  0.7d0/)
    wall(5, :) = (/r0 - 0.8d0,  0.2d0/)
    wall(6, :) = (/r0 - a,       0.d0/)
    wall(7, :) = (/r0 - 0.8d0, -0.2d0/)
    DO ir = 1, 6
       faces(ir, :) = (/ir, ir + 1/)
    ENDDO
    faces(7, :) = (/7, 1/)
    ref_nodes = (/-1.d0, 1.d0/)

    CALL geometry%init(r, z, psi, ierr, message)
    CALL assert_ok(ierr, message)
    CALL geometry%analyze(0.96d0*a*a, wall, faces, ref_nodes, ierr, message)
    CALL assert_ok(ierr, message)
    CALL assert_true(geometry%topology_kind == topology_limited, &
         'sharp-corner limiter was not classified as limited')
    CALL assert_close(geometry%psi_lcfs, a*a, 2.d-6, &
         'sharp-corner wall-derived LCFS flux')
    CALL assert_close(geometry%r_wall_contact, r0 - a, 2.d-6, &
         'sharp-corner contact R')
    CALL assert_close(geometry%z_wall_contact, 0.d0, 2.d-6, &
         'sharp-corner contact Z')
  END SUBROUTINE test_limited_sharp_corner_contact

  SUBROUTINE test_high_order_wall_and_grid_resolution()
    TYPE(equilibrium_geometry_t) :: coarse, fine
    INTEGER, PARAMETER :: nr_coarse = 17, nz_coarse = 19
    INTEGER, PARAMETER :: nr_fine = 83, nz_fine = 79
    REAL*8, PARAMETER :: r0 = 2.d0, a = 0.5d0
    REAL*8 :: rc(nr_coarse), zc(nz_coarse), psic(nz_coarse, nr_coarse)
    REAL*8 :: rf(nr_fine), zf(nz_fine), psif(nz_fine, nr_fine)
    REAL*8 :: wall(8, 2), ref_nodes(3)
    INTEGER :: faces(4, 3), ierr
    CHARACTER(LEN=256) :: message

    CALL fill_circular_flux(r0, nr_coarse, nz_coarse, rc, zc, psic)
    CALL fill_circular_flux(r0, nr_fine, nz_fine, rf, zf, psif)

    wall(1, :) = (/r0 - a - 0.2d0, -0.6d0/)
    wall(2, :) = (/r0 + 0.7d0,      -0.6d0/)
    wall(3, :) = (/r0 + 0.7d0,       0.6d0/)
    wall(4, :) = (/r0 - a - 0.2d0,  0.6d0/)
    wall(5, :) = (/r0 + 0.1d0,      -0.6d0/)
    wall(6, :) = (/r0 + 0.7d0,       0.d0/)
    wall(7, :) = (/r0 + 0.1d0,       0.6d0/)
    wall(8, :) = (/r0 - a,            0.d0/)
    faces(1, :) = (/1, 5, 2/)
    faces(2, :) = (/2, 6, 3/)
    faces(3, :) = (/3, 7, 4/)
    faces(4, :) = (/4, 8, 1/)
    ref_nodes = (/-1.d0, 0.d0, 1.d0/)

    CALL coarse%init(rc, zc, psic, ierr, message)
    CALL assert_ok(ierr, message)
    CALL coarse%analyze(1.04d0*a*a, wall, faces, ref_nodes, ierr, message)
    CALL assert_ok(ierr, message)
    CALL fine%init(rf, zf, psif, ierr, message)
    CALL assert_ok(ierr, message)
    CALL fine%analyze(0.97d0*a*a, wall, faces, ref_nodes, ierr, message)
    CALL assert_ok(ierr, message)

    CALL assert_close(coarse%psi_lcfs, a*a, 3.d-6, &
         'quadratic-wall coarse-grid LCFS flux')
    CALL assert_close(fine%psi_lcfs, a*a, 3.d-6, &
         'quadratic-wall fine-grid LCFS flux')
    CALL assert_close(coarse%psi_lcfs, fine%psi_lcfs, 2.d-7, &
         'coarse/fine equilibrium LCFS consistency')
    CALL assert_close(coarse%r_wall_contact, r0 - a, 3.d-6, &
         'quadratic high-order wall contact R')
    CALL assert_close(coarse%z_wall_contact, 0.d0, 3.d-6, &
         'quadratic high-order wall contact Z')
  END SUBROUTINE test_high_order_wall_and_grid_resolution

  SUBROUTINE test_lower_single_null_regions()
    TYPE(equilibrium_geometry_t) :: geometry
    INTEGER, PARAMETER :: nr = 81, nz = 121
    REAL*8, PARAMETER :: r0 = 2.d0, b = 0.6d0, k = 1.d0
    REAL*8, PARAMETER :: psi_x = k*b**3/6.d0
    REAL*8 :: r(nr), z(nz), psi(nz, nr)
    REAL*8 :: wall(4, 2), ref_nodes(2)
    INTEGER :: faces(4, 2), ierr, ir, iz, region
    REAL*8 :: psi_n, rho, normal(2), psi_n_core, rho_core, target_psi
    CHARACTER(LEN=256) :: message

    DO ir = 1, nr
       r(ir) = 1.2d0 + 1.6d0*REAL(ir - 1)/REAL(nr - 1)
    ENDDO
    DO iz = 1, nz
       z(iz) = -1.2d0 + 2.4d0*REAL(iz - 1)/REAL(nz - 1)
    ENDDO
    DO iz = 1, nz
       DO ir = 1, nr
          psi(iz, ir) = (r(ir) - r0)**2 + &
               k*(z(iz)**3/3.d0 + b*z(iz)**2/2.d0)
       ENDDO
    ENDDO

    wall(1, :) = (/1.35d0, -1.d0/)
    wall(2, :) = (/2.65d0, -1.d0/)
    wall(3, :) = (/2.65d0, 0.9d0/)
    wall(4, :) = (/1.35d0, 0.9d0/)
    faces(1, :) = (/1, 2/)
    faces(2, :) = (/2, 3/)
    faces(3, :) = (/3, 4/)
    faces(4, :) = (/4, 1/)
    ref_nodes = (/-1.d0, 1.d0/)

    CALL geometry%init(r, z, psi, ierr, message)
    CALL assert_ok(ierr, message)
    CALL geometry%analyze(1.20d0*psi_x, wall, faces, ref_nodes, ierr, message)
    CALL assert_ok(ierr, message)
    CALL assert_true(geometry%topology_kind == topology_lower_single_null, &
         'polynomial diverted case was not classified as lower single null')
    CALL assert_close(geometry%r_axis, r0, 2.d-3, 'diverted magnetic-axis R')
    CALL assert_close(geometry%z_axis, 0.d0, 2.d-3, 'diverted magnetic-axis Z')
    CALL assert_close(geometry%r_xpoint, r0, 2.d-3, 'lower X-point R')
    CALL assert_close(geometry%z_xpoint, -b, 2.d-3, 'lower X-point Z')
    CALL assert_close(geometry%psi_lcfs, psi_x, 4.d-4, 'X-point LCFS flux')
    CALL assert_close(geometry%a_minor, SQRT(psi_x), 2.d-3, &
         'diverted LCFS half width')
    CALL assert_true(geometry%lcfs_extrema_refined, &
         'diverted LCFS extrema were not refined')

    CALL geometry%evaluate_topology(r0, 0.15d0, psi_n, rho, normal, region)
    CALL assert_true(region == magnetic_region_core, 'diverted core classification')
    target_psi = k*(0.2d0**3/3.d0 + b*0.2d0**2/2.d0)
    CALL geometry%evaluate_topology(r0, 0.2d0, psi_n_core, rho_core, normal, region)
    CALL assert_true(region == magnetic_region_core, 'equal-rho core classification')
    CALL geometry%evaluate_topology(r0 + SQRT(target_psi), -0.9d0, psi_n, rho, normal, region)
    CALL assert_true(psi_n < 1.d0, 'private-flux fixture must have psi_n below one')
    CALL assert_close(rho, rho_core, 2.d-3, 'equal-rho core/PFR values')
    CALL assert_true(region == magnetic_region_private_flux, &
         'disconnected private-flux classification')
    CALL geometry%evaluate_topology(r0 + 0.5d0, 0.d0, psi_n, rho, normal, region)
    CALL assert_true(psi_n > 1.d0, 'main-SOL fixture must have rho above one')
    CALL assert_true(region == magnetic_region_main_sol, 'main-SOL classification')
  END SUBROUTINE test_lower_single_null_regions

  SUBROUTINE fill_circular_flux(r0, nr, nz, r, z, psi)
    REAL*8, INTENT(IN) :: r0
    INTEGER, INTENT(IN) :: nr, nz
    REAL*8, INTENT(OUT) :: r(nr), z(nz), psi(nz, nr)
    INTEGER :: ir, iz

    DO ir = 1, nr
       r(ir) = r0 - 0.9d0 + 1.8d0*REAL(ir - 1)/REAL(nr - 1)
    ENDDO
    DO iz = 1, nz
       z(iz) = -0.8d0 + 1.6d0*REAL(iz - 1)/REAL(nz - 1)
    ENDDO
    DO iz = 1, nz
       DO ir = 1, nr
          psi(iz, ir) = (r(ir) - r0)**2 + z(iz)**2
       ENDDO
    ENDDO
  END SUBROUTINE fill_circular_flux

  SUBROUTINE assert_ok(ierr, message)
    INTEGER, INTENT(IN) :: ierr
    CHARACTER(LEN=*), INTENT(IN) :: message
    IF (ierr /= 0) THEN
       WRITE (*, '(A,I0,2A)') 'FAIL(', ierr, '): ', TRIM(message)
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_ok

  SUBROUTINE assert_close(value, expected, tolerance, label)
    REAL*8, INTENT(IN) :: value, expected, tolerance
    CHARACTER(LEN=*), INTENT(IN) :: label
    IF (ABS(value - expected) > tolerance) THEN
       WRITE (*, '(3A,3ES16.7)') 'FAIL: ', TRIM(label), ': ', value, expected, tolerance
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_close

  SUBROUTINE assert_true(condition, label)
    LOGICAL, INTENT(IN) :: condition
    CHARACTER(LEN=*), INTENT(IN) :: label
    IF (.NOT. condition) THEN
       WRITE (*, '(2A)') 'FAIL: ', TRIM(label)
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_true

END PROGRAM test_magnetic_topology
