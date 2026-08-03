PROGRAM test_magnetic_geometry_state
  USE magnetic_topology
  USE magnetic_geometry_state
  IMPLICIT NONE

  CALL test_signed_field_fit()
  CALL test_exact_cache_and_refresh()
  WRITE (*, '(A)') 'magnetic geometry state checks: PASS'

CONTAINS

  SUBROUTINE test_signed_field_fit()
    TYPE(equilibrium_geometry_t) :: geometry
    TYPE(poloidal_field_fit_t) :: fit
    INTEGER, PARAMETER :: nr = 31, nz = 29
    REAL*8, PARAMETER :: r0 = 2.d0, alpha = -1.d0/(2.d0*ACOS(-1.d0))
    REAL*8 :: r(nr), z(nz), psi(nz, nr), br(nz, nr), bz(nz, nr)
    REAL*8 :: psi_value, psi_r, psi_z, perturbation
    INTEGER :: ir, iz, ierr
    CHARACTER(LEN=256) :: message

    CALL fill_circular_flux(r0, r, z, psi)
    CALL geometry%init(r, z, psi, ierr, message)
    CALL assert_ok(ierr, message)
    DO ir = 1, nr
       DO iz = 1, nz
          CALL geometry%evaluate_flux(r(ir), z(iz), psi_value, psi_r, psi_z)
          br(iz, ir) = alpha*(-psi_z/r(ir))
          bz(iz, ir) = alpha*( psi_r/r(ir))
       ENDDO
    ENDDO

    CALL fit_poloidal_field_scale(geometry, r, z, br, bz, 1.d0, fit)
    CALL assert_true(fit%is_valid, 'signed 1/(2pi) field fit is invalid')
    CALL assert_close(fit%alpha, alpha, 2.d-13, 'signed 1/(2pi) fit')
    CALL assert_close(fit%relative_rms, 0.d0, 2.d-13, 'exact field-fit residual')

    DO ir = 1, nr
       DO iz = 1, nz
          perturbation = 0.01d0*SIN(0.7d0*REAL(ir) + 0.3d0*REAL(iz))
          br(iz, ir) = br(iz, ir)*(1.d0 + perturbation)
          bz(iz, ir) = bz(iz, ir)*(1.d0 - perturbation)
       ENDDO
    ENDDO
    CALL fit_poloidal_field_scale(geometry, r, z, br, bz, 1.d0, fit)
    CALL assert_true(fit%is_valid, 'percent-perturbed field fit is invalid')
    CALL assert_close(fit%alpha, alpha, 2.d-4, 'percent-perturbed fit scale')
    CALL assert_true(fit%relative_rms > 1.d-3 .AND. fit%relative_rms < 0.02d0, &
         'percent-perturbed field residual')
  END SUBROUTINE test_signed_field_fit

  SUBROUTINE test_exact_cache_and_refresh()
    TYPE(equilibrium_geometry_t) :: geometry
    TYPE(magnetic_geometry_cache_t) :: cache
    INTEGER, PARAMETER :: nr = 41, nz = 41
    REAL*8, PARAMETER :: r0 = 2.d0, a = 0.5d0
    REAL*8 :: r(nr), z(nz), psi(nz, nr), wall(4, 2), wall_ref(2)
    REAL*8 :: nodes(4, 2), volume_shape(1, 4), face_shape(1, 2)
    INTEGER :: wall_faces(4, 2), elements(1, 4), face_nodes(4, 2)
    INTEGER :: ierr, generation
    CHARACTER(LEN=256) :: message

    CALL fill_circular_flux(r0, r, z, psi)
    wall(1, :) = (/r0 - a, -0.7d0/)
    wall(2, :) = (/r0 + 0.7d0, -0.7d0/)
    wall(3, :) = (/r0 + 0.7d0,  0.7d0/)
    wall(4, :) = (/r0 - a,  0.7d0/)
    wall_faces = RESHAPE((/1, 2, 2, 3, 3, 4, 4, 1/), SHAPE(wall_faces), ORDER=(/2, 1/))
    wall_ref = (/-1.d0, 1.d0/)
    CALL geometry%init(r, z, psi, ierr, message)
    CALL assert_ok(ierr, message)
    CALL geometry%analyze(1.02d0*a*a, wall, wall_faces, wall_ref, ierr, message)
    CALL assert_ok(ierr, message)

    nodes(1, :) = (/r0 - 0.25d0, -0.25d0/)
    nodes(2, :) = (/r0 + 0.25d0, -0.25d0/)
    nodes(3, :) = (/r0 + 0.25d0,  0.25d0/)
    nodes(4, :) = (/r0 - 0.25d0,  0.25d0/)
    elements(1, :) = (/1, 2, 3, 4/)
    face_nodes(1, :) = (/1, 2/)
    face_nodes(2, :) = (/2, 3/)
    face_nodes(3, :) = (/3, 4/)
    face_nodes(4, :) = (/4, 1/)
    volume_shape(1, :) = 0.25d0
    face_shape(1, :) = 0.5d0

    CALL cache%build(geometry, nodes, elements, volume_shape, face_nodes, face_shape)
    CALL assert_true(cache%is_initialized, 'geometry cache was not initialized')
    CALL assert_close(cache%nodal_rho(1), SQRT(0.5d0), 2.d-5, &
         'exact nodal rho')
    CALL assert_close(cache%volume_rho(1, 1), 0.d0, 2.d-8, &
         'exact volume-quadrature rho')
    CALL assert_close(cache%face_rho(1, 1, 1), 0.5d0, 2.d-5, &
         'exact face-quadrature rho')
    CALL assert_true(ALL(cache%nodal_region == magnetic_region_core), &
         'nodal topology regions')
    generation = cache%generation
    CALL cache%build(geometry, nodes, elements, volume_shape, face_nodes, face_shape)
    CALL assert_true(cache%generation == generation + 1, &
         'cache generation did not advance on refresh')
    CALL assert_true(cache%equilibrium_generation == geometry%generation, &
         'cache/equilibrium generation mismatch')
  END SUBROUTINE test_exact_cache_and_refresh

  SUBROUTINE fill_circular_flux(r0, r, z, psi)
    REAL*8, INTENT(IN) :: r0
    REAL*8, INTENT(OUT) :: r(:), z(:), psi(:, :)
    INTEGER :: ir, iz

    DO ir = 1, SIZE(r)
       r(ir) = r0 - 0.8d0 + 1.6d0*REAL(ir - 1)/REAL(SIZE(r) - 1)
    ENDDO
    DO iz = 1, SIZE(z)
       z(iz) = -0.8d0 + 1.6d0*REAL(iz - 1)/REAL(SIZE(z) - 1)
    ENDDO
    DO ir = 1, SIZE(r)
       DO iz = 1, SIZE(z)
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
       WRITE (*, '(3A,3ES16.7)') 'FAIL: ', TRIM(label), ': ', &
            value, expected, tolerance
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

END PROGRAM test_magnetic_geometry_state
