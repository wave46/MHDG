MODULE magnetic_geometry_state
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE magnetic_topology, ONLY: equilibrium_geometry_t
  IMPLICIT NONE

  PRIVATE

  REAL*8, PARAMETER :: near_null_fraction = 1.d-6
  REAL*8, PARAMETER :: restart_comparison_tolerance = 1.d-6

  TYPE, PUBLIC :: poloidal_field_fit_t
     LOGICAL :: is_valid = .FALSE.
     REAL*8 :: alpha = 0.d0
     REAL*8 :: relative_rms = HUGE(0.d0)
     INTEGER :: sample_count = 0
  END TYPE poloidal_field_fit_t

  TYPE, PUBLIC :: magnetic_geometry_cache_t
     LOGICAL :: is_initialized = .FALSE.
     INTEGER :: generation = 0
     INTEGER :: equilibrium_generation = 0
     REAL*8, ALLOCATABLE :: nodal_psi_normalized(:)
     REAL*8, ALLOCATABLE :: nodal_rho(:)
     REAL*8, ALLOCATABLE :: nodal_normal(:, :)
     INTEGER, ALLOCATABLE :: nodal_region(:)
     REAL*8, ALLOCATABLE :: volume_psi_normalized(:, :)
     REAL*8, ALLOCATABLE :: volume_rho(:, :)
     REAL*8, ALLOCATABLE :: volume_normal(:, :, :)
     INTEGER, ALLOCATABLE :: volume_region(:, :)
     REAL*8, ALLOCATABLE :: face_psi_normalized(:, :, :)
     REAL*8, ALLOCATABLE :: face_rho(:, :, :)
     REAL*8, ALLOCATABLE :: face_normal(:, :, :, :)
     INTEGER, ALLOCATABLE :: face_region(:, :, :)
   CONTAINS
     PROCEDURE :: build => magnetic_geometry_cache_build
     PROCEDURE :: clear => magnetic_geometry_cache_clear
  END TYPE magnetic_geometry_cache_t

  TYPE(equilibrium_geometry_t), SAVE, PUBLIC :: magnetic_equilibrium
  TYPE(magnetic_geometry_cache_t), SAVE, PUBLIC :: magnetic_geometry_cache
  TYPE(poloidal_field_fit_t), SAVE, PUBLIC :: poloidal_field_fit

  LOGICAL, SAVE :: restart_reference_is_present = .FALSE.
  REAL*8, SAVE :: restart_psi_lcfs = 0.d0
  REAL*8, SAVE :: restart_r_axis = 0.d0
  REAL*8, SAVE :: restart_z_axis = 0.d0
  REAL*8, SAVE :: restart_a_minor = 0.d0

  PUBLIC :: fit_poloidal_field_scale
  PUBLIC :: reset_magnetic_geometry_state
  PUBLIC :: set_restart_geometry_reference
  PUBLIC :: compare_restart_geometry

CONTAINS

  SUBROUTINE reset_magnetic_geometry_state()
    CALL magnetic_equilibrium%clear()
    CALL magnetic_geometry_cache%clear()
    poloidal_field_fit = poloidal_field_fit_t()
  END SUBROUTINE reset_magnetic_geometry_state

  SUBROUTINE fit_poloidal_field_scale(geometry, r, z, br, bz, length_scale, fit)
    TYPE(equilibrium_geometry_t), INTENT(IN) :: geometry
    REAL*8, INTENT(IN) :: r(:), z(:)
    REAL*8, INTENT(IN) :: br(:, :), bz(:, :)
    REAL*8, INTENT(IN) :: length_scale
    TYPE(poloidal_field_fit_t), INTENT(OUT) :: fit
    INTEGER :: ir, iz
    REAL*8 :: psi, psi_r, psi_z, gr, gz, gnorm, bnorm
    REAL*8 :: max_gnorm, max_bnorm, numerator, denominator
    REAL*8 :: residual, reference_norm, r_safe, length_scale_squared

    fit = poloidal_field_fit_t()
    IF (.NOT. geometry%is_initialized) RETURN
    IF (SIZE(br, 1) /= SIZE(z) .OR. SIZE(br, 2) /= SIZE(r)) RETURN
    IF (ANY(SHAPE(br) /= SHAPE(bz)) .OR. length_scale <= 0.d0) RETURN

    length_scale_squared = length_scale*length_scale
    max_gnorm = 0.d0
    max_bnorm = 0.d0
    DO ir = 1, SIZE(r)
       DO iz = 1, SIZE(z)
          CALL geometry%evaluate_flux(r(ir), z(iz), psi, psi_r, psi_z)
          r_safe = MAX(ABS(r(ir)), 1.d-12)
          gr = -psi_z/(r_safe*length_scale_squared)
          gz = psi_r/(r_safe*length_scale_squared)
          gnorm = HYPOT(gr, gz)
          bnorm = HYPOT(br(iz, ir), bz(iz, ir))
          IF (ieee_is_finite(gnorm)) max_gnorm = MAX(max_gnorm, gnorm)
          IF (ieee_is_finite(bnorm)) max_bnorm = MAX(max_bnorm, bnorm)
       ENDDO
    ENDDO
    IF (max_gnorm <= TINY(1.d0) .OR. max_bnorm <= TINY(1.d0)) RETURN

    numerator = 0.d0
    denominator = 0.d0
    reference_norm = 0.d0
    DO ir = 1, SIZE(r)
       DO iz = 1, SIZE(z)
          CALL geometry%evaluate_flux(r(ir), z(iz), psi, psi_r, psi_z)
          r_safe = MAX(ABS(r(ir)), 1.d-12)
          gr = -psi_z/(r_safe*length_scale_squared)
          gz = psi_r/(r_safe*length_scale_squared)
          gnorm = HYPOT(gr, gz)
          bnorm = HYPOT(br(iz, ir), bz(iz, ir))
          IF (.NOT. ieee_is_finite(gnorm) .OR. .NOT. ieee_is_finite(bnorm)) CYCLE
          IF (gnorm < near_null_fraction*max_gnorm .OR. &
               bnorm < near_null_fraction*max_bnorm) CYCLE
          numerator = numerator + gr*br(iz, ir) + gz*bz(iz, ir)
          denominator = denominator + gr*gr + gz*gz
          reference_norm = reference_norm + br(iz, ir)**2 + bz(iz, ir)**2
          fit%sample_count = fit%sample_count + 1
       ENDDO
    ENDDO
    IF (fit%sample_count < 4 .OR. denominator <= TINY(1.d0) .OR. &
         reference_norm <= TINY(1.d0)) RETURN

    fit%alpha = numerator/denominator
    residual = MAX(reference_norm - 2.d0*fit%alpha*numerator + &
         fit%alpha*fit%alpha*denominator, 0.d0)
    fit%relative_rms = SQRT(residual/reference_norm)
    fit%is_valid = ieee_is_finite(fit%alpha) .AND. &
         ieee_is_finite(fit%relative_rms)
  END SUBROUTINE fit_poloidal_field_scale

  SUBROUTINE magnetic_geometry_cache_build(this, geometry, node_coordinates, &
       elements, volume_shape, face_nodes, face_shape)
    CLASS(magnetic_geometry_cache_t), INTENT(INOUT) :: this
    TYPE(equilibrium_geometry_t), INTENT(IN) :: geometry
    REAL*8, INTENT(IN) :: node_coordinates(:, :)
    INTEGER, INTENT(IN) :: elements(:, :), face_nodes(:, :)
    REAL*8, INTENT(IN) :: volume_shape(:, :), face_shape(:, :)
    INTEGER :: inode, ielem, iface, igauss
    REAL*8 :: xy(2)
    REAL*8, ALLOCATABLE :: element_coordinates(:, :), face_coordinates(:, :)

    CALL this%clear()
    IF (.NOT. geometry%is_initialized) RETURN
    IF (SIZE(node_coordinates, 2) < 2) RETURN

    ALLOCATE(this%nodal_psi_normalized(SIZE(node_coordinates, 1)))
    ALLOCATE(this%nodal_rho(SIZE(node_coordinates, 1)))
    ALLOCATE(this%nodal_normal(SIZE(node_coordinates, 1), 2))
    ALLOCATE(this%nodal_region(SIZE(node_coordinates, 1)))
    DO inode = 1, SIZE(node_coordinates, 1)
       CALL geometry%evaluate_topology(node_coordinates(inode, 1), &
            node_coordinates(inode, 2), this%nodal_psi_normalized(inode), &
            this%nodal_rho(inode), this%nodal_normal(inode, :), &
            this%nodal_region(inode))
    ENDDO

    ALLOCATE(this%volume_psi_normalized(SIZE(volume_shape, 1), SIZE(elements, 1)))
    ALLOCATE(this%volume_rho(SIZE(volume_shape, 1), SIZE(elements, 1)))
    ALLOCATE(this%volume_normal(SIZE(volume_shape, 1), SIZE(elements, 1), 2))
    ALLOCATE(this%volume_region(SIZE(volume_shape, 1), SIZE(elements, 1)))
    ALLOCATE(element_coordinates(SIZE(elements, 2), 2))
    DO ielem = 1, SIZE(elements, 1)
       element_coordinates = node_coordinates(elements(ielem, :), 1:2)
       DO igauss = 1, SIZE(volume_shape, 1)
          xy = MATMUL(volume_shape(igauss, :), element_coordinates)
          CALL geometry%evaluate_topology(xy(1), xy(2), &
               this%volume_psi_normalized(igauss, ielem), &
               this%volume_rho(igauss, ielem), &
               this%volume_normal(igauss, ielem, :), &
               this%volume_region(igauss, ielem))
       ENDDO
    ENDDO
    DEALLOCATE(element_coordinates)

    ALLOCATE(this%face_psi_normalized(SIZE(face_shape, 1), &
         SIZE(face_nodes, 1), SIZE(elements, 1)))
    ALLOCATE(this%face_rho(SIZE(face_shape, 1), SIZE(face_nodes, 1), &
         SIZE(elements, 1)))
    ALLOCATE(this%face_normal(SIZE(face_shape, 1), SIZE(face_nodes, 1), &
         SIZE(elements, 1), 2))
    ALLOCATE(this%face_region(SIZE(face_shape, 1), SIZE(face_nodes, 1), &
         SIZE(elements, 1)))
    ALLOCATE(face_coordinates(SIZE(face_nodes, 2), 2))
    DO ielem = 1, SIZE(elements, 1)
       DO iface = 1, SIZE(face_nodes, 1)
          face_coordinates = node_coordinates(&
               elements(ielem, face_nodes(iface, :)), 1:2)
          DO igauss = 1, SIZE(face_shape, 1)
             xy = MATMUL(face_shape(igauss, :), face_coordinates)
             CALL geometry%evaluate_topology(xy(1), xy(2), &
                  this%face_psi_normalized(igauss, iface, ielem), &
                  this%face_rho(igauss, iface, ielem), &
                  this%face_normal(igauss, iface, ielem, :), &
                  this%face_region(igauss, iface, ielem))
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(face_coordinates)

    this%generation = this%generation + 1
    this%equilibrium_generation = geometry%generation
    this%is_initialized = .TRUE.
  END SUBROUTINE magnetic_geometry_cache_build

  SUBROUTINE magnetic_geometry_cache_clear(this)
    CLASS(magnetic_geometry_cache_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%nodal_psi_normalized)) DEALLOCATE(this%nodal_psi_normalized)
    IF (ALLOCATED(this%nodal_rho)) DEALLOCATE(this%nodal_rho)
    IF (ALLOCATED(this%nodal_normal)) DEALLOCATE(this%nodal_normal)
    IF (ALLOCATED(this%nodal_region)) DEALLOCATE(this%nodal_region)
    IF (ALLOCATED(this%volume_psi_normalized)) DEALLOCATE(this%volume_psi_normalized)
    IF (ALLOCATED(this%volume_rho)) DEALLOCATE(this%volume_rho)
    IF (ALLOCATED(this%volume_normal)) DEALLOCATE(this%volume_normal)
    IF (ALLOCATED(this%volume_region)) DEALLOCATE(this%volume_region)
    IF (ALLOCATED(this%face_psi_normalized)) DEALLOCATE(this%face_psi_normalized)
    IF (ALLOCATED(this%face_rho)) DEALLOCATE(this%face_rho)
    IF (ALLOCATED(this%face_normal)) DEALLOCATE(this%face_normal)
    IF (ALLOCATED(this%face_region)) DEALLOCATE(this%face_region)
    this%is_initialized = .FALSE.
    this%equilibrium_generation = 0
  END SUBROUTINE magnetic_geometry_cache_clear

  SUBROUTINE set_restart_geometry_reference(psi_lcfs, r_axis, z_axis, a_minor)
    REAL*8, INTENT(IN) :: psi_lcfs, r_axis, z_axis, a_minor

    restart_reference_is_present = .TRUE.
    restart_psi_lcfs = psi_lcfs
    restart_r_axis = r_axis
    restart_z_axis = z_axis
    restart_a_minor = a_minor
  END SUBROUTINE set_restart_geometry_reference

  SUBROUTINE compare_restart_geometry(length_scale, mismatch, message)
    REAL*8, INTENT(IN) :: length_scale
    LOGICAL, INTENT(OUT) :: mismatch
    CHARACTER(LEN=*), INTENT(OUT) :: message
    REAL*8 :: psi_error, axis_error, a_error, current_r, current_z, current_a

    mismatch = .FALSE.
    message = ''
    IF (.NOT. restart_reference_is_present) RETURN
    restart_reference_is_present = .FALSE.

    current_r = magnetic_equilibrium%r_axis*length_scale
    current_z = magnetic_equilibrium%z_axis*length_scale
    current_a = magnetic_equilibrium%a_minor*length_scale
    psi_error = ABS(magnetic_equilibrium%psi_lcfs - restart_psi_lcfs)/ &
         MAX(1.d0, ABS(magnetic_equilibrium%psi_lcfs), ABS(restart_psi_lcfs))
    axis_error = HYPOT(current_r - restart_r_axis, current_z - restart_z_axis)/ &
         MAX(1.d0, ABS(current_r), ABS(current_z), ABS(restart_r_axis), &
         ABS(restart_z_axis))
    a_error = ABS(current_a - restart_a_minor)/ &
         MAX(1.d0, ABS(current_a), ABS(restart_a_minor))
    mismatch = MAX(psi_error, axis_error, a_error) > restart_comparison_tolerance
    IF (mismatch) THEN
       WRITE(message, '(A,3(1X,ES10.3))') &
            'Restart magnetic geometry differs (psi, axis, a):', &
            psi_error, axis_error, a_error
    ENDIF
  END SUBROUTINE compare_restart_geometry

END MODULE magnetic_geometry_state
