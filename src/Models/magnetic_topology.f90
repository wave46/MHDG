MODULE magnetic_topology
  USE interpolation, ONLY: build_bicubic_derivatives, &
       eval_bicubic_with_derivatives, eval_bicubic_with_2nd_derivatives
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: topology_unknown = 0
  INTEGER, PARAMETER, PUBLIC :: topology_limited = 1
  INTEGER, PARAMETER, PUBLIC :: topology_lower_single_null = 2

  INTEGER, PARAMETER, PUBLIC :: lcfs_source_unknown = 0
  INTEGER, PARAMETER, PUBLIC :: lcfs_source_wall = 1
  INTEGER, PARAMETER, PUBLIC :: lcfs_source_xpoint = 2

  INTEGER, PARAMETER, PUBLIC :: magnetic_region_undefined = 0
  INTEGER, PARAMETER, PUBLIC :: magnetic_region_core = 1
  INTEGER, PARAMETER, PUBLIC :: magnetic_region_main_sol = 2
  INTEGER, PARAMETER, PUBLIC :: magnetic_region_private_flux = 3

  INTEGER, PARAMETER :: max_critical_points = 64
  INTEGER, PARAMETER :: lcfs_polygon_points = 720
  INTEGER, PARAMETER :: ray_scan_points = 512
  INTEGER, PARAMETER :: wall_samples_per_node = 6
  INTEGER, PARAMETER :: newton_max_iterations = 30
  REAL*8, PARAMETER :: geometry_tol = 1.d-12

  TYPE :: critical_point_t
     REAL*8 :: r = 0.d0
     REAL*8 :: z = 0.d0
     REAL*8 :: psi = 0.d0
     REAL*8 :: hessian_det = 0.d0
     LOGICAL :: elliptic = .FALSE.
     LOGICAL :: saddle = .FALSE.
  END TYPE critical_point_t

  TYPE, PUBLIC :: equilibrium_geometry_t
     LOGICAL :: is_initialized = .FALSE.
     INTEGER :: nr = 0
     INTEGER :: nz = 0
     INTEGER :: generation = 0
     INTEGER :: topology_kind = topology_unknown
     INTEGER :: lcfs_source = lcfs_source_unknown
     REAL*8 :: psi_sep_input = 0.d0
     REAL*8 :: psi_axis = 0.d0
     REAL*8 :: psi_lcfs = 0.d0
     REAL*8 :: r_axis = 0.d0
     REAL*8 :: z_axis = 0.d0
     REAL*8 :: r_xpoint = 0.d0
     REAL*8 :: z_xpoint = 0.d0
     REAL*8 :: r_wall_contact = 0.d0
     REAL*8 :: z_wall_contact = 0.d0
     REAL*8 :: r_lcfs_min = 0.d0
     REAL*8 :: r_lcfs_max = 0.d0
     REAL*8 :: a_minor = 0.d0
     REAL*8, ALLOCATABLE :: r(:)
     REAL*8, ALLOCATABLE :: z(:)
     REAL*8, ALLOCATABLE :: psi(:, :)
     REAL*8, ALLOCATABLE :: psi_r(:, :)
     REAL*8, ALLOCATABLE :: psi_z(:, :)
     REAL*8, ALLOCATABLE :: psi_rz(:, :)
     REAL*8, ALLOCATABLE :: lcfs_r(:)
     REAL*8, ALLOCATABLE :: lcfs_z(:)
   CONTAINS
     PROCEDURE :: init => equilibrium_geometry_init
     PROCEDURE :: analyze => equilibrium_geometry_analyze
     PROCEDURE :: evaluate_flux => equilibrium_geometry_evaluate_flux
     PROCEDURE :: evaluate_topology => equilibrium_geometry_evaluate_topology
     PROCEDURE :: classify => equilibrium_geometry_classify
     PROCEDURE :: clear => equilibrium_geometry_clear
     FINAL :: equilibrium_geometry_finalize
  END TYPE equilibrium_geometry_t

CONTAINS

  SUBROUTINE equilibrium_geometry_init(this, r, z, psi, ierr, message)
    CLASS(equilibrium_geometry_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: r(:), z(:), psi(:, :)
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    CALL this%clear()
    ierr = 0
    message = ''

    IF (SIZE(r) < 4 .OR. SIZE(z) < 4) THEN
       ierr = 1
       message = 'Equilibrium grid needs at least four points in each direction'
       RETURN
    ENDIF
    IF (SIZE(psi, 1) /= SIZE(z) .OR. SIZE(psi, 2) /= SIZE(r)) THEN
       ierr = 2
       message = 'Equilibrium psi dimensions do not match the coordinate vectors'
       RETURN
    ENDIF
    IF (ANY(r(2:) <= r(:SIZE(r) - 1)) .OR. ANY(z(2:) <= z(:SIZE(z) - 1))) THEN
       ierr = 3
       message = 'Equilibrium coordinate vectors must be strictly increasing'
       RETURN
    ENDIF

    this%nr = SIZE(r)
    this%nz = SIZE(z)
    ALLOCATE(this%r(this%nr), this%z(this%nz))
    ALLOCATE(this%psi(this%nz, this%nr))
    ALLOCATE(this%psi_r(this%nz, this%nr))
    ALLOCATE(this%psi_z(this%nz, this%nr))
    ALLOCATE(this%psi_rz(this%nz, this%nr))
    this%r = r
    this%z = z
    this%psi = psi
    CALL build_bicubic_derivatives(this%nz, this%z, this%nr, this%r, &
         this%psi, this%psi_r, this%psi_z, this%psi_rz)
    this%generation = this%generation + 1
    this%is_initialized = .TRUE.
  END SUBROUTINE equilibrium_geometry_init

  SUBROUTINE equilibrium_geometry_clear(this)
    CLASS(equilibrium_geometry_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%r)) DEALLOCATE(this%r)
    IF (ALLOCATED(this%z)) DEALLOCATE(this%z)
    IF (ALLOCATED(this%psi)) DEALLOCATE(this%psi)
    IF (ALLOCATED(this%psi_r)) DEALLOCATE(this%psi_r)
    IF (ALLOCATED(this%psi_z)) DEALLOCATE(this%psi_z)
    IF (ALLOCATED(this%psi_rz)) DEALLOCATE(this%psi_rz)
    IF (ALLOCATED(this%lcfs_r)) DEALLOCATE(this%lcfs_r)
    IF (ALLOCATED(this%lcfs_z)) DEALLOCATE(this%lcfs_z)
    this%is_initialized = .FALSE.
    this%nr = 0
    this%nz = 0
    this%topology_kind = topology_unknown
    this%lcfs_source = lcfs_source_unknown
    this%psi_sep_input = 0.d0
    this%psi_axis = 0.d0
    this%psi_lcfs = 0.d0
    this%r_axis = 0.d0
    this%z_axis = 0.d0
    this%r_xpoint = 0.d0
    this%z_xpoint = 0.d0
    this%r_wall_contact = 0.d0
    this%z_wall_contact = 0.d0
    this%r_lcfs_min = 0.d0
    this%r_lcfs_max = 0.d0
    this%a_minor = 0.d0
  END SUBROUTINE equilibrium_geometry_clear

  SUBROUTINE equilibrium_geometry_finalize(this)
    TYPE(equilibrium_geometry_t), INTENT(INOUT) :: this
    CALL this%clear()
  END SUBROUTINE equilibrium_geometry_finalize

  SUBROUTINE equilibrium_geometry_analyze(this, psi_sep_input, wall_coordinates, &
       wall_faces, wall_reference_nodes, ierr, message)
    CLASS(equilibrium_geometry_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: psi_sep_input
    REAL*8, INTENT(IN) :: wall_coordinates(:, :)
    INTEGER, INTENT(IN) :: wall_faces(:, :)
    REAL*8, INTENT(IN) :: wall_reference_nodes(:)
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    TYPE(critical_point_t) :: points(max_critical_points)
    INTEGER :: npoints, iaxis, ixpoint
    REAL*8 :: psi_wall, span_seed, outward_sign
    REAL*8 :: wall_event_delta, x_event_delta, event_tolerance
    LOGICAL :: wall_found

    ierr = 0
    message = ''
    IF (.NOT. this%is_initialized) THEN
       ierr = 10
       message = 'Equilibrium geometry has not been initialized'
       RETURN
    ENDIF
    IF (SIZE(wall_coordinates, 2) /= 2 .OR. &
         SIZE(wall_faces, 2) /= SIZE(wall_reference_nodes)) THEN
       ierr = 11
       message = 'Wall coordinates/connectivity/reference nodes are inconsistent'
       RETURN
    ENDIF
    IF (SIZE(wall_faces, 1) <= 0) THEN
       ierr = 12
       message = 'A non-periodic physical wall is required for LCFS analysis'
       RETURN
    ENDIF

    this%psi_sep_input = psi_sep_input
    CALL find_critical_points(this, wall_coordinates, wall_faces, &
         wall_reference_nodes, points, npoints, ierr, message)
    IF (ierr /= 0) RETURN

    CALL select_axis(points, npoints, iaxis, ierr, message)
    IF (ierr /= 0) RETURN
    this%r_axis = points(iaxis)%r
    this%z_axis = points(iaxis)%z
    this%psi_axis = points(iaxis)%psi

    span_seed = psi_sep_input - this%psi_axis
    IF (ABS(span_seed) <= geometry_tol*MAX(1.d0, ABS(this%psi_axis))) THEN
       ierr = 13
       message = 'Input psiSep is indistinguishable from the magnetic-axis flux'
       RETURN
    ENDIF

    outward_sign = SIGN(1.d0, span_seed)
    CALL select_lower_xpoint(this, points, npoints, span_seed, ixpoint)
    CALL find_wall_contact(this, wall_coordinates, wall_faces, wall_reference_nodes, &
         outward_sign, psi_wall, this%r_wall_contact, &
         this%z_wall_contact, wall_found)

    wall_event_delta = HUGE(0.d0)
    IF (wall_found) wall_event_delta = outward_sign*(psi_wall - this%psi_axis)
    x_event_delta = HUGE(0.d0)
    IF (ixpoint > 0) x_event_delta = outward_sign*(points(ixpoint)%psi - this%psi_axis)
    event_tolerance = 1.d-6*ABS(span_seed)

    ! Candidate events have already been checked against the axis-connected
    ! component.  Their signed flux distance can therefore order which event
    ! is reached first without confusing equal-flux private-flux components.
    IF (ixpoint > 0 .AND. x_event_delta <= wall_event_delta + event_tolerance) THEN
       this%topology_kind = topology_lower_single_null
       this%lcfs_source = lcfs_source_xpoint
       this%r_xpoint = points(ixpoint)%r
       this%z_xpoint = points(ixpoint)%z
       this%psi_lcfs = points(ixpoint)%psi
    ELSEIF (wall_found) THEN
       this%topology_kind = topology_limited
       this%lcfs_source = lcfs_source_wall
       this%psi_lcfs = psi_wall
    ELSEIF (ixpoint > 0) THEN
       this%topology_kind = topology_lower_single_null
       this%lcfs_source = lcfs_source_xpoint
       this%r_xpoint = points(ixpoint)%r
       this%z_xpoint = points(ixpoint)%z
       this%psi_lcfs = points(ixpoint)%psi
    ELSE
       ierr = 14
       message = 'Could not identify either a wall-limited LCFS or a lower X-point'
       RETURN
    ENDIF

    IF (ABS(this%psi_lcfs - this%psi_axis) <= &
         geometry_tol*MAX(1.d0, ABS(this%psi_axis))) THEN
       ierr = 15
       message = 'Detected LCFS has zero flux span from the magnetic axis'
       RETURN
    ENDIF

    CALL build_axis_connected_lcfs(this, ierr, message)
    IF (ierr /= 0) RETURN
    CALL refine_lcfs_radial_extrema(this)
    this%a_minor = 0.5d0*(this%r_lcfs_max - this%r_lcfs_min)
    IF (this%a_minor <= geometry_tol) THEN
       ierr = 16
       message = 'Detected LCFS has a non-positive minor radius'
    ENDIF
  END SUBROUTINE equilibrium_geometry_analyze

  SUBROUTINE equilibrium_geometry_evaluate_flux(this, r, z, psi, dpsi_dr, &
       dpsi_dz, d2psi_dr2, d2psi_dz2, d2psi_drdz)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: r, z
    REAL*8, INTENT(OUT) :: psi, dpsi_dr, dpsi_dz
    REAL*8, INTENT(OUT), OPTIONAL :: d2psi_dr2, d2psi_dz2, d2psi_drdz
    REAL*8 :: drr, dzz, drz

    IF (PRESENT(d2psi_dr2) .OR. PRESENT(d2psi_dz2) .OR. PRESENT(d2psi_drdz)) THEN
       CALL eval_bicubic_with_2nd_derivatives(this%nz, this%z, this%nr, this%r, &
            this%psi, this%psi_r, this%psi_z, this%psi_rz, z, r, psi, &
            dpsi_dz, dpsi_dr, dzz, drr, drz)
       IF (PRESENT(d2psi_dr2)) d2psi_dr2 = drr
       IF (PRESENT(d2psi_dz2)) d2psi_dz2 = dzz
       IF (PRESENT(d2psi_drdz)) d2psi_drdz = drz
    ELSE
       CALL eval_bicubic_with_derivatives(this%nz, this%z, this%nr, this%r, &
            this%psi, this%psi_r, this%psi_z, this%psi_rz, z, r, psi, &
            dpsi_dz, dpsi_dr)
    ENDIF
  END SUBROUTINE equilibrium_geometry_evaluate_flux

  SUBROUTINE equilibrium_geometry_evaluate_topology(this, r, z, psi_normalized, rho, &
       outward_normal, region)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: r, z
    REAL*8, INTENT(OUT) :: psi_normalized, rho, outward_normal(2)
    INTEGER, INTENT(OUT) :: region
    REAL*8 :: psi, dpsi_dr, dpsi_dz, span, norm_grad

    CALL this%evaluate_flux(r, z, psi, dpsi_dr, dpsi_dz)
    span = this%psi_lcfs - this%psi_axis
    psi_normalized = (psi - this%psi_axis)/span
    rho = SQRT(MAX(psi_normalized, 0.d0))
    outward_normal = SIGN(1.d0, span)*(/dpsi_dr, dpsi_dz/)
    norm_grad = NORM2(outward_normal)
    IF (norm_grad > geometry_tol) THEN
       outward_normal = outward_normal/norm_grad
    ELSE
       outward_normal = 0.d0
    ENDIF
    region = this%classify(r, z, psi_normalized)
  END SUBROUTINE equilibrium_geometry_evaluate_topology

  INTEGER FUNCTION equilibrium_geometry_classify(this, r, z, psi_normalized)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: r, z, psi_normalized

    equilibrium_geometry_classify = magnetic_region_undefined
    IF (.NOT. this%is_initialized .OR. .NOT. ALLOCATED(this%lcfs_r)) RETURN
    IF (point_in_polygon(r, z, this%lcfs_r, this%lcfs_z)) THEN
       equilibrium_geometry_classify = magnetic_region_core
    ELSEIF (psi_normalized < 1.d0) THEN
       equilibrium_geometry_classify = magnetic_region_private_flux
    ELSE
       equilibrium_geometry_classify = magnetic_region_main_sol
    ENDIF
  END FUNCTION equilibrium_geometry_classify

  SUBROUTINE find_critical_points(this, wall_coordinates, wall_faces, &
       wall_reference_nodes, points, npoints, ierr, message)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: wall_coordinates(:, :), wall_reference_nodes(:)
    INTEGER, INTENT(IN) :: wall_faces(:, :)
    TYPE(critical_point_t), INTENT(OUT) :: points(max_critical_points)
    INTEGER, INTENT(OUT) :: npoints, ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: ir, iz
    REAL*8 :: gr(4), gz(4), r0, z0
    TYPE(critical_point_t) :: point
    LOGICAL :: converged

    npoints = 0
    ierr = 0
    message = ''
    DO iz = 1, this%nz - 1
       DO ir = 1, this%nr - 1
          gr = (/this%psi_r(iz, ir), this%psi_r(iz, ir + 1), &
               this%psi_r(iz + 1, ir), this%psi_r(iz + 1, ir + 1)/)
          gz = (/this%psi_z(iz, ir), this%psi_z(iz, ir + 1), &
               this%psi_z(iz + 1, ir), this%psi_z(iz + 1, ir + 1)/)
          IF (MINVAL(gr) > 0.d0 .OR. MAXVAL(gr) < 0.d0) CYCLE
          IF (MINVAL(gz) > 0.d0 .OR. MAXVAL(gz) < 0.d0) CYCLE
          r0 = 0.5d0*(this%r(ir) + this%r(ir + 1))
          z0 = 0.5d0*(this%z(iz) + this%z(iz + 1))
          CALL refine_critical_point(this, r0, z0, point, converged)
          IF (.NOT. converged) CYCLE
          IF (.NOT. point_inside_wall(point%r, point%z, wall_coordinates, &
               wall_faces, wall_reference_nodes)) CYCLE
          IF (critical_point_is_duplicate(this, points, npoints, point)) CYCLE
          IF (npoints >= max_critical_points) THEN
             ierr = 20
             message = 'Too many magnetic critical points for supported topology'
             RETURN
          ENDIF
          npoints = npoints + 1
          points(npoints) = point
       ENDDO
    ENDDO
    IF (npoints == 0) THEN
       ierr = 21
       message = 'No magnetic critical point was found inside the physical wall'
    ENDIF
  END SUBROUTINE find_critical_points

  SUBROUTINE refine_critical_point(this, r, z, point, converged)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: r, z
    TYPE(critical_point_t), INTENT(OUT) :: point
    LOGICAL, INTENT(OUT) :: converged
    INTEGER :: iteration
    REAL*8 :: rr, zz, psi, pr, pz, prr, pzz, prz, det, dr, dz
    REAL*8 :: grad_scale, domain_scale

    rr = r
    zz = z
    converged = .FALSE.
    domain_scale = MAX(this%r(this%nr) - this%r(1), this%z(this%nz) - this%z(1))
    DO iteration = 1, newton_max_iterations
       CALL this%evaluate_flux(rr, zz, psi, pr, pz, prr, pzz, prz)
       det = prr*pzz - prz*prz
       IF (ABS(det) <= geometry_tol) RETURN
       dr = (-pr*pzz + pz*prz)/det
       dz = (-pz*prr + pr*prz)/det
       rr = MIN(MAX(rr + dr, this%r(1)), this%r(this%nr))
       zz = MIN(MAX(zz + dz, this%z(1)), this%z(this%nz))
       grad_scale = MAX(1.d0, ABS(psi)/MAX(domain_scale, geometry_tol))
       IF (SQRT(pr*pr + pz*pz) <= 1.d-9*grad_scale .AND. &
            SQRT(dr*dr + dz*dz) <= 1.d-9*domain_scale) THEN
          converged = .TRUE.
          EXIT
       ENDIF
    ENDDO
    IF (.NOT. converged) RETURN

    CALL this%evaluate_flux(rr, zz, psi, pr, pz, prr, pzz, prz)
    point%r = rr
    point%z = zz
    point%psi = psi
    point%hessian_det = prr*pzz - prz*prz
    point%elliptic = point%hessian_det > geometry_tol
    point%saddle = point%hessian_det < -geometry_tol
  END SUBROUTINE refine_critical_point

  LOGICAL FUNCTION critical_point_is_duplicate(this, points, npoints, point)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    TYPE(critical_point_t), INTENT(IN) :: points(max_critical_points), point
    INTEGER, INTENT(IN) :: npoints
    INTEGER :: i
    REAL*8 :: tolerance

    critical_point_is_duplicate = .FALSE.
    tolerance = 1.d-6*MAX(this%r(this%nr) - this%r(1), &
         this%z(this%nz) - this%z(1))
    DO i = 1, npoints
       IF (HYPOT(points(i)%r - point%r, points(i)%z - point%z) <= tolerance) THEN
          critical_point_is_duplicate = .TRUE.
          RETURN
       ENDIF
    ENDDO
  END FUNCTION critical_point_is_duplicate

  SUBROUTINE select_axis(points, npoints, iaxis, ierr, message)
    TYPE(critical_point_t), INTENT(IN) :: points(max_critical_points)
    INTEGER, INTENT(IN) :: npoints
    INTEGER, INTENT(OUT) :: iaxis, ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: i, count_axis

    iaxis = 0
    count_axis = 0
    ierr = 0
    message = ''
    DO i = 1, npoints
       IF (.NOT. points(i)%elliptic) CYCLE
       count_axis = count_axis + 1
       iaxis = i
    ENDDO
    IF (count_axis == 0) THEN
       ierr = 22
       message = 'No elliptic magnetic axis was found inside the physical wall'
    ELSEIF (count_axis > 1) THEN
       ierr = 23
       message = 'Multiple magnetic axes are outside the supported PR03 topology'
    ENDIF
  END SUBROUTINE select_axis

  SUBROUTINE select_lower_xpoint(this, points, npoints, span_seed, ixpoint)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    TYPE(critical_point_t), INTENT(IN) :: points(max_critical_points)
    INTEGER, INTENT(IN) :: npoints
    REAL*8, INTENT(IN) :: span_seed
    INTEGER, INTENT(OUT) :: ixpoint
    INTEGER :: i
    REAL*8 :: candidate_score, best_score, theta, radius, distance, domain_scale
    LOGICAL :: connected

    ixpoint = 0
    best_score = HUGE(0.d0)
    domain_scale = MAX(this%r(this%nr) - this%r(1), &
         this%z(this%nz) - this%z(1))
    DO i = 1, npoints
       IF (.NOT. points(i)%saddle) CYCLE
       IF (points(i)%z >= this%z_axis) CYCLE
       IF (SIGN(1.d0, span_seed)*(points(i)%psi - this%psi_axis) <= geometry_tol) CYCLE
       theta = ATAN2(points(i)%z - this%z_axis, points(i)%r - this%r_axis)
       distance = HYPOT(points(i)%r - this%r_axis, points(i)%z - this%z_axis)
       CALL first_level_crossing(this, points(i)%psi, theta, radius, connected)
       IF (.NOT. connected) CYCLE
       IF (ABS(radius - distance) > 1.d-5*domain_scale) CYCLE
       candidate_score = ABS((points(i)%psi - this%psi_sep_input)/span_seed)
       IF (candidate_score < best_score) THEN
          ixpoint = i
          best_score = candidate_score
       ENDIF
    ENDDO
  END SUBROUTINE select_lower_xpoint

  SUBROUTINE find_wall_contact(this, wall_coordinates, wall_faces, wall_reference_nodes, &
       outward_sign, psi_contact, r_contact, z_contact, found)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: wall_coordinates(:, :), wall_reference_nodes(:)
    INTEGER, INTENT(IN) :: wall_faces(:, :)
    REAL*8, INTENT(IN) :: outward_sign
    REAL*8, INTENT(OUT) :: psi_contact, r_contact, z_contact
    LOGICAL, INTENT(OUT) :: found
    INTEGER :: iface, isample, nsample
    REAL*8 :: s0, s1, s2, f0, f1, f2, smin, r, z, psi
    REAL*8 :: best_delta, delta

    found = .FALSE.
    best_delta = HUGE(0.d0)
    psi_contact = 0.d0
    r_contact = 0.d0
    z_contact = 0.d0
    nsample = MAX(12, wall_samples_per_node*SIZE(wall_reference_nodes))

    DO iface = 1, SIZE(wall_faces, 1)
       s0 = -1.d0
       CALL wall_flux(this, wall_coordinates, wall_faces(iface, :), &
            wall_reference_nodes, s0, outward_sign, f0, r, z, psi)
       CALL consider_wall_candidate(f0, r, z, psi)
       s1 = -1.d0 + 2.d0/REAL(nsample)
       CALL wall_flux(this, wall_coordinates, wall_faces(iface, :), &
            wall_reference_nodes, s1, outward_sign, f1, r, z, psi)
       CALL consider_wall_candidate(f1, r, z, psi)
       DO isample = 2, nsample
          s2 = -1.d0 + 2.d0*REAL(isample)/REAL(nsample)
          CALL wall_flux(this, wall_coordinates, wall_faces(iface, :), &
               wall_reference_nodes, s2, outward_sign, f2, r, z, psi)
          CALL consider_wall_candidate(f2, r, z, psi)
          IF (f1 <= f0 .AND. f1 <= f2) THEN
             CALL minimize_wall_interval(this, wall_coordinates, wall_faces(iface, :), &
                  wall_reference_nodes, outward_sign, s0, s2, smin)
             CALL wall_flux(this, wall_coordinates, wall_faces(iface, :), &
                  wall_reference_nodes, smin, outward_sign, delta, r, z, psi)
             CALL consider_wall_candidate(delta, r, z, psi)
          ENDIF
          s0 = s1
          f0 = f1
          s1 = s2
          f1 = f2
       ENDDO
    ENDDO

  CONTAINS
    SUBROUTINE consider_wall_candidate(candidate_delta, candidate_r, candidate_z, candidate_psi)
      REAL*8, INTENT(IN) :: candidate_delta, candidate_r, candidate_z, candidate_psi
      REAL*8 :: theta, radius, distance, domain_scale
      LOGICAL :: connected
      IF (candidate_delta <= geometry_tol) RETURN
      theta = ATAN2(candidate_z - this%z_axis, candidate_r - this%r_axis)
      distance = HYPOT(candidate_r - this%r_axis, candidate_z - this%z_axis)
      CALL first_level_crossing(this, candidate_psi, theta, radius, connected)
      IF (.NOT. connected) RETURN
      domain_scale = MAX(this%r(this%nr) - this%r(1), &
           this%z(this%nz) - this%z(1))
      IF (ABS(radius - distance) > 1.d-5*domain_scale) RETURN
      IF (candidate_delta < best_delta) THEN
         best_delta = candidate_delta
         psi_contact = candidate_psi
         r_contact = candidate_r
         z_contact = candidate_z
         found = .TRUE.
      ENDIF
    END SUBROUTINE consider_wall_candidate
  END SUBROUTINE find_wall_contact

  SUBROUTINE minimize_wall_interval(this, wall_coordinates, wall_face, &
       wall_reference_nodes, outward_sign, sa, sb, smin)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: wall_coordinates(:, :), wall_reference_nodes(:)
    INTEGER, INTENT(IN) :: wall_face(:)
    REAL*8, INTENT(IN) :: outward_sign, sa, sb
    REAL*8, INTENT(OUT) :: smin
    INTEGER :: iteration
    REAL*8, PARAMETER :: golden = 0.6180339887498948482d0
    REAL*8 :: left, right, s1, s2, f1, f2, r, z, psi

    left = sa
    right = sb
    s1 = right - golden*(right - left)
    s2 = left + golden*(right - left)
    CALL wall_flux(this, wall_coordinates, wall_face, wall_reference_nodes, &
         s1, outward_sign, f1, r, z, psi)
    CALL wall_flux(this, wall_coordinates, wall_face, wall_reference_nodes, &
         s2, outward_sign, f2, r, z, psi)
    DO iteration = 1, 50
       IF (f1 <= f2) THEN
          right = s2
          s2 = s1
          f2 = f1
          s1 = right - golden*(right - left)
          CALL wall_flux(this, wall_coordinates, wall_face, wall_reference_nodes, &
               s1, outward_sign, f1, r, z, psi)
       ELSE
          left = s1
          s1 = s2
          f1 = f2
          s2 = left + golden*(right - left)
          CALL wall_flux(this, wall_coordinates, wall_face, wall_reference_nodes, &
               s2, outward_sign, f2, r, z, psi)
       ENDIF
       IF (ABS(right - left) <= 1.d-12) EXIT
    ENDDO
    IF (f1 <= f2) THEN
       smin = s1
    ELSE
       smin = s2
    ENDIF
  END SUBROUTINE minimize_wall_interval

  SUBROUTINE wall_flux(this, wall_coordinates, wall_face, wall_reference_nodes, &
       s, outward_sign, delta, r, z, psi)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: wall_coordinates(:, :), wall_reference_nodes(:)
    INTEGER, INTENT(IN) :: wall_face(:)
    REAL*8, INTENT(IN) :: s, outward_sign
    REAL*8, INTENT(OUT) :: delta, r, z, psi
    REAL*8 :: pr, pz

    CALL evaluate_wall_position(wall_coordinates, wall_face, wall_reference_nodes, &
         s, r, z)
    CALL this%evaluate_flux(r, z, psi, pr, pz)
    delta = outward_sign*(psi - this%psi_axis)
  END SUBROUTINE wall_flux

  SUBROUTINE evaluate_wall_position(wall_coordinates, wall_face, reference_nodes, &
       s, r, z)
    REAL*8, INTENT(IN) :: wall_coordinates(:, :), reference_nodes(:), s
    INTEGER, INTENT(IN) :: wall_face(:)
    REAL*8, INTENT(OUT) :: r, z
    INTEGER :: i, j, n
    REAL*8 :: basis

    n = SIZE(reference_nodes)
    r = 0.d0
    z = 0.d0
    DO i = 1, n
       basis = 1.d0
       DO j = 1, n
          IF (j == i) CYCLE
          basis = basis*(s - reference_nodes(j))/(reference_nodes(i) - reference_nodes(j))
       ENDDO
       r = r + basis*wall_coordinates(wall_face(i), 1)
       z = z + basis*wall_coordinates(wall_face(i), 2)
    ENDDO
  END SUBROUTINE evaluate_wall_position

  LOGICAL FUNCTION point_inside_wall(r, z, wall_coordinates, wall_faces, reference_nodes)
    REAL*8, INTENT(IN) :: r, z, wall_coordinates(:, :), reference_nodes(:)
    INTEGER, INTENT(IN) :: wall_faces(:, :)
    INTEGER :: iface, isample, nsample
    REAL*8 :: s, r0, z0, r1, z1, intersection_r

    point_inside_wall = .FALSE.
    nsample = MAX(8, 3*SIZE(reference_nodes))
    DO iface = 1, SIZE(wall_faces, 1)
       CALL evaluate_wall_position(wall_coordinates, wall_faces(iface, :), &
            reference_nodes, -1.d0, r0, z0)
       DO isample = 1, nsample
          s = -1.d0 + 2.d0*REAL(isample)/REAL(nsample)
          CALL evaluate_wall_position(wall_coordinates, wall_faces(iface, :), &
               reference_nodes, s, r1, z1)
          IF ((z0 > z) .NEQV. (z1 > z)) THEN
             intersection_r = r0 + (z - z0)*(r1 - r0)/(z1 - z0)
             IF (intersection_r > r) point_inside_wall = .NOT. point_inside_wall
          ENDIF
          r0 = r1
          z0 = z1
       ENDDO
    ENDDO
  END FUNCTION point_inside_wall

  SUBROUTINE build_axis_connected_lcfs(this, ierr, message)
    CLASS(equilibrium_geometry_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: i
    REAL*8 :: theta, radius
    LOGICAL :: found

    ierr = 0
    message = ''
    IF (ALLOCATED(this%lcfs_r)) DEALLOCATE(this%lcfs_r)
    IF (ALLOCATED(this%lcfs_z)) DEALLOCATE(this%lcfs_z)
    ALLOCATE(this%lcfs_r(lcfs_polygon_points), this%lcfs_z(lcfs_polygon_points))

    DO i = 1, lcfs_polygon_points
       theta = 2.d0*ACOS(-1.d0)*REAL(i - 1)/REAL(lcfs_polygon_points)
       CALL first_lcfs_crossing(this, theta, radius, found)
       IF (.NOT. found) THEN
          ierr = 30
          WRITE(message, '(A,I0,A)') 'Axis-connected LCFS was not found on radial ray ', i, &
               '; unsupported non-star-shaped or inconsistent topology'
          RETURN
       ENDIF
       this%lcfs_r(i) = this%r_axis + radius*COS(theta)
       this%lcfs_z(i) = this%z_axis + radius*SIN(theta)
    ENDDO
  END SUBROUTINE build_axis_connected_lcfs

  SUBROUTINE first_lcfs_crossing(this, theta, radius, found)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: theta
    REAL*8, INTENT(OUT) :: radius
    LOGICAL, INTENT(OUT) :: found
    CALL first_level_crossing(this, this%psi_lcfs, theta, radius, found)
  END SUBROUTINE first_lcfs_crossing

  SUBROUTINE first_level_crossing(this, psi_level, theta, radius, found)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: psi_level, theta
    REAL*8, INTENT(OUT) :: radius
    LOGICAL, INTENT(OUT) :: found
    INTEGER :: i, iteration, imax
    REAL*8 :: c, s, rmax, t0, t1, f0, f1, tm, fm, psi, pr, pz
    REAL*8 :: prr, pzz, prz, dfdt, d2fdt2, tleft, tright, dt
    REAL*8 :: signed_span, max_f

    c = COS(theta)
    s = SIN(theta)
    rmax = ray_limit(this, c, s)
    signed_span = SIGN(1.d0, psi_level - this%psi_axis)
    t0 = 0.d0
    f0 = -ABS(psi_level - this%psi_axis)
    max_f = f0
    imax = 0
    found = .FALSE.
    radius = 0.d0

    DO i = 1, ray_scan_points
       t1 = rmax*REAL(i)/REAL(ray_scan_points)
       CALL this%evaluate_flux(this%r_axis + t1*c, this%z_axis + t1*s, psi, pr, pz)
       f1 = signed_span*(psi - psi_level)
       IF (f1 > max_f) THEN
          max_f = f1
          imax = i
       ENDIF
       IF (f1 >= 0.d0) THEN
          DO iteration = 1, 60
             tm = 0.5d0*(t0 + t1)
             CALL this%evaluate_flux(this%r_axis + tm*c, this%z_axis + tm*s, psi, pr, pz)
             fm = signed_span*(psi - psi_level)
             IF (fm >= 0.d0) THEN
                t1 = tm
             ELSE
                t0 = tm
             ENDIF
          ENDDO
          radius = 0.5d0*(t0 + t1)
          found = .TRUE.
          RETURN
       ENDIF
       t0 = t1
       f0 = f1
    ENDDO

    IF (imax <= 0) RETURN

    ! At an X point the separatrix can merely touch a radial ray, so a sign
    ! change is not guaranteed.  Refine the sampled maximum of psi along the
    ! ray before deciding whether this is that tangential LCFS contact.
    dt = rmax/REAL(ray_scan_points)
    tleft = MAX(0.d0, REAL(imax - 1)*dt)
    tright = MIN(rmax, REAL(imax + 1)*dt)
    tm = REAL(imax)*dt
    DO iteration = 1, newton_max_iterations
       CALL this%evaluate_flux(this%r_axis + tm*c, this%z_axis + tm*s, psi, &
            pr, pz, prr, pzz, prz)
       dfdt = signed_span*(pr*c + pz*s)
       d2fdt2 = signed_span*(prr*c*c + 2.d0*prz*c*s + pzz*s*s)
       IF (ABS(d2fdt2) <= geometry_tol) EXIT
       t1 = MIN(MAX(tm - dfdt/d2fdt2, tleft), tright)
       IF (ABS(t1 - tm) <= 1.d-13*MAX(1.d0, rmax)) THEN
          tm = t1
          EXIT
       ENDIF
       tm = t1
    ENDDO
    CALL this%evaluate_flux(this%r_axis + tm*c, this%z_axis + tm*s, psi, pr, pz)
    fm = signed_span*(psi - psi_level)
    IF (ABS(fm) <= 1.d-7*ABS(psi_level - this%psi_axis)) THEN
       radius = tm
       found = .TRUE.
    ENDIF
  END SUBROUTINE first_level_crossing

  REAL*8 FUNCTION ray_limit(this, c, s)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: c, s
    REAL*8 :: tr, tz

    tr = HUGE(0.d0)
    tz = HUGE(0.d0)
    IF (c > geometry_tol) tr = (this%r(this%nr) - this%r_axis)/c
    IF (c < -geometry_tol) tr = (this%r(1) - this%r_axis)/c
    IF (s > geometry_tol) tz = (this%z(this%nz) - this%z_axis)/s
    IF (s < -geometry_tol) tz = (this%z(1) - this%z_axis)/s
    ray_limit = 0.999999d0*MIN(tr, tz)
  END FUNCTION ray_limit

  SUBROUTINE refine_lcfs_radial_extrema(this)
    CLASS(equilibrium_geometry_t), INTENT(INOUT) :: this
    INTEGER :: imin(1), imax(1)
    REAL*8 :: rmin, zmin, rmax, zmax
    LOGICAL :: okmin, okmax

    imin = MINLOC(this%lcfs_r)
    imax = MAXLOC(this%lcfs_r)
    rmin = this%lcfs_r(imin(1))
    zmin = this%lcfs_z(imin(1))
    rmax = this%lcfs_r(imax(1))
    zmax = this%lcfs_z(imax(1))
    CALL refine_radial_extremum(this, rmin, zmin, okmin)
    CALL refine_radial_extremum(this, rmax, zmax, okmax)
    this%r_lcfs_min = rmin
    this%r_lcfs_max = rmax
  END SUBROUTINE refine_lcfs_radial_extrema

  SUBROUTINE refine_radial_extremum(this, r, z, converged)
    CLASS(equilibrium_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(INOUT) :: r, z
    LOGICAL, INTENT(OUT) :: converged
    INTEGER :: iteration
    REAL*8 :: psi, pr, pz, prr, pzz, prz, det, dr, dz

    converged = .FALSE.
    DO iteration = 1, newton_max_iterations
       CALL this%evaluate_flux(r, z, psi, pr, pz, prr, pzz, prz)
       det = pr*pzz - pz*prz
       IF (ABS(det) <= geometry_tol) RETURN
       dr = (-(psi - this%psi_lcfs)*pzz + pz*pz)/det
       dz = (-pr*pz + prz*(psi - this%psi_lcfs))/det
       r = MIN(MAX(r + dr, this%r(1)), this%r(this%nr))
       z = MIN(MAX(z + dz, this%z(1)), this%z(this%nz))
       IF (HYPOT(dr, dz) <= 1.d-11*MAX(1.d0, ABS(r))) THEN
          converged = .TRUE.
          RETURN
       ENDIF
    ENDDO
  END SUBROUTINE refine_radial_extremum

  LOGICAL FUNCTION point_in_polygon(r, z, polygon_r, polygon_z)
    REAL*8, INTENT(IN) :: r, z, polygon_r(:), polygon_z(:)
    INTEGER :: i, j, n
    REAL*8 :: crossing_r

    point_in_polygon = .FALSE.
    n = SIZE(polygon_r)
    j = n
    DO i = 1, n
       IF ((polygon_z(i) > z) .NEQV. (polygon_z(j) > z)) THEN
          crossing_r = polygon_r(i) + (z - polygon_z(i))* &
               (polygon_r(j) - polygon_r(i))/(polygon_z(j) - polygon_z(i))
          IF (crossing_r > r) point_in_polygon = .NOT. point_in_polygon
       ENDIF
       j = i
    ENDDO
  END FUNCTION point_in_polygon

END MODULE magnetic_topology
