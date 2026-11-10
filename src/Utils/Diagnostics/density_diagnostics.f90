! Optional chord reporting alongside balance diagnostics; feedback is unchanged.
MODULE density_diagnostics_module
  USE globals
  USE MPI_OMP
  USE interpolation, ONLY: point_in_mesh, nodesearch, integrate_line_samples
  USE magnetic_geometry_state, ONLY: magnetic_equilibrium
  USE magnetic_topology, ONLY: magnetic_region_core
  USE HDF5_io_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: measure_density_chord, report_density_diagnostics, save_density_diagnostics
  PUBLIC :: reset_density_diagnostics, density_diagnostics_available
  INTEGER, PARAMETER :: density_samples = 1000, no_owner = HUGE(1)
  REAL*8, PARAMETER :: intersection_tolerance = 1.d-12
  TYPE :: density_snapshot_type
     REAL*8 :: chord(4) = 0.d0 ! physical metres: R1,Z1,R2,Z2
     REAL*8 :: measurement_time = 0.d0
     REAL*8 :: full_integral = 0.d0
     REAL*8 :: mesh_length = 0.d0
     REAL*8 :: lcfs_length = 0.d0
  END TYPE density_snapshot_type
  TYPE(density_snapshot_type), SAVE :: snapshot
  LOGICAL, SAVE :: density_diagnostics_available = .FALSE.

CONTAINS

  SUBROUTINE reset_density_diagnostics()
    density_diagnostics_available = .FALSE.
    snapshot = density_snapshot_type()
  END SUBROUTINE reset_density_diagnostics

  SUBROUTINE density_error(message)
    CHARACTER(LEN=*), INTENT(IN) :: message
    IF (MPIvar%glob_id == 0) WRITE(6,'(2A)') 'Density chord: ', message
    STOP 1
  END SUBROUTINE density_error

  ! Geometric lengths; full integral uses the feedback's sampling convention.
  SUBROUTINE measure_density_chord(chord, full, mesh_length, lcfs_length)
    REAL*8, INTENT(IN) :: chord(4) ! physical metres: R1,Z1,R2,Z2
    REAL*8, INTENT(OUT) :: full, mesh_length, lcfs_length
    REAL*8, ALLOCATABLE :: points(:,:), segment_lengths(:), values(:)
    INTEGER, ALLOCATABLE :: owners(:), elements(:), connectivity(:,:)
    REAL*8 :: chord_length
    INTEGER :: i

    elements = [(i,i=1,Mesh%Nelems)]
#ifdef PARALL
    elements = PACK(elements,Mesh%ghostElems == 0)
#endif
    connectivity = Mesh%T(elements,:)
    CALL build_chord_points(chord,connectivity,points,segment_lengths,chord_length)
    CALL sample_owned_density(points,elements,connectivity,owners,values)
    CALL measure_chord_lengths(points(:,density_samples+1:),segment_lengths, &
         owners(density_samples+1:),chord_length,mesh_length,lcfs_length)
    full = integrate_line_samples(points(1,:density_samples),points(2,:density_samples), &
         values,owners(:density_samples) /= no_owner)*simpar%refval_density
  END SUBROUTINE measure_density_chord

  ! Density quadrature points followed by midpoints of geometric intervals.
  SUBROUTINE build_chord_points(chord,connectivity,points,segment_lengths,chord_length)
    REAL*8, INTENT(IN) :: chord(4)
    INTEGER, INTENT(IN) :: connectivity(:,:)
    REAL*8, ALLOCATABLE, INTENT(OUT) :: points(:,:),segment_lengths(:)
    REAL*8, INTENT(OUT) :: chord_length
    REAL*8, ALLOCATABLE :: cuts(:)
    REAL*8 :: origin(2),direction(2)
    INTEGER :: i,j,segments

    IF (.NOT. magnetic_equilibrium%is_initialized) CALL density_error('LCFS geometry is unavailable')
    IF (.NOT. ALLOCATED(magnetic_equilibrium%lcfs_r)) CALL density_error('LCFS polygon is unavailable')
    origin = chord(1:2)/phys%lscale
    direction = (chord(3:4)-chord(1:2))/phys%lscale
    chord_length = NORM2(direction)*phys%lscale
    IF (chord_length <= 0.d0) CALL density_error('endpoints must be distinct')
    CALL collect_chord_cuts(origin,direction,connectivity,cuts)
    segments = COUNT(cuts(2:)-cuts(:SIZE(cuts)-1) > intersection_tolerance)
    ALLOCATE(points(2,density_samples+segments),segment_lengths(segments))
    DO i = 1,density_samples
       points(:,i) = (chord(1:2)+(chord(3:4)-chord(1:2))/density_samples*(i-1))/phys%lscale
    ENDDO
    j = 0
    DO i = 1,SIZE(cuts)-1
       IF (cuts(i+1)-cuts(i) <= intersection_tolerance) CYCLE
       j = j+1
       points(:,density_samples+j) = origin+0.5d0*(cuts(i)+cuts(i+1))*direction
       segment_lengths(j) = (cuts(i+1)-cuts(i))*chord_length
    ENDDO
  END SUBROUTINE build_chord_points

  ! Sorted chord parameters at mesh/LCFS crossings, including both endpoints.
  SUBROUTINE collect_chord_cuts(origin,direction,connectivity,cuts)
    REAL*8, INTENT(IN) :: origin(2),direction(2)
    INTEGER, INTENT(IN) :: connectivity(:,:)
    REAL*8, ALLOCATABLE, INTENT(OUT) :: cuts(:)
    REAL*8 :: swap
    INTEGER :: nv,nc,i,j
#ifdef PARALL
    INTEGER :: ierr
    INTEGER, ALLOCATABLE :: counts(:),offsets(:)
    REAL*8, ALLOCATABLE :: all_cuts(:)
#endif
    nv = SIZE(magnetic_equilibrium%lcfs_r)
    ALLOCATE(cuts(nv+3*SIZE(connectivity,1)+2))
    cuts(1:2) = [0.d0,1.d0]
    nc = 2
    DO i = 1,nv
       j = MOD(i,nv)+1
       CALL append_intersection([magnetic_equilibrium%lcfs_r(i),magnetic_equilibrium%lcfs_z(i)], &
                                [magnetic_equilibrium%lcfs_r(j),magnetic_equilibrium%lcfs_z(j)])
    ENDDO
    DO i = 1,SIZE(connectivity,1)
       DO j = 1,3
          CALL append_intersection(Mesh%X(connectivity(i,j),:),Mesh%X(connectivity(i,MOD(j,3)+1),:))
       ENDDO
    ENDDO
#ifdef PARALL
    ALLOCATE(counts(MPIvar%glob_size),offsets(MPIvar%glob_size))
    CALL MPI_ALLGATHER(nc,1,MPI_INTEGER,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    offsets(1) = 0
    DO i = 2,MPIvar%glob_size
       offsets(i) = offsets(i-1)+counts(i-1)
    ENDDO
    ALLOCATE(all_cuts(SUM(counts)))
    CALL MPI_ALLGATHERV(cuts,nc,MPI_REAL8,all_cuts,counts,offsets,MPI_REAL8,MPI_COMM_WORLD,ierr)
    CALL MOVE_ALLOC(all_cuts,cuts)
    nc = SIZE(cuts)
#else
    cuts = cuts(:nc)
#endif
    DO i = 2,nc
       swap = cuts(i)
       j = i-1
       DO WHILE (j >= 1)
          IF (cuts(j) <= swap) EXIT
          cuts(j+1) = cuts(j)
          j = j-1
       ENDDO
       cuts(j+1) = swap
    ENDDO

  CONTAINS
    SUBROUTINE append_intersection(first,last)
      REAL*8, INTENT(IN) :: first(2),last(2)
      REAL*8 :: edge(2),offset(2),denominator,a,b
      edge = last-first
      offset = first-origin
      denominator = direction(1)*edge(2)-direction(2)*edge(1)
      IF (ABS(denominator) <= intersection_tolerance*NORM2(direction)*NORM2(edge)) RETURN
      a = (offset(1)*edge(2)-offset(2)*edge(1))/denominator
      b = (offset(1)*direction(2)-offset(2)*direction(1))/denominator
      IF (a <= 0.d0 .OR. a >= 1.d0 .OR. b < -intersection_tolerance .OR. b > 1.d0+intersection_tolerance) RETURN
      nc = nc+1
      cuts(nc) = a
    END SUBROUTINE append_intersection
  END SUBROUTINE collect_chord_cuts

  ! Owned elements only; shared-edge samples get one rank owner.
  SUBROUTINE sample_owned_density(points,elements,connectivity,owners,values)
    REAL*8, INTENT(IN) :: points(:,:)
    INTEGER, INTENT(IN) :: elements(:),connectivity(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: owners(:)
    REAL*8, ALLOCATABLE, INTENT(OUT) :: values(:)
    INTEGER :: i,local_element,element,node
    LOGICAL :: inside
#ifdef PARALL
    INTEGER :: ierr
#endif
    ALLOCATE(owners(SIZE(points,2)),values(density_samples))
    owners = no_owner
    values = 0.d0
    DO i = 1,SIZE(points,2)
       CALL point_in_mesh(points(1,i),points(2,i),Mesh%X,connectivity,inside,local_element)
       IF (.NOT. inside) CYCLE
       owners(i) = MPIvar%glob_id
       IF (i > density_samples) CYCLE
       element = elements(local_element)
       node = nodesearch(points(1,i),points(2,i),refElPol%Nnodes2D, &
                        Mesh%X(Mesh%T(element,:),1),Mesh%X(Mesh%T(element,:),2))
       values(i) = sol%u(((element-1)*refElPol%Nnodes2D+node-1)*phys%Neq+1)
    ENDDO
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,owners,SIZE(owners),MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
#endif
    WHERE (owners(:density_samples) /= MPIvar%glob_id) values = 0.d0
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,values,density_samples,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
  END SUBROUTINE sample_owned_density

  SUBROUTINE measure_chord_lengths(points,segment_lengths,owners,chord_length,mesh_length,lcfs_length)
    REAL*8, INTENT(IN) :: points(:,:),segment_lengths(:),chord_length
    INTEGER, INTENT(IN) :: owners(:)
    REAL*8, INTENT(OUT) :: mesh_length,lcfs_length
    INTEGER :: i

    mesh_length = 0.d0
    lcfs_length = 0.d0
    DO i = 1,SIZE(segment_lengths)
       IF (owners(i) /= no_owner) mesh_length = mesh_length+segment_lengths(i)
       IF (magnetic_equilibrium%classify(points(1,i),points(2,i),2.d0) /= magnetic_region_core) CYCLE
       IF (owners(i) == no_owner) CALL density_error('mesh does not cover the inside-LCFS chord')
       lcfs_length = lcfs_length+segment_lengths(i)
    ENDDO
    IF (lcfs_length <= intersection_tolerance*chord_length) &
         CALL density_error('chord does not cross the LCFS interior')
  END SUBROUTINE measure_chord_lengths

  ! Called beside the balance report, before the Newton update. This snapshot
  ! is intentionally not recomputed from the later HDF5 solution fields.
  SUBROUTINE report_density_diagnostics()
    IF (.NOT. utils%line_integrated_density_diagnostics) RETURN
    snapshot%chord = input%density_chord
    CALL measure_density_chord(snapshot%chord,snapshot%full_integral,snapshot%mesh_length,snapshot%lcfs_length)
    snapshot%measurement_time = time%t*simpar%refval_time
    density_diagnostics_available = .TRUE.
    IF (MPIvar%glob_id /= 0) RETURN
    WRITE(6,'(A,ES10.3,A,ES10.3,A,ES10.3,A)') &
         'Line-integrated density: raw = ',snapshot%full_integral,' [m^-2]; average (in-mesh chord) = ', &
         snapshot%full_integral/snapshot%mesh_length,' [m^-3]; average (inside-LCFS chord) = ', &
         snapshot%full_integral/snapshot%lcfs_length,' [m^-3]'
  END SUBROUTINE report_density_diagnostics

  SUBROUTINE save_density_diagnostics(parent_id)
    INTEGER(HID_T), INTENT(IN) :: parent_id
    INTEGER(HID_T) :: group_id
    IF (.NOT. utils%line_integrated_density_diagnostics .OR. .NOT. density_diagnostics_available) RETURN
    CALL HDF5_group_create('density_chord',parent_id,group_id)
    CALL HDF5_real_saving(group_id,snapshot%measurement_time,'measurement_time_seconds')
    CALL HDF5_real_saving(group_id,snapshot%full_integral,'full_integral_m-2')
    CALL HDF5_real_saving(group_id,snapshot%mesh_length,'mesh_length_m')
    CALL HDF5_real_saving(group_id,snapshot%lcfs_length,'lcfs_length_m')
    CALL HDF5_array1D_saving(group_id,snapshot%chord,4,'endpoints_m')
    CALL HDF5_group_close(group_id)
  END SUBROUTINE save_density_diagnostics
END MODULE density_diagnostics_module
