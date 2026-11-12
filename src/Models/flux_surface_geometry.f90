! Closed-surface circuit lengths and effective inverse aspect ratios.
MODULE flux_surface_geometry
  USE magnetic_topology, ONLY: equilibrium_geometry_t,lcfs_source_xpoint
  USE interpolation, ONLY: eval_bicubic_value,find_cell_and_local_coordinate,interpolate
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: flux_surface_geometry_t,core_surface_geometry
  REAL*8, PARAMETER :: two_pi=2.d0*ACOS(-1.d0)
  INTEGER, PARAMETER :: uniform_angle_count=720,xpoint_angle_count=123

  TYPE :: flux_surface_geometry_t
    INTEGER :: generation=-1
    LOGICAL :: separatrix_xpoint=.FALSE.
    REAL*8, ALLOCATABLE :: rho(:),L(:),eps(:)
  CONTAINS
    PROCEDURE :: build => build_profiles
    PROCEDURE :: interp => interpolate_profiles
  END TYPE flux_surface_geometry_t

  TYPE :: contour_point_t
    REAL*8 :: major_radius,height,flux_gradient(2),arc_per_radian
  END TYPE contour_point_t

  TYPE(flux_surface_geometry_t), SAVE :: core_surface_geometry

CONTAINS

  SUBROUTINE build_profiles(this,geometry,bt,btr,btz,btrz,alpha,length_scale,br,bz)
    CLASS(flux_surface_geometry_t), INTENT(INOUT) :: this
    TYPE(equilibrium_geometry_t), INTENT(IN) :: geometry
    REAL*8, INTENT(IN) :: bt(:,:),btr(:,:),btz(:,:),btrz(:,:),alpha,length_scale
    REAL*8, OPTIONAL, INTENT(IN) :: br(:,:),bz(:,:)
    REAL*8, ALLOCATABLE :: angles(:),weights(:)

    CALL contour_quadrature(geometry,angles,weights)
    CALL allocate_surface_grid(this)
    CALL integrate_surface_profiles(this,geometry,angles,weights,bt,btr,btz,btrz,alpha,length_scale,br,bz)
    CALL extrapolate_axis(this)
    this%separatrix_xpoint=geometry%lcfs_source==lcfs_source_xpoint
    this%generation=geometry%generation
  END SUBROUTINE build_profiles

  SUBROUTINE interpolate_profiles(this,rho,L,eps)
    CLASS(flux_surface_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(OUT) :: L,eps
    INTEGER :: lower,upper,last
    REAL*8 :: fraction

    IF (.NOT.ALLOCATED(this%rho)) ERROR STOP 'Closed-surface geometry profiles are uninitialized'
    last=SIZE(this%rho)
    IF (rho<=0.d0) THEN
      L=this%L(1)
      eps=0.d0
    ELSEIF (rho>=this%rho(last)) THEN
      L=extrapolate_edge_length(this,rho)
      eps=this%eps(last)
    ELSE
      CALL find_cell_and_local_coordinate(last,this%rho,rho,lower,fraction)
      upper=lower+1
      L=(1.d0-fraction)*this%L(lower)+fraction*this%L(upper)
      eps=(1.d0-fraction)*this%eps(lower)+fraction*this%eps(upper)
    ENDIF
  END SUBROUTINE interpolate_profiles

  SUBROUTINE integrate_surface_profiles(this,geometry,angles,weights,bt,btr,btz,btrz,alpha,length_scale,br,bz)
    CLASS(flux_surface_geometry_t), INTENT(INOUT) :: this
    TYPE(equilibrium_geometry_t), INTENT(IN) :: geometry
    REAL*8, INTENT(IN) :: angles(:),weights(:),bt(:,:),btr(:,:),btz(:,:),btrz(:,:),alpha,length_scale
    REAL*8, OPTIONAL, INTENT(IN) :: br(:,:),bz(:,:)
    REAL*8 :: flux_level
    INTEGER :: surface

    DO surface=2,SIZE(this%rho)
      flux_level=geometry%psi_axis+this%rho(surface)**2*(geometry%psi_lcfs-geometry%psi_axis)
      CALL integrate_closed_contour(geometry,flux_level,angles,weights,bt,btr,btz,btrz, &
           alpha,length_scale,this%L(surface),this%eps(surface),br,bz)
    ENDDO
  END SUBROUTINE integrate_surface_profiles

  SUBROUTINE integrate_closed_contour(geometry,flux_level,angles,weights,bt,btr,btz,btrz, &
       alpha,length_scale,parallel_length,inverse_aspect_ratio,br,bz)
    TYPE(equilibrium_geometry_t), INTENT(IN) :: geometry
    REAL*8, INTENT(IN) :: flux_level,angles(:),weights(:),bt(:,:),btr(:,:),btz(:,:),btrz(:,:)
    REAL*8, INTENT(IN) :: alpha,length_scale
    REAL*8, INTENT(OUT) :: parallel_length,inverse_aspect_ratio
    REAL*8, OPTIONAL, INTENT(IN) :: br(:,:),bz(:,:)
    TYPE(contour_point_t) :: point
    REAL*8 :: arc_length,poloidal_field,toroidal_field,poloidal_length,radius_integral
    INTEGER :: sample

    poloidal_length=0.d0
    radius_integral=0.d0
    parallel_length=0.d0
    DO sample=1,SIZE(angles)
      CALL trace_contour_point(geometry,flux_level,angles(sample),point)
      arc_length=point%arc_per_radian*weights(sample)
      poloidal_field=poloidal_field_strength(geometry,point,alpha,length_scale,br,bz)
      CALL eval_bicubic_value(geometry%nz,geometry%z,geometry%nr,geometry%r, &
           bt,btr,btz,btrz,point%height,point%major_radius,toroidal_field)
      poloidal_length=poloidal_length+arc_length
      radius_integral=radius_integral+point%major_radius*arc_length
      parallel_length=parallel_length+arc_length*SQRT(poloidal_field**2+toroidal_field**2)/poloidal_field*length_scale
    ENDDO
    inverse_aspect_ratio=poloidal_length**2/(two_pi*radius_integral)
  END SUBROUTINE integrate_closed_contour

  SUBROUTINE trace_contour_point(geometry,flux_level,angle,point)
    TYPE(equilibrium_geometry_t), INTENT(IN) :: geometry
    REAL*8, INTENT(IN) :: flux_level,angle
    TYPE(contour_point_t), INTENT(OUT) :: point
    REAL*8 :: ray_radius,cos_angle,sin_angle,flux_value
    REAL*8 :: radial_flux_derivative,angular_flux_derivative,radius_derivative
    LOGICAL :: found

    CALL geometry%closed_surface_radius(flux_level,angle,ray_radius,found)
    IF (.NOT.found) ERROR STOP 'Could not trace closed flux surface'
    cos_angle=COS(angle)
    sin_angle=SIN(angle)
    point%major_radius=geometry%r_axis+ray_radius*cos_angle
    point%height=geometry%z_axis+ray_radius*sin_angle
    CALL geometry%evaluate_flux(point%major_radius,point%height,flux_value, &
         point%flux_gradient(1),point%flux_gradient(2))
    ! Implicit differentiation of psi(ray_radius(theta),theta)=flux_level.
    radial_flux_derivative=DOT_PRODUCT(point%flux_gradient,[cos_angle,sin_angle])
    angular_flux_derivative=ray_radius*DOT_PRODUCT(point%flux_gradient,[-sin_angle,cos_angle])
    IF (ABS(radial_flux_derivative)<1.d-30) ERROR STOP 'Singular radial contour parameter'
    radius_derivative=-angular_flux_derivative/radial_flux_derivative
    point%arc_per_radian=SQRT(ray_radius**2+radius_derivative**2)
  END SUBROUTINE trace_contour_point

  REAL*8 FUNCTION poloidal_field_strength(geometry,point,alpha,length_scale,br,bz) RESULT(field)
    TYPE(equilibrium_geometry_t), INTENT(IN) :: geometry
    TYPE(contour_point_t), INTENT(IN) :: point
    REAL*8, INTENT(IN) :: alpha,length_scale
    REAL*8, OPTIONAL, INTENT(IN) :: br(:,:),bz(:,:)
    REAL*8 :: radial_field,vertical_field

    IF (PRESENT(br).AND.PRESENT(bz)) THEN
      radial_field=interpolate(geometry%nz,geometry%z,geometry%nr,geometry%r, &
           br,point%height,point%major_radius,1.d-12)
      vertical_field=interpolate(geometry%nz,geometry%z,geometry%nr,geometry%r, &
           bz,point%height,point%major_radius,1.d-12)
      field=SQRT(radial_field**2+vertical_field**2)
    ELSE
      field=ABS(alpha)*NORM2(point%flux_gradient)/(point%major_radius*length_scale**2)
    ENDIF
  END FUNCTION poloidal_field_strength

  SUBROUTINE contour_quadrature(geometry,angles,weights)
    TYPE(equilibrium_geometry_t), INTENT(IN) :: geometry
    REAL*8, ALLOCATABLE, INTENT(OUT) :: angles(:),weights(:)
    REAL*8 :: work(uniform_angle_count+xpoint_angle_count)
    INTEGER :: count

    work(:uniform_angle_count)=uniform_contour_angles()
    count=uniform_angle_count
    IF (geometry%lcfs_source==lcfs_source_xpoint) THEN
      work(count+1:)=xpoint_cluster_angles(geometry)
      count=SIZE(work)
    ENDIF
    CALL sort_angles(work(:count))
    CALL remove_repeated_angles(work,count)
    angles=work(:count)
    weights=periodic_trapezoid_weights(angles)
  END SUBROUTINE contour_quadrature

  FUNCTION uniform_contour_angles() RESULT(angles)
    REAL*8 :: angles(uniform_angle_count)
    INTEGER :: sample

    angles=two_pi*[(REAL(sample-1,8),sample=1,uniform_angle_count)]/REAL(uniform_angle_count,8)
  END FUNCTION uniform_contour_angles

  FUNCTION xpoint_cluster_angles(geometry) RESULT(angles)
    TYPE(equilibrium_geometry_t), INTENT(IN) :: geometry
    REAL*8 :: angles(xpoint_angle_count),xpoint_angle,offset
    INTEGER :: sample

    xpoint_angle=ATAN2(geometry%z_xpoint-geometry%z_axis,geometry%r_xpoint-geometry%r_axis)
    angles(1)=MODULO(xpoint_angle,two_pi)
    DO sample=0,60
      offset=10.d0**(-6.d0+5.d0*REAL(sample,8)/60.d0)
      angles(2*sample+2:2*sample+3)=MODULO([xpoint_angle-offset,xpoint_angle+offset],two_pi)
    ENDDO
  END FUNCTION xpoint_cluster_angles

  FUNCTION periodic_trapezoid_weights(angles) RESULT(weights)
    REAL*8, INTENT(IN) :: angles(:)
    REAL*8 :: weights(SIZE(angles)),previous(SIZE(angles)),following(SIZE(angles))

    previous=CSHIFT(angles,-1)
    following=CSHIFT(angles,1)
    previous(1)=previous(1)-two_pi
    following(SIZE(angles))=following(SIZE(angles))+two_pi
    weights=.5d0*(following-previous)
  END FUNCTION periodic_trapezoid_weights

  SUBROUTINE sort_angles(angles)
    REAL*8, INTENT(INOUT) :: angles(:)
    REAL*8 :: current
    INTEGER :: sample,previous

    DO sample=2,SIZE(angles)
      current=angles(sample)
      previous=sample-1
      DO WHILE (previous>=1)
        IF (angles(previous)<=current) EXIT
        angles(previous+1)=angles(previous)
        previous=previous-1
      ENDDO
      angles(previous+1)=current
    ENDDO
  END SUBROUTINE sort_angles

  SUBROUTINE remove_repeated_angles(angles,count)
    REAL*8, INTENT(INOUT) :: angles(:)
    INTEGER, INTENT(INOUT) :: count
    INTEGER :: sample,unique_count

    unique_count=1
    DO sample=2,count
      IF (angles(sample)-angles(unique_count)>1.d-12) THEN
        unique_count=unique_count+1
        angles(unique_count)=angles(sample)
      ENDIF
    ENDDO
    count=unique_count
  END SUBROUTINE remove_repeated_angles

  SUBROUTINE allocate_surface_grid(this)
    CLASS(flux_surface_geometry_t), INTENT(INOUT) :: this
    INTEGER :: surface

    IF (ALLOCATED(this%rho)) DEALLOCATE(this%rho,this%L,this%eps)
    ALLOCATE(this%rho(122),this%L(122),this%eps(122))
    this%rho(1:2)=[0.d0,.001d0]
    this%rho(3:101)=[(.01d0*REAL(surface,8),surface=1,99)]
    this%rho(102:)=[(SQRT(1.d0-10.d0**(-2.d0-.25d0*REAL(surface,8))),surface=0,20)]
  END SUBROUTINE allocate_surface_grid

  SUBROUTINE extrapolate_axis(this)
    CLASS(flux_surface_geometry_t), INTENT(INOUT) :: this
    REAL*8 :: slope_in_flux

    slope_in_flux=(this%L(3)-this%L(2))/(this%rho(3)**2-this%rho(2)**2)
    this%L(1)=this%L(2)-this%rho(2)**2*slope_in_flux
    this%eps(1)=0.d0
    IF (ANY(this%L<=0.d0).OR.ANY(this%eps<0.d0)) ERROR STOP 'Invalid closed-surface geometry profile'
  END SUBROUTINE extrapolate_axis

  REAL*8 FUNCTION extrapolate_edge_length(this,rho) RESULT(length)
    CLASS(flux_surface_geometry_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8 :: last_coordinates(2),coordinate,slope
    INTEGER :: last

    last=SIZE(this%rho)
    IF (this%separatrix_xpoint) THEN
      coordinate=-LOG(MAX(1.d0-MIN(rho,1.d0)**2,1.d-14))
      last_coordinates=-LOG(1.d0-this%rho(last-1:last)**2)
    ELSE
      coordinate=MIN(rho,1.d0)**2
      last_coordinates=this%rho(last-1:last)**2
    ENDIF
    slope=(this%L(last)-this%L(last-1))/(last_coordinates(2)-last_coordinates(1))
    length=this%L(last)+(coordinate-last_coordinates(2))*slope
  END FUNCTION extrapolate_edge_length

END MODULE flux_surface_geometry
