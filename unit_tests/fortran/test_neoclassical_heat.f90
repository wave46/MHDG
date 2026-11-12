PROGRAM test_neoclassical_heat
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys, switch, simpar
  USE physics, ONLY: parallel_heat_conduction_flux,heat_conductivity,conductive_heat_flux_result
  USE flux_surface_geometry, ONLY: flux_surface_geometry_t,core_surface_geometry
  USE magnetic_topology, ONLY: equilibrium_geometry_t
  USE interpolation, ONLY: build_bicubic_derivatives
  IMPLICIT NONE
  INTEGER :: row,i,j,k,ierr,checks
  REAL*8 :: temperature,density,epsilon,factor,slope,plus,minus,unused,step
  REAL*8 :: expected,length,aspect,radius,pi,worst
  REAL*8 :: temperatures(4),densities(3),aspects(3)
  INTEGER, PARAMETER :: n=81
  REAL*8 :: r(n),z(n),psi(n,n),bt(n,n),btr(n,n),btz(n,n),btrz(n,n),br(n,n),bz(n,n)
  TYPE(equilibrium_geometry_t) :: equilibrium
  TYPE(flux_surface_geometry_t) :: profiles,tail
  CHARACTER(256) :: message

  phys%Mref=1.d0;phys%epn=2.5d0
  phys%neoclassical_heat_conduction=.TRUE.
  switch%testcase=2 ! Check the correction independently of the temperature cap.
  simpar%refval_temperature=1.d0;simpar%refval_density=1.d19
  ALLOCATE(core_surface_geometry%rho(2),core_surface_geometry%L(2),core_surface_geometry%eps(2))
  core_surface_geometry%rho=[0.d0,1.d0]
  checks=0;worst=0.d0;pi=ACOS(-1.d0);step=1.d-5
  temperatures=[1.d0,50.d0,500.d0,1.d4]
  densities=[1.d18,1.d19,1.d21];aspects=[0.d0,.01d0,.3d0]
  DO row=3,4
    DO i=1,SIZE(temperatures)
      DO j=1,SIZE(densities)
        DO k=1,SIZE(aspects)
          temperature=temperatures(i);density=densities(j);epsilon=aspects(k)
          CALL conductivity_factor(temperature,density,25.d0,epsilon,row,factor,slope)
          IF (.NOT.ieee_is_finite(factor).OR.factor<=0.d0) ERROR STOP 'invalid heat factor or axis limit'
          CALL conductivity_factor(temperature,density*(1.d0+step),25.d0,epsilon,row,plus,unused)
          CALL conductivity_factor(temperature,density*(1.d0-step),25.d0,epsilon,row,minus,unused)
          CALL close((plus-minus)/(2.d0*step)/factor,slope,2.d-8,'density logarithmic derivative')
          CALL conductivity_factor(temperature*(1.d0+step),density,25.d0,epsilon,row,plus,unused)
          CALL conductivity_factor(temperature*(1.d0-step),density,25.d0,epsilon,row,minus,unused)
          CALL close((plus-minus)/(2.d0*step)/factor,-2.d0*slope,2.d-8,'temperature logarithmic derivative')
        ENDDO
      ENDDO
    ENDDO
    CALL conductivity_factor(.01d0,1.d22,25.d0,.2d0,row,factor,slope)
    expected=.58d0*.89d0/.52d0
    IF (row==3) expected=1.60d0*.74d0/1.17d0
    CALL close(factor,expected,1.d-7,'collisional asymptote')
  ENDDO

  CALL check_object_reuse()

  ! Circular psi, Btor=C/R: |B|/|Bpol| is constant on each circle.
  ! This independently checks full-circuit length, 2*pi normalization,
  ! effective inverse aspect ratio, axis extrapolation and physical length units.
  DO i=1,n
    r(i)=1.75d0+1.5d0*REAL(i-1)/REAL(n-1)
    z(i)=-.75d0+1.5d0*REAL(i-1)/REAL(n-1)
  ENDDO
  DO i=1,n
    psi(:,i)=(r(i)-2.5d0)**2+z**2
    bt(:,i)=2.d0/r(i)
  ENDDO
  CALL equilibrium%init(r,z,psi,ierr,message)
  IF (ierr/=0) ERROR STOP 'circular equilibrium initialization'
  equilibrium%r_axis=2.5d0;equilibrium%z_axis=0.d0
  equilibrium%psi_axis=0.d0;equilibrium%psi_lcfs=.25d0
  CALL build_bicubic_derivatives(n,z,n,r,bt,btr,btz,btrz)
  DO j=1,2
    CALL profiles%build(equilibrium,bt,btr,btz,btrz,1.d0,REAL(j,8))
    DO i=1,SIZE(profiles%rho),7
      radius=.5d0*profiles%rho(i)
      CALL profiles%interp(profiles%rho(i),length,aspect)
      expected=2.d0*pi*REAL(j,8)*SQRT(radius**2+REAL(j,8)**4)
      CALL close(length/expected,1.d0,3.d-7,'circular parallel length')
      CALL close(aspect,radius/2.5d0,3.d-8,'circular inverse aspect ratio')
    ENDDO
  ENDDO
  CALL profiles%interp(1.d0,length,aspect)
  expected=4.d0*pi*SQRT(.5d0**2+16.d0)
  CALL close(length/expected,1.d0,3.d-7,'finite limiter-edge length')
  DO i=1,n
    br(:,i)=-4.d0*z/r(i)
    bz(:,i)=4.d0*(r(i)-2.5d0)/r(i)
  ENDDO
  CALL profiles%build(equilibrium,bt,btr,btz,btrz,1.d0,1.d0,br,bz)
  CALL profiles%interp(.5d0,length,aspect)
  expected=2.d0*pi*SQRT(.25d0**2+.25d0)
  CALL close(length/expected,1.d0,3.d-4,'supplied-field circuit length')
  ALLOCATE(tail%rho(3),tail%L(3),tail%eps(3))
  tail%rho=[0.d0,.9d0,.99d0]
  tail%L=2.d0;tail%L(2:)=2.d0-3.d0*LOG(1.d0-tail%rho(2:)**2)
  tail%eps=.2d0;tail%separatrix_xpoint=.TRUE.
  CALL tail%interp(.999d0,length,aspect)
  CALL close(length,2.d0-3.d0*LOG(1.d0-.999d0**2),2.d-13,'X-point logarithmic continuation')
  IF (profiles%generation/=equilibrium%generation) ERROR STOP 'geometry provenance'
  WRITE(*,*) 'Neoclassical heat: PASS; checks, worst scaled discrepancy:',checks,worst
CONTAINS
  SUBROUTINE conductivity_factor(temperature_eV,density_m3,length_m,epsilon,row,factor,slope)
    REAL*8, INTENT(IN) :: temperature_eV,density_m3,length_m,epsilon
    INTEGER, INTENT(IN) :: row
    REAL*8, INTENT(OUT) :: factor,slope
    REAL*8 :: state(4),energy,alpha
    TYPE(parallel_heat_conduction_flux) :: heat_flux
    TYPE(heat_conductivity) :: conductivity

    core_surface_geometry%L=length_m;core_surface_geometry%eps=epsilon
    energy=1.5d0*phys%Mref*temperature_eV/simpar%refval_temperature
    state=[density_m3/simpar%refval_density,0.d0,0.d0,0.d0]
    state(row)=state(1)*energy
    CALL heat_flux%init(state,1,.5d0)
    CALL heat_flux%evaluate_conductivities()
    conductivity=heat_flux%ion_heat_flux%conductivity
    IF (row==4) conductivity=heat_flux%electron_heat_flux%conductivity
    alpha=conductivity%value;slope=conductivity%density_degree
    factor=alpha/energy**phys%epn
  END SUBROUTINE conductivity_factor

  SUBROUTINE check_object_reuse()
    TYPE(parallel_heat_conduction_flux) :: reused,fresh
    REAL*8 :: state(4),gradient(3,4),field(3),saved_state(4),saved_gradient(3,4),saved_field(3)
    INTEGER :: region

    state=[1.d0,.1d0,2.d0,3.d0];gradient=.2d0;field=[.2d0,.3d0,.4d0]
    phys%T_fluxlim_maxi=10.d0;phys%T_fluxlim_maxe=10.d0
    phys%diff_pari=.4d0;phys%diff_pare=.7d0;phys%c_fli=.6d0;phys%c_fle=.3d0
    switch%flux_limiter=.TRUE.
    core_surface_geometry%L=25.d0;core_surface_geometry%eps=.2d0
    CALL reused%init(state,1,.5d0,gradient,field)
    saved_state=state;saved_gradient=gradient;saved_field=field
    state=state*1.1d0;gradient=-2.d0*gradient;field=-field
    CALL reused%evaluate_fluxes()
    CALL fresh%init(saved_state,1,.5d0,saved_gradient,saved_field)
    CALL fresh%evaluate_fluxes()
    CALL compare_result(reused%ion_heat_flux,fresh%ion_heat_flux)
    CALL compare_result(reused%electron_heat_flux,fresh%electron_heat_flux)
    CALL reused%evaluate_fluxes()
    CALL compare_result(reused%ion_heat_flux,fresh%ion_heat_flux)
    CALL compare_result(reused%electron_heat_flux,fresh%electron_heat_flux)
    CALL reused%evaluate_conductivities()
    IF (reused%ion_heat_flux%diffusive_flux/=0.d0 .OR. reused%electron_heat_flux%limiter/=1.d0) &
         ERROR STOP 'conductivity evaluation retained previous flux results'
    CALL reused%evaluate_fluxes()
    CALL compare_result(reused%ion_heat_flux,fresh%ion_heat_flux)
    CALL compare_result(reused%electron_heat_flux,fresh%electron_heat_flux)
    state=saved_state;gradient=saved_gradient;field=saved_field
    DO region=2,1,-1
      state(3:4)=state(3:4)*1.2d0
      CALL reused%init(state,region,.5d0,gradient(:2,:),field)
      CALL reused%evaluate_fluxes()
      CALL fresh%init(state,region,.5d0,gradient(:2,:),field)
      CALL fresh%evaluate_fluxes()
      CALL compare_result(reused%ion_heat_flux,fresh%ion_heat_flux)
      CALL compare_result(reused%electron_heat_flux,fresh%electron_heat_flux)
    ENDDO
    phys%neoclassical_heat_conduction=.FALSE.
    CALL reused%init(state,2,.5d0)
    CALL reused%evaluate_conductivities()
    CALL fresh%init(state,2,.5d0)
    CALL fresh%evaluate_conductivities()
    CALL compare_result(reused%ion_heat_flux,fresh%ion_heat_flux)
    CALL compare_result(reused%electron_heat_flux,fresh%electron_heat_flux)
    IF (reused%ion_heat_flux%diffusive_flux/=0.d0 .OR. reused%electron_heat_flux%limiter/=1.d0) &
         ERROR STOP 'conductivity-only evaluation retained flux data'
    phys%neoclassical_heat_conduction=.TRUE.
  END SUBROUTINE check_object_reuse

  SUBROUTINE compare_result(actual,reference)
    TYPE(conductive_heat_flux_result), INTENT(IN) :: actual,reference
    CALL close(actual%conductivity%value,reference%conductivity%value,1.d-14,'reused conductivity')
    CALL close(MAXVAL(ABS(actual%conductivity%state_derivative-reference%conductivity%state_derivative)), &
         0.d0,1.d-14,'reused conductivity derivative')
    CALL close(actual%diffusive_flux,reference%diffusive_flux,1.d-14,'reused flux')
    CALL close(MAXVAL(ABS(actual%state_derivative-reference%state_derivative)),0.d0,1.d-14,'reused state derivative')
    CALL close(MAXVAL(ABS(actual%gradient_derivative-reference%gradient_derivative)),0.d0,1.d-14,'reused gradient derivative')
    CALL close(actual%affine,reference%affine,1.d-14,'reused affine')
    CALL close(actual%limiter,reference%limiter,1.d-14,'reused limiter')
  END SUBROUTINE compare_result

  SUBROUTINE close(actual,reference,tolerance,label)
    REAL*8, INTENT(IN) :: actual,reference,tolerance
    CHARACTER(*), INTENT(IN) :: label
    REAL*8 :: error
    error=ABS(actual-reference)/MAX(1.d0,ABS(reference))
    IF (.NOT.ieee_is_finite(error).OR.error>tolerance) THEN
      WRITE(*,*) label,actual,reference,error
      ERROR STOP 'neoclassical heat check'
    ENDIF
    checks=checks+1;worst=MAX(worst,error)
  END SUBROUTINE close
END PROGRAM test_neoclassical_heat
