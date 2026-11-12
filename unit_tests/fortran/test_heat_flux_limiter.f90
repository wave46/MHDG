PROGRAM test_heat_flux_limiter
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys,simpar,switch
  USE physics, ONLY: parallel_heat_conduction_flux,conductive_heat_flux_result,computeVi,computeVe
  USE adimensionalization_module, ONLY: ion_mass_si,electron_mass_si
  IMPLICIT NONE
  REAL*8, PARAMETER :: perturbation=2.d-7,heat_flux_floor=1.d-10
  INTEGER :: species,checks,reconstructions
  REAL*8 :: worst

  simpar%refval_temperature=1.d0
  phys%Mref=1.3d0;phys%epn=2.5d0
  phys%diff_pari=.4d0;phys%diff_pare=.7d0
  phys%c_fli=.2d0;phys%c_fle=.3d0
  phys%T_fluxlim_maxi=100.d0;phys%T_fluxlim_maxe=100.d0
  switch%testcase=50
  checks=0;reconstructions=0;worst=0.d0
  DO species=1,2
    CALL check_fixture(species,'ordinary',1.d0,2.d0,.TRUE.)
    CALL check_fixture(species,'temperature floor',1.d0,.5d-5,.TRUE.)
    CALL check_fixture(species,'temperature and flux floors',1.d-5,.5d-5,.TRUE.)
    CALL check_fixture(species,'flux floor',1.d-14,2.d0,.TRUE.)
    CALL check_fixture(species,'limiter disabled',1.d0,2.d0,.FALSE.)
  ENDDO
  WRITE(*,*) 'Heat limiter: PASS; held-phi derivatives, reconstructions, worst:',checks,reconstructions,worst

CONTAINS

  SUBROUTINE check_fixture(which,label,density,temperature,enabled)
    INTEGER, INTENT(IN) :: which
    CHARACTER(*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: density,temperature
    LOGICAL, INTENT(IN) :: enabled
    REAL*8 :: u(5),q(2,5),b(3),v(5),ju(5),jq(2,5)
    REAL*8 :: alpha,coefficient,qfs,qscale,cap,x,phi(2),expected,affine,realized
    REAL*8 :: trial(5),qtrial(2,5),scale,plus,minus
    INTEGER :: direction,j,k,energy
    TYPE(parallel_heat_conduction_flux) :: heat_flux
    TYPE(conductive_heat_flux_result) :: response

    switch%flux_limiter=enabled
    energy=which+2
    b=[.3d0,-.4d0,SQRT(.75d0)]
    u=[density,.2d0*density,density*(1.5d0*phys%Mref*temperature+.02d0), &
         density*1.5d0*phys%Mref*temperature,.1d0*density]
    CALL coefficients(u,which,v,alpha,coefficient,qfs)
    cap=MAX(ABS(qfs),heat_flux_floor)
    qscale=.2d0*cap/(coefficient*alpha*v(energy)*b(1))
    DO direction=-1,1
      q=0.d0;q(1,energy)=direction*.7d0*qscale
      x=coefficient*alpha*DOT_PRODUCT(MATMUL(q,v),b(1:2))
      CALL heat_flux%init(u,Q=q,b=b)
      CALL heat_flux%evaluate_fluxes()
      response=heat_flux%ion_heat_flux
      IF (which==2) response=heat_flux%electron_heat_flux
      phi(which)=response%limiter
      expected=1.d0
      IF (enabled) expected=1.d0/(1.d0+ABS(x)/(MERGE(phys%c_fli,phys%c_fle,which==1)*cap))
      CALL close(phi(which),expected,1.d0,label//' limiter value')
      ju=0.d0;ju(:4)=response%state_derivative
      jq=0.d0;jq(:,:4)=response%gradient_derivative(:2,:)
      affine=response%affine
      realized=DOT_PRODUCT(ju,u)+SUM(jq*q)-affine
      CALL close(realized,phi(which)*x,cap,label//' physical-flux reconstruction')
      CALL close(DOT_PRODUCT(ju,u)+SUM(jq*q),0.d0,cap,label//' joint degree-zero contraction')
      reconstructions=reconstructions+1
      DO j=1,5
        scale=MAX(ABS(u(j)),density*1.d-5)
        IF (j==3.OR.j==4) scale=density*1.5d0*phys%Mref*temperature
        trial=u;trial(j)=u(j)+perturbation*scale
        CALL held_value(trial,q,b,which,phi(which),plus)
        trial(j)=u(j)-perturbation*scale
        CALL held_value(trial,q,b,which,phi(which),minus)
        CALL close((plus-minus)/(2.d0*perturbation),ju(j)*scale,cap,label//' state tangent')
      ENDDO
      DO k=1,2
        qtrial=q;qtrial(k,energy)=q(k,energy)+perturbation*qscale
        CALL held_value(u,qtrial,b,which,phi(which),plus)
        qtrial(k,energy)=q(k,energy)-perturbation*qscale
        CALL held_value(u,qtrial,b,which,phi(which),minus)
        CALL close((plus-minus)/(2.d0*perturbation),jq(k,energy)*qscale,cap,label//' gradient tangent')
      ENDDO
    ENDDO
  END SUBROUTINE check_fixture

  SUBROUTINE coefficients(state,which,v,alpha,coefficient,qfs)
    REAL*8, INTENT(IN) :: state(5)
    INTEGER, INTENT(IN) :: which
    REAL*8, INTENT(OUT) :: v(5),alpha,coefficient,qfs
    REAL*8 :: energy,minimum_energy
    TYPE(parallel_heat_conduction_flux) :: heat_flux
    CALL heat_flux%init(state)
    CALL heat_flux%evaluate_conductivities()
    minimum_energy=1.d-5/simpar%refval_temperature*1.5d0*phys%Mref
    IF (which==1) THEN
      alpha=heat_flux%ion_heat_flux%conductivity%value
      CALL computeVi(state,v)
      coefficient=phys%diff_pari
      energy=MAX(state(3)/state(1)-.5d0*(state(2)/state(1))**2,minimum_energy)
      qfs=(2.d0/3.d0)**1.5d0*energy**1.5d0*state(1)
    ELSE
      alpha=heat_flux%electron_heat_flux%conductivity%value
      CALL computeVe(state,v)
      coefficient=phys%diff_pare
      energy=MAX(state(4)/state(1),minimum_energy)
      qfs=(2.d0/3.d0)**1.5d0*SQRT(ion_mass_si/electron_mass_si)*energy**1.5d0*state(1)
    ENDIF
    coefficient=coefficient*(2.d0/(3.d0*phys%Mref))**(1.d0+phys%epn)
  END SUBROUTINE coefficients

  SUBROUTINE held_value(state,gradient,b,which,phi,result)
    REAL*8, INTENT(IN) :: state(5),gradient(2,5),b(3),phi
    INTEGER, INTENT(IN) :: which
    REAL*8, INTENT(OUT) :: result
    LOGICAL :: enabled
    TYPE(parallel_heat_conduction_flux) :: heat_flux
    enabled=switch%flux_limiter;switch%flux_limiter=.FALSE.
    CALL heat_flux%init(state,Q=gradient,b=b)
    CALL heat_flux%evaluate_fluxes()
    result=phi*heat_flux%ion_heat_flux%diffusive_flux
    IF (which==2) result=phi*heat_flux%electron_heat_flux%diffusive_flux
    switch%flux_limiter=enabled
  END SUBROUTINE held_value

  SUBROUTINE close(actual,expected,normalization,label)
    REAL*8, INTENT(IN) :: actual,expected,normalization
    CHARACTER(*), INTENT(IN) :: label
    REAL*8 :: error
    error=ABS(actual-expected)/MAX(normalization,ABS(actual),ABS(expected))
    IF (.NOT.ieee_is_finite(error).OR.error>3.d-6) THEN
      WRITE(*,*) label,actual,expected,error
      ERROR STOP 'heat limiter consistency'
    ENDIF
    checks=checks+1;worst=MAX(worst,error)
  END SUBROUTINE close
END PROGRAM test_heat_flux_limiter
