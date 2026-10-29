PROGRAM test_stabilization_values
  USE globals, ONLY: phys,simpar,switch,numer,Mesh,refElPol,refElTor,stab_automatic
  USE physics, ONLY: set_model_layout,cons2phys,computeTauGaussPoints,jacobianMatricesFace
  IMPLICIT NONE
#ifdef TOR3D
  INTEGER, PARAMETER :: ndim=3
#else
  INTEGER, PARAMETER :: ndim=2
#endif
  REAL*8 :: u(1,5),q(5*ndim),b(3),n(ndim),xy(2),diffusion(5,5)
  REAL*8 :: tau(5,5),expected(5,5),original(5,5)
  REAL*8, ALLOCATABLE :: up(:,:)
  REAL*8 :: advective,bn,length_factor,geometry
  INTEGER :: i,checks

  CALL set_model_layout()
  ALLOCATE(up(1,phys%npv))
  simpar%Ndim=ndim;simpar%Neq=5
  phys%lscale=1.901d-3
  phys%Mref=1.3d0;phys%diff_pari=.4d0;phys%diff_pare=.7d0
  phys%T_fluxlim_maxi=.6d0;phys%T_fluxlim_maxe=.9d0
  phys%c_fli=.2d0;phys%c_fle=.3d0
  switch%testcase=50;numer%stab=stab_automatic
  ALLOCATE(Mesh%elemSize(1));Mesh%elemSize=.7d0;refElPol%ndeg=3
  b=[.2d0,.4d0,SQRT(.8d0)];n=0.d0;n(1:2)=[.6d0,.8d0]
  xy=[2.1d0,-.3d0]
  length_factor=refElPol%ndeg/Mesh%elemSize(1)
  diffusion=0.d0
  DO i=1,5
    diffusion(i,i)=.1d0*i
  ENDDO
  checks=0
  CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
  CALL check_state([1.d0,3.d0,6.5d0,.4d0,.2d0])
  CALL check_state([1.d-8,0.d0,2.d-8,1.d-8,3.d-9])
  CALL check_state([1.d0,0.d0,-1.d-6,1.d0,.1d0])
  n(1:2)=[-.6d0,-.8d0]
  CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
  n(1:2)=[2.d0,-1.d0]/SQRT(5.d0)
  CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
#ifdef TOR3D
  CALL check_toroidal_faces()
#endif
  WRITE(*,*) 'Stabilization values: PASS; dimension, checks:',ndim,checks

CONTAINS

  SUBROUTINE check_state(state)
    REAL*8, INTENT(IN) :: state(5)
    u(1,:)=state
    CALL cons2phys(u,up)
    bn=DOT_PRODUCT(b(1:ndim),n)
    advective=ABS(bn)*MAX(5.d0/3.d0*ABS(state(2)/state(1)), &
         ABS(state(2)/state(1))+SQRT(ABS(10.d0/9.d0*(state(3)+state(4))/state(1) &
         -5.d0/9.d0*(state(2)/state(1))**2)))
    expected=0.d0
    DO i=1,5
      expected(i,i)=advective+diffusion(i,i)*length_factor
    ENDDO
    geometry=ABS(bn)*NORM2(b(1:ndim))*length_factor/state(1)
    expected(3,3)=expected(3,3)+geometry*phys%diff_pari*MIN(up(1,7),phys%T_fluxlim_maxi)**2.5d0
    expected(4,4)=expected(4,4)+geometry*phys%diff_pare*MIN(up(1,8),phys%T_fluxlim_maxe)**2.5d0
    q=0.d0;switch%flux_limiter=.FALSE.
    CALL computeTauGaussPoints(up(1,:),state,q,b,n,1,0.,xy,tau,diffusion)
    CALL close(tau,expected)
    CALL check_convective_spectrum(state,tau(1,1)-diffusion(1,1)*length_factor)
    original=tau
    q=[(DBLE(i),i=1,5*ndim)]*1.d5;switch%flux_limiter=.TRUE.
    CALL computeTauGaussPoints(up(1,:),state,q,b,n,1,1.,xy,tau,diffusion)
    IF (ANY(tau/=original)) ERROR STOP 'tau depends on limiter or gradients'
    checks=checks+1
  END SUBROUTINE check_state

#ifdef TOR3D
  SUBROUTINE check_toroidal_faces()
    refElTor%ndeg=2;numer%ntor=8;numer%tmax=4.d0
    n=[0.d0,0.d0,1.d0]
    length_factor=2.d0/(xy(1)*.5d0)
    CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
    n=-n
    CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
    ! All rows, including neutrals, scale with the toroidal cell length.
    numer%ntor=16;length_factor=2.d0/(xy(1)*.25d0)
    CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
    xy(1)=2.d0*xy(1);length_factor=2.d0/(xy(1)*.25d0)
    CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
  END SUBROUTINE check_toroidal_faces
#endif

  SUBROUTINE check_convective_spectrum(state,advective_tau)
    REAL*8, INTENT(IN) :: state(5),advective_tau
    REAL*8 :: normal_jacobian(5,5),plasma(4,4),identity(4,4),polynomial(4,4)
    REAL*8 :: velocity,ion_energy,electron_energy,wave_speed,waves(4),scale
    INTEGER :: k

    velocity=state(2)/state(1)
    ion_energy=state(3)/state(1)-.5d0*velocity**2
    electron_energy=state(4)/state(1)
    IF (MIN(ion_energy,electron_energy)<=0.d0) RETURN
    wave_speed=SQRT(5.d0/3.d0*phys%Mref*(up(1,7)+up(1,8)))
    waves=bn*[velocity-wave_speed,velocity,velocity+wave_speed,5.d0/3.d0*velocity]
    IF (ABS(advective_tau-MAXVAL(ABS(waves)))>2.d-13*MAX(1.d0,ABS(advective_tau))) &
         ERROR STOP 'convective tau does not match characteristic speeds'

    CALL jacobianMatricesFace(state,bn,normal_jacobian)
    scale=MAX(1.d0,MAXVAL(ABS(normal_jacobian(1:4,1:4))),MAXVAL(ABS(waves)))
    plasma=normal_jacobian(1:4,1:4)/scale
    waves=waves/scale
    identity=0.d0
    DO k=1,4
      identity(k,k)=1.d0
    ENDDO
    ! The four proposed eigenvalues must annihilate the production flux block.
    polynomial=identity
    DO k=1,4
      polynomial=MATMUL(polynomial,plasma-waves(k)*identity)
    ENDDO
    IF (MAXVAL(ABS(polynomial))>2.d-12) ERROR STOP 'convection characteristic polynomial'
    checks=checks+2
  END SUBROUTINE check_convective_spectrum

  SUBROUTINE close(actual,reference)
    REAL*8, INTENT(IN) :: actual(:,:),reference(:,:)
    IF (MAXVAL(ABS(actual-reference)/MAX(1.d0,ABS(reference)))>2.d-13) &
         ERROR STOP 'stabilization coefficient value'
    checks=checks+1
  END SUBROUTINE close
END PROGRAM test_stabilization_values
