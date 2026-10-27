PROGRAM test_stabilization_values
  USE globals, ONLY: phys,simpar,switch,numer,Mesh,refElPol
  USE physics, ONLY: set_model_layout,cons2phys,computeTauGaussPoints
  IMPLICIT NONE
  REAL*8 :: u(1,5),q(10),b(3),n(2),diffusion(5,5),tau(5,5),expected(5,5),original(5,5)
  REAL*8, ALLOCATABLE :: up(:,:)
  REAL*8 :: advective,bn,length_factor,geometry
  INTEGER :: i,checks

  CALL set_model_layout()
  ALLOCATE(up(1,phys%npv))
  simpar%Ndim=2;simpar%Neq=5
  phys%Mref=1.3d0;phys%diff_pari=.4d0;phys%diff_pare=.7d0
  phys%T_fluxlim_maxi=.6d0;phys%T_fluxlim_maxe=.9d0
  phys%c_fli=.2d0;phys%c_fle=.3d0
  switch%testcase=50;numer%stab=5
  ALLOCATE(Mesh%elemSize(1));Mesh%elemSize=.7d0;refElPol%ndeg=3
  b=[.2d0,.4d0,SQRT(.8d0)];n=[.6d0,.8d0]
  diffusion=0.d0
  DO i=1,5
    diffusion(i,i)=.1d0*i
  ENDDO
  checks=0
  CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
  CALL check_state([1.d0,3.d0,6.5d0,.4d0,.2d0])
  CALL check_state([1.d-8,0.d0,2.d-8,1.d-8,3.d-9])
  CALL check_state([1.d0,0.d0,-1.d-6,1.d0,.1d0])
  n=[-.6d0,-.8d0]
  CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
  n=[2.d0,-1.d0]/SQRT(5.d0)
  CALL check_state([.8d0,-.15d0,1.6d0,2.2d0,.1d0])
  WRITE(*,*) 'Stabilization values: PASS; coefficient/limiter/gradient checks:',checks

CONTAINS

  SUBROUTINE check_state(state)
    REAL*8, INTENT(IN) :: state(5)
    u(1,:)=state
    CALL cons2phys(u,up)
    bn=DOT_PRODUCT(b(1:2),n)
    length_factor=refElPol%ndeg/Mesh%elemSize(1)
    advective=MAX(ABS(5.d0/3.d0*state(2)/state(1)*bn), &
         ABS(.3d0*bn*(3.d0*state(1)+SQRT(ABS(10.d0*state(1)*(state(3)+state(4))-5.d0*state(2)**2)))/state(1)))
    expected=0.d0
    DO i=1,5
      expected(i,i)=advective+diffusion(i,i)*length_factor
    ENDDO
    geometry=ABS(bn)*NORM2(b(1:2))*length_factor/state(1)
    expected(3,3)=expected(3,3)+geometry*phys%diff_pari*MIN(up(1,7),phys%T_fluxlim_maxi)**2.5d0
    expected(4,4)=expected(4,4)+geometry*phys%diff_pare*MIN(up(1,8),phys%T_fluxlim_maxe)**2.5d0
    q=0.d0;switch%flux_limiter=.FALSE.
    CALL computeTauGaussPoints(up(1,:),state,q,b,n,1,0.,[0.d0,0.d0],tau,diffusion)
    CALL close(tau,expected)
    original=tau
    q=[(DBLE(i),i=1,10)]*1.d5;switch%flux_limiter=.TRUE.
    CALL computeTauGaussPoints(up(1,:),state,q,b,n,1,1.,[0.d0,0.d0],tau,diffusion)
    IF (ANY(tau/=original)) ERROR STOP 'tau depends on limiter or gradients'
    checks=checks+1
  END SUBROUTINE check_state

  SUBROUTINE close(actual,reference)
    REAL*8, INTENT(IN) :: actual(:,:),reference(:,:)
    IF (MAXVAL(ABS(actual-reference)/MAX(1.d0,ABS(reference)))>2.d-13) &
         ERROR STOP 'stabilization coefficient value'
    checks=checks+1
  END SUBROUTINE close
END PROGRAM test_stabilization_values
