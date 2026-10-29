#ifndef NEUTRALGAMMA
#error This check requires NEUTRALGAMMA.
#endif
PROGRAM test_neutral_gamma_transport
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys, switch
  USE physics, ONLY: set_model_layout, jacobianMatrices, jacobianMatricesFace, &
       GimpMatrixN, computeVun, compute_dVun_dU, phys2cons, cons2phys
  IMPLICIT NONE
  REAL*8, ALLOCATABLE :: U(:), up(:), um(:), v(:), vp(:), vm(:), dv(:,:), a(:,:), an(:,:), g(:,:), q(:,:)
  REAL*8, ALLOCATABLE :: physical(:,:), conservative(:,:)
  REAL*8 :: scales(7), densities(3), flux, pressure, step, denominator, numerical
  REAL*8 :: worst, bn, eta, divb
  REAL*8, PARAMETER :: perturbation=1.d-7, threshold=3.d-6
  INTEGER :: inn, ign, neq, density_case, direction, j, k, row, checks, states

  CALL set_model_layout()
  phys%Mref=1.3d0
  inn=phys%idx_rhon_eq
  ign=phys%idx_gamman_eq
  neq=phys%Neq
  ALLOCATE(U(neq),up(neq),um(neq),v(neq),vp(neq),vm(neq),dv(neq,neq), &
       a(neq,neq),an(neq,neq),g(neq,neq),q(2,neq))
  ALLOCATE(physical(1,phys%npv),conservative(1,neq))
  densities=[1.d-8,1.d-7,1.d-3]
  bn=-.6d0
  eta=2.1d0
  divb=.37d0
  checks=0
  states=0
  worst=0.d0
  DO density_case=1,SIZE(densities)
    DO direction=-1,1
      U=0.d0
      U(1:4)=[2.d0,.4d0,3.d0,4.d0]
      U(inn)=densities(density_case)
      U(ign)=DBLE(direction)*.3d0*U(inn)
      IF (phys%idx_k_eq>0) U(phys%idx_k_eq)=.01d0
      scales=1.d0
      scales(1:4)=U(1:4)
      scales(inn)=U(inn)
      scales(ign)=U(inn)
      CALL jacobianMatrices(U,a)
      CALL jacobianMatricesFace(U,bn,an)
      CALL GimpMatrixN(U,divb,g)
      CALL computeVun(U,v)
      CALL compute_dVun_dU(U,dv)
      flux=momentum_flux(U)
      pressure=neutral_pressure(U)
      CALL check(DOT_PRODUCT(a(ign,:),U),flux,ABS(flux))
      CALL check(DOT_PRODUCT(an(ign,:),U),bn*flux,ABS(bn*flux))
      CALL check(DOT_PRODUCT(g(ign,:),U),divb*pressure,ABS(divb*pressure))
      CALL check(a(inn,ign),1.d0,1.d0)
      DO j=1,neq
        step=perturbation*scales(j)
        up=U
        um=U
        up(j)=up(j)+step
        um(j)=um(j)-step
        denominator=2.d0*perturbation
        numerical=(momentum_flux(up)-momentum_flux(um))/denominator
        CALL check(numerical,a(ign,j)*scales(j),ABS(flux))
        CALL check(bn*numerical,an(ign,j)*scales(j),ABS(bn*flux))
        numerical=divb*(neutral_pressure(up)-neutral_pressure(um))/denominator
        CALL check(numerical,g(ign,j)*scales(j),ABS(divb*pressure))
        numerical=(velocity(up)-velocity(um))/denominator
        CALL check(numerical,v(j)*scales(j),1.d0)
        CALL computeVun(up,vp)
        CALL computeVun(um,vm)
        DO row=inn,ign
          CALL check((vp(row)-vm(row))/denominator,dv(row,j)*scales(j),1.d0/U(inn))
        ENDDO
      ENDDO
      q=0.d0
      q(:,inn)=[.1d0,.2d0]*U(inn)
      q(:,ign)=[.04d0,-.03d0]*U(inn)
      DO k=1,2
        up=U+perturbation*q(k,:)
        um=U-perturbation*q(k,:)
        denominator=2.d0*perturbation
        numerical=eta*(velocity(up)-velocity(um))/denominator
        CALL check(eta*DOT_PRODUCT(q(k,:),v),numerical,eta)
      ENDDO
      CALL cons2phys(RESHAPE(U,[1,neq]),physical)
      CALL check(physical(1,phys%idx_un_pv),velocity(U),1.d0)
      CALL phys2cons(physical,conservative)
      CALL check(conservative(1,ign),U(ign),U(inn))
      states=states+1
    ENDDO
  ENDDO
  WRITE(*,*) 'NeutralGamma transport: PASS; states, checks, worst error:',states,checks,worst

CONTAINS

  FUNCTION velocity(state) RESULT(value)
    REAL*8, INTENT(IN) :: state(:)
    REAL*8 :: value
    value=state(ign)/state(inn)
  END FUNCTION velocity

  FUNCTION neutral_pressure(state) RESULT(value)
    REAL*8, INTENT(IN) :: state(:)
    REAL*8 :: value
    value=state(inn)*2.d0/3.d0*(state(3)/state(1)-.5d0*(state(2)/state(1))**2)
  END FUNCTION neutral_pressure

  FUNCTION momentum_flux(state) RESULT(value)
    REAL*8, INTENT(IN) :: state(:)
    REAL*8 :: value
    value=state(ign)*velocity(state)+neutral_pressure(state)
  END FUNCTION momentum_flux

  SUBROUTINE check(actual,expected,scale)
    REAL*8, INTENT(IN) :: actual,expected,scale
    REAL*8 :: error
    error=ABS(actual-expected)/MAX(scale,ABS(actual),ABS(expected),1.d-30)
    IF (.NOT. ieee_is_finite(error) .OR. error>threshold) THEN
      WRITE(*,*) 'density case, momentum direction:',density_case,direction
      WRITE(*,*) 'actual, expected, scaled error:',actual,expected,error
      ERROR STOP 'neutral gamma transport'
    ENDIF
    worst=MAX(worst,error)
    checks=checks+1
  END SUBROUTINE check

END PROGRAM test_neutral_gamma_transport
