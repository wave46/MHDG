#ifdef TEMPERATURE
! Model coefficients only; HDG owns jumps, matrix placement and boundary laws.
SUBMODULE (physics) physics_stabilization
  IMPLICIT NONE
CONTAINS

  MODULE PROCEDURE evaluate_stab5
    REAL*8 :: advective,bn,bnorm,Etan,heat_part
    REAL*8 :: dadvective(SIZE(uc)),dheat(SIZE(uc)),dneutral(SIZE(uc))
    INTEGER :: i
    LOGICAL :: differentiate

    differentiate = .FALSE.
    IF (PRESENT(dtau_du)) dtau_du = 0.d0
    ! The optional tangent is currently supported for the 2D five-equation model.
#if defined(NEUTRAL) && !defined(NEUTRALGAMMA) && !defined(KEQUATION) && !defined(TOR3D)
    differentiate = PRESENT(dtau_du)
#endif
    bn = DOT_PRODUCT(b(1:SIZE(normal)),normal)
    bnorm = NORM2(b(1:SIZE(normal)))
    CALL advective_stabilization(up,uc,bn,advective,dadvective,differentiate)
    tau = 0.d0
    DO i = 1,SIZE(uc)
      tau(i,i) = advective+diff_iso(i,i)*length_factor
      IF (differentiate) dtau_du(i,:) = dadvective
    ENDDO
    CALL thermal_stabilization(3,phys%T_fluxlim_maxi,phys%diff_pari,heat_part,dheat)
    tau(3,3) = tau(3,3)+heat_part
    IF (differentiate) dtau_du(3,:) = dtau_du(3,:)+dheat
    CALL thermal_stabilization(4,phys%T_fluxlim_maxe,phys%diff_pare,heat_part,dheat)
    tau(4,4) = tau(4,4)+heat_part
    IF (differentiate) dtau_du(4,:) = dtau_du(4,:)+dheat
#ifdef NEUTRAL
    IF (differentiate) THEN
      ! setLocalDiff supplies unlimited Dnn here; fixed additions have zero derivative.
      CALL compute_Dnn_dU(uc,dneutral)
      dtau_du(phys%idx_rhon_eq,:) = dtau_du(phys%idx_rhon_eq,:)+length_factor*dneutral
    ENDIF
#endif
#ifdef NEUTRALGAMMA
    IF (phys%idx_gamman_eq > 0 .AND. phys%idx_rhon_eq > 0) THEN
      CALL computeEtan(uc,Etan)
      tau(phys%idx_gamman_eq,phys%idx_gamman_eq) = &
           advective+Etan/uc(phys%idx_rhon_eq)*length_factor
    ENDIF
#endif

  CONTAINS

    SUBROUTINE thermal_stabilization(row,temperature_cap,parallel_diffusion,thermal_tau,derivative)
      INTEGER, INTENT(IN) :: row
      REAL*8, INTENT(IN) :: temperature_cap,parallel_diffusion
      REAL*8, INTENT(OUT) :: thermal_tau,derivative(:)
      REAL*8 :: temperature,geometry,dtemperature(SIZE(uc))

      temperature = MIN(temperature_cap,up(row+4))
      geometry = ABS(bn)*parallel_diffusion*bnorm*length_factor
      thermal_tau = ABS(bn)*parallel_diffusion*temperature**2.5d0*bnorm/uc(1)*length_factor
      derivative = 0.d0
      IF (.NOT.differentiate) RETURN
      CALL capped_temperature_derivative(uc,row,up(row+4),temperature_cap,dtemperature)
      derivative = geometry*2.5d0*temperature**1.5d0/uc(1)*dtemperature
      derivative(1) = derivative(1)-thermal_tau/uc(1)
    END SUBROUTINE thermal_stabilization

  END PROCEDURE evaluate_stab5

  SUBROUTINE advective_stabilization(up,uc,bn,advective_tau,derivative,differentiate)
    REAL*8, INTENT(IN) :: up(:),uc(:),bn
    REAL*8, INTENT(OUT) :: advective_tau,derivative(:)
    LOGICAL, INTENT(IN) :: differentiate
    REAL*8 :: first,second,chi,root,npositive,dchi(SIZE(uc)),dvelocity(SIZE(uc))

    ! Coupled plasma convection: max(|5*u/3|, |u|+c_flux), projected on the normal.
    first = 5.d0/3.d0*up(2)*bn
    chi = 10.d0*uc(3)*uc(1)+10.d0*uc(4)*uc(1)-5.d0*uc(2)**2
    root = SQRT(ABS(chi))
    second = bn*(ABS(up(2))+root/(3.d0*uc(1)))
    advective_tau = MAX(ABS(first),ABS(second))
    derivative = 0.d0
    IF (.NOT.differentiate) RETURN
    IF (ABS(first) >= ABS(second)) THEN
      npositive = MAX(uc(1),1.d-20)
      IF (first == 0.d0) RETURN
      derivative(2) = 5.d0/3.d0*bn/npositive*SIGN(1.d0,first)
      IF (uc(1) > 1.d-20) derivative(1) = -derivative(2)*uc(2)/npositive
    ELSEIF (second /= 0.d0) THEN
      npositive = MAX(uc(1),1.d-20)
      dvelocity = 0.d0
      dvelocity(2) = 1.d0/npositive
      IF (uc(1) > 1.d-20) dvelocity(1) = -uc(2)/npositive**2
      IF (up(2) /= 0.d0) derivative = bn*SIGN(1.d0,up(2))*dvelocity
      IF (chi /= 0.d0) THEN
        dchi = 0.d0
        dchi(1:4) = [10.d0*(uc(3)+uc(4)),-10.d0*uc(2),10.d0*uc(1),10.d0*uc(1)]
        derivative = derivative+bn/(3.d0*uc(1))*0.5d0*SIGN(1.d0,chi)/root*dchi
      ENDIF
      derivative(1) = derivative(1)-bn*root/(3.d0*uc(1)**2)
      derivative = SIGN(1.d0,second)*derivative
    ENDIF
    ! At chi=0 the square-root branch has no finite derivative; freeze that part.
  END SUBROUTINE advective_stabilization

  SUBROUTINE capped_temperature_derivative(uc,row,temperature,cap,derivative)
    REAL*8, INTENT(IN) :: uc(:),temperature,cap
    INTEGER, INTENT(IN) :: row
    REAL*8, INTENT(OUT) :: derivative(:)
    REAL*8 :: npositive,velocity,scale

    derivative = 0.d0
    npositive = MAX(uc(1),1.d-20)
    velocity = uc(2)/npositive
    scale = 2.d0/(3.d0*phys%Mref)
    IF (temperature <= 1.d-20 .OR. temperature > cap) RETURN
    derivative(row) = scale/npositive
    IF (uc(1) > 1.d-20) derivative(1) = -scale*uc(row)/npositive**2
    IF (row == 3) THEN
      derivative(2) = -scale*velocity/npositive
      IF (uc(1) > 1.d-20) derivative(1) = derivative(1)+scale*velocity**2/npositive
    ENDIF
  END SUBROUTINE capped_temperature_derivative

END SUBMODULE physics_stabilization
#endif
