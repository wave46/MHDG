MODULE bohm_boundary
  USE smooth_transition, ONLY: smoothstep01, smoothstep01_derivative
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: evaluate_bohm_blend, evaluate_isothermal_bohm_target

  REAL*8, PARAMETER :: mach_width=0.20d0, energy_floor=1.d-20

CONTAINS

  PURE SUBROUTINE evaluate_bohm_blend(uhat,cs,dcs,bn,angle_threshold,energy_threshold,p0,ps,pf,dp0,dps,dpf)
    REAL*8, INTENT(IN) :: uhat(:),cs,dcs(:),bn,angle_threshold,energy_threshold
    REAL*8, INTENT(OUT) :: p0,ps,pf,dp0(:),dps(:),dpf(:)
    REAL*8 :: beta,wa,wm,dbeta(SIZE(uhat)),dwm(SIZE(uhat))

    dbeta=0.d0
    CALL energy_suppression(uhat(3),energy_threshold,beta,dbeta(3))
    wa=angular_target(bn,angle_threshold)
    CALL mach_activation(uhat,cs,dcs,bn,beta,dbeta,wm,dwm)
    p0=wm*(1.d0-beta*wa)
    ps=wm*beta*wa
    pf=1.d0-wm
    dp0=(1.d0-beta*wa)*dwm-wm*wa*dbeta
    dps=beta*wa*dwm+wm*wa*dbeta
    dpf=-dwm
  END SUBROUTINE evaluate_bohm_blend

  PURE SUBROUTINE evaluate_isothermal_bohm_target(uhat,velocity,cs,bn,angle_threshold,energy_threshold, &
       delta,free_weight,target_coefficient)
    REAL*8, INTENT(IN) :: uhat(:),velocity,cs,bn,angle_threshold,energy_threshold
    REAL*8, INTENT(OUT) :: delta,free_weight,target_coefficient
    REAL*8 :: beta,sn

    CALL energy_suppression(uhat(3),energy_threshold,beta)
    sn=SIGN(1.d0,bn)
    delta=1.d0
    free_weight=0.d0
    target_coefficient=0.d0
    IF(ABS(bn)<=angle_threshold)THEN
      IF(angle_threshold>0.d0)target_coefficient=sn*cs*beta*ABS(bn)/angle_threshold
      RETURN
    ENDIF
    IF(ABS(velocity)<=cs/beta)THEN
      target_coefficient=sn*cs*beta
      RETURN
    ENDIF
    delta=0.d0
    free_weight=1.d0
  END SUBROUTINE evaluate_isothermal_bohm_target

  PURE SUBROUTINE energy_suppression(energy_density,threshold,beta,derivative)
    REAL*8, INTENT(IN) :: energy_density,threshold
    REAL*8, INTENT(OUT) :: beta
    REAL*8, OPTIONAL, INTENT(OUT) :: derivative
    REAL*8 :: energy

    energy=MAX(energy_density,energy_floor)
    beta=energy/(energy+threshold)
    IF(.NOT.PRESENT(derivative))RETURN
    derivative=0.d0
    IF(energy_density>energy_floor)derivative=beta*(1.d0-beta)/energy
  END SUBROUTINE energy_suppression

  PURE FUNCTION angular_target(bn,threshold) RESULT(weight)
    REAL*8, INTENT(IN) :: bn,threshold
    REAL*8 :: weight

    weight=0.d0
    IF(bn==0.d0)RETURN
    weight=1.d0
    IF(threshold>0.d0)weight=smoothstep01(ABS(bn)/threshold)
  END FUNCTION angular_target

  PURE SUBROUTINE mach_activation(uhat,cs,dcs,bn,beta,dbeta,weight,derivative)
    REAL*8, INTENT(IN) :: uhat(:),cs,dcs(:),bn,beta,dbeta(:)
    REAL*8, INTENT(OUT) :: weight,derivative(:)
    REAL*8 :: mach,z,sn,dmach(SIZE(uhat))

    sn=MERGE(1.d0,-1.d0,bn>=0.d0)
    mach=sn*uhat(2)/(uhat(1)*cs)
    weight=1.d0
    derivative=0.d0
    IF(beta*mach<=1.d0)RETURN
    weight=0.d0
    IF(beta*(mach-mach_width)>=1.d0)RETURN

    ! Test compact support before dividing by a possibly tiny beta.
    z=(mach-1.d0/beta)/mach_width
    dmach=-mach*dcs/cs
    dmach(1)=dmach(1)-mach/uhat(1)
    dmach(2)=dmach(2)+sn/(uhat(1)*cs)
    weight=1.d0-smoothstep01(z)
    derivative=-smoothstep01_derivative(z)/mach_width*(dmach+(dbeta/beta)/beta)
  END SUBROUTINE mach_activation

END MODULE bohm_boundary
